/*
 * Test 180: cancellation and progress reports (qaws_progress.h)
 *
 * For each long call (curve and surface batch closest points, prepared sets
 * and their queries, surface seams and shatter, 2D seams and shatter):
 *   - with a progress that never cancels, the same result as without one;
 *     reports start at stage 0, each stage's done never decreases, the last
 *     report of a stage has done equal to total;
 *   - cancelled at its k-th poll, for every k until it completes: the call
 *     returns QAWS_STATUS_CANCELLED and leaves no object behind (the
 *     runner's allocation checks see the rest).
 * The executor case: the batches through a serial executor, cancelled
 * between its waves the same way.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

typedef struct pg_probe
{
	unsigned int polls, cancel_at;   /* cancel_at 0: never */
	unsigned int reports, bad, stages_seen;
	unsigned int last_stage;
	unsigned long long last_done, last_total;
	int ended_ok;                    /* the last report had done == total */
} pg_probe;

static int pg_cancel(void* user)
{
	pg_probe* p = (pg_probe*)user;
	p->polls++;
	return p->cancel_at && p->polls >= p->cancel_at;
}

static void pg_report(qaws_progress_report const* r, void* user)
{
	pg_probe* p = (pg_probe*)user;
	if (!r->stage || r->stage_index >= r->stage_count || r->done > r->total)
		p->bad++;
	if (p->reports && r->stage_index == p->last_stage && r->done < p->last_done)
		p->bad++;
	if (p->reports && r->stage_index < p->last_stage)
		p->bad++;
	if (!p->reports || r->stage_index != p->last_stage)
		p->stages_seen++;
	p->last_stage = r->stage_index;
	p->last_done = r->done;
	p->last_total = r->total;
	p->ended_ok = r->done == r->total;
	p->reports++;
}

static qaws_progress pg_make(pg_probe* p, unsigned int cancel_at)
{
	qaws_progress g;
	memset(p, 0, sizeof *p);
	p->cancel_at = cancel_at;
	memset(&g, 0, sizeof g);
	g.cancel = pg_cancel;
	g.report = pg_report;
	g.user_data = p;
	return g;
}

/* a serial executor: the chunks in order on this thread */
static void pg_serial_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	(void)user;
	task(ctx, 0, count);
}

/* the call under test, run with progress g; returns its status (objects it made destroyed) */
typedef qaws_status (*pg_call)(void* ctx, qaws_progress const* g);

/* every k until the call completes: CANCELLED; then the clean run's reports checked */
static void pg_sweep(char const* name, pg_call call, void* ctx)
{
	pg_probe p;
	qaws_progress g;
	unsigned int k, cancelled = 0, wrong = 0, polls;
	char msg[256];
	qaws_status st;
	g = pg_make(&p, 0);
	st = call(ctx, &g);
	polls = p.polls;
	sprintf(msg, "%s: completes with a progress (status %d, %u polls, %u reports, %u stages, %u bad, ends %d)", name, (int)st, polls, p.reports,
		p.stages_seen, p.bad, p.ended_ok);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && polls > 0 && p.reports > 0 && !p.bad && p.ended_ok, msg);
	for (k = 1; k <= polls; k++)
	{
		g = pg_make(&p, k);
		st = call(ctx, &g);
		if (st == QAWS_STATUS_CANCELLED)
			cancelled++;
		else
			wrong++;
	}
	sprintf(msg, "%s: cancelled at each of its %u polls (%u not)", name, polls, wrong);
	TEST_ASSERT(cancelled == polls && !wrong, msg);
}

/* ---- curve batch ---- */

typedef struct pg_curves
{
	qaws_curve* c[3];
	qaws_scalar* pts;
	unsigned int npts;
	qaws_closest_point* ref;
	qaws_closest_point* out;
	qaws_batch_executor const* ex;
} pg_curves;

static qaws_status pg_curve_closest(void* ctx, qaws_progress const* g)
{
	pg_curves* C = (pg_curves*)ctx;
	qaws_closest_desc d;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.curves = (qaws_curve const* const*)C->c;
	d.curve_count = 3;
	d.points = C->pts;
	d.point_count = C->npts;
	d.executor = C->ex;
	d.progress = g;
	st = qaws_curve_batch_find_closest(C->out, NULL, &d);
	if (st == QAWS_STATUS_OK)
	{
		unsigned int i;
		for (i = 0; i < C->npts; i++)
			if (C->out[i].curve != C->ref[i].curve || C->out[i].parameter != C->ref[i].parameter || C->out[i].distance != C->ref[i].distance)
				st = QAWS_STATUS_INTERNAL_ERROR;   /* (not the result without progress) */
	}
	return st;
}

static qaws_status pg_curve_set(void* ctx, qaws_progress const* g)
{
	pg_curves* C = (pg_curves*)ctx;
	qaws_curve_batch_desc d;
	qaws_curve_set* set = (qaws_curve_set*)1;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.curves = (qaws_curve const* const*)C->c;
	d.curve_count = 3;
	d.progress = g;
	st = qaws_curve_set_create(&set, &d);
	if (st != QAWS_STATUS_OK)
		return set ? QAWS_STATUS_INTERNAL_ERROR : st;   /* (a failed create leaves no set) */
	qaws_curve_set_destroy(set);
	return st;
}

static qaws_status pg_curve_set_query(void* ctx, qaws_progress const* g)
{
	pg_curves* C = (pg_curves*)ctx;
	qaws_curve_batch_desc d;
	qaws_curve_set* set = NULL;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.curves = (qaws_curve const* const*)C->c;
	d.curve_count = 3;
	d.executor = C->ex;
	st = qaws_curve_set_create(&set, &d);
	if (st != QAWS_STATUS_OK)
		return QAWS_STATUS_INTERNAL_ERROR;
	st = qaws_curve_set_find_closest_ex(C->out, NULL, set, C->pts, C->npts, 0, g);
	qaws_curve_set_destroy(set);
	return st;
}

static void test_curve_batch(void)
{
	pg_curves C;
	qaws_closest_desc d;
	qaws_batch_executor ex;
	unsigned int i;
	memset(&C, 0, sizeof C);
	qaws_curve_create_ellipse_2d(&C.c[0], 0, 0, 1, 1, 0);
	qaws_curve_create_ellipse_2d(&C.c[1], 1.5, 0, 0.7, 0.4, 0.3);
	qaws_curve_create_ellipse_2d(&C.c[2], -1, 1, 0.5, 0.9, 1.1);
	C.npts = 3000;
	C.pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * C.npts);
	C.ref = (qaws_closest_point*)malloc(sizeof *C.ref * C.npts);
	C.out = (qaws_closest_point*)malloc(sizeof *C.out * C.npts);
	for (i = 0; i < C.npts; i++)
	{
		C.pts[2 * i] = (qaws_scalar)(-2.5 + 5.0 * fmod(i * 0.6180339887, 1.0));
		C.pts[2 * i + 1] = (qaws_scalar)(-2.5 + 5.0 * fmod(i * 0.4142135623, 1.0));
	}
	memset(&d, 0, sizeof d);
	d.curves = (qaws_curve const* const*)C.c;
	d.curve_count = 3;
	d.points = C.pts;
	d.point_count = C.npts;
	TEST_ASSERT_STATUS(qaws_curve_batch_find_closest(C.ref, NULL, &d));
	pg_sweep("curve batch closest", pg_curve_closest, &C);
	ex.parallel_for = pg_serial_for;
	ex.user = NULL;
	C.ex = &ex;
	pg_sweep("curve batch closest, executor", pg_curve_closest, &C);
	pg_sweep("curve set create", pg_curve_set, &C);
	pg_sweep("curve set query, executor", pg_curve_set_query, &C);
	for (i = 0; i < 3; i++)
		qaws_curve_destroy(C.c[i]);
	free(C.pts);
	free(C.ref);
	free(C.out);
}

/* ---- surface batch ---- */

typedef struct pg_surfaces
{
	qaws_surface* s[6];
	qaws_scalar* pts;
	unsigned int npts;
	qaws_surface_batch_closest* ref;
	qaws_surface_batch_closest* out;
	qaws_batch_executor const* ex;
	int query;                       /* pg_surface_set: the progress to the query, not the create */
	unsigned int seams;              /* pg_seams: the count without progress */
} pg_surfaces;

static qaws_vec3 pg_v(double x, double y, double z)
{
	qaws_vec3 v;
	v.x = (qaws_scalar)x;
	v.y = (qaws_scalar)y;
	v.z = (qaws_scalar)z;
	return v;
}

static qaws_surface* pg_bilinear(qaws_vec3 p00, qaws_vec3 p10, qaws_vec3 p01, qaws_vec3 p11)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	d.p00 = p00;
	d.p10 = p10;
	d.p01 = p01;
	d.p11 = p11;
	qaws_surface_create_bilinear(&s, &d);
	return s;
}

static qaws_status pg_surface_closest(void* ctx, qaws_progress const* g)
{
	pg_surfaces* S = (pg_surfaces*)ctx;
	qaws_surface_batch_closest_desc d;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.surfaces = (qaws_surface const* const*)S->s;
	d.surface_count = 6;
	d.points = S->pts;
	d.point_count = S->npts;
	d.executor = S->ex;
	d.progress = g;
	st = qaws_surface_batch_find_closest(S->out, NULL, &d);
	if (st == QAWS_STATUS_OK && memcmp(S->out, S->ref, sizeof *S->out * S->npts) != 0)
		st = QAWS_STATUS_INTERNAL_ERROR;
	return st;
}

static qaws_status pg_surface_set(void* ctx, qaws_progress const* g)
{
	pg_surfaces* S = (pg_surfaces*)ctx;
	qaws_surface_batch_desc d;
	qaws_surface_set* set = NULL;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.surfaces = (qaws_surface const* const*)S->s;
	d.surface_count = 6;
	d.executor = S->ex;
	d.progress = S->query ? NULL : g;
	st = qaws_surface_set_create(&set, &d);
	if (st != QAWS_STATUS_OK)
		return set ? QAWS_STATUS_INTERNAL_ERROR : st;
	st = qaws_surface_set_find_closest_ex(S->out, NULL, set, S->pts, S->npts, 0, S->query ? g : NULL);
	qaws_surface_set_destroy(set);
	if (st == QAWS_STATUS_OK && memcmp(S->out, S->ref, sizeof *S->out * S->npts) != 0)
		st = QAWS_STATUS_INTERNAL_ERROR;
	return st;
}

/* seams of a sheet and a cube, and its shatter */
static qaws_status pg_seams(void* ctx, qaws_progress const* g)
{
	pg_surfaces* S = (pg_surfaces*)ctx;
	qaws_surface* sheet = pg_bilinear(pg_v(-1, -1, 0.1), pg_v(1, -1, 0.1), pg_v(-1, 1, 0.1), pg_v(1, 1, 0.1));
	qaws_surface_seam_desc d;
	qaws_surface_seams* R = (qaws_surface_seams*)1;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.a = (qaws_surface const* const*)&sheet;
	d.a_count = 1;
	d.b = (qaws_surface const* const*)S->s;
	d.b_count = 6;
	d.kinds = QAWS_SEAM_A_B;
	d.resolution = 24;
	d.executor = S->ex;
	d.progress = g;
	st = qaws_surface_seams_create(&R, &d);
	if (st == QAWS_STATUS_OK)
	{
		if (!g)
			S->seams = qaws_surface_seams_get_count(R);
		else if (qaws_surface_seams_get_count(R) != S->seams)
			st = QAWS_STATUS_INTERNAL_ERROR;
		qaws_surface_seams_destroy(R);
	}
	else if (R)
		st = QAWS_STATUS_INTERNAL_ERROR;   /* (a failed create leaves none) */
	qaws_surface_destroy(sheet);
	return st;
}

static qaws_status pg_shatter(void* ctx, qaws_progress const* g)
{
	pg_surfaces* S = (pg_surfaces*)ctx;
	qaws_surface* sheet = pg_bilinear(pg_v(-1, -1, 0.1), pg_v(1, -1, 0.1), pg_v(-1, 1, 0.1), pg_v(1, 1, 0.1));
	qaws_surface_seam_desc d;
	qaws_surface_shatter* H = (qaws_surface_shatter*)1;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.a = (qaws_surface const* const*)&sheet;
	d.a_count = 1;
	d.b = (qaws_surface const* const*)S->s;
	d.b_count = 6;
	d.executor = S->ex;
	d.progress = g;
	st = qaws_surface_shatter_create(&H, &d, 1);
	if (st == QAWS_STATUS_OK)
	{
		if (qaws_surface_shatter_get_count(H) != 2)
			st = QAWS_STATUS_INTERNAL_ERROR;
		qaws_surface_shatter_destroy(H);
	}
	else if (H)
		st = QAWS_STATUS_INTERNAL_ERROR;
	qaws_surface_destroy(sheet);
	return st;
}

/* the seams of a sheet and the cube serially and through an executor: the same points, bit for bit */
static void pg_seams_same(pg_surfaces* S, qaws_batch_executor const* ex)
{
	qaws_surface* sheet = pg_bilinear(pg_v(-1, -1, 0.1), pg_v(1, -1, 0.1), pg_v(-1, 1, 0.1), pg_v(1, 1, 0.1));
	qaws_surface_seam_desc d;
	qaws_surface_seams *R0 = NULL, *R1 = NULL;
	unsigned int i, k, diff = 0, n = 0;
	char msg[160];
	memset(&d, 0, sizeof d);
	d.a = (qaws_surface const* const*)&sheet;
	d.a_count = 1;
	d.b = (qaws_surface const* const*)S->s;
	d.b_count = 6;
	d.kinds = QAWS_SEAM_A_B | QAWS_SEAM_B_B;
	TEST_ASSERT_STATUS(qaws_surface_seams_create(&R0, &d));
	d.executor = ex;
	TEST_ASSERT_STATUS(qaws_surface_seams_create(&R1, &d));
	if (qaws_surface_seams_get_count(R0) != qaws_surface_seams_get_count(R1))
		diff++;
	for (i = 0; !diff && i < qaws_surface_seams_get_count(R0); i++)
	{
		qaws_seam_info a, b;
		qaws_seam_point *pa, *pb;
		qaws_surface_seams_get(&a, NULL, R0, i);
		qaws_surface_seams_get(&b, NULL, R1, i);
		if (a.point_count != b.point_count || a.kind != b.kind || a.surface1 != b.surface1 || a.surface2 != b.surface2 || a.closed != b.closed)
		{
			diff++;
			break;
		}
		pa = (qaws_seam_point*)malloc(sizeof *pa * (a.point_count + 1));
		pb = (qaws_seam_point*)malloc(sizeof *pb * (b.point_count + 1));
		qaws_surface_seams_get(&a, pa, R0, i);
		qaws_surface_seams_get(&b, pb, R1, i);
		for (k = 0; k < a.point_count; k++, n++)
			if (pa[k].position.x != pb[k].position.x || pa[k].position.y != pb[k].position.y || pa[k].position.z != pb[k].position.z
				|| pa[k].u1 != pb[k].u1 || pa[k].v1 != pb[k].v1 || pa[k].u2 != pb[k].u2 || pa[k].v2 != pb[k].v2)
				diff++;
		free(pa);
		free(pb);
	}
	sprintf(msg, "surface seams through an executor: %u seams, %u points, %u differ", qaws_surface_seams_get_count(R0), n, diff);
	printf("  %s\n", msg);
	TEST_ASSERT(!diff && n > 0, msg);
	qaws_surface_seams_destroy(R0);
	qaws_surface_seams_destroy(R1);
	qaws_surface_destroy(sheet);
}

static void test_surface_batch(void)
{
	pg_surfaces S;
	qaws_surface_batch_closest_desc d;
	qaws_batch_executor ex;
	double x0 = -0.5, y0 = -0.5, z0 = -0.5, x1 = 0.5, y1 = 0.5, z1 = 0.5;
	unsigned int i;
	memset(&S, 0, sizeof S);
	S.s[0] = pg_bilinear(pg_v(x0, y0, z0), pg_v(x0, y1, z0), pg_v(x1, y0, z0), pg_v(x1, y1, z0));
	S.s[1] = pg_bilinear(pg_v(x0, y0, z1), pg_v(x1, y0, z1), pg_v(x0, y1, z1), pg_v(x1, y1, z1));
	S.s[2] = pg_bilinear(pg_v(x0, y0, z0), pg_v(x1, y0, z0), pg_v(x0, y0, z1), pg_v(x1, y0, z1));
	S.s[3] = pg_bilinear(pg_v(x0, y1, z0), pg_v(x0, y1, z1), pg_v(x1, y1, z0), pg_v(x1, y1, z1));
	S.s[4] = pg_bilinear(pg_v(x0, y0, z0), pg_v(x0, y0, z1), pg_v(x0, y1, z0), pg_v(x0, y1, z1));
	S.s[5] = pg_bilinear(pg_v(x1, y0, z0), pg_v(x1, y1, z0), pg_v(x1, y0, z1), pg_v(x1, y1, z1));
	S.npts = 2000;
	S.pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * S.npts);
	S.ref = (qaws_surface_batch_closest*)malloc(sizeof *S.ref * S.npts);
	S.out = (qaws_surface_batch_closest*)malloc(sizeof *S.out * S.npts);
	for (i = 0; i < S.npts; i++)
	{
		S.pts[3 * i] = (qaws_scalar)(-1.5 + 3.0 * fmod(i * 0.6180339887, 1.0));
		S.pts[3 * i + 1] = (qaws_scalar)(-1.5 + 3.0 * fmod(i * 0.4142135623, 1.0));
		S.pts[3 * i + 2] = (qaws_scalar)(-1.5 + 3.0 * fmod(i * 0.7320508075, 1.0));
	}
	memset(&d, 0, sizeof d);
	d.surfaces = (qaws_surface const* const*)S.s;
	d.surface_count = 6;
	d.points = S.pts;
	d.point_count = S.npts;
	TEST_ASSERT_STATUS(qaws_surface_batch_find_closest(S.ref, NULL, &d));
	pg_sweep("surface batch closest", pg_surface_closest, &S);
	ex.parallel_for = pg_serial_for;
	ex.user = NULL;
	S.ex = &ex;
	pg_sweep("surface batch closest, executor", pg_surface_closest, &S);
	pg_sweep("surface set create", pg_surface_set, &S);
	S.query = 1;
	pg_sweep("surface set query, executor", pg_surface_set, &S);
	S.ex = NULL;
	TEST_ASSERT(pg_seams(&S, NULL) == QAWS_STATUS_OK && S.seams > 0, "surface seams without progress");
	pg_sweep("surface seams", pg_seams, &S);
	pg_sweep("surface shatter", pg_shatter, &S);
	S.ex = &ex;   /* (the same seams and pieces through an executor) */
	pg_seams_same(&S, &ex);
	pg_sweep("surface seams, executor", pg_seams, &S);
	pg_sweep("surface shatter, executor", pg_shatter, &S);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(S.s[i]);
	free(S.pts);
	free(S.ref);
	free(S.out);
}

/* ---- 2D regions: seams, shatter, points located ---- */

typedef struct pg_regions
{
	qaws_curve* c[2];
	qaws_curve const* l0[1];
	qaws_curve const* l1[1];
	qaws_path_2d a, b;
	qaws_batch_executor const* ex;
	qaws_vec2* pts;                  /* the points the prepared region locates, ref the single calls' answers */
	unsigned int npts;
	qaws_point_location* ref;
	qaws_point_location* out;
} pg_regions;

static void pg_region_desc(pg_regions* R, qaws_region_seam_desc* d, qaws_progress const* g)
{
	memset(d, 0, sizeof *d);
	d->a = &R->a;
	d->a_count = 1;
	d->b = &R->b;
	d->b_count = 1;
	d->b_is_region = 1;
	d->kinds = QAWS_SEAM_A_B;
	d->executor = R->ex;
	d->progress = g;
}

static qaws_status pg_region_seams(void* ctx, qaws_progress const* g)
{
	qaws_region_seam_desc d;
	qaws_region_seam_point pts[8];
	unsigned int n = 0;
	qaws_status st;
	pg_region_desc((pg_regions*)ctx, &d, g);
	st = qaws_region_seams_2d(pts, &n, &d, 8);
	if (st == QAWS_STATUS_OK && n != 2)
		st = QAWS_STATUS_INTERNAL_ERROR;
	return st;
}

static qaws_status pg_region_shatter(void* ctx, qaws_progress const* g)
{
	qaws_region_seam_desc d;
	qaws_region_shatter* H = (qaws_region_shatter*)1;
	qaws_status st;
	pg_region_desc((pg_regions*)ctx, &d, g);
	st = qaws_region_shatter_create_2d(&H, &d);
	if (st == QAWS_STATUS_OK)
	{
		if (qaws_region_shatter_get_count(H) != 2)
			st = QAWS_STATUS_INTERNAL_ERROR;
		qaws_region_shatter_destroy(H);
	}
	else if (H)
		st = QAWS_STATUS_INTERNAL_ERROR;
	return st;
}

/* the two disks as one region (non-zero): every point located as qaws_region_locate_point_2d does */
static qaws_status pg_region_locate(void* ctx, qaws_progress const* g)
{
	pg_regions* R = (pg_regions*)ctx;
	qaws_path_2d both[2];
	qaws_region_2d* r = NULL;
	qaws_status st;
	both[0] = R->a, both[1] = R->b;
	st = qaws_region_create_2d(&r, both, 2, QAWS_FILL_NON_ZERO, 0);
	if (st != QAWS_STATUS_OK)
		return QAWS_STATUS_INTERNAL_ERROR;
	memset(R->out, 0xFF, sizeof *R->out * R->npts);
	st = qaws_region_locate_points_2d(R->out, r, g, R->pts, R->npts, R->ex);
	qaws_region_destroy_2d(r);
	if (st == QAWS_STATUS_OK && memcmp(R->out, R->ref, sizeof *R->out * R->npts) != 0)
		st = QAWS_STATUS_INTERNAL_ERROR;
	return st;
}

static void test_regions(void)
{
	pg_regions R;
	qaws_batch_executor ex;
	memset(&R, 0, sizeof R);
	qaws_curve_create_ellipse_2d(&R.c[0], 0, 0, 1, 1, 0);
	qaws_curve_create_ellipse_2d(&R.c[1], 1.2, 0, 1, 1, 0);
	R.l0[0] = R.c[0];
	R.l1[0] = R.c[1];
	R.a.curves = R.l0;
	R.a.curve_count = 1;
	R.a.closed = 1;
	R.b.curves = R.l1;
	R.b.curve_count = 1;
	R.b.closed = 1;
	{
		/* a 40 x 25 grid over both disks and the points (1, 0), (0.2, 0), (2.2, 0) on their circles */
		qaws_path_2d both[2];
		unsigned int i, on = 0, in = 0;
		R.npts = 1003;
		R.pts = (qaws_vec2*)malloc(sizeof *R.pts * R.npts);
		R.ref = (qaws_point_location*)malloc(sizeof *R.ref * R.npts);
		R.out = (qaws_point_location*)malloc(sizeof *R.out * R.npts);
		for (i = 0; i < 1000; i++)
			R.pts[i].x = (qaws_scalar)(-1.3 + 3.8 * (i % 40) / 39.0), R.pts[i].y = (qaws_scalar)(-1.2 + 2.4 * (i / 40) / 24.0);
		R.pts[1000].x = 1, R.pts[1000].y = 0;
		R.pts[1001].x = (qaws_scalar)0.2, R.pts[1001].y = 0;
		R.pts[1002].x = (qaws_scalar)2.2, R.pts[1002].y = 0;
		both[0] = R.a, both[1] = R.b;
		for (i = 0; i < R.npts; i++)
		{
			qaws_region_locate_point_2d(&R.ref[i], both, 2, QAWS_FILL_NON_ZERO, R.pts[i], 0);
			on += R.ref[i] == QAWS_POINT_ON, in += R.ref[i] == QAWS_POINT_INSIDE;
		}
		TEST_ASSERT(on >= 3 && in > 100 && in + on < R.npts, "the region's points: some on its boundary, some inside, some outside");
	}
	pg_sweep("2D seams", pg_region_seams, &R);
	pg_sweep("2D shatter", pg_region_shatter, &R);
	pg_sweep("2D region, points located", pg_region_locate, &R);
	ex.parallel_for = pg_serial_for;
	ex.user = NULL;
	R.ex = &ex;
	pg_sweep("2D seams, executor", pg_region_seams, &R);
	pg_sweep("2D shatter, executor", pg_region_shatter, &R);
	pg_sweep("2D region, points located, executor", pg_region_locate, &R);
	qaws_curve_destroy(R.c[0]);
	qaws_curve_destroy(R.c[1]);
	free(R.pts);
	free(R.ref);
	free(R.out);
}

int test_180_progress_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 180: cancellation and progress reports\n");
	test_curve_batch();
	test_surface_batch();
	test_regions();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
