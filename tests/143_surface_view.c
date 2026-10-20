/*
 * Test 143: borrowed surfaces (qaws_view.h)
 *
 *   - parity: prepared and borrowed surfaces against the owned B-spline and
 *     NURBS surfaces on the same grid (non-uniform knots with a double knot,
 *     clamped uniform default knots, weights and homogeneous points):
 *     position, du, dv, normal, duu, duv, dvv;
 *   - layouts: interleaved f32 records with a stride against packed f64 and
 *     planar f64 holding the same values;
 *   - lifetime: points rewritten in host memory between evaluations without
 *     preparing again;
 *   - range evaluation from four host threads over disjoint ranges, bit for
 *     bit the serial outputs;
 *   - allocation failures without leaks, prepare_into, validation;
 *   - operations on a borrowed surface (area, closest point) as on the owned.
 */

#include "test_common.h"
#include "test_threads.h"
#include "qaws_view.h"
#include "qaws_diff_ops.h"
#include "qaws_diff.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SV_U 6
#define SV_V 5
#define SV_UD 3
#define SV_VD 2
#define SV_N (SV_U * SV_V)

#if QAWS_SCALAR_IS_FLOAT
#define SV_TOL 1e-4
#else
#define SV_TOL 1e-11
#endif

static double const sv_uk[SV_U + SV_UD + 1] = { 0, 0, 0, 0, 0.4, 0.4, 1, 1, 1, 1 };
static double const sv_vk[SV_V + SV_VD + 1] = { 0, 0, 0, 0.3, 0.75, 1, 1, 1 };

static void sv_points(double* p, double phase)
{
	unsigned int i, j;
	for (i = 0; i < SV_U; i++)
		for (j = 0; j < SV_V; j++)
		{
			double* q = p + (i * SV_V + j) * 3;
			q[0] = i / (double)(SV_U - 1) + 0.05 * sin(1.7 * j + phase);
			q[1] = j / (double)(SV_V - 1) + 0.04 * cos(2.1 * i + phase);
			q[2] = 0.3 * sin(1.3 * i + phase) * cos(0.9 * j);
		}
}

static void sv_weights(double* w)
{
	unsigned int i;
	for (i = 0; i < SV_N; i++)
		w[i] = 0.6 + 0.5 * (double)((i * 7) % 5) / 4;
}

static qaws_array_view sv_packed(void const* data, qaws_element_type t, unsigned int comps, unsigned int n, size_t stride)
{
	qaws_array_view v;
	memset(&v, 0, sizeof v);
	v.data = data;
	v.type = t;
	v.components = comps;
	v.count = n;
	v.stride = stride;
	return v;
}

#define SV_ALL (QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DUU \
	| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV)

static double sv_diff(qaws_surface_eval_result const* a, qaws_surface_eval_result const* b)
{
	qaws_vec3 const* va[7] = { &a->position, &a->du, &a->dv, &a->normal, &a->duu, &a->duv, &a->dvv };
	qaws_vec3 const* vb[7] = { &b->position, &b->du, &b->dv, &b->normal, &b->duu, &b->duv, &b->dvv };
	double worst = 0;
	unsigned int r;
	for (r = 0; r < 7; r++)
		if (a->valid_flags & b->valid_flags & (1u << r))
		{
			double e = fabs((double)va[r]->x - vb[r]->x) + fabs((double)va[r]->y - vb[r]->y) + fabs((double)va[r]->z - vb[r]->z);
			double s = 1 + fabs((double)vb[r]->x) + fabs((double)vb[r]->y) + fabs((double)vb[r]->z);
			worst = e / s > worst ? e / s : worst;
		}
	return worst;
}

/* an owned B-spline (w NULL) or NURBS surface on the grid; knots NULL: the default ones */
static qaws_surface* sv_owned(double const* p, double const* w, int knots)
{
	qaws_vec3 cp[SV_N];
	qaws_scalar ws[SV_N], uk[SV_U + SV_UD + 1], vk[SV_V + SV_VD + 1];
	qaws_surface* s = NULL;
	unsigned int i;
	for (i = 0; i < SV_N; i++)
	{
		cp[i].x = (qaws_scalar)p[i * 3];
		cp[i].y = (qaws_scalar)p[i * 3 + 1];
		cp[i].z = (qaws_scalar)p[i * 3 + 2];
		ws[i] = w ? (qaws_scalar)w[i] : 1;
	}
	for (i = 0; i < SV_U + SV_UD + 1; i++)
		uk[i] = (qaws_scalar)sv_uk[i];
	for (i = 0; i < SV_V + SV_VD + 1; i++)
		vk[i] = (qaws_scalar)sv_vk[i];
	if (w)
	{
		qaws_surface_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = SV_UD;
		d.v_degree = SV_VD;
		d.control_points = cp;
		d.u_point_count = SV_U;
		d.v_point_count = SV_V;
		d.weights = ws;
		d.u_knots = knots ? uk : NULL;
		d.u_knot_count = knots ? SV_U + SV_UD + 1 : 0;
		d.v_knots = knots ? vk : NULL;
		d.v_knot_count = knots ? SV_V + SV_VD + 1 : 0;
		qaws_surface_create_nurbs(&s, &d);
	}
	else
	{
		qaws_surface_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = SV_UD;
		d.v_degree = SV_VD;
		d.control_points = cp;
		d.u_point_count = SV_U;
		d.v_point_count = SV_V;
		d.u_knots = knots ? uk : NULL;
		d.u_knot_count = knots ? SV_U + SV_UD + 1 : 0;
		d.v_knots = knots ? vk : NULL;
		d.v_knot_count = knots ? SV_V + SV_VD + 1 : 0;
		qaws_surface_create_bspline(&s, &d);
	}
	return s;
}

static qaws_surface_view sv_view(double const* p, unsigned int comps, double const* w, int knots)
{
	qaws_surface_view v;
	memset(&v, 0, sizeof v);
	v.u_degree = SV_UD;
	v.v_degree = SV_VD;
	v.u_count = SV_U;
	v.v_count = SV_V;
	v.points = sv_packed(p, QAWS_ELEMENT_F64, comps, SV_N, 0);
	if (w)
		v.weights = sv_packed(w, QAWS_ELEMENT_F64, 1, SV_N, 0);
	if (knots)
	{
		v.u_knots = sv_packed(sv_uk, QAWS_ELEMENT_F64, 1, SV_U + SV_UD + 1, 0);
		v.v_knots = sv_packed(sv_vk, QAWS_ELEMENT_F64, 1, SV_V + SV_VD + 1, 0);
	}
	return v;
}

/* worst difference over a 13 x 11 grid of the domain, prepared (and the borrowed surface when given) against owned */
static double sv_compare(qaws_surface_prepared const* p, qaws_surface const* borrowed, qaws_surface const* owned)
{
	qaws_range ru = qaws_surface_get_u_range(owned), rv = qaws_surface_get_v_range(owned);
	double worst = 0, e;
	unsigned int i, j;
	for (i = 0; i <= 12; i++)
		for (j = 0; j <= 10; j++)
		{
			qaws_scalar u = ru.min_value + (ru.max_value - ru.min_value) * (qaws_scalar)i / 12;
			qaws_scalar v = rv.min_value + (rv.max_value - rv.min_value) * (qaws_scalar)j / 10;
			qaws_surface_eval_result a, b;
			if (qaws_surface_evaluate(&b, owned, u, v, SV_ALL) != QAWS_STATUS_OK)
				return 1e30;
			if (qaws_surface_prepared_evaluate(&a, p, u, v, SV_ALL) != QAWS_STATUS_OK || (b.valid_flags & ~a.valid_flags))
				return 1e30;
			e = sv_diff(&a, &b);
			worst = e > worst ? e : worst;
			if (borrowed)
			{
				if (qaws_surface_evaluate(&a, borrowed, u, v, SV_ALL) != QAWS_STATUS_OK)
					return 1e30;
				e = sv_diff(&a, &b);
				worst = e > worst ? e : worst;
			}
		}
	return worst;
}

static void test_parity(void)
{
	double p[SV_N * 3], w[SV_N], h[SV_N * 4];
	double e[4];
	char const* name[4] = { "B-spline", "B-spline, default knots", "NURBS, weights", "NURBS, homogeneous points" };
	unsigned int c, i;
	char msg[300];
	sv_points(p, 0.4);
	sv_weights(w);
	for (i = 0; i < SV_N; i++)
	{
		h[i * 4 + 0] = p[i * 3 + 0] * w[i];
		h[i * 4 + 1] = p[i * 3 + 1] * w[i];
		h[i * 4 + 2] = p[i * 3 + 2] * w[i];
		h[i * 4 + 3] = w[i];
	}
	for (c = 0; c < 4; c++)
	{
		qaws_surface_view v = c == 3 ? sv_view(h, 4, NULL, 1) : sv_view(p, 3, c == 2 ? w : NULL, c != 1);
		qaws_surface* owned = sv_owned(p, c >= 2 ? w : NULL, c != 1);
		qaws_surface_prepared* pr = NULL;
		qaws_surface* borrowed = NULL;
		e[c] = 1e30;
		if (owned && qaws_surface_prepare(&pr, &v, NULL) == QAWS_STATUS_OK && qaws_surface_create_borrowed(&borrowed, &v, NULL) == QAWS_STATUS_OK)
		{
			e[c] = sv_compare(pr, borrowed, owned);
			if (c == 0)
				TEST_ASSERT(qaws_surface_get_kind(borrowed) == QAWS_SURFACE_KIND_VIEW && !qaws_surface_is_rational(borrowed)
					&& qaws_surface_prepared_dependencies(pr) == (QAWS_DEPENDS_KNOTS | QAWS_DEPENDS_LAYOUT),
					"borrowed surface: kind, polynomial, depends on knots and layout");
			if (c == 2)
				TEST_ASSERT(qaws_surface_is_rational(borrowed), "borrowed surface with weights: rational");
		}
		qaws_surface_prepared_destroy(pr);
		qaws_surface_destroy(borrowed);
		qaws_surface_destroy(owned);
	}
	sprintf(msg, "surface parity, position to second partials and normal: %s %.2g, %s %.2g, %s %.2g, %s %.2g", name[0], e[0], name[1], e[1],
		name[2], e[2], name[3], e[3]);
	TEST_ASSERT(e[0] <= SV_TOL && e[1] <= SV_TOL && e[2] <= SV_TOL && e[3] <= SV_TOL * 10, msg);
	printf("  %s\n", msg);
}

typedef struct sv_record
{
	float position[3];
	unsigned int flags;
	float uv[2];
} sv_record;

static void test_layouts_and_mutation(void)
{
	double p[SV_N * 3], pf[SV_N * 3], px[SV_N], py[SV_N], pz[SV_N];
	sv_record rec[SV_N];
	qaws_surface_view v;
	qaws_surface_prepared *a = NULL, *b = NULL, *c = NULL;
	qaws_surface* owned;
	unsigned int i, j, diff = 0;
	double e;
	char msg[200];
	sv_points(p, 1.3);
	for (i = 0; i < SV_N; i++)
	{
		memset(&rec[i], 0, sizeof rec[i]);
		for (j = 0; j < 3; j++)
		{
			rec[i].position[j] = (float)p[i * 3 + j];
			pf[i * 3 + j] = rec[i].position[j];
		}
		px[i] = pf[i * 3];
		py[i] = pf[i * 3 + 1];
		pz[i] = pf[i * 3 + 2];
	}
	v = sv_view(pf, 3, NULL, 1);
	TEST_ASSERT(qaws_surface_prepare(&b, &v, NULL) == QAWS_STATUS_OK, "prepare: packed f64");
	v.points = sv_packed(rec[0].position, QAWS_ELEMENT_F32, 3, SV_N, sizeof(sv_record));
	TEST_ASSERT(qaws_surface_prepare(&a, &v, NULL) == QAWS_STATUS_OK, "prepare: interleaved f32 records");
	memset(&v.points, 0, sizeof v.points);
	v.points.planes[0] = px;
	v.points.planes[1] = py;
	v.points.planes[2] = pz;
	v.points.type = QAWS_ELEMENT_F64;
	v.points.components = 3;
	v.points.count = SV_N;
	TEST_ASSERT(qaws_surface_prepare(&c, &v, NULL) == QAWS_STATUS_OK, "prepare: planar f64");
	if (!a || !b || !c)
		return;
	for (i = 0; i <= 20; i++)
		for (j = 0; j <= 20; j++)
		{
			qaws_surface_eval_result ra, rb, rc;
			qaws_surface_prepared_evaluate(&ra, a, (qaws_scalar)(i / 20.0), (qaws_scalar)(j / 20.0), SV_ALL);
			qaws_surface_prepared_evaluate(&rb, b, (qaws_scalar)(i / 20.0), (qaws_scalar)(j / 20.0), SV_ALL);
			qaws_surface_prepared_evaluate(&rc, c, (qaws_scalar)(i / 20.0), (qaws_scalar)(j / 20.0), SV_ALL);
			diff += memcmp(&ra, &rb, sizeof ra) != 0 || memcmp(&rb, &rc, sizeof rb) != 0;
		}
	TEST_ASSERT(diff == 0, "surface: f32 records, packed f64 and planar f64 holding the same values evaluate bit for bit alike");
	/* the host rewrites the planar points in place; c is kept */
	sv_points(p, 2.2);
	for (i = 0; i < SV_N; i++)
	{
		px[i] = p[i * 3];
		py[i] = p[i * 3 + 1];
		pz[i] = p[i * 3 + 2];
	}
	owned = sv_owned(p, NULL, 1);
	e = owned ? sv_compare(c, NULL, owned) : 1e30;
	sprintf(msg, "surface points rewritten in place without preparing again: %.2g", e);
	TEST_ASSERT(e <= SV_TOL, msg);
	qaws_surface_destroy(owned);
	qaws_surface_prepared_destroy(a);
	qaws_surface_prepared_destroy(b);
	qaws_surface_prepared_destroy(c);
}

#define SV_RANGE 2048
#define SV_THREADS 4

typedef struct sv_job
{
	qaws_surface_prepared const* p;
	qaws_array_view const* params;
	qaws_surface_outputs const* out;
	unsigned int begin, end;
	unsigned int status;
} sv_job;

static void sv_range_job(void* a)
{
	sv_job* j = (sv_job*)a;
	unsigned int r;
	for (r = 0; r < 3; r++)
		j->status |= (unsigned int)qaws_surface_prepared_evaluate_range(j->p, j->params, j->begin, j->end, j->out);
}

static void test_ranges(void)
{
	static double uv[SV_RANGE * 2], pos_s[SV_RANGE * 3], nrm_s[SV_RANGE * 3], pos[SV_RANGE * 3], nrm[SV_RANGE * 3];
	double p[SV_N * 3], w[SV_N];
	qaws_surface_view v;
	qaws_surface_prepared* pr = NULL;
	qaws_array_view pv;
	qaws_surface_outputs o;
	sv_job job[SV_THREADS];
	void* args[SV_THREADS];
	unsigned int i, bad = 0, st = 0;
	char msg[200];
	sv_points(p, 0.1);
	sv_weights(w);
	for (i = 0; i < SV_RANGE; i++)
	{
		uv[i * 2] = fmod(0.618034 * i, 1.0);
		uv[i * 2 + 1] = fmod(0.414214 * i, 1.0);
	}
	v = sv_view(p, 3, w, 1);
	if (qaws_surface_prepare(&pr, &v, NULL) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "surface ranges: prepared");
		return;
	}
	pv = sv_packed(uv, QAWS_ELEMENT_F64, 2, SV_RANGE, 0);
	memset(&o, 0, sizeof o);
	o.position.data = pos_s;
	o.position.type = QAWS_ELEMENT_F64;
	o.position.components = 3;
	o.position.count = SV_RANGE;
	o.normal = o.position;
	o.normal.data = nrm_s;
	TEST_ASSERT(qaws_surface_prepared_evaluate_range(pr, &pv, 0, SV_RANGE, &o) == QAWS_STATUS_OK, "surface range evaluation");
	for (i = 0; i < SV_RANGE; i += 5)
	{
		qaws_surface_eval_result r;
		qaws_surface_prepared_evaluate(&r, pr, (qaws_scalar)uv[i * 2], (qaws_scalar)uv[i * 2 + 1], QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
		bad += (double)r.position.x != pos_s[i * 3] || (double)r.position.z != pos_s[i * 3 + 2] || (double)r.normal.y != nrm_s[i * 3 + 1];
	}
	TEST_ASSERT(bad == 0, "surface range evaluation: the single evaluations' values");
	o.position.data = pos;
	o.normal.data = nrm;
	for (i = 0; i < SV_THREADS; i++)
	{
		job[i].p = pr;
		job[i].params = &pv;
		job[i].out = &o;
		job[i].begin = SV_RANGE * i / SV_THREADS;
		job[i].end = SV_RANGE * (i + 1) / SV_THREADS;
		job[i].status = 0;
		args[i] = &job[i];
	}
	test_run_threads(SV_THREADS, sv_range_job, args);
	for (i = 0; i < SV_THREADS; i++)
		st |= job[i].status;
	sprintf(msg, "surface: %u threads over disjoint ranges of %u parameters, bit for bit the serial outputs", SV_THREADS, SV_RANGE);
	TEST_ASSERT(st == 0 && memcmp(pos, pos_s, sizeof pos) == 0 && memcmp(nrm, nrm_s, sizeof nrm) == 0, msg);
	pv.components = 1;
	TEST_ASSERT(qaws_surface_prepared_evaluate_range(pr, &pv, 0, SV_RANGE, &o) == QAWS_STATUS_INVALID_ARGUMENT, "surface range: one-component parameters refused");
	qaws_surface_prepared_destroy(pr);
}

typedef struct sv_alloc_state
{
	int calls, fail_at, live;
} sv_alloc_state;

static void* sv_alloc(unsigned long size, void* u)
{
	sv_alloc_state* s = (sv_alloc_state*)u;
	if (s->calls++ == s->fail_at)
		return NULL;
	s->live++;
	return malloc(size);
}

static void sv_dealloc(void* ptr, void* u)
{
	if (ptr)
		((sv_alloc_state*)u)->live--;
	free(ptr);
}

static void test_allocation_and_validation(void)
{
	double p[SV_N * 3], w[SV_N];
	qaws_surface_view v;
	sv_alloc_state s;
	qaws_allocator al;
	unsigned int f, which, leaks = 0, wrong = 0, ok[3] = { 0, 0, 0 };
	qaws_surface_prepared *pi = (qaws_surface_prepared*)&v, *pa = NULL;
	unsigned char* mem;
	size_t bytes;
	char msg[256];
	sv_points(p, 0.6);
	sv_weights(w);
	v = sv_view(p, 3, w, 1);
	al.alloc = sv_alloc;
	al.dealloc = sv_dealloc;
	al.user_data = &s;
	for (which = 0; which < 3; which++)
		for (f = 0; f < 16 && !ok[which]; f++)
		{
			qaws_surface_prepared* pr = NULL;
			qaws_surface* sf = NULL;
			qaws_status st;
			s.calls = 0;
			s.fail_at = (int)f;
			s.live = 0;
			st = which == 0 ? qaws_surface_prepare(&pr, &v, &al) : which == 1 ? qaws_surface_create_borrowed(&sf, &v, &al)
				: qaws_surface_create_from_view(&sf, &v, &al);
			if (st == QAWS_STATUS_OK)
			{
				ok[which] = f + 1;
				qaws_surface_prepared_destroy(pr);
				qaws_surface_destroy(sf);
			}
			else
				wrong += st != QAWS_STATUS_ALLOCATION_FAILURE || pr || sf;
			leaks += s.live != 0;
		}
	sprintf(msg, "surface allocation failures: prepare, borrowed, copy succeed once %u, %u, %u allocations pass; %u wrong statuses, %u leaks",
		ok[0] - 1, ok[1] - 1, ok[2] - 1, wrong, leaks);
	TEST_ASSERT(ok[0] == 2 && ok[1] && ok[2] && wrong == 0 && leaks == 0, msg);
	bytes = qaws_surface_prepared_size(&v);
	mem = (unsigned char*)malloc(bytes + 32);
	{
		unsigned char* at = (unsigned char*)(((size_t)mem + 15) & ~(size_t)15);
		TEST_ASSERT(qaws_surface_prepare_into(&pi, at, &v, bytes - 1) == QAWS_STATUS_BUFFER_TOO_SMALL && !pi, "surface prepare_into: one byte short refused");
		TEST_ASSERT(qaws_surface_prepare_into(&pi, at, &v, bytes) == QAWS_STATUS_OK && qaws_surface_prepare(&pa, &v, NULL) == QAWS_STATUS_OK,
			"surface prepare_into and prepare");
		if (pi && pa)
		{
			qaws_surface_eval_result a, b;
			qaws_surface_prepared_evaluate(&a, pi, (qaws_scalar)0.37, (qaws_scalar)0.81, SV_ALL);
			qaws_surface_prepared_evaluate(&b, pa, (qaws_scalar)0.37, (qaws_scalar)0.81, SV_ALL);
			TEST_ASSERT(memcmp(&a, &b, sizeof a) == 0, "surface prepare_into evaluates as prepare");
		}
		qaws_surface_prepared_destroy(pi);
		qaws_surface_prepared_destroy(pa);
	}
	free(mem);
	/* validation */
	v.u_count = SV_U - 1;
	TEST_ASSERT(qaws_surface_view_validate(&v) == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT, "grid not matching the points refused");
	v.u_count = SV_U;
	v.v_knots.count = 7;
	TEST_ASSERT(qaws_surface_view_validate(&v) == QAWS_STATUS_INVALID_KNOT_COUNT, "v knot count refused");
	v.v_knots.count = SV_V + SV_VD + 1;
	v.u_degree = 0;
	TEST_ASSERT(qaws_surface_view_validate(&v) == QAWS_STATUS_INVALID_DEGREE, "degree 0 refused");
	v.u_degree = SV_UD;
	w[7] = -1;
	TEST_ASSERT(qaws_surface_view_validate(&v) == QAWS_STATUS_INVALID_ARGUMENT, "negative weight refused");
	w[7] = 1;
	TEST_ASSERT(qaws_surface_view_validate(&v) == QAWS_STATUS_OK, "valid surface view");
}

static void test_operations(void)
{
	double p[SV_N * 3];
	qaws_surface_view v;
	qaws_surface *owned, *borrowed = NULL;
	qaws_scalar ao = 0, ab = 0;
	qaws_vec3 q;
	qaws_surface_closest_point co, cb;
	char msg[256];
	sv_points(p, 0.8);
	v = sv_view(p, 3, NULL, 1);
	owned = sv_owned(p, NULL, 1);
	if (!owned || qaws_surface_create_borrowed(&borrowed, &v, NULL) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "surface operations: created");
		return;
	}
	TEST_ASSERT(qaws_surface_compute_area(&ao, owned) == QAWS_STATUS_OK && qaws_surface_compute_area(&ab, borrowed) == QAWS_STATUS_OK, "areas computed");
	sprintf(msg, "area: owned %.12g, borrowed %.12g", (double)ao, (double)ab);
	TEST_ASSERT(fabs((double)ao - ab) <= SV_TOL * (1 + fabs((double)ao)), msg);
	q.x = (qaws_scalar)0.42;
	q.y = (qaws_scalar)0.57;
	q.z = (qaws_scalar)0.6;
	memset(&co, 0, sizeof co);
	memset(&cb, 0, sizeof cb);
	TEST_ASSERT(qaws_surface_closest_point_tangent(&co, NULL, owned, NULL, NULL, q, NULL) == QAWS_STATUS_OK
		&& qaws_surface_closest_point_tangent(&cb, NULL, borrowed, NULL, NULL, q, NULL) == QAWS_STATUS_OK, "closest points found");
	sprintf(msg, "closest point on the borrowed surface: (%.12g, %.12g), owned (%.12g, %.12g)", (double)cb.u, (double)cb.v, (double)co.u, (double)co.v);
	TEST_ASSERT(fabs((double)co.u - cb.u) + fabs((double)co.v - cb.v) + fabs((double)co.distance - cb.distance) <= SV_TOL * 10, msg);
	qaws_surface_destroy(borrowed);
	qaws_surface_destroy(owned);
}

/* third-order jets (the closest point and curvature queries read them) against the owned surfaces' */
static void test_jets(void)
{
	double p[SV_N * 3], w[SV_N], worst = 0;
	unsigned int c, i, j, b;
	char msg[200];
	sv_points(p, 1.7);
	sv_weights(w);
	for (c = 0; c < 2; c++)
	{
		qaws_surface_view v = sv_view(p, 3, c ? w : NULL, 1);
		qaws_surface *owned = sv_owned(p, c ? w : NULL, 1), *borrowed = NULL;
		if (!owned || qaws_surface_create_borrowed(&borrowed, &v, NULL) != QAWS_STATUS_OK)
		{
			TEST_ASSERT(0, "jets: surfaces created");
			qaws_surface_destroy(owned);
			return;
		}
		for (i = 0; i <= 9; i++)
			for (j = 0; j <= 9; j++)
			{
				qaws_surface_jet a, o;
				if (qaws_surface_eval_jet(&a, borrowed, (qaws_scalar)(i / 9.0), (qaws_scalar)(j / 9.0), QAWS_SJET_ORDER3) != QAWS_STATUS_OK
					|| qaws_surface_eval_jet(&o, owned, (qaws_scalar)(i / 9.0), (qaws_scalar)(j / 9.0), QAWS_SJET_ORDER3) != QAWS_STATUS_OK)
				{
					worst = 1e30;
					continue;
				}
				for (b = 0; b < QAWS_SURFACE_JET_COUNT; b++)
				{
					double e = (fabs((double)a.d[b].x - o.d[b].x) + fabs((double)a.d[b].y - o.d[b].y) + fabs((double)a.d[b].z - o.d[b].z))
						/ (1 + fabs((double)o.d[b].x) + fabs((double)o.d[b].y) + fabs((double)o.d[b].z));
					worst = e > worst ? e : worst;
				}
			}
		qaws_surface_destroy(borrowed);
		qaws_surface_destroy(owned);
	}
	sprintf(msg, "third-order jets of borrowed B-spline and NURBS surfaces against the owned: %.2g", worst);
	TEST_ASSERT(worst <= SV_TOL * 100, msg);
	printf("  %s\n", msg);
}

int test_143_surface_view_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 143: borrowed surfaces\n");
	test_parity();
	test_layouts_and_mutation();
	test_ranges();
	test_allocation_and_validation();
	test_operations();
	test_jets();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
