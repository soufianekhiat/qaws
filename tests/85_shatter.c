/*
 * Test 85: seams and shatter of curved surfaces (qaws_shatter.h)
 *
 *   - two square sheets crossing at right angles: one open seam along their
 *     line of crossing, its length, every point on both sheets;
 *   - a saddle (a bilinear patch) and a plane above its middle: two open
 *     seams (hyperbola branches), every point on both surfaces;
 *   - A's own seams: two crossing sheets given as one operand (one A-A
 *     seam), a cube of six faces (none: they only share borders);
 *   - shatter: a square sheet by a cube (six bilinear faces, a solid): two
 *     pieces, the one inside the cube a quarter of the sheet; the saddle by
 *     the plane: three pieces, not told inside or outside;
 *   - 2D regions of curves: two circles' crossings, a bowtie's own crossing,
 *     a short buffer; a circle shattered by a circle (the lens inside, the
 *     crescent outside, areas exact), by an open cut, a ring by a line;
 *   - refusals.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#if QAWS_SCALAR_IS_FLOAT
#define SHT_TOL 2e-4
#else
#define SHT_TOL 1e-8
#endif

static qaws_vec3 sht_v(double x, double y, double z)
{
	qaws_vec3 v;
	v.x = (qaws_scalar)x;
	v.y = (qaws_scalar)y;
	v.z = (qaws_scalar)z;
	return v;
}

static qaws_surface* sht_bilinear(qaws_vec3 p00, qaws_vec3 p10, qaws_vec3 p01, qaws_vec3 p11)
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

/* the six faces of the box [lo, hi], normals outward (S_u x S_v) */
static void sht_box(double const* lo, double const* hi, qaws_surface** f)
{
	double x0 = lo[0], y0 = lo[1], z0 = lo[2], x1 = hi[0], y1 = hi[1], z1 = hi[2];
	f[0] = sht_bilinear(sht_v(x0, y0, z0), sht_v(x0, y1, z0), sht_v(x1, y0, z0), sht_v(x1, y1, z0));   /* z = z0, -z */
	f[1] = sht_bilinear(sht_v(x0, y0, z1), sht_v(x1, y0, z1), sht_v(x0, y1, z1), sht_v(x1, y1, z1));   /* z = z1, +z */
	f[2] = sht_bilinear(sht_v(x0, y0, z0), sht_v(x1, y0, z0), sht_v(x0, y0, z1), sht_v(x1, y0, z1));   /* y = y0, -y */
	f[3] = sht_bilinear(sht_v(x0, y1, z0), sht_v(x0, y1, z1), sht_v(x1, y1, z0), sht_v(x1, y1, z1));   /* y = y1, +y */
	f[4] = sht_bilinear(sht_v(x0, y0, z0), sht_v(x0, y0, z1), sht_v(x0, y1, z0), sht_v(x0, y1, z1));   /* x = x0, -x */
	f[5] = sht_bilinear(sht_v(x1, y0, z0), sht_v(x1, y1, z0), sht_v(x1, y0, z1), sht_v(x1, y1, z1));   /* x = x1, +x */
}

/* the largest gap between the two surfaces at a seam's points, and its length */
static void sht_check_seam(qaws_surface_seams const* S, unsigned int i, qaws_surface const* const* a, qaws_surface const* const* b, double* gap,
	double* len, qaws_seam_info* info)
{
	qaws_seam_point* p;
	unsigned int j;
	qaws_seam_info in;
	qaws_surface_seams_get(&in, NULL, S, i);
	p = (qaws_seam_point*)malloc(sizeof(qaws_seam_point) * (in.point_count + 1));
	qaws_surface_seams_get(NULL, p, S, i);
	*gap = 0;
	*len = 0;
	for (j = 0; j < in.point_count; j++)
	{
		qaws_surface_eval_result r1, r2;
		qaws_surface const* s1 = in.operand1 ? b[in.surface1] : a[in.surface1];
		qaws_surface const* s2 = in.operand2 ? b[in.surface2] : a[in.surface2];
		double dx, dy, dz;
		qaws_surface_evaluate(&r1, s1, p[j].u1, p[j].v1, QAWS_EVAL_FLAG_POSITION);
		qaws_surface_evaluate(&r2, s2, p[j].u2, p[j].v2, QAWS_EVAL_FLAG_POSITION);
		dx = (double)r1.position.x - (double)r2.position.x;
		dy = (double)r1.position.y - (double)r2.position.y;
		dz = (double)r1.position.z - (double)r2.position.z;
		*gap = fmax(*gap, sqrt(dx * dx + dy * dy + dz * dz));
		if (j + 1 < in.point_count + (in.closed ? 1u : 0u))
		{
			qaws_seam_point const* q = &p[(j + 1) % in.point_count];
			dx = (double)q->position.x - (double)p[j].position.x;
			dy = (double)q->position.y - (double)p[j].position.y;
			dz = (double)q->position.z - (double)p[j].position.z;
			*len += sqrt(dx * dx + dy * dy + dz * dz);
		}
	}
	free(p);
	*info = in;
}

static void test_seams(void)
{
	qaws_surface *sh = sht_bilinear(sht_v(-1, -1, 0), sht_v(1, -1, 0), sht_v(-1, 1, 0), sht_v(1, 1, 0));
	qaws_surface *up = sht_bilinear(sht_v(0.13, -1.3, -1), sht_v(0.13, 1.3, -1), sht_v(0.13, -1.3, 1), sht_v(0.13, 1.3, 1));
	qaws_surface *saddle = sht_bilinear(sht_v(-1, -1, 0.5), sht_v(1, -1, -0.5), sht_v(-1, 1, -0.5), sht_v(1, 1, 0.5));
	qaws_surface *plane = sht_bilinear(sht_v(-2, -2, 0.2), sht_v(2, -2, 0.2), sht_v(-2, 2, 0.2), sht_v(2, 2, 0.2));
	qaws_surface* cube[6];
	double const lo[3] = { -0.5, -0.5, -0.5 }, hi[3] = { 0.5, 0.5, 0.5 };
	qaws_surface_seam_desc d;
	qaws_surface_seams* S = NULL;
	qaws_seam_info in;
	double gap, len;
	unsigned int i;
	char msg[256];
	sht_box(lo, hi, cube);
	/* two sheets crossing along x = 0.13, y in [-1, 1] */
	memset(&d, 0, sizeof d);
	d.a = (qaws_surface const* const*)&sh;
	d.a_count = 1;
	d.b = (qaws_surface const* const*)&up;
	d.b_count = 1;
	d.kinds = QAWS_SEAM_A_B;
	TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
	TEST_ASSERT(qaws_surface_seams_get_count(S) == 1, "two sheets: one seam");
	if (qaws_surface_seams_get_count(S) == 1)
	{
		sht_check_seam(S, 0, d.a, d.b, &gap, &len, &in);
		sprintf(msg, "two sheets: open %d, kind %u, length %.12f (2), gap %.3g", !in.closed, in.kind, len, gap);
		printf("  %s\n", msg);
		TEST_ASSERT(!in.closed && in.kind == QAWS_SEAM_A_B && in.operand1 == 0 && in.operand2 == 1 && fabs(len - 2) < 0.02 && gap < SHT_TOL, msg);
	}
	qaws_surface_seams_destroy(S);
	/* a saddle and a plane: two hyperbola branches */
	d.a = (qaws_surface const* const*)&saddle;
	d.b = (qaws_surface const* const*)&plane;
	S = NULL;
	TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
	sprintf(msg, "saddle and plane: %u seams (2)", qaws_surface_seams_get_count(S));
	printf("  %s\n", msg);
	TEST_ASSERT(qaws_surface_seams_get_count(S) == 2, msg);
	for (i = 0; i < qaws_surface_seams_get_count(S); i++)
	{
		sht_check_seam(S, i, d.a, d.b, &gap, &len, &in);
		sprintf(msg, "saddle and plane, seam %u: %u points, open %d, length %.6f, gap %.3g", i, in.point_count, !in.closed, len, gap);
		printf("  %s\n", msg);
		TEST_ASSERT(!in.closed && in.point_count > 10 && gap < SHT_TOL, msg);
	}
	qaws_surface_seams_destroy(S);
	/* A's own: two sheets as one operand cross; a cube's faces only share borders */
	{
		qaws_surface* two[2];
		two[0] = sh;
		two[1] = up;
		memset(&d, 0, sizeof d);
		d.a = (qaws_surface const* const*)two;
		d.a_count = 2;
		d.kinds = QAWS_SEAM_A_A;
		S = NULL;
		TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
		TEST_ASSERT(qaws_surface_seams_get_count(S) == 1, "two sheets as one operand: one A-A seam");
		if (qaws_surface_seams_get_count(S) == 1)
		{
			sht_check_seam(S, 0, d.a, d.a, &gap, &len, &in);
			TEST_ASSERT(in.kind == QAWS_SEAM_A_A && fabs(len - 2) < 0.02 && gap < SHT_TOL, "two sheets as one operand: the A-A seam on both");
		}
		qaws_surface_seams_destroy(S);
		d.a = (qaws_surface const* const*)cube;
		d.a_count = 6;
		S = NULL;
		TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
		sprintf(msg, "a cube's faces: %u A-A seams (0)", qaws_surface_seams_get_count(S));
		printf("  %s\n", msg);
		TEST_ASSERT(qaws_surface_seams_get_count(S) == 0, msg);
		qaws_surface_seams_destroy(S);
	}
	/* refusals */
	S = (qaws_surface_seams*)1;
	memset(&d, 0, sizeof d);
	TEST_ASSERT(qaws_surface_seams_create(&S, &d) == QAWS_STATUS_INVALID_ARGUMENT && S == NULL, "seams: no surface refused");
	d.a = (qaws_surface const* const*)&sh;
	d.a_count = 1;
	TEST_ASSERT(qaws_surface_seams_create(&S, &d) == QAWS_STATUS_INVALID_ARGUMENT, "seams: no kind refused");
	qaws_surface_destroy(sh);
	qaws_surface_destroy(up);
	qaws_surface_destroy(saddle);
	qaws_surface_destroy(plane);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(cube[i]);
}

/* the share of a trimmed piece's (u, v) domain inside its trim, sampled on a grid */
static double sht_trim_share(qaws_surface const* piece)
{
	unsigned int i, j, in_count = 0, n = 100;
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			int in = 0;
			qaws_surface_trim_contains(&in, piece, (qaws_scalar)((i + 0.5) / n), (qaws_scalar)((j + 0.5) / n));
			in_count += in != 0;
		}
	return (double)in_count / (n * n);
}

static void test_shatter(void)
{
	qaws_surface* sh = sht_bilinear(sht_v(-1, -1, 0.05), sht_v(1, -1, 0.05), sht_v(-1, 1, 0.05), sht_v(1, 1, 0.05));
	qaws_surface* saddle = sht_bilinear(sht_v(-1, -1, 0.5), sht_v(1, -1, -0.5), sht_v(-1, 1, -0.5), sht_v(1, 1, 0.5));
	qaws_surface* plane = sht_bilinear(sht_v(-2, -2, 0.2), sht_v(2, -2, 0.2), sht_v(-2, 2, 0.2), sht_v(2, 2, 0.2));
	qaws_surface* cube[6];
	double const lo[3] = { -0.5, -0.5, -0.5 }, hi[3] = { 0.5, 0.5, 0.5 };
	qaws_surface_seam_desc d;
	qaws_surface_shatter* H = NULL;
	unsigned int i, n;
	double share_in = 0, share_out = 0, total = 0;
	int bad = 0;
	char msg[256];
	sht_box(lo, hi, cube);
	memset(&d, 0, sizeof d);
	d.a = (qaws_surface const* const*)&sh;
	d.a_count = 1;
	d.b = (qaws_surface const* const*)cube;
	d.b_count = 6;
	TEST_ASSERT_STATUS(qaws_surface_shatter_create(&H, &d, 1));
	n = qaws_surface_shatter_get_count(H);
	for (i = 0; i < n; i++)
	{
		qaws_surface const* p = NULL;
		unsigned int from = 99;
		int inside = -2;
		double s;
		qaws_surface_shatter_get(&p, &from, &inside, H, i);
		s = sht_trim_share(p);
		bad += from != 0 || inside < 0;
		if (inside == 1)
			share_in += s;
		else
			share_out += s;
		total += s;
	}
	sprintf(msg, "sheet by a cube: %u pieces (2), inside %.4f (0.25), outside %.4f (0.75)", n, share_in, share_out);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 2 && !bad && fabs(share_in - 0.25) < 0.01 && fabs(share_out - 0.75) < 0.01, msg);
	qaws_surface_shatter_destroy(H);
	/* the saddle by the plane: three pieces, not told */
	d.a = (qaws_surface const* const*)&saddle;
	d.b = (qaws_surface const* const*)&plane;
	d.b_count = 1;
	H = NULL;
	TEST_ASSERT_STATUS(qaws_surface_shatter_create(&H, &d, 0));
	n = qaws_surface_shatter_get_count(H);
	total = 0;
	bad = 0;
	for (i = 0; i < n; i++)
	{
		qaws_surface const* p = NULL;
		int inside = 0;
		qaws_surface_shatter_get(&p, NULL, &inside, H, i);
		total += sht_trim_share(p);
		bad += inside != -1;
	}
	sprintf(msg, "saddle by a plane: %u pieces (3), their shares summing to %.4f (1)", n, total);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 3 && !bad && fabs(total - 1) < 0.01, msg);
	qaws_surface_shatter_destroy(H);
	/* refusals: B empty */
	d.b_count = 0;
	H = (qaws_surface_shatter*)1;
	TEST_ASSERT(qaws_surface_shatter_create(&H, &d, 0) == QAWS_STATUS_INVALID_ARGUMENT && H == NULL, "shatter: no B refused");
	qaws_surface_destroy(sh);
	qaws_surface_destroy(saddle);
	qaws_surface_destroy(plane);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(cube[i]);
}

/* ---- 2D regions of curves ---- */

static qaws_path_2d sht_path(qaws_curve const** c, unsigned int n, int closed)
{
	qaws_path_2d p;
	p.curves = c;
	p.curve_count = n;
	p.closed = closed;
	return p;
}

static void test_regions(void)
{
	qaws_curve *c0 = NULL, *c1 = NULL, *bow = NULL, *cut = NULL, *outer = NULL, *inner = NULL, *line = NULL;
	qaws_curve const *l0[1], *l1[1], *lb[1], *lc[1], *lo[1], *li[1], *ll[1];
	qaws_path_2d A, B, Bow, Cut, Ring[2], Line;
	qaws_region_seam_desc d;
	qaws_region_seam_point pts[16];
	qaws_region_shatter* H = NULL;
	unsigned int n = 0, i, k;
	double lens = 2 * acos(0.5) - 0.5 * sqrt(3.0), area_in = 0, area_out = 0, total = 0;
	int bad = 0;
	char msg[256];
	qaws_scalar const bow_pts[8] = { -1, -1, 1, 1, 1, -1, -1, 1 };
	qaws_scalar const cut_pts[8] = { -1.5, -0.3, -0.2, 0.25, 0.4, -0.2, 1.5, 0.35 };
	qaws_scalar const line_pts[4] = { 0.05, -3, -0.05, 3 };
	qaws_curve_create_ellipse_2d(&c0, 0, 0, 1, 1, 0);
	qaws_curve_create_ellipse_2d(&c1, 1, 0, 1, 1, 0);
	qaws_curve_create_polyline_2d(&bow, bow_pts, 4, 1);
	qaws_curve_create_polyline_2d(&cut, cut_pts, 4, 0);
	qaws_curve_create_ellipse_2d(&outer, 0, 0, 2, 2, 0);
	qaws_curve_create_ellipse_2d(&inner, 0, 0, 1, 1, 0);
	qaws_curve_create_polyline_2d(&line, line_pts, 2, 0);
	l0[0] = c0; l1[0] = c1; lb[0] = bow; lc[0] = cut; lo[0] = outer; ll[0] = line;
	{
		qaws_curve* rev = NULL;
		qaws_curve_reverse(&rev, inner);
		qaws_curve_destroy(inner);
		inner = rev;
		li[0] = inner;
	}
	A = sht_path(l0, 1, 1);
	B = sht_path(l1, 1, 1);
	Bow = sht_path(lb, 1, 1);
	Cut = sht_path(lc, 1, 0);
	Ring[0] = sht_path(lo, 1, 1);
	Ring[1] = sht_path(li, 1, 1);
	Line = sht_path(ll, 1, 0);
	/* two unit circles a radius apart: crossings at (1/2, +-sqrt(3)/2) */
	memset(&d, 0, sizeof d);
	d.a = &A; d.a_count = 1; d.a_rule = QAWS_FILL_NON_ZERO;
	d.b = &B; d.b_count = 1; d.b_rule = QAWS_FILL_NON_ZERO;
	d.b_is_region = 1;
	d.kinds = QAWS_SEAM_A_B;
	TEST_ASSERT_STATUS(qaws_region_seams_2d(pts, &n, &d, 16));
	for (i = 0; i < n; i++)
		bad += fabs((double)pts[i].position.x - 0.5) > SHT_TOL || fabs(fabs((double)pts[i].position.y) - sqrt(3.0) / 2) > SHT_TOL
			|| pts[i].kind != QAWS_SEAM_A_B;
	sprintf(msg, "2D seams of two circles: %u points (2), %d off", n, bad);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 2 && !bad, msg);
	/* a bowtie crosses itself once, at the origin */
	d.a = &Bow;
	d.b_count = 0;
	d.kinds = QAWS_SEAM_A_A;
	n = 0;
	TEST_ASSERT_STATUS(qaws_region_seams_2d(pts, &n, &d, 16));
	sprintf(msg, "2D seams of a bowtie: %u points (1) at (%.3g, %.3g)", n, n ? (double)pts[0].position.x : 9.0, n ? (double)pts[0].position.y : 9.0);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 1 && fabs((double)pts[0].position.x) < SHT_TOL && fabs((double)pts[0].position.y) < SHT_TOL && pts[0].kind == QAWS_SEAM_A_A, msg);
	/* buffer too small */
	d.a = &A; d.b = &B; d.b_count = 1; d.kinds = QAWS_SEAM_A_B;
	TEST_ASSERT(qaws_region_seams_2d(pts, &n, &d, 1) == QAWS_STATUS_BUFFER_TOO_SMALL && n == 2, "2D seams: a short buffer told");
	/* the first circle shattered by the second: the lens inside, the crescent outside, areas exact */
	TEST_ASSERT_STATUS(qaws_region_shatter_create_2d(&H, &d));
	n = qaws_region_shatter_get_count(H);
	bad = 0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d const* pp = NULL;
		unsigned int np = 0;
		int inside = -2;
		qaws_scalar a = 0;
		qaws_region_shatter_get(&pp, &np, &inside, H, i);
		qaws_region_compute_area_2d(&a, pp, np);
		if (inside == 1)
			area_in += (double)a;
		else if (inside == 0)
			area_out += (double)a;
		else
			bad++;
		for (k = 0; k < np; k++)
			bad += !pp[k].closed || !pp[k].curve_count;
	}
	sprintf(msg, "2D shatter of a circle by a circle: %u pieces (2), lens %.12f (%.12f), crescent %.12f (%.12f)", n, area_in, lens, area_out,
		3.14159265358979323846 - lens);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 2 && !bad && fabs(area_in - lens) < 100 * SHT_TOL && fabs(area_out - (3.14159265358979323846 - lens)) < 100 * SHT_TOL, msg);
	qaws_region_shatter_destroy(H);
	/* a circle by an open polyline (no inside): two pieces, the circle's area */
	d.b = &Cut;
	d.b_is_region = 0;
	H = NULL;
	TEST_ASSERT_STATUS(qaws_region_shatter_create_2d(&H, &d));
	n = qaws_region_shatter_get_count(H);
	bad = 0;
	total = 0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d const* pp = NULL;
		unsigned int np = 0;
		int inside = 0;
		qaws_scalar a = 0;
		qaws_region_shatter_get(&pp, &np, &inside, H, i);
		qaws_region_compute_area_2d(&a, pp, np);
		total += (double)a;
		bad += inside != -1;
	}
	sprintf(msg, "2D shatter of a circle by an open cut: %u pieces (2), area %.12f (pi)", n, total);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 2 && !bad && fabs(total - 3.14159265358979323846) < 100 * SHT_TOL, msg);
	qaws_region_shatter_destroy(H);
	/* a ring (a hole: a clockwise inner circle) cut by a line: two half rings, the ring's area */
	d.a = Ring;
	d.a_count = 2;
	d.b = &Line;
	H = NULL;
	TEST_ASSERT_STATUS(qaws_region_shatter_create_2d(&H, &d));
	n = qaws_region_shatter_get_count(H);
	total = 0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d const* pp = NULL;
		unsigned int np = 0;
		qaws_scalar a = 0;
		qaws_region_shatter_get(&pp, &np, NULL, H, i);
		qaws_region_compute_area_2d(&a, pp, np);
		total += (double)a;
	}
	sprintf(msg, "2D shatter of a ring by a line: %u pieces (2), area %.12f (3 pi)", n, total);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 2 && fabs(total - 3 * 3.14159265358979323846) < 100 * SHT_TOL, msg);
	qaws_region_shatter_destroy(H);
	/* refusals */
	d.b_count = 0;
	H = (qaws_region_shatter*)1;
	TEST_ASSERT(qaws_region_shatter_create_2d(&H, &d) == QAWS_STATUS_INVALID_ARGUMENT && H == NULL, "2D shatter: no B refused");
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(bow);
	qaws_curve_destroy(cut);
	qaws_curve_destroy(outer);
	qaws_curve_destroy(inner);
	qaws_curve_destroy(line);
}

int test_85_shatter_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 85: seams and shatter of curved surfaces\n");
	test_seams();
	test_shatter();
	test_regions();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
