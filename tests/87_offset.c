/*
 * Test 87: Offsetting paths of curves (qaws_offset)
 *
 *   - a 10 x 10 square grown by 1 with each join type: miter 144, square
 *     132 + 8 sqrt 2, bevel 142, round 140 + pi; shrunk by 1: 64; by 6: empty
 *   - a square with a hole: the outer grows, the hole shrinks
 *   - a segment of length 10 with each end type: butt 20, square 24, round
 *     20 + pi, joined (round) 20 + pi
 *   - circles: an arc-kind circle stays exact (concentric arcs), a NURBS
 *     circle's offset is fitted within the tolerance
 *   - an open Bezier with round ends: 2 d L + pi d^2
 *   - a thin spike with miter joins stays within the miter limit
 *   - a clockwise group grows like a counter-clockwise one
 *   - the delta callback, constant, gives the plain offset
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define OST_PI 3.14159265358979323846

#if QAWS_SCALAR_IS_FLOAT
#define OST_TOL 2e-3
#else
#define OST_TOL 1e-8
#endif

static double ost_area(qaws_clip_result const* r)
{
	unsigned int i, n = qaws_clip_result_get_path_count(r);
	double sum = 0.0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d p;
		qaws_scalar a = 0;
		qaws_clip_result_get_path(r, i, &p);
		qaws_path_compute_area_2d(&p, &a);
		sum += a;
	}
	return sum;
}

static qaws_curve* ost_poly(double const* xy, unsigned int n, int closed)
{
	qaws_scalar p[64];
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	qaws_curve_create_polyline_2d(p, n, closed, &c);
	return c;
}

static double ost_run(qaws_curve const** cs, unsigned int n, int closed, double delta, qaws_join_type jt, qaws_end_type et,
	unsigned int* out_paths)
{
	qaws_path_2d p[4];
	qaws_clip_result* r = NULL;
	unsigned int i;
	double a = -1e9;
	for (i = 0; i < n; i++) { p[i].curves = &cs[i]; p[i].curve_count = 1; p[i].closed = closed; }
	if (qaws_offset_paths(p, n, (qaws_scalar)delta, jt, et, &r) == QAWS_STATUS_OK)
	{
		a = ost_area(r);
		if (out_paths) *out_paths = qaws_clip_result_get_path_count(r);
	}
	qaws_clip_result_destroy(r);
	return a;
}

static void test_square(void)
{
	static double const sq[8] = { 0, 0, 10, 0, 10, 10, 0, 10 };
	qaws_curve* c = ost_poly(sq, 4, 1);
	qaws_curve const* cc = c;
	double am = ost_run(&cc, 1, 1, 1, QAWS_JOIN_MITER, QAWS_END_POLYGON, NULL);
	double as = ost_run(&cc, 1, 1, 1, QAWS_JOIN_SQUARE, QAWS_END_POLYGON, NULL);
	double ab = ost_run(&cc, 1, 1, 1, QAWS_JOIN_BEVEL, QAWS_END_POLYGON, NULL);
	double ar = ost_run(&cc, 1, 1, 1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, NULL);
	double ai = ost_run(&cc, 1, 1, -1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, NULL);
	unsigned int ne = 9;
	double az = ost_run(&cc, 1, 1, -6, QAWS_JOIN_MITER, QAWS_END_POLYGON, &ne);
	char msg[300];
	sprintf(msg, "square grown by 1: miter 144 (%.10g), square 132 + 8 sqrt 2 (%.10g), bevel 142 (%.10g), round 140 + pi (%.10g)", am, as, ab, ar);
	TEST_ASSERT(fabs(am - 144) < OST_TOL * 100 && fabs(as - (132 + 8 * sqrt(2.0))) < OST_TOL * 100 && fabs(ab - 142) < OST_TOL * 100 &&
		fabs(ar - (140 + OST_PI)) < OST_TOL * 100, msg);
	sprintf(msg, "square shrunk by 1: 64 (%.10g); by 6: nothing (%u paths)", ai, ne);
	TEST_ASSERT(fabs(ai - 64) < OST_TOL * 100 && ne == 0, msg);
	qaws_curve_destroy(c);
}

static void test_hole_and_orientation(void)
{
	static double const o[8] = { 0, 0, 10, 0, 10, 10, 0, 10 }, h[8] = { 3, 3, 3, 7, 7, 7, 7, 3 };
	static double const cw[8] = { 0, 0, 0, 10, 10, 10, 10, 0 };
	qaws_curve *co = ost_poly(o, 4, 1), *ch = ost_poly(h, 4, 1), *ccw = ost_poly(cw, 4, 1);
	qaws_curve const* cs[2];
	unsigned int np = 0;
	double a, b;
	char msg[256];
	cs[0] = co; cs[1] = ch;
	a = ost_run(cs, 2, 1, 1, QAWS_JOIN_MITER, QAWS_END_POLYGON, &np);
	sprintf(msg, "square with a hole grown by 1 (miter): 144 - 4 (%.10g), 2 paths (%u)", a, np);
	TEST_ASSERT(fabs(a - 140) < OST_TOL * 100 && np == 2, msg);
	cs[0] = ccw;
	b = ost_run(cs, 1, 1, 1, QAWS_JOIN_MITER, QAWS_END_POLYGON, NULL);
	sprintf(msg, "a clockwise square grows too, keeping its orientation (area %.10g)", b);
	TEST_ASSERT(fabs(b + 144) < OST_TOL * 100, msg);
	qaws_curve_destroy(co); qaws_curve_destroy(ch); qaws_curve_destroy(ccw);
}

static void test_open(void)
{
	static double const seg[4] = { 0, 0, 10, 0 };
	qaws_curve* c = ost_poly(seg, 2, 0);
	qaws_curve const* cc = c;
	double ab = ost_run(&cc, 1, 0, 1, QAWS_JOIN_ROUND, QAWS_END_BUTT, NULL);
	double as = ost_run(&cc, 1, 0, 1, QAWS_JOIN_ROUND, QAWS_END_SQUARE, NULL);
	double ar = ost_run(&cc, 1, 0, 1, QAWS_JOIN_ROUND, QAWS_END_ROUND, NULL);
	double aj = ost_run(&cc, 1, 0, 1, QAWS_JOIN_ROUND, QAWS_END_JOINED, NULL);
	char msg[300];
	sprintf(msg, "segment of length 10, delta 1: butt 20 (%.10g), square 24 (%.10g), round 20 + pi (%.10g), joined round 20 + pi (%.10g)",
		ab, as, ar, aj);
	TEST_ASSERT(fabs(ab - 20) < OST_TOL * 100 && fabs(as - 24) < OST_TOL * 100 && fabs(ar - (20 + OST_PI)) < OST_TOL * 100 &&
		fabs(aj - (20 + OST_PI)) < OST_TOL * 100, msg);
	qaws_curve_destroy(c);
}

static void test_curves(void)
{
	char msg[300];
	/* an arc-kind circle of radius 5: concentric arcs, exact */
	{
		qaws_arc_segment s[2];
		qaws_arc_desc ad;
		qaws_curve* c = NULL;
		qaws_curve const* cc;
		double a;
		memset(s, 0, sizeof(s));
		s[0].radius = 5; s[0].angle_start = 0; s[0].angle_end = (qaws_scalar)OST_PI;
		s[1].radius = 5; s[1].angle_start = (qaws_scalar)OST_PI; s[1].angle_end = (qaws_scalar)(2 * OST_PI);
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = s; ad.segment_count = 2;
		qaws_curve_create_arc(&ad, &c);
		cc = c;
		a = ost_run(&cc, 1, 1, 1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, NULL);
		sprintf(msg, "arc circle r 5 grown by 1: 36 pi (%.12g vs %.12g)", a, 36 * OST_PI);
		TEST_ASSERT(fabs(a - 36 * OST_PI) < OST_TOL * 100, msg);
		qaws_curve_destroy(c);
	}
	/* a NURBS circle: fitted offset */
	{
		qaws_curve* c = NULL;
		qaws_curve const* cc;
		double a, b;
		qaws_curve_create_ellipse_2d(0, 0, 5, 5, 0, &c);
		cc = c;
		a = ost_run(&cc, 1, 1, 1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, NULL);
		b = ost_run(&cc, 1, 1, -2, QAWS_JOIN_ROUND, QAWS_END_POLYGON, NULL);
		sprintf(msg, "NURBS circle r 5: grown by 1 36 pi (%.10g), shrunk by 2 9 pi (%.10g), within the fit tolerance", a, b);
		TEST_ASSERT(fabs(a - 36 * OST_PI) < 1e-4 && fabs(b - 9 * OST_PI) < 1e-4, msg);
		qaws_curve_destroy(c);
	}
	/* an open gentle Bezier with round ends: 2 d L + pi d^2 */
	{
		static double const cp[8] = { 0, 0, 4, 3, 8, -3, 12, 0 };
		qaws_scalar p[8];
		qaws_bezier_desc d;
		qaws_curve* c = NULL;
		qaws_curve const* cc;
		qaws_path_2d path;
		qaws_scalar L = 0;
		double a, want;
		unsigned int i;
		for (i = 0; i < 8; i++) p[i] = (qaws_scalar)cp[i];
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p; d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &c);
		cc = c;
		path.curves = &cc; path.curve_count = 1; path.closed = 0;
		qaws_path_compute_length_2d(&path, &L);
		{
			/* the reference length by a dense chord sum, against the path length
			   (agreeing to 1e-9) */
			double acc = 0.0, px = 0, py = 0;
			unsigned int k;
			for (k = 0; k <= 200000; k++)
			{
				qaws_eval_result_2d ev;
				qaws_curve_evaluate_2d(c, (qaws_scalar)(k / 200000.0), QAWS_EVAL_FLAG_POSITION, &ev);
				if (k) acc += hypot(ev.position.x - px, ev.position.y - py);
				px = ev.position.x; py = ev.position.y;
			}
			sprintf(msg, "path length by Gauss-Legendre against a dense chord sum (%.12g vs %.12g)", (double)L, acc);
			TEST_ASSERT(fabs(L - acc) < 1e-8 * acc + (double)OST_TOL, msg);
			L = (qaws_scalar)acc;
		}
		a = ost_run(&cc, 1, 0, 0.5, QAWS_JOIN_ROUND, QAWS_END_ROUND, NULL);
		want = 2 * 0.5 * L + OST_PI * 0.25;
		sprintf(msg, "open Bezier, round ends, delta 0.5: 2 d L + pi d^2 (%.10g vs %.10g)", a, want);
		TEST_ASSERT(fabs(a - want) < 1e-5 * want, msg);
		qaws_curve_destroy(c);
	}
}

static qaws_scalar ost_const(void* user, unsigned int group, unsigned int path, qaws_vec2 point, qaws_vec2 normal)
{
	(void)user; (void)group; (void)path; (void)point; (void)normal;
	return 1;
}

static void test_miter_and_callback(void)
{
	/* a thin isosceles triangle (apex angle about 11 degrees) */
	static double const tri[6] = { 0, 0, 10, 0, 5, 50 };
	qaws_curve* c = ost_poly(tri, 3, 1);
	qaws_curve const* cc = c;
	qaws_path_2d p;
	qaws_offset_group g;
	qaws_offset_desc d;
	qaws_clip_result* r = NULL;
	qaws_vec2 lo, hi;
	double am, ac, ap;
	char msg[256];
	unsigned int n = 0;
	p.curves = &cc; p.curve_count = 1; p.closed = 1;
	memset(&d, 0, sizeof(d));
	g.paths = &p; g.path_count = 1; g.join_type = QAWS_JOIN_MITER; g.end_type = QAWS_END_POLYGON;
	d.groups = &g; d.group_count = 1; d.delta = 1; d.miter_limit = 2;
	qaws_offset_execute(&d, &r);
	lo.y = hi.y = 0;
	if (r && qaws_clip_result_get_path_count(r) == 1)
	{
		qaws_path_2d o;
		qaws_clip_result_get_path(r, 0, &o);
		qaws_path_compute_bounds_2d(&o, &lo, &hi);
		n = 1;
	}
	am = r ? ost_area(r) : 0;
	qaws_clip_result_destroy(r);
	sprintf(msg, "thin triangle, miter limit 2: the apex squared off within 2 deltas (top %.6g <= 52), one path (%u)", (double)hi.y, n);
	TEST_ASSERT(n == 1 && hi.y <= 52 + 1e-9 && hi.y > 50.9, msg);
	/* constant callback = plain delta */
	r = NULL;
	d.delta_fn = ost_const;
	qaws_offset_execute(&d, &r);
	ac = r ? ost_area(r) : 0;
	qaws_clip_result_destroy(r);
	ap = ost_run(&cc, 1, 1, 1, QAWS_JOIN_MITER, QAWS_END_POLYGON, NULL);
	sprintf(msg, "a constant delta callback gives the plain offset (%.10g vs %.10g, miter run %.10g)", ac, ap, am);
	TEST_ASSERT(fabs(ac - ap) < OST_TOL * 100 && fabs(am - ap) < OST_TOL * 100, msg);
	qaws_curve_destroy(c);
}

int test_87_offset_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 87: Offsetting paths of curves\n");
	test_square();
	test_hole_and_orientation();
	test_open();
	test_curves();
	test_miter_and_callback();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
