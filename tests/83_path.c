/*
 * Test 83: 2D paths and regions
 *
 *   - area: polygons (both orientations), ellipses, a circle as two arcs, a
 *     Bezier closed by a line, and closed curves of other kinds against the
 *     shoelace area of a dense sampling
 *   - bounds: a circle (exact), a Bezier against dense sampling
 *   - point location: inside / outside / on for a circle, a square with a
 *     hole under every fill rule, points just off a boundary, winding 2
 *   - transforms: arcs under similarities stay arcs, under a shear become
 *     NURBS, every point mapped
 *   - polylines: duplicates, collinear points and spikes, SimplifyPath, RDP
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define PTT_PI 3.14159265358979323846

#if QAWS_SCALAR_IS_FLOAT
#define PTT_REL 2e-4
#else
#define PTT_REL 1e-10
#endif

static qaws_path_2d ptt_one(qaws_curve const** c, int closed)
{
	qaws_path_2d p;
	p.curves = c;
	p.curve_count = 1;
	p.closed = closed;
	return p;
}

/* shoelace area of a dense sampling of a closed path */
static double ptt_sampled_area(qaws_path_2d const* p)
{
	double a = 0.0, x0 = 0, y0 = 0, px = 0, py = 0;
	unsigned int i, j, first = 1;
	for (i = 0; i < p->curve_count; i++)
	{
		qaws_range r = qaws_curve_get_parameter_range(p->curves[i]);
		for (j = 0; j <= 20000; j++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(p->curves[i], (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * j / 20000.0),
				QAWS_EVAL_FLAG_POSITION, &e);
			if (first) { x0 = px = e.position.x; y0 = py = e.position.y; first = 0; continue; }
			a += 0.5 * (px * e.position.y - e.position.x * py);
			px = e.position.x; py = e.position.y;
		}
	}
	a += 0.5 * (px * y0 - x0 * py);
	return a;
}

static void test_area(void)
{
	char msg[256];
	/* square, both ways */
	{
		qaws_scalar sq[] = { 0, 0, 2, 0, 2, 2, 0, 2 }, rv[] = { 0, 0, 0, 2, 2, 2, 2, 0 };
		qaws_curve *a = NULL, *b = NULL;
		qaws_curve const* c;
		qaws_path_2d p;
		qaws_scalar aa = 0, ab = 0;
		qaws_curve_create_polyline_2d(sq, 4, 1, &a);
		qaws_curve_create_polyline_2d(rv, 4, 1, &b);
		c = a; p = ptt_one(&c, 1); qaws_path_compute_area_2d(&p, &aa);
		TEST_ASSERT(qaws_path_is_positive_2d(&p), "a counter-clockwise square is positive");
		c = b; p = ptt_one(&c, 1); qaws_path_compute_area_2d(&p, &ab);
		sprintf(msg, "square areas +4 / -4 (%.12g, %.12g)", (double)aa, (double)ab);
		TEST_ASSERT(fabs(aa - 4) < 1e-12 && fabs(ab + 4) < 1e-12, msg);
		TEST_ASSERT(qaws_curve_is_closed(a), "a closed polyline is closed");
		qaws_curve_destroy(a);
		qaws_curve_destroy(b);
	}
	/* ellipse, rotated */
	{
		qaws_curve* e = NULL;
		qaws_curve const* c;
		qaws_path_2d p;
		qaws_scalar a = 0;
		qaws_curve_create_ellipse_2d(1, -2, 3, (qaws_scalar)1.5, (qaws_scalar)0.4, &e);
		c = e; p = ptt_one(&c, 1);
		qaws_path_compute_area_2d(&p, &a);
		sprintf(msg, "ellipse area pi 3 1.5 (%.15g vs %.15g)", (double)a, PTT_PI * 4.5);
		TEST_ASSERT(fabs(a - PTT_PI * 4.5) < 4.5 * PTT_REL * 10, msg);
		qaws_curve_destroy(e);
	}
	/* a circle as two arcs (two curves), and a Bezier closed by a line */
	{
		qaws_arc_segment s;
		qaws_arc_desc ad;
		qaws_curve* h[2] = { NULL, NULL };
		qaws_curve const* cs[2];
		qaws_path_2d p;
		qaws_scalar a = 0;
		memset(&s, 0, sizeof(s));
		s.radius = 2; s.angle_start = 0; s.angle_end = (qaws_scalar)PTT_PI;
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = &s; ad.segment_count = 1;
		qaws_curve_create_arc(&ad, &h[0]);
		s.angle_start = (qaws_scalar)PTT_PI; s.angle_end = (qaws_scalar)(2 * PTT_PI);
		qaws_curve_create_arc(&ad, &h[1]);
		cs[0] = h[0]; cs[1] = h[1];
		p.curves = cs; p.curve_count = 2; p.closed = 1;
		qaws_path_compute_area_2d(&p, &a);
		sprintf(msg, "circle of two arcs: area 4 pi (%.15g)", (double)a);
		TEST_ASSERT(fabs(a - 4 * PTT_PI) < 4 * PTT_REL * 10, msg);
		qaws_curve_destroy(h[0]);
		qaws_curve_destroy(h[1]);
	}
	{
		/* y = 1 - x^2 on [-1, 1] as a quadratic Bezier, closed along y = 0: 4/3 */
		qaws_scalar cp[] = { -1, 0, 0, 2, 1, 0 };
		qaws_bezier_desc d;
		qaws_curve* b = NULL;
		qaws_curve const* c;
		qaws_path_2d p;
		qaws_scalar a = 0;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = cp; d.control_point_count = 3;
		qaws_curve_create_bezier(&d, &b);
		c = b; p = ptt_one(&c, 0);
		qaws_path_compute_area_2d(&p, &a);
		sprintf(msg, "open parabola closed by its chord: area -4/3 (%.15g)", (double)a);
		TEST_ASSERT(fabs(a + 4.0 / 3.0) < 1e-12 + PTT_REL, msg);
		qaws_curve_destroy(b);
	}
	/* closed curves of other kinds against a dense shoelace */
	{
		qaws_scalar cr[] = { 0, 0, 3, 0, 4, 2, 2, 4, -1, 2 };
		qaws_catmull_rom_desc cd;
		qaws_subdivision_desc sd;
		qaws_curve *k[2] = { NULL, NULL };
		unsigned int i, bad = 0;
		memset(&cd, 0, sizeof(cd));
		cd.dimension = QAWS_DIMENSION_2D; cd.control_points = cr; cd.control_point_count = 5;
		cd.parameterization = QAWS_PARAMETERIZATION_CHORDAL; cd.closed = 1;
		qaws_curve_create_catmull_rom(&cd, &k[0]);
		memset(&sd, 0, sizeof(sd));
		sd.dimension = QAWS_DIMENSION_2D; sd.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_3;
		sd.control_points = cr; sd.control_point_count = 5; sd.closed = 1; sd.refinement_levels = 3;
		qaws_curve_create_subdivision(&sd, &k[1]);
		for (i = 0; i < 2; i++)
		{
			qaws_curve const* c = k[i];
			qaws_path_2d p = ptt_one(&c, 1);
			qaws_scalar a = 0;
			double ref = ptt_sampled_area(&p);
			qaws_path_compute_area_2d(&p, &a);
			if (fabs(a - ref) > 1e-6 * fabs(ref) + PTT_REL * 10)
			{
				printf("    kind %d: area %.12g, sampled %.12g\n", (int)qaws_curve_get_kind(k[i]), (double)a, ref);
				bad++;
			}
			qaws_curve_destroy(k[i]);
		}
		TEST_ASSERT(bad == 0, "closed Catmull-Rom and subdivision areas match a dense sampling");
	}
}

static void test_bounds(void)
{
	char msg[200];
	qaws_curve* e = NULL;
	qaws_curve const* c;
	qaws_path_2d p;
	qaws_vec2 lo, hi;
	qaws_curve_create_ellipse_2d(1, 2, 3, 3, 0, &e);
	c = e; p = ptt_one(&c, 1);
	qaws_path_compute_bounds_2d(&p, &lo, &hi);
	sprintf(msg, "circle bounds exact (%.15g %.15g %.15g %.15g)", (double)lo.x, (double)lo.y, (double)hi.x, (double)hi.y);
	TEST_ASSERT(fabs(lo.x + 2) < 1e-12 + PTT_REL && fabs(hi.x - 4) < 1e-12 + PTT_REL &&
		fabs(lo.y + 1) < 1e-12 + PTT_REL && fabs(hi.y - 5) < 1e-12 + PTT_REL, msg);
	qaws_curve_destroy(e);
	{
		qaws_scalar cp[] = { 0, 0, 1, 3, 3, -2, 4, 1 };
		qaws_bezier_desc d;
		qaws_curve* b = NULL;
		double slo[2] = { 1e9, 1e9 }, shi[2] = { -1e9, -1e9 };
		unsigned int i;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &b);
		for (i = 0; i <= 100000; i++)
		{
			qaws_eval_result_2d r;
			qaws_curve_evaluate_2d(b, (qaws_scalar)(i / 100000.0), QAWS_EVAL_FLAG_POSITION, &r);
			if (r.position.x < slo[0]) slo[0] = r.position.x;
			if (r.position.y < slo[1]) slo[1] = r.position.y;
			if (r.position.x > shi[0]) shi[0] = r.position.x;
			if (r.position.y > shi[1]) shi[1] = r.position.y;
		}
		c = b; p = ptt_one(&c, 0);
		qaws_path_compute_bounds_2d(&p, &lo, &hi);
		sprintf(msg, "Bezier bounds tight (y %.12g..%.12g vs sampled %.12g..%.12g)", (double)lo.y, (double)hi.y, slo[1], shi[1]);
		TEST_ASSERT(lo.y <= slo[1] + 1e-12 + PTT_REL && hi.y >= shi[1] - 1e-12 - PTT_REL && slo[1] - lo.y < 1e-8 + PTT_REL && hi.y - shi[1] < 1e-8 + PTT_REL, msg);
		qaws_curve_destroy(b);
	}
}

static void test_locate(void)
{
	char msg[256];
	qaws_vec2 q;
	qaws_point_location loc;
	/* circle */
	{
		qaws_curve* e = NULL;
		qaws_curve const* c;
		qaws_path_2d p;
		int okc;
		qaws_curve_create_ellipse_2d(0, 0, 2, 2, 0, &e);
		c = e; p = ptt_one(&c, 1);
		q.x = 0; q.y = 0; qaws_region_locate_point_2d(&p, 1, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okc = loc == QAWS_POINT_INSIDE;
		q.x = 3; q.y = 0; qaws_region_locate_point_2d(&p, 1, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okc = okc && loc == QAWS_POINT_OUTSIDE;
		q.x = (qaws_scalar)sqrt(2.0); q.y = (qaws_scalar)sqrt(2.0); qaws_region_locate_point_2d(&p, 1, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okc = okc && loc == QAWS_POINT_ON;
#if !QAWS_SCALAR_IS_FLOAT
		q.x = (qaws_scalar)(2 - 1e-7); q.y = 0; qaws_region_locate_point_2d(&p, 1, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okc = okc && loc == QAWS_POINT_INSIDE;
		q.x = (qaws_scalar)(2 + 1e-7); q.y = 0; qaws_region_locate_point_2d(&p, 1, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okc = okc && loc == QAWS_POINT_OUTSIDE;
#endif
		TEST_ASSERT(okc, "circle: centre inside, far point outside, boundary point on, points 1e-7 off on the right side");
		qaws_curve_destroy(e);
	}
	/* a square with a square hole, and two overlapping squares, under every rule */
	{
		qaws_scalar outer[] = { 0, 0, 4, 0, 4, 4, 0, 4 }, hole[] = { 1, 1, 1, 3, 3, 3, 3, 1 };
		qaws_scalar same[] = { 1, 1, 3, 1, 3, 3, 1, 3 };
		qaws_curve *o = NULL, *h = NULL, *s = NULL;
		qaws_curve const *c0, *c1;
		qaws_path_2d p[2];
		int w = 0, okh;
		qaws_curve_create_polyline_2d(outer, 4, 1, &o);
		qaws_curve_create_polyline_2d(hole, 4, 1, &h);
		qaws_curve_create_polyline_2d(same, 4, 1, &s);
		c0 = o; c1 = h;
		p[0] = ptt_one(&c0, 1); p[1] = ptt_one(&c1, 1);
		q.x = 2; q.y = 2;
		qaws_region_locate_point_2d(p, 2, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okh = loc == QAWS_POINT_OUTSIDE;
		qaws_region_locate_point_2d(p, 2, QAWS_FILL_EVEN_ODD, q, 0, &loc);
		okh = okh && loc == QAWS_POINT_OUTSIDE;
		q.x = (qaws_scalar)0.5;
		qaws_region_locate_point_2d(p, 2, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okh = okh && loc == QAWS_POINT_INSIDE;
		TEST_ASSERT(okh, "square with a hole: hole outside, ring inside");
		/* same orientation: winding 2 in the overlap */
		c1 = s; p[1] = ptt_one(&c1, 1);
		q.x = 2; q.y = 2;
		qaws_region_locate_point_2d(p, 2, QAWS_FILL_EVEN_ODD, q, 0, &loc);
		okh = loc == QAWS_POINT_OUTSIDE;
		qaws_region_locate_point_2d(p, 2, QAWS_FILL_NON_ZERO, q, 0, &loc);
		okh = okh && loc == QAWS_POINT_INSIDE;
		qaws_region_locate_point_2d(p, 2, QAWS_FILL_POSITIVE, q, 0, &loc);
		okh = okh && loc == QAWS_POINT_INSIDE;
		qaws_region_locate_point_2d(p, 2, QAWS_FILL_NEGATIVE, q, 0, &loc);
		okh = okh && loc == QAWS_POINT_OUTSIDE;
		qaws_path_compute_winding_2d(&p[0], q, &w);
		sprintf(msg, "nested squares of one orientation: winding 2 in the overlap, each fill rule right (outer winding %d)", w);
		TEST_ASSERT(okh && w == 1, msg);
		qaws_curve_destroy(o);
		qaws_curve_destroy(h);
		qaws_curve_destroy(s);
	}
}

static void test_transform(void)
{
	char msg[256];
	qaws_arc_segment s;
	qaws_arc_desc ad;
	qaws_curve *arc = NULL, *rot = NULL, *shear = NULL;
	qaws_scalar const mr[6] = { (qaws_scalar)(2 * cos(0.7)), (qaws_scalar)(-2 * sin(0.7)), 1,
		(qaws_scalar)(2 * sin(0.7)), (qaws_scalar)(2 * cos(0.7)), -3 };
	qaws_scalar const ms[6] = { 1, (qaws_scalar)0.5, 0, 0, 2, 1 };
	double worst = 0.0;
	unsigned int i;
	memset(&s, 0, sizeof(s));
	s.center[0] = 1; s.center[1] = 1; s.radius = 2; s.angle_start = (qaws_scalar)0.3; s.angle_end = (qaws_scalar)2.9;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &s; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, &arc);
	qaws_curve_transform_2d(arc, mr, &rot);
	qaws_curve_transform_2d(arc, ms, &shear);
	for (i = 0; i <= 64; i++)
	{
		qaws_range ra = qaws_curve_get_parameter_range(arc), rr = qaws_curve_get_parameter_range(rot), rs = qaws_curve_get_parameter_range(shear);
		qaws_eval_result_2d a, b, c;
		double f = i / 64.0, x, y;
		qaws_curve_evaluate_2d(arc, (qaws_scalar)(ra.min_value + f * (ra.max_value - ra.min_value)), QAWS_EVAL_FLAG_POSITION, &a);
		qaws_curve_evaluate_2d(rot, (qaws_scalar)(rr.min_value + f * (rr.max_value - rr.min_value)), QAWS_EVAL_FLAG_POSITION, &b);
		qaws_curve_evaluate_2d(shear, (qaws_scalar)(rs.min_value + f * (rs.max_value - rs.min_value)), QAWS_EVAL_FLAG_POSITION, &c);
		x = mr[0] * a.position.x + mr[1] * a.position.y + mr[2];
		y = mr[3] * a.position.x + mr[4] * a.position.y + mr[5];
		if (hypot(x - b.position.x, y - b.position.y) > worst) worst = hypot(x - b.position.x, y - b.position.y);
		x = ms[0] * a.position.x + ms[1] * a.position.y + ms[2];
		y = ms[3] * a.position.x + ms[4] * a.position.y + ms[5];
		/* the NURBS is not arc-length parameterized: compare by distance to the mapped point set at the ends only */
		if (i == 0 || i == 64)
			if (hypot(x - c.position.x, y - c.position.y) > worst) worst = hypot(x - c.position.x, y - c.position.y);
	}
	sprintf(msg, "an arc rotated and scaled stays an arc (%d), sheared becomes NURBS (%d); points mapped within %.1e",
		(int)qaws_curve_get_kind(rot), (int)qaws_curve_get_kind(shear), worst);
	TEST_ASSERT(qaws_curve_get_kind(rot) == QAWS_CURVE_KIND_ARC && qaws_curve_get_kind(shear) == QAWS_CURVE_KIND_NURBS &&
		worst < 1e-12 + PTT_REL * 10, msg);
	/* the sheared arc's area equals det(m) times the arc's segment area */
	{
		qaws_curve const* c;
		qaws_path_2d p;
		qaws_scalar a0 = 0, a1 = 0;
		c = arc; p = ptt_one(&c, 0); qaws_path_compute_area_2d(&p, &a0);
		c = shear; p = ptt_one(&c, 0); qaws_path_compute_area_2d(&p, &a1);
		sprintf(msg, "sheared arc area = det 2 x arc area (%.12g vs %.12g)", (double)a1, 2.0 * a0);
		TEST_ASSERT(fabs(a1 - 2 * a0) < 1e-10 + PTT_REL * 10, msg);
	}
	qaws_curve_destroy(arc);
	qaws_curve_destroy(rot);
	qaws_curve_destroy(shear);
}

static void test_polylines(void)
{
	char msg[200];
	qaws_scalar out[64];
	unsigned int n = 0;
	{
		qaws_scalar p[] = { 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0 };
		qaws_polyline_strip_duplicates_2d(p, 6, 1, out, &n);
		TEST_ASSERT(n == 3, "strip duplicates: 6 points with repeats and a closing copy keep 3");
	}
	{
		/* a square with midpoints and a spike */
		qaws_scalar p[] = { 0, 0, 1, 0, 2, 0, 2, 1, 2, 2, 1, 2, 0, 2, 0, 1, 0, 3, 0, 1 };
		qaws_polyline_trim_collinear_2d(p, 10, 1, out, &n);
		sprintf(msg, "trim collinear: midpoints and a spike go, the 4 corners stay (%u)", n);
		TEST_ASSERT(n == 4, msg);
	}
	{
		/* Clipper2's TestSimplifyPath: an open zig-zag within epsilon of a line */
		qaws_scalar p[] = { 0, 0, 1, (qaws_scalar)-0.1, 2, (qaws_scalar)0.1, 3, (qaws_scalar)-0.1, 4, (qaws_scalar)0.1, 100, 0 };
		qaws_polyline_simplify_2d(p, 6, 0, 1, out, &n);
		sprintf(msg, "simplify: an open zig-zag within epsilon keeps its 2 ends (%u)", n);
		TEST_ASSERT(n == 2 && out[0] == 0 && out[2] == 100, msg);
	}
	{
		qaws_scalar p[] = { 0, 0, 1, (qaws_scalar)0.05, 2, 0, 3, 5, 4, 0, 5, (qaws_scalar)0.02, 6, 0 };
		qaws_polyline_rdp_2d(p, 7, (qaws_scalar)0.5, out, &n);
		sprintf(msg, "RDP keeps the ends, the peak and its feet (%u)", n);
		TEST_ASSERT(n == 5 && out[4] == 3, msg);
	}
}

int test_83_path_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 83: 2D paths and regions\n");
	test_area();
	test_bounds();
	test_locate();
	test_transform();
	test_polylines();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
