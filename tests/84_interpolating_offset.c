/*
 * Test 84: splines solved from points, parallel curves as cubics
 *
 *   - every method passes through its points, open and closed
 *   - continuity at the joins: C2 for the uniform cubic, G1 for Hobby, the curvature's magnitude for kappa
 *   - Hobby through four points of a circle stays close to the circle
 *   - collinear points, duplicates, scale, invalid input
 *   - offsets: the distance to a sampled curve within the tolerance on a cubic, an
 *     arc, a B-spline, an inflected S; left is positive; straight lines stay straight;
 *     a corner between spans is bridged; few segments on a cubic; invalid input
 */

#include "test_common.h"
#include <math.h>

#if QAWS_SCALAR_IS_FLOAT
#define IO_EPS 2e-4
#define IO_TOL 1e-3
#else
#define IO_EPS 1e-9
#define IO_TOL 1e-4
#endif

static qaws_vec2 const io_points[6] = { { 0, 0 }, { 2, 1 }, { 4, 0.5f }, { 5, 3 }, { 3, 4.5f }, { 0.5f, 3 } };

static qaws_curve* io_make(qaws_interpolating_method m, qaws_vec2 const* p, unsigned int n, int closed)
{
	qaws_interpolating_spline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof d);
	d.points = p;
	d.point_count = n;
	d.method = m;
	d.closed = closed;
	return qaws_curve_create_interpolating_spline(&c, &d) == QAWS_STATUS_OK ? c : NULL;
}

static qaws_eval_result_2d io_at(qaws_curve const* c, double t)
{
	qaws_eval_result_2d r;
	memset(&r, 0, sizeof r);
	qaws_curve_evaluate_2d(&r, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
	return r;
}

/* the distance from p to the curve, by dense sampling refined by the segment between samples */
static double io_distance(qaws_curve const* c, double x, double y, unsigned int samples)
{
	qaws_range range = qaws_curve_get_parameter_range(c);
	double best = HUGE_VAL, px = 0, py = 0;
	unsigned int i;
	for (i = 0; i <= samples; i++)
	{
		qaws_eval_result_2d r = io_at(c, (double)range.min_value + (double)(range.max_value - range.min_value) * i / samples);
		double qx = (double)r.position.x, qy = (double)r.position.y;
		if (i)
		{
			double dx = qx - px, dy = qy - py, l = dx * dx + dy * dy, s = l > 0 ? ((x - px) * dx + (y - py) * dy) / l : 0, ex, ey, e;
			s = s < 0 ? 0 : s > 1 ? 1 : s;
			ex = px + s * dx - x;
			ey = py + s * dy - y;
			e = sqrt(ex * ex + ey * ey);
			best = e < best ? e : best;
		}
		px = qx;
		py = qy;
	}
	return best;
}

static double io_curvature(qaws_eval_result_2d const* r)
{
	double dx = (double)r->d1.x, dy = (double)r->d1.y;
	return (dx * (double)r->d2.y - dy * (double)r->d2.x) / pow(dx * dx + dy * dy, 1.5);
}

static void test_through_points(void)
{
	static qaws_interpolating_method const methods[3] = { QAWS_INTERPOLATING_CUBIC_C2, QAWS_INTERPOLATING_HOBBY, QAWS_INTERPOLATING_KAPPA };
	int m, closed;
	unsigned int i;
	for (m = 0; m < 3; m++)
		for (closed = 0; closed < 2; closed++)
		{
			qaws_curve* c = io_make(methods[m], io_points, 6, closed);
			double worst = 0;
			char msg[96];
			sprintf(msg, "method %d closed %d is made", m, closed);
			TEST_ASSERT(c != NULL, msg);
			if (!c)
				continue;
			for (i = 0; i < 6; i++)
			{
				double d = io_distance(c, (double)io_points[i].x, (double)io_points[i].y, 6000);
				worst = d > worst ? d : worst;
			}
			sprintf(msg, "method %d closed %d passes through its points (%.2e)", m, closed, worst);
			TEST_ASSERT(worst < 1e-3, msg);
			if (closed)
			{
				qaws_range r = qaws_curve_get_parameter_range(c);
				qaws_eval_result_2d a = io_at(c, (double)r.min_value), b = io_at(c, (double)r.max_value);
				sprintf(msg, "method %d closes", m);
				TEST_ASSERT(fabs((double)a.position.x - (double)b.position.x) < IO_EPS && fabs((double)a.position.y - (double)b.position.y) < IO_EPS, msg);
			}
			qaws_curve_destroy(c);
		}
}

static void test_continuity(void)
{
	double const h = QAWS_SCALAR_IS_FLOAT ? 1e-3 : 1e-7;
	qaws_curve* c2 = io_make(QAWS_INTERPOLATING_CUBIC_C2, io_points, 6, 0);
	qaws_curve* hobby = io_make(QAWS_INTERPOLATING_HOBBY, io_points, 6, 1);
	qaws_curve* kappa = io_make(QAWS_INTERPOLATING_KAPPA, io_points, 6, 1);
	double wc = 0, wh = 0, wk = 0;
	unsigned int i;
	for (i = 1; i < 5 && c2; i++)
	{
		qaws_eval_result_2d a = io_at(c2, i - h), b = io_at(c2, i + h);
		double d = fabs((double)a.d1.x - (double)b.d1.x) + fabs((double)a.d1.y - (double)b.d1.y) + fabs((double)a.d2.x - (double)b.d2.x) + fabs((double)a.d2.y - (double)b.d2.y);
		wc = d > wc ? d : wc;
	}
	TEST_ASSERT(c2 && wc < (QAWS_SCALAR_IS_FLOAT ? 0.2 : 1e-4), "the uniform cubic is C2 at its points");
	for (i = 1; i < 6 && hobby; i++)
	{
		qaws_eval_result_2d a = io_at(hobby, i - h), b = io_at(hobby, i + h);
		double cross = ((double)a.d1.x * (double)b.d1.y - (double)a.d1.y * (double)b.d1.x)
			/ sqrt(((double)a.d1.x * (double)a.d1.x + (double)a.d1.y * (double)a.d1.y) * ((double)b.d1.x * (double)b.d1.x + (double)b.d1.y * (double)b.d1.y));
		wh = fabs(cross) > wh ? fabs(cross) : wh;
	}
	TEST_ASSERT(hobby && wh < (QAWS_SCALAR_IS_FLOAT ? 0.05 : 1e-5), "Hobby's spline is G1 at its points");
	for (i = 1; i < 6 && kappa; i++)
	{
		qaws_eval_result_2d a = io_at(kappa, i - h), b = io_at(kappa, i + h);
		double ka = fabs(io_curvature(&a)), kb = fabs(io_curvature(&b)), rel = fabs(ka - kb) / (ka + kb + 1e-3);
		wk = rel > wk ? rel : wk;
	}
	TEST_ASSERT(kappa && wk < (QAWS_SCALAR_IS_FLOAT ? 0.1 : 2e-3), "the kappa spline's curvature magnitude matches at its joins");
	qaws_curve_destroy(c2);
	qaws_curve_destroy(hobby);
	qaws_curve_destroy(kappa);
}

static void test_circle(void)
{
	static qaws_vec2 const p[4] = { { 2, 0 }, { 0, 2 }, { -2, 0 }, { 0, -2 } };
	qaws_curve* c = io_make(QAWS_INTERPOLATING_HOBBY, p, 4, 1);
	double worst = 0;
	unsigned int i;
	TEST_ASSERT(c != NULL, "Hobby through a circle's points");
	for (i = 0; c && i <= 400; i++)
	{
		qaws_eval_result_2d r = io_at(c, 4.0 * i / 400);
		double e = fabs(sqrt((double)r.position.x * (double)r.position.x + (double)r.position.y * (double)r.position.y) - 2);
		worst = e > worst ? e : worst;
	}
	TEST_ASSERT(worst < 2 * 1.5e-3, "Hobby's spline through four points of a circle is close to it");
	qaws_curve_destroy(c);
}

static void test_edge_cases(void)
{
	static qaws_vec2 const line[4] = { { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 3 } };
	static qaws_vec2 const dup[5] = { { 0, 0 }, { 0, 0 }, { 1, 1 }, { 2, 0 }, { 2, 0 } };
	static qaws_vec2 const two[2] = { { 0, 0 }, { 4, 2 } };
	static qaws_vec2 const three[3] = { { 0, 0 }, { 1, 2 }, { 3, 0 } };
	static qaws_vec2 const same[3] = { { 1, 1 }, { 1, 1 }, { 1, 1 } };
	static qaws_interpolating_method const methods[3] = { QAWS_INTERPOLATING_CUBIC_C2, QAWS_INTERPOLATING_HOBBY, QAWS_INTERPOLATING_KAPPA };
	qaws_interpolating_spline_desc d;
	qaws_curve* c = NULL;
	qaws_vec2 bad[3];
	int m;
	for (m = 0; m < 3; m++)
	{
		qaws_curve* a = io_make(methods[m], line, 4, 0);
		qaws_curve* b = io_make(methods[m], dup, 5, 0);
		qaws_curve* e = io_make(methods[m], two, 2, 0);
		qaws_curve* f = io_make(methods[m], three, 3, 0);
		char msg[80];
		sprintf(msg, "method %d: collinear points, duplicates, two and three points", m);
		TEST_ASSERT(a && b && e && f, msg);
		if (a)
		{
			qaws_eval_result_2d r = io_at(a, 0.5 * (double)qaws_curve_get_parameter_range(a).max_value);
			sprintf(msg, "method %d: collinear points stay on their line", m);
			TEST_ASSERT(fabs((double)r.position.x - (double)r.position.y) < IO_EPS, msg);
		}
		if (e)
		{
			qaws_eval_result_2d r = io_at(e, 0.5);
			sprintf(msg, "method %d: two points are a straight segment", m);
			TEST_ASSERT(fabs((double)r.position.x - 2) < IO_EPS && fabs((double)r.position.y - 1) < IO_EPS, msg);
		}
		qaws_curve_destroy(a);
		qaws_curve_destroy(b);
		qaws_curve_destroy(e);
		qaws_curve_destroy(f);
		sprintf(msg, "method %d: one distinct point is refused", m);
		memset(&d, 0, sizeof d);
		d.points = same;
		d.point_count = 3;
		d.method = methods[m];
		TEST_ASSERT(qaws_curve_create_interpolating_spline(&c, &d) == QAWS_STATUS_DEGENERATE_CURVE, msg);
	}
	memset(&d, 0, sizeof d);
	d.points = two;
	d.point_count = 2;
	d.closed = 1;
	TEST_ASSERT(qaws_curve_create_interpolating_spline(&c, &d) == QAWS_STATUS_DEGENERATE_CURVE, "a closed curve needs three distinct points");
	d.closed = 0;
	d.method = (qaws_interpolating_method)99;
	TEST_ASSERT(qaws_curve_create_interpolating_spline(&c, &d) == QAWS_STATUS_INVALID_ARGUMENT, "an unknown method is refused");
	TEST_ASSERT(qaws_curve_create_interpolating_spline(&c, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "no description is refused");
	bad[0] = three[0];
	bad[1] = three[1];
	bad[2] = three[2];
	bad[1].x = (qaws_scalar)(HUGE_VAL - HUGE_VAL);
	d.method = QAWS_INTERPOLATING_HOBBY;
	d.points = bad;
	d.point_count = 3;
	TEST_ASSERT(qaws_curve_create_interpolating_spline(&c, &d) == QAWS_STATUS_INVALID_ARGUMENT, "NaN is refused");
	/* far from the origin and large: the shape scales with the points */
	{
		qaws_vec2 scaled[6], moved;
		unsigned int i;
		qaws_curve *u = io_make(QAWS_INTERPOLATING_KAPPA, io_points, 6, 1), *v;
		double worst = 0;
		double const s = QAWS_SCALAR_IS_FLOAT ? 100.0 : 1e5;
		for (i = 0; i < 6; i++)
		{
			moved.x = (qaws_scalar)((double)io_points[i].x * s + 3 * s);
			moved.y = (qaws_scalar)((double)io_points[i].y * s - 2 * s);
			scaled[i] = moved;
		}
		v = io_make(QAWS_INTERPOLATING_KAPPA, scaled, 6, 1);
		for (i = 0; u && v && i <= 60; i++)
		{
			qaws_eval_result_2d a = io_at(u, 6.0 * i / 60), b = io_at(v, 6.0 * i / 60);
			double e = fabs(((double)b.position.x - 3 * s) / s - (double)a.position.x) + fabs(((double)b.position.y + 2 * s) / s - (double)a.position.y);
			worst = e > worst ? e : worst;
		}
		TEST_ASSERT(u && v && worst < (QAWS_SCALAR_IS_FLOAT ? 5e-2 : 1e-7), "the solve does not depend on the points' scale");
		qaws_curve_destroy(u);
		qaws_curve_destroy(v);
	}
}

static qaws_curve* io_bezier(double const* xy)
{
	qaws_vec2 cp[4];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 4; i++)
	{
		cp[i].x = (qaws_scalar)xy[2 * i];
		cp[i].y = (qaws_scalar)xy[2 * i + 1];
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = cp;
	d.control_point_count = 4;
	return qaws_curve_create_bezier(&c, &d) == QAWS_STATUS_OK ? c : NULL;
}

/* every point of the offset is the distance from the input, away from cusp arcs */
static double io_offset_error(qaws_curve const* in, qaws_curve const* off, double distance)
{
	qaws_range r = qaws_curve_get_parameter_range(off);
	double worst = 0;
	unsigned int i;
	for (i = 0; i <= 400; i++)
	{
		qaws_eval_result_2d p = io_at(off, (double)r.min_value + (double)(r.max_value - r.min_value) * i / 400);
		double e = fabs(io_distance(in, (double)p.position.x, (double)p.position.y, 1500) - fabs(distance));
		worst = e > worst ? e : worst;
	}
	return worst;
}

static void test_offset(void)
{
	static double const smooth[8] = { 0, 0, 2, 3, 5, 3, 7, 0 };
	static double const s_curve[8] = { 0, 0, 3, 4, 3, -4, 6, 0 };
	static double const straight[8] = { 0, 0, 1, 0, 2, 0, 3, 0 };
	static double const wide[8] = { 0, 0, 0, 8, 8, 8, 8, 0 };
	double const* const shapes[3] = { smooth, s_curve, wide };
	static char const* const names[3] = { "a smooth cubic", "an inflected cubic", "a wide arch" };
	static double const distances[4] = { 0.3, -0.3, 0.8, -0.8 };
	int s, d;
	for (s = 0; s < 3; s++)
		for (d = 0; d < 4; d++)
		{
			qaws_curve *in = io_bezier(shapes[s]), *off = NULL;
			qaws_status st = qaws_curve_offset_cubic_2d(&off, in, (qaws_scalar)distances[d], (qaws_scalar)IO_TOL);
			char msg[96];
			sprintf(msg, "%s offset by %.1f is made", names[s], distances[d]);
			TEST_ASSERT(st == QAWS_STATUS_OK && off != NULL, msg);
			if (off)
			{
				double e = io_offset_error(in, off, distances[d]);
				sprintf(msg, "%s offset by %.1f is within the tolerance (%.2e)", names[s], distances[d], e);
				TEST_ASSERT(e < 3 * IO_TOL, msg);
				qaws_curve_destroy(off);
			}
			qaws_curve_destroy(in);
		}
	/* left is positive, like qaws_curve_offset_2d; a straight line stays one */
	{
		qaws_curve *in = io_bezier(straight), *off = NULL;
		qaws_eval_result_2d a, b;
		TEST_ASSERT(qaws_curve_offset_cubic_2d(&off, in, (qaws_scalar)0.5, (qaws_scalar)IO_TOL) == QAWS_STATUS_OK && off, "a straight line's offset is made");
		a = io_at(off, 0);
		b = io_at(off, (double)qaws_curve_get_parameter_range(off).max_value);
		TEST_ASSERT(fabs((double)a.position.y - 0.5) < IO_EPS && fabs((double)b.position.y - 0.5) < IO_EPS && fabs((double)a.position.x) < IO_EPS && fabs((double)b.position.x - 3) < IO_EPS,
			"positive distance moves a line travelling along +x to +y");
		TEST_ASSERT(qaws_curve_get_span_count(off) == 1, "a straight line's offset is one segment");
		qaws_curve_destroy(off);
		qaws_curve_destroy(in);
	}
	/* a circle made of an arc stays a circle of the offset radius */
	{
		qaws_arc_desc ad;
		qaws_arc_segment seg;
		qaws_curve *arc = NULL, *off = NULL;
		unsigned int i;
		double worst = 0;
		memset(&seg, 0, sizeof seg);
		seg.center[0] = 1;
		seg.center[1] = 2;
		seg.radius = 3;
		seg.angle_start = 0;
		seg.angle_end = (qaws_scalar)3.0;
		memset(&ad, 0, sizeof ad);
		ad.dimension = QAWS_DIMENSION_2D;
		ad.segments = &seg;
		ad.segment_count = 1;
		TEST_ASSERT(qaws_curve_create_arc(&arc, &ad) == QAWS_STATUS_OK, "an arc to offset");
		TEST_ASSERT(qaws_curve_offset_cubic_2d(&off, arc, (qaws_scalar)-0.5, (qaws_scalar)IO_TOL) == QAWS_STATUS_OK && off, "an arc's offset is made");
		for (i = 0; off && i <= 200; i++)
		{
			qaws_range r = qaws_curve_get_parameter_range(off);
			qaws_eval_result_2d p = io_at(off, (double)r.min_value + (double)(r.max_value - r.min_value) * i / 200);
			double dx = (double)p.position.x - 1, dy = (double)p.position.y - 2, e = fabs(sqrt(dx * dx + dy * dy) - 3.5);
			worst = e > worst ? e : worst;
		}
		TEST_ASSERT(worst < 3 * IO_TOL, "an arc's offset to the right of its travel has the larger radius");
		qaws_curve_destroy(off);
		qaws_curve_destroy(arc);
	}
	/* spans of a B-spline, and a polyline's corner bridged */
	{
		static qaws_vec2 const cp[7] = { { 0, 0 }, { 1, 2 }, { 3, 3 }, { 5, 1 }, { 7, 2 }, { 8, 4 }, { 10, 3 } };
		static qaws_vec2 const poly[3] = { { 0, 0 }, { 2, 0 }, { 2, 2 } };
		qaws_bspline_desc bd;
		qaws_curve *in = NULL, *off = NULL, *pl = NULL;
		memset(&bd, 0, sizeof bd);
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 3;
		bd.control_points = cp;
		bd.control_point_count = 7;
		TEST_ASSERT(qaws_curve_create_bspline(&in, &bd) == QAWS_STATUS_OK, "a B-spline to offset");
		TEST_ASSERT(qaws_curve_offset_cubic_2d(&off, in, (qaws_scalar)0.4, (qaws_scalar)IO_TOL) == QAWS_STATUS_OK && off, "a B-spline's offset is made");
		if (off)
		{
			double e = io_offset_error(in, off, 0.4);
			TEST_ASSERT(e < 3 * IO_TOL, "a B-spline's offset is within the tolerance");
			TEST_ASSERT(qaws_curve_get_span_count(off) < 400, "its spans need no refitting");
			qaws_curve_destroy(off);
			off = NULL;
		}
		bd.degree = 1;
		bd.control_points = poly;
		bd.control_point_count = 3;
		TEST_ASSERT(qaws_curve_create_bspline(&pl, &bd) == QAWS_STATUS_OK, "a polyline to offset");
		TEST_ASSERT(qaws_curve_offset_cubic_2d(&off, pl, (qaws_scalar)-0.5, (qaws_scalar)IO_TOL) == QAWS_STATUS_OK && off, "a polyline's offset is made");
		if (off)
		{
			qaws_range r = qaws_curve_get_parameter_range(off);
			qaws_eval_result_2d a = io_at(off, (double)r.min_value), b = io_at(off, (double)r.max_value);
			TEST_ASSERT(fabs((double)a.position.y + 0.5) < IO_EPS && fabs((double)b.position.x - 2.5) < IO_EPS, "a corner's two offsets are both there");
			TEST_ASSERT(qaws_curve_get_span_count(off) == 3, "the corner between them is bridged by one line");
			qaws_curve_destroy(off);
		}
		qaws_curve_destroy(pl);
		qaws_curve_destroy(in);
	}
	/* a cusp: the arc about it keeps the chain connected */
	{
		static double const cusp[8] = { 0, 0, 2, 2, 0, 2, 2, 0 };
		qaws_curve *in = io_bezier(cusp), *off = NULL;
		TEST_ASSERT(qaws_curve_offset_cubic_2d(&off, in, (qaws_scalar)0.2, (qaws_scalar)IO_TOL) == QAWS_STATUS_OK && off, "a cusped cubic's offset is made");
		if (off)
		{
			unsigned int i, n = qaws_curve_get_span_count(off), gaps = 0;
			for (i = 1; i < n; i++)
			{
				qaws_eval_result_2d a = io_at(off, i - 1e-6), b = io_at(off, i + 1e-6);
				gaps += fabs((double)a.position.x - (double)b.position.x) + fabs((double)a.position.y - (double)b.position.y) > 1e-2;
			}
			TEST_ASSERT(gaps == 0, "the offset of a cusp has no gaps");
			qaws_curve_destroy(off);
		}
		qaws_curve_destroy(in);
	}
	/* refusals */
	{
		qaws_curve *in = io_bezier(smooth), *off = NULL;
		TEST_ASSERT(qaws_curve_offset_cubic_2d(NULL, in, 1, 1) == QAWS_STATUS_INVALID_ARGUMENT, "no output is refused");
		TEST_ASSERT(qaws_curve_offset_cubic_2d(&off, NULL, 1, 1) == QAWS_STATUS_INVALID_ARGUMENT, "no curve is refused");
		TEST_ASSERT(qaws_curve_offset_cubic_2d(&off, in, 1, 0) == QAWS_STATUS_INVALID_ARGUMENT, "a zero tolerance is refused");
		qaws_curve_destroy(in);
	}
}

int test_84_interpolating_offset_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 84: splines solved from points, parallel curves as cubics\n");
	test_through_points();
	test_continuity();
	test_circle();
	test_edge_cases();
	test_offset();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
