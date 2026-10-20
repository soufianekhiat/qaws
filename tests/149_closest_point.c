/*
 * Test 149: certified closest points (qaws_curve_find_closest_point_certified_2d / 3d)
 *
 *   - a circle of radius 5 as a cubic B-spline with a narrow spike towards
 *     its center, between the 32 samples the search starts from: the circle
 *     is a ridge of equal distances where sampling and Newton stop; the
 *     certified search finds the spike, below a dense sampling's minimum;
 *   - a rational circle: the exact distance, a tight gap;
 *   - the circle's center (every point closest): within the work bound;
 *   - a borrowed view of the spiked curve: the same answer;
 *   - a 3D Bezier curve against dense sampling at several points;
 *   - other families: unsupported; qaws_curve_find_closest_parameter_* still answers.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define CP_N 200

static double cp_dense_2d(qaws_curve const* c, qaws_vec2 q, unsigned int n)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double best = 1e300;
	unsigned int i;
	for (i = 0; i <= n; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(&e, c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / n, QAWS_EVAL_FLAG_POSITION);
		best = fmin(best, hypot(e.position.x - q.x, e.position.y - q.y));
	}
	return best;
}

/* the certified gap's contract: within 1e-9 of the extent in double, QAWS_EPSILON (1e-6) of it in float; the tests allow
   ten times the double one */
#define CP_GAP (sizeof(qaws_scalar) == sizeof(double) ? 1e-8 : 1e-6)

/* the rounding of a distance between evaluated points of a curve of the given extent: the evaluation and the sum of
   squares err by a few eps of the coordinates (the double bounds are kept as they were, below this) */
static double cp_round(double extent, double dbl)
{
	return fmax(dbl, 8 * TEST_EPS * extent);
}

static void cp_spiked_points(double* p)
{
	unsigned int i;
	for (i = 0; i < CP_N; i++)
	{
		/* point 67 pulled to radius 0.5: its support, parameters [64, 68] of [0, 197], lies between samples 10 and 11 of 32 */
		double r = i == 67 ? 0.5 : 5.0, a = 2 * M_PI * i / (CP_N - 1);
		p[2 * i] = r * cos(a);
		p[2 * i + 1] = r * sin(a);
	}
}

static void test_spike(void)
{
	static double p[CP_N * 2];
	static qaws_scalar ps[CP_N * 2];   /* the control points as the scalar (the view below reads the doubles) */
	unsigned int k;
	qaws_vec2 q = { 0, 0 };
	qaws_bspline_desc d;
	qaws_curve_view v;
	qaws_curve *c = NULL, *b = NULL;
	qaws_scalar t = -1, dist = -1, gap = -1, t2 = -1, dist2 = -1, gap2 = -1, ts = -1;
	double dense;
	char msg[240];
	cp_spiked_points(p);
	for (k = 0; k < CP_N * 2; k++)
		ps[k] = (qaws_scalar)p[k];
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = ps; d.control_point_count = CP_N;
	TEST_ASSERT(qaws_curve_create_bspline(&c, &d) == QAWS_STATUS_OK, "spiked circle created");
	/* sampling and Newton alone, from the 32 samples: they stop on the circle */
	{
		qaws_range r = qaws_curve_get_parameter_range(c);
		double best = 1e300;
		unsigned int i;
		for (i = 0; i < 32; i++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(&e, c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 31, QAWS_EVAL_FLAG_POSITION);
			best = fmin(best, hypot(e.position.x, e.position.y));
		}
		snprintf(msg, sizeof msg, "the 32 samples see no spike (nearest %.3f)", best);
		printf("  %s\n", msg);
		TEST_ASSERT(best > 4.9, msg);
	}
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(&t, &dist, &gap, c, q) == QAWS_STATUS_OK, "certified search runs");
	dense = cp_dense_2d(c, q, 400000);
	snprintf(msg, sizeof msg, "the spike found: distance %.12f at t = %.6f, gap %.1e; dense sampling %.12f", (double)dist, (double)t, (double)gap, dense);
	printf("  %s\n", msg);
	TEST_ASSERT(dist < 4 && dist <= dense + cp_round(5, 1e-12) && t > 64 && t < 68, msg);
	TEST_ASSERT(gap >= 0 && gap <= CP_GAP * 10, "gap within the tolerance");
	TEST_ASSERT(qaws_curve_find_closest_parameter_2d(&ts, c, q) == QAWS_STATUS_OK && ts == t, "find_closest_parameter gives the same");

	/* the same curve borrowed */
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_2D;
	v.degree = 3;
	v.points.data = p; v.points.type = QAWS_ELEMENT_F64; v.points.components = 2; v.points.count = CP_N; v.points.stride = 2 * sizeof(double);
	TEST_ASSERT(qaws_curve_create_borrowed(&b, &v, NULL) == QAWS_STATUS_OK, "borrowed view created");
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(&t2, &dist2, &gap2, b, q) == QAWS_STATUS_OK, "certified on the view");
	snprintf(msg, sizeof msg, "the view: the same point (%.1e in t, %.1e in distance)", fabs((double)t2 - t), fabs((double)dist2 - dist));
	TEST_ASSERT(fabs((double)t2 - t) <= 1e-10 && fabs((double)dist2 - dist) <= 1e-12, msg);
	qaws_curve_destroy(b);
	qaws_curve_destroy(c);
}

static void test_circle(void)
{
	/* the nine-point rational circle of radius 2 */
	double const s = 0.5 * sqrt(2.0);
	qaws_vec2 p[9] = { { 2, 0 }, { 2, 2 }, { 0, 2 }, { -2, 2 }, { -2, 0 }, { -2, -2 }, { 0, -2 }, { 2, -2 }, { 2, 0 } };
	qaws_scalar w[9] = { 1, (qaws_scalar)s, 1, (qaws_scalar)s, 1, (qaws_scalar)s, 1, (qaws_scalar)s, 1 };
	qaws_scalar knots[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4 };
	qaws_vec2 q = { 3, 1 }, center = { 0, 0 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	qaws_scalar t = -1, dist = -1, gap = -1;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = p; d.control_point_count = 9;
	d.knots = knots; d.knot_count = 12; d.weights = w; d.weight_count = 9;
	TEST_ASSERT(qaws_curve_create_nurbs(&c, &d) == QAWS_STATUS_OK, "rational circle created");
	qaws_curve_find_closest_point_certified_2d(&t, &dist, &gap, c, q);
	snprintf(msg, sizeof msg, "rational circle: distance sqrt(10) - 2 (error %.1e, gap %.1e)", fabs(dist - (sqrt(10.0) - 2)), (double)gap);
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(dist - (sqrt(10.0) - 2)) <= cp_round(2, 1e-13) && gap <= CP_GAP, msg);
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(&t, &dist, &gap, c, center) == QAWS_STATUS_OK, "the center: it terminates");
	snprintf(msg, sizeof msg, "the center: distance 2 (error %.1e), gap %.2e within the work bound", fabs(dist - 2), (double)gap);
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(dist - 2) <= cp_round(2, 1e-13) && gap >= 0 && gap <= 2, msg);
	qaws_curve_destroy(c);
}

static void test_bezier_3d(void)
{
	qaws_vec3 p[6] = { { 0, 0, 0 }, { 1, 3, -1 }, { 3, -2, 2 }, { 4, 3, 0 }, { 6, -1, 1 }, { 7, 1, -2 } };
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	double worst = 0, gap_max = 0;
	unsigned int k, i;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D; d.degree = 5; d.control_points = p; d.control_point_count = 6;
	qaws_curve_create_bezier(&c, &d);
	for (k = 0; k < 24; k++)
	{
		qaws_vec3 q;
		qaws_scalar t, dist, gap;
		double dense = 1e300;
		q.x = (qaws_scalar)(3.5 + 4 * cos(1.7 * k)); q.y = (qaws_scalar)(2.5 * sin(2.3 * k)); q.z = (qaws_scalar)(2 * cos(0.9 * k + 1));
		qaws_curve_find_closest_point_certified_3d(&t, &dist, &gap, c, q);
		for (i = 0; i <= 100000; i++)
		{
			qaws_eval_result_3d e;
			qaws_curve_evaluate_3d(&e, c, (qaws_scalar)i / 100000, QAWS_EVAL_FLAG_POSITION);
			dense = fmin(dense, sqrt((e.position.x - q.x) * (e.position.x - q.x) + (e.position.y - q.y) * (e.position.y - q.y)
				+ (e.position.z - q.z) * (e.position.z - q.z)));
		}
		/* the dense minimum is no lower than the certified distance, and the certified one is within its gap of the truth */
		worst = fmax(worst, dist - dense);
		gap_max = fmax(gap_max, gap);
	}
	snprintf(msg, sizeof msg, "3D quintic, 24 points: never above dense sampling (worst %.1e), gaps up to %.1e", worst, gap_max);
	printf("  %s\n", msg);
	TEST_ASSERT(worst <= cp_round(7, 1e-12) && gap_max <= CP_GAP * 10, msg);
	qaws_curve_destroy(c);
}

static void test_other_families(void)
{
	qaws_vec2 pts[2] = { { 0, 0 }, { 4, 0 } }, der[2] = { { 0, 4 }, { 0, -4 } }, q = { 2, 3 };
	qaws_hermite_desc d;
	qaws_curve* c = NULL;
	qaws_scalar t = -1, dist, gap;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.points = pts; d.derivatives = der; d.point_count = 2; d.derivative_count = 2;
	TEST_ASSERT(qaws_curve_create_hermite(&c, &d) == QAWS_STATUS_OK, "Hermite curve created");
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(&t, &dist, &gap, c, q) == QAWS_STATUS_UNSUPPORTED_OPERATION, "Hermite: not certified");
	TEST_ASSERT(qaws_curve_find_closest_parameter_2d(&t, c, q) == QAWS_STATUS_OK && fabs(t - 0.5) <= 1e-6, "Hermite: still answered (t = 0.5)");
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(&t, NULL, &gap, c, q) == QAWS_STATUS_INVALID_ARGUMENT, "no output refused");
	qaws_curve_destroy(c);
}

#define CP_S 40

/* a spherical cap of radius 5 around the origin, point (14, 14) pulled to radius 0.5: its support, [11, 15] / 37 in u and
   v, lies between the 8 x 8 grid's samples 2 and 3 */
static void cp_cap_points(qaws_vec3* p)
{
	unsigned int i, j;
	for (i = 0; i < CP_S; i++)
		for (j = 0; j < CP_S; j++)
		{
			double th = 0.3 + 2.5 * i / (CP_S - 1), ph = 3.0 * j / (CP_S - 1), r = i == 14 && j == 14 ? 0.5 : 5.0;
			p[i * CP_S + j].x = (qaws_scalar)(r * sin(th) * cos(ph));
			p[i * CP_S + j].y = (qaws_scalar)(r * sin(th) * sin(ph));
			p[i * CP_S + j].z = (qaws_scalar)(r * cos(th));
		}
}

static void test_surfaces(void)
{
	static qaws_vec3 p[CP_S * CP_S];
	qaws_vec3 q = { 0, 0, 0 }, cl, cl2;
	qaws_surface_bspline_desc d;
	qaws_surface_view v;
	qaws_surface *s = NULL, *b = NULL;
	qaws_scalar u = -1, w = -1, dist = -1, gap = -1, u2, w2, dist2, gap2;
	double grid = 1e300, dense = 1e300;
	unsigned int i, j;
	char msg[260];
	cp_cap_points(p);
	memset(&d, 0, sizeof d);
	d.u_degree = 3; d.v_degree = 3; d.control_points = p; d.u_point_count = CP_S; d.v_point_count = CP_S;
	TEST_ASSERT(qaws_surface_create_bspline(&s, &d) == QAWS_STATUS_OK, "spiked cap created");
	for (i = 0; i <= 600; i++)
		for (j = 0; j <= 600; j++)
		{
			qaws_surface_eval_result e;
			double r;
			qaws_surface_evaluate(&e, s, (qaws_scalar)i / 600, (qaws_scalar)j / 600, QAWS_SURFACE_EVAL_POSITION);
			r = sqrt((double)e.position.x * e.position.x + (double)e.position.y * e.position.y + (double)e.position.z * e.position.z);
			dense = fmin(dense, r);
			if (i % 600 == 0 || (i * 7) % 600 == 0)
				if (j % 600 == 0 || (j * 7) % 600 == 0)
					grid = fmin(grid, r);
		}
	TEST_ASSERT(qaws_surface_find_closest_point_certified(&u, &w, &cl, &dist, &gap, s, q) == QAWS_STATUS_OK, "certified on the cap");
	snprintf(msg, sizeof msg, "cap: the 8 x 8 grid sees %.3f; certified %.12f at (%.4f, %.4f), gap %.1e; 601 x 601 samples %.12f",
		grid, (double)dist, (double)u, (double)w, (double)gap, dense);
	printf("  %s\n", msg);
	TEST_ASSERT(grid > 4.9 && dist < 4 && dist <= dense + cp_round(5, 1e-12) && gap <= CP_GAP * 10, msg);
	TEST_ASSERT(fabs(sqrt((double)cl.x * cl.x + (double)cl.y * cl.y + (double)cl.z * cl.z) - dist) <= cp_round(5, 1e-12), "the closest point at that distance");
	/* the same cap borrowed */
	memset(&v, 0, sizeof v);
	v.u_degree = 3; v.v_degree = 3; v.u_count = CP_S; v.v_count = CP_S;
	v.points.data = p; v.points.type = sizeof(qaws_scalar) == 8 ? QAWS_ELEMENT_F64 : QAWS_ELEMENT_F32; v.points.components = 3;
	v.points.count = CP_S * CP_S; v.points.stride = sizeof(qaws_vec3);
	TEST_ASSERT(qaws_surface_create_borrowed(&b, &v, NULL) == QAWS_STATUS_OK, "borrowed cap created");
	qaws_surface_find_closest_point_certified(&u2, &w2, &cl2, &dist2, &gap2, b, q);
	snprintf(msg, sizeof msg, "borrowed cap: the same point (%.1e in u, %.1e in v)", fabs((double)u2 - u), fabs((double)w2 - w));
	TEST_ASSERT(fabs((double)u2 - u) <= 1e-10 && fabs((double)w2 - w) <= 1e-10 && fabs((double)dist2 - dist) <= 1e-12, msg);
	qaws_surface_destroy(b);
	qaws_surface_destroy(s);

	/* a bicubic Bezier patch at points around it: never above dense sampling */
	{
		qaws_vec3 n[16];
		qaws_surface_bezier_desc bd;
		double worst = 0, gmax = 0;
		unsigned int k;
		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
			{
				n[i * 4 + j].x = (qaws_scalar)i;
				n[i * 4 + j].y = (qaws_scalar)j;
				n[i * 4 + j].z = (qaws_scalar)(1.5 * sin(1.7 * i + 0.9 * j) - 0.5 * cos(2.1 * i * j));
			}
		memset(&bd, 0, sizeof bd);
		bd.u_degree = 3; bd.v_degree = 3; bd.control_points = n; bd.u_point_count = 4; bd.v_point_count = 4;
		TEST_ASSERT(qaws_surface_create_bezier(&s, &bd) == QAWS_STATUS_OK, "bicubic patch created");
		for (k = 0; k < 12; k++)
		{
			qaws_vec3 x;
			dense = 1e300;
			x.x = (qaws_scalar)(1.5 + 2 * cos(1.3 * k)); x.y = (qaws_scalar)(1.5 + 2 * sin(2.1 * k)); x.z = (qaws_scalar)(1.5 * cos(0.7 * k + 1));
			qaws_surface_find_closest_point_certified(&u, &w, &cl, &dist, &gap, s, x);
			for (i = 0; i <= 300; i++)
				for (j = 0; j <= 300; j++)
				{
					qaws_surface_eval_result e;
					qaws_surface_evaluate(&e, s, (qaws_scalar)i / 300, (qaws_scalar)j / 300, QAWS_SURFACE_EVAL_POSITION);
					dense = fmin(dense, sqrt(((double)e.position.x - x.x) * ((double)e.position.x - x.x) + ((double)e.position.y - x.y)
						* ((double)e.position.y - x.y) + ((double)e.position.z - x.z) * ((double)e.position.z - x.z)));
				}
			worst = fmax(worst, dist - dense);
			gmax = fmax(gmax, gap);
		}
		snprintf(msg, sizeof msg, "bicubic patch, 12 points: never above dense sampling (worst %.1e), gaps up to %.1e", worst, gmax);
		printf("  %s\n", msg);
		TEST_ASSERT(worst <= cp_round(3, 1e-12) && gmax <= CP_GAP * 10, msg);
		qaws_surface_destroy(s);
	}
	/* a flat square: exact distances, inside and past its edge */
	{
		qaws_vec3 n[4] = { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 }, { 1, 1, 0 } }, x1 = { 0.3f, 0.6f, 2 }, x2 = { 1.5f, 0.5f, 1 };
		qaws_surface_bezier_desc bd;
		memset(&bd, 0, sizeof bd);
		bd.u_degree = 1; bd.v_degree = 1; bd.control_points = n; bd.u_point_count = 2; bd.v_point_count = 2;
		qaws_surface_create_bezier(&s, &bd);
		qaws_surface_find_closest_point_certified(&u, &w, &cl, &dist, &gap, s, x1);
		TEST_ASSERT(fabs(dist - 2) <= 1e-14 && fabs(u - (qaws_scalar)0.3) <= 1e-7 && fabs(w - (qaws_scalar)0.6) <= 1e-7, "square: above it");
		qaws_surface_find_closest_point_certified(&u, &w, &cl, &dist, &gap, s, x2);
		TEST_ASSERT(fabs(dist - sqrt(1.25)) <= cp_round(1, 1e-14) && u == 1 && fabs(w - (qaws_scalar)0.5) <= 1e-7, "square: past its edge");
		qaws_surface_destroy(s);
	}
}

int test_149_closest_point_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 149: certified closest points\n");
	test_spike();
	test_circle();
	test_bezier_3d();
	test_other_families();
	test_surfaces();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
