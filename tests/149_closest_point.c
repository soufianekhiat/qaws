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
		qaws_curve_evaluate_2d(c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / n, QAWS_EVAL_FLAG_POSITION, &e);
		best = fmin(best, hypot(e.position.x - q.x, e.position.y - q.y));
	}
	return best;
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
	qaws_vec2 q = { 0, 0 };
	qaws_bspline_desc d;
	qaws_curve_view v;
	qaws_curve *c = NULL, *b = NULL;
	qaws_scalar t = -1, dist = -1, gap = -1, t2 = -1, dist2 = -1, gap2 = -1, ts = -1;
	double dense;
	char msg[240];
	cp_spiked_points(p);
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p; d.control_point_count = CP_N;
	TEST_ASSERT(qaws_curve_create_bspline(&d, &c) == QAWS_STATUS_OK, "spiked circle created");
	/* sampling and Newton alone, from the 32 samples: they stop on the circle */
	{
		qaws_range r = qaws_curve_get_parameter_range(c);
		double best = 1e300;
		unsigned int i;
		for (i = 0; i < 32; i++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 31, QAWS_EVAL_FLAG_POSITION, &e);
			best = fmin(best, hypot(e.position.x, e.position.y));
		}
		snprintf(msg, sizeof msg, "the 32 samples see no spike (nearest %.3f)", best);
		printf("  %s\n", msg);
		TEST_ASSERT(best > 4.9, msg);
	}
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(c, q, &t, &dist, &gap) == QAWS_STATUS_OK, "certified search runs");
	dense = cp_dense_2d(c, q, 400000);
	snprintf(msg, sizeof msg, "the spike found: distance %.12f at t = %.6f, gap %.1e; dense sampling %.12f", (double)dist, (double)t, (double)gap, dense);
	printf("  %s\n", msg);
	TEST_ASSERT(dist < 4 && dist <= dense + 1e-12 && t > 64 && t < 68, msg);
	TEST_ASSERT(gap >= 0 && gap <= 1e-8 * 10, "gap within the tolerance");
	TEST_ASSERT(qaws_curve_find_closest_parameter_2d(c, q, &ts) == QAWS_STATUS_OK && ts == t, "find_closest_parameter gives the same");

	/* the same curve borrowed */
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_2D;
	v.degree = 3;
	v.points.data = p; v.points.type = QAWS_ELEMENT_F64; v.points.components = 2; v.points.count = CP_N; v.points.stride = 2 * sizeof(double);
	TEST_ASSERT(qaws_curve_create_borrowed(&v, NULL, &b) == QAWS_STATUS_OK, "borrowed view created");
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(b, q, &t2, &dist2, &gap2) == QAWS_STATUS_OK, "certified on the view");
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
	TEST_ASSERT(qaws_curve_create_nurbs(&d, &c) == QAWS_STATUS_OK, "rational circle created");
	qaws_curve_find_closest_point_certified_2d(c, q, &t, &dist, &gap);
	snprintf(msg, sizeof msg, "rational circle: distance sqrt(10) - 2 (error %.1e, gap %.1e)", fabs(dist - (sqrt(10.0) - 2)), (double)gap);
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(dist - (sqrt(10.0) - 2)) <= 1e-13 && gap <= 1e-8, msg);
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(c, center, &t, &dist, &gap) == QAWS_STATUS_OK, "the center: it terminates");
	snprintf(msg, sizeof msg, "the center: distance 2 (error %.1e), gap %.2e within the work bound", fabs(dist - 2), (double)gap);
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(dist - 2) <= 1e-13 && gap >= 0 && gap <= 2, msg);
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
	qaws_curve_create_bezier(&d, &c);
	for (k = 0; k < 24; k++)
	{
		qaws_vec3 q;
		qaws_scalar t, dist, gap;
		double dense = 1e300;
		q.x = (qaws_scalar)(3.5 + 4 * cos(1.7 * k)); q.y = (qaws_scalar)(2.5 * sin(2.3 * k)); q.z = (qaws_scalar)(2 * cos(0.9 * k + 1));
		qaws_curve_find_closest_point_certified_3d(c, q, &t, &dist, &gap);
		for (i = 0; i <= 100000; i++)
		{
			qaws_eval_result_3d e;
			qaws_curve_evaluate_3d(c, (qaws_scalar)i / 100000, QAWS_EVAL_FLAG_POSITION, &e);
			dense = fmin(dense, sqrt((e.position.x - q.x) * (e.position.x - q.x) + (e.position.y - q.y) * (e.position.y - q.y)
				+ (e.position.z - q.z) * (e.position.z - q.z)));
		}
		/* the dense minimum is no lower than the certified distance, and the certified one is within its gap of the truth */
		worst = fmax(worst, dist - dense);
		gap_max = fmax(gap_max, gap);
	}
	snprintf(msg, sizeof msg, "3D quintic, 24 points: never above dense sampling (worst %.1e), gaps up to %.1e", worst, gap_max);
	printf("  %s\n", msg);
	TEST_ASSERT(worst <= 1e-12 && gap_max <= 1e-7, msg);
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
	TEST_ASSERT(qaws_curve_create_hermite(&d, &c) == QAWS_STATUS_OK, "Hermite curve created");
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(c, q, &t, &dist, &gap) == QAWS_STATUS_UNSUPPORTED_OPERATION, "Hermite: not certified");
	TEST_ASSERT(qaws_curve_find_closest_parameter_2d(c, q, &t) == QAWS_STATUS_OK && fabs(t - 0.5) <= 1e-6, "Hermite: still answered (t = 0.5)");
	TEST_ASSERT(qaws_curve_find_closest_point_certified_2d(c, q, &t, NULL, &gap) == QAWS_STATUS_INVALID_ARGUMENT, "no output refused");
	qaws_curve_destroy(c);
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
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
