/*
 * Test 148: arc length and its inverse (qaws_curve_compute_arc_length,
 * qaws_curve_compute_parameter_at_arc_length)
 *
 *   - a cubic B-spline against Mathematica (tests/reference/58_arc_length.h,
 *     30 digits): its length and the parameters at three lengths;
 *   - a rational quarter circle: pi r;
 *   - a cusp (zero speed at the start) and a stationary span (zero speed over
 *     a whole span): the inverse lands where the length is reached;
 *   - from a parameter inside the range; the ends; refusals.
 */

#include "test_common.h"
#include "reference/58_arc_length.h"
#include <math.h>
#include <string.h>

static qaws_scalar al_length(qaws_curve const* c, double a, double b)
{
	qaws_scalar l = 0;
	qaws_curve_compute_arc_length(c, (qaws_scalar)a, (qaws_scalar)b, &l);
	return l;
}

static void test_reference(void)
{
	/* the curve of tests/reference/58_arc_length.wls */
	qaws_vec3 p[6] = { { 0, 0, 0 }, { 1, 2, 0.5 }, { 2.5, 2, -0.5 }, { 3.5, 0.5, 1 }, { 5, 1, 0 }, { 6, 3, 0.5 } };
	qaws_scalar knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };
	double sigma[3], total = ref_total_arc[0];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	char msg[200];
	unsigned int i;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = p; d.control_point_count = 6; d.knots = knots; d.knot_count = 10;
	TEST_ASSERT(qaws_curve_create_bspline(&d, &c) == QAWS_STATUS_OK, "reference B-spline created");
	snprintf(msg, sizeof msg, "length = Mathematica's (relative error %.1e)", fabs(al_length(c, 0, 3) - total) / total);
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(al_length(c, 0, 3) - total) <= 1e-13 * total, msg);
	/* the targets {distance, fraction}: 0.35 of the length, 1.7, the length less 0.5 */
	sigma[0] = 0.35 * total;
	sigma[1] = 1.7;
	sigma[2] = total - 0.5;
	for (i = 0; i < 3; i++)
	{
		qaws_scalar t = -1;
		qaws_status st = qaws_curve_compute_parameter_at_arc_length(c, 0, (qaws_scalar)sigma[i], &t);
		snprintf(msg, sizeof msg, "parameter at length %.4f = Mathematica's (error %.1e)", sigma[i], fabs(t - ref_t_arc[i]));
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && fabs(t - ref_t_arc[i]) <= 1e-11, msg);
	}
	/* from inside the range: the length from there is reached */
	{
		qaws_scalar t = -1;
		qaws_status st = qaws_curve_compute_parameter_at_arc_length(c, (qaws_scalar)1.3, (qaws_scalar)2.0, &t);
		TEST_ASSERT(st == QAWS_STATUS_OK && fabs(al_length(c, 1.3, t) - 2.0) <= 1e-11, "from t = 1.3, 2.0 further");
	}
	qaws_curve_destroy(c);
}

static void test_circle(void)
{
	qaws_vec2 p[3] = { { 2, 0 }, { 2, 2 }, { 0, 2 } };
	qaws_scalar knots[6] = { 0, 0, 0, 1, 1, 1 }, w[3] = { 1, (qaws_scalar)(0.5 * sqrt(2.0)), 1 }, t = -1;
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = p; d.control_point_count = 3;
	d.knots = knots; d.knot_count = 6; d.weights = w; d.weight_count = 3;
	TEST_ASSERT(qaws_curve_create_nurbs(&d, &c) == QAWS_STATUS_OK, "quarter circle created");
	snprintf(msg, sizeof msg, "quarter circle of radius 2: pi (error %.1e)", fabs(al_length(c, 0, 1) - M_PI));
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(al_length(c, 0, 1) - M_PI) <= 1e-13, msg);
	/* a third of the arc: the point at 30 degrees */
	{
		qaws_eval_result_2d e;
		qaws_curve_compute_parameter_at_arc_length(c, 0, (qaws_scalar)(M_PI / 3), &t);
		qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION, &e);
		snprintf(msg, sizeof msg, "a third of the arc: the point at 30 degrees (error %.1e)", fabs(e.position.x - sqrt(3.0)) + fabs(e.position.y - 1));
		TEST_ASSERT(fabs(e.position.x - sqrt(3.0)) + fabs(e.position.y - 1) <= 1e-11, msg);
	}
	qaws_curve_destroy(c);
}

static void test_zero_speed(void)
{
	/* a cusp at the start: P0 = P1 */
	qaws_vec2 b[4] = { { 0, 0 }, { 0, 0 }, { 1, 1 }, { 2, 0 } };
	/* a polyline with a repeated point: no motion over [1, 2] */
	qaws_vec2 q[4] = { { 0, 0 }, { 1, 0 }, { 1, 0 }, { 2, 0 } };
	qaws_scalar knots[6] = { 0, 0, 1, 2, 3, 3 }, t = -1;
	qaws_bezier_desc bd;
	qaws_bspline_desc qd;
	qaws_curve *cusp = NULL, *poly = NULL;
	double total, worst = 0;
	char msg[200];
	unsigned int i;
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 3; bd.control_points = b; bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &cusp);
	total = al_length(cusp, 0, 1);
	for (i = 0; i <= 20; i++)
	{
		double s = total * (i == 0 ? 1e-9 : i / 20.0);
		qaws_status st = qaws_curve_compute_parameter_at_arc_length(cusp, 0, (qaws_scalar)s, &t);
		double e = st == QAWS_STATUS_OK ? fabs(al_length(cusp, 0, t) - s) : 1;
		if (e > worst)
			worst = e;
	}
	snprintf(msg, sizeof msg, "cusp at the start: the lengths reached (worst %.1e of %.3f)", worst, total);
	printf("  %s\n", msg);
	TEST_ASSERT(worst <= 1e-11 * total, msg);

	memset(&qd, 0, sizeof qd);
	qd.dimension = QAWS_DIMENSION_2D; qd.degree = 1; qd.control_points = q; qd.control_point_count = 4; qd.knots = knots; qd.knot_count = 6;
	qaws_curve_create_bspline(&qd, &poly);
	TEST_ASSERT(fabs(al_length(poly, 0, 3) - 2) <= 1e-14, "polyline with a stationary span: length 2");
	qaws_curve_compute_parameter_at_arc_length(poly, 0, 1, &t);
	TEST_ASSERT(t >= 1 - 1e-9 && t <= 2 + 1e-9, "length 1: in the stationary span");
	qaws_curve_compute_parameter_at_arc_length(poly, 0, (qaws_scalar)1.5, &t);
	snprintf(msg, sizeof msg, "length 1.5: t = 2.5 past the stationary span (%.12f)", (double)t);
	TEST_ASSERT(fabs(t - 2.5) <= 1e-11, msg);
	qaws_curve_destroy(cusp);
	qaws_curve_destroy(poly);
}

static void test_ends_and_refusals(void)
{
	qaws_vec2 p[3] = { { 0, 0 }, { 1, 2 }, { 3, 1 } };
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	qaws_scalar t = -1, total;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = p; d.control_point_count = 3;
	qaws_curve_create_bezier(&d, &c);
	total = al_length(c, 0, 1);
	TEST_ASSERT(qaws_curve_compute_parameter_at_arc_length(c, (qaws_scalar)0.25, 0, &t) == QAWS_STATUS_OK && t == (qaws_scalar)0.25, "length 0: the start");
	TEST_ASSERT(qaws_curve_compute_parameter_at_arc_length(c, 0, total, &t) == QAWS_STATUS_OK && t == 1, "the whole length: the end");
	TEST_ASSERT(qaws_curve_compute_parameter_at_arc_length(c, 0, total * (qaws_scalar)1.001, &t) == QAWS_STATUS_OUT_OF_RANGE, "past the end refused");
	TEST_ASSERT(qaws_curve_compute_parameter_at_arc_length(c, 0, -1, &t) == QAWS_STATUS_OUT_OF_RANGE, "a negative length refused");
	TEST_ASSERT(qaws_curve_compute_parameter_at_arc_length(c, 2, 0, &t) == QAWS_STATUS_OUT_OF_RANGE, "a start outside the range refused");
	TEST_ASSERT(qaws_curve_compute_parameter_at_arc_length(NULL, 0, 0, &t) == QAWS_STATUS_INVALID_ARGUMENT, "no curve refused");
	TEST_ASSERT(qaws_curve_compute_parameter_at_arc_length(c, 0, 0, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "no output refused");
	qaws_curve_destroy(c);
}

int test_148_arc_length_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 148: arc length and its inverse\n");
	test_reference();
	test_circle();
	test_zero_speed();
	test_ends_and_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
