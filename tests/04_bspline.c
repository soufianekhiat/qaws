#include "test_common.h"

static void test_bspline(void)
{
	printf("test_bspline\n");

	qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
	qaws_bspline_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;
	desc.knots = NULL;
	desc.knot_count = 0;
	desc.is_uniform = 1;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bspline(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_BSPLINE, "kind == BSPLINE");

	/* Evaluate at endpoints of clamped uniform B-spline */
	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	s = qaws_curve_evaluate_2d(curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "start y");

	s = qaws_curve_evaluate_2d(curve, range.max_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)4.0), "end x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "end y");

	qaws_curve_destroy(curve);
}

static void test_bspline_custom_knots(void)
{
	printf("test_bspline_custom_knots\n");

	qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
	qaws_scalar knots[] = { 0, 0, 0, 0, 1, 1, 1, 1 }; /* clamped cubic, 4 CPs */

	qaws_bspline_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;
	desc.knots = knots;
	desc.knot_count = 8;
	desc.is_uniform = 0;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bspline(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	s = qaws_curve_evaluate_2d(curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "custom knots start x");

	s = qaws_curve_evaluate_2d(curve, range.max_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)4.0), "custom knots end x");

	qaws_curve_destroy(curve);
}

/* Unclamped knot vectors with a knot repeated at a domain end leave empty
   spans there: evaluation at the end must use the nearest non-empty span
   (the curve then passes through a control point). */
static void test_bspline_repeated_end_knots(void)
{
	static qaws_scalar const knots_start[10] = { 0, 1, 2, 3, 3, 3, 4, 5, 6, 7 };
	static qaws_scalar const knots_end[10] = { 0, 1, 2, 3, 4, 4, 4, 5, 6, 7 };
	qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0}, {6, 1}, {7, 3} };
	int which;
	printf("test_bspline_repeated_end_knots\n");
	for (which = 0; which < 2; which++)
	{
		qaws_bspline_desc desc;
		qaws_curve* curve = NULL;
		qaws_eval_result_2d at, near;
		qaws_range range;
		qaws_scalar t, inside;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = points;
		desc.control_point_count = 6;
		desc.knots = which ? knots_end : knots_start;
		desc.knot_count = 10;
		TEST_ASSERT_STATUS(qaws_curve_create_bspline(&desc, &curve));
		range = qaws_curve_get_parameter_range(curve);
		t = which ? range.max_value : range.min_value;
		inside = which ? t - (qaws_scalar)1e-4 : t + (qaws_scalar)1e-4;
		TEST_ASSERT_STATUS(qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &at));
		TEST_ASSERT_STATUS(qaws_curve_evaluate_2d(curve, inside, QAWS_EVAL_FLAG_POSITION, &near));
		TEST_ASSERT(fabs(at.position.x - near.position.x) < 1e-2 && fabs(at.position.y - near.position.y) < 1e-2,
			which ? "triple knot at the domain end: finite and continuous" : "triple knot at the domain start: finite and continuous");
		/* the triple knot interpolates a control point */
		TEST_ASSERT(approx_eq(at.position.x, which ? (qaws_scalar)4.0 : (qaws_scalar)3.0) &&
			approx_eq(at.position.y, which ? (qaws_scalar)0.0 : (qaws_scalar)2.0), "triple knot interpolates its control point (P2 at the start, P3 at the end)");
		qaws_curve_destroy(curve);
	}
}

int test_04_bspline_main(void)
{
	g_pass = 0; g_fail = 0;
	test_bspline();
	test_bspline_custom_knots();
	test_bspline_repeated_end_knots();
	return g_fail > 0 ? 1 : 0;
}
