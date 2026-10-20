#include "test_common.h"

static void test_nurbs(void)
{
	printf("test_nurbs\n");

	/* Quarter-circle NURBS: degree 2, 3 control points */
	qaws_vec2 points[] = { {1, 0}, {1, 1}, {0, 1} };
	qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };
	qaws_scalar weights[] = { 1, (qaws_scalar)0.70710678118, 1 };

	qaws_nurbs_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 2;
	desc.control_points = points;
	desc.control_point_count = 3;
	desc.knots = knots;
	desc.knot_count = 6;
	desc.weights = weights;
	desc.weight_count = 3;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_nurbs(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_NURBS, "kind == NURBS");
	TEST_ASSERT(qaws_curve_is_rational(curve) == 1, "is_rational");

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "start x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "start y");

	/* Midpoint of quarter circle should be near (cos(45), sin(45)) */
	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	qaws_scalar dist = (qaws_scalar)sqrt((double)(r.position.x * r.position.x + r.position.y * r.position.y));
	TEST_ASSERT(approx_eq(dist, (qaws_scalar)1.0), "midpoint on unit circle");

	qaws_curve_destroy(curve);
}

int test_05_nurbs_main(void)
{
	g_pass = 0; g_fail = 0;
	test_nurbs();
	return g_fail > 0 ? 1 : 0;
}
