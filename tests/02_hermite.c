#include "test_common.h"

static void test_hermite(void)
{
	printf("test_hermite\n");

	qaws_vec2 points[] = { {0, 0}, {1, 0} };
	qaws_vec2 tangents[] = { {1, 1}, {1, -1} };

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = points;
	desc.derivatives = tangents;
	desc.point_count = 2;
	desc.derivative_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_HERMITE, "kind == HERMITE");

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");

	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "end x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "end y");

	qaws_curve_destroy(curve);
}

static void test_hermite_3d(void)
{
	printf("test_hermite_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {1, 0, 0} };
	qaws_vec3 tangents[] = { {1, 0, 1}, {1, 0, -1} };

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 3;
	desc.points = points;
	desc.derivatives = tangents;
	desc.point_count = 2;
	desc.derivative_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_3d r;
	s = qaws_curve_evaluate_3d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "hermite 3d start x");
	TEST_ASSERT(approx_eq(r.position.z, (qaws_scalar)0.0), "hermite 3d start z");

	s = qaws_curve_evaluate_3d(&r, curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "hermite 3d end x");

	qaws_curve_destroy(curve);
}

int test_02_hermite_main(void)
{
	g_pass = 0; g_fail = 0;
	test_hermite();
	test_hermite_3d();
	return g_fail > 0 ? 1 : 0;
}
