#include "test_common.h"

static void test_error_handling(void)
{
	printf("test_error_handling\n");

	qaws_status s;

	/* NULL arguments */
	s = qaws_curve_create_bezier(NULL, NULL);
	TEST_ASSERT(s != QAWS_STATUS_OK, "NULL desc rejected");

	qaws_curve* curve = NULL;
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = NULL;
	desc.control_point_count = 2;
	s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT(s != QAWS_STATUS_OK, "NULL control_points rejected");

	/* Invalid degree */
	qaws_vec2 points[] = { {0, 0} };
	desc.control_points = points;
	desc.control_point_count = 1;
	desc.degree = 0;
	s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT(s != QAWS_STATUS_OK, "degree 0 rejected");

	/* Wrong CP count */
	desc.degree = 1;
	desc.control_point_count = 5;
	s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT(s != QAWS_STATUS_OK, "wrong CP count rejected");

	/* Invalid dimension */
	desc.dimension = (qaws_dimension)99;
	desc.degree = 1;
	desc.control_point_count = 2;
	qaws_vec2 pts2[] = { {0,0}, {1,1} };
	desc.control_points = pts2;
	s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT(s != QAWS_STATUS_OK, "invalid dimension rejected");

	/* Evaluate NULL curve */
	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(&r, NULL, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT(s != QAWS_STATUS_OK, "NULL curve evaluate rejected");

	/* Wrong dimension eval */
	desc.dimension = QAWS_DIMENSION_2D;
	s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_3d r3;
	s = qaws_curve_evaluate_3d(&r3, curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION, "2D curve + evaluate_3d => INVALID_DIMENSION");

	qaws_curve_destroy(curve);
}

int test_24_error_handling_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_error_handling();

	printf("24_error_handling: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
