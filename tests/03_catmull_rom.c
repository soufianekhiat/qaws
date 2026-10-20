#include "test_common.h"

static void test_catmull_rom(void)
{
	printf("test_catmull_rom\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };

	qaws_catmull_rom_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 4;
	desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_catmull_rom(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_CATMULL_ROM, "kind == CATMULL_ROM");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 1, "span_count == 1");

	/* The curve interpolates the interior points (P1 and P2) */
	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "interpolates P1 x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)1.0), "interpolates P1 y");

	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.0), "interpolates P2 x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "interpolates P2 y");

	qaws_curve_destroy(curve);
}

static void test_catmull_rom_variants(void)
{
	printf("test_catmull_rom_variants\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };

	/* Uniform */
	{
		qaws_catmull_rom_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 4;
		desc.parameterization = QAWS_PARAMETERIZATION_UNIFORM;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_catmull_rom(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "uniform interpolates P1 x");

		qaws_curve_destroy(curve);
	}

	/* Chordal */
	{
		qaws_catmull_rom_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 4;
		desc.parameterization = QAWS_PARAMETERIZATION_CHORDAL;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_catmull_rom(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "chordal interpolates P1 x");

		qaws_curve_destroy(curve);
	}
}

int test_03_catmull_rom_main(void)
{
	g_pass = 0; g_fail = 0;
	test_catmull_rom();
	test_catmull_rom_variants();
	return g_fail > 0 ? 1 : 0;
}
