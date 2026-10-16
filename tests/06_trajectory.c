#include "test_common.h"

static void test_trajectory(void)
{
	printf("test_trajectory\n");

	qaws_vec2 positions[] = { {0, 0}, {5, 0}, {5, 5} };
	qaws_scalar times[] = { 0, 1, 2 };

	qaws_trajectory_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.key_positions = positions;
	desc.key_count = 3;
	desc.key_times = times;
	desc.key_time_count = 3;
	desc.key_velocities = NULL;
	desc.key_velocity_count = 0;
	desc.key_accelerations = NULL;
	desc.key_acceleration_count = 0;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_trajectory(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_TRAJECTORY, "kind == TRAJECTORY");

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "start y");

	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "key1 x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "key1 y");

	qaws_curve_destroy(curve);
}

int test_06_trajectory_main(void)
{
	g_pass = 0; g_fail = 0;
	test_trajectory();
	return g_fail > 0 ? 1 : 0;
}
