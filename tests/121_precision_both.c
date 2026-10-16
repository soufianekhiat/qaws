/*
 * Test 121: both precisions in one binary
 *
 * The library holds every scalar-bearing function and type twice, as _f32
 * and _f64 (docs/parameter_convention_plan.md). qaws_both.h declares both
 * in one translation unit: the same cubic Bezier built and evaluated
 * through each, the two agreeing to single precision, the double one
 * exact to double; the types' sizes those of each scalar; the unsuffixed
 * names still the caller's precision.
 */

#include "qaws_both.h"
#include "test_common.h"
#include <math.h>

static void test_both(void)
{
	static double const P[4][2] = { { 0, 0 }, { 1, 2 }, { 3, 2 }, { 4, 0 } };
	double const t = 0.3, b0 = (1 - t) * (1 - t) * (1 - t), b1 = 3 * (1 - t) * (1 - t) * t, b2 = 3 * (1 - t) * t * t, b3 = t * t * t;
	double const ex = b0 * P[0][0] + b1 * P[1][0] + b2 * P[2][0] + b3 * P[3][0], ey = b0 * P[0][1] + b1 * P[1][1] + b2 * P[2][1] + b3 * P[3][1];
	qaws_vec2_f32 p32[4];
	qaws_vec2_f64 p64[4];
	qaws_bezier_desc_f32 d32;
	qaws_bezier_desc_f64 d64;
	qaws_curve_f32* c32 = NULL;
	qaws_curve_f64* c64 = NULL;
	qaws_eval_result_2d_f32 r32 = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, 0 };
	qaws_eval_result_2d_f64 r64 = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, 0 };
	qaws_status s32, s64;
	double e32, e64, agree;
	unsigned int i;
	int sizes_ok;   /* (a variable: MSVC flags a constant condition in the assertion) */
	char msg[300];
	for (i = 0; i < 4; i++)
	{
		p32[i].x = (qaws_f32)P[i][0]; p32[i].y = (qaws_f32)P[i][1];
		p64[i].x = P[i][0]; p64[i].y = P[i][1];
	}
	d32.dimension = QAWS_DIMENSION_2D; d32.degree = 3; d32.control_points = p32; d32.control_point_count = 4;
	d64.dimension = QAWS_DIMENSION_2D; d64.degree = 3; d64.control_points = p64; d64.control_point_count = 4;
	s32 = qaws_curve_create_bezier_f32(&d32, &c32);
	s64 = qaws_curve_create_bezier_f64(&d64, &c64);
	if (s32 == QAWS_STATUS_OK) s32 = qaws_curve_evaluate_2d_f32(c32, (qaws_f32)t, QAWS_EVAL_FLAG_POSITION, &r32);
	if (s64 == QAWS_STATUS_OK) s64 = qaws_curve_evaluate_2d_f64(c64, t, QAWS_EVAL_FLAG_POSITION, &r64);
	e32 = fmax(fabs((double)r32.position.x - ex), fabs((double)r32.position.y - ey));
	e64 = fmax(fabs(r64.position.x - ex), fabs(r64.position.y - ey));
	agree = fmax(fabs((double)r32.position.x - r64.position.x), fabs((double)r32.position.y - r64.position.y));
	sprintf(msg, "a cubic Bezier through both precisions (status %d, %d): f32 off %.2g, f64 off %.2g, apart %.2g", (int)s32, (int)s64, e32, e64, agree);
	TEST_ASSERT(s32 == QAWS_STATUS_OK && s64 == QAWS_STATUS_OK && e32 <= 1e-6 && e64 <= 1e-15 && agree <= 1e-6, msg);
	qaws_curve_destroy_f32(c32);
	qaws_curve_destroy_f64(c64);
	sprintf(msg, "sizes: qaws_vec2_f32 %u, qaws_vec2_f64 %u, qaws_scalar %u (the caller's, QAWS_SCALAR_IS_FLOAT %d)", (unsigned int)sizeof(qaws_vec2_f32),
		(unsigned int)sizeof(qaws_vec2_f64), (unsigned int)sizeof(qaws_scalar), (int)QAWS_SCALAR_IS_FLOAT);
	sizes_ok = sizeof(qaws_vec2_f32) == 2 * sizeof(float) && sizeof(qaws_vec2_f64) == 2 * sizeof(double) &&
		sizeof(qaws_scalar) == (QAWS_SCALAR_IS_FLOAT ? sizeof(float) : sizeof(double)) && sizeof(qaws_vec2) == 2 * sizeof(qaws_scalar);
	TEST_ASSERT(sizes_ok, msg);
}

int test_121_precision_both_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 121: both precisions in one binary\n");
	test_both();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
