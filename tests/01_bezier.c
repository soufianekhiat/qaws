#include "test_common.h"

static void test_bezier_linear(void)
{
	printf("test_bezier_linear\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(curve != NULL, "curve allocated");

	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_BEZIER, "kind == BEZIER");
	TEST_ASSERT(qaws_curve_get_degree(curve) == 1, "degree == 1");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 1, "span_count == 1");

	qaws_eval_result_2d result;
	s = qaws_curve_evaluate_2d(&result, curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(result.position.x, (qaws_scalar)0.5), "midpoint x");
	TEST_ASSERT(approx_eq(result.position.y, (qaws_scalar)0.5), "midpoint y");

	qaws_curve_destroy(curve);
}

static void test_bezier_cubic(void)
{
	printf("test_bezier_cubic\n");

	qaws_vec2 points[] = { {0, 0}, {0, 1}, {1, 1}, {1, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	/* Evaluate at endpoints */
	qaws_eval_result_2d r0, r1;
	s = qaws_curve_evaluate_2d(&r0, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r0.position.x, (qaws_scalar)0.0), "start x");
	TEST_ASSERT(approx_eq(r0.position.y, (qaws_scalar)0.0), "start y");

	s = qaws_curve_evaluate_2d(&r1, curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r1.position.x, (qaws_scalar)1.0), "end x");
	TEST_ASSERT(approx_eq(r1.position.y, (qaws_scalar)0.0), "end y");

	/* Evaluate at midpoint */
	qaws_eval_result_2d rm;
	s = qaws_curve_evaluate_2d(&rm, curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(rm.position.x, (qaws_scalar)0.5), "mid x");
	TEST_ASSERT(approx_eq(rm.position.y, (qaws_scalar)0.75), "mid y");

	qaws_curve_destroy(curve);
}

static void test_bezier_3d(void)
{
	printf("test_bezier_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_dimension(curve) == QAWS_DIMENSION_3D, "dim == 3D");

	qaws_eval_result_3d r;
	s = qaws_curve_evaluate_3d(&r, curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");

	s = qaws_curve_evaluate_3d(&r, curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.z, (qaws_scalar)1.0), "end z");

	qaws_curve_destroy(curve);
}

static void test_bezier_derivatives(void)
{
	printf("test_bezier_derivatives\n");

	/* Linear Bezier from (0,0) to (3,4): D1 should be (3,4) everywhere */
	qaws_vec2 pts_lin[] = { {0, 0}, {3, 4} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = pts_lin;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_D1);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)3.0), "linear D1.x == 3");
	TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)4.0), "linear D1.y == 4");

	qaws_curve_destroy(curve);

	/* Cubic Bezier: (0,0),(0,1),(1,1),(1,0) at t=0.5
	   D1 = 3*((0,1)-(0,0)) at t=0 => D1(0) = (0,3)
	   By symmetry, D1(0.5) should have x>0, y=0.
	   D1(t) = 3[(1-t)^2(P1-P0) + 2(1-t)t(P2-P1) + t^2(P3-P2)]
	   D1(0.5) = 3[0.25*(0,1) + 0.5*(1,0) + 0.25*(0,-1)] = 3*(0.5,0) = (1.5,0) */
	qaws_vec2 pts_cub[] = { {0, 0}, {0, 1}, {1, 1}, {1, 0} };
	desc.degree = 3;
	desc.control_points = pts_cub;
	desc.control_point_count = 4;

	s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	s = qaws_curve_evaluate_2d(&r, curve,
		(qaws_scalar)0.5, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)1.5), "cubic D1.x at 0.5");
	TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)0.0), "cubic D1.y at 0.5");
	/* D2 by symmetry at midpoint: D2.x = 0, D2.y = -6 */
	TEST_ASSERT(approx_eq(r.d2.x, (qaws_scalar)0.0), "cubic D2.x at 0.5");
	TEST_ASSERT(approx_eq(r.d2.y, (qaws_scalar)-6.0), "cubic D2.y at 0.5");
	TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D3, "D3 flag set");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Span-level evaluation                                               */
/* ------------------------------------------------------------------ */

static void test_evaluate_span(void)
{
	printf("test_evaluate_span\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_span_2d(&r, curve, 0, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.5), "span eval midpoint x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.5), "span eval midpoint y");

	/* Out of range span index */
	s = qaws_curve_evaluate_span_2d(&r, curve, 5, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT(s == QAWS_STATUS_OUT_OF_RANGE, "out-of-range span returns error");

	qaws_curve_destroy(curve);
}

int test_01_bezier_main(void)
{
	g_pass = 0; g_fail = 0;
	test_bezier_linear();
	test_bezier_cubic();
	test_bezier_3d();
	test_bezier_derivatives();
	test_evaluate_span();
	return g_fail > 0 ? 1 : 0;
}
