#include "test_common.h"

static void test_polynomial(void)
{
	printf("test_polynomial\n");

	/* Linear: P(t) = (t, 2t), t in [0, 1] */
	{
		qaws_scalar coeffs[] = { 0, 0,   1, 2 }; /* c0=(0,0), c1=(1,2) */
		qaws_polynomial_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 1;
		desc.coefficients = coeffs;
		desc.coefficient_count = 2;
		desc.t_min = 0;
		desc.t_max = 1;

		TEST_ASSERT_STATUS(qaws_curve_create_polynomial(&crv, &desc));
		TEST_ASSERT(qaws_curve_get_kind(crv) == QAWS_CURVE_KIND_POLYNOMIAL, "poly kind");
		TEST_ASSERT(qaws_curve_is_rational(crv) == 0, "poly not rational");

		qaws_curve_evaluate_2d(&r, crv, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.5) && approx_eq(r.position.y, 1),
			"poly linear pos at 0.5");
		TEST_ASSERT(approx_eq(r.d1.x, 1) && approx_eq(r.d1.y, 2),
			"poly linear d1");

		/* Endpoints */
		qaws_curve_evaluate_2d(&r, crv, 0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r.position.x, 0) && approx_eq(r.position.y, 0),
			"poly linear start");
		qaws_curve_evaluate_2d(&r, crv, 1, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r.position.x, 1) && approx_eq(r.position.y, 2),
			"poly linear end");

		qaws_curve_destroy(crv);
	}

	/* Cubic: P(t) = (t^3, t^2), t in [0, 1] */
	{
		qaws_scalar coeffs[] = { 0,0,  0,0,  0,1,  1,0 }; /* c0,c1,c2,c3 */
		qaws_polynomial_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.coefficients = coeffs;
		desc.coefficient_count = 4;
		desc.t_min = 0;
		desc.t_max = 1;

		TEST_ASSERT_STATUS(qaws_curve_create_polynomial(&crv, &desc));

		qaws_curve_evaluate_2d(&r, crv, 1, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r.position.x, 1) && approx_eq(r.position.y, 1),
			"poly cubic at t=1");

		qaws_curve_evaluate_2d(&r, crv,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.125), "poly cubic x at 0.5");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.25), "poly cubic y at 0.5");
		/* D1 of (t^3, t^2) = (3t^2, 2t), at t=0.5 => (0.75, 1.0) */
		TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)0.75), "poly cubic d1.x at 0.5");
		TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)1.0), "poly cubic d1.y at 0.5");
		/* D2 = (6t, 2), at t=0.5 => (3, 2) */
		TEST_ASSERT(approx_eq(r.d2.x, (qaws_scalar)3.0), "poly cubic d2.x at 0.5");
		TEST_ASSERT(approx_eq(r.d2.y, (qaws_scalar)2.0), "poly cubic d2.y at 0.5");
		/* D3 = (6, 0) */
		TEST_ASSERT(approx_eq(r.d3.x, (qaws_scalar)6.0), "poly cubic d3.x");
		TEST_ASSERT(approx_eq(r.d3.y, (qaws_scalar)0.0), "poly cubic d3.y");

		qaws_curve_destroy(crv);
	}

	/* Non-unit parameter range: t in [-1, 1] */
	{
		/* P(t) = (t, t^2): parabola */
		qaws_scalar coeffs[] = { 0, 0,   1, 0,   0, 1 };
		qaws_polynomial_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.coefficients = coeffs;
		desc.coefficient_count = 3;
		desc.t_min = -1;
		desc.t_max = 1;

		TEST_ASSERT_STATUS(qaws_curve_create_polynomial(&crv, &desc));

		/* At global param -1 (start): should be (-1, 1) */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv);
			qaws_curve_evaluate_2d(&r, crv, rng.min_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r.position.x, -1), "poly parabola start x");
			TEST_ASSERT(approx_eq(r.position.y, 1), "poly parabola start y");

			/* At global param 1 (end): (1, 1) */
			qaws_curve_evaluate_2d(&r, crv, rng.max_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r.position.x, 1), "poly parabola end x");
			TEST_ASSERT(approx_eq(r.position.y, 1), "poly parabola end y");
		}

		qaws_curve_destroy(crv);
	}

	/* Error: NULL desc */
	{
		qaws_curve *crv = NULL;
		TEST_ASSERT(qaws_curve_create_polynomial(&crv, NULL) != QAWS_STATUS_OK,
			"poly null desc rejected");
	}
}

int test_10_polynomial_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_polynomial();

	printf("10_polynomial: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
