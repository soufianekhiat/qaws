#include "test_common.h"

static void test_clothoid(void)
{
	printf("test_clothoid\n");

	/* Straight line clothoid: kappa0 = kappa1 = 0 */
	{
		qaws_clothoid_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;

		memset(&desc, 0, sizeof(desc));
		desc.origin_x = 0;
		desc.origin_y = 0;
		desc.start_angle = 0; /* heading along +x */
		desc.start_curvature = 0;
		desc.end_curvature = 0;
		desc.length = 2;

		TEST_ASSERT_STATUS(qaws_curve_create_clothoid(&desc, &crv));
		TEST_ASSERT(crv != NULL, "clothoid created");
		TEST_ASSERT(qaws_curve_get_kind(crv) == QAWS_CURVE_KIND_CLOTHOID, "clothoid kind");
		TEST_ASSERT(qaws_curve_is_rational(crv) == 0, "clothoid not rational");
		TEST_ASSERT(qaws_curve_get_dimension(crv) == QAWS_DIMENSION_2D, "clothoid 2d");

		/* At start: origin */
		qaws_curve_evaluate_2d(crv, 0, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, 0) && approx_eq(r.position.y, 0),
			"clothoid start");

		/* At end (s=2): should be at (2, 0) for zero curvature heading +x */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv);
			TEST_ASSERT(approx_eq(rng.max_value, 2), "clothoid range = length");
			qaws_curve_evaluate_2d(crv, rng.max_value, QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT(approx_eq(r.position.x, 2) && approx_eq(r.position.y, 0),
				"clothoid straight line end");
		}

		/* D1 for straight line: the unit tangent (cos(0), sin(0)) = (1, 0)
		   (the parameter is arc length) */
		qaws_curve_evaluate_2d(crv, 0, QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(approx_eq(r.d1.x, 1), "clothoid straight d1.x");
		TEST_ASSERT(approx_eq(r.d1.y, 0), "clothoid straight d1.y");

		qaws_curve_destroy(crv);
	}

	/* Circular arc clothoid: kappa0 = kappa1 = 1 (radius 1 circle) */
	{
		qaws_clothoid_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;
		qaws_scalar pi_half = (qaws_scalar)1.5707963268;

		memset(&desc, 0, sizeof(desc));
		desc.origin_x = 0;
		desc.origin_y = 0;
		desc.start_angle = 0;
		desc.start_curvature = 1;
		desc.end_curvature = 1;
		desc.length = pi_half; /* quarter circle */

		TEST_ASSERT_STATUS(qaws_curve_create_clothoid(&desc, &crv));

		/* End point: (sin(pi/2), 1-cos(pi/2)) = (1, 1) */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv);
			qaws_curve_evaluate_2d(crv, rng.max_value, QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT(approx_eq(r.position.x, 1), "clothoid circle x");
			TEST_ASSERT(approx_eq(r.position.y, 1), "clothoid circle y");
		}

		/* D1 at start = (cos(0), sin(0)) = (1, 0): the unit tangent */
		qaws_curve_evaluate_2d(crv, 0,
			QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);
		TEST_ASSERT(approx_eq(r.d1.x, 1), "clothoid d1.x at start");
		TEST_ASSERT(approx_eq(r.d1.y, 0), "clothoid d1.y at start");
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D2, "clothoid d2 valid");
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D3, "clothoid d3 valid");
		/* D2 at start: L^2 * kappa * (-sin(0), cos(0)) = L^2 * 1 * (0, 1) */
		TEST_ASSERT(approx_eq(r.d2.x, 0), "clothoid d2.x at start");
		TEST_ASSERT(r.d2.y > 0, "clothoid d2.y positive at start");

		qaws_curve_destroy(crv);
	}

	/* 3D clothoid should fail (2D only) */
	{
		qaws_eval_result_3d r3d;
		qaws_clothoid_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.start_curvature = 0;
		desc.end_curvature = 1;
		desc.length = 1;

		qaws_curve_create_clothoid(&desc, &crv);
		if (crv)
		{
			qaws_status s = qaws_curve_evaluate_3d(crv, 0, QAWS_EVAL_FLAG_POSITION, &r3d);
			TEST_ASSERT(s != QAWS_STATUS_OK, "clothoid 3d eval rejected");
			qaws_curve_destroy(crv);
		}
	}

	/* Error: zero length */
	{
		qaws_clothoid_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.length = 0;
		TEST_ASSERT(qaws_curve_create_clothoid(&desc, &crv) != QAWS_STATUS_OK,
			"clothoid zero length rejected");
	}
}

int test_11_clothoid_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_clothoid();

	printf("11_clothoid: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
