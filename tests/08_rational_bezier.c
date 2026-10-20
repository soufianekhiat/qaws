#include "test_common.h"

static void test_rational_bezier(void)
{
	printf("test_rational_bezier\n");

	/* Quadratic rational Bezier: exact circle arc (quarter circle) */
	{
		qaws_scalar cp[] = { 1, 0,   1, 1,   0, 1 };
		qaws_scalar w[] = { 1, (qaws_scalar)0.70710678118, 1 }; /* 1, 1/sqrt(2), 1 */
		qaws_rational_bezier_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;
		qaws_scalar dist;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = cp;
		desc.control_point_count = 3;
		desc.weights = w;
		desc.weight_count = 3;

		TEST_ASSERT_STATUS(qaws_curve_create_rational_bezier(&crv, &desc));
		TEST_ASSERT(crv != NULL, "rational bezier created");
		TEST_ASSERT(qaws_curve_get_kind(crv) == QAWS_CURVE_KIND_RATIONAL_BEZIER,
			"rational bezier kind");
		TEST_ASSERT(qaws_curve_get_degree(crv) == 2, "rational bezier degree");

		/* At t=0: should be (1,0) */
		qaws_curve_evaluate_2d(&r, crv, 0, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		TEST_ASSERT(approx_eq(r.position.x, 1) && approx_eq(r.position.y, 0),
			"rational bezier start point");
		/* D1 at start should be tangent to circle: direction (0,1) */
		TEST_ASSERT(r.d1.x < (qaws_scalar)0.01, "rbez d1.x near zero at start");
		TEST_ASSERT(r.d1.y > (qaws_scalar)0.1, "rbez d1.y positive at start");

		/* At t=1: should be (0,1) */
		qaws_curve_evaluate_2d(&r, crv, 1, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		TEST_ASSERT(approx_eq(r.position.x, 0) && approx_eq(r.position.y, 1),
			"rational bezier end point");
		/* D1 at end: tangent direction (-1, 0) */
		TEST_ASSERT(r.d1.x < (qaws_scalar)-0.1, "rbez d1.x negative at end");

		/* At t=0.5: should be on unit circle */
		qaws_curve_evaluate_2d(&r, crv,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		dist = (qaws_scalar)sqrt((double)(r.position.x * r.position.x + r.position.y * r.position.y));
		TEST_ASSERT(approx_eq(dist, 1), "rational bezier midpoint on circle");
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D2, "rbez d2 valid");

		/* Multiple points along the arc should all be on unit circle */
		{
			unsigned int i;
			int all_on_circle = 1;
			for (i = 1; i < 10; ++i)
			{
				qaws_scalar t = (qaws_scalar)i / (qaws_scalar)10;
				qaws_curve_evaluate_2d(&r, crv, t, QAWS_EVAL_FLAG_POSITION);
				dist = (qaws_scalar)sqrt((double)(r.position.x * r.position.x + r.position.y * r.position.y));
				if (!approx_eq(dist, 1)) all_on_circle = 0;
			}
			TEST_ASSERT(all_on_circle, "rbez all samples on unit circle");
		}

		/* is_rational should return 1 */
		TEST_ASSERT(qaws_curve_is_rational(crv) == 1, "rational bezier is_rational");
		TEST_ASSERT(qaws_curve_is_closed(crv) == 0, "rational bezier not closed");

		qaws_curve_destroy(crv);
	}

	/* Error handling: NULL desc */
	{
		qaws_curve *crv = NULL;
		TEST_ASSERT(qaws_curve_create_rational_bezier(&crv, NULL) != QAWS_STATUS_OK,
			"rbez null desc rejected");
	}

	/* 3D rational Bezier */
	{
		qaws_scalar cp3d[] = { 1,0,0,  1,1,0,  0,1,0 };
		qaws_scalar w3d[] = { 1, (qaws_scalar)0.70710678118, 1 };
		qaws_rational_bezier_desc desc3d;
		qaws_curve *crv3d = NULL;
		qaws_eval_result_3d r3d;

		memset(&desc3d, 0, sizeof(desc3d));
		desc3d.dimension = QAWS_DIMENSION_3D;
		desc3d.degree = 2;
		desc3d.control_points = cp3d;
		desc3d.control_point_count = 3;
		desc3d.weights = w3d;
		desc3d.weight_count = 3;

		TEST_ASSERT_STATUS(qaws_curve_create_rational_bezier(&crv3d, &desc3d));
		TEST_ASSERT(qaws_curve_get_dimension(crv3d) == QAWS_DIMENSION_3D,
			"rbez 3d dimension");

		qaws_curve_evaluate_3d(&r3d, crv3d,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		TEST_ASSERT(r3d.valid_flags & QAWS_EVAL_FLAG_POSITION, "rational bezier 3d pos");
		TEST_ASSERT(r3d.valid_flags & QAWS_EVAL_FLAG_D1, "rational bezier 3d d1");
		TEST_ASSERT(r3d.valid_flags & QAWS_EVAL_FLAG_D2, "rational bezier 3d d2");
		/* z should be 0 since all points are in z=0 plane */
		TEST_ASSERT(approx_eq(r3d.position.z, 0), "rbez 3d z=0");

		qaws_curve_destroy(crv3d);
	}
}

int test_08_rational_bezier_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_rational_bezier();

	printf("08_rational_bezier: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
