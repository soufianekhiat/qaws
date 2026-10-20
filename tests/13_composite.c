#include "test_common.h"

static void test_composite(void)
{
	printf("test_composite\n");

	/* Chain two Bezier segments into a composite */
	{
		qaws_scalar cp1[] = { 0,0, 1,1 };
		qaws_scalar cp2[] = { 1,1, 2,0 };
		qaws_bezier_desc bd1, bd2;
		qaws_curve *seg1 = NULL, *seg2 = NULL, *comp = NULL;
		qaws_composite_desc cd;
		qaws_curve *segs[2];
		qaws_eval_result_2d r;

		memset(&bd1, 0, sizeof(bd1));
		bd1.dimension = QAWS_DIMENSION_2D;
		bd1.degree = 1;
		bd1.control_points = cp1;
		bd1.control_point_count = 2;

		memset(&bd2, 0, sizeof(bd2));
		bd2.dimension = QAWS_DIMENSION_2D;
		bd2.degree = 1;
		bd2.control_points = cp2;
		bd2.control_point_count = 2;

		TEST_ASSERT_STATUS(qaws_curve_create_bezier(&seg1, &bd1));
		TEST_ASSERT_STATUS(qaws_curve_create_bezier(&seg2, &bd2));

		segs[0] = seg1;
		segs[1] = seg2;

		memset(&cd, 0, sizeof(cd));
		cd.dimension = QAWS_DIMENSION_2D;
		cd.segments = segs;
		cd.segment_count = 2;

		TEST_ASSERT_STATUS(qaws_curve_create_composite(&comp, &cd));
		TEST_ASSERT(comp != NULL, "composite created");
		TEST_ASSERT(qaws_curve_get_kind(comp) == QAWS_CURVE_KIND_COMPOSITE,
			"composite kind");
		TEST_ASSERT(qaws_curve_get_span_count(comp) == 2, "composite span_count");

		/* Range should be [0, 2] */
		{
			qaws_range rng = qaws_curve_get_parameter_range(comp);
			TEST_ASSERT(approx_eq(rng.min_value, 0), "composite range min");
			TEST_ASSERT(approx_eq(rng.max_value, 2), "composite range max");
		}

		/* t=0: start of first segment = (0,0) */
		qaws_curve_evaluate_2d(&r, comp, 0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r.position.x, 0) && approx_eq(r.position.y, 0),
			"composite start");

		/* t=1: end of first segment / start of second = (1,1) */
		qaws_curve_evaluate_2d(&r, comp, 1, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r.position.x, 1) && approx_eq(r.position.y, 1),
			"composite junction");

		/* t=2: end = (2,0) */
		qaws_curve_evaluate_2d(&r, comp, 2, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r.position.x, 2) && approx_eq(r.position.y, 0),
			"composite end");

		/* Derivatives: linear segments have constant D1 */
		qaws_curve_evaluate_2d(&r, comp,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D1, "composite d1 valid");
		/* D1 of first segment (0,0)->(1,1) should be (1,1) scaled by chain rule */
		TEST_ASSERT(approx_eq(r.d1.x, 1) && approx_eq(r.d1.y, 1),
			"composite d1 first seg");

		qaws_curve_evaluate_2d(&r, comp,
			(qaws_scalar)1.5, QAWS_EVAL_FLAG_D1);
		/* D1 of second segment (1,1)->(2,0) should be (1,-1) */
		TEST_ASSERT(approx_eq(r.d1.x, 1) && approx_eq(r.d1.y, -1),
			"composite d1 second seg");

		/* Not closed (start != end) */
		TEST_ASSERT(qaws_curve_get_continuity(comp) == QAWS_CONTINUITY_C0,
			"composite continuity C0");

		qaws_curve_destroy(comp); /* destroys seg1, seg2 too */
	}

	/* Error handling: NULL desc */
	{
		qaws_curve *crv = NULL;
		TEST_ASSERT(qaws_curve_create_composite(&crv, NULL) != QAWS_STATUS_OK,
			"composite null desc rejected");
	}
}

int test_13_composite_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_composite();

	printf("13_composite: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
