#include "test_common.h"

static void test_subdivision(void)
{
	printf("test_subdivision\n");

	/* Open Chaikin subdivision */
	{
		qaws_scalar cp[] = { 0,0,  1,2,  3,1,  4,0 };
		qaws_subdivision_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.scheme = QAWS_SUBDIVISION_CHAIKIN;
		desc.control_points = cp;
		desc.control_point_count = 4;
		desc.closed = 0;
		desc.refinement_levels = 5;

		TEST_ASSERT_STATUS(qaws_curve_create_subdivision(&crv, &desc));
		TEST_ASSERT(crv != NULL, "subdivision chaikin created");
		TEST_ASSERT(qaws_curve_get_kind(crv) == QAWS_CURVE_KIND_SUBDIVISION,
			"subdivision kind");
		TEST_ASSERT(qaws_curve_get_continuity(crv) == QAWS_CONTINUITY_C1,
			"chaikin continuity C1");
		TEST_ASSERT(qaws_curve_is_closed(crv) == 0, "subdivision open not closed");

		/* Endpoints should be near original endpoints */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv);
			qaws_curve_evaluate_2d(&r, crv, rng.min_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(r.position.x >= -1 && r.position.x <= 1, "subdivision start x plausible");

			qaws_curve_evaluate_2d(&r, crv, rng.max_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(r.position.x >= 3 && r.position.x <= 5, "subdivision end x plausible");
		}

		/* Should support all derivatives */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv);
			qaws_scalar mid = (rng.min_value + rng.max_value) * (qaws_scalar)0.5;
			qaws_curve_evaluate_2d(&r, crv,
				mid, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3);
			TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D1, "subdivision d1 valid");
			TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D2, "subdivision d2 valid");
			TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D3, "subdivision d3 valid");
			/* D1 should point generally in +x direction */
			TEST_ASSERT(r.d1.x > 0, "subdivision d1.x positive at mid");
		}

		/* Smoothness: sample the curve and check it has no big jumps */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv);
			qaws_eval_result_2d prev, curr;
			unsigned int i;
			int smooth = 1;
			qaws_curve_evaluate_2d(&prev, crv, rng.min_value, QAWS_EVAL_FLAG_POSITION);
			for (i = 1; i <= 50; ++i)
			{
				qaws_scalar t = rng.min_value + (rng.max_value - rng.min_value) * (qaws_scalar)i / (qaws_scalar)50;
				qaws_scalar dx, dy, jump;
				qaws_curve_evaluate_2d(&curr, crv, t, QAWS_EVAL_FLAG_POSITION);
				dx = curr.position.x - prev.position.x;
				dy = curr.position.y - prev.position.y;
				jump = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));
				if (jump > (qaws_scalar)0.5) smooth = 0;
				prev = curr;
			}
			TEST_ASSERT(smooth, "subdivision chaikin smooth curve");
		}

		qaws_curve_destroy(crv);
	}

	/* Closed Lane-Riesenfeld (cubic) */
	{
		qaws_scalar cp[] = { 0,0,  2,0,  2,2,  0,2 };
		qaws_subdivision_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r_start, r_end;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		desc.closed = 1;
		desc.refinement_levels = 5;

		TEST_ASSERT_STATUS(qaws_curve_create_subdivision(&crv, &desc));
		TEST_ASSERT(qaws_curve_is_closed(crv) == 1, "subdivision closed");
		TEST_ASSERT(qaws_curve_get_continuity(crv) == QAWS_CONTINUITY_C2,
			"subdivision lr3 continuity");

		/* Closed: start and end should match */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv);
			qaws_curve_evaluate_2d(&r_start, crv, rng.min_value, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&r_end, crv, rng.max_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(
				approx_eq(r_start.position.x, r_end.position.x) &&
				approx_eq(r_start.position.y, r_end.position.y),
				"subdivision closed start==end");
		}

		qaws_curve_destroy(crv);
	}

	/* Lane-Riesenfeld order 4 */
	{
		qaws_scalar cp[] = { 0,0,  1,2,  3,1,  4,0 };
		qaws_subdivision_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_4;
		desc.control_points = cp;
		desc.control_point_count = 4;
		desc.closed = 0;
		desc.refinement_levels = 5;

		TEST_ASSERT_STATUS(qaws_curve_create_subdivision(&crv, &desc));
		TEST_ASSERT(qaws_curve_get_continuity(crv) == QAWS_CONTINUITY_C3,
			"subdivision lr4 continuity C3");
		qaws_curve_destroy(crv);
	}

	/* 3D subdivision */
	{
		qaws_scalar cp3d[] = { 0,0,0,  1,1,1,  2,0,2,  3,1,0 };
		qaws_subdivision_desc desc3d;
		qaws_curve *crv3d = NULL;
		qaws_eval_result_3d r3d;

		memset(&desc3d, 0, sizeof(desc3d));
		desc3d.dimension = QAWS_DIMENSION_3D;
		desc3d.scheme = QAWS_SUBDIVISION_CHAIKIN;
		desc3d.control_points = cp3d;
		desc3d.control_point_count = 4;
		desc3d.closed = 0;
		desc3d.refinement_levels = 4;

		TEST_ASSERT_STATUS(qaws_curve_create_subdivision(&crv3d, &desc3d));
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv3d);
			qaws_scalar mid = (rng.min_value + rng.max_value) * (qaws_scalar)0.5;
			qaws_curve_evaluate_3d(&r3d, crv3d,
				mid, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
			TEST_ASSERT(r3d.valid_flags & QAWS_EVAL_FLAG_POSITION, "subdivision 3d pos");
			TEST_ASSERT(r3d.valid_flags & QAWS_EVAL_FLAG_D1, "subdivision 3d d1");
		}

		qaws_curve_destroy(crv3d);
	}
}

int test_12_subdivision_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_subdivision();

	printf("12_subdivision: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
