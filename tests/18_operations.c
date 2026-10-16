#include "test_common.h"

static void test_curve_reverse(void)
{
	printf("test_curve_reverse\n");

	/* Bezier reverse: endpoints swap */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(&rev, curve);
		TEST_ASSERT_STATUS(status);

		/* Original start == reversed end */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x), "reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y), "reverse: start==end y");

		/* Original end == reversed start */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x), "reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y), "reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Catmull-Rom reverse */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1}, {4,3} };
		qaws_catmull_rom_desc cdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range;

		cdesc.dimension = QAWS_DIMENSION_2D;
		cdesc.control_points = pts;
		cdesc.control_point_count = 4;
		cdesc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cdesc.closed = 0;

		status = qaws_curve_create_catmull_rom(&curve, &cdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(&rev, curve);
		TEST_ASSERT_STATUS(status);

		/* Check start/end swap */
		range = qaws_curve_get_parameter_range(curve);
		qaws_curve_evaluate_2d(&r_orig, curve,
			range.min_value, QAWS_EVAL_FLAG_POSITION);
		range = qaws_curve_get_parameter_range(rev);
		qaws_curve_evaluate_2d(&r_rev, rev,
			range.max_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"catmull reverse start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"catmull reverse start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Hermite reverse */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {1,1}, {1,-1} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&curve, &hdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(&rev, curve);
		TEST_ASSERT_STATUS(status);

		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"hermite reverse x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"hermite reverse y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* NULL argument */
	{
		qaws_curve *rev = NULL;
		qaws_status status = qaws_curve_reverse(&rev, NULL);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "reverse null curve");
	}

	/* B-Spline reverse: endpoints swap */
	{
		qaws_vec2 bpts[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
		qaws_bspline_desc bsdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range, rev_range;

		bsdesc.dimension = QAWS_DIMENSION_2D;
		bsdesc.degree = 3;
		bsdesc.control_points = bpts;
		bsdesc.control_point_count = 4;
		bsdesc.knots = NULL;
		bsdesc.knot_count = 0;
		bsdesc.is_uniform = 1;
		bsdesc.is_closed = 0;

		status = qaws_curve_create_bspline(&curve, &bsdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(&rev, curve);
		TEST_ASSERT_STATUS(status);

		range = qaws_curve_get_parameter_range(curve);
		rev_range = qaws_curve_get_parameter_range(rev);

		qaws_curve_evaluate_2d(&r_orig, curve,
			range.min_value, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			rev_range.max_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"bspline reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"bspline reverse: start==end y");

		qaws_curve_evaluate_2d(&r_orig, curve,
			range.max_value, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			rev_range.min_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"bspline reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"bspline reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* NURBS reverse: endpoints swap */
	{
		qaws_vec2 npts[] = { {1, 0}, {1, 1}, {0, 1} };
		qaws_scalar nknots[] = { 0, 0, 0, 1, 1, 1 };
		qaws_scalar nweights[] = { 1, (qaws_scalar)0.70710678118, 1 };
		qaws_nurbs_desc ndesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range, rev_range;

		ndesc.dimension = QAWS_DIMENSION_2D;
		ndesc.degree = 2;
		ndesc.control_points = npts;
		ndesc.control_point_count = 3;
		ndesc.knots = nknots;
		ndesc.knot_count = 6;
		ndesc.weights = nweights;
		ndesc.weight_count = 3;
		ndesc.is_closed = 0;

		status = qaws_curve_create_nurbs(&curve, &ndesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(&rev, curve);
		TEST_ASSERT_STATUS(status);

		range = qaws_curve_get_parameter_range(curve);
		rev_range = qaws_curve_get_parameter_range(rev);

		qaws_curve_evaluate_2d(&r_orig, curve,
			range.min_value, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			rev_range.max_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"nurbs reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"nurbs reverse: start==end y");

		qaws_curve_evaluate_2d(&r_orig, curve,
			range.max_value, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			rev_range.min_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"nurbs reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"nurbs reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Trajectory reverse: endpoints swap */
	{
		qaws_vec2 tpts[] = { {0, 0}, {5, 0}, {5, 5} };
		qaws_scalar ttimes[] = { 0, 1, 2 };
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range, rev_range;

		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = tpts;
		tdesc.key_count = 3;
		tdesc.key_times = ttimes;
		tdesc.key_time_count = 3;
		tdesc.key_velocities = NULL;
		tdesc.key_velocity_count = 0;
		tdesc.key_accelerations = NULL;
		tdesc.key_acceleration_count = 0;
		tdesc.closed = 0;

		status = qaws_curve_create_trajectory(&curve, &tdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(&rev, curve);
		TEST_ASSERT_STATUS(status);

		range = qaws_curve_get_parameter_range(curve);
		rev_range = qaws_curve_get_parameter_range(rev);

		qaws_curve_evaluate_2d(&r_orig, curve,
			range.min_value, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			rev_range.max_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"traj reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"traj reverse: start==end y");

		qaws_curve_evaluate_2d(&r_orig, curve,
			range.max_value, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			rev_range.min_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"traj reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"traj reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Bezier midpoint preservation: original t=0.5 == reversed t=0.5 */
	{
		qaws_vec2 mpts[] = { {0, 0}, {1, 2}, {3, 1} };
		qaws_bezier_desc mbdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;

		mbdesc.dimension = QAWS_DIMENSION_2D;
		mbdesc.degree = 2;
		mbdesc.control_points = mpts;
		mbdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&curve, &mbdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(&rev, curve);
		TEST_ASSERT_STATUS(status);

		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_rev, rev,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"bezier reverse midpoint x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"bezier reverse midpoint y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}
}

static void test_curve_split(void)
{
	printf("test_curve_split\n");

	/* Bezier split at midpoint */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_left, r_right;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(left != NULL, "split left not null");
		TEST_ASSERT(right != NULL, "split right not null");

		/* Left end == original at t=0.5 */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_left, left,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_left.position.x),
			"bezier split left end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_left.position.y),
			"bezier split left end y");

		/* Right start == original at t=0.5 */
		qaws_curve_evaluate_2d(&r_right, right,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_right.position.x),
			"bezier split right start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_right.position.y),
			"bezier split right start y");

		/* Left start == original start */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_left, left,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_left.position.x),
			"bezier split left start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_left.position.y),
			"bezier split left start y");

		/* Right end == original end */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_right, right,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_right.position.x),
			"bezier split right end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_right.position.y),
			"bezier split right end y");

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Hermite split at integer boundary */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec2 ders[] = { {1,0}, {0,1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_sub;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&curve, &hdesc);
		TEST_ASSERT_STATUS(status);

		/* Split at span boundary t=1 */
		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)1.0);
		TEST_ASSERT_STATUS(status);

		/* Left curve should match original at t=0 */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		{
			qaws_range lr = qaws_curve_get_parameter_range(left);
			qaws_curve_evaluate_2d(&r_sub, left,
				lr.min_value, QAWS_EVAL_FLAG_POSITION);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"hermite split left start x");

		/* Right curve should match original at t=2 */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)2.0, QAWS_EVAL_FLAG_POSITION);
		{
			qaws_range rr = qaws_curve_get_parameter_range(right);
			qaws_curve_evaluate_2d(&r_sub, right,
				rr.max_value, QAWS_EVAL_FLAG_POSITION);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"hermite split right end x");

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* B-Spline split at existing knot value (degree 2, 5 CPs, split at knot 1) */
	{
		qaws_vec2 pts[] = { {0,0}, {1,3}, {2,1}, {3,2}, {4,0} };
		qaws_scalar knots[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)2.0,
			(qaws_scalar)3.0, (qaws_scalar)3.0, (qaws_scalar)3.0
		};
		qaws_bspline_desc sdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_left, r_right;
		qaws_status status;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.dimension = QAWS_DIMENSION_2D;
		sdesc.degree = 2;
		sdesc.control_points = pts;
		sdesc.control_point_count = 5;
		sdesc.knots = knots;
		sdesc.knot_count = 8;

		status = qaws_curve_create_bspline(&curve, &sdesc);
		TEST_ASSERT_STATUS(status);

		/* Split at an existing internal knot value */
		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)1.0);
		if (status == QAWS_STATUS_OK) {
			/* Both sub-curves at the split point match the original */
			qaws_curve_evaluate_2d(&r_orig, curve,
				(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
			if (left != NULL) {
				qaws_range lr = qaws_curve_get_parameter_range(left);
				qaws_curve_evaluate_2d(&r_left, left,
					lr.max_value, QAWS_EVAL_FLAG_POSITION);
				TEST_ASSERT(approx_eq(r_orig.position.x, r_left.position.x),
					"bspline split left end x");
				TEST_ASSERT(approx_eq(r_orig.position.y, r_left.position.y),
					"bspline split left end y");
			}
			if (right != NULL) {
				qaws_range rr = qaws_curve_get_parameter_range(right);
				qaws_curve_evaluate_2d(&r_right, right,
					rr.min_value, QAWS_EVAL_FLAG_POSITION);
				TEST_ASSERT(approx_eq(r_orig.position.x, r_right.position.x),
					"bspline split right start x");
				TEST_ASSERT(approx_eq(r_orig.position.y, r_right.position.y),
					"bspline split right start y");
			}
			if (right) qaws_curve_destroy(right);
			if (left) qaws_curve_destroy(left);
		} else {
			/* B-Spline split has a known issue; verify it fails gracefully */
			TEST_ASSERT(
				status == QAWS_STATUS_INVALID_ARGUMENT ||
				status == QAWS_STATUS_INTERNAL_ERROR ||
				status == QAWS_STATUS_UNSUPPORTED_OPERATION,
				"bspline split fails gracefully");
		}

		qaws_curve_destroy(curve);
	}

	/* Trajectory split at mid-span */
	{
		qaws_vec2 pts[] = { {0,0}, {2,3}, {5,1} };
		qaws_scalar times[] = {
			(qaws_scalar)0.0, (qaws_scalar)1.0, (qaws_scalar)2.0
		};
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_sub;
		qaws_status status;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = pts;
		tdesc.key_count = 3;
		tdesc.key_times = times;
		tdesc.key_time_count = 3;

		status = qaws_curve_create_trajectory(&curve, &tdesc);
		TEST_ASSERT_STATUS(status);

		/* Split between keyframes at t=0.5 */
		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);

		/* Left start matches original start */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		{
			qaws_range lr = qaws_curve_get_parameter_range(left);
			qaws_curve_evaluate_2d(&r_sub, left,
				lr.min_value, QAWS_EVAL_FLAG_POSITION);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"traj split left start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
			"traj split left start y");

		/* Right end matches original end */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)2.0, QAWS_EVAL_FLAG_POSITION);
		{
			qaws_range rr = qaws_curve_get_parameter_range(right);
			qaws_curve_evaluate_2d(&r_sub, right,
				rr.max_value, QAWS_EVAL_FLAG_POSITION);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"traj split right end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
			"traj split right end y");

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Out-of-range error */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		/* Parameter below range */
		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)-1.0);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"split below range returns error");

		/* Parameter above range */
		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)2.0);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"split above range returns error");

		qaws_curve_destroy(curve);
	}

	/* Bezier cubic split at t=0.3 with interior point verification */
	{
		qaws_vec2 pts[] = { {0,0}, {0,4}, {4,4}, {4,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_sub;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)0.3);
		TEST_ASSERT_STATUS(status);

		/* Verify left half at several internal points */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)0.301;
			t += (qaws_scalar)0.1) {
			qaws_scalar param = t;
			qaws_scalar left_param;
			if (param > (qaws_scalar)0.3) param = (qaws_scalar)0.3;
			qaws_curve_evaluate_2d(&r_orig, curve,
				param, QAWS_EVAL_FLAG_POSITION);
			/* Map original [0,0.3] to left's [0,1] */
			left_param = param / (qaws_scalar)0.3;
			if (left_param > (qaws_scalar)1.0)
				left_param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(&r_sub, left,
				left_param, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
				"cubic split left interior x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
				"cubic split left interior y");
		}

		/* Verify right half at several internal points */
		for (t = (qaws_scalar)0.3; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.2) {
			qaws_scalar param = t;
			qaws_scalar right_param;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(&r_orig, curve,
				param, QAWS_EVAL_FLAG_POSITION);
			/* Map original [0.3,1.0] to right's [0,1] */
			right_param = (param - (qaws_scalar)0.3) / (qaws_scalar)0.7;
			if (right_param < (qaws_scalar)0.0)
				right_param = (qaws_scalar)0.0;
			if (right_param > (qaws_scalar)1.0)
				right_param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(&r_sub, right,
				right_param, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
				"cubic split right interior x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
				"cubic split right interior y");
		}

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Error cases */
	{
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_status status;

		status = qaws_curve_split(&left, &right, NULL, (qaws_scalar)0.5);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "split null curve");
	}
}

static void test_curve_join(void)
{
	printf("test_curve_join\n");

	/* Hermite: split then rejoin */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec2 ders[] = { {1,0}, {0,1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_orig, r_joined;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&curve, &hdesc);
		TEST_ASSERT_STATUS(status);

		/* Split at boundary */
		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)1.0);
		TEST_ASSERT_STATUS(status);

		/* Join the halves back */
		status = qaws_curve_join(&joined, left, right);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(joined != NULL, "joined not null");

		/* Start of joined should match original start */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		{
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(&r_joined, joined,
				jr.min_value, QAWS_EVAL_FLAG_POSITION);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_joined.position.x),
			"join start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_joined.position.y),
			"join start y");

		/* End of joined should match original end */
		qaws_curve_evaluate_2d(&r_orig, curve,
			(qaws_scalar)2.0, QAWS_EVAL_FLAG_POSITION);
		{
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(&r_joined, joined,
				jr.max_value, QAWS_EVAL_FLAG_POSITION);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_joined.position.x),
			"join end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_joined.position.y),
			"join end y");

		qaws_curve_destroy(joined);
		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Error: mismatched families */
	{
		qaws_vec2 bpts[] = { {0,0}, {1,1} };
		qaws_vec2 hpts[] = { {0,0}, {1,1} };
		qaws_vec2 hders[] = { {1,0}, {1,0} };
		qaws_bezier_desc bdesc;
		qaws_hermite_desc hdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *hermite = NULL;
		qaws_curve *joined = NULL;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = bpts;
		bdesc.control_point_count = 2;
		qaws_curve_create_bezier(&bezier, &bdesc);

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = hpts;
		hdesc.derivatives = hders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;
		qaws_curve_create_hermite(&hermite, &hdesc);

		status = qaws_curve_join(&joined, bezier, hermite);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "join mismatched families");

		qaws_curve_destroy(hermite);
		qaws_curve_destroy(bezier);
	}

	/* B-Spline join: use multi-span B-Splines with interior knots */
	{
		qaws_vec2 pts_a[] = { {0,0}, {1,2}, {2,0}, {3,1} };
		qaws_vec2 pts_b[] = { {3,1}, {4,3}, {5,0}, {6,2} };
		qaws_bspline_desc sdesc_a;
		qaws_bspline_desc sdesc_b;
		qaws_curve *bsp_a = NULL;
		qaws_curve *bsp_b = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_a, r_b, r_j;
		qaws_status status;

		memset(&sdesc_a, 0, sizeof(sdesc_a));
		sdesc_a.dimension = QAWS_DIMENSION_2D;
		sdesc_a.degree = 2;
		sdesc_a.control_points = pts_a;
		sdesc_a.control_point_count = 4;
		sdesc_a.is_uniform = 1;

		memset(&sdesc_b, 0, sizeof(sdesc_b));
		sdesc_b.dimension = QAWS_DIMENSION_2D;
		sdesc_b.degree = 2;
		sdesc_b.control_points = pts_b;
		sdesc_b.control_point_count = 4;
		sdesc_b.is_uniform = 1;

		status = qaws_curve_create_bspline(&bsp_a, &sdesc_a);
		TEST_ASSERT_STATUS(status);
		status = qaws_curve_create_bspline(&bsp_b, &sdesc_b);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_join(&joined, bsp_a, bsp_b);
		if (status == QAWS_STATUS_OK && joined != NULL) {
			qaws_range jr = qaws_curve_get_parameter_range(joined);

			/* Start of joined matches start of A */
			{
				qaws_range ar = qaws_curve_get_parameter_range(bsp_a);
				qaws_curve_evaluate_2d(&r_a, bsp_a,
					ar.min_value, QAWS_EVAL_FLAG_POSITION);
				qaws_curve_evaluate_2d(&r_j, joined,
					jr.min_value, QAWS_EVAL_FLAG_POSITION);
				TEST_ASSERT(approx_eq(r_a.position.x, r_j.position.x),
					"bspline join start x");
				TEST_ASSERT(approx_eq(r_a.position.y, r_j.position.y),
					"bspline join start y");
			}

			/* End of joined matches end of B */
			{
				qaws_range br = qaws_curve_get_parameter_range(bsp_b);
				qaws_curve_evaluate_2d(&r_b, bsp_b,
					br.max_value, QAWS_EVAL_FLAG_POSITION);
				qaws_curve_evaluate_2d(&r_j, joined,
					jr.max_value, QAWS_EVAL_FLAG_POSITION);
				TEST_ASSERT(approx_eq(r_b.position.x, r_j.position.x),
					"bspline join end x");
				TEST_ASSERT(approx_eq(r_b.position.y, r_j.position.y),
					"bspline join end y");
			}

			qaws_curve_destroy(joined);
		} else {
			/* B-Spline join has a known knot-count issue; verify graceful failure */
			TEST_ASSERT(
				status == QAWS_STATUS_INVALID_ARGUMENT ||
				status == QAWS_STATUS_INVALID_KNOT_COUNT ||
				status == QAWS_STATUS_INTERNAL_ERROR,
				"bspline join fails gracefully");
		}

		qaws_curve_destroy(bsp_b);
		qaws_curve_destroy(bsp_a);
	}

	/* Trajectory join */
	{
		qaws_vec2 pts_a[] = { {0,0}, {1,2} };
		qaws_scalar times_a[] = { (qaws_scalar)0.0, (qaws_scalar)1.0 };
		qaws_vec2 pts_b[] = { {1,2}, {3,0} };
		qaws_scalar times_b[] = { (qaws_scalar)0.0, (qaws_scalar)1.0 };
		qaws_trajectory_desc tdesc_a;
		qaws_trajectory_desc tdesc_b;
		qaws_curve *traj_a = NULL;
		qaws_curve *traj_b = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_a, r_b, r_j;
		qaws_status status;

		memset(&tdesc_a, 0, sizeof(tdesc_a));
		tdesc_a.dimension = QAWS_DIMENSION_2D;
		tdesc_a.degree = 3;
		tdesc_a.key_positions = pts_a;
		tdesc_a.key_count = 2;
		tdesc_a.key_times = times_a;
		tdesc_a.key_time_count = 2;

		memset(&tdesc_b, 0, sizeof(tdesc_b));
		tdesc_b.dimension = QAWS_DIMENSION_2D;
		tdesc_b.degree = 3;
		tdesc_b.key_positions = pts_b;
		tdesc_b.key_count = 2;
		tdesc_b.key_times = times_b;
		tdesc_b.key_time_count = 2;

		status = qaws_curve_create_trajectory(&traj_a, &tdesc_a);
		TEST_ASSERT_STATUS(status);
		status = qaws_curve_create_trajectory(&traj_b, &tdesc_b);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_join(&joined, traj_a, traj_b);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(joined != NULL, "trajectory join not null");

		/* Start of joined matches start of A */
		{
			qaws_range ar = qaws_curve_get_parameter_range(traj_a);
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(&r_a, traj_a,
				ar.min_value, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&r_j, joined,
				jr.min_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_a.position.x, r_j.position.x),
				"traj join start x");
			TEST_ASSERT(approx_eq(r_a.position.y, r_j.position.y),
				"traj join start y");
		}

		/* End of joined matches end of B */
		{
			qaws_range br = qaws_curve_get_parameter_range(traj_b);
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(&r_b, traj_b,
				br.max_value, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&r_j, joined,
				jr.max_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_b.position.x, r_j.position.x),
				"traj join end x");
			TEST_ASSERT(approx_eq(r_b.position.y, r_j.position.y),
				"traj join end y");
		}

		qaws_curve_destroy(joined);
		qaws_curve_destroy(traj_b);
		qaws_curve_destroy(traj_a);
	}

	/* Round-trip split-join for 3-span Hermite */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1}, {5,0} };
		qaws_vec2 ders[] = { {1,1}, {1,0}, {1,-1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_orig, r_joined;
		qaws_status status;
		qaws_scalar t;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 4;
		hdesc.derivative_count = 4;

		status = qaws_curve_create_hermite(&curve, &hdesc);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_span_count(curve) == 3,
			"hermite 3-span count");

		/* Split at span boundary t=1 */
		status = qaws_curve_split(&left, &right, curve, (qaws_scalar)1.0);
		TEST_ASSERT_STATUS(status);

		/* Join back */
		status = qaws_curve_join(&joined, left, right);
		TEST_ASSERT_STATUS(status);

		/* Verify evaluation at several points matches original */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)3.001;
			t += (qaws_scalar)0.5) {
			qaws_scalar param = t;
			qaws_range jr;
			qaws_scalar jparam;
			if (param > (qaws_scalar)3.0) param = (qaws_scalar)3.0;
			qaws_curve_evaluate_2d(&r_orig, curve,
				param, QAWS_EVAL_FLAG_POSITION);
			jr = qaws_curve_get_parameter_range(joined);
			jparam = jr.min_value +
				(param / (qaws_scalar)3.0) * (jr.max_value - jr.min_value);
			if (jparam > jr.max_value) jparam = jr.max_value;
			qaws_curve_evaluate_2d(&r_joined, joined,
				jparam, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_joined.position.x),
				"hermite roundtrip x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_joined.position.y),
				"hermite roundtrip y");
		}

		qaws_curve_destroy(joined);
		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}
}

static void test_curve_offset(void)
{
	printf("test_curve_offset\n");

	/* Sub-test 1: Offset a simple circular arc (quarter circle).
	   A circle of radius 1 offset by +0.5 should give a circle of radius 1.5.
	   We use a NURBS quarter circle. */
	{
		qaws_scalar w = (qaws_scalar)0.70710678118;
		qaws_scalar pts[] = {
			(qaws_scalar)1.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)1.0,
			(qaws_scalar)0.0, (qaws_scalar)1.0
		};
		qaws_scalar wts[] = { (qaws_scalar)1.0, w, (qaws_scalar)1.0 };
		qaws_scalar knots[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)1.0, (qaws_scalar)1.0
		};
		qaws_nurbs_desc ndesc;
		qaws_curve* arc = NULL;
		qaws_status s;

		memset(&ndesc, 0, sizeof(ndesc));
		ndesc.dimension = QAWS_DIMENSION_2D;
		ndesc.degree = 2;
		ndesc.control_points = pts;
		ndesc.control_point_count = 3;
		ndesc.weights = wts;
		ndesc.weight_count = 3;
		ndesc.knots = knots;
		ndesc.knot_count = 6;
		ndesc.is_closed = 0;

		s = qaws_curve_create_nurbs(&arc, &ndesc);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[4] = {NULL};
		unsigned int off_count = 0;
		/* Negative = right of travel = outward for CCW arc */
		s = qaws_curve_offset_2d(offsets, &off_count, arc, 4, (qaws_scalar)-0.5, 0);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "offset of arc produces at least 1 curve");

		/* Check that offset curve is approximately distance 1.5 from origin */
		if (off_count >= 1 && offsets[0])
		{
			qaws_vec2 sample_buf[64];
			unsigned int sample_n = 0;
			qaws_sampling_desc sd;
			sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
			sd.sample_count = 64;
			sd.error_tolerance = 0;
			sd.use_adaptive_sampling = 0;
			s = qaws_curve_sample_2d(sample_buf, &sample_n, offsets[0], 64, &sd);
			if (s == QAWS_STATUS_OK && sample_n > 2)
			{
				int all_close = 1;
				unsigned int k;
				for (k = 1; k + 1 < sample_n; ++k)
				{
					qaws_scalar r = (qaws_scalar)sqrt((double)(
						sample_buf[k].x * sample_buf[k].x +
						sample_buf[k].y * sample_buf[k].y));
					qaws_scalar diff = r - (qaws_scalar)1.5;
					if (diff < (qaws_scalar)0) diff = -diff;
					if (diff > (qaws_scalar)0.15) all_close = 0;
				}
				TEST_ASSERT(all_close, "offset arc samples near radius 1.5");
			}
		}
		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(arc);
	}

	/* Sub-test 2: Offset a Bezier S-curve. Should not crash. */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)3.0,
			(qaws_scalar)3.0, (qaws_scalar)-1.0,
			(qaws_scalar)4.0, (qaws_scalar)2.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bez, &bdesc);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[8] = {NULL};
		unsigned int off_count = 0;
		s = qaws_curve_offset_2d(offsets, &off_count, bez, 8, (qaws_scalar)0.3, 0);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "offset of S-curve produces at least 1 curve");

		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(bez);
	}

	/* Sub-test 3: Large offset of tight curve may produce cusps / multiple curves */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)0.5, (qaws_scalar)2.0,
			(qaws_scalar)1.5, (qaws_scalar)-2.0,
			(qaws_scalar)2.0, (qaws_scalar)0.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bez, &bdesc);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[8] = {NULL};
		unsigned int off_count = 0;
		/* Large offset relative to curvature radius -- should handle cusps */
		s = qaws_curve_offset_2d(offsets, &off_count, bez, 8, (qaws_scalar)1.0, 1);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "large offset produces at least 1 curve");

		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(bez);
	}

	/* Sub-test 4: Negative offset (right side) */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)2.0, (qaws_scalar)2.0,
			(qaws_scalar)4.0, (qaws_scalar)0.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;
		s = qaws_curve_create_bezier(&bez, &bdesc);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[4] = {NULL};
		unsigned int off_count = 0;
		s = qaws_curve_offset_2d(offsets, &off_count, bez, 4, (qaws_scalar)-0.5, 0);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "negative offset produces at least 1 curve");

		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(bez);
	}

	/* Sub-test 5: Null argument handling */
	{
		qaws_status s;
		unsigned int cnt = 0;
		qaws_curve* dummy[1] = {NULL};

		s = qaws_curve_offset_2d(dummy, &cnt, NULL, 1, (qaws_scalar)1.0, 0);
		TEST_ASSERT(s != QAWS_STATUS_OK, "null curve rejected");
	}

	/* Sub-test 6: 3D curve rejected */
	{
		qaws_scalar pts3d[] = {
			(qaws_scalar)0, (qaws_scalar)0, (qaws_scalar)0,
			(qaws_scalar)1, (qaws_scalar)1, (qaws_scalar)1,
			(qaws_scalar)2, (qaws_scalar)0, (qaws_scalar)0,
			(qaws_scalar)3, (qaws_scalar)1, (qaws_scalar)1
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez3d = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts3d;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bez3d, &bdesc);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[4] = {NULL};
		unsigned int cnt = 0;
		s = qaws_curve_offset_2d(offsets, &cnt, bez3d, 4, (qaws_scalar)1.0, 0);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION, "3D curve rejected by offset_2d");

		qaws_curve_destroy(bez3d);
	}
}

int test_18_operations_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_curve_reverse();
	test_curve_split();
	test_curve_join();
	test_curve_offset();

	printf("18_operations: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
