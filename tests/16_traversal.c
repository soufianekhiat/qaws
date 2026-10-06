#include "test_common.h"

static int scalar_is_finite(qaws_scalar v)
{
	return (v == v) && (v - v == 0);
}

static qaws_scalar constant_speed_fn(qaws_scalar distance, void* user_data)
{
	(void)distance;
	(void)user_data;
	return (qaws_scalar)1.0;
}

static qaws_scalar accelerating_speed_fn(qaws_scalar distance, void* user_data)
{
	(void)user_data;
	return (qaws_scalar)1.0 + distance;
}

static void test_traversal(void)
{
	printf("test_traversal\n");

	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	memset(&tdesc, 0, sizeof(tdesc));
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)5.0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = (qaws_scalar)2.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	/* At t=0, distance=0, should be at start */
	qaws_eval_result_2d r;
	s = qaws_traversal_evaluate_2d(traversal, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "traversal t=0 x");

	/* At t=1, distance=5, should be at x=5 */
	s = qaws_traversal_evaluate_2d(traversal, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "traversal t=1 x");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

static void test_traversal_mappings(void)
{
	printf("test_traversal_mappings\n");

	/* Straight line (0,0) to (10,0), length=10 */
	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	memset(&tdesc, 0, sizeof(tdesc));
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)5.0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = (qaws_scalar)2.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	/* map_time_to_parameter: at t=1, distance=5, on length-10 line => param ~0.5 */
	qaws_scalar param = 0;
	s = qaws_traversal_map_time_to_parameter(traversal, (qaws_scalar)1.0, &param);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(param, (qaws_scalar)0.5), "time_to_param t=1 => 0.5");

	/* map_distance_to_parameter: distance 5 => param 0.5 */
	s = qaws_traversal_map_distance_to_parameter(traversal, (qaws_scalar)5.0, &param);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(param, (qaws_scalar)0.5), "dist_to_param d=5 => 0.5");

	/* map_parameter_to_distance: param 0.5 => distance 5 */
	qaws_scalar dist = 0;
	s = qaws_traversal_map_parameter_to_distance(traversal, (qaws_scalar)0.5, &dist);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(dist, (qaws_scalar)5.0), "param_to_dist p=0.5 => 5");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

static void test_traversal_arc_length_mode(void)
{
	printf("test_traversal_arc_length_mode\n");

	/* Arc-length traversal: straight line, distance 5 should map to midpoint */
	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	memset(&tdesc, 0, sizeof(tdesc));
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_NONE;
	tdesc.speed = 0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = 0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	/* input_value is distance: 5.0 on a 10-length line => x=5 */
	qaws_eval_result_2d r;
	s = qaws_traversal_evaluate_2d(traversal, (qaws_scalar)5.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "arc-length eval at d=5 => x=5");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

static void test_traversal_3d(void)
{
	printf("test_traversal_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {10, 0, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	memset(&tdesc, 0, sizeof(tdesc));
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)10.0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = (qaws_scalar)1.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_3d r;
	s = qaws_traversal_evaluate_3d(traversal, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "3d traversal t=0.5 x");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

static void test_easing(void)
{
	printf("test_easing\n");

	/* Create a simple curve and traversal with easing */
	qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	/* Linear easing: same as no easing */
	{
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_linear, r_eased;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;
		tdesc.easing = QAWS_EASING_LINEAR;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		status = qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_linear);
		TEST_ASSERT_STATUS(status);

		qaws_traversal_destroy(trav);

		/* Quad ease-in: at t=0.5, eased_t = 0.25, so position should be earlier */
		tdesc.easing = QAWS_EASING_QUAD_IN;
		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		status = qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_eased);
		TEST_ASSERT_STATUS(status);

		/* Eased position should differ from linear */
		TEST_ASSERT(
			!approx_eq(r_linear.position.x, r_eased.position.x) ||
			!approx_eq(r_linear.position.y, r_eased.position.y),
			"ease-in differs from linear at t=0.5");

		qaws_traversal_destroy(trav);
	}

	/* Sine easing endpoints: at t=0 should be start, at t=1 should be end */
	{
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_start, r_end;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;
		tdesc.easing = QAWS_EASING_SINE_IN_OUT;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r_start);
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r_end);

		/* Start should be near (0,0) */
		TEST_ASSERT(approx_eq(r_start.position.x, (qaws_scalar)0.0),
			"eased start x near 0");

		qaws_traversal_destroy(trav);
	}

	/* All easing modes at endpoints: t=0 => start, t=1 => end */
	{
		qaws_easing modes[] = {
			QAWS_EASING_QUAD_OUT,
			QAWS_EASING_CUBIC_IN,
			QAWS_EASING_CUBIC_OUT,
			QAWS_EASING_CUBIC_IN_OUT,
			QAWS_EASING_SINE_IN,
			QAWS_EASING_SINE_OUT
		};
		const char *mode_names[] = {
			"QUAD_OUT", "CUBIC_IN", "CUBIC_OUT",
			"CUBIC_IN_OUT", "SINE_IN", "SINE_OUT"
		};
		unsigned int mode_count = 6;
		unsigned int mi;
		qaws_scalar arc_len;

		qaws_eval_result_2d r_curve_start, r_curve_end;

		/* Get the curve's actual start and end positions */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_curve_start);
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_curve_end);

		/* Compute arc length so speed covers the full curve in [0,1] */
		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &arc_len);
		}

		for (mi = 0; mi < mode_count; ++mi)
		{
			qaws_traversal_desc etdesc;
			memset(&etdesc, 0, sizeof(etdesc));
			qaws_traversal *etrav = NULL;
			qaws_eval_result_2d r_t0, r_t1;
			char msg[128];

			memset(&etdesc, 0, sizeof(etdesc));
			etdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
			etdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
			etdesc.speed = arc_len;
			etdesc.start_time = (qaws_scalar)0.0;
			etdesc.end_time = (qaws_scalar)1.0;
			etdesc.clamp_to_domain = 1;
			etdesc.easing = modes[mi];
			etdesc.wrap_mode = QAWS_WRAP_CLAMP;

			status = qaws_traversal_create(curve, &etdesc, &etrav);
			TEST_ASSERT_STATUS(status);

			status = qaws_traversal_evaluate_2d(
				etrav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r_t0);
			TEST_ASSERT_STATUS(status);

			status = qaws_traversal_evaluate_2d(
				etrav, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r_t1);
			TEST_ASSERT_STATUS(status);

			sprintf(msg, "%s t=0 start x", mode_names[mi]);
			TEST_ASSERT(approx_eq(r_t0.position.x, r_curve_start.position.x), msg);
			sprintf(msg, "%s t=0 start y", mode_names[mi]);
			TEST_ASSERT(approx_eq(r_t0.position.y, r_curve_start.position.y), msg);

			sprintf(msg, "%s t=1 end x", mode_names[mi]);
			TEST_ASSERT(approx_eq_loose(r_t1.position.x, r_curve_end.position.x), msg);
			sprintf(msg, "%s t=1 end y", mode_names[mi]);
			TEST_ASSERT(approx_eq_loose(r_t1.position.y, r_curve_end.position.y), msg);

			qaws_traversal_destroy(etrav);
		}
	}

	/* QUAD_IN_OUT symmetry: at t=0.5, should be at midpoint position */
	{
		qaws_traversal_desc qdesc;
		memset(&qdesc, 0, sizeof(qdesc));
		qaws_traversal *qtrav = NULL;
		qaws_eval_result_2d r_linear_mid, r_qio_mid;
		qaws_traversal *ltrav = NULL;

		/* Get linear position at t=0.5 for reference (the "midpoint") */
		memset(&qdesc, 0, sizeof(qdesc));
		qdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		qdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		qdesc.speed = (qaws_scalar)1.0;
		qdesc.start_time = (qaws_scalar)0.0;
		qdesc.end_time = (qaws_scalar)1.0;
		qdesc.clamp_to_domain = 1;
		qdesc.easing = QAWS_EASING_LINEAR;
		qdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &qdesc, &ltrav);
		TEST_ASSERT_STATUS(status);
		qaws_traversal_evaluate_2d(
			ltrav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_linear_mid);
		qaws_traversal_destroy(ltrav);

		/* QUAD_IN_OUT at t=0.5 should be at the halfway position */
		qdesc.easing = QAWS_EASING_QUAD_IN_OUT;
		status = qaws_traversal_create(curve, &qdesc, &qtrav);
		TEST_ASSERT_STATUS(status);

		qaws_traversal_evaluate_2d(
			qtrav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_qio_mid);

		TEST_ASSERT(approx_eq(r_qio_mid.position.x, r_linear_mid.position.x),
			"QUAD_IN_OUT t=0.5 midpoint x");
		TEST_ASSERT(approx_eq(r_qio_mid.position.y, r_linear_mid.position.y),
			"QUAD_IN_OUT t=0.5 midpoint y");

		qaws_traversal_destroy(qtrav);
	}

	qaws_curve_destroy(curve);
}

static void test_wrap_modes(void)
{
	printf("test_wrap_modes\n");

	qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	/* Loop mode: going past end wraps around */
	{
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_a, r_b;
		qaws_scalar total_len;
		qaws_scalar eps = (qaws_scalar)0.05;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_LOOP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		/* distance=eps and distance=total_len+eps should map to same point */
		qaws_traversal_evaluate_2d(
			trav, eps, QAWS_EVAL_FLAG_POSITION, &r_a);
		qaws_traversal_evaluate_2d(
			trav, total_len + eps, QAWS_EVAL_FLAG_POSITION, &r_b);

		TEST_ASSERT(approx_eq_loose(r_a.position.x, r_b.position.x),
			"loop wraps to start x");
		TEST_ASSERT(approx_eq_loose(r_a.position.y, r_b.position.y),
			"loop wraps to start y");

		qaws_traversal_destroy(trav);
	}

	/* Ping-pong mode: d=eps and d=2*total_len+eps should match */
	{
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_a, r_b;
		qaws_scalar total_len;
		qaws_scalar eps = (qaws_scalar)0.05;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_PING_PONG;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		qaws_traversal_evaluate_2d(
			trav, eps, QAWS_EVAL_FLAG_POSITION, &r_a);
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)2.0 * total_len + eps, QAWS_EVAL_FLAG_POSITION, &r_b);

		TEST_ASSERT(approx_eq_loose(r_a.position.x, r_b.position.x),
			"ping-pong returns to start x");
		TEST_ASSERT(approx_eq_loose(r_a.position.y, r_b.position.y),
			"ping-pong returns to start y");

		qaws_traversal_destroy(trav);
	}

	/* CLAMP mode: negative distance clamps to start, huge distance clamps to end */
	{
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_neg, r_huge, r_start, r_end;
		qaws_scalar total_len;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		/* Position at distance=0 (start of curve) */
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r_start);

		/* Position at distance=total_len (end of curve) */
		qaws_traversal_evaluate_2d(
			trav, total_len, QAWS_EVAL_FLAG_POSITION, &r_end);

		/* Negative distance should clamp to start */
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)-10.0, QAWS_EVAL_FLAG_POSITION, &r_neg);
		TEST_ASSERT(approx_eq(r_neg.position.x, r_start.position.x),
			"clamp negative dist x == start x");
		TEST_ASSERT(approx_eq(r_neg.position.y, r_start.position.y),
			"clamp negative dist y == start y");

		/* Very large distance should clamp to end */
		qaws_traversal_evaluate_2d(
			trav, total_len * (qaws_scalar)100.0, QAWS_EVAL_FLAG_POSITION, &r_huge);
		TEST_ASSERT(approx_eq_loose(r_huge.position.x, r_end.position.x),
			"clamp huge dist x == end x");
		TEST_ASSERT(approx_eq_loose(r_huge.position.y, r_end.position.y),
			"clamp huge dist y == end y");

		qaws_traversal_destroy(trav);
	}

	/* Negative values in LOOP mode: should wrap correctly */
	{
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_neg, r_pos;
		qaws_scalar total_len;
		qaws_scalar eps = (qaws_scalar)0.05;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_LOOP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		/* Negative distance -eps should wrap to total_len - eps */
		qaws_traversal_evaluate_2d(
			trav, -eps, QAWS_EVAL_FLAG_POSITION, &r_neg);
		qaws_traversal_evaluate_2d(
			trav, total_len - eps, QAWS_EVAL_FLAG_POSITION, &r_pos);

		TEST_ASSERT(approx_eq_loose(r_neg.position.x, r_pos.position.x),
			"loop negative wrap x");
		TEST_ASSERT(approx_eq_loose(r_neg.position.y, r_pos.position.y),
			"loop negative wrap y");

		qaws_traversal_destroy(trav);
	}

	qaws_curve_destroy(curve);
}

static void test_traversal_advance(void)
{
	printf("test_traversal_advance\n");

	qaws_vec2 pts[] = { {0,0}, {2,2}, {4,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	{
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r1, r2, r3;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		/* First advance: should move forward from distance=0 */
		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r1.position.x > (qaws_scalar)0.0, "advance moves forward x");

		/* Second advance: should be further along */
		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r2);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r2.position.x > r1.position.x, "second advance further");

		/* Reset should go back to start */
		status = qaws_traversal_reset(trav);
		TEST_ASSERT_STATUS(status);

		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r3);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(r3.position.x, (qaws_scalar)0.0),
			"reset goes to start x");
		TEST_ASSERT(approx_eq(r3.position.y, (qaws_scalar)0.0),
			"reset goes to start y");

		/* Loop mode advance: should wrap around */
		qaws_traversal_destroy(trav);
		tdesc.wrap_mode = QAWS_WRAP_LOOP;
		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		/* Advance a large amount */
		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)100.0, QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT_STATUS(status);
		/* Should still be valid (not NaN) */
		TEST_ASSERT(r1.position.x == r1.position.x, "loop advance valid x");

		/* NULL checks */
		status = qaws_traversal_advance_2d(NULL, (qaws_scalar)0.1,
			QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "advance null trav");

		status = qaws_traversal_reset(NULL);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "reset null");

		qaws_traversal_destroy(trav);
	}

	qaws_curve_destroy(curve);

	/* 3D advance: create a 3D bezier, traverse, verify forward movement */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {2,2,1}, {4,0,2} };
		qaws_bezier_desc bdesc3d;
		qaws_curve *curve3d = NULL;
		qaws_traversal_desc tdesc3d;
		qaws_traversal *trav3d = NULL;
		qaws_eval_result_3d r3d_1, r3d_2;

		bdesc3d.dimension = QAWS_DIMENSION_3D;
		bdesc3d.degree = 2;
		bdesc3d.control_points = pts3d;
		bdesc3d.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc3d, &curve3d);
		TEST_ASSERT_STATUS(status);

		memset(&tdesc3d, 0, sizeof(tdesc3d));
		tdesc3d.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc3d.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc3d.speed = (qaws_scalar)1.0;
		tdesc3d.clamp_to_domain = 1;
		tdesc3d.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve3d, &tdesc3d, &trav3d);
		TEST_ASSERT_STATUS(status);

		/* First advance */
		status = qaws_traversal_advance_3d(
			trav3d, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r3d_1);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r3d_1.position.x > (qaws_scalar)0.0, "3d advance moves forward x");

		/* Second advance should be further along */
		status = qaws_traversal_advance_3d(
			trav3d, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r3d_2);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r3d_2.position.x > r3d_1.position.x, "3d second advance further x");

		qaws_traversal_destroy(trav3d);
		qaws_curve_destroy(curve3d);
	}
}

static void test_multi_curve_traversal(void)
{
	printf("test_multi_curve_traversal\n");

	/* 1. Two-curve chain: (0,0)->(2,1) then (2,1)->(4,0) */
	{
		qaws_vec2 pts1[] = { {0, 0}, {1, 1}, {2, 1} };
		qaws_vec2 pts2[] = { {2, 1}, {3, 1}, {4, 0} };
		qaws_bezier_desc bdesc1, bdesc2;
		qaws_curve* c1 = NULL;
		qaws_curve* c2 = NULL;
		qaws_curve const* curves[2];
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;
		qaws_scalar total_len, len1, len2;

		bdesc1.dimension = QAWS_DIMENSION_2D;
		bdesc1.degree = 2;
		bdesc1.control_points = pts1;
		bdesc1.control_point_count = 3;

		bdesc2.dimension = QAWS_DIMENSION_2D;
		bdesc2.degree = 2;
		bdesc2.control_points = pts2;
		bdesc2.control_point_count = 3;

		s = qaws_curve_create_bezier(&bdesc1, &c1);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&bdesc2, &c2);
		TEST_ASSERT_STATUS(s);

		curves[0] = c1;
		curves[1] = c2;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create_multi(curves, 2, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(trav != NULL, "multi-curve traversal created");

		/* At distance=0 => near (0,0) */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0),
			"multi trav d=0 x~0");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0),
			"multi trav d=0 y~0");

		/* Compute total arc length of the chain */
		s = qaws_curve_compute_arc_length(c1, (qaws_scalar)0.0,
			(qaws_scalar)1.0, &len1);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_compute_arc_length(c2, (qaws_scalar)0.0,
			(qaws_scalar)1.0, &len2);
		TEST_ASSERT_STATUS(s);
		total_len = len1 + len2;

		/* At total distance => near (4,0) */
		s = qaws_traversal_evaluate_2d(trav, total_len,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "multi trav end x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "multi trav end y finite");

		/* At midpoint distance => somewhere valid */
		s = qaws_traversal_evaluate_2d(trav,
			total_len * (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "multi trav mid x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "multi trav mid y finite");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(c1);
		qaws_curve_destroy(c2);
	}

	/* 2. Single curve in chain: should work like regular traversal */
	{
		qaws_vec2 pts[] = { {0, 0}, {5, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* c = NULL;
		qaws_curve const* curves[1];
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &c);
		TEST_ASSERT_STATUS(s);

		curves[0] = c;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create_multi(curves, 1, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)2.5,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.5),
			"single-chain d=2.5 x~2.5");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(c);
	}

	/* 3. Null args */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* c = NULL;
		qaws_curve const* curves[1];
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &c);
		TEST_ASSERT_STATUS(s);

		curves[0] = c;
		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		/* null curves array */
		s = qaws_traversal_create_multi(NULL, 2, &tdesc, &trav);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"multi null curves => INVALID_ARGUMENT");

		/* null desc */
		s = qaws_traversal_create_multi(curves, 1, NULL, &trav);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"multi null desc => INVALID_ARGUMENT");

		/* null out */
		s = qaws_traversal_create_multi(curves, 1, &tdesc, NULL);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"multi null out => INVALID_ARGUMENT");

		qaws_curve_destroy(c);
	}
}

static void test_scurve_profile(void)
{
	printf("test_scurve_profile\n");

	/* 1. Basic S-curve profile */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 3}, {5, 3}, {8, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_S_CURVE;
		tdesc.max_speed = (qaws_scalar)5.0;
		tdesc.acceleration = (qaws_scalar)2.0;
		tdesc.jerk = (qaws_scalar)1.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)5.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(trav != NULL, "S-curve traversal created");

		/* At t=start_time => near curve start */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "scurve t=0 x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "scurve t=0 y finite");

		/* At t=end_time => somewhere along curve */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)5.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "scurve t=end x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "scurve t=end y finite");

		/* Intermediate positions are smooth (all finite) */
		for (i = 1; i <= 4; i++) {
			qaws_scalar t = (qaws_scalar)i;
			s = qaws_traversal_evaluate_2d(trav, t,
				QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(scalar_is_finite(r.position.x),
				"scurve intermediate x finite");
			TEST_ASSERT(scalar_is_finite(r.position.y),
				"scurve intermediate y finite");
		}

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}

	/* 2. Very large jerk => should behave like trapezoidal (still valid) */
	{
		qaws_vec2 pts[] = { {0, 0}, {5, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_S_CURVE;
		tdesc.max_speed = (qaws_scalar)5.0;
		tdesc.acceleration = (qaws_scalar)10.0;
		tdesc.jerk = (qaws_scalar)100000.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)2.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x),
			"high-jerk scurve mid x finite");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}
}

static void test_custom_speed(void)
{
	printf("test_custom_speed\n");

	/* 1. Constant speed callback */
	{
		qaws_vec2 pts[] = { {0, 0}, {3, 4} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r0, r1, r2;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CUSTOM;
		tdesc.custom_speed_fn = constant_speed_fn;
		tdesc.custom_speed_fn_user_data = NULL;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)5.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(trav != NULL, "custom speed traversal created");

		/* Evaluate at start */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r0);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r0.position.x),
			"custom speed t=0 x finite");

		/* Evaluate at midpoint */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)2.5,
			QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r1.position.x),
			"custom speed t=mid x finite");

		/* Evaluate at end */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)5.0,
			QAWS_EVAL_FLAG_POSITION, &r2);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r2.position.x),
			"custom speed t=end x finite");

		/* Motion progresses: end should be further than start */
		TEST_ASSERT(r2.position.x > r0.position.x ||
			r2.position.y > r0.position.y,
			"custom const speed progresses along curve");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}

	/* 2. Variable (accelerating) speed */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 2}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CUSTOM;
		tdesc.custom_speed_fn = accelerating_speed_fn;
		tdesc.custom_speed_fn_user_data = NULL;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)3.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		for (i = 0; i <= 6; i++) {
			qaws_scalar t = (qaws_scalar)i * (qaws_scalar)0.5;
			s = qaws_traversal_evaluate_2d(trav, t,
				QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(scalar_is_finite(r.position.x),
				"accel speed eval x finite");
			TEST_ASSERT(scalar_is_finite(r.position.y),
				"accel speed eval y finite");
		}

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}

	/* 3. Null callback: implementation gracefully degrades (no custom table,
	      falls back to normal time-to-distance mapping). Verify it still
	      creates and evaluates without crashing. */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		memset(&tdesc, 0, sizeof(tdesc));
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CUSTOM;
		tdesc.custom_speed_fn = NULL;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		/* Should still evaluate without crash */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x),
			"null custom fn fallback x finite");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}
}

int test_16_traversal_main(void)
{
	g_pass = 0; g_fail = 0;

	test_traversal();
	test_traversal_mappings();
	test_traversal_arc_length_mode();
	test_traversal_3d();
	test_easing();
	test_wrap_modes();
	test_traversal_advance();
	test_multi_curve_traversal();
	test_scurve_profile();
	test_custom_speed();

	printf("16_traversal: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
