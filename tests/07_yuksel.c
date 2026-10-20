#include "test_common.h"

static int scalar_is_finite(qaws_scalar v)
{
	return (v == v) && (v - v == 0);
}

static void test_yuksel_bezier_mode(void)
{
	printf("test_yuksel_bezier_mode\n");

	/* Triangle-ish open curve: 5 points */
	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, 0}, {3, 2}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(curve != NULL, "yuksel bezier curve allocated");
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_YUKSEL, "kind == YUKSEL");
	TEST_ASSERT(qaws_curve_get_continuity(curve) == QAWS_CONTINUITY_C2, "continuity == C2");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 4, "span_count == 4");
	TEST_ASSERT(qaws_curve_is_closed(curve) == 0, "not closed");

	/* Evaluate at parameter 0 (first span start) */
	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(&r, curve,
		(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_POSITION, "position valid");
	TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D1, "d1 valid");

	/* Evaluate at each integer parameter to check interpolation */
	{
		unsigned int pi;
		qaws_range range = qaws_curve_get_parameter_range(curve);
		/* Sample throughout to check no NaN or crash */
		for (pi = 0; pi < 40; pi++) {
			qaws_scalar t = range.min_value +
				(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)39;
			s = qaws_curve_evaluate_2d(&r, curve, t, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(s == QAWS_STATUS_OK, "eval ok across range");
		}
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_closed(void)
{
	printf("test_yuksel_closed\n");

	/* Pentagon */
	qaws_vec2 points[5];
	{
		unsigned int pi;
		for (pi = 0; pi < 5; pi++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979323846) * (qaws_scalar)pi / (qaws_scalar)5;
			points[pi].x = (qaws_scalar)cos((double)angle);
			points[pi].y = (qaws_scalar)sin((double)angle);
		}
	}

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 1;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_is_closed(curve) == 1, "closed");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 5, "5 spans for closed pentagon");

	/* Closed curve: first and last points should match */
	qaws_eval_result_2d r0, r1;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	s = qaws_curve_evaluate_2d(&r0, curve, range.min_value, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	s = qaws_curve_evaluate_2d(&r1, curve, range.max_value, QAWS_EVAL_FLAG_POSITION);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r0.position.x, r1.position.x), "closed: start.x == end.x");
	TEST_ASSERT(approx_eq(r0.position.y, r1.position.y), "closed: start.y == end.y");

	qaws_curve_destroy(curve);
}

static void test_yuksel_circular_mode(void)
{
	printf("test_yuksel_circular_mode\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 1}, {2, 0}, {3, 1}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	/* Sample throughout */
	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_2d(&r, curve, t, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(s == QAWS_STATUS_OK, "circular eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_elliptical_mode(void)
{
	printf("test_yuksel_elliptical_mode\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, 0}, {3, 2}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_ELLIPTICAL;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_2d(&r, curve, t, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(s == QAWS_STATUS_OK, "elliptical eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_hybrid_mode(void)
{
	printf("test_yuksel_hybrid_mode\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, 0}, {3, 2}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_HYBRID;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_2d(&r, curve, t, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(s == QAWS_STATUS_OK, "hybrid eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_3d(void)
{
	printf("test_yuksel_3d\n");

	qaws_vec3 points[] = {
		{0, 0, 0}, {1, 1, 1}, {2, 0, 2}, {3, 1, 1}, {4, 0, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_dimension(curve) == QAWS_DIMENSION_3D, "dim == 3D");

	qaws_eval_result_3d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_3d(&r, curve, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		TEST_ASSERT(s == QAWS_STATUS_OK, "3d eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_error_handling(void)
{
	printf("test_yuksel_error_handling\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 2; /* too few */
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&curve, &desc);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT, "too few points");
	TEST_ASSERT(curve == NULL, "curve is NULL on error");

	/* 3D with circular mode should now succeed */
	{
		qaws_vec3 pts3[] = { {0,0,0}, {1,1,1}, {2,0,0} };
		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts3;
		desc.control_point_count = 3;
		desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		s = qaws_curve_create_yuksel(&curve, &desc);
		TEST_ASSERT(s == QAWS_STATUS_OK, "circular 3d accepted");
		if (curve) {
			/* Verify we can evaluate the curve */
			qaws_eval_result_3d res3d;
			qaws_status es = qaws_curve_evaluate_3d(&res3d, curve,
				(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(es == QAWS_STATUS_OK, "circular 3d eval ok");
			qaws_curve_destroy(curve);
			curve = NULL;
		}
	}

	/* NULL desc */
	s = qaws_curve_create_yuksel(&curve, NULL);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "null desc rejected");
}

static void test_yuksel_3d_circular(void)
{
	printf("test_yuksel_3d_circular\n");

	/* 1. 3D circular Yuksel: 5 points in 3D */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 2, 1}, {2, 0, 2}, {3, 2, 1}, {4, 0, 0}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_eval_result_3d r;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&curve, &desc);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(curve != NULL, "3D circular yuksel created");

		range = qaws_curve_get_parameter_range(curve);

		/* Evaluate at start */
		s = qaws_curve_evaluate_3d(&r, curve,
			range.min_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d circ start x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "3d circ start y finite");
		TEST_ASSERT(scalar_is_finite(r.position.z), "3d circ start z finite");

		/* Evaluate at midpoint */
		{
			qaws_scalar tmid = (range.min_value + range.max_value)
				* (qaws_scalar)0.5;
			s = qaws_curve_evaluate_3d(&r, curve,
				tmid, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(scalar_is_finite(r.position.x), "3d circ mid x finite");
		}

		/* Evaluate at end */
		s = qaws_curve_evaluate_3d(&r, curve,
			range.max_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d circ end x finite");

		qaws_curve_destroy(curve);
	}

	/* 2. 3D elliptical Yuksel */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 2, 1}, {2, 0, 2}, {3, 2, 1}, {4, 0, 0}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_eval_result_3d r;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_ELLIPTICAL;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&curve, &desc);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(curve != NULL, "3D elliptical yuksel created");

		range = qaws_curve_get_parameter_range(curve);

		s = qaws_curve_evaluate_3d(&r, curve,
			range.min_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d ellip start finite");

		s = qaws_curve_evaluate_3d(&r, curve,
			range.max_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d ellip end finite");

		qaws_curve_destroy(curve);
	}

	/* 3. 3D hybrid Yuksel */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 2, 1}, {2, 0, 2}, {3, 2, 1}, {4, 0, 0}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_eval_result_3d r;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_HYBRID;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&curve, &desc);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(curve != NULL, "3D hybrid yuksel created");

		range = qaws_curve_get_parameter_range(curve);

		s = qaws_curve_evaluate_3d(&r, curve,
			range.min_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d hybrid start finite");

		s = qaws_curve_evaluate_3d(&r, curve,
			range.max_value, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d hybrid end finite");

		qaws_curve_destroy(curve);
	}

	/* 4. Collinear 3D points: should fall back gracefully */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&curve, &desc);
		/* Should succeed (fallback to bezier) or return degenerate */
		TEST_ASSERT(s == QAWS_STATUS_OK || s == QAWS_STATUS_DEGENERATE_CURVE,
			"collinear 3D yuksel handled");

		if (s == QAWS_STATUS_OK && curve != NULL) {
			qaws_eval_result_3d r;
			qaws_range range = qaws_curve_get_parameter_range(curve);
			unsigned int pi;

			/* Evaluate across the range - should not crash */
			for (pi = 0; pi < 10; pi++) {
				qaws_scalar t = range.min_value +
					(qaws_scalar)pi * (range.max_value - range.min_value)
					/ (qaws_scalar)9;
				qaws_status es = qaws_curve_evaluate_3d(&r, curve,
					t, QAWS_EVAL_FLAG_POSITION);
				TEST_ASSERT(es == QAWS_STATUS_OK, "collinear 3d eval ok");
			}

			qaws_curve_destroy(curve);
		}
	}
}

int test_07_yuksel_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_yuksel_bezier_mode();
	test_yuksel_closed();
	test_yuksel_circular_mode();
	test_yuksel_elliptical_mode();
	test_yuksel_hybrid_mode();
	test_yuksel_3d();
	test_yuksel_error_handling();
	test_yuksel_3d_circular();

	printf("07_yuksel: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
