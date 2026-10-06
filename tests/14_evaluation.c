#include "test_common.h"

static void test_arc_length(void)
{
	printf("test_arc_length\n");

	/* Straight line Bezier: arc length should equal distance */
	qaws_vec2 points[] = { {0, 0}, {3, 4} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_scalar length = 0;
	s = qaws_curve_compute_arc_length(curve, (qaws_scalar)0.0, (qaws_scalar)1.0, &length);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(length, (qaws_scalar)5.0), "straight line length == 5");

	qaws_curve_destroy(curve);
}

static void test_inspection_flags(void)
{
	printf("test_inspection_flags\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	TEST_ASSERT(qaws_curve_is_closed(curve) == 0, "bezier is_closed == 0");
	TEST_ASSERT(qaws_curve_is_periodic(curve) == 0, "bezier is_periodic == 0");
	TEST_ASSERT(qaws_curve_is_rational(curve) == 0, "bezier is_rational == 0");
	TEST_ASSERT(qaws_curve_get_continuity(curve) >= QAWS_CONTINUITY_C0, "bezier continuity >= C0");

	qaws_curve_destroy(curve);
}

static void test_bounds_3d(void)
{
	printf("test_bounds_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {1, 2, 3}, {4, 1, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 2;
	desc.control_points = points;
	desc.control_point_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec3 bmin, bmax;
	s = qaws_curve_compute_bounds_3d(curve, &bmin, &bmax);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(bmin.x >= (qaws_scalar)-0.1, "3d bounds min x");
	TEST_ASSERT(bmax.x <= (qaws_scalar)4.1, "3d bounds max x");
	TEST_ASSERT(bmin.z >= (qaws_scalar)-0.1, "3d bounds min z");
	TEST_ASSERT(bmax.z <= (qaws_scalar)3.1, "3d bounds max z");

	qaws_curve_destroy(curve);
}

static void test_precomputed_coefficients(void)
{
	printf("test_precomputed_coefficients\n");

	/* Hermite: verify D1 derivatives still correct with precomputed coeffs */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {0,2}, {0,-2} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* At t=0, D1 should be the tangent (0, 2) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "hermite coeffs pos x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "hermite coeffs pos y");
		TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)0.0), "hermite coeffs d1 x");
		TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)2.0), "hermite coeffs d1 y");

		/* At t=1, position should be (2,0), D1 should be (0,-2) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.0), "hermite coeffs end pos x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "hermite coeffs end pos y");
		TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)0.0), "hermite coeffs end d1 x");
		TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)-2.0), "hermite coeffs end d1 y");

		qaws_curve_destroy(curve);
	}

	/* Trajectory: verify positions and derivatives */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_scalar times[] = {
			(qaws_scalar)0.0, (qaws_scalar)1.0, (qaws_scalar)2.0 };
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = pts;
		tdesc.key_count = 3;
		tdesc.key_times = times;
		tdesc.key_time_count = 3;

		status = qaws_curve_create_trajectory(&tdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Start should be (0,0) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "traj coeffs start x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "traj coeffs start y");

		/* At t=1, should be (1,1) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "traj coeffs mid x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)1.0), "traj coeffs mid y");

		/* At t=2, should be (2,0) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)2.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.0), "traj coeffs end x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "traj coeffs end y");

		qaws_curve_destroy(curve);
	}

	/* Hermite D2/D3 check */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {0,2}, {0,-2} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* At t=0, evaluate all derivatives */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 |
			QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);

		/* D2 should be finite (not NaN) */
		TEST_ASSERT(r.d2.x == r.d2.x, "hermite d2 x finite at t=0");
		TEST_ASSERT(r.d2.y == r.d2.y, "hermite d2 y finite at t=0");

		/* D3 should be finite */
		TEST_ASSERT(r.d3.x == r.d3.x, "hermite d3 x finite at t=0");
		TEST_ASSERT(r.d3.y == r.d3.y, "hermite d3 y finite at t=0");

		/* At t=0.5, D2 and D3 should also be finite */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);
		TEST_ASSERT(r.d2.x == r.d2.x, "hermite d2 x finite at t=0.5");
		TEST_ASSERT(r.d2.y == r.d2.y, "hermite d2 y finite at t=0.5");
		TEST_ASSERT(r.d3.x == r.d3.x, "hermite d3 x finite at t=0.5");
		TEST_ASSERT(r.d3.y == r.d3.y, "hermite d3 y finite at t=0.5");

		/*
		 * For a cubic Hermite, D3 is constant.
		 * Verify D3 at t=0 and t=1 are the same.
		 */
		{
			qaws_eval_result_2d r0, r1;
			qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
				QAWS_EVAL_FLAG_D3, &r0);
			qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
				QAWS_EVAL_FLAG_D3, &r1);
			TEST_ASSERT(approx_eq(r0.d3.x, r1.d3.x),
				"hermite d3 constant x");
			TEST_ASSERT(approx_eq(r0.d3.y, r1.d3.y),
				"hermite d3 constant y");
		}

		qaws_curve_destroy(curve);
	}

	/* Trajectory D1 check: velocity at start matches expected direction */
	{
		qaws_vec2 pts[] = { {0,0}, {3,4}, {6,0} };
		qaws_scalar times[] = {
			(qaws_scalar)0.0, (qaws_scalar)1.0, (qaws_scalar)2.0
		};
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = pts;
		tdesc.key_count = 3;
		tdesc.key_times = times;
		tdesc.key_time_count = 3;

		status = qaws_curve_create_trajectory(&tdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* D1 at t=0 should point roughly toward (3,4) => positive x, positive y */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(r.d1.x > (qaws_scalar)0.0, "traj d1 start x positive");
		TEST_ASSERT(r.d1.y > (qaws_scalar)0.0, "traj d1 start y positive");

		/* D1 should be finite everywhere */
		TEST_ASSERT(r.d1.x == r.d1.x, "traj d1 start x finite");
		TEST_ASSERT(r.d1.y == r.d1.y, "traj d1 start y finite");

		qaws_curve_destroy(curve);
	}

	/* Multi-span Hermite: verify positions at all span boundaries */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1} };
		qaws_vec2 ders[] = { {1,1}, {1,0}, {1,-1} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_span_count(curve) == 2,
			"multi-span hermite span count");

		/* Boundary t=0 -> point 0 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0),
			"multi hermite boundary 0 x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0),
			"multi hermite boundary 0 y");

		/* Boundary t=1 -> point 1 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0),
			"multi hermite boundary 1 x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)2.0),
			"multi hermite boundary 1 y");

		/* Boundary t=2 -> point 2 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)2.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)3.0),
			"multi hermite boundary 2 x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)1.0),
			"multi hermite boundary 2 y");

		qaws_curve_destroy(curve);
	}
}

/* D1, D2 and D3 of every curve family against central differences of the
   order below, at parameters away from span boundaries: the derivatives are
   with respect to the curve parameter (families evaluate spans on a
   normalized local parameter and must scale). */
static int derivative_consistent(qaws_curve const* c, char const* name, int max_order)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	static double const fr[5] = { 0.137, 0.291, 0.457, 0.613, 0.871 };
	double h = (r.max_value - r.min_value) * (QAWS_SCALAR_IS_FLOAT ? 1e-3 : 1e-6);
	double tol = QAWS_SCALAR_IS_FLOAT ? 3e-2 : 1e-5;
	int k, ok = 1;
	unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;
	for (k = 0; k < 5; k++)
	{
		qaws_scalar t = (qaws_scalar)(r.min_value + fr[k] * (r.max_value - r.min_value));
		qaws_eval_result_2d e, ep, em;
		double scale;
		qaws_curve_evaluate_2d(c, t, flags, &e);
		qaws_curve_evaluate_2d(c, (qaws_scalar)(t + h), flags, &ep);
		qaws_curve_evaluate_2d(c, (qaws_scalar)(t - h), flags, &em);
		scale = 1 + sqrt((double)e.d1.x * e.d1.x + (double)e.d1.y * e.d1.y);
		if (fabs(e.d1.x - (ep.position.x - em.position.x) / (2 * h)) > tol * scale ||
		    fabs(e.d1.y - (ep.position.y - em.position.y) / (2 * h)) > tol * scale)
		{
			printf("  %s: D1 at t=%g (%g %g) vs differences (%g %g)\n", name, (double)t, (double)e.d1.x, (double)e.d1.y,
				(ep.position.x - em.position.x) / (2 * h), (ep.position.y - em.position.y) / (2 * h));
			ok = 0;
		}
		scale = 1 + sqrt((double)e.d2.x * e.d2.x + (double)e.d2.y * e.d2.y);
		if (max_order >= 2 && (e.valid_flags & QAWS_EVAL_FLAG_D2) &&
		    (fabs(e.d2.x - (ep.d1.x - em.d1.x) / (2 * h)) > tol * scale || fabs(e.d2.y - (ep.d1.y - em.d1.y) / (2 * h)) > tol * scale))
		{
			printf("  %s: D2 at t=%g (%g %g) vs differences (%g %g)\n", name, (double)t, (double)e.d2.x, (double)e.d2.y,
				(ep.d1.x - em.d1.x) / (2 * h), (ep.d1.y - em.d1.y) / (2 * h));
			ok = 0;
		}
		scale = 1 + sqrt((double)e.d3.x * e.d3.x + (double)e.d3.y * e.d3.y);
		if (max_order >= 3 && (e.valid_flags & QAWS_EVAL_FLAG_D3) &&
		    (fabs(e.d3.x - (ep.d2.x - em.d2.x) / (2 * h)) > tol * scale || fabs(e.d3.y - (ep.d2.y - em.d2.y) / (2 * h)) > tol * scale))
		{
			printf("  %s: D3 at t=%g (%g %g) vs differences (%g %g)\n", name, (double)t, (double)e.d3.x, (double)e.d3.y,
				(ep.d2.x - em.d2.x) / (2 * h), (ep.d2.y - em.d2.y) / (2 * h));
			ok = 0;
		}
	}
	return ok;
}

static void test_family_derivatives(void)
{
	static qaws_scalar const p4[8] = { 0, 0, 1, 2, 3, 2.5f, 4, 0.5f };
	static qaws_scalar const p6[12] = { 0, 0, 1, 1.5f, 2.2f, -0.4f, 3.1f, 1.2f, 4.4f, 0.3f, 5.2f, 1.6f };
	static qaws_scalar const knots[10] = { 0, 0, 0, 0, 0.3f, 1.1f, 2, 2, 2, 2 };
	static qaws_scalar const w6[6] = { 1, 1.4f, 0.7f, 1.2f, 0.9f, 1 };
	char msg[128];
	qaws_curve* c;

#define CHECK_FAMILY(curve_expr, name, order) \
	do { c = NULL; curve_expr; if (c) { sprintf(msg, "%s: derivatives match differences", name); \
		TEST_ASSERT(derivative_consistent(c, name, order), msg); qaws_curve_destroy(c); } \
		else { sprintf(msg, "%s: created", name); TEST_ASSERT(0, msg); } } while (0)

	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p4; d.control_point_count = 4;
		CHECK_FAMILY(qaws_curve_create_bezier(&d, &c), "bezier", 3);
	}
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p6; d.control_point_count = 6;
		d.knots = knots; d.knot_count = 10;
		CHECK_FAMILY(qaws_curve_create_bspline(&d, &c), "bspline (non-uniform knots)", 3);
	}
	{
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p6; d.control_point_count = 6;
		d.knots = knots; d.knot_count = 10; d.weights = w6; d.weight_count = 6;
		CHECK_FAMILY(qaws_curve_create_nurbs(&d, &c), "nurbs", 3);
	}
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p4; d.control_point_count = 4;
		d.weights = w6; d.weight_count = 4;
		CHECK_FAMILY(qaws_curve_create_rational_bezier(&d, &c), "rational bezier", 3);
	}
	{
		static qaws_scalar const der[8] = { 1, 0, 0.5f, 1, 1, -1, 0, 1 };
		qaws_hermite_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.points = p4; d.derivatives = der; d.point_count = 4; d.derivative_count = 4;
		CHECK_FAMILY(qaws_curve_create_hermite(&d, &c), "hermite", 3);
	}
	{
		qaws_catmull_rom_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.control_points = p6; d.control_point_count = 6;
		d.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		CHECK_FAMILY(qaws_curve_create_catmull_rom(&d, &c), "catmull-rom (centripetal)", 3);
	}
	{
		qaws_yuksel_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.control_points = p6; d.control_point_count = 6; d.mode = QAWS_YUKSEL_MODE_BEZIER;
		CHECK_FAMILY(qaws_curve_create_yuksel(&d, &c), "yuksel (bezier)", 2);
	}
	{
		qaws_yuksel_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.control_points = p6; d.control_point_count = 6; d.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		CHECK_FAMILY(qaws_curve_create_yuksel(&d, &c), "yuksel (circular)", 2);
	}
	{
		static qaws_scalar const times[4] = { 0, 0.7f, 1.9f, 3.0f };
		qaws_trajectory_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.key_positions = p4; d.key_count = 4; d.key_times = times; d.key_time_count = 4;
		CHECK_FAMILY(qaws_curve_create_trajectory(&d, &c), "trajectory (non-uniform times)", 3);
	}
	{
		static qaws_scalar const coef[8] = { 0.2f, -0.1f, 1.1f, 0.4f, -0.3f, 0.9f, 0.05f, -0.2f };
		qaws_polynomial_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.coefficients = coef; d.coefficient_count = 4;
		d.t_min = (qaws_scalar)0.5; d.t_max = (qaws_scalar)2.5;
		CHECK_FAMILY(qaws_curve_create_polynomial(&d, &c), "polynomial on [0.5, 2.5]", 3);
	}
	{
		qaws_clothoid_desc d;
		memset(&d, 0, sizeof(d));
		d.start_angle = (qaws_scalar)0.3; d.start_curvature = (qaws_scalar)0.2; d.end_curvature = (qaws_scalar)1.5;
		d.length = (qaws_scalar)3.0;
		CHECK_FAMILY(qaws_curve_create_clothoid(&d, &c), "clothoid", 3);
	}
	{
		qaws_arc_segment seg[2];
		qaws_arc_desc d;
		memset(seg, 0, sizeof(seg));
		seg[0].radius = (qaws_scalar)1.5; seg[0].angle_start = (qaws_scalar)0.2; seg[0].angle_end = (qaws_scalar)1.4;
		seg[1].center[0] = (qaws_scalar)2; seg[1].radius = (qaws_scalar)0.6; seg[1].angle_start = (qaws_scalar)3.0;
		seg[1].angle_end = (qaws_scalar)1.0;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = seg; d.segment_count = 2;
		CHECK_FAMILY(qaws_curve_create_arc(&d, &c), "arc (two segments)", 3);
	}
	{
		qaws_subdivision_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_3; d.control_points = p6; d.control_point_count = 6;
		d.refinement_levels = 6;
		CHECK_FAMILY(qaws_curve_create_subdivision(&d, &c), "subdivision (cubic)", 1);
	}
	{
		qaws_bezier_desc b;
		qaws_arc_segment seg;
		qaws_arc_desc a;
		qaws_curve* parts[2] = { NULL, NULL };
		qaws_composite_desc d;
		memset(&b, 0, sizeof(b));
		b.dimension = QAWS_DIMENSION_2D; b.degree = 3; b.control_points = p4; b.control_point_count = 4;
		qaws_curve_create_bezier(&b, &parts[0]);
		memset(&seg, 0, sizeof(seg));
		seg.center[0] = (qaws_scalar)4; seg.center[1] = (qaws_scalar)1.5; seg.radius = (qaws_scalar)1.0;
		seg.angle_start = (qaws_scalar)-1.6; seg.angle_end = (qaws_scalar)0.4;
		memset(&a, 0, sizeof(a));
		a.dimension = QAWS_DIMENSION_2D; a.segments = &seg; a.segment_count = 1;
		qaws_curve_create_arc(&a, &parts[1]);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = parts; d.segment_count = 2;
		CHECK_FAMILY(qaws_curve_create_composite(&d, &c), "composite (bezier + arc)", 3);
	}
#undef CHECK_FAMILY
}

int test_14_evaluation_main(void)
{
	g_pass = 0; g_fail = 0;

	test_arc_length();
	test_inspection_flags();
	test_bounds_3d();
	test_precomputed_coefficients();
	test_family_derivatives();

	printf("14_evaluation: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
