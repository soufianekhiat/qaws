#include "test_common.h"

static int scalar_is_finite(qaws_scalar v)
{
	return (v == v) && (v - v == 0);
}

static void test_inspect(void)
{
	printf("test_inspect\n");

	qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 2;
	desc.control_points = points;
	desc.control_point_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	/* Bounds */
	qaws_vec2 bmin, bmax;
	s = qaws_curve_compute_bounds_2d(&bmin, &bmax, curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(bmin.x >= (qaws_scalar)-0.1, "bounds min x reasonable");
	TEST_ASSERT(bmax.x <= (qaws_scalar)3.1, "bounds max x reasonable");

	/* Get control points back */
	qaws_vec2 out_pts[3];
	unsigned int out_count = 0;
	s = qaws_bezier_get_control_points(out_pts, &out_count, curve, 3);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(out_count == 3, "got 3 control points");
	TEST_ASSERT(approx_eq(out_pts[0].x, (qaws_scalar)0.0), "cp0 x");

	qaws_curve_destroy(curve);
}

static void test_closest_point(void)
{
	printf("test_closest_point\n");

	/* Straight line from (0,0) to (10,0): closest to (5,3) should be t=0.5 */
	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 target = { 5, 3 };
	qaws_scalar param = 0;
	s = qaws_curve_find_closest_parameter_2d(&param, curve, target);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(param, (qaws_scalar)0.5), "closest param == 0.5");

	qaws_curve_destroy(curve);
}

static void test_span_continuity(void)
{
	printf("test_span_continuity\n");

	/* Multi-span Hermite: 3 points => 2 spans => 1 interior boundary */
	qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0} };
	qaws_vec2 tans[] = { {1, 0}, {0, 0}, {1, 0} };

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = pts;
	desc.derivatives = tans;
	desc.point_count = 3;
	desc.derivative_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&curve, &desc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 2, "hermite 3pt span_count == 2");

	qaws_continuity cont;
	s = qaws_curve_get_span_continuity(&cont, curve, 0);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(cont >= QAWS_CONTINUITY_C0, "span continuity >= C0");

	qaws_curve_destroy(curve);
}

static void test_family_inspection(void)
{
	printf("test_family_inspection\n");

	/* bspline_get_knots */
	{
		qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
		qaws_bspline_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = points;
		desc.control_point_count = 4;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_uniform = 1;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_bspline(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		qaws_scalar knots[16];
		unsigned int knot_count = 0;
		s = qaws_bspline_get_knots(knots, &knot_count, curve, 16);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(knot_count > 0, "bspline knot_count > 0");
		/* Clamped: first degree+1 knots should be equal */
		TEST_ASSERT(approx_eq(knots[0], knots[1]), "clamped knots start equal");

		/* Wrong family: bezier_get_control_points on bspline */
		qaws_vec2 dummy[4];
		unsigned int dummy_count = 0;
		s = qaws_bezier_get_control_points(dummy, &dummy_count, curve, 4);
		TEST_ASSERT(s == QAWS_STATUS_UNSUPPORTED_OPERATION, "bezier_get_cp on bspline => UNSUPPORTED");

		qaws_curve_destroy(curve);
	}

	/* nurbs_get_weights */
	{
		qaws_vec2 points[] = { {1, 0}, {1, 1}, {0, 1} };
		qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };
		qaws_scalar weights[] = { 1, (qaws_scalar)0.707, 1 };

		qaws_nurbs_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = points;
		desc.control_point_count = 3;
		desc.knots = knots;
		desc.knot_count = 6;
		desc.weights = weights;
		desc.weight_count = 3;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_nurbs(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		qaws_scalar out_w[3];
		unsigned int w_count = 0;
		s = qaws_nurbs_get_weights(out_w, &w_count, curve, 3);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(w_count == 3, "nurbs weight_count == 3");
		TEST_ASSERT(approx_eq(out_w[0], (qaws_scalar)1.0), "weight 0 == 1");

		qaws_curve_destroy(curve);
	}
}

static void test_frenet_frame(void)
{
	printf("test_frenet_frame\n");

	/* 2D normal: perpendicular to tangent */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec2 normal;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		/* Horizontal line: tangent is (1,0), normal should be (0,1) */
		status = qaws_curve_compute_normal_2d(&normal, curve, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(normal.x, (qaws_scalar)0.0), "2D normal x=0 for horizontal");
		TEST_ASSERT(approx_eq(normal.y, (qaws_scalar)1.0), "2D normal y=1 for horizontal");

		qaws_curve_destroy(curve);
	}

	/* 2D normal for vertical line */
	{
		qaws_vec2 pts[] = { {0, 0}, {0, 2} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec2 normal;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		/* Vertical line going up: tangent (0,1), normal (-1,0) */
		status = qaws_curve_compute_normal_2d(&normal, curve, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(normal.x, (qaws_scalar)-1.0), "2D normal x=-1 for vertical");
		TEST_ASSERT(approx_eq(normal.y, (qaws_scalar)0.0), "2D normal y=0 for vertical");

		qaws_curve_destroy(curve);
	}

	/* 3D Frenet frame: planar curve in XY */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,2,0}, {2,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec3 T, N, B;
		qaws_scalar dot_TN, dot_TB, dot_NB;
		qaws_scalar T_len, N_len, B_len;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_compute_frenet_frame_3d(
			&T, &N, &B, curve, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);

		/* T, N, B should be orthonormal */
		T_len = (qaws_scalar)sqrt((double)(T.x*T.x + T.y*T.y + T.z*T.z));
		N_len = (qaws_scalar)sqrt((double)(N.x*N.x + N.y*N.y + N.z*N.z));
		B_len = (qaws_scalar)sqrt((double)(B.x*B.x + B.y*B.y + B.z*B.z));
		TEST_ASSERT(approx_eq(T_len, (qaws_scalar)1.0), "T is unit length");
		TEST_ASSERT(approx_eq(N_len, (qaws_scalar)1.0), "N is unit length");
		TEST_ASSERT(approx_eq(B_len, (qaws_scalar)1.0), "B is unit length");

		dot_TN = T.x*N.x + T.y*N.y + T.z*N.z;
		dot_TB = T.x*B.x + T.y*B.y + T.z*B.z;
		dot_NB = N.x*B.x + N.y*B.y + N.z*B.z;
		TEST_ASSERT(approx_eq(dot_TN, (qaws_scalar)0.0), "T perpendicular to N");
		TEST_ASSERT(approx_eq(dot_TB, (qaws_scalar)0.0), "T perpendicular to B");
		TEST_ASSERT(approx_eq(dot_NB, (qaws_scalar)0.0), "N perpendicular to B");

		/* For planar XY curve, binormal should point along Z */
		TEST_ASSERT(B.z != (qaws_scalar)0.0, "B has Z component for XY curve");

		qaws_curve_destroy(curve);
	}

	/* 3D Frenet frame: straight line (degenerate) */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec3 T, N, B;
		qaws_scalar dot_TN;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		/* Straight line: D1xD2=0, should still produce valid frame */
		status = qaws_curve_compute_frenet_frame_3d(
			&T, &N, &B, curve, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(T.x, (qaws_scalar)1.0), "straight line T along X");

		dot_TN = T.x*N.x + T.y*N.y + T.z*N.z;
		TEST_ASSERT(approx_eq(dot_TN, (qaws_scalar)0.0), "T perp N for straight line");

		qaws_curve_destroy(curve);
	}

	/* Dimension mismatch */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec3 T, N, B;
		qaws_vec2 normal2d;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;
		status = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_compute_frenet_frame_3d(
			&T, &N, &B, curve, (qaws_scalar)0.5);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION, "3D frame rejects 2D curve");

		status = qaws_curve_compute_normal_2d(&normal2d, NULL, (qaws_scalar)0.5);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "null curve");

		qaws_curve_destroy(curve);
	}
}

static void test_geometric_helpers(void)
{
	printf("test_geometric_helpers\n");

	/* Create a known curve: quadratic Bezier in 2D */
	qaws_vec2 pts2d[] = { {0, 0}, {1, 2}, {2, 0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve2d = NULL;
	qaws_status status;
	qaws_vec2 tangent2d;
	qaws_scalar curvature2d, speed_val;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts2d;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&curve2d, &bdesc);
	TEST_ASSERT_STATUS(status);

	/* Tangent at t=0 should point roughly toward (1,2), normalized */
	status = qaws_curve_compute_tangent_2d(&tangent2d, curve2d, (qaws_scalar)0.0);
	TEST_ASSERT_STATUS(status);
	{
		qaws_scalar len = (qaws_scalar)sqrt(
			(double)(tangent2d.x * tangent2d.x + tangent2d.y * tangent2d.y));
		TEST_ASSERT(approx_eq(len, (qaws_scalar)1.0), "tangent is unit length");
	}
	TEST_ASSERT(tangent2d.x > (qaws_scalar)0.0, "tangent2d x > 0 at t=0");
	TEST_ASSERT(tangent2d.y > (qaws_scalar)0.0, "tangent2d y > 0 at t=0");

	/* Tangent at t=0.5 should be horizontal (symmetric parabola) */
	status = qaws_curve_compute_tangent_2d(&tangent2d, curve2d, (qaws_scalar)0.5);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(approx_eq(tangent2d.y, (qaws_scalar)0.0), "tangent horizontal at midpoint");

	/* Curvature should be negative (concave down parabola) */
	status = qaws_curve_compute_curvature_2d(
		&curvature2d, curve2d, (qaws_scalar)0.5);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(curvature2d < (qaws_scalar)0.0, "negative curvature for concave-down");

	/* Speed at t=0 */
	status = qaws_curve_compute_speed(&speed_val, curve2d, (qaws_scalar)0.0);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(speed_val > (qaws_scalar)0.0, "speed > 0");

	/* Dimension mismatch errors */
	{
		qaws_vec3 tangent3d;
		qaws_scalar curvature3d, torsion;

		status = qaws_curve_compute_tangent_3d(
			&tangent3d, curve2d, (qaws_scalar)0.5);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION,
			"tangent_3d rejects 2D curve");

		status = qaws_curve_compute_curvature_3d(
			&curvature3d, curve2d, (qaws_scalar)0.5);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION,
			"curvature_3d rejects 2D curve");

		status = qaws_curve_compute_torsion_3d(
			&torsion, curve2d, (qaws_scalar)0.5);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION,
			"torsion_3d rejects 2D curve");
	}

	qaws_curve_destroy(curve2d);

	/* 3D curve tests */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {1,2,0}, {2,0,0} };
		qaws_curve *curve3d = NULL;
		qaws_vec3 tangent3d;
		qaws_scalar curvature3d, torsion, speed3d;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 2;
		bdesc.control_points = pts3d;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&curve3d, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_compute_tangent_3d(
			&tangent3d, curve3d, (qaws_scalar)0.0);
		TEST_ASSERT_STATUS(status);
		{
			qaws_scalar len = (qaws_scalar)sqrt((double)(
				tangent3d.x * tangent3d.x +
				tangent3d.y * tangent3d.y +
				tangent3d.z * tangent3d.z));
			TEST_ASSERT(approx_eq(len, (qaws_scalar)1.0), "3D tangent is unit length");
		}

		status = qaws_curve_compute_curvature_3d(
			&curvature3d, curve3d, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(curvature3d > (qaws_scalar)0.0, "3D curvature > 0");

		/* Planar curve has zero torsion */
		status = qaws_curve_compute_torsion_3d(
			&torsion, curve3d, (qaws_scalar)0.5);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(torsion, (qaws_scalar)0.0),
			"planar curve has zero torsion");

		status = qaws_curve_compute_speed(
			&speed3d, curve3d, (qaws_scalar)0.0);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(speed3d > (qaws_scalar)0.0, "3D speed > 0");

		qaws_curve_destroy(curve3d);
	}

	/* NULL argument tests */
	status = qaws_curve_compute_tangent_2d(&tangent2d, NULL, (qaws_scalar)0.0);
	TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "tangent null curve");

	status = qaws_curve_compute_speed(&speed_val, NULL, (qaws_scalar)0.0);
	TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "speed null curve");
}

static void test_inflection_points(void)
{
	printf("test_inflection_points\n");

	/* 1. Cubic Bezier S-curve: CPs (0,0),(1,4),(2,-2),(3,0)
	 *    inflection at t = 5/9 ~ 0.556 (avoids exact sample alignment) */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 4}, {2, -2}, {3, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(params, &count, curve, 4);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 1, "S-curve has 1 inflection point");

		if (count >= 1) {
			qaws_scalar diff = params[0] - (qaws_scalar)(5.0 / 9.0);
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.1, "inflection near t=5/9");

			/* Curvature at inflection should be near 0 */
			{
				qaws_scalar curv = 0;
				s = qaws_curve_compute_curvature_2d(&curv, curve, params[0]);
				TEST_ASSERT_STATUS(s);
				if (curv < 0) curv = -curv;
				TEST_ASSERT(curv < (qaws_scalar)0.1, "curvature ~0 at inflection");
			}
		}

		qaws_curve_destroy(curve);
	}

	/* 2. Straight line: collinear points, no inflection */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 2}, {3, 3} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 99;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(params, &count, curve, 4);
		TEST_ASSERT(s == QAWS_STATUS_OK || s == QAWS_STATUS_DEGENERATE_CURVE,
			"collinear inflection handled gracefully");
		if (s == QAWS_STATUS_OK) {
			TEST_ASSERT(count == 0, "collinear has 0 inflections");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Quadratic Bezier: constant curvature sign, 0 inflections */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 99;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(params, &count, curve, 4);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 0, "quadratic has 0 inflections");

		qaws_curve_destroy(curve);
	}

	/* 4. Null out_parameters: requires non-NULL => INVALID_ARGUMENT */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 4}, {2, -2}, {3, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(NULL, &count, curve, 0);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null out_parameters => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}

	/* 5. Null args */
	{
		qaws_status s;
		unsigned int count = 0;

		s = qaws_curve_find_inflection_points(NULL, &count, NULL, 0);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "null curve => INVALID_ARGUMENT");

		/* Null out_count */
		{
			qaws_vec2 pts[] = { {0, 0}, {1, 3}, {2, -3}, {3, 0} };
			qaws_bezier_desc desc;
			qaws_curve* curve = NULL;
			qaws_scalar params[4];

			desc.dimension = QAWS_DIMENSION_2D;
			desc.degree = 3;
			desc.control_points = pts;
			desc.control_point_count = 4;

			s = qaws_curve_create_bezier(&curve, &desc);
			TEST_ASSERT_STATUS(s);

			s = qaws_curve_find_inflection_points(params, NULL, curve, 4);
			TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
				"null out_count => INVALID_ARGUMENT");

			qaws_curve_destroy(curve);
		}
	}
}

static void test_extrema(void)
{
	printf("test_extrema\n");

	/* 1. Cubic Bezier: (0,0),(1,4),(2,1),(3,0) - Y extremum near t~0.377
	 *    (asymmetric to avoid landing on exact sample grid) */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 4}, {2, 1}, {3, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(params, &count, curve, 4, 1);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 1, "cubic has >=1 Y-extremum");

		if (count >= 1) {
			qaws_scalar diff = params[0] - (qaws_scalar)0.377;
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.1,
				"Y-extremum near t~0.377");

			/* D1.y should be ~0 at the extremum */
			{
				qaws_eval_result_2d r;
				s = qaws_curve_evaluate_2d(&r, curve,
					params[0], QAWS_EVAL_FLAG_D1);
				TEST_ASSERT_STATUS(s);
				{
					qaws_scalar d1y = r.d1.y;
					if (d1y < 0) d1y = -d1y;
					TEST_ASSERT(d1y < (qaws_scalar)0.1,
						"D1.y ~0 at Y-extremum");
				}
			}
		}

		qaws_curve_destroy(curve);
	}

	/* 2. X-extrema of S-curve */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 1}, {-1, 2}, {1, 3} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(params, &count, curve, 4, 0);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 1, "S-curve has X-extrema");

		qaws_curve_destroy(curve);
	}

	/* 3. Invalid axis: axis=2 for 2D curve */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(params, &count, curve, 4, 2);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT
			|| s == QAWS_STATUS_OUT_OF_RANGE,
			"axis=2 on 2D curve => error");

		qaws_curve_destroy(curve);
	}

	/* 4. 3D curve: find z-axis extrema (axis=2) using cubic
	 *    (0,0,0),(1,1,4),(2,0,1),(3,0,0) -> z extremum at t~0.377 */
	{
		qaws_vec3 pts[] = { {0, 0, 0}, {1, 1, 4}, {2, 0, 1}, {3, 0, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(params, &count, curve, 4, 2);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 1, "3D curve has z-extremum");

		if (count >= 1) {
			qaws_scalar diff = params[0] - (qaws_scalar)0.377;
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.1,
				"z-extremum near t~0.377");
		}

		qaws_curve_destroy(curve);
	}

	/* 5. Null args */
	{
		qaws_status s;
		unsigned int count = 0;

		s = qaws_curve_find_extrema(NULL, &count, NULL, 0, 0);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve extrema => INVALID_ARGUMENT");
	}
}

static void test_curvature_comb(void)
{
	printf("test_curvature_comb\n");

	/* 1. 2D curvature comb on quadratic parabola */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_2d samples[10];
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_2d(samples, curve, 10, 10);
		TEST_ASSERT_STATUS(s);

		for (i = 0; i < 10; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].position.x),
				"2d comb position.x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].position.y),
				"2d comb position.y finite");
			TEST_ASSERT(scalar_is_finite(samples[i].curvature),
				"2d comb curvature finite");
			{
				qaws_scalar nx = samples[i].normal.x;
				qaws_scalar ny = samples[i].normal.y;
				qaws_scalar nlen = (qaws_scalar)sqrt(
					(double)(nx * nx + ny * ny));
				TEST_ASSERT(approx_eq(nlen, (qaws_scalar)1.0),
					"2d comb normal is unit length");
			}
		}

		qaws_curve_destroy(curve);
	}

	/* 2. 3D curvature comb */
	{
		qaws_vec3 pts[] = { {0, 0, 0}, {1, 2, 1}, {2, 0, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_3d samples[10];
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_3d(samples, curve, 10, 10);
		TEST_ASSERT_STATUS(s);

		for (i = 0; i < 10; i++) {
			qaws_scalar nx = samples[i].normal.x;
			qaws_scalar ny = samples[i].normal.y;
			qaws_scalar nz = samples[i].normal.z;
			qaws_scalar nlen = (qaws_scalar)sqrt(
				(double)(nx * nx + ny * ny + nz * nz));
			TEST_ASSERT(approx_eq(nlen, (qaws_scalar)1.0),
				"3d comb normal is unit length");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Buffer too small: capacity < sample_count */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_2d samples[3];

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_2d(samples, curve, 10, 3);
		TEST_ASSERT(s == QAWS_STATUS_BUFFER_TOO_SMALL,
			"comb buffer too small");

		qaws_curve_destroy(curve);
	}

	/* 4. sample_count < 2 => INVALID_ARGUMENT */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_2d samples[2];

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_2d(samples, curve, 1, 2);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"sample_count < 2 => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}

	/* 5. Null args */
	{
		qaws_status s;
		qaws_curvature_sample_2d samples[10];

		s = qaws_curve_compute_curvature_comb_2d(samples, NULL, 10, 10);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve comb => INVALID_ARGUMENT");
	}
}

static void test_winding_number(void)
{
	printf("test_winding_number\n");

	/* 1. Closed Catmull-Rom circle: inside and outside points */
	{
		qaws_vec2 points[8];
		qaws_catmull_rom_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		int winding = 0;
		unsigned int pi;

		/* 8 points around a circle of radius 2 */
		for (pi = 0; pi < 8; pi++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979323846)
				* (qaws_scalar)pi / (qaws_scalar)8;
			points[pi].x = (qaws_scalar)(2.0 * cos((double)angle));
			points[pi].y = (qaws_scalar)(2.0 * sin((double)angle));
		}

		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 8;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 1;

		s = qaws_curve_create_catmull_rom(&curve, &desc);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(qaws_curve_is_closed(curve) == 1, "circle is closed");

		/* Inside point (origin): winding number should be +/-1 */
		{
			qaws_vec2 inside = { 0, 0 };
			s = qaws_curve_compute_winding_number_2d(&winding, curve, inside);
			TEST_ASSERT_STATUS(s);
			{
				int abs_w = winding < 0 ? -winding : winding;
				TEST_ASSERT(abs_w == 1, "inside winding = +/-1");
			}
		}

		/* Outside point: winding number should be 0 */
		{
			qaws_vec2 outside = { 100, 100 };
			s = qaws_curve_compute_winding_number_2d(&winding, curve, outside);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(winding == 0, "outside winding = 0");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. Non-closed curve => error */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_catmull_rom_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		int winding = 0;
		qaws_vec2 pt = { 1, 0 };

		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = pts;
		desc.control_point_count = 4;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 0;

		s = qaws_curve_create_catmull_rom(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_winding_number_2d(&winding, curve, pt);
		TEST_ASSERT(s != QAWS_STATUS_OK,
			"open curve winding => error");

		qaws_curve_destroy(curve);
	}

	/* 3. Dimension check: 3D curve => INVALID_DIMENSION */
	{
		qaws_vec3 pts[] = { {0, 0, 0}, {1, 1, 0}, {2, 0, 0}, {0, 0, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		int winding = 0;
		qaws_vec2 pt = { 1, 0 };

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_winding_number_2d(&winding, curve, pt);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION,
			"3D curve winding => INVALID_DIMENSION");

		qaws_curve_destroy(curve);
	}

	/* 4. Null args */
	{
		qaws_status s;
		int winding = 0;
		qaws_vec2 pt = { 0, 0 };

		s = qaws_curve_compute_winding_number_2d(&winding, NULL, pt);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve winding => INVALID_ARGUMENT");
	}
}

int test_17_inspection_main(void)
{
	g_pass = 0; g_fail = 0;

	test_inspect();
	test_closest_point();
	test_span_continuity();
	test_family_inspection();
	test_frenet_frame();
	test_geometric_helpers();
	test_inflection_points();
	test_extrema();
	test_curvature_comb();
	test_winding_number();

	printf("17_inspection: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
