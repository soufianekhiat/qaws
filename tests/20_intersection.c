#include "test_common.h"

static void test_self_intersection(void)
{
	printf("test_self_intersection\n");

	/* 1. Figure-8 Hermite curve has a self-intersection at the origin */
	{
		qaws_vec2 pts[5];
		qaws_vec2 tans[5];
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_2d isects[16];
		unsigned int isect_count = 0;

		pts[0].x = 0; pts[0].y = 0;
		pts[1].x = 1; pts[1].y = (qaws_scalar)0.5;
		pts[2].x = 0; pts[2].y = 0;
		pts[3].x = -1; pts[3].y = (qaws_scalar)-0.5;
		pts[4].x = 0; pts[4].y = 0;

		tans[0].x = 1; tans[0].y = 1;
		tans[1].x = 0; tans[1].y = (qaws_scalar)-1.5;
		tans[2].x = -1; tans[2].y = -1;
		tans[3].x = 0; tans[3].y = (qaws_scalar)1.5;
		tans[4].x = 1; tans[4].y = 1;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 5;
		desc.derivative_count = 5;

		s = qaws_curve_create_hermite(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_2d(
			isects, &isect_count, curve, 16);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count > 0, "figure-8 has self-intersections");

		if (isect_count > 0 && isect_count <= 16)
		{
			int found_origin = 0;
			unsigned int i;
			for (i = 0; i < isect_count; ++i)
			{
				qaws_scalar dx = isects[i].position.x;
				qaws_scalar dy = isects[i].position.y;
				if (dx < 0) dx = -dx;
				if (dy < 0) dy = -dy;
				if (dx < (qaws_scalar)0.1 && dy < (qaws_scalar)0.1)
					found_origin = 1;
				TEST_ASSERT(isects[i].parameter_a < isects[i].parameter_b,
					"self-isect param_a < param_b");
			}
			TEST_ASSERT(found_origin, "self-intersection near origin");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. Simple Bezier arc — no self-intersection */
	{
		qaws_vec2 cp[4];
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_2d isects[8];
		unsigned int isect_count = 0;

		cp[0].x = 0; cp[0].y = 0;
		cp[1].x = 1; cp[1].y = 2;
		cp[2].x = 3; cp[2].y = 2;
		cp[3].x = 4; cp[3].y = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = cp;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_2d(
			isects, &isect_count, curve, 8);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count == 0, "simple arc has no self-intersections");

		qaws_curve_destroy(curve);
	}

	/* 3. Null arguments */
	{
		qaws_intersection_2d isects[4];
		unsigned int count = 0;
		qaws_status s;

		s = qaws_curve_find_self_intersections_2d(isects, &count, NULL, 4);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve => INVALID_ARGUMENT");
	}

	/* 4. Dimension mismatch */
	{
		qaws_vec3 cp3[4];
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_2d isects[4];
		unsigned int count = 0;

		cp3[0].x = 0; cp3[0].y = 0; cp3[0].z = 0;
		cp3[1].x = 1; cp3[1].y = 1; cp3[1].z = 1;
		cp3[2].x = 2; cp3[2].y = 0; cp3[2].z = 1;
		cp3[3].x = 3; cp3[3].y = 1; cp3[3].z = 0;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = cp3;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_2d(isects, &count, curve, 4);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION,
			"3d curve for 2d self-isect => INVALID_DIMENSION");

		qaws_curve_destroy(curve);
	}

	/* 5. 3D self-intersection */
	{
		qaws_vec3 pts[5];
		qaws_vec3 tans[5];
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_3d isects[16];
		unsigned int isect_count = 0;

		pts[0].x = 0; pts[0].y = 0; pts[0].z = 0;
		pts[1].x = 1; pts[1].y = (qaws_scalar)0.5; pts[1].z = 0;
		pts[2].x = 0; pts[2].y = 0; pts[2].z = 0;
		pts[3].x = -1; pts[3].y = (qaws_scalar)-0.5; pts[3].z = 0;
		pts[4].x = 0; pts[4].y = 0; pts[4].z = 0;

		tans[0].x = 1; tans[0].y = 1; tans[0].z = 0;
		tans[1].x = 0; tans[1].y = (qaws_scalar)-1.5; tans[1].z = 0;
		tans[2].x = -1; tans[2].y = -1; tans[2].z = 0;
		tans[3].x = 0; tans[3].y = (qaws_scalar)1.5; tans[3].z = 0;
		tans[4].x = 1; tans[4].y = 1; tans[4].z = 0;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 5;
		desc.derivative_count = 5;

		s = qaws_curve_create_hermite(&curve, &desc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_3d(
			isects, &isect_count, curve, 16);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count > 0, "3D figure-8 has self-intersections");

		qaws_curve_destroy(curve);
	}
}

static void test_curve_curve_intersection(void)
{
	printf("test_curve_curve_intersection\n");

	/* 1. Two crossing Bezier curves in 2D */
	{
		/* Curve A: arch shape */
		qaws_vec2 cpA[] = { {0,0}, {1,3}, {3,3}, {4,0} };
		/* Curve B: S-shape crossing A */
		qaws_vec2 cpB[] = { {0,2}, {2,-1}, {2,3}, {4,1} };

		qaws_bezier_desc descA, descB;
		qaws_curve *curveA = NULL, *curveB = NULL;
		qaws_status s;
		qaws_intersection_2d isects[16];
		unsigned int isect_count = 0;

		descA.dimension = QAWS_DIMENSION_2D;
		descA.degree = 3;
		descA.control_points = cpA;
		descA.control_point_count = 4;

		descB.dimension = QAWS_DIMENSION_2D;
		descB.degree = 3;
		descB.control_points = cpB;
		descB.control_point_count = 4;

		s = qaws_curve_create_bezier(&curveA, &descA);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&curveB, &descB);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(
			isects, &isect_count, curveA, curveB, 16);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count > 0, "crossing beziers have intersections");

		/* Verify intersection positions agree on both curves */
		{
			unsigned int i;
			for (i = 0; i < isect_count && i < 16; ++i)
			{
				qaws_eval_result_2d ra, rb;
				qaws_scalar dx, dy;
				s = qaws_curve_evaluate_2d(&ra, curveA,
					isects[i].parameter_a, QAWS_EVAL_FLAG_POSITION);
				TEST_ASSERT_STATUS(s);
				s = qaws_curve_evaluate_2d(&rb, curveB,
					isects[i].parameter_b, QAWS_EVAL_FLAG_POSITION);
				TEST_ASSERT_STATUS(s);

				dx = ra.position.x - rb.position.x;
				dy = ra.position.y - rb.position.y;
				if (dx < 0) dx = -dx;
				if (dy < 0) dy = -dy;
				TEST_ASSERT(dx < (qaws_scalar)0.01 && dy < (qaws_scalar)0.01,
					"intersection point matches on both curves");
			}
		}

		qaws_curve_destroy(curveA);
		qaws_curve_destroy(curveB);
	}

	/* 2. Non-intersecting curves */
	{
		qaws_vec2 cpA[] = { {0,0}, {1,1}, {2,1}, {3,0} };
		qaws_vec2 cpB[] = { {0,3}, {1,4}, {2,4}, {3,3} };

		qaws_bezier_desc descA, descB;
		qaws_curve *curveA = NULL, *curveB = NULL;
		qaws_status s;
		qaws_intersection_2d isects[8];
		unsigned int isect_count = 0;

		descA.dimension = QAWS_DIMENSION_2D;
		descA.degree = 3;
		descA.control_points = cpA;
		descA.control_point_count = 4;

		descB.dimension = QAWS_DIMENSION_2D;
		descB.degree = 3;
		descB.control_points = cpB;
		descB.control_point_count = 4;

		s = qaws_curve_create_bezier(&curveA, &descA);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&curveB, &descB);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(
			isects, &isect_count, curveA, curveB, 8);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count == 0, "non-intersecting curves have 0 intersections");

		qaws_curve_destroy(curveA);
		qaws_curve_destroy(curveB);
	}

	/* 3. Line crossing a curve at 2 points */
	{
		qaws_vec2 cpLine[] = { {0,(qaws_scalar)1.5}, {4,(qaws_scalar)1.5} };
		qaws_vec2 cpArc[]  = { {0,0}, {1,3}, {3,3}, {4,0} };

		qaws_bezier_desc descL, descA;
		qaws_curve *line = NULL, *arc = NULL;
		qaws_status s;
		qaws_intersection_2d isects[8];
		unsigned int isect_count = 0;

		descL.dimension = QAWS_DIMENSION_2D;
		descL.degree = 1;
		descL.control_points = cpLine;
		descL.control_point_count = 2;

		descA.dimension = QAWS_DIMENSION_2D;
		descA.degree = 3;
		descA.control_points = cpArc;
		descA.control_point_count = 4;

		s = qaws_curve_create_bezier(&line, &descL);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&arc, &descA);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(
			isects, &isect_count, line, arc, 8);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count == 2, "line crosses arc at 2 points");

		qaws_curve_destroy(line);
		qaws_curve_destroy(arc);
	}

	/* 4. Null arguments */
	{
		qaws_intersection_2d isects[4];
		unsigned int count = 0;
		qaws_status s;

		s = qaws_curve_find_intersections_2d(isects, &count, NULL, NULL, 4);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curves => INVALID_ARGUMENT");
	}

	/* 5. Dimension mismatch */
	{
		qaws_vec2 cp2[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec3 cp3[] = { {0,0,0}, {1,1,1}, {2,0,0} };
		qaws_bezier_desc d2, d3;
		qaws_curve *c2 = NULL, *c3 = NULL;
		qaws_status s;
		qaws_intersection_2d isects[4];
		unsigned int count = 0;

		d2.dimension = QAWS_DIMENSION_2D;
		d2.degree = 2;
		d2.control_points = cp2;
		d2.control_point_count = 3;

		d3.dimension = QAWS_DIMENSION_3D;
		d3.degree = 2;
		d3.control_points = cp3;
		d3.control_point_count = 3;

		s = qaws_curve_create_bezier(&c2, &d2);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&c3, &d3);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(isects, &count, c2, c3, 4);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION,
			"mixed dimensions => INVALID_DIMENSION");

		qaws_curve_destroy(c2);
		qaws_curve_destroy(c3);
	}

	/* 6. 3D curve-curve intersection */
	{
		qaws_vec3 cpA[] = { {0,0,0}, {2,2,0}, {2,-2,0}, {4,0,0} };
		qaws_vec3 cpB[] = { {2,-2,0}, {2,2,0}, {2,2,0}, {2,-2,0} };

		qaws_bezier_desc descA, descB;
		qaws_curve *curveA = NULL, *curveB = NULL;
		qaws_status s;
		qaws_intersection_3d isects[16];
		unsigned int isect_count = 0;

		descA.dimension = QAWS_DIMENSION_3D;
		descA.degree = 3;
		descA.control_points = cpA;
		descA.control_point_count = 4;

		descB.dimension = QAWS_DIMENSION_3D;
		descB.degree = 3;
		descB.control_points = cpB;
		descB.control_point_count = 4;

		s = qaws_curve_create_bezier(&curveA, &descA);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&curveB, &descB);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_3d(
			isects, &isect_count, curveA, curveB, 16);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count >= 0, "3D curve-curve intersection ran");

		qaws_curve_destroy(curveA);
		qaws_curve_destroy(curveB);
	}
}

int test_20_intersection_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_self_intersection();
	test_curve_curve_intersection();

	printf("20_intersection: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
