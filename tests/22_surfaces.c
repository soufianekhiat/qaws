#include "test_common.h"

static void test_surface_bezier(void)
{
	/* Bilinear Bezier patch: 4 corner points, degree (1,1).
	   S(u,v) = (1-u)(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11
	   Corners: (0,0,0), (2,0,0), (0,3,0), (2,3,1) */
	qaws_surface* surf = NULL;
	qaws_surface_bezier_desc desc;
	qaws_surface_eval_result r;
	qaws_status s;
	qaws_vec3 cp[4];
	qaws_range ur, vr;

	printf("test_surface_bezier\n");

	cp[0].x = 0; cp[0].y = 0; cp[0].z = 0;  /* P[0][0] */
	cp[1].x = 0; cp[1].y = 3; cp[1].z = 0;  /* P[0][1] */
	cp[2].x = 2; cp[2].y = 0; cp[2].z = 0;  /* P[1][0] */
	cp[3].x = 2; cp[3].y = 3; cp[3].z = 1;  /* P[1][1] */

	desc.u_degree = 1;
	desc.v_degree = 1;
	desc.control_points = cp;
	desc.u_point_count = 2;
	desc.v_point_count = 2;

	s = qaws_surface_create_bezier(&desc, &surf);
	TEST_ASSERT(s == QAWS_STATUS_OK, "bezier surface create");
	TEST_ASSERT(surf != NULL, "bezier surface not null");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_BEZIER, "bezier kind");
	TEST_ASSERT(qaws_surface_get_u_degree(surf) == 1, "bezier u degree");
	TEST_ASSERT(qaws_surface_get_v_degree(surf) == 1, "bezier v degree");
	TEST_ASSERT(qaws_surface_is_rational(surf) == 0, "bezier not rational");

	ur = qaws_surface_get_u_range(surf);
	vr = qaws_surface_get_v_range(surf);
	TEST_ASSERT(ur.min_value == 0 && ur.max_value == 1, "bezier u range");
	TEST_ASSERT(vr.min_value == 0 && vr.max_value == 1, "bezier v range");

	/* Evaluate at corner (0,0) -> P00 = (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0, QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "bezier eval (0,0)");
	TEST_ASSERT(fabs(r.position.x) < 0.01f && fabs(r.position.y) < 0.01f
		&& fabs(r.position.z) < 0.01f, "bezier pos (0,0) = origin");

	/* Evaluate at corner (1,1) -> P11 = (2,3,1) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1, QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "bezier eval (1,1)");
	TEST_ASSERT(fabs(r.position.x - 2) < 0.01f && fabs(r.position.y - 3) < 0.01f
		&& fabs(r.position.z - 1) < 0.01f, "bezier pos (1,1) = (2,3,1)");

	/* Evaluate at center (0.5,0.5) -> average = (1, 1.5, 0.25) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "bezier eval center");
	TEST_ASSERT(fabs(r.position.x - 1) < 0.01f && fabs(r.position.y - 1.5f) < 0.01f,
		"bezier center position");
	TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_NORMAL, "bezier normal computed");

	/* Sample grid */
	{
		qaws_vec3 pts[25];
		unsigned int count = 0;
		s = qaws_surface_sample(surf, 5, 5, pts, 25, &count);
		TEST_ASSERT(s == QAWS_STATUS_OK, "bezier sample ok");
		TEST_ASSERT(count == 25, "bezier sample count 25");
	}

	qaws_surface_destroy(surf);
}

static void test_surface_bspline(void)
{
	/* Bicubic B-spline surface: 5x5 control points, degree (3,3).
	   Create a simple hill shape. */
	qaws_surface* surf = NULL;
	qaws_surface_bspline_desc desc;
	qaws_surface_eval_result r;
	qaws_status s;
	qaws_vec3 cp[25]; /* 5x5 */
	unsigned int i, j;

	printf("test_surface_bspline\n");

	for (i = 0; i < 5; i++)
	{
		for (j = 0; j < 5; j++)
		{
			unsigned int idx = i * 5 + j;
			cp[idx].x = (qaws_scalar)i;
			cp[idx].y = (qaws_scalar)j;
			/* Smooth hill: z peaks at center */
			cp[idx].z = (qaws_scalar)(2.0 * exp(-0.5 * ((i - 2.0) * (i - 2.0) + (j - 2.0) * (j - 2.0))));
		}
	}

	memset(&desc, 0, sizeof(desc));
	desc.u_degree = 3;
	desc.v_degree = 3;
	desc.control_points = cp;
	desc.u_point_count = 5;
	desc.v_point_count = 5;
	desc.u_knot_count = 0; /* auto uniform */
	desc.v_knot_count = 0;

	s = qaws_surface_create_bspline(&desc, &surf);
	TEST_ASSERT(s == QAWS_STATUS_OK, "bspline surface create");
	TEST_ASSERT(surf != NULL, "bspline surface not null");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_BSPLINE, "bspline kind");
	TEST_ASSERT(qaws_surface_get_u_degree(surf) == 3, "bspline u degree 3");
	TEST_ASSERT(qaws_surface_is_rational(surf) == 0, "bspline not rational");

	/* Evaluate position + derivatives at domain center */
	{
		qaws_range ur = qaws_surface_get_u_range(surf);
		qaws_range vr = qaws_surface_get_v_range(surf);
		qaws_scalar um = (ur.min_value + ur.max_value) * (qaws_scalar)0.5;
		qaws_scalar vm = (vr.min_value + vr.max_value) * (qaws_scalar)0.5;

		memset(&r, 0, sizeof(r));
		s = qaws_surface_evaluate(surf, um, vm,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU
			| QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL, &r);
		TEST_ASSERT(s == QAWS_STATUS_OK, "bspline eval center ok");
		TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_POSITION, "bspline pos valid");
		TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_NORMAL, "bspline normal valid");
		/* Z should be positive at the hill */
		TEST_ASSERT(r.position.z > 0, "bspline z > 0 at hill");
	}

	/* Sample grid */
	{
		qaws_vec3 pts[100];
		unsigned int count = 0;
		s = qaws_surface_sample(surf, 10, 10, pts, 100, &count);
		TEST_ASSERT(s == QAWS_STATUS_OK, "bspline sample ok");
		TEST_ASSERT(count == 100, "bspline sample count 100");
	}

	qaws_surface_destroy(surf);
}

static void test_surface_nurbs(void)
{
	/* NURBS sphere octant: degree (2,2) with 3x3 CPs and weights.
	   This tests rational surface evaluation. */
	qaws_surface* surf = NULL;
	qaws_surface_nurbs_desc desc;
	qaws_surface_eval_result r;
	qaws_status s;
	qaws_scalar w = (qaws_scalar)0.707106781f; /* 1/sqrt(2) */

	/* 3x3 control points for a quarter sphere octant */
	qaws_vec3 cp[9] = {
		{0, 0, 1}, {0, 0, w}, {0, 0, 0},    /* row 0 */
		{0, 0, w}, {0, 0, 0.5f}, {0, w, 0},  /* row 1 */
		{1, 0, 0}, {w, 0, 0}, {0, 1, 0}      /* row 2 */
	};
	qaws_scalar weights[9] = {
		1,  w,  1,
		w,  0.5f, w,
		1,  w,  1
	};

	printf("test_surface_nurbs\n");

	memset(&desc, 0, sizeof(desc));
	desc.u_degree = 2;
	desc.v_degree = 2;
	desc.control_points = cp;
	desc.u_point_count = 3;
	desc.v_point_count = 3;
	desc.weights = weights;
	desc.u_knot_count = 0; /* auto uniform */
	desc.v_knot_count = 0;

	s = qaws_surface_create_nurbs(&desc, &surf);
	TEST_ASSERT(s == QAWS_STATUS_OK, "nurbs surface create");
	TEST_ASSERT(surf != NULL, "nurbs surface not null");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_NURBS, "nurbs kind");
	TEST_ASSERT(qaws_surface_is_rational(surf) == 1, "nurbs is rational");

	/* Evaluate at domain center */
	{
		qaws_range ur = qaws_surface_get_u_range(surf);
		qaws_range vr = qaws_surface_get_v_range(surf);
		qaws_scalar um = (ur.min_value + ur.max_value) * (qaws_scalar)0.5;
		qaws_scalar vm = (vr.min_value + vr.max_value) * (qaws_scalar)0.5;

		memset(&r, 0, sizeof(r));
		s = qaws_surface_evaluate(surf, um, vm,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
		TEST_ASSERT(s == QAWS_STATUS_OK, "nurbs eval center");
		TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_POSITION, "nurbs pos valid");
		TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_NORMAL, "nurbs normal valid");
	}

	/* Evaluate at corners */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0, QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "nurbs eval corner (0,0)");

	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1, QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "nurbs eval corner (1,1)");

	/* Sample grid */
	{
		qaws_vec3 pts[64];
		unsigned int count = 0;
		s = qaws_surface_sample(surf, 8, 8, pts, 64, &count);
		TEST_ASSERT(s == QAWS_STATUS_OK, "nurbs sample ok");
		TEST_ASSERT(count == 64, "nurbs sample count 64");
	}

	qaws_surface_destroy(surf);
}

static void test_surface_error_handling(void)
{
	qaws_surface* surf = NULL;
	qaws_surface_bezier_desc bdesc;
	qaws_surface_bspline_desc sdesc;
	qaws_surface_nurbs_desc ndesc;
	qaws_status s;

	printf("test_surface_error_handling\n");

	/* Null desc */
	s = qaws_surface_create_bezier(NULL, &surf);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "bezier null desc");

	/* Null output */
	memset(&bdesc, 0, sizeof(bdesc));
	s = qaws_surface_create_bezier(&bdesc, NULL);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "bezier null out");

	/* Bad degree */
	bdesc.u_degree = 0;
	bdesc.v_degree = 1;
	s = qaws_surface_create_bezier(&bdesc, &surf);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_DEGREE, "bezier bad degree");

	/* Bad point count */
	{
		qaws_vec3 cp[4] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0}};
		bdesc.u_degree = 1;
		bdesc.v_degree = 1;
		bdesc.control_points = cp;
		bdesc.u_point_count = 3; /* wrong: should be 2 for degree 1 */
		bdesc.v_point_count = 2;
		s = qaws_surface_create_bezier(&bdesc, &surf);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT, "bezier bad cp count");
	}

	/* B-spline null desc */
	s = qaws_surface_create_bspline(NULL, &surf);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "bspline null desc");

	/* B-spline too few CPs */
	memset(&sdesc, 0, sizeof(sdesc));
	{
		qaws_vec3 cp2[4] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0}};
		sdesc.u_degree = 3;
		sdesc.v_degree = 3;
		sdesc.control_points = cp2;
		sdesc.u_point_count = 2;
		sdesc.v_point_count = 2;
		s = qaws_surface_create_bspline(&sdesc, &surf);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT, "bspline too few CPs");
	}

	/* NURBS null weights */
	memset(&ndesc, 0, sizeof(ndesc));
	{
		qaws_vec3 cp3[4] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0}};
		ndesc.u_degree = 1;
		ndesc.v_degree = 1;
		ndesc.control_points = cp3;
		ndesc.u_point_count = 2;
		ndesc.v_point_count = 2;
		ndesc.weights = NULL;
		s = qaws_surface_create_nurbs(&ndesc, &surf);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "nurbs null weights");
	}

	/* Null surface queries */
	TEST_ASSERT(qaws_surface_get_kind(NULL) == QAWS_SURFACE_KIND_INVALID, "null kind");
	TEST_ASSERT(qaws_surface_get_u_degree(NULL) == 0, "null u_degree");
	TEST_ASSERT(qaws_surface_get_v_degree(NULL) == 0, "null v_degree");
	TEST_ASSERT(qaws_surface_is_rational(NULL) == 0, "null rational");
}

static void test_surface_ruled(void)
{
	/* Ruled surface between two 3D Bezier curves */
	qaws_curve* ca = NULL;
	qaws_curve* cb = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_ruled_desc rdesc;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_surface_ruled\n");

	/* Curve A: straight line from (0,0,0) to (4,0,0) */
	{
		qaws_scalar pts_a[] = {0,0,0, 4,0,0};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 1;
		d.control_points = pts_a;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&ca, &d);
	}

	/* Curve B: arc from (0,0,2) through (2,2,2) to (4,0,2) */
	{
		qaws_scalar pts_b[] = {0,0,2, 1,2,2, 3,2,2, 4,0,2};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = pts_b;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&cb, &d);
	}

	rdesc.curve_a = ca;
	rdesc.curve_b = cb;
	s = qaws_surface_create_ruled(&rdesc, &surf);
	TEST_ASSERT(s == QAWS_STATUS_OK, "ruled create ok");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_RULED, "ruled kind");

	/* At v=0 should be on curve A */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0, QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "ruled eval (0,0)");
	TEST_ASSERT(fabs(r.position.x) < 0.01f && fabs(r.position.z) < 0.01f,
		"ruled v=0 on curve A start");

	/* At v=1 should be on curve B */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 1, QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "ruled eval (0,1)");
	TEST_ASSERT(fabs(r.position.z - 2) < 0.01f, "ruled v=1 on curve B");

	/* dv should be B-A direction (upward) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU
		| QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "ruled eval center");
	TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_NORMAL, "ruled normal valid");
	/* Position z should be ~1 at v=0.5 */
	TEST_ASSERT(fabs(r.position.z - 1) < 0.2f, "ruled midpoint z ~1");

	/* dvv should be 0 (linear in v) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_DVV, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "ruled dvv eval");
	TEST_ASSERT(fabs(r.dvv.x) < 0.01f && fabs(r.dvv.y) < 0.01f
		&& fabs(r.dvv.z) < 0.01f, "ruled dvv = 0");

	/* Sample grid */
	{
		qaws_vec3 pts[100];
		unsigned int count = 0;
		s = qaws_surface_sample(surf, 10, 10, pts, 100, &count);
		TEST_ASSERT(s == QAWS_STATUS_OK, "ruled sample ok");
		TEST_ASSERT(count == 100, "ruled sample count 100");
	}

	qaws_surface_destroy(surf);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
}

static void test_surface_swept(void)
{
	/* Swept surface: circle cross-section swept along a 3D Bezier path */
	qaws_curve* path = NULL;
	qaws_curve* profile = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_swept_desc sdesc;
	qaws_surface_eval_result r;
	qaws_status s;
	unsigned int i;

	printf("test_surface_swept\n");

	/* Path: cubic Bezier curve going from (0,0,0) to (4,0,3) with some curvature */
	{
		qaws_scalar pts[] = {0,0,0, 1,2,1, 3,2,2, 4,0,3};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = pts;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&path, &d);
	}

	/* Profile: circle approximation using 2D cubic Bezier */
	{
		qaws_scalar pts[26]; /* 13 points x 2 */
		qaws_bezier_desc d;
		for (i = 0; i < 13; i++)
		{
			double angle = 2.0 * M_PI * (double)i / 12.0;
			pts[i * 2 + 0] = (qaws_scalar)cos(angle);
			pts[i * 2 + 1] = (qaws_scalar)sin(angle);
		}
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 12;
		d.control_points = pts;
		d.control_point_count = 13;
		qaws_curve_create_bezier(&profile, &d);
	}

	memset(&sdesc, 0, sizeof(sdesc));
	sdesc.path = path;
	sdesc.profile = profile;
	sdesc.scale = (qaws_scalar)0.3;

	s = qaws_surface_create_swept(&sdesc, &surf);
	TEST_ASSERT(s == QAWS_STATUS_OK, "swept create ok");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_SWEPT, "swept kind");

	/* Evaluate at path start, profile center-ish */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK, "swept eval ok");
	TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_POSITION, "swept pos valid");

	/* Sample grid */
	{
		qaws_vec3 pts[200];
		unsigned int count = 0;
		s = qaws_surface_sample(surf, 20, 10, pts, 200, &count);
		TEST_ASSERT(s == QAWS_STATUS_OK, "swept sample ok");
		TEST_ASSERT(count == 200, "swept sample count 200");
	}

	qaws_surface_destroy(surf);
	qaws_curve_destroy(path);
	qaws_curve_destroy(profile);
}

int test_22_surfaces_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_surface_bezier();
	test_surface_bspline();
	test_surface_nurbs();
	test_surface_error_handling();
	test_surface_ruled();
	test_surface_swept();

	printf("22_surfaces: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
