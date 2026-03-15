#include "test_common.h"

/* ------------------------------------------------------------------ */
/*  Helper: make a 3D line curve from A to B                           */
/* ------------------------------------------------------------------ */
static qaws_curve* make_line_3d(qaws_vec3 a, qaws_vec3 b)
{
	qaws_curve* c = NULL;
	qaws_scalar pts[6];
	qaws_bezier_desc d;
	pts[0] = a.x; pts[1] = a.y; pts[2] = a.z;
	pts[3] = b.x; pts[4] = b.y; pts[5] = b.z;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 1;
	d.control_points = pts;
	d.control_point_count = 2;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

/* ------------------------------------------------------------------ */
/*  Helper: make a cubic 3D Bezier curve (4 control points)            */
/* ------------------------------------------------------------------ */
static qaws_curve* make_cubic_3d(qaws_vec3 a, qaws_vec3 b,
	qaws_vec3 c_pt, qaws_vec3 d_pt)
{
	qaws_curve* crv = NULL;
	qaws_scalar pts[12];
	qaws_bezier_desc d;
	pts[0]  = a.x;    pts[1]  = a.y;    pts[2]  = a.z;
	pts[3]  = b.x;    pts[4]  = b.y;    pts[5]  = b.z;
	pts[6]  = c_pt.x; pts[7]  = c_pt.y; pts[8]  = c_pt.z;
	pts[9]  = d_pt.x; pts[10] = d_pt.y; pts[11] = d_pt.z;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = pts;
	d.control_point_count = 4;
	qaws_curve_create_bezier(&d, &crv);
	return crv;
}

/* ------------------------------------------------------------------ */
/*  Test: Circular pipe                                                */
/* ------------------------------------------------------------------ */
static void test_circular_pipe(void)
{
	qaws_curve* path = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r0, r1;
	qaws_status s;
	qaws_scalar dist;

	printf("test_circular_pipe\n");

	/* 3D cubic Bezier helix-like path */
	{
		qaws_vec3 p0 = {0, 0, 0};
		qaws_vec3 p1 = {2, 2, 1};
		qaws_vec3 p2 = {4, 0, 2};
		qaws_vec3 p3 = {6, 2, 3};
		path = make_cubic_3d(p0, p1, p2, p3);
	}
	TEST_ASSERT(path != NULL, "pipe path created");

	/* Circular pipe: radius_x = 0.5, radius_y = 0 (defaults to circular) */
	{
		qaws_surface_pipe_desc desc;
		desc.path = path;
		desc.radius_x = (qaws_scalar)0.5;
		desc.radius_y = 0;
		s = qaws_surface_create_pipe(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "pipe surface created");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_PIPE,
		"pipe kind correct");

	/* Evaluate at u=0, v=0 and u=0, v=0.5: these should be at radius 0.5
	   from the path start, in opposite directions */
	memset(&r0, 0, sizeof(r0));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r0);
	TEST_ASSERT_STATUS(s);

	memset(&r1, 0, sizeof(r1));
	s = qaws_surface_evaluate(surf, 0, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r1);
	TEST_ASSERT_STATUS(s);

	/* Distance between these two opposite points should be ~2*radius = 1.0 */
	{
		qaws_scalar dx = r1.position.x - r0.position.x;
		qaws_scalar dy = r1.position.y - r0.position.y;
		qaws_scalar dz = r1.position.z - r0.position.z;
		dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy + dz * dz));
	}
	TEST_ASSERT(approx_eq_loose(dist, (qaws_scalar)1.0), "pipe opposite pts dist ~1.0");

	/* Null/invalid args */
	{
		qaws_surface* tmp = NULL;
		s = qaws_surface_create_pipe(NULL, &tmp);
		TEST_ASSERT(s != QAWS_STATUS_OK, "pipe null desc rejected");
	}

	qaws_surface_destroy(surf);
	qaws_curve_destroy(path);
}

/* ------------------------------------------------------------------ */
/*  Test: Elliptical pipe                                              */
/* ------------------------------------------------------------------ */
static void test_elliptical_pipe(void)
{
	qaws_curve* path = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_elliptical_pipe\n");

	/* Same helix-like path */
	{
		qaws_vec3 p0 = {0, 0, 0};
		qaws_vec3 p1 = {2, 2, 1};
		qaws_vec3 p2 = {4, 0, 2};
		qaws_vec3 p3 = {6, 2, 3};
		path = make_cubic_3d(p0, p1, p2, p3);
	}
	TEST_ASSERT(path != NULL, "elliptical pipe path created");

	/* Elliptical pipe: radius_x=0.5, radius_y=0.3 */
	{
		qaws_surface_pipe_desc desc;
		desc.path = path;
		desc.radius_x = (qaws_scalar)0.5;
		desc.radius_y = (qaws_scalar)0.3;
		s = qaws_surface_create_pipe(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "elliptical pipe surface created");

	/* Check that eval at u=0.5, v=0.25 produces a valid point */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.25,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	/* Point should be somewhere in the vicinity of the path midpoint */
	TEST_ASSERT(r.position.x > (qaws_scalar)-1.0 && r.position.x < (qaws_scalar)8.0,
		"elliptical pipe eval x reasonable");
	TEST_ASSERT(r.position.y > (qaws_scalar)-2.0 && r.position.y < (qaws_scalar)4.0,
		"elliptical pipe eval y reasonable");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(path);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ circular pipe                                          */
/* ------------------------------------------------------------------ */
static void visual_obj_circular_pipe(void)
{
	obj_writer w;
	qaws_curve* path = NULL;
	qaws_surface* surf = NULL;

	printf("visual_obj_circular_pipe\n");
	svg_ensure_output_dir();

	/* Cubic path */
	{
		qaws_vec3 p0 = {0, 0, 0};
		qaws_vec3 p1 = {3, 3, 1};
		qaws_vec3 p2 = {6, 0, 2};
		qaws_vec3 p3 = {9, 3, 3};
		path = make_cubic_3d(p0, p1, p2, p3);
	}
	if (!path) return;

	{
		qaws_surface_pipe_desc desc;
		desc.path = path;
		desc.radius_x = (qaws_scalar)0.4;
		desc.radius_y = 0;
		qaws_surface_create_pipe(&desc, &surf);
	}
	if (!surf) { qaws_curve_destroy(path); return; }

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/35_circular_pipe.obj",
		OBJ_OUTPUT_DIR "/35_circular_pipe.mtl"))
	{
		qaws_surface_destroy(surf); qaws_curve_destroy(path); return;
	}

	obj_material(&w, "pipe", 0.9, 0.6, 0.2);
	obj_group(&w, "circular_pipe");
	obj_use_material(&w, "pipe");
	obj_surface_mesh(&w, surf, 48, 16);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/35_circular_pipe.obj\n");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(path);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ elliptical pipe                                        */
/* ------------------------------------------------------------------ */
static void visual_obj_elliptical_pipe(void)
{
	obj_writer w;
	qaws_curve* path = NULL;
	qaws_surface* surf = NULL;

	printf("visual_obj_elliptical_pipe\n");
	svg_ensure_output_dir();

	/* Same path as circular */
	{
		qaws_vec3 p0 = {0, 0, 0};
		qaws_vec3 p1 = {3, 3, 1};
		qaws_vec3 p2 = {6, 0, 2};
		qaws_vec3 p3 = {9, 3, 3};
		path = make_cubic_3d(p0, p1, p2, p3);
	}
	if (!path) return;

	{
		qaws_surface_pipe_desc desc;
		desc.path = path;
		desc.radius_x = (qaws_scalar)0.6;
		desc.radius_y = (qaws_scalar)0.3;
		qaws_surface_create_pipe(&desc, &surf);
	}
	if (!surf) { qaws_curve_destroy(path); return; }

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/35_elliptical_pipe.obj",
		OBJ_OUTPUT_DIR "/35_elliptical_pipe.mtl"))
	{
		qaws_surface_destroy(surf); qaws_curve_destroy(path); return;
	}

	obj_material(&w, "epipe", 0.3, 0.7, 0.8);
	obj_group(&w, "elliptical_pipe");
	obj_use_material(&w, "epipe");
	obj_surface_mesh(&w, surf, 48, 16);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/35_elliptical_pipe.obj\n");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(path);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_35_surface_pipe_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_circular_pipe();
	test_elliptical_pipe();

	/* Visual output */
	visual_obj_circular_pipe();
	visual_obj_elliptical_pipe();

	printf("35_surface_pipe: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
