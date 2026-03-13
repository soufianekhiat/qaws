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
/*  Helper: make a quadratic 3D Bezier curve                           */
/* ------------------------------------------------------------------ */
static qaws_curve* make_quad_3d(qaws_vec3 a, qaws_vec3 b, qaws_vec3 c_pt)
{
	qaws_curve* crv = NULL;
	qaws_scalar pts[9];
	qaws_bezier_desc d;
	pts[0] = a.x; pts[1] = a.y; pts[2] = a.z;
	pts[3] = b.x; pts[4] = b.y; pts[5] = b.z;
	pts[6] = c_pt.x; pts[7] = c_pt.y; pts[8] = c_pt.z;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 2;
	d.control_points = pts;
	d.control_point_count = 3;
	qaws_curve_create_bezier(&d, &crv);
	return crv;
}

/* ------------------------------------------------------------------ */
/*  Helper: create a dome-shaped Coons patch                           */
/* ------------------------------------------------------------------ */
static qaws_surface* make_dome_coons(
	qaws_curve** out_c0, qaws_curve** out_c1,
	qaws_curve** out_d0, qaws_curve** out_d1)
{
	qaws_surface* surf = NULL;
	qaws_vec3 p00 = {-1,0,0}, p10 = {1,0,0}, p01 = {-1,1,0}, p11 = {1,1,0};
	qaws_vec3 m_bot = {0,0,(qaws_scalar)0.5};
	qaws_vec3 m_top = {0,1,(qaws_scalar)0.5};
	qaws_vec3 m_left = {-1,(qaws_scalar)0.5,(qaws_scalar)0.3};
	qaws_vec3 m_right = {1,(qaws_scalar)0.5,(qaws_scalar)0.3};
	qaws_surface_coons_desc desc;

	*out_c0 = make_quad_3d(p00, m_bot, p10);
	*out_c1 = make_quad_3d(p01, m_top, p11);
	*out_d0 = make_quad_3d(p00, m_left, p01);
	*out_d1 = make_quad_3d(p10, m_right, p11);

	desc.c0 = *out_c0; desc.c1 = *out_c1;
	desc.d0 = *out_d0; desc.d1 = *out_d1;
	qaws_surface_create_coons(&desc, &surf);
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Test: plane-plane intersection                                     */
/* ------------------------------------------------------------------ */
static void test_plane_plane_intersection(void)
{
	qaws_surface* surf_a = NULL;
	qaws_surface* surf_b = NULL;
	qaws_ssi_desc desc;
	qaws_ssi_curve out_curves[4];
	qaws_ssi_point point_buffer[500];
	unsigned int curve_count = 0;
	qaws_status s;
	unsigned int ci, pi;

	printf("test_plane_plane_intersection\n");

	/* Surface A: flat in XY plane at z=0 */
	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = -2; bd.p00.y = -2; bd.p00.z = 0;
		bd.p10.x =  2; bd.p10.y = -2; bd.p10.z = 0;
		bd.p01.x = -2; bd.p01.y =  2; bd.p01.z = 0;
		bd.p11.x =  2; bd.p11.y =  2; bd.p11.z = 0;
		s = qaws_surface_create_bilinear(&bd, &surf_a);
		TEST_ASSERT_STATUS(s);
	}

	/* Surface B: flat in XZ plane at y=0 */
	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = -2; bd.p00.y = 0; bd.p00.z = -2;
		bd.p10.x =  2; bd.p10.y = 0; bd.p10.z = -2;
		bd.p01.x = -2; bd.p01.y = 0; bd.p01.z =  2;
		bd.p11.x =  2; bd.p11.y = 0; bd.p11.z =  2;
		s = qaws_surface_create_bilinear(&bd, &surf_b);
		TEST_ASSERT_STATUS(s);
	}

	TEST_ASSERT(surf_a != NULL, "plane A created");
	TEST_ASSERT(surf_b != NULL, "plane B created");

	/* Intersect: the two planes meet along the X axis (y=0, z=0) */
	memset(&desc, 0, sizeof(desc));
	desc.surface_a = surf_a;
	desc.surface_b = surf_b;
	desc.tolerance = (qaws_scalar)1e-4;
	desc.grid_samples = 10;

	memset(out_curves, 0, sizeof(out_curves));
	memset(point_buffer, 0, sizeof(point_buffer));

	s = qaws_surface_intersect(&desc, out_curves, 4, &curve_count,
		point_buffer, 500);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(curve_count >= 1, "plane-plane: at least 1 intersection curve");

	/* Verify intersection points lie on y~0, z~0 */
	for (ci = 0; ci < curve_count; ci++)
	{
		for (pi = 0; pi < out_curves[ci].point_count; pi++)
		{
			qaws_ssi_point* pt = &out_curves[ci].points[pi];
			TEST_ASSERT(fabs(pt->position.y) < (qaws_scalar)0.05,
				"plane-plane isect y~0");
			TEST_ASSERT(fabs(pt->position.z) < (qaws_scalar)0.05,
				"plane-plane isect z~0");
		}
	}

	qaws_surface_destroy(surf_a);
	qaws_surface_destroy(surf_b);
}

/* ------------------------------------------------------------------ */
/*  Test: dome-plane intersection                                      */
/* ------------------------------------------------------------------ */
static void test_sphere_plane_intersection(void)
{
	qaws_curve* dc0 = NULL, *dc1 = NULL, *dd0 = NULL, *dd1 = NULL;
	qaws_surface* surf_a = NULL;
	qaws_surface* surf_b = NULL;
	qaws_ssi_desc desc;
	qaws_ssi_curve out_curves[4];
	qaws_ssi_point point_buffer[500];
	unsigned int curve_count = 0;
	qaws_status s;
	unsigned int ci, pi;

	printf("test_sphere_plane_intersection\n");

	/* Surface A: dome Coons patch */
	surf_a = make_dome_coons(&dc0, &dc1, &dd0, &dd1);
	TEST_ASSERT(surf_a != NULL, "dome created for intersection");

	/* Surface B: flat bilinear at z=0.3 */
	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = -2; bd.p00.y = -1; bd.p00.z = (qaws_scalar)0.3;
		bd.p10.x =  2; bd.p10.y = -1; bd.p10.z = (qaws_scalar)0.3;
		bd.p01.x = -2; bd.p01.y =  2; bd.p01.z = (qaws_scalar)0.3;
		bd.p11.x =  2; bd.p11.y =  2; bd.p11.z = (qaws_scalar)0.3;
		s = qaws_surface_create_bilinear(&bd, &surf_b);
		TEST_ASSERT_STATUS(s);
	}

	TEST_ASSERT(surf_b != NULL, "flat plane at z=0.3 created");

	/* Intersect */
	memset(&desc, 0, sizeof(desc));
	desc.surface_a = surf_a;
	desc.surface_b = surf_b;
	desc.tolerance = (qaws_scalar)1e-4;
	desc.grid_samples = 10;

	memset(out_curves, 0, sizeof(out_curves));
	memset(point_buffer, 0, sizeof(point_buffer));

	s = qaws_surface_intersect(&desc, out_curves, 4, &curve_count,
		point_buffer, 500);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(curve_count >= 1, "dome-plane: at least 1 intersection curve");

	/* All intersection points should have z ~ 0.3 */
	for (ci = 0; ci < curve_count; ci++)
	{
		for (pi = 0; pi < out_curves[ci].point_count; pi++)
		{
			qaws_ssi_point* pt = &out_curves[ci].points[pi];
			TEST_ASSERT(fabs(pt->position.z - (qaws_scalar)0.3) < (qaws_scalar)0.05,
				"dome-plane isect z~0.3");
		}
	}

	qaws_surface_destroy(surf_a);
	qaws_surface_destroy(surf_b);
	qaws_curve_destroy(dc0);
	qaws_curve_destroy(dc1);
	qaws_curve_destroy(dd0);
	qaws_curve_destroy(dd1);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ surface-surface intersection                           */
/* ------------------------------------------------------------------ */
static void visual_obj_intersection(void)
{
	obj_writer w;
	qaws_curve* dc0 = NULL, *dc1 = NULL, *dd0 = NULL, *dd1 = NULL;
	qaws_surface* surf_a = NULL;
	qaws_surface* surf_b = NULL;
	qaws_ssi_desc desc;
	qaws_ssi_curve out_curves[8];
	qaws_ssi_point point_buffer[1000];
	unsigned int curve_count = 0;
	qaws_status s;
	unsigned int ci, pi;

	printf("visual_obj_intersection\n");
	svg_ensure_output_dir();

	/* Surface A: tall dome Coons patch (-1..1 in x, 0..1 in y, peak z~0.5) */
	surf_a = make_dome_coons(&dc0, &dc1, &dd0, &dd1);
	if (!surf_a) goto cleanup_isect;

	/* Surface B: horizontal plane at z=0.25, clearly slicing through the dome.
	   Extends well beyond the dome so the overlap is visible. */
	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = (qaws_scalar)-1.5; bd.p00.y = (qaws_scalar)-0.5; bd.p00.z = (qaws_scalar)0.25;
		bd.p10.x = (qaws_scalar) 1.5; bd.p10.y = (qaws_scalar)-0.5; bd.p10.z = (qaws_scalar)0.25;
		bd.p01.x = (qaws_scalar)-1.5; bd.p01.y = (qaws_scalar) 1.5; bd.p01.z = (qaws_scalar)0.25;
		bd.p11.x = (qaws_scalar) 1.5; bd.p11.y = (qaws_scalar) 1.5; bd.p11.z = (qaws_scalar)0.25;
		s = qaws_surface_create_bilinear(&bd, &surf_b);
		if (s != QAWS_STATUS_OK || !surf_b) goto cleanup_isect;
	}

	/* Intersect */
	memset(&desc, 0, sizeof(desc));
	desc.surface_a = surf_a;
	desc.surface_b = surf_b;
	desc.tolerance = (qaws_scalar)1e-4;
	desc.grid_samples = 20;

	memset(out_curves, 0, sizeof(out_curves));
	memset(point_buffer, 0, sizeof(point_buffer));

	s = qaws_surface_intersect(&desc, out_curves, 8, &curve_count,
		point_buffer, 1000);
	if (s != QAWS_STATUS_OK) goto cleanup_isect;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/surface_intersection.obj",
		OBJ_OUTPUT_DIR "/surface_intersection.mtl")) goto cleanup_isect;

	/* Surface A: dome */
	obj_material(&w, "surf_a", 0.3, 0.6, 0.9);
	obj_group(&w, "surface_a");
	obj_use_material(&w, "surf_a");
	obj_surface_mesh(&w, surf_a, 32, 32);

	/* Surface B: tilted plane */
	obj_material(&w, "surf_b", 0.9, 0.8, 0.3);
	obj_group(&w, "surface_b");
	obj_use_material(&w, "surf_b");
	obj_surface_mesh(&w, surf_b, 16, 16);

	/* Intersection points as red spheres */
	obj_material(&w, "isect", 1.0, 0.0, 0.0);
	obj_group(&w, "intersection");
	obj_use_material(&w, "isect");
	for (ci = 0; ci < curve_count; ci++)
	{
		for (pi = 0; pi < out_curves[ci].point_count; pi++)
		{
			obj_sphere(&w, out_curves[ci].points[pi].position, (qaws_scalar)0.03);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/surface_intersection.obj\n");

cleanup_isect:
	qaws_surface_destroy(surf_a);
	qaws_surface_destroy(surf_b);
	qaws_curve_destroy(dc0);
	qaws_curve_destroy(dc1);
	qaws_curve_destroy(dd0);
	qaws_curve_destroy(dd1);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_40_surface_intersect_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_plane_plane_intersection();
	test_sphere_plane_intersection();

	/* Visual output */
	visual_obj_intersection();

	printf("40_surface_intersect: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
