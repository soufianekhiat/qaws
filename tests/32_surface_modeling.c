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
/*  Helper: make a 2D circle profile (quadratic Bezier approximation)  */
/* ------------------------------------------------------------------ */
static qaws_curve* make_circle_profile_2d(qaws_scalar radius)
{
	qaws_curve* crv = NULL;
	qaws_vec2 pts[3];
	qaws_bezier_desc d;
	pts[0].x = radius; pts[0].y = 0;
	pts[1].x = radius; pts[1].y = radius;
	pts[2].x = 0;      pts[2].y = radius;
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = pts;
	d.control_point_count = 3;
	qaws_curve_create_bezier(&d, &crv);
	return crv;
}

/* ------------------------------------------------------------------ */
/*  Test: Coons patch                                                  */
/* ------------------------------------------------------------------ */
static void test_coons_patch(void)
{
	qaws_curve* c0 = NULL; /* bottom: v=0 */
	qaws_curve* c1 = NULL; /* top: v=1 */
	qaws_curve* d0 = NULL; /* left: u=0 */
	qaws_curve* d1 = NULL; /* right: u=1 */
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_coons_patch\n");

	/* Flat quad: corners at (0,0,0), (4,0,0), (0,3,0), (4,3,0) */
	{
		qaws_vec3 p00 = {0,0,0}, p10 = {4,0,0}, p01 = {0,3,0}, p11 = {4,3,0};
		c0 = make_line_3d(p00, p10); /* bottom */
		c1 = make_line_3d(p01, p11); /* top */
		d0 = make_line_3d(p00, p01); /* left */
		d1 = make_line_3d(p10, p11); /* right */
	}

	{
		qaws_surface_coons_desc desc;
		desc.c0 = c0; desc.c1 = c1;
		desc.d0 = d0; desc.d1 = d1;
		s = qaws_surface_create_coons(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "coons surface created");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_COONS,
		"coons kind correct");

	/* Corner (0,0) should be (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "coons (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "coons (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "coons (0,0) z");

	/* Corner (1,1) should be (4,3,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "coons (1,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "coons (1,1) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "coons (1,1) z");

	/* Center (0.5,0.5) of flat quad should be (2, 1.5, 0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)2.0), "coons center x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.5), "coons center y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "coons center z");

	/* Normal should point in z direction for flat patch */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.normal.z) > (qaws_scalar)0.9, "coons flat normal ~z");

	/* Null args */
	s = qaws_surface_create_coons(NULL, &surf);
	TEST_ASSERT(s != QAWS_STATUS_OK, "coons null desc rejected");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);

	/* Coons with curved boundaries */
	{
		qaws_vec3 p00 = {0,0,0}, p10 = {4,0,0}, p01 = {0,3,0}, p11 = {4,3,0};
		qaws_vec3 m_bot = {2,0,1};  /* raised bottom mid */
		qaws_vec3 m_top = {2,3,1};  /* raised top mid */
		qaws_surface_coons_desc desc;

		c0 = make_quad_3d(p00, m_bot, p10);
		c1 = make_quad_3d(p01, m_top, p11);
		d0 = make_line_3d(p00, p01);
		d1 = make_line_3d(p10, p11);

		desc.c0 = c0; desc.c1 = c1;
		desc.d0 = d0; desc.d1 = d1;
		surf = NULL;
		s = qaws_surface_create_coons(&desc, &surf);
		TEST_ASSERT_STATUS(s);

		/* Center should be raised in z due to curved boundaries */
		memset(&r, 0, sizeof(r));
		qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
			QAWS_SURFACE_EVAL_POSITION, &r);
		TEST_ASSERT(r.position.z > (qaws_scalar)0.1, "curved coons center raised");

		qaws_surface_destroy(surf);
		qaws_curve_destroy(c0);
		qaws_curve_destroy(c1);
		qaws_curve_destroy(d0);
		qaws_curve_destroy(d1);
	}
}

/* ------------------------------------------------------------------ */
/*  Test: Extrusion surface                                            */
/* ------------------------------------------------------------------ */
static void test_extrusion_surface(void)
{
	qaws_curve* profile = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_extrusion_surface\n");

	/* 2D quadratic profile extruded in Z */
	{
		qaws_vec2 pts[] = { {0,0}, {2,1}, {4,0} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 2;
		d.control_points = pts;
		d.control_point_count = 3;
		qaws_curve_create_bezier(&d, &profile);
	}

	{
		qaws_surface_extrusion_desc desc;
		desc.profile = profile;
		desc.direction.x = 0; desc.direction.y = 0; desc.direction.z = 1;
		desc.length = (qaws_scalar)5.0;
		s = qaws_surface_create_extrusion(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "extrusion surface created");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_EXTRUSION,
		"extrusion kind correct");

	/* At v=0, should be the profile at z=0 */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "extr (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "extr (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "extr (0,0) z");

	/* At v=1, should be the profile at z=5 */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "extr (0,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)5.0), "extr (0,1) z=5");

	/* At u=1, v=0, should be end of profile (4,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "extr (1,0) x=4");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "extr (1,0) z=0");

	/* dv should be the direction vector */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.dv.z, (qaws_scalar)5.0), "extr dv.z = length");
	TEST_ASSERT(fabs(r.dv.x) < TOLERANCE_LOOSE, "extr dv.x = 0");
	TEST_ASSERT(fabs(r.dv.y) < TOLERANCE_LOOSE, "extr dv.y = 0");

	/* Null args */
	{
		qaws_surface* tmp = NULL;
		s = qaws_surface_create_extrusion(NULL, &tmp);
		TEST_ASSERT(s != QAWS_STATUS_OK, "extrusion null desc rejected");
	}

	qaws_surface_destroy(surf);
	qaws_curve_destroy(profile);
}

/* ------------------------------------------------------------------ */
/*  Test: Revolution surface                                           */
/* ------------------------------------------------------------------ */
static void test_revolution_surface(void)
{
	qaws_curve* profile = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;
	qaws_scalar pi = (qaws_scalar)(3.14159265358979323846);

	printf("test_revolution_surface\n");

	/* Profile in XZ plane: a line from (1,0) to (1,2) -- cylinder radius 1, height 2 */
	{
		qaws_vec2 pts[] = { {1, 0}, {1, 2} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 1;
		d.control_points = pts;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &profile);
	}

	{
		qaws_surface_revolution_desc desc;
		desc.profile = profile;
		desc.axis_origin.x = 0; desc.axis_origin.y = 0; desc.axis_origin.z = 0;
		desc.axis_direction.x = 0; desc.axis_direction.y = 0; desc.axis_direction.z = 1;
		desc.angle = 2 * pi;  /* full revolution */
		s = qaws_surface_create_revolution(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "revolution surface created");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_REVOLUTION,
		"revolution kind correct");

	/* At u=0, v=0: profile start at angle 0 -> (1, 0, 0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)1.0), "rev (0,0) x=1");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "rev (0,0) y=0");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "rev (0,0) z=0");

	/* At u=0.25 (90 degrees), v=0: should be (0, 1, 0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.25, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "rev 90deg x~0");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.0), "rev 90deg y=1");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "rev 90deg z=0");

	/* At u=0, v=1: should be at height 2 -> (1, 0, 2) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)1.0), "rev top x=1");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "rev top y=0");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)2.0), "rev top z=2");

	/* Points should all be at radius 1 from the z axis */
	{
		unsigned int i;
		for (i = 0; i < 8; i++)
		{
			qaws_scalar u = (qaws_scalar)i / (qaws_scalar)8;
			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surf, u, 0, QAWS_SURFACE_EVAL_POSITION, &r);
			{
				qaws_scalar rad = (qaws_scalar)sqrt((double)(r.position.x * r.position.x +
					r.position.y * r.position.y));
				TEST_ASSERT(approx_eq_loose(rad, (qaws_scalar)1.0), "rev radius=1");
			}
		}
	}

	/* Normal should point outward (away from z axis) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
	TEST_ASSERT_STATUS(s);
	/* At angle 0, normal should have significant x component */
	TEST_ASSERT(fabs(r.normal.x) > (qaws_scalar)0.5 ||
		fabs(r.normal.y) > (qaws_scalar)0.5,
		"rev normal points outward");

	/* Null args */
	{
		qaws_surface* tmp = NULL;
		s = qaws_surface_create_revolution(NULL, &tmp);
		TEST_ASSERT(s != QAWS_STATUS_OK, "revolution null desc rejected");
	}

	qaws_surface_destroy(surf);
	qaws_curve_destroy(profile);

	/* Half revolution (semicircle cross-section) */
	{
		qaws_vec2 pts[] = { {2, 0}, {2, 3} };
		qaws_bezier_desc d;
		qaws_surface_revolution_desc rd;

		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 1;
		d.control_points = pts;
		d.control_point_count = 2;
		profile = NULL;
		qaws_curve_create_bezier(&d, &profile);

		rd.profile = profile;
		rd.axis_origin.x = 0; rd.axis_origin.y = 0; rd.axis_origin.z = 0;
		rd.axis_direction.x = 0; rd.axis_direction.y = 0; rd.axis_direction.z = 1;
		rd.angle = pi;  /* half revolution */
		surf = NULL;
		s = qaws_surface_create_revolution(&rd, &surf);
		TEST_ASSERT_STATUS(s);

		/* At u=0.5 (90 degrees), v=0: should be (0, 2, 0) */
		memset(&r, 0, sizeof(r));
		qaws_surface_evaluate(surf, (qaws_scalar)0.5, 0,
			QAWS_SURFACE_EVAL_POSITION, &r);
		TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "half rev 90 x~0");
		/* At u=1 (180 degrees), should be (-2, 0, 0) */
		memset(&r, 0, sizeof(r));
		qaws_surface_evaluate(surf, 1, 0,
			QAWS_SURFACE_EVAL_POSITION, &r);
		TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)-2.0), "half rev 180 x=-2");

		qaws_surface_destroy(surf);
		qaws_curve_destroy(profile);
	}
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ Coons patch                                            */
/* ------------------------------------------------------------------ */
static void visual_obj_coons(void)
{
	obj_writer w;
	qaws_curve* c0 = NULL, *c1 = NULL, *d0 = NULL, *d1 = NULL;
	qaws_surface* surf = NULL;
	qaws_vec3 p00 = {0,0,0}, p10 = {4,0,0}, p01 = {0,3,0}, p11 = {4,3,0};
	qaws_vec3 m_bot = {2,0,(qaws_scalar)1.5}, m_top = {2,3,1};
	qaws_vec3 m_left = {0,(qaws_scalar)1.5,(qaws_scalar)0.5};
	qaws_vec3 m_right = {4,(qaws_scalar)1.5,(qaws_scalar)0.5};

	printf("visual_obj_coons\n");
	svg_ensure_output_dir();

	c0 = make_quad_3d(p00, m_bot, p10);
	c1 = make_quad_3d(p01, m_top, p11);
	d0 = make_quad_3d(p00, m_left, p01);
	d1 = make_quad_3d(p10, m_right, p11);

	{
		qaws_surface_coons_desc desc;
		desc.c0 = c0; desc.c1 = c1; desc.d0 = d0; desc.d1 = d1;
		qaws_surface_create_coons(&desc, &surf);
	}
	if (!surf) goto cleanup_coons;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/coons_patch.obj",
		OBJ_OUTPUT_DIR "/coons_patch.mtl")) goto cleanup_coons;

	obj_material(&w, "coons", 0.2, 0.7, 0.9);
	obj_group(&w, "coons_surface");
	obj_use_material(&w, "coons");
	obj_surface_mesh(&w, surf, 32, 32);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/coons_patch.obj\n");

cleanup_coons:
	qaws_surface_destroy(surf);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ Extrusion surface                                      */
/* ------------------------------------------------------------------ */
static void visual_obj_extrusion(void)
{
	obj_writer w;
	qaws_curve* profile = NULL;
	qaws_surface* surf = NULL;

	printf("visual_obj_extrusion\n");
	svg_ensure_output_dir();

	/* S-curve profile extruded in Z */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,-1}, {4,1} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 3;
		d.control_points = pts;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &profile);
	}
	if (!profile) return;

	{
		qaws_surface_extrusion_desc desc;
		desc.profile = profile;
		desc.direction.x = 0; desc.direction.y = 0; desc.direction.z = 1;
		desc.length = (qaws_scalar)3.0;
		qaws_surface_create_extrusion(&desc, &surf);
	}
	if (!surf) { qaws_curve_destroy(profile); return; }

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/extrusion_surface.obj",
		OBJ_OUTPUT_DIR "/extrusion_surface.mtl"))
	{
		qaws_surface_destroy(surf); qaws_curve_destroy(profile); return;
	}

	obj_material(&w, "extr", 0.9, 0.5, 0.2);
	obj_group(&w, "extrusion");
	obj_use_material(&w, "extr");
	obj_surface_mesh(&w, surf, 32, 8);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/extrusion_surface.obj\n");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(profile);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ Revolution surface (vase)                              */
/* ------------------------------------------------------------------ */
static void visual_obj_revolution(void)
{
	obj_writer w;
	qaws_curve* profile = NULL;
	qaws_surface* surf = NULL;
	qaws_scalar pi = (qaws_scalar)(3.14159265358979323846);

	printf("visual_obj_revolution\n");
	svg_ensure_output_dir();

	/* Vase-shaped profile */
	{
		qaws_vec2 pts[] = {
			{(qaws_scalar)1.5, 0},
			{(qaws_scalar)0.5, 1},
			{(qaws_scalar)1.2, 2},
			{(qaws_scalar)0.8, 3}
		};
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 3;
		d.control_points = pts;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &profile);
	}
	if (!profile) return;

	{
		qaws_surface_revolution_desc desc;
		desc.profile = profile;
		desc.axis_origin.x = 0; desc.axis_origin.y = 0; desc.axis_origin.z = 0;
		desc.axis_direction.x = 0; desc.axis_direction.y = 0; desc.axis_direction.z = 1;
		desc.angle = 2 * pi;
		qaws_surface_create_revolution(&desc, &surf);
	}
	if (!surf) { qaws_curve_destroy(profile); return; }

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/revolution_vase.obj",
		OBJ_OUTPUT_DIR "/revolution_vase.mtl"))
	{
		qaws_surface_destroy(surf); qaws_curve_destroy(profile); return;
	}

	obj_material(&w, "vase", 0.8, 0.3, 0.6);
	obj_group(&w, "vase");
	obj_use_material(&w, "vase");
	obj_surface_mesh(&w, surf, 48, 24);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/revolution_vase.obj\n");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(profile);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_32_surface_modeling_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_coons_patch();
	test_extrusion_surface();
	test_revolution_surface();

	/* Visual output */
	visual_obj_coons();
	visual_obj_extrusion();
	visual_obj_revolution();

	printf("32_surface_modeling: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
