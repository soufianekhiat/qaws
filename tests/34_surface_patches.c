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
/*  Test: Bilinear patch                                               */
/* ------------------------------------------------------------------ */
static void test_bilinear_patch(void)
{
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_bilinear_patch\n");

	/* Flat quad: corners at (0,0,0), (4,0,0), (0,3,0), (4,3,0) */
	{
		qaws_surface_bilinear_desc desc;
		desc.p00.x = 0; desc.p00.y = 0; desc.p00.z = 0;
		desc.p10.x = 4; desc.p10.y = 0; desc.p10.z = 0;
		desc.p01.x = 0; desc.p01.y = 3; desc.p01.z = 0;
		desc.p11.x = 4; desc.p11.y = 3; desc.p11.z = 0;
		s = qaws_surface_create_bilinear(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "bilinear surface created");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_BILINEAR,
		"bilinear kind correct");

	/* Corner (0,0) should be (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "bilinear (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "bilinear (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "bilinear (0,0) z");

	/* Corner (1,1) should be (4,3,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "bilinear (1,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "bilinear (1,1) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "bilinear (1,1) z");

	/* Center (0.5,0.5) of flat quad should be (2, 1.5, 0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)2.0), "bilinear center x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.5), "bilinear center y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "bilinear center z");

	/* Normal should point in z direction for flat patch */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.normal.z) > (qaws_scalar)0.9, "bilinear flat normal ~z");

	/* Null args */
	{
		qaws_surface* tmp = NULL;
		s = qaws_surface_create_bilinear(NULL, &tmp);
		TEST_ASSERT(s != QAWS_STATUS_OK, "bilinear null desc rejected");
	}

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: Biquadratic patch                                            */
/* ------------------------------------------------------------------ */
static void test_biquadratic_patch(void)
{
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_biquadratic_patch\n");

	/* Dome-like shape: 3x3 grid with center raised in z */
	{
		qaws_surface_biquadratic_desc desc;
		/* Row 0 (v=0): bottom */
		desc.control_points[0].x = 0; desc.control_points[0].y = 0; desc.control_points[0].z = 0;
		desc.control_points[1].x = 2; desc.control_points[1].y = 0; desc.control_points[1].z = 0;
		desc.control_points[2].x = 4; desc.control_points[2].y = 0; desc.control_points[2].z = 0;
		/* Row 1 (v=0.5): middle, center raised */
		desc.control_points[3].x = 0;    desc.control_points[3].y = (qaws_scalar)1.5; desc.control_points[3].z = 0;
		desc.control_points[4].x = 2;    desc.control_points[4].y = (qaws_scalar)1.5; desc.control_points[4].z = 1;
		desc.control_points[5].x = 4;    desc.control_points[5].y = (qaws_scalar)1.5; desc.control_points[5].z = 0;
		/* Row 2 (v=1): top */
		desc.control_points[6].x = 0; desc.control_points[6].y = 3; desc.control_points[6].z = 0;
		desc.control_points[7].x = 2; desc.control_points[7].y = 3; desc.control_points[7].z = 0;
		desc.control_points[8].x = 4; desc.control_points[8].y = 3; desc.control_points[8].z = 0;
		s = qaws_surface_create_biquadratic(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "biquadratic surface created");
	TEST_ASSERT(qaws_surface_get_kind(surf) == QAWS_SURFACE_KIND_BIQUADRATIC,
		"biquadratic kind correct");

	/* Corner (0,0) should be (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "biquad (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "biquad (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "biquad (0,0) z");

	/* Corner (1,0) should be (4,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "biquad (1,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "biquad (1,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "biquad (1,0) z");

	/* Corner (0,1) should be (0,3,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "biquad (0,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "biquad (0,1) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "biquad (0,1) z");

	/* Corner (1,1) should be (4,3,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "biquad (1,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "biquad (1,1) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "biquad (1,1) z");

	/* Center eval(0.5,0.5) should be raised in z (z > 0.5) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(r.position.z > (qaws_scalar)0.1, "biquad center raised in z");

	/* Null args */
	{
		qaws_surface* tmp = NULL;
		s = qaws_surface_create_biquadratic(NULL, &tmp);
		TEST_ASSERT(s != QAWS_STATUS_OK, "biquadratic null desc rejected");
	}

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ bilinear patch                                         */
/* ------------------------------------------------------------------ */
static void visual_obj_bilinear(void)
{
	obj_writer w;
	qaws_surface* surf = NULL;

	printf("visual_obj_bilinear\n");
	svg_ensure_output_dir();

	/* Twisted saddle bilinear patch — corners at very different z heights */
	{
		qaws_surface_bilinear_desc desc;
		desc.p00.x = 0; desc.p00.y = 0; desc.p00.z = 0;
		desc.p10.x = 4; desc.p10.y = 0; desc.p10.z = 3;
		desc.p01.x = 0; desc.p01.y = 3; desc.p01.z = 3;
		desc.p11.x = 4; desc.p11.y = 3; desc.p11.z = 0;
		qaws_surface_create_bilinear(&desc, &surf);
	}
	if (!surf) return;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/34_bilinear_patch.obj",
		OBJ_OUTPUT_DIR "/34_bilinear_patch.mtl"))
	{
		qaws_surface_destroy(surf); return;
	}

	obj_material(&w, "bilinear", 0.3, 0.8, 0.5);
	obj_group(&w, "bilinear_surface");
	obj_use_material(&w, "bilinear");
	obj_surface_mesh(&w, surf, 16, 16);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/34_bilinear_patch.obj\n");

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ biquadratic patch                                      */
/* ------------------------------------------------------------------ */
static void visual_obj_biquadratic(void)
{
	obj_writer w;
	qaws_surface* surf = NULL;

	printf("visual_obj_biquadratic\n");
	svg_ensure_output_dir();

	/* Tall dome biquadratic patch — center raised high, edges pulled down */
	{
		qaws_surface_biquadratic_desc desc;
		desc.control_points[0].x = 0; desc.control_points[0].y = 0; desc.control_points[0].z = -1;
		desc.control_points[1].x = 2; desc.control_points[1].y = 0; desc.control_points[1].z = 1;
		desc.control_points[2].x = 4; desc.control_points[2].y = 0; desc.control_points[2].z = -1;
		desc.control_points[3].x = 0;    desc.control_points[3].y = (qaws_scalar)1.5; desc.control_points[3].z = 1;
		desc.control_points[4].x = 2;    desc.control_points[4].y = (qaws_scalar)1.5; desc.control_points[4].z = 4;
		desc.control_points[5].x = 4;    desc.control_points[5].y = (qaws_scalar)1.5; desc.control_points[5].z = 1;
		desc.control_points[6].x = 0; desc.control_points[6].y = 3; desc.control_points[6].z = -1;
		desc.control_points[7].x = 2; desc.control_points[7].y = 3; desc.control_points[7].z = 1;
		desc.control_points[8].x = 4; desc.control_points[8].y = 3; desc.control_points[8].z = -1;
		qaws_surface_create_biquadratic(&desc, &surf);
	}
	if (!surf) return;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/34_biquadratic_patch.obj",
		OBJ_OUTPUT_DIR "/34_biquadratic_patch.mtl"))
	{
		qaws_surface_destroy(surf); return;
	}

	obj_material(&w, "biquad", 0.7, 0.4, 0.9);
	obj_group(&w, "biquadratic_surface");
	obj_use_material(&w, "biquad");
	obj_surface_mesh(&w, surf, 32, 32);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/34_biquadratic_patch.obj\n");

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_34_surface_patches_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_bilinear_patch();
	test_biquadratic_patch();

	/* Visual output */
	visual_obj_bilinear();
	visual_obj_biquadratic();

	printf("34_surface_patches: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
