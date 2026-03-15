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
/*  Test: Loft with 2 sections                                         */
/* ------------------------------------------------------------------ */
static void test_loft_two_sections(void)
{
	qaws_curve* s0 = NULL;
	qaws_curve* s1 = NULL;
	qaws_curve const* sections[2];
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_loft_two_sections\n");

	/* Section 0: line from (0,0,0) to (4,0,0) */
	{
		qaws_vec3 a = {0, 0, 0}, b = {4, 0, 0};
		s0 = make_line_3d(a, b);
	}
	/* Section 1: line from (0,3,2) to (4,3,2) */
	{
		qaws_vec3 a = {0, 3, 2}, b = {4, 3, 2};
		s1 = make_line_3d(a, b);
	}
	TEST_ASSERT(s0 != NULL && s1 != NULL, "loft sections created");

	sections[0] = s0;
	sections[1] = s1;

	{
		qaws_surface_loft_desc desc;
		desc.sections = sections;
		desc.section_count = 2;
		desc.v_parameters = NULL; /* uniform */
		s = qaws_surface_create_loft(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "loft 2-section surface created");

	/* eval(0,0) ~ (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "loft2 (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "loft2 (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "loft2 (0,0) z");

	/* eval(0,1) ~ (0,3,2) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "loft2 (0,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "loft2 (0,1) y");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)2.0), "loft2 (0,1) z");

	/* eval(0.5, 0.5) ~ (2, 1.5, 1) midpoint */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)2.0), "loft2 mid x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.5), "loft2 mid y");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)1.0), "loft2 mid z");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(s0);
	qaws_curve_destroy(s1);
}

/* ------------------------------------------------------------------ */
/*  Test: Loft with 4 sections                                         */
/* ------------------------------------------------------------------ */
static void test_loft_four_sections(void)
{
	qaws_curve* s0 = NULL;
	qaws_curve* s1 = NULL;
	qaws_curve* s2 = NULL;
	qaws_curve* s3 = NULL;
	qaws_curve const* sections[4];
	qaws_scalar v_params[4];
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_loft_four_sections\n");

	/* Section 0 at z=0: line (0,0,0)-(4,0,0) */
	{
		qaws_vec3 a = {0, 0, 0}, b = {4, 0, 0};
		s0 = make_line_3d(a, b);
	}
	/* Section 1: quadratic with raised middle (0,1,0),(2,1,2),(4,1,0) */
	{
		qaws_vec3 a = {0, 1, 0}, b = {2, 1, 2}, c = {4, 1, 0};
		s1 = make_quad_3d(a, b, c);
	}
	/* Section 2: quadratic with raised middle (0,2,0),(2,2,1),(4,2,0) */
	{
		qaws_vec3 a = {0, 2, 0}, b = {2, 2, 1}, c = {4, 2, 0};
		s2 = make_quad_3d(a, b, c);
	}
	/* Section 3: line (0,3,0)-(4,3,0) */
	{
		qaws_vec3 a = {0, 3, 0}, b = {4, 3, 0};
		s3 = make_line_3d(a, b);
	}
	TEST_ASSERT(s0 != NULL && s1 != NULL && s2 != NULL && s3 != NULL,
		"loft4 sections created");

	sections[0] = s0;
	sections[1] = s1;
	sections[2] = s2;
	sections[3] = s3;

	v_params[0] = 0;
	v_params[1] = (qaws_scalar)0.33;
	v_params[2] = (qaws_scalar)0.67;
	v_params[3] = (qaws_scalar)1.0;

	{
		qaws_surface_loft_desc desc;
		desc.sections = sections;
		desc.section_count = 4;
		desc.v_parameters = v_params;
		s = qaws_surface_create_loft(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "loft 4-section surface created");

	/* Corner (0,0) ~ (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "loft4 (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "loft4 (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "loft4 (0,0) z");

	/* Corner (1,1) ~ (4,3,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "loft4 (1,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "loft4 (1,1) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "loft4 (1,1) z");

	/* Center should have a raised z due to middle sections */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(r.position.z > (qaws_scalar)0.1, "loft4 center raised in z");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(s0);
	qaws_curve_destroy(s1);
	qaws_curve_destroy(s2);
	qaws_curve_destroy(s3);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ lofted surface                                         */
/* ------------------------------------------------------------------ */
static void visual_obj_loft(void)
{
	obj_writer w;
	qaws_curve* s0 = NULL;
	qaws_curve* s1 = NULL;
	qaws_curve* s2 = NULL;
	qaws_curve* s3 = NULL;
	qaws_curve const* sections[4];
	qaws_scalar v_params[4];
	qaws_surface* surf = NULL;

	printf("visual_obj_loft\n");
	svg_ensure_output_dir();

	/* s0: line from (-2,0,0) to (2,0,0) */
	{
		qaws_vec3 a = {-2, 0, 0}, b = {2, 0, 0};
		s0 = make_line_3d(a, b);
	}
	/* s1: quadratic from (-1,0,1) through (0,1,1) to (1,0,1) */
	{
		qaws_vec3 a = {-1, 0, 1}, b = {0, 1, 1}, c = {1, 0, 1};
		s1 = make_quad_3d(a, b, c);
	}
	/* s2: quadratic from (-1.5,0,2) through (0,1.5,2) to (1.5,0,2) */
	{
		qaws_vec3 a = {(qaws_scalar)-1.5, 0, 2}, b = {0, (qaws_scalar)1.5, 2}, c = {(qaws_scalar)1.5, 0, 2};
		s2 = make_quad_3d(a, b, c);
	}
	/* s3: line from (-1,0,3) to (1,0,3) */
	{
		qaws_vec3 a = {-1, 0, 3}, b = {1, 0, 3};
		s3 = make_line_3d(a, b);
	}

	if (!s0 || !s1 || !s2 || !s3) goto cleanup_loft;

	sections[0] = s0;
	sections[1] = s1;
	sections[2] = s2;
	sections[3] = s3;

	v_params[0] = 0;
	v_params[1] = (qaws_scalar)0.33;
	v_params[2] = (qaws_scalar)0.67;
	v_params[3] = (qaws_scalar)1.0;

	{
		qaws_surface_loft_desc desc;
		desc.sections = sections;
		desc.section_count = 4;
		desc.v_parameters = v_params;
		qaws_surface_create_loft(&desc, &surf);
	}
	if (!surf) goto cleanup_loft;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/36_lofted_surface.obj",
		OBJ_OUTPUT_DIR "/36_lofted_surface.mtl")) goto cleanup_loft;

	obj_material(&w, "loft", 0.5, 0.8, 0.3);
	obj_group(&w, "lofted_surface");
	obj_use_material(&w, "loft");
	obj_surface_mesh(&w, surf, 32, 32);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/36_lofted_surface.obj\n");

cleanup_loft:
	qaws_surface_destroy(surf);
	qaws_curve_destroy(s0);
	qaws_curve_destroy(s1);
	qaws_curve_destroy(s2);
	qaws_curve_destroy(s3);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_36_surface_loft_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_loft_two_sections();
	test_loft_four_sections();

	/* Visual output */
	visual_obj_loft();

	printf("36_surface_loft: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
