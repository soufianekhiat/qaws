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
/*  Test: Gordon 2x2 (equivalent to Coons patch)                       */
/* ------------------------------------------------------------------ */
static void test_gordon_2x2(void)
{
	qaws_curve* uc0 = NULL; /* u-curve at v=0 (bottom) */
	qaws_curve* uc1 = NULL; /* u-curve at v=1 (top) */
	qaws_curve* vc0 = NULL; /* v-curve at u=0 (left) */
	qaws_curve* vc1 = NULL; /* v-curve at u=1 (right) */
	qaws_curve const* u_curves[2];
	qaws_curve const* v_curves[2];
	qaws_scalar v_params[2];
	qaws_scalar u_params[2];
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_gordon_2x2\n");

	/* Flat quad: (0,0,0), (4,0,0), (0,3,0), (4,3,0) */
	{
		qaws_vec3 p00 = {0, 0, 0}, p10 = {4, 0, 0};
		qaws_vec3 p01 = {0, 3, 0}, p11 = {4, 3, 0};

		uc0 = make_line_3d(p00, p10); /* bottom: v=0 */
		uc1 = make_line_3d(p01, p11); /* top: v=1 */
		vc0 = make_line_3d(p00, p01); /* left: u=0 */
		vc1 = make_line_3d(p10, p11); /* right: u=1 */
	}
	TEST_ASSERT(uc0 != NULL && uc1 != NULL && vc0 != NULL && vc1 != NULL,
		"gordon 2x2 curves created");

	u_curves[0] = uc0;
	u_curves[1] = uc1;
	v_curves[0] = vc0;
	v_curves[1] = vc1;
	v_params[0] = 0;
	v_params[1] = 1;
	u_params[0] = 0;
	u_params[1] = 1;

	{
		qaws_surface_gordon_desc desc;
		desc.u_curves = u_curves;
		desc.u_curve_count = 2;
		desc.v_params = v_params;
		desc.v_curves = v_curves;
		desc.v_curve_count = 2;
		desc.u_params = u_params;
		s = qaws_surface_create_gordon(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "gordon 2x2 surface created");

	/* Corner (0,0) ~ (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "gordon2 (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "gordon2 (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "gordon2 (0,0) z");

	/* Corner (1,1) ~ (4,3,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "gordon2 (1,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "gordon2 (1,1) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "gordon2 (1,1) z");

	/* Center (0.5,0.5) ~ (2,1.5,0) -- same as Coons on flat quad */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)2.0), "gordon2 center x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.5), "gordon2 center y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "gordon2 center z");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(uc0);
	qaws_curve_destroy(uc1);
	qaws_curve_destroy(vc0);
	qaws_curve_destroy(vc1);
}

/* ------------------------------------------------------------------ */
/*  Test: Gordon 3x3                                                   */
/* ------------------------------------------------------------------ */
static void test_gordon_3x3(void)
{
	qaws_curve* uc0 = NULL; /* u-curve at v=0 */
	qaws_curve* uc1 = NULL; /* u-curve at v=0.5 */
	qaws_curve* uc2 = NULL; /* u-curve at v=1 */
	qaws_curve* vc0 = NULL; /* v-curve at u=0 */
	qaws_curve* vc1 = NULL; /* v-curve at u=0.5 */
	qaws_curve* vc2 = NULL; /* v-curve at u=1 */
	qaws_curve const* u_curves[3];
	qaws_curve const* v_curves[3];
	qaws_scalar v_params[3];
	qaws_scalar u_params[3];
	qaws_surface* surf = NULL;
	qaws_surface_eval_result r;
	qaws_status s;

	printf("test_gordon_3x3\n");

	/*
	 * 3x3 grid with raised center point at (2, 1.5, 1):
	 *
	 *   (0,3,0) --- (2,3,0) --- (4,3,0)
	 *      |           |           |
	 *   (0,1.5,0) - (2,1.5,1) - (4,1.5,0)
	 *      |           |           |
	 *   (0,0,0) --- (2,0,0) --- (4,0,0)
	 */

	/* u-curves (horizontal, parameterized in u) */
	{
		qaws_vec3 a = {0, 0, 0}, b = {2, 0, 0}, c = {4, 0, 0};
		uc0 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {0, (qaws_scalar)1.5, 0}, b = {2, (qaws_scalar)1.5, 1}, c = {4, (qaws_scalar)1.5, 0};
		uc1 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {0, 3, 0}, b = {2, 3, 0}, c = {4, 3, 0};
		uc2 = make_quad_3d(a, b, c);
	}

	/* v-curves (vertical, parameterized in v) */
	{
		qaws_vec3 a = {0, 0, 0}, b = {0, (qaws_scalar)1.5, 0}, c = {0, 3, 0};
		vc0 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {2, 0, 0}, b = {2, (qaws_scalar)1.5, 1}, c = {2, 3, 0};
		vc1 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {4, 0, 0}, b = {4, (qaws_scalar)1.5, 0}, c = {4, 3, 0};
		vc2 = make_quad_3d(a, b, c);
	}

	TEST_ASSERT(uc0 && uc1 && uc2 && vc0 && vc1 && vc2,
		"gordon 3x3 curves created");

	u_curves[0] = uc0;
	u_curves[1] = uc1;
	u_curves[2] = uc2;
	v_curves[0] = vc0;
	v_curves[1] = vc1;
	v_curves[2] = vc2;

	v_params[0] = 0;
	v_params[1] = (qaws_scalar)0.5;
	v_params[2] = 1;
	u_params[0] = 0;
	u_params[1] = (qaws_scalar)0.5;
	u_params[2] = 1;

	{
		qaws_surface_gordon_desc desc;
		desc.u_curves = u_curves;
		desc.u_curve_count = 3;
		desc.v_params = v_params;
		desc.v_curves = v_curves;
		desc.v_curve_count = 3;
		desc.u_params = u_params;
		s = qaws_surface_create_gordon(&desc, &surf);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(surf != NULL, "gordon 3x3 surface created");

	/* Verify interpolation at grid intersections */
	/* Corner (0,0) ~ (0,0,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "gordon3 (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "gordon3 (0,0) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "gordon3 (0,0) z");

	/* Center (0.5,0.5) ~ (2, 1.5, 1) -- interpolation at grid center */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)2.0), "gordon3 center x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.5), "gordon3 center y");
	TEST_ASSERT(r.position.z > (qaws_scalar)0.3, "gordon3 center raised z");

	/* Corner (1,1) ~ (4,3,0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(surf, 1, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "gordon3 (1,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "gordon3 (1,1) y");
	TEST_ASSERT(fabs(r.position.z) < TOLERANCE_LOOSE, "gordon3 (1,1) z");

	qaws_surface_destroy(surf);
	qaws_curve_destroy(uc0);
	qaws_curve_destroy(uc1);
	qaws_curve_destroy(uc2);
	qaws_curve_destroy(vc0);
	qaws_curve_destroy(vc1);
	qaws_curve_destroy(vc2);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ Gordon surface                                         */
/* ------------------------------------------------------------------ */
static void visual_obj_gordon(void)
{
	obj_writer w;
	qaws_curve* uc0 = NULL;
	qaws_curve* uc1 = NULL;
	qaws_curve* uc2 = NULL;
	qaws_curve* vc0 = NULL;
	qaws_curve* vc1 = NULL;
	qaws_curve* vc2 = NULL;
	qaws_curve const* u_curves[3];
	qaws_curve const* v_curves[3];
	qaws_scalar v_params[3];
	qaws_scalar u_params[3];
	qaws_surface* surf = NULL;

	printf("visual_obj_gordon\n");
	svg_ensure_output_dir();

	/* u-curves with curved boundaries */
	{
		qaws_vec3 a = {0, 0, 0}, b = {2, 0, (qaws_scalar)0.5}, c = {4, 0, 0};
		uc0 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {0, (qaws_scalar)1.5, 0}, b = {2, (qaws_scalar)1.5, (qaws_scalar)1.5}, c = {4, (qaws_scalar)1.5, 0};
		uc1 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {0, 3, 0}, b = {2, 3, (qaws_scalar)0.5}, c = {4, 3, 0};
		uc2 = make_quad_3d(a, b, c);
	}

	/* v-curves with curved boundaries */
	{
		qaws_vec3 a = {0, 0, 0}, b = {0, (qaws_scalar)1.5, 0}, c = {0, 3, 0};
		vc0 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {2, 0, (qaws_scalar)0.5}, b = {2, (qaws_scalar)1.5, (qaws_scalar)1.5}, c = {2, 3, (qaws_scalar)0.5};
		vc1 = make_quad_3d(a, b, c);
	}
	{
		qaws_vec3 a = {4, 0, 0}, b = {4, (qaws_scalar)1.5, 0}, c = {4, 3, 0};
		vc2 = make_quad_3d(a, b, c);
	}

	if (!uc0 || !uc1 || !uc2 || !vc0 || !vc1 || !vc2) goto cleanup_gordon;

	u_curves[0] = uc0;
	u_curves[1] = uc1;
	u_curves[2] = uc2;
	v_curves[0] = vc0;
	v_curves[1] = vc1;
	v_curves[2] = vc2;

	v_params[0] = 0;
	v_params[1] = (qaws_scalar)0.5;
	v_params[2] = 1;
	u_params[0] = 0;
	u_params[1] = (qaws_scalar)0.5;
	u_params[2] = 1;

	{
		qaws_surface_gordon_desc desc;
		desc.u_curves = u_curves;
		desc.u_curve_count = 3;
		desc.v_params = v_params;
		desc.v_curves = v_curves;
		desc.v_curve_count = 3;
		desc.u_params = u_params;
		qaws_surface_create_gordon(&desc, &surf);
	}
	if (!surf) goto cleanup_gordon;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/37_gordon_surface.obj",
		OBJ_OUTPUT_DIR "/37_gordon_surface.mtl")) goto cleanup_gordon;

	obj_material(&w, "gordon", 0.8, 0.5, 0.3);
	obj_group(&w, "gordon_surface");
	obj_use_material(&w, "gordon");
	obj_surface_mesh(&w, surf, 32, 32);
	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/37_gordon_surface.obj\n");

cleanup_gordon:
	qaws_surface_destroy(surf);
	qaws_curve_destroy(uc0);
	qaws_curve_destroy(uc1);
	qaws_curve_destroy(uc2);
	qaws_curve_destroy(vc0);
	qaws_curve_destroy(vc1);
	qaws_curve_destroy(vc2);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_37_surface_gordon_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_gordon_2x2();
	test_gordon_3x3();

	/* Visual output */
	visual_obj_gordon();

	printf("37_surface_gordon: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
