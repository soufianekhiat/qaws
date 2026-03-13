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
/*  Helper: create a flat Coons patch (rectangle in XY plane)          */
/* ------------------------------------------------------------------ */
static qaws_surface* make_flat_coons(
	qaws_curve** out_c0, qaws_curve** out_c1,
	qaws_curve** out_d0, qaws_curve** out_d1)
{
	qaws_surface* surf = NULL;
	qaws_vec3 p00 = {0,0,0}, p10 = {4,0,0}, p01 = {0,3,0}, p11 = {4,3,0};
	qaws_surface_coons_desc desc;

	*out_c0 = make_line_3d(p00, p10);
	*out_c1 = make_line_3d(p01, p11);
	*out_d0 = make_line_3d(p00, p01);
	*out_d1 = make_line_3d(p10, p11);

	desc.c0 = *out_c0; desc.c1 = *out_c1;
	desc.d0 = *out_d0; desc.d1 = *out_d1;
	qaws_surface_create_coons(&desc, &surf);
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Test: surface offset (positive)                                    */
/* ------------------------------------------------------------------ */
static void test_surface_offset(void)
{
	qaws_curve* c0 = NULL, *c1 = NULL, *d0 = NULL, *d1 = NULL;
	qaws_surface* base = NULL;
	qaws_surface* offset = NULL;
	qaws_surface_eval_result r;
	qaws_status s;
	qaws_surface_offset_desc odesc;

	printf("test_surface_offset\n");

	base = make_flat_coons(&c0, &c1, &d0, &d1);
	TEST_ASSERT(base != NULL, "base coons created for offset");

	/* Create offset surface, distance +1.0 */
	odesc.base = base;
	odesc.distance = (qaws_scalar)1.0;
	s = qaws_surface_create_offset(&odesc, &offset);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(offset != NULL, "offset surface created");
	TEST_ASSERT(qaws_surface_get_kind(offset) == QAWS_SURFACE_KIND_OFFSET,
		"offset kind correct");

	/* Center (0.5, 0.5): flat patch normal is +Z, so offset pushes z up by 1.0 */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(offset, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)2.0), "offset center x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.5), "offset center y");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)1.0), "offset center z");

	/* Corner (0,0): should be offset by 1.0 in Z */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(offset, 0, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(r.position.x) < TOLERANCE_LOOSE, "offset (0,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "offset (0,0) y");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)1.0), "offset (0,0) z");

	/* Corner (1,1): should be (4, 3, 1.0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(offset, 1, 1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "offset (1,1) x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)3.0), "offset (1,1) y");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)1.0), "offset (1,1) z");

	/* Corner (1,0): should be (4, 0, 1.0) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(offset, 1, 0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)4.0), "offset (1,0) x");
	TEST_ASSERT(fabs(r.position.y) < TOLERANCE_LOOSE, "offset (1,0) y");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)1.0), "offset (1,0) z");

	/* Null args rejected */
	{
		qaws_surface* tmp = NULL;
		s = qaws_surface_create_offset(NULL, &tmp);
		TEST_ASSERT(s != QAWS_STATUS_OK, "offset null desc rejected");
	}

	qaws_surface_destroy(offset);
	qaws_surface_destroy(base);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);
}

/* ------------------------------------------------------------------ */
/*  Test: surface offset (negative)                                    */
/* ------------------------------------------------------------------ */
static void test_surface_offset_negative(void)
{
	qaws_curve* c0 = NULL, *c1 = NULL, *d0 = NULL, *d1 = NULL;
	qaws_surface* base = NULL;
	qaws_surface* offset = NULL;
	qaws_surface_eval_result r;
	qaws_status s;
	qaws_surface_offset_desc odesc;

	printf("test_surface_offset_negative\n");

	base = make_flat_coons(&c0, &c1, &d0, &d1);
	TEST_ASSERT(base != NULL, "base coons created for neg offset");

	/* Offset by -0.5 */
	odesc.base = base;
	odesc.distance = (qaws_scalar)-0.5;
	s = qaws_surface_create_offset(&odesc, &offset);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(offset != NULL, "negative offset surface created");

	/* Center (0.5, 0.5): should be (2, 1.5, -0.5) */
	memset(&r, 0, sizeof(r));
	s = qaws_surface_evaluate(offset, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(r.position.x, (qaws_scalar)2.0), "neg offset center x");
	TEST_ASSERT(approx_eq_loose(r.position.y, (qaws_scalar)1.5), "neg offset center y");
	TEST_ASSERT(approx_eq_loose(r.position.z, (qaws_scalar)-0.5), "neg offset center z");

	qaws_surface_destroy(offset);
	qaws_surface_destroy(base);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);
}

/* ------------------------------------------------------------------ */
/*  Helper: create a hemisphere surface of revolution                   */
/* ------------------------------------------------------------------ */
static qaws_surface* make_hemisphere(qaws_curve** out_profile)
{
	qaws_curve* prof = NULL;
	qaws_surface* surf = NULL;
	qaws_scalar pts[6];
	qaws_bezier_desc bd;
	qaws_surface_revolution_desc rd;

	/* Quadratic Bezier approximation of quarter circle (r=2):
	   (2, 0) -> (2, 2) -> (0, 2) in XY plane.
	   Profile x = radius, y = height along axis. */
	pts[0] = 2; pts[1] = 0;
	pts[2] = 2; pts[3] = 2;
	pts[4] = 0; pts[5] = 2;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D;
	bd.degree = 2;
	bd.control_points = pts;
	bd.control_point_count = 3;
	qaws_curve_create_bezier(&bd, &prof);
	if (!prof) return NULL;

	memset(&rd, 0, sizeof(rd));
	rd.profile = prof;
	rd.axis_origin.x = 0; rd.axis_origin.y = 0; rd.axis_origin.z = 0;
	rd.axis_direction.x = 0; rd.axis_direction.y = 0; rd.axis_direction.z = 1;
	rd.angle = 0; /* full 2*PI */
	qaws_surface_create_revolution(&rd, &surf);
	*out_profile = prof;
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ hemisphere with positive and negative offsets           */
/* ------------------------------------------------------------------ */
static void visual_obj_offset(void)
{
	obj_writer w;
	qaws_curve* prof = NULL;
	qaws_surface* base = NULL;
	qaws_surface* off_pos = NULL;
	qaws_surface* off_neg = NULL;
	qaws_surface_offset_desc odesc;

	printf("visual_obj_offset\n");
	svg_ensure_output_dir();

	base = make_hemisphere(&prof);
	if (!base) goto cleanup_offset;

	/* Positive offset (+0.5) — larger shell outside */
	odesc.base = base;
	odesc.distance = (qaws_scalar)0.5;
	qaws_surface_create_offset(&odesc, &off_pos);

	/* Negative offset (-0.5) — smaller shell inside */
	odesc.distance = (qaws_scalar)-0.5;
	qaws_surface_create_offset(&odesc, &off_neg);

	if (!off_pos || !off_neg) goto cleanup_offset;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/surface_offset.obj",
		OBJ_OUTPUT_DIR "/surface_offset.mtl")) goto cleanup_offset;

	obj_material(&w, "base", 0.5, 0.5, 0.8);
	obj_material(&w, "off_pos", 0.9, 0.3, 0.3);
	obj_material(&w, "off_neg", 0.3, 0.9, 0.3);

	obj_group(&w, "base_hemisphere");
	obj_use_material(&w, "base");
	obj_surface_mesh(&w, base, 32, 32);

	obj_group(&w, "positive_offset");
	obj_use_material(&w, "off_pos");
	obj_surface_mesh(&w, off_pos, 32, 32);

	obj_group(&w, "negative_offset");
	obj_use_material(&w, "off_neg");
	obj_surface_mesh(&w, off_neg, 32, 32);

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/surface_offset.obj\n");

cleanup_offset:
	qaws_surface_destroy(off_neg);
	qaws_surface_destroy(off_pos);
	qaws_surface_destroy(base);
	qaws_curve_destroy(prof);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ self-intersecting offset                                */
/* ------------------------------------------------------------------ */
static void visual_obj_offset_selfintersect(void)
{
	obj_writer w;
	qaws_curve* c0 = NULL, *c1 = NULL, *d0 = NULL, *d1 = NULL;
	qaws_surface* base = NULL;
	qaws_surface* offset = NULL;
	qaws_surface_offset_desc odesc;

	/* Sharp ridge Coons patch — high curvature along center ridge
	   so a large offset will self-intersect */
	qaws_vec3 p00 = {0,0,0}, p10 = {4,0,0}, p01 = {0,3,0}, p11 = {4,3,0};
	qaws_vec3 m_bot = {2,0,4};   /* tall ridge at center */
	qaws_vec3 m_top = {2,3,4};
	qaws_vec3 m_left = {0,(qaws_scalar)1.5,0};
	qaws_vec3 m_right = {4,(qaws_scalar)1.5,0};

	printf("visual_obj_offset_selfintersect\n");
	svg_ensure_output_dir();

	c0 = make_quad_3d(p00, m_bot, p10);
	c1 = make_quad_3d(p01, m_top, p11);
	d0 = make_quad_3d(p00, m_left, p01);
	d1 = make_quad_3d(p10, m_right, p11);

	{
		qaws_surface_coons_desc desc;
		desc.c0 = c0; desc.c1 = c1; desc.d0 = d0; desc.d1 = d1;
		qaws_surface_create_coons(&desc, &base);
	}
	if (!base) goto cleanup_si;

	/* Large offset exceeds radius of curvature -> self-intersection */
	odesc.base = base;
	odesc.distance = (qaws_scalar)2.0;
	qaws_surface_create_offset(&odesc, &offset);
	if (!offset) goto cleanup_si;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/offset_selfintersect.obj",
		OBJ_OUTPUT_DIR "/offset_selfintersect.mtl")) goto cleanup_si;

	obj_material(&w, "ridge", 0.5, 0.5, 0.8);
	obj_material(&w, "offset", 0.9, 0.3, 0.3);

	obj_group(&w, "base_ridge");
	obj_use_material(&w, "ridge");
	obj_surface_mesh(&w, base, 32, 32);

	obj_group(&w, "offset_surface");
	obj_use_material(&w, "offset");
	obj_surface_mesh(&w, offset, 48, 48);

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/offset_selfintersect.obj\n");

cleanup_si:
	qaws_surface_destroy(offset);
	qaws_surface_destroy(base);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_38_surface_offset_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_surface_offset();
	test_surface_offset_negative();

	/* Visual output */
	visual_obj_offset();
	visual_obj_offset_selfintersect();

	printf("38_surface_offset: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
