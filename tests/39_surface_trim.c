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
/*  Helper: make a quadratic 2D Bezier curve (for trim loops)          */
/* ------------------------------------------------------------------ */
static qaws_curve* make_quad_2d(qaws_vec2 a, qaws_vec2 b, qaws_vec2 c_pt)
{
	qaws_curve* crv = NULL;
	qaws_scalar pts[6];
	qaws_bezier_desc d;
	pts[0] = a.x; pts[1] = a.y;
	pts[2] = b.x; pts[3] = b.y;
	pts[4] = c_pt.x; pts[5] = c_pt.y;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
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
/*  Helper: create 4 quadratic 2D arcs forming a circular trim loop    */
/*  centered at (cx,cy) with given radius.                             */
/*  arc_curves[4] is filled with the 4 arcs.                           */
/* ------------------------------------------------------------------ */
static void make_circle_trim_arcs(
	qaws_scalar cx, qaws_scalar cy, qaws_scalar radius,
	qaws_curve* arc_curves[4])
{
	qaws_vec2 a, b, c;

	/* Arc 1: (cx+r, cy) -> (cx+r, cy+r) -> (cx, cy+r) */
	a.x = cx + radius; a.y = cy;
	b.x = cx + radius; b.y = cy + radius;
	c.x = cx;          c.y = cy + radius;
	arc_curves[0] = make_quad_2d(a, b, c);

	/* Arc 2: (cx, cy+r) -> (cx-r, cy+r) -> (cx-r, cy) */
	a.x = cx;          a.y = cy + radius;
	b.x = cx - radius; b.y = cy + radius;
	c.x = cx - radius; c.y = cy;
	arc_curves[1] = make_quad_2d(a, b, c);

	/* Arc 3: (cx-r, cy) -> (cx-r, cy-r) -> (cx, cy-r) */
	a.x = cx - radius; a.y = cy;
	b.x = cx - radius; b.y = cy - radius;
	c.x = cx;          c.y = cy - radius;
	arc_curves[2] = make_quad_2d(a, b, c);

	/* Arc 4: (cx, cy-r) -> (cx+r, cy-r) -> (cx+r, cy) */
	a.x = cx;          a.y = cy - radius;
	b.x = cx + radius; b.y = cy - radius;
	c.x = cx + radius; c.y = cy;
	arc_curves[3] = make_quad_2d(a, b, c);
}

/* ------------------------------------------------------------------ */
/*  Test: trimmed surface with outer boundary (contains)               */
/* ------------------------------------------------------------------ */
static void test_surface_trim_contains(void)
{
	qaws_curve* c0 = NULL, *c1 = NULL, *d0 = NULL, *d1 = NULL;
	qaws_surface* base = NULL;
	qaws_surface* trimmed = NULL;
	qaws_curve* arcs[4];
	qaws_curve const* arc_ptrs[4];
	qaws_trim_loop loop;
	qaws_surface_trim_desc tdesc;
	qaws_surface_eval_result r_base, r_trim;
	qaws_status s;
	int inside;
	unsigned int i;

	printf("test_surface_trim_contains\n");

	base = make_flat_coons(&c0, &c1, &d0, &d1);
	TEST_ASSERT(base != NULL, "base coons created for trim");

	/* Create circular trim loop centered at (0.5, 0.5) radius 0.3 */
	make_circle_trim_arcs((qaws_scalar)0.5, (qaws_scalar)0.5, (qaws_scalar)0.3, arcs);
	for (i = 0; i < 4; i++)
		arc_ptrs[i] = arcs[i];

	loop.curves = arc_ptrs;
	loop.curve_count = 4;
	loop.is_outer = 1;

	tdesc.base = base;
	tdesc.loops = &loop;
	tdesc.loop_count = 1;
	s = qaws_surface_create_trimmed(&tdesc, &trimmed);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(trimmed != NULL, "trimmed surface created");
	TEST_ASSERT(qaws_surface_get_kind(trimmed) == QAWS_SURFACE_KIND_TRIMMED,
		"trimmed kind correct");

	/* Center (0.5, 0.5) should be inside the trim */
	inside = -1;
	s = qaws_surface_trim_contains(trimmed,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &inside);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(inside == 1, "trim center is inside");

	/* Corner (0.0, 0.0) should be outside the trim */
	inside = -1;
	s = qaws_surface_trim_contains(trimmed, 0, 0, &inside);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(inside == 0, "trim corner is outside");

	/* Eval at center should match base surface eval */
	memset(&r_base, 0, sizeof(r_base));
	memset(&r_trim, 0, sizeof(r_trim));
	qaws_surface_evaluate(base, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r_base);
	qaws_surface_evaluate(trimmed, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r_trim);
	TEST_ASSERT(approx_eq_loose(r_base.position.x, r_trim.position.x),
		"trim eval matches base x");
	TEST_ASSERT(approx_eq_loose(r_base.position.y, r_trim.position.y),
		"trim eval matches base y");
	TEST_ASSERT(approx_eq_loose(r_base.position.z, r_trim.position.z),
		"trim eval matches base z");

	qaws_surface_destroy(trimmed);
	qaws_surface_destroy(base);
	for (i = 0; i < 4; i++)
		qaws_curve_destroy(arcs[i]);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);
}

/* ------------------------------------------------------------------ */
/*  Test: trimmed surface with hole                                    */
/* ------------------------------------------------------------------ */
static void test_surface_trim_hole(void)
{
	qaws_curve* c0 = NULL, *c1 = NULL, *d0 = NULL, *d1 = NULL;
	qaws_surface* base = NULL;
	qaws_surface* trimmed = NULL;
	qaws_curve* arcs[4];
	qaws_curve const* arc_ptrs[4];
	qaws_trim_loop loop;
	qaws_surface_trim_desc tdesc;
	qaws_status s;
	int inside;
	unsigned int i;

	printf("test_surface_trim_hole\n");

	base = make_flat_coons(&c0, &c1, &d0, &d1);
	TEST_ASSERT(base != NULL, "base coons created for hole");

	/* Same circle, but as a hole (is_outer = 0) */
	make_circle_trim_arcs((qaws_scalar)0.5, (qaws_scalar)0.5, (qaws_scalar)0.3, arcs);
	for (i = 0; i < 4; i++)
		arc_ptrs[i] = arcs[i];

	loop.curves = arc_ptrs;
	loop.curve_count = 4;
	loop.is_outer = 0;

	tdesc.base = base;
	tdesc.loops = &loop;
	tdesc.loop_count = 1;
	s = qaws_surface_create_trimmed(&tdesc, &trimmed);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(trimmed != NULL, "trimmed surface with hole created");

	/* Center (0.5, 0.5) is inside the hole, so OUTSIDE the trimmed region */
	inside = -1;
	s = qaws_surface_trim_contains(trimmed,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &inside);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(inside == 0, "hole center is outside trimmed region");

	/* Corner (0.1, 0.1) is outside the hole, so INSIDE the trimmed region */
	inside = -1;
	s = qaws_surface_trim_contains(trimmed,
		(qaws_scalar)0.1, (qaws_scalar)0.1, &inside);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(inside == 1, "hole corner is inside trimmed region");

	qaws_surface_destroy(trimmed);
	qaws_surface_destroy(base);
	for (i = 0; i < 4; i++)
		qaws_curve_destroy(arcs[i]);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ trimmed surface                                        */
/* ------------------------------------------------------------------ */
static void visual_obj_trimmed(void)
{
	obj_writer w;
	qaws_curve* c0 = NULL, *c1 = NULL, *d0 = NULL, *d1 = NULL;
	qaws_surface* base = NULL;
	qaws_surface* trimmed = NULL;
	qaws_curve* arcs[4];
	qaws_curve const* arc_ptrs[4];
	qaws_trim_loop loop;
	qaws_surface_trim_desc tdesc;
	qaws_vec3 p00 = {0,0,0}, p10 = {4,0,0}, p01 = {0,3,0}, p11 = {4,3,0};
	qaws_vec3 m_bot = {2,0,(qaws_scalar)1.5};
	qaws_vec3 m_top = {2,3,1};
	qaws_vec3 m_left = {0,(qaws_scalar)1.5,(qaws_scalar)0.5};
	qaws_vec3 m_right = {4,(qaws_scalar)1.5,(qaws_scalar)0.5};
	unsigned int res = 48;
	unsigned int ui, vi, i;
	unsigned int first_v, first_n;

	printf("visual_obj_trimmed\n");
	svg_ensure_output_dir();

	/* Curved Coons base */
	c0 = make_quad_3d(p00, m_bot, p10);
	c1 = make_quad_3d(p01, m_top, p11);
	d0 = make_quad_3d(p00, m_left, p01);
	d1 = make_quad_3d(p10, m_right, p11);

	{
		qaws_surface_coons_desc desc;
		desc.c0 = c0; desc.c1 = c1; desc.d0 = d0; desc.d1 = d1;
		qaws_surface_create_coons(&desc, &base);
	}
	if (!base) goto cleanup_trim;

	/* Circular trim centered at (0.5, 0.5) radius 0.4 */
	make_circle_trim_arcs((qaws_scalar)0.5, (qaws_scalar)0.5, (qaws_scalar)0.4, arcs);
	for (i = 0; i < 4; i++)
		arc_ptrs[i] = arcs[i];

	loop.curves = arc_ptrs;
	loop.curve_count = 4;
	loop.is_outer = 1;

	tdesc.base = base;
	tdesc.loops = &loop;
	tdesc.loop_count = 1;
	qaws_surface_create_trimmed(&tdesc, &trimmed);
	if (!trimmed) goto cleanup_trim;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/trimmed_surface.obj",
		OBJ_OUTPUT_DIR "/trimmed_surface.mtl")) goto cleanup_trim;

	obj_material(&w, "trimmed", 0.3, 0.8, 0.6);
	obj_group(&w, "trimmed_surface");
	obj_use_material(&w, "trimmed");

	/* Custom trimmed mesh: emit vertices for the full grid, then only
	   emit faces where all 4 corners are inside the trim region. */
	first_v = w.vertex_count + 1;
	first_n = w.normal_count + 1;
	for (ui = 0; ui < res; ui++)
	{
		qaws_scalar u = (qaws_scalar)ui / (qaws_scalar)(res - 1);
		for (vi = 0; vi < res; vi++)
		{
			qaws_scalar v = (qaws_scalar)vi / (qaws_scalar)(res - 1);
			qaws_surface_eval_result r;
			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(base, u, v,
				QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
			obj_vertex(&w, r.position);
			obj_normal(&w, r.normal);
		}
	}

	for (ui = 0; ui < res - 1; ui++)
	{
		for (vi = 0; vi < res - 1; vi++)
		{
			qaws_scalar u0 = (qaws_scalar)ui / (qaws_scalar)(res - 1);
			qaws_scalar u1 = (qaws_scalar)(ui + 1) / (qaws_scalar)(res - 1);
			qaws_scalar v0 = (qaws_scalar)vi / (qaws_scalar)(res - 1);
			qaws_scalar v1 = (qaws_scalar)(vi + 1) / (qaws_scalar)(res - 1);
			int in00 = 0, in01 = 0, in10 = 0, in11 = 0;

			qaws_surface_trim_contains(trimmed, u0, v0, &in00);
			qaws_surface_trim_contains(trimmed, u0, v1, &in01);
			qaws_surface_trim_contains(trimmed, u1, v0, &in10);
			qaws_surface_trim_contains(trimmed, u1, v1, &in11);

			if (in00 && in01 && in10 && in11)
			{
				unsigned int v00_idx = first_v + ui * res + vi;
				unsigned int v01_idx = v00_idx + 1;
				unsigned int v10_idx = v00_idx + res;
				unsigned int v11_idx = v10_idx + 1;
				unsigned int n00_idx = first_n + ui * res + vi;
				unsigned int n01_idx = n00_idx + 1;
				unsigned int n10_idx = n00_idx + res;
				unsigned int n11_idx = n10_idx + 1;
				fprintf(w.fp, "f %u//%u %u//%u %u//%u\n",
					v00_idx, n00_idx,
					v10_idx, n10_idx,
					v01_idx, n01_idx);
				fprintf(w.fp, "f %u//%u %u//%u %u//%u\n",
					v10_idx, n10_idx,
					v11_idx, n11_idx,
					v01_idx, n01_idx);
			}
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/trimmed_surface.obj\n");

cleanup_trim:
	qaws_surface_destroy(trimmed);
	qaws_surface_destroy(base);
	for (i = 0; i < 4; i++)
		qaws_curve_destroy(arcs[i]);
	qaws_curve_destroy(c0);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(d0);
	qaws_curve_destroy(d1);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_39_surface_trim_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_surface_trim_contains();
	test_surface_trim_hole();

	/* Visual output */
	visual_obj_trimmed();

	printf("39_surface_trim: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
