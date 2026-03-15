#include "test_common.h"
#include "qaws_platform.h"

/* ------------------------------------------------------------------ */
/*  Test: curve length matching                                        */
/* ------------------------------------------------------------------ */
static void test_curve_length_matching(void)
{
	qaws_curve* curve_a = NULL;
	qaws_curve* curve_b = NULL;
	qaws_curve* matched = NULL;
	qaws_scalar len_a, len_b;
	qaws_status s;

	printf("test_curve_length_matching\n");

	/* Curve A: short line from (0,0) to (1,0) */
	{
		qaws_scalar pts[4] = {0, 0, 1, 0};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 1;
		bd.control_points = pts;
		bd.control_point_count = 2;
		s = qaws_curve_create_bezier(&bd, &curve_a);
		TEST_ASSERT_STATUS(s);
	}

	/* Curve B: longer line from (0,0) to (3,0) -- arc length = 3 */
	{
		qaws_scalar pts[4] = {0, 0, 3, 0};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 1;
		bd.control_points = pts;
		bd.control_point_count = 2;
		s = qaws_curve_create_bezier(&bd, &curve_b);
		TEST_ASSERT_STATUS(s);
	}

	qaws_curve_compute_arc_length(curve_a, (qaws_scalar)0, (qaws_scalar)1, &len_a);
	qaws_curve_compute_arc_length(curve_b, (qaws_scalar)0, (qaws_scalar)1, &len_b);

	TEST_ASSERT(QAWS_FABS(len_a - (qaws_scalar)1.0) < (qaws_scalar)0.01, "curve A length ~ 1");
	TEST_ASSERT(QAWS_FABS(len_b - (qaws_scalar)3.0) < (qaws_scalar)0.01, "curve B length ~ 3");

	/* Match B to A's length */
	s = qaws_curve_match_arc_length(curve_a, curve_b, 256, &matched);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(matched != NULL, "length match produced a curve");

	/* The matched curve's domain should be [0, len_a] */
	{
		qaws_range range = qaws_curve_get_parameter_range(matched);
		TEST_ASSERT(QAWS_FABS(range.min_value) < (qaws_scalar)0.01,
			"matched curve starts at 0");
		TEST_ASSERT(QAWS_FABS(range.max_value - len_a) < (qaws_scalar)0.1,
			"matched curve ends at arc_length_a");
	}

	/* Evaluate matched curve at midpoint: should give halfway along B's geometry */
	{
		qaws_eval_result_2d er;
		qaws_range range = qaws_curve_get_parameter_range(matched);
		qaws_scalar mid = (range.min_value + range.max_value) * (qaws_scalar)0.5;
		memset(&er, 0, sizeof(er));
		qaws_curve_evaluate_2d(matched, mid, QAWS_EVAL_FLAG_POSITION, &er);
		/* At 50% of A's length, we should be at 50% of B's geometry = (1.5, 0) */
		TEST_ASSERT(QAWS_FABS(er.position.x - (qaws_scalar)1.5) < (qaws_scalar)0.2,
			"matched midpoint x ~ 1.5");
		TEST_ASSERT(QAWS_FABS(er.position.y) < (qaws_scalar)0.1,
			"matched midpoint y ~ 0");
	}

	qaws_curve_destroy(matched);
	qaws_curve_destroy(curve_a);
	qaws_curve_destroy(curve_b);
}

/* ------------------------------------------------------------------ */
/*  Test: 3D curve offsetting                                          */
/* ------------------------------------------------------------------ */
static void test_3d_offset(void)
{
	qaws_curve* curve = NULL;
	qaws_curve* offset_curve = NULL;
	qaws_status s;
	qaws_vec3 direction;

	printf("test_3d_offset\n");

	/* 3D line along X axis */
	{
		qaws_scalar pts[6] = {0, 0, 0, 2, 0, 0};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_3D;
		bd.degree = 1;
		bd.control_points = pts;
		bd.control_point_count = 2;
		s = qaws_curve_create_bezier(&bd, &curve);
		TEST_ASSERT_STATUS(s);
	}

	/* Offset in +Z direction by 0.5 */
	direction.x = 0; direction.y = 0; direction.z = 1;
	s = qaws_curve_offset_3d(curve, (qaws_scalar)0.5, 0, &direction, NULL, 64, &offset_curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(offset_curve != NULL, "3D offset produced a curve");

	/* Check that offset points are at z ~ 0.5 */
	{
		qaws_range range = qaws_curve_get_parameter_range(offset_curve);
		unsigned int i;
		for (i = 0; i <= 10; i++)
		{
			qaws_eval_result_3d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)10;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(offset_curve, t, QAWS_EVAL_FLAG_POSITION, &er);
			TEST_ASSERT(QAWS_FABS(er.position.z - (qaws_scalar)0.5) < (qaws_scalar)0.1,
				"offset z ~ 0.5");
			TEST_ASSERT(QAWS_FABS(er.position.y) < (qaws_scalar)0.1,
				"offset y ~ 0");
		}
	}

	qaws_curve_destroy(offset_curve);
	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/*  Test: approximate curve merging                                    */
/* ------------------------------------------------------------------ */
static void test_curve_merge(void)
{
	qaws_curve* seg1 = NULL;
	qaws_curve* seg2 = NULL;
	qaws_curve* seg3 = NULL;
	qaws_curve const* chain[3];
	qaws_curve* merged[4];
	unsigned int merged_count = 0;
	qaws_status s;

	printf("test_curve_merge\n");

	/* Three connected line segments forming a polyline */
	{
		qaws_scalar pts1[4] = {0, 0, 1, 1};
		qaws_scalar pts2[4] = {1, 1, 2, 0};
		qaws_scalar pts3[4] = {2, 0, 3, 1};
		qaws_bezier_desc bd;

		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 1;
		bd.control_point_count = 2;

		bd.control_points = pts1;
		s = qaws_curve_create_bezier(&bd, &seg1);
		TEST_ASSERT_STATUS(s);

		bd.control_points = pts2;
		s = qaws_curve_create_bezier(&bd, &seg2);
		TEST_ASSERT_STATUS(s);

		bd.control_points = pts3;
		s = qaws_curve_create_bezier(&bd, &seg3);
		TEST_ASSERT_STATUS(s);
	}

	chain[0] = seg1;
	chain[1] = seg2;
	chain[2] = seg3;

	memset(merged, 0, sizeof(merged));
	s = qaws_curve_merge_chain(chain, 3, 3, (qaws_scalar)0.5, merged, 4, &merged_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(merged_count >= 1, "merge produced at least 1 curve");
	TEST_ASSERT(merged[0] != NULL, "merged curve is non-null");

	/* Verify the merged curve passes near the junction points */
	{
		qaws_range range = qaws_curve_get_parameter_range(merged[0]);
		qaws_eval_result_2d er;
		qaws_scalar t_third = range.min_value +
			(range.max_value - range.min_value) / (qaws_scalar)3.0;

		memset(&er, 0, sizeof(er));
		qaws_curve_evaluate_2d(merged[0], t_third, QAWS_EVAL_FLAG_POSITION, &er);
		/* Near (1,1) junction */
		TEST_ASSERT(QAWS_FABS(er.position.x - (qaws_scalar)1.0) < (qaws_scalar)0.5,
			"merged 1/3 point x near junction");
	}

	{
		unsigned int mi;
		for (mi = 0; mi < merged_count; mi++)
			qaws_curve_destroy(merged[mi]);
	}
	qaws_curve_destroy(seg1);
	qaws_curve_destroy(seg2);
	qaws_curve_destroy(seg3);
}

/* ------------------------------------------------------------------ */
/*  Test: fillet on 2D composite curve                                 */
/* ------------------------------------------------------------------ */
static void test_fillet_2d(void)
{
	qaws_curve* seg1 = NULL;
	qaws_curve* seg2 = NULL;
	qaws_curve** segs = NULL;
	qaws_composite_desc cdesc;
	qaws_curve* composite = NULL;
	qaws_curve* filleted = NULL;
	qaws_status s;

	printf("test_fillet_2d\n");

	/* Two line segments meeting at a 90-degree corner at (1,0) */
	{
		qaws_scalar pts1[4] = {0, 0, 1, 0};
		qaws_scalar pts2[4] = {1, 0, 1, 1};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 1;
		bd.control_point_count = 2;

		bd.control_points = pts1;
		s = qaws_curve_create_bezier(&bd, &seg1);
		TEST_ASSERT_STATUS(s);

		bd.control_points = pts2;
		s = qaws_curve_create_bezier(&bd, &seg2);
		TEST_ASSERT_STATUS(s);
	}

	/* Create composite */
	segs = (qaws_curve**)malloc(sizeof(qaws_curve*) * 2);
	segs[0] = seg1;
	segs[1] = seg2;
	memset(&cdesc, 0, sizeof(cdesc));
	cdesc.dimension = QAWS_DIMENSION_2D;
	cdesc.segments = segs;
	cdesc.segment_count = 2;
	s = qaws_curve_create_composite(&cdesc, &composite);
	TEST_ASSERT_STATUS(s);
	free(segs);

	/* Apply fillet */
	s = qaws_curve_fillet_2d(composite, (qaws_scalar)0.2, &filleted);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(filleted != NULL, "fillet produced a curve");

	/* The filleted curve should have more spans than the original */
	TEST_ASSERT(qaws_curve_get_span_count(filleted) >= 2,
		"filleted curve has at least 2 spans");

	qaws_curve_destroy(filleted);
	qaws_curve_destroy(composite);
}

/* ------------------------------------------------------------------ */
/*  Test: chamfer on 2D composite curve                                */
/* ------------------------------------------------------------------ */
static void test_chamfer_2d(void)
{
	qaws_curve* seg1 = NULL;
	qaws_curve* seg2 = NULL;
	qaws_curve** segs = NULL;
	qaws_composite_desc cdesc;
	qaws_curve* composite = NULL;
	qaws_curve* chamfered = NULL;
	qaws_status s;

	printf("test_chamfer_2d\n");

	/* Two line segments with a corner */
	{
		qaws_scalar pts1[4] = {0, 0, 2, 0};
		qaws_scalar pts2[4] = {2, 0, 2, 2};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 1;
		bd.control_point_count = 2;

		bd.control_points = pts1;
		s = qaws_curve_create_bezier(&bd, &seg1);
		TEST_ASSERT_STATUS(s);

		bd.control_points = pts2;
		s = qaws_curve_create_bezier(&bd, &seg2);
		TEST_ASSERT_STATUS(s);
	}

	segs = (qaws_curve**)malloc(sizeof(qaws_curve*) * 2);
	segs[0] = seg1;
	segs[1] = seg2;
	memset(&cdesc, 0, sizeof(cdesc));
	cdesc.dimension = QAWS_DIMENSION_2D;
	cdesc.segments = segs;
	cdesc.segment_count = 2;
	s = qaws_curve_create_composite(&cdesc, &composite);
	TEST_ASSERT_STATUS(s);
	free(segs);

	s = qaws_curve_chamfer_2d(composite, (qaws_scalar)0.3, &chamfered);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(chamfered != NULL, "chamfer produced a curve");
	TEST_ASSERT(qaws_curve_get_span_count(chamfered) >= 2,
		"chamfered curve has at least 2 spans");

	qaws_curve_destroy(chamfered);
	qaws_curve_destroy(composite);
}

/* ------------------------------------------------------------------ */
/*  Visual: 3D curve offset                                            */
/* ------------------------------------------------------------------ */
static void visual_obj_3d_offset(void)
{
	obj_writer w;
	qaws_curve* helix = NULL;
	qaws_curve* offset_curve = NULL;
	qaws_status s;
	qaws_vec3 direction;
	unsigned int i;

	printf("visual_obj_3d_offset\n");
	svg_ensure_output_dir();

	/* Create a 3D helix curve via B-spline fitting */
	{
		unsigned int n = 64;
		qaws_scalar* coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n * 3);
		qaws_scalar* params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n);
		qaws_bspline_fit_desc fd;
		unsigned int hi;

		for (hi = 0; hi < n; hi++)
		{
			qaws_scalar t = (qaws_scalar)hi / (qaws_scalar)(n - 1);
			qaws_scalar angle = t * (qaws_scalar)6.283185;
			params[hi] = t;
			coords[hi * 3 + 0] = QAWS_COS(angle);
			coords[hi * 3 + 1] = QAWS_SIN(angle);
			coords[hi * 3 + 2] = t * (qaws_scalar)2.0;
		}
		memset(&fd, 0, sizeof(fd));
		fd.dimension = QAWS_DIMENSION_3D;
		fd.data_points = coords;
		fd.data_point_count = n;
		fd.degree = 3;
		fd.control_point_count = 20;
		fd.parameters = params;
		qaws_curve_fit_bspline(&fd, &helix);
		free(params); free(coords);
	}
	if (!helix) return;

	/* Offset in Z direction */
	direction.x = 0; direction.y = 0; direction.z = 1;
	s = qaws_curve_offset_3d(helix, (qaws_scalar)0.3, 0, &direction, NULL, 128, &offset_curve);
	if (s != QAWS_STATUS_OK || !offset_curve)
	{
		qaws_curve_destroy(helix);
		return;
	}

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/42_curve_offset_3d.obj",
		OBJ_OUTPUT_DIR "/42_curve_offset_3d.mtl"))
	{
		qaws_curve_destroy(helix); qaws_curve_destroy(offset_curve);
		return;
	}

	obj_material(&w, "original", 0.3, 0.6, 0.9);
	obj_material(&w, "offset", 1.0, 0.3, 0.0);

	obj_group(&w, "original");
	obj_use_material(&w, "original");
	{
		qaws_range range = qaws_curve_get_parameter_range(helix);
		for (i = 0; i <= 100; i++)
		{
			qaws_eval_result_3d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)100;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(helix, t, QAWS_EVAL_FLAG_POSITION, &er);
			obj_sphere(&w, er.position, (qaws_scalar)0.02);
		}
	}

	obj_group(&w, "offset_curve");
	obj_use_material(&w, "offset");
	{
		qaws_range range = qaws_curve_get_parameter_range(offset_curve);
		for (i = 0; i <= 100; i++)
		{
			qaws_eval_result_3d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)100;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(offset_curve, t, QAWS_EVAL_FLAG_POSITION, &er);
			obj_sphere(&w, er.position, (qaws_scalar)0.02);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/42_curve_offset_3d.obj\n");

	qaws_curve_destroy(offset_curve);
	qaws_curve_destroy(helix);
}

/* ------------------------------------------------------------------ */
/*  Visual: fillet and chamfer                                         */
/* ------------------------------------------------------------------ */
static void visual_svg_fillet_chamfer(void)
{
	svg_writer svg;
	qaws_curve* seg1 = NULL;
	qaws_curve* seg2 = NULL;
	qaws_curve* seg3 = NULL;
	qaws_curve** segs = NULL;
	qaws_composite_desc cdesc;
	qaws_curve* composite = NULL;
	qaws_curve* filleted = NULL;
	qaws_curve* chamfered = NULL;
	qaws_status s;

	printf("visual_svg_fillet_chamfer\n");
	svg_ensure_output_dir();

	/* Three-segment L-shape */
	{
		qaws_scalar pts1[4] = {0, 0, 2, 0};
		qaws_scalar pts2[4] = {2, 0, 2, 2};
		qaws_scalar pts3[4] = {2, 2, 4, 2};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 1;
		bd.control_point_count = 2;

		bd.control_points = pts1;
		qaws_curve_create_bezier(&bd, &seg1);
		bd.control_points = pts2;
		qaws_curve_create_bezier(&bd, &seg2);
		bd.control_points = pts3;
		qaws_curve_create_bezier(&bd, &seg3);
	}

	segs = (qaws_curve**)malloc(sizeof(qaws_curve*) * 3);
	segs[0] = seg1; segs[1] = seg2; segs[2] = seg3;
	memset(&cdesc, 0, sizeof(cdesc));
	cdesc.dimension = QAWS_DIMENSION_2D;
	cdesc.segments = segs;
	cdesc.segment_count = 3;
	s = qaws_curve_create_composite(&cdesc, &composite);
	free(segs);
	if (s != QAWS_STATUS_OK) return;

	qaws_curve_fillet_2d(composite, (qaws_scalar)0.4, &filleted);
	qaws_curve_chamfer_2d(composite, (qaws_scalar)0.4, &chamfered);

	if (!svg_open(&svg, SVG_OUTPUT_DIR "/42_fillet_chamfer.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)3.5,
		(qaws_scalar)600, (qaws_scalar)420))
		goto cleanup_fc;

	/* Draw original (grey) */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(composite, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#666666", (qaws_scalar)3.0);
	}

	/* Draw filleted (red) */
	if (filleted)
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(filleted, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.0);
	}

	/* Draw chamfered (blue) */
	if (chamfered)
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(chamfered, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)2.0);
	}

	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)3.0, "fillet (red), chamfer (blue)", "#aaaaaa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/42_fillet_chamfer.svg\n");

cleanup_fc:
	qaws_curve_destroy(filleted);
	qaws_curve_destroy(chamfered);
	qaws_curve_destroy(composite);
}

/* ------------------------------------------------------------------ */
/*  Visual: curve length matching                                      */
/* ------------------------------------------------------------------ */
static void visual_svg_length_match(void)
{
	svg_writer svg;
	qaws_curve* curve_a = NULL;
	qaws_curve* curve_b = NULL;
	qaws_curve* matched = NULL;
	qaws_status s;
	unsigned int i;

	printf("visual_svg_length_match\n");
	svg_ensure_output_dir();

	/* Curve A: short S-curve */
	{
		qaws_scalar pts[8] = {0, 0, (qaws_scalar)0.5, (qaws_scalar)1.0,
			(qaws_scalar)1.5, (qaws_scalar)-1.0, 2, 0};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 3;
		bd.control_points = pts;
		bd.control_point_count = 4;
		qaws_curve_create_bezier(&bd, &curve_a);
	}

	/* Curve B: longer S-curve */
	{
		qaws_scalar pts[8] = {0, (qaws_scalar)-2.5, (qaws_scalar)2.0, (qaws_scalar)-1.5,
			(qaws_scalar)4.0, (qaws_scalar)-3.5, (qaws_scalar)6.0, (qaws_scalar)-2.5};
		qaws_bezier_desc bd;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 3;
		bd.control_points = pts;
		bd.control_point_count = 4;
		qaws_curve_create_bezier(&bd, &curve_b);
	}

	if (!curve_a || !curve_b) goto cleanup_lm;

	s = qaws_curve_match_arc_length(curve_a, curve_b, 256, &matched);
	if (s != QAWS_STATUS_OK) goto cleanup_lm;

	if (!svg_open(&svg, SVG_OUTPUT_DIR "/42_length_match.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-4.0, (qaws_scalar)7.0, (qaws_scalar)6.0,
		(qaws_scalar)600, (qaws_scalar)500))
		goto cleanup_lm;

	/* Draw curve A (green) */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve_a, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)2.5);
	}
	svg_label(&svg, (qaws_scalar)0.1, (qaws_scalar)1.2, "A (short)", "#50fa7b");

	/* Draw curve B original (grey dashed) */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve_b, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#666666", (qaws_scalar)1.5);
	}
	svg_label(&svg, (qaws_scalar)0.1, (qaws_scalar)-1.8, "B (original)", "#666666");

	/* Draw matched curve (red) */
	if (matched)
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(matched, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.0);
		svg_label(&svg, (qaws_scalar)0.1, (qaws_scalar)-2.3, "B matched", "#e94560");

		/* Draw synchronized markers */
		{
			qaws_range range_a = qaws_curve_get_parameter_range(curve_a);
			qaws_range range_m = qaws_curve_get_parameter_range(matched);
			for (i = 0; i <= 10; i++)
			{
				qaws_eval_result_2d er_a, er_m;
				qaws_scalar frac = (qaws_scalar)i / (qaws_scalar)10;
				qaws_scalar t_a = range_a.min_value + frac * (range_a.max_value - range_a.min_value);
				qaws_scalar t_m = range_m.min_value + frac * (range_m.max_value - range_m.min_value);
				memset(&er_a, 0, sizeof(er_a));
				memset(&er_m, 0, sizeof(er_m));
				qaws_curve_evaluate_2d(curve_a, t_a, QAWS_EVAL_FLAG_POSITION, &er_a);
				qaws_curve_evaluate_2d(matched, t_m, QAWS_EVAL_FLAG_POSITION, &er_m);
				svg_line(&svg, er_a.position.x, er_a.position.y,
					er_m.position.x, er_m.position.y, "#e9456080", (qaws_scalar)1.0);
				{
					qaws_vec2 dot_a = er_a.position;
					qaws_vec2 dot_m = er_m.position;
					svg_dots(&svg, &dot_a, 1, "#50fa7b", (qaws_scalar)3.0);
					svg_dots(&svg, &dot_m, 1, "#e94560", (qaws_scalar)3.0);
				}
			}
		}
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/42_length_match.svg\n");

cleanup_lm:
	qaws_curve_destroy(matched);
	qaws_curve_destroy(curve_a);
	qaws_curve_destroy(curve_b);
}

/* ------------------------------------------------------------------ */
/*  Visual: curve merging                                              */
/* ------------------------------------------------------------------ */
static void visual_svg_merge(void)
{
	svg_writer svg;
	qaws_curve* seg1 = NULL;
	qaws_curve* seg2 = NULL;
	qaws_curve* seg3 = NULL;
	qaws_curve* seg4 = NULL;
	qaws_curve const* chain[4];
	qaws_curve* merged[4];
	unsigned int merged_count = 0;
	qaws_status s;
	unsigned int ci;

	printf("visual_svg_merge\n");
	svg_ensure_output_dir();

	/* Four connected cubic S-curves */
	{
		qaws_scalar pts1[8] = {0, 0, (qaws_scalar)0.2, (qaws_scalar)1.8, (qaws_scalar)0.8, (qaws_scalar)-0.3, 1, (qaws_scalar)1.5};
		qaws_scalar pts2[8] = {1, (qaws_scalar)1.5, (qaws_scalar)1.2, (qaws_scalar)2.5, (qaws_scalar)1.8, (qaws_scalar)-0.5, 2, (qaws_scalar)0.8};
		qaws_scalar pts3[8] = {2, (qaws_scalar)0.8, (qaws_scalar)2.15, (qaws_scalar)1.8, (qaws_scalar)2.85, (qaws_scalar)-0.2, 3, (qaws_scalar)1.2};
		qaws_scalar pts4[8] = {3, (qaws_scalar)1.2, (qaws_scalar)3.15, (qaws_scalar)2.2, (qaws_scalar)3.85, (qaws_scalar)0.5, 4, 0};
		qaws_bezier_desc bd;

		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 3;
		bd.control_point_count = 4;

		bd.control_points = pts1; qaws_curve_create_bezier(&bd, &seg1);
		bd.control_points = pts2; qaws_curve_create_bezier(&bd, &seg2);
		bd.control_points = pts3; qaws_curve_create_bezier(&bd, &seg3);
		bd.control_points = pts4; qaws_curve_create_bezier(&bd, &seg4);
	}

	if (!seg1 || !seg2 || !seg3 || !seg4) goto cleanup_merge;

	chain[0] = seg1; chain[1] = seg2; chain[2] = seg3; chain[3] = seg4;
	memset(merged, 0, sizeof(merged));
	s = qaws_curve_merge_chain(chain, 4, 3, (qaws_scalar)0.5, merged, 4, &merged_count);

	if (!svg_open(&svg, SVG_OUTPUT_DIR "/42_curve_merge.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.0, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
		goto cleanup_merge;

	/* Draw original segments */
	{
		char const* colors[4] = {"#888888", "#666666", "#888888", "#666666"};
		for (ci = 0; ci < 4; ci++)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(chain[ci], buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, colors[ci], (qaws_scalar)3.0);
		}
	}

	/* Draw junction dots */
	{
		qaws_vec2 junctions[3];
		junctions[0].x = 1; junctions[0].y = (qaws_scalar)1.5;
		junctions[1].x = 2; junctions[1].y = (qaws_scalar)0.8;
		junctions[2].x = 3; junctions[2].y = (qaws_scalar)1.2;
		svg_dots(&svg, junctions, 3, "#ffffff", (qaws_scalar)4.0);
	}

	/* Draw merged result */
	if (s == QAWS_STATUS_OK && merged_count > 0 && merged[0])
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(merged[0], buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.0);
	}

	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)2.0, "segments (grey), merged (red)", "#aaaaaa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/42_curve_merge.svg\n");

cleanup_merge:
	for (ci = 0; ci < merged_count; ci++)
		qaws_curve_destroy(merged[ci]);
	qaws_curve_destroy(seg1);
	qaws_curve_destroy(seg2);
	qaws_curve_destroy(seg3);
	qaws_curve_destroy(seg4);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_42_curve_operations_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_curve_length_matching();
	test_3d_offset();
	test_curve_merge();
	test_fillet_2d();
	test_chamfer_2d();

	/* Visual output */
	visual_obj_3d_offset();
	visual_svg_fillet_chamfer();
	visual_svg_length_match();
	visual_svg_merge();

	printf("42_curve_operations: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
