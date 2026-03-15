#include "test_common.h"

static void test_svg_inflection_points(void)
{
	printf("test_svg_inflection_points\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_inflection_points.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.0, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)500, (qaws_scalar)600))
		return;

	qaws_vec2 cp[] = { {0, 0}, {1, 4}, {2, -2}, {3, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Find and draw inflection points */
	qaws_scalar infl_params[8];
	unsigned int infl_count = 0;
	s = qaws_curve_find_inflection_points(curve, infl_params, 8, &infl_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(infl_count > 0, "inflection points found");

	for (unsigned int i = 0; i < infl_count; i++) {
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, infl_params[i],
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		svg_dot(&svg, r.position.x, r.position.y, "#e94560", (qaws_scalar)6);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.5, "Inflection Points", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/27_inflection_points.svg\n");
}

static void test_svg_extrema(void)
{
	printf("test_svg_extrema\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_extrema.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-5.0, (qaws_scalar)5.5, (qaws_scalar)12.0,
		(qaws_scalar)500, (qaws_scalar)800))
		return;

	/* --- Row 1: Wavy cubic Bezier --- */
	{
		qaws_vec2 cp[] = { {0,4}, {1,8}, {3,0}, {4,6} };
		qaws_bezier_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		if (qaws_curve_create_bezier(&desc, &crv) == QAWS_STATUS_OK) {
			svg_draw_extrema(&svg, crv, (qaws_scalar)3, (qaws_scalar)7,
				"#888888", "S-curve (cubic Bezier)");
			svg_polyline(&svg, cp, 4, "#ffffff20", (qaws_scalar)0.8);
			svg_dots(&svg, cp, 4, "#ffffff50", (qaws_scalar)2);
			qaws_curve_destroy(crv);
		}
	}

	/* --- Row 2: Tight loop Bezier --- */
	{
		qaws_vec2 cp[] = { {0,0}, {2,4}, {2,-2}, {4,2} };
		qaws_bezier_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		if (qaws_curve_create_bezier(&desc, &crv) == QAWS_STATUS_OK) {
			svg_draw_extrema(&svg, crv, (qaws_scalar)-2, (qaws_scalar)3,
				"#ab9df2", "Loop crossing (cubic Bezier)");
			svg_polyline(&svg, cp, 4, "#ffffff20", (qaws_scalar)0.8);
			svg_dots(&svg, cp, 4, "#ffffff50", (qaws_scalar)2);
			qaws_curve_destroy(crv);
		}
	}

	/* --- Row 3: Hermite 3-point wave --- */
	{
		qaws_scalar pts[] = { 0,-3,  2,-2,  4,-3 };
		qaws_scalar derivs[] = { 2,3,  0,-4,  2,3 };
		qaws_hermite_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = derivs;
		desc.point_count = 3;
		desc.derivative_count = 3;
		if (qaws_curve_create_hermite(&desc, &crv) == QAWS_STATUS_OK) {
			svg_draw_extrema(&svg, crv, (qaws_scalar)-5, (qaws_scalar)-1.5,
				"#f5a623", "W-shaped (Hermite)");
			qaws_curve_destroy(crv);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)7.2,
		"Extrema Detection", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)-4.5,
		"Blue dots=X extrema  Green dots=Y extrema", "#666688");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/27_extrema.svg\n");
}

static void test_svg_curvature_comb(void)
{
	printf("test_svg_curvature_comb\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_curvature_comb.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-2.0, (qaws_scalar)5.0, (qaws_scalar)5.0,
		(qaws_scalar)600, (qaws_scalar)500))
		return;

	/* Smooth S-shaped cubic Bezier — single inflection, clean curvature profile */
	qaws_vec2 cp[] = { {0, 0}, {1, 3}, {3, -1}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Compute curvature comb */
	{
		unsigned int comb_count = 48;
		qaws_curvature_sample_2d comb[48];
		s = qaws_curve_compute_curvature_comb_2d(curve, comb_count, comb, 48);
		TEST_ASSERT_STATUS(s);

		/* Find max |curvature| to normalize tooth length */
		qaws_scalar max_k = (qaws_scalar)0.001;
		for (unsigned int i = 0; i < comb_count; i++) {
			qaws_scalar ak = comb[i].curvature;
			if (ak < 0) ak = -ak;
			if (ak > max_k) max_k = ak;
		}

		/* Scale so the tallest tooth is about 0.8 world units */
		qaws_scalar scale = (qaws_scalar)0.8 / max_k;
		qaws_vec2 tips[48];

		for (unsigned int i = 0; i < comb_count; i++) {
			qaws_scalar tip_x = comb[i].position.x + comb[i].normal.x * comb[i].curvature * scale;
			qaws_scalar tip_y = comb[i].position.y + comb[i].normal.y * comb[i].curvature * scale;
			tips[i].x = tip_x;
			tips[i].y = tip_y;

			/* Draw comb tooth from curve to tip */
			svg_line(&svg, comb[i].position.x, comb[i].position.y,
				tip_x, tip_y, "#e9456088", (qaws_scalar)0.8);
		}

		/* Connect comb tips as an envelope */
		svg_polyline(&svg, tips, comb_count, "#e94560", (qaws_scalar)1.2);
	}

	/* Draw control polygon faintly */
	svg_polyline(&svg, cp, 4, "#ffffff30", (qaws_scalar)0.8);
	svg_dots(&svg, cp, 4, "#ffffff60", (qaws_scalar)2.5);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.8, "Curvature Comb (teeth = signed curvature)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/27_curvature_comb.svg\n");
}

static void test_svg_winding_number(void)
{
	printf("test_svg_winding_number\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_winding_number.svg",
		(qaws_scalar)-1.0, (qaws_scalar)-1.0, (qaws_scalar)12.0, (qaws_scalar)8.0,
		(qaws_scalar)750, (qaws_scalar)500))
		return;

	qaws_status s;

	/* --- Shape 1: Simple closed loop (CCW) --- */
	{
		qaws_vec2 pts[8];
		unsigned int i;
		qaws_scalar cx = (qaws_scalar)2.0, cy = (qaws_scalar)3.5, rad = (qaws_scalar)1.5;
		qaws_catmull_rom_desc desc;
		qaws_curve *crv = NULL;

		for (i = 0; i < 8; i++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979 * (double)i / 8.0);
			pts[i].x = cx + rad * (qaws_scalar)cos((double)angle);
			pts[i].y = cy + rad * (qaws_scalar)sin((double)angle);
		}
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = pts;
		desc.control_point_count = 8;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 1;

		if (qaws_curve_create_catmull_rom(&desc, &crv) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)2.5);

			/* Grid of test points */
			{
				int gx, gy;
				for (gx = 0; gx <= 8; gx++) {
					for (gy = 0; gy <= 8; gy++) {
						qaws_vec2 p;
						int w = 0;
						p.x = cx - rad - (qaws_scalar)0.3 + (qaws_scalar)gx * (rad * 2 + (qaws_scalar)0.6) / (qaws_scalar)8;
						p.y = cy - rad - (qaws_scalar)0.3 + (qaws_scalar)gy * (rad * 2 + (qaws_scalar)0.6) / (qaws_scalar)8;
						s = qaws_curve_compute_winding_number_2d(crv, p, &w);
						if (s == QAWS_STATUS_OK) {
							char const* col = (w != 0) ? "#50fa7b" : "#e94560";
							svg_dot(&svg, p.x, p.y, col, (qaws_scalar)2.5);
						}
					}
				}
			}
			svg_label(&svg, (qaws_scalar)0.2, (qaws_scalar)5.5, "Simple loop (w=1)", "#8888aa");
			qaws_curve_destroy(crv);
		}
	}

	/* --- Shape 2: Figure-eight (self-crossing, winding 0 in center) --- */
	{
		qaws_vec2 cp[] = { {6,4}, {8,7}, {9,1}, {7,4}, {5,7}, {6,1}, {8,4} };
		qaws_catmull_rom_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = cp;
		desc.control_point_count = 7;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 1;

		if (qaws_curve_create_catmull_rom(&desc, &crv) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#ab9df2", (qaws_scalar)2.5);

			/* Test points across the figure-eight */
			{
				qaws_vec2 test_pts[] = {
					{(qaws_scalar)6.5, (qaws_scalar)5.0},
					{(qaws_scalar)7.0, (qaws_scalar)4.0},
					{(qaws_scalar)7.5, (qaws_scalar)3.0},
					{(qaws_scalar)8.0, (qaws_scalar)2.0},
					{(qaws_scalar)5.5, (qaws_scalar)4.0},
					{(qaws_scalar)9.0, (qaws_scalar)4.0}
				};
				unsigned int ti;
				for (ti = 0; ti < 6; ti++) {
					int w = 0;
					s = qaws_curve_compute_winding_number_2d(crv, test_pts[ti], &w);
					if (s == QAWS_STATUS_OK) {
						char lb[16];
						char const* col = (w > 0) ? "#50fa7b" : (w < 0 ? "#f5a623" : "#e94560");
						sprintf(lb, "w=%d", w);
						svg_dot(&svg, test_pts[ti].x, test_pts[ti].y, col, (qaws_scalar)4);
						svg_label(&svg, test_pts[ti].x + (qaws_scalar)0.15,
							test_pts[ti].y + (qaws_scalar)0.15, lb, col);
					}
				}
			}
			svg_label(&svg, (qaws_scalar)5.5, (qaws_scalar)5.5, "Figure-eight", "#8888aa");
			qaws_curve_destroy(crv);
		}
	}

	/* --- Shape 3: Star shape (concave) --- */
	{
		qaws_vec2 pts[10];
		unsigned int i;
		qaws_scalar cx = (qaws_scalar)3.0, cy = (qaws_scalar)0.0;
		qaws_catmull_rom_desc desc;
		qaws_curve *crv = NULL;

		for (i = 0; i < 10; i++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979 * (double)i / 10.0);
			qaws_scalar r = (i % 2 == 0) ? (qaws_scalar)1.5 : (qaws_scalar)0.6;
			pts[i].x = cx + r * (qaws_scalar)cos((double)angle);
			pts[i].y = cy + r * (qaws_scalar)sin((double)angle);
		}
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = pts;
		desc.control_point_count = 10;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 1;

		if (qaws_curve_create_catmull_rom(&desc, &crv) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)2.5);

			/* Test center and tips */
			{
				qaws_vec2 test_pts[] = {
					{cx, cy},
					{cx + (qaws_scalar)1.0, cy},
					{cx, cy + (qaws_scalar)1.0},
					{cx - (qaws_scalar)2.0, cy}
				};
				unsigned int ti;
				for (ti = 0; ti < 4; ti++) {
					int w = 0;
					s = qaws_curve_compute_winding_number_2d(crv, test_pts[ti], &w);
					if (s == QAWS_STATUS_OK) {
						char lb[16];
						char const* col = (w != 0) ? "#50fa7b" : "#e94560";
						sprintf(lb, "w=%d", w);
						svg_dot(&svg, test_pts[ti].x, test_pts[ti].y, col, (qaws_scalar)4);
						svg_label(&svg, test_pts[ti].x + (qaws_scalar)0.15,
							test_pts[ti].y + (qaws_scalar)0.15, lb, col);
					}
				}
			}
			svg_label(&svg, (qaws_scalar)1.5, (qaws_scalar)-1.5, "Star (concave)", "#8888aa");
			qaws_curve_destroy(crv);
		}
	}

	/* --- Shape 4: Tiny ellipse-like for variety --- */
	{
		qaws_vec2 pts[6];
		unsigned int i;
		qaws_scalar cx = (qaws_scalar)8.5, cy = (qaws_scalar)0.0;
		qaws_catmull_rom_desc desc;
		qaws_curve *crv = NULL;

		for (i = 0; i < 6; i++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979 * (double)i / 6.0);
			pts[i].x = cx + (qaws_scalar)2.0 * (qaws_scalar)cos((double)angle);
			pts[i].y = cy + (qaws_scalar)0.8 * (qaws_scalar)sin((double)angle);
		}
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = pts;
		desc.control_point_count = 6;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 1;

		if (qaws_curve_create_catmull_rom(&desc, &crv) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#78dce8", (qaws_scalar)2.5);

			{
				qaws_vec2 test_pts[] = {
					{cx, cy},
					{cx + (qaws_scalar)1.5, cy},
					{cx, cy + (qaws_scalar)1.2}
				};
				unsigned int ti;
				for (ti = 0; ti < 3; ti++) {
					int w = 0;
					s = qaws_curve_compute_winding_number_2d(crv, test_pts[ti], &w);
					if (s == QAWS_STATUS_OK) {
						char lb[16];
						char const* col = (w != 0) ? "#50fa7b" : "#e94560";
						sprintf(lb, "w=%d", w);
						svg_dot(&svg, test_pts[ti].x, test_pts[ti].y, col, (qaws_scalar)4);
						svg_label(&svg, test_pts[ti].x + (qaws_scalar)0.15,
							test_pts[ti].y + (qaws_scalar)0.15, lb, col);
					}
				}
			}
			svg_label(&svg, (qaws_scalar)6.5, (qaws_scalar)-1.5, "Ellipse-like", "#8888aa");
			qaws_curve_destroy(crv);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.5, (qaws_scalar)7.0,
		"Winding Number", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.5, (qaws_scalar)6.5,
		"Green=inside (w!=0)  Red=outside (w=0)", "#666688");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/27_winding_number.svg\n");
}

static void test_svg_curvature_weighted(void)
{
	printf("test_svg_curvature_weighted\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_curvature_weighted.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.0, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)500, (qaws_scalar)600))
		return;

	qaws_vec2 cp[] = { {0, 0}, {(qaws_scalar)0.5, 4}, {(qaws_scalar)3.5, -2}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Uniform sampling: 20 dots in blue */
	{
		qaws_vec2 uniform_pts[20];
		qaws_sampling_desc sd;
		sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sd.sample_count = 20;
		sd.error_tolerance = 0;
		sd.use_adaptive_sampling = 0;
		unsigned int ucount = 0;
		s = qaws_curve_sample_2d(curve, &sd, uniform_pts, 20, &ucount);
		TEST_ASSERT_STATUS(s);
		svg_dots(&svg, uniform_pts, ucount, "#6272a4", (qaws_scalar)4);
	}

	/* Curvature-weighted sampling: 20 dots in red */
	{
		qaws_vec2 cw_pts[20];
		unsigned int cw_count = 0;
		s = qaws_curve_sample_curvature_weighted_2d(curve, 20, cw_pts, 20, &cw_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(cw_count > 0, "curvature-weighted samples produced");
		svg_dots(&svg, cw_pts, cw_count, "#e94560", (qaws_scalar)4);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.5,
		"Curvature-weighted (red) vs Uniform (blue)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/27_curvature_weighted.svg\n");
}

static void test_svg_feature_preserving(void)
{
	printf("test_svg_feature_preserving\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_feature_preserving.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-5.0, (qaws_scalar)5.5, (qaws_scalar)11.0,
		(qaws_scalar)500, (qaws_scalar)600))
		return;

	qaws_vec2 cp[] = { {0, 0}, {(qaws_scalar)0.5, 4}, {(qaws_scalar)3.5, -2}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_scalar offset_y = (qaws_scalar)4.0;

	/* --- Top row: Uniform sampling (blue) --- */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		for (unsigned int i = 0; i < n; i++) buf[i].y += offset_y;
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
	}
	{
		qaws_vec2 uniform_pts[8];
		qaws_sampling_desc sd;
		sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sd.sample_count = 8;
		sd.error_tolerance = 0;
		sd.use_adaptive_sampling = 0;
		unsigned int ucount = 0;
		s = qaws_curve_sample_2d(curve, &sd, uniform_pts, 8, &ucount);
		TEST_ASSERT_STATUS(s);
		for (unsigned int i = 0; i < ucount; i++) uniform_pts[i].y += offset_y;
		svg_dots(&svg, uniform_pts, ucount, "#6272a4", (qaws_scalar)4);
	}

	/* --- Bottom row: Feature-preserving sampling (green) --- */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
	}
	{
		qaws_vec2 fp_pts[32];
		unsigned int fp_count = 0;
		s = qaws_curve_sample_feature_preserving_2d(curve, 8, fp_pts, 32, &fp_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(fp_count >= 8, "feature-preserving has at least base_count samples");
		svg_dots(&svg, fp_pts, fp_count, "#50fa7b", (qaws_scalar)4);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)5.8,
		"Uniform 8 pts (blue)", "#6272a4");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)1.8,
		"Feature-preserving (green, adds inflections+extrema)", "#50fa7b");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/27_feature_preserving.svg\n");
}

static void test_svg_multi_traversal(void)
{
	printf("test_svg_multi_traversal\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_multi_traversal.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.5, (qaws_scalar)3.0,
		(qaws_scalar)600, (qaws_scalar)300))
		return;

	/* 3 connected cubic Bezier curves */
	qaws_vec2 cp1[] = { {0, 0}, {(qaws_scalar)0.5, 2}, {1, 0}, {(qaws_scalar)1.5, 0} };
	qaws_vec2 cp2[] = { {(qaws_scalar)1.5, 0}, {2, 0}, {(qaws_scalar)2.5, 2}, {3, 0} };
	qaws_vec2 cp3[] = { {3, 0}, {(qaws_scalar)3.5, 2}, {4, 0}, {(qaws_scalar)4.5, 0} };

	qaws_bezier_desc d1, d2, d3;
	d1.dimension = QAWS_DIMENSION_2D; d1.degree = 3;
	d1.control_points = cp1; d1.control_point_count = 4;
	d2.dimension = QAWS_DIMENSION_2D; d2.degree = 3;
	d2.control_points = cp2; d2.control_point_count = 4;
	d3.dimension = QAWS_DIMENSION_2D; d3.degree = 3;
	d3.control_points = cp3; d3.control_point_count = 4;

	qaws_curve* c1 = NULL;
	qaws_curve* c2 = NULL;
	qaws_curve* c3 = NULL;
	qaws_status s;
	s = qaws_curve_create_bezier(&d1, &c1); TEST_ASSERT_STATUS(s);
	s = qaws_curve_create_bezier(&d2, &c2); TEST_ASSERT_STATUS(s);
	s = qaws_curve_create_bezier(&d3, &c3); TEST_ASSERT_STATUS(s);

	/* Draw each segment in a different color */
	{
		qaws_vec2 seg[SVG_SAMPLES];
		unsigned int sn;
		sn = svg_sample_curve(c1, seg, SVG_SAMPLES);
		svg_polyline(&svg, seg, sn, "#e94560", (qaws_scalar)2.5);
		sn = svg_sample_curve(c2, seg, SVG_SAMPLES);
		svg_polyline(&svg, seg, sn, "#50fa7b", (qaws_scalar)2.5);
		sn = svg_sample_curve(c3, seg, SVG_SAMPLES);
		svg_polyline(&svg, seg, sn, "#6272a4", (qaws_scalar)2.5);
	}

	/* Create multi-curve traversal with constant speed + arc-length */
	qaws_curve const* curves[3];
	curves[0] = c1; curves[1] = c2; curves[2] = c3;

	qaws_traversal_desc tdesc;
	memset(&tdesc, 0, sizeof(tdesc));
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)1.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* trav = NULL;
	s = qaws_traversal_create_multi(curves, 3, &tdesc, &trav);
	TEST_ASSERT_STATUS(s);

	/* Compute total arc length */
	qaws_scalar len1, len2, len3, total_len;
	qaws_curve_compute_arc_length(c1, (qaws_scalar)0.0, (qaws_scalar)1.0, &len1);
	qaws_curve_compute_arc_length(c2, (qaws_scalar)0.0, (qaws_scalar)1.0, &len2);
	qaws_curve_compute_arc_length(c3, (qaws_scalar)0.0, (qaws_scalar)1.0, &len3);
	total_len = len1 + len2 + len3;

	/* Evaluate at 32 uniform distances and draw dots */
	for (int i = 0; i < 32; i++) {
		qaws_scalar dist = total_len * (qaws_scalar)i / (qaws_scalar)31;
		qaws_eval_result_2d r;
		s = qaws_traversal_evaluate_2d(trav, dist,
			QAWS_EVAL_FLAG_POSITION, &r);
		if (s == QAWS_STATUS_OK) {
			svg_dot(&svg, r.position.x, r.position.y, "#ffffff", (qaws_scalar)3);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.5,
		"Multi-curve traversal (uniform arc-length)", "#8888aa");
	svg_close(&svg);
	qaws_traversal_destroy(trav);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(c2);
	qaws_curve_destroy(c3);
	printf("  -> " SVG_OUTPUT_DIR "/27_multi_traversal.svg\n");
}

static void test_svg_frenet_frame(void)
{
	printf("test_svg_frenet_frame\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_frenet_frame.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.5, (qaws_scalar)5.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	qaws_vec2 cp[] = { {0, 0}, {1, 3}, {3, -1}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* At 8 evenly spaced parameters, draw tangent and normal vectors */
	qaws_scalar arrow_len = (qaws_scalar)0.3;
	for (int i = 0; i < 8; i++) {
		qaws_scalar t = (qaws_scalar)i / (qaws_scalar)7;
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) continue;

		qaws_vec2 tangent, normal;
		s = qaws_curve_compute_tangent_2d(curve, t, &tangent);
		if (s != QAWS_STATUS_OK) continue;
		s = qaws_curve_compute_normal_2d(curve, t, &normal);
		if (s != QAWS_STATUS_OK) continue;

		/* Draw tangent arrow (green) */
		svg_line(&svg, r.position.x, r.position.y,
			r.position.x + tangent.x * arrow_len,
			r.position.y + tangent.y * arrow_len,
			"#50fa7b", (qaws_scalar)1.5);

		/* Draw normal arrow (red) */
		svg_line(&svg, r.position.x, r.position.y,
			r.position.x + normal.x * arrow_len,
			r.position.y + normal.y * arrow_len,
			"#e94560", (qaws_scalar)1.5);

		/* Draw a small dot at the sample point */
		svg_dot(&svg, r.position.x, r.position.y, "#ffffff", (qaws_scalar)2.5);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5,
		"Frenet Frame (green=T, red=N)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/27_frenet_frame.svg\n");
}

static void test_svg_self_intersection(void)
{
	printf("test_svg_self_intersection\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_self_intersection.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.0, (qaws_scalar)3.0, (qaws_scalar)2.0,
		(qaws_scalar)600, (qaws_scalar)400))
		return;

	/* Figure-8 Hermite */
	qaws_vec2 pts[] = {
		{(qaws_scalar)0.0,  (qaws_scalar)0.0},
		{(qaws_scalar)1.0,  (qaws_scalar)0.5},
		{(qaws_scalar)0.0,  (qaws_scalar)0.0},
		{(qaws_scalar)-1.0, (qaws_scalar)-0.5},
		{(qaws_scalar)0.0,  (qaws_scalar)0.0}
	};
	qaws_vec2 tans[] = {
		{(qaws_scalar)1.0,  (qaws_scalar)1.0},
		{(qaws_scalar)0.0,  (qaws_scalar)-1.5},
		{(qaws_scalar)-1.0, (qaws_scalar)-1.0},
		{(qaws_scalar)0.0,  (qaws_scalar)1.5},
		{(qaws_scalar)1.0,  (qaws_scalar)1.0}
	};

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = pts;
	desc.derivatives = tans;
	desc.point_count = 5;
	desc.derivative_count = 5;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw curve */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Find self-intersections */
	qaws_intersection_2d isects[16];
	unsigned int isect_count = 0;
	s = qaws_curve_find_self_intersections_2d(curve, isects, 16, &isect_count);
	TEST_ASSERT_STATUS(s);

	/* Draw intersection points in red */
	for (unsigned int i = 0; i < isect_count && i < 16; ++i)
	{
		svg_dot(&svg, isects[i].position.x, isects[i].position.y,
			"#e94560", (qaws_scalar)6);
	}

	/* Draw key points in green */
	svg_dots(&svg, pts, 5, "#50fa7b", (qaws_scalar)3);

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)0.8,
		"Self-intersection (red dots)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/27_self_intersection.svg\n");
}

static void test_svg_curve_curve_intersection(void)
{
	printf("test_svg_curve_curve_intersection\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_curve_curve_intersection.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-5.0, (qaws_scalar)10.5, (qaws_scalar)11.0,
		(qaws_scalar)700, (qaws_scalar)700))
		return;

	qaws_status s;

	/* --- Example 1 (top-left): Arch vs S-curve (2 intersections) --- */
	{
		qaws_vec2 cpA[] = { {0,3}, {1,6}, {3,6}, {4,3} };
		qaws_vec2 cpB[] = { {0,5}, {2,2}, {2,6}, {4,4} };
		qaws_bezier_desc dA, dB;
		qaws_curve *cA = NULL, *cB = NULL;

		memset(&dA, 0, sizeof(dA));
		dA.dimension = QAWS_DIMENSION_2D; dA.degree = 3;
		dA.control_points = cpA; dA.control_point_count = 4;
		memset(&dB, 0, sizeof(dB));
		dB.dimension = QAWS_DIMENSION_2D; dB.degree = 3;
		dB.control_points = cpB; dB.control_point_count = 4;

		if (qaws_curve_create_bezier(&dA, &cA) == QAWS_STATUS_OK &&
			qaws_curve_create_bezier(&dB, &cB) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n;
			qaws_intersection_2d isects[16];
			unsigned int ic = 0;
			char lbl[64];

			n = svg_sample_curve(cA, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)2.5);
			n = svg_sample_curve(cB, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)2.5);

			s = qaws_curve_find_intersections_2d(cA, cB, isects, 16, &ic);
			if (s == QAWS_STATUS_OK) {
				unsigned int i;
				for (i = 0; i < ic; i++)
					svg_ring(&svg, isects[i].position.x, isects[i].position.y,
						"#e94560", (qaws_scalar)6, (qaws_scalar)2);
			}
			sprintf(lbl, "Arch vs S-curve (%u pts)", ic);
			svg_label(&svg, (qaws_scalar)0, (qaws_scalar)6.2, lbl, "#8888aa");
		}
		qaws_curve_destroy(cA); qaws_curve_destroy(cB);
	}

	/* --- Example 2 (top-right): Line vs loop (3 intersections) --- */
	{
		qaws_vec2 cpA[] = { {5,3}, {9,6} };
		qaws_vec2 cpB[] = { {5,4}, {7,7}, {7,2}, {9,5} };
		qaws_bezier_desc dA, dB;
		qaws_curve *cA = NULL, *cB = NULL;

		memset(&dA, 0, sizeof(dA));
		dA.dimension = QAWS_DIMENSION_2D; dA.degree = 1;
		dA.control_points = cpA; dA.control_point_count = 2;
		memset(&dB, 0, sizeof(dB));
		dB.dimension = QAWS_DIMENSION_2D; dB.degree = 3;
		dB.control_points = cpB; dB.control_point_count = 4;

		if (qaws_curve_create_bezier(&dA, &cA) == QAWS_STATUS_OK &&
			qaws_curve_create_bezier(&dB, &cB) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n;
			qaws_intersection_2d isects[16];
			unsigned int ic = 0;
			char lbl[64];

			n = svg_sample_curve(cA, buf, SVG_SAMPLES);
			svg_polyline_dashed(&svg, buf, n, "#f5a623", (qaws_scalar)1.5, "6,4");
			n = svg_sample_curve(cB, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#78dce8", (qaws_scalar)2.5);

			s = qaws_curve_find_intersections_2d(cA, cB, isects, 16, &ic);
			if (s == QAWS_STATUS_OK) {
				unsigned int i;
				for (i = 0; i < ic; i++)
					svg_ring(&svg, isects[i].position.x, isects[i].position.y,
						"#e94560", (qaws_scalar)6, (qaws_scalar)2);
			}
			sprintf(lbl, "Line vs loop (%u pts)", ic);
			svg_label(&svg, (qaws_scalar)5, (qaws_scalar)6.2, lbl, "#8888aa");
		}
		qaws_curve_destroy(cA); qaws_curve_destroy(cB);
	}

	/* --- Example 3 (bottom-left): Self-intersecting Hermite figure-8 --- */
	{
		qaws_scalar pts[] = { 2,-1,  3,-0.5,  2,-1,  1,-0.5,  2,-1 };
		qaws_scalar derivs[] = { 1,1,  0,-1.5,  -1,-1,  0,1.5,  1,1 };
		qaws_hermite_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = derivs;
		desc.point_count = 5;
		desc.derivative_count = 5;

		if (qaws_curve_create_hermite(&desc, &crv) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n;
			qaws_intersection_2d isects[16];
			unsigned int ic = 0;
			char lbl[64];

			n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#ab9df2", (qaws_scalar)2.5);

			s = qaws_curve_find_self_intersections_2d(crv, isects, 16, &ic);
			if (s == QAWS_STATUS_OK) {
				unsigned int i;
				for (i = 0; i < ic; i++) {
					svg_ring(&svg, isects[i].position.x, isects[i].position.y,
						"#e94560", (qaws_scalar)7, (qaws_scalar)2);
					svg_dot(&svg, isects[i].position.x, isects[i].position.y,
						"#e94560", (qaws_scalar)3);
				}
			}
			sprintf(lbl, "Self-intersection (%u pts)", ic);
			svg_label(&svg, (qaws_scalar)0, (qaws_scalar)-0.5, lbl, "#8888aa");
			qaws_curve_destroy(crv);
		}
	}

	/* --- Example 4 (bottom-right): Self-intersecting Hermite vs Catmull-Rom --- */
	{
		/* Hermite figure-8 centered at (7, -1.5) */
		qaws_scalar h_pts[] = { 7,-1.5,  8,-1,  7,-1.5,  6,-1,  7,-1.5 };
		qaws_scalar h_der[] = { 1,1,  0,-1.5,  -1,-1,  0,1.5,  1,1 };
		qaws_hermite_desc hd;
		qaws_curve *cA = NULL;

		/* Catmull-Rom wave through the figure-8 */
		qaws_vec2 cr_pts[] = { {5,-3}, {6,-1}, {7,-2.5}, {8,-1}, {9,-3} };
		qaws_catmull_rom_desc cd;
		qaws_curve *cB = NULL;

		memset(&hd, 0, sizeof(hd));
		hd.dimension = QAWS_DIMENSION_2D;
		hd.degree = 3;
		hd.points = h_pts;
		hd.derivatives = h_der;
		hd.point_count = 5;
		hd.derivative_count = 5;

		memset(&cd, 0, sizeof(cd));
		cd.dimension = QAWS_DIMENSION_2D;
		cd.control_points = cr_pts;
		cd.control_point_count = 5;
		cd.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cd.closed = 0;

		if (qaws_curve_create_hermite(&hd, &cA) == QAWS_STATUS_OK &&
			qaws_curve_create_catmull_rom(&cd, &cB) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n;
			qaws_intersection_2d isects[16];
			unsigned int ic = 0;
			char lbl[64];

			n = svg_sample_curve(cA, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#ff6188", (qaws_scalar)2.5);
			n = svg_sample_curve(cB, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)2.5);

			s = qaws_curve_find_intersections_2d(cA, cB, isects, 16, &ic);
			if (s == QAWS_STATUS_OK) {
				unsigned int i;
				for (i = 0; i < ic; i++)
					svg_ring(&svg, isects[i].position.x, isects[i].position.y,
						"#e94560", (qaws_scalar)6, (qaws_scalar)2);
			}
			sprintf(lbl, "Hermite fig-8 vs Catmull-Rom (%u pts)", ic);
			svg_label(&svg, (qaws_scalar)5, (qaws_scalar)-0.5, lbl, "#8888aa");
		}
		qaws_curve_destroy(cA); qaws_curve_destroy(cB);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)6.5,
		"Curve Intersections", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)-4.5,
		"Red rings = intersection points", "#666688");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/27_curve_curve_intersection.svg\n");
}

static void test_svg_offset(void)
{
	printf("test_svg_offset\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_curve_offset.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-7.0, (qaws_scalar)7.0, (qaws_scalar)11.5,
		(qaws_scalar)700, (qaws_scalar)1150))
		return;

	/* --- Row 1: Bezier S-curve with left/right offsets --- */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)3.0,
			(qaws_scalar)3.0, (qaws_scalar)-1.0,
			(qaws_scalar)4.0, (qaws_scalar)2.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;
		unsigned int k;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		if (s != QAWS_STATUS_OK) goto row1_end;

		/* Draw original curve (gray) */
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(bez, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

		/* Left offset (blue) */
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.4, 0, offs, 4, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)1.5);
				qaws_curve_destroy(offs[k]);
			}
		}

		/* Right offset (green) */
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)-0.4, 0, offs, 4, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)1.5);
				qaws_curve_destroy(offs[k]);
			}
		}

		svg_label(&svg, (qaws_scalar)-1.0, (qaws_scalar)3.5,
			"S-curve: left (blue) + right (green)", "#8888aa");
		qaws_curve_destroy(bez);
row1_end:;
	}

	/* --- Row 2: Tight S-curve with large offset (cusp region) --- */
	{
		/* Shifted down by -4 */
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)-2.5,
			(qaws_scalar)0.5, (qaws_scalar)-0.5,
			(qaws_scalar)1.5, (qaws_scalar)-4.5,
			(qaws_scalar)2.0, (qaws_scalar)-2.5
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;
		unsigned int k;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		if (s != QAWS_STATUS_OK) goto row2_end;

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(bez, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

		/* Large left offset WITHOUT trim (raw, dashed orange) */
		{
			qaws_curve* offs[8] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.8, 0, offs, 8, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#f0a030", (qaws_scalar)3.0);
				qaws_curve_destroy(offs[k]);
			}
		}

		/* Large left offset WITH trim (solid red) */
		{
			qaws_curve* offs[8] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.8, 1, offs, 8, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)1.0);
				qaws_curve_destroy(offs[k]);
			}
		}

		/* Smaller offset for comparison (no cusps expected) */
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.15, 0, offs, 4, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)1.5);
				qaws_curve_destroy(offs[k]);
			}
		}

		svg_label(&svg, (qaws_scalar)-1.0, (qaws_scalar)-0.3,
			"Tight S: offset=0.8 raw (orange) vs trimmed (red), small=0.15 (blue)", "#8888aa");
		qaws_curve_destroy(bez);
row2_end:;
	}

	/* --- Row 3: "Stroke" visualization — paired left+right offsets at multiple widths --- */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)-5.0,
			(qaws_scalar)1.5, (qaws_scalar)-3.5,
			(qaws_scalar)2.5, (qaws_scalar)-6.5,
			(qaws_scalar)4.0, (qaws_scalar)-5.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;
		unsigned int k;
		static const qaws_scalar widths[] = {
			(qaws_scalar)0.15, (qaws_scalar)0.3, (qaws_scalar)0.6
		};
		static const char* colors_l[] = { "#6272a4", "#bd93f9", "#ff79c6" };
		static const char* colors_r[] = { "#6272a4", "#bd93f9", "#ff79c6" };

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		if (s != QAWS_STATUS_OK) goto row3_end;

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(bez, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

		for (k = 0; k < 3; ++k)
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0, m;
			/* Left */
			s = qaws_curve_offset_2d(bez, widths[k], 0, offs, 4, &off_cnt);
			for (m = 0; m < off_cnt; ++m) {
				if (!offs[m]) continue;
				n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, colors_l[k], (qaws_scalar)1.0);
				qaws_curve_destroy(offs[m]);
			}
			/* Right */
			off_cnt = 0;
			memset(offs, 0, sizeof(offs));
			s = qaws_curve_offset_2d(bez, -widths[k], 0, offs, 4, &off_cnt);
			for (m = 0; m < off_cnt; ++m) {
				if (!offs[m]) continue;
				n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, colors_r[k], (qaws_scalar)1.0);
				qaws_curve_destroy(offs[m]);
			}
		}

		svg_label(&svg, (qaws_scalar)-1.0, (qaws_scalar)-2.8,
			"Stroke widths: 0.15 (blue) 0.3 (purple) 0.6 (pink)", "#8888aa");
		qaws_curve_destroy(bez);
row3_end:;
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/27_curve_offset.svg\n");
}

static void test_svg_offset_closed(void)
{
	printf("test_svg_offset_closed\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_offset_closed.svg",
		(qaws_scalar)-4.5, (qaws_scalar)-3.5, (qaws_scalar)8.0, (qaws_scalar)21.0,
		(qaws_scalar)600, (qaws_scalar)1600))
		return;

	/* Cayley oval: locus where harmonic mean of distances to two
	 * foci equals b.  1/r + 1/r' = 2/b, foci at (+-a, 0).
	 * b = a gives figure-8; b slightly > a gives an oval that
	 * is very close to self-intersecting (pinched waist). */
	{
		#define CAYLEY_N 64
		qaws_scalar pts[CAYLEY_N * 2];
		qaws_catmull_rom_desc cr_desc;
		qaws_curve* crv = NULL;
		qaws_status s;
		unsigned int k;
		double ca = 1.0, cb = 1.08; /* close to self-intersecting */

		for (k = 0; k < CAYLEY_N; ++k)
		{
			double theta = 2.0 * 3.14159265358979323846 * (double)k / (double)CAYLEY_N;
			double ct = cos(theta), st = sin(theta);
			/* Bisect for R: 1/r + 1/r' = 2/b where
			 * r = sqrt(R^2 + 2aR*cos(t) + a^2)
			 * r' = sqrt(R^2 - 2aR*cos(t) + a^2) */
			double lo = 0.001, hi = 4.0, mid, R;
			int iter;
			for (iter = 0; iter < 60; ++iter)
			{
				double r1, r2, val;
				mid = 0.5 * (lo + hi);
				r1 = sqrt(mid * mid + 2.0 * ca * mid * ct + ca * ca);
				r2 = sqrt(mid * mid - 2.0 * ca * mid * ct + ca * ca);
				val = 1.0 / r1 + 1.0 / r2;
				if (val > 2.0 / cb) lo = mid;
				else hi = mid;
			}
			R = 0.5 * (lo + hi);
			pts[k * 2 + 0] = (qaws_scalar)(R * ct);
			pts[k * 2 + 1] = (qaws_scalar)(R * st);
		}

		memset(&cr_desc, 0, sizeof(cr_desc));
		cr_desc.dimension = QAWS_DIMENSION_2D;
		cr_desc.control_points = pts;
		cr_desc.control_point_count = CAYLEY_N;
		cr_desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cr_desc.closed = 1;
		s = qaws_curve_create_catmull_rom(&cr_desc, &crv);
		if (s != QAWS_STATUS_OK) goto closed_end;

		/* Three rows stacked vertically: trim=0, trim=1, trim=2 */
		{
			static const qaws_scalar dists_out[] = {
				(qaws_scalar)0.3, (qaws_scalar)0.6, (qaws_scalar)1.0
			};
			static const qaws_scalar dists_in[] = {
				(qaws_scalar)-0.3, (qaws_scalar)-0.6, (qaws_scalar)-1.0
			};
			static const char* colors_out[] = { "#50fa7b", "#8be9fd", "#ff79c6" };
			static const char* colors_in[]  = { "#6272a4", "#bd93f9", "#e94560" };
			static const int trim_vals[] = { 0, 1, 2 };
			static const char* labels[] = {
				"trim=0 (raw)", "trim=1 (distance)", "trim=2 (cleanup)"
			};
			int row;
			for (row = 0; row < 3; ++row)
			{
				int do_trim = trim_vals[row];
				qaws_scalar yoff = (qaws_scalar)row * (qaws_scalar)7.0;

				/* Draw base curve shifted */
				{
					qaws_vec2 buf[SVG_SAMPLES];
					unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES), m;
					for (m = 0; m < n; ++m) buf[m].y += yoff;
					svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
				}

				/* Outward offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_out[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_out[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				/* Inward offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_in[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_in[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				svg_label(&svg, (qaws_scalar)-4.0, (qaws_scalar)(-3.0) + yoff,
					labels[row], "#8888aa");
			}
		}
		qaws_curve_destroy(crv);
		#undef CAYLEY_N
closed_end:;
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/27_offset_closed.svg\n");
}

static void test_svg_offset_selfintersect(void)
{
	printf("test_svg_offset_selfintersect\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/27_offset_selfintersect.svg",
		(qaws_scalar)-4.0, (qaws_scalar)-2.5, (qaws_scalar)8.0, (qaws_scalar)18.0,
		(qaws_scalar)600, (qaws_scalar)1350))
		return;

	/* Lemniscate of Bernoulli (figure-8, self-intersects at origin)
	 * x(t) = a*cos(t) / (1 + sin^2(t))
	 * y(t) = a*sin(t)*cos(t) / (1 + sin^2(t))
	 * t in [0, 2*pi] traces the full curve. */
	{
		#define LEMN_OFF_N 48
		qaws_scalar pts[LEMN_OFF_N * 2];
		qaws_catmull_rom_desc cr_desc;
		qaws_curve* crv = NULL;
		qaws_status s;
		unsigned int k;
		double a = 2.5;

		for (k = 0; k < LEMN_OFF_N; ++k)
		{
			double t = 2.0 * 3.14159265358979323846 * (double)k / (double)LEMN_OFF_N;
			double st = sin(t);
			double ct = cos(t);
			double denom = 1.0 + st * st;
			pts[k * 2 + 0] = (qaws_scalar)(a * ct / denom);
			pts[k * 2 + 1] = (qaws_scalar)(a * st * ct / denom);
		}

		memset(&cr_desc, 0, sizeof(cr_desc));
		cr_desc.dimension = QAWS_DIMENSION_2D;
		cr_desc.control_points = pts;
		cr_desc.control_point_count = LEMN_OFF_N;
		cr_desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cr_desc.closed = 1;
		s = qaws_curve_create_catmull_rom(&cr_desc, &crv);
		if (s != QAWS_STATUS_OK) goto si_end;

		/* Three rows stacked vertically: trim=0, trim=1, trim=2 */
		{
			static const qaws_scalar dists_l[] = {
				(qaws_scalar)0.2, (qaws_scalar)0.5, (qaws_scalar)0.8
			};
			static const qaws_scalar dists_r[] = {
				(qaws_scalar)-0.2, (qaws_scalar)-0.5, (qaws_scalar)-0.8
			};
			static const char* colors_l[] = { "#50fa7b", "#8be9fd", "#ff79c6" };
			static const char* colors_r[] = { "#6272a4", "#bd93f9", "#e94560" };
			static const int trim_vals[] = { 0, 1, 2 };
			static const char* labels[] = {
				"trim=0 (raw)", "trim=1 (distance)", "trim=2 (cleanup)"
			};
			int row;
			for (row = 0; row < 3; ++row)
			{
				int do_trim = trim_vals[row];
				qaws_scalar yoff = (qaws_scalar)row * (qaws_scalar)6.0;

				/* Draw base curve shifted */
				{
					qaws_vec2 buf[SVG_SAMPLES];
					unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES), m;
					for (m = 0; m < n; ++m) buf[m].y += yoff;
					svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
				}

				/* Left offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_l[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_l[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				/* Right offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_r[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_r[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				svg_label(&svg, (qaws_scalar)-3.5, (qaws_scalar)(-2.0) + yoff,
					labels[row], "#8888aa");
			}
		}
		qaws_curve_destroy(crv);
		#undef LEMN_OFF_N
si_end:;
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/27_offset_selfintersect.svg\n");
}

int test_27_svg_analysis_main(void)
{
	svg_ensure_output_dir();
	test_svg_inflection_points();
	test_svg_extrema();
	test_svg_curvature_comb();
	test_svg_winding_number();
	test_svg_curvature_weighted();
	test_svg_feature_preserving();
	test_svg_multi_traversal();
	test_svg_frenet_frame();
	test_svg_self_intersection();
	test_svg_curve_curve_intersection();
	test_svg_offset();
	test_svg_offset_closed();
	test_svg_offset_selfintersect();
	return 0;
}
