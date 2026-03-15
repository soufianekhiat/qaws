#include "test_common.h"

static void test_svg_curve_reverse(void)
{
	printf("test_svg_curve_reverse\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_curve_reverse.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.5, (qaws_scalar)8.0, (qaws_scalar)8.0,
		(qaws_scalar)600, (qaws_scalar)600))
		return;

	/* --- Row 1: S-curve Bezier --- */
	{
		qaws_vec2 cp[] = { {0,2}, {0,5}, {4,-1}, {4,2} };
		qaws_bezier_desc desc;
		qaws_curve *crv = NULL, *rev = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;

		if (qaws_curve_create_bezier(&desc, &crv) == QAWS_STATUS_OK &&
			qaws_curve_reverse(crv, &rev) == QAWS_STATUS_OK) {

			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n;
			n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)3.0);
			svg_arrow(&svg, crv, (qaws_scalar)0.3, "#e94560", (qaws_scalar)0.3);
			svg_arrow(&svg, crv, (qaws_scalar)0.7, "#e94560", (qaws_scalar)0.3);

			n = svg_sample_curve(rev, buf, SVG_SAMPLES);
			svg_polyline_dashed(&svg, buf, n, "#78dce8", (qaws_scalar)2.0, "6,4");
			svg_arrow(&svg, rev, (qaws_scalar)0.3, "#78dce8", (qaws_scalar)0.3);
			svg_arrow(&svg, rev, (qaws_scalar)0.7, "#78dce8", (qaws_scalar)0.3);

			/* Start/end markers */
			svg_ring(&svg, cp[0].x, cp[0].y, "#e94560", (qaws_scalar)5, (qaws_scalar)1.5);
			svg_dot(&svg, cp[0].x, cp[0].y, "#e94560", (qaws_scalar)2);
			svg_ring(&svg, cp[3].x, cp[3].y, "#78dce8", (qaws_scalar)5, (qaws_scalar)1.5);
			svg_dot(&svg, cp[3].x, cp[3].y, "#78dce8", (qaws_scalar)2);

			svg_label(&svg, (qaws_scalar)0.5, (qaws_scalar)4.2, "Cubic Bezier S-curve", "#8888aa");
			qaws_curve_destroy(rev);
			qaws_curve_destroy(crv);
		}
	}

	/* --- Row 2: Hermite curve --- */
	{
		qaws_scalar pts[] = { 0,-2,  2,0,  4,-2 };
		qaws_scalar derivs[] = { 2,2,  0,-3,  2,2 };
		qaws_hermite_desc desc;
		qaws_curve *crv = NULL, *rev = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = derivs;
		desc.point_count = 3;
		desc.derivative_count = 3;

		if (qaws_curve_create_hermite(&desc, &crv) == QAWS_STATUS_OK &&
			qaws_curve_reverse(crv, &rev) == QAWS_STATUS_OK) {

			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n;
			n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)3.0);
			svg_arrow(&svg, crv, (qaws_scalar)0.25, "#f5a623", (qaws_scalar)0.3);
			svg_arrow(&svg, crv, (qaws_scalar)0.75, "#f5a623", (qaws_scalar)0.3);

			n = svg_sample_curve(rev, buf, SVG_SAMPLES);
			svg_polyline_dashed(&svg, buf, n, "#ab9df2", (qaws_scalar)2.0, "6,4");
			svg_arrow(&svg, rev, (qaws_scalar)0.25, "#ab9df2", (qaws_scalar)0.3);
			svg_arrow(&svg, rev, (qaws_scalar)0.75, "#ab9df2", (qaws_scalar)0.3);

			svg_label(&svg, (qaws_scalar)0.5, (qaws_scalar)-0.3, "Hermite 2-span", "#8888aa");
			qaws_curve_destroy(rev);
			qaws_curve_destroy(crv);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.6, "Curve Reverse", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)-3.0,
		"Solid=original  Dashed=reversed  Arrows=direction", "#666688");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_curve_reverse.svg\n");
}

/* SVG: curve split - split cubic Bezier at t=0.5 */
static void test_svg_curve_split(void)
{
	printf("test_svg_curve_split\n");

	qaws_vec2 cp[] = { {0,0}, {1,4}, {3,-2}, {4,2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	qaws_curve* left = NULL;
	qaws_curve* right = NULL;
	s = qaws_curve_split(curve, (qaws_scalar)0.5, &left, &right);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: split failed\n");
		qaws_curve_destroy(curve);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_curve_split.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.0, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)600, (qaws_scalar)500))
	{
		qaws_curve_destroy(curve);
		qaws_curve_destroy(left);
		qaws_curve_destroy(right);
		return;
	}

	/* Original curve faintly */
	qaws_vec2 samples_full[SVG_SAMPLES];
	unsigned int n_full = svg_sample_curve(curve, samples_full, SVG_SAMPLES);
	svg_polyline(&svg, samples_full, n_full, "#666666", (qaws_scalar)1.0);

	/* Left half in red */
	qaws_vec2 samples_left[SVG_SAMPLES];
	unsigned int n_left = svg_sample_curve(left, samples_left, SVG_SAMPLES);
	svg_polyline(&svg, samples_left, n_left, "#e94560", (qaws_scalar)2.5);

	/* Right half in green */
	qaws_vec2 samples_right[SVG_SAMPLES];
	unsigned int n_right = svg_sample_curve(right, samples_right, SVG_SAMPLES);
	svg_polyline(&svg, samples_right, n_right, "#50fa7b", (qaws_scalar)2.5);

	/* Split point as large yellow dot */
	{
		qaws_eval_result_2d r;
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
		svg_dots(&svg, &r.position, 1, "#f1fa8c", (qaws_scalar)7);
	}

	/* Control points of original in white */
	svg_dots(&svg, cp, 4, "#ffffff", (qaws_scalar)4);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.5, "Curve Split at t=0.5", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.2, "Left (red)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.9, "Right (green)", "#50fa7b");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_curve_split.svg\n");
	qaws_curve_destroy(curve);
	qaws_curve_destroy(left);
	qaws_curve_destroy(right);
}

/* SVG: curve join - two Hermite curves joined end-to-end */
static void test_svg_curve_join(void)
{
	printf("test_svg_curve_join\n");

	/* Curve A: points (0,0),(2,2), tangents (2,1),(1,-1) */
	qaws_vec2 pts_a[] = { {0,0}, {2,2} };
	qaws_vec2 tan_a[] = { {2,1}, {1,-1} };
	qaws_hermite_desc desc_a;
	desc_a.dimension = QAWS_DIMENSION_2D;
	desc_a.degree = 3;
	desc_a.points = pts_a;
	desc_a.derivatives = tan_a;
	desc_a.point_count = 2;
	desc_a.derivative_count = 2;

	qaws_curve* curve_a = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc_a, &curve_a);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create A failed\n"); return; }

	/* Curve B: points (2,2),(4,0), tangents (1,-1),(2,1) */
	qaws_vec2 pts_b[] = { {2,2}, {4,0} };
	qaws_vec2 tan_b[] = { {1,-1}, {2,1} };
	qaws_hermite_desc desc_b;
	desc_b.dimension = QAWS_DIMENSION_2D;
	desc_b.degree = 3;
	desc_b.points = pts_b;
	desc_b.derivatives = tan_b;
	desc_b.point_count = 2;
	desc_b.derivative_count = 2;

	qaws_curve* curve_b = NULL;
	s = qaws_curve_create_hermite(&desc_b, &curve_b);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: create B failed\n");
		qaws_curve_destroy(curve_a);
		return;
	}

	qaws_curve* joined = NULL;
	s = qaws_curve_join(curve_a, curve_b, &joined);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: join failed\n");
		qaws_curve_destroy(curve_a);
		qaws_curve_destroy(curve_b);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_curve_join.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.0, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{
		qaws_curve_destroy(curve_a);
		qaws_curve_destroy(curve_b);
		qaws_curve_destroy(joined);
		return;
	}

	/* Draw joined curve in green (thicker, underneath) */
	qaws_vec2 samples_joined[SVG_SAMPLES];
	unsigned int n_joined = svg_sample_curve(joined, samples_joined, SVG_SAMPLES);
	svg_polyline(&svg, samples_joined, n_joined, "#50fa7b", (qaws_scalar)4.0);

	/* Draw curve A in red */
	qaws_vec2 samples_a[SVG_SAMPLES];
	unsigned int n_a = svg_sample_curve(curve_a, samples_a, SVG_SAMPLES);
	svg_polyline(&svg, samples_a, n_a, "#e94560", (qaws_scalar)2.0);

	/* Draw curve B in blue */
	qaws_vec2 samples_b[SVG_SAMPLES];
	unsigned int n_b = svg_sample_curve(curve_b, samples_b, SVG_SAMPLES);
	svg_polyline(&svg, samples_b, n_b, "#6272a4", (qaws_scalar)2.0);

	/* Key points */
	qaws_vec2 key_pts[] = { {0,0}, {2,2}, {4,0} };
	svg_dots(&svg, key_pts, 3, "#ffffff", (qaws_scalar)5);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Curve Join", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Curve A (red)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.6, "Curve B (blue)", "#6272a4");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.3, "Joined (green)", "#50fa7b");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_curve_join.svg\n");
	qaws_curve_destroy(curve_a);
	qaws_curve_destroy(curve_b);
	qaws_curve_destroy(joined);
}

/* SVG: degree elevation - quadratic elevated to cubic should overlap */
static void test_svg_degree_elevation(void)
{
	printf("test_svg_degree_elevation\n");

	qaws_vec2 cp[] = { {0,0}, {2,3}, {4,0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 2;
	desc.control_points = cp;
	desc.control_point_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	qaws_curve* elevated = NULL;
	s = qaws_curve_elevate_degree(curve, &elevated);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: elevate failed\n");
		qaws_curve_destroy(curve);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_degree_elevation.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{
		qaws_curve_destroy(curve);
		qaws_curve_destroy(elevated);
		return;
	}

	/* Original quadratic in red */
	qaws_vec2 samples_orig[SVG_SAMPLES];
	unsigned int n_orig = svg_sample_curve(curve, samples_orig, SVG_SAMPLES);
	svg_polyline(&svg, samples_orig, n_orig, "#e94560", (qaws_scalar)3.0);

	/* Elevated cubic in blue (slightly thinner) */
	qaws_vec2 samples_elev[SVG_SAMPLES];
	unsigned int n_elev = svg_sample_curve(elevated, samples_elev, SVG_SAMPLES);
	svg_polyline(&svg, samples_elev, n_elev, "#6272a4", (qaws_scalar)1.5);

	/* Original CPs as red dots */
	svg_dots(&svg, cp, 3, "#e94560", (qaws_scalar)5);

	/* Elevated CPs as blue dots */
	{
		qaws_vec2 elev_cp[8];
		unsigned int cp_count = 0;
		s = qaws_bezier_get_control_points(elevated, elev_cp, 8, &cp_count);
		if (s == QAWS_STATUS_OK && cp_count > 0) {
			svg_dots(&svg, elev_cp, cp_count, "#6272a4", (qaws_scalar)5);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5, "Degree Elevation", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Quadratic (red, 3 CPs)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Cubic (blue, 4 CPs)", "#6272a4");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_degree_elevation.svg\n");
	qaws_curve_destroy(curve);
	qaws_curve_destroy(elevated);
}

/* SVG: family conversion - Hermite to Bezier should overlap exactly */
static void test_svg_family_conversion(void)
{
	printf("test_svg_family_conversion\n");

	qaws_vec2 pts[] = { {0,0}, {4,2} };
	qaws_vec2 tans[] = { {2,4}, {2,-4} };
	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = pts;
	desc.derivatives = tans;
	desc.point_count = 2;
	desc.derivative_count = 2;

	qaws_curve* hermite = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc, &hermite);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create hermite failed\n"); return; }

	qaws_curve* bezier = NULL;
	s = qaws_curve_convert_hermite_to_bezier(hermite, 0, &bezier);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: convert failed\n");
		qaws_curve_destroy(hermite);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_family_conversion.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.0, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{
		qaws_curve_destroy(hermite);
		qaws_curve_destroy(bezier);
		return;
	}

	/* Original Hermite in red */
	qaws_vec2 samples_herm[SVG_SAMPLES];
	unsigned int n_herm = svg_sample_curve(hermite, samples_herm, SVG_SAMPLES);
	svg_polyline(&svg, samples_herm, n_herm, "#e94560", (qaws_scalar)3.0);

	/* Converted Bezier in blue */
	qaws_vec2 samples_bez[SVG_SAMPLES];
	unsigned int n_bez = svg_sample_curve(bezier, samples_bez, SVG_SAMPLES);
	svg_polyline(&svg, samples_bez, n_bez, "#6272a4", (qaws_scalar)1.5);

	/* Hermite points as red dots */
	svg_dots(&svg, pts, 2, "#e94560", (qaws_scalar)5);

	/* Bezier CPs as blue dots */
	{
		qaws_vec2 bez_cp[8];
		unsigned int cp_count = 0;
		s = qaws_bezier_get_control_points(bezier, bez_cp, 8, &cp_count);
		if (s == QAWS_STATUS_OK && cp_count > 0) {
			svg_dots(&svg, bez_cp, cp_count, "#6272a4", (qaws_scalar)5);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Family Conversion", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Hermite (red)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.6, "Bezier (blue)", "#6272a4");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_family_conversion.svg\n");
	qaws_curve_destroy(hermite);
	qaws_curve_destroy(bezier);
}

/* SVG: easing curves - show different easing modes as dot distributions along a curve */
static void test_svg_easing_curves(void)
{
	printf("test_svg_easing_curves\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_easing_curves.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.5, (qaws_scalar)7.5,
		(qaws_scalar)550, (qaws_scalar)750))
		return;

	/* Each easing gets its own horizontal line at a different Y level */
	qaws_easing modes[] = {
		QAWS_EASING_LINEAR,
		QAWS_EASING_QUAD_IN,
		QAWS_EASING_QUAD_OUT,
		QAWS_EASING_QUAD_IN_OUT,
		QAWS_EASING_CUBIC_IN,
		QAWS_EASING_CUBIC_OUT,
		QAWS_EASING_CUBIC_IN_OUT,
		QAWS_EASING_SINE_IN_OUT
	};
	char const* names[] = {
		"Linear", "Quad In", "Quad Out", "Quad InOut",
		"Cubic In", "Cubic Out", "Cubic InOut", "Sine InOut"
	};
	char const* colors[] = {
		"#ffffff", "#e94560", "#50fa7b", "#f5a623",
		"#78dce8", "#ab9df2", "#ff6188", "#6272a4"
	};
	unsigned int mode_count = 8;

	/* Use a simple straight line so we can see the easing effect purely */
	qaws_vec2 line_cp[] = { {0,0}, {4,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *line = NULL;
	qaws_scalar arc_len = 0;
	qaws_status s;
	unsigned int mi;

	memset(&bdesc, 0, sizeof(bdesc));
	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 1;
	bdesc.control_points = line_cp;
	bdesc.control_point_count = 2;
	s = qaws_curve_create_bezier(&bdesc, &line);
	if (s != QAWS_STATUS_OK) return;

	{
		qaws_range pr = qaws_curve_get_parameter_range(line);
		qaws_curve_compute_arc_length(line, pr.min_value, pr.max_value, &arc_len);
	}

	#define EASE_DOTS 24
	for (mi = 0; mi < mode_count; mi++) {
		qaws_scalar row_y = (qaws_scalar)(mode_count - 1 - mi) * (qaws_scalar)0.85;
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		unsigned int di;

		/* Draw baseline */
		svg_line(&svg, (qaws_scalar)0, row_y, (qaws_scalar)4, row_y,
			"#2a2a4a", (qaws_scalar)0.5);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = arc_len;
		tdesc.start_time = (qaws_scalar)0;
		tdesc.end_time = (qaws_scalar)1;
		tdesc.clamp_to_domain = 1;
		tdesc.easing = modes[mi];
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		s = qaws_traversal_create(line, &tdesc, &trav);
		if (s != QAWS_STATUS_OK) continue;

		for (di = 0; di < EASE_DOTS; di++) {
			qaws_scalar t = (qaws_scalar)di / (qaws_scalar)(EASE_DOTS - 1);
			qaws_eval_result_2d r;
			qaws_traversal_evaluate_2d(trav, t, QAWS_EVAL_FLAG_POSITION, &r);
			svg_dot(&svg, r.position.x, row_y, colors[mi], (qaws_scalar)3.5);
		}

		svg_label(&svg, (qaws_scalar)4.2, row_y + (qaws_scalar)0.1,
			names[mi], colors[mi]);

		qaws_traversal_destroy(trav);
	}
	#undef EASE_DOTS

	svg_label(&svg, (qaws_scalar)0.5, (qaws_scalar)7.0,
		"Easing Functions (dot spacing = speed)", "#8888aa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_easing_curves.svg\n");
	qaws_curve_destroy(line);
}

/* SVG: adaptive vs uniform sampling on a high-curvature curve */
static void test_svg_adaptive_sampling(void)
{
	printf("test_svg_adaptive_sampling\n");

	/* Cubic Bezier with tight curvature */
	qaws_vec2 cp[] = { {0,0}, {0,3}, {1,-1}, {4,2} };
	qaws_bezier_desc bdesc;
	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 3;
	bdesc.control_points = cp;
	bdesc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&bdesc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_adaptive_sampling.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.0, (qaws_scalar)4.5,
		(qaws_scalar)600, (qaws_scalar)400))
	{ qaws_curve_destroy(curve); return; }

	/* Draw the curve in gray */
	{
		qaws_vec2 full_samples[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, full_samples, SVG_SAMPLES);
		svg_polyline(&svg, full_samples, n, "#444466", (qaws_scalar)1.5);
	}

	/* Uniform sampling with 8 points — few enough that curvature gaps are visible */
	{
		qaws_sampling_desc sdesc;
		qaws_vec2 uniform_pts[64];
		unsigned int uniform_count = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 0;
		sdesc.error_tolerance = 0;

		s = qaws_curve_sample_2d(curve, &sdesc, uniform_pts, 64, &uniform_count);
		if (s == QAWS_STATUS_OK) {
			svg_dots(&svg, uniform_pts, uniform_count, "#e94560", (qaws_scalar)5);
		}
	}

	/* Adaptive sampling — starts with 8 segments, subdivides where chord error > tolerance */
	{
		qaws_sampling_desc sdesc;
		qaws_vec2 adaptive_pts[512];
		unsigned int adaptive_count = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.05;

		s = qaws_curve_sample_2d(curve, &sdesc, adaptive_pts, 512, &adaptive_count);
		if (s == QAWS_STATUS_OK) {
			svg_dots(&svg, adaptive_pts, adaptive_count, "#6272a4", (qaws_scalar)3);
		}
	}

	/* Control points */
	svg_dots(&svg, cp, 4, "#ffffff", (qaws_scalar)3);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Adaptive vs Uniform Sampling", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Uniform 8 pts (red, large)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.6, "Adaptive 8+subdivisions (blue, small)", "#6272a4");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_adaptive_sampling.svg\n");
	qaws_curve_destroy(curve);
}

static void test_svg_bspline_fitting(void)
{
	printf("test_svg_bspline_fitting\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_bspline_fitting.svg",
		(qaws_scalar)-1.8, (qaws_scalar)-2.0, (qaws_scalar)8.5, (qaws_scalar)5.5,
		(qaws_scalar)800, (qaws_scalar)520))
		return;

	/* --- Row 1: Semicircle fit --- */
	{
		qaws_vec2 pts[21];
		unsigned int i;
		qaws_bspline_fit_desc desc;
		qaws_curve *crv = NULL;

		for (i = 0; i < 21; ++i) {
			qaws_scalar angle = (qaws_scalar)3.14159265358979 * (qaws_scalar)i / (qaws_scalar)20;
			pts[i].x = (qaws_scalar)cos((double)angle);
			pts[i].y = (qaws_scalar)sin((double)angle) + (qaws_scalar)1.5;
		}

		/* Ideal semicircle reference */
		{
			qaws_vec2 ref[101];
			unsigned int ri;
			for (ri = 0; ri <= 100; ++ri) {
				qaws_scalar a = (qaws_scalar)3.14159265358979 * (qaws_scalar)ri / (qaws_scalar)100;
				ref[ri].x = (qaws_scalar)cos((double)a);
				ref[ri].y = (qaws_scalar)sin((double)a) + (qaws_scalar)1.5;
			}
			svg_polyline_dashed(&svg, ref, 101, "#444466", (qaws_scalar)1.0, "4,4");
		}

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.data_points = pts;
		desc.data_point_count = 21;
		desc.degree = 3;
		desc.control_point_count = 8;

		if (qaws_curve_fit_bspline(&desc, &crv) == QAWS_STATUS_OK) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)2.5);
			qaws_curve_destroy(crv);
		}

		for (i = 0; i < 21; ++i) {
			svg_ring(&svg, pts[i].x, pts[i].y, "#f5a623", (qaws_scalar)5, (qaws_scalar)1.5);
			svg_dot(&svg, pts[i].x, pts[i].y, "#f5a623", (qaws_scalar)2);
		}
		svg_label(&svg, (qaws_scalar)-1.5, (qaws_scalar)3.0,
			"Semicircle (21 pts -> 8 CPs)", "#8888aa");
	}

	/* --- Row 2: Noisy sine wave with varying CP count --- */
	{
		#define FIT_PTS 31
		qaws_vec2 pts[FIT_PTS];
		unsigned int i;

		for (i = 0; i < FIT_PTS; ++i) {
			qaws_scalar t = (qaws_scalar)i / (qaws_scalar)(FIT_PTS - 1);
			pts[i].x = t * (qaws_scalar)6.0;
			pts[i].y = (qaws_scalar)sin(2.5 * 3.14159265358979 * (double)t) * (qaws_scalar)0.8
				- (qaws_scalar)1.0;
			/* Add slight noise */
			pts[i].y += ((qaws_scalar)((int)(i * 7 + 3) % 11) / (qaws_scalar)11 - (qaws_scalar)0.5) * (qaws_scalar)0.15;
		}

		/* Data points as small white dots */
		for (i = 0; i < FIT_PTS; ++i)
			svg_dot(&svg, pts[i].x, pts[i].y, "#ffffff", (qaws_scalar)2.5);

		/* Fit with 6 CPs (underfitting) */
		{
			qaws_bspline_fit_desc desc;
			qaws_curve *crv = NULL;
			memset(&desc, 0, sizeof(desc));
			desc.dimension = QAWS_DIMENSION_2D;
			desc.data_points = pts;
			desc.data_point_count = FIT_PTS;
			desc.degree = 3;
			desc.control_point_count = 6;
			if (qaws_curve_fit_bspline(&desc, &crv) == QAWS_STATUS_OK) {
				qaws_vec2 buf[SVG_SAMPLES];
				unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
				svg_polyline_dashed(&svg, buf, n, "#6272a4", (qaws_scalar)2.0, "6,4");
				qaws_curve_destroy(crv);
			}
		}

		/* Fit with 12 CPs (good fit) */
		{
			qaws_bspline_fit_desc desc;
			qaws_curve *crv = NULL;
			memset(&desc, 0, sizeof(desc));
			desc.dimension = QAWS_DIMENSION_2D;
			desc.data_points = pts;
			desc.data_point_count = FIT_PTS;
			desc.degree = 3;
			desc.control_point_count = 12;
			if (qaws_curve_fit_bspline(&desc, &crv) == QAWS_STATUS_OK) {
				qaws_vec2 buf[SVG_SAMPLES];
				unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.5);
				qaws_curve_destroy(crv);
			}
		}

		svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)-1.6,
			"Noisy sine: dashed blue=6 CPs, red=12 CPs, white=data", "#666688");

		#undef FIT_PTS
	}

	svg_label(&svg, (qaws_scalar)-1.5, (qaws_scalar)3.3,
		"B-spline least-squares fitting", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_bspline_fitting.svg\n");
}

static void test_svg_arc_length_reparam(void)
{
	printf("test_svg_arc_length_reparam\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_arc_length_reparam.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)3.0, (qaws_scalar)3.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	/* Create a cubic Bezier with non-uniform speed */
	{
		qaws_vec2 cp[] = { {0,0}, {0,2}, {2,2}, {2,0} };
		qaws_bezier_desc desc;
		qaws_curve *src = NULL, *rep = NULL;
		qaws_range rng;
		unsigned int i;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &src);

		/* Draw source curve */
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(src, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#444466", (qaws_scalar)1.5);
		}

		/* Draw uniform parametric samples (non-uniform spacing) */
		for (i = 0; i <= 10; ++i) {
			qaws_eval_result_2d r;
			qaws_scalar t = (qaws_scalar)i / (qaws_scalar)10;
			qaws_curve_evaluate_2d(src, t, QAWS_EVAL_FLAG_POSITION, &r);
			svg_dot(&svg, r.position.x, r.position.y, "#e94560", (qaws_scalar)4);
		}

		/* Reparameterize and draw uniform arc-length samples */
		if (qaws_curve_reparameterize_arc_length(src, 512, &rep) == QAWS_STATUS_OK) {
			rng = qaws_curve_get_parameter_range(rep);
			for (i = 0; i <= 10; ++i) {
				qaws_eval_result_2d r;
				qaws_scalar s = rng.min_value + (rng.max_value - rng.min_value) * (qaws_scalar)i / (qaws_scalar)10;
				qaws_curve_evaluate_2d(rep, s, QAWS_EVAL_FLAG_POSITION, &r);
				svg_dot(&svg, r.position.x, r.position.y, "#50fa7b", (qaws_scalar)4);
			}
			qaws_curve_destroy(rep);
		}

		svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.6,
			"Arc-length reparam", "#8888aa");
		svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.3,
			"Red=parametric  Green=arc-length", "#666688");

		svg_dots(&svg, cp, 4, "#ab9df2", (qaws_scalar)3);
		qaws_curve_destroy(src);
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_arc_length_reparam.svg\n");
}

static void test_svg_degree_reduction(void)
{
	printf("test_svg_degree_reduction\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/26_degree_reduction.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.5, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)600, (qaws_scalar)700))
		return;

	/* Two examples stacked: S-curve (top) and loop (bottom) */

	/* --- Top: S-curve cubic --- */
	{
		qaws_vec2 cp[] = { {0,2}, {0,4}, {4,0}, {4,2} };
		qaws_bezier_desc desc;
		qaws_curve *src = NULL, *elevated = NULL, *reduced = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &src);

		if (qaws_curve_elevate_degree(src, &elevated) == QAWS_STATUS_OK)
			qaws_curve_reduce_degree(elevated, &reduced);

		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(src, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)4.0);
		}
		if (elevated) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(elevated, buf, SVG_SAMPLES);
			svg_polyline_dashed(&svg, buf, n, "#78dce8", (qaws_scalar)2.5, "8,4");
		}
		if (reduced) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(reduced, buf, SVG_SAMPLES);
			svg_polyline_dashed(&svg, buf, n, "#ff6188", (qaws_scalar)2.0, "3,5");
			qaws_curve_destroy(reduced);
		}
		if (elevated) qaws_curve_destroy(elevated);

		svg_polyline(&svg, cp, 4, "#ffffff20", (qaws_scalar)0.8);
		svg_dots(&svg, cp, 4, "#50fa7b", (qaws_scalar)3);
		svg_label(&svg, (qaws_scalar)0.3, (qaws_scalar)3.7, "S-curve (deg 3->4->3)", "#8888aa");
		qaws_curve_destroy(src);
	}

	/* --- Bottom: quintic loop, reduce twice (5->4->3) --- */
	{
		qaws_vec2 cp[] = { {0,-2}, {1,1}, {2,-3}, {3,1}, {4,-3}, {4.5,-1} };
		qaws_bezier_desc desc;
		qaws_curve *src = NULL, *r1 = NULL, *r2 = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 5;
		desc.control_points = cp;
		desc.control_point_count = 6;
		qaws_curve_create_bezier(&desc, &src);

		qaws_curve_reduce_degree(src, &r1);
		if (r1) qaws_curve_reduce_degree(r1, &r2);

		/* Original quintic: thick green */
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(src, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)4.0);
		}
		/* Reduced to quartic: dashed orange */
		if (r1) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(r1, buf, SVG_SAMPLES);
			svg_polyline_dashed(&svg, buf, n, "#f5a623", (qaws_scalar)2.5, "8,4");
		}
		/* Reduced to cubic: dotted pink */
		if (r2) {
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(r2, buf, SVG_SAMPLES);
			svg_polyline_dashed(&svg, buf, n, "#ff6188", (qaws_scalar)2.0, "3,5");
			qaws_curve_destroy(r2);
		}
		if (r1) qaws_curve_destroy(r1);

		svg_polyline(&svg, cp, 6, "#ffffff20", (qaws_scalar)0.8);
		svg_dots(&svg, cp, 6, "#50fa7b", (qaws_scalar)3);
		svg_label(&svg, (qaws_scalar)0.3, (qaws_scalar)-0.3, "Complex curve (deg 5->4->3)", "#8888aa");
		qaws_curve_destroy(src);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.0,
		"Degree Reduction", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)-3.0,
		"Green=original  Cyan/Orange=1st reduce  Pink=2nd reduce", "#666688");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/26_degree_reduction.svg\n");
}

int test_26_svg_operations_main(void)
{
	svg_ensure_output_dir();
	test_svg_curve_reverse();
	test_svg_curve_split();
	test_svg_curve_join();
	test_svg_degree_elevation();
	test_svg_family_conversion();
	test_svg_easing_curves();
	test_svg_adaptive_sampling();
	test_svg_bspline_fitting();
	test_svg_arc_length_reparam();
	test_svg_degree_reduction();
	return 0;
}
