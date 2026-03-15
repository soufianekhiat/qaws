#include "test_common.h"

static void test_svg_circle(void)
{
	printf("test_svg_circle\n");

	/* Full circle via 4 quarter-arc NURBS patches joined together.
	   We'll draw each quarter as a separate NURBS. */
	qaws_scalar w = (qaws_scalar)0.70710678118; /* 1/sqrt(2) */
	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_circle.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.5, (qaws_scalar)3.0, (qaws_scalar)3.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	char const* colors[] = { "#e94560", "#f5a623", "#50fa7b", "#6272a4" };

	/* Quarter arcs: Q0 (1,0)->(0,1), Q1 (0,1)->(-1,0), Q2 (-1,0)->(0,-1), Q3 (0,-1)->(1,0) */
	qaws_vec2 q_cp[4][3] = {
		{ {1, 0}, {1, 1}, {0, 1} },
		{ {0, 1}, {-1, 1}, {-1, 0} },
		{ {-1, 0}, {-1, -1}, {0, -1} },
		{ {0, -1}, {1, -1}, {1, 0} }
	};
	qaws_scalar q_w[] = { 1, w, 1 };
	qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };

	qaws_vec2 all_pts[SVG_SAMPLES * 4];
	unsigned int total = 0;

	for (int qi = 0; qi < 4; qi++) {
		qaws_nurbs_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = q_cp[qi];
		desc.control_point_count = 3;
		desc.knots = knots;
		desc.knot_count = 6;
		desc.weights = q_w;
		desc.weight_count = 3;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_nurbs(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, colors[qi], (qaws_scalar)2.5);
		for (unsigned int i = 0; i < n; i++) all_pts[total++] = buf[i];
		svg_dots(&svg, q_cp[qi], 3, "#ffffff", (qaws_scalar)3);
		qaws_curve_destroy(curve);
	}

	/* Verify: check that sampled points are on the unit circle */
	int on_circle = 1;
	for (unsigned int i = 0; i < total; i++) {
		qaws_scalar r = (qaws_scalar)sqrt(
			(double)(all_pts[i].x * all_pts[i].x + all_pts[i].y * all_pts[i].y));
		if (r < (qaws_scalar)0.99 || r > (qaws_scalar)1.01) { on_circle = 0; break; }
	}
	TEST_ASSERT(on_circle, "circle: all samples on unit circle");

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)1.3, "Circle (4x NURBS quarter-arcs)", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/25_circle.svg\n");
}

static void test_svg_lemniscate(void)
{
	printf("test_svg_lemniscate\n");

	/* Lemniscate of Bernoulli (infinity symbol) approximated as a Bezier.
	   x(t) = cos(t) / (1 + sin^2(t))
	   y(t) = sin(t)*cos(t) / (1 + sin^2(t))
	   We sample the parametric form directly, then fit Bezier control points
	   for two lobes and draw them. Simpler: just build a Bezier curve that
	   approximates each half, using sampled points as guide. */

	/* Approach: we'll manually sample the lemniscate formula and create
	   a piecewise Bezier (degree-3) approximation by using Hermite
	   interpolation on several key points per lobe. */
	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_lemniscate.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-0.8, (qaws_scalar)3.0, (qaws_scalar)1.6,
		(qaws_scalar)600, (qaws_scalar)320))
		return;

	/* Parametric lemniscate: sample N points over t=[0, 2*pi] */
	#define LEMN_N 512
	qaws_vec2 lemn_pts[LEMN_N];
	{
		double pi2 = 6.283185307179586;
		for (int i = 0; i < LEMN_N; i++) {
			double t = pi2 * (double)i / (double)(LEMN_N - 1);
			double s = sin(t);
			double c = cos(t);
			double denom = 1.0 + s * s;
			lemn_pts[i].x = (qaws_scalar)(c / denom);
			lemn_pts[i].y = (qaws_scalar)(s * c / denom);
		}
	}

	/* Draw the raw parametric lemniscate as a polyline */
	svg_polyline(&svg, lemn_pts, LEMN_N, "#50fa7b", (qaws_scalar)1.0);

	/* Now create a Bezier approximation of the right lobe using Hermite.
	   Right lobe: from (1,0) up through top, back to (0,0), down to (1,0).
	   We'll use two Hermite spans. */
	{
		/* Right lobe: 3 key points with tangents */
		qaws_vec2 pts[] = {
			{(qaws_scalar)0.0, (qaws_scalar)0.0},
			{(qaws_scalar)1.0, (qaws_scalar)0.0},  /* rightmost */
			{(qaws_scalar)0.0, (qaws_scalar)0.0}
		};
		qaws_vec2 tans[] = {
			{(qaws_scalar)1.0, (qaws_scalar)1.0},
			{(qaws_scalar)0.0, (qaws_scalar)-1.5},
			{(qaws_scalar)-1.0, (qaws_scalar)1.0}
		};

		qaws_hermite_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 3;
		desc.derivative_count = 3;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_hermite(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.5);
		svg_dots(&svg, pts, 3, "#ffffff", (qaws_scalar)3);
		qaws_curve_destroy(curve);
	}

	/* Left lobe */
	{
		qaws_vec2 pts[] = {
			{(qaws_scalar)0.0, (qaws_scalar)0.0},
			{(qaws_scalar)-1.0, (qaws_scalar)0.0},
			{(qaws_scalar)0.0, (qaws_scalar)0.0}
		};
		qaws_vec2 tans[] = {
			{(qaws_scalar)-1.0, (qaws_scalar)-1.0},
			{(qaws_scalar)0.0, (qaws_scalar)1.5},
			{(qaws_scalar)1.0, (qaws_scalar)-1.0}
		};

		qaws_hermite_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 3;
		desc.derivative_count = 3;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_hermite(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)2.5);
		qaws_curve_destroy(curve);
	}

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)0.65,
		"Lemniscate (Hermite lobes + parametric ref)", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/25_lemniscate.svg\n");
}

static void test_svg_bezier_s(void)
{
	printf("test_svg_bezier_s\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_bezier_s_curve.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)500, (qaws_scalar)400))
		return;

	/* S-shaped cubic Bezier */
	qaws_vec2 cp[] = { {0, 0}, {0, 3}, {4, 0}, {4, 3} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);

	/* Draw control polygon (dashed) */
	svg_polyline(&svg, cp, 4, "#444466", (qaws_scalar)1);
	/* Draw curve */
	svg_polyline(&svg, buf, n, "#ff6188", (qaws_scalar)2.5);
	/* Draw control points */
	svg_dots(&svg, cp, 4, "#a9dc76", (qaws_scalar)4);

	/* Verify endpoints */
	TEST_ASSERT(approx_eq(buf[0].x, (qaws_scalar)0.0), "s-curve start x");
	TEST_ASSERT(approx_eq(buf[n - 1].x, (qaws_scalar)4.0), "s-curve end x");

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5, "Cubic Bezier S-curve", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/25_bezier_s_curve.svg\n");
}

static void test_svg_catmull_rom_wave(void)
{
	printf("test_svg_catmull_rom_wave\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_catmull_rom_wave.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-2.0, (qaws_scalar)8.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)300))
		return;

	/* Wave-like Catmull-Rom through 6 control points */
	qaws_vec2 all_cp[] = {
		{0, 0}, {1, (qaws_scalar)1.5}, {2, (qaws_scalar)-1.5},
		{3, (qaws_scalar)1.5}, {4, (qaws_scalar)-1.5}, {5, (qaws_scalar)1.5},
		{6, 0}
	};

	/* Catmull-Rom needs at least 4 points for 1 interior span.
	   With 7 control points we get 5 interior spans (P1..P5). */
	qaws_catmull_rom_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = all_cp;
	desc.control_point_count = 7;
	desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_catmull_rom(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#78dce8", (qaws_scalar)2.5);

	/* Show interpolated points (interior: indices 1..5) */
	svg_dots(&svg, all_cp + 1, 5, "#ff6188", (qaws_scalar)4);
	/* Show ghost points */
	svg_dots(&svg, all_cp, 1, "#666688", (qaws_scalar)3);
	svg_dots(&svg, all_cp + 6, 1, "#666688", (qaws_scalar)3);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)1.7, "Catmull-Rom wave (centripetal)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/25_catmull_rom_wave.svg\n");
}

static void test_svg_bspline(void)
{
	printf("test_svg_bspline\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_bspline.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)500, (qaws_scalar)400))
		return;

	/* Cubic B-Spline with 6 control points */
	qaws_vec2 cp[] = {
		{0, 0}, {1, 3}, {2, 0}, {3, 3}, {4, 0}, {4, 3}
	};
	qaws_bspline_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 6;
	desc.knots = NULL;
	desc.knot_count = 0;
	desc.is_uniform = 1;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bspline(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);

	/* Control polygon */
	svg_polyline(&svg, cp, 6, "#444466", (qaws_scalar)1);
	/* Curve */
	svg_polyline(&svg, buf, n, "#ab9df2", (qaws_scalar)2.5);
	/* Control points */
	svg_dots(&svg, cp, 6, "#a9dc76", (qaws_scalar)4);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5, "Cubic B-Spline (6 CPs, uniform clamped)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/25_bspline.svg\n");
}

static void test_svg_nurbs_ellipse(void)
{
	printf("test_svg_nurbs_ellipse\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_nurbs_ellipse.svg",
		(qaws_scalar)-3.0, (qaws_scalar)-2.0, (qaws_scalar)6.0, (qaws_scalar)4.0,
		(qaws_scalar)500, (qaws_scalar)340))
		return;

	/* Ellipse: scale unit circle x by 2 */
	qaws_scalar w = (qaws_scalar)0.70710678118;
	char const* colors[] = { "#e94560", "#f5a623", "#50fa7b", "#6272a4" };

	qaws_vec2 q_cp[4][3] = {
		{ {2, 0}, {2, 1}, {0, 1} },
		{ {0, 1}, {-2, 1}, {-2, 0} },
		{ {-2, 0}, {-2, -1}, {0, -1} },
		{ {0, -1}, {2, -1}, {2, 0} }
	};
	qaws_scalar q_w[] = { 1, w, 1 };
	qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };

	for (int qi = 0; qi < 4; qi++) {
		qaws_nurbs_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = q_cp[qi];
		desc.control_point_count = 3;
		desc.knots = knots;
		desc.knot_count = 6;
		desc.weights = q_w;
		desc.weight_count = 3;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_nurbs(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, colors[qi], (qaws_scalar)2.5);
		svg_dots(&svg, q_cp[qi], 3, "#ffffff80", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	svg_label(&svg, (qaws_scalar)-2.8, (qaws_scalar)1.7, "Ellipse (4x NURBS, a=2, b=1)", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/25_nurbs_ellipse.svg\n");
}

static void test_svg_hermite_figure8(void)
{
	printf("test_svg_hermite_figure8\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_hermite_figure8.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.5, (qaws_scalar)3.0, (qaws_scalar)3.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	/* Figure-8 using 5-point Hermite loop:
	   origin -> right lobe top -> origin -> left lobe bottom -> origin */
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

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#ff6188", (qaws_scalar)2.5);
	svg_dots(&svg, pts, 5, "#a9dc76", (qaws_scalar)4);

	/* Label each control point with its index.
	   Points 0, 2, 4 all sit at the origin — offset labels so all 5 are visible. */
	{
		qaws_scalar lx[] = {
			(qaws_scalar) 0.08, (qaws_scalar) 0.08,
			(qaws_scalar)-0.20, (qaws_scalar) 0.08,
			(qaws_scalar) 0.08
		};
		qaws_scalar ly[] = {
			(qaws_scalar) 0.12, (qaws_scalar) 0.10,
			(qaws_scalar)-0.06, (qaws_scalar) 0.10,
			(qaws_scalar)-0.12
		};
		for (int i = 0; i < 5; i++) {
			char lb[4];
			sprintf(lb, "%d", i);
			svg_label(&svg, pts[i].x + lx[i], pts[i].y + ly[i], lb, "#a9dc76");
		}
	}

	/* Verify: starts and ends at origin */
	TEST_ASSERT(approx_eq(buf[0].x, (qaws_scalar)0.0), "figure8 start x");
	TEST_ASSERT(approx_eq(buf[n - 1].x, (qaws_scalar)0.0), "figure8 end x");

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)1.3, "Figure-8 (Hermite, 5 key points: 0,2,4 at origin)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/25_hermite_figure8.svg\n");
}

static void test_svg_all_families(void)
{
	printf("test_svg_all_families\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_all_families.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.0, (qaws_scalar)5.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n;

	/* 1. Bezier cubic */
	{
		qaws_vec2 cp[] = { {0, 0}, {0, 3}, {4, 0}, {4, 3} };
		qaws_bezier_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;

		qaws_curve* curve = NULL;
		qaws_curve_create_bezier(&desc, &curve);
		n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#ff6188", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	/* 2. Hermite */
	{
		qaws_vec2 pts[] = { {0, (qaws_scalar)-0.5}, {4, (qaws_scalar)-0.5} };
		qaws_vec2 tans[] = { {2, 2}, {2, -2} };
		qaws_hermite_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 2;
		desc.derivative_count = 2;

		qaws_curve* curve = NULL;
		qaws_curve_create_hermite(&desc, &curve);
		n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	/* 3. B-Spline */
	{
		qaws_vec2 cp[] = { {0, (qaws_scalar)-1}, {1, 2}, {2, (qaws_scalar)-1}, {3, 2}, {4, (qaws_scalar)-1} };
		qaws_bspline_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 5;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_uniform = 1;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_curve_create_bspline(&desc, &curve);
		n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#78dce8", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	/* Legend */
	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)3.3, "Bezier", "#ff6188");
	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)3.0, "Hermite", "#f5a623");
	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)2.7, "B-Spline", "#78dce8");

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.7, "All families comparison", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/25_all_families.svg\n");
}

static void test_svg_yuksel_bezier(void)
{
	printf("test_svg_yuksel_bezier\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, -1}, {3, 2}, {4, -1}, {5, 2}, {6, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 7;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_curve_create_yuksel(&desc, &curve);

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_yuksel_bezier.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-2.0, (qaws_scalar)7.0, (qaws_scalar)5.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{ qaws_curve_destroy(curve); return; }

	qaws_vec2 samples[SVG_SAMPLES];
	unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
	svg_polyline(&svg, samples, count, "#e94560", (qaws_scalar)2);
	svg_dots(&svg, points, 7, "#50fa7b", (qaws_scalar)5);
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.8, "Yuksel Bezier", "#e94560");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/25_yuksel_bezier.svg\n");
	qaws_curve_destroy(curve);
}

/* SVG: Yuksel closed pentagon - all 4 modes comparison */
static void test_svg_yuksel_modes(void)
{
	printf("test_svg_yuksel_modes\n");

	qaws_vec2 points[6];
	{
		unsigned int pi;
		for (pi = 0; pi < 6; pi++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979323846) * (qaws_scalar)pi / (qaws_scalar)6;
			points[pi].x = (qaws_scalar)(2.0 * cos((double)angle));
			points[pi].y = (qaws_scalar)(2.0 * sin((double)angle));
		}
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_yuksel_modes.svg",
		(qaws_scalar)-3.0, (qaws_scalar)-3.0, (qaws_scalar)6.0, (qaws_scalar)6.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	char const* colors[] = { "#e94560", "#f5a623", "#50fa7b", "#78dce8" };
	char const* names[] = { "Bezier", "Circular", "Elliptical", "Hybrid" };
	int modes[] = {
		QAWS_YUKSEL_MODE_BEZIER,
		QAWS_YUKSEL_MODE_CIRCULAR,
		QAWS_YUKSEL_MODE_ELLIPTICAL,
		QAWS_YUKSEL_MODE_HYBRID
	};

	unsigned int m;
	for (m = 0; m < 4; m++) {
		qaws_yuksel_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 6;
		desc.mode = (qaws_yuksel_mode)modes[m];
		desc.closed = 1;

		qaws_curve* curve = NULL;
		if (qaws_curve_create_yuksel(&desc, &curve) == QAWS_STATUS_OK) {
			qaws_vec2 samples[SVG_SAMPLES];
			unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
			svg_polyline(&svg, samples, count, colors[m], (qaws_scalar)1.5);
			qaws_curve_destroy(curve);
		}

		svg_label(&svg, (qaws_scalar)-2.8, (qaws_scalar)2.8 - (qaws_scalar)m * (qaws_scalar)0.35,
			names[m], colors[m]);
	}

	svg_dots(&svg, points, 6, "#ffffff", (qaws_scalar)4);
	svg_label(&svg, (qaws_scalar)-2.8, (qaws_scalar)3.2, "Yuksel modes (hexagon)", "#8888aa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/25_yuksel_modes.svg\n");
}

/* SVG: Yuksel vs Catmull-Rom comparison */
static void test_svg_yuksel_vs_catmull(void)
{
	printf("test_svg_yuksel_vs_catmull\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, -0.5}, {3, 1.5}, {4, 0}
	};

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/25_yuksel_vs_catmull.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.5, (qaws_scalar)4.5,
		(qaws_scalar)600, (qaws_scalar)400))
		return;

	/* Yuksel Bezier */
	{
		qaws_yuksel_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_BEZIER;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		if (qaws_curve_create_yuksel(&desc, &curve) == QAWS_STATUS_OK) {
			qaws_vec2 samples[SVG_SAMPLES];
			unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
			svg_polyline(&svg, samples, count, "#e94560", (qaws_scalar)2);
			qaws_curve_destroy(curve);
		}
	}

	/* Catmull-Rom centripetal (needs 2 extra ghost points) */
	{
		qaws_vec2 cr_pts[] = {
			{-1, -2}, {0, 0}, {1, 2}, {2, -0.5}, {3, 1.5}, {4, 0}, {5, -2}
		};
		qaws_catmull_rom_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = cr_pts;
		desc.control_point_count = 7;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		if (qaws_curve_create_catmull_rom(&desc, &curve) == QAWS_STATUS_OK) {
			qaws_vec2 samples[SVG_SAMPLES];
			unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
			svg_polyline(&svg, samples, count, "#50fa7b", (qaws_scalar)2);
			qaws_curve_destroy(curve);
		}
	}

	svg_dots(&svg, points, 5, "#ffffff", (qaws_scalar)5);
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.8, "Yuksel (red) vs Catmull-Rom (green)", "#8888aa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/25_yuksel_vs_catmull.svg\n");
}

int test_25_svg_curves_main(void)
{
	svg_ensure_output_dir();
	test_svg_circle();
	test_svg_lemniscate();
	test_svg_bezier_s();
	test_svg_catmull_rom_wave();
	test_svg_bspline();
	test_svg_nurbs_ellipse();
	test_svg_hermite_figure8();
	test_svg_all_families();
	test_svg_yuksel_bezier();
	test_svg_yuksel_modes();
	test_svg_yuksel_vs_catmull();
	return 0;
}
