#include "test_common.h"
#include "qaws_boolean_2d.h"

/* ------------------------------------------------------------------ */
/*  Helper: create a circle curve in 2D                                */
/* ------------------------------------------------------------------ */
static qaws_curve* make_circle_2d(qaws_scalar cx, qaws_scalar cy, qaws_scalar r)
{
	qaws_arc_segment seg;
	qaws_arc_desc desc;
	qaws_curve* curve = NULL;

	memset(&seg, 0, sizeof(seg));
	seg.center[0] = cx;
	seg.center[1] = cy;
	seg.radius = r;
	seg.angle_start = (qaws_scalar)0.0;
	seg.angle_end = (qaws_scalar)6.283185307;

	memset(&desc, 0, sizeof(desc));
	desc.dimension = QAWS_DIMENSION_2D;
	desc.segments = &seg;
	desc.segment_count = 1;

	qaws_curve_create_arc(&desc, &curve);
	return curve;
}

/* ------------------------------------------------------------------ */
/*  Test: boolean union of two overlapping circles                     */
/* ------------------------------------------------------------------ */
static void test_boolean_union(void)
{
	qaws_curve* circle_a = NULL;
	qaws_curve* circle_b = NULL;
	qaws_curve* result[4];
	unsigned int result_count = 0;
	qaws_status s;

	printf("test_boolean_union\n");

	circle_a = make_circle_2d(0, 0, (qaws_scalar)1.0);
	circle_b = make_circle_2d((qaws_scalar)0.8, 0, (qaws_scalar)1.0);

	TEST_ASSERT(circle_a != NULL, "circle A created");
	TEST_ASSERT(circle_b != NULL, "circle B created");

	memset(result, 0, sizeof(result));
	s = qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_UNION,
		result, 4, &result_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(result_count >= 1, "union produced at least 1 boundary");
	TEST_ASSERT(result[0] != NULL, "union boundary is non-null");

	{
		unsigned int ri;
		for (ri = 0; ri < result_count; ri++)
			qaws_curve_destroy(result[ri]);
	}
	qaws_curve_destroy(circle_a);
	qaws_curve_destroy(circle_b);
}

/* ------------------------------------------------------------------ */
/*  Test: boolean intersection of two overlapping circles              */
/* ------------------------------------------------------------------ */
static void test_boolean_intersection(void)
{
	qaws_curve* circle_a = NULL;
	qaws_curve* circle_b = NULL;
	qaws_curve* result[4];
	unsigned int result_count = 0;
	qaws_status s;

	printf("test_boolean_intersection\n");

	circle_a = make_circle_2d(0, 0, (qaws_scalar)1.0);
	circle_b = make_circle_2d((qaws_scalar)0.8, 0, (qaws_scalar)1.0);

	memset(result, 0, sizeof(result));
	s = qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_INTERSECTION,
		result, 4, &result_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(result_count >= 1, "intersection produced at least 1 boundary");

	{
		unsigned int ri;
		for (ri = 0; ri < result_count; ri++)
			qaws_curve_destroy(result[ri]);
	}
	qaws_curve_destroy(circle_a);
	qaws_curve_destroy(circle_b);
}

/* ------------------------------------------------------------------ */
/*  Test: boolean difference                                           */
/* ------------------------------------------------------------------ */
static void test_boolean_difference(void)
{
	qaws_curve* circle_a = NULL;
	qaws_curve* circle_b = NULL;
	qaws_curve* result[4];
	unsigned int result_count = 0;
	qaws_status s;

	printf("test_boolean_difference\n");

	circle_a = make_circle_2d(0, 0, (qaws_scalar)1.0);
	circle_b = make_circle_2d((qaws_scalar)0.8, 0, (qaws_scalar)1.0);

	memset(result, 0, sizeof(result));
	s = qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_DIFFERENCE,
		result, 4, &result_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(result_count >= 1, "difference produced at least 1 boundary");

	{
		unsigned int ri;
		for (ri = 0; ri < result_count; ri++)
			qaws_curve_destroy(result[ri]);
	}
	qaws_curve_destroy(circle_a);
	qaws_curve_destroy(circle_b);
}

/* ------------------------------------------------------------------ */
/*  Test: disjoint circles (no intersection)                           */
/* ------------------------------------------------------------------ */
static void test_boolean_disjoint(void)
{
	qaws_curve* circle_a = NULL;
	qaws_curve* circle_b = NULL;
	qaws_curve* result[4];
	unsigned int result_count = 0;
	qaws_status s;

	printf("test_boolean_disjoint\n");

	circle_a = make_circle_2d(0, 0, (qaws_scalar)0.5);
	circle_b = make_circle_2d((qaws_scalar)3.0, 0, (qaws_scalar)0.5);

	/* Union of disjoint circles: should get 2 boundaries */
	memset(result, 0, sizeof(result));
	s = qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_UNION,
		result, 4, &result_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(result_count == 2, "disjoint union: 2 boundaries");

	{
		unsigned int ri;
		for (ri = 0; ri < result_count; ri++)
			qaws_curve_destroy(result[ri]);
	}

	/* Intersection of disjoint circles: empty */
	result_count = 0;
	memset(result, 0, sizeof(result));
	s = qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_INTERSECTION,
		result, 4, &result_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(result_count == 0, "disjoint intersection: empty");

	qaws_curve_destroy(circle_a);
	qaws_curve_destroy(circle_b);
}

/* ------------------------------------------------------------------ */
/*  Internal: draw a boolean result into the svg at a given offset     */
/* ------------------------------------------------------------------ */
static void draw_bool_result_svg(svg_writer* svg,
	qaws_curve* const* results, unsigned int count,
	char const* color, qaws_scalar x_off, qaws_scalar y_off)
{
	unsigned int ri;
	for (ri = 0; ri < count; ri++)
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n, pi;
		if (!results[ri]) continue;
		n = svg_sample_curve(results[ri], buf, SVG_SAMPLES);
		/* Apply offset */
		for (pi = 0; pi < n; pi++)
		{
			buf[pi].x += x_off;
			buf[pi].y += y_off;
		}
		svg_polyline(svg, buf, n, color, (qaws_scalar)2.0);
	}
}

/* ------------------------------------------------------------------ */
/*  Internal: draw two input circles as grey outlines with offset      */
/* ------------------------------------------------------------------ */
static void draw_input_circles(svg_writer* svg,
	qaws_curve* ca, qaws_curve* cb,
	qaws_scalar x_off, qaws_scalar y_off)
{
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n, pi;

	n = svg_sample_curve(ca, buf, SVG_SAMPLES);
	for (pi = 0; pi < n; pi++) { buf[pi].x += x_off; buf[pi].y += y_off; }
	svg_polyline(svg, buf, n, "#555555", (qaws_scalar)1.0);

	n = svg_sample_curve(cb, buf, SVG_SAMPLES);
	for (pi = 0; pi < n; pi++) { buf[pi].x += x_off; buf[pi].y += y_off; }
	svg_polyline(svg, buf, n, "#555555", (qaws_scalar)1.0);
}

/* ------------------------------------------------------------------ */
/*  Visual: SVG boolean operations (union, intersection, difference)   */
/* ------------------------------------------------------------------ */
static void visual_svg_boolean(void)
{
	svg_writer svg;
	qaws_curve* circle_a = NULL;
	qaws_curve* circle_b = NULL;
	qaws_curve* result[4];
	unsigned int result_count = 0;
	unsigned int ri;

	printf("visual_svg_boolean\n");
	svg_ensure_output_dir();

	circle_a = make_circle_2d(0, 0, (qaws_scalar)1.0);
	circle_b = make_circle_2d((qaws_scalar)0.8, 0, (qaws_scalar)1.0);
	if (!circle_a || !circle_b) goto cleanup_vis;

	if (!svg_open(&svg, SVG_OUTPUT_DIR "/43_boolean_2d.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.5, (qaws_scalar)11.0, (qaws_scalar)3.0,
		(qaws_scalar)900, (qaws_scalar)250))
		goto cleanup_vis;

	/* Column 1: Union */
	draw_input_circles(&svg, circle_a, circle_b, 0, 0);
	memset(result, 0, sizeof(result));
	result_count = 0;
	qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_UNION, result, 4, &result_count);
	draw_bool_result_svg(&svg, result, result_count, "#50fa7b", 0, 0);
	for (ri = 0; ri < result_count; ri++) qaws_curve_destroy(result[ri]);
	svg_label(&svg, (qaws_scalar)-0.2, (qaws_scalar)-1.3, "Union", "#50fa7b");

	/* Column 2: Intersection */
	draw_input_circles(&svg, circle_a, circle_b, (qaws_scalar)3.5, 0);
	memset(result, 0, sizeof(result));
	result_count = 0;
	qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_INTERSECTION, result, 4, &result_count);
	draw_bool_result_svg(&svg, result, result_count, "#e94560", (qaws_scalar)3.5, 0);
	for (ri = 0; ri < result_count; ri++) qaws_curve_destroy(result[ri]);
	svg_label(&svg, (qaws_scalar)3.3, (qaws_scalar)-1.3, "Intersection", "#e94560");

	/* Column 3: Difference */
	draw_input_circles(&svg, circle_a, circle_b, (qaws_scalar)7.0, 0);
	memset(result, 0, sizeof(result));
	result_count = 0;
	qaws_boolean_2d(circle_a, circle_b, QAWS_BOOLEAN_DIFFERENCE, result, 4, &result_count);
	draw_bool_result_svg(&svg, result, result_count, "#6272a4", (qaws_scalar)7.0, 0);
	for (ri = 0; ri < result_count; ri++) qaws_curve_destroy(result[ri]);
	svg_label(&svg, (qaws_scalar)6.8, (qaws_scalar)-1.3, "Difference", "#6272a4");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/43_boolean_2d.svg\n");

cleanup_vis:
	qaws_curve_destroy(circle_a);
	qaws_curve_destroy(circle_b);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_43_boolean_2d_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_boolean_union();
	test_boolean_intersection();
	test_boolean_difference();
	test_boolean_disjoint();

	/* Visual output */
	visual_svg_boolean();

	printf("43_boolean_2d: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
