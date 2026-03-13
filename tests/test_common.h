/*
 * test_common.h - Shared test infrastructure for qaws tests
 *
 * Include guard: QAWS_TEST_COMMON_H
 * All functions are declared static for inclusion in multiple translation units.
 */

#ifndef QAWS_TEST_COMMON_H
#define QAWS_TEST_COMMON_H

#include "qaws.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path) mkdir(path, 0755)
#endif

#define SVG_OUTPUT_DIR "tests_output"
#define OBJ_OUTPUT_DIR "tests_output"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SVG_SAMPLES 256

static int g_pass = 0;
static int g_fail = 0;

#define TEST_ASSERT(cond, msg) \
	do { \
		if (cond) { g_pass++; } \
		else { g_fail++; printf("  FAIL: %s (line %d)\n", msg, __LINE__); } \
	} while (0)

#define TEST_ASSERT_STATUS(status) \
	TEST_ASSERT((status) == QAWS_STATUS_OK, #status " == OK")

#if QAWS_SCALAR_IS_FLOAT
#define TOLERANCE 1e-4f
#define TOLERANCE_LOOSE 1e-3f
#else
#define TOLERANCE 1e-8
#define TOLERANCE_LOOSE 1e-4
#endif

static int approx_eq_loose(qaws_scalar a, qaws_scalar b)
{
	qaws_scalar diff = a - b;
	if (diff < 0) diff = -diff;
	return diff < TOLERANCE_LOOSE;
}

static int approx_eq(qaws_scalar a, qaws_scalar b)
{
	qaws_scalar diff = a - b;
	if (diff < 0) diff = -diff;
	return diff < TOLERANCE;
}

/* ================================================================== */
/* Trivial SVG writer                                                  */
/* ================================================================== */

typedef struct svg_writer
{
	FILE* fp;
	qaws_scalar view_x, view_y, view_w, view_h;
	qaws_scalar svg_w, svg_h;
	qaws_scalar padding;
} svg_writer;

static int svg_open(svg_writer* w, char const* path,
	qaws_scalar view_x, qaws_scalar view_y,
	qaws_scalar view_w, qaws_scalar view_h,
	qaws_scalar svg_w, qaws_scalar svg_h)
{
	w->fp = fopen(path, "w");
	if (!w->fp) return 0;
	w->view_x = view_x; w->view_y = view_y;
	w->view_w = view_w; w->view_h = view_h;
	w->svg_w = svg_w;   w->svg_h = svg_h;
	w->padding = (qaws_scalar)20;
	fprintf(w->fp,
		"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		"<svg xmlns=\"http://www.w3.org/2000/svg\" "
		"width=\"%.0f\" height=\"%.0f\" "
		"viewBox=\"0 0 %.0f %.0f\" "
		"style=\"background:#1a1a2e\">\n",
		(double)svg_w, (double)svg_h, (double)svg_w, (double)svg_h);
	/* Grid */
	fprintf(w->fp,
		"<defs><pattern id=\"grid\" width=\"%.4f\" height=\"%.4f\" patternUnits=\"userSpaceOnUse\">"
		"<path d=\"M %.4f 0 L 0 0 0 %.4f\" fill=\"none\" stroke=\"#2a2a4a\" stroke-width=\"0.5\"/>"
		"</pattern></defs>\n"
		"<rect width=\"100%%\" height=\"100%%\" fill=\"url(#grid)\"/>\n",
		(double)((svg_w - 2 * w->padding) / view_w),
		(double)((svg_h - 2 * w->padding) / view_h),
		(double)((svg_w - 2 * w->padding) / view_w),
		(double)((svg_h - 2 * w->padding) / view_h));
	/* Axes */
	{
		qaws_scalar scale_x = (svg_w - 2 * w->padding) / view_w;
		qaws_scalar scale_y = (svg_h - 2 * w->padding) / view_h;
		qaws_scalar ox = w->padding + (0 - view_x) * scale_x;
		qaws_scalar oy = svg_h - w->padding - (0 - view_y) * scale_y;
		fprintf(w->fp,
			"<line x1=\"%.2f\" y1=\"0\" x2=\"%.2f\" y2=\"%.0f\" stroke=\"#3a3a5a\" stroke-width=\"1\"/>\n"
			"<line x1=\"0\" y1=\"%.2f\" x2=\"%.0f\" y2=\"%.2f\" stroke=\"#3a3a5a\" stroke-width=\"1\"/>\n",
			(double)ox, (double)ox, (double)svg_h,
			(double)oy, (double)svg_w, (double)oy);
	}
	return 1;
}

static void svg_map(svg_writer const* w, qaws_scalar wx, qaws_scalar wy,
	qaws_scalar* sx, qaws_scalar* sy)
{
	qaws_scalar scale_x = (w->svg_w - 2 * w->padding) / w->view_w;
	qaws_scalar scale_y = (w->svg_h - 2 * w->padding) / w->view_h;
	*sx = w->padding + (wx - w->view_x) * scale_x;
	*sy = w->svg_h - w->padding - (wy - w->view_y) * scale_y;
}

static void svg_polyline(svg_writer* w, qaws_vec2 const* pts, unsigned int count,
	char const* stroke, qaws_scalar stroke_width)
{
	if (count < 2) return;
	fprintf(w->fp, "<polyline points=\"");
	for (unsigned int i = 0; i < count; i++) {
		qaws_scalar sx, sy;
		svg_map(w, pts[i].x, pts[i].y, &sx, &sy);
		fprintf(w->fp, "%.4f,%.4f ", (double)sx, (double)sy);
	}
	fprintf(w->fp, "\" fill=\"none\" stroke=\"%s\" stroke-width=\"%.1f\" "
		"stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n",
		stroke, (double)stroke_width);
}

static void svg_dots(svg_writer* w, qaws_vec2 const* pts, unsigned int count,
	char const* fill, qaws_scalar radius)
{
	for (unsigned int i = 0; i < count; i++) {
		qaws_scalar sx, sy;
		svg_map(w, pts[i].x, pts[i].y, &sx, &sy);
		fprintf(w->fp, "<circle cx=\"%.4f\" cy=\"%.4f\" r=\"%.1f\" fill=\"%s\"/>\n",
			(double)sx, (double)sy, (double)radius, fill);
	}
}

static void svg_label(svg_writer* w, qaws_scalar wx, qaws_scalar wy,
	char const* text, char const* fill)
{
	qaws_scalar sx, sy;
	svg_map(w, wx, wy, &sx, &sy);
	fprintf(w->fp,
		"<text x=\"%.2f\" y=\"%.2f\" font-family=\"monospace\" font-size=\"13\" "
		"fill=\"%s\">%s</text>\n",
		(double)sx, (double)(sy - 8), fill, text);
}

static void svg_close(svg_writer* w)
{
	fprintf(w->fp, "</svg>\n");
	fclose(w->fp);
	w->fp = NULL;
}

static void svg_line(svg_writer* w, qaws_scalar x0, qaws_scalar y0,
	qaws_scalar x1, qaws_scalar y1, char const* stroke, qaws_scalar width)
{
	qaws_scalar sx0, sy0, sx1, sy1;
	svg_map(w, x0, y0, &sx0, &sy0);
	svg_map(w, x1, y1, &sx1, &sy1);
	fprintf(w->fp, "<line x1=\"%.4f\" y1=\"%.4f\" x2=\"%.4f\" y2=\"%.4f\" "
		"stroke=\"%s\" stroke-width=\"%.1f\"/>\n",
		(double)sx0, (double)sy0, (double)sx1, (double)sy1,
		stroke, (double)width);
}

static void svg_polyline_dashed(svg_writer* w, qaws_vec2 const* pts, unsigned int count,
	char const* stroke, qaws_scalar stroke_width, char const* dasharray)
{
	if (count < 2) return;
	fprintf(w->fp, "<polyline points=\"");
	for (unsigned int i = 0; i < count; i++) {
		qaws_scalar sx, sy;
		svg_map(w, pts[i].x, pts[i].y, &sx, &sy);
		fprintf(w->fp, "%.4f,%.4f ", (double)sx, (double)sy);
	}
	fprintf(w->fp, "\" fill=\"none\" stroke=\"%s\" stroke-width=\"%.1f\" "
		"stroke-dasharray=\"%s\" "
		"stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n",
		stroke, (double)stroke_width, dasharray);
}

static void svg_ring(svg_writer* w, qaws_scalar wx, qaws_scalar wy,
	char const* stroke, qaws_scalar r, qaws_scalar stroke_width)
{
	qaws_scalar sx, sy;
	svg_map(w, wx, wy, &sx, &sy);
	fprintf(w->fp, "<circle cx=\"%.4f\" cy=\"%.4f\" r=\"%.1f\" fill=\"none\" "
		"stroke=\"%s\" stroke-width=\"%.1f\"/>\n",
		(double)sx, (double)sy, (double)r, stroke, (double)stroke_width);
}

static void svg_dot(svg_writer* w, qaws_scalar wx, qaws_scalar wy,
	char const* fill, qaws_scalar r)
{
	qaws_scalar sx, sy;
	svg_map(w, wx, wy, &sx, &sy);
	fprintf(w->fp, "<circle cx=\"%.4f\" cy=\"%.4f\" r=\"%.1f\" fill=\"%s\"/>\n",
		(double)sx, (double)sy, (double)r, fill);
}

/* Helper: sample a 2D curve into a buffer */
static unsigned int svg_sample_curve(qaws_curve const* curve,
	qaws_vec2* buf, unsigned int capacity)
{
	qaws_sampling_desc sd;
	sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
	sd.sample_count = capacity;
	sd.error_tolerance = 0;
	sd.use_adaptive_sampling = 0;
	unsigned int count = 0;
	qaws_curve_sample_2d(curve, &sd, buf, capacity, &count);
	return count;
}

static void svg_ensure_output_dir(void)
{
	//MKDIR("working_dir");
	MKDIR(SVG_OUTPUT_DIR);
}

static void svg_arrow(svg_writer* w, qaws_curve const* crv, qaws_scalar frac,
	char const* color, qaws_scalar size)
{
	/* Draw small direction arrow at parameter fraction along curve */
	qaws_eval_result_2d r;
	qaws_range rng = qaws_curve_get_parameter_range(crv);
	qaws_scalar t = rng.min_value + (rng.max_value - rng.min_value) * frac;
	qaws_scalar len, nx, ny;
	qaws_curve_evaluate_2d(crv, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
	len = (qaws_scalar)sqrt((double)(r.d1.x * r.d1.x + r.d1.y * r.d1.y));
	if (len < (qaws_scalar)1e-10) return;
	nx = r.d1.x / len; ny = r.d1.y / len;
	/* Arrowhead: two lines from tip backward */
	svg_line(w, r.position.x, r.position.y,
		r.position.x - nx * size + ny * size * (qaws_scalar)0.4,
		r.position.y - ny * size - nx * size * (qaws_scalar)0.4,
		color, (qaws_scalar)1.5);
	svg_line(w, r.position.x, r.position.y,
		r.position.x - nx * size - ny * size * (qaws_scalar)0.4,
		r.position.y - ny * size + nx * size * (qaws_scalar)0.4,
		color, (qaws_scalar)1.5);
}

static void svg_draw_extrema(svg_writer* svg, qaws_curve* curve, qaws_scalar y_lo, qaws_scalar y_hi,
	char const* curve_color, char const* label)
{
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(svg, buf, n, curve_color, (qaws_scalar)2.5);

	/* X-extrema: blue dots + vertical lines */
	{
		qaws_scalar params[8];
		unsigned int count = 0, i;
		if (qaws_curve_find_extrema(curve, 0, params, 8, &count) == QAWS_STATUS_OK) {
			for (i = 0; i < count; i++) {
				qaws_eval_result_2d r;
				qaws_curve_evaluate_2d(curve, params[i], QAWS_EVAL_FLAG_POSITION, &r);
				svg_dot(svg, r.position.x, r.position.y, "#6272a4", (qaws_scalar)5);
				svg_line(svg, r.position.x, y_lo, r.position.x, y_hi,
					"#6272a450", (qaws_scalar)0.5);
			}
		}
	}
	/* Y-extrema: green dots + horizontal lines */
	{
		qaws_scalar params[8];
		unsigned int count = 0;
		qaws_range rng = qaws_curve_get_parameter_range(curve);
		qaws_eval_result_2d r0, r1;
		qaws_scalar x0, x1;
		qaws_curve_evaluate_2d(curve, rng.min_value, QAWS_EVAL_FLAG_POSITION, &r0);
		qaws_curve_evaluate_2d(curve, rng.max_value, QAWS_EVAL_FLAG_POSITION, &r1);
		x0 = r0.position.x < r1.position.x ? r0.position.x - (qaws_scalar)0.3 : r1.position.x - (qaws_scalar)0.3;
		x1 = r0.position.x > r1.position.x ? r0.position.x + (qaws_scalar)0.3 : r1.position.x + (qaws_scalar)0.3;
		if (qaws_curve_find_extrema(curve, 1, params, 8, &count) == QAWS_STATUS_OK) {
			unsigned int yi;
			for (yi = 0; yi < count; yi++) {
				qaws_eval_result_2d r;
				qaws_curve_evaluate_2d(curve, params[yi], QAWS_EVAL_FLAG_POSITION, &r);
				svg_dot(svg, r.position.x, r.position.y, "#50fa7b", (qaws_scalar)5);
				svg_line(svg, x0, r.position.y, x1, r.position.y,
					"#50fa7b50", (qaws_scalar)0.5);
			}
		}
	}
	svg_label(svg, (qaws_scalar)-0.3, y_hi - (qaws_scalar)0.3, label, "#8888aa");
}

/* ================================================================== */
/* OBJ writer helpers for 3D visualization                             */
/* ================================================================== */

typedef struct obj_writer
{
	FILE* fp;
	FILE* mtl_fp;
	unsigned int vertex_count;   /* running total, OBJ is 1-indexed */
	unsigned int normal_count;
	char mtl_name[256];
} obj_writer;

static int obj_open(obj_writer* w, char const* obj_path, char const* mtl_path)
{
	w->fp = fopen(obj_path, "w");
	if (!w->fp) return 0;
	w->mtl_fp = NULL;
	w->vertex_count = 0;
	w->normal_count = 0;
	if (mtl_path)
	{
		char const* slash;
		char const* bslash;
		char const* last;
		w->mtl_fp = fopen(mtl_path, "w");
		if (!w->mtl_fp) { fclose(w->fp); return 0; }
		/* Write mtllib with just the filename */
		slash = strrchr(mtl_path, '/');
		bslash = strrchr(mtl_path, '\\');
		last = slash > bslash ? slash : bslash;
		if (last) last++; else last = mtl_path;
		strncpy(w->mtl_name, last, sizeof(w->mtl_name) - 1);
		w->mtl_name[sizeof(w->mtl_name) - 1] = '\0';
		fprintf(w->fp, "mtllib %s\n", w->mtl_name);
	}
	fprintf(w->fp, "# Qaws 3D test output\n\n");
	return 1;
}

static void obj_close(obj_writer* w)
{
	if (w->fp) { fclose(w->fp); w->fp = NULL; }
	if (w->mtl_fp) { fclose(w->mtl_fp); w->mtl_fp = NULL; }
}

static void obj_material(obj_writer* w, char const* name,
	double r, double g, double b)
{
	if (!w->mtl_fp) return;
	fprintf(w->mtl_fp,
		"newmtl %s\n"
		"Ka %.3f %.3f %.3f\n"
		"Kd %.3f %.3f %.3f\n"
		"Ks 0.3 0.3 0.3\n"
		"Ns 50.0\n"
		"d 1.0\n\n",
		name, r * 0.2, g * 0.2, b * 0.2, r, g, b);
}

static void obj_group(obj_writer* w, char const* name)
{
	fprintf(w->fp, "\ng %s\n", name);
}

static void obj_use_material(obj_writer* w, char const* name)
{
	fprintf(w->fp, "usemtl %s\n", name);
}

static void obj_comment(obj_writer* w, char const* text)
{
	fprintf(w->fp, "# %s\n", text);
}

/* Write a single vertex, return its 1-based index */
static unsigned int obj_vertex(obj_writer* w, qaws_vec3 p)
{
	fprintf(w->fp, "v %.6f %.6f %.6f\n", (double)p.x, (double)p.y, (double)p.z);
	return ++w->vertex_count;
}

/* Write a vertex normal, return its 1-based index */
static unsigned int obj_normal(obj_writer* w, qaws_vec3 n)
{
	fprintf(w->fp, "vn %.6f %.6f %.6f\n", (double)n.x, (double)n.y, (double)n.z);
	return ++w->normal_count;
}

/* Sample a 3D curve and write as polyline. Returns first vertex index. */
static unsigned int obj_curve_polyline(obj_writer* w, qaws_curve const* curve,
	unsigned int sample_count)
{
	unsigned int first_vi, count, ci;
	qaws_vec3 buf[512];
	qaws_sampling_desc sd;

	if (sample_count > 512) sample_count = 512;
	memset(&sd, 0, sizeof(sd));
	sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
	sd.sample_count = sample_count;
	count = 0;
	qaws_curve_sample_3d(curve, &sd, buf, 512, &count);
	if (count < 2) return 0;

	first_vi = w->vertex_count + 1;
	for (ci = 0; ci < count; ci++)
		obj_vertex(w, buf[ci]);

	/* Write line segments (pairs, not one long polyline — better viewer compat) */
	for (ci = 0; ci + 1 < count; ci++)
		fprintf(w->fp, "l %u %u\n", first_vi + ci, first_vi + ci + 1);
	return first_vi;
}

/* Write a UV sphere at a point (8 longitude x 6 latitude segments) */
static void obj_sphere(obj_writer* w, qaws_vec3 center, qaws_scalar radius)
{
	unsigned int slices = 8, stacks = 6;
	unsigned int first_vi = w->vertex_count + 1;
	unsigned int si, ti;
	qaws_vec3 v;

	/* Top pole */
	v.x = center.x; v.y = center.y + radius; v.z = center.z;
	obj_vertex(w, v);
	/* Body vertices: (stacks-1) rings of (slices) vertices */
	for (ti = 1; ti < stacks; ti++)
	{
		double phi = M_PI * (double)ti / (double)stacks;
		double sp = sin(phi), cp = cos(phi);
		for (si = 0; si < slices; si++)
		{
			double theta = 2.0 * M_PI * (double)si / (double)slices;
			v.x = center.x + radius * (qaws_scalar)(sp * cos(theta));
			v.y = center.y + radius * (qaws_scalar)(cp);
			v.z = center.z + radius * (qaws_scalar)(sp * sin(theta));
			obj_vertex(w, v);
		}
	}
	/* Bottom pole */
	v.x = center.x; v.y = center.y - radius; v.z = center.z;
	obj_vertex(w, v);

	/* Top cap: triangles from pole to first ring */
	for (si = 0; si < slices; si++)
	{
		unsigned int next = (si + 1) % slices;
		fprintf(w->fp, "f %u %u %u\n",
			first_vi, first_vi + 1 + si, first_vi + 1 + next);
	}
	/* Body: quad strips between adjacent rings */
	for (ti = 0; ti < stacks - 2; ti++)
	{
		unsigned int ring0 = first_vi + 1 + ti * slices;
		unsigned int ring1 = first_vi + 1 + (ti + 1) * slices;
		for (si = 0; si < slices; si++)
		{
			unsigned int next = (si + 1) % slices;
			fprintf(w->fp, "f %u %u %u %u\n",
				ring0 + si, ring0 + next, ring1 + next, ring1 + si);
		}
	}
	/* Bottom cap: triangles from last ring to bottom pole */
	{
		unsigned int last_ring = first_vi + 1 + (stacks - 2) * slices;
		unsigned int bot_pole = first_vi + 1 + (stacks - 1) * slices;
		for (si = 0; si < slices; si++)
		{
			unsigned int next = (si + 1) % slices;
			fprintf(w->fp, "f %u %u %u\n",
				last_ring + si, bot_pole, last_ring + next);
		}
	}
}

/* Backward-compatible alias */
static void obj_marker(obj_writer* w, qaws_vec3 center, qaws_scalar size)
{
	obj_sphere(w, center, size);
}

/* Write a line segment between two points */
static void obj_line_segment(obj_writer* w, qaws_vec3 a, qaws_vec3 b)
{
	unsigned int va = obj_vertex(w, a);
	unsigned int vb = obj_vertex(w, b);
	fprintf(w->fp, "l %u %u\n", va, vb);
}

/* Write a Frenet frame as 3 line segments: T(red-ish), N(green-ish), B(blue-ish) */
static void obj_frenet_frame(obj_writer* w, qaws_vec3 pos,
	qaws_vec3 T, qaws_vec3 N, qaws_vec3 B, qaws_scalar length)
{
	qaws_vec3 tip;
	tip.x = pos.x + T.x * length;
	tip.y = pos.y + T.y * length;
	tip.z = pos.z + T.z * length;
	obj_line_segment(w, pos, tip);

	tip.x = pos.x + N.x * length;
	tip.y = pos.y + N.y * length;
	tip.z = pos.z + N.z * length;
	obj_line_segment(w, pos, tip);

	tip.x = pos.x + B.x * length;
	tip.y = pos.y + B.y * length;
	tip.z = pos.z + B.z * length;
	obj_line_segment(w, pos, tip);
}

/* Helper: normalize a vec3 in-place, return length */
static double v3_normalize(qaws_vec3* v)
{
	double len = sqrt((double)(v->x*v->x + v->y*v->y + v->z*v->z));
	if (len > 1e-12) { v->x /= (qaws_scalar)len; v->y /= (qaws_scalar)len; v->z /= (qaws_scalar)len; }
	return len;
}

static double v3_dot(qaws_vec3 a, qaws_vec3 b)
{
	return (double)a.x*(double)b.x + (double)a.y*(double)b.y + (double)a.z*(double)b.z;
}

static qaws_vec3 v3_cross(qaws_vec3 a, qaws_vec3 b)
{
	qaws_vec3 r;
	r.x = a.y*b.z - a.z*b.y;
	r.y = a.z*b.x - a.x*b.z;
	r.z = a.x*b.y - a.y*b.x;
	return r;
}

/* Rotate vector v around axis by angle (Rodrigues' formula) */
static qaws_vec3 v3_rotate(qaws_vec3 v, qaws_vec3 axis, double cos_a, double sin_a)
{
	/* v*cos + (axis x v)*sin + axis*(axis.v)*(1-cos) */
	double adv = v3_dot(axis, v);
	qaws_vec3 axv = v3_cross(axis, v);
	qaws_vec3 r;
	r.x = (qaws_scalar)((double)v.x * cos_a + (double)axv.x * sin_a + (double)axis.x * adv * (1.0 - cos_a));
	r.y = (qaws_scalar)((double)v.y * cos_a + (double)axv.y * sin_a + (double)axis.y * adv * (1.0 - cos_a));
	r.z = (qaws_scalar)((double)v.z * cos_a + (double)axv.z * sin_a + (double)axis.z * adv * (1.0 - cos_a));
	return r;
}

/* Build an initial perpendicular frame for a tangent vector */
static void build_initial_frame(qaws_vec3 T_vec, qaws_vec3* N_out, qaws_vec3* B_out)
{
	qaws_vec3 up = { 0, 1, 0 };
	*N_out = v3_cross(up, T_vec);
	if (v3_normalize(N_out) < 1e-8)
	{
		up.x = 1; up.y = 0; up.z = 0;
		*N_out = v3_cross(up, T_vec);
		v3_normalize(N_out);
	}
	*B_out = v3_cross(T_vec, *N_out);
	v3_normalize(B_out);
}

/* Generate a tube mesh along a 3D curve using parallel transport frame
   (Bishop frame) to avoid Frenet flipping artifacts */
static void obj_tube(obj_writer* w, qaws_curve const* curve,
	unsigned int length_segments, unsigned int circle_segments,
	qaws_scalar radius)
{
	unsigned int si, ci;
	unsigned int first_vi = w->vertex_count + 1;
	unsigned int first_ni = w->normal_count + 1;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	qaws_scalar t_min = range.min_value;
	qaws_scalar t_len = range.max_value - range.min_value;
	qaws_vec3 prev_T = {0,0,0}, prev_N = {0,0,0}, prev_B = {0,0,0};
	int frame_initialized = 0;
	unsigned int rings_emitted = 0;

	for (si = 0; si <= length_segments; si++)
	{
		qaws_scalar t = t_min + t_len * (qaws_scalar)si / (qaws_scalar)length_segments;
		qaws_vec3 T_vec, N_vec, B_vec;
		qaws_eval_result_3d r;

		if (qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION, &r) != QAWS_STATUS_OK)
			continue;
		if (qaws_curve_compute_tangent_3d(curve, t, &T_vec) != QAWS_STATUS_OK)
			continue;
		v3_normalize(&T_vec);

		/* Ensure consistent tangent direction: if tangent flips vs previous,
		   negate it. This handles span boundary discontinuities. */
		if (frame_initialized && v3_dot(prev_T, T_vec) < 0.0)
		{
			T_vec.x = -T_vec.x;
			T_vec.y = -T_vec.y;
			T_vec.z = -T_vec.z;
		}

		if (!frame_initialized)
		{
			/* First point: build arbitrary perpendicular frame */
			build_initial_frame(T_vec, &N_vec, &B_vec);
			frame_initialized = 1;
		}
		else
		{
			/* Parallel transport: rotate previous N,B to align with new tangent.
			   Rotation axis = prev_T x new_T, angle = acos(prev_T . new_T) */
			double d = v3_dot(prev_T, T_vec);
			if (d > 1.0) d = 1.0;
			if (d < -1.0) d = -1.0;

			if (d > 0.9999999)
			{
				/* Nearly parallel: keep previous frame */
				N_vec = prev_N;
				B_vec = prev_B;
			}
			else
			{
				qaws_vec3 rot_axis = v3_cross(prev_T, T_vec);
				double sin_a, cos_a;
				v3_normalize(&rot_axis);
				cos_a = d;
				sin_a = sqrt(1.0 - d * d);
				N_vec = v3_rotate(prev_N, rot_axis, cos_a, sin_a);
				B_vec = v3_rotate(prev_B, rot_axis, cos_a, sin_a);
				v3_normalize(&N_vec);
				v3_normalize(&B_vec);
			}
		}

		prev_T = T_vec;
		prev_N = N_vec;
		prev_B = B_vec;

		for (ci = 0; ci < circle_segments; ci++)
		{
			double angle = 2.0 * M_PI * (double)ci / (double)circle_segments;
			double ca = cos(angle), sa = sin(angle);
			qaws_vec3 pos, nrm;
			pos.x = r.position.x + radius * (qaws_scalar)(ca * (double)N_vec.x + sa * (double)B_vec.x);
			pos.y = r.position.y + radius * (qaws_scalar)(ca * (double)N_vec.y + sa * (double)B_vec.y);
			pos.z = r.position.z + radius * (qaws_scalar)(ca * (double)N_vec.z + sa * (double)B_vec.z);
			nrm.x = (qaws_scalar)(ca * (double)N_vec.x + sa * (double)B_vec.x);
			nrm.y = (qaws_scalar)(ca * (double)N_vec.y + sa * (double)B_vec.y);
			nrm.z = (qaws_scalar)(ca * (double)N_vec.z + sa * (double)B_vec.z);
			obj_vertex(w, pos);
			obj_normal(w, nrm);
		}
		rings_emitted++;
	}

	/* Generate quad faces connecting adjacent rings */
	if (rings_emitted > 1)
	{
		for (si = 0; si < rings_emitted - 1; si++)
		{
			for (ci = 0; ci < circle_segments; ci++)
			{
				unsigned int ci_next = (ci + 1) % circle_segments;
				unsigned int v00 = first_vi + si * circle_segments + ci;
				unsigned int v01 = first_vi + si * circle_segments + ci_next;
				unsigned int v10 = first_vi + (si + 1) * circle_segments + ci;
				unsigned int v11 = first_vi + (si + 1) * circle_segments + ci_next;
				unsigned int n00 = first_ni + si * circle_segments + ci;
				unsigned int n01 = first_ni + si * circle_segments + ci_next;
				unsigned int n10 = first_ni + (si + 1) * circle_segments + ci;
				unsigned int n11 = first_ni + (si + 1) * circle_segments + ci_next;
				fprintf(w->fp, "f %u//%u %u//%u %u//%u %u//%u\n",
					v00, n00, v01, n01, v11, n11, v10, n10);
			}
		}
	}
}

/* Write a curvature comb for a 3D curve as lines from curve to comb tips */
static void obj_curvature_comb(obj_writer* w, qaws_curve const* curve,
	unsigned int sample_count, qaws_scalar scale)
{
	qaws_curvature_sample_3d samples[256];
	unsigned int ci;
	if (sample_count > 256) sample_count = 256;
	if (qaws_curve_compute_curvature_comb_3d(curve, sample_count, samples, 256)
		!= QAWS_STATUS_OK)
		return;
	for (ci = 0; ci < sample_count; ci++)
	{
		qaws_vec3 tip;
		tip.x = samples[ci].position.x + samples[ci].normal.x * samples[ci].curvature * scale;
		tip.y = samples[ci].position.y + samples[ci].normal.y * samples[ci].curvature * scale;
		tip.z = samples[ci].position.z + samples[ci].normal.z * samples[ci].curvature * scale;
		obj_line_segment(w, samples[ci].position, tip);
	}
}

/* Write a curvature comb as a visible triangle-strip ribbon surface
   between the curve and the comb tips.  Each comb spine is also drawn
   as a chain of small spheres so it renders in all OBJ viewers. */
static void obj_curvature_comb_surface(obj_writer* w, qaws_curve const* curve,
	unsigned int sample_count, qaws_scalar scale)
{
	qaws_curvature_sample_3d samples[256];
	unsigned int first_v, ci;
	if (sample_count > 256) sample_count = 256;
	if (qaws_curve_compute_curvature_comb_3d(curve, sample_count, samples, 256)
		!= QAWS_STATUS_OK)
		return;

	/* Emit pairs of vertices: curve point and comb tip */
	first_v = w->vertex_count + 1;
	for (ci = 0; ci < sample_count; ci++)
	{
		qaws_vec3 tip;
		tip.x = samples[ci].position.x + samples[ci].normal.x * samples[ci].curvature * scale;
		tip.y = samples[ci].position.y + samples[ci].normal.y * samples[ci].curvature * scale;
		tip.z = samples[ci].position.z + samples[ci].normal.z * samples[ci].curvature * scale;
		obj_vertex(w, samples[ci].position); /* even index */
		obj_vertex(w, tip);                  /* odd  index */
	}

	/* Triangle strip: two triangles per segment connecting consecutive spines */
	for (ci = 0; ci + 1 < sample_count; ci++)
	{
		unsigned int base = first_v + ci * 2;
		unsigned int next = base + 2;
		/* base+0 = curve[ci], base+1 = tip[ci]
		   next+0 = curve[ci+1], next+1 = tip[ci+1] */
		fprintf(w->fp, "f %u %u %u\n", base, base + 1, next + 1);
		fprintf(w->fp, "f %u %u %u\n", base, next + 1, next);
	}

	/* Comb spine spheres (every 4th sample) for visibility */
	for (ci = 0; ci < sample_count; ci += 4)
	{
		qaws_vec3 tip;
		unsigned int si;
		tip.x = samples[ci].position.x + samples[ci].normal.x * samples[ci].curvature * scale;
		tip.y = samples[ci].position.y + samples[ci].normal.y * samples[ci].curvature * scale;
		tip.z = samples[ci].position.z + samples[ci].normal.z * samples[ci].curvature * scale;
		for (si = 0; si <= 4; si++)
		{
			qaws_scalar f = (qaws_scalar)si / (qaws_scalar)4;
			qaws_vec3 p;
			p.x = samples[ci].position.x + f * (tip.x - samples[ci].position.x);
			p.y = samples[ci].position.y + f * (tip.y - samples[ci].position.y);
			p.z = samples[ci].position.z + f * (tip.z - samples[ci].position.z);
			obj_sphere(w, p, (qaws_scalar)0.012);
		}
	}

	/* Ridge line: connect comb tips as a tube-like chain */
	for (ci = 0; ci < sample_count; ci += 2)
	{
		qaws_vec3 tip;
		tip.x = samples[ci].position.x + samples[ci].normal.x * samples[ci].curvature * scale;
		tip.y = samples[ci].position.y + samples[ci].normal.y * samples[ci].curvature * scale;
		tip.z = samples[ci].position.z + samples[ci].normal.z * samples[ci].curvature * scale;
		obj_sphere(w, tip, (qaws_scalar)0.015);
	}
}

/* Write a triangle mesh for a sampled surface grid.
   Auto-detects normal orientation using two heuristics:
   1) V-closure: if v=0 and v=1 map to the same position (tube/pipe),
      du x dv points inward, so flip.
   2) Centroid: if the normal at domain center points toward the centroid
      of the four corners, the normal points into the concavity, so flip.
   Both vertex normals and face winding are flipped together. */
static void obj_surface_mesh(
	obj_writer* w, qaws_surface const* surf,
	unsigned int u_samples, unsigned int v_samples)
{
	qaws_range ur = qaws_surface_get_u_range(surf);
	qaws_range vr = qaws_surface_get_v_range(surf);
	qaws_scalar u_min = ur.min_value, u_len = ur.max_value - ur.min_value;
	qaws_scalar v_min = vr.min_value, v_len = vr.max_value - vr.min_value;
	unsigned int first_v = w->vertex_count + 1;
	unsigned int first_n = w->normal_count + 1;
	unsigned int ui, vi;
	int flip = 0;

	{
		qaws_surface_eval_result rc, r00, r10, r01, r11, rv0, rv1;
		qaws_vec3 centroid;
		qaws_scalar dx, dy, dz, dot_val, gap, diag;

		memset(&rc, 0, sizeof(rc));
		memset(&r00, 0, sizeof(r00));
		memset(&r10, 0, sizeof(r10));
		memset(&r01, 0, sizeof(r01));
		memset(&r11, 0, sizeof(r11));
		memset(&rv0, 0, sizeof(rv0));
		memset(&rv1, 0, sizeof(rv1));

		qaws_surface_evaluate(surf,
			u_min + u_len * (qaws_scalar)0.5,
			v_min + v_len * (qaws_scalar)0.5,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &rc);
		qaws_surface_evaluate(surf, u_min, v_min,
			QAWS_SURFACE_EVAL_POSITION, &r00);
		qaws_surface_evaluate(surf, u_min + u_len, v_min,
			QAWS_SURFACE_EVAL_POSITION, &r10);
		qaws_surface_evaluate(surf, u_min, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &r01);
		qaws_surface_evaluate(surf, u_min + u_len, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &r11);

		/* Heuristic 1: V-closure detection.
		   If S(u_mid, v_min) ~ S(u_mid, v_max), the surface wraps
		   around in v (pipe, swept). For these, du x dv points inward. */
		qaws_surface_evaluate(surf,
			u_min + u_len * (qaws_scalar)0.5, v_min,
			QAWS_SURFACE_EVAL_POSITION, &rv0);
		qaws_surface_evaluate(surf,
			u_min + u_len * (qaws_scalar)0.5, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &rv1);

		dx = rv1.position.x - rv0.position.x;
		dy = rv1.position.y - rv0.position.y;
		dz = rv1.position.z - rv0.position.z;
		gap = dx * dx + dy * dy + dz * dz;

		dx = r11.position.x - r00.position.x;
		dy = r11.position.y - r00.position.y;
		dz = r11.position.z - r00.position.z;
		diag = dx * dx + dy * dy + dz * dz;

		/* Compare squared distances: gap^2 < (0.01 * diag)^2 = 0.0001 * diag^2 */
		if (gap < diag * (qaws_scalar)0.0001 || gap < (qaws_scalar)1e-10)
		{
			flip = 1;
		}

		/* Heuristic 2 (for open surfaces): centroid check.
		   If the normal at center points toward the corner centroid,
		   the surface faces inward toward its concavity. */
		if (!flip)
		{
			centroid.x = (r00.position.x + r10.position.x + r01.position.x + r11.position.x) * (qaws_scalar)0.25;
			centroid.y = (r00.position.y + r10.position.y + r01.position.y + r11.position.y) * (qaws_scalar)0.25;
			centroid.z = (r00.position.z + r10.position.z + r01.position.z + r11.position.z) * (qaws_scalar)0.25;

			dx = centroid.x - rc.position.x;
			dy = centroid.y - rc.position.y;
			dz = centroid.z - rc.position.z;

			dot_val = rc.normal.x * dx + rc.normal.y * dy + rc.normal.z * dz;

			if (dot_val > (qaws_scalar)0.01)
				flip = 1;
		}
	}

	/* Emit vertices and normals */
	for (ui = 0; ui < u_samples; ui++)
	{
		qaws_scalar u = u_min + u_len * (qaws_scalar)ui / (qaws_scalar)(u_samples - 1);
		for (vi = 0; vi < v_samples; vi++)
		{
			qaws_scalar v = v_min + v_len * (qaws_scalar)vi / (qaws_scalar)(v_samples - 1);
			qaws_surface_eval_result r;
			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surf, u, v,
				QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
			obj_vertex(w, r.position);
			if (flip)
			{
				qaws_vec3 fn;
				fn.x = -r.normal.x;
				fn.y = -r.normal.y;
				fn.z = -r.normal.z;
				obj_normal(w, fn);
			}
			else
			{
				obj_normal(w, r.normal);
			}
		}
	}

	/* Emit quad faces (as two triangles each).
	   Default winding matches du x dv; flipped winding for inverted normals. */
	for (ui = 0; ui < u_samples - 1; ui++)
	{
		for (vi = 0; vi < v_samples - 1; vi++)
		{
			unsigned int v00 = first_v + ui * v_samples + vi;
			unsigned int v01 = v00 + 1;
			unsigned int v10 = v00 + v_samples;
			unsigned int v11 = v10 + 1;
			unsigned int n00 = first_n + ui * v_samples + vi;
			unsigned int n01 = n00 + 1;
			unsigned int n10 = n00 + v_samples;
			unsigned int n11 = n10 + 1;
			if (flip)
			{
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v00, n00, v01, n01, v10, n10);
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v01, n01, v11, n11, v10, n10);
			}
			else
			{
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v00, n00, v10, n10, v01, n01);
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v10, n10, v11, n11, v01, n01);
			}
		}
	}
}

/* Write a single vertex with per-vertex RGB color (OBJ extension), return 1-based index */
static unsigned int obj_vertex_color(obj_writer* w, qaws_vec3 p,
	double r, double g, double b)
{
	fprintf(w->fp, "v %.6f %.6f %.6f %.4f %.4f %.4f\n",
		(double)p.x, (double)p.y, (double)p.z, r, g, b);
	return ++w->vertex_count;
}

/* Map normalized curvature t in [-1,1] to blue-white-red diverging color */
static void curvature_color(double t, double* r, double* g, double* b)
{
	if (t < -1.0) t = -1.0;
	if (t > 1.0) t = 1.0;
	if (t < 0.0)
	{
		*r = 1.0 + t;
		*g = 1.0 + t;
		*b = 1.0;
	}
	else
	{
		*r = 1.0;
		*g = 1.0 - t;
		*b = 1.0 - t;
	}
}

/* Surface mesh with per-vertex Gaussian curvature coloring.
   Uses the OBJ per-vertex color extension (v x y z r g b).
   Blue = negative K (saddle), white = zero, red = positive K (dome). */
static void obj_surface_mesh_curvature(
	obj_writer* w, qaws_surface const* surf,
	unsigned int u_samples, unsigned int v_samples)
{
	qaws_range ur = qaws_surface_get_u_range(surf);
	qaws_range vr = qaws_surface_get_v_range(surf);
	qaws_scalar u_min = ur.min_value, u_len = ur.max_value - ur.min_value;
	qaws_scalar v_min = vr.min_value, v_len = vr.max_value - vr.min_value;
	unsigned int first_v, first_n;
	unsigned int ui, vi, total;
	int flip = 0;
	qaws_scalar* kappa = NULL;
	qaws_scalar max_abs_k = 0;

	total = u_samples * v_samples;
	kappa = (qaws_scalar*)malloc(total * sizeof(qaws_scalar));
	if (!kappa) return;

	/* First pass: compute Gaussian curvature at each sample */
	for (ui = 0; ui < u_samples; ui++)
	{
		qaws_scalar u = u_min + u_len * (qaws_scalar)ui / (qaws_scalar)(u_samples - 1);
		for (vi = 0; vi < v_samples; vi++)
		{
			qaws_scalar v = v_min + v_len * (qaws_scalar)vi / (qaws_scalar)(v_samples - 1);
			qaws_scalar K = 0;
			unsigned int idx = ui * v_samples + vi;
			qaws_surface_compute_gaussian_curvature(surf, u, v, &K);
			kappa[idx] = K;
			if (fabs((double)K) > (double)max_abs_k)
				max_abs_k = (qaws_scalar)fabs((double)K);
		}
	}
	if ((double)max_abs_k < 1e-10) max_abs_k = (qaws_scalar)1.0;

	/* Auto-flip detection (same heuristics as obj_surface_mesh) */
	{
		qaws_surface_eval_result rc, r00, r10, r01, r11, rv0, rv1;
		qaws_vec3 centroid;
		qaws_scalar dx, dy, dz, dot_val, gap, diag;

		memset(&rc, 0, sizeof(rc));
		memset(&r00, 0, sizeof(r00));
		memset(&r10, 0, sizeof(r10));
		memset(&r01, 0, sizeof(r01));
		memset(&r11, 0, sizeof(r11));
		memset(&rv0, 0, sizeof(rv0));
		memset(&rv1, 0, sizeof(rv1));

		qaws_surface_evaluate(surf,
			u_min + u_len * (qaws_scalar)0.5,
			v_min + v_len * (qaws_scalar)0.5,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &rc);
		qaws_surface_evaluate(surf, u_min, v_min,
			QAWS_SURFACE_EVAL_POSITION, &r00);
		qaws_surface_evaluate(surf, u_min + u_len, v_min,
			QAWS_SURFACE_EVAL_POSITION, &r10);
		qaws_surface_evaluate(surf, u_min, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &r01);
		qaws_surface_evaluate(surf, u_min + u_len, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &r11);

		qaws_surface_evaluate(surf,
			u_min + u_len * (qaws_scalar)0.5, v_min,
			QAWS_SURFACE_EVAL_POSITION, &rv0);
		qaws_surface_evaluate(surf,
			u_min + u_len * (qaws_scalar)0.5, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &rv1);

		dx = rv1.position.x - rv0.position.x;
		dy = rv1.position.y - rv0.position.y;
		dz = rv1.position.z - rv0.position.z;
		gap = dx * dx + dy * dy + dz * dz;

		dx = r11.position.x - r00.position.x;
		dy = r11.position.y - r00.position.y;
		dz = r11.position.z - r00.position.z;
		diag = dx * dx + dy * dy + dz * dz;

		if (gap < diag * (qaws_scalar)0.0001 || gap < (qaws_scalar)1e-10)
			flip = 1;

		if (!flip)
		{
			centroid.x = (r00.position.x + r10.position.x + r01.position.x + r11.position.x) * (qaws_scalar)0.25;
			centroid.y = (r00.position.y + r10.position.y + r01.position.y + r11.position.y) * (qaws_scalar)0.25;
			centroid.z = (r00.position.z + r10.position.z + r01.position.z + r11.position.z) * (qaws_scalar)0.25;

			dx = centroid.x - rc.position.x;
			dy = centroid.y - rc.position.y;
			dz = centroid.z - rc.position.z;

			dot_val = rc.normal.x * dx + rc.normal.y * dy + rc.normal.z * dz;
			if (dot_val > (qaws_scalar)0.01)
				flip = 1;
		}
	}

	/* Second pass: write colored vertices and normals */
	first_v = w->vertex_count + 1;
	first_n = w->normal_count + 1;
	for (ui = 0; ui < u_samples; ui++)
	{
		qaws_scalar u = u_min + u_len * (qaws_scalar)ui / (qaws_scalar)(u_samples - 1);
		for (vi = 0; vi < v_samples; vi++)
		{
			qaws_scalar v = v_min + v_len * (qaws_scalar)vi / (qaws_scalar)(v_samples - 1);
			qaws_surface_eval_result r;
			unsigned int idx = ui * v_samples + vi;
			double cr, cg, cb, t;

			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surf, u, v,
				QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);

			t = (double)kappa[idx] / (double)max_abs_k;
			curvature_color(t, &cr, &cg, &cb);
			obj_vertex_color(w, r.position, cr, cg, cb);

			if (flip)
			{
				qaws_vec3 fn;
				fn.x = -r.normal.x;
				fn.y = -r.normal.y;
				fn.z = -r.normal.z;
				obj_normal(w, fn);
			}
			else
			{
				obj_normal(w, r.normal);
			}
		}
	}

	/* Emit faces */
	for (ui = 0; ui < u_samples - 1; ui++)
	{
		for (vi = 0; vi < v_samples - 1; vi++)
		{
			unsigned int v00 = first_v + ui * v_samples + vi;
			unsigned int v01 = v00 + 1;
			unsigned int v10 = v00 + v_samples;
			unsigned int v11 = v10 + 1;
			unsigned int n00 = first_n + ui * v_samples + vi;
			unsigned int n01 = n00 + 1;
			unsigned int n10 = n00 + v_samples;
			unsigned int n11 = n10 + 1;
			if (flip)
			{
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v00, n00, v01, n01, v10, n10);
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v01, n01, v11, n11, v10, n10);
			}
			else
			{
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v00, n00, v10, n10, v01, n01);
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n", v10, n10, v11, n11, v01, n01);
			}
		}
	}

	free(kappa);
}

/* Draw control point grid as markers and wireframe */
static void obj_surface_cp_grid(
	obj_writer* w, qaws_vec3 const* cp,
	unsigned int u_count, unsigned int v_count,
	qaws_scalar marker_size)
{
	unsigned int i, j;
	unsigned int first_v;

	/* Draw markers at each CP */
	for (i = 0; i < u_count; i++)
		for (j = 0; j < v_count; j++)
			obj_sphere(w, cp[i * v_count + j], marker_size);

	/* Draw grid lines */
	first_v = w->vertex_count + 1;
	for (i = 0; i < u_count; i++)
		for (j = 0; j < v_count; j++)
			obj_vertex(w, cp[i * v_count + j]);

	/* U-direction lines (pairs for viewer compatibility) */
	for (i = 0; i < u_count; i++)
		for (j = 0; j + 1 < v_count; j++)
			fprintf(w->fp, "l %u %u\n",
				first_v + i * v_count + j,
				first_v + i * v_count + j + 1);
	/* V-direction lines */
	for (j = 0; j < v_count; j++)
		for (i = 0; i + 1 < u_count; i++)
			fprintf(w->fp, "l %u %u\n",
				first_v + i * v_count + j,
				first_v + (i + 1) * v_count + j);
}

#endif /* QAWS_TEST_COMMON_H */
