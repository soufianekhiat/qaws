#ifndef EXAMPLE_SVG_H
#define EXAMPLE_SVG_H

/* Small SVG writer, plotting and optimizer helpers shared by the examples. */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

static qaws_vec3 v3(qaws_scalar x, qaws_scalar y, qaws_scalar z)
{
	qaws_vec3 r;
	r.x = x;
	r.y = y;
	r.z = z;
	return r;
}

/* ================================================================== */
/*  Minimal SVG writer                                                */
/* ================================================================== */

typedef struct svg
{
	FILE* f;
} svg;

typedef struct viewport
{
	double x0, y0, w, h;         /* screen rectangle */
	double xmin, xmax, ymin, ymax; /* world rectangle */
} viewport;

static double vx(viewport const* v, double x) { return v->x0 + (x - v->xmin) / (v->xmax - v->xmin) * v->w; }
static double vy(viewport const* v, double y) { return v->y0 + v->h - (y - v->ymin) / (v->ymax - v->ymin) * v->h; }

static int svg_open(svg* s, char const* path, int w, int h, char const* title, char const* subtitle)
{
	s->f = fopen(path, "w");
	if (!s->f)
		return 0;
	fprintf(s->f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\" "
		"font-family=\"Segoe UI, Helvetica, Arial, sans-serif\">\n", w, h, w, h);
	fprintf(s->f, "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n", w, h);
	fprintf(s->f, "<text x=\"24\" y=\"36\" font-size=\"22\" font-weight=\"600\" fill=\"#1b1f24\">%s</text>\n", title);
	fprintf(s->f, "<text x=\"24\" y=\"60\" font-size=\"14\" fill=\"#57606a\">%s</text>\n", subtitle);
	return 1;
}

static void svg_close(svg* s)
{
	fprintf(s->f, "</svg>\n");
	fclose(s->f);
}

static void svg_panel(svg* s, viewport const* v, char const* label)
{
	fprintf(s->f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#f6f8fa\" stroke=\"#d0d7de\"/>\n",
		v->x0, v->y0, v->w, v->h);
	if (label)
		fprintf(s->f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">%s</text>\n",
			v->x0 + 10, v->y0 + 20, label);
}

static void svg_polyline(svg* s, double const* xy, int n, char const* color, double width, double opacity, int dashed)
{
	int i;
	fprintf(s->f, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-opacity=\"%.2f\" "
		"stroke-linejoin=\"round\" stroke-linecap=\"round\"%s points=\"", color, width, opacity,
		dashed ? " stroke-dasharray=\"5,4\"" : "");
	for (i = 0; i < n; i++)
		fprintf(s->f, "%.2f,%.2f ", xy[2 * i], xy[2 * i + 1]);
	fprintf(s->f, "\"/>\n");
}

static void svg_line(svg* s, double x0, double y0, double x1, double y1, char const* color, double width, double opacity)
{
	fprintf(s->f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"%.2f\" "
		"stroke-opacity=\"%.2f\" stroke-linecap=\"round\"/>\n", x0, y0, x1, y1, color, width, opacity);
}

static void svg_circle(svg* s, double x, double y, double r, char const* fill, char const* stroke)
{
	fprintf(s->f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"1\"/>\n",
		x, y, r, fill, stroke);
}

static void svg_text(svg* s, double x, double y, int size, char const* color, char const* anchor, char const* text)
{
	fprintf(s->f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"%d\" fill=\"%s\" text-anchor=\"%s\">%s</text>\n",
		x, y, size, color, anchor, text);
}

/* Perceptual heat ramp (dark blue -> teal -> yellow -> red). */
static void heat(double t, char* out)
{
	static double const stops[5][3] = {
		{ 0.18, 0.20, 0.55 }, { 0.13, 0.55, 0.65 }, { 0.36, 0.78, 0.40 }, { 0.98, 0.80, 0.18 }, { 0.86, 0.20, 0.15 } };
	double x;
	int i;
	if (t < 0) t = 0;
	if (t > 1) t = 1;
	x = t * 4.0;
	i = (int)x;
	if (i > 3) i = 3;
	x -= i;
	sprintf(out, "rgb(%d,%d,%d)",
		(int)(255 * (stops[i][0] + (stops[i + 1][0] - stops[i][0]) * x)),
		(int)(255 * (stops[i][1] + (stops[i + 1][1] - stops[i][1]) * x)),
		(int)(255 * (stops[i][2] + (stops[i + 1][2] - stops[i][2]) * x)));
}

static void svg_colorbar(svg* s, double x, double y, double w, double h, char const* lo, char const* hi)
{
	int i;
	char c[32];
	for (i = 0; i < 40; i++)
	{
		heat(i / 39.0, c);
		fprintf(s->f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.2f\" height=\"%.1f\" fill=\"%s\"/>\n",
			x + w * i / 40.0, y, w / 40.0 + 0.5, h, c);
	}
	svg_text(s, x, y + h + 14, 11, "#57606a", "start", lo);
	svg_text(s, x + w, y + h + 14, 11, "#57606a", "end", hi);
}

/* Log-scale loss plot inside a viewport. */
static void svg_loss_plot(svg* s, viewport const* v, double const* loss, int n, char const* color, char const* label)
{
	double lo = 1e300, hi = -1e300, xy[2 * 1024];
	int i, m = n > 1024 ? 1024 : n;
	char buf[96];
	for (i = 0; i < m; i++)
	{
		double l = log10(loss[i] > 1e-300 ? loss[i] : 1e-300);
		if (l < lo) lo = l;
		if (l > hi) hi = l;
	}
	if (hi - lo < 1e-9) hi = lo + 1;
	svg_panel(s, v, label);
	for (i = 0; i < m; i++)
	{
		double l = log10(loss[i] > 1e-300 ? loss[i] : 1e-300);
		xy[2 * i] = v->x0 + 12 + (v->w - 24) * i / (double)(m - 1);
		xy[2 * i + 1] = v->y0 + 32 + (v->h - 50) * (hi - l) / (hi - lo);
	}
	svg_polyline(s, xy, m, color, 2.2, 1.0, 0);
	sprintf(buf, "start %.3g", loss[0]);
	svg_text(s, v->x0 + 12, v->y0 + v->h - 6, 11, "#57606a", "start", buf);
	sprintf(buf, "end %.3g (%d iterations)", loss[m - 1], m - 1);
	svg_text(s, v->x0 + v->w - 12, v->y0 + v->h - 6, 11, "#57606a", "end", buf);
}

/* Adam on a flat parameter vector. */
typedef struct adam
{
	double m[256], v[256];
	int t;
} adam;

static void adam_step(adam* a, qaws_scalar* x, qaws_scalar const* g, int n, double lr)
{
	int i;
	a->t++;
	for (i = 0; i < n; i++)
	{
		double mh, vh;
		a->m[i] = 0.9 * a->m[i] + 0.1 * g[i];
		a->v[i] = 0.999 * a->v[i] + 0.001 * g[i] * g[i];
		mh = a->m[i] / (1 - pow(0.9, a->t));
		vh = a->v[i] / (1 - pow(0.999, a->t));
		x[i] -= (qaws_scalar)(lr * mh / (sqrt(vh) + 1e-12));
	}
}

static qaws_diff_views one_field(qaws_field_view* fv, qaws_diff_field field, qaws_scalar* data, unsigned int count, unsigned int comps)
{
	qaws_diff_views views;
	*fv = qaws_field_view_make(field, data, count, comps);
	views.fields = fv;
	views.field_count = 1;
	views.children = NULL;
	views.child_count = 0;
	return views;
}

static qaws_curve* bspline_2d(qaws_scalar const* cps, unsigned int n)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = n;
	d.is_uniform = 1;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static void curve_polyline(qaws_curve const* c, viewport const* v, double* xy, int n)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	int i;
	for (i = 0; i < n; i++)
	{
		qaws_eval_result_2d e;
		qaws_scalar t = r.min_value + (r.max_value - r.min_value) * i / (qaws_scalar)(n - 1);
		qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION, &e);
		xy[2 * i] = vx(v, e.position.x);
		xy[2 * i + 1] = vy(v, e.position.y);
	}
}

static void control_polygon(svg* s, viewport const* v, qaws_scalar const* cps, unsigned int n, char const* color)
{
	double xy[2 * 64];
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		xy[2 * i] = vx(v, cps[2 * i]);
		xy[2 * i + 1] = vy(v, cps[2 * i + 1]);
	}
	svg_polyline(s, xy, (int)n, color, 1.0, 0.6, 1);
	for (i = 0; i < n; i++)
		svg_circle(s, xy[2 * i], xy[2 * i + 1], 3.5, "#ffffff", color);
}

/* ================================================================== */
/*  Shared projection for 3D figures                                  */
/* ================================================================== */

typedef struct projection
{
	double cx, cy, scale, zscale;
} projection;

static void project(projection const* p, double x, double y, double z, double* sx, double* sy)
{
	*sx = p->cx + (x - y) * 0.866 * p->scale;
	*sy = p->cy + ((x + y) * 0.5 - z * p->zscale) * p->scale;
}

typedef double (*quad_value_fn)(void const* user, qaws_scalar u, qaws_scalar v);

/* Draws a (u,v) grid of quads, back to front, colored by value(u,v). */
static void draw_quads(svg* s, projection const* pr, qaws_surface const* surf, int nu, int nv,
	quad_value_fn value, void const* user, double lo, double hi)
{
	int i, j, k, n = nu * nv, *order;
	double* depth;
	char col[32];

	order = (int*)malloc(sizeof(int) * (size_t)n);
	depth = (double*)malloc(sizeof(double) * (size_t)n);
	for (k = 0; k < n; k++)
	{
		qaws_surface_eval_result r;
		i = k / nv;
		j = k % nv;
		qaws_surface_evaluate(surf, (qaws_scalar)((i + 0.5) / nu), (qaws_scalar)((j + 0.5) / nv), QAWS_SURFACE_EVAL_POSITION, &r);
		depth[k] = r.position.x + r.position.y;
		order[k] = k;
	}
	/* insertion sort by depth: far (small x+y) first */
	for (k = 1; k < n; k++)
	{
		int key = order[k], m = k - 1;
		while (m >= 0 && depth[order[m]] > depth[key])
		{
			order[m + 1] = order[m];
			m--;
		}
		order[m + 1] = key;
	}
	for (k = 0; k < n; k++)
	{
		qaws_surface_eval_result r[4];
		double x[4], y[4];
		int c;
		i = order[k] / nv;
		j = order[k] % nv;
		qaws_surface_evaluate(surf, (qaws_scalar)(i / (double)nu), (qaws_scalar)(j / (double)nv), QAWS_SURFACE_EVAL_POSITION, &r[0]);
		qaws_surface_evaluate(surf, (qaws_scalar)((i + 1) / (double)nu), (qaws_scalar)(j / (double)nv), QAWS_SURFACE_EVAL_POSITION, &r[1]);
		qaws_surface_evaluate(surf, (qaws_scalar)((i + 1) / (double)nu), (qaws_scalar)((j + 1) / (double)nv), QAWS_SURFACE_EVAL_POSITION, &r[2]);
		qaws_surface_evaluate(surf, (qaws_scalar)(i / (double)nu), (qaws_scalar)((j + 1) / (double)nv), QAWS_SURFACE_EVAL_POSITION, &r[3]);
		for (c = 0; c < 4; c++)
			project(pr, r[c].position.x, r[c].position.y, r[c].position.z, &x[c], &y[c]);
		heat((value(user, (qaws_scalar)((i + 0.5) / nu), (qaws_scalar)((j + 0.5) / nv)) - lo) / (hi - lo), col);
		fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"#ffffff\" stroke-width=\"0.35\" stroke-opacity=\"0.55\"/>\n",
			x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], col);
	}
	free(order);
	free(depth);
}

#endif /* EXAMPLE_SVG_H */
