/*
 * diff_showcase.c - Applications of qaws differentiability
 *
 * Writes SVG figures into ./showcase/:
 *   1_curve_fit.svg        B-spline fitted to data by adjoint gradient descent
 *   2_sensitivity.svg      "Affected by" (adjoint) and "Affects" (tangent) maps
 *   3_fairing.svg          curvature fairing through the geometry kernel
 *   4_nurbs_weight.svg     learning a NURBS weight: a parabola becomes a circle
 *   5_surface.svg          height-field fit of a B-spline surface (z only, by
 *                          component mask) and Gaussian curvature sensitivity
 *   6_vase.svg             surface of revolution fitted to scan points through
 *                          its profile curve (child views)
 *   7_coons.svg            Coons patch faired by editing two boundary curves
 *
 * Only the public qaws API is used. No finite differences anywhere: every
 * gradient comes from the library's adjoint rules.
 */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

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
/*  1. Fitting a B-spline to data with adjoints                      */
/* ================================================================== */

#define FIT_CP 9
#define FIT_SAMPLES 80
#define FIT_ITERS 600

static double fit_target_y(double x)
{
	return 0.55 * sin(2 * PI * x) + 0.25 * cos(5 * PI * x) * x;
}

static void demo_curve_fit(void)
{
	qaws_scalar cps[FIT_CP * 2], grad[FIT_CP * 2];
	qaws_scalar ts[FIT_SAMPLES];
	double tx[FIT_SAMPLES], ty[FIT_SAMPLES], loss[FIT_ITERS + 1];
	qaws_scalar snapshots[4][FIT_CP * 2];
	int snap_iter[4] = { 0, 25, 100, FIT_ITERS };
	int snap = 0, it, i;
	adam opt;
	svg s;
	viewport v = { 30, 80, 760, 470, -0.05, 1.05, -1.0, 1.0 };
	viewport lv = { 810, 80, 360, 220, 0, 0, 0, 0 };
	char buf[160];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < FIT_CP; i++)
	{
		cps[2 * i] = (qaws_scalar)i / (FIT_CP - 1);
		cps[2 * i + 1] = 0;
	}
	for (i = 0; i < FIT_SAMPLES; i++)
	{
		double x = i / (double)(FIT_SAMPLES - 1);
		tx[i] = x;
		ty[i] = fit_target_y(x) + 0.03 * sin(37.0 * x * x);
		ts[i] = (qaws_scalar)((FIT_CP - 3) * x);
	}

	for (it = 0; it <= FIT_ITERS; it++)
	{
		qaws_curve* c = bspline_2d(cps, FIT_CP);
		qaws_curve_jet_2d primal[FIT_SAMPLES], tangent[FIT_SAMPLES], ybar[FIT_SAMPLES];
		qaws_field_view fv;
		qaws_diff_views views;
		double l = 0;

		/* Forward: primal positions (no parameter tangent). */
		qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, FIT_SAMPLES, QAWS_EVAL_FLAG_POSITION, NULL, primal, tangent);
		memset(ybar, 0, sizeof(ybar));
		for (i = 0; i < FIT_SAMPLES; i++)
		{
			double dx = primal[i].d[0].x - tx[i], dy = primal[i].d[0].y - ty[i];
			l += (dx * dx + dy * dy) / FIT_SAMPLES;
			ybar[i].d[0].x = (qaws_scalar)(2 * dx / FIT_SAMPLES);
			ybar[i].d[0].y = (qaws_scalar)(2 * dy / FIT_SAMPLES);
		}
		loss[it] = l;
		if (snap < 4 && it == snap_iter[snap])
			memcpy(snapshots[snap++], cps, sizeof(cps));

		/* Reverse: one batch adjoint gives dL/dP for every control point. */
		memset(grad, 0, sizeof(grad));
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, FIT_CP, 2);
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, FIT_SAMPLES, QAWS_EVAL_FLAG_POSITION, ybar, &views, NULL);
		qaws_curve_destroy(c);
		if (it < FIT_ITERS)
			adam_step(&opt, cps, grad, FIT_CP * 2, 0.02);
	}

	svg_open(&s, "showcase/1_curve_fit.svg", 1200, 600, "Curve fitting with adjoints",
		"Cubic B-spline (9 control points) fitted to 80 samples: loss gradient = one batch adjoint per iteration, Adam updates");
	svg_panel(&s, &v, "data, intermediate fits and final curve");
	for (i = 0; i < FIT_SAMPLES; i++)
		svg_circle(&s, vx(&v, tx[i]), vy(&v, ty[i]), 2.6, "#d0d7de", "#8c959f");
	{
		static char const* colors[4] = { "#afb8c1", "#d29922", "#bf8700", "#0969da" };
		double xy[2 * 300];
		for (i = 0; i < 4; i++)
		{
			qaws_curve* c = bspline_2d(snapshots[i], FIT_CP);
			curve_polyline(c, &v, xy, 300);
			svg_polyline(&s, xy, 300, colors[i], i == 3 ? 3.0 : 1.6, i == 3 ? 1.0 : 0.8, i == 0);
			qaws_curve_destroy(c);
		}
		control_polygon(&s, &v, snapshots[3], FIT_CP, "#0969da");
	}
	svg_loss_plot(&s, &lv, loss, FIT_ITERS + 1, "#cf222e", "mean squared error (log scale)");
	sprintf(buf, "iterations shown: 0 (dashed), 25, 100, %d (blue)", FIT_ITERS);
	svg_text(&s, 810, 330, 13, "#24292f", "start", buf);
	svg_text(&s, 810, 352, 13, "#24292f", "start", "API: qaws_curve_eval_batch_tangent_2d (primal)");
	svg_text(&s, 810, 372, 13, "#24292f", "start", "     qaws_curve_eval_batch_adjoint_2d (dL/dP)");
	sprintf(buf, "final loss %.3e", loss[FIT_ITERS]);
	svg_text(&s, 810, 400, 15, "#0969da", "start", buf);
	svg_close(&s);
	printf("1_curve_fit: loss %.4e -> %.4e\n", loss[0], loss[FIT_ITERS]);
}

/* ================================================================== */
/*  2. Sensitivity: "affected by" and "affects"                       */
/* ================================================================== */

static void demo_sensitivity(void)
{
	static qaws_scalar const cps[10 * 2] = {
		0.0, 0.2, 0.8, 1.2, 1.6, 0.1, 2.4, 1.4, 3.2, 0.6, 4.0, 1.6, 4.8, 0.3, 5.6, 1.1, 6.4, 0.2, 7.2, 0.9 };
	qaws_curve* c = bspline_2d(cps, 10);
	qaws_range r = qaws_curve_get_parameter_range(c);
	viewport a = { 30, 80, 560, 440, -0.4, 7.6, -0.6, 2.0 };
	viewport b = { 610, 80, 560, 440, -0.4, 7.6, -0.6, 2.0 };
	qaws_scalar t_star = (qaws_scalar)2.6;
	unsigned int k_star = 5;
	double xy[2 * 400], wmax = 0, infl[10];
	qaws_scalar pbar[20];
	qaws_field_view fv;
	qaws_diff_views views;
	qaws_curve_jet_2d ybar, primal, tg;
	char col[32], buf[128];
	unsigned int i, j;
	svg s;

	/* Affected by: adjoint of the position at t*, seeded along x then y. */
	for (i = 0; i < 10; i++)
		infl[i] = 0;
	for (j = 0; j < 2; j++)
	{
		memset(pbar, 0, sizeof(pbar));
		memset(&ybar, 0, sizeof(ybar));
		if (j == 0) ybar.d[0].x = 1; else ybar.d[0].y = 1;
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, pbar, 10, 2);
		qaws_curve_eval_adjoint_2d(NULL, c, t_star, QAWS_EVAL_FLAG_POSITION, &ybar, &views, NULL);
		for (i = 0; i < 10; i++)
			infl[i] += pbar[2 * i] * pbar[2 * i] + pbar[2 * i + 1] * pbar[2 * i + 1];
	}
	for (i = 0; i < 10; i++)
	{
		infl[i] = sqrt(infl[i] / 2);
		if (infl[i] > wmax) wmax = infl[i];
	}

	svg_open(&s, "showcase/2_sensitivity.svg", 1200, 600, "Sensitivity: affected by / affects",
		"Left: adjoint of one curve point -> influence of each control point.  Right: tangent of one control point -> displacement along the curve.");

	svg_panel(&s, &a, "Affected by: which control points move the marked point?");
	curve_polyline(c, &a, xy, 400);
	svg_polyline(&s, xy, 400, "#57606a", 2.4, 1, 0);
	{
		double cxy[20];
		for (i = 0; i < 10; i++) { cxy[2 * i] = vx(&a, cps[2 * i]); cxy[2 * i + 1] = vy(&a, cps[2 * i + 1]); }
		svg_polyline(&s, cxy, 10, "#8c959f", 1, 0.6, 1);
		for (i = 0; i < 10; i++)
		{
			heat(infl[i] / wmax, col);
			svg_circle(&s, cxy[2 * i], cxy[2 * i + 1], 4 + 16 * infl[i] / wmax, infl[i] > 0 ? col : "#ffffff", "#24292f");
			sprintf(buf, "%.2f", infl[i]);
			svg_text(&s, cxy[2 * i], cxy[2 * i + 1] - 22, 11, "#24292f", "middle", buf);
		}
	}
	{
		qaws_curve_jet_2d p, t;
		qaws_curve_eval_tangent_2d(NULL, c, t_star, 0, QAWS_EVAL_FLAG_POSITION, NULL, &p, &t);
		svg_circle(&s, vx(&a, p.d[0].x), vy(&a, p.d[0].y), 7, "#cf222e", "#ffffff");
		sprintf(buf, "C(t = %.1f)", (double)t_star);
		svg_text(&s, vx(&a, p.d[0].x) + 10, vy(&a, p.d[0].y) + 22, 13, "#cf222e", "start", buf);
	}
	svg_text(&s, a.x0 + 12, a.y0 + a.h - 30, 12, "#57606a", "start",
		"only the 4 control points in the local support respond (numbers = basis weight)");
	svg_colorbar(&s, a.x0 + 12, a.y0 + a.h - 22, 200, 8, "0", "max");

	/* Affects: tangent of every curve point when control point k moves up. */
	svg_panel(&s, &b, "Affects: where does moving one control point act?");
	{
		qaws_scalar dir[20];
		double mags[400], mmax = 0;
		memset(dir, 0, sizeof(dir));
		dir[2 * k_star + 1] = 1;
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, dir, 10, 2);
		for (i = 0; i < 400; i++)
		{
			qaws_scalar t = r.min_value + (r.max_value - r.min_value) * i / (qaws_scalar)399;
			qaws_curve_eval_tangent_2d(NULL, c, t, 0, QAWS_EVAL_FLAG_POSITION, &views, &primal, &tg);
			xy[2 * i] = vx(&b, primal.d[0].x);
			xy[2 * i + 1] = vy(&b, primal.d[0].y);
			mags[i] = sqrt(tg.d[0].x * tg.d[0].x + tg.d[0].y * tg.d[0].y);
			if (mags[i] > mmax) mmax = mags[i];
		}
		for (i = 0; i + 1 < 400; i++)
		{
			heat(mags[i] / mmax, col);
			svg_line(&s, xy[2 * i], xy[2 * i + 1], xy[2 * i + 2], xy[2 * i + 3], col, 5, 1);
		}
		for (i = 0; i < 400; i += 12)
			if (mags[i] > 0.02 * mmax)
				svg_line(&s, xy[2 * i], xy[2 * i + 1], xy[2 * i], xy[2 * i + 1] - 60 * mags[i], "#cf222e", 1.5, 0.8);
		{
			double cxy[20];
			for (i = 0; i < 10; i++) { cxy[2 * i] = vx(&b, cps[2 * i]); cxy[2 * i + 1] = vy(&b, cps[2 * i + 1]); }
			svg_polyline(&s, cxy, 10, "#8c959f", 1, 0.6, 1);
			for (i = 0; i < 10; i++)
				svg_circle(&s, cxy[2 * i], cxy[2 * i + 1], i == k_star ? 8 : 3.5, i == k_star ? "#cf222e" : "#ffffff", "#24292f");
			svg_line(&s, cxy[2 * k_star], cxy[2 * k_star + 1], cxy[2 * k_star], cxy[2 * k_star + 1] - 50, "#cf222e", 2.5, 1);
		}
		svg_text(&s, b.x0 + 12, b.y0 + b.h - 30, 12, "#57606a", "start",
			"color and arrows = |dC/dP_k| for an upward move of the red control point");
		svg_colorbar(&s, b.x0 + 12, b.y0 + b.h - 22, 200, 8, "0", "max");
	}
	svg_close(&s);
	qaws_curve_destroy(c);
	printf("2_sensitivity: influence weights at t*=%.2f written\n", (double)t_star);
}

/* ================================================================== */
/*  3. Curvature fairing through the geometry kernel                  */
/* ================================================================== */

#define FAIR_CP 16
#define FAIR_SAMPLES 160
#define FAIR_ITERS 500

static void fairing_energy(qaws_scalar const* cps, qaws_scalar const* ref, qaws_scalar const* ts,
	double lambda, double* out_energy, qaws_scalar* grad, double* kappa_out)
{
	qaws_curve* c = bspline_2d(cps, FAIR_CP);
	qaws_curve_jet_2d primal[FAIR_SAMPLES], tangent[FAIR_SAMPLES], ybar[FAIR_SAMPLES];
	qaws_field_view fv;
	qaws_diff_views views;
	double e = 0;
	int i;

	qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, FAIR_SAMPLES, 0x7, NULL, primal, tangent);
	memset(ybar, 0, sizeof(ybar));
	for (i = 0; i < FAIR_SAMPLES; i++)
	{
		qaws_curve_geometry_2d g, gbar;
		qaws_curve_geometry_eval_2d(&primal[i], NULL, NULL, &g, NULL, NULL, NULL);
		if (kappa_out)
			kappa_out[i] = g.curvature;
		/* E = mean(kappa^2 * speed) + lambda * mean |C - C_ref|^2 */
		e += g.curvature * g.curvature * g.speed / FAIR_SAMPLES;
		memset(&gbar, 0, sizeof(gbar));
		gbar.curvature = (qaws_scalar)(2 * g.curvature * g.speed / FAIR_SAMPLES);
		gbar.speed = (qaws_scalar)(g.curvature * g.curvature / FAIR_SAMPLES);
		qaws_curve_geometry_adjoint_2d(&primal[i], &gbar, &ybar[i], NULL);
	}
	{
		qaws_curve* r = bspline_2d(ref, FAIR_CP);
		for (i = 0; i < FAIR_SAMPLES; i++)
		{
			qaws_eval_result_2d er;
			double dx, dy;
			qaws_curve_evaluate_2d(r, ts[i], QAWS_EVAL_FLAG_POSITION, &er);
			dx = primal[i].d[0].x - er.position.x;
			dy = primal[i].d[0].y - er.position.y;
			e += lambda * (dx * dx + dy * dy) / FAIR_SAMPLES;
			ybar[i].d[0].x += (qaws_scalar)(2 * lambda * dx / FAIR_SAMPLES);
			ybar[i].d[0].y += (qaws_scalar)(2 * lambda * dy / FAIR_SAMPLES);
			ybar[i].channels |= QAWS_EVAL_FLAG_POSITION;
		}
		qaws_curve_destroy(r);
	}
	memset(grad, 0, sizeof(qaws_scalar) * FAIR_CP * 2);
	views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, FAIR_CP, 2);
	qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, FAIR_SAMPLES, 0x7, ybar, &views, NULL);
	qaws_curve_destroy(c);
	*out_energy = e;
}

static void draw_comb(svg* s, viewport const* v, qaws_scalar const* cps, qaws_scalar const* ts, char const* curve_color)
{
	qaws_curve* c = bspline_2d(cps, FAIR_CP);
	double xy[2 * FAIR_SAMPLES];
	int i;
	for (i = 0; i < FAIR_SAMPLES; i++)
	{
		qaws_curve_jet_2d p, t;
		qaws_curve_geometry_2d g;
		double tipx, tipy;
		qaws_curve_eval_tangent_2d(NULL, c, ts[i], 0, 0x7, NULL, &p, &t);
		qaws_curve_geometry_eval_2d(&p, NULL, NULL, &g, NULL, NULL, NULL);
		xy[2 * i] = vx(v, p.d[0].x);
		xy[2 * i + 1] = vy(v, p.d[0].y);
		tipx = vx(v, p.d[0].x - 0.12 * g.curvature * g.normal.x);
		tipy = vy(v, p.d[0].y - 0.12 * g.curvature * g.normal.y);
		svg_line(s, xy[2 * i], xy[2 * i + 1], tipx, tipy, "#bf8700", 1, 0.55);
	}
	svg_polyline(s, xy, FAIR_SAMPLES, curve_color, 2.6, 1, 0);
	qaws_curve_destroy(c);
}

static void demo_fairing(void)
{
	qaws_scalar cps[FAIR_CP * 2], ref[FAIR_CP * 2], grad[FAIR_CP * 2], ts[FAIR_SAMPLES];
	double energy[FAIR_ITERS + 1], k0[FAIR_SAMPLES], k1[FAIR_SAMPLES];
	viewport a = { 30, 80, 560, 300, -0.3, 7.8, -1.6, 1.6 };
	viewport b = { 610, 80, 560, 300, -0.3, 7.8, -1.6, 1.6 };
	viewport kv = { 30, 400, 760, 180, 0, 0, 0, 0 };
	viewport lv = { 810, 400, 360, 180, 0, 0, 0, 0 };
	adam opt;
	int i, it;
	svg s;
	char buf[128];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < FAIR_CP; i++)
	{
		double x = 7.5 * i / (FAIR_CP - 1);
		cps[2 * i] = (qaws_scalar)(x + 0.12 * sin(9.1 * i));
		cps[2 * i + 1] = (qaws_scalar)(0.8 * sin(0.9 * x) + 0.35 * sin(5.3 * i + 0.4));
	}
	memcpy(ref, cps, sizeof(cps));
	for (i = 0; i < FAIR_SAMPLES; i++)
		ts[i] = (qaws_scalar)((FAIR_CP - 3) * (i + 0.5) / FAIR_SAMPLES);

	for (it = 0; it <= FAIR_ITERS; it++)
	{
		fairing_energy(cps, ref, ts, 0.6, &energy[it], grad, it == 0 ? k0 : (it == FAIR_ITERS ? k1 : NULL));
		if (it < FAIR_ITERS)
			adam_step(&opt, cps, grad, FAIR_CP * 2, 0.01);
	}

	svg_open(&s, "showcase/3_fairing.svg", 1200, 600, "Curvature fairing through the geometry kernel",
		"Minimize mean(kappa^2 |C'|) + 0.6 mean|C - C0|^2: curvature adjoint -> jet adjoint (D1, D2) -> control points. Orange = curvature comb.");
	svg_panel(&s, &a, "before");
	draw_comb(&s, &a, ref, ts, "#57606a");
	control_polygon(&s, &a, ref, FAIR_CP, "#8c959f");
	svg_panel(&s, &b, "after fairing");
	draw_comb(&s, &b, cps, ts, "#0969da");
	control_polygon(&s, &b, cps, FAIR_CP, "#0969da");

	svg_panel(&s, &kv, "curvature along the curve (gray before, blue after)");
	{
		double xy0[2 * FAIR_SAMPLES], xy1[2 * FAIR_SAMPLES], km = 0;
		for (i = 0; i < FAIR_SAMPLES; i++)
			if (fabs(k0[i]) > km) km = fabs(k0[i]);
		for (i = 0; i < FAIR_SAMPLES; i++)
		{
			double x = kv.x0 + 12 + (kv.w - 24) * i / (double)(FAIR_SAMPLES - 1);
			xy0[2 * i] = xy1[2 * i] = x;
			xy0[2 * i + 1] = kv.y0 + kv.h / 2 + 10 - (kv.h / 2 - 20) * k0[i] / km;
			xy1[2 * i + 1] = kv.y0 + kv.h / 2 + 10 - (kv.h / 2 - 20) * k1[i] / km;
		}
		svg_line(&s, kv.x0 + 12, kv.y0 + kv.h / 2 + 10, kv.x0 + kv.w - 12, kv.y0 + kv.h / 2 + 10, "#d0d7de", 1, 1);
		svg_polyline(&s, xy0, FAIR_SAMPLES, "#8c959f", 1.8, 1, 0);
		svg_polyline(&s, xy1, FAIR_SAMPLES, "#0969da", 2.2, 1, 0);
	}
	svg_loss_plot(&s, &lv, energy, FAIR_ITERS + 1, "#cf222e", "fairing energy (log scale)");
	svg_close(&s);
	sprintf(buf, "%.4e -> %.4e", energy[0], energy[FAIR_ITERS]);
	printf("3_fairing: energy %s\n", buf);
}

/* ================================================================== */
/*  4. Learning a NURBS weight: parabola -> exact circle              */
/* ================================================================== */

#define W_ITERS 80

static qaws_curve* quarter(qaws_scalar const* w)
{
	static qaws_scalar const cps[6] = { 1, 0, 1, 1, 0, 1 };
	qaws_rational_bezier_desc d;
	qaws_curve* c = NULL;
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = cps;
	d.control_point_count = 3;
	d.weights = w;
	d.weight_count = 3;
	qaws_curve_create_rational_bezier(&d, &c);
	return c;
}

static void demo_nurbs_weight(void)
{
	qaws_scalar w[3] = { 1, 1, 1 };
	qaws_scalar ts[64];
	double loss[W_ITERS + 1], wh[W_ITERS + 1];
	unsigned char active[3] = { 0, 1, 0 };
	viewport a = { 30, 80, 520, 500, -0.15, 1.25, -0.15, 1.25 };
	viewport lv = { 580, 80, 590, 220, 0, 0, 0, 0 };
	viewport wv = { 580, 320, 590, 220, 0, W_ITERS, 0.6, 1.05 };
	int it, i;
	svg s;
	char buf[512];

	for (i = 0; i < 64; i++)
		ts[i] = (qaws_scalar)(i / 63.0);

	for (it = 0; it <= W_ITERS; it++)
	{
		qaws_curve* c = quarter(w);
		qaws_curve_jet_2d p[64], t[64], ybar[64];
		qaws_scalar wbar[3] = { 0, 0, 0 };
		qaws_field_view fv;
		qaws_diff_views views;
		double l = 0;

		qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, 64, QAWS_EVAL_FLAG_POSITION, NULL, p, t);
		memset(ybar, 0, sizeof(ybar));
		for (i = 0; i < 64; i++)
		{
			double r2 = p[i].d[0].x * p[i].d[0].x + p[i].d[0].y * p[i].d[0].y;
			double e = r2 - 1;
			l += e * e / 64;
			ybar[i].d[0].x = (qaws_scalar)(4 * e * p[i].d[0].x / 64);
			ybar[i].d[0].y = (qaws_scalar)(4 * e * p[i].d[0].y / 64);
		}
		loss[it] = l > 1e-30 ? l : 1e-30;
		wh[it] = w[1];

		/* Only the middle weight is active: the mask keeps the end weights fixed. */
		views = one_field(&fv, QAWS_FIELD_WEIGHTS, wbar, 3, 1);
		fv.active = active;
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, 64, QAWS_EVAL_FLAG_POSITION, ybar, &views, NULL);
		qaws_curve_destroy(c);
		if (it < W_ITERS)
			w[1] -= (qaws_scalar)(1.5 * wbar[1]);
	}

	sprintf(buf, "Rational quadratic, end weights fixed by an activity mask; gradient of mean(|C|^2 - 1)^2 w.r.t. the middle weight. Learned w = %.6f (exact: sqrt(2)/2 = 0.707107)", (double)w[1]);
	svg_open(&s, "showcase/4_nurbs_weight.svg", 1200, 600, "Learning a NURBS weight: a parabola becomes a circle", buf);
	svg_panel(&s, &a, "quarter arc");
	{
		double xy[2 * 200];
		qaws_scalar w0[3] = { 1, 1, 1 };
		qaws_curve* c0 = quarter(w0);
		qaws_curve* c1 = quarter(w);
		for (i = 0; i < 200; i++)
		{
			double th = 0.5 * PI * i / 199.0;
			xy[2 * i] = vx(&a, cos(th));
			xy[2 * i + 1] = vy(&a, sin(th));
		}
		svg_polyline(&s, xy, 200, "#d0d7de", 9, 1, 0);
		curve_polyline(c0, &a, xy, 200);
		svg_polyline(&s, xy, 200, "#8c959f", 2, 1, 1);
		curve_polyline(c1, &a, xy, 200);
		svg_polyline(&s, xy, 200, "#0969da", 2.6, 1, 0);
		svg_circle(&s, vx(&a, 1), vy(&a, 0), 4, "#ffffff", "#24292f");
		svg_circle(&s, vx(&a, 1), vy(&a, 1), 6, "#cf222e", "#24292f");
		svg_circle(&s, vx(&a, 0), vy(&a, 1), 4, "#ffffff", "#24292f");
		svg_text(&s, vx(&a, 1) - 8, vy(&a, 1) - 10, 12, "#cf222e", "end", "P1 (learned weight)");
		svg_text(&s, a.x0 + 12, a.y0 + a.h - 30, 12, "#57606a", "start", "wide gray: unit circle   dashed: w = 1 (parabola)");
		svg_text(&s, a.x0 + 12, a.y0 + a.h - 12, 12, "#0969da", "start", "blue: learned rational curve");
		qaws_curve_destroy(c0);
		qaws_curve_destroy(c1);
	}
	svg_loss_plot(&s, &lv, loss, W_ITERS + 1, "#cf222e", "loss (log scale)");
	svg_panel(&s, &wv, "middle weight over iterations");
	{
		double xy[2 * (W_ITERS + 1)];
		for (i = 0; i <= W_ITERS; i++)
		{
			xy[2 * i] = wv.x0 + 12 + (wv.w - 24) * i / (double)W_ITERS;
			xy[2 * i + 1] = wv.y0 + 30 + (wv.h - 50) * (1.05 - wh[i]) / 0.45;
		}
		svg_line(&s, wv.x0 + 12, wv.y0 + 30 + (wv.h - 50) * (1.05 - 0.7071068) / 0.45,
			wv.x0 + wv.w - 12, wv.y0 + 30 + (wv.h - 50) * (1.05 - 0.7071068) / 0.45, "#2da44e", 1.5, 1);
		svg_text(&s, wv.x0 + wv.w - 14, wv.y0 + 26 + (wv.h - 50) * (1.05 - 0.7071068) / 0.45, 11, "#2da44e", "end", "sqrt(2)/2");
		svg_polyline(&s, xy, W_ITERS + 1, "#0969da", 2.4, 1, 0);
	}
	svg_close(&s);
	printf("4_nurbs_weight: w = %.7f (target 0.7071068)\n", (double)w[1]);
}

/* ================================================================== */
/*  5. Surface: masked height-field fit + curvature sensitivity       */
/* ================================================================== */

#define SN 6
#define SG 18
#define S_ITERS 400

static double target_height(double x, double y)
{
	double r2 = (x - 0.45) * (x - 0.45) + (y - 0.55) * (y - 0.55);
	return 0.55 * exp(-r2 / 0.08) - 0.25 * exp(-((x - 0.8) * (x - 0.8) + (y - 0.2) * (y - 0.2)) / 0.05);
}

static qaws_surface* bspline_surface(qaws_scalar const* cps)
{
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = SN;
	d.v_point_count = SN;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

/* Isometric projection of a world point into a viewport. */
static void iso(viewport const* v, double x, double y, double z, double* sx, double* sy)
{
	double px = (x - y) * 0.866, py = (x + y) * 0.5 - z * 1.4;
	*sx = v->x0 + v->w / 2 + px * v->w * 0.5;
	*sy = v->y0 + v->h * 0.34 + py * v->h * 0.42;
}

static void draw_surface(svg* s, viewport const* v, qaws_surface const* surf, double zlo, double zhi)
{
	int i, j;
	char col[32];
	for (i = 0; i < 24; i++)
		for (j = 0; j < 24; j++)
		{
			qaws_scalar u0 = (qaws_scalar)(i / 24.0), u1 = (qaws_scalar)((i + 1) / 24.0);
			qaws_scalar v0 = (qaws_scalar)(j / 24.0), v1 = (qaws_scalar)((j + 1) / 24.0);
			qaws_surface_eval_result r[4];
			double x[4], y[4], zc = 0;
			int k;
			qaws_surface_evaluate(surf, u0, v0, QAWS_SURFACE_EVAL_POSITION, &r[0]);
			qaws_surface_evaluate(surf, u1, v0, QAWS_SURFACE_EVAL_POSITION, &r[1]);
			qaws_surface_evaluate(surf, u1, v1, QAWS_SURFACE_EVAL_POSITION, &r[2]);
			qaws_surface_evaluate(surf, u0, v1, QAWS_SURFACE_EVAL_POSITION, &r[3]);
			for (k = 0; k < 4; k++)
			{
				iso(v, r[k].position.x, r[k].position.y, r[k].position.z, &x[k], &y[k]);
				zc += r[k].position.z / 4;
			}
			heat((zc - zlo) / (zhi - zlo), col);
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"#ffffff\" stroke-width=\"0.4\" stroke-opacity=\"0.6\"/>\n",
				x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], col);
		}
}

static void demo_surface(void)
{
	qaws_scalar cps[SN * SN * 3], grad[SN * SN * 3];
	qaws_scalar us[SG * SG], vs[SG * SG];
	double tz[SG * SG], loss[S_ITERS + 1];
	viewport a = { 30, 80, 370, 420, 0, 0, 0, 0 };
	viewport b = { 415, 80, 370, 420, 0, 0, 0, 0 };
	viewport c = { 800, 80, 370, 420, 0, 0, 0, 0 };
	viewport lv = { 30, 515, 755, 75, 0, 0, 0, 0 };
	adam opt;
	int i, j, it;
	svg s;
	char buf[200];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < SN; i++)
		for (j = 0; j < SN; j++)
		{
			qaws_scalar* p = &cps[(i * SN + j) * 3];
			p[0] = (qaws_scalar)i / (SN - 1);
			p[1] = (qaws_scalar)j / (SN - 1);
			p[2] = 0;
		}
	for (i = 0; i < SG; i++)
		for (j = 0; j < SG; j++)
		{
			us[i * SG + j] = (qaws_scalar)(i / (SG - 1.0));
			vs[i * SG + j] = (qaws_scalar)(j / (SG - 1.0));
		}

	svg_open(&s, "showcase/5_surface.svg", 1200, 600, "Surfaces: masked height-field fit and curvature sensitivity",
		"Bicubic B-spline (6x6 net). Fit: only z is active (component mask), x/y stay on the grid. Right: d(Gaussian curvature at the red point)/d(control point z).");
	svg_panel(&s, &a, "initial (flat) + target samples");
	{
		qaws_surface* s0 = bspline_surface(cps);
		draw_surface(&s, &a, s0, -0.3, 0.6);
		qaws_surface_destroy(s0);
	}

	for (it = 0; it <= S_ITERS; it++)
	{
		qaws_surface* sf = bspline_surface(cps);
		static qaws_surface_jet primal[SG * SG], tangent[SG * SG], ybar[SG * SG];
		qaws_field_view fv;
		qaws_diff_views views;
		double l = 0;

		qaws_surface_eval_batch_tangent(NULL, sf, us, vs, NULL, NULL, SG * SG, QAWS_SJET_P, NULL, primal, tangent);
		memset(ybar, 0, sizeof(ybar));
		for (i = 0; i < SG * SG; i++)
		{
			double dz;
			if (it == 0)
				tz[i] = target_height(primal[i].d[0].x, primal[i].d[0].y);
			dz = primal[i].d[0].z - tz[i];
			l += dz * dz / (SG * SG);
			ybar[i].d[0].z = (qaws_scalar)(2 * dz / (SG * SG));
		}
		loss[it] = l;

		memset(grad, 0, sizeof(grad));
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, SN * SN, 3);
		fv.component_mask = 1u << 2; /* z only */
		qaws_surface_eval_batch_adjoint(NULL, sf, us, vs, SG * SG, QAWS_SJET_P, ybar, &views, NULL, NULL);
		qaws_surface_destroy(sf);
		if (it < S_ITERS)
			adam_step(&opt, cps, grad, SN * SN * 3, 0.01);
	}

	for (i = 0; i < SG * SG; i += 3)
	{
		double x, y;
		iso(&a, us[i], vs[i], tz[i], &x, &y);
		svg_circle(&s, x, y, 1.6, "#24292f", "none");
	}

	svg_panel(&s, &b, "after fitting (z only)");
	{
		qaws_surface* s1 = bspline_surface(cps);
		draw_surface(&s, &b, s1, -0.3, 0.6);
		qaws_surface_destroy(s1);
	}

	/* Gaussian curvature sensitivity: K at one point -> geometry adjoint ->
	   surface adjoint on the control net. */
	svg_panel(&s, &c, "dK / d(z_ij) at the red point");
	{
		qaws_surface* s1 = bspline_surface(cps);
		qaws_scalar u0 = (qaws_scalar)0.30, v0 = (qaws_scalar)0.70;
		qaws_surface_jet p, ybar;
		qaws_surface_geometry g, gbar;
		qaws_field_view fv;
		qaws_diff_views views;
		double gmax = 0, x, y;
		char col[32];

		qaws_surface_eval_jet(s1, u0, v0, QAWS_SJET_ORDER2, &p);
		qaws_surface_geometry_eval(&p, NULL, NULL, &g, NULL, NULL, NULL);
		memset(&gbar, 0, sizeof(gbar));
		gbar.gaussian = 1;
		memset(&ybar, 0, sizeof(ybar));
		qaws_surface_geometry_adjoint(&p, &gbar, &ybar, NULL);
		memset(grad, 0, sizeof(grad));
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, SN * SN, 3);
		qaws_surface_eval_adjoint(NULL, s1, u0, v0, ybar.channels, &ybar, &views, NULL, NULL);
		for (i = 0; i < SN * SN; i++)
			if (fabs(grad[i * 3 + 2]) > gmax) gmax = fabs(grad[i * 3 + 2]);

		/* Top view of the control net: one cell per control point, colored
		   by dK/dz (diverging: blue lowers K, red raises it). */
		{
			double gx0 = c.x0 + 45, gy0 = c.y0 + 45, cell = (c.w - 90) / SN;
			for (i = 0; i < SN; i++)
				for (j = 0; j < SN; j++)
				{
					double sv = grad[(i * SN + j) * 3 + 2] / gmax;
					x = gx0 + i * cell;
					y = gy0 + (SN - 1 - j) * cell;
					heat(0.5 + 0.5 * sv, col);
					fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\" stroke=\"#ffffff\" stroke-width=\"2\"/>\n",
						x, y, cell, cell, sv == 0 ? "#eaeef2" : col);
					sprintf(buf, "%+.2f", sv);
					svg_text(&s, x + cell / 2, y + cell / 2 + 4, 11, "#1b1f24", "middle", sv == 0 ? "0" : buf);
				}
			x = gx0 + p.d[0].x * (SN - 1) * cell + cell / 2;
			y = gy0 + (1 - p.d[0].y) * (SN - 1) * cell + cell / 2;
			svg_circle(&s, x, y, 7, "#cf222e", "#ffffff");
			svg_text(&s, gx0, gy0 - 8, 11, "#57606a", "start", "top view of the 6x6 control net (normalized dK/dz)");
		}
		sprintf(buf, "K = %.3f, mean = %.3f", (double)g.gaussian, (double)g.mean);
		svg_text(&s, c.x0 + 12, c.y0 + c.h - 34, 12, "#24292f", "start", buf);
		svg_colorbar(&s, c.x0 + 12, c.y0 + c.h - 24, 200, 8, "lowers K", "raises K");
		qaws_surface_destroy(s1);
	}
	svg_loss_plot(&s, &lv, loss, S_ITERS + 1, "#cf222e", "height-fit MSE (log)");
	svg_close(&s);
	printf("5_surface: fit MSE %.4e -> %.4e\n", loss[0], loss[S_ITERS]);
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

/* ================================================================== */
/*  6. Vase from scanned points: fitting a revolution profile         */
/* ================================================================== */

#define VASE_CP 7
#define VASE_RINGS 22
#define VASE_SPOKES 14
#define VASE_N (VASE_RINGS * VASE_SPOKES)
#define VASE_ITERS 500

static double vase_radius(double z)
{
	return 0.85 + 0.38 * sin(2.1 * z + 0.3) - 0.08 * z;
}

static qaws_surface* vase_surface(qaws_curve* profile)
{
	qaws_surface_revolution_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.profile = profile;
	d.axis_origin = v3(0, 0, 0);
	d.axis_direction = v3(0, 0, 1);
	d.angle = 0; /* full turn */
	qaws_surface_create_revolution(&d, &s);
	return s;
}

static double vase_height_value(void const* user, qaws_scalar u, qaws_scalar v)
{
	(void)user;
	(void)u;
	return v;
}

static void demo_vase(void)
{
	qaws_scalar cps[VASE_CP * 2], grad[VASE_CP * 2];
	qaws_scalar us[VASE_N], vs[VASE_N];
	qaws_vec3 targets[VASE_N];
	double loss[VASE_ITERS + 1];
	qaws_scalar initial[VASE_CP * 2];
	adam opt;
	int i, j, it;
	svg s;
	char buf[200];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < VASE_CP; i++)
	{
		cps[2 * i] = 1.0f;                                  /* radius */
		cps[2 * i + 1] = (qaws_scalar)(3.0 * i / (VASE_CP - 1)); /* height */
	}
	memcpy(initial, cps, sizeof(cps));
	for (i = 0; i < VASE_RINGS; i++)
		for (j = 0; j < VASE_SPOKES; j++)
		{
			int k = i * VASE_SPOKES + j;
			double v = (i + 0.5) / VASE_RINGS, u = (j + 0.25 * (i % 2)) / VASE_SPOKES;
			double z = 3.0 * v, r = vase_radius(z) * (1 + 0.015 * sin(13.0 * k));
			us[k] = (qaws_scalar)u;
			vs[k] = (qaws_scalar)v;
			targets[k] = v3((qaws_scalar)(r * cos(2 * PI * u)), (qaws_scalar)(r * sin(2 * PI * u)), (qaws_scalar)z);
		}

	for (it = 0; it <= VASE_ITERS; it++)
	{
		qaws_curve* profile = bspline_2d(cps, VASE_CP);
		qaws_surface* vase = vase_surface(profile);
		static qaws_surface_jet primal[VASE_N], tangent[VASE_N], ybar[VASE_N];
		qaws_field_view fv;
		qaws_diff_views child, views;
		double l = 0;

		qaws_surface_eval_batch_tangent(NULL, vase, us, vs, NULL, NULL, VASE_N, QAWS_SJET_P, NULL, primal, tangent);
		memset(ybar, 0, sizeof(ybar));
		for (i = 0; i < VASE_N; i++)
		{
			qaws_vec3 d = v3(primal[i].d[0].x - targets[i].x, primal[i].d[0].y - targets[i].y, primal[i].d[0].z - targets[i].z);
			l += (d.x * d.x + d.y * d.y + d.z * d.z) / VASE_N;
			ybar[i].d[0] = v3((qaws_scalar)(2 * d.x / VASE_N), (qaws_scalar)(2 * d.y / VASE_N), (qaws_scalar)(2 * d.z / VASE_N));
		}
		loss[it] = l;

		/* The revolution has no own fields in this fit; its profile is
		   child 0, so the profile control point adjoints live there. */
		memset(grad, 0, sizeof(grad));
		child = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, VASE_CP, 2);
		views.fields = NULL;
		views.field_count = 0;
		views.children = &child;
		views.child_count = 1;
		qaws_surface_eval_batch_adjoint(NULL, vase, us, vs, VASE_N, QAWS_SJET_P, ybar, &views, NULL, NULL);
		qaws_surface_destroy(vase);
		qaws_curve_destroy(profile);
		if (it < VASE_ITERS)
			adam_step(&opt, cps, grad, VASE_CP * 2, 0.01);
	}

	svg_open(&s, "showcase/6_vase.svg", 1200, 600, "Vase from scanned points: fitting a surface of revolution",
		"308 noisy points; the gradient flows surface -> revolution rule -> child[0] (profile B-spline, 7 control points). Color = height parameter.");
	{
		viewport a = { 30, 80, 370, 500, 0, 0, 0, 0 };
		viewport b = { 415, 80, 370, 500, 0, 0, 0, 0 };
		viewport c = { 800, 80, 370, 280, -0.2, 1.8, -0.2, 3.2 };
		viewport lv = { 800, 380, 370, 200, 0, 0, 0, 0 };
		projection pa = { 215, 470, 95, 1.0 }, pb = { 600, 470, 95, 1.0 };
		qaws_curve* p0 = bspline_2d(initial, VASE_CP);
		qaws_curve* p1 = bspline_2d(cps, VASE_CP);
		qaws_surface* s0 = vase_surface(p0);
		qaws_surface* s1 = vase_surface(p1);
		double xy[2 * 200];

		svg_panel(&s, &a, "initial profile (cylinder) + scan points");
		draw_quads(&s, &pa, s0, 28, 16, vase_height_value, NULL, 0, 1);
		for (i = 0; i < VASE_N; i++)
		{
			double x, y;
			project(&pa, targets[i].x, targets[i].y, targets[i].z, &x, &y);
			svg_circle(&s, x, y, 1.5, "#24292f", "none");
		}
		svg_panel(&s, &b, "fitted vase");
		draw_quads(&s, &pb, s1, 28, 16, vase_height_value, NULL, 0, 1);

		/* Profile plot: (radius, height). */
		svg_panel(&s, &c, "profile r(z): dashed initial, blue fitted, green truth");
		for (i = 0; i < 200; i++)
		{
			double z = 3.0 * i / 199.0;
			xy[2 * i] = vx(&c, vase_radius(z));
			xy[2 * i + 1] = vy(&c, z);
		}
		svg_polyline(&s, xy, 200, "#2da44e", 5, 0.45, 0);
		curve_polyline(p0, &c, xy, 200);
		svg_polyline(&s, xy, 200, "#8c959f", 1.6, 1, 1);
		curve_polyline(p1, &c, xy, 200);
		svg_polyline(&s, xy, 200, "#0969da", 2.4, 1, 0);
		control_polygon(&s, &c, cps, VASE_CP, "#0969da");
		svg_loss_plot(&s, &lv, loss, VASE_ITERS + 1, "#cf222e", "point-to-surface MSE (log)");

		qaws_surface_destroy(s0);
		qaws_surface_destroy(s1);
		qaws_curve_destroy(p0);
		qaws_curve_destroy(p1);
	}
	svg_close(&s);
	sprintf(buf, "%.4e -> %.4e", loss[0], loss[VASE_ITERS]);
	printf("6_vase: MSE %s\n", buf);
}

/* ================================================================== */
/*  7. Coons patch fairing by editing two boundaries                  */
/* ================================================================== */

#define CB_CP 6
#define CB_GRID 9
#define CB_ITERS 400

static qaws_scalar const g_cb_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_curve* cb_curve(qaws_scalar const* p)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = CB_CP;
	d.knots = g_cb_knots;
	d.knot_count = 10;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

typedef struct cb_patch
{
	qaws_curve* curves[4];
	qaws_surface* surface;
} cb_patch;

static void cb_build(qaws_scalar (*p)[CB_CP * 3], cb_patch* out)
{
	qaws_surface_coons_desc d;
	int i;
	for (i = 0; i < 4; i++)
		out->curves[i] = cb_curve(p[i]);
	memset(&d, 0, sizeof(d));
	d.c0 = out->curves[0];
	d.c1 = out->curves[1];
	d.d0 = out->curves[2];
	d.d1 = out->curves[3];
	out->surface = NULL;
	qaws_surface_create_coons(&d, &out->surface);
}

static void cb_destroy(cb_patch* p)
{
	int i;
	qaws_surface_destroy(p->surface);
	for (i = 0; i < 4; i++)
		qaws_curve_destroy(p->curves[i]);
}

static double cb_abs_mean(void const* user, qaws_scalar u, qaws_scalar v)
{
	qaws_surface_jet j;
	qaws_surface_geometry g;
	qaws_surface_eval_jet((qaws_surface const*)user, u, v, QAWS_SJET_ORDER2, &j);
	qaws_surface_geometry_eval(&j, NULL, NULL, &g, NULL, NULL, NULL);
	return fabs(g.mean);
}

/* Mean-curvature energy on a grid and its adjoint into boundary views. */
static double cb_energy(cb_patch const* patch, qaws_diff_views* views)
{
	double e = 0;
	int i, j;
	for (i = 0; i < CB_GRID; i++)
		for (j = 0; j < CB_GRID; j++)
		{
			qaws_scalar u = (qaws_scalar)((i + 0.5) / CB_GRID), v = (qaws_scalar)((j + 0.5) / CB_GRID);
			qaws_surface_jet p, ybar;
			qaws_surface_geometry g, gbar;
			qaws_surface_eval_jet(patch->surface, u, v, QAWS_SJET_ORDER2, &p);
			qaws_surface_geometry_eval(&p, NULL, NULL, &g, NULL, NULL, NULL);
			e += g.mean * g.mean / (CB_GRID * CB_GRID);
			if (views)
			{
				memset(&gbar, 0, sizeof(gbar));
				gbar.mean = (qaws_scalar)(2 * g.mean / (CB_GRID * CB_GRID));
				memset(&ybar, 0, sizeof(ybar));
				qaws_surface_geometry_adjoint(&p, &gbar, &ybar, NULL);
				qaws_surface_eval_adjoint(NULL, patch->surface, u, v, ybar.channels, &ybar, views, NULL, NULL);
			}
		}
	return e;
}

static void cb_draw_boundaries(svg* s, projection const* pr, cb_patch const* patch)
{
	int i, k;
	for (i = 0; i < 4; i++)
	{
		double xy[2 * 80];
		qaws_range r = qaws_curve_get_parameter_range(patch->curves[i]);
		for (k = 0; k < 80; k++)
		{
			qaws_eval_result_3d e;
			qaws_curve_evaluate_3d(patch->curves[i], r.min_value + (r.max_value - r.min_value) * k / (qaws_scalar)79,
				QAWS_EVAL_FLAG_POSITION, &e);
			project(pr, e.position.x, e.position.y, e.position.z, &xy[2 * k], &xy[2 * k + 1]);
		}
		svg_polyline(s, xy, 80, (i == 1 || i == 3) ? "#cf222e" : "#24292f", (i == 1 || i == 3) ? 3 : 2, 1, 0);
	}
}

static void demo_coons(void)
{
	qaws_scalar p[4][CB_CP * 3], initial[4][CB_CP * 3];
	qaws_scalar grad[4][CB_CP * 3];
	unsigned char interior[CB_CP] = { 0, 1, 1, 1, 1, 0 };
	double energy[CB_ITERS + 1];
	adam opt;
	int i, n, it;
	svg s;
	char buf[200];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < 4; i++)
		for (n = 0; n < CB_CP; n++)
		{
			qaws_scalar t = (qaws_scalar)n * (qaws_scalar)0.6;
			qaws_scalar* q = &p[i][3 * n];
			qaws_scalar amp = (i == 1 || i == 3) ? (qaws_scalar)0.7 : (qaws_scalar)0.12;
			qaws_scalar wave = (n == 0 || n == CB_CP - 1) ? (qaws_scalar)0 : (qaws_scalar)(amp * sin(2.3 * n + i));
			if (i == 0) { q[0] = t; q[1] = 0; }
			else if (i == 1) { q[0] = t; q[1] = 3; }
			else if (i == 2) { q[0] = 0; q[1] = t; }
			else { q[0] = 3; q[1] = t; }
			q[2] = wave;
		}
	memcpy(initial, p, sizeof(p));

	for (it = 0; it <= CB_ITERS; it++)
	{
		cb_patch patch;
		qaws_field_view fv[4];
		qaws_diff_views child[4], views;
		cb_build(p, &patch);
		memset(grad, 0, sizeof(grad));
		for (i = 0; i < 4; i++)
		{
			child[i] = one_field(&fv[i], QAWS_FIELD_CONTROL_POINTS, grad[i], CB_CP, 3);
			child[i].fields = &fv[i];
			fv[i].active = interior;          /* end points stay on the corners */
			fv[i].component_mask = 1u << 2;   /* only heights move */
		}
		views.fields = NULL;
		views.field_count = 0;
		views.children = child;
		views.child_count = 4;
		/* Only boundaries c1 and d1 (children 1 and 3) are editable. */
		child[0].field_count = 0;
		child[2].field_count = 0;
		energy[it] = cb_energy(&patch, &views);
		cb_destroy(&patch);
		if (it < CB_ITERS)
			adam_step(&opt, &p[0][0], &grad[0][0], 4 * CB_CP * 3, 0.01);
	}

	svg_open(&s, "showcase/7_coons.svg", 1200, 600, "Coons patch fairing by editing two boundaries",
		"Minimize mean(H^2) over the patch; adjoint: mean curvature -> geometry -> Coons rule -> child curves. Editable: red boundaries, interior heights only (masks).");
	{
		viewport a = { 30, 80, 560, 400, 0, 0, 0, 0 };
		viewport b = { 610, 80, 560, 400, 0, 0, 0, 0 };
		viewport lv = { 30, 495, 1140, 95, 0, 0, 0, 0 };
		projection pa = { 310, 250, 72, 1.6 }, pb = { 890, 250, 72, 1.6 };
		cb_patch before, after;
		cb_build(initial, &before);
		cb_build(p, &after);
		svg_panel(&s, &a, "before: |mean curvature|");
		draw_quads(&s, &pa, before.surface, 24, 24, cb_abs_mean, before.surface, 0, 0.8);
		cb_draw_boundaries(&s, &pa, &before);
		svg_panel(&s, &b, "after: two boundaries edited by the optimizer");
		draw_quads(&s, &pb, after.surface, 24, 24, cb_abs_mean, after.surface, 0, 0.8);
		cb_draw_boundaries(&s, &pb, &after);
		svg_colorbar(&s, b.x0 + b.w - 230, b.y0 + b.h - 26, 200, 8, "|H| = 0", "0.8");
		svg_loss_plot(&s, &lv, energy, CB_ITERS + 1, "#cf222e", "mean(H^2) (log)");
		cb_destroy(&before);
		cb_destroy(&after);
	}
	svg_close(&s);
	sprintf(buf, "%.4e -> %.4e", energy[0], energy[CB_ITERS]);
	printf("7_coons: mean(H^2) %s\n", buf);
}

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_curve_fit();
	demo_sensitivity();
	demo_fairing();
	demo_nurbs_weight();
	demo_surface();
	demo_vase();
	demo_coons();
	return 0;
}
