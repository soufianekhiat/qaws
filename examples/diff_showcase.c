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
 *   14_arc_length.svg      constant-speed samples: first and second order
 *                          tangents, and a fit by Newton-CG on exact HVPs
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

#include "example_svg.h"

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

/* ================================================================== */
/*  8. Curve networks: Gordon "affects" and loft "affected by"         */
/* ================================================================== */

static qaws_scalar const g_net_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_curve* net_curve(qaws_scalar const* p)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 6;
	d.knots = g_net_knots;
	d.knot_count = 10;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static double net_height(double x, double y)
{
	return 0.35 * sin(1.1 * x) * cos(0.9 * y) + 0.15 * x;
}

typedef struct affects_ctx
{
	qaws_surface const* surface;
	qaws_diff_views const* views;
} affects_ctx;

static double affects_value(void const* user, qaws_scalar u, qaws_scalar v)
{
	affects_ctx const* a = (affects_ctx const*)user;
	qaws_surface_jet p, t;
	qaws_surface_eval_tangent(NULL, a->surface, u, v, 0, 0, QAWS_SJET_P, a->views, &p, &t);
	return sqrt(t.d[0].x * t.d[0].x + t.d[0].y * t.d[0].y + t.d[0].z * t.d[0].z);
}

static void draw_curve3(svg* s, projection const* pr, qaws_curve const* c, char const* color, double width)
{
	double xy[2 * 80];
	qaws_range r = qaws_curve_get_parameter_range(c);
	int k;
	for (k = 0; k < 80; k++)
	{
		qaws_eval_result_3d e;
		qaws_curve_evaluate_3d(c, r.min_value + (r.max_value - r.min_value) * k / (qaws_scalar)79, QAWS_EVAL_FLAG_POSITION, &e);
		project(pr, e.position.x, e.position.y, e.position.z, &xy[2 * k], &xy[2 * k + 1]);
	}
	svg_polyline(s, xy, 80, color, width, 1, 0);
}

static void demo_networks(void)
{
	static qaws_scalar const params[3] = { 0, 0.5f, 1 };
	qaws_scalar net[6][18];
	qaws_curve* curves[6];
	qaws_curve const* uc[3];
	qaws_curve const* vc[3];
	qaws_surface* gordon = NULL;
	qaws_surface_gordon_desc gd;
	svg s;
	int i, n;
	unsigned int sel_curve = 4, sel_cp = 3;  /* middle v-curve, interior control point */
	viewport a = { 30, 80, 560, 500, 0, 0, 0, 0 };
	viewport b = { 610, 80, 560, 500, 0, 0, 0, 0 };
	projection pa = { 310, 230, 95, 1.6 };

	for (i = 0; i < 6; i++)
	{
		for (n = 0; n < 6; n++)
		{
			double sv = n * 0.6, fixed = (i % 3) * 1.5;
			qaws_scalar* p = &net[i][3 * n];
			if (i < 3) { p[0] = (qaws_scalar)sv; p[1] = (qaws_scalar)fixed; }
			else { p[0] = (qaws_scalar)fixed; p[1] = (qaws_scalar)sv; }
			p[2] = (qaws_scalar)net_height(p[0], p[1]);
		}
		curves[i] = net_curve(net[i]);
	}
	for (i = 0; i < 3; i++)
	{
		uc[i] = curves[i];
		vc[i] = curves[3 + i];
	}
	gd.u_curves = uc;
	gd.u_curve_count = 3;
	gd.v_params = params;
	gd.v_curves = vc;
	gd.v_curve_count = 3;
	gd.u_params = params;
	qaws_surface_create_gordon(&gd, &gordon);

	svg_open(&s, "showcase/8_networks.svg", 1200, 600, "Curve networks: Gordon affects / loft affected by",
		"Left: tangent of one control point of the middle v-curve -> |dS| over a Gordon surface. Right: adjoint of one loft point -> influence per section.");

	/* Gordon: affects map through child views (child 4 = middle v-curve). */
	svg_panel(&s, &a, "Gordon: where does the red control point act?");
	{
		qaws_scalar dir[6][18];
		qaws_field_view fv[6];
		qaws_diff_views child[6], views;
		affects_ctx actx;
		double sx, sy;
		memset(dir, 0, sizeof(dir));
		dir[sel_curve][3 * sel_cp + 2] = 1;  /* move up */
		for (i = 0; i < 6; i++)
			child[i] = one_field(&fv[i], QAWS_FIELD_CONTROL_POINTS, dir[i], 6, 3);
		for (i = 0; i < 6; i++)
			child[i].fields = &fv[i];
		views.fields = NULL;
		views.field_count = 0;
		views.children = child;
		views.child_count = 6;
		actx.surface = gordon;
		actx.views = &views;
		draw_quads(&s, &pa, gordon, 26, 26, affects_value, &actx, 0, 1);
		for (i = 0; i < 6; i++)
			draw_curve3(&s, &pa, curves[i], i == (int)sel_curve ? "#cf222e" : "#24292f", i == (int)sel_curve ? 3 : 1.6);
		project(&pa, net[sel_curve][3 * sel_cp], net[sel_curve][3 * sel_cp + 1], net[sel_curve][3 * sel_cp + 2], &sx, &sy);
		svg_line(&s, sx, sy, sx, sy - 45, "#cf222e", 2.5, 1);
		svg_circle(&s, sx, sy, 6, "#cf222e", "#ffffff");
		svg_colorbar(&s, a.x0 + 12, a.y0 + a.h - 24, 200, 8, "|dS| = 0", "1");
		svg_text(&s, a.x0 + 12, a.y0 + a.h - 34, 12, "#57606a", "start",
			"the edit spreads along the v-curve and fades with the Catmull-Rom blend in u");
	}

	/* Loft: affected-by for one surface point through the section curves. */
	svg_panel(&s, &b, "Loft: which sections move the marked point?");
	{
		qaws_scalar sec[5][18], bars[5][18];
		qaws_curve* sections[5];
		qaws_curve const* sp[5];
		qaws_surface* loft = NULL;
		qaws_surface_loft_desc ld;
		qaws_field_view fv[5];
		qaws_diff_views child[5], views;
		qaws_surface_jet ybar, p, t;
		qaws_scalar u0 = (qaws_scalar)0.5, v0 = (qaws_scalar)0.42;
		double infl[5], imax = 0, sx, sy;
		projection pb = { 890, 180, 60, 1.6 };
		char buf[64], col[32];

		for (i = 0; i < 5; i++)
		{
			for (n = 0; n < 6; n++)
			{
				qaws_scalar* q = &sec[i][3 * n];
				q[0] = (qaws_scalar)(n * 0.6);
				q[1] = (qaws_scalar)(i * 1.0);
				q[2] = (qaws_scalar)(0.5 * sin(1.3 * n + i) * (0.4 + 0.15 * i));
			}
			sections[i] = net_curve(sec[i]);
			sp[i] = sections[i];
		}
		ld.sections = sp;
		ld.section_count = 5;
		ld.v_parameters = NULL;
		qaws_surface_create_loft(&ld, &loft);

		memset(bars, 0, sizeof(bars));
		for (i = 0; i < 5; i++)
			child[i] = one_field(&fv[i], QAWS_FIELD_CONTROL_POINTS, bars[i], 6, 3);
		for (i = 0; i < 5; i++)
			child[i].fields = &fv[i];
		views.fields = NULL;
		views.field_count = 0;
		views.children = child;
		views.child_count = 5;
		memset(&ybar, 0, sizeof(ybar));
		ybar.d[0] = v3(0, 0, 1);  /* height of the marked point */
		qaws_surface_eval_adjoint(NULL, loft, u0, v0, QAWS_SJET_P, &ybar, &views, NULL, NULL);
		for (i = 0; i < 5; i++)
		{
			infl[i] = 0;
			for (n = 0; n < 6; n++)
				infl[i] += fabs(bars[i][3 * n + 2]);
			if (infl[i] > imax) imax = infl[i];
		}

		draw_quads(&s, &pb, loft, 26, 26, vase_height_value, NULL, 0, 1);
		for (i = 0; i < 5; i++)
		{
			heat(infl[i] / imax, col);
			draw_curve3(&s, &pb, sections[i], col, 2 + 4 * infl[i] / imax);
		}
		qaws_surface_eval_tangent(NULL, loft, u0, v0, 0, 0, QAWS_SJET_P, NULL, &p, &t);
		project(&pb, p.d[0].x, p.d[0].y, p.d[0].z, &sx, &sy);
		svg_circle(&s, sx, sy, 7, "#cf222e", "#ffffff");

		/* bar chart of influence per section (signed sum shown as magnitude) */
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"175\" fill=\"#ffffff\" stroke=\"#d0d7de\"/>\n", b.x0 + 10, b.y0 + b.h - 185, b.w - 20);
		svg_text(&s, b.x0 + 20, b.y0 + b.h - 168, 12, "#57606a", "start", "sum over each section of |d height / d control point z|");
		for (i = 0; i < 5; i++)
		{
			double bx = b.x0 + 40 + i * 100, bh = 95 * infl[i] / imax;
			heat(infl[i] / imax, col);
			fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"60\" height=\"%.1f\" fill=\"%s\"/>\n", bx, b.y0 + b.h - 40 - bh, bh, col);
			sprintf(buf, "section %d", i);
			svg_text(&s, bx + 30, b.y0 + b.h - 24, 11, "#24292f", "middle", buf);
			sprintf(buf, "%.2f", infl[i]);
			svg_text(&s, bx + 30, b.y0 + b.h - 46 - bh, 11, "#24292f", "middle", buf);
		}
		qaws_surface_destroy(loft);
		for (i = 0; i < 5; i++)
			qaws_curve_destroy(sections[i]);
	}

	svg_close(&s);
	qaws_surface_destroy(gordon);
	for (i = 0; i < 6; i++)
		qaws_curve_destroy(curves[i]);
	printf("8_networks: written\n");
}

/* ================================================================== */
/*  9. Projection-based fitting and the validity map of closest point */
/* ================================================================== */

#define PRJ_CP 10
#define PRJ_N 70
#define PRJ_ITERS 400

static void prj_target(int i, double* x, double* y)
{
	double s = i / (double)(PRJ_N - 1);
	*x = 0.7 + 4.4 * s;
	*y = 2.1 + 1.25 * sin(2 * PI * 1.25 * s) + 0.5 * s;
}

static void demo_projection(void)
{
	qaws_scalar cps[PRJ_CP * 2], grad[PRJ_CP * 2];
	qaws_vec3 targets[PRJ_N];
	double loss[PRJ_ITERS + 1];
	adam opt;
	int i, it;
	svg s;
	viewport a = { 30, 80, 560, 500, 0.2, 5.8, -0.3, 4.6 };
	viewport b = { 610, 80, 560, 500, 0.2, 5.8, -0.3, 4.6 };
	viewport lv = { 360, 470, 220, 100, 0, 0, 0, 0 };
	char buf[160];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < PRJ_N; i++)
	{
		double x, y;
		prj_target(i, &x, &y);
		targets[i] = v3((qaws_scalar)x, (qaws_scalar)y, 0);
	}
	for (i = 0; i < PRJ_CP; i++)
	{
		cps[2 * i] = (qaws_scalar)(1.0 + 3.6 * i / (PRJ_CP - 1));
		cps[2 * i + 1] = (qaws_scalar)2.1;
	}

	for (it = 0; it <= PRJ_ITERS; it++)
	{
		qaws_curve* c = bspline_2d(cps, PRJ_CP);
		qaws_field_view fv;
		qaws_diff_views views;
		double l = 0;
		memset(grad, 0, sizeof(grad));
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, PRJ_CP, 2);
		for (i = 0; i < PRJ_N; i++)
		{
			/* d^2 to the curve: no correspondences, the foot point is
			   re-solved and differentiated through its optimality. */
			qaws_curve_closest_point val, bar;
			qaws_curve_closest_point_tangent(NULL, c, targets[i], NULL, NULL, &val, NULL);
			l += val.distance * val.distance / PRJ_N;
			memset(&bar, 0, sizeof(bar));
			bar.distance = (qaws_scalar)(2 * val.distance / PRJ_N);
			qaws_curve_closest_point_adjoint(NULL, c, targets[i], &bar, &views, NULL);
		}
		loss[it] = l;
		qaws_curve_destroy(c);
		if (it < PRJ_ITERS)
			adam_step(&opt, cps, grad, PRJ_CP * 2, 0.02);
	}

	svg_open(&s, "showcase/9_projection.svg", 1200, 600, "Closest point through the implicit function theorem",
		"Left: fit without correspondences, minimizing squared point-to-curve distance (closest-point adjoint). Right: validity reported for a grid of queries.");

	svg_panel(&s, &a, "fit by distance to the curve");
	{
		qaws_curve* c = bspline_2d(cps, PRJ_CP);
		double xy[2 * 300];
		curve_polyline(c, &a, xy, 300);
		for (i = 0; i < PRJ_N; i++)
		{
			qaws_curve_closest_point val;
			qaws_curve_closest_point_tangent(NULL, c, targets[i], NULL, NULL, &val, NULL);
			svg_line(&s, vx(&a, targets[i].x), vy(&a, targets[i].y), vx(&a, val.position.x), vy(&a, val.position.y), "#cf222e", 1, 0.7);
			svg_circle(&s, vx(&a, targets[i].x), vy(&a, targets[i].y), 2.6, "#d0d7de", "#57606a");
		}
		svg_polyline(&s, xy, 300, "#0969da", 2.6, 1, 0);
		control_polygon(&s, &a, cps, PRJ_CP, "#0969da");
		{
			qaws_scalar init[PRJ_CP * 2];
			qaws_curve* c0;
			for (i = 0; i < PRJ_CP; i++)
			{
				init[2 * i] = (qaws_scalar)(1.0 + 3.6 * i / (PRJ_CP - 1));
				init[2 * i + 1] = (qaws_scalar)2.1;
			}
			c0 = bspline_2d(init, PRJ_CP);
			curve_polyline(c0, &a, xy, 300);
			svg_polyline(&s, xy, 300, "#8c959f", 1.6, 1, 1);
			qaws_curve_destroy(c0);
		}
		svg_loss_plot(&s, &lv, loss, PRJ_ITERS + 1, "#cf222e", "mean d^2");
		qaws_curve_destroy(c);
	}

	/* Validity map: the medial axis shows up as the ambiguous region. */
	svg_panel(&s, &b, "validity of d(foot point)/d(query) over the plane");
	{
		qaws_curve* c = bspline_2d(cps, PRJ_CP);
		double xy[2 * 300];
		int gx, gy, nx = 112, ny = 100, counts[3] = { 0, 0, 0 };
		double cw = b.w / nx, ch = b.h / ny;
		for (gx = 0; gx < nx; gx++)
			for (gy = 0; gy < ny; gy++)
			{
				qaws_diff_context ctx;
				qaws_diff_report report;
				qaws_curve_closest_point val;
				double wx = b.xmin + (b.xmax - b.xmin) * (gx + 0.5) / nx;
				double wy = b.ymax - (b.ymax - b.ymin) * (gy + 0.5) / ny;
				char col[40];
				double rel;
				qaws_diff_context_init(&ctx);
				qaws_diff_report_reset(&report);
				ctx.report = &report;
				qaws_curve_closest_point_tangent(&ctx, c, v3((qaws_scalar)wx, (qaws_scalar)wy, 0), NULL, NULL, &val, NULL);
				/* Shade by the reported gap to the best competing foot point,
				   relative to the distance: the medial axis is where it vanishes. */
				rel = report.branch_gap / (val.distance + 1e-9);
				if (rel > 0.3) rel = 0.3;
				rel /= 0.3;
				if (report.validity == QAWS_DIFF_AMBIGUOUS || report.validity == QAWS_DIFF_ILL_CONDITIONED)
				{
					sprintf(col, "#cf222e");
					counts[1]++;
				}
				else if (report.validity == QAWS_DIFF_VALID_LOCALLY)
				{
					sprintf(col, "rgb(%d,%d,%d)", (int)(150 + 70 * rel), (int)(200 + 40 * rel), 255);
					counts[2]++;
				}
				else
				{
					sprintf(col, "rgb(%d,%d,%d)", (int)(235 - 17 * rel), (int)(120 + 131 * rel), (int)(110 + 115 * rel));
					counts[0]++;
				}
				fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\"/>\n",
					b.x0 + gx * cw, b.y0 + gy * ch, cw + 0.3, ch + 0.3, col);
			}
		curve_polyline(c, &b, xy, 300);
		svg_polyline(&s, xy, 300, "#0969da", 2.6, 1, 0);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"330\" height=\"64\" fill=\"#ffffff\" stroke=\"#d0d7de\"/>\n",
			b.x0 + 10, b.y0 + b.h - 74);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"12\" height=\"12\" fill=\"#dafbe1\" stroke=\"#8c959f\"/>\n", b.x0 + 20, b.y0 + b.h - 66);
		svg_text(&s, b.x0 + 38, b.y0 + b.h - 56, 12, "#24292f", "start", "valid; pinker = smaller gap to a competing foot point");
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"12\" height=\"12\" fill=\"#b6e3ff\" stroke=\"#8c959f\"/>\n", b.x0 + 20, b.y0 + b.h - 48);
		svg_text(&s, b.x0 + 38, b.y0 + b.h - 38, 12, "#24292f", "start", "valid locally: endpoint (active set frozen)");
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"12\" height=\"12\" fill=\"#cf222e\"/>\n", b.x0 + 20, b.y0 + b.h - 30);
		svg_text(&s, b.x0 + 38, b.y0 + b.h - 20, 12, "#24292f", "start", "ambiguous: gap &lt; 0.5% of the distance (medial axis)");
		qaws_curve_destroy(c);
		sprintf(buf, "%d valid, %d ambiguous, %d endpoint", counts[0], counts[1], counts[2]);
		printf("9_projection: fit %.4e -> %.4e, map: %s\n", loss[0], loss[PRJ_ITERS], buf);
	}
	svg_close(&s);
}

/* ================================================================== */
/*  10. Soap film: area minimization with gradient descent vs Newton  */
/*      steps built from Hessian-vector products                      */
/* ================================================================== */

#define SF_N 10
#define SF_CP (SF_N * SF_N)
#define SF_GD_ITERS 60
#define SF_NEWTON_ITERS 12

static qaws_surface* film_surface(qaws_scalar const* cps)
{
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = SF_N;
	d.v_point_count = SF_N;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

static void film_init(qaws_scalar* cps, unsigned char* active)
{
	int i, j;
	for (i = 0; i < SF_N; i++)
		for (j = 0; j < SF_N; j++)
		{
			double x = i / (double)(SF_N - 1), y = j / (double)(SF_N - 1);
			qaws_scalar* p = &cps[(i * SF_N + j) * 3];
			int boundary = i == 0 || j == 0 || i == SF_N - 1 || j == SF_N - 1;
			p[0] = (qaws_scalar)x;
			p[1] = (qaws_scalar)y;
			/* wavy frame; interior starts as a tilted bump far from minimal */
			p[2] = boundary ? (qaws_scalar)(0.3 * sin(2 * PI * x) * cos(PI * y) + 0.2 * sin(3 * PI * y) * (x - 0.5))
			                : (qaws_scalar)(0.45 * sin(PI * x) * sin(PI * y));
			active[i * SF_N + j] = (unsigned char)!boundary;
		}
}

static double film_area(qaws_scalar const* cps)
{
	qaws_surface* s = film_surface(cps);
	qaws_scalar a = 0;
	qaws_surface_functional_eval(NULL, s, QAWS_FUNCTIONAL_AREA, 0, NULL, &a, NULL, NULL);
	qaws_surface_destroy(s);
	return a;
}

/* Masked view: interior control points, z component only. */
static qaws_diff_views film_view(qaws_field_view* fv, qaws_scalar* data, unsigned char const* active)
{
	qaws_diff_views v = one_field(fv, QAWS_FIELD_CONTROL_POINTS, data, SF_CP, 3);
	fv->active = active;
	fv->component_mask = 1u << 2;
	return v;
}

static void film_gradient(qaws_scalar const* cps, unsigned char const* active, qaws_scalar* g)
{
	qaws_surface* s = film_surface(cps);
	qaws_field_view fv;
	qaws_diff_views v;
	memset(g, 0, sizeof(qaws_scalar) * SF_CP * 3);
	v = film_view(&fv, g, active);
	qaws_surface_functional_gradient(NULL, s, QAWS_FUNCTIONAL_AREA, 0, &v, NULL);
	qaws_surface_destroy(s);
}

static void film_hvp(qaws_surface const* s, unsigned char const* active, qaws_scalar* dir, qaws_scalar* out)
{
	qaws_field_view fd, fo;
	qaws_diff_views vd = film_view(&fd, dir, active), vo;
	memset(out, 0, sizeof(qaws_scalar) * SF_CP * 3);
	vo = film_view(&fo, out, active);
	qaws_surface_functional_hvp(NULL, s, QAWS_FUNCTIONAL_AREA, 0, &vd, &vo);
}

static double vdot(qaws_scalar const* a, qaws_scalar const* b, int n)
{
	double s = 0;
	int i;
	for (i = 0; i < n; i++)
		s += (double)a[i] * b[i];
	return s;
}

static double film_mean_abs(void const* user, qaws_scalar u, qaws_scalar v)
{
	qaws_surface_jet j;
	qaws_surface_geometry g;
	qaws_surface_eval_jet((qaws_surface const*)user, u, v, QAWS_SJET_ORDER2, &j);
	qaws_surface_geometry_eval(&j, NULL, NULL, &g, NULL, NULL, NULL);
	return fabs(g.mean);
}

static void demo_soap_film(void)
{
	enum { N3 = SF_CP * 3 };
	qaws_scalar init[N3], gd[N3], nt[N3], g[N3];
	unsigned char active[SF_CP];
	double e_gd[SF_GD_ITERS + 1], e_nt[SF_GD_ITERS + 1];
	int it, i, hvp_calls = 0;
	svg s;
	char buf[220];

	film_init(init, active);
	memcpy(gd, init, sizeof(init));
	memcpy(nt, init, sizeof(init));

	/* Plain gradient descent with a fixed step. */
	for (it = 0; it <= SF_GD_ITERS; it++)
	{
		e_gd[it] = film_area(gd);
		if (it == SF_GD_ITERS)
			break;
		film_gradient(gd, active, g);
		for (i = 0; i < N3; i++)
			gd[i] -= (qaws_scalar)(1.5 * g[i]);
	}

	/* Newton-CG: solve H p = -g with conjugate gradients on HVPs, then a
	   backtracking line search on the area. */
	for (it = 0; it <= SF_NEWTON_ITERS; it++)
	{
		qaws_scalar p[N3], r[N3], dd[N3], hd[N3];
		qaws_surface* sf;
		double rr, e0;
		int k;
		e_nt[it] = e0 = film_area(nt);
		if (it == SF_NEWTON_ITERS)
			break;
		film_gradient(nt, active, g);
		sf = film_surface(nt);
		memset(p, 0, sizeof(p));
		for (i = 0; i < N3; i++)
			r[i] = -g[i];
		memcpy(dd, r, sizeof(r));
		rr = vdot(r, r, N3);
		for (k = 0; k < 25 && rr > 1e-24; k++)
		{
			double dhd, alpha, rr_new;
			film_hvp(sf, active, dd, hd);
			hvp_calls++;
			dhd = vdot(dd, hd, N3);
			if (dhd <= 0)
			{
				if (k == 0)
					memcpy(p, r, sizeof(p));
				break;
			}
			alpha = rr / dhd;
			for (i = 0; i < N3; i++)
			{
				p[i] += (qaws_scalar)(alpha * dd[i]);
				r[i] -= (qaws_scalar)(alpha * hd[i]);
			}
			rr_new = vdot(r, r, N3);
			for (i = 0; i < N3; i++)
				dd[i] = (qaws_scalar)(r[i] + rr_new / rr * dd[i]);
			rr = rr_new;
		}
		qaws_surface_destroy(sf);
		{
			double step = 1;
			qaws_scalar trial[N3];
			for (k = 0; k < 20; k++)
			{
				for (i = 0; i < N3; i++)
					trial[i] = (qaws_scalar)(nt[i] + step * p[i]);
				if (film_area(trial) < e0)
					break;
				step *= 0.5;
			}
			if (k < 20)
				memcpy(nt, trial, sizeof(trial));
		}
	}

	svg_open(&s, "showcase/10_soap_film.svg", 1200, 600, "Soap film: area minimization with Hessian-vector products",
		"Bicubic B-spline, boundary fixed, interior heights free (masks). Newton steps solve H p = -g by conjugate gradients on the library's direct HVPs.");
	{
		viewport a = { 30, 80, 370, 500, 0, 0, 0, 0 };
		viewport b = { 415, 80, 370, 500, 0, 0, 0, 0 };
		viewport lv = { 800, 80, 370, 300, 0, 0, 0, 0 };
		projection pa = { 215, 300, 200, 1.0 }, pb = { 600, 300, 200, 1.0 };
		qaws_surface* s0 = film_surface(init);
		qaws_surface* s1 = film_surface(nt);
		double emin = e_nt[SF_NEWTON_ITERS], xy0[2 * (SF_GD_ITERS + 1)], xy1[2 * (SF_NEWTON_ITERS + 1)];
		double lo = 1e300, hi = -1e300;

		svg_panel(&s, &a, "initial interior: |mean curvature|");
		draw_quads(&s, &pa, s0, 26, 26, film_mean_abs, s0, 0, 1.5);
		svg_panel(&s, &b, "after Newton: mean curvature driven toward 0");
		draw_quads(&s, &pb, s1, 26, 26, film_mean_abs, s1, 0, 1.5);
		svg_colorbar(&s, b.x0 + 12, b.y0 + b.h - 24, 200, 8, "|H| = 0", "1.5");

		/* Convergence: log10(E - E_min) per iteration. */
		svg_panel(&s, &lv, "log10(area - final area) per iteration");
		for (i = 0; i <= SF_GD_ITERS; i++)
		{
			double l = log10(e_gd[i] - emin + 1e-16);
			if (l < lo) lo = l;
			if (l > hi) hi = l;
		}
		for (i = 0; i <= SF_NEWTON_ITERS; i++)
		{
			double l = log10(e_nt[i] - emin + 1e-16);
			if (l < lo) lo = l;
			if (l > hi) hi = l;
		}
		if (lo < -14) lo = -14;
		for (i = 0; i <= SF_GD_ITERS; i++)
		{
			double l = log10(e_gd[i] - emin + 1e-16);
			if (l < lo) l = lo;
			xy0[2 * i] = lv.x0 + 14 + (lv.w - 28) * i / (double)SF_GD_ITERS;
			xy0[2 * i + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - l) / (hi - lo);
		}
		for (i = 0; i <= SF_NEWTON_ITERS; i++)
		{
			double l = log10(e_nt[i] - emin + 1e-16);
			if (l < lo) l = lo;
			xy1[2 * i] = lv.x0 + 14 + (lv.w - 28) * i / (double)SF_GD_ITERS;
			xy1[2 * i + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - l) / (hi - lo);
		}
		svg_polyline(&s, xy0, SF_GD_ITERS + 1, "#8c959f", 2, 1, 0);
		svg_polyline(&s, xy1, SF_NEWTON_ITERS + 1, "#0969da", 2.6, 1, 0);
		for (i = 0; i <= SF_NEWTON_ITERS; i++)
			svg_circle(&s, xy1[2 * i], xy1[2 * i + 1], 3, "#0969da", "#ffffff");
		svg_text(&s, lv.x0 + lv.w - 14, xy0[2 * SF_GD_ITERS + 1] - 8, 12, "#57606a", "end", "gradient descent");
		svg_text(&s, lv.x0 + 120, lv.y0 + 60, 12, "#0969da", "start", "Newton-CG on HVPs");

		sprintf(buf, "initial area %.6f", e_gd[0]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 26, 13, "#24292f", "start", buf);
		sprintf(buf, "gradient descent, %d steps: %.8f", SF_GD_ITERS, e_gd[SF_GD_ITERS]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 48, 13, "#57606a", "start", buf);
		sprintf(buf, "Newton-CG, %d steps (%d HVPs): %.8f", SF_NEWTON_ITERS, hvp_calls, e_nt[SF_NEWTON_ITERS]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 70, 13, "#0969da", "start", buf);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 98, 12, "#57606a", "start", "API: qaws_surface_functional_gradient / _hvp");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 116, 12, "#57606a", "start", "masks: interior elements, z component");
		qaws_surface_destroy(s0);
		qaws_surface_destroy(s1);
	}
	svg_close(&s);
	printf("10_soap_film: area %.8f -> GD %.8f, Newton %.8f (%d HVPs)\n", e_gd[0], e_gd[SF_GD_ITERS], e_nt[SF_NEWTON_ITERS], hvp_calls);
}


/* ================================================================== */
/*  11. Parameter correction: optimizing sample parameters through     */
/*      the differentiable least-squares fit                          */
/* ================================================================== */

#define PC_M 36
#define PC_N 9
#define PC_KC (PC_N + 3 + 1)
#define PC_ITERS 300
#define PC_LAMBDA 2e-7
#define PC_SWITCH 30

static void pc_data(qaws_scalar* d)
{
	int i;
	for (i = 0; i < PC_M; i++)
	{
		/* samples crowd toward the start of the arc */
		double w = pow(i / (double)(PC_M - 1), 1.9);
		double s = 1.55 * PI * w;
		double r = 1.0 + 0.28 * cos(3.0 * s);
		d[2 * i] = (qaws_scalar)(r * cos(s));
		d[2 * i + 1] = (qaws_scalar)(r * sin(s));
	}
}

static qaws_curve* pc_fit(qaws_scalar const* data, qaws_scalar const* t, qaws_diff_map** map)
{
	qaws_bspline_fit_desc desc;
	qaws_curve* c = NULL;
	memset(&desc, 0, sizeof(desc));
	desc.dimension = QAWS_DIMENSION_2D;
	desc.data_points = data;
	desc.data_point_count = PC_M;
	desc.degree = 3;
	desc.control_point_count = PC_N;
	desc.parameters = t;
	if (map)
		qaws_curve_fit_bspline_diff(&desc, &c, map);
	else
		qaws_curve_fit_bspline(&desc, &c);
	return c;
}

/* Mean squared distance between data and the fitted curve at the samples. */
static double pc_error(qaws_curve const* c, qaws_scalar const* data, qaws_scalar const* t)
{
	double e = 0;
	int k;
	for (k = 0; k < PC_M; k++)
	{
		qaws_eval_result_2d r;
		double dx, dy;
		qaws_curve_evaluate_2d(c, t[k], QAWS_EVAL_FLAG_POSITION, &r);
		dx = r.position.x - data[2 * k];
		dy = r.position.y - data[2 * k + 1];
		e += dx * dx + dy * dy;
	}
	return e / PC_M;
}

static double pc_error_of(qaws_scalar const* data, qaws_scalar const* t)
{
	qaws_curve* c = pc_fit(data, t, NULL);
	double e = pc_error(c, data, t);
	qaws_curve_destroy(c);
	return e;
}

/* Keep parameters strictly increasing with fixed ends. */
static void pc_project(qaws_scalar* t)
{
	int k;
	t[0] = 0;
	t[PC_M - 1] = 1;
	for (k = 1; k < PC_M - 1; k++)
	{
		if (t[k] < t[k - 1] + (qaws_scalar)1e-4) t[k] = t[k - 1] + (qaws_scalar)1e-4;
		if (t[k] > (qaws_scalar)(1 - 1e-4 * (PC_M - 1 - k))) t[k] = (qaws_scalar)(1 - 1e-4 * (PC_M - 1 - k));
	}
}

/* Exact gradient of the fit error with respect to the sample parameters:
   the samples move along the curve (t adjoint) and the curve itself is
   refit (control points and knots pulled back through the fit map). */
static double pc_gradient(qaws_scalar const* data, qaws_scalar const* t, qaws_scalar* g, double* out_bending)
{
	qaws_diff_map* map = NULL;
	qaws_curve* c = pc_fit(data, t, &map);
	qaws_curve_jet_2d ybar[PC_M];
	qaws_scalar cbar[PC_N * 2], kbar[PC_KC], tbar[PC_M], pbar[PC_M];
	qaws_field_view out_f[2], in_f;
	qaws_diff_views out_v, in_v;
	qaws_diff_views const* outs[1];
	qaws_diff_views* ins[2];
	double e = 0;
	int k;

	memset(ybar, 0, sizeof(ybar));
	for (k = 0; k < PC_M; k++)
	{
		qaws_eval_result_2d r;
		double dx, dy;
		qaws_curve_evaluate_2d(c, t[k], QAWS_EVAL_FLAG_POSITION, &r);
		dx = r.position.x - data[2 * k];
		dy = r.position.y - data[2 * k + 1];
		e += dx * dx + dy * dy;
		ybar[k].d[0].x = (qaws_scalar)(2 * dx / PC_M);
		ybar[k].d[0].y = (qaws_scalar)(2 * dy / PC_M);
		ybar[k].channels = 1;
	}
	memset(cbar, 0, sizeof(cbar));
	memset(kbar, 0, sizeof(kbar));
	memset(tbar, 0, sizeof(tbar));
	memset(pbar, 0, sizeof(pbar));
	out_f[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, cbar, PC_N, 2);
	out_f[1] = qaws_field_view_make(QAWS_FIELD_KNOTS, kbar, PC_KC, 1);
	out_v.fields = out_f;
	out_v.field_count = 2;
	out_v.children = NULL;
	out_v.child_count = 0;
	qaws_curve_eval_batch_adjoint_2d(NULL, c, t, PC_M, 1, ybar, &out_v, tbar);

	/* Fairness keeps the refit from looping between sparse samples: the
	   bending gradient lands on control points and knots alike. */
	{
		qaws_scalar cb[PC_N * 2], kb[PC_KC], bend = 0;
		qaws_field_view bf[2];
		qaws_diff_views bv;
		int i;
		memset(cb, 0, sizeof(cb));
		memset(kb, 0, sizeof(kb));
		bf[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, cb, PC_N, 2);
		bf[1] = qaws_field_view_make(QAWS_FIELD_KNOTS, kb, PC_KC, 1);
		bv = out_v;
		bv.fields = bf;
		qaws_curve_functional_gradient(NULL, c, QAWS_FUNCTIONAL_BENDING, 0, &bv, &bend);
		for (i = 0; i < PC_N * 2; i++)
			cbar[i] += (qaws_scalar)(PC_LAMBDA * cb[i]);
		for (i = 0; i < PC_KC; i++)
			kbar[i] += (qaws_scalar)(PC_LAMBDA * kb[i]);
		if (out_bending)
			*out_bending = bend;
	}

	in_v = one_field(&in_f, QAWS_FIELD_PARAMETER, pbar, PC_M, 1);
	outs[0] = &out_v;
	ins[0] = NULL;
	ins[1] = &in_v;
	qaws_diff_map_adjoint(map, NULL, outs, 1, ins, 2);

	for (k = 0; k < PC_M; k++)
		g[k] = (k == 0 || k == PC_M - 1) ? 0 : tbar[k] + pbar[k];
	qaws_diff_map_destroy(map);
	qaws_curve_destroy(c);
	return e / PC_M;
}

/* Classic parameter correction: project every data point onto the current
   fit (Newton on the foot point), then refit. */
static void pc_hoschek_step(qaws_scalar const* data, qaws_scalar* t)
{
	qaws_curve* c = pc_fit(data, t, NULL);
	int k, it;
	for (k = 1; k < PC_M - 1; k++)
		for (it = 0; it < 4; it++)
		{
			qaws_eval_result_2d r;
			double dx, dy, num, den;
			qaws_curve_evaluate_2d(c, t[k], QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &r);
			dx = r.position.x - data[2 * k];
			dy = r.position.y - data[2 * k + 1];
			num = dx * r.d1.x + dy * r.d1.y;
			den = r.d1.x * r.d1.x + r.d1.y * r.d1.y + dx * r.d2.x + dy * r.d2.y;
			if (den > 1e-12)
			{
				/* safeguarded: stay between the neighbours, accept only descent */
				double lo = t[k - 1] + 1e-4, hi = t[k + 1] - 1e-4, step = -num / den;
				int h;
				for (h = 0; h < 6; h++, step *= 0.5)
				{
					qaws_eval_result_2d q;
					double tn = t[k] + step, ex, ey;
					if (tn < lo) tn = lo;
					if (tn > hi) tn = hi;
					qaws_curve_evaluate_2d(c, (qaws_scalar)tn, QAWS_EVAL_FLAG_POSITION, &q);
					ex = q.position.x - data[2 * k];
					ey = q.position.y - data[2 * k + 1];
					if (ex * ex + ey * ey < dx * dx + dy * dy)
					{
						t[k] = (qaws_scalar)tn;
						break;
					}
				}
			}
		}
	pc_project(t);
	qaws_curve_destroy(c);
}

static void pc_panel(svg* s, viewport const* v, char const* label, qaws_scalar const* data, qaws_scalar const* t,
	char const* color, char* buf)
{
	qaws_curve* c = pc_fit(data, t, NULL);
	double xy[2 * 400];
	int k;
	svg_panel(s, v, label);
	curve_polyline(c, v, xy, 400);
	svg_polyline(s, xy, 400, color, 2.4, 1, 0);
	for (k = 0; k < PC_M; k++)
	{
		qaws_eval_result_2d r;
		qaws_curve_evaluate_2d(c, t[k], QAWS_EVAL_FLAG_POSITION, &r);
		svg_line(s, vx(v, data[2 * k]), vy(v, data[2 * k + 1]), vx(v, r.position.x), vy(v, r.position.y), "#cf222e", 1.6, 0.9);
		svg_circle(s, vx(v, r.position.x), vy(v, r.position.y), 2.2, color, color);
	}
	for (k = 0; k < PC_M; k++)
		svg_circle(s, vx(v, data[2 * k]), vy(v, data[2 * k + 1]), 3.2, "#ffffff", "#24292f");
	sprintf(buf, "RMS error %.2e", sqrt(pc_error(c, data, t)));
	svg_text(s, v->x0 + v->w - 10, v->y0 + v->h - 10, 13, "#24292f", "end", buf);
	qaws_curve_destroy(c);
}

/* Fit error plus weighted bending: the objective foot-point correction
   cannot see, since it only moves samples to their closest points. */
static double pc_objective(qaws_scalar const* data, qaws_scalar const* t, double* mse, double* bending)
{
	qaws_curve* c = pc_fit(data, t, NULL);
	qaws_scalar b = 0;
	qaws_curve_functional_eval(NULL, c, QAWS_FUNCTIONAL_BENDING, 0, NULL, &b, NULL, NULL);
	*mse = pc_error(c, data, t);
	*bending = b;
	qaws_curve_destroy(c);
	return *mse + PC_LAMBDA * b;
}

static void demo_parameter_correction(void)
{
	qaws_scalar data[2 * PC_M], tc[PC_M], th[PC_M], to[PC_M], g[PC_M];
	double j_opt[PC_ITERS + 1], j_hos[PC_ITERS + 1], j_chord, mse, bend;
	double mse_c, bend_c, mse_h, bend_h, mse_o, bend_o;
	adam opt;
	int k, it;
	svg s;
	char buf[256];

	pc_data(data);
	{
		/* chord length, as qaws_curve_fit_bspline does without parameters */
		double total = 0, acc = 0;
		for (k = 1; k < PC_M; k++)
			total += hypot(data[2 * k] - data[2 * k - 2], data[2 * k + 1] - data[2 * k - 1]);
		tc[0] = 0;
		for (k = 1; k < PC_M; k++)
		{
			acc += hypot(data[2 * k] - data[2 * k - 2], data[2 * k + 1] - data[2 * k - 1]);
			tc[k] = (qaws_scalar)(acc / total);
		}
		tc[PC_M - 1] = 1;
	}
	j_chord = pc_objective(data, tc, &mse_c, &bend_c);

	/* Both start from chord length. Foot-point correction keeps projecting;
	   from PC_SWITCH on, the exact gradient of the objective takes over. */
	memcpy(th, tc, sizeof(tc));
	memset(&opt, 0, sizeof(opt));
	for (it = 0; it <= PC_ITERS; it++)
	{
		j_hos[it] = pc_objective(data, th, &mse, &bend);
		if (it <= PC_SWITCH)
		{
			memcpy(to, th, sizeof(th));
			j_opt[it] = j_hos[it];
		}
		else
		{
			mse = pc_gradient(data, to, g, &bend);
			j_opt[it] = mse + PC_LAMBDA * bend;
		}
		if (it == PC_ITERS)
			break;
		if (it >= PC_SWITCH)
		{
			if (it == PC_SWITCH)
				pc_gradient(data, to, g, &bend);
			adam_step(&opt, to, g, PC_M, 1.5e-3 * (1.0 - 0.9 * (it - PC_SWITCH) / (double)(PC_ITERS - PC_SWITCH)));
			pc_project(to);
		}
		pc_hoschek_step(data, th);
	}
	pc_objective(data, th, &mse_h, &bend_h);
	pc_objective(data, to, &mse_o, &bend_o);

	svg_open(&s, "showcase/11_parameter_correction.svg", 1500, 560,
		"Parameter correction through the differentiable least-squares fit",
		"Sample parameters minimize fit error + fairness. The gradient runs through the whole fit: chord, knot averaging, "
		"basis and normal equations (qaws_curve_fit_bspline_diff) plus knot derivatives.");
	{
		viewport a = { 24, 80, 350, 440, -1.45, 1.45, -1.55, 1.95 };
		viewport b = a, c = a, lv = { 1140, 80, 340, 300, 0, 0, 0, 0 };
		double lo = 1e300, hi = -1e300, xy0[2 * (PC_ITERS + 1)], xy1[2 * (PC_ITERS + 1)], ys;
		b.x0 = 390;
		c.x0 = 756;
		pc_panel(&s, &a, "chord length (fit default)", data, tc, "#9a6700", buf);
		pc_panel(&s, &b, "foot-point correction", data, th, "#8c959f", buf);
		pc_panel(&s, &c, "exact gradient through the fit map", data, to, "#0969da", buf);
		sprintf(buf, "bending %.0f", bend_c);
		svg_text(&s, a.x0 + a.w - 10, a.y0 + a.h - 28, 12, "#57606a", "end", buf);
		sprintf(buf, "bending %.0f", bend_h);
		svg_text(&s, b.x0 + b.w - 10, b.y0 + b.h - 28, 12, "#57606a", "end", buf);
		sprintf(buf, "bending %.0f", bend_o);
		svg_text(&s, c.x0 + c.w - 10, c.y0 + c.h - 28, 12, "#57606a", "end", buf);

		for (it = 0; it <= PC_ITERS; it++)
		{
			double l0 = log10(j_opt[it]), l1 = log10(j_hos[it]);
			if (l0 < lo) lo = l0;
			if (l1 < lo) lo = l1;
			if (l0 > hi) hi = l0;
			if (l1 > hi) hi = l1;
		}
		lo -= 0.05;
		hi += 0.05;
		svg_panel(&s, &lv, "log10 objective per iteration");
		for (it = 0; it <= PC_ITERS; it++)
		{
			xy0[2 * it] = lv.x0 + 14 + (lv.w - 28) * it / (double)PC_ITERS;
			xy0[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(j_opt[it])) / (hi - lo);
			xy1[2 * it] = xy0[2 * it];
			xy1[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(j_hos[it])) / (hi - lo);
		}
		ys = lv.x0 + 14 + (lv.w - 28) * PC_SWITCH / (double)PC_ITERS;
		svg_line(&s, ys, lv.y0 + 30, ys, lv.y0 + lv.h - 12, "#d0d7de", 1.2, 1);
		svg_text(&s, ys + 4, lv.y0 + lv.h - 14, 11, "#57606a", "start", "gradient starts");
		svg_polyline(&s, xy1, PC_ITERS + 1, "#8c959f", 2, 1, 0);
		svg_polyline(&s, xy0, PC_ITERS + 1, "#0969da", 2.6, 1, 0);
		svg_text(&s, lv.x0 + lv.w - 14, xy1[2 * PC_ITERS + 1] - 8, 11, "#57606a", "end", "foot-point correction");
		svg_text(&s, lv.x0 + lv.w - 14, xy0[2 * PC_ITERS + 1] + 16, 11, "#0969da", "end", "exact gradient");

		sprintf(buf, "objective = mean squared error + %.0e x bending", PC_LAMBDA);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 26, 13, "#24292f", "start", buf);
		sprintf(buf, "chord %.3e   foot-point %.3e", j_chord, j_hos[PC_ITERS]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 48, 13, "#57606a", "start", buf);
		sprintf(buf, "exact gradient %.3e", j_opt[PC_ITERS]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 70, 13, "#0969da", "start", buf);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 98, 12, "#57606a", "start", "red: data to fitted sample, dots: C(t_k)");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 116, 12, "#57606a", "start", "gradient = t adjoint + fit map adjoint");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 134, 12, "#57606a", "start", "of control points and knots");
	}
	svg_close(&s);
	printf("11_parameter_correction: objective chord %.4e (rms %.3e, bend %.0f), foot-point %.4e (rms %.3e, bend %.0f), "
		"gradient %.4e (rms %.3e, bend %.0f)\n", j_chord, sqrt(mse_c), bend_c, j_hos[PC_ITERS], sqrt(mse_h), bend_h,
		j_opt[PC_ITERS], sqrt(mse_o), bend_o);
}


/* ================================================================== */
/*  12. Arch design from implicit features: intersections with a deck */
/*      and the apex (extremum) placed by adjoints of their solves     */
/* ================================================================== */

#define AR_N 8
#define AR_ITERS 160
#define AR_DECK 1.0
#define AR_APEX_X 3.0
#define AR_APEX_Y 2.4
#define AR_X0 1.1
#define AR_X1 4.9

static qaws_curve* arch_curve(qaws_scalar const* cps)
{
	return bspline_2d(cps, AR_N);
}

static qaws_curve* arch_deck(void)
{
	static qaws_scalar const d[8] = { -0.5, AR_DECK, 2.0, AR_DECK, 4.0, AR_DECK, 6.5, AR_DECK };
	return bspline_2d(d, 4);
}

typedef struct arch_features
{
	int ok;
	qaws_scalar t_apex, t_cross[2], tb_cross[2];
	qaws_vec3 apex, cross[2];
} arch_features;

/* Finds the features with the discrete finders (the active set), then
   accumulates the gradient of the design loss through the adjoints of
   the implicit solves. Returns the loss. */
static double arch_loss(qaws_curve const* arch, qaws_curve const* deck, qaws_scalar* grad, arch_features* f)
{
	qaws_intersection_2d hits[4];
	qaws_scalar ext[8];
	unsigned int nh = 0, ne = 0, i;
	double loss = 0;
	qaws_field_view fv;
	qaws_diff_views gv = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, AR_N, 2);
	qaws_vec3 up = v3(0, 1, 0);

	memset(f, 0, sizeof(*f));
	if (grad)
		memset(grad, 0, sizeof(qaws_scalar) * AR_N * 2);
	qaws_curve_find_intersections_2d(arch, deck, hits, 4, &nh);
	qaws_curve_find_extrema(arch, 1, ext, 8, &ne);
	if (nh < 2 || ne < 1)
		return 1e9;
	/* apex: the highest extremum; crossings: left-most and right-most */
	{
		qaws_scalar best = -1e30f;
		for (i = 0; i < ne; i++)
		{
			qaws_eval_result_2d r;
			qaws_curve_evaluate_2d(arch, ext[i], QAWS_EVAL_FLAG_POSITION, &r);
			if (r.position.y > best)
			{
				best = r.position.y;
				f->t_apex = ext[i];
			}
		}
		f->t_cross[0] = f->t_cross[1] = hits[0].parameter_a;
		f->tb_cross[0] = f->tb_cross[1] = hits[0].parameter_b;
		for (i = 1; i < nh; i++)
		{
			if (hits[i].parameter_a < f->t_cross[0]) { f->t_cross[0] = hits[i].parameter_a; f->tb_cross[0] = hits[i].parameter_b; }
			if (hits[i].parameter_a > f->t_cross[1]) { f->t_cross[1] = hits[i].parameter_a; f->tb_cross[1] = hits[i].parameter_b; }
		}
	}
	{
		qaws_curve_extremum e, eb;
		qaws_curve_extremum_tangent(NULL, arch, up, f->t_apex, NULL, NULL, &e, NULL);
		f->apex = e.position;
		loss += (e.position.x - AR_APEX_X) * (e.position.x - AR_APEX_X) + (e.value - AR_APEX_Y) * (e.value - AR_APEX_Y);
		if (grad)
		{
			memset(&eb, 0, sizeof(eb));
			eb.position.x = (qaws_scalar)(2 * (e.position.x - AR_APEX_X));
			eb.value = (qaws_scalar)(2 * (e.value - AR_APEX_Y));
			qaws_curve_extremum_adjoint(NULL, arch, up, f->t_apex, &eb, &gv, NULL);
		}
	}
	for (i = 0; i < 2; i++)
	{
		qaws_curve_pair_point p, pb;
		double target = i == 0 ? AR_X0 : AR_X1;
		qaws_curve_pair_point_tangent(NULL, arch, deck, f->t_cross[i], f->tb_cross[i], NULL, NULL, &p, NULL);
		f->cross[i] = p.position_a;
		loss += (p.position_a.x - target) * (p.position_a.x - target);
		if (grad)
		{
			memset(&pb, 0, sizeof(pb));
			pb.position_a.x = (qaws_scalar)(2 * (p.position_a.x - target));
			qaws_curve_pair_point_adjoint(NULL, arch, deck, f->t_cross[i], f->tb_cross[i], &pb, &gv, NULL);
		}
	}
	/* a little fairness */
	{
		qaws_scalar bg[AR_N * 2], bend = 0;
		qaws_field_view bf;
		qaws_diff_views bv = one_field(&bf, QAWS_FIELD_CONTROL_POINTS, bg, AR_N, 2);
		memset(bg, 0, sizeof(bg));
		qaws_curve_functional_gradient(NULL, arch, QAWS_FUNCTIONAL_BENDING, 0, &bv, &bend);
		loss += 1e-3 * bend;
		if (grad)
			for (i = 0; i < AR_N * 2; i++)
				grad[i] += (qaws_scalar)(1e-3 * bg[i]);
	}
	f->ok = 1;
	return loss;
}

static void demo_arch(void)
{
	qaws_scalar cps[AR_N * 2], g[AR_N * 2], init[AR_N * 2];
	qaws_scalar ghosts[12][AR_N * 2];
	arch_features f0, f1;
	double loss[AR_ITERS + 1];
	qaws_curve* deck = arch_deck();
	adam opt;
	int it, i, ng = 0;
	svg s;
	char buf[200];
	static double const start[AR_N][2] = {
		{ 0.0, 0.0 }, { 0.3, 0.9 }, { 0.9, 1.9 }, { 2.4, 1.6 }, { 3.6, 3.0 }, { 4.8, 2.2 }, { 5.6, 1.0 }, { 6.0, 0.0 } };

	for (i = 0; i < AR_N; i++)
	{
		cps[2 * i] = (qaws_scalar)start[i][0];
		cps[2 * i + 1] = (qaws_scalar)start[i][1];
	}
	memcpy(init, cps, sizeof(cps));
	memset(&opt, 0, sizeof(opt));
	for (it = 0; it <= AR_ITERS; it++)
	{
		qaws_curve* c = arch_curve(cps);
		arch_features f;
		loss[it] = arch_loss(c, deck, g, &f);
		qaws_curve_destroy(c);
		if (it % 16 == 0 && ng < 12)
			memcpy(ghosts[ng++], cps, sizeof(cps));
		if (it == AR_ITERS)
			break;
		/* feet stay on the ground */
		g[0] = g[1] = g[2 * AR_N - 2] = g[2 * AR_N - 1] = 0;
		adam_step(&opt, cps, g, AR_N * 2, 0.04 * (1.0 - 0.8 * it / (double)AR_ITERS));
	}

	svg_open(&s, "showcase/12_arch_features.svg", 1200, 600, "Arch design from implicit features",
		"Deck crossings (curve pair solves) and the apex (extremum solve) are steered to targets by adjoints of the implicit "
		"function theorem; the finders only pick which solutions exist.");
	{
		viewport v = { 24, 80, 760, 500, -0.6, 6.6, -0.5, 3.6 };
		viewport lv = { 800, 80, 376, 300, 0, 0, 0, 0 };
		double xy[2 * 300];
		qaws_curve* c;

		svg_panel(&s, &v, "initial (grey), iterations (light blue), final (blue)");
		/* ground and deck */
		svg_line(&s, vx(&v, -0.6), vy(&v, 0), vx(&v, 6.6), vy(&v, 0), "#8c959f", 2, 1);
		svg_line(&s, vx(&v, -0.5), vy(&v, AR_DECK), vx(&v, 6.5), vy(&v, AR_DECK), "#57606a", 5, 0.85);
		svg_text(&s, vx(&v, 6.45), vy(&v, AR_DECK) - 8, 12, "#57606a", "end", "deck");
		for (i = 1; i < ng - 1; i++)
		{
			c = arch_curve(ghosts[i]);
			curve_polyline(c, &v, xy, 300);
			svg_polyline(&s, xy, 300, "#54aeff", 1.4, 0.45, 0);
			qaws_curve_destroy(c);
		}
		c = arch_curve(init);
		arch_loss(c, deck, NULL, &f0);
		curve_polyline(c, &v, xy, 300);
		svg_polyline(&s, xy, 300, "#8c959f", 2, 1, 1);
		qaws_curve_destroy(c);
		c = arch_curve(cps);
		arch_loss(c, deck, NULL, &f1);
		curve_polyline(c, &v, xy, 300);
		svg_polyline(&s, xy, 300, "#0969da", 3.2, 1, 0);
		control_polygon(&s, &v, cps, AR_N, "#0969da");
		qaws_curve_destroy(c);

		/* targets */
		svg_line(&s, vx(&v, AR_X0), vy(&v, AR_DECK) - 18, vx(&v, AR_X0), vy(&v, AR_DECK) + 18, "#cf222e", 2, 0.9);
		svg_line(&s, vx(&v, AR_X1), vy(&v, AR_DECK) - 18, vx(&v, AR_X1), vy(&v, AR_DECK) + 18, "#cf222e", 2, 0.9);
		svg_line(&s, vx(&v, AR_APEX_X) - 12, vy(&v, AR_APEX_Y), vx(&v, AR_APEX_X) + 12, vy(&v, AR_APEX_Y), "#cf222e", 2, 0.9);
		svg_line(&s, vx(&v, AR_APEX_X), vy(&v, AR_APEX_Y) - 12, vx(&v, AR_APEX_X), vy(&v, AR_APEX_Y) + 12, "#cf222e", 2, 0.9);
		svg_text(&s, vx(&v, AR_APEX_X) + 16, vy(&v, AR_APEX_Y) - 8, 12, "#cf222e", "start", "apex target");
		svg_text(&s, vx(&v, AR_X0), vy(&v, AR_DECK) + 34, 12, "#cf222e", "middle", "crossing target");
		svg_text(&s, vx(&v, AR_X1), vy(&v, AR_DECK) + 34, 12, "#cf222e", "middle", "crossing target");
		/* features: initial (grey) and final (blue) */
		svg_circle(&s, vx(&v, f0.apex.x), vy(&v, f0.apex.y), 5, "#ffffff", "#57606a");
		svg_circle(&s, vx(&v, f0.cross[0].x), vy(&v, f0.cross[0].y), 5, "#ffffff", "#57606a");
		svg_circle(&s, vx(&v, f0.cross[1].x), vy(&v, f0.cross[1].y), 5, "#ffffff", "#57606a");
		svg_circle(&s, vx(&v, f1.apex.x), vy(&v, f1.apex.y), 6, "#0969da", "#ffffff");
		svg_circle(&s, vx(&v, f1.cross[0].x), vy(&v, f1.cross[0].y), 6, "#0969da", "#ffffff");
		svg_circle(&s, vx(&v, f1.cross[1].x), vy(&v, f1.cross[1].y), 6, "#0969da", "#ffffff");

		svg_loss_plot(&s, &lv, loss, AR_ITERS + 1, "#0969da", "design loss (log scale)");
		sprintf(buf, "apex (%.3f, %.3f)   target (%.1f, %.1f)", f1.apex.x, f1.apex.y, AR_APEX_X, AR_APEX_Y);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 26, 13, "#24292f", "start", buf);
		sprintf(buf, "crossings x = %.3f, %.3f   targets %.1f, %.1f", f1.cross[0].x, f1.cross[1].x, AR_X0, AR_X1);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 48, 13, "#24292f", "start", buf);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 76, 12, "#57606a", "start", "qaws_curve_pair_point_adjoint (deck crossings)");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 94, 12, "#57606a", "start", "qaws_curve_extremum_adjoint (apex height, position)");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 112, 12, "#57606a", "start", "+ 1e-3 bending (functional gradient), feet masked");
	}
	svg_close(&s);
	qaws_curve_destroy(deck);
	printf("12_arch_features: loss %.4e -> %.4e, apex (%.4f, %.4f), crossings %.4f %.4f\n", loss[0], loss[AR_ITERS],
		f1.apex.x, f1.apex.y, f1.cross[0].x, f1.cross[1].x);
}


/* ================================================================== */
/*  13. Knot placement: fitting a surface with a crease by moving its  */
/*      knots together with its control points                        */
/* ================================================================== */

#define KP_NU 8
#define KP_NV 6
#define KP_UK (KP_NU + 4)
#define KP_VK (KP_NV + 4)
#define KP_SU 48
#define KP_SV 36
#define KP_ITERS 500

static double kp_target(double u, double v)
{
	/* a crease along u = 0.68 and a gentle wave across v */
	return 0.35 * tanh(18.0 * (u - 0.68)) + 0.12 * sin(2.0 * PI * v) * (1.0 - u);
}

static qaws_surface* kp_surface(qaws_scalar const* cps, qaws_scalar const* uk, qaws_scalar const* vk)
{
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = KP_NU;
	d.v_point_count = KP_NV;
	d.u_knots = uk;
	d.u_knot_count = KP_UK;
	d.v_knots = vk;
	d.v_knot_count = KP_VK;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

static void kp_samples(qaws_scalar* us, qaws_scalar* vs)
{
	int i, j;
	for (i = 0; i < KP_SU; i++)
		for (j = 0; j < KP_SV; j++)
		{
			us[i * KP_SV + j] = (qaws_scalar)((i + 0.5) / KP_SU);
			vs[i * KP_SV + j] = (qaws_scalar)((j + 0.5) / KP_SV);
		}
}

/* Mean squared height error and its gradient (z of control points, and
   interior knots when move_knots). */
static double kp_loss(qaws_scalar const* cps, qaws_scalar const* uk, qaws_scalar const* vk, int move_knots,
	qaws_scalar* g_cp, qaws_scalar* g_uk, qaws_scalar* g_vk)
{
	static qaws_scalar us[KP_SU * KP_SV], vs[KP_SU * KP_SV];
	static qaws_surface_jet jets[KP_SU * KP_SV], bars[KP_SU * KP_SV];
	static int init = 0;
	unsigned char uact[KP_UK], vact[KP_VK];
	qaws_surface* s = kp_surface(cps, uk, vk);
	qaws_field_view fv[3];
	qaws_diff_views views;
	double loss = 0;
	int k, n = KP_SU * KP_SV;

	if (!init)
	{
		kp_samples(us, vs);
		init = 1;
	}
	for (k = 0; k < n; k++)
	{
		qaws_surface_jet j;
		double e;
		qaws_surface_eval_jet(s, us[k], vs[k], QAWS_SJET_P, &j);
		jets[k] = j;
		e = j.d[0].z - kp_target(us[k], vs[k]);
		loss += e * e;
		memset(&bars[k], 0, sizeof(bars[k]));
		bars[k].d[0].z = (qaws_scalar)(2 * e / n);
		bars[k].channels = QAWS_SJET_P;
	}
	if (g_cp)
	{
		memset(g_cp, 0, sizeof(qaws_scalar) * KP_NU * KP_NV * 3);
		memset(g_uk, 0, sizeof(qaws_scalar) * KP_UK);
		memset(g_vk, 0, sizeof(qaws_scalar) * KP_VK);
		for (k = 0; k < KP_UK; k++) uact[k] = (unsigned char)(k > 3 && k < KP_UK - 4);
		for (k = 0; k < KP_VK; k++) vact[k] = (unsigned char)(k > 3 && k < KP_VK - 4);
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, g_cp, KP_NU * KP_NV, 3);
		fv[0].component_mask = 1u << 2;
		fv[1] = qaws_field_view_make(QAWS_FIELD_U_KNOTS, g_uk, KP_UK, 1);
		fv[1].active = uact;
		fv[2] = qaws_field_view_make(QAWS_FIELD_V_KNOTS, g_vk, KP_VK, 1);
		fv[2].active = vact;
		views.fields = fv;
		views.field_count = move_knots ? 3u : 1u;
		views.children = NULL;
		views.child_count = 0;
		qaws_surface_eval_batch_adjoint(NULL, s, us, vs, (unsigned int)n, QAWS_SJET_P, bars, &views, NULL, NULL);
	}
	qaws_surface_destroy(s);
	return loss / n;
}

/* Interior knots stay ordered with a minimum gap. */
static void kp_project(qaws_scalar* k, int count)
{
	int i, lo = 4, hi = count - 5;
	for (i = lo; i <= hi; i++)
	{
		double min = (i == lo ? 0.0 : k[i - 1]) + 0.025;
		if (k[i] < min) k[i] = (qaws_scalar)min;
	}
	for (i = hi; i >= lo; i--)
	{
		double max = (i == hi ? 1.0 : k[i + 1]) - 0.025;
		if (k[i] > max) k[i] = (qaws_scalar)max;
	}
}

static void kp_fit(qaws_scalar* cps, qaws_scalar* uk, qaws_scalar* vk, int move_knots, double* loss)
{
	static qaws_scalar g_cp[KP_NU * KP_NV * 3], g_uk[KP_UK], g_vk[KP_VK];
	qaws_scalar z[KP_NU * KP_NV], gz[KP_NU * KP_NV];
	adam a_cp, a_k;
	int it, i;
	memset(&a_cp, 0, sizeof(a_cp));
	memset(&a_k, 0, sizeof(a_k));
	for (it = 0; it <= KP_ITERS; it++)
	{
		loss[it] = kp_loss(cps, uk, vk, move_knots, g_cp, g_uk, g_vk);
		if (it == KP_ITERS)
			break;
		for (i = 0; i < KP_NU * KP_NV; i++)
		{
			z[i] = cps[i * 3 + 2];
			gz[i] = g_cp[i * 3 + 2];
		}
		adam_step(&a_cp, z, gz, KP_NU * KP_NV, 0.02 * (1.0 - 0.7 * it / (double)KP_ITERS));
		for (i = 0; i < KP_NU * KP_NV; i++)
			cps[i * 3 + 2] = z[i];
		if (move_knots)
		{
			qaws_scalar kk[KP_UK + KP_VK], gk[KP_UK + KP_VK];
			memcpy(kk, uk, sizeof(qaws_scalar) * KP_UK);
			memcpy(kk + KP_UK, vk, sizeof(qaws_scalar) * KP_VK);
			memcpy(gk, g_uk, sizeof(qaws_scalar) * KP_UK);
			memcpy(gk + KP_UK, g_vk, sizeof(qaws_scalar) * KP_VK);
			adam_step(&a_k, kk, gk, KP_UK + KP_VK, 0.004 * (1.0 - 0.7 * it / (double)KP_ITERS));
			memcpy(uk, kk, sizeof(qaws_scalar) * KP_UK);
			memcpy(vk, kk + KP_UK, sizeof(qaws_scalar) * KP_VK);
			kp_project(uk, KP_UK);
			kp_project(vk, KP_VK);
		}
	}
}

static void kp_error_map(svg* s, viewport const* v, qaws_scalar const* cps, qaws_scalar const* uk, qaws_scalar const* vk,
	double vmax, char* buf)
{
	qaws_surface* sf = kp_surface(cps, uk, vk);
	int i, j;
	double cw = v->w / 60.0, ch = (v->h - 30) / 45.0;
	for (i = 0; i < 60; i++)
		for (j = 0; j < 45; j++)
		{
			qaws_surface_jet jt;
			double u = (i + 0.5) / 60.0, w = (j + 0.5) / 45.0, e;
			char col[32];
			qaws_surface_eval_jet(sf, (qaws_scalar)u, (qaws_scalar)w, QAWS_SJET_P, &jt);
			e = fabs(jt.d[0].z - kp_target(u, w));
			heat(e / vmax, col);
			fprintf(s->f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\"/>\n",
				v->x0 + i * cw, v->y0 + 30 + (44 - j) * ch, cw + 0.4, ch + 0.4, col);
		}
	for (i = 4; i < KP_UK - 4; i++)
		svg_line(s, v->x0 + uk[i] * v->w, v->y0 + 30, v->x0 + uk[i] * v->w, v->y0 + v->h, "#ffffff", 2, 0.95);
	for (i = 4; i < KP_VK - 4; i++)
		svg_line(s, v->x0, v->y0 + v->h - vk[i] * (v->h - 30), v->x0 + v->w, v->y0 + v->h - vk[i] * (v->h - 30), "#ffffff", 2, 0.95);
	svg_line(s, v->x0 + 0.68 * v->w, v->y0 + 30, v->x0 + 0.68 * v->w, v->y0 + v->h, "#1b1f24", 1.2, 0.6);
	qaws_surface_destroy(sf);
	(void)buf;
}

static void demo_knot_placement(void)
{
	static qaws_scalar cps0[KP_NU * KP_NV * 3], cps1[KP_NU * KP_NV * 3];
	qaws_scalar uk0[KP_UK], vk0[KP_VK], uk1[KP_UK], vk1[KP_VK];
	static double l0[KP_ITERS + 1], l1[KP_ITERS + 1];
	int i, j;
	svg s;
	char buf[256];

	for (i = 0; i < KP_NU; i++)
		for (j = 0; j < KP_NV; j++)
		{
			qaws_scalar* p = &cps0[(i * KP_NV + j) * 3];
			p[0] = (qaws_scalar)(i / (double)(KP_NU - 1));
			p[1] = (qaws_scalar)(j / (double)(KP_NV - 1));
			p[2] = 0;
		}
	for (i = 0; i < KP_UK; i++)
		uk0[i] = (qaws_scalar)(i < 4 ? 0.0 : (i >= KP_UK - 4 ? 1.0 : (i - 3) / (double)(KP_NU - 3)));
	for (i = 0; i < KP_VK; i++)
		vk0[i] = (qaws_scalar)(i < 4 ? 0.0 : (i >= KP_VK - 4 ? 1.0 : (i - 3) / (double)(KP_NV - 3)));
	memcpy(cps1, cps0, sizeof(cps0));
	memcpy(uk1, uk0, sizeof(uk0));
	memcpy(vk1, vk0, sizeof(vk0));

	kp_fit(cps0, uk0, vk0, 0, l0);
	kp_fit(cps1, uk1, vk1, 1, l1);

	svg_open(&s, "showcase/13_knot_placement.svg", 1200, 600, "Knot placement: fitting a crease by moving knots",
		"Bicubic B-spline (8 x 6) fit to a height field with a crease at u = 0.68 (dark line). Interior U and V knots "
		"(white lines) move with the control points through knot adjoints.");
	{
		viewport a = { 24, 80, 370, 470, 0, 1, 0, 1 }, b = { 410, 80, 370, 470, 0, 1, 0, 1 };
		viewport lv = { 800, 80, 376, 300, 0, 0, 0, 0 };
		double vmax = 0.06, xy0[2 * (KP_ITERS + 1)], xy1[2 * (KP_ITERS + 1)], lo = 1e300, hi = -1e300;
		int it;
		svg_panel(&s, &a, "control points only, uniform knots: |error|");
		kp_error_map(&s, &a, cps0, uk0, vk0, vmax, buf);
		svg_panel(&s, &b, "control points + knots: |error|");
		kp_error_map(&s, &b, cps1, uk1, vk1, vmax, buf);
		svg_colorbar(&s, b.x0 + 12, b.y0 + b.h + 8, 200, 8, "0", "0.06");

		for (it = 0; it <= KP_ITERS; it++)
		{
			double a0 = log10(sqrt(l0[it])), a1 = log10(sqrt(l1[it]));
			if (a0 < lo) lo = a0;
			if (a1 < lo) lo = a1;
			if (a0 > hi) hi = a0;
			if (a1 > hi) hi = a1;
		}
		svg_panel(&s, &lv, "log10 RMS error per iteration");
		for (it = 0; it <= KP_ITERS; it++)
		{
			xy0[2 * it] = xy1[2 * it] = lv.x0 + 14 + (lv.w - 28) * it / (double)KP_ITERS;
			xy0[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(sqrt(l0[it]))) / (hi - lo);
			xy1[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(sqrt(l1[it]))) / (hi - lo);
		}
		svg_polyline(&s, xy0, KP_ITERS + 1, "#8c959f", 2, 1, 0);
		svg_polyline(&s, xy1, KP_ITERS + 1, "#0969da", 2.6, 1, 0);
		svg_text(&s, lv.x0 + lv.w - 14, xy0[2 * KP_ITERS + 1] - 8, 11, "#57606a", "end", "control points only");
		svg_text(&s, lv.x0 + lv.w - 14, xy1[2 * KP_ITERS + 1] + 16, 11, "#0969da", "end", "control points + knots");

		sprintf(buf, "RMS: uniform knots %.4f   moved knots %.4f", sqrt(l0[KP_ITERS]), sqrt(l1[KP_ITERS]));
		svg_text(&s, lv.x0, lv.y0 + lv.h + 26, 13, "#24292f", "start", buf);
		sprintf(buf, "U knots: %.3f %.3f %.3f %.3f", uk1[4], uk1[5], uk1[6], uk1[7]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 48, 13, "#0969da", "start", buf);
		sprintf(buf, "V knots: %.3f %.3f", vk1[4], vk1[5]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 70, 13, "#0969da", "start", buf);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 98, 12, "#57606a", "start", "qaws_surface_eval_batch_adjoint with views:");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 116, 12, "#57606a", "start", "control points (z mask), U_KNOTS, V_KNOTS");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 134, 12, "#57606a", "start", "(interior knots active, ordered by projection)");
	}
	svg_close(&s);
	printf("13_knot_placement: RMS uniform %.5f, moved knots %.5f, u knots %.3f %.3f %.3f %.3f\n",
		sqrt(l0[KP_ITERS]), sqrt(l1[KP_ITERS]), uk1[4], uk1[5], uk1[6], uk1[7]);
}

/* ================================================================== */
/*  14. Constant-speed sampling: first and second order derivatives  */
/* ================================================================== */

#define AL_CP 10
#define AL_N 32
#define AL_P (2 * AL_CP)

/* Target: AL_N points at equal arc length along a figure-like curve. */
static void al_targets(double* q)
{
	enum { M = 4000 };
	static double px[M + 1], py[M + 1], cum[M + 1];
	int i, j = 0;
	for (i = 0; i <= M; i++)
	{
		double u = (double)i / M;
		px[i] = 0.4 + 6.4 * u;
		py[i] = 0.8 + 0.55 * sin(2.4 * PI * u) + 0.35 * sin(7.0 * PI * u * u);
		cum[i] = i ? cum[i - 1] + hypot(px[i] - px[i - 1], py[i] - py[i - 1]) : 0;
	}
	for (i = 0; i < AL_N; i++)
	{
		double s = cum[M] * i / (AL_N - 1), f;
		while (j + 1 < M && cum[j + 1] < s)
			j++;
		f = (s - cum[j]) / (cum[j + 1] - cum[j] + 1e-300);
		q[2 * i] = px[j] + f * (px[j + 1] - px[j]);
		q[2 * i + 1] = py[j] + f * (py[j + 1] - py[j]);
	}
}

static qaws_arc_length_target g_al_targets[AL_N];

/* E = sum |p_i - q_i|^2 over the constant-speed samples; gradient through
   the sample adjoint (p_bar = 2 (p - q)) when g is given. */
static double al_energy(qaws_scalar const* cps, double const* q, qaws_scalar* g, qaws_arc_length_sample* adj_out)
{
	qaws_curve* c = bspline_2d(cps, AL_CP);
	qaws_arc_length_sample s[AL_N], adj[AL_N];
	double e = 0;
	int i;
	qaws_curve_arc_length_sample_tangent(NULL, c, g_al_targets, NULL, AL_N, 0, NULL, s, NULL, NULL, NULL);
	for (i = 0; i < AL_N; i++)
	{
		double dx = s[i].position.x - q[2 * i], dy = s[i].position.y - q[2 * i + 1];
		e += dx * dx + dy * dy;
		adj[i].t = 0;
		adj[i].position = v3((qaws_scalar)(2 * dx), (qaws_scalar)(2 * dy), 0);
	}
	if (g)
	{
		qaws_field_view fv;
		qaws_diff_views views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, g, AL_CP, 2);
		memset(g, 0, sizeof(qaws_scalar) * AL_P);
		qaws_curve_arc_length_sample_adjoint(NULL, c, g_al_targets, AL_N, 0, adj, &views, NULL);
	}
	if (adj_out)
		memcpy(adj_out, adj, sizeof(adj));
	qaws_curve_destroy(c);
	return e;
}

/* Exact Hessian-vector product of E: J^T (2 J v) (first order forward then
   backward) plus the second order term of the samples (sample HVP with
   p_bar = 2 (p - q)). */
static void al_hvp(qaws_scalar const* cps, qaws_arc_length_sample const* adj, qaws_scalar const* v, qaws_scalar* out)
{
	qaws_curve* c = bspline_2d(cps, AL_CP);
	qaws_arc_length_sample jv[AL_N], w[AL_N];
	qaws_scalar vv[AL_P];
	qaws_field_view fv, fo;
	qaws_diff_views vin, vout;
	int i;
	memcpy(vv, v, sizeof(vv));
	vin = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, vv, AL_CP, 2);
	vout = one_field(&fo, QAWS_FIELD_CONTROL_POINTS, out, AL_CP, 2);
	memset(out, 0, sizeof(qaws_scalar) * AL_P);
	qaws_curve_arc_length_sample_tangent(NULL, c, g_al_targets, NULL, AL_N, 0, &vin, NULL, jv, NULL, NULL);
	for (i = 0; i < AL_N; i++)
	{
		w[i].t = 0;
		w[i].position = v3(2 * jv[i].position.x, 2 * jv[i].position.y, 0);
	}
	qaws_curve_arc_length_sample_adjoint(NULL, c, g_al_targets, AL_N, 0, w, &vout, NULL);
	qaws_curve_arc_length_sample_hvp(NULL, c, g_al_targets, AL_N, 0, adj, &vin, &vout);
	qaws_curve_destroy(c);
}

static double al_dot(qaws_scalar const* a, qaws_scalar const* b)
{
	double s = 0;
	int i;
	for (i = 0; i < AL_P; i++)
		s += (double)a[i] * b[i];
	return s;
}

/* Backtracking line search along d from x (energy e0, slope gd < 0). */
static double al_line(qaws_scalar* x, qaws_scalar const* d, double e0, double gd, double const* q, double step)
{
	qaws_scalar trial[AL_P];
	int i, k;
	for (k = 0; k < 40; k++, step *= 0.5)
	{
		double e;
		for (i = 0; i < AL_P; i++)
			trial[i] = (qaws_scalar)(x[i] + step * d[i]);
		e = al_energy(trial, q, NULL, NULL);
		if (e <= e0 + 1e-4 * step * gd)
		{
			memcpy(x, trial, sizeof(trial));
			return e;
		}
	}
	return e0;
}

static void demo_arc_length(void)
{
	static qaws_scalar const wiggle[AL_P] = { 0.0, 0.3, 0.6, 1.5, 0.9, 0.2, 1.3, 1.6, 2.6, 1.7, 3.4, 0.1, 4.6, 0.4, 5.2, 1.7, 6.4, 1.2, 7.2, 0.4 };
	viewport a = { 30, 80, 560, 500, -0.4, 7.6, -0.8, 2.6 };
	viewport b = { 610, 80, 560, 330, -0.2, 7.4, -0.4, 2.2 };
	viewport lp = { 610, 420, 560, 160, 0, 1, 0, 1 };
	double q[2 * AL_N], xy[2 * 400], loss_gd[61], loss_nt[61];
	qaws_scalar gd[AL_P], nt[AL_P];
	char buf[200];
	int i, it, n_gd = 60, n_nt = 12;
	svg s;

	for (i = 0; i < AL_N; i++)
	{
		g_al_targets[i].distance = 0;
		g_al_targets[i].fraction = (qaws_scalar)((double)i / (AL_N - 1));
	}
	svg_open(&s, "showcase/14_arc_length.svg", 1200, 620, "Constant-speed sampling: first and second order, forward and backward",
		"Samples at equal arc length, t solving L(t) = s: tangents and second tangents (forward), adjoints and Hessian-vector products (backward).");

	/* --- left: parameter-uniform vs constant-speed samples, and their tangents */
	svg_panel(&s, &a, "Constant speed vs uniform parameter; sample tangents for one moving control point");
	{
		qaws_curve* c = bspline_2d(wiggle, AL_CP);
		qaws_range r = qaws_curve_get_parameter_range(c);
		qaws_arc_length_sample val[AL_N], t1[AL_N], t2[AL_N];
		qaws_scalar dir[AL_P];
		qaws_field_view fv;
		qaws_diff_views views;
		unsigned int k_star = 4;
		double cxy[AL_P];
		curve_polyline(c, &a, xy, 400);
		svg_polyline(&s, xy, 400, "#57606a", 2.2, 1, 0);
		for (i = 0; i < AL_CP; i++)
		{
			cxy[2 * i] = vx(&a, wiggle[2 * i]);
			cxy[2 * i + 1] = vy(&a, wiggle[2 * i + 1]);
		}
		svg_polyline(&s, cxy, AL_CP, "#8c959f", 1, 0.6, 1);
		for (i = 0; i < AL_CP; i++)
			svg_circle(&s, cxy[2 * i], cxy[2 * i + 1], i == (int)k_star ? 7 : 3, i == (int)k_star ? "#cf222e" : "#ffffff", "#24292f");
		svg_line(&s, cxy[2 * k_star], cxy[2 * k_star + 1], cxy[2 * k_star] + 40, cxy[2 * k_star + 1], "#cf222e", 2.5, 1);
		/* parameter-uniform samples */
		for (i = 0; i < AL_N; i++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * i / (AL_N - 1)), QAWS_EVAL_FLAG_POSITION, &e);
			svg_circle(&s, vx(&a, e.position.x), vy(&a, e.position.y) + 26, 3.2, "#ffffff", "#8c959f");
		}
		/* constant-speed samples with first and second tangents for the control point moving right */
		memset(dir, 0, sizeof(dir));
		dir[2 * k_star] = 1;
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, dir, AL_CP, 2);
		qaws_curve_arc_length_sample_tangent(NULL, c, g_al_targets, NULL, AL_N, 0, &views, val, t1, t2, NULL);
		for (i = 0; i < AL_N; i++)
		{
			double x0 = vx(&a, val[i].position.x), y0 = vy(&a, val[i].position.y);
			double k1 = 0.9 * a.w / (a.xmax - a.xmin);
			svg_line(&s, x0, y0, x0 + k1 * t1[i].position.x, y0 - k1 * t1[i].position.y, "#0969da", 1.8, 0.9);
			svg_line(&s, x0, y0, x0 + 2.5 * k1 * t2[i].position.x, y0 - 2.5 * k1 * t2[i].position.y, "#8250df", 1.4, 0.9);
			svg_circle(&s, x0, y0, 3.6, "#cf222e", "#ffffff");
		}
		svg_text(&s, a.x0 + 12, a.y0 + a.h - 48, 12, "#57606a", "start",
			"hollow (shifted down): uniform parameter, bunched where the curve is slow;  red: constant speed");
		svg_text(&s, a.x0 + 12, a.y0 + a.h - 30, 12, "#0969da", "start",
			"blue: dp/dP for the red control point moving right (first order; samples slide along the curve)");
		svg_text(&s, a.x0 + 12, a.y0 + a.h - 12, 12, "#8250df", "start", "violet: d2p/dP2 (second order, x2.5)");
		qaws_curve_destroy(c);
	}

	/* --- right: fit the constant-speed samples to targets, GD vs Newton-CG */
	al_targets(q);
	for (i = 0; i < AL_CP; i++)
	{
		gd[2 * i] = (qaws_scalar)(0.4 + 6.4 * i / (AL_CP - 1));
		gd[2 * i + 1] = (qaws_scalar)0.8;
	}
	memcpy(nt, gd, sizeof(gd));
	{
		qaws_scalar g[AL_P];
		double e = al_energy(gd, q, g, NULL), step = 0.05;
		loss_gd[0] = e;
		for (it = 1; it <= n_gd; it++)
		{
			qaws_scalar d[AL_P];
			double gg = 0;
			for (i = 0; i < AL_P; i++)
			{
				d[i] = -g[i];
				gg += (double)g[i] * g[i];
			}
			e = al_line(gd, d, e, -gg, q, 4 * step);
			e = al_energy(gd, q, g, NULL);
			loss_gd[it] = e;
		}
	}
	{
		/* Newton-CG with the exact HVP (Steihaug truncation on negative curvature) */
		qaws_scalar g[AL_P];
		qaws_arc_length_sample adj[AL_N];
		double e = al_energy(nt, q, g, adj);
		loss_nt[0] = e;
		for (it = 1; it <= n_nt; it++)
		{
			qaws_scalar x[AL_P], r[AL_P], p[AL_P], hp[AL_P];
			double rr, gd_dot;
			int k;
			memset(x, 0, sizeof(x));
			for (i = 0; i < AL_P; i++)
			{
				r[i] = -g[i];
				p[i] = r[i];
			}
			rr = al_dot(r, r);
			for (k = 0; k < 40 && rr > 1e-24; k++)
			{
				double php, alpha, rr2, beta;
				al_hvp(nt, adj, p, hp);
				php = al_dot(p, hp);
				if (php <= 1e-14 * al_dot(p, p))
				{
					if (k == 0)
						memcpy(x, r, sizeof(x));
					break;
				}
				alpha = rr / php;
				for (i = 0; i < AL_P; i++)
				{
					x[i] += (qaws_scalar)(alpha * p[i]);
					r[i] -= (qaws_scalar)(alpha * hp[i]);
				}
				rr2 = al_dot(r, r);
				beta = rr2 / rr;
				rr = rr2;
				for (i = 0; i < AL_P; i++)
					p[i] = (qaws_scalar)(r[i] + beta * p[i]);
			}
			gd_dot = al_dot(g, x);
			if (gd_dot >= 0)
			{
				for (i = 0; i < AL_P; i++)
					x[i] = -g[i];
				gd_dot = -al_dot(g, g);
			}
			e = al_line(nt, x, e, gd_dot, q, 1.0);
			e = al_energy(nt, q, g, adj);
			loss_nt[it] = e;
		}
	}
	svg_panel(&s, &b, "Fit: constant-speed samples of a 10-point B-spline onto 32 targets");
	{
		qaws_curve* cg = bspline_2d(gd, AL_CP);
		qaws_curve* cn = bspline_2d(nt, AL_CP);
		qaws_arc_length_sample sn[AL_N];
		curve_polyline(cg, &b, xy, 400);
		svg_polyline(&s, xy, 400, "#d4a72c", 2, 0.9, 1);
		curve_polyline(cn, &b, xy, 400);
		svg_polyline(&s, xy, 400, "#1a7f37", 2.4, 1, 0);
		qaws_curve_arc_length_sample_tangent(NULL, cn, g_al_targets, NULL, AL_N, 0, NULL, sn, NULL, NULL, NULL);
		for (i = 0; i < AL_N; i++)
		{
			svg_line(&s, vx(&b, q[2 * i]), vy(&b, q[2 * i + 1]), vx(&b, sn[i].position.x), vy(&b, sn[i].position.y), "#cf222e", 1, 0.8);
			svg_circle(&s, vx(&b, q[2 * i]), vy(&b, q[2 * i + 1]), 3.2, "#ffffff", "#cf222e");
			svg_circle(&s, vx(&b, sn[i].position.x), vy(&b, sn[i].position.y), 2.4, "#1a7f37", "#1a7f37");
		}
		control_polygon(&s, &b, nt, AL_CP, "#1a7f37");
		qaws_curve_destroy(cg);
		qaws_curve_destroy(cn);
	}
	sprintf(buf, "green: Newton-CG, %d it., exact HVP (J^T J v + second order)   dashed: gradient descent, %d it.", n_nt, n_gd);
	svg_text(&s, b.x0 + 12, b.y0 + b.h - 10, 11, "#57606a", "start", buf);
	{
		/* both losses on one log axis */
		double lo = 1e300, hi = -1e300, pts[2 * 61];
		int k;
		svg_panel(&s, &lp, "Energy (log10) per iteration");
		for (k = 0; k <= n_gd; k++) { double l = log10(loss_gd[k] + 1e-300); if (l < lo) lo = l; if (l > hi) hi = l; }
		for (k = 0; k <= n_nt; k++) { double l = log10(loss_nt[k] + 1e-300); if (l < lo) lo = l; if (l > hi) hi = l; }
		for (k = 0; k <= n_gd; k++)
		{
			pts[2 * k] = lp.x0 + 16 + (lp.w - 32) * k / (double)n_gd;
			pts[2 * k + 1] = lp.y0 + 30 + (lp.h - 46) * (hi - log10(loss_gd[k] + 1e-300)) / (hi - lo);
		}
		svg_polyline(&s, pts, n_gd + 1, "#d4a72c", 2, 1, 0);
		for (k = 0; k <= n_nt; k++)
		{
			pts[2 * k] = lp.x0 + 16 + (lp.w - 32) * k / (double)n_gd;
			pts[2 * k + 1] = lp.y0 + 30 + (lp.h - 46) * (hi - log10(loss_nt[k] + 1e-300)) / (hi - lo);
		}
		svg_polyline(&s, pts, n_nt + 1, "#1a7f37", 2.4, 1, 0);
		sprintf(buf, "start %.3g   gradient descent %.3g   Newton-CG %.3g", loss_gd[0], loss_gd[n_gd], loss_nt[n_nt]);
		svg_text(&s, lp.x0 + lp.w - 12, lp.y0 + lp.h - 6, 11, "#57606a", "end", buf);
	}
	svg_close(&s);
	printf("14_arc_length: start %.4g, gradient descent (%d it) %.4g, Newton-CG (%d it) %.4g\n", loss_gd[0], n_gd, loss_gd[n_gd],
		n_nt, loss_nt[n_nt]);
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
	demo_networks();
	demo_projection();
	demo_soap_film();
	demo_parameter_correction();
	demo_arch();
	demo_knot_placement();
	demo_arc_length();
	return 0;
}
