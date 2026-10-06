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

