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

