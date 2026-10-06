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

