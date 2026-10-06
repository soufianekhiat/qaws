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

