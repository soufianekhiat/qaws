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


