/* ================================================================== */
/*  18. Knots as sampling parameters: toward arc-length parameters    */
/* ================================================================== */

#define KS_N 31
#define KS_CP 10
#define KS_KNOTS 14

static qaws_curve* ks_curve(qaws_scalar const* cps, qaws_scalar const* knots)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = KS_CP;
	d.knots = knots;
	d.knot_count = KS_KNOTS;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

/* Loss sum (t_i - tau_i)^2 / N of the arc-length samples against uniform
   parameters, with its knot gradient (grad may be NULL). */
static double ks_loss(qaws_curve const* c, qaws_cdf_target const* tg, qaws_scalar* grad)
{
	qaws_cdf_sample val[KS_N], adj[KS_N];
	qaws_range r = qaws_curve_get_parameter_range(c);
	double loss = 0;
	int i;
	qaws_curve_cdf_sample_tangent(NULL, c, NULL, tg, NULL, NULL, KS_N, 0, NULL, val, NULL, NULL, NULL);
	memset(adj, 0, sizeof(adj));
	for (i = 0; i < KS_N; i++)
	{
		double tau = r.min_value + (r.max_value - r.min_value) * i / (KS_N - 1.0), e = val[i].t - tau;
		loss += e * e / KS_N;
		adj[i].t = (qaws_scalar)(2 * e / KS_N);
	}
	if (grad)
	{
		qaws_field_view fv;
		qaws_diff_views views = one_field(&fv, QAWS_FIELD_KNOTS, grad, KS_KNOTS, 1);
		memset(grad, 0, sizeof(qaws_scalar) * KS_KNOTS);
		qaws_curve_cdf_sample_adjoint(NULL, c, NULL, tg, KS_N, 0, adj, &views, NULL);
	}
	return loss;
}

/* Speed ratio max |C'| / min |C'| over the domain. */
static double ks_speed_ratio(qaws_curve const* c)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double lo = 1e300, hi = 0;
	int i;
	for (i = 0; i <= 400; i++)
	{
		qaws_eval_result_2d e;
		double s;
		qaws_curve_evaluate_2d(c, r.min_value + (r.max_value - r.min_value) * i / (qaws_scalar)400, QAWS_EVAL_FLAG_D1, &e);
		s = sqrt((double)e.d1.x * e.d1.x + (double)e.d1.y * e.d1.y);
		if (s < lo) lo = s;
		if (s > hi) hi = s;
	}
	return hi / lo;
}

/* Curve, its uniform-parameter points (blue), arc-length samples (red) and
   interior knots (green diamonds). */
static void ks_panel(svg* s, viewport const* v, qaws_curve const* c, qaws_scalar const* knots, qaws_cdf_target const* tg,
	char const* label, double loss)
{
	qaws_cdf_sample val[KS_N];
	qaws_range r = qaws_curve_get_parameter_range(c);
	double xy[2 * 400];
	char buf[128];
	int i;
	svg_panel(s, v, label);
	curve_polyline(c, v, xy, 400);
	svg_polyline(s, xy, 400, "#57606a", 2, 1, 0);
	qaws_curve_cdf_sample_tangent(NULL, c, NULL, tg, NULL, NULL, KS_N, 0, NULL, val, NULL, NULL, NULL);
	for (i = 0; i < KS_N; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, r.min_value + (r.max_value - r.min_value) * i / (qaws_scalar)(KS_N - 1), QAWS_EVAL_FLAG_POSITION, &e);
		svg_line(s, vx(v, e.position.x), vy(v, e.position.y), vx(v, val[i].position.x), vy(v, val[i].position.y), "#8250df", 1, 0.6);
		svg_circle(s, vx(v, e.position.x), vy(v, e.position.y), 2.6, "#0969da", "#ffffff");
		svg_circle(s, vx(v, val[i].position.x), vy(v, val[i].position.y), 2.6, "#cf222e", "#ffffff");
	}
	for (i = 4; i < KS_KNOTS - 4; i++)
	{
		qaws_eval_result_2d e;
		double x0, y0;
		qaws_curve_evaluate_2d(c, knots[i], QAWS_EVAL_FLAG_POSITION, &e);
		x0 = vx(v, e.position.x);
		y0 = vy(v, e.position.y);
		fprintf(s->f, "<path d=\"M%.1f %.1fL%.1f %.1fL%.1f %.1fL%.1f %.1fZ\" fill=\"#1a7f37\" stroke=\"#ffffff\"/>\n",
			x0, y0 - 6, x0 + 6, y0, x0, y0 + 6, x0 - 6, y0);
	}
	sprintf(buf, "loss %.3g, speed ratio %.2f", loss, ks_speed_ratio(c));
	svg_text(s, v->x0 + 12, v->y0 + v->h - 12, 12, "#57606a", "start", buf);
}

static void demo_knot_sampling(void)
{
	static qaws_scalar const cps[2 * KS_CP] = { 0.0, 0.0, 0.25, 0.7, 0.6, 1.1, 1.0, 0.9, 1.3, 0.4,
		2.8, 0.0, 4.4, 0.9, 5.5, 2.1, 6.8, 1.7, 7.3, 0.2 };
	static qaws_scalar const knots0[KS_KNOTS] = { 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7 };
	qaws_scalar knots[KS_KNOTS], grad[KS_KNOTS], kd[KS_KNOTS];
	qaws_cdf_target tg[KS_N];
	qaws_curve* c;
	double losses[301];
	adam opt;
	int i, it;
	svg s;
	viewport va = { 20, 80, 395, 440, -0.4, 7.7, -0.5, 2.6 };
	viewport vb = { 425, 80, 395, 440, -0.4, 7.7, -0.5, 2.6 };
	viewport vc = { 830, 80, 395, 440, -0.4, 7.7, -0.5, 2.6 };
	viewport vl = { 1235, 80, 245, 440, 0, 1, 0, 1 };
	for (i = 0; i < KS_N; i++)
	{
		tg[i].distance = 0;
		tg[i].fraction = (qaws_scalar)((double)i / (KS_N - 1));
	}
	memcpy(knots, knots0, sizeof(knots));
	svg_open(&s, "showcase/18_knot_sampling.svg", 1500, 560, "Knots as parameters of inverse-CDF sampling",
		"Red: 31 arc-length samples; blue: uniform parameters t = 7 i/30; green: interior knots. The knot adjoint of the samples drives the knots toward an arc-length parameterization.");

	/* a: tangents of the arc-length samples for knot 7 moving right */
	{
		qaws_cdf_sample val[KS_N], t1[KS_N];
		qaws_field_view fv;
		qaws_diff_views views;
		double xy[2 * 400], k1 = 0.6 * va.w / (va.xmax - va.xmin);
		qaws_eval_result_2d e;
		c = ks_curve(cps, knots);
		memset(kd, 0, sizeof(kd));
		kd[7] = 1;
		views = one_field(&fv, QAWS_FIELD_KNOTS, kd, KS_KNOTS, 1);
		svg_panel(&s, &va, "d samples / d knot u7 (first order tangents)");
		curve_polyline(c, &va, xy, 400);
		svg_polyline(&s, xy, 400, "#57606a", 2, 1, 0);
		qaws_curve_cdf_sample_tangent(NULL, c, NULL, tg, NULL, NULL, KS_N, 0, &views, val, t1, NULL, NULL);
		{
			/* longest arrow 50 px */
			double mx = 1e-12;
			for (i = 0; i < KS_N; i++)
			{
				double l = sqrt((double)t1[i].position.x * t1[i].position.x + (double)t1[i].position.y * t1[i].position.y);
				if (l > mx) mx = l;
			}
			k1 = 50 / mx;
		}
		for (i = 0; i < KS_N; i++)
		{
			double x0 = vx(&va, val[i].position.x), y0 = vy(&va, val[i].position.y);
			svg_line(&s, x0, y0, x0 + k1 * t1[i].position.x, y0 - k1 * t1[i].position.y, "#0969da", 1.6, 0.9);
			svg_circle(&s, x0, y0, 2.8, "#cf222e", "#ffffff");
		}
		qaws_curve_evaluate_2d(c, knots[7], QAWS_EVAL_FLAG_POSITION, &e);
		svg_circle(&s, vx(&va, e.position.x), vy(&va, e.position.y), 6, "#1a7f37", "#24292f");
		svg_text(&s, va.x0 + 12, va.y0 + va.h - 12, 12, "#57606a", "start", "arrows: sample tangents (longest 50 px)");
	}

	losses[0] = ks_loss(c, tg, NULL);
	ks_panel(&s, &vb, c, knots, tg, "uniform knots: parameters bunch where the curve is slow", losses[0]);
	qaws_curve_destroy(c);

	/* adjoint descent on the 6 interior knots, keeping them ordered */
	memset(&opt, 0, sizeof(opt));
	for (it = 0; it < 300; it++)
	{
		qaws_scalar x[6], g6[6];
		c = ks_curve(cps, knots);
		losses[it] = ks_loss(c, tg, grad);
		qaws_curve_destroy(c);
		for (i = 0; i < 6; i++)
		{
			x[i] = knots[4 + i];
			g6[i] = grad[4 + i];
		}
		adam_step(&opt, x, g6, 6, 0.02);
		for (i = 0; i < 6; i++)
		{
			double lo = (i == 0 ? 0.0 : (double)x[i - 1]) + 0.05, hi = 6.95 - 0.05 * (5 - i);
			knots[4 + i] = (qaws_scalar)(x[i] < lo ? lo : (x[i] > hi ? hi : x[i]));
		}
	}
	c = ks_curve(cps, knots);
	losses[300] = ks_loss(c, tg, NULL);
	ks_panel(&s, &vc, c, knots, tg, "300 adjoint steps on the interior knots (shape follows)", losses[300]);
	qaws_curve_destroy(c);
	svg_loss_plot(&s, &vl, losses, 301, "#8250df", "loss sum (t_i - tau_i)^2 / N");
	svg_close(&s);
	printf("18_knot_sampling: loss %.4g -> %.4g, knots %.3f %.3f %.3f %.3f %.3f %.3f\n", losses[0], losses[300],
		(double)knots[4], (double)knots[5], (double)knots[6], (double)knots[7], (double)knots[8], (double)knots[9]);
}
