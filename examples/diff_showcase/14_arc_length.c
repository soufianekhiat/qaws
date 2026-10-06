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

static qaws_cdf_target g_al_targets[AL_N];

/* E = sum |p_i - q_i|^2 over the constant-speed samples; gradient through
   the sample adjoint (p_bar = 2 (p - q)) when g is given. */
static double al_energy(qaws_scalar const* cps, double const* q, qaws_scalar* g, qaws_cdf_sample* adj_out)
{
	qaws_curve* c = bspline_2d(cps, AL_CP);
	qaws_cdf_sample s[AL_N], adj[AL_N];
	double e = 0;
	int i;
	qaws_curve_cdf_sample_tangent(NULL, c, NULL, g_al_targets, NULL, NULL, AL_N, 0, NULL, s, NULL, NULL, NULL);
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
		qaws_curve_cdf_sample_adjoint(NULL, c, NULL, g_al_targets, AL_N, 0, adj, &views, NULL);
	}
	if (adj_out)
		memcpy(adj_out, adj, sizeof(adj));
	qaws_curve_destroy(c);
	return e;
}

/* Exact Hessian-vector product of E: J^T (2 J v) (first order forward then
   backward) plus the second order term of the samples (sample HVP with
   p_bar = 2 (p - q)). */
static void al_hvp(qaws_scalar const* cps, qaws_cdf_sample const* adj, qaws_scalar const* v, qaws_scalar* out)
{
	qaws_curve* c = bspline_2d(cps, AL_CP);
	qaws_cdf_sample jv[AL_N], w[AL_N];
	qaws_scalar vv[AL_P];
	qaws_field_view fv, fo;
	qaws_diff_views vin, vout;
	int i;
	memcpy(vv, v, sizeof(vv));
	vin = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, vv, AL_CP, 2);
	vout = one_field(&fo, QAWS_FIELD_CONTROL_POINTS, out, AL_CP, 2);
	memset(out, 0, sizeof(qaws_scalar) * AL_P);
	qaws_curve_cdf_sample_tangent(NULL, c, NULL, g_al_targets, NULL, NULL, AL_N, 0, &vin, NULL, jv, NULL, NULL);
	for (i = 0; i < AL_N; i++)
	{
		w[i].t = 0;
		w[i].position = v3(2 * jv[i].position.x, 2 * jv[i].position.y, 0);
	}
	qaws_curve_cdf_sample_adjoint(NULL, c, NULL, g_al_targets, AL_N, 0, w, &vout, NULL);
	qaws_curve_cdf_sample_hvp(NULL, c, NULL, g_al_targets, AL_N, 0, adj, &vin, &vout);
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
		qaws_cdf_sample val[AL_N], t1[AL_N], t2[AL_N];
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
		qaws_curve_cdf_sample_tangent(NULL, c, NULL, g_al_targets, NULL, NULL, AL_N, 0, &views, val, t1, t2, NULL);
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
		qaws_cdf_sample adj[AL_N];
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
		qaws_cdf_sample sn[AL_N];
		curve_polyline(cg, &b, xy, 400);
		svg_polyline(&s, xy, 400, "#d4a72c", 2, 0.9, 1);
		curve_polyline(cn, &b, xy, 400);
		svg_polyline(&s, xy, 400, "#1a7f37", 2.4, 1, 0);
		qaws_curve_cdf_sample_tangent(NULL, cn, NULL, g_al_targets, NULL, NULL, AL_N, 0, NULL, sn, NULL, NULL, NULL);
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

