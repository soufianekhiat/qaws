/* ================================================================== */
/*  2. Exact rational Bezier evaluation against f32 and f64           */
/* ================================================================== */

#define EB_DEG 16
#define EB_SAMPLES 400

/* Naive single-precision rational De Casteljau (homogeneous, projected). */
static void eb_f32(float const* h, int n, float t, float* out)
{
	float q[(EB_DEG + 1) * 3];
	int r, i, c;
	memcpy(q, h, sizeof(float) * (size_t)(n + 1) * 3);
	for (r = 1; r <= n; r++)
		for (i = 0; i + r <= n; i++)
			for (c = 0; c < 3; c++)
				q[i * 3 + c] = (1.0f - t) * q[i * 3 + c] + t * q[(i + 1) * 3 + c];
	out[0] = q[0] / q[2];
	out[1] = q[1] / q[2];
}

static double eb_rel(double const* a, double const* exact)
{
	double d = hypot(a[0] - exact[0], a[1] - exact[1]), m = hypot(exact[0], exact[1]);
	double r = d / (m > 0 ? m : 1);
	return r > 1e-18 ? r : 1e-18;
}

static void eb_plot(svg* s, viewport const* v, double const* err, int n, char const* color, double width)
{
	double xy[2 * EB_SAMPLES];
	int i;
	for (i = 0; i < n; i++)
	{
		xy[2 * i] = vx(v, (double)i / (n - 1));
		xy[2 * i + 1] = vy(v, log10(err[i]));
	}
	svg_polyline(s, xy, n, color, width, 0.9, 0);
}

static void eb_axes(svg* s, viewport const* v, char const* title)
{
	int e;
	char buf[32];
	svg_panel(s, v, title);
	for (e = (int)v->ymin; e <= (int)v->ymax; e += 2)
	{
		svg_line(s, v->x0, vy(v, e), v->x0 + v->w, vy(v, e), "#d0d7de", 0.6, 1);
		sprintf(buf, "1e%d", e);
		svg_text(s, v->x0 + 4, vy(v, e) - 2, 10, "#8c959f", "start", buf);
	}
}

static void demo_exact_bezier(void)
{
	qaws_scalar cps[(EB_DEG + 1) * 2], ws[EB_DEG + 1];
	float h32[(EB_DEG + 1) * 3];
	double e32[EB_SAMPLES], e0[EB_SAMPLES], e1[EB_SAMPLES], e2[EB_SAMPLES], e3[EB_SAMPLES], worst[5] = { 0, 0, 0, 0, 0 };
	qaws_rational_bezier_desc d;
	qaws_curve* c = NULL;
	qaws_exact_curve* ec = NULL;
	qaws_exact_desc desc;
	qaws_exact_report rep;
	viewport vc = { 20, 80, 380, 380, -9, 9, -9, 9 };
	viewport vp = { 420, 80, 400, 380, 0, 1, -18, 0 };
	viewport vd = { 840, 80, 400, 380, 0, 1, -18, 0 };
	char buf[300];
	int i;
	svg s;
	/* a looping degree-16 curve, weights from 1 to 10^4, every control
	   point on the 2^-20 lattice */
	for (i = 0; i <= EB_DEG; i++)
	{
		double a = 2 * 3.14159265358979 * i / EB_DEG * 1.35, r = 3 + 4 * fabs(sin(1.7 * i));
		cps[2 * i] = (qaws_scalar)ldexp(nearbyint(ldexp(r * cos(a), 20)), -20);
		cps[2 * i + 1] = (qaws_scalar)ldexp(nearbyint(ldexp(r * sin(a), 20)), -20);
		ws[i] = (qaws_scalar)nearbyint(pow(10.0, 4.0 * (0.5 + 0.5 * sin(2.3 * i))));
		h32[3 * i] = (float)(cps[2 * i] * ws[i]);
		h32[3 * i + 1] = (float)(cps[2 * i + 1] * ws[i]);
		h32[3 * i + 2] = (float)ws[i];
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = EB_DEG;
	d.control_points = cps;
	d.control_point_count = EB_DEG + 1;
	d.weights = ws;
	d.weight_count = EB_DEG + 1;
	qaws_curve_create_rational_bezier(&d, &c);
	qaws_exact_desc_default(&desc);
	qaws_exact_curve_prepare(&desc, c, &ec, &rep);
	for (i = 0; i < EB_SAMPLES; i++)
	{
		/* parameters on the 2^-24 lattice, so every evaluator sees the same t */
		double t = ldexp(nearbyint(ldexp((i + 0.5) / EB_SAMPLES, 24)), -24), ex[8], f[2];
		float g[2];
		qaws_eval_result_2d r;
		qaws_exact_curve_evaluate(ec, t, 3, ex, NULL);
		qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);
		eb_f32(h32, EB_DEG, (float)t, g);
		f[0] = g[0];
		f[1] = g[1];
		e32[i] = eb_rel(f, ex);
		f[0] = r.position.x; f[1] = r.position.y; e0[i] = eb_rel(f, ex);
		f[0] = r.d1.x; f[1] = r.d1.y; e1[i] = eb_rel(f, ex + 2);
		f[0] = r.d2.x; f[1] = r.d2.y; e2[i] = eb_rel(f, ex + 4);
		f[0] = r.d3.x; f[1] = r.d3.y; e3[i] = eb_rel(f, ex + 6);
		worst[0] = e32[i] > worst[0] ? e32[i] : worst[0];
		worst[1] = e0[i] > worst[1] ? e0[i] : worst[1];
		worst[2] = e1[i] > worst[2] ? e1[i] : worst[2];
		worst[3] = e2[i] > worst[3] ? e2[i] : worst[3];
		worst[4] = e3[i] > worst[4] ? e3[i] : worst[4];
	}
	svg_open(&s, "showcase/exact2_bezier.svg", 1260, 520, "Exact rational Bezier evaluation: the reference that f32 and f64 are measured against",
		"Degree 16, weights from 1 to 10^4 (circle size), control points on the 2^-20 lattice. The exact kernel evaluates in multi-limb integers and rounds once; errors are relative to it.");
	svg_panel(&s, &vc, "the curve and its control polygon");
	{
		double xy[2 * 401], pg[2 * (EB_DEG + 1)];
		for (i = 0; i <= 400; i++)
		{
			double ex[2];
			qaws_exact_curve_evaluate(ec, i / 400.0, 0, ex, NULL);
			xy[2 * i] = vx(&vc, ex[0]);
			xy[2 * i + 1] = vy(&vc, ex[1]);
		}
		for (i = 0; i <= EB_DEG; i++)
		{
			pg[2 * i] = vx(&vc, cps[2 * i]);
			pg[2 * i + 1] = vy(&vc, cps[2 * i + 1]);
		}
		svg_polyline(&s, pg, EB_DEG + 1, "#8c959f", 0.8, 0.8, 1);
		for (i = 0; i <= EB_DEG; i++)
			svg_circle(&s, pg[2 * i], pg[2 * i + 1], 1.5 + 1.6 * log10(ws[i]), "#ffffff", "#8c959f");
		svg_polyline(&s, xy, 401, "#24292f", 2.2, 1, 0);
	}
	eb_axes(&s, &vp, "position: relative error against exact");
	eb_plot(&s, &vp, e32, EB_SAMPLES, "#bc4c00", 1.4);
	eb_plot(&s, &vp, e0, EB_SAMPLES, "#0969da", 1.4);
	svg_text(&s, vp.x0 + vp.w - 10, vp.y0 + 46, 12, "#bc4c00", "end", "naive f32 De Casteljau");
	svg_text(&s, vp.x0 + vp.w - 10, vp.y0 + 62, 12, "#0969da", "end", "qaws f64 evaluator");
	eb_axes(&s, &vd, "f64 derivatives: relative error against exact");
	eb_plot(&s, &vd, e1, EB_SAMPLES, "#1a7f37", 1.3);
	eb_plot(&s, &vd, e2, EB_SAMPLES, "#8250df", 1.3);
	eb_plot(&s, &vd, e3, EB_SAMPLES, "#cf222e", 1.3);
	svg_text(&s, vd.x0 + vd.w - 10, vd.y0 + 46, 12, "#1a7f37", "end", "C'");
	svg_text(&s, vd.x0 + vd.w - 10, vd.y0 + 62, 12, "#8250df", "end", "C''");
	svg_text(&s, vd.x0 + vd.w - 10, vd.y0 + 78, 12, "#cf222e", "end", "C'''");
	sprintf(buf, "worst relative error: f32 position %.1e, f64 position %.1e, C' %.1e, C'' %.1e, C''' %.1e; exact path: correctly rounded, storage %u bits",
		worst[0], worst[1], worst[2], worst[3], worst[4], rep.storage_bits);
	svg_text(&s, 20, 495, 13, "#57606a", "start", buf);
	svg_close(&s);
	printf("exact2_bezier: worst f32 %.2e, f64 %.2e / %.2e / %.2e / %.2e\n", worst[0], worst[1], worst[2], worst[3], worst[4]);
	qaws_exact_curve_destroy(ec);
	qaws_curve_destroy(c);
}
