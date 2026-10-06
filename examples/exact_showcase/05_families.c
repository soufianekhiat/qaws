/* ================================================================== */
/*  5. Exact Catmull-Rom spans and polynomials far from the origin    */
/* ================================================================== */

#define EF_SAMPLES 1025

/* Naive single-precision Horner on the monomial coefficients. */
static float ef_horner_f32(double const* a, int n, float t)
{
	float r = (float)a[n];
	int k;
	for (k = n - 1; k >= 0; k--)
		r = r * t + (float)a[k];
	return r;
}

static void demo_families(void)
{
	static char const* const span_colors[7] = { "#0969da", "#1a7f37", "#bc4c00", "#8250df", "#cf222e", "#1b7c83", "#9a6700" };
	/* Chebyshev T6(u) = 32u^6 - 48u^4 + 18u^2 - 1, and u = t - 100 */
	static double const cheb[7] = { -1, 0, 18, 0, -48, 0, 32 };
	qaws_scalar crp[7 * 2], co[7 * 2];
	double a[7], e32[EF_SAMPLES], e64[EF_SAMPLES], worst32 = 0, worst64 = 0;
	qaws_catmull_rom_desc cd;
	qaws_polynomial_desc pd;
	qaws_curve* cr = NULL;
	qaws_curve* poly = NULL;
	qaws_exact_curve* ecr = NULL;
	qaws_exact_curve* epoly = NULL;
	qaws_exact_desc desc;
	qaws_exact_report rcr, rpoly;
	viewport vc = { 20, 80, 380, 380, -5.2, 5.2, -5.2, 5.2 };
	viewport vp = { 420, 80, 400, 380, -1.05, 1.05, -1.7, 1.7 };
	viewport ve = { 840, 80, 400, 380, -1, 1, -18, 9 };
	unsigned int s, i, j, k;
	char buf[300];
	svg sv;
	qaws_exact_desc_default(&desc);

	/* closed uniform Catmull-Rom through 7 lattice points */
	for (i = 0; i < 7; i++)
	{
		double ang = 2 * 3.14159265358979 * i / 7, r = (i % 2) ? 2.2 : 4.4;
		crp[2 * i] = (qaws_scalar)ldexp(nearbyint(ldexp(r * cos(ang), 20)), -20);
		crp[2 * i + 1] = (qaws_scalar)ldexp(nearbyint(ldexp(r * sin(ang), 20)), -20);
	}
	memset(&cd, 0, sizeof(cd));
	cd.dimension = QAWS_DIMENSION_2D;
	cd.control_points = crp;
	cd.control_point_count = 7;
	cd.parameterization = QAWS_PARAMETERIZATION_UNIFORM;
	cd.closed = 1;
	qaws_curve_create_catmull_rom(&cd, &cr);
	qaws_exact_curve_prepare(&desc, cr, &ecr, &rcr);

	/* the polynomial in t: a_k = sum_j cheb_j C(j, k) (-100)^(j - k), integers below 2^53 */
	for (k = 0; k <= 6; k++)
	{
		a[k] = 0;
		for (j = k; j <= 6; j++)
		{
			double b = 1, p = 1;
			for (i = 0; i < k; i++)
				b = b * (j - i) / (i + 1);
			for (i = 0; i < j - k; i++)
				p *= -100;
			a[k] += cheb[j] * b * p;
		}
		co[2 * k] = (qaws_scalar)(k == 0 ? -100 : (k == 1 ? 1 : 0));   /* x = t - 100 */
		co[2 * k + 1] = (qaws_scalar)a[k];
	}
	memset(&pd, 0, sizeof(pd));
	pd.dimension = QAWS_DIMENSION_2D;
	pd.degree = 6;
	pd.coefficients = co;
	pd.coefficient_count = 7;
	pd.t_min = 99;
	pd.t_max = 101;
	qaws_curve_create_polynomial(&pd, &poly);
	qaws_exact_curve_prepare(&desc, poly, &epoly, &rpoly);

	svg_open(&sv, "showcase/exact5_families.svg", 1260, 520, "Exact Hermite, Catmull-Rom and polynomial curves",
		"Left: a closed uniform Catmull-Rom curve converted exactly into one integer cubic Bezier per span (common factor 6). Middle, right: T6(t - 100) as monomials in t on [99, 101].");

	svg_panel(&sv, &vc, "closed uniform Catmull-Rom: exact Bezier spans");
	{
		double pg[2 * 8];
		for (i = 0; i <= 7; i++)
		{
			pg[2 * i] = vx(&vc, crp[2 * (i % 7)]);
			pg[2 * i + 1] = vy(&vc, crp[2 * (i % 7) + 1]);
		}
		svg_polyline(&sv, pg, 8, "#d0d7de", 1.0, 1, 1);
	}
	for (s = 0; s < qaws_exact_curve_span_count(ecr); s++)
	{
		double pts[2 * 4], t0, t1, xy[2 * 61], pg[2 * 4];
		unsigned int deg;
		qaws_exact_curve_span_bezier(ecr, s, &deg, &t0, &t1, pts, NULL);
		for (i = 0; i <= deg; i++)
		{
			pg[2 * i] = vx(&vc, pts[2 * i]);
			pg[2 * i + 1] = vy(&vc, pts[2 * i + 1]);
		}
		svg_polyline(&sv, pg, (int)deg + 1, span_colors[s % 7], 1.0, 0.8, 1);
		for (i = 1; i < deg; i++)
			svg_circle(&sv, pg[2 * i], pg[2 * i + 1], 2.5, "#ffffff", span_colors[s % 7]);
		for (i = 0; i <= 60; i++)
		{
			double out[2];
			qaws_exact_curve_evaluate(ecr, t0 + (t1 - t0) * i / 60.0, 0, out, NULL);
			xy[2 * i] = vx(&vc, out[0]);
			xy[2 * i + 1] = vy(&vc, out[1]);
		}
		svg_polyline(&sv, xy, 61, span_colors[s % 7], 3.0, 1, 0);
	}
	for (i = 0; i < 7; i++)
		svg_circle(&sv, vx(&vc, crp[2 * i]), vy(&vc, crp[2 * i + 1]), 4, "#24292f", "#ffffff");
	sprintf(buf, "7 spans, integers up to %u bits, nothing rounded", rcr.storage_bits);
	svg_text(&sv, vc.x0 + 12, vc.y0 + vc.h - 12, 12, "#57606a", "start", buf);

	/* the polynomial: exact curve, f32 Horner samples, errors */
	svg_panel(&sv, &vp, "T6(t - 100) from monomial coefficients");
	{
		double xy[2 * EF_SAMPLES];
		unsigned int nxy = 0;
		for (i = 0; i < EF_SAMPLES; i++)
		{
			double t = 99 + i / 512.0, ex[2];
			float g = ef_horner_f32(a, 6, (float)t);
			qaws_eval_result_2d r;
			qaws_exact_curve_evaluate(epoly, t, 0, ex, NULL);
			qaws_curve_evaluate_2d(poly, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION, &r);
			e32[i] = fabs(g - ex[1]);
			e64[i] = fabs(r.position.y - ex[1]);
			if (e32[i] > worst32) worst32 = e32[i];
			if (e64[i] > worst64) worst64 = e64[i];
			e32[i] = e32[i] > 1e-17 ? e32[i] : 1e-17;
			e64[i] = e64[i] > 1e-17 ? e64[i] : 1e-17;
			if (i % 4 == 0 && g > vp.ymin && g < vp.ymax)
				svg_circle(&sv, vx(&vp, ex[0]), vy(&vp, g), 1.6, "#bc4c00", "none");
			xy[2 * nxy] = vx(&vp, ex[0]);
			xy[2 * nxy + 1] = vy(&vp, ex[1]);
			nxy++;
		}
		svg_polyline(&sv, xy, (int)nxy, "#24292f", 2.2, 1, 0);
		svg_text(&sv, vp.x0 + vp.w - 10, vp.y0 + 46, 12, "#24292f", "end", "exact (rounded once)");
		svg_text(&sv, vp.x0 + vp.w - 10, vp.y0 + 62, 12, "#bc4c00", "end", "f32 Horner: most samples off the chart");
	}
	eb_axes(&sv, &ve, "absolute error against exact (log10)");
	{
		/* envelope: the largest error over each block of 8 samples */
		double xy[2 * EF_SAMPLES];
		unsigned int nb = (EF_SAMPLES - 1) / 8, b;
		for (k = 0; k < 2; k++)
		{
			double const* e = k == 0 ? e32 : e64;
			for (b = 0; b < nb; b++)
			{
				double m = 0;
				for (i = 8 * b; i <= 8 * b + 8; i++)
					m = e[i] > m ? e[i] : m;
				xy[2 * b] = vx(&ve, -1 + (8 * b + 4) / 512.0);
				xy[2 * b + 1] = vy(&ve, log10(m));
			}
			svg_polyline(&sv, xy, (int)nb, k == 0 ? "#bc4c00" : "#0969da", 1.6, 0.95, 0);
		}
		svg_line(&sv, ve.x0, vy(&ve, log10(ldexp(1, -53))), ve.x0 + ve.w, vy(&ve, log10(ldexp(1, -53))), "#1a7f37", 1.2, 1);
		svg_text(&sv, ve.x0 + ve.w - 10, vy(&ve, 5.6), 12, "#bc4c00", "end", "naive f32 Horner");
		svg_text(&sv, ve.x0 + ve.w - 10, vy(&ve, -3.4), 12, "#0969da", "end", "qaws f64 polynomial evaluator");
		svg_text(&sv, ve.x0 + ve.w - 10, vy(&ve, log10(ldexp(1, -53))) - 6, 12, "#1a7f37", "end", "exact path: below half an ulp");
	}
	sprintf(buf, "worst absolute error on values in [-1, 1]: f32 %.1e, f64 %.1e; exact path: Bernstein form from the dyadic coefficients, %u-bit integers, correctly rounded",
		worst32, worst64, rpoly.storage_bits);
	svg_text(&sv, 20, 495, 13, "#57606a", "start", buf);
	svg_close(&sv);
	printf("exact5_families: CR %u bits; polynomial %u bits, worst f32 %.2e, f64 %.2e\n", rcr.storage_bits, rpoly.storage_bits, worst32, worst64);
	qaws_exact_curve_destroy(ecr);
	qaws_exact_curve_destroy(epoly);
	qaws_curve_destroy(cr);
	qaws_curve_destroy(poly);
}
