/* ================================================================== */
/*  4. Exact NURBS: local Bezier extraction and the integer budget    */
/* ================================================================== */

static void demo_nurbs_extraction(void)
{
	static char const* const span_colors[6] = { "#0969da", "#1a7f37", "#bc4c00", "#8250df", "#cf222e", "#1b7c83" };
	qaws_scalar cps[9 * 2] = { 0, 0, 1.0f, 2.5f, 3.0f, 3.0f, 4.0f, 1.0f, 5.5f, 0.0f, 6.0f, 2.5f, 7.5f, 3.2f, 9.0f, 1.0f, 10.0f, 3.0f };
	qaws_scalar ws[9] = { 1, 4, 1, 2, 1, 6, 1, 3, 1 };
	/* degree 3, a triple interior knot at 2 (a corner), uneven spacing */
	qaws_scalar knots[13] = { 0, 0, 0, 0, 1.0f, 2.0f, 2.0f, 2.0f, 2.75f, 4, 4, 4, 4 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_desc desc;
	qaws_exact_report rep;
	viewport vc = { 20, 80, 640, 400, -0.5, 10.5, -0.6, 4.4 };
	viewport vb = { 680, 80, 580, 400, 0.5, 16.5, 0, 1700 };
	unsigned int s, i, p;
	char buf[300];
	svg sv;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = 9;
	d.knots = knots;
	d.knot_count = 13;
	d.weights = ws;
	d.weight_count = 9;
	qaws_curve_create_nurbs(&d, &c);
	qaws_exact_desc_default(&desc);
	qaws_exact_curve_prepare(&desc, c, &e, &rep);
	svg_open(&sv, "showcase/exact4_nurbs.svg", 1280, 540, "Exact NURBS: local Bezier extraction and the integer budget",
		"Left: a cubic NURBS (weights 1..6, a triple knot making a corner) split exactly into one integer rational Bezier per span (colored control polygons). Right: storage width at full precision.");
	svg_panel(&sv, &vc, "NURBS control polygon (grey) and the exact Bezier spans");
	{
		double pg[2 * 9];
		for (i = 0; i < 9; i++)
		{
			pg[2 * i] = vx(&vc, cps[2 * i]);
			pg[2 * i + 1] = vy(&vc, cps[2 * i + 1]);
		}
		svg_polyline(&sv, pg, 9, "#8c959f", 1.0, 0.8, 1);
		for (i = 0; i < 9; i++)
			svg_circle(&sv, pg[2 * i], pg[2 * i + 1], 2 + ws[i], "#ffffff", "#8c959f");
	}
	for (s = 0; s < qaws_exact_curve_span_count(e); s++)
	{
		double pts[2 * 17], w[17], t0, t1, xy[2 * 101], pg[2 * 17];
		unsigned int deg;
		qaws_exact_curve_span_bezier(e, s, &deg, &t0, &t1, pts, w);
		for (i = 0; i <= deg; i++)
		{
			pg[2 * i] = vx(&vc, pts[2 * i]);
			pg[2 * i + 1] = vy(&vc, pts[2 * i + 1]);
		}
		svg_polyline(&sv, pg, (int)deg + 1, span_colors[s % 6], 1.4, 0.9, 0);
		for (i = 0; i <= deg; i++)
			svg_circle(&sv, pg[2 * i], pg[2 * i + 1], 2.5, span_colors[s % 6], "#ffffff");
		for (i = 0; i <= 100; i++)
		{
			double out[2];
			qaws_exact_curve_evaluate(e, t0 + (t1 - t0) * i / 100.0, 0, out, NULL);
			xy[2 * i] = vx(&vc, out[0]);
			xy[2 * i + 1] = vy(&vc, out[1]);
		}
		svg_polyline(&sv, xy, 101, span_colors[s % 6], 3.0, 1, 0);
	}
	sprintf(buf, "%u spans; homogeneous integers up to %u bits; exact equality with Mathematica in test 67", qaws_exact_curve_span_count(e),
		rep.storage_bits);
	svg_text(&sv, vc.x0 + 12, vc.y0 + vc.h - 12, 12, "#57606a", "start", buf);

	/* the budget: random full-precision NURBS of every degree */
	svg_panel(&sv, &vb, "storage bits by degree; color: highest derivative within 2048 bits");
	{
		static char const* const order_colors[5] = { "#d0d7de", "#bc4c00", "#f2a900", "#8250df", "#1a7f37" };
		uint64_t state = 0x853C49E6748FEA9Bull;
		for (p = 1; p <= 16; p++)
		{
			qaws_scalar rc[24 * 2], rw[24], rk[48];
			unsigned int n = p + 4, nk = n + p + 1, order = 0;
			qaws_nurbs_desc rd;
			qaws_curve* rcu = NULL;
			qaws_exact_curve* re = NULL;
			qaws_exact_report rr;
			for (i = 0; i < n * 2; i++)
			{
				state ^= state << 13; state ^= state >> 7; state ^= state << 17;
				rc[i] = (qaws_scalar)ldexp((double)(int64_t)(state % (1u << 26)) - (1 << 25), -20);
			}
			for (i = 0; i < n; i++)
			{
				state ^= state << 13; state ^= state >> 7; state ^= state << 17;
				rw[i] = (qaws_scalar)(1 + (double)(state % 16000000u));
			}
			for (i = 0; i < nk; i++)
			{
				state ^= state << 13; state ^= state >> 7; state ^= state << 17;
				rk[i] = i <= p ? 0 : (i >= n ? 1 : (qaws_scalar)ldexp((double)(1 + state % ((1u << 24) - 2)), -24));
			}
			for (i = p + 1; i < n; i++)
			{
				unsigned int j;
				for (j = i; j > p + 1 && rk[j - 1] > rk[j]; j--)
				{
					qaws_scalar tmp = rk[j];
					rk[j] = rk[j - 1];
					rk[j - 1] = tmp;
				}
			}
			memset(&rd, 0, sizeof(rd));
			rd.dimension = QAWS_DIMENSION_2D;
			rd.degree = p;
			rd.control_points = rc;
			rd.control_point_count = n;
			rd.knots = rk;
			rd.knot_count = nk;
			rd.weights = rw;
			rd.weight_count = n;
			qaws_curve_create_nurbs(&rd, &rcu);
			if (qaws_exact_curve_prepare(&desc, rcu, &re, &rr) == QAWS_STATUS_OK)
			{
				unsigned int k;
				double out[8], x0 = vx(&vb, p - 0.38), x1 = vx(&vb, p + 0.38);
				for (k = 0; k <= 3; k++)
				{
					if (qaws_exact_curve_evaluate(re, 0.37, k, out, NULL) != QAWS_STATUS_OK)
						break;
					order = k + 1;
				}
				fprintf(sv.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n", x0, vy(&vb, rr.storage_bits), x1 - x0,
					vy(&vb, 0) - vy(&vb, rr.storage_bits), order_colors[order]);
				sprintf(buf, "%u", p);
				svg_text(&sv, (x0 + x1) / 2, vb.y0 + vb.h - 8, 11, "#57606a", "middle", buf);
			}
			qaws_exact_curve_destroy(re);
			qaws_curve_destroy(rcu);
		}
		svg_text(&sv, vb.x0 + 12, vb.y0 + 46, 12, "#1a7f37", "start", "C''' fits");
		svg_text(&sv, vb.x0 + 12, vb.y0 + 62, 12, "#8250df", "start", "C'' fits");
		svg_text(&sv, vb.x0 + 12, vb.y0 + 78, 12, "#f2a900", "start", "C' fits");
		svg_text(&sv, vb.x0 + 12, vb.y0 + 94, 12, "#bc4c00", "start", "position only (higher orders: EXACT_RANGE_EXCEEDED)");
	}
	svg_close(&sv);
	printf("exact4_nurbs: %u spans, %u storage bits\n", qaws_exact_curve_span_count(e), rep.storage_bits);
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
}
