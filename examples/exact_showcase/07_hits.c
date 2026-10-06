/* ================================================================== */
/*  7. Certified intersections near a tangency                        */
/* ================================================================== */

#define EH_K 48

static void demo_hits(void)
{
	/* y(t) = 12 t^2 (1 - t)^2 (max 3/4 at t = 1/2), x(t) = 4 t */
	qaws_scalar cps[10] = { 0, 0, 1, 0, 2, 2, 3, 0, 4, 0 };
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_desc desc;
	viewport vc = { 20, 80, 380, 380, -0.2, 4.2, -0.1, 0.9 };
	viewport vn = { 440, 80, 380, 380, 0.5, EH_K + 0.5, -0.3, 3.3 };
	viewport vs = { 860, 80, 380, 380, 0.5, EH_K + 0.5, -20, 2 };
	int exact_n[EH_K + 1], f64_n[EH_K + 1], samp_n[EH_K + 1], k, i, f64_wrong = 0, samp_wrong = 0;
	double sep[EH_K + 1], wid[EH_K + 1];
	unsigned int tangent_count = 0;
	qaws_exact_hit th[4];
	char buf[300];
	svg sv;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 4;
	d.control_points = cps;
	d.control_point_count = 5;
	qaws_curve_create_bezier(&d, &c);
	qaws_exact_desc_default(&desc);
	qaws_exact_curve_prepare(&desc, c, &e, NULL);
	for (k = 1; k <= EH_K; k++)
	{
		double h = 0.75 - ldexp(1.0, -k), p0[2] = { -1, h }, p1[2] = { 5, h };
		qaws_exact_hit hits[8];
		unsigned int n = 0;
		qaws_intersection_2d ix[16];
		unsigned int nf = 0;
		qaws_scalar seg[4] = { -1, (qaws_scalar)h, 5, (qaws_scalar)h };
		qaws_bezier_desc ld;
		qaws_curve* line = NULL;
		double prev = 0;
		qaws_exact_curve_line_hits(e, p0, p1, 0, hits, 8, &n);
		exact_n[k] = (int)n;
		sep[k] = n == 2 ? hits[1].t_lo - hits[0].t_hi : 0;
		wid[k] = n == 2 ? (hits[0].t_hi - hits[0].t_lo > hits[1].t_hi - hits[1].t_lo ? hits[0].t_hi - hits[0].t_lo : hits[1].t_hi - hits[1].t_lo) : 0;
		memset(&ld, 0, sizeof(ld));
		ld.dimension = QAWS_DIMENSION_2D;
		ld.degree = 1;
		ld.control_points = seg;
		ld.control_point_count = 2;
		qaws_curve_create_bezier(&ld, &line);
		qaws_curve_find_intersections_2d(c, line, ix, 16, &nf);
		f64_n[k] = (int)nf;
		qaws_curve_destroy(line);
		/* sign changes of y(t) - h over 100000 samples */
		samp_n[k] = 0;
		for (i = 0; i <= 100000; i++)
		{
			double t = (i + 0.37) / 100000.0, y = 12 * t * t * (1 - t) * (1 - t) - h;
			if (i > 0 && ((prev < 0) != (y < 0)))
				samp_n[k]++;
			prev = y;
		}
		f64_wrong += f64_n[k] != exact_n[k];
		samp_wrong += samp_n[k] != exact_n[k];
	}
	{
		double p0[2] = { -1, 0.75 }, p1[2] = { 5, 0.75 };
		qaws_exact_curve_line_hits(e, p0, p1, 0, th, 4, &tangent_count);
	}

	svg_open(&sv, "showcase/exact7_hits.svg", 1260, 520, "Certified intersections near a tangency",
		"A quartic bump y = 12 t^2 (1 - t)^2 (top 3/4) cut by the lines y = 3/4 - 2^-k: always two crossings, closing in like 2^(-k/2).");
	svg_panel(&sv, &vc, "the curve and lines y = 3/4 - 2^-k (k = 1..6)");
	{
		double xy[2 * 201];
		for (i = 0; i <= 200; i++)
		{
			double out[2];
			qaws_exact_curve_evaluate(e, i / 200.0, 0, out, NULL);
			xy[2 * i] = vx(&vc, out[0]);
			xy[2 * i + 1] = vy(&vc, out[1]);
		}
		for (k = 1; k <= 6; k++)
			svg_line(&sv, vx(&vc, -0.2), vy(&vc, 0.75 - ldexp(1.0, -k)), vx(&vc, 4.2), vy(&vc, 0.75 - ldexp(1.0, -k)), "#8c959f", 0.8, 1);
		svg_line(&sv, vx(&vc, -0.2), vy(&vc, 0.75), vx(&vc, 4.2), vy(&vc, 0.75), "#1a7f37", 1.2, 1);
		svg_polyline(&sv, xy, 201, "#24292f", 2.4, 1, 0);
		sprintf(buf, "y = 3/4: exact tangent point at t = %g (%u hit)", th[0].t_lo, tangent_count);
		svg_text(&sv, vc.x0 + 12, vc.y0 + vc.h - 12, 12, "#1a7f37", "start", buf);
	}
	svg_panel(&sv, &vn, "intersections found, by k");
	for (k = 0; k <= 3; k++)
	{
		sprintf(buf, "%d", k);
		svg_line(&sv, vn.x0, vy(&vn, k), vn.x0 + vn.w, vy(&vn, k), "#d0d7de", 0.6, 1);
		svg_text(&sv, vn.x0 + 4, vy(&vn, k) - 2, 10, "#8c959f", "start", buf);
	}
	for (k = 1; k <= EH_K; k++)
	{
		svg_circle(&sv, vx(&vn, k), vy(&vn, samp_n[k] + 0.12), 3.0, "#bc4c00", "none");
		svg_circle(&sv, vx(&vn, k), vy(&vn, f64_n[k] - 0.12), 3.0, "#0969da", "none");
		svg_circle(&sv, vx(&vn, k), vy(&vn, exact_n[k]), 3.0, "#1a7f37", "none");
	}
	for (k = 8; k <= EH_K; k += 8)
	{
		sprintf(buf, "k = %d", k);
		svg_text(&sv, vx(&vn, k), vn.y0 + vn.h - 6, 10, "#8c959f", "middle", buf);
	}
	svg_text(&sv, vn.x0 + vn.w - 10, vn.y0 + 46, 12, "#1a7f37", "end", "exact: 2 for every k");
	svg_text(&sv, vn.x0 + vn.w - 10, vn.y0 + 62, 12, "#0969da", "end", "qaws f64 curve/curve intersection");
	svg_text(&sv, vn.x0 + vn.w - 10, vn.y0 + 78, 12, "#bc4c00", "end", "f64 sign changes, 10^5 samples");

	svg_panel(&sv, &vs, "log10 in t: root separation and enclosure width");
	for (k = -20; k <= 0; k += 4)
	{
		sprintf(buf, "1e%d", k);
		svg_line(&sv, vs.x0, vy(&vs, k), vs.x0 + vs.w, vy(&vs, k), "#d0d7de", 0.6, 1);
		svg_text(&sv, vs.x0 + 4, vy(&vs, k) - 2, 10, "#8c959f", "start", buf);
	}
	{
		double a[2 * EH_K], b[2 * EH_K];
		int na = 0, nb = 0;
		for (k = 1; k <= EH_K; k++)
		{
			if (sep[k] > 0)
			{
				a[2 * na] = vx(&vs, k);
				a[2 * na + 1] = vy(&vs, log10(sep[k]));
				na++;
			}
			if (wid[k] > 0)
			{
				b[2 * nb] = vx(&vs, k);
				b[2 * nb + 1] = vy(&vs, log10(wid[k]));
				nb++;
			}
		}
		svg_polyline(&sv, a, na, "#24292f", 2.0, 1, 0);
		svg_polyline(&sv, b, nb, "#1a7f37", 2.0, 1, 0);
	}
	svg_text(&sv, vs.x0 + vs.w - 10, vy(&vs, -9.5), 12, "#24292f", "end", "gap between the two crossings");
	svg_text(&sv, vs.x0 + vs.w - 10, vy(&vs, -14.8), 12, "#1a7f37", "end", "certified interval width (2 ulps)");
	sprintf(buf, "k = 1..%d: f64 curve/curve intersection wrong for %d lines, f64 sampling for %d; exact: every pair isolated and enclosed", EH_K, f64_wrong,
		samp_wrong);
	svg_text(&sv, 20, 495, 13, "#57606a", "start", buf);
	svg_close(&sv);
	printf("exact7_hits: f64 wrong %d, sampling wrong %d of %d; tangent hits %u at %.17g\n", f64_wrong, samp_wrong, EH_K, tangent_count, th[0].t_lo);
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
}
