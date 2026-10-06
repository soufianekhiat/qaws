/* ================================================================== */
/*  8. Certified curve / curve intersections near a tangency          */
/* ================================================================== */

#define EC_K 28

static void demo_curve_hits(void)
{
	/* A: (4 s, 12 s^2 (1 - s)^2), top 3/4 at s = 1/2.
	   B: (4 r, 3/2 - 2^-k - 12 r^2 (1 - r)^2), bottom 3/4 - 2^-k at r = 1/2. */
	qaws_scalar acp[10] = { 0, 0, 1, 0, 2, 2, 3, 0, 4, 0 };
	qaws_bezier_desc d;
	qaws_curve* ca = NULL;
	qaws_exact_curve* ea = NULL;
	qaws_exact_desc desc;
	viewport vc = { 20, 80, 380, 380, -0.2, 4.2, -0.1, 1.6 };
	viewport vn = { 440, 80, 380, 380, 0.5, EC_K + 0.5, -0.3, 3.3 };
	viewport vz = { 860, 80, 380, 380, 1.99, 2.01, 0.7499, 0.75005 };
	int exact_n[EC_K + 1], f64_n[EC_K + 1], k, i, f64_wrong = 0, failed = 0;
	char buf[300];
	svg sv;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 4;
	d.control_points = acp;
	d.control_point_count = 5;
	qaws_curve_create_bezier(&d, &ca);
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -29;   /* every 2^-k, k <= 28, on the lattice (x up to 4) */
	desc.coord_bits = 32;
	qaws_exact_curve_prepare(&desc, ca, &ea, NULL);
	svg_open(&sv, "showcase/exact8_curve_hits.svg", 1260, 520, "Certified curve / curve intersections near a tangency",
		"Two quartic bumps, the upper one dipping 2^-k below the top of the lower one: always two crossings, 2^(-k/2) apart (coordinates on the 2^-29 lattice).");
	for (k = 1; k <= EC_K; k++)
	{
		double h = 1.5 - ldexp(1.0, -k);
		qaws_scalar bcp[10];
		qaws_curve* cb = NULL;
		qaws_exact_curve* eb = NULL;
		qaws_exact_pair hits[8];
		qaws_intersection_2d ix[16];
		unsigned int n = 0, nf = 0;
		for (i = 0; i < 5; i++)
		{
			bcp[2 * i] = acp[2 * i];
			bcp[2 * i + 1] = (qaws_scalar)(h - acp[2 * i + 1]);
		}
		d.control_points = bcp;
		qaws_curve_create_bezier(&d, &cb);
		qaws_exact_curve_prepare(&desc, cb, &eb, NULL);
		if (qaws_exact_curve_curve_hits(ea, eb, hits, 8, &n) != QAWS_STATUS_OK)
			failed++;
		exact_n[k] = (int)n;
		qaws_curve_find_intersections_2d(ca, cb, ix, 16, &nf);
		f64_n[k] = (int)nf;
		f64_wrong += f64_n[k] != exact_n[k];
		if (k <= 5)
		{
			double xy[2 * 201];
			for (i = 0; i <= 200; i++)
			{
				double out[2];
				qaws_exact_curve_evaluate(eb, i / 200.0, 0, out, NULL);
				xy[2 * i] = vx(&vc, out[0]);
				xy[2 * i + 1] = vy(&vc, out[1]);
			}
			if (k == 1)
				svg_panel(&sv, &vc, "the lower bump and the upper ones for k = 1..5");
			svg_polyline(&sv, xy, 201, "#8c959f", 1.2, 1, 0);
		}
		if (k == 16)
		{
			/* zoom: both curves near the touching point and the certified crossings */
			double xa[2 * 201], xb[2 * 201];
			unsigned int j;
			svg_panel(&sv, &vz, "k = 16, zoomed: the two certified crossings");
			for (i = 0; i <= 200; i++)
			{
				double t = 0.4975 + 0.005 * i / 200.0, pa[2], pb[2];
				qaws_exact_curve_evaluate(ea, t, 0, pa, NULL);
				qaws_exact_curve_evaluate(eb, t, 0, pb, NULL);
				xa[2 * i] = vx(&vz, pa[0]);
				xa[2 * i + 1] = vy(&vz, pa[1]);
				xb[2 * i] = vx(&vz, pb[0]);
				xb[2 * i + 1] = vy(&vz, pb[1]);
			}
			svg_polyline(&sv, xa, 201, "#24292f", 2.0, 1, 0);
			svg_polyline(&sv, xb, 201, "#0969da", 2.0, 1, 0);
			for (j = 0; j < n; j++)
			{
				double p[2];
				qaws_exact_curve_evaluate(ea, 0.5 * (hits[j].a_lo + hits[j].a_hi), 0, p, NULL);
				svg_circle(&sv, vx(&vz, p[0]), vy(&vz, p[1]), 6, "none", "#1a7f37");
				sprintf(buf, "s in [%.15f, %.15f]", hits[j].a_lo, hits[j].a_hi);
				svg_text(&sv, vz.x0 + 12, vz.y0 + vz.h - 30 + 16 * j, 11, "#1a7f37", "start", buf);
			}
		}
		qaws_exact_curve_destroy(eb);
		qaws_curve_destroy(cb);
	}
	{
		double xy[2 * 201];
		for (i = 0; i <= 200; i++)
		{
			double out[2];
			qaws_exact_curve_evaluate(ea, i / 200.0, 0, out, NULL);
			xy[2 * i] = vx(&vc, out[0]);
			xy[2 * i + 1] = vy(&vc, out[1]);
		}
		svg_polyline(&sv, xy, 201, "#24292f", 2.4, 1, 0);
	}
	svg_panel(&sv, &vn, "intersections found, by k");
	for (k = 0; k <= 3; k++)
	{
		sprintf(buf, "%d", k);
		svg_line(&sv, vn.x0, vy(&vn, k), vn.x0 + vn.w, vy(&vn, k), "#d0d7de", 0.6, 1);
		svg_text(&sv, vn.x0 + 4, vy(&vn, k) - 2, 10, "#8c959f", "start", buf);
	}
	for (k = 1; k <= EC_K; k++)
	{
		svg_circle(&sv, vx(&vn, k), vy(&vn, f64_n[k] - 0.12), 3.0, "#0969da", "none");
		svg_circle(&sv, vx(&vn, k), vy(&vn, exact_n[k]), 3.0, "#1a7f37", "none");
	}
	for (k = 8; k <= EC_K; k += 8)
	{
		sprintf(buf, "k = %d", k);
		svg_text(&sv, vx(&vn, k), vn.y0 + vn.h - 6, 10, "#8c959f", "middle", buf);
	}
	svg_text(&sv, vn.x0 + vn.w - 10, vn.y0 + 46, 12, "#1a7f37", "end", "exact: 2 for every k");
	svg_text(&sv, vn.x0 + vn.w - 10, vn.y0 + 62, 12, "#0969da", "end", "qaws f64 curve/curve intersection");
	sprintf(buf, "k = 1..%d: f64 intersection count wrong for %d pairs; exact: implicitization, root isolation, %d uncertified", EC_K, f64_wrong, failed);
	svg_text(&sv, 20, 495, 13, "#57606a", "start", buf);
	svg_close(&sv);
	printf("exact8_curve_hits: f64 wrong %d of %d, exact failures %d\n", f64_wrong, EC_K, failed);
	qaws_exact_curve_destroy(ea);
	qaws_curve_destroy(ca);
}
