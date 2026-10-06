/* ================================================================== */
/*  12. Certified curve / surface intersections                       */
/* ================================================================== */

static void demo_curve_surface(void)
{
	qaws_vec3 cps[16];
	qaws_scalar ws[16], kn[7] = { 0, 0, 0, 1, 2, 2, 2 }, ccp[13 * 3], ckn[17];
	qaws_surface_nurbs_desc d;
	qaws_bspline_desc cd;
	qaws_surface* s = NULL;
	qaws_curve* c = NULL;
	qaws_exact_surface* es = NULL;
	qaws_exact_curve* ec = NULL;
	qaws_exact_curve_surface_hit hits[64];
	unsigned int nh = 0, i, j, k;
	qaws_status st;
	projection pr = { 330, 250, 70, 1.0 };
	char buf[300];
	clock_t t0;
	double secs;
	svg sv;
	/* the wave of figure 10 */
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			cps[i * 4 + j].x = (qaws_scalar)i;
			cps[i * 4 + j].y = (qaws_scalar)j;
			cps[i * 4 + j].z = (qaws_scalar)(((i + j) % 2) ? 1.25 : -0.75) * (qaws_scalar)(i == 0 || i == 3 ? 0.5 : 1);
			ws[i * 4 + j] = (qaws_scalar)((i == 1 && j == 2) ? 4 : 1);
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cps;
	d.u_point_count = 4;
	d.v_point_count = 4;
	d.weights = ws;
	d.u_knots = kn;
	d.u_knot_count = 7;
	d.v_knots = kn;
	d.v_knot_count = 7;
	qaws_surface_create_nurbs(&d, &s);
	/* a spiral weaving up and down through it (a cubic B-spline, 10 spans) */
	for (k = 0; k < 13; k++)
	{
		double th = 3.14159265358979 * k / 4;
		ccp[3 * k] = (qaws_scalar)ldexp(nearbyint(ldexp(1.5 + 1.15 * cos(th), 12)), -12);
		ccp[3 * k + 1] = (qaws_scalar)ldexp(nearbyint(ldexp(1.5 + 1.15 * sin(th), 12)), -12);
		ccp[3 * k + 2] = (qaws_scalar)ldexp(nearbyint(ldexp(1.6 * cos(1.5 * th), 12)), -12);
	}
	for (k = 0; k < 17; k++)
		ckn[k] = (qaws_scalar)(k < 4 ? 0 : (k > 12 ? 10 : k - 3));
	memset(&cd, 0, sizeof(cd));
	cd.dimension = QAWS_DIMENSION_3D;
	cd.degree = 3;
	cd.control_points = ccp;
	cd.control_point_count = 13;
	cd.knots = ckn;
	cd.knot_count = 17;
	qaws_curve_create_bspline(&cd, &c);
	qaws_exact_surface_prepare(NULL, s, &es, NULL);
	qaws_exact_curve_prepare(NULL, c, &ec, NULL);
	t0 = clock();
	st = qaws_exact_curve_surface_hits(ec, es, hits, 64, &nh);
	secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
	/* hits by curve parameter */
	for (i = 1; i < nh; i++)
	{
		qaws_exact_curve_surface_hit key = hits[i];
		int m = (int)i - 1;
		while (m >= 0 && hits[m].t_lo > key.t_lo)
		{
			hits[m + 1] = hits[m];
			m--;
		}
		hits[m + 1] = key;
	}

	svg_open(&sv, "showcase/exact12_curve_surface.svg", 1180, 560, "Certified curve / surface intersections",
		"A cubic B-spline spiral through a rational NURBS wave: exact Bernstein subdivision, Miranda existence and an exact uniqueness test per crossing.");
	svg_text(&sv, 20, 96, 14, "#24292f", "start", "the curve alternates sides at each certified crossing");
	/* the surface: iso-parameter lines */
	for (k = 0; k <= 16; k++)
	{
		double xy[2 * 41], xy2[2 * 41];
		for (i = 0; i <= 40; i++)
		{
			double p[3], q[3];
			qaws_exact_surface_evaluate(es, 2.0 * k / 16, 2.0 * i / 40, 0, p, NULL);
			qaws_exact_surface_evaluate(es, 2.0 * i / 40, 2.0 * k / 16, 0, q, NULL);
			project(&pr, p[0], p[1], p[2], &xy[2 * i], &xy[2 * i + 1]);
			project(&pr, q[0], q[1], q[2], &xy2[2 * i], &xy2[2 * i + 1]);
		}
		svg_polyline(&sv, xy, 41, "#8c959f", 0.7, 0.8, 0);
		svg_polyline(&sv, xy2, 41, "#8c959f", 0.7, 0.8, 0);
	}
	/* the curve, its pieces between crossings in alternating colors */
	{
		double tprev = 0;
		for (k = 0; k <= nh; k++)
		{
			double tnext = k < nh ? 0.5 * (hits[k].t_lo + hits[k].t_hi) : 10, xy[2 * 201];
			for (i = 0; i <= 200; i++)
			{
				double p[3];
				qaws_exact_curve_evaluate(ec, tprev + (tnext - tprev) * i / 200.0, 0, p, NULL);
				project(&pr, p[0], p[1], p[2], &xy[2 * i], &xy[2 * i + 1]);
			}
			svg_polyline(&sv, xy, 201, (k % 2) ? "#bc4c00" : "#0969da", 2.6, 1, 0);
			tprev = tnext;
		}
	}
	for (k = 0; k < nh; k++)
	{
		double p[3], x, y;
		qaws_exact_curve_evaluate(ec, 0.5 * (hits[k].t_lo + hits[k].t_hi), 0, p, NULL);
		project(&pr, p[0], p[1], p[2], &x, &y);
		svg_circle(&sv, x, y, 5, "#ffffff", "#1a7f37");
	}
	/* the enclosures */
	svg_text(&sv, 700, 96, 14, "#24292f", "start", "certified enclosures (curve t; surface u, v)");
	for (k = 0; k < nh && k < 18; k++)
	{
		sprintf(buf, "t %.12f  u %.12f  v %.12f   width %.0e", 0.5 * (hits[k].t_lo + hits[k].t_hi), 0.5 * (hits[k].u_lo + hits[k].u_hi),
			0.5 * (hits[k].v_lo + hits[k].v_hi), hits[k].t_hi - hits[k].t_lo);
		svg_text(&sv, 700, 124 + 20 * k, 12, "#57606a", "start", buf);
	}
	sprintf(buf, "%u crossings certified (status %d) in %.1f s; 10 curve spans x 4 patches", nh, (int)st, secs);
	svg_text(&sv, 20, 540, 13, "#57606a", "start", buf);
	svg_close(&sv);
	printf("exact12_curve_surface: %u hits, status %d, %.1f s\n", nh, (int)st, secs);
	qaws_exact_surface_destroy(es);
	qaws_exact_curve_destroy(ec);
	qaws_surface_destroy(s);
	qaws_curve_destroy(c);
}
