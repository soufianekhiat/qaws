/* ================================================================== */
/*  9. Certified self-intersections: a star and a knot                */
/* ================================================================== */

static void demo_self_hits(void)
{
	qaws_scalar star[6 * 2], knot[15 * 3], kn_star[9], kn_knot[19];
	qaws_bspline_desc d;
	qaws_curve* cs = NULL;
	qaws_curve* ck = NULL;
	qaws_exact_curve* es = NULL;
	qaws_exact_curve* ek = NULL;
	qaws_exact_desc desc;
	qaws_exact_pair hs[32], hk[32];
	qaws_intersection_2d f2[64];
	qaws_intersection_3d f3[64];
	unsigned int ns = 0, nk = 0, nf2 = 0, nf3 = 0, i;
	qaws_status ss, sk;
	viewport va = { 20, 80, 400, 400, -5.5, 5.5, -5.5, 5.5 };
	viewport vb = { 460, 80, 400, 400, -4.2, 4.2, -4.2, 4.2 };
	char buf[300];
	svg sv;
	qaws_exact_desc_default(&desc);
	/* a pentagram traced by a quadratic B-spline (open, clamped) */
	for (i = 0; i < 5; i++)
	{
		double a = 3.14159265358979 / 2 + 4 * 3.14159265358979 * (i % 5) / 5;
		star[2 * i] = (qaws_scalar)ldexp(nearbyint(ldexp(4.8 * cos(a), 12)), -12);
		star[2 * i + 1] = (qaws_scalar)ldexp(nearbyint(ldexp(4.8 * sin(a), 12)), -12);
	}
	star[10] = star[0] * 0.9f;
	star[11] = star[1] * 0.9f;
	for (i = 0; i < 9; i++)
		kn_star[i] = (qaws_scalar)(i < 3 ? 0 : (i > 5 ? 4 : i - 2));
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = star;
	d.control_point_count = 6;
	d.knots = kn_star;
	d.knot_count = 9;
	qaws_curve_create_bspline(&d, &cs);
	qaws_exact_curve_prepare(&desc, cs, &es, NULL);
	ss = qaws_exact_curve_self_hits(es, hs, 32, &ns);
	qaws_curve_find_self_intersections_2d(cs, f2, 64, &nf2);
	/* a trefoil knot: a closed (periodic) cubic B-spline in space */
	for (i = 0; i < 12; i++)
	{
		double t = 2 * 3.14159265358979 * i / 12;
		knot[3 * i] = (qaws_scalar)ldexp(nearbyint(ldexp(sin(t) + 2 * sin(2 * t), 12)), -12);
		knot[3 * i + 1] = (qaws_scalar)ldexp(nearbyint(ldexp(cos(t) - 2 * cos(2 * t), 12)), -12);
		knot[3 * i + 2] = (qaws_scalar)ldexp(nearbyint(ldexp(-1.5 * sin(3 * t), 12)), -12);
	}
	for (i = 12; i < 15; i++)
	{
		knot[3 * i] = knot[3 * (i - 12)];
		knot[3 * i + 1] = knot[3 * (i - 12) + 1];
		knot[3 * i + 2] = knot[3 * (i - 12) + 2];
	}
	for (i = 0; i < 19; i++)
		kn_knot[i] = (qaws_scalar)i;
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = knot;
	d.control_point_count = 15;
	d.knots = kn_knot;
	d.knot_count = 19;
	qaws_curve_create_bspline(&d, &ck);
	qaws_exact_curve_prepare(&desc, ck, &ek, NULL);
	sk = qaws_exact_curve_self_hits(ek, hk, 32, &nk);
	qaws_curve_find_self_intersections_3d(ck, f3, 64, &nf3);

	svg_open(&sv, "showcase/exact9_self_hits.svg", 880, 540, "Certified self-intersections",
		"Left: a B-spline star, crossings isolated exactly. Right: a trefoil knot, whose shadow crosses itself while the curve never does.");
	svg_panel(&sv, &va, "2D: crossings of a quadratic B-spline star");
	{
		double xy[2 * 701], t0 = 0, t1 = 4;
		for (i = 0; i <= 700; i++)
		{
			double out[2];
			qaws_exact_curve_evaluate(es, t0 + (t1 - t0) * i / 700.0, 0, out, NULL);
			xy[2 * i] = vx(&va, out[0]);
			xy[2 * i + 1] = vy(&va, out[1]);
		}
		svg_polyline(&sv, xy, 701, "#24292f", 2.0, 1, 0);
		for (i = 0; i < ns; i++)
		{
			double p[2];
			qaws_exact_curve_evaluate(es, 0.5 * (hs[i].a_lo + hs[i].a_hi), 0, p, NULL);
			svg_circle(&sv, vx(&va, p[0]), vy(&va, p[1]), 6, "none", "#1a7f37");
		}
		sprintf(buf, "exact: %u crossings (status %d); qaws f64 finder: %u; Mathematica: 5", ns, (int)ss, nf2);
		svg_text(&sv, va.x0 + 12, va.y0 + va.h - 12, 12, "#57606a", "start", buf);
	}
	svg_panel(&sv, &vb, "3D: a trefoil knot (xy shadow, depth shaded)");
	{
		double t0 = 3, t1 = 15;
		for (i = 0; i < 600; i++)
		{
			double p[3], q[3];
			char col[32];
			qaws_exact_curve_evaluate(ek, t0 + (t1 - t0) * i / 600.0, 0, p, NULL);
			qaws_exact_curve_evaluate(ek, t0 + (t1 - t0) * (i + 1) / 600.0, 0, q, NULL);
			heat(0.5 + 0.3 * p[2], col);
			svg_line(&sv, vx(&vb, p[0]), vy(&vb, p[1]), vx(&vb, q[0]), vy(&vb, q[1]), col, 3.0, 1);
		}
		sprintf(buf, "exact: %u self-intersections (status %d), the 3 shadow crossings rejected", nk, (int)sk);
		svg_text(&sv, vb.x0 + 12, vb.y0 + vb.h - 28, 12, "#57606a", "start", buf);
		sprintf(buf, "qaws f64 finder: %u (Mathematica: 0)", nf3);
		svg_text(&sv, vb.x0 + 12, vb.y0 + vb.h - 12, 12, "#57606a", "start", buf);
	}
	svg_close(&sv);
	printf("exact9_self_hits: star %u (status %d, f64 %u), knot %u (status %d, f64 %u)\n", ns, (int)ss, nf2, nk, (int)sk, nf3);
	qaws_exact_curve_destroy(es);
	qaws_exact_curve_destroy(ek);
	qaws_curve_destroy(cs);
	qaws_curve_destroy(ck);
}
