/* ================================================================== */
/*  11. Certified 2D Booleans from the source curves                  */
/* ================================================================== */

static qaws_curve* bo_conic_curve(double dx, double dy, double sx)
{
	qaws_scalar cps[14] = { 2, 0, 2, 4, -1, 2, -4, 0, -1, -2, 2, -4, 2, 0 };
	qaws_scalar ws[7] = { 2, 1, 2, 1, 2, 1, 2 }, kn[10] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 7; i++)
	{
		cps[2 * i] = (qaws_scalar)(cps[2 * i] * sx + dx);
		cps[2 * i + 1] = (qaws_scalar)(cps[2 * i + 1] + dy);
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = cps;
	d.control_point_count = 7;
	d.weights = ws;
	d.weight_count = 7;
	d.knots = kn;
	d.knot_count = 10;
	qaws_curve_create_nurbs(&d, &c);
	return c;
}

static void bo_draw_piece(svg* sv, viewport const* v, qaws_exact_curve const* c, qaws_exact_piece const* p, char const* color)
{
	double ts, te, span, t0 = 0.5 * (p->t0_lo + p->t0_hi), t1 = 0.5 * (p->t1_lo + p->t1_hi), xy[2 * 301];
	unsigned int sc = qaws_exact_curve_span_count(c), i;
	qaws_exact_curve_span_bezier(c, 0, NULL, &ts, NULL, NULL, NULL);
	qaws_exact_curve_span_bezier(c, sc - 1, NULL, NULL, &te, NULL, NULL);
	span = te - ts;
	if (!p->reversed && t1 <= t0) t1 += span;
	if (p->reversed && t1 >= t0) t1 -= span;
	for (i = 0; i <= 300; i++)
	{
		double t = t0 + (t1 - t0) * i / 300.0, q[2];
		while (t > te) t -= span;
		while (t < ts) t += span;
		qaws_exact_curve_evaluate(c, t, 0, q, NULL);
		xy[2 * i] = vx(v, q[0]);
		xy[2 * i + 1] = vy(v, q[1]);
	}
	svg_polyline(sv, xy, 301, color, 3.0, 1, 0);
}

static void demo_boolean(void)
{
	static char const* const names[3] = { "union", "intersection", "a minus b" };
	qaws_curve* ca = bo_conic_curve(0, 0, 1);
	qaws_curve* cb = bo_conic_curve(1.5, 0.75, -0.75);
	qaws_exact_curve* ea = NULL;
	qaws_exact_curve* eb = NULL;
	qaws_exact_curve const* curves[2];
	qaws_exact_pair x[16];
	unsigned int nx = 0, op, k;
	char buf[300];
	svg sv;
	qaws_exact_curve_prepare(NULL, ca, &ea, NULL);
	qaws_exact_curve_prepare(NULL, cb, &eb, NULL);
	curves[0] = ea;
	curves[1] = eb;
	qaws_exact_curve_curve_hits(ea, eb, x, 16, &nx);
	svg_open(&sv, "showcase/exact11_boolean.svg", 1260, 480, "Certified 2D Booleans from the source curves",
		"Two NURBS conic regions (a in blue, b mirrored, in orange). Each result boundary is made of pieces of the inputs, cut at certified crossings.");
	for (op = 0; op < 3; op++)
	{
		viewport v = { 20 + 413.0 * op, 80, 393, 360, -4.6, 5.4, -4.4, 4.8 };
		qaws_exact_piece pieces[32];
		qaws_exact_loop loops[8];
		unsigned int np = 0, nl = 0;
		qaws_status st = qaws_exact_boolean_2d(ea, eb, op, pieces, 32, &np, loops, 8, &nl);
		svg_panel(&sv, &v, names[op]);
		/* the inputs, faint */
		for (k = 0; k < 2; k++)
		{
			double xy[2 * 301];
			unsigned int i;
			for (i = 0; i <= 300; i++)
			{
				double q[2];
				qaws_exact_curve_evaluate(curves[k], 3.0 * i / 300.0, 0, q, NULL);
				xy[2 * i] = vx(&v, q[0]);
				xy[2 * i + 1] = vy(&v, q[1]);
			}
			svg_polyline(&sv, xy, 301, "#d0d7de", 1.2, 1, 1);
		}
		for (k = 0; k < np; k++)
			bo_draw_piece(&sv, &v, curves[pieces[k].region], &pieces[k], pieces[k].region == 0 ? "#0969da" : "#bc4c00");
		for (k = 0; k < nx; k++)
		{
			double q[2];
			qaws_exact_curve_evaluate(ea, 0.5 * (x[k].a_lo + x[k].a_hi), 0, q, NULL);
			svg_circle(&sv, vx(&v, q[0]), vy(&v, q[1]), 5, "#ffffff", "#1a7f37");
		}
		sprintf(buf, "%u loop(s), %u pieces (status %d)", nl, np, (int)st);
		svg_text(&sv, v.x0 + 12, v.y0 + v.h - 12, 12, "#57606a", "start", buf);
	}
	sprintf(buf, "%u certified crossings (green), parameters enclosed to about 1e-16; every piece classified by an exact winding number", nx);
	svg_text(&sv, 20, 465, 13, "#57606a", "start", buf);
	svg_close(&sv);
	printf("exact11_boolean: %u crossings\n", nx);
	qaws_exact_curve_destroy(ea);
	qaws_exact_curve_destroy(eb);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
}
