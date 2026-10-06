/* ================================================================== */
/*  13. Certified surface / surface intersection curves               */
/* ================================================================== */

static qaws_surface* ssi_tube_surface(double k, int along_x, double h)
{
	static double const ex[7] = { 2, 2, -1, -4, -1, 2, 2 }, ey[7] = { 0, 4, 2, 0, -2, -4, 0 }, ew[7] = { 2, 1, 2, 1, 2, 1, 2 };
	qaws_vec3 cps[14];
	qaws_scalar ws[14], ku[10] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 }, kv[4] = { 0, 0, 1, 1 };
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 7; i++)
		for (j = 0; j < 2; j++)
		{
			double a = k * ex[i], b = k * ey[i], t = j ? h : -h;
			qaws_vec3* p = &cps[i * 2 + j];
			if (along_x)
			{
				p->x = (qaws_scalar)t; p->y = (qaws_scalar)a; p->z = (qaws_scalar)b;
			}
			else
			{
				p->x = (qaws_scalar)a; p->y = (qaws_scalar)b; p->z = (qaws_scalar)t;
			}
			ws[i * 2 + j] = (qaws_scalar)ew[i];
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 2;
	d.v_degree = 1;
	d.control_points = cps;
	d.u_point_count = 7;
	d.v_point_count = 2;
	d.weights = ws;
	d.u_knots = ku;
	d.u_knot_count = 10;
	d.v_knots = kv;
	d.v_knot_count = 4;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

static void demo_ssi(void)
{
	static qaws_exact_ssi_point pts[4096];
	static qaws_ssi_point fpts[20000];
	static char const* const loop_colors[4] = { "#cf222e", "#1a7f37", "#8250df", "#bc4c00" };
	qaws_exact_ssi_branch br[16];
	qaws_ssi_curve fc[16];
	qaws_ssi_desc fd;
	qaws_surface* sa = ssi_tube_surface(1, 0, 4);
	qaws_surface* sb = ssi_tube_surface(0.5, 1, 6);
	qaws_exact_surface* a = NULL;
	qaws_exact_surface* b = NULL;
	unsigned int np = 0, nb = 0, nf = 0, k, i, j;
	qaws_status st;
	projection pr = { 300, 300, 46, 1.0 };
	viewport vp = { 640, 80, 520, 400, -0.1, 3.1, -0.05, 1.05 };
	char buf[300];
	clock_t t0;
	double secs;
	svg sv;
	qaws_exact_surface_prepare(NULL, sa, &a, NULL);
	qaws_exact_surface_prepare(NULL, sb, &b, NULL);
	t0 = clock();
	st = qaws_exact_surface_surface_hits(a, b, 3, pts, 4096, &np, br, 16, &nb);
	secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
	memset(&fd, 0, sizeof(fd));
	fd.surface_a = sa;
	fd.surface_b = sb;
	qaws_surface_intersect(&fd, fc, 16, &nf, fpts, 20000);

	svg_open(&sv, "showcase/exact13_ssi.svg", 1180, 540, "Certified surface / surface intersection curves",
		"Two elliptic tubes (rational NURBS conics, extruded): the thin one pierces the wide one. Each loop is a chain of certified points; between two of them lies exactly one smooth arc.");
	svg_text(&sv, 20, 96, 14, "#24292f", "start", "the tubes and their certified intersection loops");
	/* the tubes: rulings and cross-sections */
	for (k = 0; k < 2; k++)
	{
		qaws_exact_surface const* s = k == 0 ? a : b;
		for (i = 0; i <= 24; i++)
		{
			double xy[2 * 2], p[3];
			for (j = 0; j < 2; j++)
			{
				qaws_exact_surface_evaluate(s, 3.0 * i / 24, (double)j, 0, p, NULL);
				project(&pr, p[0], p[1], p[2], &xy[2 * j], &xy[2 * j + 1]);
			}
			svg_polyline(&sv, xy, 2, k == 0 ? "#8c959f" : "#0969da", 0.6, 0.6, 0);
		}
		for (j = 0; j <= 4; j++)
		{
			double xy[2 * 97], p[3];
			for (i = 0; i <= 96; i++)
			{
				qaws_exact_surface_evaluate(s, 3.0 * i / 96, j / 4.0, 0, p, NULL);
				project(&pr, p[0], p[1], p[2], &xy[2 * i], &xy[2 * i + 1]);
			}
			svg_polyline(&sv, xy, 97, k == 0 ? "#8c959f" : "#0969da", 0.8, 0.7, 0);
		}
	}
	for (k = 0; k < nb; k++)
	{
		double xy[2 * 1025];
		unsigned int n = br[k].count + (br[k].closed ? 1 : 0), m;
		for (m = 0; m < n && m < 1024; m++)
		{
			qaws_exact_ssi_point const* q = &pts[br[k].first + (m % br[k].count)];
			double p[3];
			qaws_exact_surface_evaluate(a, 0.5 * (q->u1_lo + q->u1_hi), 0.5 * (q->v1_lo + q->v1_hi), 0, p, NULL);
			project(&pr, p[0], p[1], p[2], &xy[2 * m], &xy[2 * m + 1]);
		}
		svg_polyline(&sv, xy, (int)m, loop_colors[k % 4], 2.6, 1, 0);
		for (m = 0; m < br[k].count; m++)
			svg_circle(&sv, xy[2 * m], xy[2 * m + 1], 1.8, loop_colors[k % 4], "none");
	}
	/* parameter space of the wide tube: the loops cross its seam u1 = 0 = 3 */
	svg_panel(&sv, &vp, "the loops in the wide tube's parameters (u1, v1)");
	svg_line(&sv, vx(&vp, 0), vy(&vp, 0), vx(&vp, 0), vy(&vp, 1), "#57606a", 1.0, 1);
	svg_line(&sv, vx(&vp, 3), vy(&vp, 0), vx(&vp, 3), vy(&vp, 1), "#57606a", 1.0, 1);
	svg_text(&sv, vx(&vp, 3) - 4, vy(&vp, 1) + 14, 11, "#57606a", "end", "seam u1 = 3 (= 0)");
	for (k = 0; k < nb; k++)
		for (i = 0; i < br[k].count; i++)
		{
			qaws_exact_ssi_point const* q = &pts[br[k].first + i];
			svg_circle(&sv, vx(&vp, 0.5 * (q->u1_lo + q->u1_hi)), vy(&vp, 0.5 * (q->v1_lo + q->v1_hi)), 2.2, loop_colors[k % 4], "none");
		}
	sprintf(buf, "certified: %u closed loops, %u points (status %d, %.1f s); qaws f64 marching: %u curves", nb, np, (int)st, secs, nf);
	svg_text(&sv, 20, 520, 13, "#57606a", "start", buf);
	svg_close(&sv);
	printf("exact13_ssi: %u branches (%u closed first), %u points, status %d, %.1f s; f64 %u curves\n", nb, nb ? (unsigned int)br[0].closed : 0, np, (int)st,
		secs, nf);
	qaws_exact_surface_destroy(a);
	qaws_exact_surface_destroy(b);
	qaws_surface_destroy(sa);
	qaws_surface_destroy(sb);
}
