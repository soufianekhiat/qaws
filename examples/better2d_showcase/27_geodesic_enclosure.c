/* Figure 27: validated enclosures of geodesics (qaws_surface_geodesic_enclose,
   qaws_geodesic_enclosure_cross). On a bicubic B-spline bump, three
   geodesics, each cut into steps with a box in (u, v) proven to hold it
   (interval arithmetic, Lohner's moving frame); two of them proven to cross
   exactly once, the third proven to meet neither. Zoomed in (u, v): the
   step boxes round the crossing with the certified crossing box, and the
   end box of one geodesic with the float geodesic's end inside it. */

static qaws_surface* s27_bump(void)
{
	qaws_vec3 cp[36];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
		{
			double x = i / 5.0, y = j / 5.0;
			cp[6 * i + j].x = (qaws_scalar)(2 * x - 1); cp[6 * i + j].y = (qaws_scalar)(2 * y - 1);
			cp[6 * i + j].z = (qaws_scalar)(1.0 * exp(-8 * ((x - 0.5) * (x - 0.5) + (y - 0.45) * (y - 0.45))));
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 3; d.v_degree = 3; d.control_points = cp; d.u_point_count = 6; d.v_point_count = 6;
	qaws_surface_create_bspline(&s, &d);
	return s;
}

/* a (u, v) box as a rectangle on the surface */
static void s27_box_3d(ssd_scene* sc, qaws_surface const* f, qaws_vec2 lo, qaws_vec2 hi, int const* rgb)
{
	double uv[10];
	uv[0] = lo.x; uv[1] = lo.y; uv[2] = hi.x; uv[3] = lo.y; uv[4] = hi.x; uv[5] = hi.y; uv[6] = lo.x; uv[7] = hi.y; uv[8] = lo.x; uv[9] = lo.y;
	ssd_uv_line(sc, f, uv, 5, rgb, 0.8);
}

static void s27_box_2d(svg* s, viewport const* v, double ulo, double vlo, double uhi, double vhi, char const* color, double width)
{
	double xy[10];
	/* (clipped to the window) */
	if (ulo < v->xmin) ulo = v->xmin;
	if (uhi > v->xmax) uhi = v->xmax;
	if (vlo < v->ymin) vlo = v->ymin;
	if (vhi > v->ymax) vhi = v->ymax;
	if (ulo > uhi || vlo > vhi)
		return;
	xy[0] = vx(v, ulo); xy[1] = vy(v, vlo); xy[2] = vx(v, uhi); xy[3] = vy(v, vlo); xy[4] = vx(v, uhi); xy[5] = vy(v, vhi);
	xy[6] = vx(v, ulo); xy[7] = vy(v, vhi); xy[8] = xy[0]; xy[9] = xy[1];
	svg_polyline(s, xy, 5, color, width, 0.9, 0);
}

/* the float geodesic in (u, v) within the window, as a polyline */
static void s27_float_2d(svg* s, viewport const* v, qaws_surface const* f, double const* g, double ta, double tb, char const* color)
{
	enum { N = 400 };
	static double xy[2 * (N + 1)];
	int i, n = 0;
	for (i = 0; i <= N; i++)
	{
		qaws_vec2 q, p, d;
		p.x = (qaws_scalar)g[0]; p.y = (qaws_scalar)g[1]; d.x = (qaws_scalar)g[2]; d.y = (qaws_scalar)g[3];
		qaws_surface_exp(&q, NULL, f, p, d, (qaws_scalar)(ta + (tb - ta) * i / N), 0);
		if (q.x < v->xmin || q.x > v->xmax || q.y < v->ymin || q.y > v->ymax)
		{
			if (n > 1) svg_polyline(s, xy, n, color, 1.6, 0.9, 0);
			n = 0;
			continue;
		}
		xy[2 * n] = vx(v, q.x); xy[2 * n + 1] = vy(v, q.y);
		n++;
	}
	if (n > 1) svg_polyline(s, xy, n, color, 1.6, 0.9, 0);
}

/* the step boxes within the window; the stretch of length they cover */
static unsigned int s27_steps_2d(svg* s, viewport const* v, qaws_geodesic_enclosure const* e, char const* color, double* ta, double* tb)
{
	unsigned int i, n = 0;
	*ta = 1e300; *tb = -1e300;
	for (i = 0; i < qaws_geodesic_enclosure_get_step_count(e); i++)
	{
		qaws_vec2 lo, hi;
		qaws_scalar t0, t1;
		qaws_geodesic_enclosure_get_step(&t0, &t1, &lo, &hi, e, i);
		if (hi.x < v->xmin || lo.x > v->xmax || hi.y < v->ymin || lo.y > v->ymax)
			continue;
		s27_box_2d(s, v, lo.x, lo.y, hi.x, hi.y, color, 0.7);
		if (t0 < *ta) *ta = t0;
		if (t1 > *tb) *tb = t1;
		n++;
	}
	return n;
}

static void demo_geodesic_enclosure(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 420, ph = 420;
	static int const blue[3] = { 9, 105, 218 }, red[3] = { 207, 34, 46 }, green[3] = { 26, 127, 55 }, black[3] = { 27, 31, 36 };
	static double const G[3][5] = { { 0.1, 0.2, 0.7, 0.5, 1.6 }, { 0.2, 0.75, 0.6, -0.6, 1.4 }, { 0.1, 0.9, 1, 0.05, 1.3 } };
	static int const* const col[3] = { blue, red, green };
	static char const* const hex[3] = { "#0969da", "#cf222e", "#1a7f37" };
	qaws_surface* f[1];
	qaws_geodesic_enclosure* e[3] = { NULL, NULL, NULL };
	qaws_vec2 clo[2], chi[2];
	unsigned int i, k, n01 = 0, n02 = 0, steps = 0;
	qaws_status x01 = QAWS_STATUS_INTERNAL_ERROR, x02 = QAWS_STATUS_INTERNAL_ERROR;
	char label[240];
	if (!svg_open(&s, "showcase/b2d27_geodesic_enclosure.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + (ph + 44)),
		"Validated geodesics",
		"qaws_surface_geodesic_enclose: boxes proven to hold each geodesic (intervals, Lohner's frame); crossings certified"))
		return;
	f[0] = s27_bump();
	for (k = 0; k < 3; k++)
	{
		qaws_vec2 p, d;
		p.x = (qaws_scalar)G[k][0]; p.y = (qaws_scalar)G[k][1]; d.x = (qaws_scalar)G[k][2]; d.y = (qaws_scalar)G[k][3];
		if (qaws_surface_geodesic_enclose(&e[k], f[0], 0, p, d, (qaws_scalar)G[k][4]) == QAWS_STATUS_OK)
			steps += qaws_geodesic_enclosure_get_step_count(e[k]);
	}
	if (e[0] && e[1]) x01 = qaws_geodesic_enclosure_cross(&clo[0], &chi[0], &n01, e[0], e[1], 1);
	if (e[0] && e[2]) x02 = qaws_geodesic_enclosure_cross(&clo[1], &chi[1], &n02, e[0], e[2], 1);
	/* the bump: the geodesics, every 40th step box, the crossing */
	{
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		static double const box[1][4] = { { 0, 1, 0, 1 } };
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24; v3.y0 = 104; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -0.7, 0.75, 150, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		sc.view = &view; sc.box = box; sc.res = 70;
		sprintf(label, "three geodesics on a B-spline bump: %u steps", steps);
		ssd_panel(&s, &v3, label, x01 == QAWS_STATUS_OK && x02 == QAWS_STATUS_OK ?
			(n01 == 1 && n02 == 0 ? "blue and red cross once (proven); green meets neither (proven)" : "crossings certified") : "crossings not certified");
		ssd_surfaces(&sc, f, 1, NULL, 0);
		for (k = 0; k < 3; k++)
		{
			qaws_surface_curve* c = NULL;
			qaws_vec2 p, d;
			p.x = (qaws_scalar)G[k][0]; p.y = (qaws_scalar)G[k][1]; d.x = (qaws_scalar)G[k][2]; d.y = (qaws_scalar)G[k][3];
			if (qaws_surface_curve_create_geodesic_ray(&c, f[0], p, d, (qaws_scalar)G[k][4], 0) == QAWS_STATUS_OK)
				ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(c), col[k], 1.6);
			qaws_surface_curve_destroy(c);
			for (i = 0; e[k] && i < qaws_geodesic_enclosure_get_step_count(e[k]); i += 40)
			{
				qaws_vec2 lo, hi;
				qaws_geodesic_enclosure_get_step(NULL, NULL, &lo, &hi, e[k], i);
				s27_box_3d(&sc, f[0], lo, hi, col[k]);
			}
			ssd_dot(&sc, f[0], G[k][0], G[k][1], col[k], 3.2);
		}
		if (n01 == 1)
			ssd_dot(&sc, f[0], 0.5 * (clo[0].x + chi[0].x), 0.5 * (clo[0].y + chi[0].y), black, 4.0);
		ssd_draw(&s, &sc);
	}
	/* (u, v) round the crossing: the step boxes and the crossing's box */
	if (n01 == 1)
	{
		viewport v2;
		double cu = 0.5 * (clo[0].x + chi[0].x), cv = 0.5 * (clo[0].y + chi[0].y), r = 2.5 * fmax(chi[0].x - clo[0].x, chi[0].y - clo[0].y);
		unsigned int na, nb;
		double ta, tb, sa, sb;
		v2.x0 = 24 + (pw + 12); v2.y0 = 104; v2.w = pw; v2.h = ph;
		v2.xmin = cu - r; v2.xmax = cu + r; v2.ymin = cv - r; v2.ymax = cv + r;
		sprintf(label, "the crossing in (u, v), %.1e across", fmax(chi[0].x - clo[0].x, chi[0].y - clo[0].y));
		ssd_panel(&s, &v2, label, "step boxes of blue and red; black: the box holding their one crossing");
		na = s27_steps_2d(&s, &v2, e[0], hex[0], &ta, &tb);
		nb = s27_steps_2d(&s, &v2, e[1], hex[1], &sa, &sb);
		s27_float_2d(&s, &v2, f[0], G[0], ta, tb, hex[0]);
		s27_float_2d(&s, &v2, f[0], G[1], sa, sb, hex[1]);
		s27_box_2d(&s, &v2, clo[0].x, clo[0].y, chi[0].x, chi[0].y, "#1b1f24", 2.2);
		sprintf(label, "%u and %u step boxes here; transversal, each passing the other", na, nb);
		svg_text(&s, v2.x0 + 10, v2.y0 + v2.h - 12, 11, "#57606a", "start", label);
	}
	/* the end of blue: its box and the float geodesic's end */
	if (e[0])
	{
		viewport v2;
		qaws_vec2 lo, hi, q, p, d;
		double cu, cv, r;
		p.x = (qaws_scalar)G[0][0]; p.y = (qaws_scalar)G[0][1]; d.x = (qaws_scalar)G[0][2]; d.y = (qaws_scalar)G[0][3];
		qaws_geodesic_enclosure_get_end(&lo, &hi, NULL, NULL, e[0]);
		qaws_surface_exp(&q, NULL, f[0], p, d, (qaws_scalar)G[0][4], (qaws_scalar)1e-12);
		cu = 0.5 * (lo.x + hi.x); cv = 0.5 * (lo.y + hi.y); r = 1.5 * fmax(hi.x - lo.x, hi.y - lo.y);
		v2.x0 = 24 + 2 * (pw + 12); v2.y0 = 104; v2.w = pw; v2.h = ph;
		v2.xmin = cu - r; v2.xmax = cu + r; v2.ymin = cv - r; v2.ymax = cv + r;
		sprintf(label, "blue's end box, %.1e x %.1e", (double)(hi.x - lo.x), (double)(hi.y - lo.y));
		ssd_panel(&s, &v2, label, "proven to hold the end; the dot: the float geodesic's end (tolerance 1e-12)");
		s27_box_2d(&s, &v2, lo.x, lo.y, hi.x, hi.y, hex[0], 2.0);
		svg_circle(&s, vx(&v2, q.x), vy(&v2, q.y), 4, "#1b1f24", "#ffffff");
	}
	svg_close(&s);
	for (k = 0; k < 3; k++) qaws_geodesic_enclosure_destroy(e[k]);
	qaws_surface_destroy(f[0]);
	printf("figure 27: validated geodesics, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
