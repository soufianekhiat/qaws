/* Figure 28: geodesic queries (qaws_surface_curve_closest_point,
   qaws_surface_curve_intersect, qaws_surface_curve_closest_points,
   qaws_surface_region_distance) on figure 27's bump: points' closest
   points on a geodesic, joined by the shortest geodesics (square to it);
   two geodesics crossing and two apart with their closest pair; points'
   signed distances to a geodesic disk. */

/* the geodesic from a to b drawn */
static void s28_join(ssd_scene* sc, qaws_surface* f, qaws_vec2 a, qaws_vec2 b, int const* rgb, double width)
{
	qaws_surface_curve* g = NULL;
	if (fabs(a.x - b.x) + fabs(a.y - b.y) > 1e-9 && qaws_surface_curve_create_geodesic(&g, f, a, b, 0) == QAWS_STATUS_OK)
		ssd_curve_line(sc, f, qaws_surface_curve_get_uv(g), rgb, width);
	qaws_surface_curve_destroy(g);
}

static void s28_scene(ssd_scene* sc, sgf_view* view, viewport* v3, unsigned int k, double pw, double ph)
{
	static double const box[1][4] = { { 0, 1, 0, 1 } };
	memset(sc, 0, sizeof(*sc));
	v3->x0 = 24 + k * (pw + 12); v3->y0 = 104; v3->w = pw; v3->h = ph;
	sgf_view_init(view, -0.7, 0.75, 150, v3->x0 + pw / 2, v3->y0 + ph / 2 + 10);
	sc->view = view; sc->box = box; sc->res = 70;
}

static void demo_geodesic_query(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 420, ph = 420, pi = 3.14159265358979323846;
	static int const blue[3] = { 9, 105, 218 }, red[3] = { 207, 34, 46 }, green[3] = { 26, 127, 55 }, black[3] = { 27, 31, 36 }, purple[3] = { 130, 80, 223 };
	static int const fblue[3] = { 84, 174, 255 };
	qaws_surface* f[1];
	viewport v3;
	sgf_view view;
	ssd_scene sc;
	char label[240];
	unsigned int i;
	if (!svg_open(&s, "showcase/b2d28_geodesic_query.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + (ph + 44)),
		"Geodesic queries",
		"distances, closest points and crossings along the surface: qaws_surface_distance, _curve_closest_point, _curve_intersect, _curve_closest_points, _region_distance"))
		return;
	f[0] = s27_bump();
	/* points to a geodesic */
	{
		static double const pts[4][2] = { { 0.3, 0.75 }, { 0.55, 0.85 }, { 0.75, 0.3 }, { 0.45, 0.15 } };
		qaws_surface_curve* g = NULL;
		double worst = 0;
		qaws_surface_curve_create_geodesic(&g, f[0], sgf_uv(0.1, 0.25), sgf_uv(0.9, 0.7), 0);
		s28_scene(&sc, &view, &v3, 0, pw, ph);
		ssd_surfaces(&sc, f, 1, NULL, 0);
		if (g) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(g), blue, 1.8);
		for (i = 0; g && i < 4; i++)
		{
			qaws_vec2 q, p = sgf_uv(pts[i][0], pts[i][1]);
			qaws_scalar d = 0;
			if (qaws_surface_curve_closest_point(NULL, &q, &d, g, p, 0) != QAWS_STATUS_OK)
				continue;
			s28_join(&sc, f[0], p, q, red, 1.3);
			ssd_dot(&sc, f[0], pts[i][0], pts[i][1], black, 3.2);
			ssd_dot(&sc, f[0], q.x, q.y, red, 2.6);
			if ((double)d > worst) worst = (double)d;
		}
		sprintf(label, "points to a geodesic: the farthest %.4f away", worst);
		ssd_panel(&s, &v3, label, "red: the shortest geodesic to the curve, square to it where it lands");
		ssd_draw(&s, &sc);
		qaws_surface_curve_destroy(g);
	}
	/* two geodesics crossing; two apart and their closest pair */
	{
		qaws_surface_curve *a = NULL, *b = NULL, *c = NULL, *e = NULL;
		qaws_surface_crossing x[2];
		unsigned int nx = 0;
		qaws_vec2 pa, pb;
		qaws_scalar d = 0;
		qaws_surface_curve_create_geodesic(&a, f[0], sgf_uv(0.15, 0.3), sgf_uv(0.8, 0.65), 0);
		qaws_surface_curve_create_geodesic(&b, f[0], sgf_uv(0.2, 0.7), sgf_uv(0.8, 0.3), 0);
		qaws_surface_curve_create_geodesic(&c, f[0], sgf_uv(0.05, 0.92), sgf_uv(0.45, 0.95), 0);
		qaws_surface_curve_create_geodesic(&e, f[0], sgf_uv(0.62, 0.9), sgf_uv(0.95, 0.5), 0);
		if (a && b) qaws_surface_curve_intersect(x, &nx, a, b, 2);
		if (c && e) qaws_surface_curve_closest_points(NULL, NULL, &pa, &pb, &d, c, e, 0);
		s28_scene(&sc, &view, &v3, 1, pw, ph);
		ssd_surfaces(&sc, f, 1, NULL, 0);
		if (a) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(a), blue, 1.8);
		if (b) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(b), red, 1.8);
		if (c) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(c), green, 1.8);
		if (e) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(e), purple, 1.8);
		for (i = 0; i < nx && i < 2; i++) ssd_dot(&sc, f[0], x[i].point.x, x[i].point.y, black, 4);
		s28_join(&sc, f[0], pa, pb, black, 1.2);
		ssd_dot(&sc, f[0], pa.x, pa.y, green, 3); ssd_dot(&sc, f[0], pb.x, pb.y, purple, 3);
		sprintf(label, "crossings and closest points: %u crossing; apart by %.4f", nx, (double)d);
		ssd_panel(&s, &v3, label, "blue and red cross (black); green and purple apart, joined at their closest pair");
		ssd_draw(&s, &sc);
		qaws_surface_curve_destroy(a); qaws_surface_curve_destroy(b); qaws_surface_curve_destroy(c); qaws_surface_curve_destroy(e);
	}
	/* signed distances to a geodesic disk */
	{
		static double const pts[5][2] = { { 0.2, 0.8 }, { 0.85, 0.75 }, { 0.8, 0.15 }, { 0.5, 0.5 }, { 0.2, 0.3 } };
		qaws_surface_curve* cc = NULL;
		qaws_surface_curve const* pc[1];
		qaws_surface_path pp;
		qaws_surface_region* rg = NULL;
		qaws_surface_curve_create_geodesic_arc(&cc, f[0], sgf_uv(0.45, 0.45), (qaws_scalar)0.55, 0, (qaws_scalar)(2 * pi), 0);
		pc[0] = cc;
		pp.curves = pc; pp.curve_count = 1; pp.closed = 1;
		qaws_surface_clip(&rg, f[0], &pp, NULL, 1, 0, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO);
		s28_scene(&sc, &view, &v3, 2, pw, ph);
		ssd_surfaces(&sc, f, 1, NULL, 0);
		if (cc) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(cc), blue, 1.8);
		label[0] = 0;
		for (i = 0; rg && i < 5; i++)
		{
			qaws_vec2 q, p = sgf_uv(pts[i][0], pts[i][1]);
			qaws_scalar d = 0;
			char one[40];
			if (qaws_surface_region_distance(&d, &q, rg, p, 0) != QAWS_STATUS_OK)
				continue;
			s28_join(&sc, f[0], p, q, d < 0 ? green : red, 1.3);
			ssd_dot(&sc, f[0], pts[i][0], pts[i][1], d < 0 ? green : red, 3.2);
			sprintf(one, "%s%+.3f", i ? ", " : "", (double)d);
			strcat(label, one);
		}
		{
			char title[300];
			sprintf(title, "signed distances to a disk: %s", label);
			ssd_panel(&s, &v3, title, "red outside, green inside (negative), each to its nearest boundary point");
		}
		(void)fblue;
		ssd_draw(&s, &sc);
		qaws_surface_region_destroy(rg);
		qaws_surface_curve_destroy(cc);
	}
	/* shortest paths round two disks */
	{
		static double const to[5][2] = { { 0.85, 0.2 }, { 0.88, 0.4 }, { 0.88, 0.6 }, { 0.85, 0.8 }, { 0.7, 0.9 } };
		qaws_surface_curve *c1 = NULL, *c2 = NULL;
		qaws_surface_curve const *p1[1], *p2[1];
		qaws_surface_path pp[2];
		qaws_surface_region* rg = NULL;
		double total = 0;
		unsigned int ok = 0;
		qaws_surface_curve_create_geodesic_arc(&c1, f[0], sgf_uv(0.45, 0.35), (qaws_scalar)0.25, 0, (qaws_scalar)(2 * pi), 0);
		qaws_surface_curve_create_geodesic_arc(&c2, f[0], sgf_uv(0.55, 0.7), (qaws_scalar)0.2, 0, (qaws_scalar)(2 * pi), 0);
		p1[0] = c1; p2[0] = c2;
		pp[0].curves = p1; pp[0].curve_count = 1; pp[0].closed = 1;
		pp[1].curves = p2; pp[1].curve_count = 1; pp[1].closed = 1;
		qaws_surface_clip(&rg, f[0], &pp[0], &pp[1], 1, 1, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
		s28_scene(&sc, &view, &v3, 3, pw, ph);
		ssd_surfaces(&sc, f, 1, NULL, 0);
		if (c1) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(c1), black, 1.4);
		if (c2) ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(c2), black, 1.4);
		for (i = 0; rg && i < 5; i++)
		{
			qaws_surface_curve* path = NULL;
			qaws_scalar L = 0;
			if (qaws_surface_curve_create_shortest_path(&path, f[0], rg, sgf_uv(0.12, 0.5), sgf_uv(to[i][0], to[i][1]), 0) == QAWS_STATUS_OK)
			{
				ssd_curve_line(&sc, f[0], qaws_surface_curve_get_uv(path), red, 1.5);
				qaws_surface_curve_get_length(&L, path);
				total += (double)L;
				ok++;
			}
			qaws_surface_curve_destroy(path);
			ssd_dot(&sc, f[0], to[i][0], to[i][1], black, 2.6);
		}
		ssd_dot(&sc, f[0], 0.12, 0.5, black, 3.5);
		sprintf(label, "shortest paths round two disks: %u paths, %.4f in all", ok, total);
		ssd_panel(&s, &v3, label, "qaws_surface_curve_create_shortest_path: the disks' union as the obstacle");
		ssd_draw(&s, &sc);
		qaws_surface_region_destroy(rg);
		qaws_surface_curve_destroy(c1); qaws_surface_curve_destroy(c2);
	}
	svg_close(&s);
	qaws_surface_destroy(f[0]);
	printf("figure 28: geodesic queries, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
