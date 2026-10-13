/* Figure 26: shortest paths round holes (qaws_shell_curve_create_shortest_path)
   and rays against thin holes. A plate with four holes: paths from one
   point to a row of points beyond them, each the geodesics touching the
   holes joined along their rims (a graph of them searched by Dijkstra). A
   cylinder with two holes: the same on a curved face. A plate with holes
   far smaller than a trace step (0.008 to 0.12 across against 0.2): rays
   through them stop at their rims.
   Drawn with figure 24's shaded renderer. */

static double const s26_hx[4] = { -0.5, 0.2, 0.5, -0.1 }, s26_hy[4] = { 0.15, -0.45, 0.4, 0.85 }, s26_hr[4] = { 0.3, 0.25, 0.2, 0.18 };

/* the cylinder's two holes in (u, v) */
static double const s26_cu[2] = { 0.5, 0.5 }, s26_cv[2] = { 0.38, 0.6 }, s26_cr = 0.06;

static int s26_keep_cylinder(unsigned int face, double u, double v)
{
	unsigned int i;
	(void)face;
	for (i = 0; i < 2; i++)
		if (hypot(u - s26_cu[i], v - s26_cv[i]) < s26_cr)
			return 0;
	return 1;
}

static void demo_shell_paths(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 420, ph = 420, pi = 3.14159265358979323846;
	static double const P00[3] = { -2, -2, 0 }, P10[3] = { 2, -2, 0 }, P01[3] = { -2, 2, 0 }, P11[3] = { 2, 2, 0 };
	static int const red[3] = { 207, 34, 46 }, black[3] = { 27, 31, 36 };
	unsigned int i;
	char label[220];
	if (!svg_open(&s, "showcase/b2d26_shell_paths.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + (ph + 44)),
		"Shortest paths round holes",
		"qaws_shell_curve_create_shortest_path: geodesics touching the holes joined along their rims; thin holes stop the rays through them"))
		return;
	/* a plate with four holes */
	{
		static double const box[1][4] = { { 0.125, 0.875, 0.125, 0.875 } };
		static double const to[7][2] = { { 1.3, -1.2 }, { 1.3, -0.8 }, { 1.3, -0.4 }, { 1.3, 0.0 }, { 1.3, 0.4 }, { 1.3, 0.8 }, { 1.3, 1.2 } };
		qaws_surface* f[1];
		qaws_curve *oc[4], *hc[4];
		qaws_curve const *ol[4], *hl[4];
		qaws_path_2d pa[5];
		qaws_shell_face_desc fd[1];
		qaws_surface_region* tr[1] = { NULL };
		qaws_shell_atlas* a = NULL;
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		double tu[7][2], total = 0;
		f[0] = shf_face(P00, P10, P01, P11);
		oc[0] = s25_line(0.125, 0.125, 0.875, 0.125); oc[1] = s25_line(0.875, 0.125, 0.875, 0.875);
		oc[2] = s25_line(0.875, 0.875, 0.125, 0.875); oc[3] = s25_line(0.125, 0.875, 0.125, 0.125);
		for (i = 0; i < 4; i++)
		{
			ol[i] = oc[i];
			hc[i] = sgf_arc((s26_hx[i] + 2) / 4, (s26_hy[i] + 2) / 4, s26_hr[i] / 4, 0, 2 * pi);
			hl[i] = hc[i];
			pa[1 + i].curves = &hl[i]; pa[1 + i].curve_count = 1; pa[1 + i].closed = 1;
		}
		pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
		fd[0].surface = f[0]; fd[0].trims = pa; fd[0].trim_count = 5;
		qaws_shell_atlas_create_trimmed(fd, 1, 0, &a);
		tr[0] = s25_trim_region(f[0], pa, 5);
		for (i = 0; i < 7; i++) { tu[i][0] = (to[i][0] + 2) / 4; tu[i][1] = (to[i][1] + 2) / 4; }
		for (i = 0; i < 7; i++)
		{
			qaws_shell_curve* g = NULL;
			qaws_scalar L = 0;
			if (qaws_shell_curve_create_shortest_path(a, shf_pt(0, (-1.3 + 2) / 4, 0.5), shf_pt(0, tu[i][0], tu[i][1]), 0, &g) == QAWS_STATUS_OK)
				qaws_shell_curve_get_length(g, &L);
			total += (double)L;
			qaws_shell_curve_destroy(g);
		}
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24; v3.y0 = 104; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -0.6, 0.85, 92, v3.x0 + pw / 2, v3.y0 + ph / 2 + 20);
		sc.view = &view; sc.trims = (qaws_surface_region const* const*)tr; sc.box = box; sc.res = 90;
		sprintf(label, "a plate with four holes: seven paths, %.6f in all", total);
		ssd_panel(&s, &v3, label, "each: geodesics touching the holes, joined along their rims");
		ssd_surfaces(&sc, f, 1, NULL, 0);
		s25_paths(&sc, a, f, shf_pt(0, (-1.3 + 2) / 4, 0.5), (double const (*)[2])tu, 7, red);
		s25_trim_edges(&sc, f[0], pa, 5);
		ssd_dot(&sc, f[0], (-1.3 + 2) / 4, 0.5, black, 3.5);
		ssd_draw(&s, &sc);
		qaws_surface_region_destroy(tr[0]);
		qaws_shell_atlas_destroy(a);
		qaws_surface_destroy(f[0]);
		for (i = 0; i < 4; i++) { qaws_curve_destroy(oc[i]); qaws_curve_destroy(hc[i]); }
	}
	/* a cylinder trimmed by its rims, with two holes */
	{
		static double const box[1][4] = { { 0, 1, 0, 1 } };
		static double const to[7][2] = { { 0.72, 0.28 }, { 0.72, 0.36 }, { 0.72, 0.43 }, { 0.72, 0.49 }, { 0.72, 0.55 }, { 0.72, 0.62 }, { 0.72, 0.72 } };
		static int const in_[1] = { -1 };
		qaws_curve *prof = s25_line(0.5, 0, 0.5, -pi), *bot = s25_line(0, 0, 1, 0), *top = s25_line(1, 1, 0, 1), *hc[2];
		qaws_curve const *l1[1], *l2[1], *lh[2];
		qaws_surface* f[1];
		qaws_path_2d pa[4];
		qaws_shell_face_desc fd[1];
		qaws_shell_atlas* a = NULL;
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		f[0] = sgf_revolve(prof);
		l1[0] = bot; l2[0] = top;
		pa[0].curves = l1; pa[0].curve_count = 1; pa[0].closed = 1;
		pa[1].curves = l2; pa[1].curve_count = 1; pa[1].closed = 1;
		for (i = 0; i < 2; i++)
		{
			hc[i] = sgf_arc(s26_cu[i], s26_cv[i], s26_cr, 0, 2 * pi);
			lh[i] = hc[i];
			pa[2 + i].curves = &lh[i]; pa[2 + i].curve_count = 1; pa[2 + i].closed = 1;
		}
		fd[0].surface = f[0]; fd[0].trims = pa; fd[0].trim_count = 4;
		qaws_shell_atlas_create_trimmed(fd, 1, 0, &a);
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24 + (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, 0.5 * pi + 0.25, 0.25, 112, v3.x0 + pw / 2, v3.y0 + ph / 2 + 15);
		view.o[2] = -0.5 * pi;
		sc.view = &view; sc.keep = s26_keep_cylinder; sc.box = box; sc.res = 90; sc.inner = in_;
		ssd_panel(&s, &v3, "a cylinder with two holes", "the paths wind between and round them, on the curved face");
		ssd_surfaces(&sc, f, 1, NULL, 0);
		s25_paths(&sc, a, f, shf_pt(0, 0.28, 0.5), to, 7, red);
		s25_trim_edges(&sc, f[0], pa + 2, 2);
		ssd_dot(&sc, f[0], 0.28, 0.5, black, 3.5);
		ssd_draw(&s, &sc);
		qaws_shell_atlas_destroy(a);
		qaws_surface_destroy(f[0]);
		qaws_curve_destroy(prof); qaws_curve_destroy(bot); qaws_curve_destroy(top);
		for (i = 0; i < 2; i++) qaws_curve_destroy(hc[i]);
	}
	/* holes far smaller than a trace step */
	{
		static double const box[1][4] = { { 0.125, 0.875, 0.125, 0.875 } };
		static double const hy[5] = { -0.8, -0.4, 0.0, 0.4, 0.8 }, hr[5] = { 0.004, 0.008, 0.016, 0.03, 0.06 };
		qaws_surface* f[1];
		qaws_curve *oc[4], *hc[5];
		qaws_curve const *ol[4], *hl[5];
		qaws_path_2d pa[6];
		qaws_shell_face_desc fd[1];
		qaws_surface_region* tr[1] = { NULL };
		qaws_shell_atlas* a = NULL;
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		double dn[9][2], cx = -1.3, cy = 0.0, hx = 0.6;
		f[0] = shf_face(P00, P10, P01, P11);
		oc[0] = s25_line(0.125, 0.125, 0.875, 0.125); oc[1] = s25_line(0.875, 0.125, 0.875, 0.875);
		oc[2] = s25_line(0.875, 0.875, 0.125, 0.875); oc[3] = s25_line(0.125, 0.875, 0.125, 0.125);
		for (i = 0; i < 4; i++) ol[i] = oc[i];
		pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
		for (i = 0; i < 5; i++)
		{
			hc[i] = sgf_arc((hx + 2) / 4, (hy[i] + 2) / 4, hr[i] / 4, 0, 2 * pi);
			hl[i] = hc[i];
			pa[1 + i].curves = &hl[i]; pa[1 + i].curve_count = 1; pa[1 + i].closed = 1;
		}
		fd[0].surface = f[0]; fd[0].trims = pa; fd[0].trim_count = 6;
		qaws_shell_atlas_create_trimmed(fd, 1, 0, &a);
		tr[0] = s25_trim_region(f[0], pa, 6);
		/* rays at each hole's centre, and between them */
		for (i = 0; i < 9; i++)
		{
			double ty = i % 2 == 0 ? hy[i / 2] : 0.5 * (hy[i / 2] + hy[i / 2 + 1]), nn;
			dn[i][0] = hx - cx; dn[i][1] = ty - cy;
			nn = hypot(dn[i][0], dn[i][1]);
			dn[i][0] /= nn; dn[i][1] /= nn;
		}
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24 + 2 * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -0.6, 0.85, 92, v3.x0 + pw / 2, v3.y0 + ph / 2 + 20);
		sc.view = &view; sc.trims = (qaws_surface_region const* const*)tr; sc.box = box; sc.res = 90;
		ssd_panel(&s, &v3, "thin holes, 0.008 to 0.12 across (a trace step 0.2)", "rays at them stop at their rims; between them they pass");
		ssd_surfaces(&sc, f, 1, NULL, 0);
		s25_rays(&sc, a, f, shf_pt(0, (cx + 2) / 4, (cy + 2) / 4), (double const (*)[2])dn, 9, 2.6, red);
		s25_trim_edges(&sc, f[0], pa, 6);
		for (i = 0; i < 5; i++) ssd_dot(&sc, f[0], (hx + 2) / 4, (hy[i] + 2) / 4, black, 1.2);
		ssd_dot(&sc, f[0], (cx + 2) / 4, (cy + 2) / 4, black, 3.5);
		ssd_draw(&s, &sc);
		qaws_surface_region_destroy(tr[0]);
		qaws_shell_atlas_destroy(a);
		qaws_surface_destroy(f[0]);
		for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
		for (i = 0; i < 5; i++) qaws_curve_destroy(hc[i]);
	}
	svg_close(&s);
	printf("figure 26: shortest paths round holes, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
