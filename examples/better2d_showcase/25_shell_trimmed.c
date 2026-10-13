/* Figure 25: shells of trimmed faces (qaws_shell_atlas_create_trimmed).
   Faces bounded by curves of their (u, v), meeting along those curves: a
   plane split by a circle into a square with a hole and a disk (each with
   its own surface); a plate with a round hole and a tube hanging from it
   (the hole's trim meeting the tube's top); a cube of larger planes
   trimmed to its faces. Straightest rays and geodesic disks run over the
   trim edges. Drawn with figure 24's shaded renderer. */

/* a face's trim loops as a region of it (to draw only its trimmed part) */
static qaws_surface_region* s25_trim_region(qaws_surface const* s, qaws_path_2d const* loops, unsigned int nl)
{
	qaws_surface_path P[4];
	qaws_surface_curve* sc[16];
	qaws_surface_curve const* cc[16];
	qaws_surface_region* r = NULL;
	unsigned int l, k, n = 0;
	for (l = 0; l < nl; l++)
	{
		P[l].curves = cc + n; P[l].curve_count = loops[l].curve_count; P[l].closed = 1;
		for (k = 0; k < loops[l].curve_count; k++, n++)
		{
			sc[n] = NULL;
			qaws_surface_curve_create_uv(s, loops[l].curves[k], &sc[n]);
			cc[n] = sc[n];
		}
	}
	qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_EVEN_ODD, P, nl, NULL, 0, &r);
	for (k = 0; k < n; k++) qaws_surface_curve_destroy(sc[k]);
	return r;
}

/* a face's trim loops, dark */
static void s25_trim_edges(ssd_scene* sc, qaws_surface const* f, qaws_path_2d const* loops, unsigned int nl)
{
	static int const ink[3] = { 87, 96, 106 };
	unsigned int l, k;
	for (l = 0; l < nl; l++)
		for (k = 0; k < loops[l].curve_count; k++)
			ssd_curve_line(sc, f, loops[l].curves[k], ink, 1.1);
}

static qaws_curve* s25_line(double x0, double y0, double x1, double y1)
{
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_curve* c = NULL;
	cp[0] = (qaws_scalar)x0; cp[1] = (qaws_scalar)y0; cp[2] = (qaws_scalar)x1; cp[3] = (qaws_scalar)y1;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	qaws_curve_create_bezier(&bd, &c);
	return c;
}

/* rays from c at the given (u, v) directions over len */
static void s25_rays(ssd_scene* sc, qaws_shell_atlas const* a, qaws_surface** f, qaws_shell_point c, double const (*dirs)[2], unsigned int n, double len, int const* rgb)
{
	unsigned int i, w;
	for (i = 0; i < n; i++)
	{
		qaws_shell_curve* g = NULL;
		qaws_vec2 d;
		d.x = (qaws_scalar)dirs[i][0]; d.y = (qaws_scalar)dirs[i][1];
		if (qaws_shell_curve_create_geodesic_ray(a, c, d, (qaws_scalar)len, 0, &g) == QAWS_STATUS_OK)
			for (w = 0; w < qaws_shell_curve_get_run_count(g); w++)
			{
				unsigned int fw;
				qaws_surface_curve const* run;
				qaws_shell_curve_get_run(g, w, &fw, &run);
				ssd_curve_line(sc, f[fw], qaws_surface_curve_get_uv(run), rgb, 1.4);
			}
		qaws_shell_curve_destroy(g);
	}
}

static void demo_shell_trimmed(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 420, ph = 420, pi = 3.14159265358979323846;
	static double const P00[3] = { -2, -2, 0 }, P10[3] = { 2, -2, 0 }, P01[3] = { -2, 2, 0 }, P11[3] = { 2, 2, 0 };
	static int const blue[3] = { 9, 105, 218 }, red[3] = { 207, 34, 46 }, black[3] = { 27, 31, 36 };
	static int const fblue[3] = { 84, 174, 255 }, fred[3] = { 255, 129, 130 }, fpurple[3] = { 170, 120, 250 };
	unsigned int k, i, j;
	char label[220];
	if (!svg_open(&s, "showcase/b2d25_shell_trimmed.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Shells of trimmed faces",
		"qaws_shell_atlas_create_trimmed: faces bounded by curves of their (u, v), edges matched along those curves in space; rays and disks cross them"))
		return;
	/* the plane split by a circle */
	{
		static double const Q00[3] = { 1, -1, 0 }, Q10[3] = { 1, 1, 0 }, Q01[3] = { -1, -1, 0 }, Q11[3] = { -1, 1, 0 };
		static double const box[2][4] = { { 0.125, 0.875, 0.125, 0.875 }, { 0.1, 0.9, 0.1, 0.9 } };
		qaws_surface* f[2];
		qaws_curve *oc[4], *hole = sgf_arc(0.5, 0.5, 0.2, 0, 2 * pi), *disk = sgf_arc(0.5, 0.5, 0.4, 1.0, 1.0 + 2 * pi);
		qaws_curve const *ol[4], *hl[1], *dl[1];
		qaws_path_2d pa[2], pb[1];
		qaws_shell_face_desc fd[2];
		qaws_surface_region* tr[2];
		qaws_shell_atlas* a = NULL;
		qaws_shell_region* g = NULL;
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		ssd_paint paint[1];
		double A = 0;
		f[0] = shf_face(P00, P10, P01, P11);
		f[1] = shf_face(Q00, Q10, Q01, Q11);
		oc[0] = s25_line(0.125, 0.125, 0.875, 0.125); oc[1] = s25_line(0.875, 0.125, 0.875, 0.875);
		oc[2] = s25_line(0.875, 0.875, 0.125, 0.875); oc[3] = s25_line(0.125, 0.875, 0.125, 0.125);
		for (i = 0; i < 4; i++) ol[i] = oc[i];
		hl[0] = hole; dl[0] = disk;
		pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
		pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
		pb[0].curves = dl; pb[0].curve_count = 1; pb[0].closed = 1;
		fd[0].surface = f[0]; fd[0].trims = pa; fd[0].trim_count = 2;
		fd[1].surface = f[1]; fd[1].trims = pb; fd[1].trim_count = 1;
		qaws_shell_atlas_create_trimmed(fd, 2, 0, &a);
		tr[0] = s25_trim_region(f[0], pa, 2);
		tr[1] = s25_trim_region(f[1], pb, 1);
		g = ssd_disk(a, shf_pt(0, (0.9 + 2) / 4, (0.2 + 2) / 4), 0.5, &A);
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24; v3.y0 = 104; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -0.6, 0.85, 120, v3.x0 + pw / 2, v3.y0 + ph / 2 + 15);
		sc.view = &view; sc.trims = (qaws_surface_region const* const*)tr; sc.box = box; sc.res = 90;
		sprintf(label, "a plane split by a circle: disk %.9f, pi r^2 %.9f", A, pi * 0.25);
		ssd_panel(&s, &v3, label, "two trimmed faces of two surfaces meet along the circle; rays cross it straight");
		paint[0].r = g; memcpy(paint[0].rgb, fblue, sizeof(fblue));
		ssd_surfaces(&sc, f, 2, paint, 1);
		ssd_outline(&sc, g, f, 2, blue);
		{
			static double const dirs[5][2] = { { 1, 0.3 }, { 1, 0.15 }, { 1, 0 }, { 1, -0.15 }, { 1, -0.3 } };
			s25_rays(&sc, a, f, shf_pt(0, (-1.3 + 2) / 4, (0.1 + 2) / 4), dirs, 5, 2.5, red);
		}
		s25_trim_edges(&sc, f[0], pa, 2);
		s25_trim_edges(&sc, f[1], pb, 1);
		ssd_dot(&sc, f[0], (0.9 + 2) / 4, (0.2 + 2) / 4, black, 3.5);
		ssd_draw(&s, &sc);
		qaws_shell_region_destroy(g);
		qaws_surface_region_destroy(tr[0]); qaws_surface_region_destroy(tr[1]);
		qaws_shell_atlas_destroy(a);
		qaws_surface_destroy(f[0]); qaws_surface_destroy(f[1]);
		for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
		qaws_curve_destroy(hole); qaws_curve_destroy(disk);
	}
	/* the plate with a hole and the tube */
	for (k = 1; k < 4; k++)
	{
		static double const box[2][4] = { { 0.125, 0.875, 0.125, 0.875 }, { 0, 1, 0, 1 } };
		qaws_surface* f[2];
		qaws_curve *oc[4], *hole = sgf_arc(0.5, 0.5, 0.125, 0, 2 * pi), *prof = s25_line(0.5, 0, 0.5, -1.5);
		qaws_curve const *ol[4], *hl[1];
		qaws_path_2d pa[2];
		qaws_shell_face_desc fd[2];
		qaws_surface_region* tr[2];
		qaws_shell_atlas* a = NULL;
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		ssd_paint paint[2];
		unsigned int np = 0;
		qaws_shell_region *g = NULL, *h = NULL, *x = NULL;
		double A = 0, B = 0;
		qaws_scalar X = 0;
		f[0] = shf_face(P00, P10, P01, P11);
		f[1] = sgf_revolve(prof);
		oc[0] = s25_line(0.125, 0.125, 0.875, 0.125); oc[1] = s25_line(0.875, 0.125, 0.875, 0.875);
		oc[2] = s25_line(0.875, 0.875, 0.125, 0.875); oc[3] = s25_line(0.125, 0.875, 0.125, 0.125);
		for (i = 0; i < 4; i++) ol[i] = oc[i];
		hl[0] = hole;
		pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
		pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
		fd[0].surface = f[0]; fd[0].trims = pa; fd[0].trim_count = 2;
		fd[1].surface = f[1]; fd[1].trims = NULL; fd[1].trim_count = 0;
		qaws_shell_atlas_create_trimmed(fd, 2, 0, &a);
		tr[0] = s25_trim_region(f[0], pa, 2);
		tr[1] = NULL;
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24 + (k % 3) * (pw + 12); v3.y0 = 104 + (k / 3) * (ph + 44); v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -0.6, k == 1 ? 0.75 : -0.42, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 + (k == 1 ? 40 : -30));
		sc.view = &view; sc.trims = (qaws_surface_region const* const*)tr; sc.box = box; sc.res = 80;
		if (k == 1)
		{
			static double const dirs[7][2] = { { -1, -0.5 }, { -1, -0.33 }, { -1, -0.17 }, { -1, 0 }, { -1, 0.12 }, { -1, 0.24 }, { -1, 0.36 } };
			double dn[7][2];
			for (i = 0; i < 7; i++) { double nn = hypot(dirs[i][0], dirs[i][1]); dn[i][0] = dirs[i][0] / nn; dn[i][1] = dirs[i][1] / nn; }
			ssd_panel(&s, &v3, "a plate with a hole, a tube hanging from it", "straightest rays over the rim: helices down the tube (seen through the hole)");
			ssd_surfaces(&sc, f, 2, NULL, 0);
			s25_rays(&sc, a, f, shf_pt(0, (1.2 + 2) / 4, (0.3 + 2) / 4), (double const (*)[2])dn, 7, 2.4, red);
			ssd_dot(&sc, f[0], (1.2 + 2) / 4, (0.3 + 2) / 4, black, 3.5);
		}
		else if (k == 2)
		{
			g = ssd_disk(a, shf_pt(0, (0.75 + 2) / 4, (0.1 + 2) / 4), 0.5, &A);
			sprintf(label, "a disk over the rim: area %.6f", A);
			ssd_panel(&s, &v3, label, "from below: the part on the plate and the part down the tube");
			paint[0].r = g; memcpy(paint[0].rgb, fblue, sizeof(fblue)); np = 1;
			ssd_surfaces(&sc, f, 2, paint, np);
			ssd_outline(&sc, g, f, 2, blue);
			ssd_dot(&sc, f[0], (0.75 + 2) / 4, (0.1 + 2) / 4, black, 3.5);
		}
		else
		{
			g = ssd_disk(a, shf_pt(0, (0.75 + 2) / 4, (0.1 + 2) / 4), 0.5, &A);
			h = ssd_disk(a, shf_pt(1, 0.06, 0.25), 0.4, &B);
			if (g && h) qaws_shell_region_clip(g, h, QAWS_CLIP_UNION, &x);
			if (x) qaws_shell_region_get_area(x, &X);
			sprintf(label, "united with a disk on the tube: area %.6f", (double)X);
			ssd_panel(&s, &v3, label, "A on the plate (blue), B inside the tube (red), their union filled");
			paint[0].r = x; memcpy(paint[0].rgb, fpurple, sizeof(fpurple)); np = 1;
			ssd_surfaces(&sc, f, 2, paint, np);
			ssd_outline(&sc, g, f, 2, blue);
			ssd_outline(&sc, h, f, 2, red);
			ssd_dot(&sc, f[0], (0.75 + 2) / 4, (0.1 + 2) / 4, black, 3.5);
			ssd_dot(&sc, f[1], 0.06, 0.25, red, 3.5);
		}
		s25_trim_edges(&sc, f[0], pa, 2);
		{
			/* the tube's bottom rim */
			static int const ink[3] = { 87, 96, 106 };
			double uv[2 * 81];
			for (j = 0; j <= 80; j++) { uv[2 * j] = j / 80.0; uv[2 * j + 1] = 1; }
			ssd_uv_line(&sc, f[1], uv, 81, ink, 1.1);
		}
		ssd_draw(&s, &sc);
		qaws_shell_region_destroy(g); qaws_shell_region_destroy(h); qaws_shell_region_destroy(x);
		qaws_surface_region_destroy(tr[0]);
		qaws_shell_atlas_destroy(a);
		qaws_surface_destroy(f[0]); qaws_surface_destroy(f[1]);
		for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
		qaws_curve_destroy(hole); qaws_curve_destroy(prof);
	}
	/* a cube of planes trimmed to its faces */
	for (k = 4; k < 6; k++)
	{
		static double const P[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { 1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } };
		static int const C[6][4] = { { 4, 5, 6, 7 }, { 0, 2, 1, 3 }, { 1, 3, 5, 7 }, { 0, 4, 2, 6 }, { 2, 6, 3, 7 }, { 0, 1, 4, 5 } };
		double const lo = 1.0 / 6, span = 2.0 / 3;
		double box[6][4];
		qaws_surface* f[6];
		qaws_curve* sq[6][4];
		qaws_curve const* sl[6][4];
		qaws_path_2d loops[6];
		qaws_shell_face_desc fd[6];
		qaws_surface_region* tr[6];
		qaws_shell_atlas* a = NULL;
		qaws_shell_region* g = NULL;
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		ssd_paint paint[1];
		double A = 0;
		for (i = 0; i < 6; i++)
		{
			double q[4][3], c[3];
			unsigned int m;
			for (m = 0; m < 3; m++) c[m] = 0.25 * (P[C[i][0]][m] + P[C[i][1]][m] + P[C[i][2]][m] + P[C[i][3]][m]);
			for (j = 0; j < 4; j++)
				for (m = 0; m < 3; m++) q[j][m] = c[m] + 1.5 * (P[C[i][j]][m] - c[m]);
			f[i] = shf_face(q[0], q[1], q[2], q[3]);
			sq[i][0] = s25_line(lo, lo, 1 - lo, lo); sq[i][1] = s25_line(1 - lo, lo, 1 - lo, 1 - lo);
			sq[i][2] = s25_line(1 - lo, 1 - lo, lo, 1 - lo); sq[i][3] = s25_line(lo, 1 - lo, lo, lo);
			for (j = 0; j < 4; j++) sl[i][j] = sq[i][j];
			loops[i].curves = sl[i]; loops[i].curve_count = 4; loops[i].closed = 1;
			fd[i].surface = f[i]; fd[i].trims = &loops[i]; fd[i].trim_count = 1;
			tr[i] = NULL;
			box[i][0] = lo; box[i][1] = 1 - lo; box[i][2] = lo; box[i][3] = 1 - lo;
		}
		qaws_shell_atlas_create_trimmed(fd, 6, 0, &a);
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24 + (k % 3) * (pw + 12); v3.y0 = 104 + (k / 3) * (ph + 44); v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -0.75 * pi, 0.6, 85, v3.x0 + pw / 2, v3.y0 + ph / 2 + 15);
		sc.view = &view; sc.trims = (qaws_surface_region const* const*)tr; sc.box = (double const (*)[4])box; sc.res = 40;
		if (k == 4)
		{
			qaws_shell_curve* c = NULL;
			qaws_shell_curve_create_geodesic_circle(a, shf_pt(0, lo + span * 0.8, lo + span * 0.8), (qaws_scalar)0.5, 0, &c);
			if (c) qaws_shell_region_create(a, (qaws_shell_curve const* const*)&c, 1, &g);
			if (g) { qaws_scalar Ar = 0; qaws_shell_region_get_area(g, &Ar); A = (double)Ar; }
			qaws_shell_curve_destroy(c);
			sprintf(label, "a cube of trimmed planes: circle %.9f", A);
			ssd_panel(&s, &v3, label, "over two trim edges by a vertex, pi r^2 = 0.785398163");
		}
		else
		{
			g = ssd_disk(a, shf_pt(0, lo + span * 0.85, lo + span * 0.85), 0.75, &A);
			sprintf(label, "a disk round a trimmed corner: area %.6f", A);
			ssd_panel(&s, &v3, label, "the same as on the untrimmed cube (test 108)");
		}
		paint[0].r = g; memcpy(paint[0].rgb, fblue, sizeof(fblue));
		ssd_surfaces(&sc, f, 6, paint, 1);
		ssd_outline(&sc, g, f, 6, blue);
		for (i = 0; i < 6; i++) s25_trim_edges(&sc, f[i], &loops[i], 1);
		ssd_draw(&s, &sc);
		qaws_shell_region_destroy(g);
		qaws_shell_atlas_destroy(a);
		for (i = 0; i < 6; i++)
		{
			qaws_surface_destroy(f[i]);
			for (j = 0; j < 4; j++) qaws_curve_destroy(sq[i][j]);
		}
	}
	svg_close(&s);
	printf("figure 25: shells of trimmed faces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
