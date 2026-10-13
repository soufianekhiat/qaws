/* Figure 24: geodesic disks through saddle vertices (qaws_shell_region.h).
   Faces round a vertex with more than 2 pi of angle: five squares (5 pi/2)
   and six (the monkey saddle, 3 pi), flat, then as bicubic patches bent by
   a smooth map (curved faces and curved edges). Rays from a centre part at
   the vertex; a disk through it is the plane disk and, behind the vertex,
   the disk about it of the radius left. On the curved monkey saddle, two
   disks through the vertex combined by every Boolean. Uses figure 23's
   helpers. */

/* faces round the origin, face i spanned by axis directions d[i] and
   d[i + 1] (scaled by size): flat bilinear squares (bend 0), or bicubic
   patches whose control grid is the flat one bent (shared edges share
   their control points, so the faces still meet along them) */
static void ssd_bend(double const* p, double k, qaws_vec3* q)
{
	q->x = (qaws_scalar)(p[0] + k * 0.35 * sin(2.2 * p[1] + 0.7 * p[2]));
	q->y = (qaws_scalar)(p[1] + k * 0.35 * sin(2.2 * p[2] + 0.7 * p[0]));
	q->z = (qaws_scalar)(p[2] + k * 0.45 * sin(2.0 * p[0] + 1.3 * p[1]) * cos(1.1 * p[2]));
}

static void ssd_fan(qaws_surface** f, double const (*d)[3], unsigned int n, double size, double bend)
{
	static double const O[3] = { 0, 0, 0 };
	unsigned int i, j, q, m;
	for (i = 0; i < n; i++)
	{
		double const* a = d[i];
		double const* b = d[(i + 1) % n];
		if (bend == 0)
		{
			double pa[3], pb[3], pc[3];
			for (m = 0; m < 3; m++) { pa[m] = size * a[m]; pb[m] = size * b[m]; pc[m] = pa[m] + pb[m]; }
			f[i] = shf_face(O, pa, pb, pc);
		}
		else
		{
			qaws_vec3 cp[16];
			qaws_surface_bezier_desc bd;
			for (j = 0; j < 4; j++)
				for (q = 0; q < 4; q++)
				{
					double p[3];
					for (m = 0; m < 3; m++) p[m] = size * (a[m] * j / 3.0 + b[m] * q / 3.0);
					ssd_bend(p, bend, &cp[4 * j + q]);
				}
			memset(&bd, 0, sizeof(bd));
			bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
			f[i] = NULL;
			qaws_surface_create_bezier(&bd, &f[i]);
		}
	}
}

static void ssd_faces(svg* s, qaws_surface** f, unsigned int n, sgf_view const* view, int curved)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		if (curved) sgf_mesh(s, f[i], view, 8, 0);
		shf_face_box(s, f[i], view);
	}
}

static qaws_shell_region* ssd_disk(qaws_shell_atlas const* a, qaws_shell_point c, double r, double* area)
{
	qaws_shell_region* g = NULL;
	qaws_scalar A = 0;
	qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)r, 0, &g);
	if (g) qaws_shell_region_get_area(g, &A);
	if (area) *area = (double)A;
	return g;
}

static void ssd_fill(svg* s, qaws_shell_region const* g, qaws_shell_atlas const* a, sgf_view const* view, char const* fill)
{
	if (g && fill) shf_region(s, g, a, view, fill);
}

static void ssd_outline(svg* s, qaws_shell_region const* g, qaws_surface** f, unsigned int n, sgf_view const* view, char const* line)
{
	unsigned int k;
	if (!g)
		return;
	for (k = 0; k < n; k++)
		if (qaws_shell_region_get_face(g, k))
			sof_outline(s, qaws_shell_region_get_face(g, k), f[k], view, NULL, 0, 0, line, 1.4);
}

static void ssd_rays(svg* s, qaws_shell_atlas const* a, qaws_shell_point c, double len, sgf_view const* view)
{
	unsigned int i;
	double const pi = 3.14159265358979323846;
	for (i = 0; i < 30; i++)
	{
		qaws_shell_curve* g = NULL;
		qaws_vec2 d;
		double th = 1.25 * pi + 0.6 * pi * ((i + 0.5) / 30.0 - 0.5);
		d.x = (qaws_scalar)cos(th); d.y = (qaws_scalar)sin(th);
		if (qaws_shell_curve_create_geodesic_ray(a, c, d, (qaws_scalar)len, 0, &g) == QAWS_STATUS_OK)
			shf_curve(s, g, view, a, 0, i < 15 ? "#cf222e" : "#8250df");
		qaws_shell_curve_destroy(g);
	}
}

static void demo_shell_saddle(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 300, pi = 3.14159265358979323846, rho0 = 0.3 * sqrt(2.0);
	static double const D5[5][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 } };
	static double const D6[6][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 } };
	qaws_shell_point c = shf_pt(0, 0.3, 0.3);
	unsigned int k, i;
	char label[220];
	if (!svg_open(&s, "showcase/b2d24_shell_saddle.svg", (int)(48 + 6 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Geodesic disks through saddle vertices, flat and curved",
		"qaws_shell_region_create_geodesic_disk: faces round a vertex with more than 2 pi of angle; rays part there; bottom: two disks on a curved monkey saddle under every Boolean"))
		return;
	/* the top row: flat saddles, then bent ones */
	for (k = 0; k < 6; k++)
	{
		int six = k == 1 || k == 2 || k == 4 || k == 5, curved = k >= 3;
		unsigned int n = six ? 6 : 5;
		double A = 0;
		qaws_surface* f[6];
		qaws_shell_atlas* a = NULL;
		qaws_shell_point cc = curved ? shf_pt(0, 0.2, 0.2) : c;
		viewport v3;
		sgf_view view;
		ssd_fan(f, six ? D6 : D5, n, curved ? 1.6 : 1, curved ? 1 : 0);
		qaws_shell_atlas_create((qaws_surface const* const*)f, n, 0, &a);
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		if (six)
			sgf_view_init(&view, -0.25 * pi, 0.615, curved ? 68 : 95, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		else
			sgf_view_init(&view, 0.25 * pi, 0.5, curved ? 72 : 105, v3.x0 + pw / 2, v3.y0 + ph / 2 - 10);
		if (k == 0)
		{
			/* growing disks: short of the vertex, through it, well past it */
			static double const R[3] = { 0.75, 0.55, 0.35 };
			static char const* const F[3] = { "#b6e3ff", "#80ccff", "#54aeff" };
			svg_panel(&s, &v3, "five squares (5 pi/2): disks r 0.35, 0.55, 0.75");
			ssd_faces(&s, f, n, &view, 0);
			for (i = 0; i < 3; i++)
			{
				qaws_shell_region* g = ssd_disk(a, c, R[i], NULL);
				ssd_fill(&s, g, a, &view, F[i]);
				ssd_outline(&s, g, f, n, &view, "#0969da");
				qaws_shell_region_destroy(g);
			}
		}
		else if (k == 1 || k == 4)
		{
			svg_panel(&s, &v3, curved ? "curved monkey saddle: rays part at the vertex" : "monkey saddle (3 pi): rays part at the vertex");
			ssd_faces(&s, f, n, &view, curved);
			ssd_rays(&s, a, cc, curved ? 1.5 : 1.3, &view);
		}
		else
		{
			qaws_shell_region* g = ssd_disk(a, cc, 0.6, &A);
			if (k == 2)
				sprintf(label, "monkey saddle: %.6f, closed form %.6f", A, pi * 0.36 + 0.5 * (n * pi / 2 - 2 * pi) * (0.6 - rho0) * (0.6 - rho0));
			else
				sprintf(label, "curved %s saddle: disk r 0.6, area %.6f", six ? "monkey" : "five-face", A);
			svg_panel(&s, &v3, label);
			ssd_faces(&s, f, n, &view, curved);
			ssd_fill(&s, g, a, &view, "#54aeff");
			ssd_outline(&s, g, f, n, &view, "#0969da");
			qaws_shell_region_destroy(g);
		}
		qaws_shell_atlas_destroy(a);
		for (i = 0; i < n; i++) qaws_surface_destroy(f[i]);
	}
	/* the bottom row: two disks through the vertex of the curved monkey
	   saddle, A (blue) and B (red), and their Booleans */
	{
		static qaws_clip_type const T[4] = { QAWS_CLIP_UNION, QAWS_CLIP_INTERSECTION, QAWS_CLIP_DIFFERENCE, QAWS_CLIP_XOR };
		static char const* const TN[4] = { "A union B", "A intersect B", "A minus B", "A xor B" };
		static char const* const TF[4] = { "#c297ff", "#8250df", "#54aeff", "#ffb3b8" };
		qaws_surface* f[6];
		qaws_shell_atlas* a = NULL;
		qaws_shell_region *ga, *gb;
		double A = 0, B = 0;
		ssd_fan(f, D6, 6, 1.6, 1);
		qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &a);
		ga = ssd_disk(a, shf_pt(0, 0.2, 0.2), 0.6, &A);
		gb = ssd_disk(a, shf_pt(2, 0.22, 0.19), 0.55, &B);
		for (k = 0; k < 6; k++)
		{
			viewport v3;
			sgf_view view;
			qaws_shell_region* x = NULL;
			qaws_scalar X = 0;
			v3.x0 = 24 + k * (pw + 12); v3.y0 = 104 + ph + 44; v3.w = pw; v3.h = ph;
			sgf_view_init(&view, -0.25 * pi, 0.615, 68, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
			if (k >= 1 && k <= 4 && ga && gb)
			{
				qaws_shell_region_clip(ga, gb, T[k - 1], &x);
				if (x) qaws_shell_region_get_area(x, &X);
				sprintf(label, "curved monkey: %s, area %.6f", TN[k - 1], (double)X);
			}
			else if (k == 0)
				sprintf(label, "curved monkey: disks A %.6f, B %.6f", A, B);
			else
			{
				/* B minus A, from the other side */
				if (ga && gb) qaws_shell_region_clip(gb, ga, QAWS_CLIP_DIFFERENCE, &x);
				if (x) qaws_shell_region_get_area(x, &X);
				sprintf(label, "curved monkey: B minus A, area %.6f", (double)X);
			}
			svg_panel(&s, &v3, label);
			ssd_faces(&s, f, 6, &view, 1);
			if (k == 0)
			{
				ssd_fill(&s, ga, a, &view, "#54aeff");
				ssd_fill(&s, gb, a, &view, "#ff8182");
			}
			else
				ssd_fill(&s, x, a, &view, k == 5 ? "#ff8182" : TF[k - 1]);
			ssd_outline(&s, ga, f, 6, &view, "#0969da");
			ssd_outline(&s, gb, f, 6, &view, "#cf222e");
			qaws_shell_region_destroy(x);
		}
		qaws_shell_region_destroy(ga);
		qaws_shell_region_destroy(gb);
		qaws_shell_atlas_destroy(a);
		for (i = 0; i < 6; i++) qaws_surface_destroy(f[i]);
	}
	svg_close(&s);
	printf("figure 24: geodesic disks through saddle vertices, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
