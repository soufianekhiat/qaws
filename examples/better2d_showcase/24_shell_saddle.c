/* Figure 24: geodesic disks through saddle vertices (qaws_shell_region.h).
   Faces round a vertex with more than 2 pi of angle: five squares (5 pi/2)
   and six (the monkey saddle, 3 pi), flat and as hypars (each face's far
   corner pushed off its plane, the faces curved). Rays from a centre part
   at the vertex; a disk through it is the plane disk and, behind the
   vertex, the disk about it of the radius left. Uses figure 23's helpers. */

/* faces round the origin, face i spanned by axis directions d[i] and
   d[i + 1], its far corner pushed by w along its normal */
static void ssd_fan(qaws_surface** f, double const (*d)[3], unsigned int n, double w)
{
	static double const O[3] = { 0, 0, 0 };
	unsigned int i, k;
	for (i = 0; i < n; i++)
	{
		double const* a = d[i];
		double const* b = d[(i + 1) % n];
		double nn[3] = { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] }, p11[3];
		for (k = 0; k < 3; k++) p11[k] = a[k] + b[k] + w * nn[k];
		f[i] = shf_face(O, a, b, p11);
	}
}

static void ssd_faces(svg* s, qaws_surface** f, unsigned int n, sgf_view const* view, int curved)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		if (curved) sgf_mesh(s, f[i], view, 6, 0);
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

static void ssd_draw(svg* s, qaws_shell_region* g, qaws_shell_atlas const* a, qaws_surface** f, unsigned int n, sgf_view const* view,
	char const* fill, char const* line)
{
	unsigned int k;
	if (!g)
		return;
	if (fill) shf_region(s, g, a, view, fill);
	for (k = 0; k < n; k++)
		if (qaws_shell_region_get_face(g, k))
			sof_outline(s, qaws_shell_region_get_face(g, k), f[k], view, NULL, 0, 0, line, 1.4);
	qaws_shell_region_destroy(g);
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
	if (!svg_open(&s, "showcase/b2d24_shell_saddle.svg", (int)(48 + 6 * (pw + 12)), (int)(110 + (ph + 44)),
		"Geodesic disks through saddle vertices",
		"qaws_shell_region_create_geodesic_disk: faces round a vertex with more than 2 pi of angle; rays part there, the disk behind it is the one about the vertex of the radius left"))
		return;
	for (k = 0; k < 6; k++)
	{
		int six = k == 1 || k == 2 || k == 3 || k == 5, curved = k >= 3;
		unsigned int n = six ? 6 : 5;
		double w = curved ? 0.4 : 0, theta = n * pi / 2, A = 0;
		qaws_surface* f[6];
		qaws_shell_atlas* a = NULL;
		viewport v3;
		sgf_view view;
		ssd_fan(f, six ? D6 : D5, n, w);
		qaws_shell_atlas_create((qaws_surface const* const*)f, n, 0, &a);
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		if (six)
			sgf_view_init(&view, -0.75 * pi + 0.45, 0.85, 95, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		else
			sgf_view_init(&view, 0.25 * pi, 0.5, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 - 10);
		if (k == 0)
		{
			/* growing disks: short of the vertex, through it, well past it */
			static double const R[3] = { 0.75, 0.55, 0.35 };
			static char const* const F[3] = { "#b6e3ff", "#80ccff", "#54aeff" };
			svg_panel(&s, &v3, "five squares (5 pi/2): disks r 0.35, 0.55, 0.75");
			ssd_faces(&s, f, n, &view, 0);
			for (i = 0; i < 3; i++) ssd_draw(&s, ssd_disk(a, c, R[i], NULL), a, f, n, &view, F[i], "#0969da");
		}
		else if (k == 1)
		{
			/* a fan of straightest rays: those either side of the vertex part */
			svg_panel(&s, &v3, "monkey saddle (3 pi): rays part at the vertex");
			ssd_faces(&s, f, n, &view, 0);
			for (i = 0; i < 30; i++)
			{
				qaws_shell_curve* g = NULL;
				qaws_vec2 d;
				double th = 1.25 * pi + 0.6 * pi * ((i + 0.5) / 30.0 - 0.5);
				d.x = (qaws_scalar)cos(th); d.y = (qaws_scalar)sin(th);
				if (qaws_shell_curve_create_geodesic_ray(a, c, d, (qaws_scalar)1.3, 0, &g) == QAWS_STATUS_OK)
					shf_curve(&s, g, &view, a, 0, i < 15 ? "#cf222e" : "#8250df");
				qaws_shell_curve_destroy(g);
			}
		}
		else if (k == 4)
		{
			qaws_shell_region* g = ssd_disk(a, c, 0.6, &A);
			sprintf(label, "five hypars: disk r 0.6, area %.6f", A);
			svg_panel(&s, &v3, label);
			ssd_faces(&s, f, n, &view, 1);
			ssd_draw(&s, g, a, f, n, &view, "#54aeff", "#0969da");
		}
		else if (k == 5)
		{
			/* two disks on opposite sides of the vertex, united */
			qaws_shell_region *g1 = NULL, *g2 = NULL, *gu = NULL;
			qaws_scalar U = 0;
			qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)0.55, 0, &g1);
			qaws_shell_region_create_geodesic_disk(a, shf_pt(3, 0.35, 0.25), (qaws_scalar)0.5, 0, &g2);
			if (g1 && g2) qaws_shell_region_clip(g1, g2, QAWS_CLIP_UNION, &gu);
			if (gu) qaws_shell_region_get_area(gu, &U);
			sprintf(label, "monkey hypars: two disks united, area %.6f", (double)U);
			svg_panel(&s, &v3, label);
			ssd_faces(&s, f, n, &view, 1);
			if (gu) shf_region(&s, gu, a, &view, "#c297ff");
			for (i = 0; i < n; i++)
			{
				if (g1 && qaws_shell_region_get_face(g1, i))
					sof_outline(&s, qaws_shell_region_get_face(g1, i), f[i], &view, NULL, 0, 0, "#0969da", 1.4);
				if (g2 && qaws_shell_region_get_face(g2, i))
					sof_outline(&s, qaws_shell_region_get_face(g2, i), f[i], &view, NULL, 0, 0, "#cf222e", 1.4);
			}
			qaws_shell_region_destroy(g1); qaws_shell_region_destroy(g2); qaws_shell_region_destroy(gu);
		}
		else
		{
			qaws_shell_region* g = ssd_disk(a, c, 0.6, &A);
			if (curved)
				sprintf(label, "monkey hypars: disk r 0.6, area %.6f", A);
			else
				sprintf(label, "monkey saddle: %.6f, closed form %.6f", A, pi * 0.36 + 0.5 * (theta - 2 * pi) * (0.6 - rho0) * (0.6 - rho0));
			svg_panel(&s, &v3, label);
			ssd_faces(&s, f, n, &view, curved);
			ssd_draw(&s, g, a, f, n, &view, "#54aeff", "#0969da");
		}
		qaws_shell_atlas_destroy(a);
		for (i = 0; i < n; i++) qaws_surface_destroy(f[i]);
	}
	svg_close(&s);
	printf("figure 24: geodesic disks through saddle vertices, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
