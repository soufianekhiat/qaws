/* Figure 23: curves and regions across the faces of a shell
   (qaws_shell_region.h). A cube of six flat faces: a fan of straightest
   geodesics from one point over the edges, a geodesic circle over two
   edges by a vertex, two circles over edges united; a sphere as two
   hemispheres: circles over the equator. Uses figures 16 to 19's
   helpers. */

static qaws_surface* shf_face(double const* p00, double const* p10, double const* p01, double const* p11)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.p00.x = (qaws_scalar)p00[0]; d.p00.y = (qaws_scalar)p00[1]; d.p00.z = (qaws_scalar)p00[2];
	d.p10.x = (qaws_scalar)p10[0]; d.p10.y = (qaws_scalar)p10[1]; d.p10.z = (qaws_scalar)p10[2];
	d.p01.x = (qaws_scalar)p01[0]; d.p01.y = (qaws_scalar)p01[1]; d.p01.z = (qaws_scalar)p01[2];
	d.p11.x = (qaws_scalar)p11[0]; d.p11.y = (qaws_scalar)p11[1]; d.p11.z = (qaws_scalar)p11[2];
	qaws_surface_create_bilinear(&d, &s);
	return s;
}

/* a face's outline (its four sides), faded where it faces away */
static void shf_face_box(svg* s, qaws_surface const* f, sgf_view const* view)
{
	double uv[2 * 5] = { 0, 0, 1, 0, 1, 1, 0, 1, 0, 0 };
	double pts[2 * 41];
	unsigned int k, q;
	for (k = 0; k < 4; k++)
	{
		for (q = 0; q <= 40; q++)
		{
			pts[2 * q] = uv[2 * k] + (uv[2 * k + 2] - uv[2 * k]) * q / 40.0;
			pts[2 * q + 1] = uv[2 * k + 1] + (uv[2 * k + 3] - uv[2 * k + 1]) * q / 40.0;
		}
		sgf_curve3d(s, f, view, pts, 41, 0, "#57606a", 1.0, 0);
	}
}

static void shf_curve(svg* s, qaws_shell_curve const* c, sgf_view const* view, qaws_shell_atlas const* a, double pu, char const* color)
{
	unsigned int i;
	for (i = 0; i < qaws_shell_curve_get_run_count(c); i++)
	{
		unsigned int f;
		qaws_surface_curve const* sc;
		qaws_shell_curve_get_run(c, i, &f, &sc);
		ssf_curve(s, sc, qaws_shell_atlas_get_face(a, f), view, NULL, pu, 0, color);
	}
}

static void shf_region(svg* s, qaws_shell_region const* r, qaws_shell_atlas const* a, sgf_view const* view, char const* fill)
{
	unsigned int f;
	for (f = 0; f < qaws_shell_atlas_get_face_count(a); f++)
		if (qaws_shell_region_get_face(r, f))
			ssf_fill3d(s, qaws_shell_region_get_face(r, f), qaws_shell_atlas_get_face(a, f), view, fill);
}

static qaws_shell_point shf_pt(unsigned int face, double u, double v)
{
	qaws_shell_point p;
	p.face = face; p.uv.x = (qaws_scalar)u; p.uv.y = (qaws_scalar)v;
	return p;
}

static void demo_shell_region(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 300, pi = 3.14159265358979323846;
	static double const P[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { 1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } };
	qaws_surface* f[6];
	qaws_shell_atlas* a = NULL;
	unsigned int i, k;
	char label[220];
	if (!svg_open(&s, "showcase/b2d23_shell_region.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + (ph + 44)),
		"Curves and regions across the faces of a shell",
		"qaws_shell_region: straightest geodesics continued over edges, geodesic circles across faces, regions and Booleans face by face"))
		return;
	f[0] = shf_face(P[4], P[5], P[6], P[7]);
	f[1] = shf_face(P[0], P[2], P[1], P[3]);
	f[2] = shf_face(P[1], P[3], P[5], P[7]);
	f[3] = shf_face(P[0], P[4], P[2], P[6]);
	f[4] = shf_face(P[2], P[6], P[3], P[7]);
	f[5] = shf_face(P[0], P[1], P[4], P[5]);
	qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &a);
	for (k = 0; k < 3; k++)
	{
		viewport v3;
		sgf_view view;
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		/* from the corner (1, 1, 1): the top, +x and +y faces in view */
		sgf_view_init(&view, -0.75 * pi, 0.6, 62, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		if (k == 0)
		{
			/* a fan of straightest rays from a point of the top face */
			svg_panel(&s, &v3, "a fan of straightest geodesics over edges");
			for (i = 0; i < 6; i++) shf_face_box(&s, f[i], &view);
			for (i = 0; i < 24; i++)
			{
				qaws_shell_curve* c = NULL;
				qaws_vec2 d;
				double th = 2 * pi * (i + 0.5) / 24;
				d.x = (qaws_scalar)cos(th); d.y = (qaws_scalar)sin(th);
				/* run until a vertex or 3.2 */
				if (qaws_shell_curve_create_geodesic_ray(a, shf_pt(0, 0.6, 0.45), d, (qaws_scalar)3.2, 0, &c) == QAWS_STATUS_OK)
					shf_curve(&s, c, &view, a, 0, "#cf222e");
				qaws_shell_curve_destroy(c);
			}
		}
		else if (k == 1)
		{
			qaws_shell_curve* c = NULL;
			qaws_shell_region* r = NULL;
			qaws_scalar A = 0;
			qaws_shell_curve_create_geodesic_circle(a, shf_pt(0, 0.8, 0.8), (qaws_scalar)0.55, 0, &c);
			if (c) qaws_shell_region_create(a, (qaws_shell_curve const* const*)&c, 1, &r);
			if (r) qaws_shell_region_get_area(r, &A);
			sprintf(label, "a circle over two edges: %.6f, pi r^2 %.6f", (double)A, pi * 0.55 * 0.55);
			svg_panel(&s, &v3, label);
			for (i = 0; i < 6; i++) shf_face_box(&s, f[i], &view);
			if (r) shf_region(&s, r, a, &view, "#54aeff");
			if (c) shf_curve(&s, c, &view, a, 0, "#0969da");
			qaws_shell_region_destroy(r);
			qaws_shell_curve_destroy(c);
		}
		else
		{
			qaws_shell_curve *c1 = NULL, *c2 = NULL;
			qaws_shell_region *r1 = NULL, *r2 = NULL, *ru = NULL;
			qaws_scalar A = 0;
			qaws_shell_curve_create_geodesic_circle(a, shf_pt(0, 0.75, 0.5), (qaws_scalar)0.7, 0, &c1);
			qaws_shell_curve_create_geodesic_circle(a, shf_pt(2, 0.6, 0.8), (qaws_scalar)0.6, 0, &c2);
			if (c1) qaws_shell_region_create(a, (qaws_shell_curve const* const*)&c1, 1, &r1);
			if (c2) qaws_shell_region_create(a, (qaws_shell_curve const* const*)&c2, 1, &r2);
			if (r1 && r2) qaws_shell_region_clip(r1, r2, QAWS_CLIP_UNION, &ru);
			if (ru) qaws_shell_region_get_area(ru, &A);
			sprintf(label, "two circles over edges, united: area %.6f", (double)A);
			svg_panel(&s, &v3, label);
			for (i = 0; i < 6; i++) shf_face_box(&s, f[i], &view);
			if (ru) shf_region(&s, ru, a, &view, "#8250df");
			if (c1) shf_curve(&s, c1, &view, a, 0, "#0969da");
			if (c2) shf_curve(&s, c2, &view, a, 0, "#cf222e");
			qaws_shell_region_destroy(r1); qaws_shell_region_destroy(r2); qaws_shell_region_destroy(ru);
			qaws_shell_curve_destroy(c1); qaws_shell_curve_destroy(c2);
		}
	}
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 6; i++) qaws_surface_destroy(f[i]);
	/* two hemispheres */
	{
		qaws_curve *p1 = sgf_arc(0, 0, 1, -pi / 2, 0), *p2 = sgf_arc(0, 0, 1, 0, pi / 2);
		qaws_surface* h[2];
		qaws_shell_atlas* ah = NULL;
		viewport v3;
		sgf_view view;
		qaws_shell_curve *c1 = NULL, *c2 = NULL;
		qaws_shell_region *r1 = NULL, *r2 = NULL, *rx = NULL;
		qaws_scalar A = 0;
		h[0] = sgf_revolve(p1);
		h[1] = sgf_revolve(p2);
		qaws_shell_atlas_create((qaws_surface const* const*)h, 2, 0, &ah);
		v3.x0 = 24 + 3 * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -pi / 2 - 2 * pi * 0.3 + 0.2, 0.3, 120, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		qaws_shell_curve_create_geodesic_circle(ah, shf_pt(1, 0.3, 0.08), (qaws_scalar)0.55, 0, &c1);
		qaws_shell_curve_create_geodesic_circle(ah, shf_pt(0, 0.36, 0.85), (qaws_scalar)0.45, 0, &c2);
		if (c1) qaws_shell_region_create(ah, (qaws_shell_curve const* const*)&c1, 1, &r1);
		if (c2) qaws_shell_region_create(ah, (qaws_shell_curve const* const*)&c2, 1, &r2);
		if (r1 && r2) qaws_shell_region_clip(r1, r2, QAWS_CLIP_XOR, &rx);
		if (rx) qaws_shell_region_get_area(rx, &A);
		sprintf(label, "hemispheres: circles over the equator, xor");
		svg_panel(&s, &v3, label);
		sgf_mesh(&s, h[0], &view, 12, 1);
		sgf_mesh(&s, h[1], &view, 12, 1);
		if (rx) shf_region(&s, rx, ah, &view, "#54aeff");
		if (c1) shf_curve(&s, c1, &view, ah, 1, "#0969da");
		if (c2) shf_curve(&s, c2, &view, ah, 1, "#cf222e");
		qaws_shell_region_destroy(r1); qaws_shell_region_destroy(r2); qaws_shell_region_destroy(rx);
		qaws_shell_curve_destroy(c1); qaws_shell_curve_destroy(c2);
		qaws_shell_atlas_destroy(ah);
		qaws_surface_destroy(h[0]); qaws_surface_destroy(h[1]);
		qaws_curve_destroy(p1); qaws_curve_destroy(p2);
	}
	svg_close(&s);
	printf("figure 23: curves and regions across the faces of a shell, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
