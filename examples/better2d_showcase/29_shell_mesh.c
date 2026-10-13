/* Figure 29: sampling and meshing across a shell's faces
   (qaws_shell_region_triangulate, qaws_shell_region_sample). A geodesic
   disk round a cube's corner spans three faces: its mesh is one, conformal
   along the cube's edges (each face's boundary edges split at the other
   faces' vertices there, no T-junctions), and its samples are shared out by
   the faces' areas. */

/* a segment in space */
static void s29_seg(ssd_scene* sc, qaws_vec3 a, qaws_vec3 b, int const* rgb, double width)
{
	double x0, y0, d0, x1, y1, d1;
	ssd_item* t = ssd_new(sc);
	if (!t)
		return;
	sgf_project(sc->view, a, &x0, &y0, &d0);
	sgf_project(sc->view, b, &x1, &y1, &d1);
	t->n = 2;
	t->xy[0] = x0; t->xy[1] = y0; t->xy[2] = x1; t->xy[3] = y1;
	t->depth = 0.5 * (d0 + d1) - 0.02;
	t->rgb[0] = rgb[0]; t->rgb[1] = rgb[1]; t->rgb[2] = rgb[2];
	t->width = width;
}

static void demo_shell_mesh(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 420, ph = 420, pi = 3.14159265358979323846;
	static double const P[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { 1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } };
	static int const C[6][4] = { { 4, 5, 6, 7 }, { 0, 2, 1, 3 }, { 1, 3, 5, 7 }, { 0, 4, 2, 6 }, { 2, 6, 3, 7 }, { 0, 1, 4, 5 } };
	static int const blue[3] = { 9, 105, 218 }, red[3] = { 207, 34, 46 }, dark[3] = { 87, 96, 106 };
	static int const fblue[3] = { 84, 174, 255 };
	qaws_surface* f[6];
	qaws_shell_atlas* a = NULL;
	qaws_shell_region* g = NULL;
	qaws_shell_mesh* m = NULL;
	qaws_shell_point c;
	qaws_surface_mesh_desc md;
	unsigned int i, k;
	char label[240];
	if (!svg_open(&s, "showcase/b2d29_shell_mesh.svg", (int)(48 + 2 * (pw + 12)), (int)(110 + (ph + 44)),
		"Sampling and meshing on shells",
		"qaws_shell_region_triangulate: one mesh across the faces, conformal at their edges; qaws_shell_region_sample: shared out by area"))
		return;
	for (i = 0; i < 6; i++) f[i] = shf_face(P[C[i][0]], P[C[i][1]], P[C[i][2]], P[C[i][3]]);
	qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &a);
	c.face = 0; c.uv.x = (qaws_scalar)0.85; c.uv.y = (qaws_scalar)0.85;
	qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)0.75, 0, &g);
	memset(&md, 0, sizeof(md));
	md.edge_length = (qaws_scalar)0.2;
	if (g) qaws_shell_region_triangulate(g, &md, &m);
	for (k = 0; k < 2; k++)
	{
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		ssd_paint paint[1];
		memset(&sc, 0, sizeof(sc));
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -0.75 * pi, 0.6, 110, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		sc.view = &view; sc.res = 40;
		paint[0].r = g; memcpy(paint[0].rgb, fblue, sizeof(fblue));
		ssd_surfaces(&sc, f, 6, k == 1 ? paint : NULL, k == 1 ? 1 : 0);
		if (k == 0 && m)
		{
			unsigned int nt = qaws_shell_mesh_get_triangle_count(m);
			for (i = 0; i < nt; i++)
			{
				unsigned int t[3], e;
				qaws_vec3 x[3];
				qaws_shell_mesh_get_triangle(m, i, t);
				for (e = 0; e < 3; e++) qaws_shell_mesh_get_vertex(m, t[e], NULL, &x[e]);
				for (e = 0; e < 3; e++)
					if (t[e] < t[(e + 1) % 3])
						s29_seg(&sc, x[e], x[(e + 1) % 3], blue, 0.6);
			}
			sprintf(label, "one mesh round the corner: %u triangles, %u vertices", nt, qaws_shell_mesh_get_vertex_count(m));
			ssd_panel(&s, &v3, label, "conformal at the cube's edges: no T-junctions where faces meet");
		}
		else
		{
			enum { N = 600 };
			static qaws_shell_point pts[N];
			unsigned int per[6] = { 0, 0, 0, 0, 0, 0 };
			if (g && qaws_shell_region_sample(g, NULL, N, pts, NULL) == QAWS_STATUS_OK)
				for (i = 0; i < N; i++)
				{
					ssd_dot(&sc, f[pts[i].face], pts[i].uv.x, pts[i].uv.y, red, 1.3);
					per[pts[i].face]++;
				}
			sprintf(label, "%u samples: %u, %u and %u on the three faces", N, per[0], per[2], per[4]);
			ssd_panel(&s, &v3, label, "shared out by the faces' areas (largest remainders), each face's part sampled");
		}
		ssd_outline(&sc, g, f, 6, dark);
		ssd_edges(&sc, f, 6);
		ssd_draw(&s, &sc);
	}
	svg_close(&s);
	qaws_shell_mesh_destroy(m);
	qaws_shell_region_destroy(g);
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 6; i++) qaws_surface_destroy(f[i]);
	printf("figure 29: sampling and meshing on shells, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
