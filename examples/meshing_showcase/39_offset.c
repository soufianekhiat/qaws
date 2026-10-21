/*
 * Figure 39: offsets of meshes (qaws_mesh_3d_offset).
 *   - a U offset outward by contouring the distance: rounded about its
 *     edges, the slot narrowed (the input inside, solid);
 *   - the U moved outward by more than half its slot (mitred): flat faces
 *     and sharp edges kept, the slot closed, the crossing faces resolved;
 *   - a sphere's shell (contoured), cut open to show the wall;
 *   - a wavy sheet thickened by moving it both ways, the rim joined.
 */

static void ms_of_box(double* x, unsigned int* t, double const* lo, double const* hi, unsigned int base)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (i == 1 || i == 2 || i == 5 || i == 6) ? hi[0] : lo[0];
		x[3 * i + 1] = (i == 2 || i == 3 || i == 6 || i == 7) ? hi[1] : lo[1];
		x[3 * i + 2] = i >= 4 ? hi[2] : lo[2];
	}
	for (i = 0; i < 36; i++) t[i] = f[i] + base;
}

/* a U: a base and two prongs (one soup), made one by the self-union */
static qaws_mesh_3d* ms_of_u(void)
{
	static double const b[6][3] = { { 0, 0, 0 }, { 1, 0.3, 0.3 }, { 0, 0.25, 0 }, { 0.3, 1, 0.3 }, { 0.7, 0.25, 0 }, { 1, 1, 0.3 } };
	double x[72];
	qaws_scalar xs[72];
	unsigned int t[108], i;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* u = NULL;
	for (i = 0; i < 3; i++) ms_of_box(x + 24 * i, t + 36 * i, b[2 * i], b[2 * i + 1], 8 * i);
	for (i = 0; i < 72; i++) xs[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = xs;
	d.vertex_count = 24;
	d.face_vertices = t;
	d.face_count = 36;
	qaws_mesh_boolean(&u, NULL, &d, NULL, QAWS_BOOLEAN_UNION, NULL);
	return u;
}

/* a mesh's triangles whose centroids lie past y = cut_y (the side toward the viewer cut away), flat shaded, outlined */
static void ms_of_draw_cut(er_canvas* c, qaws_mesh_3d const* m, double const* rgb, double const* wire, double cut_y)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	double* x = ms_positions(m);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		double const *a = x + 3 * t[3 * i], *b = x + 3 * t[3 * i + 1], *d = x + 3 * t[3 * i + 2];
		if ((a[1] + b[1] + d[1]) / 3 < cut_y)
			continue;
		er_triangle(c, a, b, d, NULL, NULL, NULL, rgb, rgb, rgb);
		er_line(c, a, b, wire, 0.35);
		er_line(c, b, d, wire, 0.35);
		er_line(c, d, a, wire, 0.35);
	}
	free(x);
	free(t);
}

static void demo_offset(void)
{
	static double const solid[3] = { 120, 150, 210 }, shell[3] = { 230, 160, 100 }, wire[3] = { 60, 66, 80 }, light[3] = { 150, 156, 170 },
		green[3] = { 130, 190, 140 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	char sub[200];
	qaws_mesh_offset_desc d;
	qaws_mesh_offset_report r;
	qaws_mesh_3d *u = ms_of_u(), *o = NULL;
	if (!c || !u)
		return;
	/* contoured outward */
	memset(&d, 0, sizeof d);
	d.kind = QAWS_MESH_OFFSET_SURFACE;
	d.method = QAWS_MESH_OFFSET_CONTOUR;
	d.distance = (qaws_scalar)0.12;
	d.resolution = 56;
	qaws_mesh_3d_offset(&o, &r, u, &d);
	if (o)
	{
		sprintf(sub, "0.12 out: %.3f to %.3f away", (double)r.min_distance, (double)r.max_distance);
		ms_panel(c, 0, 0, 450, 640, o, 0.6, 0.55, "contoured", sub);
		ms_draw(c, u, solid, NULL, 0, wire, 0.6);
		ms_draw(c, o, NULL, NULL, 0, light, 0.5);
		printf("  contoured: %s\n", sub);
		qaws_mesh_3d_destroy(o);
		o = NULL;
	}
	/* moved: the slot closed */
	d.method = QAWS_MESH_OFFSET_MOVE;
	d.distance = (qaws_scalar)0.25;
	qaws_mesh_3d_offset(&o, &r, u, &d);
	if (o)
	{
		sprintf(sub, "0.25 out: slot shut, %u crossings", r.intersect.pairs);
		ms_panel(c, 450, 0, 450, 640, o, 0.6, 0.55, "moved, mitred", sub);
		ms_draw(c, o, green, NULL, 0, wire, 0.6);
		printf("  moved: %s\n", sub);
		qaws_mesh_3d_destroy(o);
		o = NULL;
	}
	{
		/* the shell of a sphere, cut open */
		qaws_mesh_3d* s = ms_rm_sphere(48, 24);
		d.kind = QAWS_MESH_OFFSET_SHELL;
		d.method = QAWS_MESH_OFFSET_CONTOUR;
		d.distance = (qaws_scalar)0.25;
		d.resolution = 56;
		qaws_mesh_3d_offset(&o, &r, s, &d);
		if (o)
		{
			sprintf(sub, "0.25 thick: inside %.3f to %.3f in", (double)r.min_distance, (double)r.max_distance);
			ms_panel(c, 900, 0, 450, 640, s, 0.5, 0.35, "a shell, cut open", sub);
			ms_of_draw_cut(c, o, shell, wire, -0.2);
			printf("  shell: %s\n", sub);
			qaws_mesh_3d_destroy(o);
			o = NULL;
		}
		qaws_mesh_3d_destroy(s);
	}
	{
		/* a wavy sheet thickened */
		unsigned int n = 28, i, j, k = 0, nv = (n + 1) * (n + 1);
		double* x = (double*)malloc(sizeof(double) * 3 * nv);
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
		qaws_mesh_3d* s;
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				double* p = x + 3 * (j * (n + 1) + i), a = (double)i / n, b = (double)j / n;
				p[0] = a; p[1] = b; p[2] = 0.12 * sin(2 * MS_PI * a) * cos(1.5 * MS_PI * b);
			}
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int v = j * (n + 1) + i;
				t[k++] = v; t[k++] = v + 1; t[k++] = v + n + 2;
				t[k++] = v; t[k++] = v + n + 2; t[k++] = v + n + 1;
			}
		s = ms_make(x, nv, t, NULL, 2 * n * n);
		d.kind = QAWS_MESH_OFFSET_THICKEN;
		d.method = QAWS_MESH_OFFSET_MOVE;
		d.distance = (qaws_scalar)0.08;
		qaws_mesh_3d_offset(&o, &r, s, &d);
		if (o)
		{
			sprintf(sub, "0.08 thick, the rim joined: closed %d", r.closed);
			ms_panel(c, 1350, 0, 450, 640, o, 0.7, 0.45, "a sheet thickened", sub);
			ms_draw(c, o, solid, NULL, 0, wire, 0.4);
			printf("  thickened: %s\n", sub);
			qaws_mesh_3d_destroy(o);
		}
		qaws_mesh_3d_destroy(s);
		free(x);
		free(t);
	}
	qaws_mesh_3d_destroy(u);
	er_save(c, "showcase/mesh39_offset.png");
	printf("mesh39_offset.png\n");
	er_destroy(c);
}
