/*
 * Figure 36: the mesh connector, finished (qaws_mesh_2d_connect, chains,
 * attributes, qaws_volume_mesh_connect).
 *   - two planar meshes, a strip between facing sides of 8 and 3 edges: the
 *     counts graded, triangles zipped;
 *   - two tubes carrying a vertex attribute (0 on one, 1 on the other): the
 *     bridge blends it, coloured;
 *   - a block of 3 x 3 hexahedra and a turned copy above it, their facing
 *     patches joined by layers of hexahedra (the volume's border drawn).
 */

/* a planar grid of nx x ny cells over [x0, x0 + w] x [y0, y0 + h] */
static qaws_mesh_2d* ms_c2_grid(unsigned int nx, unsigned int ny, double x0, double y0, double w, double h)
{
	unsigned int i, j, k = 0, nv = (nx + 1) * (ny + 1);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nx * ny);
	qaws_mesh_2d_desc d;
	qaws_mesh_2d* m = NULL;
	for (j = 0; j <= ny; j++)
		for (i = 0; i <= nx; i++)
		{
			x[2 * (j * (nx + 1) + i)] = (qaws_scalar)(x0 + w * i / nx);
			x[2 * (j * (nx + 1) + i) + 1] = (qaws_scalar)(y0 + h * j / ny);
		}
	for (j = 0; j < ny; j++)
		for (i = 0; i < nx; i++)
		{
			unsigned int v = j * (nx + 1) + i;
			t[k++] = v; t[k++] = v + 1; t[k++] = v + nx + 2;
			t[k++] = v; t[k++] = v + nx + 2; t[k++] = v + nx + 1;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = 2 * nx * ny;
	qaws_mesh_2d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

/* a block of n x n x 1 hexahedra of side s, turned by ang about x, at c */
static qaws_volume_mesh* ms_c2_block(unsigned int n, double s, double ang, double const* c)
{
	unsigned int q = n + 1, i, j, k, cc = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * q * q * 2);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 8 * n * n);
	qaws_cell_type* ty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * n * n);
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* m = NULL;
	for (k = 0; k < 2; k++)
		for (j = 0; j < q; j++)
			for (i = 0; i < q; i++)
			{
				unsigned int v = (k * q + j) * q + i;
				double px = s * ((double)i - 0.5 * n), py = s * ((double)j - 0.5 * n), pz = s * k;
				x[3 * v] = (qaws_scalar)(c[0] + px);
				x[3 * v + 1] = (qaws_scalar)(c[1] + cos(ang) * py - sin(ang) * pz);
				x[3 * v + 2] = (qaws_scalar)(c[2] + sin(ang) * py + cos(ang) * pz);
			}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++, cc++)
		{
			unsigned int b0 = j * q + i, t0 = q * q + b0;
			unsigned int const h[8] = { b0, b0 + 1, b0 + q + 1, b0 + q, t0, t0 + 1, t0 + q + 1, t0 + q };
			memcpy(cv + 8 * cc, h, sizeof h);
			ty[cc] = QAWS_CELL_HEXAHEDRON;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 2 * q * q;
	d.cell_vertices = cv;
	d.cell_types = ty;
	d.cell_count = n * n;
	qaws_volume_mesh_create(&m, &d);
	free(x); free(cv); free(ty);
	return m;
}

static void demo_connect_more(void)
{
	static double const col[3] = { 180, 190, 214 }, tint[3] = { 236, 160, 90 }, wire[3] = { 60, 66, 80 };
	er_canvas* c = er_create(1800, 720, 3, 252, 252, 250);
	char sub[200];
	if (!c)
		return;
	{
		/* planar: a 8 x 8 grid and a 3 x 3 one up and to the right, their facing sides joined by a strip */
		qaws_mesh_2d *a = ms_c2_grid(8, 8, 0, 0, 1, 1), *b = ms_c2_grid(3, 3, 1.8, 0.5, 0.6, 0.6), *m = NULL;
		qaws_connect_desc d;
		qaws_connect_report r;
		memset(&d, 0, sizeof d);
		d.flags = QAWS_CONNECT_CHAINS;
		d.a_end = 8 * 9 + 8;
		d.b_end = 3 * 4;
		qaws_mesh_2d_connect(&m, &r, a, 8, b, 0, &d);
		if (m)
		{
			qaws_mesh_3d const* m3 = qaws_mesh_2d_get_mesh(m);
			sprintf(sub, "8 to 3 edges in %u rows, smallest angle %.0f deg", r.rows, (double)r.min_angle * 180 / MS_PI);
			ms_panel(c, 0, 0, 600, 720, m3, 0.0, 1.5707, "planar meshes, a strip", sub);
			ms_cn_draw(c, m3, col, tint, wire);
			printf("  planar strip: %s\n", sub);
		}
		qaws_mesh_2d_destroy(m);
		qaws_mesh_2d_destroy(a);
		qaws_mesh_2d_destroy(b);
	}
	{
		/* tubes carrying an attribute */
		double const ca[3] = { 0, 0, -1.2 }, da[3] = { 0, 0, 1 }, cb[3] = { 1.4, 0, 1.2 }, db[3] = { 0.7071067811865476, 0, 0.7071067811865476 };
		qaws_mesh_3d *a = ms_cn_tube(ca, da, 0.3, 1.2, 24, 6), *b = ms_cn_tube(cb, db, 0.3, 1.2, 24, 6), *m = NULL;
		qaws_attribute_desc ad;
		qaws_connect_desc d;
		qaws_connect_report r;
		unsigned int ia = 0, i;
		qaws_scalar v[24 * 7];
		for (i = 0; i < 24 * 7; i++)
			v[i] = 0;
		memset(&ad, 0, sizeof ad);
		ad.name = "heat";
		ad.element = QAWS_MESH_VERTEX;
		ad.components = 1;
		qaws_mesh_3d_add_attribute(&ia, a, v, &ad);
		for (i = 0; i < 24 * 7; i++)
			v[i] = 1;
		qaws_mesh_3d_add_attribute(&ia, b, v, &ad);
		memset(&d, 0, sizeof d);
		d.quads = 1;
		qaws_mesh_3d_connect(&m, &r, a, 24 * 6, b, 0, &d);
		if (m)
		{
			unsigned int nv = qaws_mesh_3d_get_vertex_count(m), io = 0;
			qaws_scalar* h = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
			double* rgb = (double*)malloc(sizeof(double) * 3 * nv);
			qaws_mesh_3d_find_attribute(&io, m, "heat");
			qaws_mesh_3d_get_attribute(h, m, io);
			for (i = 0; i < nv; i++)
				ms_colormap(h[i], rgb + 3 * i);
			sprintf(sub, "an attribute 0 and 1 blended over %u rows", r.rows);
			ms_panel(c, 600, 0, 600, 720, m, 0.3, 0.2, "attributes carried", sub);
			ms_draw(c, m, NULL, rgb, 1, wire, 0.3);
			printf("  attribute: %s\n", sub);
			free(h);
			free(rgb);
		}
		qaws_mesh_3d_destroy(m);
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
	}
	{
		/* volumes: a block and a turned copy above it */
		double const c0[3] = { 0, 0, 0 }, c1[3] = { 0.6, 0, 1.6 };
		qaws_volume_mesh *a = ms_c2_block(3, 0.4, 0, c0), *b = ms_c2_block(3, 0.4, 0.5, c1), *m = NULL;
		unsigned int top[9], bot[9], i;
		qaws_connect_desc d;
		qaws_connect_report r;
		for (i = 0; i < 9; i++)
		{
			top[i] = 6 * i + 1;
			bot[i] = 6 * i;
		}
		memset(&d, 0, sizeof d);
		d.rows = 4;
		qaws_volume_mesh_connect(&m, &r, a, top, b, bot, 9, 9, &d);
		if (m)
		{
			qaws_mesh_3d* bd = NULL;
			qaws_volume_mesh_get_boundary(&bd, m);
			sprintf(sub, "%u hexahedra between, scaled Jacobian %.2f and up", r.bridge_faces, (double)r.min_scaled_jacobian);
			ms_panel(c, 1200, 0, 600, 720, bd, 0.5, 0.25, "volumes: hexahedra", sub);
			ms_draw(c, bd, col, NULL, 0, wire, 0.8);
			printf("  volumes: %s\n", sub);
			qaws_mesh_3d_destroy(bd);
		}
		qaws_volume_mesh_destroy(m);
		qaws_volume_mesh_destroy(a);
		qaws_volume_mesh_destroy(b);
	}
	er_save(c, "showcase/mesh36_connect_more.png");
	printf("mesh36_connect_more.png\n");
	er_destroy(c);
}
