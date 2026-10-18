/*
 * Figure 1: the corner store.
 *   - a triangle torus and its dual, read from the stored face adjacency:
 *     a vertex per face, a hexagon per vertex;
 *   - a cube of quads fanned into triangles, then refined 1-to-4 three times
 *     by the public edits alone (every edge split, the new diagonals
 *     flipped): its 12 creases and 8 corners kept, every edit checked;
 *   - an icosahedron refined the same way four times with a linear attribute
 *     interpolated by every edit: the colours are the attribute, exact.
 */

static qaws_mesh_3d* ms_torus(unsigned int nu, unsigned int nv, double R, double r)
{
	unsigned int i, j, k = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nu * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv);
	qaws_mesh_3d* m;
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * MS_PI * i / nu, v = 2 * MS_PI * (j + 0.5 * (i & 1)) / nv;
			x[3 * (j * nu + i)] = (R + r * cos(v)) * cos(u);
			x[3 * (j * nu + i) + 1] = (R + r * cos(v)) * sin(u);
			x[3 * (j * nu + i) + 2] = r * sin(v);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, d = ((j + 1) % nv) * nu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	m = ms_make(x, nu * nv, f, NULL, 2 * nu * nv);
	free(x);
	free(f);
	return m;
}

/* One 1-to-4 step through the public edits: every edge split at its
   midpoint, then every new edge from a new vertex to an old one that is not
   a half of the new vertex's edge flipped (the old triangle's middle one
   formed). Counts the edits. */
static void ms_subdivide(qaws_mesh_3d* m, unsigned int* edits)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), nv0 = qaws_mesh_3d_get_vertex_count(m), i, v;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	unsigned int* src = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(m, e);
	for (i = 0; i < ne; i++)
		if (qaws_mesh_3d_split_edge(m, e[2 * i], e[2 * i + 1], NULL, NULL, &v) == QAWS_STATUS_OK)
		{
			src[2 * (v - nv0)] = e[2 * i];
			src[2 * (v - nv0) + 1] = e[2 * i + 1];
			(*edits)++;
		}
	free(e);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(m, e);
	for (i = 0; i < ne; i++)
	{
		unsigned int p = e[2 * i], q = e[2 * i + 1], nw, old;
		if ((p >= nv0) == (q >= nv0))
			continue;
		nw = p >= nv0 ? p : q;
		old = p >= nv0 ? q : p;
		if (old != src[2 * (nw - nv0)] && old != src[2 * (nw - nv0) + 1] && qaws_mesh_3d_can_flip_edge(m, p, q)
			&& qaws_mesh_3d_flip_edge(m, p, q) == QAWS_STATUS_OK)
			(*edits)++;
	}
	free(e);
	free(src);
}

static void demo_core(void)
{
	static double const wire[3] = { 55, 60, 72 }, red[3] = { 214, 39, 40 }, dark[3] = { 20, 20, 20 };
	static double const tint[3] = { 178, 192, 212 }, warm[3] = { 226, 200, 160 };
	er_canvas* c = er_create(2200, 640, 3, 252, 252, 250);
	qaws_mesh_3d *torus = ms_torus(30, 14, 2, 0.75), *dual = NULL, *cube = NULL, *ico = NULL;
	qaws_mesh_topology t;
	char sub[200];
	unsigned int i, edits, nv, ok;
	if (!c || !torus)
		return;
	/* the torus and its dual */
	qaws_mesh_3d_get_topology(torus, &t);
	sprintf(sub, "%u triangles, chi %d, genus %u", qaws_mesh_3d_get_face_count(torus), t.euler_characteristic, t.genus);
	ms_panel(c, 0, 0, 550, 640, torus, 0.5, 0.75, "a triangle torus", sub);
	ms_draw(c, torus, tint, NULL, 1, wire, 1.1);
	qaws_mesh_3d_create_dual(torus, &dual);
	if (dual)
	{
		qaws_mesh_quality q;
		qaws_mesh_3d_measure_quality(dual, &q);
		sprintf(sub, "%u faces, %u hexagons: the stored adjacency", qaws_mesh_3d_get_face_count(dual), q.polygons);
		ms_panel(c, 550, 0, 550, 640, dual, 0.5, 0.75, "its dual", sub);
		ms_draw(c, dual, warm, NULL, 1, wire, 1.4);
	}
	/* the cube: fanned, detected, edited */
	{
		static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
		static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
		static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
		unsigned int nfe, corners = 0, *fe;
		double* p;
		cube = ms_make(x, 8, f, fs, 6);
		for (i = 0; i < 6; i++)
			qaws_mesh_3d_split_face(cube, i, NULL, NULL, NULL);
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)(30 * MS_PI / 180), 0, 0);
		edits = 0;
		ms_subdivide(cube, &edits);
		ms_subdivide(cube, &edits);
		ms_subdivide(cube, &edits);
		ok = qaws_mesh_3d_validate(cube) == QAWS_STATUS_OK;
		nfe = qaws_mesh_3d_get_feature_edge_count(cube);
		nv = qaws_mesh_3d_get_vertex_count(cube);
		for (i = 0; i < nv; i++)
			corners += qaws_mesh_3d_is_feature_corner(cube, i);
		sprintf(sub, "%u edits: %u crease edges, %u corners, valid %d", edits, nfe, corners, ok);
		ms_panel(c, 1100, 0, 550, 640, cube, 0.62, 0.55, "1-to-4 by splits and flips", sub);
		ms_draw(c, cube, tint, NULL, 0, wire, 1.0);
		fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe ? nfe : 1));
		qaws_mesh_3d_get_feature_edges(cube, fe);
		p = ms_positions(cube);
		for (i = 0; i < nfe; i++)
			er_line(c, p + 3 * fe[2 * i], p + 3 * fe[2 * i + 1], red, 3.2);
		for (i = 0; i < nv; i++)
			if (qaws_mesh_3d_is_feature_corner(cube, i))
				er_point(c, p + 3 * i, dark, 6);
		free(fe);
		free(p);
	}
	/* the icosahedron: a linear attribute through the edits */
	{
		static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
			3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
		double g = (1 + sqrt(5.0)) / 2, x[36], *p, *col, err = 0, lo = 1e300, hi = -1e300;
		double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
		qaws_scalar vals[12], *got;
		qaws_attribute_desc d;
		unsigned int a = 0;
		for (i = 0; i < 36; i++)
			x[i] = v[i];
		ico = ms_make(x, 12, f, NULL, 20);
		for (i = 0; i < 12; i++)
			vals[i] = (qaws_scalar)(x[3 * i] - 0.3 * x[3 * i + 1] + 0.2 * x[3 * i + 2]);
		memset(&d, 0, sizeof d);
		d.name = "f";
		d.element = QAWS_MESH_VERTEX;
		d.components = 1;
		d.rule = QAWS_ATTRIBUTE_LINEAR;
		qaws_mesh_3d_add_attribute(ico, &d, vals, &a);
		edits = 0;
		for (i = 0; i < 4; i++)
			ms_subdivide(ico, &edits);
		nv = qaws_mesh_3d_get_vertex_count(ico);
		p = ms_positions(ico);
		got = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		col = (double*)malloc(sizeof(double) * 3 * nv);
		qaws_mesh_3d_get_attribute(ico, a, got);
		for (i = 0; i < nv; i++)
		{
			double want = p[3 * i] - 0.3 * p[3 * i + 1] + 0.2 * p[3 * i + 2];
			err = fmax(err, fabs(got[i] - want));
			lo = fmin(lo, got[i]);
			hi = fmax(hi, got[i]);
		}
		for (i = 0; i < nv; i++)
			ms_colormap((got[i] - lo) / (hi - lo), col + 3 * i);
		sprintf(sub, "%u edits, %u vertices: error %.1e", edits, nv, err);
		ms_panel(c, 1650, 0, 550, 640, ico, 0.3, 0.45, "an attribute kept exact", sub);
		ms_draw(c, ico, NULL, col, 0, wire, 0.9);
		free(p); free(got); free(col);
	}
	er_save(c, "showcase/mesh1_core.png");
	printf("mesh1_core.png: torus %u faces, dual %u faces, cube %u vertices, icosahedron %u vertices\n", qaws_mesh_3d_get_face_count(torus),
		dual ? qaws_mesh_3d_get_face_count(dual) : 0, qaws_mesh_3d_get_vertex_count(cube), qaws_mesh_3d_get_vertex_count(ico));
	qaws_mesh_3d_destroy(torus);
	qaws_mesh_3d_destroy(dual);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(ico);
	er_destroy(c);
}
