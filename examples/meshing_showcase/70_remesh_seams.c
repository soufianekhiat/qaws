/*
 * Figure 70: attributes through remeshing (qaws_mesh_remesh, the input's own
 * reference). A cube with a colour per triangle (a seam along every side's
 * diagonal, where nothing creases), per-side UV charts (a seam on every
 * edge) and a value per edge (thick, one hue per value):
 *   - the input;
 *   - remeshed: the diagonal seams kept as lines (no face across one), the
 *     edges' values on their pieces;
 *   - its UVs as colours: each side's chart exact to its edges.
 * Printed: faces, features, faces across a seam.
 */

/* a hue for a value (0: none) */
static void ms_rs_hue(double v, double* rgb)
{
	double h = 2.4 * v;
	rgb[0] = 255 * (0.55 + 0.4 * cos(h));
	rgb[1] = 255 * (0.55 + 0.4 * cos(h - 2.0944));
	rgb[2] = 255 * (0.55 + 0.4 * cos(h + 2.0944));
}

/* the edges with a value, thick in its hue */
static void ms_rs_edges(er_canvas* c, qaws_mesh_3d const* m)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), i, el;
	unsigned int* ev = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
	double* v = ms_at_get(m, "e", &el);
	double bias = c->depth_bias;
	if (!v)
	{
		free(ev);
		return;
	}
	qaws_mesh_3d_get_edges(ev, m);
	c->depth_bias = 1e30;
	for (i = 0; i < ne; i++)
		if (v[i] > 0.5)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, ev[2 * i]), b = qaws_mesh_3d_get_vertex(m, ev[2 * i + 1]);
			double pa[3] = { a.x, a.y, a.z }, pb[3] = { b.x, b.y, b.z }, rgb[3];
			ms_rs_hue(v[i], rgb);
			er_line(c, pa, pb, rgb, 5);
		}
	c->depth_bias = bias;
	free(ev); free(v);
}

/* per-face colours from the per-triangle value, per-corner colours from the UVs */
static double* ms_rs_face_colours(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_face_count(m), i, el;
	double *t = ms_at_get(m, "tri", &el), *o = (double*)malloc(sizeof(double) * 3 * (nt + 1));
	for (i = 0; i < nt; i++)
		ms_rs_hue(t ? 1 + t[i] : 0, o + 3 * i);
	free(t);
	return o;
}

static double* ms_rs_uv_colours(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_face_count(m), i, el;
	double *uv = ms_at_get(m, "uv", &el), *o = (double*)malloc(sizeof(double) * 9 * (nt + 1));
	for (i = 0; uv && i < 3 * nt; i++)
	{
		/* (u's side and place as hue, v as brightness) */
		double u = uv[2 * i], v = uv[2 * i + 1], rgb[3];
		unsigned int k;
		ms_rs_hue(u, rgb);
		for (k = 0; k < 3; k++)
			o[3 * i + k] = rgb[k] * (0.55 + 0.45 * v);
	}
	free(uv);
	return o;
}

static void demo_remesh_seams(void)
{
	static unsigned int const F[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	static double const wire[3] = { 70, 70, 80 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	double x[24], uv[72], tri[12], *e;
	unsigned int i, k, ne;
	qaws_mesh_3d *m, *o = NULL;
	qaws_remesh_desc rd;
	char sub[200];
	if (!c)
		return;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (i == 1 || i == 2 || i == 5 || i == 6) ? 1 : 0;
		x[3 * i + 1] = (i == 2 || i == 3 || i == 6 || i == 7) ? 1 : 0;
		x[3 * i + 2] = i >= 4 ? 1 : 0;
	}
	m = ms_make(x, 8, F, NULL, 12);
	if (!m)
		return;
	qaws_mesh_3d_detect_features(m, (qaws_scalar)0.5, 0, 0);
	for (i = 0; i < 12; i++)
	{
		unsigned int ax = i < 4 ? 2 : i < 8 ? 1 : 0;
		tri[i] = i;
		for (k = 0; k < 3; k++)
		{
			unsigned int v = F[3 * i + k];
			uv[6 * i + 2 * k] = (ax == 0 ? x[3 * v + 1] : x[3 * v]) + 2 * (i / 2);
			uv[6 * i + 2 * k + 1] = ax == 2 ? x[3 * v + 1] : x[3 * v + 2];
		}
	}
	ms_at_add(m, "uv", QAWS_MESH_CORNER, 2, QAWS_ATTRIBUTE_LINEAR, uv, 36);
	ms_at_add(m, "tri", QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_NEAREST, tri, 12);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (double*)malloc(sizeof(double) * ne);
	{
		unsigned int* ev = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
		qaws_mesh_3d_get_edges(ev, m);
		for (i = 0; i < ne; i++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, ev[2 * i]), b = qaws_mesh_3d_get_vertex(m, ev[2 * i + 1]);
			int diff = (a.x != b.x) + (a.y != b.y) + (a.z != b.z);
			e[i] = diff == 1 ? 1 + (ev[2 * i] + ev[2 * i + 1]) % 7 : 0;
		}
		free(ev);
	}
	ms_at_add(m, "e", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, e, ne);
	free(e);
	{
		double* fc = ms_rs_face_colours(m);
		ms_panel(c, 0, 0, 600, 640, m, -0.6, 0.6, "the input", "a colour per triangle, a value per edge");
		ms_at_draw(c, m, NULL, fc, NULL, wire, 1.0);
		ms_rs_edges(c, m);
		free(fc);
	}
	memset(&rd, 0, sizeof rd);
	rd.target_length = (qaws_scalar)0.2;
	if (qaws_mesh_remesh(&o, NULL, m, &rd) == QAWS_STATUS_OK && o)
	{
		double *fc = ms_rs_face_colours(o), *cc = ms_rs_uv_colours(o);
		sprintf(sub, "%u faces, %u features (the cube's edges)", qaws_mesh_3d_get_face_count(o), qaws_mesh_3d_get_feature_edge_count(o));
		ms_panel(c, 600, 0, 600, 640, o, -0.6, 0.6, "remeshed: seams kept", sub);
		ms_at_draw(c, o, NULL, fc, NULL, wire, 0.6);
		ms_rs_edges(c, o);
		printf("  remeshed: %s\n", sub);
		ms_panel(c, 1200, 0, 600, 640, o, -0.6, 0.6, "its UV charts", "each side's exact to its edges");
		ms_at_draw(c, o, cc, NULL, NULL, wire, 0.6);
		free(fc); free(cc);
	}
	er_save(c, "showcase/mesh70_remesh_seams.png");
	printf("mesh70_remesh_seams.png\n");
	qaws_mesh_3d_destroy(m); qaws_mesh_3d_destroy(o);
	er_destroy(c);
}
