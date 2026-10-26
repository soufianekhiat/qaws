/*
 * Figure 62: more attributes (qaws_mesh.h's rules, qaws_mesh_3d_boolean,
 * qaws_mesh_sweep_boolean).
 *   - one quad's corners red, green, blue and yellow, subdivided five times
 *     (Catmull-Clark) under three rules: LINEAR blends the stored sRGB
 *     values (dark, muddy midpoints), SRGB blends in linear light, OKLAB in
 *     Oklab (perceptually even);
 *   - the CSG emblem in colour: a cube and a ball intersected, three
 *     cylinders taken away, each solid's colour and normals (flat faces,
 *     smooth ball and cylinders) carried through five Booleans; the feature
 *     edges drawn (the solids' creases and the cuts);
 *   - painting by sweeps: three strokes of coloured ball brushes united with
 *     a slab (qaws_mesh_sweep_boolean, UNION), each keeping its colour;
 *   - edge attributes: a box's edges weighted 0 to 1 (drawn from blue to
 *     red) less a ball: each edge's pieces keep its weight, the cut curves
 *     (weight 0) become features (black).
 * Printed: faces and features.
 */

/* a corner attribute filled per corner from its position and its face's normal */
typedef void (*ms_at2_fn)(double const* p, double const* fn, void const* user, double* out);

static void ms_at2_corners(qaws_mesh_3d* m, char const* name, unsigned int n, qaws_attribute_rule rule, ms_at2_fn fn, void const* user)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), f, j, q = 0;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	double* v = (double*)malloc(sizeof(double) * (n * nc + 1));
	qaws_mesh_3d_get_faces(fv, fs, m);
	for (f = 0; f < nf; q += fs[f++])
	{
		double nn[3] = { 0, 0, 0 }, l;
		for (j = 0; j < fs[f]; j++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, fv[q + j]), b = qaws_mesh_3d_get_vertex(m, fv[q + (j + 1) % fs[f]]);
			nn[0] += ((double)a.y - b.y) * ((double)a.z + b.z);
			nn[1] += ((double)a.z - b.z) * ((double)a.x + b.x);
			nn[2] += ((double)a.x - b.x) * ((double)a.y + b.y);
		}
		l = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
		for (j = 0; j < 3; j++) nn[j] = l > 0 ? nn[j] / l : 0;
		for (j = 0; j < fs[f]; j++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, fv[q + j]);
			double p[3] = { (double)a.x, (double)a.y, (double)a.z };
			fn(p, nn, user, v + (size_t)n * (q + j));
		}
	}
	ms_at_add(m, name, QAWS_MESH_CORNER, n, rule, v, nc);
	free(fv); free(fs); free(v);
}

static void ms_at2_constant(double const* p, double const* fn, void const* user, double* out)
{
	(void)p; (void)fn;
	memcpy(out, user, 3 * sizeof(double));
}

static void ms_at2_flat(double const* p, double const* fn, void const* user, double* out)
{
	(void)p; (void)user;
	memcpy(out, fn, 3 * sizeof(double));
}

/* radial from a centre (a ball) */
static void ms_at2_radial(double const* p, double const* fn, void const* user, double* out)
{
	double const* c = (double const*)user;
	double d[3] = { p[0] - c[0], p[1] - c[1], p[2] - c[2] }, l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
	unsigned int k;
	(void)fn;
	for (k = 0; k < 3; k++) out[k] = d[k] / l;
}

/* a cylinder along axis (user: the axis as a double): radial on the side, the face's on the caps */
static void ms_at2_cylinder(double const* p, double const* fn, void const* user, double* out)
{
	int axis = (int)*(double const*)user;
	double l;
	if (fabs(fn[axis]) > 0.9)
	{
		memcpy(out, fn, 3 * sizeof(double));
		return;
	}
	out[axis] = 0;
	out[(axis + 1) % 3] = p[(axis + 1) % 3];
	out[(axis + 2) % 3] = p[(axis + 2) % 3];
	l = sqrt(out[0] * out[0] + out[1] * out[1] + out[2] * out[2]);
	out[0] /= l; out[1] /= l; out[2] /= l;
}

/* a mesh's colour and normal per corner: constant colour, normals by fn */
static void ms_at2_paint(qaws_mesh_3d* m, double const* rgb, ms_at2_fn nfn, void const* user)
{
	ms_at2_corners(m, "color", 3, QAWS_ATTRIBUTE_SRGB, ms_at2_constant, rgb);
	ms_at2_corners(m, "normal", 3, QAWS_ATTRIBUTE_NORMAL, nfn, user);
}

static void ms_at2_features(er_canvas* c, qaws_mesh_3d const* m, double const* col, double width)
{
	unsigned int n = qaws_mesh_3d_get_feature_edge_count(m), i;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (n + 1));
	double* x = ms_positions(m);
	qaws_mesh_3d_get_feature_edges(e, m);
	for (i = 0; i < n; i++)
		er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], col, width);
	free(e); free(x);
}

static qaws_mesh_3d* ms_at2_cylinder_mesh(int axis, double rad, double len)
{
	qaws_scalar x[3 * 66];
	unsigned int t[3 * 128], nv, nt, i;
	double xd[3 * 66];
	ms_bo_cylinder(axis, rad, len, 32, x, t, &nv, &nt);
	for (i = 0; i < 3 * nv; i++) xd[i] = (double)x[i];
	return ms_make(xd, nv, t, NULL, nt);
}

static void demo_attributes_more(void)
{
	er_canvas* c = er_create(2400, 1280, 3, 252, 252, 250);
	char sub[200];
	unsigned int k;
	if (!c)
		return;
	{
		/* one quad's corner colours subdivided under three rules */
		static double const x[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
		static unsigned int const f[4] = { 0, 1, 2, 3 }, fs[1] = { 4 };
		static double const col[12] = { 1, 0, 0, 0, 0.8, 0, 0.1, 0.2, 1, 1, 0.9, 0 };
		static qaws_attribute_rule const rule[3] = { QAWS_ATTRIBUTE_LINEAR, QAWS_ATTRIBUTE_SRGB, QAWS_ATTRIBUTE_OKLAB };
		static char const* const title[3] = { "LINEAR: the stored values", "SRGB: in linear light", "OKLAB: perceptual" };
		for (k = 0; k < 3; k++)
		{
			qaws_mesh_3d *q = ms_make(x, 4, f, fs, 1), *s = NULL;
			qaws_mesh_subdivide_desc d;
			unsigned int el = 0;
			ms_at_add(q, "color", QAWS_MESH_VERTEX, 3, rule[k], col, 4);
			memset(&d, 0, sizeof d);
			d.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
			d.levels = 5;
			if (q && qaws_mesh_3d_subdivide(&s, q, &d) == QAWS_STATUS_OK && s)
			{
				/* (per vertex: each triangle's corners from its vertices) */
				unsigned int nt = qaws_mesh_3d_get_triangle_count(s), i, j;
				unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
				double *v = ms_at_get(s, "color", &el), *cc = (double*)malloc(sizeof(double) * 9 * nt), *xs = ms_positions(s);
				qaws_mesh_3d_get_triangles(t, s);
				for (i = 0; i < nt; i++)
					for (j = 0; j < 9; j++)
						cc[9 * i + j] = 255 * fmin(1, fmax(0, v[3 * t[3 * i + j / 3] + j % 3]));
				sprintf(sub, "%u faces, corners red green blue yellow", qaws_mesh_3d_get_face_count(s));
				ms_panel(c, 800 * (int)k, 0, 800, 640, s, 0, 1.5707, title[k], sub);
				for (i = 0; i < nt; i++)
				{
					static double const up[3] = { 0, 0, 1 };
					er_triangle(c, xs + 3 * t[3 * i], xs + 3 * t[3 * i + 1], xs + 3 * t[3 * i + 2], up, up, up, cc + 9 * i, cc + 9 * i + 3, cc + 9 * i + 6);
				}
				printf("  %s: %s\n", title[k], sub);
				free(t); free(v); free(cc); free(xs);
			}
			qaws_mesh_3d_destroy(q);
			qaws_mesh_3d_destroy(s);
		}
	}
	{
		/* the CSG emblem in colour */
		static double const blue[3] = { 0.25, 0.45, 0.85 }, gold[3] = { 0.95, 0.7, 0.2 }, cyl[3][3] = { { 0.85, 0.2, 0.25 }, { 0.2, 0.7, 0.35 }, { 0.6, 0.3, 0.8 } };
		static double const centre[3] = { 0, 0, 0 }, black[3] = { 30, 30, 36 };
		qaws_mesh_3d *cube = ms_sw_box(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5), *ball = ms_rm_sphere(48, 24), *cur = NULL;
		unsigned int nv = ball ? qaws_mesh_3d_get_vertex_count(ball) : 0, i;
		if (ball)
		{
			qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			qaws_mesh_3d_get_positions(x, ball);
			for (i = 0; i < 3 * nv; i++) x[i] = (qaws_scalar)(0.68 * (double)x[i]);
			qaws_mesh_3d_set_positions(ball, x);
			free(x);
		}
		ms_at2_paint(cube, blue, ms_at2_flat, NULL);
		ms_at2_paint(ball, gold, ms_at2_radial, centre);
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)0.5, 0, 0);
		qaws_mesh_3d_boolean(&cur, NULL, cube, ball, QAWS_BOOLEAN_INTERSECTION, NULL);
		for (k = 0; cur && k < 3; k++)
		{
			qaws_mesh_3d *cy = ms_at2_cylinder_mesh((int)k, 0.3, 0.8), *next = NULL;
			double axis = (double)k;
			ms_at2_paint(cy, cyl[k], ms_at2_cylinder, &axis);
			qaws_mesh_3d_detect_features(cy, (qaws_scalar)0.5, 0, 0);
			qaws_mesh_3d_boolean(&next, NULL, cur, cy, QAWS_BOOLEAN_DIFFERENCE, NULL);
			qaws_mesh_3d_destroy(cur);
			qaws_mesh_3d_destroy(cy);
			cur = next;
		}
		if (cur)
		{
			unsigned int ec = 0, en = 0;
			double *cc = ms_at_colours(cur, &ec), *cn = ms_at_get(cur, "normal", &en);
			sprintf(sub, "five Booleans: %u faces, %u feature edges", qaws_mesh_3d_get_face_count(cur), qaws_mesh_3d_get_feature_edge_count(cur));
			ms_panel(c, 0, 640, 800, 640, cur, 0.7, 0.55, "the CSG emblem, in colour", sub);
			ms_at_draw(c, cur, cc, NULL, cn, NULL, 0);
			ms_at2_features(c, cur, black, 0.9);
			printf("  emblem: %s\n", sub);
			free(cc); free(cn);
		}
		qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(ball); qaws_mesh_3d_destroy(cur);
	}
	{
		/* painting: three strokes of ball brushes united with a slab */
		static double const slab_c[3] = { 0.92, 0.9, 0.84 }, brush_c[3][3] = { { 0.9, 0.25, 0.3 }, { 0.2, 0.55, 0.9 }, { 0.25, 0.75, 0.35 } };
		static double const path[3][6] = { { -0.3, 0.6, 0.25, 4.3, 0.8, 0.25 }, { 0.2, 2.6, 0.25, 3.8, 1.2, 0.25 }, { 2.0, -0.2, 0.25, 2.2, 3.2, 0.25 } };
		static double const turn[3][2] = { { 0, 3 }, { 2, -3 }, { -3, 3 } };
		qaws_mesh_3d *cur = ms_sw_box(0, 4, 0, 3, -0.3, 0);
		double zero[3] = { 0, 0, 0 };
		ms_at2_paint(cur, slab_c, ms_at2_flat, NULL);
		for (k = 0; cur && k < 3; k++)
		{
			qaws_mesh_3d *brush = ms_rm_sphere(24, 12), *next = NULL;
			qaws_motion_pose s, e;
			qaws_motion* m = NULL;
			qaws_motion const* cm[1];
			qaws_sweep_desc d;
			unsigned int nv = brush ? qaws_mesh_3d_get_vertex_count(brush) : 0, i;
			qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
			qaws_mesh_3d_get_positions(x, brush);
			for (i = 0; i < 3 * nv; i++) x[i] = (qaws_scalar)(0.22 * (double)x[i]);
			qaws_mesh_3d_set_positions(brush, x);
			free(x);
			ms_at2_paint(brush, brush_c[k], ms_at2_radial, zero);
			memset(&s, 0, sizeof s);
			memset(&e, 0, sizeof e);
			s.rotation[0] = e.rotation[0] = 1;
			s.position.x = (qaws_scalar)path[k][0]; s.position.y = (qaws_scalar)path[k][1]; s.position.z = (qaws_scalar)path[k][2];
			e.position.x = (qaws_scalar)path[k][3]; e.position.y = (qaws_scalar)path[k][4]; e.position.z = (qaws_scalar)path[k][5];
			s.velocity.x = (qaws_scalar)(path[k][3] - path[k][0]); s.velocity.y = (qaws_scalar)turn[k][0];
			e.velocity.x = (qaws_scalar)(path[k][3] - path[k][0]); e.velocity.y = (qaws_scalar)turn[k][1];
			if (k == 2)
			{
				s.velocity.x = (qaws_scalar)turn[k][0]; s.velocity.y = (qaws_scalar)(path[k][4] - path[k][1]);
				e.velocity.x = (qaws_scalar)turn[k][1]; e.velocity.y = (qaws_scalar)(path[k][4] - path[k][1]);
			}
			qaws_motion_create_poses(&m, &s, &e, 0);
			cm[0] = m;
			memset(&d, 0, sizeof d);
			d.tolerance = (qaws_scalar)0.01;
			if (m && qaws_mesh_sweep_boolean(&next, NULL, brush, cm, cur, 1, QAWS_BOOLEAN_UNION, &d) == QAWS_STATUS_OK && next)
			{
				qaws_mesh_3d_destroy(cur);
				cur = next;
			}
			qaws_motion_destroy(m);
			qaws_mesh_3d_destroy(brush);
		}
		if (cur)
		{
			unsigned int ec = 0, en = 0;
			double *cc = ms_at_colours(cur, &ec), *cn = ms_at_get(cur, "normal", &en);
			sprintf(sub, "three strokes united: %u faces", qaws_mesh_3d_get_face_count(cur));
			ms_panel(c, 800, 640, 800, 640, cur, -0.35, 0.95, "painting by sweeps", sub);
			ms_at_draw(c, cur, cc, NULL, cn, NULL, 0);
			printf("  painting: %s\n", sub);
			free(cc); free(cn);
		}
		qaws_mesh_3d_destroy(cur);
	}
	{
		/* edge weights through a bite */
		static double const grey[3] = { 225, 225, 228 }, black[3] = { 20, 20, 24 };
		qaws_mesh_3d *box = ms_sw_box(-1, 1, -0.7, 0.7, -0.5, 0.5), *ball = ms_rm_sphere(48, 24), *o = NULL;
		unsigned int nv = ball ? qaws_mesh_3d_get_vertex_count(ball) : 0, i;
		{
			qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			qaws_mesh_3d_get_positions(x, ball);
			for (i = 0; i < nv; i++)
			{
				x[3 * i] = (qaws_scalar)(0.75 * (double)x[3 * i] + 0.9);
				x[3 * i + 1] = (qaws_scalar)(0.75 * (double)x[3 * i + 1] - 0.55);
				x[3 * i + 2] = (qaws_scalar)(0.75 * (double)x[3 * i + 2] + 0.4);
			}
			qaws_mesh_3d_set_positions(ball, x);
			free(x);
		}
		{
			/* each box edge weighted by its direction and place (diagonals 0) */
			unsigned int ne = qaws_mesh_3d_get_edge_count(box), j = 0;
			unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
			double* w = (double*)malloc(sizeof(double) * ne);
			qaws_mesh_3d_get_edges(e, box);
			for (i = 0; i < ne; i++)
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(box, e[2 * i]), q = qaws_mesh_3d_get_vertex(box, e[2 * i + 1]);
				int same = (p.x == q.x) + (p.y == q.y) + (p.z == q.z);
				w[i] = same == 2 ? 0.1 + 0.9 * ((j++ * 5) % 12) / 11.0 : 0;
			}
			ms_at_add(box, "crease", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, w, ne);
			free(e); free(w);
		}
		qaws_mesh_3d_detect_features(box, (qaws_scalar)0.5, 0, 0);
		if (box && ball && qaws_mesh_3d_boolean(&o, NULL, box, ball, QAWS_BOOLEAN_DIFFERENCE, NULL) == QAWS_STATUS_OK && o)
		{
			unsigned int el = 0, ne = qaws_mesh_3d_get_edge_count(o);
			unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
			double *w = ms_at_get(o, "crease", &el), *x = ms_positions(o);
			sprintf(sub, "%u faces; weights per %s, %u feature edges", qaws_mesh_3d_get_face_count(o), el == QAWS_MESH_EDGE ? "edge" : "?", qaws_mesh_3d_get_feature_edge_count(o));
			ms_panel(c, 1600, 640, 800, 640, o, -0.6, 0.6, "edge weights through a bite", sub);
			ms_draw(c, o, grey, NULL, 0, NULL, 0);
			ms_at2_features(c, o, black, 1.2);
			qaws_mesh_3d_get_edges(e, o);
			for (i = 0; w && i < ne; i++)
				if (w[i] > 0)
				{
					double col[3] = { 255 * w[i], 60, 255 * (1 - w[i]) };
					er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], col, 3.2);
				}
			printf("  edge weights: %s\n", sub);
			free(e); free(w); free(x);
		}
		qaws_mesh_3d_destroy(box); qaws_mesh_3d_destroy(ball); qaws_mesh_3d_destroy(o);
	}
	er_save(c, "showcase/mesh62_attributes_more.png");
	printf("mesh62_attributes_more.png\n");
	er_destroy(c);
}
