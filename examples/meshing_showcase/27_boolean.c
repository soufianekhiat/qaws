/*
 * Figure 27: mesh Booleans (qaws_mesh_boolean).
 *   - two spheres: their union, intersection and difference (the
 *     difference cut open to show the bite, its faces from B turned);
 *   - the CSG emblem: a cube and a sphere intersected, three cylinders
 *     taken away, each step a Boolean of the last result.
 */

static void ms_bo_cylinder(int axis, double rad, double len, unsigned int n, qaws_scalar* x, unsigned int* t, unsigned int* nv, unsigned int* nt)
{
	unsigned int i, k;
	for (i = 0; i < n; i++)
		for (k = 0; k < 2; k++)
		{
			double p[3], a = 2 * MS_PI * i / n;
			p[axis] = k ? len : -len;
			p[(axis + 1) % 3] = rad * cos(a);
			p[(axis + 2) % 3] = rad * sin(a);
			x[3 * (2 * i + k)] = (qaws_scalar)p[0];
			x[3 * (2 * i + k) + 1] = (qaws_scalar)p[1];
			x[3 * (2 * i + k) + 2] = (qaws_scalar)p[2];
		}
	/* the centres of the caps */
	for (k = 0; k < 2; k++)
	{
		double p[3] = { 0, 0, 0 };
		p[axis] = k ? len : -len;
		x[3 * (2 * n + k)] = (qaws_scalar)p[0];
		x[3 * (2 * n + k) + 1] = (qaws_scalar)p[1];
		x[3 * (2 * n + k) + 2] = (qaws_scalar)p[2];
	}
	*nt = 0;
	for (i = 0; i < n; i++)
	{
		unsigned int a0 = 2 * i, a1 = 2 * i + 1, b0 = 2 * ((i + 1) % n), b1 = b0 + 1;
		unsigned int* f = t + 3 * *nt;
		f[0] = a0; f[1] = b0; f[2] = b1;
		f[3] = a0; f[4] = b1; f[5] = a1;
		f[6] = 2 * n; f[7] = b0; f[8] = a0;
		f[9] = 2 * n + 1; f[10] = a1; f[11] = b1;
		*nt += 4;
	}
	*nv = 2 * n + 2;
}

/* a mesh's triangles as a desc over arrays it owns (freed by the caller) */
static qaws_mesh_3d_desc ms_bo_desc(qaws_mesh_3d const* m, qaws_scalar** x, unsigned int** t)
{
	qaws_mesh_3d_desc d;
	memset(&d, 0, sizeof d);
	d.vertex_count = qaws_mesh_3d_get_vertex_count(m);
	d.face_count = qaws_mesh_3d_get_triangle_count(m);
	*x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (d.vertex_count + 1));
	*t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (d.face_count + 1));
	qaws_mesh_3d_get_positions(*x, m);
	qaws_mesh_3d_get_triangles(*t, m);
	d.positions = *x;
	d.face_vertices = *t;
	return d;
}

static int ms_bo_front(double const* g, int part)
{
	(void)part;
	return g[1] < 0.15;   /* the difference cut open: the half facing the viewer removed */
}

static void demo_boolean(void)
{
	static double const col[3][3] = { { 150, 178, 220 }, { 236, 176, 110 }, { 150, 205, 150 } }, wire[3] = { 50, 55, 70 };
	static char const* const title[3] = { "union", "intersection", "difference, cut open" };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	unsigned int na = 48, nb = 24, nv = 0, nt = 0, nv0, nt0, i, k, split[2];
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (2 + na * nb));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 12 * na * nb);
	qaws_mesh_3d_desc A, B;
	qaws_boolean_report r;
	char sub[200];
	if (!c || !x || !t)
		return;
	ms_ix_sphere(ms_ix_c[0], ms_ix_c[0][3], na, nb, x, &nv, t, &nt);
	nv0 = nv;
	nt0 = nt;
	ms_ix_sphere(ms_ix_c[1], ms_ix_c[1][3], na, nb, x, &nv, t, &nt);
	for (i = 3 * nt0; i < 3 * nt; i++)
		t[i] -= nv0;
	memset(&A, 0, sizeof A);
	A.positions = x;
	A.vertex_count = nv0;
	A.face_vertices = t;
	A.face_count = nt0;
	B = A;
	B.positions = x + 3 * nv0;
	B.vertex_count = nv - nv0;
	B.face_vertices = t + 3 * nt0;
	B.face_count = nt - nt0;
	split[0] = nt0;
	split[1] = 0;
	for (k = 0; k < 3; k++)
	{
		qaws_mesh_3d *m = NULL, *pick;
		double* rgb = NULL;
		qaws_mesh_boolean(&m, &r, &A, &B, (qaws_boolean_op)k, NULL);
		if (!m)
			continue;
		pick = ms_ix_pick(m, k == 2 ? ms_bo_front : NULL, split, col, &rgb);
		if (pick)
		{
			sprintf(sub, "%u faces, closed %d", qaws_mesh_3d_get_triangle_count(m), r.closed);
			ms_panel(c, 550 * (int)k, 0, 550, 700, pick, -1.2, 0.5, title[k], sub);
			ms_draw_faces(c, pick, rgb);
			ms_draw(c, pick, NULL, NULL, 0, wire, 0.3);
			printf("  spheres, %s: %s\n", title[k], sub);
		}
		qaws_mesh_3d_destroy(pick);
		qaws_mesh_3d_destroy(m);
		free(rgb);
	}
	/* the CSG emblem */
	{
		static double const lo[3] = { -0.5, -0.5, -0.5 }, zero[3] = { 0, 0, 0 };
		qaws_scalar bx[24], cx[3 * 66];
		unsigned int bt[36], ct[3 * 128], bn = 0, bf = 0, cn, cf, step;
		qaws_mesh_3d *cur = NULL, *next = NULL;
		qaws_mesh_3d_desc D, E;
		ms_ix_box(lo, ms_ix_c[0], 0, bx, &bn, bt, &bf);   /* (the axis is unused at angle 0) */
		nv = 0;
		nt = 0;
		ms_ix_sphere(zero, 0.68, na, nb, x, &nv, t, &nt);
		memset(&D, 0, sizeof D);
		D.positions = bx;
		D.vertex_count = bn;
		D.face_vertices = bt;
		D.face_count = bf;
		E = D;
		E.positions = x;
		E.vertex_count = nv;
		E.face_vertices = t;
		E.face_count = nt;
		qaws_mesh_boolean(&cur, &r, &D, &E, QAWS_BOOLEAN_INTERSECTION, NULL);
		for (step = 0; cur && step < 3; step++)
		{
			qaws_scalar* px;
			unsigned int* pt;
			D = ms_bo_desc(cur, &px, &pt);
			ms_bo_cylinder((int)step, 0.3, 0.8, 32, cx, ct, &cn, &cf);
			memset(&E, 0, sizeof E);
			E.positions = cx;
			E.vertex_count = cn;
			E.face_vertices = ct;
			E.face_count = cf;
			next = NULL;
			qaws_mesh_boolean(&next, &r, &D, &E, QAWS_BOOLEAN_DIFFERENCE, NULL);
			qaws_mesh_3d_destroy(cur);
			cur = next;
			free(px);
			free(pt);
		}
		if (cur)
		{
			static double const warm[3] = { 226, 196, 150 };
			sprintf(sub, "five Booleans: %u faces, closed %d", qaws_mesh_3d_get_triangle_count(cur), r.closed);
			ms_panel(c, 1650, 0, 550, 700, cur, 0.7, 0.55, "the CSG emblem", sub);
			ms_draw(c, cur, warm, NULL, 0, wire, 0.3);
			printf("  emblem: %s, %u crossing left\n", sub, r.intersect.remaining);
		}
		qaws_mesh_3d_destroy(cur);
	}
	er_save(c, "showcase/mesh27_boolean.png");
	printf("mesh27_boolean.png\n");
	free(x);
	free(t);
	er_destroy(c);
}
