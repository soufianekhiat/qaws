/*
 * Figure 26: self-intersections resolved (qaws_mesh_resolve_intersections).
 *   - two spheres crossing, as given (a self-intersecting soup);
 *   - resolved, the faces outside the other sphere kept: their union, the
 *     seam where they were cut conforming on both sides;
 *   - the faces inside the other kept: their intersection, the lens;
 *   - three cubes crossing at angles, resolved: their cuts meeting at
 *     triple points, coloured by cube.
 */

static void ms_ix_sphere(double const* c, double rad, unsigned int na, unsigned int nb, qaws_scalar* x, unsigned int* nv, unsigned int* t, unsigned int* nt)
{
	unsigned int i, j, v0 = *nv;
	for (i = 0; i < 2; i++)
	{
		x[3 * (v0 + i)] = (qaws_scalar)c[0];
		x[3 * (v0 + i) + 1] = (qaws_scalar)c[1];
		x[3 * (v0 + i) + 2] = (qaws_scalar)(c[2] + (i ? -rad : rad));
	}
	for (j = 1; j < nb; j++)
		for (i = 0; i < na; i++)
		{
			double th = MS_PI * j / nb, ph = 2 * MS_PI * i / na;
			unsigned int v = v0 + 2 + (j - 1) * na + i;
			x[3 * v] = (qaws_scalar)(c[0] + rad * sin(th) * cos(ph));
			x[3 * v + 1] = (qaws_scalar)(c[1] + rad * sin(th) * sin(ph));
			x[3 * v + 2] = (qaws_scalar)(c[2] + rad * cos(th));
		}
	for (i = 0; i < na; i++)
	{
		unsigned int* f = t + 3 * *nt;
		f[0] = v0; f[1] = v0 + 2 + i; f[2] = v0 + 2 + (i + 1) % na;
		f[3] = v0 + 1; f[4] = v0 + 2 + (nb - 2) * na + (i + 1) % na; f[5] = v0 + 2 + (nb - 2) * na + i;
		*nt += 2;
	}
	for (j = 0; j + 2 < nb; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = v0 + 2 + j * na + i, b = v0 + 2 + j * na + (i + 1) % na, cc = a + na, d = b + na;
			unsigned int* f = t + 3 * *nt;
			f[0] = a; f[1] = cc; f[2] = d;
			f[3] = a; f[4] = d; f[5] = b;
			*nt += 2;
		}
	*nv += 2 + na * (nb - 1);
}

static void ms_ix_box(double const* lo, double const* axis, double angle, qaws_scalar* x, unsigned int* nv, unsigned int* t, unsigned int* nt)
{
	static unsigned int const bt[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5 };
	static int const bx[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
	double co = cos(angle), si = sin(angle);
	unsigned int i, k;
	for (i = 0; i < 8; i++)
	{
		double p[3], q[3], d;
		for (k = 0; k < 3; k++)
			p[k] = bx[3 * i + k] - 0.5;
		d = axis[0] * p[0] + axis[1] * p[1] + axis[2] * p[2];
		q[0] = p[0] * co + (axis[1] * p[2] - axis[2] * p[1]) * si + axis[0] * d * (1 - co);
		q[1] = p[1] * co + (axis[2] * p[0] - axis[0] * p[2]) * si + axis[1] * d * (1 - co);
		q[2] = p[2] * co + (axis[0] * p[1] - axis[1] * p[0]) * si + axis[2] * d * (1 - co);
		for (k = 0; k < 3; k++)
			x[3 * (*nv + i) + k] = (qaws_scalar)(q[k] + lo[k] + 0.5);
	}
	for (i = 0; i < 36; i++)
		t[3 * *nt + i] = *nv + bt[i];
	*nv += 8;
	*nt += 12;
}

/* the resolved faces whose centroids pass keep(), coloured by their source's part (faces below split come from part 0) */
static qaws_mesh_3d* ms_ix_pick(qaws_mesh_3d const* m, int (*keep)(double const*, int), unsigned int const* split, double const (*rgb)[3],
	double** out_rgb)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), a = 0, i, k, n = 0;
	double* x = ms_positions(m);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	unsigned int* kt = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	qaws_scalar* src = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nt + 1));
	double* col = (double*)malloc(sizeof(double) * 3 * (nt + 1));
	qaws_mesh_3d* r;
	qaws_mesh_3d_get_triangles(t, m);
	qaws_mesh_3d_find_attribute(&a, m, "source");
	qaws_mesh_3d_get_attribute(src, m, a);
	for (i = 0; i < nt; i++)
	{
		double g[3] = { 0, 0, 0 };
		int part = 0;
		while (split[part] != 0 && (unsigned int)src[i] >= split[part])
			part++;
		for (k = 0; k < 3; k++)
		{
			g[0] += x[3 * t[3 * i + k]] / 3;
			g[1] += x[3 * t[3 * i + k] + 1] / 3;
			g[2] += x[3 * t[3 * i + k] + 2] / 3;
		}
		if (keep && !keep(g, part))
			continue;
		memcpy(kt + 3 * n, t + 3 * i, sizeof(unsigned int) * 3);
		memcpy(col + 3 * n, rgb[part], sizeof(double) * 3);
		n++;
	}
	r = ms_make(x, nv, kt, NULL, n);
	*out_rgb = col;
	free(x);
	free(t);
	free(kt);
	free(src);
	return r;
}

static double const ms_ix_c[2][4] = { { 0, 0, 0, 1 }, { 0.75, 0.3, 0.2, 0.8 } };

static int ms_ix_outside(double const* g, int part)
{
	double const* o = ms_ix_c[1 - part];
	return (g[0] - o[0]) * (g[0] - o[0]) + (g[1] - o[1]) * (g[1] - o[1]) + (g[2] - o[2]) * (g[2] - o[2]) > o[3] * o[3];
}

static int ms_ix_inside(double const* g, int part)
{
	return !ms_ix_outside(g, part);
}

static void demo_intersect(void)
{
	static double const col[3][3] = { { 150, 178, 220 }, { 236, 176, 110 }, { 150, 205, 150 } }, wire[3] = { 50, 55, 70 };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	unsigned int na = 48, nb = 24, nv = 0, nt = 0, i, split[3];
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (2 + na * nb));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 12 * na * nb);
	qaws_intersect_report r;
	qaws_mesh_3d *res = NULL, *pick;
	double* rgb = NULL;
	char sub[200];
	if (!c || !x || !t)
		return;
	ms_ix_sphere(ms_ix_c[0], ms_ix_c[0][3], na, nb, x, &nv, t, &nt);
	split[0] = nt;
	split[1] = 0;
	ms_ix_sphere(ms_ix_c[1], ms_ix_c[1][3], na, nb, x, &nv, t, &nt);
	{
		double* xd = (double*)malloc(sizeof(double) * 3 * nv);
		double* fc = (double*)malloc(sizeof(double) * 3 * nt);
		qaws_mesh_3d* in;
		for (i = 0; i < 3 * nv; i++)
			xd[i] = x[i];
		for (i = 0; i < nt; i++)
			memcpy(fc + 3 * i, col[i >= split[0]], sizeof(double) * 3);
		in = ms_make(xd, nv, t, NULL, nt);
		if (in)
		{
			sprintf(sub, "%u triangles crossing each other", nt);
			ms_panel(c, 0, 0, 550, 700, in, -1.2, 0.5, "two spheres crossing", sub);
			ms_draw_faces(c, in, fc);
			ms_draw(c, in, NULL, NULL, 0, wire, 0.4);
		}
		qaws_mesh_3d_destroy(in);
		free(xd);
		free(fc);
	}
	qaws_mesh_resolve_intersections(&res, &r, x, t, NULL, nv, nt, NULL);
	if (res)
	{
		printf("  spheres: %u pairs, %u points added, %u faces cut, %u faces, %u crossing left, %u unresolved\n", r.pairs, r.points, r.split_faces,
			qaws_mesh_3d_get_triangle_count(res), r.remaining, r.unresolved);
		pick = ms_ix_pick(res, ms_ix_outside, split, col, &rgb);
		if (pick)
		{
			sprintf(sub, "%u faces cut, %u crossing left", r.split_faces, r.remaining);
			ms_panel(c, 550, 0, 550, 700, pick, -1.2, 0.5, "resolved: outer faces kept", sub);
			ms_draw_faces(c, pick, rgb);
			ms_draw(c, pick, NULL, NULL, 0, wire, 0.4);
		}
		qaws_mesh_3d_destroy(pick);
		free(rgb);
		pick = ms_ix_pick(res, ms_ix_inside, split, col, &rgb);
		if (pick)
		{
			ms_panel(c, 1100, 0, 550, 700, pick, -1.2, 0.5, "inner faces kept: the lens", "the seam matched on both sides");
			ms_draw_faces(c, pick, rgb);
			ms_draw(c, pick, NULL, NULL, 0, wire, 0.4);
		}
		qaws_mesh_3d_destroy(pick);
		free(rgb);
	}
	qaws_mesh_3d_destroy(res);
	res = NULL;
	/* three cubes at angles */
	{
		static double const lo[3][3] = { { 0, 0, 0 }, { 0.45, 0.35, 0.25 }, { -0.3, 0.4, 0.5 } };
		static double const ax[3][3] = { { 0, 0, 1 }, { 0.57735026918962573, 0.57735026918962573, 0.57735026918962573 }, { 0.6, 0, 0.8 } };
		static double const ang[3] = { 0, 0.3, 0.45 };
		qaws_scalar bx[72];
		unsigned int bt[108], bn = 0, bf = 0, b;
		for (b = 0; b < 3; b++)
			ms_ix_box(lo[b], ax[b], ang[b], bx, &bn, bt, &bf);
		split[0] = 12;
		split[1] = 24;
		split[2] = 0;
		qaws_mesh_resolve_intersections(&res, &r, bx, bt, NULL, bn, bf, NULL);
		if (res)
		{
			pick = ms_ix_pick(res, NULL, split, col, &rgb);
			if (pick)
			{
				sprintf(sub, "%u pairs, %u faces, %u crossing left", r.pairs, qaws_mesh_3d_get_triangle_count(res), r.remaining);
				ms_panel(c, 1650, 0, 550, 700, pick, 0.75, 0.45, "three cubes crossing", sub);
				ms_draw_faces(c, pick, rgb);
				ms_draw(c, pick, NULL, NULL, 0, wire, 0.6);
				printf("  cubes: %s\n", sub);
			}
			qaws_mesh_3d_destroy(pick);
			free(rgb);
		}
		qaws_mesh_3d_destroy(res);
	}
	er_save(c, "showcase/mesh26_intersect.png");
	printf("mesh26_intersect.png\n");
	free(x);
	free(t);
	er_destroy(c);
}
