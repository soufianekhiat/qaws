/*
 * Figure 28: geodesics on meshes (qaws_mesh_geodesic, the heat method).
 *   - a torus: the distance from one point, its isolines, and paths back
 *     to the point from round the torus (some over the top, some through
 *     the hole);
 *   - a bumpy sphere with three sources: each vertex's distance to the
 *     nearest, the isolines meeting where two sources are as near, the
 *     regions' borders in red, placed inside the triangles from each
 *     source's exact distances (qaws_mesh_3d_geodesic_borders);
 *   - a sheet with a hill, the source at its left side: paths from the right
 *     side go round the hill, not over it.
 * Printed: the factoring and query times.
 */

/* isolines of d (per vertex) every step, drawn on the triangles */
static void ms_gd_isolines(er_canvas* c, qaws_mesh_3d const* m, qaws_scalar const* d, double step, double const* col)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), f, k;
	double* x = ms_positions(m);
	for (f = 0; f < nt; f++)
	{
		unsigned int v[3];
		double lo, hi, L;
		qaws_mesh_3d_get_triangle(v, m, f);
		lo = fmin(fmin(d[v[0]], d[v[1]]), d[v[2]]);
		hi = fmax(fmax(d[v[0]], d[v[1]]), d[v[2]]);
		for (L = ceil(lo / step) * step; L < hi; L += step)
		{
			double p[2][3];
			int n = 0;
			for (k = 0; k < 3 && n < 2; k++)
			{
				double a = d[v[k]], b = d[v[(k + 1) % 3]];
				if ((a - L) * (b - L) < 0)
				{
					double t = (L - a) / (b - a);
					int j;
					for (j = 0; j < 3; j++)
						p[n][j] = x[3 * v[k] + j] + t * (x[3 * v[(k + 1) % 3] + j] - x[3 * v[k] + j]);
					n++;
				}
			}
			if (n == 2)
				er_line(c, p[0], p[1], col, 1.3);
		}
	}
	free(x);
}

/* the mesh coloured by d over [0, top], its isolines, the paths from the given locations, the sources */
static void ms_gd_draw(er_canvas* c, qaws_mesh_3d const* m, qaws_mesh_geodesic const* g, qaws_scalar const* d, double top,
	qaws_mesh_location const* src, unsigned int ns, qaws_mesh_location const* from, unsigned int nfrom)
{
	static double const iso[3] = { 40, 40, 50 }, path[3] = { 20, 20, 20 }, dot[3] = { 250, 220, 40 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), v, i, k;
	double* rgb = (double*)malloc(sizeof(double) * 3 * nv);
	qaws_scalar* pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 4096);
	for (v = 0; v < nv; v++)
		ms_colormap(d[v] < (qaws_scalar)HUGE_VAL ? 0.15 + 0.7 * d[v] / top : 1, rgb + 3 * v);
	ms_draw(c, m, NULL, rgb, 1, NULL, 0);
	ms_gd_isolines(c, m, d, top / 16, iso);
	for (i = 0; i < nfrom; i++)
	{
		unsigned int count = 0;
		if (qaws_mesh_geodesic_path(pts, NULL, &count, g, d, src, ns, from[i], 4096) != QAWS_STATUS_OK)
			continue;
		for (k = 1; k < count; k++)
		{
			double a[3] = { pts[3 * k - 3], pts[3 * k - 2], pts[3 * k - 1] }, b[3] = { pts[3 * k], pts[3 * k + 1], pts[3 * k + 2] };
			er_line(c, a, b, path, 3.2);
		}
	}
	for (i = 0; i < ns; i++)
	{
		unsigned int t[3];
		double p[3];
		qaws_mesh_3d_get_triangle(t, m, src[i].patch);
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, t[0]), b = qaws_mesh_3d_get_vertex(m, t[1]), e = qaws_mesh_3d_get_vertex(m, t[2]);
			double const A[3] = { a.x, a.y, a.z }, B[3] = { b.x, b.y, b.z }, E[3] = { e.x, e.y, e.z };
			p[k] = (1 - src[i].a - src[i].b) * A[k] + src[i].a * B[k] + src[i].b * E[k];
		}
		er_point(c, p, dot, 9);
	}
	free(rgb);
	free(pts);
}

/* the largest finite distance */
static double ms_gd_top(qaws_scalar const* d, unsigned int n)
{
	double top = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
		if (d[i] < (qaws_scalar)HUGE_VAL && d[i] > top)
			top = d[i];
	return top > 0 ? top : 1;
}

/* the geodesic object, the distances from the sources, timed into sub */
static qaws_mesh_geodesic* ms_gd_solve(qaws_mesh_3d const* m, qaws_mesh_location const* src, unsigned int ns, qaws_scalar* d, char* sub)
{
	qaws_mesh_geodesic* g = NULL;
	clock_t t0 = clock(), t1;
	qaws_mesh_geodesic_create(&g, m, NULL);
	t1 = clock();
	if (g)
		qaws_mesh_geodesic_distance(d, g, src, ns);
	sprintf(sub, "%u vertices: factored in %.0f ms, a query %.0f ms", qaws_mesh_3d_get_vertex_count(m), 1000.0 * (double)(t1 - t0) / CLOCKS_PER_SEC,
		1000.0 * (double)(clock() - t1) / CLOCKS_PER_SEC);
	return g;
}

static void demo_geodesic(void)
{
	er_canvas* c = er_create(2250, 720, 3, 252, 252, 250);
	char sub[200];
	if (!c)
		return;
	{
		/* the torus, the source on its outer equator, paths from round it */
		enum { NU = 160, NV = 60 };
		double* x = (double*)malloc(sizeof(double) * 3 * NU * NV);
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * NU * NV), i;
		qaws_mesh_3d* m;
		qaws_mesh_location src, from[9];
		qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * NU * NV);
		qaws_mesh_geodesic* g;
		ms_torus_grid(NU, NV, x, t);
		m = ms_make(x, NU * NV, t, NULL, 2 * NU * NV);
		src.patch = 2 * (0 * NV + 0);
		src.a = (qaws_scalar)0.3;
		src.b = (qaws_scalar)0.3;
		for (i = 0; i < 9; i++)
		{
			unsigned int u = NU / 2 + (i % 3) * NU / 10 - NU / 10, w = (i / 3) * NV / 3 + NV / 6;
			from[i].patch = 2 * (u * NV + w);
			from[i].a = (qaws_scalar)0.3;
			from[i].b = (qaws_scalar)0.3;
		}
		g = ms_gd_solve(m, &src, 1, d, sub);
		ms_panel(c, 0, 0, 750, 720, m, -0.5, 0.95, "a torus: the distance from one point", sub);
		ms_gd_draw(c, m, g, d, ms_gd_top(d, NU * NV), &src, 1, from, 9);
		printf("  torus: %s\n", sub);
		qaws_mesh_geodesic_destroy(g);
		qaws_mesh_3d_destroy(m);
		free(x); free(t); free(d);
	}
	{
		/* a bumpy sphere, three sources */
		unsigned int na = 160, nb = 80, nv = 0, nt = 0, i;
		double const zero[3] = { 0, 0, 0 };
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (2 + na * nb));
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nb);
		qaws_mesh_3d_desc desc;
		qaws_mesh_3d* m = NULL;
		qaws_mesh_location src[3];
		qaws_mesh_geodesic* g;
		qaws_scalar* d;
		ms_ix_sphere(zero, 1, na, nb, x, &nv, t, &nt);
		for (i = 0; i < nv; i++)
		{
			double p[3] = { x[3 * i], x[3 * i + 1], x[3 * i + 2] }, r = 1 + 0.08 * sin(5 * p[0]) * sin(5 * p[1]) * sin(5 * p[2] + 1);
			x[3 * i] = (qaws_scalar)(r * p[0]);
			x[3 * i + 1] = (qaws_scalar)(r * p[1]);
			x[3 * i + 2] = (qaws_scalar)(r * p[2]);
		}
		memset(&desc, 0, sizeof desc);
		desc.positions = x;
		desc.vertex_count = nv;
		desc.face_vertices = t;
		desc.face_count = nt;
		qaws_mesh_3d_create(&m, &desc);
		d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		{
			/* the first source facing the viewer, the others about 70 degrees from it on either side: the triangles nearest
			   those directions */
			double dir[3][3] = { { 0.35, -0.8, -0.45 }, { -0.75, -0.65, 0.1 }, { 0.85, -0.45, 0.3 } };
			for (i = 0; i < 3; i++)
			{
				unsigned int f, best = 0;
				double bd = -2;
				for (f = 0; f < nt; f++)
				{
					double cc[3] = { 0, 0, 0 }, l, s;
					unsigned int k;
					for (k = 0; k < 3; k++)
					{
						cc[0] += x[3 * t[3 * f + k]];
						cc[1] += x[3 * t[3 * f + k] + 1];
						cc[2] += x[3 * t[3 * f + k] + 2];
					}
					l = sqrt(cc[0] * cc[0] + cc[1] * cc[1] + cc[2] * cc[2]);
					s = (cc[0] * dir[i][0] + cc[1] * dir[i][1] + cc[2] * dir[i][2]) / (l * sqrt(dir[i][0] * dir[i][0] + dir[i][1] * dir[i][1] + dir[i][2] * dir[i][2]));
					if (s > bd)
					{
						bd = s;
						best = f;
					}
				}
				src[i].patch = best;
				src[i].a = (qaws_scalar)0.3;
				src[i].b = (qaws_scalar)0.3;
			}
		}
		g = ms_gd_solve(m, src, 3, d, sub);
		ms_panel(c, 750, 0, 750, 720, m, 0.4, 0.35, "a bumpy sphere: the nearest of three", sub);
		ms_gd_draw(c, m, g, d, ms_gd_top(d, nv), src, 3, NULL, 0);
		{
			/* the regions' borders inside the triangles, from each source's exact distances */
			static double const red[3] = { 210, 40, 40 };
			qaws_scalar* ds = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			unsigned int cnt = 0, cap = 8192, k, s, ok = ds != NULL;
			qaws_mesh_location* ends = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * 2 * cap);
			clock_t t0 = clock();
			for (s = 0; s < 3 && ok; s++)
				ok = qaws_mesh_3d_geodesic_exact(ds + (size_t)s * nv, m, &src[s], 1, (qaws_scalar)HUGE_VAL) == QAWS_STATUS_OK;
			if (ok && ends && qaws_mesh_3d_geodesic_borders(ends, NULL, &cnt, m, ds, 3, cap) == QAWS_STATUS_OK)
				for (k = 0; k < cnt; k++)
				{
					double pa[3], pb[3];
					unsigned int z;
					for (z = 0; z < 2; z++)
					{
						unsigned int tv[3], j;
						qaws_mesh_location l = ends[2 * k + z];
						double* q = z ? pb : pa;
						qaws_mesh_3d_get_triangle(tv, m, l.patch);
						for (j = 0; j < 3; j++)
							q[j] = (1 - (double)l.a - (double)l.b) * x[3 * tv[0] + j] + (double)l.a * x[3 * tv[1] + j] + (double)l.b * x[3 * tv[2] + j];
					}
					er_line(c, pa, pb, red, 2.4);
				}
			printf("  bumpy sphere: %u border pieces, exact distances and borders in %.0f ms\n", cnt, 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC);
			free(ds);
			free(ends);
		}
		printf("  bumpy sphere: %s\n", sub);
		qaws_mesh_geodesic_destroy(g);
		qaws_mesh_3d_destroy(m);
		free(x); free(t); free(d);
	}
	{
		/* a sheet with a hill: paths round it */
		enum { N = 120 };
		double* x = (double*)malloc(sizeof(double) * 3 * (N + 1) * (N + 1));
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * N * N), i, j, k = 0;
		qaws_mesh_3d* m;
		qaws_mesh_location src, from[7];
		qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (N + 1) * (N + 1));
		qaws_mesh_geodesic* g;
		for (j = 0; j <= N; j++)
			for (i = 0; i <= N; i++)
			{
				double u = (double)i / N, v = (double)j / N, r2 = (u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5);
				x[3 * (j * (N + 1) + i)] = 2 * u - 1;
				x[3 * (j * (N + 1) + i) + 1] = 2 * v - 1;
				x[3 * (j * (N + 1) + i) + 2] = 0.9 * exp(-r2 / 0.02);
			}
		for (j = 0; j < N; j++)
			for (i = 0; i < N; i++)
			{
				unsigned int a = j * (N + 1) + i;
				t[k++] = a; t[k++] = a + 1; t[k++] = a + N + 2;
				t[k++] = a; t[k++] = a + N + 2; t[k++] = a + N + 1;
			}
		m = ms_make(x, (N + 1) * (N + 1), t, NULL, 2 * N * N);
		src.patch = 2 * ((N / 2) * N + 2);
		src.a = (qaws_scalar)0.3;
		src.b = (qaws_scalar)0.3;
		for (i = 0; i < 7; i++)
		{
			from[i].patch = 2 * ((N / 8 + i * N / 8) * N + N - 3);
			from[i].a = (qaws_scalar)0.3;
			from[i].b = (qaws_scalar)0.3;
		}
		g = ms_gd_solve(m, &src, 1, d, sub);
		ms_panel(c, 1500, 0, 750, 720, m, -0.35, 0.75, "a hill: paths go round it", sub);
		ms_gd_draw(c, m, g, d, ms_gd_top(d, (N + 1) * (N + 1)), &src, 1, from, 7);
		printf("  hill: %s\n", sub);
		qaws_mesh_geodesic_destroy(g);
		qaws_mesh_3d_destroy(m);
		free(x); free(t); free(d);
	}
	er_save(c, "showcase/mesh28_geodesic.png");
	printf("mesh28_geodesic.png\n");
	er_destroy(c);
}
