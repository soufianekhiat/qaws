/*
 * Figure 49: exact polyhedral geodesics (qaws_mesh_3d_geodesic_exact, ICH
 * windows) beside the heat method.
 *   - a cube's surface from a corner: the isolines circles about the corner
 *     on its three faces, unfolded across the edges beyond;
 *   - a sheet with a hill, the source at its left side: the heat method's
 *     distance, then the exact one (the isolines bending round the hill);
 * Printed: the times and the heat method's largest error against the exact.
 */

/* the cube [-1, 1]^3's surface, n x n quads a side cut in two, vertices shared */
static qaws_mesh_3d* ms_gx_cube(unsigned int n)
{
	unsigned int N = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * N * N * N), nv = 0, nf = 0, a, s, u, v, k;
	double* x = (double*)malloc(sizeof(double) * 3 * N * N * N);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
	qaws_mesh_3d* m;
	for (k = 0; k < N * N * N; k++)
	{
		unsigned int q[3];
		q[0] = k % N;
		q[1] = (k / N) % N;
		q[2] = k / (N * N);
		id[k] = 0xFFFFFFFFu;
		if (q[0] % n && q[1] % n && q[2] % n)
			continue;
		id[k] = nv;
		for (a = 0; a < 3; a++) x[3 * nv + a] = 2.0 * q[a] / n - 1;
		nv++;
	}
	for (a = 0; a < 3; a++)
		for (s = 0; s < 2; s++)
			for (u = 0; u < n; u++)
				for (v = 0; v < n; v++)
				{
					unsigned int q4[4], c, b = (a + 1) % 3, e = (a + 2) % 3;
					for (c = 0; c < 4; c++)
					{
						unsigned int q[3];
						q[a] = s * n;
						q[b] = u + (c == 1 || c == 2);
						q[e] = v + (c >= 2);
						q4[c] = id[q[0] + N * (q[1] + N * q[2])];
					}
					if (!s)
					{
						c = q4[1];
						q4[1] = q4[3];
						q4[3] = c;
					}
					f[3 * nf] = q4[0]; f[3 * nf + 1] = q4[1]; f[3 * nf + 2] = q4[2]; nf++;
					f[3 * nf] = q4[0]; f[3 * nf + 1] = q4[2]; f[3 * nf + 2] = q4[3]; nf++;
				}
	m = ms_make(x, nv, f, NULL, nf);
	free(id);
	free(x);
	free(f);
	return m;
}

static void demo_geodesic_exact(void)
{
	er_canvas* c = er_create(2250, 720, 3, 252, 252, 250);
	char sub[200];
	if (!c)
		return;
	{
		/* the cube from its corner (-1, -1, -1), the view onto the far faces' meeting */
		qaws_mesh_3d* m = ms_gx_cube(24);
		unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), f, k;
		qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_mesh_location src;
		clock_t t0;
		src.patch = 0;
		src.a = src.b = 0;
		for (f = 0; f < nt; f++)
		{
			unsigned int v[3];
			qaws_mesh_3d_get_triangle(v, m, f);
			for (k = 0; k < 3; k++)
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(m, v[k]);
				if (p.x == -1 && p.y == -1 && p.z == -1)
				{
					src.patch = f;
					src.a = (qaws_scalar)(k == 1);
					src.b = (qaws_scalar)(k == 2);
				}
			}
		}
		t0 = clock();
		qaws_mesh_3d_geodesic_exact(d, m, &src, 1);
		sprintf(sub, "%u vertices in %.0f ms; far corner %.6f (2 sqrt 5 %.6f)", nv, 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC,
			(double)d[nv - 1], 2 * sqrt(5.0));
		ms_panel(c, 0, 0, 750, 720, m, -0.6, 0.5, "a cube from a corner, exact", sub);
		ms_gd_draw(c, m, NULL, d, ms_gd_top(d, nv), &src, 1, NULL, 0);
		printf("  cube: %s\n", sub);
		qaws_mesh_3d_destroy(m);
		free(d);
	}
	{
		/* a sheet with a hill: the heat method, then exact */
		enum { N = 72 };
		double* x = (double*)malloc(sizeof(double) * 3 * (N + 1) * (N + 1));
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * N * N), i, j, k = 0, nv = (N + 1) * (N + 1);
		qaws_mesh_3d* m;
		qaws_mesh_location src;
		qaws_scalar* dh = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_scalar* de = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_mesh_geodesic* g;
		double top, err = 0, ms;
		clock_t t0;
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
		m = ms_make(x, nv, t, NULL, 2 * N * N);
		src.patch = 2 * ((N / 2) * N + 2);
		src.a = (qaws_scalar)0.3;
		src.b = (qaws_scalar)0.3;
		g = ms_gd_solve(m, &src, 1, dh, sub);
		t0 = clock();
		qaws_mesh_3d_geodesic_exact(de, m, &src, 1);
		ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
		top = ms_gd_top(de, nv);
		for (i = 0; i < nv; i++) err = fmax(err, fabs((double)dh[i] - (double)de[i]));
		ms_panel(c, 750, 0, 750, 720, m, -0.35, 0.75, "a hill: the heat method", sub);
		ms_gd_draw(c, m, g, dh, top, &src, 1, NULL, 0);
		printf("  hill, heat: %s\n", sub);
		sprintf(sub, "exact in %.0f ms; the heat method off by %.3f at most", ms, err);
		ms_panel(c, 1500, 0, 750, 720, m, -0.35, 0.75, "a hill: exact", sub);
		ms_gd_draw(c, m, NULL, de, top, &src, 1, NULL, 0);
		printf("  hill, exact: %s\n", sub);
		qaws_mesh_geodesic_destroy(g);
		qaws_mesh_3d_destroy(m);
		free(x); free(t); free(dh); free(de);
	}
	er_save(c, "showcase/mesh49_geodesic_exact.png");
	printf("mesh49_geodesic_exact.png\n");
	er_destroy(c);
}
