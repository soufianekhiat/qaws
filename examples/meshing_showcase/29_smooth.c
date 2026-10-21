/*
 * Figure 29: smoothing and fairing (qaws_mesh_3d_smooth, qaws_mesh_3d_fair).
 *   - a noisy torus, then smoothed by 30 Laplacian steps (it shrinks), 30
 *     Taubin steps (it does not), one bi-Laplacian step;
 *   - a noisy roof smoothed by one long implicit step without features (the crease rounds) and with
 *     the crease as a feature line (it stays sharp);
 *   - a sphere's cap thrown about, then faired as a membrane and as a thin
 *     plate (the sphere comes back).
 * Printed: the volumes.
 */

static double ms_sm_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		v += ((double)a.x * (b.y * e.z - b.z * e.y) + (double)a.y * (b.z * e.x - b.x * e.z) + (double)a.z * (b.x * e.y - b.y * e.x)) / 6;
	}
	return v;
}

/* a torus with noise along its normals */
static qaws_mesh_3d* ms_sm_torus(void)
{
	enum { NU = 96, NV = 36 };
	double* x = (double*)malloc(sizeof(double) * 3 * NU * NV);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * NU * NV), i, j, seed = 77;
	qaws_mesh_3d* m;
	ms_torus_grid(NU, NV, x, t);
	for (i = 0; i < NU; i++)
		for (j = 0; j < NV; j++)
		{
			double u = 2 * MS_PI * i / NU, v = 2 * MS_PI * j / NV, a = 0.035 * (ms_uniform(&seed) - 0.5);
			double* p = x + 3 * (i * NV + j);
			p[0] += a * cos(v) * cos(u);
			p[1] += a * cos(v) * sin(u);
			p[2] += a * sin(v);
		}
	m = ms_make(x, NU * NV, t, NULL, 2 * NU * NV);
	free(x);
	free(t);
	return m;
}

/* a roof: two planes meeting along a crease at x = 0, noise in z */
static qaws_mesh_3d* ms_sm_roof(int feature)
{
	enum { N = 40 };
	double* x = (double*)malloc(sizeof(double) * 3 * (N + 1) * (N + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * N * N), i, j, k = 0, seed = 5;
	unsigned int e[2 * N];
	qaws_mesh_3d* m;
	for (j = 0; j <= N; j++)
		for (i = 0; i <= N; i++)
		{
			double u = 2.0 * i / N - 1, v = 2.0 * j / N - 1;
			double* p = x + 3 * (j * (N + 1) + i);
			p[0] = u;
			p[1] = v;
			p[2] = 0.6 - 0.6 * fabs(u) + (i && j && i < N && j < N ? 0.03 * (ms_uniform(&seed) - 0.5) : 0);
		}
	for (j = 0; j < N; j++)
		for (i = 0; i < N; i++)
		{
			unsigned int v = j * (N + 1) + i;
			int fl = i < N / 2;
			t[k++] = v; t[k++] = v + 1; t[k++] = fl ? v + N + 1 : v + N + 2;
			t[k++] = fl ? v + 1 : v; t[k++] = v + N + 2; t[k++] = v + N + 1;
		}
	m = ms_make(x, (N + 1) * (N + 1), t, NULL, 2 * N * N);
	for (j = 0; j < N; j++)
	{
		e[2 * j] = j * (N + 1) + N / 2;
		e[2 * j + 1] = (j + 1) * (N + 1) + N / 2;
	}
	if (feature)
		qaws_mesh_3d_set_feature_edges(m, e, N);
	free(x);
	free(t);
	return m;
}

static void demo_smooth(void)
{
	static double const col[3] = { 180, 190, 214 }, warm[3] = { 226, 196, 150 }, wire[3] = { 70, 75, 90 }, red[3] = { 200, 40, 40 };
	static char const* const tname[4] = { "a noisy torus", "30 Laplacian steps", "30 Taubin steps", "one bi-Laplacian step" };
	er_canvas* c = er_create(2400, 1400, 3, 252, 252, 250);
	char sub[200];
	unsigned int k;
	double v0 = 0;
	if (!c)
		return;
	for (k = 0; k < 4; k++)
	{
		qaws_mesh_3d* m = ms_sm_torus();
		qaws_smooth_desc d;
		memset(&d, 0, sizeof d);
		d.method = k == 3 ? QAWS_SMOOTH_BILAPLACIAN : k == 2 ? QAWS_SMOOTH_TAUBIN : QAWS_SMOOTH_LAPLACIAN;
		d.iterations = k == 3 ? 1 : 30;
		if (k)
			qaws_mesh_3d_smooth(m, &d);
		if (!k)
			v0 = ms_sm_volume(m);
		sprintf(sub, "volume %.4f (%+.1f%%)", ms_sm_volume(m), 100 * (ms_sm_volume(m) / v0 - 1));
		/* (the same frame for all four: fitted to the noisy one) */
		{
			qaws_mesh_3d* frame = ms_sm_torus();
			ms_panel(c, 600 * (int)k, 0, 600, 700, frame, -0.4, 0.8, tname[k], sub);
			qaws_mesh_3d_destroy(frame);
		}
		ms_draw(c, m, col, NULL, 0, NULL, 0);
		printf("  torus, %s: %s\n", tname[k], sub);
		qaws_mesh_3d_destroy(m);
	}
	for (k = 0; k < 2; k++)
	{
		qaws_mesh_3d* m = ms_sm_roof((int)k);
		qaws_smooth_desc d;
		unsigned int e[80], ne, i;
		memset(&d, 0, sizeof d);
		d.method = QAWS_SMOOTH_IMPLICIT;
		d.step = 40;
		qaws_mesh_3d_smooth(m, &d);
		sprintf(sub, "implicit, t = 40 h^2%s", k ? ", along the crease only" : "");
		ms_panel(c, 600 * (int)k, 700, 600, 700, m, 0.0, 0.3, k ? "the crease a feature" : "a noisy roof, no features", sub);
		ms_draw(c, m, warm, NULL, 1, wire, 0.25);
		ne = qaws_mesh_3d_get_feature_edge_count(m);
		if (ne && ne <= 40)
		{
			qaws_mesh_3d_get_feature_edges(e, m);
			for (i = 0; i < ne; i++)
			{
				qaws_vec3 a = qaws_mesh_3d_get_vertex(m, e[2 * i]), b = qaws_mesh_3d_get_vertex(m, e[2 * i + 1]);
				double A[3] = { a.x, a.y, a.z }, B[3] = { b.x, b.y, b.z };
				er_line(c, A, B, red, 3);
			}
		}
		qaws_mesh_3d_destroy(m);
	}
	for (k = 0; k < 2; k++)
	{
		/* a sphere, its cap thrown about, faired */
		unsigned int na = 64, nb = 32, nv = 0, nt = 0, i, seed = 3;
		double const zero[3] = { 0, 0, 0 };
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (2 + na * nb));
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nb);
		unsigned char* fr;
		qaws_mesh_3d_desc desc;
		qaws_mesh_3d* m = NULL;
		qaws_fair_desc f;
		double* rgb;
		ms_ix_sphere(zero, 1, na, nb, x, &nv, t, &nt);
		fr = (unsigned char*)calloc(nv, 1);
		rgb = (double*)malloc(sizeof(double) * 3 * nv);
		for (i = 0; i < nv; i++)
		{
			fr[i] = x[3 * i + 1] < -0.5;
			if (fr[i])
			{
				x[3 * i] += (qaws_scalar)(0.15 * (ms_uniform(&seed) - 0.5));
				x[3 * i + 1] += (qaws_scalar)(0.15 * (ms_uniform(&seed) - 0.5));
				x[3 * i + 2] += (qaws_scalar)(0.15 * (ms_uniform(&seed) - 0.5));
			}
			memcpy(rgb + 3 * i, fr[i] ? warm : col, sizeof col);
		}
		memset(&desc, 0, sizeof desc);
		desc.positions = x;
		desc.vertex_count = nv;
		desc.face_vertices = t;
		desc.face_count = nt;
		qaws_mesh_3d_create(&m, &desc);
		memset(&f, 0, sizeof f);
		f.order = k ? 2 : 1;
		qaws_mesh_3d_fair(m, fr, &f);
		sprintf(sub, "the warm cap's %s", k ? "thin plate: tangent to the sphere round it" : "membrane: the least area");
		ms_panel(c, 1200 + 600 * (int)k, 700, 600, 700, m, 0.3, -0.45, k ? "a cap faired, order 2" : "a cap faired, order 1", sub);
		ms_draw(c, m, NULL, rgb, 1, wire, 0.2);
		qaws_mesh_3d_destroy(m);
		free(x); free(t); free(fr); free(rgb);
	}
	er_save(c, "showcase/mesh29_smooth.png");
	printf("mesh29_smooth.png\n");
	er_destroy(c);
}
