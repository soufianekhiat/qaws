/*
 * Figure 30: deformation and parameterization (qaws_mesh_3d_deform,
 * qaws_mesh_3d_parameterize).
 *   - a bar held at its foot, its top ring turned a quarter and carried
 *     over: bent as rigidly as possible (ARAP), and by Laplacian editing
 *     (its rings shear instead of turning);
 *   - a sheet with a hill flattened as rigidly as possible: its (u, v) grid
 *     drawn on it, and the flat layout.
 * Printed: the distortions.
 */

/* lines of constant u (or v) every step on the triangles, from the corners' (u, v) */
static void ms_df_uvlines(er_canvas* c, qaws_mesh_3d const* m, qaws_scalar const* uv, unsigned int comp, double step, double const* col)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), f, k;
	for (f = 0; f < nt; f++)
	{
		unsigned int v[3];
		double x[3][3], d[3], lo, hi, L;
		qaws_mesh_3d_get_triangle(v, m, f);
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, v[k]);
			x[k][0] = p.x;
			x[k][1] = p.y;
			x[k][2] = p.z;
			d[k] = uv[2 * (3 * f + k) + comp];
		}
		lo = fmin(fmin(d[0], d[1]), d[2]);
		hi = fmax(fmax(d[0], d[1]), d[2]);
		for (L = ceil(lo / step) * step; L < hi; L += step)
		{
			double p[2][3];
			int n = 0, j;
			for (k = 0; k < 3 && n < 2; k++)
				if ((d[k] - L) * (d[(k + 1) % 3] - L) < 0)
				{
					double t = (L - d[k]) / (d[(k + 1) % 3] - d[k]);
					for (j = 0; j < 3; j++)
						p[n][j] = x[k][j] + t * (x[(k + 1) % 3][j] - x[k][j]);
					n++;
				}
			if (n == 2)
				er_line(c, p[0], p[1], col, 1.4);
		}
	}
}

/* a bar of radius r and height h, nu round and nv along, open at both ends */
static qaws_mesh_3d* ms_df_bar(unsigned int nu, unsigned int nv)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * (nv + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j <= nv; j++)
		for (i = 0; i < nu; i++)
		{
			double a = 2 * MS_PI * i / nu;
			x[3 * (j * nu + i)] = 0.15 * cos(a);
			x[3 * (j * nu + i) + 1] = 0.15 * sin(a);
			x[3 * (j * nu + i) + 2] = 1.6 * j / nv;
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b; t[k++] = b + nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = a + nu;
		}
	m = ms_make(x, nu * (nv + 1), t, NULL, 2 * nu * nv);
	free(x);
	free(t);
	return m;
}

static void demo_deform(void)
{
	static double const col[3] = { 180, 190, 214 }, warm[3] = { 226, 196, 150 }, wire[3] = { 70, 75, 90 }, red[3] = { 200, 50, 40 };
	static double const blue[3] = { 40, 80, 200 };
	er_canvas* c = er_create(2400, 720, 3, 252, 252, 250);
	char sub[200];
	unsigned int k;
	if (!c)
		return;
	for (k = 0; k < 2; k++)
	{
		enum { NU = 24, NV = 40 };
		qaws_mesh_3d* m = ms_df_bar(NU, NV);
		unsigned int h[2 * NU], i;
		qaws_scalar tg[6 * NU];
		qaws_deform_desc d;
		double r = 1.6 / (MS_PI / 2);
		for (i = 0; i < NU; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i), q = qaws_mesh_3d_get_vertex(m, NV * NU + i);
			h[i] = i;
			tg[3 * i] = p.x;
			tg[3 * i + 1] = p.y;
			tg[3 * i + 2] = p.z;
			/* the top ring turned a quarter about y and carried to the end of a quarter circle */
			h[NU + i] = NV * NU + i;
			tg[3 * (NU + i)] = (qaws_scalar)(r + (q.z - 1.6));
			tg[3 * (NU + i) + 1] = q.y;
			tg[3 * (NU + i) + 2] = (qaws_scalar)(r - q.x);
		}
		memset(&d, 0, sizeof d);
		d.method = k ? QAWS_DEFORM_LAPLACIAN : QAWS_DEFORM_ARAP;
		d.iterations = 0;
		qaws_mesh_3d_deform(m, h, tg, 2 * NU, &d);
		sprintf(sub, "%s", k ? "the rings shear instead of turning" : "local-global rounds: the rings turn with it");
		ms_panel(c, 600 * (int)k, 0, 600, 720, m, 0.0, 0.25, k ? "Laplacian editing" : "ARAP: a bar bent", sub);
		ms_draw(c, m, col, NULL, 1, wire, 0.4);
		qaws_mesh_3d_destroy(m);
	}
	{
		/* a sheet with a hill, flattened (ARAP) */
		enum { N = 30 };
		double* x = (double*)malloc(sizeof(double) * 3 * (N + 1) * (N + 1));
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * N * N), i, j, kk = 0, nt;
		qaws_mesh_3d *m, *flat;
		qaws_scalar* uv;
		qaws_param_desc d;
		qaws_param_report r;
		for (j = 0; j <= N; j++)
			for (i = 0; i <= N; i++)
			{
				double u = 2.0 * i / N - 1, v = 2.0 * j / N - 1;
				x[3 * (j * (N + 1) + i)] = u;
				x[3 * (j * (N + 1) + i) + 1] = v;
				x[3 * (j * (N + 1) + i) + 2] = 0.7 * exp(-(u * u + v * v) / 0.18);
			}
		for (j = 0; j < N; j++)
			for (i = 0; i < N; i++)
			{
				unsigned int a = j * (N + 1) + i;
				t[kk++] = a; t[kk++] = a + 1; t[kk++] = a + N + 2;
				t[kk++] = a; t[kk++] = a + N + 2; t[kk++] = a + N + 1;
			}
		m = ms_make(x, (N + 1) * (N + 1), t, NULL, 2 * N * N);
		nt = qaws_mesh_3d_get_triangle_count(m);
		uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nt);
		memset(&d, 0, sizeof d);
		d.method = QAWS_PARAM_ARAP;
		qaws_mesh_3d_parameterize(uv, &r, m, NULL, 0, &d);
		sprintf(sub, "angles %.3f rad, areas %.3f, %u turned", (double)r.angle_distortion, (double)r.area_distortion, r.flipped);
		ms_panel(c, 1200, 0, 600, 720, m, -0.6, 0.7, "a hill flattened (ARAP)", sub);
		ms_draw(c, m, warm, NULL, 1, NULL, 0);
		ms_df_uvlines(c, m, uv, 0, 0.1, red);
		ms_df_uvlines(c, m, uv, 1, 0.1, blue);
		printf("  hill, ARAP: %s\n", sub);
		{
			/* the layout: a soup of the triangles at their (u, v) */
			double* fx = (double*)malloc(sizeof(double) * 9 * nt);
			unsigned int* ft = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
			for (i = 0; i < 3 * nt; i++)
			{
				fx[3 * i] = uv[2 * i];
				fx[3 * i + 1] = uv[2 * i + 1];
				fx[3 * i + 2] = 0;
				ft[i] = i;
			}
			flat = ms_make(fx, 3 * nt, ft, NULL, nt);
			ms_panel(c, 1800, 0, 600, 720, flat, 0.0, 1.5707, "the layout in the plane", "each triangle near its own shape");
			ms_draw(c, flat, warm, NULL, 0, wire, 0.3);
			qaws_mesh_3d_destroy(flat);
			free(fx);
			free(ft);
		}
		qaws_mesh_3d_destroy(m);
		free(x); free(t); free(uv);
	}
	er_save(c, "showcase/mesh30_deform.png");
	printf("mesh30_deform.png\n");
	er_destroy(c);
}
