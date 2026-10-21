/*
 * Figure 20: curved elements (qaws_mesh_generate_nodes_surface).
 *   - a coarse NURBS torus in quads, the elements straight (order 1);
 *   - the same quads at order 2, each drawn through its nine nodes on the
 *     torus (the edge and inner nodes as dots);
 *   - the wavy B-spline patch in cubic triangles (ten nodes each).
 * Elements are drawn finely through their Lagrange maps; the subtitles give
 * the largest distance from those maps to the exact surface.
 */

/* the 1D Lagrange basis on o + 1 equally spaced nodes at t */
static double ms_lagrange(unsigned int o, unsigned int i, double t)
{
	double r = 1;
	unsigned int j;
	for (j = 0; j <= o; j++)
		if (j != i)
			r *= (t * o - j) / ((double)i - j);
	return r;
}

/* the point at (s, t) of element en (k corners, order o; nodes in qaws_mesh_generate_node_counts' order): a quad's tensor
   map over (s, t) along its sides 0-1 and 0-3, a triangle's on its lattice at barycentric (1 - s - t, s, t) */
static void ms_curved_point(double* out, double const* x, unsigned int const* en, unsigned int k, unsigned int o, double s, double t)
{
	unsigned int i, j, c, idx[4][4];
	double l[3];
	out[0] = out[1] = out[2] = 0;
	if (k == 4)
	{
		idx[0][0] = en[0]; idx[o][0] = en[1]; idx[o][o] = en[2]; idx[0][o] = en[3];
		for (j = 1; j < o; j++)
		{
			idx[j][0] = en[4 + j - 1];
			idx[o][j] = en[4 + (o - 1) + j - 1];
			idx[o - j][o] = en[4 + 2 * (o - 1) + j - 1];
			idx[0][o - j] = en[4 + 3 * (o - 1) + j - 1];
		}
		for (j = 1; j < o; j++)
			for (i = 1; i < o; i++)
				idx[i][j] = en[4 * o + (j - 1) * (o - 1) + i - 1];
		for (i = 0; i <= o; i++)
			for (j = 0; j <= o; j++)
			{
				double w = ms_lagrange(o, i, s) * ms_lagrange(o, j, t);
				for (c = 0; c < 3; c++)
					out[c] += w * x[3 * (size_t)idx[i][j] + c];
			}
		return;
	}
	/* a triangle's node at lattice (a, b, d), a + b + d = o: phi = prod over the three of the 1D factors (o l - m) / (m + 1) */
	l[0] = 1 - s - t;
	l[1] = s;
	l[2] = t;
	for (i = 0; i < 3 * o + (o == 3); i++)
	{
		unsigned int lat[3] = { 0, 0, 0 }, m, q;
		double w = 1;
		if (i < 3)
			lat[i] = o;
		else if (i < 3 * o)
		{
			unsigned int e = (i - 3) / (o - 1), jj = (i - 3) % (o - 1) + 1;
			lat[e] = o - jj;
			lat[(e + 1) % 3] = jj;
		}
		else
			lat[0] = lat[1] = lat[2] = 1;
		for (q = 0; q < 3; q++)
			for (m = 0; m < lat[q]; m++)
				w *= (o * l[q] - m) / (m + 1);
		for (c = 0; c < 3; c++)
			out[c] += w * x[3 * (size_t)en[i] + c];
	}
}

/* every element drawn through its map on an n x n grid (smooth normals by differences), its edges as curves,
   its new nodes as dots; returns the largest distance of the grid's points to the surface by `off` (NULL: none) */
static double ms_draw_curved(er_canvas* c, double const* x, unsigned int nn, unsigned int nv, unsigned int const* en, unsigned int const* fs,
	unsigned int nf, unsigned int o, double const* rgb, double const* wire, double const* dot, double (*off)(double const*))
{
	unsigned int f, e = 0, n = 10, i, j;
	double worst = 0;
	for (f = 0; f < nf; e += fs[f] * o + (fs[f] == 4 ? (o - 1) * (o - 1) : fs[f] == 3 && o == 3), f++)
	{
		unsigned int k = fs[f];
		double P[11][11][3], N[11][11][3];
		for (i = 0; i <= n; i++)
			for (j = 0; j <= n; j++)
			{
				double s = (double)i / n, t = (double)j / n, a[3], b[3], du[3], dv[3], h = 1e-4, len;
				if (k == 3 && i + j > n)
					continue;
				ms_curved_point(P[i][j], x, en + e, k, o, s, t);
				ms_curved_point(a, x, en + e, k, o, s + h, t);
				ms_curved_point(b, x, en + e, k, o, s, t + h);
				du[0] = a[0] - P[i][j][0]; du[1] = a[1] - P[i][j][1]; du[2] = a[2] - P[i][j][2];
				dv[0] = b[0] - P[i][j][0]; dv[1] = b[1] - P[i][j][1]; dv[2] = b[2] - P[i][j][2];
				N[i][j][0] = du[1] * dv[2] - du[2] * dv[1];
				N[i][j][1] = du[2] * dv[0] - du[0] * dv[2];
				N[i][j][2] = du[0] * dv[1] - du[1] * dv[0];
				len = sqrt(N[i][j][0] * N[i][j][0] + N[i][j][1] * N[i][j][1] + N[i][j][2] * N[i][j][2]);
				if (len > 0)
					N[i][j][0] /= len, N[i][j][1] /= len, N[i][j][2] /= len;
				if (off)
					worst = fmax(worst, off(P[i][j]));
			}
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (k == 4)
				{
					er_triangle(c, P[i][j], P[i + 1][j], P[i + 1][j + 1], N[i][j], N[i + 1][j], N[i + 1][j + 1], rgb, rgb, rgb);
					er_triangle(c, P[i][j], P[i + 1][j + 1], P[i][j + 1], N[i][j], N[i + 1][j + 1], N[i][j + 1], rgb, rgb, rgb);
				}
				else if (i + j < n)
				{
					er_triangle(c, P[i][j], P[i + 1][j], P[i][j + 1], N[i][j], N[i + 1][j], N[i][j + 1], rgb, rgb, rgb);
					if (i + j + 1 < n)
						er_triangle(c, P[i + 1][j], P[i + 1][j + 1], P[i][j + 1], N[i + 1][j], N[i + 1][j + 1], N[i][j + 1], rgb, rgb, rgb);
				}
		/* the sides: s along 0-1, then 1-2, 2-3 (or 2-0), 3-0 */
		for (i = 0; i < n; i++)
		{
			if (k == 4)
			{
				er_line(c, P[i][0], P[i + 1][0], wire, 1.6);
				er_line(c, P[n][i], P[n][i + 1], wire, 1.6);
				er_line(c, P[i][n], P[i + 1][n], wire, 1.6);
				er_line(c, P[0][i], P[0][i + 1], wire, 1.6);
			}
			else
			{
				er_line(c, P[i][0], P[i + 1][0], wire, 1.6);
				er_line(c, P[n - i][i], P[n - i - 1][i + 1], wire, 1.6);
				er_line(c, P[0][i], P[0][i + 1], wire, 1.6);
			}
		}
	}
	for (i = nv; dot && i < nn; i++)
		er_point(c, x + 3 * (size_t)i, dot, 2.6);
	return worst;
}

static double ms_off_torus(double const* p)
{
	double rho = sqrt(p[0] * p[0] + p[1] * p[1]);
	return fabs(sqrt((rho - 2) * (rho - 2) + p[2] * p[2]) - 0.7);
}

/* a generated mesh's nodes at order o, as doubles; the face sizes and element nodes */
static double* ms_nodes(qaws_mesh_3d const* m, qaws_surface const* s, unsigned int o, unsigned int* nn, unsigned int** en, unsigned int** fs)
{
	unsigned int ne = 0, nf = qaws_mesh_3d_get_face_count(m), i;
	qaws_scalar* p;
	double* x;
	unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(m));
	qaws_mesh_generate_node_counts(nn, &ne, m, o);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * *nn);
	x = (double*)malloc(sizeof(double) * 3 * *nn);
	*en = (unsigned int*)malloc(sizeof(unsigned int) * ne);
	*fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	qaws_mesh_3d_get_faces(fv, *fs, m);
	free(fv);
	if (qaws_mesh_generate_nodes_surface(p, NULL, *en, m, s, o) != QAWS_STATUS_OK)
		printf("  nodes failed\n");
	for (i = 0; i < 3 * *nn; i++)
		x[i] = (double)p[i];
	free(p);
	return x;
}

/* a panel of a generated mesh's elements at order o through its nodes on s */
static void ms_curved_panel(er_canvas* c, int x, qaws_mesh_3d const* m, qaws_surface const* s, unsigned int o, double yaw, double pitch,
	double (*off)(double const*), char const* title)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, red[3] = { 200, 60, 50 };
	unsigned int nn, *en, *fs, nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_face_count(m);
	double* nx = ms_nodes(m, s, o, &nn, &en, &fs);
	double d;
	char sub[200];
	er_panel(c, x, 72, 733, 628, nx, nn, yaw, pitch, 0.06);
	d = ms_draw_curved(c, nx, nn, nv, en, fs, nf, o, tint, wire, o > 1 ? red : NULL, off);
	if (off)
		sprintf(sub, "%u elements, %u nodes, %.1e off the surface", nf, nn, d);
	else
		sprintf(sub, "%u elements, %u nodes", nf, nn);
	printf("  %s: %s\n", title, sub);
	er_text(c, x + 14, 12, title, 3, 30, 30, 30);
	er_text(c, x + 14, 46, sub, 2, 95, 95, 95);
	free(nx); free(en); free(fs);
}

static void demo_curved(void)
{
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	qaws_surface *torus = ms_nurbs_torus(2, 0.7), *wavy = ms_wavy();
	qaws_mesh_3d *tq = NULL, *wt = NULL;
	if (!c || !torus || !wavy)
		return;
	tq = ms_generate(torus, NULL, QAWS_ELEMENTS_QUADS, 0.6, 0, NULL);
	wt = ms_generate(wavy, NULL, QAWS_ELEMENTS_TRIANGLES, 0.3, 0, NULL);
	if (tq)
	{
		ms_curved_panel(c, 0, tq, torus, 1, 0.5, 0.75, ms_off_torus, "a torus in quads, straight");
		ms_curved_panel(c, 733, tq, torus, 2, 0.5, 0.75, ms_off_torus, "the same quads at order 2");
	}
	if (wt)
		ms_curved_panel(c, 1466, wt, wavy, 3, 0.6, 0.9, NULL, "a wavy patch in cubic triangles");
	er_save(c, "showcase/mesh20_curved.png");
	printf("mesh20_curved.png\n");
	qaws_mesh_3d_destroy(tq);
	qaws_mesh_3d_destroy(wt);
	qaws_surface_destroy(torus);
	qaws_surface_destroy(wavy);
	er_destroy(c);
}
