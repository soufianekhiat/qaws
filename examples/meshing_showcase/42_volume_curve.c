/*
 * Figure 42: curved volume elements (qaws_volume_mesh_quadratic_nodes),
 * cut open, each quadratic face drawn finely.
 *   - a sphere's tetrahedra, straight (their quadratic nodes at the
 *     middles);
 *   - the same curved on a PN reference of a finer sphere;
 *   - a ball of hexahedra (an O-grid: a cube and six blocks out to the
 *     sphere, one border face per border cell), straight;
 *   - the same curved: triquadratic, the border on the sphere.
 */

static unsigned char const ms_vc_tf[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
static unsigned char const ms_vc_te[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
static unsigned char const ms_vc_hf[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
static unsigned char const ms_vc_he[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
static unsigned char const ms_vc_hvf[6][4] = { { 0, 3, 7, 4 }, { 1, 2, 6, 5 }, { 0, 1, 5, 4 }, { 3, 2, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };

static unsigned int ms_vc_slot(int tet, unsigned int i, unsigned int j)
{
	unsigned int e;
	for (e = 0; e < (tet ? 6u : 12u); e++)
	{
		unsigned int a = tet ? ms_vc_te[e][0] : ms_vc_he[e][0], b = tet ? ms_vc_te[e][1] : ms_vc_he[e][1];
		if ((a == i && b == j) || (a == j && b == i))
			return (tet ? 4 : 8) + e;
	}
	return 0;
}

/* a quadratic face drawn as n x n pieces: a triangle from its corners and edges' nodes (p[0..5]: corners, then edges 01 12 20),
   or a quad (p[0..8]: corners, edges 01 12 23 30, centre) */
static void ms_vc_face(er_canvas* c, double const* const* p, int quad, unsigned int n, double const* rgb, double const* wire)
{
	unsigned int i, j, k;
	double g[13][13][3];
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			double u = (double)i / n, v = (double)j / n;
			for (k = 0; k < 3; k++)
				if (quad)
				{
					double s = 2 * u - 1, t = 2 * v - 1, a[3] = { 0.5 * s * (s - 1), 1 - s * s, 0.5 * s * (s + 1) }, b[3] = { 0.5 * t * (t - 1), 1 - t * t, 0.5 * t * (t + 1) };
					/* (i, j) in {0, 1, 2}^2 to the node */
					static unsigned char const at[3][3] = { { 0, 4, 1 }, { 7, 8, 5 }, { 3, 6, 2 } };
					unsigned int x, y;
					double val = 0;
					for (y = 0; y < 3; y++)
						for (x = 0; x < 3; x++) val += a[x] * b[y] * p[at[y][x]][k];
					g[j][i][k] = val;
				}
				else if (i + j <= n)
				{
					double l1 = u, l2 = v, l0 = 1 - u - v;
					g[j][i][k] = p[0][k] * l0 * (2 * l0 - 1) + p[1][k] * l1 * (2 * l1 - 1) + p[2][k] * l2 * (2 * l2 - 1) + 4 * (l0 * l1 * p[3][k] + l1 * l2 * p[4][k] + l2 * l0 * p[5][k]);
				}
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			if (quad || i + j + 1 <= n)
				er_triangle(c, g[j][i], g[j][i + 1], g[j + 1][i], NULL, NULL, NULL, rgb, rgb, rgb);
			if (quad || i + j + 2 <= n)
				er_triangle(c, g[j][i + 1], g[j + 1][i + 1], g[j + 1][i], NULL, NULL, NULL, rgb, rgb, rgb);
		}
	/* the face's sides */
	for (i = 0; i < n; i++)
	{
		er_line(c, g[0][i], g[0][i + 1], wire, 0.7);
		if (quad)
		{
			er_line(c, g[n][i], g[n][i + 1], wire, 0.7);
			er_line(c, g[i][0], g[i + 1][0], wire, 0.7);
			er_line(c, g[i][n], g[i + 1][n], wire, 0.7);
		}
		else
		{
			er_line(c, g[i][0], g[i + 1][0], wire, 0.7);
			er_line(c, g[i][n - i], g[i + 1][n - i - 1], wire, 0.7);
		}
	}
}

/* the cells with their centroid at y > 0, their faces not shared with another kept cell drawn (border ones in rgb, cut ones in cut) */
static void ms_vc_draw(er_canvas* c, qaws_volume_mesh const* v, qaws_scalar const* pos, unsigned int const* en, double const* rgb, double const* cut,
	double const* wire)
{
	unsigned int nc = qaws_volume_mesh_get_cell_count(v), nf = qaws_volume_mesh_get_facet_count(v), cc, f, k, q = 0, fq = 0;
	qaws_cell_type* ty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * nc);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 8 * nc);
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	unsigned int *e0 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *f0 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
	unsigned char* kept = (unsigned char*)malloc(nc);
	double* x = (double*)malloc(sizeof(double) * 3 * 27 * nc);
	qaws_volume_mesh_get_cells(cv, ty, v);
	qaws_volume_mesh_get_cell_adjacency(adj, v);
	for (cc = 0; cc < nc; cc++)
	{
		int tet = ty[cc] == QAWS_CELL_TETRAHEDRON;
		double y = 0;
		e0[cc] = q;
		f0[cc] = fq;
		for (k = 0; k < (tet ? 4u : 8u); k++) y += (double)pos[3 * en[q + k] + 1];
		kept[cc] = (unsigned char)(y > 0);
		q += tet ? 10 : 27;
		fq += tet ? 4 : 6;
	}
	e0[nc] = q;
	f0[nc] = fq;
	for (k = 0; k < 3 * q; k++) x[k] = 0;
	for (cc = 0; cc < nc; cc++)
	{
		int tet = ty[cc] == QAWS_CELL_TETRAHEDRON;
		unsigned int const* e = en + e0[cc];
		if (!kept[cc])
			continue;
		for (f = 0; f < (tet ? 4u : 6u); f++)
		{
			unsigned int o = adj[f0[cc] + f], m = tet ? 3 : 4, i;
			double pt[9][3];
			double const* pp[9];
			if (o != QAWS_MESH_NONE && kept[o])
				continue;
			for (i = 0; i < m; i++)
			{
				unsigned int a = tet ? ms_vc_tf[f][i] : ms_vc_hf[f][i], b = tet ? ms_vc_tf[f][(i + 1) % m] : ms_vc_hf[f][(i + 1) % m];
				for (k = 0; k < 3; k++)
				{
					pt[i][k] = (double)pos[3 * e[a] + k];
					pt[m + i][k] = (double)pos[3 * e[ms_vc_slot(tet, a, b)] + k];
				}
			}
			if (!tet)
			{
				/* the face's centre node: the VTK face with the same corners */
				unsigned int s, hit = 20;
				for (s = 0; s < 6; s++)
				{
					unsigned int z, w, all = 1;
					for (z = 0; z < 4; z++)
					{
						int found = 0;
						for (w = 0; w < 4; w++) found |= ms_vc_hvf[s][w] == ms_vc_hf[f][z];
						all &= (unsigned int)found;
					}
					if (all)
						hit = 20 + s;
				}
				for (k = 0; k < 3; k++) pt[8][k] = (double)pos[3 * e[hit] + k];
			}
			for (i = 0; i < 9; i++) pp[i] = pt[i];
			ms_vc_face(c, pp, !tet, tet ? 6 : 8, o == QAWS_MESH_NONE ? rgb : cut, wire);
		}
	}
	free(ty); free(cv); free(adj); free(e0); free(f0); free(kept); free(x);
}

/* the index of point x among the n kept (appended if new) */
static unsigned int ms_vc_point(double* pts, unsigned int* n, double const* x)
{
	unsigned int i;
	for (i = 0; i < *n; i++)
		if (fabs(pts[3 * i] - x[0]) + fabs(pts[3 * i + 1] - x[1]) + fabs(pts[3 * i + 2] - x[2]) < 1e-9)
			return i;
	memcpy(pts + 3 * *n, x, sizeof(double) * 3);
	return (*n)++;
}

/* A ball of hexahedra as an O-grid: an inner cube of half-size a (n^3 cells) and six blocks of n x n x m cells out to the
   unit sphere, each point between the cube's face and its projection; every border cell has one face on the sphere. */
static qaws_volume_mesh* ms_vc_ball(unsigned int n, unsigned int m, double a)
{
	unsigned int cap = 8 * (n * n * n + 6 * n * n * m), nh = 0, np = 0, f, i, j, k, l;
	double* pts = (double*)malloc(sizeof(double) * 3 * cap);
	unsigned int* h = (unsigned int*)malloc(sizeof(unsigned int) * cap);
	qaws_cell_type* ty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * cap / 8);
	qaws_scalar* xs;
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* v = NULL;
	for (k = 0; k < n; k++)
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int z;
				for (z = 0; z < 8; z++)
				{
					unsigned int di = z == 1 || z == 2 || z == 5 || z == 6, dj = z == 2 || z == 3 || z == 6 || z == 7, dk = z >= 4;
					double x[3] = { a * (2.0 * (i + di) / n - 1), a * (2.0 * (j + dj) / n - 1), a * (2.0 * (k + dk) / n - 1) };
					h[8 * nh + z] = ms_vc_point(pts, &np, x);
				}
				nh++;
			}
	for (f = 0; f < 6; f++)
	{
		unsigned int ax = f / 2, u = (ax + 1) % 3, w = (ax + 2) % 3;
		double sg = f % 2 ? 1 : -1;
		for (l = 0; l < m; l++)
			for (j = 0; j < n; j++)
				for (i = 0; i < n; i++)
				{
					unsigned int z;
					for (z = 0; z < 8; z++)
					{
						unsigned int di = z == 1 || z == 2 || z == 5 || z == 6, dj = z == 2 || z == 3 || z == 6 || z == 7, dl = z >= 4, q;
						double c[3], s[3], x[3], len, t = (double)(l + dl) / m;
						c[ax] = sg;
						c[u] = 2.0 * (i + di) / n - 1;
						c[w] = 2.0 * (j + dj) / n - 1;
						len = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
						for (q = 0; q < 3; q++)
						{
							s[q] = c[q] / len;
							x[q] = (1 - t) * a * c[q] + t * s[q];
						}
						h[8 * nh + z] = ms_vc_point(pts, &np, x);
					}
					nh++;
				}
	}
	/* each cell turned positive (its base mirrored when not) */
	for (i = 0; i < nh; i++)
	{
		unsigned int* c8 = h + 8 * i, s;
		double e1[3], e3[3], e4[3], dt;
		for (k = 0; k < 3; k++)
		{
			e1[k] = pts[3 * c8[1] + k] - pts[3 * c8[0] + k];
			e3[k] = pts[3 * c8[3] + k] - pts[3 * c8[0] + k];
			e4[k] = pts[3 * c8[4] + k] - pts[3 * c8[0] + k];
		}
		dt = e1[0] * (e3[1] * e4[2] - e3[2] * e4[1]) - e1[1] * (e3[0] * e4[2] - e3[2] * e4[0]) + e1[2] * (e3[0] * e4[1] - e3[1] * e4[0]);
		if (dt < 0)
		{
			s = c8[1]; c8[1] = c8[3]; c8[3] = s;
			s = c8[5]; c8[5] = c8[7]; c8[7] = s;
		}
		ty[i] = QAWS_CELL_HEXAHEDRON;
	}
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * np);
	for (i = 0; i < 3 * np; i++) xs[i] = (qaws_scalar)pts[i];
	memset(&d, 0, sizeof d);
	d.positions = xs;
	d.vertex_count = np;
	d.cell_vertices = h;
	d.cell_types = ty;
	d.cell_count = nh;
	qaws_volume_mesh_create(&v, &d);
	free(pts);
	free(h);
	free(ty);
	free(xs);
	return v;
}

/* a mesh's quadratic nodes on the reference (NULL: straight); drawn in a panel */
static void ms_vc_panel(er_canvas* c, int px, qaws_volume_mesh const* v, qaws_mesh_reference const* ref, char const* title, double const* rgb,
	double const* cut, double const* wire)
{
	unsigned int nn = 0, nen = 0, i;
	qaws_scalar* pos;
	unsigned int* en;
	double* pd;
	qaws_volume_curve_report r;
	char sub[200];
	if (qaws_volume_mesh_quadratic_node_counts(&nn, &nen, v) != QAWS_STATUS_OK)
		return;
	pos = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn);
	en = (unsigned int*)malloc(sizeof(unsigned int) * nen);
	pd = (double*)malloc(sizeof(double) * 3 * nn);
	memset(&r, 0, sizeof r);
	qaws_volume_mesh_quadratic_nodes(pos, en, &r, v, ref, NULL);
	for (i = 0; i < 3 * nn; i++) pd[i] = (double)pos[i];
	if (ref)
		sprintf(sub, "%u cells, Jacobian >= %.2f, %u relaxed", qaws_volume_mesh_get_cell_count(v), (double)r.min_jacobian, r.relaxed);
	else
		sprintf(sub, "%u cells, sagitta to the sphere", qaws_volume_mesh_get_cell_count(v));
	er_panel(c, px, 72, 450, 628, pd, nn, 0.55, 0.45, 0.06);
	er_text(c, px + 14, 12, title, 3, 30, 30, 30);
	er_text(c, px + 14, 46, sub, 2, 95, 95, 95);
	ms_vc_draw(c, v, pos, en, rgb, cut, wire);
	printf("  %s: %s, sagitta %.4f\n", title, sub, (double)r.max_sagitta);
	free(pos);
	free(en);
	free(pd);
}

static void demo_volume_curve(void)
{
	static double const wire[3] = { 60, 66, 80 }, skin[3] = { 120, 150, 210 }, cut[3] = { 232, 196, 120 }, skin2[3] = { 130, 190, 140 };
	er_canvas* c = er_create(1800, 700, 3, 252, 252, 250);
	qaws_mesh_3d *s = ms_rm_sphere(16, 8), *fine = ms_rm_sphere(64, 32);
	qaws_mesh_reference* ref = NULL;
	qaws_volume_mesh *t = NULL, *b = ms_vc_ball(3, 2, 0.45);
	qaws_tetrahedralize_desc td;
	if (!c || !s || !fine || !b)
		return;
	memset(&td, 0, sizeof td);
	qaws_volume_mesh_tetrahedralize(&t, NULL, s, &td);
	qaws_mesh_reference_create(&ref, fine, NULL);
	if (t && ref)
	{
		ms_vc_panel(c, 0, t, NULL, "straight", skin, cut, wire);
		ms_vc_panel(c, 450, t, ref, "curved P2", skin, cut, wire);
	}
	if (ref)
	{
		ms_vc_panel(c, 900, b, NULL, "a ball of hexahedra", skin2, cut, wire);
		ms_vc_panel(c, 1350, b, ref, "curved Q2", skin2, cut, wire);
	}
	qaws_volume_mesh_destroy(t);
	qaws_volume_mesh_destroy(b);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(fine);
	er_save(c, "showcase/mesh42_volume_curve.png");
	printf("mesh42_volume_curve.png\n");
	er_destroy(c);
}
