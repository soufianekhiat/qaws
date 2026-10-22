/*
 * Test 174: curved volume elements (qaws_volume_mesh_quadratic_nodes),
 * small meshes
 *
 *   - a sphere's tetrahedra (16 x 8 rings) curved on a PN reference of a
 *     finer sphere: the border's edge nodes' distance to the unit sphere
 *     well under the straight middles' but for a few flat cells with two
 *     border facets (relaxed: their corners' Jacobians hang on border nodes
 *     only), every element's Jacobian above the bound,
 *     the counts (vertices plus edges; 10 per tetrahedron);
 *   - a cube's one hexahedron curved on the sphere through its corners: 27
 *     nodes, the edges' and faces' on the sphere, the centre at the centre,
 *     its Jacobian's Bernstein bound under the sampled determinant;
 *   - a tetrahedron on a small sphere inside it: its nodes pulled in until
 *     the element folds, relaxed back (its Jacobian kept above the bound);
 *   - no reference: the straight middles, Jacobian 1;
 *   - cubic: a cube's hexahedron with a spherical top, its top's nodes on
 *     the cap, its bound under the sampled determinant; two cubes sharing a face labelled differently,
 *     every node at its lattice point; the sphere's tetrahedra;
 *   - untangling: a 7-cell O-grid ball made cubic: corner cells under the
 *     bound brought over it by moving their inside nodes, the border kept;
 *   - refusals: prisms.
 */

#include "test_common.h"
#include "qaws_volume_curve.h"
#include "qaws_tetrahedralize.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define VC_PI 3.14159265358979323846

/* a sphere of radius r: rings between single poles */
static qaws_mesh_3d* vc_sphere(unsigned int nu, unsigned int nv, double r)
{
	unsigned int n = nu * (nv - 1) + 2, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * (nv - 1));
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = VC_PI * j / nv, ph = 2 * VC_PI * i / nu;
			qaws_scalar* p = x + 3 * ((j - 1) * nu + i);
			p[0] = (qaws_scalar)(r * sin(th) * cos(ph)); p[1] = (qaws_scalar)(r * sin(th) * sin(ph)); p[2] = (qaws_scalar)(r * cos(th));
		}
	x[3 * (n - 2)] = 0; x[3 * (n - 2) + 1] = 0; x[3 * (n - 2) + 2] = (qaws_scalar)r;
	x[3 * (n - 1)] = 0; x[3 * (n - 1) + 1] = 0; x[3 * (n - 1) + 2] = (qaws_scalar)-r;
	for (i = 0; i < nu; i++)
	{
		unsigned int a = i, b = (i + 1) % nu, o = (nv - 2) * nu;
		t[k++] = n - 2; t[k++] = a; t[k++] = b;
		t[k++] = n - 1; t[k++] = o + b; t[k++] = o + a;
	}
	for (j = 0; j + 2 < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = b;
			t[k++] = a; t[k++] = a + nu; t[k++] = b + nu;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = n;
	d.face_vertices = t;
	d.face_count = k / 3;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

static qaws_volume_mesh* vc_cells(qaws_scalar const* x, unsigned int nv, unsigned int const* cv, qaws_cell_type const* types, unsigned int nc)
{
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.cell_vertices = cv;
	d.cell_types = types;
	d.cell_count = nc;
	qaws_volume_mesh_create(&m, &d);
	return m;
}

static double vc_norm(qaws_scalar const* p)
{
	return sqrt((double)p[0] * p[0] + (double)p[1] * p[1] + (double)p[2] * p[2]);
}

static void test_sphere(void)
{
	qaws_mesh_3d *s = vc_sphere(16, 8, 1), *fine = vc_sphere(48, 24, 1);
	qaws_mesh_reference* ref = NULL;
	qaws_volume_mesh* v = NULL;
	qaws_tetrahedralize_desc td;
	qaws_volume_curve_report r;
	unsigned int nn = 0, nen = 0, nv, nc, ne;
	qaws_scalar *p = NULL, *q = NULL;
	unsigned int *en = NULL, *en2 = NULL;
	double curved = 0, straight = 0;
	unsigned int on = 0, total = 0;
	qaws_status st = QAWS_STATUS_INVALID_ARGUMENT;
	char msg[400];
	memset(&td, 0, sizeof td);
	memset(&r, 0, sizeof r);
	qaws_volume_mesh_tetrahedralize(&v, NULL, s, &td);
	qaws_mesh_reference_create(&ref, fine, NULL);
	if (v && ref)
	{
		unsigned int i;
		qaws_mesh_3d* b = NULL;
		qaws_volume_mesh_quadratic_node_counts(&nn, &nen, v);
		nv = qaws_volume_mesh_get_vertex_count(v);
		nc = qaws_volume_mesh_get_cell_count(v);
		p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn);
		q = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn);
		en = (unsigned int*)malloc(sizeof(unsigned int) * nen);
		en2 = (unsigned int*)malloc(sizeof(unsigned int) * nen);
		st = qaws_volume_mesh_quadratic_nodes(p, en, &r, v, ref, NULL);
		qaws_volume_mesh_quadratic_nodes(q, en2, NULL, v, NULL, NULL);
		/* the border's edge nodes: those of border facets (the boundary mesh's edges) */
		qaws_volume_mesh_get_boundary(&b, v);
		if (b)
		{
			unsigned int nbe = qaws_mesh_3d_get_edge_count(b), j, att;
			unsigned int* be = (unsigned int*)malloc(sizeof(unsigned int) * 2 * nbe);
			qaws_scalar* vid = (qaws_scalar*)malloc(sizeof(qaws_scalar) * qaws_mesh_3d_get_vertex_count(b));
			qaws_mesh_3d_get_edges(be, b);
			qaws_mesh_3d_find_attribute(&att, b, "vertex");
			qaws_mesh_3d_get_attribute(vid, b, att);
			for (j = 0; j < nbe; j++)
			{
				unsigned int a = (unsigned int)vid[be[2 * j]], c = (unsigned int)vid[be[2 * j + 1]], cc;
				/* the element node of edge (a, c): search the elements */
				for (cc = 0; cc < nc; cc++)
				{
					unsigned int const* e = en + 10 * cc;
					unsigned int k;
					static unsigned char const ed[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
					for (k = 0; k < 6; k++)
						if ((e[ed[k][0]] == a && e[ed[k][1]] == c) || (e[ed[k][0]] == c && e[ed[k][1]] == a))
							break;
					if (k < 6)
					{
						curved = fabs(vc_norm(p + 3 * e[4 + k]) - 1);
						on += curved < 2e-3;
						total++;
						straight = fmax(straight, fabs(vc_norm(q + 3 * e[4 + k]) - 1));
						break;
					}
				}
			}
			free(be);
			free(vid);
			qaws_mesh_3d_destroy(b);
		}
		for (i = 0, ne = 0; i < nen; i++) ne += en[i] != en2[i];
		sprintf(msg, "a sphere's %u tetrahedra curved: %u nodes (%u vertices + edges), %u border nodes; border edge nodes within 2e-3 of the sphere %u of %u (straight: off by up to %.2e), sagitta %.3f, min Jacobian %.3f, relaxed %u; numbering the same without a reference %d",
			nc, nn, nv, r.border_nodes, on, total, straight, (double)r.max_sagitta, (double)r.min_jacobian, r.relaxed, ne == 0);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && nen == 10 * nc && nn > nv && on >= total - total / 10 && straight > 0.03 && r.min_jacobian >= 0.1 && r.relaxed <= nc / 50 &&
			ne == 0, msg);
	}
	else
		TEST_ASSERT(0, "the sphere's tetrahedra and reference");
	free(p);
	free(q);
	free(en);
	free(en2);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(fine);
}

/* a triquadratic hexahedron's Jacobian determinant at (s, t, u) (nodes x: 27, VTK's order) */
static double vc_hex_det(double const (*x)[3], double s, double t, double u)
{
	static signed char const at[27][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { 1, 1, 1 }, { -1, 1, 1 },
		{ 0, -1, -1 }, { 1, 0, -1 }, { 0, 1, -1 }, { -1, 0, -1 }, { 0, -1, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { -1, 0, 1 }, { -1, -1, 0 }, { 1, -1, 0 }, { 1, 1, 0 },
		{ -1, 1, 0 }, { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 }, { 0, 0, 0 } };
	double q[3] = { s, t, u }, v[3][3], dv[3][3], col[3][3];
	unsigned int n, a, d;
	for (a = 0; a < 3; a++)
	{
		v[a][0] = 0.5 * q[a] * (q[a] - 1); v[a][1] = 1 - q[a] * q[a]; v[a][2] = 0.5 * q[a] * (q[a] + 1);
		dv[a][0] = q[a] - 0.5; dv[a][1] = -2 * q[a]; dv[a][2] = q[a] + 0.5;
	}
	memset(col, 0, sizeof col);
	for (n = 0; n < 27; n++)
	{
		unsigned int i = (unsigned int)(at[n][0] + 1), j = (unsigned int)(at[n][1] + 1), k = (unsigned int)(at[n][2] + 1);
		double w[3];
		w[0] = dv[0][i] * v[1][j] * v[2][k];
		w[1] = v[0][i] * dv[1][j] * v[2][k];
		w[2] = v[0][i] * v[1][j] * dv[2][k];
		for (a = 0; a < 3; a++)
			for (d = 0; d < 3; d++) col[a][d] += w[a] * x[n][d];
	}
	return col[0][0] * (col[1][1] * col[2][2] - col[1][2] * col[2][1]) - col[0][1] * (col[1][0] * col[2][2] - col[1][2] * col[2][0]) +
		col[0][2] * (col[1][0] * col[2][1] - col[1][1] * col[2][0]);
}

static void test_hex(void)
{
	static unsigned int const h[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	qaws_scalar x[24], p[27 * 3];
	qaws_cell_type ty = QAWS_CELL_HEXAHEDRON;
	qaws_volume_mesh* v;
	qaws_mesh_3d* sph = vc_sphere(64, 32, sqrt(3.0));
	qaws_mesh_reference* ref = NULL;
	qaws_mesh_reference_desc rd;
	qaws_volume_curve_report r;
	unsigned int en[27], nn = 0, nen = 0, i;
	double off = 0, worst = HUGE_VAL;
	qaws_status st;
	char msg[300];
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)(i == 1 || i == 2 || i == 5 || i == 6 ? 1 : -1);
		x[3 * i + 1] = (qaws_scalar)(i == 2 || i == 3 || i == 6 || i == 7 ? 1 : -1);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? 1 : -1);
	}
	v = vc_cells(x, 8, h, &ty, 1);
	memset(&rd, 0, sizeof rd);
	rd.interpolant = QAWS_INTERPOLANT_PN;
	qaws_mesh_reference_create(&ref, sph, &rd);
	memset(&r, 0, sizeof r);
	qaws_volume_mesh_quadratic_node_counts(&nn, &nen, v);
	st = qaws_volume_mesh_quadratic_nodes(p, en, &r, v, ref, NULL);
	for (i = 8; i < 26; i++) off = fmax(off, fabs(vc_norm(p + 3 * en[i]) - sqrt(3.0)));
	{
		/* the bound under the determinant (the straight cube's is 1): sampled at 11^3 points, never below it, nor far above */
		double xn[27][3];
		unsigned int a, b, c;
		for (i = 0; i < 27; i++)
			for (a = 0; a < 3; a++) xn[i][a] = (double)p[3 * en[i] + a];
		for (a = 0; a <= 10; a++)
			for (b = 0; b <= 10; b++)
				for (c = 0; c <= 10; c++) worst = fmin(worst, vc_hex_det((double const(*)[3])xn, -1 + 0.2 * a, -1 + 0.2 * b, -1 + 0.2 * c));
	}
	sprintf(msg, "a cube's hexahedron curved on the sphere through its corners: %u nodes, %u border nodes, edges' and faces' off the sphere %.1e, centre at %.1e, min Jacobian bounded at %.3f (sampled %.3f)",
		nn, r.border_nodes, off, vc_norm(p + 3 * en[26]), (double)r.min_jacobian, worst);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && nn == 27 && nen == 27 && r.border_nodes == 18 && off < 1e-3 && vc_norm(p + 3 * en[26]) < 1e-6 &&
		r.min_jacobian > 0.1 && r.relaxed == 0 && worst >= (double)r.min_jacobian - 1e-12 && (double)r.min_jacobian > 0.5 * worst, msg);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(sph);
}

/* a quadratic tetrahedron's Jacobian determinant at barycentric l (nodes x: 10, VTK's order) */
static double vc_tet_det(double const (*x)[3], double const* l)
{
	static unsigned char const ed[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
	double g[10][4], col[3][3];
	unsigned int n, k, d;
	memset(g, 0, sizeof g);
	for (k = 0; k < 4; k++) g[k][k] = 4 * l[k] - 1;
	for (n = 0; n < 6; n++)
	{
		g[4 + n][ed[n][0]] = 4 * l[ed[n][1]];
		g[4 + n][ed[n][1]] = 4 * l[ed[n][0]];
	}
	for (k = 0; k < 3; k++)
		for (d = 0; d < 3; d++)
		{
			double s = 0;
			for (n = 0; n < 10; n++) s += x[n][d] * (g[n][k + 1] - g[n][0]);
			col[k][d] = s;
		}
	return col[0][0] * (col[1][1] * col[2][2] - col[1][2] * col[2][1]) - col[0][1] * (col[1][0] * col[2][2] - col[1][2] * col[2][0]) +
		col[0][2] * (col[1][0] * col[2][1] - col[1][1] * col[2][0]);
}

/* the tricubic hexahedron's 64 nodes as documented: (i, j, k) in [0, 3]^3 per slot */
static void vc_hex3_lattice(unsigned char (*lat)[3])
{
	static unsigned char const ed[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
	static unsigned char const fa[6][4] = { { 0, 3, 7, 4 }, { 1, 2, 6, 5 }, { 0, 1, 5, 4 }, { 3, 2, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };
	unsigned char c[8][3];
	unsigned int q, m, a, b, i, j, k, n = 0, d;
	for (q = 0; q < 8; q++)
	{
		c[q][0] = (unsigned char)((q == 1 || q == 2 || q == 5 || q == 6) ? 3 : 0);
		c[q][1] = (unsigned char)((q == 2 || q == 3 || q == 6 || q == 7) ? 3 : 0);
		c[q][2] = (unsigned char)(q >= 4 ? 3 : 0);
		memcpy(lat[n++], c[q], 3);
	}
	for (q = 0; q < 12; q++)
		for (m = 1; m < 3; m++, n++)
			for (d = 0; d < 3; d++) lat[n][d] = (unsigned char)((c[ed[q][0]][d] * (3 - m) + c[ed[q][1]][d] * m) / 3);
	for (q = 0; q < 6; q++)
		for (b = 1; b < 3; b++)
			for (a = 1; a < 3; a++, n++)
				for (d = 0; d < 3; d++)
					lat[n][d] = (unsigned char)((3 * c[fa[q][0]][d] + (c[fa[q][1]][d] - c[fa[q][0]][d]) * (int)a + (c[fa[q][3]][d] - c[fa[q][0]][d]) * (int)b) / 3);
	for (k = 1; k < 3; k++)
		for (j = 1; j < 3; j++)
			for (i = 1; i < 3; i++, n++)
			{
				lat[n][0] = (unsigned char)i;
				lat[n][1] = (unsigned char)j;
				lat[n][2] = (unsigned char)k;
			}
}

/* a tricubic hexahedron's Jacobian determinant at (s, t, u) */
static double vc_hex3_det(double const (*x)[3], unsigned char const (*lat)[3], double s, double t, double u)
{
	double q[3] = { s, t, u }, v[3][4], dv[3][4], col[3][3];
	unsigned int n, a, m, r, w, d;
	for (a = 0; a < 3; a++)
		for (m = 0; m < 4; m++)
		{
			double tm = -1 + 2.0 * m / 3, val = 1, der = 0;
			for (r = 0; r < 4; r++)
				if (r != m)
				{
					double tr = -1 + 2.0 * r / 3, term = 1 / (tm - tr);
					val *= (q[a] - tr) / (tm - tr);
					for (w = 0; w < 4; w++)
						if (w != m && w != r)
							term *= (q[a] - (-1 + 2.0 * w / 3)) / (tm - (-1 + 2.0 * w / 3));
					der += term;
				}
			v[a][m] = val;
			dv[a][m] = der;
		}
	memset(col, 0, sizeof col);
	for (n = 0; n < 64; n++)
	{
		double wt[3];
		wt[0] = dv[0][lat[n][0]] * v[1][lat[n][1]] * v[2][lat[n][2]];
		wt[1] = v[0][lat[n][0]] * dv[1][lat[n][1]] * v[2][lat[n][2]];
		wt[2] = v[0][lat[n][0]] * v[1][lat[n][1]] * dv[2][lat[n][2]];
		for (a = 0; a < 3; a++)
			for (d = 0; d < 3; d++) col[a][d] += wt[a] * x[n][d];
	}
	return col[0][0] * (col[1][1] * col[2][2] - col[1][2] * col[2][1]) - col[0][1] * (col[1][0] * col[2][2] - col[1][2] * col[2][0]) +
		col[0][2] * (col[1][0] * col[2][1] - col[1][1] * col[2][0]);
}

/* the pillow: the cube [-1, 1]^3's surface, n x n quads a side, its top a cap of the sphere of radius 3 through the top
   corners (centre (0, 0, 1 - sqrt 7)), the sides lifted to meet it (z from -1 to the cap, linearly) */
#define VC_CAP_H (1 - 2.6457513110645907)
static double vc_cap(double x, double y)
{
	return VC_CAP_H + sqrt(9 - x * x - y * y);
}

static qaws_mesh_3d* vc_pillow(unsigned int n)
{
	unsigned int N = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * N * N * N), nv = 0, nf = 0, a, s, u, v, k;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * N * N * N);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (k = 0; k < N * N * N; k++)
	{
		unsigned int q[3];
		double c[3];
		q[0] = k % N;
		q[1] = (k / N) % N;
		q[2] = k / (N * N);
		id[k] = 0xFFFFFFFFu;
		if (q[0] % n && q[1] % n && q[2] % n)
			continue;
		for (a = 0; a < 3; a++) c[a] = 2.0 * q[a] / n - 1;
		c[2] = -1 + (c[2] + 1) / 2 * (vc_cap(c[0], c[1]) + 1);
		id[k] = nv;
		for (a = 0; a < 3; a++) x[3 * nv + a] = (qaws_scalar)c[a];
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
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(id);
	free(x);
	free(f);
	return m;
}

/* cubic nodes: a cube's hexahedron with a curved top (the pillow; a hexahedron whose every face lies on one sphere has its
   faces tangent along its edges, its Jacobian 0 there: degenerate), its top nodes on the cap, its sides' on their planes,
   the bound under the sampled determinant (the nodes in the documented order); two hexahedra sharing a face labelled
   differently (every node at its lattice point: the shared nodes agree); the sphere's tetrahedra */
static void test_cubic(void)
{
	static unsigned int const h[8] = { 0, 1, 2, 3, 4, 5, 6, 7 }, h2[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 9, 2, 1, 8, 11, 6, 5, 10 };
	unsigned char lat[64][3];
	qaws_scalar x[24], x2[36], p[64 * 3], *p2 = NULL;
	qaws_cell_type ty[2] = { QAWS_CELL_HEXAHEDRON, QAWS_CELL_HEXAHEDRON };
	qaws_volume_mesh *v = NULL, *v2 = NULL;
	qaws_mesh_3d* sph = vc_pillow(24);
	qaws_mesh_reference* ref = NULL;
	qaws_volume_curve_report r;
	unsigned int en[64], en2[128], nn = 0, nen = 0, nn2 = 0, nen2 = 0, i;
	double off = 0, side = 0, worst = HUGE_VAL, mis = 0;
	qaws_status st;
	char msg[400];
	vc_hex3_lattice(lat);
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)(i == 1 || i == 2 || i == 5 || i == 6 ? 1 : -1);
		x[3 * i + 1] = (qaws_scalar)(i == 2 || i == 3 || i == 6 || i == 7 ? 1 : -1);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? 1 : -1);
	}
	v = vc_cells(x, 8, h, ty, 1);
	qaws_mesh_reference_create(&ref, sph, NULL);
	memset(&r, 0, sizeof r);
	qaws_volume_mesh_cubic_node_counts(&nn, &nen, v);
	st = qaws_volume_mesh_cubic_nodes(p, en, &r, v, ref, NULL);
	for (i = 8; i < 56; i++)
	{
		/* the top face's own nodes on the cap (its edges' lie on the lifted sides already), the sides' on their planes */
		qaws_scalar const* q = p + 3 * en[i];
		if (lat[i][2] == 3 && lat[i][0] % 3 && lat[i][1] % 3)
			off = fmax(off, fabs(sqrt((double)q[0] * q[0] + (double)q[1] * q[1] + ((double)q[2] - VC_CAP_H) * ((double)q[2] - VC_CAP_H)) - 3));
		if (lat[i][0] == 0 || lat[i][0] == 3)
			side = fmax(side, fabs(fabs((double)q[0]) - 1));
		if (lat[i][1] == 0 || lat[i][1] == 3)
			side = fmax(side, fabs(fabs((double)q[1]) - 1));
	}
	{
		double xn[64][3];
		unsigned int a, b, c;
		for (i = 0; i < 64; i++)
			for (a = 0; a < 3; a++) xn[i][a] = (double)p[3 * en[i] + a];
		for (a = 0; a <= 10; a++)
			for (b = 0; b <= 10; b++)
				for (c = 0; c <= 10; c++)
					worst = fmin(worst, vc_hex3_det((double const(*)[3])xn, (unsigned char const(*)[3])lat, -1 + 0.2 * a, -1 + 0.2 * b, -1 + 0.2 * c));
	}
	/* two cubes along x, the second labelled from its own corner: every node at its lattice point without a reference */
	for (i = 0; i < 12; i++)
	{
		x2[3 * i] = (qaws_scalar)(i < 8 ? (i == 1 || i == 2 || i == 5 || i == 6 ? 1 : 0) : 2);
		x2[3 * i + 1] = (qaws_scalar)(i < 8 ? (i == 2 || i == 3 || i == 6 || i == 7 ? 1 : 0) : (i == 9 || i == 11 ? 1 : 0));
		x2[3 * i + 2] = (qaws_scalar)(i < 8 ? (i >= 4 ? 1 : 0) : (i >= 10 ? 1 : 0));
	}
	/* (cell 2: from its corner (2, 1, 0), turned half about z, so the shared face has another frame in it) */
	v2 = vc_cells(x2, 12, h2, ty, 2);
	qaws_volume_mesh_cubic_node_counts(&nn2, &nen2, v2);
	p2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nn2 ? nn2 : 1));
	if (v2 && p2 && nen2 == 128 && qaws_volume_mesh_cubic_nodes(p2, en2, NULL, v2, NULL, NULL) == QAWS_STATUS_OK)
	{
		unsigned int c, n, d, q;
		for (c = 0; c < 2; c++)
			for (n = 0; n < 64; n++)
			{
				/* the trilinear point of the lattice node from the cell's corners */
				double want[3] = { 0, 0, 0 };
				for (q = 0; q < 8; q++)
				{
					double wv = (lat[q][0] ? lat[n][0] / 3.0 : 1 - lat[n][0] / 3.0) * (lat[q][1] ? lat[n][1] / 3.0 : 1 - lat[n][1] / 3.0) *
						(lat[q][2] ? lat[n][2] / 3.0 : 1 - lat[n][2] / 3.0);
					for (d = 0; d < 3; d++) want[d] += wv * (double)x2[3 * h2[8 * c + q] + d];
				}
				for (d = 0; d < 3; d++) mis = fmax(mis, fabs(want[d] - (double)p2[3 * en2[64 * c + n] + d]));
			}
	}
	else
		mis = 1;
	sprintf(msg, "cubic: a cube's hexahedron with a spherical top: %u nodes, %u border nodes, the top face's off the cap %.1e, the sides' off their planes %.1e, min Jacobian bounded at %.3f (sampled %.3f); two cubes sharing a face: %u nodes, every node at its lattice point to %.1e",
		nn, r.border_nodes, off, side, (double)r.min_jacobian, worst, nn2, mis);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && nn == 64 && nen == 64 && r.border_nodes == 48 && off < 1e-3 && side < 1e-9 && r.min_jacobian > 0.1 && r.relaxed == 0 &&
		worst >= (double)r.min_jacobian - 1e-12 && (double)r.min_jacobian > 0.5 * worst && nn2 == 112 && mis < 1e-12, msg);
	{
		/* the sphere's tetrahedra, cubic */
		qaws_mesh_3d *s = vc_sphere(16, 8, 1), *fine = vc_sphere(48, 24, 1);
		qaws_mesh_reference* rf = NULL;
		qaws_volume_mesh* vt = NULL;
		qaws_tetrahedralize_desc td;
		unsigned int tn = 0, ten = 0;
		qaws_scalar* pt = NULL;
		unsigned int* et = NULL;
		memset(&td, 0, sizeof td);
		memset(&r, 0, sizeof r);
		qaws_volume_mesh_tetrahedralize(&vt, NULL, s, &td);
		qaws_mesh_reference_create(&rf, fine, NULL);
		st = QAWS_STATUS_INVALID_ARGUMENT;
		if (vt && rf && qaws_volume_mesh_cubic_node_counts(&tn, &ten, vt) == QAWS_STATUS_OK)
		{
			pt = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * tn);
			et = (unsigned int*)malloc(sizeof(unsigned int) * ten);
			st = qaws_volume_mesh_cubic_nodes(pt, et, &r, vt, rf, NULL);
		}
		sprintf(msg, "cubic: the sphere's %u tetrahedra: %u nodes, %u border nodes, sagitta %.3f, min Jacobian %.3f, relaxed %u",
			vt ? qaws_volume_mesh_get_cell_count(vt) : 0, tn, r.border_nodes, (double)r.max_sagitta, (double)r.min_jacobian, r.relaxed);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && vt && ten == 20 * qaws_volume_mesh_get_cell_count(vt) && r.min_jacobian >= 0.1 &&
			r.relaxed <= qaws_volume_mesh_get_cell_count(vt) / 50, msg);
		free(pt);
		free(et);
		qaws_volume_mesh_destroy(vt);
		qaws_mesh_reference_destroy(rf);
		qaws_mesh_3d_destroy(s);
		qaws_mesh_3d_destroy(fine);
	}
	free(p2);
	qaws_volume_mesh_destroy(v);
	qaws_volume_mesh_destroy(v2);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(sph);
}

/* an O-grid ball of 7 hexahedra: a cube of half-size 0.5 and six blocks out to the unit sphere (each point between the
   cube's face and its projection), each cell turned positive */
static qaws_volume_mesh* vc_ograph(void)
{
	qaws_scalar x[3 * 16];
	unsigned int cv[56], i, k, f, nh = 1;
	qaws_cell_type ty[7];
	for (i = 0; i < 8; i++)
	{
		double c[3];
		c[0] = i == 1 || i == 2 || i == 5 || i == 6 ? 1 : -1;
		c[1] = i == 2 || i == 3 || i == 6 || i == 7 ? 1 : -1;
		c[2] = i >= 4 ? 1 : -1;
		for (k = 0; k < 3; k++)
		{
			x[3 * i + k] = (qaws_scalar)(0.5 * c[k]);
			x[3 * (8 + i) + k] = (qaws_scalar)(c[k] / sqrt(3.0));
		}
		cv[i] = i;
	}
	for (f = 0; f < 6; f++)
	{
		/* the cube's face f (corners on side f) and their projections */
		static unsigned char const fc[6][4] = { { 0, 3, 7, 4 }, { 1, 2, 6, 5 }, { 0, 1, 5, 4 }, { 3, 2, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };
		unsigned int* c8 = cv + 8 * nh, s;
		double e1[3], e3[3], e4[3], dt;
		for (k = 0; k < 4; k++)
		{
			c8[k] = fc[f][k];
			c8[4 + k] = 8 + fc[f][k];
		}
		for (k = 0; k < 3; k++)
		{
			e1[k] = (double)x[3 * c8[1] + k] - x[3 * c8[0] + k];
			e3[k] = (double)x[3 * c8[3] + k] - x[3 * c8[0] + k];
			e4[k] = (double)x[3 * c8[4] + k] - x[3 * c8[0] + k];
		}
		dt = e1[0] * (e3[1] * e4[2] - e3[2] * e4[1]) - e1[1] * (e3[0] * e4[2] - e3[2] * e4[0]) + e1[2] * (e3[0] * e4[1] - e3[1] * e4[0]);
		if (dt < 0)
		{
			s = c8[1]; c8[1] = c8[3]; c8[3] = s;
			s = c8[5]; c8[5] = c8[7]; c8[7] = s;
		}
		nh++;
	}
	for (i = 0; i < 7; i++) ty[i] = QAWS_CELL_HEXAHEDRON;
	return vc_cells(x, 16, cv, ty, 7);
}

/* untangling: the 7-cell O-grid ball, cubic on the sphere */
static void test_untangle(void)
{
	qaws_volume_mesh* v = vc_ograph();
	qaws_mesh_3d* sph = vc_sphere(96, 48, 1);
	qaws_mesh_reference* ref = NULL;
	qaws_volume_curve_report r;
	unsigned int nn = 0, nen = 0;
	qaws_scalar* p = NULL;
	unsigned int* en = NULL;
	qaws_status st = QAWS_STATUS_INVALID_ARGUMENT;
	char msg[300];
	memset(&r, 0, sizeof r);
	qaws_mesh_reference_create(&ref, sph, NULL);
	if (v && ref && qaws_volume_mesh_cubic_node_counts(&nn, &nen, v) == QAWS_STATUS_OK)
	{
		p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn);
		en = (unsigned int*)malloc(sizeof(unsigned int) * nen);
		st = qaws_volume_mesh_cubic_nodes(p, en, &r, v, ref, NULL);
	}
	sprintf(msg, "untangling: a 7-cell O-grid ball, cubic: %u untangled, %u relaxed, Jacobian bound %.3f, border placed %.1e off the reference",
		r.untangled, r.relaxed, (double)r.min_jacobian, (double)r.max_distance);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.untangled >= 1 && r.min_jacobian >= 0.1, msg);
	free(p);
	free(en);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(sph);
}

static void test_relax(void)
{
	static unsigned int const t[4] = { 0, 1, 2, 3 };
	static qaws_scalar const x[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };
	qaws_volume_mesh* v = vc_cells(x, 4, t, NULL, 1);
	qaws_mesh_3d* sph = vc_sphere(32, 16, 0.12);
	qaws_mesh_reference* ref = NULL;
	qaws_volume_curve_report r, r0;
	qaws_scalar p[30], q[30];
	unsigned int en[10], i, k;
	double mid = 0, sampled = 0;
	int bound_ok = 0;
	qaws_status st, s0;
	char msg[300];
	/* the small sphere at the tetrahedron's centroid */
	{
		unsigned int n = qaws_mesh_3d_get_vertex_count(sph);
		qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		qaws_mesh_3d_get_positions(ps, sph);
		for (i = 0; i < 3 * n; i++) ps[i] += (qaws_scalar)0.25;
		qaws_mesh_3d_set_positions(sph, ps);
		free(ps);
	}
	qaws_mesh_reference_create(&ref, sph, NULL);
	memset(&r, 0, sizeof r);
	st = qaws_volume_mesh_quadratic_nodes(p, en, &r, v, ref, NULL);
	s0 = qaws_volume_mesh_quadratic_nodes(q, en, &r0, v, NULL, NULL);
	for (i = 0; i < 6; i++)
	{
		static unsigned char const ed[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
		for (k = 0; k < 3; k++) mid = fmax(mid, fabs((double)q[3 * en[4 + i] + k] - 0.5 * ((double)x[3 * ed[i][0] + k] + (double)x[3 * ed[i][1] + k])));
	}
	{
		/* the bound under the determinant: sampled at 2000 points, never below it */
		double xc[10][3], straight, worst = HUGE_VAL;
		unsigned int s, seed = 7;
		static double const cen[4] = { 0.25, 0.25, 0.25, 0.25 };
		for (i = 0; i < 10; i++)
			for (k = 0; k < 3; k++) xc[i][k] = (double)q[3 * en[i] + k];
		straight = vc_tet_det((double const(*)[3])xc, cen);
		for (i = 0; i < 10; i++)
			for (k = 0; k < 3; k++) xc[i][k] = (double)p[3 * en[i] + k];
		for (s = 0; s < 2000; s++)
		{
			double l[4], tot = 0;
			for (k = 0; k < 4; k++)
			{
				seed = seed * 1664525u + 1013904223u;
				l[k] = -log(((seed >> 8) + 1) / 16777217.0);
				tot += l[k];
			}
			for (k = 0; k < 4; k++) l[k] /= tot;
			worst = fmin(worst, vc_tet_det((double const(*)[3])xc, l) / straight);
		}
		bound_ok = worst >= (double)r.min_jacobian - 1e-12;
		sampled = worst;
	}
	sprintf(msg, "a tetrahedron on a small sphere inside it: relaxed %u, min Jacobian bound %.3f (sampled at 2000 points: %.3f), its nodes left %.3f off the reference (straight %.3f); no reference: middles off %.1e, Jacobian %.3f",
		r.relaxed, (double)r.min_jacobian, sampled, (double)r.max_distance, (double)r.max_sagitta, mid, (double)r0.min_jacobian);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && bound_ok && r.relaxed == 1 && r.min_jacobian >= 0.1 && r.max_distance > 0 && r.max_distance < r.max_sagitta && s0 == QAWS_STATUS_OK &&
		mid < 1e-15 && fabs((double)r0.min_jacobian - 1) < 1e-12, msg);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(sph);
}

static void test_refusals(void)
{
	static unsigned int const pr[6] = { 0, 1, 2, 3, 4, 5 };
	static qaws_scalar const x[18] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1 };
	qaws_cell_type ty = QAWS_CELL_PRISM;
	qaws_volume_mesh* v = vc_cells(x, 6, pr, &ty, 1);
	unsigned int nn, nen;
	qaws_status st = qaws_volume_mesh_quadratic_node_counts(&nn, &nen, v);
	char msg[200];
	sprintf(msg, "refused: a prism (%d)", (int)st);
	printf("  %s\n", msg);
	TEST_ASSERT(v && st == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	qaws_volume_mesh_destroy(v);
}

int test_174_volume_curve_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 174: curved volume elements\n");
	test_sphere();
	test_hex();
	test_cubic();
	test_untangle();
	test_relax();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
