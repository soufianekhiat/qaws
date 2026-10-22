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
 *     nodes, the edges' and faces' on the sphere, the centre at the centre;
 *   - a tetrahedron on a small sphere inside it: its nodes pulled in until
 *     the element folds, relaxed back (its Jacobian kept above the bound);
 *   - no reference: the straight middles, Jacobian 1;
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
	double off = 0;
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
	sprintf(msg, "a cube's hexahedron curved on the sphere through its corners: %u nodes, %u border nodes, edges' and faces' off the sphere %.1e, centre at %.1e, min Jacobian %.3f",
		nn, r.border_nodes, off, vc_norm(p + 3 * en[26]), (double)r.min_jacobian);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && nn == 27 && nen == 27 && r.border_nodes == 18 && off < 1e-3 && vc_norm(p + 3 * en[26]) < 1e-6 &&
		r.min_jacobian > 0.1 && r.relaxed == 0, msg);
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
	test_relax();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
