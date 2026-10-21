/*
 * Test 157: curved elements (qaws_mesh_generate_nodes_surface/_shell)
 *
 *   - a NURBS torus in triangles and in quads at orders 1, 2, 3: node
 *     counts (vertices, order - 1 per edge, a quad's (order - 1)^2 inside,
 *     a cubic triangle's centroid), every node on the torus, each edge's
 *     nodes shared by its two faces; the curved elements' distance to the
 *     torus between their nodes (edges by Lagrange interpolation, quads
 *     over their tensor grid) a fortieth of the straight elements' or
 *     less at both orders (Lagrange interpolation's radial error is O(h^4)
 *     at order 2 and 3 alike);
 *   - a box of six bilinear faces as a shell, Q2: every node on its
 *     element's face (an inner node located on it), the nodes of edges
 *     along the box's edges on them;
 *   - refusals: order 0 and 4, a shell mesh without its "face" attribute.
 */

#include "test_common.h"
#include "qaws_mesh_generate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* the torus of radii R, r: the unit circle's 9 rational points times the profile circle's */
static qaws_surface* tc_torus(double R, double r)
{
	static double const cx[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, cy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	static qaws_scalar const kn[12] = { 0, 0, 0, 0.25f, 0.25f, 0.5f, 0.5f, 0.75f, 0.75f, 1, 1, 1 };
	qaws_vec3 cp[81];
	qaws_scalar w[81];
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 9; i++)
		for (j = 0; j < 9; j++)
		{
			double rad = R + r * cx[j];
			cp[i * 9 + j].x = (qaws_scalar)(rad * cx[i]);
			cp[i * 9 + j].y = (qaws_scalar)(rad * cy[i]);
			cp[i * 9 + j].z = (qaws_scalar)(r * cy[j]);
			w[i * 9 + j] = (qaws_scalar)((i & 1 ? sqrt(0.5) : 1) * (j & 1 ? sqrt(0.5) : 1));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cp;
	d.u_point_count = 9;
	d.v_point_count = 9;
	d.weights = w;
	d.u_knots = kn;
	d.u_knot_count = 12;
	d.v_knots = kn;
	d.v_knot_count = 12;
	qaws_surface_create_nurbs(&s, &d);
	return s;
}

static double tc_off_torus(double const* p)
{
	double rho = sqrt(p[0] * p[0] + p[1] * p[1]);
	return fabs(sqrt((rho - 2) * (rho - 2) + p[2] * p[2]) - 0.7);
}

/* the 1D Lagrange basis on o + 1 equally spaced nodes at t */
static void tc_lagrange(double* out, unsigned int o, double t)
{
	unsigned int i, j;
	for (i = 0; i <= o; i++)
	{
		out[i] = 1;
		for (j = 0; j <= o; j++)
			if (j != i)
				out[i] *= (t * o - j) / ((double)i - j);
	}
}

/* a generated mesh's nodes at order o */
typedef struct tc_nodes
{
	unsigned int nv, nf, nn, ne;
	unsigned int *fs, *en;
	qaws_scalar* x;
	qaws_mesh_location* loc;
	qaws_status st;
} tc_nodes;

static void tc_make(tc_nodes* t, qaws_mesh_3d const* m, qaws_surface const* one, qaws_shell_atlas const* shell, unsigned int o)
{
	memset(t, 0, sizeof *t);
	t->nv = qaws_mesh_3d_get_vertex_count(m);
	t->nf = qaws_mesh_3d_get_face_count(m);
	t->st = qaws_mesh_generate_node_counts(&t->nn, &t->ne, m, o);
	if (t->st != QAWS_STATUS_OK)
		return;
	t->fs = (unsigned int*)malloc(sizeof(unsigned int) * (t->nf + 1));
	t->en = (unsigned int*)malloc(sizeof(unsigned int) * (t->ne + 1));
	t->x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (t->nn + 1));
	t->loc = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (t->nn + 1));
	{
		unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1));
		qaws_mesh_3d_get_faces(fv, t->fs, m);
		free(fv);
	}
	t->st = one ? qaws_mesh_generate_nodes_surface(t->x, t->loc, t->en, m, one, o) : qaws_mesh_generate_nodes_shell(t->x, t->loc, t->en, m, shell, o);
}

static void tc_free(tc_nodes* t)
{
	free(t->fs); free(t->en); free(t->x); free(t->loc);
}

static void tc_point(double* out, tc_nodes const* t, unsigned int n)
{
	out[0] = (double)t->x[3 * (size_t)n];
	out[1] = (double)t->x[3 * (size_t)n + 1];
	out[2] = (double)t->x[3 * (size_t)n + 2];
}

/* the torus element by element at order o: node count, nodes off the torus, edges whose two faces disagree
   on their nodes, the elements' distance to the torus between nodes */
static void tc_check_torus(qaws_mesh_3d const* m, qaws_surface const* s, unsigned int o, int quads, double* out_between, int* out_ok, char* msg)
{
	tc_nodes t;
	unsigned int f, k, j, e = 0, mismatched = 0, ne = 0, nfe = 0, expect;
	double off = 0, between = 0;
	qaws_mesh_topology top;
	tc_make(&t, m, s, NULL, o);
	qaws_mesh_3d_get_topology(&top, m);
	*out_ok = 0;
	if (t.st != QAWS_STATUS_OK)
	{
		sprintf(msg, "torus order %u: %s", o, qaws_status_to_string(t.st));
		tc_free(&t);
		return;
	}
	for (f = 0; f < t.nf; f++)
		nfe += t.fs[f];
	ne = nfe / 2;   /* closed */
	expect = t.nv + ne * (o - 1) + (quads ? t.nf * (o - 1) * (o - 1) : o == 3 ? t.nf : 0);
	for (k = 0; k < t.nn; k++)
	{
		double p[3];
		tc_point(p, &t, k);
		off = fmax(off, tc_off_torus(p));
	}
	/* each edge's nodes against its twin's, reversed: searched by its end vertices */
	{
		unsigned int* fo = (unsigned int*)malloc(sizeof(unsigned int) * (t.nf + 1));
		for (f = 0; f < t.nf; f++)
		{
			fo[f] = e;
			e += t.fs[f] * o + (t.fs[f] == 4 ? (o - 1) * (o - 1) : t.fs[f] == 3 && o == 3 ? 1 : 0);
		}
		for (f = 0; f < t.nf && o > 1; f++)
			for (k = 0; k < t.fs[f]; k++)
			{
				unsigned int a = t.en[fo[f] + k], b = t.en[fo[f] + (k + 1) % t.fs[f]], g, l, found = 0;
				for (g = 0; g < t.nf && !found; g++)
					for (l = 0; l < t.fs[g] && !found; l++)
						if (t.en[fo[g] + l] == b && t.en[fo[g] + (l + 1) % t.fs[g]] == a)
						{
							found = 1;
							for (j = 1; j < o; j++)
								mismatched += t.en[fo[f] + t.fs[f] + k * (o - 1) + j - 1] != t.en[fo[g] + t.fs[g] + l * (o - 1) + (o - 1 - j)];
						}
				mismatched += !found;
			}
		/* between the nodes: each edge's curve at the midpoints of its node intervals; a quad's tensor grid at its cells' centres */
		for (f = 0; f < t.nf; f++)
		{
			unsigned int n = t.fs[f];
			unsigned int const* en = t.en + fo[f];
			for (k = 0; k < n; k++)
			{
				unsigned int line[4];
				line[0] = en[k];
				line[o] = en[(k + 1) % n];
				for (j = 1; j < o; j++)
					line[j] = en[n + k * (o - 1) + j - 1];
				for (j = 0; j < o; j++)
				{
					double L[4], p[3] = { 0, 0, 0 }, q[3];
					unsigned int i, c;
					tc_lagrange(L, o, (j + 0.5) / o);
					for (i = 0; i <= o; i++)
					{
						tc_point(q, &t, line[i]);
						for (c = 0; c < 3; c++)
							p[c] += L[i] * q[c];
					}
					between = fmax(between, tc_off_torus(p));
				}
			}
			if (n == 4)
			{
				unsigned int G[4][4], i, a, b, c;
				/* G[i][j]: the node at (i/o, j/o) */
				G[0][0] = en[0]; G[o][0] = en[1]; G[o][o] = en[2]; G[0][o] = en[3];
				for (j = 1; j < o; j++)
				{
					G[j][0] = en[4 + 0 * (o - 1) + j - 1];
					G[o][j] = en[4 + 1 * (o - 1) + j - 1];
					G[o - j][o] = en[4 + 2 * (o - 1) + j - 1];
					G[0][o - j] = en[4 + 3 * (o - 1) + j - 1];
				}
				for (j = 1; j < o; j++)
					for (i = 1; i < o; i++)
						G[i][j] = en[4 * o + (j - 1) * (o - 1) + (i - 1)];
				for (a = 0; a < o; a++)
					for (b = 0; b < o; b++)
					{
						double Ls[4], Lt[4], p[3] = { 0, 0, 0 }, q[3];
						tc_lagrange(Ls, o, (a + 0.5) / o);
						tc_lagrange(Lt, o, (b + 0.5) / o);
						for (i = 0; i <= o; i++)
							for (j = 0; j <= o; j++)
							{
								tc_point(q, &t, G[i][j]);
								for (c = 0; c < 3; c++)
									p[c] += Ls[i] * Lt[j] * q[c];
							}
						between = fmax(between, tc_off_torus(p));
					}
			}
		}
		free(fo);
	}
	*out_between = between;
	sprintf(msg, "torus %s order %u: %u faces, %u nodes (expected %u), %u element nodes, %u edge nodes mismatched, nodes off the torus %.2e, "
		"between nodes %.3e", quads ? "quads" : "triangles", o, t.nf, t.nn, expect, t.ne, mismatched, off, between);
	*out_ok = top.closed && t.nn == expect && e == t.ne && mismatched == 0 && off < TEST_BOUND(1e-9, 32);
	tc_free(&t);
}

static void test_torus(void)
{
	qaws_surface* s = tc_torus(2, 0.7);
	qaws_generate_desc d;
	unsigned int el;
	char msg[400];
	for (el = 0; el < 2; el++)
	{
		qaws_mesh_3d* m = NULL;
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.elements = el ? QAWS_ELEMENTS_QUADS : QAWS_ELEMENTS_TRIANGLES;
		d.target_length = (qaws_scalar)0.4;
		st = qaws_mesh_generate_surface(&m, NULL, s, &d);
		if (st != QAWS_STATUS_OK || !m)
			TEST_ASSERT(0, "torus: generated");
		else
		{
			double between[4];
			unsigned int o;
			for (o = 1; o <= 3; o++)
			{
				int ok;
				tc_check_torus(m, s, o, (int)el, &between[o], &ok, msg);
				TEST_ASSERT(ok, msg);
			}
			sprintf(msg, "torus %s: distance between nodes %.3e (order 1), %.3e (2), %.3e (3)", el ? "quads" : "triangles", between[1], between[2], between[3]);
			TEST_ASSERT(between[2] < between[1] / 40 && between[3] < between[1] / 40, msg);
		}
		qaws_mesh_3d_destroy(m);
	}
	qaws_surface_destroy(s);
}

/* the box [-h0, h0] x [-h1, h1] x [-h2, h2]: face 2k at -h_k, 2k + 1 at +h_k, in the order bottom, top, front, back, left, right */
static int const tc_axis[6] = { 2, 2, 1, 1, 0, 0 };
static double const tc_side[6] = { -1, 1, -1, 1, -1, 1 };

static void tc_box(qaws_surface** f, double const* h)
{
	static int const q[6][4][3] = {
		{ { -1, -1, -1 }, { -1, 1, -1 }, { 1, -1, -1 }, { 1, 1, -1 } }, { { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { 1, -1, -1 }, { -1, -1, 1 }, { 1, -1, 1 } }, { { -1, 1, -1 }, { -1, 1, 1 }, { 1, 1, -1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { -1, -1, 1 }, { -1, 1, -1 }, { -1, 1, 1 } }, { { 1, -1, -1 }, { 1, 1, -1 }, { 1, -1, 1 }, { 1, 1, 1 } } };
	unsigned int i, k;
	for (i = 0; i < 6; i++)
	{
		qaws_surface_bilinear_desc d;
		qaws_vec3* p[4];
		p[0] = &d.p00; p[1] = &d.p10; p[2] = &d.p01; p[3] = &d.p11;
		for (k = 0; k < 4; k++)
		{
			p[k]->x = (qaws_scalar)(h[0] * q[i][k][0]);
			p[k]->y = (qaws_scalar)(h[1] * q[i][k][1]);
			p[k]->z = (qaws_scalar)(h[2] * q[i][k][2]);
		}
		f[i] = NULL;
		qaws_surface_create_bilinear(&f[i], &d);
	}
}

/* coordinates of p at the box's extremes */
static unsigned int tc_on(double const* p, double const* h)
{
	unsigned int k, on = 0;
	for (k = 0; k < 3; k++)
		on += h[k] - fabs(p[k]) < TEST_BOUND(1e-9, 32);
	return on;
}

static void test_box(void)
{
	static double const h[3] = { 1, 0.7, 0.5 };
	qaws_surface* box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_generate_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_status st;
	char msg[400];
	unsigned int i;
	tc_box(box, h);
	st = qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	memset(&d, 0, sizeof d);
	d.elements = QAWS_ELEMENTS_QUADS;
	d.target_length = (qaws_scalar)0.25;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_generate_shell(&m, NULL, atlas, &d);
	if (st != QAWS_STATUS_OK || !m)
		TEST_ASSERT(0, "box: generated");
	else
	{
		tc_nodes t;
		unsigned int fa = QAWS_MESH_NONE, f, k, e = 0, offface = 0, edge_nodes = 0, off_edge = 0, badloc = 0, skipped = 0;
		qaws_scalar* fid = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (qaws_mesh_3d_get_face_count(m) + 1));
		tc_make(&t, m, NULL, atlas, 2);
		qaws_mesh_3d_find_attribute(&fa, m, "face");
		qaws_mesh_3d_get_attribute(fid, m, fa);
		for (f = 0; f < t.nf && t.st == QAWS_STATUS_OK; f++)
		{
			unsigned int n = t.fs[f], F = (unsigned int)fid[f];
			unsigned int const* en = t.en + e;
			unsigned int straddle = 0;
			/* an element whose corners are not all on its face (the float build's quads can straddle an edge) is the mesh's, skipped */
			for (k = 0; k < n; k++)
			{
				double p[3];
				tc_point(p, &t, en[k]);
				straddle += F >= 6 || fabs(p[tc_axis[F]] - tc_side[F] * h[tc_axis[F]]) > TEST_BOUND(1e-9, 32);
			}
			e += n * 2 + (n == 4);
			if (straddle)
			{
				skipped++;
				continue;
			}
			for (k = n; k < n * 2 + (n == 4); k++)
			{
				double p[3];
				tc_point(p, &t, en[k]);
				offface += fabs(p[tc_axis[F]] - tc_side[F] * h[tc_axis[F]]) > TEST_BOUND(1e-9, 32);
				badloc += k >= 2 * n && t.loc[en[k]].patch != F;
			}
			/* an edge along the box's edges (both ends on two extremes): its node on it too */
			for (k = 0; k < n; k++)
			{
				double a[3], b[3], p[3];
				tc_point(a, &t, en[k]);
				tc_point(b, &t, en[(k + 1) % n]);
				tc_point(p, &t, en[n + k]);
				if (tc_on(a, h) >= 2 && tc_on(b, h) >= 2)
				{
					unsigned int c, same = 0;
					for (c = 0; c < 3; c++)
						same += fabs(a[c] - b[c]) < 1e-9;
					if (same == 2)
					{
						edge_nodes++;
						off_edge += tc_on(p, h) < 2;
					}
				}
			}
		}
		sprintf(msg, "box Q2: %s, %u faces, %u nodes, %u elements skipped (corners off their face), %u new nodes off their face, %u inner nodes located on another face, "
			"%u nodes of edges along the box's edges, %u off them", qaws_status_to_string(t.st), t.nf, t.nn, skipped, offface, badloc, edge_nodes, off_edge);
		TEST_ASSERT(t.st == QAWS_STATUS_OK && t.nn > t.nv && skipped <= t.nf / 50 && offface == 0 && badloc == 0 && edge_nodes > 40 && off_edge == 0, msg);
		free(fid);
		tc_free(&t);
	}
	qaws_mesh_3d_destroy(m);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
}

static void test_refusals(void)
{
	static qaws_scalar const x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	static unsigned int const fv[3] = { 0, 1, 2 };
	qaws_surface* box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_3d_desc md;
	qaws_scalar out[30];
	unsigned int en[30], nn = 0, i;
	static double const h[3] = { 1, 1, 1 };
	qaws_status s0, s4, s5, s1;
	char msg[300];
	tc_box(box, h);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	memset(&md, 0, sizeof md);
	md.positions = x;
	md.vertex_count = 3;
	md.face_vertices = fv;
	md.face_count = 1;
	qaws_mesh_3d_create(&m, &md);
	s0 = qaws_mesh_generate_node_counts(&nn, NULL, m, 0);
	s4 = qaws_mesh_generate_nodes_surface(out, NULL, en, m, box[0], 4);
	s5 = qaws_mesh_generate_nodes_shell(out, NULL, en, m, atlas, 2);
	s1 = qaws_mesh_generate_nodes_surface(out, NULL, en, m, box[0], 1);
	sprintf(msg, "refusals: order 0 %s, order 4 %s, shell mesh without \"face\" %s; order 1 on a plain mesh %s", qaws_status_to_string(s0),
		qaws_status_to_string(s4), qaws_status_to_string(s5), qaws_status_to_string(s1));
	TEST_ASSERT(s0 == QAWS_STATUS_INVALID_ARGUMENT && s4 == QAWS_STATUS_INVALID_ARGUMENT && s5 == QAWS_STATUS_INVALID_ARGUMENT
		&& s1 == QAWS_STATUS_OK && en[0] == 0 && en[2] == 2, msg);
	qaws_mesh_3d_destroy(m);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
}

int test_157_curved_elements_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 157: curved elements\n");
	test_torus();
	test_box();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
