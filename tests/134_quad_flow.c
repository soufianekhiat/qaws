/*
 * Test 134: quads by min-cost flow (QuadriFlow's consistency)
 *
 *   - the flow module: random small networks against exhaustive
 *     enumeration (the cost the least there is), an infeasible supply
 *     reported, two runs alike, an augmentation reaching the nearer of two
 *     targets and then, its arc full, the other;
 *   - a planar grid of steps with a pair of opposite dislocations five
 *     faces apart: closed by five unit changes; one dislocation three faces
 *     from the border: out through it in three;
 *   - two cones off lattice points (their sums odd) made even;
 *   - the blob quad-dominant: no face but quads; the blob, the can and a
 *     torus all quads at the target spacing, their irregular vertices
 *     counted; a wavy B-spline patch generated in quads;
 *   - the same quads twice; consistency timed on a dense blob.
 */

#include "test_common.h"
#include "qaws_quadrangulate.h"
#include "qaws_remesh.h"
#include "qaws_mesh_generate.h"
#include "mesh/qaws_quad_flow.h"
#include "internal/qaws_internal_flow.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

static unsigned int tf_rand(unsigned int* s)
{
	*s = *s * 1664525u + 1013904223u;
	return *s >> 8;
}

/* ---- the flow module ---------------------------------------------------- */

static void test_flow_module(void)
{
	unsigned int seed = 134, g, worse = 0, unrouted = 0, broken = 0, differ = 0;
	char msg[300];
	for (g = 0; g < 50; g++)
	{
		unsigned int n = 3 + tf_rand(&seed) % 4, m = 1 + tf_rand(&seed) % 8, i, a, from[8], to[8];
		int cap[8], b[6] = { 0, 0, 0, 0, 0, 0 }, x[8], fl[8], left = 0, fl2[8];
		long long cost[8], best = -1, got = 0, total = 1;
		qaws_flow *f = NULL, *f2 = NULL;
		for (a = 0; a < m; a++)
		{
			from[a] = tf_rand(&seed) % n;
			to[a] = (from[a] + 1 + tf_rand(&seed) % (n - 1)) % n;
			cap[a] = 1 + (int)(tf_rand(&seed) % 2);
			cost[a] = tf_rand(&seed) % 10;
			x[a] = (int)(tf_rand(&seed) % (unsigned int)(cap[a] + 1));   /* (a feasible flow gives the supplies) */
			b[from[a]] += x[a];
			b[to[a]] -= x[a];
			total *= cap[a] + 1;
		}
		/* every flow within the capacities that meets the supplies */
		for (i = 0; i < (unsigned int)total; i++)
		{
			unsigned int r = i, v;
			int net[6] = { 0, 0, 0, 0, 0, 0 }, ok = 1;
			long long c = 0;
			for (a = 0; a < m; a++)
			{
				x[a] = (int)(r % (unsigned int)(cap[a] + 1));
				r /= (unsigned int)(cap[a] + 1);
				net[from[a]] += x[a];
				net[to[a]] -= x[a];
				c += x[a] * cost[a];
			}
			for (v = 0; v < n; v++)
				ok &= net[v] == b[v];
			if (ok && (best < 0 || c < best))
				best = c;
		}
		qaws_flow_create(n, m, &f);
		qaws_flow_create(n, m, &f2);
		for (a = 0; a < m; a++)
		{
			qaws_flow_arc(f, from[a], to[a], cap[a], cost[a]);
			qaws_flow_arc(f2, from[a], to[a], cap[a], cost[a]);
		}
		for (i = 0; i < n; i++)
		{
			qaws_flow_supply(f, i, b[i]);
			qaws_flow_supply(f2, i, b[i]);
		}
		qaws_flow_solve(f, &left, &got);
		qaws_flow_solve(f2, NULL, NULL);
		{
			int net[6] = { 0, 0, 0, 0, 0, 0 };
			for (a = 0; a < m; a++)
			{
				fl[a] = qaws_flow_get(f, a);
				fl2[a] = qaws_flow_get(f2, a);
				differ += fl[a] != fl2[a];
				broken += fl[a] < 0 || fl[a] > cap[a];
				net[from[a]] += fl[a];
				net[to[a]] -= fl[a];
			}
			for (i = 0; i < n; i++)
				broken += net[i] != b[i];
		}
		worse += got != best;
		unrouted += left != 0;
		qaws_flow_destroy(f);
		qaws_flow_destroy(f2);
	}
	sprintf(msg, "flow: 50 random networks, %u not the least cost, %u not all routed, %u breaking a capacity or a supply, %u arcs differing between runs",
		worse, unrouted, broken, differ);
	printf("  %s\n", msg);
	TEST_ASSERT(worse == 0 && unrouted == 0 && broken == 0 && differ == 0, msg);
	{
		/* more supply than room: reported, not failed */
		qaws_flow* f = NULL;
		int left = 0;
		qaws_flow_create(2, 1, &f);
		qaws_flow_arc(f, 0, 1, 1, 3);
		qaws_flow_supply(f, 0, 2);
		qaws_flow_supply(f, 1, -2);
		qaws_flow_solve(f, &left, NULL);
		sprintf(msg, "flow: 2 units over an arc of 1, %d unrouted", left);
		TEST_ASSERT(left == 1 && qaws_flow_get(f, 0) == 1, msg);
		qaws_flow_destroy(f);
	}
	{
		/* 0 -> 1 -> 2 costs 5, 0 -> 3 -> 4 costs 3: the nearer first, then (its arcs full) the other */
		qaws_flow* f = NULL;
		unsigned int targets[2] = { 2, 4 }, r1, r2, r3;
		long long c1, c2;
		qaws_flow_create(5, 4, &f);
		qaws_flow_arc(f, 0, 1, 1, 4);
		qaws_flow_arc(f, 1, 2, 1, 1);
		qaws_flow_arc(f, 0, 3, 1, 2);
		qaws_flow_arc(f, 3, 4, 1, 1);
		qaws_flow_augment(f, 0, targets, 2, &r1, &c1);
		qaws_flow_augment(f, 0, targets, 2, &r2, &c2);
		qaws_flow_augment(f, 0, targets, 2, &r3, NULL);
		sprintf(msg, "flow: augmenting reaches %u at %lld, then %u at %lld, then none (%u)", r1, c1, r2, c2, r3);
		TEST_ASSERT(r1 == 4 && c1 == 3 && r2 == 2 && c2 == 5 && r3 == QAWS_FLOW_NONE, msg);
		qaws_flow_destroy(f);
	}
}

/* ---- synthetic steps on a planar grid ------------------------------------ */

#define TG_N 20   /* squares per side */

typedef struct tg_grid
{
	qaws_hm m;
	mf_level L;
	int *T, *T0;
	unsigned char* K;
	double* Tc;
	double x[3 * (TG_N + 1) * (TG_N + 1)], nrm[3 * (TG_N + 1) * (TG_N + 1)], q[3 * (TG_N + 1) * (TG_N + 1)];
} tg_grid;

static unsigned int tg_v(unsigned int i, unsigned int j)
{
	return j * (TG_N + 1) + i;
}

/* square (i, j): A = (i, j) (i + 1, j) (i + 1, j + 1), B = (i, j) (i + 1, j + 1) (i, j + 1) */
static unsigned int tg_face(unsigned int i, unsigned int j, int b)
{
	return 2 * (j * TG_N + i) + (unsigned int)b;
}

static void tg_build(tg_grid* g)
{
	unsigned int fv[6 * TG_N * TG_N], i, j, k = 0, c, nv = (TG_N + 1) * (TG_N + 1);
	for (j = 0; j <= TG_N; j++)
		for (i = 0; i <= TG_N; i++)
		{
			unsigned int v = tg_v(i, j);
			g->x[3 * v] = i;
			g->x[3 * v + 1] = j;
			g->x[3 * v + 2] = 0;
			g->nrm[3 * v] = 0; g->nrm[3 * v + 1] = 0; g->nrm[3 * v + 2] = 1;
			g->q[3 * v] = 1; g->q[3 * v + 1] = 0; g->q[3 * v + 2] = 0;
		}
	for (j = 0; j < TG_N; j++)
		for (i = 0; i < TG_N; i++)
		{
			fv[k++] = tg_v(i, j); fv[k++] = tg_v(i + 1, j); fv[k++] = tg_v(i + 1, j + 1);
			fv[k++] = tg_v(i, j); fv[k++] = tg_v(i + 1, j + 1); fv[k++] = tg_v(i, j + 1);
		}
	hm_build(&g->m, 3, g->x, nv, fv, NULL, 2 * TG_N * TG_N);
	memset(&g->L, 0, sizeof g->L);
	g->L.n = nv;
	g->L.x = g->x;
	g->L.nrm = g->nrm;
	g->L.q = g->q;
	g->L.con = (unsigned char*)calloc(nv, 1);
	g->L.cq = (double*)calloc(3 * (size_t)nv, sizeof(double));
	g->T = (int*)calloc(2 * (size_t)g->m.nc, sizeof(int));
	g->T0 = (int*)calloc(2 * (size_t)g->m.nc, sizeof(int));
	g->K = (unsigned char*)calloc(g->m.nc, 1);
	g->Tc = (double*)calloc(4 * (size_t)g->m.nc, sizeof(double));
	/* the dense grid its own lattice: steps the coordinate differences */
	for (c = 0; c < g->m.nc; c++)
	{
		unsigned int a = g->m.cv[c], b = hm_to(&g->m, c);
		g->T[2 * c] = (int)(g->x[3 * b] - g->x[3 * a]);
		g->T[2 * c + 1] = (int)(g->x[3 * b + 1] - g->x[3 * a + 1]);
	}
}

/* each end's place in its cell: a on its lattice point, b a little off it (every change near a whole step dear) */
static void tg_rounding(tg_grid* g)
{
	unsigned int c;
	for (c = 0; c < 2 * g->m.nc; c++)
		g->Tc[4 * (c >> 1) + 2 + (c & 1)] = 0.05 * sin(1.3 * c + 0.4);
	memcpy(g->T0, g->T, sizeof(int) * 2 * g->m.nc);
}

static void tg_free(tg_grid* g)
{
	hm_free(&g->m);
	free(g->T); free(g->T0); free(g->K); free(g->Tc); free(g->L.con); free(g->L.cq);
}

/* edge u - w crossed from face f0 into the next: +d onto f0's corner, -d onto the other's (K = 0 there) */
static void tg_cross(tg_grid* g, unsigned int u, unsigned int w, unsigned int f0, int dx, int dy)
{
	unsigned int c1 = hm_find(&g->m, u, w), c2 = hm_find(&g->m, w, u);
	unsigned int c = c1 != HM_NONE && g->m.cf[c1] == f0 ? c1 : c2, o = c == c1 ? c2 : c1;
	g->T[2 * c] += dx;
	g->T[2 * c + 1] += dy;
	if (o != HM_NONE)
	{
		g->T[2 * o] -= dx;
		g->T[2 * o + 1] -= dy;
	}
}

/* the units the steps changed by, over the edges */
static unsigned int tg_change(tg_grid const* g)
{
	unsigned int c, n = 0;
	for (c = 0; c < g->m.nc; c++)
		if (g->m.co[c] == HM_NONE || c < g->m.co[c])
			n += (unsigned int)(abs(g->T[2 * c] - g->T0[2 * c]) + abs(g->T[2 * c + 1] - g->T0[2 * c + 1]));
	return n;
}

static void test_synthetic(void)
{
	char msg[300];
	tg_grid* g = (tg_grid*)malloc(sizeof(tg_grid));
	unsigned int left, odd, i0 = 8, j0 = 10, i;
	qf_report rep;
	/* a pair five faces apart along row j0: B(8) A(8) B(9) A(9) B(10) A(10) */
	tg_build(g);
	for (i = 0; i < 3; i++)
	{
		tg_cross(g, tg_v(i0 + i, j0), tg_v(i0 + i + 1, j0 + 1), tg_face(i0 + i, j0, 1), 1, 0);
		if (i < 2)
			tg_cross(g, tg_v(i0 + i + 1, j0), tg_v(i0 + i + 1, j0 + 1), tg_face(i0 + i, j0, 0), 1, 0);
	}
	tg_rounding(g);
	left = qf_consistent(&g->m, &g->L, g->T, g->K, g->Tc, &odd, &rep);
	sprintf(msg, "grid: a pair of dislocations five faces apart (supply %u): %u left, changed by %u (5), cost %.3f", rep.supply, left,
		tg_change(g), rep.cost * 1e-6);
	printf("  %s\n", msg);
	TEST_ASSERT(rep.supply == 2 && left == 0 && tg_change(g) == 5, msg);
	tg_free(g);
	/* one dislocation three faces from the border: A(8, 1) -> B(8, 0) -> A(8, 0) -> out */
	tg_build(g);
	tg_cross(g, tg_v(i0, 1), tg_v(i0 + 1, 1), tg_face(i0, 1, 0), 1, 0);
	tg_cross(g, tg_v(i0, 0), tg_v(i0 + 1, 1), tg_face(i0, 0, 1), 1, 0);
	tg_cross(g, tg_v(i0, 0), tg_v(i0 + 1, 0), tg_face(i0, 0, 0), 1, 0);
	tg_rounding(g);
	left = qf_consistent(&g->m, &g->L, g->T, g->K, g->Tc, &odd, &rep);
	sprintf(msg, "grid: a dislocation three faces from the border (supply %u): %u left, changed by %u (3)", rep.supply, left, tg_change(g));
	printf("  %s\n", msg);
	TEST_ASSERT(rep.supply == 1 && left == 0 && tg_change(g) == 3, msg);
	tg_free(g);
	/* Two cones: the frames turn by a quarter across a cut of three edges along row j0, from B(8) to A(9), the lattice
	   beyond it seen through g(P) = R^-1 (P - p) + p + d (p the cut's middle, d = (1, 0) odd: the cones between lattice
	   points). Corners crossing from the row's top (L) to its bottom (R): K = 1, T = g(P_b) - P_a; back: K = 3,
	   T = g^-1(P_b) - P_a. */
	tg_build(g);
	{
		int p[2] = { (int)i0 + 1, (int)j0 }, d[2] = { 1, 0 };
		unsigned int e, c, s1 = tg_face(i0, j0, 1), s2 = tg_face(i0 + 1, j0, 0), f, others = 0, before = 0, after = 0;
		unsigned char* fr = (unsigned char*)malloc(g->m.nc);
		for (e = 0; e < 3; e++)
		{
			/* even: the diagonal of square i0 + e / 2; odd: the side between it and the next */
			unsigned int k = i0 + (e + 1) / 2, bot = tg_v(e & 1 ? k : i0 + e / 2, j0), top = tg_v(e & 1 ? k : i0 + e / 2 + 1, j0 + 1), side;
			for (side = 0; side < 2; side++)
			{
				unsigned int a = side ? bot : top, b = side ? top : bot;
				int Pa[2], Pb[2], w[2], r[2];
				c = hm_find(&g->m, a, b);
				Pa[0] = (int)g->x[3 * a]; Pa[1] = (int)g->x[3 * a + 1];
				Pb[0] = (int)g->x[3 * b]; Pb[1] = (int)g->x[3 * b + 1];
				if (!side)
				{
					/* top (L) -> bottom (R): g(P_b) = R^-1 (P_b - p) + p + d */
					w[0] = Pb[0] - p[0]; w[1] = Pb[1] - p[1];
					qf_turn(w, 3, r);
					g->T[2 * c] = r[0] + p[0] + d[0] - Pa[0];
					g->T[2 * c + 1] = r[1] + p[1] + d[1] - Pa[1];
					g->K[c] = 1;
				}
				else
				{
					/* bottom (R) -> top (L): g^-1(P_b) = R (P_b - p - d) + p */
					w[0] = Pb[0] - p[0] - d[0]; w[1] = Pb[1] - p[1] - d[1];
					qf_turn(w, 1, r);
					g->T[2 * c] = r[0] + p[0] - Pa[0];
					g->T[2 * c + 1] = r[1] + p[1] - Pa[1];
					g->K[c] = 3;
				}
			}
		}
		for (f = 0; f < g->m.nf; f++)
		{
			int S[2], t = qf_closure(&g->m, f, g->T, g->K, fr, S);
			if (f == s1 || f == s2)
				before += (t & 1) && ((S[0] + S[1]) & 1);
			else
				others += t != 0 || S[0] || S[1];
		}
		tg_rounding(g);
		left = qf_consistent(&g->m, &g->L, g->T, g->K, g->Tc, &odd, &rep);
		for (f = 0; f < g->m.nf; f++)
		{
			int S[2], t = qf_closure(&g->m, f, g->T, g->K, fr, S);
			if (f == s1 || f == s2)
				after += (t & 1) && ((S[0] + S[1]) & 1);
		}
		sprintf(msg, "grid: two cones off lattice points (%u odd, %u other faces open): %u dislocations left, %u odd after (%u reported), "
			"%u path(s), changed by %u, cut %u", before, others, left, after, odd, rep.paths, tg_change(g), rep.cut);
		printf("  %s\n", msg);
		TEST_ASSERT(before == 2 && others == 0 && left == 0 && after == 0 && odd == 0 && rep.paths == 1, msg);
		free(fr);
	}
	tg_free(g);
	free(g);
}

/* ---- shapes -------------------------------------------------------------- */

static qaws_mesh_3d* tq_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int i;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_sizes = fs;
	d.face_count = nf;
	qaws_mesh_3d_create(&d, &m);
	free(p);
	return m;
}

/* the showcase's blob: an icosahedron split once 1-to-4 by public edits, pushed in and out smoothly */
static qaws_mesh_3d* tq_blob(void)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[36];
	double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	qaws_mesh_3d* m;
	unsigned int i, n, ne, nv0, w;
	unsigned int *e, *src;
	qaws_scalar* p;
	for (i = 0; i < 36; i++)
		x[i] = v[i] / sqrt(1 + g * g);
	m = tq_make(x, 12, f, NULL, 20);
	ne = qaws_mesh_3d_get_edge_count(m);
	nv0 = qaws_mesh_3d_get_vertex_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	src = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(m, e);
	for (i = 0; i < ne; i++)
		if (qaws_mesh_3d_split_edge(m, e[2 * i], e[2 * i + 1], NULL, NULL, &w) == QAWS_STATUS_OK)
		{
			src[2 * (w - nv0)] = e[2 * i];
			src[2 * (w - nv0) + 1] = e[2 * i + 1];
		}
	free(e);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(m, e);
	for (i = 0; i < ne; i++)
	{
		unsigned int a = e[2 * i], b = e[2 * i + 1], nw, old;
		if ((a >= nv0) == (b >= nv0))
			continue;
		nw = a >= nv0 ? a : b;
		old = a >= nv0 ? b : a;
		if (old != src[2 * (nw - nv0)] && old != src[2 * (nw - nv0) + 1] && qaws_mesh_3d_can_flip_edge(m, a, b))
			qaws_mesh_3d_flip_edge(m, a, b);
	}
	free(e);
	free(src);
	n = qaws_mesh_3d_get_vertex_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	qaws_mesh_3d_get_positions(m, p);
	for (i = 0; i < n; i++)
	{
		double q0 = p[3 * i], q1 = p[3 * i + 1], q2 = p[3 * i + 2], l = sqrt(q0 * q0 + q1 * q1 + q2 * q2), s;
		q0 /= l; q1 /= l; q2 /= l;
		s = 1 + 0.32 * sin(2.2 * q0 + 0.4) * cos(1.7 * q1) + 0.2 * q2 * q2;
		p[3 * i] = (qaws_scalar)(q0 * 1.25 * s);
		p[3 * i + 1] = (qaws_scalar)(q1 * s);
		p[3 * i + 2] = (qaws_scalar)(q2 * 0.9 * s);
	}
	qaws_mesh_3d_set_positions(m, p);
	free(p);
	return m;
}

static qaws_mesh_3d* tq_can(void)
{
	unsigned int n = 8, i, k = 0;
	double x[3 * 18];
	unsigned int fv[3 * 32];
	qaws_mesh_3d* m;
	for (i = 0; i < n; i++)
	{
		double a = 2 * M_PI * i / n;
		x[3 * i] = cos(a); x[3 * i + 1] = sin(a); x[3 * i + 2] = -0.8;
		x[3 * (n + i)] = cos(a); x[3 * (n + i) + 1] = sin(a); x[3 * (n + i) + 2] = 0.8;
	}
	x[6 * n] = 0; x[6 * n + 1] = 0; x[6 * n + 2] = -0.8;
	x[6 * n + 3] = 0; x[6 * n + 4] = 0; x[6 * n + 5] = 0.8;
	for (i = 0; i < n; i++)
	{
		unsigned int a = i, b = (i + 1) % n;
		fv[k++] = a; fv[k++] = b; fv[k++] = n + b;
		fv[k++] = a; fv[k++] = n + b; fv[k++] = n + a;
		fv[k++] = 2 * n; fv[k++] = b; fv[k++] = a;
		fv[k++] = 2 * n + 1; fv[k++] = n + a; fv[k++] = n + b;
	}
	m = tq_make(x, 2 * n + 2, fv, NULL, 4 * n);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(50 * M_PI / 180), 0, 0);
	return m;
}

static qaws_mesh_3d* tq_torus(unsigned int nu, unsigned int nv, double R, double r)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * M_PI * i / nu, v = 2 * M_PI * (j + 0.5 * (i & 1)) / nv;
			x[3 * (j * nu + i)] = (R + r * cos(v)) * cos(u);
			x[3 * (j * nu + i) + 1] = (R + r * cos(v)) * sin(u);
			x[3 * (j * nu + i) + 2] = r * sin(v);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, d = ((j + 1) % nv) * nu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	m = tq_make(x, nu * nv, f, NULL, 2 * nu * nv);
	free(x);
	free(f);
	return m;
}

static qaws_mesh_3d* tq_quads(qaws_mesh_3d const* m, double L, int pure, qaws_mesh_quality* q, double* seconds)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* out = NULL;
	clock_t t0 = clock();
	memset(&d, 0, sizeof d);
	d.target_length = (qaws_scalar)L;
	d.pure = pure;
	if (qaws_mesh_quadrangulate(m, &d, &out, NULL) != QAWS_STATUS_OK)
		return NULL;
	if (seconds)
		*seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
	qaws_mesh_3d_measure_quality(out, q);
	return out;
}

static void test_blob(void)
{
	qaws_mesh_3d *blob = tq_blob(), *q1, *q2, *q3;
	qaws_mesh_quality a, b, c, in;
	char msg[400];
	double t1 = 0, t2 = 0, area;
	unsigned int nf;
	qaws_mesh_3d_measure_quality(blob, &in);
	area = (double)in.area;
	q1 = tq_quads(blob, 0.16, 0, &a, &t1);
	sprintf(msg, "blob quad-dominant at 0.16: %u quads, %u triangles, %u polygons (0; floor 2; was 5 and 1), irregular %u (10), %.2f s",
		q1 ? a.quads : 0, q1 ? a.triangles : 0, q1 ? a.polygons : 0, q1 ? a.irregular : 0, t1);
	printf("  %s\n", msg);
	TEST_ASSERT(q1 && a.triangles + a.polygons == 0 && a.irregular <= 10, msg);
	q2 = tq_quads(blob, 0.11, 1, &b, &t2);
	nf = q2 ? qaws_mesh_3d_get_face_count(q2) : 0;
	sprintf(msg, "blob pure at 0.11: %u quads of %u faces (area / 0.11^2 = %.0f, +-25%%), irregular %u (8; floor 10; was 12), mean edge %.3f, %.2f s",
		q2 ? b.quads : 0, nf, area / (0.11 * 0.11), q2 ? b.irregular : 0, q2 ? (double)b.mean_edge : 0, t2);
	printf("  %s\n", msg);
	/* (the floor: measured 10, the 8 cones and a pair where a cone could not be moved onto a vertex of its face without
	   steps of two; the fallback at twice the target gives 8 but not at this spacing) */
	TEST_ASSERT(q2 && b.quads == nf && fabs(nf - area / (0.11 * 0.11)) < 0.25 * area / (0.11 * 0.11) && b.irregular <= 10, msg);
	/* the same twice */
	q3 = tq_quads(blob, 0.11, 1, &c, NULL);
	{
		int same = q2 && q3 && qaws_mesh_3d_get_vertex_count(q2) == qaws_mesh_3d_get_vertex_count(q3) && qaws_mesh_3d_get_face_count(q2) == nf;
		if (same)
		{
			unsigned int nv = qaws_mesh_3d_get_vertex_count(q2), nc = qaws_mesh_3d_get_corner_count(q2);
			qaws_scalar *p2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *p3 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			unsigned int *f2 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + nf)), *f3 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + nf));
			qaws_mesh_3d_get_positions(q2, p2);
			qaws_mesh_3d_get_positions(q3, p3);
			qaws_mesh_3d_get_faces(q2, f2, f2 + nc);
			qaws_mesh_3d_get_faces(q3, f3, f3 + nc);
			same = !memcmp(p2, p3, sizeof(qaws_scalar) * 3 * nv) && !memcmp(f2, f3, sizeof(unsigned int) * (nc + nf));
			free(p2); free(p3); free(f2); free(f3);
		}
		sprintf(msg, "blob pure twice: identical arrays %d", same);
		TEST_ASSERT(same, msg);
	}
	qaws_mesh_3d_destroy(q1);
	qaws_mesh_3d_destroy(q2);
	qaws_mesh_3d_destroy(q3);
	qaws_mesh_3d_destroy(blob);
}

static void test_can(void)
{
	qaws_mesh_3d *can = tq_can(), *q;
	qaws_mesh_quality a;
	qaws_mesh_topology top;
	char msg[400];
	unsigned int nf, nfe, f, c0, bad = 0, onrim = 1, i;
	double t = 0;
	q = tq_quads(can, 0.14, 1, &a, &t);
	if (!q)
	{
		TEST_ASSERT(0, "can: quadrangulated");
		qaws_mesh_3d_destroy(can);
		return;
	}
	qaws_mesh_3d_get_topology(q, &top);
	nf = qaws_mesh_3d_get_face_count(q);
	nfe = qaws_mesh_3d_get_feature_edge_count(q);
	{
		unsigned int nvq = qaws_mesh_3d_get_vertex_count(q), *ff = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(q));
		unsigned int *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nvq);
		qaws_mesh_3d_get_positions(q, p);
		qaws_mesh_3d_get_faces(q, ff, fs);
		for (f = 0, c0 = 0; f < nf; c0 += fs[f], f++)
		{
			double cen[3] = { 0, 0, 0 }, nn[3] = { 0, 0, 0 }, out[3];
			unsigned int j;
			for (j = 0; j < fs[f]; j++)
			{
				qaws_scalar const *u = p + 3 * ff[c0 + j], *w = p + 3 * ff[c0 + (j + 1) % fs[f]];
				nn[0] += (u[1] - w[1]) * (u[2] + w[2]); nn[1] += (u[2] - w[2]) * (u[0] + w[0]); nn[2] += (u[0] - w[0]) * (u[1] + w[1]);
				cen[0] += u[0] / fs[f]; cen[1] += u[1] / fs[f]; cen[2] += u[2] / fs[f];
			}
			if (fabs(cen[2]) > 0.79) { out[0] = 0; out[1] = 0; out[2] = cen[2] > 0 ? 1 : -1; }
			else { out[0] = cen[0]; out[1] = cen[1]; out[2] = 0; }
			bad += nn[0] * out[0] + nn[1] * out[1] + nn[2] * out[2] <= 0;
		}
		qaws_mesh_3d_get_feature_edges(q, fe);
		for (i = 0; i < 2 * nfe; i++)
			onrim &= fabs(fabs(p[3 * fe[i] + 2]) - 0.8) < 1e-6;
		free(ff); free(fs); free(fe); free(p);
	}
	sprintf(msg, "can pure at 0.14: %u quads of %u, closed %d, %u turned in, %u crease edges on the rims %u, irregular %u (12; was 22), %.2f s",
		a.quads, nf, top.closed, bad, nfe, onrim, a.irregular, t);
	printf("  %s\n", msg);
	TEST_ASSERT(a.quads == nf && qaws_mesh_3d_validate(q) == QAWS_STATUS_OK && top.closed && bad == 0 && nfe > 40 && onrim && a.irregular <= 12, msg);
	qaws_mesh_3d_destroy(q);
	qaws_mesh_3d_destroy(can);
}

static void test_torus(void)
{
	qaws_mesh_3d *t = tq_torus(30, 12, 2, 0.75), *q;
	qaws_mesh_quality a;
	char msg[300];
	double area = 4 * M_PI * M_PI * 2 * 0.75;
	unsigned int nf;
	q = tq_quads(t, 0.3, 1, &a, NULL);
	nf = q ? qaws_mesh_3d_get_face_count(q) : 0;
	sprintf(msg, "torus pure at 0.3: %u quads of %u (area / 0.09 = %.0f, +-15%%), irregular %u (0)", q ? a.quads : 0, nf, area / 0.09, q ? a.irregular : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(q && a.quads == nf && a.irregular == 0 && fabs(a.quads - area / 0.09) < 0.15 * area / 0.09, msg);
	qaws_mesh_3d_destroy(q);
	qaws_mesh_3d_destroy(t);
}

static void test_patch(void)
{
	qaws_vec3 cp[49];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	qaws_generate_desc gd;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_quality a;
	unsigned int i, j;
	char msg[300];
	for (i = 0; i < 7; i++)
		for (j = 0; j < 7; j++)
		{
			cp[i * 7 + j].x = (qaws_scalar)(i / 3.0 - 1);
			cp[i * 7 + j].y = (qaws_scalar)(j / 3.0 - 1);
			cp[i * 7 + j].z = (qaws_scalar)(0.55 * sin(1.7 * i) * cos(1.3 * j));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = cp;
	d.u_point_count = 7;
	d.v_point_count = 7;
	qaws_surface_create_bspline(&d, &s);
	memset(&gd, 0, sizeof gd);
	gd.target_length = (qaws_scalar)0.15;
	gd.elements = QAWS_ELEMENTS_QUADS;
	qaws_mesh_generate_surface(s, &gd, &m, NULL);
	if (m)
		qaws_mesh_3d_measure_quality(m, &a);
	sprintf(msg, "wavy patch in quads at 0.15: %u quads of %u, irregular %u (2; was 4)", m ? a.quads : 0, m ? qaws_mesh_3d_get_face_count(m) : 0,
		m ? a.irregular : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(m && a.quads == qaws_mesh_3d_get_face_count(m) && a.irregular <= 2, msg);
	qaws_mesh_3d_destroy(m);
	qaws_surface_destroy(s);
}

/* consistency alone on a dense remesh of the blob (about 50k faces) */
static void test_timing(void)
{
	qaws_mesh_3d *blob = tq_blob(), *dense = NULL;
	qaws_remesh_desc rd;
	mf_hier h;
	int* T;
	unsigned char* K;
	double* Tc;
	unsigned int left = HM_NONE, odd = 0;
	qf_report rep;
	char msg[400];
	clock_t t0;
	double sec = 0, rf = 0.09;
	memset(&rd, 0, sizeof rd);
	memset(&h, 0, sizeof h);
	memset(&rep, 0, sizeof rep);
	rd.target_length = (qaws_scalar)(rf / 3);
	qaws_mesh_remesh(blob, &rd, &dense, NULL);
	if (!dense || mf_build(&h, &dense->hm, NULL) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "timing: dense blob");
		qaws_mesh_3d_destroy(dense);
		qaws_mesh_3d_destroy(blob);
		return;
	}
	mf_orientations(&h, 8, NULL, 0);
	mf_positions(&h, rf, 8);
	T = (int*)malloc(sizeof(int) * 2 * dense->hm.nc);
	K = (unsigned char*)malloc(dense->hm.nc);
	Tc = (double*)malloc(sizeof(double) * 4 * dense->hm.nc);
	qf_steps(&dense->hm, &h.lv[0], rf, T, K, Tc);
	t0 = clock();
	left = qf_consistent(&dense->hm, &h.lv[0], T, K, Tc, &odd, &rep);
	sec = (double)(clock() - t0) / CLOCKS_PER_SEC;
	sprintf(msg, "consistency on %u dense faces: supply %u, %u components cut, %u unrouted, changed by %u, cost %.1f, %u parity paths, "
		"%u dislocations and %u odd cones left, %u steps of two; %.3f s (0.2; under 2)", dense->hm.nf, rep.supply, rep.cut, rep.unrouted,
		rep.changed, rep.cost * 1e-6, rep.paths, left, odd, rep.twos, sec);
	printf("  %s\n", msg);
	TEST_ASSERT(dense->hm.nf > 40000 && left == 0 && odd == 0 && sec < 2, msg);
	free(T); free(K); free(Tc);
	mf_free(&h);
	qaws_mesh_3d_destroy(dense);
	qaws_mesh_3d_destroy(blob);
}

int test_134_quad_flow_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 134: quads by min-cost flow\n");
	test_flow_module();
	test_synthetic();
	test_timing();
	test_blob();
	test_can();
	test_torus();
	test_patch();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
