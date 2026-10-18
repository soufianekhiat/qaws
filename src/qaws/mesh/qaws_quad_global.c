/*
 * GLOBAL quads: the integer-grid map of INSTANT's dense mesh from its
 * consistent steps, solved once, and its quads read off by tracing the
 * map's unit isolines (qaws_quad_global.h).
 *
 * Per dense vertex v, d_v: its place relative to its lattice point, in
 * steps of its frame (q, n x q). A corner c (a -> b) maps b's frame into
 * a's: x_a = T_c + turn(x_b, -K_c). Each triangle of the map takes the
 * frame of one of its vertices as its chart, its corners' uv turn(d, r) +
 * t with whole r and t; across an edge the two charts differ by a quarter
 * turn and a whole shift (a transition), read off a vertex both triangles
 * see through the same steps.
 */

#include "qaws_quad_global.h"
#include "qaws_quad_flow.h"
#include "../internal/qaws_internal_ldlt.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define QG_TOL 1e-9         /* a barycentric within this of 0: on that side */
#define QG_TILT 1e-6        /* the isolines' nudge to their left, in steps */

/* a triangle of the map: a dense face, or a third of a singular one round its cone */
typedef struct qg_tri
{
	double U[6], X[9];
	unsigned int v[3];      /* dense vertices; HM_NONE: the cone of a split face */
	int r[3], t[6];         /* per corner: uv = turn(d_v, r) + t (the cone: t its uv) */
	unsigned int nb[3];     /* across edge k (corner k to k + 1): the neighbour * 3 + its edge; HM_NONE on the border */
	int tr[9];              /* edge k's transition: uv there = turn(uv, tr[3k]) + (tr[3k + 1], tr[3k + 2]) */
	unsigned char line;     /* bit k: edge k on a feature or the border */
} qg_tri;

/* a lattice point seen from a triangle (its coordinates there) and the half-edges leaving it along each axis there */
typedef struct qg_reg
{
	unsigned int t, p, he[4];
	int i, j;
} qg_reg;

typedef struct qg_he
{
	unsigned int reg, to, end, next;
	int a, Q[2], aend;
} qg_he;

typedef struct qg_ctx
{
	qg_tri* tri;
	unsigned int nt;
	qg_reg* reg;
	unsigned int nreg, capreg;
	unsigned int* hash;
	unsigned int mask;
	unsigned int* uf;           /* per point: union-find parent */
	unsigned int* home;         /* per point: the triangle it was found in */
	double* bary;               /* and its barycentrics there */
	unsigned int np, capp;
} qg_ctx;

void qg_map_free(qg_map* m)
{
	free(m->x);
	free(m->uv);
	memset(m, 0, sizeof *m);
}

static void qg_turn(double const* p, int r, double* o)
{
	double x = p[0], y = p[1];
	switch (r & 3)
	{
	case 0: o[0] = x; o[1] = y; break;
	case 1: o[0] = -y; o[1] = x; break;
	case 2: o[0] = -x; o[1] = -y; break;
	default: o[0] = y; o[1] = -x; break;
	}
}

static double qg_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void qg_cross3(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double qg_cross(double const* o, double const* a, double const* b)
{
	return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
}

static void qg_bary(qg_tri const* t, double const* p, double* b)
{
	double A = qg_cross(t->U, t->U + 2, t->U + 4);
	b[0] = qg_cross(p, t->U + 2, t->U + 4) / A;
	b[1] = qg_cross(p, t->U + 4, t->U) / A;
	b[2] = qg_cross(p, t->U, t->U + 2) / A;
}

/* p (uv) moved by edge k's transition */
static void qg_across(qg_tri const* t, unsigned int k, double* p)
{
	qg_turn(p, t->tr[3 * k], p);
	p[0] += t->tr[3 * k + 1];
	p[1] += t->tr[3 * k + 2];
}

static void qg_iacross(qg_tri const* t, unsigned int k, int* Q)
{
	qf_turn(Q, t->tr[3 * k], Q);
	Q[0] += t->tr[3 * k + 1];
	Q[1] += t->tr[3 * k + 2];
}

/* whether t's closure holds p and its sector at p the axis a nudged to its left (turn((1, tilt), a)): sides p is on must
   grow along it */
static int qg_owns(qg_tri const* t, double const* p, int a)
{
	static double const e0[2] = { 1, QG_TILT };
	double b[3], bb[3], e[2], pe[2];
	int k;
	qg_bary(t, p, b);
	qg_turn(e0, a, e);
	pe[0] = p[0] + e[0];
	pe[1] = p[1] + e[1];
	qg_bary(t, pe, bb);
	for (k = 0; k < 3; k++)
		if (b[k] < -QG_TOL || (b[k] <= QG_TOL && bb[k] <= b[k]))
			return 0;
	return 1;
}

static unsigned int qg_root(unsigned int* uf, unsigned int i)
{
	while (uf[i] != i)
	{
		uf[i] = uf[uf[i]];
		i = uf[i];
	}
	return i;
}

static unsigned int qg_slot(unsigned int t, int i, int j, unsigned int mask)
{
	return ((t * 2654435761u) ^ ((unsigned int)i * 40503u) ^ ((unsigned int)j * 2246822519u)) & mask;
}

static unsigned int qg_find(qg_ctx const* c, unsigned int t, int i, int j)
{
	unsigned int s;
	for (s = qg_slot(t, i, j, c->mask); c->hash[s] != HM_NONE; s = (s + 1) & c->mask)
	{
		qg_reg const* g = &c->reg[c->hash[s]];
		if (g->t == t && g->i == i && g->j == j)
			return c->hash[s];
	}
	return HM_NONE;
}

/* point p seen from t at (i, j): a second point there merged with it; 0 on allocation failure */
static int qg_register(qg_ctx* c, unsigned int t, int i, int j, unsigned int p)
{
	unsigned int r = qg_find(c, t, i, j), s;
	if (r != HM_NONE)
	{
		unsigned int a = qg_root(c->uf, c->reg[r].p), b = qg_root(c->uf, p);
		if (a != b)
			c->uf[b] = a;
		return 1;
	}
	if (2 * (c->nreg + 1) > c->mask + 1)
	{
		unsigned int n = 2 * (c->mask + 1), k;
		unsigned int* h = (unsigned int*)malloc(sizeof(unsigned int) * n);
		qg_reg* g = (qg_reg*)realloc(c->reg, sizeof(qg_reg) * (n / 2 + 1));
		if (g)
			c->reg = g;
		if (!h || !g)
		{
			free(h);
			return 0;
		}
		for (k = 0; k < n; k++)
			h[k] = HM_NONE;
		free(c->hash);
		c->hash = h;
		c->mask = n - 1;
		for (k = 0; k < c->nreg; k++)
		{
			for (s = qg_slot(c->reg[k].t, c->reg[k].i, c->reg[k].j, c->mask); h[s] != HM_NONE; s = (s + 1) & c->mask) {}
			h[s] = k;
		}
	}
	for (s = qg_slot(t, i, j, c->mask); c->hash[s] != HM_NONE; s = (s + 1) & c->mask) {}
	c->hash[s] = c->nreg;
	c->reg[c->nreg].t = t;
	c->reg[c->nreg].i = i;
	c->reg[c->nreg].j = j;
	c->reg[c->nreg].p = p;
	c->reg[c->nreg].he[0] = c->reg[c->nreg].he[1] = c->reg[c->nreg].he[2] = c->reg[c->nreg].he[3] = HM_NONE;
	c->nreg++;
	return 1;
}

/* The lattice points of the map: each triangle's (its closure's), each registered in every triangle round it, those
   found twice merged. */
static qaws_status qg_points(qg_ctx* c)
{
	unsigned int t, k;
	for (t = 0; t < c->nt; t++)
	{
		qg_tri const* T = &c->tri[t];
		double lo[2], hi[2];
		int i, j;
		for (k = 0; k < 2; k++)
		{
			lo[k] = fmin(T->U[k], fmin(T->U[2 + k], T->U[4 + k]));
			hi[k] = fmax(T->U[k], fmax(T->U[2 + k], T->U[4 + k]));
		}
		for (i = (int)ceil(lo[0] - 1e-7); i <= (int)floor(hi[0] + 1e-7); i++)
			for (j = (int)ceil(lo[1] - 1e-7); j <= (int)floor(hi[1] + 1e-7); j++)
			{
				double P[2], b[3];
				unsigned int p, z = 0, m = 0, w;
				P[0] = i;
				P[1] = j;
				qg_bary(T, P, b);
				if (b[0] < -QG_TOL || b[1] < -QG_TOL || b[2] < -QG_TOL || qg_find(c, t, i, j) != HM_NONE)
					continue;
				if (c->np == c->capp)
				{
					unsigned int n = 2 * c->capp + 64;
					unsigned int* uf = (unsigned int*)realloc(c->uf, sizeof(unsigned int) * n);
					unsigned int* home = uf ? (unsigned int*)realloc(c->home, sizeof(unsigned int) * n) : NULL;
					double* bary = home ? (double*)realloc(c->bary, sizeof(double) * 3 * n) : NULL;
					if (uf)
						c->uf = uf;
					if (home)
						c->home = home;
					if (!bary)
						return QAWS_STATUS_ALLOCATION_FAILURE;
					c->bary = bary;
					c->capp = n;
				}
				p = c->np++;
				c->uf[p] = p;
				c->home[p] = t;
				memcpy(c->bary + 3 * (size_t)p, b, sizeof b);
				if (!qg_register(c, t, i, j, p))
					return QAWS_STATUS_ALLOCATION_FAILURE;
				for (k = 0; k < 3; k++)
					if (b[k] <= QG_TOL)
						z++;
					else
						m = k;
				if (z == 1)
				{
					/* on an edge: the triangle across */
					for (k = 0; b[k] > QG_TOL; k++) {}
					w = T->nb[(k + 1) % 3];
					if (w != HM_NONE)
					{
						int Q[2];
						Q[0] = i;
						Q[1] = j;
						qg_iacross(T, (k + 1) % 3, Q);
						if (!qg_register(c, w / 3, Q[0], Q[1], p))
							return QAWS_STATUS_ALLOCATION_FAILURE;
					}
				}
				else if (z == 2)
				{
					/* at a corner: every triangle round it, clockwise (across the corner's own edge) then
					   counter-clockwise (across the one before) until the border or back */
					unsigned int dir, round = 0;
					for (dir = 0; dir < 2 && !round; dir++)
					{
						unsigned int cur = t, cm = m, n;
						int Q[2];
						Q[0] = i;
						Q[1] = j;
						for (n = 0; n < 64; n++)
						{
							unsigned int e = dir ? (cm + 2) % 3 : cm, nx = c->tri[cur].nb[e];
							if (nx == HM_NONE)
								break;
							qg_iacross(&c->tri[cur], e, Q);
							cur = nx / 3;
							cm = dir ? nx % 3 : (nx % 3 + 1) % 3;
							if (cur == t)
							{
								round = 1;
								break;
							}
							if (!qg_register(c, cur, Q[0], Q[1], p))
								return QAWS_STATUS_ALLOCATION_FAILURE;
						}
					}
				}
			}
	}
	return QAWS_STATUS_OK;
}

/* the registration in or round t of the lattice point at Q (t's coordinates) whose sector holds axis *a (updated into
   that triangle); HM_NONE when none does (it points off the border) */
static unsigned int qg_owner(qg_ctx const* c, unsigned int t, int const* Q0, int* a)
{
	double P[2], b[3];
	unsigned int k, z = 0, m = 0, dir;
	int Q[2];
	P[0] = Q0[0];
	P[1] = Q0[1];
	if (qg_owns(&c->tri[t], P, *a))
		return qg_find(c, t, Q0[0], Q0[1]);
	qg_bary(&c->tri[t], P, b);
	for (k = 0; k < 3; k++)
		if (fabs(b[k]) <= QG_TOL)
			z++;
		else
			m = k;
	for (dir = 0; dir < 2 && z; dir++)
	{
		unsigned int cur = t, cm = m, n;
		int ax = *a;
		Q[0] = Q0[0];
		Q[1] = Q0[1];
		for (n = 0; n < (z == 1 ? 1u : 4u); n++)
		{
			unsigned int e, nx;
			if (z == 1)
			{
				for (k = 0; fabs(b[k]) > QG_TOL; k++) {}
				e = (k + 1) % 3;
			}
			else
				e = dir ? cm : (cm + 2) % 3;
			nx = c->tri[cur].nb[e];
			if (nx == HM_NONE)
				break;
			qg_iacross(&c->tri[cur], e, Q);
			ax = (ax + c->tri[cur].tr[3 * e]) & 3;
			cur = nx / 3;
			cm = dir ? (nx % 3 + 1) % 3 : nx % 3;
			if (cur == t)
				break;
			P[0] = Q[0];
			P[1] = Q[1];
			if (qg_owns(&c->tri[cur], P, ax))
			{
				*a = ax;
				return qg_find(c, cur, Q[0], Q[1]);
			}
		}
		if (z == 1)
			break;
	}
	return HM_NONE;
}

/* the isoline from registration r's point one step along axis a, nudged to its left: the registration of the point it
   reaches (HM_NONE: lost), where it arrives (the triangle, the point's and the axis's coordinates there) */
static unsigned int qg_trace(qg_ctx const* c, unsigned int r, int a, unsigned int* out_t, int* out_Q, int* out_a)
{
	static double const e0[2] = { 1, QG_TILT };
	static int const i0[2] = { 1, 0 };
	qg_reg const* g = &c->reg[r];
	unsigned int t = g->t, step, entry = 3;
	double s[2], q[2], e[2];
	int Q[2], d[2], ax = a;
	s[0] = g->i;
	s[1] = g->j;
	qf_turn(i0, a, d);
	Q[0] = g->i + d[0];
	Q[1] = g->j + d[1];
	qg_turn(e0, a, e);
	q[0] = s[0] + e[0];
	q[1] = s[1] + e[1];
	for (step = 0;; step++)
	{
		qg_tri const* T = &c->tri[t];
		double B[3], S[3], lam = HUGE_VAL;
		unsigned int k, best = 3, ex, nx;
		if (step > 4096)
			return HM_NONE;
		qg_bary(T, q, B);
		if (B[0] >= -1e-12 && B[1] >= -1e-12 && B[2] >= -1e-12)
		{
			/* the nudged end reached; the point itself maybe a hair across an edge: on to it */
			if (qg_find(c, t, Q[0], Q[1]) != HM_NONE || (q[0] == Q[0] && q[1] == Q[1]))
				break;
			q[0] = Q[0];
			q[1] = Q[1];
			entry = 3;
			continue;
		}
		qg_bary(T, s, S);
		for (k = 0; k < 3; k++)
			if (B[k] < 0 && (k + 1) % 3 != entry)
			{
				double l = S[k] > 0 ? S[k] / (S[k] - B[k]) : 0;
				if (l < lam)
				{
					lam = l;
					best = k;
				}
			}
		if (best == 3)
			return HM_NONE;
		ex = (best + 1) % 3;
		nx = T->nb[ex];
		if (nx == HM_NONE)
			return HM_NONE;
		s[0] += lam * (q[0] - s[0]);
		s[1] += lam * (q[1] - s[1]);
		qg_across(T, ex, s);
		qg_across(T, ex, q);
		qg_iacross(T, ex, Q);
		ax = (ax + T->tr[3 * ex]) & 3;
		entry = nx % 3;
		t = nx / 3;
	}
	*out_t = t;
	out_Q[0] = Q[0];
	out_Q[1] = Q[1];
	*out_a = ax;
	return qg_find(c, t, Q[0], Q[1]);
}

/* The quads: per registration, a half-edge along each axis its sector holds, traced to the next point; the next
   half-edge of a face a quarter turn left at its end; the faces the cycles, all quads or none. */
static qaws_status qg_read(qg_ctx* c, qaws_hm const* dm, mf_level const* L, qd_out* o)
{
	qg_he* he = NULL;
	unsigned int nhe = 0, r, a, h, *seen = NULL, *out = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (r = 0; r < c->nreg; r++)
		for (a = 0; a < 4; a++)
		{
			double P[2];
			P[0] = c->reg[r].i;
			P[1] = c->reg[r].j;
			if (qg_owns(&c->tri[c->reg[r].t], P, (int)a))
				c->reg[r].he[a] = nhe++;
		}
	he = (qg_he*)malloc(sizeof(qg_he) * (nhe + 1));
	seen = (unsigned int*)malloc(sizeof(unsigned int) * (nhe + 1));
	out = (unsigned int*)malloc(sizeof(unsigned int) * (c->np + 1));
	o->fv = (unsigned int*)malloc(sizeof(unsigned int) * (nhe + 1));
	o->fs = (unsigned int*)malloc(sizeof(unsigned int) * (nhe / 4 + 1));
	if (!he || !seen || !out || !o->fv || !o->fs)
		goto done;
	st = QAWS_STATUS_INVALID_TOPOLOGY;
	for (r = 0; r < c->nreg; r++)
		for (a = 0; a < 4; a++)
			if ((h = c->reg[r].he[a]) != HM_NONE)
			{
				unsigned int to = qg_trace(c, r, (int)a, &he[h].end, he[h].Q, &he[h].aend);
				if (to == HM_NONE)
					goto done;
				he[h].reg = r;
				he[h].a = (int)a;
				he[h].to = qg_root(c->uf, c->reg[to].p);
				seen[h] = HM_NONE;
			}
	/* each half-edge's next, the next of exactly one */
	for (h = 0; h < nhe; h++)
	{
		int ax = (he[h].aend + 1) & 3;
		unsigned int g = qg_owner(c, he[h].end, he[h].Q, &ax);
		he[h].next = g == HM_NONE ? HM_NONE : c->reg[g].he[ax];
		if (he[h].next == HM_NONE || seen[he[h].next] != HM_NONE)
			goto done;
		seen[he[h].next] = h;
	}
	for (h = 0; h < c->np; h++)
		out[h] = HM_NONE;
	for (h = 0; h < nhe; h++)
		seen[h] = 0;
	for (h = 0; h < nhe; h++)
	{
		unsigned int g = h, n = 0, k, *q = o->fv + o->nc;
		if (seen[h])
			continue;
		do
		{
			if (n == 4 || seen[g])
				goto done;
			seen[g] = 1;
			q[n++] = qg_root(c->uf, c->reg[he[g].reg].p);
			g = he[g].next;
		} while (g != h);
		if (n != 4 || q[0] == q[2] || q[1] == q[3])
			goto done;
		for (k = 0; k < 4; k++)
			if (out[q[k]] == HM_NONE)
				out[q[k]] = o->nv++;
		o->fs[o->nf++] = 4;
		o->nc += 4;
	}
	if (!o->nf)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	o->x = (double*)malloc(sizeof(double) * 3 * (o->nv + 1));
	o->foot = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (o->nv + 1));
	o->kind = (unsigned char*)malloc(o->nv + 1);
	if (!o->x || !o->foot || !o->kind)
		goto done;
	/* the vertices: at their barycentrics in the triangle they were found in; a corner where a pinned vertex is, a point on
	   a line vertex or a line edge on the line */
	for (h = 0; h < c->np; h++)
	{
		unsigned int q = out[h], k, best = 0, z = 0, m = 0, t = c->home[h];
		double const* b = c->bary + 3 * (size_t)h;
		qg_tri const* T = &c->tri[t];
		if (q == HM_NONE)
			continue;
		for (k = 0; k < 3; k++)
		{
			o->x[3 * (size_t)q + k] = b[0] * T->X[k] + b[1] * T->X[3 + k] + b[2] * T->X[6 + k];
			if (b[k] <= QG_TOL)
				z++;
			else
				m = k;
			if (T->v[k] != HM_NONE && (T->v[best] == HM_NONE || b[k] > b[best]))
				best = k;
		}
		o->foot[q] = dm->foot[T->v[best]];
		o->kind[q] = 0;
		if (z == 2 && T->v[m] != HM_NONE && (L->con[T->v[m]] & MF_PINNED))
		{
			memcpy(o->x + 3 * (size_t)q, dm->x + 3 * (size_t)T->v[m], sizeof(double) * 3);
			o->foot[q] = dm->foot[T->v[m]];
			o->kind[q] = 2;
		}
		else if (z == 2 && T->v[m] != HM_NONE && (L->con[T->v[m]] & MF_ON_LINE))
			o->kind[q] = 1;
		else if (z == 1)
		{
			for (k = 0; b[k] > QG_TOL; k++) {}
			o->kind[q] = (unsigned char)((T->line >> ((k + 1) % 3)) & 1);
		}
	}
	for (h = 0; h < o->nc; h++)
		o->fv[h] = out[o->fv[h]];
	st = QAWS_STATUS_OK;
done:
	free(he);
	free(seen);
	free(out);
	return st;
}

/* vertex b's q turned into the plane of normal n (the rotation taking its normal onto n: unfolded across a crease, where
   a line vertex's frame lies in one of its faces), unit; 0 when it vanishes there */
static int qg_axis(mf_level const* L, unsigned int b, double const* n, double* X)
{
	double const *q = L->q + 3 * (size_t)b, *m = L->nrm + 3 * (size_t)b;
	double k[3], kq[3], kkq[3], c = qg_dot(m, n), l;
	unsigned int i;
	qg_cross3(m, n, k);
	qg_cross3(k, q, kq);
	qg_cross3(k, kq, kkq);
	for (i = 0; i < 3; i++)
		X[i] = q[i] + kq[i] + (c > -0.99 ? kkq[i] / (1 + c) : 0);
	c = qg_dot(X, n);
	for (i = 0; i < 3; i++)
		X[i] -= c * n[i];
	l = sqrt(qg_dot(X, X));
	if (!(l > 1e-9))
		return 0;
	for (i = 0; i < 3; i++)
		X[i] /= l;
	return 1;
}

/* The places d (2 per dense vertex) minimizing the map's energy, the steps fixed: per triangle A |J - F|^2 with J the
   map's gradient and F the field of its chart's vertex in the triangle's plane, in steps; held components (fix, bits 1
   and 2) kept at 0; the rest pulled toward d by a millionth of a typical weight (the torus's free shift). */
static qaws_status qg_solve(qaws_hm const* dm, mf_level const* L, double rho, qg_tri const* tri, unsigned int nt, unsigned char const* fix,
	double* d)
{
	unsigned int n = dm->nv, v, t, k, l, *nb0 = NULL, *nb = NULL, *rp = NULL, *col = NULL, ring[64], tot = 0;
	double *val = NULL, *rhs = NULL, *xyz = NULL, mean = 0;
	qaws_internal_ldlt f;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&f, 0, sizeof f);
	nb0 = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	if (!nb0)
		goto done;
	for (v = 0; v < n; v++)
	{
		nb0[v] = tot;
		tot += 1 + (dm->vc[v] == HM_NONE ? 0 : hm_ring(dm, v, ring, 64));
	}
	nb0[n] = tot;
	nb = (unsigned int*)malloc(sizeof(unsigned int) * (tot + 1));
	rp = (unsigned int*)malloc(sizeof(unsigned int) * (2 * (size_t)n + 1));
	col = (unsigned int*)malloc(sizeof(unsigned int) * (4 * (size_t)tot + 1));
	val = (double*)calloc(4 * (size_t)tot + 1, sizeof(double));
	rhs = (double*)calloc(2 * (size_t)n + 1, sizeof(double));
	xyz = (double*)malloc(sizeof(double) * (6 * (size_t)n + 1));
	if (!nb || !rp || !col || !val || !rhs || !xyz)
		goto done;
	/* each vertex's neighbours and itself, sorted; its two rows the same columns (a 2 x 2 block per neighbour) */
	for (v = 0; v < n; v++)
	{
		unsigned int m = nb0[v + 1] - nb0[v], *w = nb + nb0[v], i, j;
		w[0] = v;
		if (m > 1)
			hm_ring(dm, v, w + 1, m - 1);
		for (i = 1; i < m; i++)
			for (j = i; j > 0 && w[j - 1] > w[j]; j--)
			{
				unsigned int x = w[j];
				w[j] = w[j - 1];
				w[j - 1] = x;
			}
		rp[2 * v] = 4 * nb0[v];
		rp[2 * v + 1] = 4 * nb0[v] + 2 * m;
		for (i = 0; i < m; i++)
		{
			col[rp[2 * v] + 2 * i] = col[rp[2 * v + 1] + 2 * i] = 2 * w[i];
			col[rp[2 * v] + 2 * i + 1] = col[rp[2 * v + 1] + 2 * i + 1] = 2 * w[i] + 1;
		}
		for (k = 0; k < 3; k++)
			xyz[6 * (size_t)v + k] = xyz[6 * (size_t)v + 3 + k] = dm->x[3 * (size_t)v + k];
	}
	rp[2 * n] = 4 * tot;
	for (t = 0; t < nt; t++)
	{
		qg_tri const* T = &tri[t];
		double e[3][3], nrm[3], A, g[3][3], X[3], Y[3], B[3][2], M[3][3];
		unsigned int b = T->v[0];
		for (k = 0; k < 3; k++)
			for (l = 0; l < 3; l++)
				e[k][l] = (T->X[3 * ((k + 2) % 3) + l] - T->X[3 * ((k + 1) % 3) + l]) / rho;   /* (the edge opposite corner k) */
		qg_cross3(e[2], e[0], nrm);
		A = sqrt(qg_dot(nrm, nrm));
		if (A <= 0)
			continue;
		for (l = 0; l < 3; l++)
			nrm[l] /= A;
		A *= 0.5;
		/* grad of corner k's hat: n x e_k / 2A */
		for (k = 0; k < 3; k++)
		{
			qg_cross3(nrm, e[k], g[k]);
			for (l = 0; l < 3; l++)
				g[k][l] /= 2 * A;
		}
		/* the chart's axes: its first vertex's q and n x q, turned into the triangle's plane */
		if (!qg_axis(L, b, nrm, X))
			continue;
		qg_cross3(nrm, X, Y);
		for (k = 0; k < 3; k++)
		{
			B[k][0] = A * qg_dot(g[k], X);
			B[k][1] = A * qg_dot(g[k], Y);
			for (l = 0; l < 3; l++)
				M[k][l] = A * qg_dot(g[k], g[l]);
		}
		/* uv_k = R_k d_k + c_k: blocks M_kl R_k^T R_l, right-hand sides R_k^T (B_k - sum_l M_kl c_l) */
		for (k = 0; k < 3; k++)
		{
			unsigned int vk = T->v[k], i;
			double s[2];
			if (vk == HM_NONE)
				continue;
			s[0] = B[k][0];
			s[1] = B[k][1];
			for (l = 0; l < 3; l++)
			{
				s[0] -= M[k][l] * T->t[2 * l];
				s[1] -= M[k][l] * T->t[2 * l + 1];
			}
			qg_turn(s, -T->r[k], s);
			rhs[2 * vk] += s[0];
			rhs[2 * vk + 1] += s[1];
			for (l = 0; l < 3; l++)
			{
				unsigned int vl = T->v[l], lo = nb0[vk], hi = nb0[vk + 1] - 1, mid;
				double c0[2] = { 1, 0 }, c1[2] = { 0, 1 }, R0[2], R1[2];
				if (vl == HM_NONE)
					continue;
				while (lo < hi)
				{
					mid = (lo + hi) / 2;
					if (nb[mid] < vl)
						lo = mid + 1;
					else
						hi = mid;
				}
				i = lo - nb0[vk];
				/* R(r_l - r_k)'s columns */
				qg_turn(c0, T->r[l] - T->r[k], R0);
				qg_turn(c1, T->r[l] - T->r[k], R1);
				val[rp[2 * vk] + 2 * i] += M[k][l] * R0[0];
				val[rp[2 * vk] + 2 * i + 1] += M[k][l] * R1[0];
				val[rp[2 * vk + 1] + 2 * i] += M[k][l] * R0[1];
				val[rp[2 * vk + 1] + 2 * i + 1] += M[k][l] * R1[1];
			}
		}
	}
	for (v = 0, k = 0; v < 2 * n; v++)
	{
		unsigned int p;
		if (fix[v / 2] & (1u << (v & 1)))
			continue;
		for (p = rp[v]; p < rp[v + 1]; p++)
			if (col[p] == v)
				mean += val[p];
		k++;
	}
	mean = k ? 1e-6 * mean / k : 1e-6;
	for (v = 0; v < 2 * n; v++)
	{
		unsigned int p;
		int held = (fix[v / 2] >> (v & 1)) & 1;
		for (p = rp[v]; p < rp[v + 1]; p++)
		{
			int ch = (fix[col[p] / 2] >> (col[p] & 1)) & 1;
			if (held || ch)
				val[p] = col[p] == v ? 1 : 0;
			else if (col[p] == v)
				val[p] += mean;
		}
		rhs[v] = held ? 0 : rhs[v] + mean * d[v];
	}
	st = qaws_internal_ldlt_factor(&f, 2 * n, rp, col, val, xyz, 0);
	if (st == QAWS_STATUS_OK)
		qaws_internal_ldlt_solve(&f, rhs, d);
done:
	qaws_internal_ldlt_free(&f);
	free(nb0); free(nb); free(rp); free(col); free(val); free(rhs); free(xyz);
	return st;
}

static void qg_uv(qg_tri* t, double const* d)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		double s[2] = { 0, 0 };
		if (t->v[k] != HM_NONE)
			qg_turn(d + 2 * (size_t)t->v[k], t->r[k], s);
		t->U[2 * k] = s[0] + t->t[2 * k];
		t->U[2 * k + 1] = s[1] + t->t[2 * k + 1];
	}
}

static int qg_folded(qg_tri const* t)
{
	return !(qg_cross(t->U, t->U + 2, t->U + 4) > 1e-9);
}

#define QG_RING 32

static double qg_margin(double const (*nl)[3], unsigned int n, double const* p)
{
	double s = HUGE_VAL;
	unsigned int i;
	for (i = 0; i < n; i++)
		s = fmin(s, nl[i][0] * p[0] + nl[i][1] * p[1] - nl[i][2]);
	return s;
}

/* p moved to where its smallest distance to the lines nl (n . p - c) is largest, within a step of where it was, its held
   components (fix bits) kept: where three of them are equal (two with one component held); 1 when it gained */
static int qg_centre(double const (*nl)[3], unsigned int n, unsigned int fix, double* p)
{
	double best = qg_margin(nl, n, p), q[2], p0[2];
	unsigned int i, j, k;
	int moved = 0;
	p0[0] = p[0];
	p0[1] = p[1];
	for (i = 0; i < n; i++)
		for (j = i + 1; j < n; j++)
		{
			if (fix)
			{
				/* along the free component f: equal at t */
				unsigned int f = (fix & 1) ? 1 : 0;
				double a = nl[i][f] - nl[j][f], b = (nl[j][0] * p0[0] + nl[j][1] * p0[1] - nl[j][2]) - (nl[i][0] * p0[0] + nl[i][1] * p0[1] - nl[i][2]), s;
				if (fabs(a) < 1e-12)
					continue;
				q[0] = p0[0];
				q[1] = p0[1];
				q[f] += b / a;
				if (fabs(q[f] - p0[f]) > 1)
					continue;
				s = qg_margin(nl, n, q);
				if (s > best + 1e-12)
				{
					best = s;
					p[0] = q[0];
					p[1] = q[1];
					moved = 1;
				}
				continue;
			}
			for (k = j + 1; k < n; k++)
			{
				double a00 = nl[i][0] - nl[j][0], a01 = nl[i][1] - nl[j][1], a10 = nl[i][0] - nl[k][0], a11 = nl[i][1] - nl[k][1];
				double b0 = nl[i][2] - nl[j][2], b1 = nl[i][2] - nl[k][2], det = a00 * a11 - a01 * a10, s;
				if (fabs(det) < 1e-12)
					continue;
				q[0] = (b0 * a11 - a01 * b1) / det;
				q[1] = (a00 * b1 - b0 * a10) / det;
				if (fabs(q[0] - p0[0]) > 1 || fabs(q[1] - p0[1]) > 1)
					continue;
				s = qg_margin(nl, n, q);
				if (s > best + 1e-12)
				{
					best = s;
					p[0] = q[0];
					p[1] = q[1];
					moved = 1;
				}
			}
		}
	return moved;
}

/* Folds undone after the solve:the free vertices of folded triangles and their rings (one more ring every 8 rounds),
   each moved to the mean of its neighbours seen from its own frame (a Tutte step, held components kept), at most 64
   rounds. Returns the folds left (HM_NONE on allocation failure). */
static unsigned int qg_untangle(qg_tri* tri, unsigned int nt, unsigned int nv, unsigned char const* fix, double* d)
{
	unsigned int *at = (unsigned int*)calloc((size_t)nv + 2, sizeof(unsigned int)), *ct = (unsigned int*)malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
	unsigned char* mark = (unsigned char*)malloc((size_t)nv + 1);
	unsigned int t, k, v, round, folds = HM_NONE;
	if (!at || !ct || !mark)
		goto done;
	for (t = 0; t < nt; t++)
		for (k = 0; k < 3; k++)
			if (tri[t].v[k] != HM_NONE)
				at[tri[t].v[k] + 2]++;
	for (v = 0; v < nv; v++)
		at[v + 2] += at[v + 1];
	for (t = 0; t < nt; t++)
		for (k = 0; k < 3; k++)
			if (tri[t].v[k] != HM_NONE)
				ct[at[tri[t].v[k] + 1]++] = 3 * t + k;
	for (round = 0;; round++)
	{
		unsigned int g;
		for (t = 0, folds = 0; t < nt; t++)
			folds += qg_folded(&tri[t]);
		if (!folds || round == 64)
			break;
		memset(mark, 0, nv);
		for (t = 0; t < nt; t++)
			if (qg_folded(&tri[t]))
				for (k = 0; k < 3; k++)
					if (tri[t].v[k] != HM_NONE)
						mark[tri[t].v[k]] = 1;
		for (g = 1; g <= 1 + round / 8; g++)
			for (v = 0; v < nv; v++)
				if (mark[v] == g)
				{
					unsigned int i;
					for (i = at[v]; i < at[v + 1]; i++)
						for (k = 0; k < 3; k++)
						{
							unsigned int w = tri[ct[i] / 3].v[k];
							if (w != HM_NONE && !mark[w])
								mark[w] = (unsigned char)(g + 1);
						}
				}
		for (v = 0; v < nv; v++)
		{
			double nl[QG_RING][3], p[2];
			unsigned int i, n = 0;
			if (!mark[v] || fix[v] == 3)
				continue;
			/* each triangle's far side as a line, its signed distance n . p - c positive inside, in v's frame */
			for (i = at[v]; i < at[v + 1] && n < QG_RING; i++)
			{
				qg_tri const* T = &tri[ct[i] / 3];
				unsigned int c = ct[i] % 3, j;
				double a[2], b[2], l;
				for (j = 0; j < 2; j++)
				{
					double* q = j ? b : a;
					unsigned int w = (c + 1 + j) % 3;
					q[0] = T->U[2 * w] - T->t[2 * c];
					q[1] = T->U[2 * w + 1] - T->t[2 * c + 1];
					qg_turn(q, -T->r[c], q);
				}
				l = hypot(b[0] - a[0], b[1] - a[1]);
				if (!(l > 0))
					continue;
				nl[n][0] = -(b[1] - a[1]) / l;
				nl[n][1] = (b[0] - a[0]) / l;
				nl[n][2] = nl[n][0] * a[0] + nl[n][1] * a[1];
				n++;
			}
			p[0] = d[2 * (size_t)v];
			p[1] = d[2 * (size_t)v + 1];
			if (n && qg_centre(nl, n, fix[v], p))
			{
				d[2 * (size_t)v] = p[0];
				d[2 * (size_t)v + 1] = p[1];
				for (i = at[v]; i < at[v + 1]; i++)
					qg_uv(&tri[ct[i] / 3], d);
			}
		}
	}
done:
	free(at);
	free(ct);
	free(mark);
	return folds;
}

/* corner k of a triangle charted from dense corner c's vertex: the vertex reached by the face's steps, k = 1 along c,
   k = 2 back along the corner before c */
static void qg_chart(qaws_hm const* dm, int const* T, unsigned char const* K, unsigned int c, qg_tri* t)
{
	unsigned int c1 = hm_next(dm, c), c2 = hm_next(dm, c1);
	int m[2];
	t->v[0] = dm->cv[c];
	t->v[1] = dm->cv[c1];
	t->v[2] = dm->cv[c2];
	t->r[0] = 0;
	t->t[0] = t->t[1] = 0;
	t->r[1] = (4 - K[c]) & 3;
	t->t[2] = T[2 * (size_t)c];
	t->t[3] = T[2 * (size_t)c + 1];
	t->r[2] = K[c2];
	qf_turn(T + 2 * (size_t)c2, K[c2], m);
	t->t[4] = -m[0];
	t->t[5] = -m[1];
}

/* The steps of an edge whose ends' frames lie in planes far apart (two line vertices on different creases by a corner:
   the rotation taking one normal onto the other is about neither crease) read again in its faces' plane: both ends' q
   turned into it, K the quarter turns between them, T the step between their lattice points there; kept when its two
   faces disagree (a crease between them) or when the edge is a feature. */
static void qg_mend(qaws_hm const* dm, mf_level const* L, double rho, int* T, unsigned char* K)
{
	unsigned int c, k;
	for (c = 0; c < dm->nc; c++)
	{
		unsigned int o = dm->co[c], a, b, s, fc[2];
		int kg[2], m[2];
		double Xa[3], Xb[3], Ya[3], n[2][3], x[3];
		if (dm->cf[c] == HM_NONE || (o != HM_NONE && o < c) || (dm->cflag[c] & HM_FEATURE) || (o != HM_NONE && (dm->cflag[o] & HM_FEATURE)))
			continue;
		a = dm->cv[c];
		b = hm_to(dm, c);
		if (qg_dot(L->nrm + 3 * (size_t)a, L->nrm + 3 * (size_t)b) > 0.9 || !(L->con[a] & L->con[b] & MF_ON_LINE))
			continue;
		fc[0] = dm->cf[c];
		fc[1] = o == HM_NONE ? fc[0] : dm->cf[o];
		for (s = 0; s < 2; s++)
		{
			double l;
			hm_face_normal(dm, fc[s], n[s]);
			l = sqrt(qg_dot(n[s], n[s]));
			if (!(l > 0) || !qg_axis(L, a, n[s], Xa) || !qg_axis(L, b, n[s], Xb))
				break;
			for (k = 0; k < 3; k++)
				n[s][k] /= l;
			qg_cross3(Xb, Xa, x);
			kg[s] = (int)floor(atan2(qg_dot(x, n[s]), qg_dot(Xb, Xa)) / 1.5707963267948966 + 0.5) & 3;
		}
		if (s < 2 || kg[0] != kg[1] || kg[0] == K[c])
			continue;
		qg_axis(L, a, n[0], Xa);
		qg_cross3(n[0], Xa, Ya);
		for (k = 0; k < 3; k++)
			x[k] = L->o[3 * (size_t)b + k] - L->o[3 * (size_t)a + k];
		K[c] = (unsigned char)kg[0];
		T[2 * (size_t)c] = (int)floor(qg_dot(x, Xa) / rho + 0.5);
		T[2 * (size_t)c + 1] = (int)floor(qg_dot(x, Ya) / rho + 0.5);
		if (o != HM_NONE)
		{
			m[0] = -T[2 * (size_t)c];
			m[1] = -T[2 * (size_t)c + 1];
			K[o] = (unsigned char)((4 - K[c]) & 3);
			qf_turn(m, K[c], T + 2 * (size_t)o);
		}
	}
}

qaws_status qg_quads(qaws_hm const* dm, mf_level const* L, double rho, int const* T0, unsigned char const* K0, qd_out* o, qg_map* map)
{
	qg_ctx c;
	unsigned int nf = dm->nf, f, i, k, *ct = NULL, cap;
	unsigned char *fix = NULL, *cone = NULL, *K = NULL;
	double* d = NULL;
	int* T = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&c, 0, sizeof c);
	memset(o, 0, sizeof *o);
	T = (int*)malloc(sizeof(int) * (2 * (size_t)dm->nc + 2));
	K = (unsigned char*)malloc((size_t)dm->nc + 1);
	if (!T || !K)
	{
		free(T);
		free(K);
		return st;
	}
	memcpy(T, T0, sizeof(int) * 2 * dm->nc);
	memcpy(K, K0, dm->nc);
	qg_mend(dm, L, rho, T, K);
	cap = 3 * nf + 1;
	c.tri = (qg_tri*)calloc(cap, sizeof(qg_tri));
	ct = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)dm->nc + 1));
	fix = (unsigned char*)calloc((size_t)dm->nv + 1, 1);
	cone = (unsigned char*)calloc((size_t)dm->nv + 1, 1);
	d = (double*)malloc(sizeof(double) * (2 * (size_t)dm->nv + 2));
	c.mask = 255;
	c.hash = (unsigned int*)malloc(sizeof(unsigned int) * 256);
	c.reg = (qg_reg*)malloc(sizeof(qg_reg) * 129);
	if (!c.tri || !ct || !fix || !cone || !d || !c.hash || !c.reg)
		goto done;
	for (k = 0; k < 256; k++)
		c.hash[k] = HM_NONE;
	for (k = 0; k < dm->nc; k++)
		ct[k] = HM_NONE;
	st = QAWS_STATUS_INVALID_TOPOLOGY;
	/* the triangles: a face whose frames do not turn charted from its first vertex (its steps must close); a singular one
	   round its cone, the fixed point of its steps' loop (a lattice point, or no map): from the vertex after the cone's
	   when the cone is a vertex's lattice point (that vertex then the cone), else split in three round it */
	for (f = 0; f < nf; f++)
	{
		unsigned int c0 = dm->fc[f], cc[3], v[3];
		int ks = 0, h[2] = { 0, 0 }, p[3][2], at = -1;
		for (i = 0; i < 3; i++)
		{
			cc[i] = c0 + i;
			v[i] = dm->cv[cc[i]];
			ks += K[cc[i]];
		}
		for (i = 0; i < 3; i++)
		{
			h[0] -= T[2 * (size_t)cc[i]];
			h[1] -= T[2 * (size_t)cc[i] + 1];
			qf_turn(h, K[cc[i]], h);
		}
		ks &= 3;
		for (i = 0, at = 0; i < 3; i++)
			at += (L->con[v[i]] & MF_PINNED) != 0;
		if (at)
		{
			/* a corner's frame lies in none of its faces (nor do its steps): the face charted from a vertex that is no
			   corner, each corner at its lattice point (itself) as that vertex sees it, a cone of the map */
			unsigned int b = 0, j;
			qg_tri* t = &c.tri[c.nt];
			double tq[3];
			if (at == 3)
				goto done;
			while (L->con[v[b]] & MF_PINNED)
				b++;
			qg_chart(dm, T, K, cc[b], t);
			qg_cross3(L->nrm + 3 * (size_t)v[b], L->q + 3 * (size_t)v[b], tq);
			for (j = 1; j < 3; j++)
				if (L->con[t->v[j]] & MF_PINNED)
				{
					double x[3];
					for (k = 0; k < 3; k++)
						x[k] = dm->x[3 * (size_t)t->v[j] + k] - L->o[3 * (size_t)v[b] + k];
					t->r[j] = 0;
					t->t[2 * j] = (int)floor(qg_dot(x, L->q + 3 * (size_t)v[b]) / rho + 0.5);
					t->t[2 * j + 1] = (int)floor(qg_dot(x, tq) / rho + 0.5);
				}
			for (i = 0; i < 3; i++)
			{
				ct[c0 + (b + i) % 3] = 3 * c.nt + i;
				if (L->con[v[i]] & MF_PINNED)
					cone[v[i]] = 1;
			}
			c.nt++;
			continue;
		}
		at = -1;
		if (!ks)
		{
			if (h[0] || h[1])
				goto done;   /* (a dislocation) */
			qg_chart(dm, T, K, c0, &c.tri[c.nt]);
			for (i = 0; i < 3; i++)
				ct[cc[i]] = 3 * c.nt + i;
			c.nt++;
			continue;
		}
		/* p - turn(p, ks) = h */
		if (ks == 2)
		{
			if ((h[0] & 1) || (h[1] & 1))
				goto done;
			p[0][0] = h[0] / 2;
			p[0][1] = h[1] / 2;
		}
		else
		{
			if ((h[0] + h[1]) & 1)
				goto done;
			p[0][0] = ks == 1 ? (h[0] - h[1]) / 2 : (h[0] + h[1]) / 2;
			p[0][1] = ks == 1 ? (h[0] + h[1]) / 2 : (h[1] - h[0]) / 2;
		}
		for (i = 1; i < 3; i++)
		{
			int s[2];
			s[0] = p[i - 1][0] - T[2 * (size_t)cc[i - 1]];
			s[1] = p[i - 1][1] - T[2 * (size_t)cc[i - 1] + 1];
			qf_turn(s, K[cc[i - 1]], p[i]);
		}
		for (i = 0; i < 3 && at < 0; i++)
			if (!p[i][0] && !p[i][1])
				at = (int)i;
		if (at >= 0)
		{
			unsigned int b = cc[(at + 1) % 3];
			cone[v[at]] = 1;
			fix[v[at]] = 3;
			qg_chart(dm, T, K, b, &c.tri[c.nt]);
			for (i = 0; i < 3; i++)
				ct[c0 + (b - c0 + i) % 3] = 3 * c.nt + i;
			c.nt++;
			continue;
		}
		{
			/* the cone at the face's centroid (its lattice point in space may lie outside the face; the solve and the
			   untangling place the face's corners round it in the map) */
			double x[3];
			unsigned int j, t0 = c.nt;
			for (k = 0; k < 3; k++)
				x[k] = (dm->x[3 * (size_t)v[0] + k] + dm->x[3 * (size_t)v[1] + k] + dm->x[3 * (size_t)v[2] + k]) / 3;
			for (i = 0; i < 3; i++)
			{
				qg_tri* t = &c.tri[c.nt];
				qg_chart(dm, T, K, cc[i], t);
				t->v[2] = HM_NONE;
				t->r[2] = 0;
				t->t[4] = p[i][0];
				t->t[5] = p[i][1];
				memcpy(t->X + 6, x, sizeof x);
				ct[cc[i]] = 3 * c.nt;
				/* edge 1 (to the cone) against the next third's edge 2 (from it) */
				t->nb[1] = 3 * (t0 + (i + 1) % 3) + 2;
				t->nb[2] = 3 * (t0 + (i + 2) % 3) + 1;
				t->line = 0;
				c.nt++;
			}
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
				{
					c.tri[t0 + j].X[k] = dm->x[3 * (size_t)v[j] + k];
					c.tri[t0 + j].X[3 + k] = dm->x[3 * (size_t)v[(j + 1) % 3] + k];
				}
		}
	}
	/* positions and neighbours across the dense edges */
	for (f = 0; f < nf; f++)
		for (i = 0; i < 3; i++)
		{
			unsigned int cc = dm->fc[f] + i, w = ct[cc], o2 = dm->co[cc];
			qg_tri* t = &c.tri[w / 3];
			if (t->v[2] != HM_NONE)
				for (k = 0; k < 3; k++)
					memcpy(t->X + 3 * k, dm->x + 3 * (size_t)t->v[k], sizeof(double) * 3);
			t->nb[w % 3] = o2 == HM_NONE ? HM_NONE : ct[o2];
			if (o2 == HM_NONE || ((dm->cflag[cc] | dm->cflag[o2]) & HM_FEATURE))
				t->line |= (unsigned char)(1u << (w % 3));
		}
	/* the held components: corners and cones at their lattice points, line vertices on their line's isoline */
	for (i = 0; i < dm->nv; i++)
	{
		double x[3], tq[3];
		unsigned int con = L->con[i];
		qg_cross3(L->nrm + 3 * (size_t)i, L->q + 3 * (size_t)i, tq);
		for (k = 0; k < 3; k++)
			x[k] = dm->x[3 * (size_t)i + k] - L->o[3 * (size_t)i + k];
		d[2 * i] = qg_dot(x, L->q + 3 * (size_t)i) / rho;
		d[2 * i + 1] = qg_dot(x, tq) / rho;
		if (dm->vc[i] == HM_NONE || (con & MF_PINNED))
			fix[i] = 3;
		else if (con & MF_ON_LINE)
			fix[i] |= fabs(qg_dot(L->cq + 3 * (size_t)i, L->q + 3 * (size_t)i)) >= fabs(qg_dot(L->cq + 3 * (size_t)i, tq)) ? 2 : 1;
		for (k = 0; k < 2; k++)
			if (fix[i] & (1u << k))
				d[2 * i + k] = 0;
	}
	st = qg_solve(dm, L, rho, c.tri, c.nt, fix, d);
	if (st != QAWS_STATUS_OK)
	{
		if (st != QAWS_STATUS_ALLOCATION_FAILURE)
			st = QAWS_STATUS_INVALID_TOPOLOGY;
		goto done;
	}
	/* the map from the places, its folds undone; a fold left: no map.
	   shortcut: the integers are the flow's, never re-quantized against the map, and a fold the untangling leaves (or a
	   cone wound by a full turn, caught at extraction) falls back on INSTANT; quantize on the map (greedy rounding with
	   re-solves, or Bi-MDF) when the fallbacks matter */
	for (i = 0; i < c.nt; i++)
		qg_uv(&c.tri[i], d);
	k = qg_untangle(c.tri, c.nt, dm->nv, fix, d);
	if (k == HM_NONE)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	st = QAWS_STATUS_INVALID_TOPOLOGY;
	if (k)
		goto done;
	/* the transitions, from a corner both triangles chart through the same steps (a vertex that is no cone), else from the
	   edge itself (both ends on lattice points); checked at both ends */
	for (i = 0; i < c.nt; i++)
		for (k = 0; k < 3; k++)
		{
			qg_tri* t = &c.tri[i];
			unsigned int nx = t->nb[k], j, a = k, b = (k + 1) % 3, ka, kb, r;
			qg_tri const* s;
			double ea[2], eb[2], err = 0;
			int m[2];
			if (nx == HM_NONE)
				continue;
			s = &c.tri[nx / 3];
			j = nx % 3;
			ka = (j + 1) % 3;
			kb = j;
			if (t->v[a] != HM_NONE && !cone[t->v[a]])
				r = (unsigned int)(s->r[ka] - t->r[a]) & 3;
			else if (t->v[b] != HM_NONE && !cone[t->v[b]])
			{
				r = (unsigned int)(s->r[kb] - t->r[b]) & 3;
				a = b;
				ka = kb;
			}
			else
			{
				double best = HUGE_VAL, de[2], dt[2];
				unsigned int q;
				dt[0] = t->U[2 * b] - t->U[2 * a];
				dt[1] = t->U[2 * b + 1] - t->U[2 * a + 1];
				for (r = 0, q = 0; q < 4; q++)
				{
					double e;
					qg_turn(dt, (int)q, de);
					e = fabs(de[0] - (s->U[2 * kb] - s->U[2 * ka])) + fabs(de[1] - (s->U[2 * kb + 1] - s->U[2 * ka + 1]));
					if (e < best)
					{
						best = e;
						r = q;
					}
				}
			}
			if (t->v[a] != HM_NONE && !cone[t->v[a]])
			{
				qf_turn(t->t + 2 * a, (int)r, m);
				t->tr[3 * k + 1] = s->t[2 * ka] - m[0];
				t->tr[3 * k + 2] = s->t[2 * ka + 1] - m[1];
			}
			else
			{
				qg_turn(t->U + 2 * a, (int)r, ea);
				t->tr[3 * k + 1] = (int)floor(s->U[2 * ka] - ea[0] + 0.5);
				t->tr[3 * k + 2] = (int)floor(s->U[2 * ka + 1] - ea[1] + 0.5);
			}
			t->tr[3 * k] = (int)r;
			a = k;
			b = (k + 1) % 3;
			memcpy(ea, t->U + 2 * a, sizeof ea);
			memcpy(eb, t->U + 2 * b, sizeof eb);
			qg_across(t, k, ea);
			qg_across(t, k, eb);
			err = fabs(ea[0] - s->U[2 * ((j + 1) % 3)]) + fabs(ea[1] - s->U[2 * ((j + 1) % 3) + 1]) + fabs(eb[0] - s->U[2 * j])
				+ fabs(eb[1] - s->U[2 * j + 1]);
			if (!(err < 1e-6))
				goto done;
		}
	st = qg_points(&c);
	if (st == QAWS_STATUS_OK)
		st = qg_read(&c, dm, L, o);
	if (st == QAWS_STATUS_OK && map)
	{
		map->x = (double*)malloc(sizeof(double) * 9 * ((size_t)c.nt + 1));
		map->uv = (double*)malloc(sizeof(double) * 6 * ((size_t)c.nt + 1));
		if (!map->x || !map->uv)
		{
			qg_map_free(map);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		map->nt = c.nt;
		for (i = 0; i < c.nt; i++)
		{
			memcpy(map->x + 9 * (size_t)i, c.tri[i].X, sizeof(double) * 9);
			memcpy(map->uv + 6 * (size_t)i, c.tri[i].U, sizeof(double) * 6);
		}
	}
done:
	if (st != QAWS_STATUS_OK)
		qd_out_free(o);
	free(c.tri); free(c.reg); free(c.hash); free(c.uf); free(c.home); free(c.bary);
	free(ct); free(fix); free(cone); free(d); free(T); free(K);
	return st;
}
