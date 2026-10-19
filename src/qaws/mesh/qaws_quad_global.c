/*
 * GLOBAL quads: the integer-grid map of INSTANT's dense mesh from its
 * consistent steps (or, when they give none, from steps quantized on the
 * map), solved once, and its quads read off by tracing the map's unit
 * isolines (qaws_quad_global.h).
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
#include "../internal/qaws_internal_memory.h"
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
	qaws_internal_free(m->x);
	qaws_internal_free(m->uv);
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
		unsigned int* h = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * n);
		qg_reg* g = (qg_reg*)qaws_internal_realloc(c->reg, sizeof(qg_reg) * (n / 2 + 1));
		if (g)
			c->reg = g;
		if (!h || !g)
		{
			qaws_internal_free(h);
			return 0;
		}
		for (k = 0; k < n; k++)
			h[k] = HM_NONE;
		qaws_internal_free(c->hash);
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
					unsigned int* uf = (unsigned int*)qaws_internal_realloc(c->uf, sizeof(unsigned int) * n);
					unsigned int* home = uf ? (unsigned int*)qaws_internal_realloc(c->home, sizeof(unsigned int) * n) : NULL;
					double* bary = home ? (double*)qaws_internal_realloc(c->bary, sizeof(double) * 3 * n) : NULL;
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
	he = (qg_he*)qaws_internal_malloc(sizeof(qg_he) * (nhe + 1));
	seen = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (nhe + 1));
	out = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (c->np + 1));
	o->fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (nhe + 1));
	o->fs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (nhe / 4 + 1));
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
	o->x = (double*)qaws_internal_malloc(sizeof(double) * 3 * (o->nv + 1));
	o->foot = (qaws_mesh_location*)qaws_internal_malloc(sizeof(qaws_mesh_location) * (o->nv + 1));
	o->kind = (unsigned char*)qaws_internal_malloc(o->nv + 1);
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
	qaws_internal_free(he);
	qaws_internal_free(seen);
	qaws_internal_free(out);
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

/* Triangle T's energy A |J - F|^2 (J the map's gradient, F the field of its chart's vertex turned into its plane, in
   steps) as sum over u and v of U^T M U - 2 B^T U + const, U its corners' coordinate; 0 when it is degenerate. */
static int qg_mb(mf_level const* L, double rho, qg_tri const* T, double (*M)[3], double (*B)[2])
{
	double e[3][3], nrm[3], A, g[3][3], X[3], Y[3];
	unsigned int k, l;
	for (k = 0; k < 3; k++)
		for (l = 0; l < 3; l++)
			e[k][l] = (T->X[3 * ((k + 2) % 3) + l] - T->X[3 * ((k + 1) % 3) + l]) / rho;   /* (the edge opposite corner k) */
	qg_cross3(e[2], e[0], nrm);
	A = sqrt(qg_dot(nrm, nrm));
	if (A <= 0)
		return 0;
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
	if (!qg_axis(L, T->v[0], nrm, X))
		return 0;
	qg_cross3(nrm, X, Y);
	for (k = 0; k < 3; k++)
	{
		B[k][0] = A * qg_dot(g[k], X);
		B[k][1] = A * qg_dot(g[k], Y);
		for (l = 0; l < 3; l++)
			M[k][l] = A * qg_dot(g[k], g[l]);
	}
	return 1;
}

/* The places d (2 per dense vertex) minimizing the map's energy (qg_mb), the steps fixed; held components (fix, bits 1
   and 2) kept at 0; the rest pulled toward d by a millionth of a typical weight (the torus's free shift). keep (may be
   NULL): the factor kept, for more right-hand sides. */
static qaws_status qg_solve(qaws_hm const* dm, mf_level const* L, double rho, qg_tri const* tri, unsigned int nt, unsigned char const* fix,
	double* d, qaws_internal_ldlt* keep)
{
	unsigned int n = dm->nv, v, t, k, l, *nb0 = NULL, *nb = NULL, *rp = NULL, *col = NULL, ring[64], tot = 0;
	double *val = NULL, *rhs = NULL, *xyz = NULL, mean = 0;
	qaws_internal_ldlt f;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&f, 0, sizeof f);
	nb0 = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (n + 1));
	if (!nb0)
		goto done;
	for (v = 0; v < n; v++)
	{
		nb0[v] = tot;
		tot += 1 + (dm->vc[v] == HM_NONE ? 0 : hm_ring(dm, v, ring, 64));
	}
	nb0[n] = tot;
	nb = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (tot + 1));
	rp = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (2 * (size_t)n + 1));
	col = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (4 * (size_t)tot + 1));
	val = (double*)qaws_internal_zalloc(4 * (size_t)tot + 1, sizeof(double));
	rhs = (double*)qaws_internal_zalloc(2 * (size_t)n + 1, sizeof(double));
	xyz = (double*)qaws_internal_malloc(sizeof(double) * (6 * (size_t)n + 1));
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
		double B[3][2], M[3][3];
		if (!qg_mb(L, rho, T, M, B))
			continue;
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
	if (st == QAWS_STATUS_OK && keep)
	{
		*keep = f;
		memset(&f, 0, sizeof f);
	}
done:
	qaws_internal_ldlt_free(&f);
	qaws_internal_free(nb0); qaws_internal_free(nb); qaws_internal_free(rp); qaws_internal_free(col); qaws_internal_free(val); qaws_internal_free(rhs); qaws_internal_free(xyz);
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
	unsigned int *at = (unsigned int*)qaws_internal_zalloc((size_t)nv + 2, sizeof(unsigned int)), *ct = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
	unsigned char* mark = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
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
	qaws_internal_free(at);
	qaws_internal_free(ct);
	qaws_internal_free(mark);
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
			if (!(l > 0))
				break;
			for (k = 0; k < 3; k++)
				n[s][k] /= l;
			if (!qg_axis(L, a, n[s], Xa) || !qg_axis(L, b, n[s], Xb))
				break;
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

/* the lattice point P (d's: 0 in its own frame) solving P = turn(P, R) + t; 0 when it is none (t of the wrong parity, or
   R = 0 and t not 0) */
static int qg_fixed(int R, int const* t, int* P)
{
	switch (R & 3)
	{
	case 0: return !t[0] && !t[1];
	case 2:
		if ((t[0] & 1) || (t[1] & 1))
			return 0;
		P[0] = t[0] / 2;
		P[1] = t[1] / 2;
		return 1;
	default:
		if ((t[0] + t[1]) & 1)
			return 0;
		P[0] = (R & 3) == 1 ? (t[0] - t[1]) / 2 : (t[0] + t[1]) / 2;
		P[1] = (R & 3) == 1 ? (t[0] + t[1]) / 2 : (t[1] - t[0]) / 2;
		return 1;
	}
}

/* A corner's place in each triangle round it: carried round from triangle to triangle through the vertex they share (its
   chart in both), from the fixed point of the turn all the way round (inside) or the first triangle's rounding (on the
   border); without that, each triangle's own rounding may disagree with its neighbour's by a step. A neighbour that is a
   cone or a corner itself leaves the roundings. sR (may be NULL): per vertex the turn round its corner when it goes all
   the way round, -1 for none; sbad: 1 where that turn fixes no lattice point (the roundings kept, a step off somewhere). */
static void qg_stars(qg_ctx* c, qaws_hm const* dm, mf_level const* L, unsigned int const* ct, unsigned char const* cone, signed char* sR,
	unsigned char* sbad)
{
	unsigned int v;
	for (v = 0; v < dm->nv; v++)
	{
		unsigned int t, k, n, round = 0, path[64][2], m = 0, i, ok = 1;
		int R = 0, tau[2] = { 0, 0 }, P[2];
		if (sR)
		{
			sR[v] = -1;
			sbad[v] = 0;
		}
		if (!(L->con[v] & MF_PINNED) || dm->vc[v] == HM_NONE || ct[dm->vc[v]] == HM_NONE)
			continue;
		/* back to the first triangle (across the edge into the corner) */
		t = ct[dm->vc[v]] / 3;
		k = ct[dm->vc[v]] % 3;
		for (n = 0; n < 64; n++)
		{
			unsigned int nx = c->tri[t].nb[(k + 2) % 3];
			if (nx == HM_NONE)
				break;
			t = nx / 3;
			k = nx % 3;
			if (t == ct[dm->vc[v]] / 3)
			{
				round = 1;
				break;
			}
		}
		/* forward across the edges out of the corner, each triangle's place as turn(P0, R) + tau */
		for (n = 0; ok; n++)
		{
			qg_tri const* T = &c->tri[t];
			unsigned int nx = T->nb[k], s, j, w;
			int q[2];
			path[m][0] = t;
			path[m][1] = k;
			m++;
			if (nx == HM_NONE)
				break;
			s = nx / 3;
			j = nx % 3;
			w = T->v[(k + 1) % 3];
			if (w == HM_NONE || cone[w] || c->tri[s].v[j] != w || m == 64)
			{
				ok = 0;
				break;
			}
			/* uv in s = turn(uv in t - w's t there, its r in s - in t) + w's t in s */
			q[0] = tau[0] - T->t[2 * ((k + 1) % 3)];
			q[1] = tau[1] - T->t[2 * ((k + 1) % 3) + 1];
			qf_turn(q, c->tri[s].r[j] - T->r[(k + 1) % 3], tau);
			tau[0] += c->tri[s].t[2 * j];
			tau[1] += c->tri[s].t[2 * j + 1];
			R = (R + c->tri[s].r[j] - T->r[(k + 1) % 3]) & 3;
			t = s;
			k = (j + 1) % 3;
			if (round && t == path[0][0])
				break;
		}
		if (!ok)
			continue;
		P[0] = c->tri[path[0][0]].t[2 * path[0][1]];
		P[1] = c->tri[path[0][0]].t[2 * path[0][1] + 1];
		if (round && sR)
			sR[v] = (signed char)R;
		if (round && !qg_fixed(R, tau, P))
		{
			if (sR)
				sbad[v] = 1;
			continue;
		}
		/* carried round again, from P */
		for (i = 0; i < m; i++)
		{
			qg_tri* T = &c->tri[path[i][0]];
			unsigned int a = path[i][1];
			T->t[2 * a] = P[0];
			T->t[2 * a + 1] = P[1];
			if (i + 1 < m)
			{
				unsigned int b = (a + 1) % 3, s = path[i + 1][0], j = (path[i + 1][1] + 2) % 3;
				int q[2];
				q[0] = P[0] - T->t[2 * b];
				q[1] = P[1] - T->t[2 * b + 1];
				qf_turn(q, c->tri[s].r[j] - T->r[b], P);
				P[0] += c->tri[s].t[2 * j];
				P[1] += c->tri[s].t[2 * j + 1];
			}
		}
	}
}

/* what a build reports for the quantization (qg_build in split mode) */
typedef struct qg_info
{
	unsigned char* fbad;        /* per face: a singular one whose cone is no lattice point */
	signed char* sR;            /* per vertex: the turn round a corner going all the way round (qg_stars), -1 for none */
	unsigned char* sbad;        /* per vertex: no lattice point fixed round that corner */
} qg_info;

/* The triangles of the map from the steps T, K into c (its triangles allocated for 3 per face), ct (per dense corner its
   triangle * 3 + its corner there), fix, cone and d (the places from INSTANT's lattice points). A face whose frames do not
   turn is charted from its first vertex (its steps must close); a singular one round its cone, the fixed point of its
   steps' loop (a lattice point, or no map): from the vertex after the cone's when the cone is a vertex's lattice point
   (that vertex then the cone), else split in three round it. info (split mode, for the quantization): every singular face
   split, the triangles then the same whatever T, a cone off the lattice marked in info->fbad instead of failing. */
static qaws_status qg_build(qaws_hm const* dm, mf_level const* L, double rho, int const* T, unsigned char const* K, qg_info const* info, qg_ctx* c,
	unsigned int* ct, unsigned char* fix, unsigned char* cone, double* d)
{
	unsigned int nf = dm->nf, f, i, k;
	c->nt = 0;
	memset(c->tri, 0, sizeof(qg_tri) * (3 * (size_t)nf + 1));
	memset(fix, 0, dm->nv);
	memset(cone, 0, dm->nv);
	for (k = 0; k < dm->nc; k++)
		ct[k] = HM_NONE;
	for (f = 0; f < nf; f++)
	{
		unsigned int c0 = dm->fc[f], cc[3], v[3];
		int ks = 0, h[2] = { 0, 0 }, p[3][2], at = -1;
		if (info)
			info->fbad[f] = 0;
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
			   corner, each corner at its lattice point (itself) as that vertex sees it, a cone of the map (qg_stars makes
			   the triangles round it agree) */
			unsigned int b = 0, j;
			qg_tri* t = &c->tri[c->nt];
			double tq[3];
			if (at == 3)
				return QAWS_STATUS_INVALID_TOPOLOGY;
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
				ct[c0 + (b + i) % 3] = 3 * c->nt + i;
				if (L->con[v[i]] & MF_PINNED)
					cone[v[i]] = 1;
			}
			c->nt++;
			continue;
		}
		at = -1;
		if (!ks)
		{
			if (h[0] || h[1])
				return QAWS_STATUS_INVALID_TOPOLOGY;   /* (a dislocation) */
			qg_chart(dm, T, K, c0, &c->tri[c->nt]);
			for (i = 0; i < 3; i++)
				ct[cc[i]] = 3 * c->nt + i;
			c->nt++;
			continue;
		}
		/* p - turn(p, ks) = h */
		if (!qg_fixed(ks, h, p[0]))
		{
			if (!info)
				return QAWS_STATUS_INVALID_TOPOLOGY;
			info->fbad[f] = 1;
			p[0][0] = h[0] / 2;
			p[0][1] = h[1] / 2;
		}
		for (i = 1; i < 3; i++)
		{
			int s[2];
			s[0] = p[i - 1][0] - T[2 * (size_t)cc[i - 1]];
			s[1] = p[i - 1][1] - T[2 * (size_t)cc[i - 1] + 1];
			qf_turn(s, K[cc[i - 1]], p[i]);
		}
		for (i = 0; i < 3 && at < 0 && !info; i++)
			if (!p[i][0] && !p[i][1])
				at = (int)i;
		if (at >= 0)
		{
			unsigned int b = cc[(at + 1) % 3];
			cone[v[at]] = 1;
			fix[v[at]] = 3;
			qg_chart(dm, T, K, b, &c->tri[c->nt]);
			for (i = 0; i < 3; i++)
				ct[c0 + (b - c0 + i) % 3] = 3 * c->nt + i;
			c->nt++;
			continue;
		}
		{
			/* the cone at the face's centroid (its lattice point in space may lie outside the face; the solve and the
			   untangling place the face's corners round it in the map) */
			double x[3];
			unsigned int j, t0 = c->nt;
			for (k = 0; k < 3; k++)
				x[k] = (dm->x[3 * (size_t)v[0] + k] + dm->x[3 * (size_t)v[1] + k] + dm->x[3 * (size_t)v[2] + k]) / 3;
			for (i = 0; i < 3; i++)
			{
				qg_tri* t = &c->tri[c->nt];
				qg_chart(dm, T, K, cc[i], t);
				t->v[2] = HM_NONE;
				t->r[2] = 0;
				t->t[4] = p[i][0];
				t->t[5] = p[i][1];
				memcpy(t->X + 6, x, sizeof x);
				ct[cc[i]] = 3 * c->nt;
				/* edge 1 (to the cone) against the next third's edge 2 (from it) */
				t->nb[1] = 3 * (t0 + (i + 1) % 3) + 2;
				t->nb[2] = 3 * (t0 + (i + 2) % 3) + 1;
				t->line = 0;
				c->nt++;
			}
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
				{
					c->tri[t0 + j].X[k] = dm->x[3 * (size_t)v[j] + k];
					c->tri[t0 + j].X[3 + k] = dm->x[3 * (size_t)v[(j + 1) % 3] + k];
				}
		}
	}
	/* positions and neighbours across the dense edges */
	for (f = 0; f < nf; f++)
		for (i = 0; i < 3; i++)
		{
			unsigned int cc = dm->fc[f] + i, w = ct[cc], o2 = dm->co[cc];
			qg_tri* t = &c->tri[w / 3];
			if (t->v[2] != HM_NONE)
				for (k = 0; k < 3; k++)
					memcpy(t->X + 3 * k, dm->x + 3 * (size_t)t->v[k], sizeof(double) * 3);
			t->nb[w % 3] = o2 == HM_NONE ? HM_NONE : ct[o2];
			if (o2 == HM_NONE || ((dm->cflag[cc] | dm->cflag[o2]) & HM_FEATURE))
				t->line |= (unsigned char)(1u << (w % 3));
		}
	qg_stars(c, dm, L, ct, cone, info ? info->sR : NULL, info ? info->sbad : NULL);
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
	return QAWS_STATUS_OK;
}

/* Each singular face's cone onto one of its vertices, as MIQ has its cones at vertices: the free vertex (on no line, no
   corner, no other cone) nearest it in the map with every singular face split, its lattice point moved onto the cone
   (T changed round it by that whole shift, every face kept closed), so that qg_build makes it the cone. */
static qaws_status qg_cones_onto(qaws_hm const* dm, mf_level const* L, double rho, int* T, unsigned char const* K)
{
	qg_ctx c;
	qg_info info;
	unsigned int f, i, k, *ct = NULL;
	unsigned char *fix = NULL, *cone = NULL, *used = NULL;
	double* d = NULL;
	int* s = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&c, 0, sizeof c);
	memset(&info, 0, sizeof info);
	c.tri = (qg_tri*)qaws_internal_malloc(sizeof(qg_tri) * (3 * (size_t)dm->nf + 1));
	ct = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)dm->nc + 1));
	fix = (unsigned char*)qaws_internal_malloc((size_t)dm->nv + 1);
	cone = (unsigned char*)qaws_internal_malloc((size_t)dm->nv + 1);
	used = (unsigned char*)qaws_internal_zalloc((size_t)dm->nv + 1, 1);
	d = (double*)qaws_internal_malloc(sizeof(double) * (2 * (size_t)dm->nv + 2));
	s = (int*)qaws_internal_zalloc(2 * (size_t)dm->nv + 2, sizeof(int));
	info.fbad = (unsigned char*)qaws_internal_malloc((size_t)dm->nf + 1);
	info.sR = (signed char*)qaws_internal_malloc((size_t)dm->nv + 1);
	info.sbad = (unsigned char*)qaws_internal_malloc((size_t)dm->nv + 1);
	if (!c.tri || !ct || !fix || !cone || !used || !d || !s || !info.fbad || !info.sR || !info.sbad)
		goto done;
	st = qg_build(dm, L, rho, T, K, &info, &c, ct, fix, cone, d);
	if (st == QAWS_STATUS_OK)
		st = qg_solve(dm, L, rho, c.tri, c.nt, fix, d, NULL);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (f = 0; f < dm->nf; f++)
	{
		unsigned int c0 = dm->fc[f], best = 3;
		int ks = (K[c0] + K[c0 + 1] + K[c0 + 2]) & 3, h[2] = { 0, 0 }, p[3][2];
		double dist = HUGE_VAL;
		if (!ks || info.fbad[f])
			continue;
		for (i = 0; i < 3; i++)
			if (L->con[dm->cv[c0 + i]] & MF_PINNED)
				break;
		if (i < 3)
			continue;
		for (i = 0; i < 3; i++)
		{
			h[0] -= T[2 * (size_t)(c0 + i)];
			h[1] -= T[2 * (size_t)(c0 + i) + 1];
			qf_turn(h, K[c0 + i], h);
		}
		qg_fixed(ks, h, p[0]);
		for (i = 1; i < 3; i++)
		{
			int q[2];
			q[0] = p[i - 1][0] - T[2 * (size_t)(c0 + i - 1)];
			q[1] = p[i - 1][1] - T[2 * (size_t)(c0 + i - 1) + 1];
			qf_turn(q, K[c0 + i - 1], p[i]);
		}
		for (i = 0; i < 3; i++)
		{
			unsigned int v = dm->cv[c0 + i];
			double e = hypot(p[i][0] - d[2 * (size_t)v], p[i][1] - d[2 * (size_t)v + 1]);
			if (!fix[v] && !used[v] && e < dist)
			{
				dist = e;
				best = i;
			}
		}
		if (best == 3)
			continue;
		used[dm->cv[c0 + best]] = 1;
		s[2 * (size_t)dm->cv[c0 + best]] = p[best][0];
		s[2 * (size_t)dm->cv[c0 + best] + 1] = p[best][1];
	}
	/* a vertex's lattice point moved by s (its frame): the steps out of it less s, those into it plus s seen from their start */
	for (k = 0; k < dm->nc; k++)
	{
		unsigned int a, b;
		int m[2];
		if (dm->cf[k] == HM_NONE)
			continue;
		a = dm->cv[k];
		b = hm_to(dm, k);
		qf_turn(s + 2 * (size_t)b, -(int)K[k], m);
		T[2 * (size_t)k] += m[0] - s[2 * (size_t)a];
		T[2 * (size_t)k + 1] += m[1] - s[2 * (size_t)a + 1];
	}
done:
	qaws_internal_free(c.tri); qaws_internal_free(ct); qaws_internal_free(fix); qaws_internal_free(cone); qaws_internal_free(used); qaws_internal_free(d); qaws_internal_free(s); qaws_internal_free(info.fbad); qaws_internal_free(info.sR); qaws_internal_free(info.sbad);
	return st;
}

/* the map of the steps T, K and its quads (qg_quads) */
static qaws_status qg_map_of(qaws_hm const* dm, mf_level const* L, double rho, int const* T, unsigned char const* K, qd_out* o, qg_map* map)
{
	qg_ctx c;
	unsigned int i, k, *ct = NULL;
	unsigned char *fix = NULL, *cone = NULL;
	double* d = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&c, 0, sizeof c);
	memset(o, 0, sizeof *o);
	c.tri = (qg_tri*)qaws_internal_malloc(sizeof(qg_tri) * (3 * (size_t)dm->nf + 1));
	ct = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)dm->nc + 1));
	fix = (unsigned char*)qaws_internal_malloc((size_t)dm->nv + 1);
	cone = (unsigned char*)qaws_internal_malloc((size_t)dm->nv + 1);
	d = (double*)qaws_internal_malloc(sizeof(double) * (2 * (size_t)dm->nv + 2));
	c.mask = 255;
	c.hash = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * 256);
	c.reg = (qg_reg*)qaws_internal_malloc(sizeof(qg_reg) * 129);
	if (!c.tri || !ct || !fix || !cone || !d || !c.hash || !c.reg)
		goto done;
	for (k = 0; k < 256; k++)
		c.hash[k] = HM_NONE;
	st = qg_build(dm, L, rho, T, K, NULL, &c, ct, fix, cone, d);
	if (st != QAWS_STATUS_OK)
		goto done;
	st = qg_solve(dm, L, rho, c.tri, c.nt, fix, d, NULL);
	if (st != QAWS_STATUS_OK)
	{
		if (st != QAWS_STATUS_ALLOCATION_FAILURE)
			st = QAWS_STATUS_INVALID_TOPOLOGY;
		goto done;
	}
	/* the map from the places, its folds undone; a fold left: no map */
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
	/* every free vertex inside wound once round: one wound twice (its triangles all positive, beside a cone whose map
	   winds a full turn short) leaves no map */
	{
		double* ang = (double*)qaws_internal_zalloc((size_t)dm->nv + 1, sizeof(double));
		unsigned int w = HM_NONE;
		if (!ang)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (i = 0; i < c.nt; i++)
			for (k = 0; k < 3; k++)
				if (c.tri[i].v[k] != HM_NONE)
				{
					double const* u = c.tri[i].U;
					unsigned int a = (k + 1) % 3, b = (k + 2) % 3;
					double e1[2], e2[2];
					e1[0] = u[2 * a] - u[2 * k];
					e1[1] = u[2 * a + 1] - u[2 * k + 1];
					e2[0] = u[2 * b] - u[2 * k];
					e2[1] = u[2 * b + 1] - u[2 * k + 1];
					ang[c.tri[i].v[k]] += atan2(e1[0] * e2[1] - e1[1] * e2[0], e1[0] * e2[0] + e1[1] * e2[1]);
				}
		for (i = 0; i < dm->nv && w == HM_NONE; i++)
			if (!fix[i] && !cone[i] && dm->vc[i] != HM_NONE && fabs(ang[i] - 6.283185307179586) > 1.5707963267948966)
				w = i;
		qaws_internal_free(ang);
		if (w != HM_NONE)
			goto done;
	}
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
		map->x = (double*)qaws_internal_malloc(sizeof(double) * 9 * ((size_t)c.nt + 1));
		map->uv = (double*)qaws_internal_malloc(sizeof(double) * 6 * ((size_t)c.nt + 1));
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
	qaws_internal_free(c.tri); qaws_internal_free(c.reg); qaws_internal_free(c.hash); qaws_internal_free(c.uf); qaws_internal_free(c.home); qaws_internal_free(c.bary);
	qaws_internal_free(ct); qaws_internal_free(fix); qaws_internal_free(cone); qaws_internal_free(d);
	return st;
}

/* ---- the quantization ----------------------------------------------------- */

#define QG_WALL (HM_NONE - 1)   /* a face no path crosses (round a corner on the border, or one qg_stars leaves) */
#define QG_VARS 192             /* the most integer unknowns quantized */

/* K's quarter turns from corner c to the end of its face (the loop of qg_build's h from c on) */
static int qg_suffix(qaws_hm const* dm, unsigned char const* K, unsigned int c)
{
	unsigned int e = dm->fc[dm->cf[c]] + 3;
	int s = 0;
	for (; c < e; c++)
		s += K[c];
	return s;
}

/* T changed by delta on corner c's edge (in its first vertex's frame) and on the way back to the node its path leads to
   (par: per face the corner toward it), each face between kept closed: the change into a face cancelled by the one out of
   it toward its parent. c on the border (a path from the border, delta along it) changes only its own face on the way.
   T NULL: nothing changed; last, lastD (may be NULL): the node's corner changed and by what. */
static void qg_path(qaws_hm const* dm, unsigned char const* K, unsigned int c, unsigned int const* par, unsigned int const* nodeof, int const* delta,
	int* T, unsigned int* last, int* lastD)
{
	int D[2], m[2];
	D[0] = delta[0];
	D[1] = delta[1];
	for (;;)
	{
		unsigned int o = dm->co[c], f, x;
		if (T)
		{
			T[2 * (size_t)c] += D[0];
			T[2 * (size_t)c + 1] += D[1];
		}
		if (o != HM_NONE)
		{
			m[0] = -D[0];
			m[1] = -D[1];
			qf_turn(m, K[c], D);
			if (T)
			{
				T[2 * (size_t)o] += D[0];
				T[2 * (size_t)o + 1] += D[1];
			}
		}
		else
			o = c;
		f = dm->cf[o];
		if (nodeof[f] != HM_NONE)
		{
			if (last)
			{
				*last = o;
				lastD[0] = D[0];
				lastD[1] = D[1];
			}
			return;
		}
		/* turn(dT_o, S_o) + turn(dT_x, S_x) = 0 keeps the face's loop */
		x = par[f];
		qf_turn(D, qg_suffix(dm, K, o) - qg_suffix(dm, K, x), m);
		D[0] = -m[0];
		D[1] = -m[1];
		c = x;
	}
}

/* A x = b by Cholesky, A n x n dense symmetric positive (overwritten) */
static int qg_chol(double* A, unsigned int n, double* b)
{
	unsigned int i, j, k;
	for (j = 0; j < n; j++)
	{
		double s = A[j * n + j];
		for (k = 0; k < j; k++)
			s -= A[j * n + k] * A[j * n + k];
		if (!(s > 0))
			return 0;
		A[j * n + j] = sqrt(s);
		for (i = j + 1; i < n; i++)
		{
			double t = A[i * n + j];
			for (k = 0; k < j; k++)
				t -= A[i * n + k] * A[j * n + k];
			A[i * n + j] = t / A[j * n + j];
		}
	}
	for (i = 0; i < n; i++)
	{
		for (k = 0; k < i; k++)
			b[i] -= A[i * n + k] * b[k];
		b[i] /= A[i * n + i];
	}
	for (i = n; i-- > 0;)
	{
		for (k = i + 1; k < n; k++)
			b[i] -= A[k * n + i] * b[k];
		b[i] /= A[i * n + i];
	}
	return 1;
}

/* one integer unknown: a whole step of its node's path, the triangles' translations it changes (dt, 6 per triangle) */
typedef struct qg_var
{
	unsigned int node, first, count;    /* its triangles: at first in the lists, count of them */
	unsigned int entry;                 /* the corner its path starts from */
	int delta[2];
} qg_var;

typedef struct qg_node
{
	unsigned int entry, parent;         /* the corner its path leaves by (in one of its faces), the node it leads to (HM_NONE: a root) */
	unsigned int f0, nf;                /* its faces (in the CSR of faces) */
	unsigned int star;                  /* the corner it is round, or HM_NONE (a singular face) */
	unsigned int bnd[2];                /* border edges its paths to the border start from (HM_NONE: none) */
	int R;
} qg_node;

/* The integers quantized on the map (MIQ's greedy rounding, Bommes, Zimmer, Kobbelt 2009, on a reduced system). The
   cones (singular faces, and the corners going all the way round) are the nodes of a spanning tree over the faces, its
   paths crossing no feature or border edge. The unknowns: a whole step of a node's path (two per node but the roots;
   steps that keep both ends on lattice points) moves the node against its parent, each face between kept closed; a
   singular face's path to the nearest border edge (two steps along it) moves it against the border, one unknown per
   direction. The map's energy E(z) over them: the triangles' translations change along the paths only, the places by
   one solve each on the base map's factor, so that Q_ij = dt_i^T M U_j and g_i = dt_i^T (M U_0 - B) over each path's
   triangles (the places stationary); z rounded greedily (the closest to whole first, with those within a tenth, the
   rest solved again), then whole steps down E from the rounding or from the flow's integers, whichever is lower. T
   changed in place; INVALID_TOPOLOGY (T as given) when it cannot: no cone, a dislocation, a cone off the lattice, more
   than QG_VARS unknowns.
   shortcut: the isovalues of feature lines (creases) and the torus's periods stay the flow's, paths cross no crease;
   quantize those too (a node per line, paths along it) when parallel creases come out a step too close or too far */
static qaws_status qg_quantize(qaws_hm const* dm, mf_level const* L, double rho, int* T, unsigned char const* K)
{
	unsigned int nf = dm->nf, nv = dm->nv, f, i, j, k, l, nn = 0, nz = 0, nl = 0, capl = 0, head, tail;
	unsigned int* dep = NULL;
	unsigned int *nodeof = NULL, *par = NULL, *own = NULL, *queue = NULL, *faces = NULL, *star = NULL, *order = NULL, *ct = NULL, *ct2 = NULL, *at = NULL;
	unsigned char *fix = NULL, *cone = NULL, *fix2 = NULL, *cone2 = NULL, *vis = NULL;
	int *Tw = NULL, (*dt)[6] = NULL;
	unsigned int* lt = NULL;
	double *d = NULL, *d2 = NULL, *D = NULL, *Q = NULL, *g = NULL, *z = NULL, *A = NULL, *b = NULL, *rhs = NULL;
	unsigned char* fixed = NULL;
	unsigned int* fr = NULL;
	qg_node* node = NULL;
	qg_var* var = NULL;
	qg_ctx c, c2;
	qg_info info;
	qaws_internal_ldlt fac;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&c, 0, sizeof c);
	memset(&c2, 0, sizeof c2);
	memset(&fac, 0, sizeof fac);
	memset(&info, 0, sizeof info);
	c.tri = (qg_tri*)qaws_internal_malloc(sizeof(qg_tri) * (3 * (size_t)nf + 1));
	c2.tri = (qg_tri*)qaws_internal_malloc(sizeof(qg_tri) * (3 * (size_t)nf + 1));
	ct = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)dm->nc + 1));
	ct2 = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)dm->nc + 1));
	fix = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
	cone = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
	fix2 = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
	cone2 = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
	d = (double*)qaws_internal_malloc(sizeof(double) * (2 * (size_t)nv + 2));
	d2 = (double*)qaws_internal_malloc(sizeof(double) * (2 * (size_t)nv + 2));
	info.fbad = (unsigned char*)qaws_internal_malloc((size_t)nf + 1);
	info.sR = (signed char*)qaws_internal_malloc((size_t)nv + 1);
	info.sbad = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
	Tw = (int*)qaws_internal_malloc(sizeof(int) * (2 * (size_t)dm->nc + 2));
	nodeof = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	par = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	own = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	queue = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	faces = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	vis = (unsigned char*)qaws_internal_zalloc((size_t)nf + 1, 1);
	dep = (unsigned int*)qaws_internal_zalloc((size_t)nf + 1, sizeof(unsigned int));
	star = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nv + 1));
	if (!c.tri || !c2.tri || !ct || !ct2 || !fix || !cone || !fix2 || !cone2 || !d || !d2 || !info.fbad || !info.sR || !info.sbad || !Tw || !nodeof || !par
		|| !own || !queue || !faces || !vis || !star || !dep)
		goto done;
	st = qg_build(dm, L, rho, T, K, &info, &c, ct, fix, cone, d);
	if (st != QAWS_STATUS_OK)
		goto done;
	/* the nodes: corners going all the way round (their faces), singular faces; faces round other corners walls */
	for (i = 0; i < nv; i++)
		star[i] = HM_NONE;
	for (f = 0; f < nf; f++)
	{
		unsigned int c0 = dm->fc[f], p = HM_NONE;
		int ks = K[c0] + K[c0 + 1] + K[c0 + 2];
		for (i = 0; i < 3 && p == HM_NONE; i++)
			if (L->con[dm->cv[c0 + i]] & MF_PINNED)
				p = dm->cv[c0 + i];
		nodeof[f] = HM_NONE;
		if (p != HM_NONE && info.sR[p] < 0)
			nodeof[f] = QG_WALL;
		else if (p != HM_NONE)
		{
			if (star[p] == HM_NONE)
				star[p] = nn++;
			nodeof[f] = star[p];
		}
		else if (ks & 3)
			nodeof[f] = nn++;
	}
	st = QAWS_STATUS_INVALID_TOPOLOGY;
	if (!nn || 4 * nn > QG_VARS)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	node = (qg_node*)qaws_internal_zalloc(nn, sizeof(qg_node));
	order = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * nn);
	if (!node || !order)
		goto done;
	for (i = 0; i < nn; i++)
		node[i].star = HM_NONE;
	for (f = 0; f < nf; f++)
		if (nodeof[f] < nn)
		{
			unsigned int c0 = dm->fc[f], x = nodeof[f];
			node[x].nf++;
			node[x].R = (K[c0] + K[c0 + 1] + K[c0 + 2]) & 3;
			for (i = 0; i < 3; i++)
				if (L->con[dm->cv[c0 + i]] & MF_PINNED)
				{
					node[x].star = dm->cv[c0 + i];
					node[x].R = info.sR[dm->cv[c0 + i]];
					break;
				}
		}
	for (i = 0, k = 0; i < nn; i++)
	{
		node[i].f0 = k;
		k += node[i].nf;
		node[i].nf = 0;
	}
	for (f = 0; f < nf; f++)
		if (nodeof[f] < nn)
			faces[node[nodeof[f]].f0 + node[nodeof[f]].nf++] = f;
	/* the spanning tree: breadth first from a root (a node of odd turn first), across edges that are no line, a node's path
	   to the node whose faces reached it */
	for (f = 0; f < nf; f++)
		par[f] = own[f] = HM_NONE;
	for (l = 0, k = 0; l < 2; l++)
		for (i = 0; i < nn; i++)
		{
			unsigned int r = i, q;
			if (vis[faces[node[r].f0]] || (unsigned int)((node[r].R & 1) == 0) != l)
				continue;
			node[r].parent = HM_NONE;
			node[r].entry = HM_NONE;
			order[k++] = r;
			head = tail = 0;
			for (q = 0; q < node[r].nf; q++)
			{
				vis[faces[node[r].f0 + q]] = 1;
				own[faces[node[r].f0 + q]] = r;
				queue[tail++] = faces[node[r].f0 + q];
			}
			while (head < tail)
			{
				unsigned int fh = queue[head++], e;
				for (e = dm->fc[fh]; e < dm->fc[fh] + 3; e++)
				{
					unsigned int o = dm->co[e], gf, x;
					if (o == HM_NONE || ((dm->cflag[e] | dm->cflag[o]) & HM_FEATURE))
						continue;
					gf = dm->cf[o];
					if (vis[gf] || nodeof[gf] == QG_WALL)
						continue;
					if (nodeof[gf] == HM_NONE)
					{
						vis[gf] = 1;
						par[gf] = o;
						dep[gf] = dep[fh] + 1;
						own[gf] = own[fh];
						queue[tail++] = gf;
						continue;
					}
					x = nodeof[gf];
					node[x].parent = own[fh];
					node[x].entry = o;
					order[k++] = x;
					for (q = 0; q < node[x].nf; q++)
					{
						vis[faces[node[x].f0 + q]] = 1;
						own[faces[node[x].f0 + q]] = x;
						queue[tail++] = faces[node[x].f0 + q];
					}
				}
			}
		}
	/* every node on the lattice, else no quantization.
	   shortcut: a node off the lattice (a corner whose turn fixes no lattice point) gives up; step it along its path by an
	   odd step (leaves first) when one shows up */
	st = QAWS_STATUS_INVALID_TOPOLOGY;
	for (i = 0; i < nn; i++)
		if (node[i].star != HM_NONE ? info.sbad[node[i].star] : info.fbad[faces[node[i].f0]])
			goto done;
	st = qg_solve(dm, L, rho, c.tri, c.nt, fix, d, &fac);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < c.nt; i++)
		qg_uv(&c.tri[i], d);
	/* the unknowns and the translations each changes */
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	var = (qg_var*)qaws_internal_malloc(sizeof(qg_var) * (4 * (size_t)nn + 1));
	at = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)c.nt + 1));
	if (!var || !at)
		goto done;
	for (i = 0; i < c.nt; i++)
		at[i] = HM_NONE;
	/* a singular face's paths to the border (two steps along it at the end, the cone then on the same parity): the
	   nearest border edge for each direction it moves the cone, from the faces its tree reached */
	for (i = 0; i < nn; i++)
		node[i].bnd[0] = node[i].bnd[1] = HM_NONE;
	for (f = 0; f < nf; f++)
		if (vis[f] && nodeof[f] == HM_NONE && node[own[f]].star == HM_NONE)
			for (k = dm->fc[f]; k < dm->fc[f] + 3; k++)
			{
				unsigned int x = own[f], a = dm->cv[k], last, cl;
				int e[2], h[2];
				if (dm->co[k] != HM_NONE || (fix[a] != 1 && fix[a] != 2))
					continue;
				e[0] = fix[a] == 2 ? 2 : 0;
				e[1] = fix[a] == 1 ? 2 : 0;
				qg_path(dm, K, k, par, nodeof, e, NULL, &last, h);
				qf_turn(h, qg_suffix(dm, K, last), e);
				cl = e[0] != 0 ? 0 : 1;
				if (node[x].bnd[cl] == HM_NONE || dep[dm->cf[node[x].bnd[cl]]] > dep[f])
					node[x].bnd[cl] = k;
			}
	for (i = 0; i < nn; i++)
	{
		unsigned int x = order[i], s;
		int even = node[x].parent != HM_NONE && (node[x].R == 2 || node[node[x].parent].R == 2);
		for (s = 0; s < 4; s++)
		{
			qg_var* w = &var[nz];
			if (s < 2 && node[x].parent == HM_NONE)
				continue;
			if (s >= 2 && node[x].bnd[s - 2] == HM_NONE)
				continue;
			w->node = x;
			if (s < 2)
			{
				w->entry = node[x].entry;
				w->delta[0] = even ? 2 * (1 - (int)s) : 1;
				w->delta[1] = even ? 2 * (int)s : (s ? -1 : 1);
			}
			else
			{
				unsigned int a = dm->cv[node[x].bnd[s - 2]];
				w->entry = node[x].bnd[s - 2];
				w->delta[0] = fix[a] == 2 ? 2 : 0;
				w->delta[1] = fix[a] == 1 ? 2 : 0;
			}
			w->first = nl;
			memcpy(Tw, T, sizeof(int) * 2 * dm->nc);
			qg_path(dm, K, w->entry, par, nodeof, w->delta, Tw, NULL, NULL);
			st = qg_build(dm, L, rho, Tw, K, &info, &c2, ct2, fix2, cone2, d2);
			if (st != QAWS_STATUS_OK)
				goto done;
			st = QAWS_STATUS_INVALID_TOPOLOGY;
			if (c2.nt != c.nt)
				goto done;
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (j = 0; j < c.nt; j++)
			{
				int e[6], any = 0;
				for (k = 0; k < 6; k++)
					any |= e[k] = c2.tri[j].t[k] - c.tri[j].t[k];
				if (!any)
					continue;
				if (nl == capl)
				{
					unsigned int nc2 = 2 * capl + 256;
					int (*ndt)[6] = (int (*)[6])qaws_internal_realloc(dt, sizeof(int) * 6 * nc2);
					unsigned int* nlt = ndt ? (unsigned int*)qaws_internal_realloc(lt, sizeof(unsigned int) * nc2) : NULL;
					if (ndt)
						dt = ndt;
					if (!nlt)
						goto done;
					lt = nlt;
					capl = nc2;
				}
				memcpy(dt[nl], e, sizeof e);
				lt[nl++] = j;
			}
			w->count = nl - w->first;
			nz++;
		}
	}
	st = QAWS_STATUS_INVALID_TOPOLOGY;
	if (!nz)
		goto done;
	/* E(z): per unknown its places' response (one solve on the factor), Q's column; g */
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	D = (double*)qaws_internal_malloc(sizeof(double) * (2 * (size_t)nv + 2));
	rhs = (double*)qaws_internal_malloc(sizeof(double) * (2 * (size_t)nv + 2));
	Q = (double*)qaws_internal_zalloc((size_t)nz * nz, sizeof(double));
	A = (double*)qaws_internal_malloc(sizeof(double) * (size_t)nz * nz);
	g = (double*)qaws_internal_zalloc(nz, sizeof(double));
	z = (double*)qaws_internal_zalloc(nz, sizeof(double));
	b = (double*)qaws_internal_malloc(sizeof(double) * nz);
	fixed = (unsigned char*)qaws_internal_zalloc(nz, 1);
	fr = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * nz);
	if (!D || !rhs || !Q || !A || !g || !z || !b || !fixed || !fr)
		goto done;
	for (i = 0; i < nz; i++)
		for (l = var[i].first; l < var[i].first + var[i].count; l++)
		{
			qg_tri const* t = &c.tri[lt[l]];
			double M[3][3], B[3][2];
			if (!qg_mb(L, rho, t, M, B))
				continue;
			for (k = 0; k < 3; k++)
				for (j = 0; j < 2; j++)
					g[i] += dt[l][2 * k + j] * (M[k][0] * t->U[j] + M[k][1] * t->U[2 + j] + M[k][2] * t->U[4 + j] - B[k][j]);
		}
	for (j = 0; j < nz; j++)
	{
		memset(rhs, 0, sizeof(double) * 2 * nv);
		for (l = var[j].first; l < var[j].first + var[j].count; l++)
		{
			qg_tri const* t = &c.tri[lt[l]];
			double M[3][3], B[3][2];
			at[lt[l]] = l;
			if (!qg_mb(L, rho, t, M, B))
				continue;
			for (k = 0; k < 3; k++)
				if (t->v[k] != HM_NONE)
				{
					double s[2];
					unsigned int q;
					for (q = 0; q < 2; q++)
						s[q] = -(M[k][0] * dt[l][q] + M[k][1] * dt[l][2 + q] + M[k][2] * dt[l][4 + q]);
					qg_turn(s, -t->r[k], s);
					rhs[2 * (size_t)t->v[k]] += s[0];
					rhs[2 * (size_t)t->v[k] + 1] += s[1];
				}
		}
		for (k = 0; k < 2 * nv; k++)
			if ((fix[k / 2] >> (k & 1)) & 1)
				rhs[k] = 0;
		qaws_internal_ldlt_solve(&fac, rhs, D);
		for (i = 0; i < nz; i++)
			for (l = var[i].first; l < var[i].first + var[i].count; l++)
			{
				qg_tri const* t = &c.tri[lt[l]];
				double M[3][3], B[3][2], U[6];
				unsigned int q;
				if (!qg_mb(L, rho, t, M, B))
					continue;
				for (k = 0; k < 3; k++)
				{
					U[2 * k] = U[2 * k + 1] = 0;
					if (t->v[k] != HM_NONE)
						qg_turn(D + 2 * (size_t)t->v[k], t->r[k], U + 2 * k);
					if (at[lt[l]] != HM_NONE)
					{
						U[2 * k] += dt[at[lt[l]]][2 * k];
						U[2 * k + 1] += dt[at[lt[l]]][2 * k + 1];
					}
				}
				for (k = 0; k < 3; k++)
					for (q = 0; q < 2; q++)
						Q[(size_t)i * nz + j] += dt[l][2 * k + q] * (M[k][0] * U[q] + M[k][1] * U[2 + q] + M[k][2] * U[4 + q]);
			}
		for (l = var[j].first; l < var[j].first + var[j].count; l++)
			at[lt[l]] = HM_NONE;
	}
	for (i = 0; i < nz; i++)
		for (j = 0; j < i; j++)
			Q[(size_t)i * nz + j] = Q[(size_t)j * nz + i] = 0.5 * (Q[(size_t)i * nz + j] + Q[(size_t)j * nz + i]);
	/* greedy rounding: the free unknowns solved, the closest to whole fixed (with those within a tenth of it), again */
	{
		double ridge = 0;
		st = QAWS_STATUS_INVALID_TOPOLOGY;
		for (i = 0; i < nz; i++)
			ridge = fmax(ridge, Q[(size_t)i * nz + i]);
		ridge *= 1e-9;
		for (;;)
		{
			unsigned int n = 0, a, a2;
			double best = HUGE_VAL;
			for (i = 0; i < nz; i++)
				if (!fixed[i])
					fr[n++] = i;
			if (!n)
				break;
			for (a = 0; a < n; a++)
			{
				b[a] = -g[fr[a]];
				for (i = 0; i < nz; i++)
					if (fixed[i])
						b[a] -= Q[(size_t)fr[a] * nz + i] * z[i];
				for (a2 = 0; a2 < n; a2++)
					A[(size_t)a * n + a2] = Q[(size_t)fr[a] * nz + fr[a2]] + (a == a2 ? ridge : 0);
			}
			if (!qg_chol(A, n, b))
				goto done;
			for (a = 0; a < n; a++)
			{
				z[fr[a]] = b[a];
				best = fmin(best, fabs(b[a] - floor(b[a] + 0.5)));
			}
			for (a = 0; a < n; a++)
				if (fabs(b[a] - floor(b[a] + 0.5)) <= fmax(best, 0.1))
				{
					fixed[fr[a]] = 1;
					z[fr[a]] = floor(b[a] + 0.5);
				}
		}
	}
	/* then whole steps down E(z) = z^T Q z / 2 + g^T z one unknown at a time, from the rounding or from the flow's
	   integers (z = 0), whichever is lower */
	{
		double e = 0, gain;
		unsigned int round;
		for (i = 0; i < nz; i++)
		{
			b[i] = g[i];
			for (j = 0; j < nz; j++)
				b[i] += 0.5 * Q[(size_t)i * nz + j] * z[j];
			e += z[i] * b[i];
		}
		if (e >= 0)
			memset(z, 0, sizeof(double) * nz);
		/* b: the gradient Q z + g */
		for (i = 0; i < nz; i++)
		{
			b[i] = g[i];
			for (j = 0; j < nz; j++)
				b[i] += Q[(size_t)i * nz + j] * z[j];
		}
		for (round = 0; round < 64; round++)
		{
			unsigned int best = HM_NONE;
			double s = 0;
			gain = 0;
			for (i = 0; i < nz; i++)
			{
				double q = Q[(size_t)i * nz + i], up = 0.5 * q + b[i], down = 0.5 * q - b[i];
				if (up < gain - 1e-12 * q)
				{
					gain = up;
					best = i;
					s = 1;
				}
				if (down < gain - 1e-12 * q)
				{
					gain = down;
					best = i;
					s = -1;
				}
			}
			if (best == HM_NONE)
				break;
			z[best] += s;
			for (i = 0; i < nz; i++)
				b[i] += s * Q[(size_t)i * nz + best];
		}
	}
	for (i = 0; i < nz; i++)
		if (z[i] != 0)
		{
			int dl[2];
			dl[0] = (int)z[i] * var[i].delta[0];
			dl[1] = (int)z[i] * var[i].delta[1];
			qg_path(dm, K, var[i].entry, par, nodeof, dl, T, NULL, NULL);
		}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_ldlt_free(&fac);
	qaws_internal_free(c.tri); qaws_internal_free(c2.tri); qaws_internal_free(ct); qaws_internal_free(ct2); qaws_internal_free(fix); qaws_internal_free(cone); qaws_internal_free(fix2); qaws_internal_free(cone2); qaws_internal_free(d); qaws_internal_free(d2);
	qaws_internal_free(info.fbad); qaws_internal_free(info.sR); qaws_internal_free(info.sbad); qaws_internal_free(Tw); qaws_internal_free(nodeof); qaws_internal_free(par); qaws_internal_free(own); qaws_internal_free(queue); qaws_internal_free(faces); qaws_internal_free(vis); qaws_internal_free(star);
	qaws_internal_free(node); qaws_internal_free(order); qaws_internal_free(dep); qaws_internal_free(var); qaws_internal_free(at); qaws_internal_free(dt); qaws_internal_free(lt); qaws_internal_free(D); qaws_internal_free(rhs); qaws_internal_free(Q); qaws_internal_free(A); qaws_internal_free(g); qaws_internal_free(z); qaws_internal_free(b); qaws_internal_free(fixed); qaws_internal_free(fr);
	return st;
}

qaws_status qg_quads(qaws_hm const* dm, mf_level const* L, double rho, int const* T0, unsigned char const* K0, qd_out* o, qg_map* map)
{
	int *T = (int*)qaws_internal_malloc(sizeof(int) * (2 * (size_t)dm->nc + 2)), *Tq = (int*)qaws_internal_malloc(sizeof(int) * (2 * (size_t)dm->nc + 2));
	unsigned char* K = (unsigned char*)qaws_internal_malloc((size_t)dm->nc + 1);
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(o, 0, sizeof *o);
	if (!T || !Tq || !K)
		goto done;
	memcpy(T, T0, sizeof(int) * 2 * dm->nc);
	memcpy(K, K0, dm->nc);
	qg_mend(dm, L, rho, T, K);
	/* the flow's integers first (where they give a map its quality is theirs), else the quantized ones with the cones on
	   vertices */
	st = qg_map_of(dm, L, rho, T, K, o, map);
	if (st != QAWS_STATUS_INVALID_TOPOLOGY)
		goto done;
	memcpy(Tq, T, sizeof(int) * 2 * dm->nc);
	st = qg_quantize(dm, L, rho, Tq, K);
	if (st == QAWS_STATUS_OK)
		st = qg_cones_onto(dm, L, rho, Tq, K);
	if (st == QAWS_STATUS_OK)
		st = qg_map_of(dm, L, rho, Tq, K, o, map);
done:
	qaws_internal_free(T);
	qaws_internal_free(Tq);
	qaws_internal_free(K);
	return st;
}
