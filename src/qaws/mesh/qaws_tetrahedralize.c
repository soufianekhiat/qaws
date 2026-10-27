/*
 * Tetrahedralization, method DELAUNAY (qaws_tetrahedralize.h): Bowyer-Watson
 * on an integer lattice with exact orientation and insphere predicates
 * (doubles first with an error bound, 384-bit integers when it does not
 * decide), each cavity grown until every facet of its border faces the new
 * point (so degenerate inputs never leave flat cells), the border recovered,
 * the inside flooded from it, then filled and refined with the border as
 * walls, then improved: 3-2 and 2-3 flips and the smoothing of Freitag and
 * Ollivier-Gooch, each step kept when it raises the worst dihedral sine.
 *
 * shortcut: improvement by flips and smoothing only (Klingner-Shewchuk's
 * edge removal, multi-face removal and vertex insertion would lift the
 * worst angles further) and no protection of sharp angles between border
 * triangles (the recovery may not converge there: NOT_CONVERGED); a
 * constrained Delaunay recovery (Diazzi et al. 2023) to come.
 */

#include "qaws_mesh_internal.h"
#include "../qaws_volume_improve.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_tetrahedralize.h"
#include "../qaws_mesh_bvh.h"
#include "../internal/qaws_internal_wide.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define TT_DEAD 16u
#define TT_IN 32u
#define TT_SPAN 1073741824.0            /* the lattice: 2^30 steps over the box */
#define TT_FAR 281474976710656LL        /* the super cell's scale: 2^48 */

/* facet i (opposite vertex i) of a cell with det(v1 - v0, v2 - v0, v3 - v0) > 0, outward */
static unsigned char const tt_face[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };

/* ======================================================================== */
/*  Exact predicates on lattice points                                      */
/* ======================================================================== */

static double tt_det3(double const* x, double const* y, double const* z, double* perm)
{
	*perm = fabs(x[0]) * (fabs(y[1] * z[2]) + fabs(y[2] * z[1])) + fabs(x[1]) * (fabs(y[2] * z[0]) + fabs(y[0] * z[2]))
		+ fabs(x[2]) * (fabs(y[0] * z[1]) + fabs(y[1] * z[0]));
	return x[0] * (y[1] * z[2] - y[2] * z[1]) + x[1] * (y[2] * z[0] - y[0] * z[2]) + x[2] * (y[0] * z[1] - y[1] * z[0]);
}

static qaws_wide tt_wide_det3(int64_t const* x, int64_t const* y, int64_t const* z)
{
	qaws_wide c1 = qaws_wide_cross(y[1], y[2], z[1], z[2]), c2 = qaws_wide_cross(y[2], y[0], z[2], z[0]), c3 = qaws_wide_cross(y[0], y[1], z[0], z[1]);
	qaws_wide x0 = qaws_wide_from_i64(x[0]), x1 = qaws_wide_from_i64(x[1]), x2 = qaws_wide_from_i64(x[2]);
	qaws_wide s = qaws_wide_mul(&x0, &c1), t = qaws_wide_mul(&x1, &c2), u = qaws_wide_mul(&x2, &c3);
	s = qaws_wide_add(s, &t);
	return qaws_wide_add(s, &u);
}

static qaws_wide tt_wide_len2(int64_t const* x)
{
	qaws_wide s = qaws_wide_dot(x[0], x[1], x[0], x[1]), z = qaws_wide_from_i64(x[2]), t = qaws_wide_mul(&z, &z);
	return qaws_wide_add(s, &t);
}

/* sign of det(b - a, c - a, d - a): positive when d is on the side (b - a) x (c - a) points to */
static int tt_orient(int64_t const* a, int64_t const* b, int64_t const* c, int64_t const* d)
{
	int64_t u[3], v[3], w[3];
	double fu[3], fv[3], fw[3], det, perm;
	int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = b[k] - a[k];
		v[k] = c[k] - a[k];
		w[k] = d[k] - a[k];
		fu[k] = (double)u[k];
		fv[k] = (double)v[k];
		fw[k] = (double)w[k];
	}
	det = tt_det3(fu, fv, fw, &perm);
	if (fabs(det) > 1e-14 * perm)
		return det > 0 ? 1 : -1;
	{
		qaws_wide s = tt_wide_det3(u, v, w);
		return qaws_wide_sign(&s);
	}
}

/* positive when e is inside the sphere through a, b, c, d (det(b - a, c - a, d - a) > 0) */
static int tt_insphere(int64_t const* a, int64_t const* b, int64_t const* c, int64_t const* d, int64_t const* e)
{
	int64_t A[3], B[3], C[3], D[3];
	double fa[3], fb[3], fc[3], fd[3], la, lb, lc, ld, p1, p2, p3, p4, d1, d2, d3, d4, det, perm;
	int k;
	for (k = 0; k < 3; k++)
	{
		A[k] = a[k] - e[k];
		B[k] = b[k] - e[k];
		C[k] = c[k] - e[k];
		D[k] = d[k] - e[k];
		fa[k] = (double)A[k];
		fb[k] = (double)B[k];
		fc[k] = (double)C[k];
		fd[k] = (double)D[k];
	}
	la = fa[0] * fa[0] + fa[1] * fa[1] + fa[2] * fa[2];
	lb = fb[0] * fb[0] + fb[1] * fb[1] + fb[2] * fb[2];
	lc = fc[0] * fc[0] + fc[1] * fc[1] + fc[2] * fc[2];
	ld = fd[0] * fd[0] + fd[1] * fd[1] + fd[2] * fd[2];
	d1 = tt_det3(fb, fc, fd, &p1);
	d2 = tt_det3(fa, fc, fd, &p2);
	d3 = tt_det3(fa, fb, fd, &p3);
	d4 = tt_det3(fa, fb, fc, &p4);
	det = la * d1 - lb * d2 + lc * d3 - ld * d4;
	perm = la * p1 + lb * p2 + lc * p3 + ld * p4;
	if (fabs(det) > 1e-13 * perm)
		return det > 0 ? 1 : -1;
	{
		qaws_wide LA = tt_wide_len2(A), LB = tt_wide_len2(B), LC = tt_wide_len2(C), LD = tt_wide_len2(D);
		qaws_wide W1 = tt_wide_det3(B, C, D), W2 = tt_wide_det3(A, C, D), W3 = tt_wide_det3(A, B, D), W4 = tt_wide_det3(A, B, C);
		qaws_wide t1 = qaws_wide_mul(&LA, &W1), t2 = qaws_wide_mul(&LB, &W2), t3 = qaws_wide_mul(&LC, &W3), t4 = qaws_wide_mul(&LD, &W4);
		t1 = qaws_wide_sub(t1, &t2);
		t1 = qaws_wide_add(t1, &t3);
		t1 = qaws_wide_sub(t1, &t4);
		return qaws_wide_sign(&t1);
	}
}

/* ======================================================================== */
/*  The triangulation                                                       */
/* ======================================================================== */

/* directed edges to values, open addressing */
typedef struct tt_emap
{
	unsigned int *ka, *kb, *val;
	unsigned int mask, cap;
} tt_emap;

static int tt_emap_reset(tt_emap* e, unsigned int n)
{
	unsigned int cap = qaws_internal_grow(16, qaws_internal_mul2(2, (size_t)n + 1), 16);
	if (!cap || (cap & (cap - 1)))
		return 0;   /* (a power of two: the probe masks with cap - 1) */
	if (cap > e->cap)
	{
		unsigned int *a = (unsigned int*)qaws_internal_realloc(e->ka, qaws_internal_mul2(sizeof(unsigned int), cap)), *b, *v;
		if (!a)
			return 0;
		e->ka = a;
		b = (unsigned int*)qaws_internal_realloc(e->kb, qaws_internal_mul2(sizeof(unsigned int), cap));
		if (!b)
			return 0;
		e->kb = b;
		v = (unsigned int*)qaws_internal_realloc(e->val, qaws_internal_mul2(sizeof(unsigned int), cap));
		if (!v)
			return 0;
		e->val = v;
		e->cap = cap;
	}
	e->mask = cap - 1;
	memset(e->ka, 0xFF, sizeof(unsigned int) * cap);
	return 1;
}

static void tt_emap_put(tt_emap* e, unsigned int a, unsigned int b, unsigned int v)
{
	unsigned int h = (a * 2654435761u ^ b * 40503u) & e->mask;
	while (e->ka[h] != HM_NONE)
		h = (h + 1) & e->mask;
	e->ka[h] = a;
	e->kb[h] = b;
	e->val[h] = v;
}

static unsigned int tt_emap_get(tt_emap const* e, unsigned int a, unsigned int b)
{
	unsigned int h = (a * 2654435761u ^ b * 40503u) & e->mask;
	while (e->ka[h] != HM_NONE)
	{
		if (e->ka[h] == a && e->kb[h] == b)
			return e->val[h];
		h = (h + 1) & e->mask;
	}
	return HM_NONE;
}

static void tt_emap_free(tt_emap* e)
{
	qaws_internal_free(e->ka);
	qaws_internal_free(e->kb);
	qaws_internal_free(e->val);
	memset(e, 0, sizeof *e);
}

typedef struct tt_mesh
{
	int64_t* q;                         /* lattice points, 3 per vertex */
	double* x;                          /* positions */
	unsigned int *vst, *v2t;            /* per vertex: a stamp; a live cell */
	unsigned int nv, capv;
	unsigned int *tv, *tn;              /* 4 per cell: vertices; the cell across facet i, HM_NONE past the super cell */
	unsigned char* tf;                  /* per cell: bit i facet i a wall; TT_DEAD; TT_IN */
	unsigned int* stamp;                /* per cell */
	unsigned int nt, capt, gen;
	unsigned int *freel, nfree;
	unsigned int last;
	uint32_t rng;
	unsigned int *cav, ncav, capcav;    /* a cavity's or a star's cells */
	unsigned int *bc, *bf, nb, capb, capbf;  /* a cavity's border facets: cell, facet */
	unsigned int *made, nmade, capmade; /* the cells the last insertion made */
	tt_emap link;
} tt_mesh;

#define TT_Q(m, v) ((m)->q + 3 * (size_t)(v))

static int tt_grow(void** p, unsigned int* cap, size_t need, size_t size)
{
	void* g;
	unsigned int c;
	if (need <= *cap)
		return 1;
	c = qaws_internal_grow(*cap, need, 64);
	if (!c)
		return 0;
	g = qaws_internal_realloc(*p, qaws_internal_mul2(size, c));
	if (!g)
		return 0;
	*p = g;
	*cap = c;
	return 1;
}

static int tt_push(unsigned int** a, unsigned int* n, unsigned int* cap, unsigned int v)
{
	if (!tt_grow((void**)a, cap, (size_t)*n + 1, sizeof(unsigned int)))
		return 0;
	(*a)[(*n)++] = v;
	return 1;
}

static void tt_free(tt_mesh* m)
{
	qaws_internal_free(m->q); qaws_internal_free(m->x); qaws_internal_free(m->vst); qaws_internal_free(m->v2t);
	qaws_internal_free(m->tv); qaws_internal_free(m->tn); qaws_internal_free(m->tf); qaws_internal_free(m->stamp); qaws_internal_free(m->freel);
	qaws_internal_free(m->cav); qaws_internal_free(m->bc); qaws_internal_free(m->bf); qaws_internal_free(m->made);
	tt_emap_free(&m->link);
	memset(m, 0, sizeof *m);
}

static unsigned int tt_add_vertex(tt_mesh* m, int64_t const* q, double const* x)
{
	unsigned int c = m->capv;
	if (m->nv + 1 > m->capv)
	{
		if (!tt_grow((void**)&m->q, &c, (size_t)m->nv + 1, sizeof(int64_t) * 3))
			return HM_NONE;
		c = m->capv;
		if (!tt_grow((void**)&m->x, &c, (size_t)m->nv + 1, sizeof(double) * 3))
			return HM_NONE;
		c = m->capv;
		if (!tt_grow((void**)&m->vst, &c, (size_t)m->nv + 1, sizeof(unsigned int)))
			return HM_NONE;
		c = m->capv;
		if (!tt_grow((void**)&m->v2t, &c, (size_t)m->nv + 1, sizeof(unsigned int)))
			return HM_NONE;
		m->capv = c;
	}
	memcpy(m->q + 3 * (size_t)m->nv, q, sizeof(int64_t) * 3);
	memcpy(m->x + 3 * (size_t)m->nv, x, sizeof(double) * 3);
	m->vst[m->nv] = 0;
	m->v2t[m->nv] = HM_NONE;
	return m->nv++;
}

static unsigned int tt_new_cell(tt_mesh* m)
{
	unsigned int t;
	if (m->nfree)
		t = m->freel[--m->nfree];
	else
	{
		unsigned int c = m->capt;
		if (m->nt + 1 > m->capt)
		{
			if (!tt_grow((void**)&m->tv, &c, (size_t)m->nt + 1, sizeof(unsigned int) * 4))
				return HM_NONE;
			c = m->capt;
			if (!tt_grow((void**)&m->tn, &c, (size_t)m->nt + 1, sizeof(unsigned int) * 4))
				return HM_NONE;
			c = m->capt;
			if (!tt_grow((void**)&m->tf, &c, (size_t)m->nt + 1, 1))
				return HM_NONE;
			c = m->capt;
			if (!tt_grow((void**)&m->stamp, &c, (size_t)m->nt + 1, sizeof(unsigned int)))
				return HM_NONE;
			c = m->capt;
			if (!tt_grow((void**)&m->freel, &c, (size_t)m->nt + 1, sizeof(unsigned int)))
				return HM_NONE;
			m->capt = c;
		}
		t = m->nt++;
	}
	m->tf[t] = 0;
	m->stamp[t] = 0;
	return t;
}

/* a live cell holding p (a walk, the facet order at random), or HM_NONE outside the super cell */
static unsigned int tt_locate(tt_mesh* m, int64_t const* p)
{
	unsigned int t = m->last, guard;
	if (t >= m->nt || (m->tf[t] & TT_DEAD))
		for (t = 0; t < m->nt && (m->tf[t] & TT_DEAD); t++) {}
	for (guard = 0; guard < 4 * m->nt + 64; guard++)
	{
		unsigned int k, r, moved = 0;
		m->rng = m->rng * 1664525u + 1013904223u;
		r = m->rng >> 30;
		for (k = 0; k < 4 && !moved; k++)
		{
			unsigned int f = (r + k) & 3;
			unsigned int const* v = m->tv + 4 * (size_t)t;
			if (tt_orient(TT_Q(m, v[tt_face[f][0]]), TT_Q(m, v[tt_face[f][1]]), TT_Q(m, v[tt_face[f][2]]), p) > 0)
			{
				if (m->tn[4 * (size_t)t + f] == HM_NONE)
					return HM_NONE;
				t = m->tn[4 * (size_t)t + f];
				moved = 1;
			}
		}
		if (!moved)
			return t;
	}
	return HM_NONE;
}

/* p inserted (Bowyer-Watson); constrained: the cavity kept off the walls, the point inside. The new vertex, or
   HM_NONE with *why: 0 memory, 1 a vertex already there, 2 refused (a wall faces away, or a vertex would be swallowed),
   3 outside */
static unsigned int tt_insert(tt_mesh* m, int64_t const* p, double const* x, int constrained, unsigned int* why)
{
	unsigned int t, i, k, pv, vg;
	t = tt_locate(m, p);
	*why = 3;
	if (t == HM_NONE || (constrained && !(m->tf[t] & TT_IN)))
		return HM_NONE;
	*why = 1;
	for (k = 0; k < 4; k++)
		if (!memcmp(TT_Q(m, m->tv[4 * (size_t)t + k]), p, sizeof(int64_t) * 3))
			return HM_NONE;
	*why = 0;
	/* the cavity: the cells whose circumsphere holds p, reached across facets (walls not crossed when constrained) */
	m->gen++;
	m->ncav = 0;
	if (!tt_push(&m->cav, &m->ncav, &m->capcav, t))
		return HM_NONE;
	m->stamp[t] = m->gen;
	for (i = 0; i < m->ncav; i++)
	{
		unsigned int c = m->cav[i], f;
		for (f = 0; f < 4; f++)
		{
			unsigned int u = m->tn[4 * (size_t)c + f];
			unsigned int const* v;
			if (u == HM_NONE || m->stamp[u] == m->gen || (constrained && (m->tf[c] & (1u << f))))
				continue;
			v = m->tv + 4 * (size_t)u;
			if (tt_insphere(TT_Q(m, v[0]), TT_Q(m, v[1]), TT_Q(m, v[2]), TT_Q(m, v[3]), p) > 0)
			{
				m->stamp[u] = m->gen;
				if (!tt_push(&m->cav, &m->ncav, &m->capcav, u))
					return HM_NONE;
			}
		}
	}
	/* grown until every facet of its border faces p */
	for (;;)
	{
		unsigned int bad = HM_NONE, badf = 0, u;
		m->nb = 0;
		for (i = 0; i < m->ncav && bad == HM_NONE; i++)
		{
			unsigned int c = m->cav[i], f;
			for (f = 0; f < 4; f++)
			{
				unsigned int const* v = m->tv + 4 * (size_t)c;
				u = m->tn[4 * (size_t)c + f];
				if (u != HM_NONE && m->stamp[u] == m->gen)
					continue;
				if (tt_orient(TT_Q(m, v[tt_face[f][0]]), TT_Q(m, v[tt_face[f][1]]), TT_Q(m, v[tt_face[f][2]]), p) >= 0)
				{
					bad = c;
					badf = f;
					break;
				}
				if (!tt_push(&m->bc, &m->nb, &m->capb, c))
					return HM_NONE;
				m->nb--;
				if (!tt_push(&m->bf, &m->nb, &m->capbf, f))
					return HM_NONE;
			}
		}
		if (bad == HM_NONE)
			break;
		u = m->tn[4 * (size_t)bad + badf];
		if (u == HM_NONE || (constrained && (m->tf[bad] & (1u << badf))))
		{
			*why = 2;
			return HM_NONE;
		}
		m->stamp[u] = m->gen;
		if (!tt_push(&m->cav, &m->ncav, &m->capcav, u))
			return HM_NONE;
	}
	/* no vertex swallowed: every cavity cell's vertices on its border */
	vg = ++m->gen;
	for (i = 0; i < m->nb; i++)
		for (k = 0; k < 3; k++)
			m->vst[m->tv[4 * (size_t)m->bc[i] + tt_face[m->bf[i]][k]]] = vg;
	for (i = 0; i < m->ncav; i++)
		for (k = 0; k < 4; k++)
			if (m->vst[m->tv[4 * (size_t)m->cav[i] + k]] != vg)
			{
				*why = 2;
				return HM_NONE;
			}
	/* the cells of the border's facets and p */
	pv = tt_add_vertex(m, p, x);
	if (pv == HM_NONE || !tt_emap_reset(&m->link, 3 * m->nb))
		return HM_NONE;
	m->nmade = 0;
	for (i = 0; i < m->nb; i++)
	{
		unsigned int c = m->bc[i], f = m->bf[i], a, b, d, u, n, g;
		a = m->tv[4 * (size_t)c + tt_face[f][0]];
		b = m->tv[4 * (size_t)c + tt_face[f][1]];
		d = m->tv[4 * (size_t)c + tt_face[f][2]];
		u = m->tn[4 * (size_t)c + f];
		n = tt_new_cell(m);
		if (n == HM_NONE || !tt_push(&m->made, &m->nmade, &m->capmade, n))
			return HM_NONE;
		m->tv[4 * (size_t)n] = a;
		m->tv[4 * (size_t)n + 1] = d;
		m->tv[4 * (size_t)n + 2] = b;
		m->tv[4 * (size_t)n + 3] = pv;
		m->tn[4 * (size_t)n] = m->tn[4 * (size_t)n + 1] = m->tn[4 * (size_t)n + 2] = HM_NONE;
		m->tn[4 * (size_t)n + 3] = u;
		m->tf[n] = (unsigned char)((((m->tf[c] >> f) & 1u) << 3) | (m->tf[c] & TT_IN));
		if (u != HM_NONE)
			for (g = 0; g < 4; g++)
				if (m->tn[4 * (size_t)u + g] == c)
				{
					m->tn[4 * (size_t)u + g] = n;
					break;
				}
		/* the base's edges: facet 0 holds b -> d, facet 1 a -> b, facet 2 d -> a */
		tt_emap_put(&m->link, b, d, 4 * n);
		tt_emap_put(&m->link, a, b, 4 * n + 1);
		tt_emap_put(&m->link, d, a, 4 * n + 2);
	}
	for (i = 0; i < m->nmade; i++)
	{
		unsigned int n = m->made[i], f;
		for (f = 0; f < 3; f++)
		{
			unsigned int const* v = m->tv + 4 * (size_t)n;
			unsigned int from = f == 0 ? v[2] : f == 1 ? v[0] : v[1], to = f == 0 ? v[1] : f == 1 ? v[2] : v[0];
			unsigned int o = tt_emap_get(&m->link, to, from);
			if (o == HM_NONE)
			{
				*why = 0;
				return HM_NONE;   /* (a cavity whose border is not a closed surface: cannot happen) */
			}
			m->tn[4 * (size_t)n + f] = o / 4;
		}
	}
	for (i = 0; i < m->ncav; i++)
	{
		m->tf[m->cav[i]] |= TT_DEAD;
		m->freel[m->nfree++] = m->cav[i];
	}
	for (i = 0; i < m->nmade; i++)
		for (k = 0; k < 4; k++)
			m->v2t[m->tv[4 * (size_t)m->made[i] + k]] = m->made[i];
	m->last = m->made[0];
	return pv;
}

/* the live cells holding vertex a, into cav */
static int tt_star(tt_mesh* m, unsigned int a)
{
	unsigned int i;
	m->gen++;
	m->ncav = 0;
	if (m->v2t[a] == HM_NONE || !tt_push(&m->cav, &m->ncav, &m->capcav, m->v2t[a]))
		return 0;
	m->stamp[m->v2t[a]] = m->gen;
	for (i = 0; i < m->ncav; i++)
	{
		unsigned int c = m->cav[i], k, ia = 4;
		for (k = 0; k < 4; k++)
			if (m->tv[4 * (size_t)c + k] == a)
				ia = k;
		for (k = 0; k < 4; k++)
		{
			unsigned int u = m->tn[4 * (size_t)c + k];
			if (k == ia || u == HM_NONE || m->stamp[u] == m->gen)
				continue;
			m->stamp[u] = m->gen;
			if (!tt_push(&m->cav, &m->ncav, &m->capcav, u))
				return 0;
		}
	}
	return 1;
}

/* a cell with face {a, b, c} and the facet it is; 0 when the face is not in the triangulation */
static int tt_find_face(tt_mesh* m, unsigned int a, unsigned int b, unsigned int c, unsigned int* oc, unsigned int* of)
{
	unsigned int i;
	if (!tt_star(m, a))
		return 0;
	for (i = 0; i < m->ncav; i++)
	{
		unsigned int t = m->cav[i], k, ia = 4, ib = 4, ic = 4;
		for (k = 0; k < 4; k++)
		{
			unsigned int v = m->tv[4 * (size_t)t + k];
			ia = v == a ? k : ia;
			ib = v == b ? k : ib;
			ic = v == c ? k : ic;
		}
		if (ib < 4 && ic < 4)
		{
			*oc = t;
			*of = 6 - ia - ib - ic;
			return 1;
		}
	}
	return 0;
}

static int tt_has_edge(tt_mesh* m, unsigned int a, unsigned int b)
{
	unsigned int i, k;
	if (!tt_star(m, a))
		return 0;
	for (i = 0; i < m->ncav; i++)
		for (k = 0; k < 4; k++)
			if (m->tv[4 * (size_t)m->cav[i] + k] == b)
				return 1;
	return 0;
}

/* ======================================================================== */
/*  Geometry in doubles                                                     */
/* ======================================================================== */

static double tt_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void tt_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

/* the squared distance from p to triangle abc (Ericson's regions) */
static double tt_tri_d2(double const* p, double const* a, double const* b, double const* c)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, v, w, q[3], denom;
	int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = b[k] - a[k];
		ac[k] = c[k] - a[k];
		ap[k] = p[k] - a[k];
		bp[k] = p[k] - b[k];
		cp[k] = p[k] - c[k];
	}
	d1 = tt_dot(ab, ap);
	d2 = tt_dot(ac, ap);
	d3 = tt_dot(ab, bp);
	d4 = tt_dot(ac, bp);
	d5 = tt_dot(ab, cp);
	d6 = tt_dot(ac, cp);
	vc = d1 * d4 - d3 * d2;
	vb = d5 * d2 - d1 * d6;
	va = d3 * d6 - d5 * d4;
	if (d1 <= 0 && d2 <= 0)
		v = 0, w = 0;
	else if (d3 >= 0 && d4 <= d3)
		v = 1, w = 0;
	else if (d6 >= 0 && d5 <= d6)
		v = 0, w = 1;
	else if (vc <= 0 && d1 >= 0 && d3 <= 0)
		v = d1 / (d1 - d3), w = 0;
	else if (vb <= 0 && d2 >= 0 && d6 <= 0)
		v = 0, w = d2 / (d2 - d6);
	else if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
	{
		w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		v = 1 - w;
	}
	else
	{
		denom = 1 / (va + vb + vc);
		v = vb * denom;
		w = vc * denom;
	}
	for (k = 0; k < 3; k++)
		q[k] = a[k] + ab[k] * v + ac[k] * w - p[k];
	return tt_dot(q, q);
}

/* the border's triangles near a point */
typedef struct tt_near
{
	double const* x;
	unsigned int const* bt;
	double const* p;
	double best;
} tt_near;

static void tt_near_visit(void* user, unsigned int item)
{
	tt_near* w = (tt_near*)user;
	unsigned int const* t = w->bt + 3 * (size_t)item;
	double d = tt_tri_d2(w->p, w->x + 3 * (size_t)t[0], w->x + 3 * (size_t)t[1], w->x + 3 * (size_t)t[2]);
	if (d < w->best)
		w->best = d;
}

/* is p within r of the border */
static int tt_near_border(qaws_bp_grid const* g, double const* x, unsigned int const* bt, double const* p, double r)
{
	tt_near w;
	unsigned int ring;
	w.x = x;
	w.bt = bt;
	w.p = p;
	w.best = HUGE_VAL;
	for (ring = 0;; ring++)
	{
		double lb = qaws_internal_grid_ring(g, p, ring, tt_near_visit, &w);
		if (w.best < r * r)
			return 1;
		if (lb == HUGE_VAL || lb > r)
			return 0;
	}
}

/* a cell's quality: the smallest sine of its dihedral angles (sin = 3 V |e| / (2 A A'), negative when inverted) */
static double tt_quality(double const* a, double const* b, double const* c, double const* d)
{
	double const* p[4];
	double area[4], e[3], vol, u[3], w[3], cr[3], q = HUGE_VAL;
	unsigned int i, j, k;
	p[0] = a; p[1] = b; p[2] = c; p[3] = d;
	for (k = 0; k < 3; k++)
	{
		u[k] = b[k] - a[k];
		w[k] = c[k] - a[k];
		e[k] = d[k] - a[k];
	}
	tt_cross(u, w, cr);
	vol = tt_dot(cr, e) / 6;
	for (i = 0; i < 4; i++)
	{
		double const *x0 = p[tt_face[i][0]], *x1 = p[tt_face[i][1]], *x2 = p[tt_face[i][2]];
		for (k = 0; k < 3; k++)
		{
			u[k] = x1[k] - x0[k];
			w[k] = x2[k] - x0[k];
		}
		tt_cross(u, w, cr);
		area[i] = 0.5 * sqrt(tt_dot(cr, cr));
	}
	/* the edge between vertices i and j joins the facets opposite the other two */
	for (i = 0; i < 4; i++)
		for (j = i + 1; j < 4; j++)
		{
			unsigned int f0 = 4, f1 = 4;
			double l;
			for (k = 0; k < 4; k++)
				if (k != i && k != j)
				{
					if (f0 == 4)
						f0 = k;
					else
						f1 = k;
				}
			for (k = 0; k < 3; k++)
				u[k] = p[j][k] - p[i][k];
			l = sqrt(tt_dot(u, u));
			q = fmin(q, area[f0] > 0 && area[f1] > 0 ? 1.5 * vol * l / (area[f0] * area[f1]) : -1);
		}
	return q;
}

/* The inside vertices moved, each to the best of a few positions (its neighbours' centre, halfway there) by the worst
   quality round it, never turning a cell over: Freitag and Ollivier-Gooch's smoothing, a few sweeps. */
static int tt_smooth(double* x, unsigned int nv, unsigned int const* cells, unsigned int nc, unsigned char const* fixed, unsigned int sweeps)
{
	unsigned int *v0 = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nv + 1, sizeof(unsigned int)), *vt = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)nc + 1)));
	unsigned int s, v, i, k;
	if (!v0 || !vt)
	{
		qaws_internal_free(v0);
		qaws_internal_free(vt);
		return 0;
	}
	for (i = 0; i < 4 * nc; i++)
		v0[cells[i] + 1]++;
	for (v = 0; v < nv; v++)
		v0[v + 1] += v0[v];
	{
		unsigned int* fill = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nv + 1)));
		if (!fill)
		{
			qaws_internal_free(v0);
			qaws_internal_free(vt);
			return 0;
		}
		memcpy(fill, v0, sizeof(unsigned int) * nv);
		for (i = 0; i < 4 * nc; i++)
			vt[fill[cells[i]]++] = i / 4;
		qaws_internal_free(fill);
	}
	for (s = 0; s < sweeps; s++)
		for (v = 0; v < nv; v++)
		{
			double keep[3], cen[3] = { 0, 0, 0 }, cand[3][3], q0 = HUGE_VAL, best;
			unsigned int nn = 0, bi = 3, c;
			if (fixed[v] || v0[v] == v0[v + 1])
				continue;
			memcpy(keep, x + 3 * (size_t)v, sizeof keep);
			for (i = v0[v]; i < v0[v + 1]; i++)
			{
				unsigned int const* t = cells + 4 * (size_t)vt[i];
				q0 = fmin(q0, tt_quality(x + 3 * (size_t)t[0], x + 3 * (size_t)t[1], x + 3 * (size_t)t[2], x + 3 * (size_t)t[3]));
				for (k = 0; k < 4; k++)
					if (t[k] != v)
					{
						cen[0] += x[3 * (size_t)t[k]];
						cen[1] += x[3 * (size_t)t[k] + 1];
						cen[2] += x[3 * (size_t)t[k] + 2];
						nn++;
					}
			}
			if (q0 > 0.5)
				continue;   /* (round it already good: sin 30 degrees) */
			for (k = 0; k < 3; k++)
			{
				cand[0][k] = cen[k] / nn;
				cand[1][k] = 0.5 * (keep[k] + cand[0][k]);
				cand[2][k] = 0.75 * keep[k] + 0.25 * cand[0][k];
			}
			best = q0;
			for (c = 0; c < 3; c++)
			{
				double q = HUGE_VAL;
				memcpy(x + 3 * (size_t)v, cand[c], sizeof keep);
				for (i = v0[v]; i < v0[v + 1] && q > best; i++)
				{
					unsigned int const* t = cells + 4 * (size_t)vt[i];
					q = fmin(q, tt_quality(x + 3 * (size_t)t[0], x + 3 * (size_t)t[1], x + 3 * (size_t)t[2], x + 3 * (size_t)t[3]));
				}
				if (q > best)
				{
					best = q;
					bi = c;
				}
			}
			memcpy(x + 3 * (size_t)v, bi < 3 ? cand[bi] : keep, sizeof keep);
		}
	qaws_internal_free(v0);
	qaws_internal_free(vt);
	return 1;
}

/* ======================================================================== */
/*  Flips                                                                   */
/* ======================================================================== */

static double tt_cell_quality(tt_mesh const* m, unsigned int const* v)
{
	return tt_quality(m->x + 3 * (size_t)v[0], m->x + 3 * (size_t)v[1], m->x + 3 * (size_t)v[2], m->x + 3 * (size_t)v[3]);
}

static void tt_sort3(unsigned int* s)
{
	unsigned int t;
	if (s[0] > s[1]) { t = s[0]; s[0] = s[1]; s[1] = t; }
	if (s[1] > s[2]) { t = s[1]; s[1] = s[2]; s[2] = t; }
	if (s[0] > s[1]) { t = s[0]; s[0] = s[1]; s[1] = t; }
}

static void tt_facet_key(unsigned int const* v, unsigned int f, unsigned int* s)
{
	s[0] = v[tt_face[f][0]];
	s[1] = v[tt_face[f][1]];
	s[2] = v[tt_face[f][2]];
	tt_sort3(s);
}

/* the cells old (a group inside, at most 4) replaced by the cells nv (at most 4, each 4 vertices, positive): the
   group's outer facets matched to the new cells' (walls kept); 0 when they do not match */
static int tt_replace(tt_mesh* m, unsigned int const* old, unsigned int nold, unsigned int const (*nv)[4], unsigned int nnew)
{
	unsigned int ext[16][6], next = 0, made[4], i, f, k, j, g, h;
	unsigned char inside = (unsigned char)(m->tf[old[0]] & TT_IN);
	for (i = 0; i < nold; i++)
		for (f = 0; f < 4; f++)
		{
			unsigned int u = m->tn[4 * (size_t)old[i] + f];
			int inner = 0;
			for (k = 0; k < nold; k++)
				inner |= u == old[k];
			if (inner)
				continue;
			tt_facet_key(m->tv + 4 * (size_t)old[i], f, ext[next]);
			ext[next][3] = old[i];
			ext[next][4] = f;
			ext[next][5] = u;
			next++;
		}
	for (i = 0; i < nnew; i++)
	{
		made[i] = tt_new_cell(m);
		if (made[i] == HM_NONE)
			return 0;
		memcpy(m->tv + 4 * (size_t)made[i], nv[i], sizeof(unsigned int) * 4);
		m->tf[made[i]] = inside;
	}
	for (i = 0; i < nnew; i++)
		for (g = 0; g < 4; g++)
		{
			unsigned int s[3], s2[3];
			int linked = 0;
			tt_facet_key(nv[i], g, s);
			for (j = 0; j < nnew && !linked; j++)
				for (h = 0; h < 4 && j != i && !linked; h++)
				{
					tt_facet_key(nv[j], h, s2);
					if (s[0] == s2[0] && s[1] == s2[1] && s[2] == s2[2])
					{
						m->tn[4 * (size_t)made[i] + g] = made[j];
						linked = 1;
					}
				}
			for (k = 0; k < next && !linked; k++)
				if (s[0] == ext[k][0] && s[1] == ext[k][1] && s[2] == ext[k][2])
				{
					unsigned int u = ext[k][5];
					m->tn[4 * (size_t)made[i] + g] = u;
					if (m->tf[ext[k][3]] & (1u << ext[k][4]))
						m->tf[made[i]] |= (unsigned char)(1u << g);
					if (u != HM_NONE)
						for (h = 0; h < 4; h++)
							if (m->tn[4 * (size_t)u + h] == ext[k][3])
							{
								m->tn[4 * (size_t)u + h] = made[i];
								break;
							}
					linked = 1;
				}
			if (!linked)
				return 0;
		}
	for (i = 0; i < nold; i++)
	{
		m->tf[old[i]] |= TT_DEAD;
		m->freel[m->nfree++] = old[i];
	}
	for (i = 0; i < nnew; i++)
		for (k = 0; k < 4; k++)
			m->v2t[nv[i][k]] = made[i];
	m->last = made[0];
	return 1;
}

/* the cells round edge {d, e} from cell t, in turn (at most cap); 0 when it meets a wall, the outside or more cells */
static unsigned int tt_edge_ring(tt_mesh const* m, unsigned int t, unsigned int d, unsigned int e, unsigned int* ring, unsigned int cap)
{
	unsigned int n = 0, cur = t, prev = HM_NONE;
	for (;;)
	{
		unsigned int f, nx = HM_NONE;
		if (n == cap)
			return 0;
		ring[n++] = cur;
		for (f = 0; f < 4; f++)
		{
			unsigned int v = m->tv[4 * (size_t)cur + f], u;
			if (v == d || v == e)
				continue;
			u = m->tn[4 * (size_t)cur + f];
			if ((m->tf[cur] & (1u << f)) || u == HM_NONE || !(m->tf[u] & TT_IN))
				return 0;
			if (u != prev && nx == HM_NONE)
				nx = u;
		}
		prev = cur;
		cur = nx;
		if (cur == t)
			return n;
	}
}

/* 3-2 and 2-3 flips round the inside cells worse than `below`, each kept when it raises the worst quality */
static int tt_flip_pass(tt_mesh* m, double below)
{
	unsigned int c, nt = m->nt;
	for (c = 0; c < nt; c++)
	{
		unsigned int i, j, done = 0, v[4];
		double q0;
		if ((m->tf[c] & TT_DEAD) || !(m->tf[c] & TT_IN))
			continue;
		memcpy(v, m->tv + 4 * (size_t)c, sizeof v);
		q0 = tt_cell_quality(m, v);
		if (q0 >= below)
			continue;
		/* 3-2: an edge with three cells round it becomes their ring's triangle */
		for (i = 0; i < 4 && !done; i++)
			for (j = i + 1; j < 4 && !done; j++)
			{
				unsigned int d = v[i], e = v[j], ring[3], r[3], nr = 0, k, h, nvv[2][4];
				double qo = HUGE_VAL, qn;
				if (tt_edge_ring(m, c, d, e, ring, 3) != 3)
					continue;
				for (k = 0; k < 3; k++)
				{
					qo = fmin(qo, tt_cell_quality(m, m->tv + 4 * (size_t)ring[k]));
					for (h = 0; h < 4; h++)
					{
						unsigned int w = m->tv[4 * (size_t)ring[k] + h], seen = w == d || w == e, z;
						for (z = 0; z < nr; z++)
							seen |= r[z] == w;
						if (!seen && nr < 3)
							r[nr++] = w;
					}
				}
				if (nr != 3)
					continue;
				nvv[0][0] = r[0]; nvv[0][1] = r[1]; nvv[0][2] = r[2]; nvv[0][3] = e;
				nvv[1][0] = r[0]; nvv[1][1] = r[2]; nvv[1][2] = r[1]; nvv[1][3] = d;
				if (tt_cell_quality(m, nvv[0]) < 0)
				{
					nvv[0][1] = r[2]; nvv[0][2] = r[1];
					nvv[1][1] = r[1]; nvv[1][2] = r[2];
				}
				qn = fmin(tt_cell_quality(m, nvv[0]), tt_cell_quality(m, nvv[1]));
				if (qn > qo + 1e-9)
				{
					if (!tt_replace(m, ring, 3, (unsigned int const(*)[4])nvv, 2))
						return 0;
					done = 1;
				}
			}
		/* 2-3: a facet shared with an inside cell becomes three cells round the edge between their apexes */
		for (i = 0; i < 4 && !done; i++)
		{
			unsigned int u = m->tn[4 * (size_t)c + i], a, b, cc, d = v[i], e = HM_NONE, k, old[2], nvv[3][4];
			double qn;
			if ((m->tf[c] & (1u << i)) || u == HM_NONE || !(m->tf[u] & TT_IN))
				continue;
			a = v[tt_face[i][0]];
			b = v[tt_face[i][1]];
			cc = v[tt_face[i][2]];
			for (k = 0; k < 4; k++)
			{
				unsigned int w = m->tv[4 * (size_t)u + k];
				if (w != a && w != b && w != cc)
					e = w;
			}
			nvv[0][0] = a; nvv[0][1] = b; nvv[0][2] = d; nvv[0][3] = e;
			nvv[1][0] = b; nvv[1][1] = cc; nvv[1][2] = d; nvv[1][3] = e;
			nvv[2][0] = cc; nvv[2][1] = a; nvv[2][2] = d; nvv[2][3] = e;
			qn = fmin(tt_cell_quality(m, nvv[0]), fmin(tt_cell_quality(m, nvv[1]), tt_cell_quality(m, nvv[2])));
			if (qn > fmin(q0, tt_cell_quality(m, m->tv + 4 * (size_t)u)) + 1e-9)
			{
				old[0] = c;
				old[1] = u;
				if (!tt_replace(m, old, 2, (unsigned int const(*)[4])nvv, 3))
					return 0;
				done = 1;
			}
		}
	}
	return 1;
}

/* ======================================================================== */
/*  The operator                                                            */
/* ======================================================================== */

typedef struct tt_ctx
{
	tt_mesh m;
	double lo[3], step;
	unsigned int *bt, nbt, capbt;       /* the border's triangles (outward) */
	unsigned char* bdead;
	unsigned int capbd;
	tt_emap bmap;                       /* its directed edges to triangles */
} tt_ctx;

static void tt_lattice(tt_ctx const* t, double const* x, int64_t* q)
{
	int k;
	for (k = 0; k < 3; k++)
		q[k] = (int64_t)floor((x[k] - t->lo[k]) / t->step + 0.5);
}

static int tt_add_tri(tt_ctx* t, unsigned int a, unsigned int b, unsigned int c)
{
	unsigned int cap = t->capbt;
	if (!tt_grow((void**)&t->bt, &cap, (size_t)t->nbt + 1, sizeof(unsigned int) * 3))
		return 0;
	t->capbt = cap;
	cap = t->capbd;
	if (!tt_grow((void**)&t->bdead, &cap, (size_t)t->nbt + 1, 1))
		return 0;
	t->capbd = cap;
	t->bt[3 * (size_t)t->nbt] = a;
	t->bt[3 * (size_t)t->nbt + 1] = b;
	t->bt[3 * (size_t)t->nbt + 2] = c;
	t->bdead[t->nbt++] = 0;
	return 1;
}

/* a point added to the border at x: into the triangulation; its vertex or HM_NONE */
static unsigned int tt_border_point(tt_ctx* t, double const* x)
{
	int64_t q[3];
	unsigned int why;
	tt_lattice(t, x, q);
	return tt_insert(&t->m, q, x, 0, &why);
}

/* integer 2D orientation in the plane that drops axis `drop` */
static int tt_orient2(int64_t const* a, int64_t const* b, int64_t const* c, int drop)
{
	int i = drop == 0 ? 1 : 0, j = drop == 2 ? 1 : 2;
	qaws_wide w = qaws_wide_cross(b[i] - a[i], b[j] - a[j], c[i] - a[i], c[j] - a[j]);
	return qaws_wide_sign(&w);
}

/* the border's directed edges to its live triangles; 0 memory */
static int tt_border_map(tt_ctx* t)
{
	unsigned int i;
	if (!tt_emap_reset(&t->bmap, 3 * t->nbt))
		return 0;
	for (i = 0; i < t->nbt; i++)
		if (!t->bdead[i])
		{
			tt_emap_put(&t->bmap, t->bt[3 * i], t->bt[3 * i + 1], i);
			tt_emap_put(&t->bmap, t->bt[3 * i + 1], t->bt[3 * i + 2], i);
			tt_emap_put(&t->bmap, t->bt[3 * i + 2], t->bt[3 * i], i);
		}
	return 1;
}

/* the vertex of border triangle j that is neither x nor y */
static unsigned int tt_third(tt_ctx const* t, unsigned int j, unsigned int x, unsigned int y)
{
	unsigned int const* T = t->bt + 3 * (size_t)j;
	return T[0] != x && T[0] != y ? T[0] : T[1] != x && T[1] != y ? T[1] : T[2];
}

/* x y z and y x w coplanar, their union a convex quad (so the diagonal x y may become z w) */
static int tt_flippable(tt_mesh const* m, unsigned int x, unsigned int y, unsigned int z, unsigned int w)
{
	double nrm[3], u1[3], u2[3];
	unsigned int k, drop;
	if (tt_orient(TT_Q(m, x), TT_Q(m, y), TT_Q(m, z), TT_Q(m, w)) != 0)
		return 0;
	for (k = 0; k < 3; k++)
	{
		u1[k] = (double)(TT_Q(m, y)[k] - TT_Q(m, x)[k]);
		u2[k] = (double)(TT_Q(m, z)[k] - TT_Q(m, x)[k]);
	}
	tt_cross(u1, u2, nrm);
	drop = fabs(nrm[0]) >= fabs(nrm[1]) && fabs(nrm[0]) >= fabs(nrm[2]) ? 0 : fabs(nrm[1]) >= fabs(nrm[2]) ? 1 : 2;
	return tt_orient2(TT_Q(m, z), TT_Q(m, w), TT_Q(m, x), (int)drop) * tt_orient2(TT_Q(m, z), TT_Q(m, w), TT_Q(m, y), (int)drop) < 0;
}

/* is w (in the plane of x y z) strictly inside their circumcircle: the sphere through them and a point off the plane */
static int tt_incircle(tt_mesh const* m, unsigned int x, unsigned int y, unsigned int z, unsigned int w)
{
	int64_t const *a = TT_Q(m, x), *b = TT_Q(m, y), *c = TT_Q(m, z);
	int64_t u[3], v[3], n[3], e[3], big = 0;
	int k, s = 0;
	for (k = 0; k < 3; k++)
	{
		u[k] = b[k] - a[k];
		v[k] = c[k] - a[k];
	}
	n[0] = u[1] * v[2] - u[2] * v[1];
	n[1] = u[2] * v[0] - u[0] * v[2];
	n[2] = u[0] * v[1] - u[1] * v[0];
	for (k = 0; k < 3; k++)
		big = n[k] > big ? n[k] : -n[k] > big ? -n[k] : big;
	while ((big >> s) > ((int64_t)1 << 30))
		s++;
	for (k = 0; k < 3; k++)
		e[k] = a[k] + (n[k] >= 0 ? n[k] >> s : -((-n[k]) >> s));   /* (any point on the normal's side will do) */
	return tt_insphere(a, b, c, e, TT_Q(m, w)) > 0;
}

/* the border's coplanar neighbours flipped until every pair is Delaunay in its plane (a convex polyhedron's border
   then lies in its Delaunay triangulation, but for ties); 0 memory */
static int tt_lawson(tt_ctx* t, unsigned int* flips)
{
	tt_mesh* m = &t->m;
	unsigned char* mark = NULL;
	unsigned int capmark = 0, pass, changed = 1;
	for (pass = 0; pass < 100 && changed; pass++)
	{
		unsigned int i, n0 = t->nbt;
		changed = 0;
		if (!tt_border_map(t) || !tt_grow((void**)&mark, &capmark, (size_t)t->nbt + 1, 1))
		{
			qaws_internal_free(mark);
			return 0;
		}
		memset(mark, 0, t->nbt);
		for (i = 0; i < n0; i++)
		{
			unsigned int e;
			if (t->bdead[i] || mark[i])
				continue;
			for (e = 0; e < 3; e++)
			{
				unsigned int x = t->bt[3 * (size_t)i + e], y = t->bt[3 * (size_t)i + (e + 1) % 3], z = t->bt[3 * (size_t)i + (e + 2) % 3], j, w;
				j = tt_emap_get(&t->bmap, y, x);
				if (j == HM_NONE || mark[j])
					continue;
				w = tt_third(t, j, x, y);
				if (!tt_flippable(m, x, y, z, w) || !tt_incircle(m, x, y, z, w))
					continue;
				mark[i] = mark[j] = 1;
				t->bdead[i] = t->bdead[j] = 1;
				if (!tt_add_tri(t, z, x, w) || !tt_add_tri(t, w, y, z))
				{
					qaws_internal_free(mark);
					return 0;
				}
				(*flips)++;
				changed = 1;
				break;
			}
		}
	}
	qaws_internal_free(mark);
	return 1;
}

/* The border recovered: every triangle a face of the triangulation. 1 done, 0 memory, -1 not converged */
static int tt_recover(tt_ctx* t, unsigned int* passes, unsigned int* splits, unsigned int* flips)
{
	tt_mesh* m = &t->m;
	unsigned char* mark = NULL;
	unsigned int capmark = 0, pass, limit = 4 * m->nv + 1000;
	int result = -1;
	for (pass = 0; pass < 200; pass++)
	{
		unsigned int i, n0, missing = 0;
		if (!tt_lawson(t, flips) || !tt_border_map(t) || !tt_grow((void**)&mark, &capmark, (size_t)t->nbt + 1, 1))
		{
			result = 0;
			break;
		}
		n0 = t->nbt;
		memset(mark, 0, t->nbt);
		for (i = 0; i < n0; i++)
		{
			unsigned int const* T;
			unsigned int a, b, c, oc, of, e, done = 0;
			if (t->bdead[i] || mark[i])
				continue;
			T = t->bt + 3 * (size_t)i;
			a = T[0];
			b = T[1];
			c = T[2];
			if (tt_find_face(m, a, b, c, &oc, &of))
				continue;
			missing++;
			/* a coplanar neighbour whose shared diagonal the triangulation turned the other way (cospherical ties): flipped */
			for (e = 0; e < 3 && !done; e++)
			{
				unsigned int x = T[e], y = T[(e + 1) % 3], z = T[(e + 2) % 3], j, w, c1, f1;
				j = tt_emap_get(&t->bmap, y, x);
				if (j == HM_NONE || mark[j])
					continue;
				w = tt_third(t, j, x, y);
				if (!tt_flippable(m, x, y, z, w) || !tt_find_face(m, z, x, w, &c1, &f1) || !tt_find_face(m, w, y, z, &c1, &f1))
					continue;
				mark[i] = mark[j] = 1;
				t->bdead[i] = t->bdead[j] = 1;
				if (!tt_add_tri(t, z, x, w) || !tt_add_tri(t, w, y, z))
				{
					result = 0;
					goto out;
				}
				(*flips)++;
				done = 1;
			}
			if (done)
				continue;
			/* an edge missing: split at its midpoint (both its triangles) */
			for (e = 0; e < 3 && !done; e++)
			{
				unsigned int x = T[e], y = T[(e + 1) % 3], z = T[(e + 2) % 3], j, w, v;
				double mid[3];
				if (tt_has_edge(m, x, y))
					continue;
				j = tt_emap_get(&t->bmap, y, x);
				if (j == HM_NONE || mark[j])
				{
					done = 1;
					break;
				}
				w = tt_third(t, j, x, y);
				for (v = 0; v < 3; v++)
					mid[v] = 0.5 * (m->x[3 * (size_t)x + v] + m->x[3 * (size_t)y + v]);
				v = tt_border_point(t, mid);
				if (v == HM_NONE)
					goto out;
				mark[i] = mark[j] = 1;
				t->bdead[i] = t->bdead[j] = 1;
				if (!tt_add_tri(t, x, v, z) || !tt_add_tri(t, v, y, z) || !tt_add_tri(t, y, v, w) || !tt_add_tri(t, v, x, w))
				{
					result = 0;
					goto out;
				}
				(*splits)++;
				done = 1;
			}
			if (done)
				continue;
			/* else: split at the circumcentre, or the longest edge when the circumcentre is outside the triangle */
			{
				double const *A = m->x + 3 * (size_t)a, *B = m->x + 3 * (size_t)b, *C = m->x + 3 * (size_t)c;
				double ab[3], ac[3], n[3], t1[3], t2[3], cc[3], nn, l2b, l2c, bary[3], tot;
				unsigned int k, v;
				for (k = 0; k < 3; k++)
				{
					ab[k] = B[k] - A[k];
					ac[k] = C[k] - A[k];
				}
				tt_cross(ab, ac, n);
				nn = tt_dot(n, n);
				l2b = tt_dot(ab, ab);
				l2c = tt_dot(ac, ac);
				tt_cross(n, ab, t1);
				tt_cross(ac, n, t2);
				for (k = 0; k < 3; k++)
					cc[k] = A[k] + (l2c * t1[k] + l2b * t2[k]) / (2 * nn);
				/* its barycentrics */
				{
					double pa[3], pb[3], pc[3], cr[3];
					for (k = 0; k < 3; k++)
					{
						pa[k] = A[k] - cc[k];
						pb[k] = B[k] - cc[k];
						pc[k] = C[k] - cc[k];
					}
					tt_cross(pb, pc, cr);
					bary[0] = tt_dot(cr, n);
					tt_cross(pc, pa, cr);
					bary[1] = tt_dot(cr, n);
					tt_cross(pa, pb, cr);
					bary[2] = tt_dot(cr, n);
					tot = bary[0] + bary[1] + bary[2];
				}
				if (nn > 0 && bary[0] > 0.02 * tot && bary[1] > 0.02 * tot && bary[2] > 0.02 * tot)
				{
					v = tt_border_point(t, cc);
					if (v == HM_NONE)
						goto out;
					mark[i] = 1;
					t->bdead[i] = 1;
					if (!tt_add_tri(t, a, b, v) || !tt_add_tri(t, b, c, v) || !tt_add_tri(t, c, a, v))
					{
						result = 0;
						goto out;
					}
					(*splits)++;
				}
				else
				{
					/* the longest edge */
					unsigned int best = 0, x, y, z, j, w;
					double lmax = -1, mid[3];
					for (e = 0; e < 3; e++)
					{
						double d2 = 0;
						for (k = 0; k < 3; k++)
							d2 += (m->x[3 * (size_t)T[(e + 1) % 3] + k] - m->x[3 * (size_t)T[e] + k]) * (m->x[3 * (size_t)T[(e + 1) % 3] + k] - m->x[3 * (size_t)T[e] + k]);
						if (d2 > lmax)
						{
							lmax = d2;
							best = e;
						}
					}
					x = T[best];
					y = T[(best + 1) % 3];
					z = T[(best + 2) % 3];
					j = tt_emap_get(&t->bmap, y, x);
					if (j == HM_NONE || mark[j])
						continue;
					w = tt_third(t, j, x, y);
					for (k = 0; k < 3; k++)
						mid[k] = 0.5 * (m->x[3 * (size_t)x + k] + m->x[3 * (size_t)y + k]);
					v = tt_border_point(t, mid);
					if (v == HM_NONE)
						goto out;
					mark[i] = mark[j] = 1;
					t->bdead[i] = t->bdead[j] = 1;
					if (!tt_add_tri(t, x, v, z) || !tt_add_tri(t, v, y, z) || !tt_add_tri(t, y, v, w) || !tt_add_tri(t, v, x, w))
					{
						result = 0;
						goto out;
					}
					(*splits)++;
				}
			}
		}
		*passes = pass + 1;
		if (!missing)
		{
			result = 1;
			break;
		}
		if (*splits > limit)
			break;
	}
out:
	qaws_internal_free(mark);
	return result;
}

/* ======================================================================== */
/*  STUFFING: isosurface stuffing (Labelle and Shewchuk 2007)               */
/* ======================================================================== */

static int tt_u64_cmp(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

/* the index of edge {a, b} among the sorted edges (present by construction) */
static size_t tt_find_edge(unsigned long long const* e, size_t ne, unsigned int a, unsigned int b)
{
	unsigned long long key = a < b ? (unsigned long long)a << 32 | b : (unsigned long long)b << 32 | a;
	size_t l = 0, h = ne;
	while (l < h)
	{
		size_t m = (l + h) / 2;
		if (e[m] < key)
			l = m + 1;
		else
			h = m;
	}
	return l;
}

/* a prism (bottom v[0..2], top v[3..5], v[k + 3] over v[k]) in three tetrahedra by Dompierre et al.'s rule: each quad
   cut through its smallest vertex, so prisms sharing a quad cut it alike; out: 12 vertices */
static void tt_prism(unsigned int const* v, unsigned int* out)
{
	unsigned int w[6], best = 0, r, fl, k;
	for (r = 0; r < 3; r++)
		for (fl = 0; fl < 2; fl++)
		{
			unsigned int at = (fl ? 3 : 0) + r;
			if (v[at] < v[best])
				best = at;
		}
	r = best % 3;
	fl = best >= 3;
	for (k = 0; k < 3; k++)
	{
		w[k] = v[(fl ? 3 : 0) + (k + r) % 3];
		w[k + 3] = v[(fl ? 0 : 3) + (k + r) % 3];
	}
	if ((w[1] < w[5] ? w[1] : w[5]) < (w[2] < w[4] ? w[2] : w[4]))
	{
		unsigned int t[12] = { w[0], w[1], w[2], w[5], w[0], w[1], w[5], w[4], w[0], w[4], w[5], w[3] };
		memcpy(out, t, sizeof t);
	}
	else
	{
		unsigned int t[12] = { w[0], w[1], w[2], w[4], w[0], w[4], w[2], w[5], w[0], w[4], w[5], w[3] };
		memcpy(out, t, sizeof t);
	}
}

/* STUFFING: a body-centred cubic lattice of spacing h over the border's box (cube corners and centres; each pair of
   neighbouring centres and the four edges of their shared square make four tetrahedra), the signed distance to the
   border at its vertices (negative inside: the winding number, so soups and small holes do), each edge crossing the
   surface cut where the distance vanishes (secant steps); a vertex closer to a cut than 0.24999 of a cube edge or 0.41189
   of a diagonal (Labelle and Shewchuk's alphas) warped onto it (its cuts gone); then each tetrahedron's inside: whole,
   a corner tetrahedron, or a prism in three (tt_prism, conforming across neighbours). The border is the surface's
   approximation by the cuts and warped vertices, not the input's triangles. */
static qaws_status tt_stuffing(qaws_mesh_3d const* boundary, qaws_tetrahedralize_desc const* desc, qaws_volume_mesh** out_mesh, qaws_mesh_report* rep)
{
	qaws_hm const* hb = &boundary->hm;
	qaws_mesh_bvh* bvh = NULL;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, h, *X = NULL, *phi = NULL, *cx = NULL;
	unsigned int n[3], nc, nce, nv, i, j, k, a, *tets = NULL, ntet = 0, cap = 0, *out = NULL, nout = 0, *used = NULL, nused = 0;
	unsigned long long* edges = NULL;
	size_t ne = 0, ecap = 0, q;
	double* t = NULL;
	qaws_scalar *pts = NULL, *dist = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < hb->nv; i++)
		if (hb->vc[i] != HM_NONE)
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], hb->x[3 * (size_t)i + k]);
				hi[k] = fmax(hi[k], hb->x[3 * (size_t)i + k]);
			}
	if (!(hi[0] >= lo[0]))
		return QAWS_STATUS_INVALID_TOPOLOGY;
	h = desc->target_length > 0 ? (double)desc->target_length : 0;
	if (h == 0)
	{
		qaws_mesh_quality qq;
		qaws_mesh_3d_measure_quality(&qq, boundary);
		h = qq.mean_edge > 0 ? (double)qq.mean_edge : 0.1 * (hi[0] - lo[0] + hi[1] - lo[1] + hi[2] - lo[2]);
	}
	/* (the cubes: a cube's margin round the box) */
	for (k = 0; k < 3; k++)
	{
		n[k] = (unsigned int)ceil((hi[k] - lo[k]) / h) + 2;
		lo[k] -= h;
	}
	if ((double)(n[0] + 1) * (n[1] + 1) * (n[2] + 1) * 2 > 4e7)
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (the lattice past 40 million vertices: the size too small for the box) */
	nc = (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
	nce = n[0] * n[1] * n[2];
	nv = nc + nce;
	X = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv));
	phi = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)nv));
	pts = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv));
	dist = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)nv));
	if (!X || !phi || !pts || !dist)
		goto done;
	for (k = 0; k <= n[2]; k++)
		for (j = 0; j <= n[1]; j++)
			for (i = 0; i <= n[0]; i++)
			{
				size_t v = i + (size_t)(n[0] + 1) * (j + (size_t)(n[1] + 1) * k);
				X[3 * v] = lo[0] + h * i;
				X[3 * v + 1] = lo[1] + h * j;
				X[3 * v + 2] = lo[2] + h * k;
			}
	for (k = 0; k < n[2]; k++)
		for (j = 0; j < n[1]; j++)
			for (i = 0; i < n[0]; i++)
			{
				size_t v = nc + i + (size_t)n[0] * (j + (size_t)n[1] * k);
				X[3 * v] = lo[0] + h * (i + 0.5);
				X[3 * v + 1] = lo[1] + h * (j + 0.5);
				X[3 * v + 2] = lo[2] + h * (k + 0.5);
			}
	if ((st = qaws_mesh_bvh_create_mesh_3d(&bvh, boundary, NULL)) != QAWS_STATUS_OK)
		goto done;
	for (q = 0; q < 3 * (size_t)nv; q++)
		pts[q] = (qaws_scalar)X[q];
	if ((st = qaws_mesh_bvh_signed_distance(dist, NULL, bvh, pts, nv, NULL)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (q = 0; q < nv; q++)
		phi[q] = (double)dist[q];
	/* the lattice's tetrahedra not wholly outside, and their edges crossing the surface */
	for (a = 0; a < 3; a++)
	{
		unsigned int b = (a + 1) % 3, c = (a + 2) % 3, m[3];
		for (m[2] = 0; m[2] < n[2]; m[2]++)
			for (m[1] = 0; m[1] < n[1]; m[1]++)
				for (m[0] = 0; m[0] < n[0]; m[0]++)
				{
					unsigned int m2[3], f[4][3], e, z, c1, c2, fv[4];
					if (m[a] + 1 >= n[a])
						continue;
					m2[0] = m[0]; m2[1] = m[1]; m2[2] = m[2];
					m2[a]++;
					c1 = nc + m[0] + n[0] * (m[1] + n[1] * m[2]);
					c2 = nc + m2[0] + n[0] * (m2[1] + n[1] * m2[2]);
					/* (the shared square at m[a] + 1 along a: its corners in order round it) */
					for (z = 0; z < 4; z++)
					{
						f[z][a] = m[a] + 1;
						f[z][b] = m[b] + (z == 1 || z == 2);
						f[z][c] = m[c] + (z >= 2);
						fv[z] = f[z][0] + (n[0] + 1) * (f[z][1] + (n[1] + 1) * f[z][2]);
					}
					for (e = 0; e < 4; e++)
					{
						unsigned int tv[4] = { c1, c2, fv[e], fv[(e + 1) % 4] }, y, inside = 0;
						for (y = 0; y < 4; y++)
							inside += phi[tv[y]] <= 0;
						if (!inside)
							continue;
						if (ntet == cap)
						{
							unsigned int nc2 = 2 * cap + 1024, *g = (unsigned int*)qaws_internal_realloc(tets, qaws_internal_mul3(sizeof(unsigned int), 4, nc2));
							if (!g)
								goto done;
							tets = g;
							cap = nc2;
						}
						memcpy(tets + 4 * (size_t)ntet++, tv, sizeof tv);
						for (y = 0; y < 4; y++)
						{
							unsigned int z2;
							for (z2 = y + 1; z2 < 4; z2++)
							{
								unsigned int p = tv[y], r = tv[z2];
								if ((phi[p] < 0 && phi[r] > 0) || (phi[p] > 0 && phi[r] < 0))
								{
									if (ne == ecap)
									{
										size_t n2 = 2 * ecap + 1024;
										unsigned long long* g = (unsigned long long*)qaws_internal_realloc(edges, qaws_internal_mul2(sizeof(unsigned long long), n2));
										if (!g)
											goto done;
										edges = g;
										ecap = n2;
									}
									edges[ne++] = p < r ? (unsigned long long)p << 32 | r : (unsigned long long)r << 32 | p;
								}
							}
						}
					}
				}
	}
	qsort(edges, ne, sizeof *edges, tt_u64_cmp);
	{
		size_t u = 0;
		for (q = 0; q < ne; q++)
			if (!u || edges[u - 1] != edges[q])
				edges[u++] = edges[q];
		ne = u;
	}
	/* the cuts: secant steps on the signed distance, all edges at once */
	t = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 4, ne + 1));
	cx = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, ne + 1));
	{
		qaws_scalar* p2 = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, ne + 1));
		qaws_scalar* d2 = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), ne + 1));
		unsigned int it;
		if (!t || !cx || !p2 || !d2)
		{
			qaws_internal_free(p2);
			qaws_internal_free(d2);
			goto done;
		}
		/* (per edge: its bracket [t0, t1] and their values) */
		for (q = 0; q < ne; q++)
		{
			unsigned int p = (unsigned int)(edges[q] >> 32), r = (unsigned int)(edges[q] & 0xffffffffu);
			t[4 * q] = 0;
			t[4 * q + 1] = 1;
			t[4 * q + 2] = phi[p];
			t[4 * q + 3] = phi[r];
		}
		for (it = 0; it < 6; it++)
		{
			for (q = 0; q < ne; q++)
			{
				unsigned int p = (unsigned int)(edges[q] >> 32), r = (unsigned int)(edges[q] & 0xffffffffu);
				double s = t[4 * q] + (t[4 * q + 1] - t[4 * q]) * t[4 * q + 2] / (t[4 * q + 2] - t[4 * q + 3]);
				for (k = 0; k < 3; k++)
				{
					cx[3 * q + k] = X[3 * (size_t)p + k] + s * (X[3 * (size_t)r + k] - X[3 * (size_t)p + k]);
					p2[3 * q + k] = (qaws_scalar)cx[3 * q + k];
				}
			}
			if (it == 5)
				break;
			if ((st = qaws_mesh_bvh_signed_distance(d2, NULL, bvh, p2, (unsigned int)ne, NULL)) != QAWS_STATUS_OK)
			{
				qaws_internal_free(p2);
				qaws_internal_free(d2);
				goto done;
			}
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (q = 0; q < ne; q++)
			{
				double s = t[4 * q] + (t[4 * q + 1] - t[4 * q]) * t[4 * q + 2] / (t[4 * q + 2] - t[4 * q + 3]), d = (double)d2[q];
				if (d == 0)
				{
					t[4 * q] = t[4 * q + 1] = s;
					t[4 * q + 2] = 1;
					t[4 * q + 3] = -1;
				}
				else if ((d < 0) == (t[4 * q + 2] < 0))
				{
					t[4 * q] = s;
					t[4 * q + 2] = d;
				}
				else
				{
					t[4 * q + 1] = s;
					t[4 * q + 3] = d;
				}
			}
		}
		qaws_internal_free(p2);
		qaws_internal_free(d2);
	}
	/* warping: each vertex onto its nearest cut within alpha of its edge's length from it */
	{
		double* best = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)nv));
		unsigned int* to = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nv));
		if (!best || !to)
		{
			qaws_internal_free(best);
			qaws_internal_free(to);
			goto done;
		}
		for (q = 0; q < nv; q++)
		{
			best[q] = HUGE_VAL;
			to[q] = HM_NONE;
		}
		for (q = 0; q < ne; q++)
		{
			unsigned int p = (unsigned int)(edges[q] >> 32), r = (unsigned int)(edges[q] & 0xffffffffu), s;
			double len = 0, alpha;
			for (k = 0; k < 3; k++)
				len += (X[3 * (size_t)p + k] - X[3 * (size_t)r + k]) * (X[3 * (size_t)p + k] - X[3 * (size_t)r + k]);
			len = sqrt(len);
			alpha = len > 0.95 * h ? 0.24999 : 0.41189;
			for (s = 0; s < 2; s++)
			{
				unsigned int v = s ? r : p;
				double d = 0;
				for (k = 0; k < 3; k++)
					d += (cx[3 * q + k] - X[3 * (size_t)v + k]) * (cx[3 * q + k] - X[3 * (size_t)v + k]);
				d = sqrt(d);
				if (d < alpha * len && d < best[v])
				{
					best[v] = d;
					to[v] = (unsigned int)q;
				}
			}
		}
		for (q = 0; q < nv; q++)
			if (to[q] != HM_NONE)
			{
				memcpy(X + 3 * q, cx + 3 * (size_t)to[q], sizeof(double) * 3);
				phi[q] = 0;
			}
		qaws_internal_free(best);
		qaws_internal_free(to);
	}
	/* each tetrahedron's inside (cut ids nv + edge's index; a zero vertex is on the surface: no cut beside it) */
	{
		size_t ocap = 0;
		for (q = 0; q < ntet; q++)
		{
			unsigned int const* tv = tets + 4 * q;
			unsigned int P[4], N[4], np = 0, nn = 0, y, piece[12], npiece = 0, z;
			int sg[4];
			for (y = 0; y < 4; y++)
			{
				sg[y] = phi[tv[y]] < 0 ? -1 : phi[tv[y]] > 0 ? 1 : 0;
				if (sg[y] > 0)
					P[np++] = y;
				else
					N[nn++] = y;
			}
#define TT_CUT(x, yv) (sg[x] == 0 ? tv[x] : (unsigned int)(nv + tt_find_edge(edges, ne, tv[x], tv[yv])))
			if (np == 0)
			{
				if (sg[0] || sg[1] || sg[2] || sg[3])
				{
					memcpy(piece, tv, sizeof(unsigned int) * 4);
					npiece = 1;
				}
			}
			else if (np == 3)
			{
				unsigned int i0 = N[0];
				if (sg[i0] < 0)
				{
					piece[0] = tv[i0];
					piece[1] = TT_CUT(i0, P[0]);
					piece[2] = TT_CUT(i0, P[1]);
					piece[3] = TT_CUT(i0, P[2]);
					npiece = 1;
				}
			}
			else if (np == 1 || np == 2)
			{
				unsigned int pr[6];
				if (np == 1)
				{
					/* (bottom: the three inside, top: their cuts toward the outside one) */
					for (y = 0; y < 3; y++)
					{
						pr[y] = tv[N[y]];
						pr[y + 3] = TT_CUT(N[y], P[0]);
					}
				}
				else
				{
					/* (the wedge: each inside one with its cuts toward the two outside) */
					pr[0] = tv[N[0]];
					pr[1] = TT_CUT(N[0], P[0]);
					pr[2] = TT_CUT(N[0], P[1]);
					pr[3] = tv[N[1]];
					pr[4] = TT_CUT(N[1], P[0]);
					pr[5] = TT_CUT(N[1], P[1]);
				}
				tt_prism(pr, piece);
				npiece = 3;
			}
#undef TT_CUT
			for (z = 0; z < npiece; z++)
			{
				unsigned int* tz = piece + 4 * z;
				if (tz[0] == tz[1] || tz[0] == tz[2] || tz[0] == tz[3] || tz[1] == tz[2] || tz[1] == tz[3] || tz[2] == tz[3])
					continue;
				if (4 * (size_t)nout + 4 > ocap)
				{
					unsigned int* g;
					ocap = 2 * ocap + 4096;
					g = (unsigned int*)qaws_internal_realloc(out, qaws_internal_mul2(sizeof(unsigned int), ocap));
					if (!g)
						goto done;
					out = g;
				}
				memcpy(out + 4 * (size_t)nout++, tz, sizeof(unsigned int) * 4);
			}
		}
	}
	/* the vertices in use (lattice ones, then cuts), the cells turned positive, the flat ones dropped */
	used = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nv + ne + 1));
	if (!used)
		goto done;
	for (q = 0; q < nv + ne; q++)
		used[q] = HM_NONE;
	{
		qaws_scalar* pos = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + ne + 1));
		unsigned int* cells = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, (size_t)nout + 1));
		unsigned int kept = 0;
		qaws_volume_mesh_desc vd;
		if (!pos || !cells)
		{
			qaws_internal_free(pos);
			qaws_internal_free(cells);
			goto done;
		}
		for (q = 0; q < nout; q++)
		{
			unsigned int* tz = out + 4 * q, y;
			double const* P4[4];
			double e1[3], e2[3], e3[3], vol;
			for (y = 0; y < 4; y++)
				P4[y] = tz[y] < nv ? X + 3 * (size_t)tz[y] : cx + 3 * (size_t)(tz[y] - nv);
			for (k = 0; k < 3; k++)
			{
				e1[k] = P4[1][k] - P4[0][k];
				e2[k] = P4[2][k] - P4[0][k];
				e3[k] = P4[3][k] - P4[0][k];
			}
			vol = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
			if (fabs(vol) < 1e-9 * h * h * h)
				continue;
			if (vol < 0)
			{
				unsigned int sw = tz[2];
				tz[2] = tz[3];
				tz[3] = sw;
			}
			for (y = 0; y < 4; y++)
			{
				if (used[tz[y]] == HM_NONE)
				{
					double const* pp = tz[y] < nv ? X + 3 * (size_t)tz[y] : cx + 3 * (size_t)(tz[y] - nv);
					used[tz[y]] = nused;
					pos[3 * (size_t)nused] = (qaws_scalar)pp[0];
					pos[3 * (size_t)nused + 1] = (qaws_scalar)pp[1];
					pos[3 * (size_t)nused + 2] = (qaws_scalar)pp[2];
					nused++;
				}
				cells[4 * (size_t)kept + y] = used[tz[y]];
			}
			kept++;
		}
		memset(&vd, 0, sizeof vd);
		vd.positions = pos;
		vd.vertex_count = nused;
		vd.cell_vertices = cells;
		vd.cell_count = kept;
		st = kept ? qaws_volume_mesh_create(out_mesh, &vd) : QAWS_STATUS_INVALID_TOPOLOGY;
		memset(rep, 0, sizeof *rep);
		rep->moves = nused;
		qaws_internal_free(pos);
		qaws_internal_free(cells);
	}
done:
	qaws_mesh_bvh_destroy(bvh);
	qaws_internal_free(X); qaws_internal_free(phi); qaws_internal_free(pts); qaws_internal_free(dist); qaws_internal_free(cx); qaws_internal_free(t);
	qaws_internal_free(tets); qaws_internal_free(edges); qaws_internal_free(out); qaws_internal_free(used);
	return st;
}

qaws_status qaws_volume_mesh_tetrahedralize(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* boundary,
	qaws_tetrahedralize_desc const* desc)
{
	tt_ctx t;
	tt_mesh* m = &t.m;
	qaws_hm const* h;
	qaws_mesh_report rep;
	qaws_bp_grid* grid = NULL;
	qaws_bp_box* boxes = NULL;
	unsigned int *map = NULL, *cells = NULL, *vmap = NULL, i, k, f, c, nin = 0, ninv = 0;
	unsigned char* ctype = NULL;
	double hi[3], ext = 0, L, bound, *xo = NULL, sumlen = 0;
	unsigned int nedge = 0;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!boundary || !desc || (unsigned int)desc->method > QAWS_TETS_STUFFING || desc->target_length < 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->method == QAWS_TETS_STUFFING)
	{
		qaws_mesh_report sr;
		st = tt_stuffing(boundary, desc, out_mesh, &sr);
		if (st == QAWS_STATUS_OK && out_report)
			*out_report = sr;
		return st;
	}
	if (desc->method != QAWS_TETS_DELAUNAY)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	h = &boundary->hm;
	memset(&t, 0, sizeof t);
	memset(&rep, 0, sizeof rep);
	/* closed, every face's edges shared */
	for (c = 0; c < h->nc; c++)
		if (h->cf[c] != HM_NONE && h->co[c] == HM_NONE)
			return QAWS_STATUS_INVALID_TOPOLOGY;
	if (!h->nf)
		return QAWS_STATUS_INVALID_TOPOLOGY;
	/* the lattice over the box */
	for (k = 0; k < 3; k++)
	{
		t.lo[k] = HUGE_VAL;
		hi[k] = -HUGE_VAL;
	}
	for (i = 0; i < h->nv; i++)
		if (h->vc[i] != HM_NONE)
			for (k = 0; k < 3; k++)
			{
				t.lo[k] = fmin(t.lo[k], h->x[3 * (size_t)i + k]);
				hi[k] = fmax(hi[k], h->x[3 * (size_t)i + k]);
			}
	for (k = 0; k < 3; k++)
		ext = fmax(ext, hi[k] - t.lo[k]);
	if (!(ext > 0) || ext == HUGE_VAL)
		return QAWS_STATUS_INVALID_ARGUMENT;
	t.step = ext / TT_SPAN;
	m->rng = 12345u;
	/* the super cell */
	{
		static int const corner[4][3] = { { -1, -1, -1 }, { 7, -1, -1 }, { -1, 7, -1 }, { -1, -1, 7 } };
		unsigned int s;
		for (k = 0; k < 4; k++)
		{
			int64_t q[3];
			double x[3];
			for (f = 0; f < 3; f++)
			{
				q[f] = (int64_t)corner[k][f] * TT_FAR;
				x[f] = t.lo[f] + (double)q[f] * t.step;
			}
			if (tt_add_vertex(m, q, x) == HM_NONE)
				goto done;
		}
		s = tt_new_cell(m);
		if (s == HM_NONE)
			goto done;
		for (k = 0; k < 4; k++)
		{
			m->tv[4 * (size_t)s + k] = k;
			m->tn[4 * (size_t)s + k] = HM_NONE;
			m->v2t[k] = s;
		}
		m->last = s;
	}
	/* the border's vertices, then its triangles (polygons fanned) */
	map = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)h->nv + 1)));
	if (!map)
		goto done;
	for (i = 0; i < h->nv; i++)
	{
		int64_t q[3];
		unsigned int why;
		map[i] = HM_NONE;
		if (h->vc[i] == HM_NONE)
			continue;
		tt_lattice(&t, h->x + 3 * (size_t)i, q);
		map[i] = tt_insert(m, q, h->x + 3 * (size_t)i, 0, &why);
		if (map[i] == HM_NONE)
		{
			st = why == 1 ? QAWS_STATUS_INVALID_ARGUMENT : why == 0 ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_INTERNAL_ERROR;
			goto done;
		}
	}
	for (f = 0; f < h->nf; f++)
	{
		unsigned int c0 = h->fc[f];
		for (k = 1; k + 1 < h->fn[f]; k++)
			if (!tt_add_tri(&t, map[h->cv[c0]], map[h->cv[c0 + k]], map[h->cv[c0 + k + 1]]))
				goto done;
		for (k = 0; k < h->fn[f]; k++)
		{
			double const *a = h->x + 3 * (size_t)h->cv[c0 + k], *b = h->x + 3 * (size_t)h->cv[c0 + (k + 1) % h->fn[f]];
			sumlen += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
			nedge++;
		}
	}
	L = desc->target_length > 0 ? (double)desc->target_length : sumlen / (nedge ? nedge : 1);
	/* the border recovered */
	{
		int r = tt_recover(&t, &rep.iterations, &rep.splits, &rep.flips);
		if (r <= 0)
		{
			st = r == 0 ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_NOT_CONVERGED;
			goto done;
		}
	}
	/* walls on both sides of every border triangle; the inside flooded from the cells within them */
	m->ncav = 0;
	for (i = 0; i < t.nbt; i++)
	{
		unsigned int const* T = t.bt + 3 * (size_t)i;
		unsigned int oc, of, u, g, same;
		if (t.bdead[i])
			continue;
		if (!tt_find_face(m, T[0], T[1], T[2], &oc, &of))
		{
			st = QAWS_STATUS_INTERNAL_ERROR;
			goto done;
		}
		u = m->tn[4 * (size_t)oc + of];
		/* the cell whose facet turns as the triangle does is inside */
		same = 0;
		for (k = 0; k < 3; k++)
			same |= m->tv[4 * (size_t)oc + tt_face[of][k]] == T[0] && m->tv[4 * (size_t)oc + tt_face[of][(k + 1) % 3]] == T[1];
		if (u == HM_NONE)
		{
			st = QAWS_STATUS_INTERNAL_ERROR;
			goto done;
		}
		for (g = 0; g < 4 && m->tn[4 * (size_t)u + g] != oc; g++) {}
		m->tf[oc] |= (unsigned char)(1u << of);
		m->tf[u] |= (unsigned char)(1u << g);
		m->tf[same ? oc : u] |= TT_IN;
	}
	{
		unsigned int *queue = NULL, nq = 0, capq = 0;
		for (c = 0; c < m->nt; c++)
			if (!(m->tf[c] & TT_DEAD) && (m->tf[c] & TT_IN) && !tt_push(&queue, &nq, &capq, c))
			{
				qaws_internal_free(queue);
				goto done;
			}
		for (i = 0; i < nq; i++)
		{
			unsigned int cc = queue[i];
			for (k = 0; k < 4; k++)
			{
				unsigned int u = m->tn[4 * (size_t)cc + k];
				if (m->tv[4 * (size_t)cc + k] < 4)
				{
					qaws_internal_free(queue);
					st = QAWS_STATUS_INVALID_TOPOLOGY;   /* (the inside reaches the super cell: not closed) */
					goto done;
				}
				if ((m->tf[cc] & (1u << k)) || u == HM_NONE || (m->tf[u] & TT_IN))
					continue;
				m->tf[u] |= TT_IN;
				if (!tt_push(&queue, &nq, &capq, u))
				{
					qaws_internal_free(queue);
					goto done;
				}
			}
		}
		qaws_internal_free(queue);
		/* consistent: no wall with the inside on both sides */
		for (c = 0; c < m->nt; c++)
		{
			if ((m->tf[c] & TT_DEAD) || !(m->tf[c] & TT_IN))
				continue;
			for (k = 0; k < 4; k++)
				if ((m->tf[c] & (1u << k)) && m->tn[4 * (size_t)c + k] != HM_NONE && (m->tf[m->tn[4 * (size_t)c + k]] & TT_IN))
				{
					st = QAWS_STATUS_INVALID_TOPOLOGY;
					goto done;
				}
		}
	}
	/* the inside filled: a body-centred cubic lattice of spacing 1.1 L, its points farther than L / 2 from the border */
	boxes = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), ((size_t)t.nbt + 1)));
	if (!boxes)
		goto done;
	{
		unsigned int nlive = 0;
		unsigned int* live = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, ((size_t)t.nbt + 1)));
		if (!live)
			goto done;
		for (i = 0; i < t.nbt; i++)
			if (!t.bdead[i])
			{
				for (k = 0; k < 3; k++)
				{
					unsigned int v = t.bt[3 * (size_t)i + k];
					live[3 * nlive + k] = v;
					for (f = 0; f < 3; f++)
					{
						boxes[nlive].lo[f] = k ? fmin(boxes[nlive].lo[f], m->x[3 * (size_t)v + f]) : m->x[3 * (size_t)v + f];
						boxes[nlive].hi[f] = k ? fmax(boxes[nlive].hi[f], m->x[3 * (size_t)v + f]) : m->x[3 * (size_t)v + f];
					}
				}
				nlive++;
			}
		qaws_internal_free(t.bt);
		t.bt = live;
		t.nbt = nlive;
		t.capbt = nlive + 1;
		memset(t.bdead, 0, nlive);
		st = qaws_internal_grid_create(&grid, boxes, nlive, 3);
		if (st != QAWS_STATUS_OK)
			goto done;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (!desc->border_only)
	{
		double hstep = 1.1 * L;
		unsigned int n[3], a, b, d, s;
		for (k = 0; k < 3; k++)
			n[k] = (unsigned int)fmin(4096, floor((hi[k] - t.lo[k]) / hstep) + 1);
		for (s = 0; s < 2; s++)
			for (a = 0; a <= n[0]; a++)
				for (b = 0; b <= n[1]; b++)
					for (d = 0; d <= n[2]; d++)
					{
						double p[3];
						int64_t q[3];
						unsigned int why;
						p[0] = t.lo[0] + (a + 0.5 * s) * hstep;
						p[1] = t.lo[1] + (b + 0.5 * s) * hstep;
						p[2] = t.lo[2] + (d + 0.5 * s) * hstep;
						if (p[0] > hi[0] || p[1] > hi[1] || p[2] > hi[2] || tt_near_border(grid, m->x, t.bt, p, 0.5 * L))
							continue;
						tt_lattice(&t, p, q);
						if (tt_insert(m, q, p, 1, &why) != HM_NONE)
							rep.moves++;
						else if (why == 0)
							goto done;
					}
		/* refined: the circumcentres of cells too long for their shortest edge, away from the border */
		bound = desc->radius_edge == 0 ? 2 : (double)desc->radius_edge;
		if (bound > 0)
		{
			unsigned int *queue = NULL, nq = 0, capq = 0, cap = 2 * m->nv + 1000, added = 0;
			for (c = 0; c < m->nt; c++)
				if (!(m->tf[c] & TT_DEAD) && (m->tf[c] & TT_IN) && !tt_push(&queue, &nq, &capq, c))
				{
					qaws_internal_free(queue);
					goto done;
				}
			while (nq && added < cap)
			{
				unsigned int cc = queue[--nq], why;
				int64_t const* q0;
				double dv[3][3], l[3], cr[3], det, ccl[3], r2, e2 = HUGE_VAL, p[3];
				int64_t q[3];
				if ((m->tf[cc] & TT_DEAD) || !(m->tf[cc] & TT_IN))
					continue;
				q0 = TT_Q(m, m->tv[4 * (size_t)cc]);
				for (k = 0; k < 3; k++)
				{
					for (f = 0; f < 3; f++)
						dv[k][f] = (double)(TT_Q(m, m->tv[4 * (size_t)cc + k + 1])[f] - q0[f]);
					l[k] = tt_dot(dv[k], dv[k]);
				}
				tt_cross(dv[1], dv[2], cr);
				det = 2 * tt_dot(dv[0], cr);
				if (!(det > 0))
					continue;
				{
					double c23[3], c31[3], c12[3];
					tt_cross(dv[1], dv[2], c23);
					tt_cross(dv[2], dv[0], c31);
					tt_cross(dv[0], dv[1], c12);
					for (f = 0; f < 3; f++)
						ccl[f] = (l[0] * c23[f] + l[1] * c31[f] + l[2] * c12[f]) / det;
				}
				r2 = tt_dot(ccl, ccl);
				for (k = 0; k < 4; k++)
					for (f = k + 1; f < 4; f++)
					{
						int64_t const *pa = TT_Q(m, m->tv[4 * (size_t)cc + k]), *pb = TT_Q(m, m->tv[4 * (size_t)cc + f]);
						double dd[3];
						unsigned int z;
						for (z = 0; z < 3; z++)
							dd[z] = (double)(pa[z] - pb[z]);
						e2 = fmin(e2, tt_dot(dd, dd));
					}
				if (r2 <= bound * bound * e2)
					continue;
				for (f = 0; f < 3; f++)
				{
					q[f] = q0[f] + (int64_t)floor(ccl[f] + 0.5);
					p[f] = t.lo[f] + (double)q[f] * t.step;
				}
				if (tt_near_border(grid, m->x, t.bt, p, 0.4 * L))
					continue;
				if (tt_insert(m, q, p, 1, &why) == HM_NONE)
				{
					if (why == 0)
					{
						qaws_internal_free(queue);
						goto done;
					}
					continue;
				}
				added++;
				rep.moves++;
				for (i = 0; i < m->nmade; i++)
					if (!tt_push(&queue, &nq, &capq, m->made[i]))
					{
						qaws_internal_free(queue);
						goto done;
					}
			}
			qaws_internal_free(queue);
		}
	}
	{
		/* improved: flips, then the inside vertices smoothed (the border's fixed), three rounds */
		unsigned char* fixed = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, 1);
		unsigned int* arr = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)m->nt + 1)));
		unsigned int round, n;
		int ok = fixed && arr;
		for (i = 0; ok && i < 3 * t.nbt; i++)
			fixed[t.bt[i]] = 1;
		for (round = 0; ok && round < 3; round++)
		{
			ok = tt_flip_pass(m, 0.42);
			for (c = 0, n = 0; ok && c < m->nt; c++)
				if (!(m->tf[c] & TT_DEAD) && (m->tf[c] & TT_IN))
					memcpy(arr + 4 * (size_t)n++, m->tv + 4 * (size_t)c, sizeof(unsigned int) * 4);
			ok = ok && tt_smooth(m->x, m->nv, arr, n, fixed, 2);
		}
		qaws_internal_free(fixed);
		qaws_internal_free(arr);
		if (!ok)
			goto done;
	}
	/* the cells inside, their vertices compacted */
	for (c = 0; c < m->nt; c++)
		nin += !(m->tf[c] & TT_DEAD) && (m->tf[c] & TT_IN);
	vmap = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->nv + 1)));
	cells = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)nin + 1)));
	ctype = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nin + 1, 1);
	xo = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)m->nv + 1)));
	if (!vmap || !cells || !ctype || !xo)
		goto done;
	for (i = 0; i < m->nv; i++)
		vmap[i] = HM_NONE;
	for (c = 0, nin = 0; c < m->nt; c++)
	{
		if ((m->tf[c] & TT_DEAD) || !(m->tf[c] & TT_IN))
			continue;
		for (k = 0; k < 4; k++)
		{
			unsigned int v = m->tv[4 * (size_t)c + k];
			if (vmap[v] == HM_NONE)
			{
				memcpy(xo + 3 * (size_t)ninv, m->x + 3 * (size_t)v, sizeof(double) * 3);
				vmap[v] = ninv++;
			}
			cells[4 * (size_t)nin + k] = vmap[v];
		}
		nin++;
	}
	st = vm_create(xo, ninv, cells, ctype, nin, out_mesh);
	if (st == QAWS_STATUS_OK && desc->keep_border && rep.splits)
	{
		/* the recovery's border points (on the border, not the input's) moved inside, then the cells improved */
		unsigned char *in = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, 1), *steiner = (unsigned char*)qaws_internal_scratch_zalloc((size_t)ninv + 1, 1);
		qaws_volume_mesh *u = NULL, *w = NULL;
		unsigned int left = 0;
		if (!in || !steiner)
		{
			qaws_internal_free(in);
			qaws_internal_free(steiner);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (i = 0; i < h->nv; i++)
			if (map[i] != HM_NONE)
				in[map[i]] = 1;
		for (i = 0; i < 3 * t.nbt; i++)
			if (!in[t.bt[i]] && vmap[t.bt[i]] != HM_NONE)
				steiner[vmap[t.bt[i]]] = 1;
		st = vim_unsplit_border(&u, &rep.collapses, &left, *out_mesh, steiner);
		if (st == QAWS_STATUS_OK)
			st = qaws_volume_mesh_improve(&w, NULL, u, NULL);
		qaws_volume_mesh_destroy(u);
		qaws_internal_free(in);
		qaws_internal_free(steiner);
		if (st == QAWS_STATUS_OK)
		{
			qaws_volume_mesh_destroy(*out_mesh);
			*out_mesh = w;
		}
	}
	if (st == QAWS_STATUS_OK && qaws_mesh_3d_get_attribute_count(boundary))
		st = vm_attr_from_surface(*out_mesh, boundary);
	if (st == QAWS_STATUS_OK)
	{
		qaws_mesh_3d* b = NULL;
		if (qaws_volume_mesh_get_boundary(&b, *out_mesh) == QAWS_STATUS_OK)
			qaws_mesh_3d_measure_quality(&rep.quality, b);
		qaws_mesh_3d_destroy(b);
	}
done:
	qaws_internal_free(map); qaws_internal_free(cells); qaws_internal_free(vmap); qaws_internal_free(ctype); qaws_internal_free(xo); qaws_internal_free(boxes);
	qaws_internal_free(t.bt); qaws_internal_free(t.bdead);
	tt_emap_free(&t.bmap);
	qaws_internal_grid_destroy(grid);
	tt_free(m);
	if (st != QAWS_STATUS_OK)
	{
		qaws_volume_mesh_destroy(*out_mesh);
		*out_mesh = NULL;
		return st;
	}
	if (out_report)
		*out_report = rep;
	return QAWS_STATUS_OK;
}
