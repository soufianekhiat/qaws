/*
 * Repair for farq (qaws_mesh_repair.h). The stages work on plain arrays of
 * triangles (the corner store splits non-manifold input, which a soup is
 * until welded); the result goes into the store at the end. Edges are
 * found by sorting the half-edges by their two vertices.
 */

#include "qaws_mesh_core.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_repair.h"
#include "../internal/qaws_internal_broadphase.h"
#include "../internal/qaws_internal_ldlt.h"
#include "../qaws_exact.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RP_PI 3.14159265358979323846
#define RP_DP_MAX 500u          /* loops filled by the dynamic programme (n^3 / 6 steps, 7 n^2 doubles); longer ones by a fan */

typedef struct rp_mesh
{
	double* x;                  /* 3 per vertex */
	unsigned int nv, capv;
	unsigned int* t;            /* 3 per triangle */
	int* src;                   /* per triangle: the given face, -1 for a fill */
	unsigned char* fl;          /* per triangle: turned from the given orientation */
	unsigned int nt, capt;
} rp_mesh;

typedef struct rp_work
{
	rp_mesh m;
	qaws_mesh_repair_report* r;
	qaws_status st;
} rp_work;

/* ---- vectors ------------------------------------------------------------ */

static void mrp_sub(double const* a, double const* b, double* o)
{
	o[0] = a[0] - b[0];
	o[1] = a[1] - b[1];
	o[2] = a[2] - b[2];
}

static double mrp_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void rp_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double rp_dist(double const* a, double const* b)
{
	double d[3];
	mrp_sub(a, b, d);
	return sqrt(mrp_dot(d, d));
}

/* the normal of (a, b, c), as long as twice the area */
static void rp_normal(double const* a, double const* b, double const* c, double* n)
{
	double u[3], v[3];
	mrp_sub(b, a, u);
	mrp_sub(c, a, v);
	rp_cross(u, v, n);
}

/* n made unit (zero when it has no length); its length */
static double rp_unit(double* n)
{
	double l = sqrt(mrp_dot(n, n));
	n[0] = l > 0 ? n[0] / l : 0;
	n[1] = l > 0 ? n[1] / l : 0;
	n[2] = l > 0 ? n[2] / l : 0;
	return l;
}

/* the angle between u and v */
static double rp_angle(double const* u, double const* v)
{
	double c[3];
	rp_cross(u, v, c);
	return atan2(sqrt(mrp_dot(c, c)), mrp_dot(u, v));
}

static int rp_finite(double v)
{
	return v - v == 0;
}

#define RP_X(m, v) ((m)->x + 3 * (size_t)(v))

/* ---- storage ------------------------------------------------------------ */

static unsigned int rp_add_v(rp_mesh* m, double const* p)
{
	if (m->nv == m->capv)
	{
		unsigned int cap = qaws_internal_grow(m->capv, (size_t)m->nv + 1, 64);
		double* x = cap ? (double*)qaws_internal_realloc(m->x, qaws_internal_mul3(sizeof(double), 3, cap)) : NULL;
		if (!x)
			return HM_NONE;
		m->x = x;
		m->capv = cap;
	}
	memcpy(RP_X(m, m->nv), p, sizeof(double) * 3);
	return m->nv++;
}

static int rp_add_t(rp_mesh* m, unsigned int a, unsigned int b, unsigned int c, int src, unsigned char fl)
{
	if (m->nt == m->capt)
	{
		unsigned int cap = qaws_internal_grow(m->capt, (size_t)m->nt + 1, 64);
		unsigned int* t = cap ? (unsigned int*)qaws_internal_realloc(m->t, qaws_internal_mul3(sizeof(unsigned int), 3, cap)) : NULL;
		int* s;
		unsigned char* f;
		if (!t)
			return 0;
		m->t = t;
		s = (int*)qaws_internal_realloc(m->src, qaws_internal_mul2(sizeof(int), cap));
		if (!s)
			return 0;
		m->src = s;
		f = (unsigned char*)qaws_internal_realloc(m->fl, cap);
		if (!f)
			return 0;
		m->fl = f;
		m->capt = cap;
	}
	m->t[3 * (size_t)m->nt] = a;
	m->t[3 * (size_t)m->nt + 1] = b;
	m->t[3 * (size_t)m->nt + 2] = c;
	m->src[m->nt] = src;
	m->fl[m->nt++] = fl;
	return 1;
}

/* the triangles below n with del set removed, the others kept in order */
static void rp_compact(rp_mesh* m, unsigned char const* del, unsigned int n)
{
	unsigned int f, k = 0;
	for (f = 0; f < m->nt; f++)
	{
		if (f < n && del[f])
			continue;
		memmove(m->t + 3 * (size_t)k, m->t + 3 * (size_t)f, sizeof(unsigned int) * 3);
		m->src[k] = m->src[f];
		m->fl[k++] = m->fl[f];
	}
	m->nt = k;
}

static void rp_flip(rp_mesh* m, unsigned int f)
{
	unsigned int s = m->t[3 * (size_t)f + 1];
	m->t[3 * (size_t)f + 1] = m->t[3 * (size_t)f + 2];
	m->t[3 * (size_t)f + 2] = s;
	m->fl[f] ^= 1;
}

/* ---- edges: half-edges sorted by their vertices --------------------------- */

typedef struct rp_he
{
	unsigned long long key;     /* lower vertex << 32 | higher */
	unsigned int f, k;          /* the triangle, its corner (the edge from corner k to k + 1) */
} rp_he;

static unsigned long long rp_key(unsigned int a, unsigned int b)
{
	return a < b ? (unsigned long long)a << 32 | b : (unsigned long long)b << 32 | a;
}

static int rp_he_cmp(void const* p, void const* q)
{
	rp_he const* a = (rp_he const*)p;
	rp_he const* b = (rp_he const*)q;
	if (a->key != b->key)
		return a->key < b->key ? -1 : 1;
	if (a->f != b->f)
		return a->f < b->f ? -1 : 1;
	return a->k < b->k ? -1 : a->k > b->k;
}

/* the half-edges of triangles [t0, t1), sorted; NULL on allocation failure */
static rp_he* rp_edges(rp_mesh const* m, unsigned int t0, unsigned int t1, unsigned int* out_n)
{
	unsigned int n = 3 * (t1 - t0), f, k;
	rp_he* E = (rp_he*)qaws_internal_malloc(qaws_internal_mul2(sizeof(rp_he), (n ? n : 1)));
	*out_n = n;
	if (!E)
		return NULL;
	for (f = t0; f < t1; f++)
		for (k = 0; k < 3; k++)
		{
			rp_he* e = E + 3 * (size_t)(f - t0) + k;
			e->key = rp_key(m->t[3 * (size_t)f + k], m->t[3 * (size_t)f + (k + 1) % 3]);
			e->f = f;
			e->k = k;
		}
	qsort(E, n, sizeof(rp_he), rp_he_cmp);
	return E;
}

/* the first record of the edge (n when it has none) */
static unsigned int rp_lower(rp_he const* E, unsigned int n, unsigned long long key)
{
	unsigned int lo = 0, hi = n;
	while (lo < hi)
	{
		unsigned int mid = lo + (hi - lo) / 2;
		if (E[mid].key < key)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo;
}

static int rp_has_edge(rp_he const* E, unsigned int n, unsigned int a, unsigned int b)
{
	unsigned long long key = rp_key(a, b);
	unsigned int i = rp_lower(E, n, key);
	return i < n && E[i].key == key;
}

/* 1 when the half-edge runs from its lower vertex to its higher */
static int rp_dir(rp_mesh const* m, rp_he const* e)
{
	return m->t[3 * (size_t)e->f + e->k] < m->t[3 * (size_t)e->f + (e->k + 1) % 3];
}

/* the end of the group of records r is in */
static unsigned int rp_group_end(rp_he const* E, unsigned int n, unsigned int r)
{
	unsigned int e = r + 1;
	while (e < n && E[e].key == E[r].key)
		e++;
	return e;
}

/* ---- merging ------------------------------------------------------------ */

static unsigned int rp_find(unsigned int* uf, unsigned int i)
{
	while (uf[i] != i)
	{
		uf[i] = uf[uf[i]];
		i = uf[i];
	}
	return i;
}

/* the lower root kept: a group's root is its least vertex whatever the order of the unions */
static void rp_union(unsigned int* uf, unsigned int a, unsigned int b)
{
	a = rp_find(uf, a);
	b = rp_find(uf, b);
	if (a < b)
		uf[b] = a;
	else if (b < a)
		uf[a] = b;
}

/* each group of uf moved to its mean, held by its root; the triangles' vertices mapped. Returns the vertices merged away. */
static unsigned int rp_merge(rp_work* w, unsigned int* uf)
{
	rp_mesh* m = &w->m;
	double* s = (double*)qaws_internal_zalloc(qaws_internal_mul2(4, ((size_t)(m->nv ? m->nv : 1))), sizeof(double));
	unsigned int v, n = 0;
	size_t i;
	if (!s)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return 0;
	}
	for (v = 0; v < m->nv; v++)
	{
		unsigned int r = rp_find(uf, v);
		double* g = s + 4 * (size_t)r;
		g[0] += RP_X(m, v)[0];
		g[1] += RP_X(m, v)[1];
		g[2] += RP_X(m, v)[2];
		g[3] += 1;
		n += r != v;
	}
	for (v = 0; v < m->nv; v++)
		if (uf[v] == v && s[4 * (size_t)v + 3] > 1)
		{
			RP_X(m, v)[0] = s[4 * (size_t)v] / s[4 * (size_t)v + 3];
			RP_X(m, v)[1] = s[4 * (size_t)v + 1] / s[4 * (size_t)v + 3];
			RP_X(m, v)[2] = s[4 * (size_t)v + 2] / s[4 * (size_t)v + 3];
		}
	for (i = 0; i < 3 * (size_t)m->nt; i++)
		m->t[i] = rp_find(uf, m->t[i]);
	qaws_internal_free(s);
	return n;
}

typedef struct rp_tri_key
{
	unsigned int v[3], f;
	unsigned int odd;           /* turned against its sorted order */
} rp_tri_key;

static int rp_tri_cmp(void const* p, void const* q)
{
	rp_tri_key const* a = (rp_tri_key const*)p;
	rp_tri_key const* b = (rp_tri_key const*)q;
	unsigned int k;
	for (k = 0; k < 3; k++)
		if (a->v[k] != b->v[k])
			return a->v[k] < b->v[k] ? -1 : 1;
	return a->f < b->f ? -1 : a->f > b->f;
}

/* triangles with a vertex twice dropped; triangles on the same three vertices kept once (the majority's turn, none on a tie) */
static void rp_clean(rp_work* w)
{
	rp_mesh* m = &w->m;
	unsigned int f, n = 0, i, j;
	rp_tri_key* K;
	unsigned char* del;
	for (f = 0; f < m->nt; f++)
	{
		unsigned int a = m->t[3 * (size_t)f], b = m->t[3 * (size_t)f + 1], c = m->t[3 * (size_t)f + 2];
		if (a == b || b == c || c == a)
		{
			w->r->dropped_faces++;
			continue;
		}
		m->t[3 * (size_t)n] = a;
		m->t[3 * (size_t)n + 1] = b;
		m->t[3 * (size_t)n + 2] = c;
		m->src[n] = m->src[f];
		m->fl[n++] = m->fl[f];
	}
	m->nt = n;
	K = (rp_tri_key*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(rp_tri_key), (n ? n : 1)));
	del = (unsigned char*)qaws_internal_scratch_zalloc(n ? n : 1, 1);
	if (!K || !del)
	{
		qaws_internal_free(K);
		qaws_internal_free(del);
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	for (f = 0; f < n; f++)
	{
		unsigned int const* t = m->t + 3 * (size_t)f;
		unsigned int lo = t[0] < t[1] ? (t[0] < t[2] ? 0 : 2) : (t[1] < t[2] ? 1 : 2);
		unsigned int a = t[lo], b = t[(lo + 1) % 3], c = t[(lo + 2) % 3];
		K[f].v[0] = a;
		K[f].v[1] = b < c ? b : c;
		K[f].v[2] = b < c ? c : b;
		K[f].odd = b > c;
		K[f].f = f;
	}
	qsort(K, n, sizeof(rp_tri_key), rp_tri_cmp);
	for (i = 0; i < n; i = j)
	{
		unsigned int odd = 0, keep = HM_NONE, k;
		for (j = i + 1; j < n && !memcmp(K[j].v, K[i].v, sizeof K[i].v); j++) {}
		if (j - i == 1)
			continue;
		for (k = i; k < j; k++)
			odd += K[k].odd;
		if (2 * odd != j - i)
			for (k = i; k < j && keep == HM_NONE; k++)
				if (K[k].odd == (unsigned int)(2 * odd > j - i))
					keep = K[k].f;
		for (k = i; k < j; k++)
			if (K[k].f != keep)
			{
				del[K[k].f] = 1;
				w->r->duplicate_faces++;
			}
	}
	rp_compact(m, del, n);
	qaws_internal_free(K);
	qaws_internal_free(del);
}

/* ---- welding ------------------------------------------------------------ */

typedef struct rp_weld_ctx
{
	double const* x;
	unsigned int const* id;     /* box: its vertex */
	unsigned int* uf;
	double tol2;
} rp_weld_ctx;

static qaws_status rp_weld_visit(void* user, unsigned int i, unsigned int j)
{
	rp_weld_ctx* c = (rp_weld_ctx*)user;
	double d[3];
	mrp_sub(c->x + 3 * (size_t)c->id[i], c->x + 3 * (size_t)c->id[j], d);
	if (mrp_dot(d, d) <= c->tol2)
		rp_union(c->uf, c->id[i], c->id[j]);
	return QAWS_STATUS_OK;
}

/* the used points within tol of each other merged; the union-find returned (NULL on failure) */
static unsigned int* rp_weld(rp_work* w, unsigned char const* used, double tol)
{
	rp_mesh* m = &w->m;
	unsigned int* uf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv ? m->nv : 1)));
	unsigned int* id = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv ? m->nv : 1)));
	qaws_bp_box* box = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (m->nv ? m->nv : 1)));
	unsigned int v, n = 0, k;
	rp_weld_ctx c;
	if (!uf || !id || !box)
		goto fail;
	for (v = 0; v < m->nv; v++)
	{
		uf[v] = v;
		if (!used[v])
			continue;
		for (k = 0; k < 3; k++)
		{
			box[n].lo[k] = RP_X(m, v)[k] - 0.5 * tol;
			box[n].hi[k] = RP_X(m, v)[k] + 0.5 * tol;
		}
		id[n++] = v;
	}
	c.x = m->x;
	c.id = id;
	c.uf = uf;
	c.tol2 = tol * tol;
	if (n > 1 && qaws_internal_broadphase(NULL, box, n, 3, NULL, rp_weld_visit, &c) != QAWS_STATUS_OK)
		goto fail;
	w->r->welded_vertices += rp_merge(w, uf);
	qaws_internal_free(id);
	qaws_internal_free(box);
	return uf;
fail:
	qaws_internal_free(uf);
	qaws_internal_free(id);
	qaws_internal_free(box);
	w->st = QAWS_STATUS_ALLOCATION_FAILURE;
	return NULL;
}

/* ---- polygons ----------------------------------------------------------- */

static double rp_orient2(double const* a, double const* b, double const* c)
{
	return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

/* is corner e of the cnt corners left (id) an ear: convex, no other corner in its triangle */
static int rp_is_ear(double const* uv, unsigned int const* id, unsigned int cnt, unsigned int e)
{
	unsigned int p = id[(e + cnt - 1) % cnt], c = id[e], q = id[(e + 1) % cnt], j;
	double const *A = uv + 2 * (size_t)p, *B = uv + 2 * (size_t)c, *C = uv + 2 * (size_t)q;
	if (rp_orient2(A, B, C) <= 0)
		return 0;
	for (j = 0; j < cnt; j++)
	{
		unsigned int o = id[j];
		double const* X = uv + 2 * (size_t)o;
		if (o != p && o != c && o != q && rp_orient2(A, B, X) >= 0 && rp_orient2(B, C, X) >= 0 && rp_orient2(C, A, X) >= 0)
			return 0;
	}
	return 1;
}

/* polygon P (n > 3 vertices) cut into triangles by ear clipping in the plane of its Newell normal; where no ear is found
   (a polygon crossing itself) the first corner is cut anyway.
   shortcut: n^2 per ear search, slow past thousands of corners; upgrade to a monotone decomposition if such faces come up */
static void rp_ear_clip(rp_work* w, unsigned int const* P, unsigned int n, int src)
{
	rp_mesh* m = &w->m;
	double N[3] = { 0, 0, 0 };
	double* uv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 2, n));
	unsigned int* id = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), n));
	unsigned int i, cnt = n, ax = 0, k;
	if (!uv || !id)
	{
		qaws_internal_free(uv);
		qaws_internal_free(id);
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	for (i = 0; i < n; i++)
	{
		double const *a = RP_X(m, P[i]), *b = RP_X(m, P[(i + 1) % n]);
		N[0] += (a[1] - b[1]) * (a[2] + b[2]);
		N[1] += (a[2] - b[2]) * (a[0] + b[0]);
		N[2] += (a[0] - b[0]) * (a[1] + b[1]);
	}
	for (k = 1; k < 3; k++)
		if (fabs(N[k]) > fabs(N[ax]))
			ax = k;
	for (i = 0; i < n; i++)
	{
		double const* a = RP_X(m, P[i]);
		uv[2 * i] = a[(ax + 1) % 3];
		uv[2 * i + 1] = N[ax] < 0 ? -a[(ax + 2) % 3] : a[(ax + 2) % 3];
		id[i] = i;
	}
	while (cnt > 3)
	{
		unsigned int e;
		for (e = 0; e < cnt && !rp_is_ear(uv, id, cnt, e); e++) {}
		if (e == cnt)
			e = 0;
		if (!rp_add_t(m, P[id[(e + cnt - 1) % cnt]], P[id[e]], P[id[(e + 1) % cnt]], src, 0))
			break;
		memmove(id + e, id + e + 1, sizeof(unsigned int) * (cnt - e - 1));
		cnt--;
	}
	if (cnt == 3 && !rp_add_t(m, P[id[0]], P[id[1]], P[id[2]], src, 0))
		cnt = 0;
	if (cnt != 3)
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
	qaws_internal_free(uv);
	qaws_internal_free(id);
}

/* ---- degenerate faces --------------------------------------------------- */

/* faces lower than h: a needle's short edge collapsed; a cap removed and the faces across its long edge split at its corner */
static void rp_degenerate(rp_work* w, double h)
{
	rp_mesh* m = &w->m;
	unsigned int pass;
	for (pass = 0; pass < 16 && h > 0 && w->st == QAWS_STATUS_OK; pass++)
	{
		unsigned int ne, nt0 = m->nt, f, changed = 0, v;
		rp_he* E = rp_edges(m, 0, nt0, &ne);
		unsigned int* uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv ? m->nv : 1)));
		unsigned char* mark = (unsigned char*)qaws_internal_zalloc(nt0 ? nt0 : 1, 1);   /* 1 touched this pass, 2 removed */
		if (!E || !uf || !mark)
		{
			qaws_internal_free(E);
			qaws_internal_free(uf);
			qaws_internal_free(mark);
			w->st = QAWS_STATUS_ALLOCATION_FAILURE;
			return;
		}
		for (v = 0; v < m->nv; v++)
			uf[v] = v;
		for (f = 0; f < nt0 && w->st == QAWS_STATUS_OK; f++)
		{
			double l[3], n[3], lmax = 0, lmin = HUGE_VAL;
			unsigned int k, kmax = 0, kmin = 0, a, b, p, r, r1;
			int ok = 1;
			if (mark[f])
				continue;
			for (k = 0; k < 3; k++)
			{
				l[k] = rp_dist(RP_X(m, m->t[3 * (size_t)f + k]), RP_X(m, m->t[3 * (size_t)f + (k + 1) % 3]));
				if (l[k] > lmax)
				{
					lmax = l[k];
					kmax = k;
				}
				if (l[k] < lmin)
				{
					lmin = l[k];
					kmin = k;
				}
			}
			rp_normal(RP_X(m, m->t[3 * (size_t)f]), RP_X(m, m->t[3 * (size_t)f + 1]), RP_X(m, m->t[3 * (size_t)f + 2]), n);
			if (lmax > 0 && sqrt(mrp_dot(n, n)) >= h * lmax)
				continue;
			if (lmin < h)
			{
				rp_union(uf, m->t[3 * (size_t)f + kmin], m->t[3 * (size_t)f + (kmin + 1) % 3]);
				mark[f] = 1;
				changed++;
				w->r->needles++;
				continue;
			}
			a = m->t[3 * (size_t)f + kmax];
			b = m->t[3 * (size_t)f + (kmax + 1) % 3];
			p = m->t[3 * (size_t)f + (kmax + 2) % 3];
			r = rp_lower(E, ne, rp_key(a, b));
			r1 = rp_group_end(E, ne, r);
			for (k = r; k < r1; k++)
				if (E[k].f != f && mark[E[k].f])
					ok = 0;
			if (!ok)
				continue;
			for (; r < r1; r++)
			{
				unsigned int g = E[r].f, kg = E[r].k, x, y, q;
				if (g == f)
					continue;
				x = m->t[3 * (size_t)g + kg];
				y = m->t[3 * (size_t)g + (kg + 1) % 3];
				q = m->t[3 * (size_t)g + (kg + 2) % 3];
				m->t[3 * (size_t)g] = x;
				m->t[3 * (size_t)g + 1] = p;
				m->t[3 * (size_t)g + 2] = q;
				if (!rp_add_t(m, p, y, q, m->src[g], m->fl[g]))
					w->st = QAWS_STATUS_ALLOCATION_FAILURE;
				mark[g] = 1;
			}
			mark[f] = 2;
			changed++;
			w->r->caps++;
		}
		if (changed && w->st == QAWS_STATUS_OK)
		{
			for (f = 0; f < nt0; f++)
				mark[f] = mark[f] == 2;
			rp_compact(m, mark, nt0);
			rp_merge(w, uf);
			if (w->st == QAWS_STATUS_OK)
				rp_clean(w);
		}
		qaws_internal_free(E);
		qaws_internal_free(uf);
		qaws_internal_free(mark);
		if (!changed)
			break;
	}
}

/* ---- gaps --------------------------------------------------------------- */

typedef struct rp_stitch_ctx
{
	rp_mesh const* m;
	rp_he const* E;
	unsigned int ne, nbv;
	unsigned int const* bv;     /* the boundary vertices (boxes [0, nbv)) */
	unsigned int const* be;     /* the boundary edges' records (boxes from nbv) */
	unsigned int* uf;
	unsigned char* welded;      /* per vertex: merged with another boundary vertex */
	unsigned int* best;         /* per boundary vertex: the boundary edge it snaps onto, HM_NONE */
	double *bs, *bd;            /* its parameter on the edge, its distance */
	double tol;
} rp_stitch_ctx;

static int rp_stitch_accept(void* user, unsigned int i, unsigned int j)
{
	(void)j;
	return i < ((rp_stitch_ctx const*)user)->nbv;
}

static qaws_status rp_stitch_visit(void* user, unsigned int i, unsigned int j)
{
	rp_stitch_ctx* c = (rp_stitch_ctx*)user;
	rp_mesh const* m = c->m;
	unsigned int v = c->bv[i];
	double const* P = RP_X(m, v);
	if (j < c->nbv)
	{
		unsigned int u = c->bv[j];
		if (rp_dist(P, RP_X(m, u)) <= c->tol && !rp_has_edge(c->E, c->ne, v, u))
		{
			rp_union(c->uf, v, u);
			c->welded[v] = c->welded[u] = 1;
		}
	}
	else
	{
		rp_he const* e = c->E + c->be[j - c->nbv];
		unsigned int a = m->t[3 * (size_t)e->f + e->k], b = m->t[3 * (size_t)e->f + (e->k + 1) % 3], q = m->t[3 * (size_t)e->f + (e->k + 2) % 3];
		double const *A = RP_X(m, a), *B = RP_X(m, b);
		double ab[3], ap[3], F[3], L2, s, d;
		if (v == a || v == b || v == q || rp_has_edge(c->E, c->ne, v, a) || rp_has_edge(c->E, c->ne, v, b)
			|| rp_dist(P, A) <= c->tol || rp_dist(P, B) <= c->tol)
			return QAWS_STATUS_OK;   /* (near an end: merged with it instead) */
		mrp_sub(B, A, ab);
		mrp_sub(P, A, ap);
		L2 = mrp_dot(ab, ab);
		s = L2 > 0 ? mrp_dot(ap, ab) / L2 : 0;
		if (s <= 0 || s >= 1)
			return QAWS_STATUS_OK;
		F[0] = A[0] + s * ab[0];
		F[1] = A[1] + s * ab[1];
		F[2] = A[2] + s * ab[2];
		d = rp_dist(P, F);
		if (d <= c->tol && (c->best[i] == HM_NONE || d < c->bd[i] || (d == c->bd[i] && j - c->nbv < c->best[i])))
		{
			c->best[i] = j - c->nbv;
			c->bs[i] = s;
			c->bd[i] = d;
		}
	}
	return QAWS_STATUS_OK;
}

typedef struct rp_snap
{
	unsigned int e, v;
	double s;
} rp_snap;

static int rp_snap_cmp(void const* p, void const* q)
{
	rp_snap const* a = (rp_snap const*)p;
	rp_snap const* b = (rp_snap const*)q;
	if (a->e != b->e)
		return a->e < b->e ? -1 : 1;
	if (a->s != b->s)
		return a->s < b->s ? -1 : 1;
	return a->v < b->v ? -1 : a->v > b->v;
}

/* boundary vertices within tol of another boundary vertex merged, or of a boundary edge: its face split there; until none */
static void rp_stitch(rp_work* w, double tol)
{
	rp_mesh* m = &w->m;
	unsigned int pass;
	for (pass = 0; pass < 8 && w->st == QAWS_STATUS_OK; pass++)
	{
		unsigned int ne, nbe = 0, nbv = 0, r, v, i, nsn = 0, snapped = 0, merged, k, nt0 = m->nt;
		rp_he* E = rp_edges(m, 0, nt0, &ne);
		unsigned int *be = NULL, *bv = NULL, *uf = NULL, *best = NULL;
		unsigned char *isb = NULL, *welded = NULL, *used = NULL;
		double *bs = NULL, *bd = NULL;
		qaws_bp_box* box = NULL;
		rp_snap* sn = NULL;
		rp_stitch_ctx c;
		if (!E)
			goto oom;
		be = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (ne ? ne : 1)));
		isb = (unsigned char*)qaws_internal_scratch_zalloc(m->nv ? m->nv : 1, 1);
		if (!be || !isb)
			goto oom;
		for (r = 0; r < ne; r = k)
		{
			k = rp_group_end(E, ne, r);
			if (k - r == 1)
			{
				be[nbe++] = r;
				isb[(unsigned int)(E[r].key >> 32)] = isb[(unsigned int)E[r].key] = 1;
			}
		}
		if (!nbe)
		{
			qaws_internal_free(E);
			qaws_internal_free(be);
			qaws_internal_free(isb);
			break;
		}
		for (v = 0; v < m->nv; v++)
			nbv += isb[v];
		bv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), nbv));
		uf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), m->nv));
		welded = (unsigned char*)qaws_internal_zalloc(m->nv, 1);
		best = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), nbv));
		bs = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nbv));
		bd = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nbv));
		box = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (nbv + nbe)));
		sn = (rp_snap*)qaws_internal_malloc(qaws_internal_mul2(sizeof(rp_snap), nbv));
		used = (unsigned char*)qaws_internal_scratch_zalloc(nt0, 1);
		if (!bv || !uf || !welded || !best || !bs || !bd || !box || !sn || !used)
			goto oom;
		for (v = 0, i = 0; v < m->nv; v++)
		{
			uf[v] = v;
			if (!isb[v])
				continue;
			for (k = 0; k < 3; k++)
			{
				box[i].lo[k] = RP_X(m, v)[k] - 0.5 * tol;
				box[i].hi[k] = RP_X(m, v)[k] + 0.5 * tol;
			}
			best[i] = HM_NONE;
			bv[i++] = v;
		}
		for (i = 0; i < nbe; i++)
		{
			double const* A = RP_X(m, (unsigned int)(E[be[i]].key >> 32));
			double const* B = RP_X(m, (unsigned int)E[be[i]].key);
			for (k = 0; k < 3; k++)
			{
				box[nbv + i].lo[k] = fmin(A[k], B[k]) - 0.5 * tol;
				box[nbv + i].hi[k] = fmax(A[k], B[k]) + 0.5 * tol;
			}
		}
		c.m = m;
		c.E = E;
		c.ne = ne;
		c.nbv = nbv;
		c.bv = bv;
		c.be = be;
		c.uf = uf;
		c.welded = welded;
		c.best = best;
		c.bs = bs;
		c.bd = bd;
		c.tol = tol;
		if (qaws_internal_broadphase(NULL, box, nbv + nbe, 3, rp_stitch_accept, rp_stitch_visit, &c) != QAWS_STATUS_OK)
			goto oom;
		for (i = 0; i < nbv; i++)
			if (best[i] != HM_NONE && !welded[bv[i]])
			{
				sn[nsn].e = best[i];
				sn[nsn].v = bv[i];
				sn[nsn++].s = bs[i];
			}
		qsort(sn, nsn, sizeof(rp_snap), rp_snap_cmp);
		/* each face split along one edge a pass: (a, v1, q), (v1, v2, q), ..., (vk, b, q) */
		for (i = 0; i < nsn; i = k)
		{
			rp_he const* e = E + be[sn[i].e];
			unsigned int f = e->f, a, b, q, j;
			for (k = i + 1; k < nsn && sn[k].e == sn[i].e; k++) {}
			if (used[f])
				continue;
			used[f] = 1;
			a = m->t[3 * (size_t)f + e->k];
			b = m->t[3 * (size_t)f + (e->k + 1) % 3];
			q = m->t[3 * (size_t)f + (e->k + 2) % 3];
			m->t[3 * (size_t)f] = a;
			m->t[3 * (size_t)f + 1] = sn[i].v;
			m->t[3 * (size_t)f + 2] = q;
			for (j = i; j < k && w->st == QAWS_STATUS_OK; j++)
				if (!rp_add_t(m, sn[j].v, j + 1 < k ? sn[j + 1].v : b, q, m->src[f], m->fl[f]))
					w->st = QAWS_STATUS_ALLOCATION_FAILURE;
			snapped += k - i;
		}
		merged = w->st == QAWS_STATUS_OK ? rp_merge(w, uf) : 0;
		w->r->stitched_vertices += snapped + merged;
		if (w->st == QAWS_STATUS_OK)
			rp_clean(w);
		qaws_internal_free(E); qaws_internal_free(be); qaws_internal_free(isb); qaws_internal_free(bv); qaws_internal_free(uf); qaws_internal_free(welded); qaws_internal_free(best); qaws_internal_free(bs); qaws_internal_free(bd); qaws_internal_free(box); qaws_internal_free(sn); qaws_internal_free(used);
		if (!snapped && !merged)
			break;
		continue;
	oom:
		qaws_internal_free(E); qaws_internal_free(be); qaws_internal_free(isb); qaws_internal_free(bv); qaws_internal_free(uf); qaws_internal_free(welded); qaws_internal_free(best); qaws_internal_free(bs); qaws_internal_free(bd); qaws_internal_free(box); qaws_internal_free(sn); qaws_internal_free(used);
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
}

/* ---- small components --------------------------------------------------- */

/* parts (faces sharing vertices) with less area than min_area or fewer faces than min_faces removed */
static void rp_components(rp_work* w, double min_area, unsigned int min_faces)
{
	rp_mesh* m = &w->m;
	unsigned int* uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv ? m->nv : 1)));
	unsigned int* cnt = (unsigned int*)qaws_internal_scratch_zalloc(m->nv ? m->nv : 1, sizeof(unsigned int));
	double* area = (double*)qaws_internal_scratch_zalloc(m->nv ? m->nv : 1, sizeof(double));
	unsigned char* del = (unsigned char*)qaws_internal_scratch_zalloc(m->nt ? m->nt : 1, 1);
	unsigned int v, f;
	if (!uf || !cnt || !area || !del)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (v = 0; v < m->nv; v++)
		uf[v] = v;
	for (f = 0; f < m->nt; f++)
	{
		rp_union(uf, m->t[3 * (size_t)f], m->t[3 * (size_t)f + 1]);
		rp_union(uf, m->t[3 * (size_t)f], m->t[3 * (size_t)f + 2]);
	}
	for (f = 0; f < m->nt; f++)
	{
		unsigned int r = rp_find(uf, m->t[3 * (size_t)f]);
		double n[3];
		rp_normal(RP_X(m, m->t[3 * (size_t)f]), RP_X(m, m->t[3 * (size_t)f + 1]), RP_X(m, m->t[3 * (size_t)f + 2]), n);
		area[r] += 0.5 * sqrt(mrp_dot(n, n));
		cnt[r]++;
	}
	for (v = 0; v < m->nv; v++)
		if (cnt[v] && (area[v] < min_area || cnt[v] < min_faces))
		{
			w->r->removed_components++;
			w->r->removed_faces += cnt[v];
		}
	for (f = 0; f < m->nt; f++)
	{
		unsigned int r = rp_find(uf, m->t[3 * (size_t)f]);
		del[f] = area[r] < min_area || cnt[r] < min_faces;
	}
	rp_compact(m, del, m->nt);
done:
	qaws_internal_free(uf);
	qaws_internal_free(cnt);
	qaws_internal_free(area);
	qaws_internal_free(del);
}

/* ---- orientation -------------------------------------------------------- */

/* the parts across manifold edges (shell per face), each made consistent by a breadth-first walk turning faces; bad per part:
   no consistent orientation exists (a conflict met). Returns the parts' count. */
static unsigned int rp_consistent(rp_work* w, unsigned int* shell, unsigned char** out_bad)
{
	rp_mesh* m = &w->m;
	unsigned int ne, r, k, f, ns = 0;
	rp_he* E = rp_edges(m, 0, m->nt, &ne);
	unsigned int* pos = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (ne ? ne : 1)));
	unsigned int* mate = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (ne ? ne : 1)));
	unsigned int* queue = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nt ? m->nt : 1)));
	unsigned char* o = (unsigned char*)qaws_internal_scratch_zalloc(m->nt ? m->nt : 1, 1);
	unsigned char* bad = (unsigned char*)qaws_internal_zalloc(m->nt ? m->nt : 1, 1);
	*out_bad = NULL;
	if (!E || !pos || !mate || !queue || !o || !bad)
	{
		qaws_internal_free(bad);
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (r = 0; r < ne; r = k)
	{
		unsigned int i;
		k = rp_group_end(E, ne, r);
		for (i = r; i < k; i++)
		{
			pos[3 * E[i].f + E[i].k] = i;
			mate[i] = k - r == 2 ? (i == r ? r + 1 : r) : HM_NONE;
		}
	}
	for (f = 0; f < m->nt; f++)
		shell[f] = HM_NONE;
	for (f = 0; f < m->nt; f++)
	{
		unsigned int qn = 1, qi;
		if (shell[f] != HM_NONE)
			continue;
		shell[f] = ns;
		queue[0] = f;
		for (qi = 0; qi < qn; qi++)
		{
			unsigned int g = queue[qi];
			for (k = 0; k < 3; k++)
			{
				unsigned int rg = pos[3 * g + k], rh = mate[rg], h;
				int dg, dh;
				if (rh == HM_NONE)
					continue;
				h = E[rh].f;
				dg = rp_dir(m, E + rg) ^ o[g];
				dh = rp_dir(m, E + rh);
				if (shell[h] == HM_NONE)
				{
					shell[h] = ns;
					o[h] = (unsigned char)(dg == dh);
					queue[qn++] = h;
				}
				else if (dg == (dh ^ o[h]))
					bad[ns] = 1;
			}
		}
		ns++;
	}
	for (f = 0; f < m->nt; f++)
		if (o[f])
			rp_flip(m, f);
	*out_bad = bad;
done:
	qaws_internal_free(E);
	qaws_internal_free(pos);
	qaws_internal_free(mate);
	qaws_internal_free(queue);
	qaws_internal_free(o);
	return ns;
}

/* the solid angle of triangle (a, b, c) seen from the origin (Van Oosterom and Strackee) */
static double rp_solid_angle(double const* a, double const* b, double const* c)
{
	double la = sqrt(mrp_dot(a, a)), lb = sqrt(mrp_dot(b, b)), lc = sqrt(mrp_dot(c, c)), bc[3];
	rp_cross(b, c, bc);
	return 2 * atan2(mrp_dot(a, bc), la * lb * lc + mrp_dot(a, b) * lc + mrp_dot(b, c) * la + mrp_dot(c, a) * lb);
}

/* each part turned so its signed volume is positive, or negative at an odd depth inside other closed parts (removed instead
   with remove_nested); the parts with no consistent orientation counted */
static void rp_outward(rp_work* w, int remove_nested)
{
	rp_mesh* m = &w->m;
	unsigned int nt = m->nt, ns, s, f, r, k, ne = 0, i;
	unsigned int *shell = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt ? nt : 1))), *off = NULL, *list = NULL, *depth = NULL, *stamp = NULL;
	unsigned char *bad = NULL, *open = NULL, *del = NULL;
	double *vol = NULL, *box = NULL;
	rp_he* E = NULL;
	if (!shell)
		goto oom;
	ns = rp_consistent(w, shell, &bad);
	if (w->st != QAWS_STATUS_OK)
		goto done;
	E = rp_edges(m, 0, nt, &ne);
	off = (unsigned int*)qaws_internal_zalloc(ns + 1, sizeof(unsigned int));
	list = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt ? nt : 1)));
	depth = (unsigned int*)qaws_internal_zalloc(ns ? ns : 1, sizeof(unsigned int));
	stamp = (unsigned int*)qaws_internal_zalloc(m->nv ? m->nv : 1, sizeof(unsigned int));
	open = (unsigned char*)qaws_internal_scratch_zalloc(ns ? ns : 1, 1);
	del = (unsigned char*)qaws_internal_scratch_zalloc(nt ? nt : 1, 1);
	vol = (double*)qaws_internal_scratch_zalloc(ns ? ns : 1, sizeof(double));
	box = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 6, (ns ? ns : 1)));
	if (!E || !off || !list || !depth || !stamp || !open || !del || !vol || !box)
		goto oom;
	for (s = 0; s < ns; s++)
	{
		w->r->non_orientable += bad[s];
		box[6 * s] = box[6 * s + 1] = box[6 * s + 2] = HUGE_VAL;
		box[6 * s + 3] = box[6 * s + 4] = box[6 * s + 5] = -HUGE_VAL;
	}
	/* a part is closed when its own half-edges cancel on every edge */
	for (r = 0; r < ne; r = k)
	{
		unsigned int a, b;
		k = rp_group_end(E, ne, r);
		for (a = r; a < k; a++)
		{
			int net = 0;
			for (b = r; b < k; b++)
				if (shell[E[b].f] == shell[E[a].f])
					net += rp_dir(m, E + b) ? 1 : -1;
			if (net)
				open[shell[E[a].f]] = 1;
		}
	}
	for (f = 0; f < nt; f++)
		off[shell[f] + 1]++;
	for (s = 0; s < ns; s++)
		off[s + 1] += off[s];
	for (f = 0; f < nt; f++)
	{
		double* bx = box + 6 * shell[f];
		list[off[shell[f]] + depth[shell[f]]++] = f;
		for (k = 0; k < 3; k++)
			for (i = 0; i < 3; i++)
			{
				double c = RP_X(m, m->t[3 * (size_t)f + k])[i];
				bx[i] = fmin(bx[i], c);
				bx[3 + i] = fmax(bx[3 + i], c);
			}
	}
	/* signed volumes about each part's box centre (the same for any centre when closed) */
	for (f = 0; f < nt; f++)
	{
		double const* bx = box + 6 * shell[f];
		double c[3], a[3], b[3], d[3], n[3];
		for (i = 0; i < 3; i++)
			c[i] = 0.5 * (bx[i] + bx[3 + i]);
		mrp_sub(RP_X(m, m->t[3 * (size_t)f]), c, a);
		mrp_sub(RP_X(m, m->t[3 * (size_t)f + 1]), c, b);
		mrp_sub(RP_X(m, m->t[3 * (size_t)f + 2]), c, d);
		rp_cross(b, d, n);
		vol[shell[f]] += mrp_dot(a, n) / 6;
	}
	memset(depth, 0, sizeof(unsigned int) * (ns ? ns : 1));
	/* depth: the closed parts T round a part S, judged at up to 16 of S's vertices by T's winding number
	   shortcut: 16 samples per pair, every closed part tried against every part inside its box; upgrade to a tree of parts
	   if inputs with thousands of closed parts come up */
	for (s = 0; s < ns; s++)
	{
		double const* bt = box + 6 * s;
		unsigned int q;
		if (open[s] || bad[s])
			continue;
		for (i = off[s]; i < off[s + 1]; i++)
			for (k = 0; k < 3; k++)
				stamp[m->t[3 * (size_t)list[i] + k]] = s + 1;
		for (q = 0; q < ns; q++)
		{
			double const* bq = box + 6 * q;
			unsigned int n = off[q + 1] - off[q], j, tested = 0, inside = 0;
			if (q == s || bq[0] < bt[0] || bq[1] < bt[1] || bq[2] < bt[2] || bq[3] > bt[3] || bq[4] > bt[4] || bq[5] > bt[5])
				continue;
			for (j = 0; j < (n < 16 ? n : 16); j++)
			{
				unsigned int v = m->t[3 * (size_t)list[off[q] + (size_t)j * n / (n < 16 ? n : 16)]];
				double const* p = RP_X(m, v);
				double wn = 0;
				if (stamp[v] == s + 1)
					continue;
				for (i = off[s]; i < off[s + 1]; i++)
				{
					double a[3], b[3], d[3];
					mrp_sub(RP_X(m, m->t[3 * (size_t)list[i]]), p, a);
					mrp_sub(RP_X(m, m->t[3 * (size_t)list[i] + 1]), p, b);
					mrp_sub(RP_X(m, m->t[3 * (size_t)list[i] + 2]), p, d);
					wn += rp_solid_angle(a, b, d);
				}
				tested++;
				inside += fabs(wn) > 2 * RP_PI;
			}
			if (tested && inside == tested)
				depth[q]++;
		}
	}
	for (s = 0; s < ns; s++)
	{
		int want = depth[s] % 2 ? -1 : 1;
		if (depth[s])
			w->r->nested_shells++;
		for (i = off[s]; i < off[s + 1]; i++)
		{
			if (remove_nested && depth[s])
				del[list[i]] = 1;
			else if (vol[s] * want < 0)
				rp_flip(m, list[i]);
		}
	}
	rp_compact(m, del, nt);
	goto done;
oom:
	w->st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(shell); qaws_internal_free(off); qaws_internal_free(list); qaws_internal_free(depth); qaws_internal_free(stamp); qaws_internal_free(bad); qaws_internal_free(open); qaws_internal_free(del); qaws_internal_free(vol); qaws_internal_free(box); qaws_internal_free(E);
}

/* ---- holes -------------------------------------------------------------- */

typedef struct rp_fill_ctx
{
	rp_he const* E;             /* the half-edges before the fills */
	unsigned int ne;
	double perimeter;
	double* sigma;              /* per vertex: the length its edges should have */
	unsigned int scap;
	unsigned int* O;            /* scratch: per loop edge the third vertex of its face */
	unsigned char* lock;        /* per vertex (as sigma): on a crease path across a fill, its edges never flipped, never faired */
	unsigned int* crease;       /* per vertex before the fills: the far end of its one crease edge (HM_NONE: none or several) */
	unsigned char* crease_used; /* per vertex before the fills: its crease continued already */
	unsigned int ncrease;
} rp_fill_ctx;

static int rp_sigma_room(rp_fill_ctx* c, unsigned int nv)
{
	if (nv > c->scap)
	{
		unsigned int cap = qaws_internal_grow(c->scap, nv, 16);
		double* s = cap ? (double*)qaws_internal_realloc(c->sigma, qaws_internal_mul2(sizeof(double), cap)) : NULL;
		unsigned char* l;
		if (!s)
			return 0;
		c->sigma = s;
		l = (unsigned char*)qaws_internal_realloc(c->lock, cap);
		if (!l)
			return 0;
		memset(l + c->scap, 0, cap - c->scap);
		c->lock = l;
		c->scap = cap;
	}
	return 1;
}

/* the cosine of the angle between unit normals; -1 (the worst) when one is degenerate */
static double rp_cos(double const* a, double const* b)
{
	double d;
	if (b[0] == HUGE_VAL)
		return 1;
	d = mrp_dot(a, b);
	return mrp_dot(a, a) > 0 && mrp_dot(b, b) > 0 ? d : -1;
}

/* loop L (its edges L[i] -> L[i + 1] those of faces whose third vertices are O[i]) filled by Liepa's programme: the
   triangulation largest in its smallest dihedral cosine (with the faces round the loop and between its own triangles), then
   least in area; turned against the loop. A diagonal that is an edge already is taken only when the loop leaves no other. */
static void rp_fill_dp(rp_work* w, rp_fill_ctx const* c, unsigned int const* L, unsigned int const* O, unsigned int n)
{
	rp_mesh* m = &w->m;
	size_t nn = (size_t)n * n;
	double* A = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), nn));          /* the triangulation of [i, k]: its smallest cosine */
	double* S = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), nn));          /* its area */
	double* N = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, nn));      /* its triangle on edge (i, k): the unit normal */
	double* Ne = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, n));      /* the faces round the loop */
	unsigned int* K = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), nn));
	unsigned int* stack = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, n));
	unsigned int i, len, top = 0;
	if (!A || !S || !N || !Ne || !K || !stack)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < n; i++)
	{
		/* an edge with no face beside it (a crease path, a halving diagonal): its neighbour judged level with the new face */
		if (O[i] == HM_NONE)
			Ne[3 * i] = Ne[3 * i + 1] = Ne[3 * i + 2] = HUGE_VAL;
		else
		{
			rp_normal(RP_X(m, L[i]), RP_X(m, L[(i + 1) % n]), RP_X(m, O[i]), Ne + 3 * i);
			rp_unit(Ne + 3 * i);
		}
		if (i + 1 < n)
		{
			A[i * (size_t)n + i + 1] = 1;
			S[i * (size_t)n + i + 1] = 0;
		}
	}
	for (len = 2; len < n; len++)
		for (i = 0; i + len < n; i++)
		{
			unsigned int k = i + len, j, best = i + 1;
			double ba = -HUGE_VAL, bs = HUGE_VAL, bn[3] = { 0, 0, 0 };
			for (j = i + 1; j < k; j++)
			{
				double nr[3], a, s;
				rp_normal(RP_X(m, L[k]), RP_X(m, L[j]), RP_X(m, L[i]), nr);
				s = S[i * (size_t)n + j] + S[j * (size_t)n + k] + 0.5 * rp_unit(nr);
				a = fmin(A[i * (size_t)n + j], A[j * (size_t)n + k]);
				a = fmin(a, rp_cos(nr, j == i + 1 ? Ne + 3 * i : N + 3 * (i * (size_t)n + j)));
				a = fmin(a, rp_cos(nr, k == j + 1 ? Ne + 3 * j : N + 3 * (j * (size_t)n + k)));
				if (i == 0 && k == n - 1)
					a = fmin(a, rp_cos(nr, Ne + 3 * (n - 1)));
				if (a > ba + 1e-12 || (a >= ba - 1e-12 && s < bs))
				{
					ba = a;
					bs = s;
					best = j;
					memcpy(bn, nr, sizeof bn);
				}
			}
			if (len < n - 1 && rp_has_edge(c->E, c->ne, L[i], L[k]))
				ba -= 4;
			A[i * (size_t)n + k] = ba;
			S[i * (size_t)n + k] = bs;
			K[i * (size_t)n + k] = best;
			memcpy(N + 3 * (i * (size_t)n + k), bn, sizeof bn);
		}
	stack[top++] = 0;
	stack[top++] = n - 1;
	while (top && w->st == QAWS_STATUS_OK)
	{
		unsigned int k = stack[--top], j;
		i = stack[--top];
		if (k - i < 2)
			continue;
		j = K[i * (size_t)n + k];
		if (!rp_add_t(m, L[k], L[j], L[i], -1, 0))
			w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		stack[top++] = i;
		stack[top++] = j;
		stack[top++] = j;
		stack[top++] = k;
	}
done:
	qaws_internal_free(A);
	qaws_internal_free(S);
	qaws_internal_free(N);
	qaws_internal_free(Ne);
	qaws_internal_free(K);
	qaws_internal_free(stack);
}

/* loop L (third vertices O) filled: by the programme up to RP_DP_MAX vertices, else halved by its shortest diagonal between
   vertices a third to a half of the loop apart (sampled), the diagonal an edge with no face beside it in each half,
   recursively */
static void rp_fill_long(rp_work* w, rp_fill_ctx const* c, unsigned int const* L, unsigned int const* O, unsigned int n)
{
	rp_mesh* m = &w->m;
	unsigned int step = n / 256 + 1, i, j, bi = 0, bj = n / 2, na, nb, k;
	double best = HUGE_VAL;
	unsigned int *A, *B;
	if (n <= RP_DP_MAX)
	{
		rp_fill_dp(w, c, L, O, n);
		return;
	}
	for (i = 0; i < n; i += step)
		for (j = i + n / 3; j <= i + n / 2; j += step)
		{
			double d = rp_dist(RP_X(m, L[i]), RP_X(m, L[j % n]));
			if (d < best)
			{
				best = d;
				bi = i;
				bj = j % n;
			}
		}
	if (bj < bi)
	{
		unsigned int t = bi;
		bi = bj;
		bj = t;
	}
	na = bj - bi + 1;
	nb = n - na + 2;
	A = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)na));
	B = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nb));
	if (!A || !B)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		qaws_internal_free(A);
		qaws_internal_free(B);
		return;
	}
	/* A: L[bi..bj], closed by bj -> bi; B: L[bj..bi + n], closed by bi -> bj */
	for (k = 0; k < na; k++)
	{
		A[k] = L[bi + k];
		A[na + k] = O[bi + k];
	}
	A[na + na - 1] = HM_NONE;
	for (k = 0; k < nb; k++)
	{
		B[k] = L[(bj + k) % n];
		B[nb + k] = O[(bj + k) % n];
	}
	B[nb + nb - 1] = HM_NONE;
	rp_fill_long(w, c, A, A + na, na);
	if (w->st == QAWS_STATUS_OK)
		rp_fill_long(w, c, B, B + nb, nb);
	qaws_internal_free(A);
	qaws_internal_free(B);
}

/* a new vertex on a crease path: its length wanted, locked */
static unsigned int rp_path_vertex(rp_work* w, rp_fill_ctx* c, double const* p, double sigma)
{
	unsigned int v = rp_add_v(&w->m, p);
	if (v == HM_NONE || !rp_sigma_room(c, w->m.nv))
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return HM_NONE;
	}
	c->sigma[v] = sigma;
	c->lock[v] = 1;
	return v;
}

/* the inner vertices of a path from vertex a to point q (b: its vertex, or HM_NONE for a point not made yet): a cubic Hermite
   leaving a along da and arriving along -db (straight when da or db is NULL), one vertex per wanted length; into P, the count
   returned (at most cap) */
static unsigned int rp_path(rp_work* w, rp_fill_ctx* c, unsigned int a, double const* q, double sq, double const* da, double const* db, unsigned int* P,
	unsigned int cap)
{
	double p[3], len, sg, x[3];
	unsigned int s, i, k;
	memcpy(p, RP_X(&w->m, a), sizeof p);
	len = rp_dist(p, q);
	sg = 0.5 * (c->sigma[a] + sq);
	s = sg > 0 ? (unsigned int)floor(len / sg + 0.5) : 1;
	s = s < 1 ? 1 : s > cap + 1 ? cap + 1 : s;
	for (i = 1; i < s; i++)
	{
		double t = (double)i / s, h00 = 2 * t * t * t - 3 * t * t + 1, h10 = t * t * t - 2 * t * t + t, h01 = -2 * t * t * t + 3 * t * t, h11 = t * t * t - t * t;
		for (k = 0; k < 3; k++)
			x[k] = da && db ? h00 * p[k] + h10 * len * da[k] + h01 * q[k] - h11 * len * db[k] : p[k] + t * (q[k] - p[k]);
		P[i - 1] = rp_path_vertex(w, c, x, sg);
		if (P[i - 1] == HM_NONE)
			return 0;
	}
	return s - 1;
}

/* a sub-loop: L[i0..i1] along the loop (indices mod n) with its third vertices, then the k vertices of Q, its edges after the
   loop's free */
static void rp_sub_loop(unsigned int const* L, unsigned int const* O, unsigned int n, unsigned int i0, unsigned int i1, unsigned int const* Q,
	unsigned int k, unsigned int* S, unsigned int* SO, unsigned int* out_n)
{
	unsigned int j = 0, i;
	for (i = i0;; i = (i + 1) % n)
	{
		S[j] = L[i];
		SO[j++] = i == i1 ? HM_NONE : O[i];
		if (i == i1)
			break;
	}
	for (i = 0; i < k; i++)
	{
		S[j] = Q[i];
		SO[j++] = HM_NONE;
	}
	*out_n = j;
}

static void rp_fill_feature(rp_work* w, rp_fill_ctx* c, unsigned int const* L, unsigned int const* O, unsigned int n, unsigned int depth);

/* the loop split along crease paths inside it: two crease ends facing each other joined, or three and more meeting at a point
   (a corner), and each part filled on its own; otherwise filled whole */
static void rp_fill_feature(rp_work* w, rp_fill_ctx* c, unsigned int const* L, unsigned int const* O, unsigned int n, unsigned int depth)
{
	rp_mesh* m = &w->m;
	unsigned int ent[16], ne = 0, i, j, k, ba = 0, bb = 0, *S = NULL, *SO = NULL, *P = NULL;
	double dir[16][3], best = HUGE_VAL;
	for (i = 0; i < n && ne < 16 && depth < 8; i++)
		if (L[i] < c->ncrease && c->crease[L[i]] != HM_NONE && !c->crease_used[L[i]])
		{
			mrp_sub(RP_X(m, L[i]), RP_X(m, c->crease[L[i]]), dir[ne]);
			if (rp_unit(dir[ne]) > 0)
				ent[ne++] = i;
		}
	if (ne < 2)
	{
		rp_fill_long(w, c, L, O, n);
		return;
	}
	/* the pair facing each other best (each end's crease within 37 degrees of the chord), at least two loop steps apart both ways */
	for (i = 0; i < ne; i++)
		for (j = i + 1; j < ne; j++)
		{
			double ch[3], l, ci, cj;
			unsigned int gap = ent[j] - ent[i];
			if (gap < 2 || n - gap < 2)
				continue;
			mrp_sub(RP_X(m, L[ent[j]]), RP_X(m, L[ent[i]]), ch);
			l = rp_unit(ch);
			if (!(l > 0))
				continue;
			ci = mrp_dot(dir[i], ch);
			cj = -mrp_dot(dir[j], ch);
			if (ci > 0.8 && cj > 0.8 && 2 - ci - cj < best)
			{
				best = 2 - ci - cj;
				ba = i;
				bb = j;
			}
		}
	S = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)n + 4096));
	SO = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)n + 4096));
	P = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), 4096));
	if (!S || !SO || !P)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if (best < HUGE_VAL)
	{
		unsigned int a = L[ent[ba]], b = L[ent[bb]], np, ns;
		double q[3];
		memcpy(q, RP_X(m, b), sizeof q);
		c->crease_used[a] = c->crease_used[b] = 1;
		c->lock[a] = c->lock[b] = 1;
		np = rp_path(w, c, a, q, c->sigma[b], dir[ba], dir[bb], P, 2048);
		if (w->st != QAWS_STATUS_OK)
			goto done;
		/* the side from a to b along the loop, back over the path reversed; then the other side, forward over the path */
		for (k = 0; k < np; k++)
			P[2048 + k] = P[np - 1 - k];
		rp_sub_loop(L, O, n, ent[ba], ent[bb], P + 2048, np, S, SO, &ns);
		rp_fill_feature(w, c, S, SO, ns, depth + 1);
		if (w->st == QAWS_STATUS_OK)
		{
			rp_sub_loop(L, O, n, ent[bb], ent[ba], P, np, S, SO, &ns);
			rp_fill_feature(w, c, S, SO, ns, depth + 1);
		}
		goto done;
	}
	if (ne >= 3)
	{
		/* a corner: the point nearest every crease line, ahead of each end and within a quarter of its distance of each line */
		double A[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, rhs[3] = { 0, 0, 0 }, C[3], det, sg = 0;
		int ok = 1;
		unsigned int cv, np[16], off[16], tot = 0, ns;
		for (i = 0; i < ne; i++)
		{
			double const* p = RP_X(m, L[ent[i]]);
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
				{
					double e = (j == k) - dir[i][j] * dir[i][k];
					A[3 * j + k] += e;
					rhs[j] += e * p[k];
				}
			sg += c->sigma[L[ent[i]]] / ne;
		}
		det = A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) + A[2] * (A[3] * A[7] - A[4] * A[6]);
		if (fabs(det) < 1e-12)
			ok = 0;
		else
			for (j = 0; j < 3; j++)
			{
				double M[9];
				memcpy(M, A, sizeof M);
				M[j] = rhs[0];
				M[3 + j] = rhs[1];
				M[6 + j] = rhs[2];
				C[j] = (M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) + M[2] * (M[3] * M[7] - M[4] * M[6])) / det;
			}
		for (i = 0; i < ne && ok; i++)
		{
			double v[3], t, off2[3];
			mrp_sub(C, RP_X(m, L[ent[i]]), v);
			t = mrp_dot(v, dir[i]);
			off2[0] = v[0] - t * dir[i][0];
			off2[1] = v[1] - t * dir[i][1];
			off2[2] = v[2] - t * dir[i][2];
			ok = t > 0 && sqrt(mrp_dot(off2, off2)) < 0.25 * t;
		}
		if (ok)
		{
			cv = rp_path_vertex(w, c, C, sg);
			if (cv == HM_NONE)
				goto done;
			for (i = 0; i < ne; i++)
			{
				unsigned int e = L[ent[i]];
				c->crease_used[e] = 1;
				c->lock[e] = 1;
				off[i] = tot;
				np[i] = rp_path(w, c, e, C, sg, NULL, NULL, P + tot, (2048 - tot) / 2);
				if (w->st != QAWS_STATUS_OK)
					goto done;
				tot += np[i];
			}
			/* each sector: the loop from one end to the next, in over the next end's path, the corner, out over this end's path */
			for (i = 0; i < ne && w->st == QAWS_STATUS_OK; i++)
			{
				unsigned int nx = (i + 1) % ne, q = 0;
				for (k = 0; k < np[nx]; k++)
					P[2048 + q++] = P[off[nx] + k];
				P[2048 + q++] = cv;
				for (k = np[i]; k-- > 0;)
					P[2048 + q++] = P[off[i] + k];
				rp_sub_loop(L, O, n, ent[i], ent[nx], P + 2048, q, S, SO, &ns);
				rp_fill_feature(w, c, S, SO, ns, depth + 1);
			}
			goto done;
		}
	}
	rp_fill_long(w, c, L, O, n);
done:
	qaws_internal_free(S);
	qaws_internal_free(SO);
	qaws_internal_free(P);
}

/* one boundary loop (pv: its vertices, pr: its boundary records' indices into be) filled, or counted as left */
static void rp_fill_loop(rp_work* w, rp_fill_ctx* c, unsigned int const* be, unsigned int const* pv, unsigned int const* pr, unsigned int n)
{
	rp_mesh* m = &w->m;
	double perimeter = 0;
	unsigned int i, skip = n < 3;
	for (i = 0; i < n; i++)
	{
		rp_he const* e = c->E + be[pr[i]];
		c->O[i] = m->t[3 * (size_t)e->f + (e->k + 2) % 3];
		perimeter += rp_dist(RP_X(m, pv[i]), RP_X(m, pv[(i + 1) % n]));
	}
	if (skip || perimeter > c->perimeter)
	{
		w->r->holes_left++;
		return;
	}
	rp_fill_feature(w, c, pv, c->O, n, 0);
	w->r->holes_filled++;
}

/* the fills' inner edges (triangles from t0) flipped toward Delaunay (opposite angles summing past a half turn), no edge
   doubled, no face turned over */
static void rp_relax(rp_work* w, rp_fill_ctx const* c, unsigned int t0)
{
	rp_mesh* m = &w->m;
	unsigned int sweep;
	for (sweep = 0; sweep < 50 && w->st == QAWS_STATUS_OK; sweep++)
	{
		unsigned int nf, r, flips = 0;
		rp_he* F = rp_edges(m, t0, m->nt, &nf);
		unsigned char* touched = (unsigned char*)qaws_internal_scratch_zalloc(m->nt - t0 + 1, 1);
		if (!F || !touched)
		{
			qaws_internal_free(F);
			qaws_internal_free(touched);
			w->st = QAWS_STATUS_ALLOCATION_FAILURE;
			return;
		}
		for (r = 0; r + 1 < nf; r++)
		{
			unsigned int f = F[r].f, g = F[r + 1].f, a, b, p, q;
			double u[3], v[3], n0[3], n1[3], n2[3], nf2[3], ap, aq;
			if (F[r].key != F[r + 1].key || (r + 2 < nf && F[r + 2].key == F[r].key) || touched[f - t0] || touched[g - t0])
				continue;
			a = m->t[3 * (size_t)f + F[r].k];
			b = m->t[3 * (size_t)f + (F[r].k + 1) % 3];
			p = m->t[3 * (size_t)f + (F[r].k + 2) % 3];
			if (m->t[3 * (size_t)g + F[r + 1].k] != b)
				continue;
			q = m->t[3 * (size_t)g + (F[r + 1].k + 2) % 3];
			if (p == q || rp_has_edge(F, nf, p, q) || rp_has_edge(c->E, c->ne, p, q) || (c->lock && c->lock[a] && c->lock[b]))
				continue;
			mrp_sub(RP_X(m, a), RP_X(m, p), u);
			mrp_sub(RP_X(m, b), RP_X(m, p), v);
			ap = rp_angle(u, v);
			mrp_sub(RP_X(m, a), RP_X(m, q), u);
			mrp_sub(RP_X(m, b), RP_X(m, q), v);
			aq = rp_angle(u, v);
			if (ap + aq <= RP_PI + 1e-9)
				continue;
			rp_normal(RP_X(m, a), RP_X(m, b), RP_X(m, p), n0);
			rp_normal(RP_X(m, b), RP_X(m, a), RP_X(m, q), nf2);
			n0[0] += nf2[0];
			n0[1] += nf2[1];
			n0[2] += nf2[2];
			rp_normal(RP_X(m, p), RP_X(m, a), RP_X(m, q), n1);
			rp_normal(RP_X(m, q), RP_X(m, b), RP_X(m, p), n2);
			if (mrp_dot(n1, n0) <= 0 || mrp_dot(n2, n0) <= 0)
				continue;
			m->t[3 * (size_t)f] = p;
			m->t[3 * (size_t)f + 1] = a;
			m->t[3 * (size_t)f + 2] = q;
			m->t[3 * (size_t)g] = q;
			m->t[3 * (size_t)g + 1] = b;
			m->t[3 * (size_t)g + 2] = p;
			touched[f - t0] = touched[g - t0] = 1;
			flips++;
		}
		qaws_internal_free(F);
		qaws_internal_free(touched);
		if (!flips)
			break;
	}
}

/* Liepa's refinement of the fills (triangles from t0): a triangle split at its centroid c when sqrt(2) |c - x| exceeds both
   the length wanted at c (its corners' mean) and at each corner x; flips between rounds */
static void rp_refine(rp_work* w, rp_fill_ctx* c, unsigned int t0)
{
	rp_mesh* m = &w->m;
	unsigned int round;
	for (round = 0; round < 32 && w->st == QAWS_STATUS_OK; round++)
	{
		unsigned int nt = m->nt, f, splits = 0;
		for (f = t0; f < nt && w->st == QAWS_STATUS_OK; f++)
		{
			unsigned int a = m->t[3 * (size_t)f], b = m->t[3 * (size_t)f + 1], d = m->t[3 * (size_t)f + 2], k, v;
			unsigned int abd[3];
			double C[3], sc = (c->sigma[a] + c->sigma[b] + c->sigma[d]) / 3;
			int ok = 1;
			abd[0] = a;
			abd[1] = b;
			abd[2] = d;
			for (k = 0; k < 3; k++)
				C[k] = (RP_X(m, a)[k] + RP_X(m, b)[k] + RP_X(m, d)[k]) / 3;
			for (k = 0; k < 3; k++)
			{
				double l = sqrt(2.0) * rp_dist(C, RP_X(m, abd[k]));
				if (l <= sc || l <= c->sigma[abd[k]])
					ok = 0;
			}
			if (!ok)
				continue;
			v = rp_add_v(m, C);
			if (v == HM_NONE || !rp_sigma_room(c, m->nv) || !rp_add_t(m, b, d, v, -1, 0) || !rp_add_t(m, d, a, v, -1, 0))
			{
				w->st = QAWS_STATUS_ALLOCATION_FAILURE;
				return;
			}
			c->sigma[v] = sc;
			m->t[3 * (size_t)f + 2] = v;
			splits++;
		}
		if (!splits)
			break;
		rp_relax(w, c, t0);
	}
}

typedef struct rp_trip
{
	unsigned int i, j;
	double v;
} rp_trip;

static int rp_trip_cmp(void const* p, void const* q)
{
	rp_trip const* a = (rp_trip const*)p;
	rp_trip const* b = (rp_trip const*)q;
	if (a->i != b->i)
		return a->i < b->i ? -1 : 1;
	return a->j < b->j ? -1 : a->j > b->j;
}

/* the fills' new vertices (from v0) faired: the Laplacian (Liepa's scale-dependent umbrella, neighbours weighted by
   1 / length) at them and at the vertices next to them made least in the sum of squares (the bi-Laplacian's normal
   equations, one LDL^T), every other vertex held */
#define RP_FREE(v) ((v) >= v0 && fi[(v) - v0] != HM_NONE)
static void rp_fair(rp_work* w, unsigned int v0, unsigned char const* lock)
{
	rp_mesh* m = &w->m;
	unsigned int nv = m->nv, nfree = 0, ne, r, k, v, i, j, maxd = 0, nt = 0, ct = 0;
	unsigned int* fi = NULL;   /* per vertex from v0: its unknown, or HM_NONE when locked */
	rp_he* E = NULL;
	unsigned int *deg = NULL, *off = NULL, *adj = NULL, *rp = NULL, *col = NULL, *fc = NULL;
	double *val = NULL, *b = NULL, *xyz = NULL, *sol = NULL, *fw = NULL;
	rp_trip* T = NULL;
	qaws_internal_ldlt L;
	qaws_status st;
	memset(&L, 0, sizeof L);
	if (nv <= v0)
		return;
	fi = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), nv - v0));
	if (!fi)
		goto oom;
	for (v = v0; v < nv; v++)
		fi[v - v0] = lock && lock[v] ? HM_NONE : nfree++;
	if (!nfree)
		goto done;
	E = rp_edges(m, 0, m->nt, &ne);
	deg = (unsigned int*)qaws_internal_scratch_zalloc(nv, sizeof(unsigned int));
	off = (unsigned int*)qaws_internal_scratch_zalloc(nv + 1, sizeof(unsigned int));
	if (!E || !deg || !off)
		goto oom;
	for (r = 0; r < ne; r = rp_group_end(E, ne, r))
	{
		deg[(unsigned int)(E[r].key >> 32)]++;
		deg[(unsigned int)E[r].key]++;
	}
	for (v = 0; v < nv; v++)
	{
		off[v + 1] = off[v] + deg[v];
		maxd = deg[v] > maxd ? deg[v] : maxd;
		deg[v] = 0;
	}
	adj = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (off[nv] ? off[nv] : 1)));
	fc = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (maxd + 1)));
	fw = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (maxd + 1)));
	b = (double*)qaws_internal_scratch_zalloc(qaws_internal_mul2(3, (size_t)nfree), sizeof(double));
	if (!adj || !fc || !fw || !b)
		goto oom;
	for (r = 0; r < ne; r = rp_group_end(E, ne, r))
	{
		unsigned int a = (unsigned int)(E[r].key >> 32), c = (unsigned int)E[r].key;
		adj[off[a] + deg[a]++] = c;
		adj[off[c] + deg[c]++] = a;
	}
	/* rows: the Laplacian at each free vertex and at each vertex next to one; its free entries go into the normal matrix, its
	   held ones into the right-hand side */
	for (v = 0; v < nv; v++)
	{
		unsigned int nfc = 0;
		double rh[3] = { 0, 0, 0 }, inv;
		int row = RP_FREE(v);
		for (i = off[v]; i < off[v + 1] && !row; i++)
			row = RP_FREE(adj[i]);
		/* (no row at a crease path's vertex: its neighbours lie on both sides) */
		if (!row || !deg[v] || (lock && lock[v]))
			continue;
		for (i = off[v], inv = 0; i < off[v + 1]; i++)
			inv += 1 / fmax(rp_dist(RP_X(m, v), RP_X(m, adj[i])), 1e-300);
		if (RP_FREE(v))
		{
			fc[nfc] = fi[v - v0];
			fw[nfc++] = -1;
		}
		else
			for (k = 0; k < 3; k++)
				rh[k] -= RP_X(m, v)[k];
		for (i = off[v]; i < off[v + 1]; i++)
			if (RP_FREE(adj[i]))
			{
				fc[nfc] = fi[adj[i] - v0];
				fw[nfc++] = 1 / (fmax(rp_dist(RP_X(m, v), RP_X(m, adj[i])), 1e-300) * inv);
			}
			else
				for (k = 0; k < 3; k++)
					rh[k] += RP_X(m, adj[i])[k] / (fmax(rp_dist(RP_X(m, v), RP_X(m, adj[i])), 1e-300) * inv);
		if (nt + (size_t)nfc * nfc > ct)
		{
			unsigned int cap = qaws_internal_grow(ct, nt + (size_t)nfc * nfc, 16);
			rp_trip* g = cap ? (rp_trip*)qaws_internal_realloc(T, qaws_internal_mul2(sizeof(rp_trip), cap)) : NULL;
			if (!g)
				goto oom;
			T = g;
			ct = cap;
		}
		for (i = 0; i < nfc; i++)
		{
			for (k = 0; k < 3; k++)
				b[3 * (size_t)fc[i] + k] -= fw[i] * rh[k];
			for (j = 0; j < nfc; j++)
			{
				T[nt].i = fc[i];
				T[nt].j = fc[j];
				T[nt++].v = fw[i] * fw[j];
			}
		}
	}
	qsort(T, nt, sizeof(rp_trip), rp_trip_cmp);
	rp = (unsigned int*)qaws_internal_scratch_zalloc(nfree + 1, sizeof(unsigned int));
	col = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt ? nt : 1)));
	val = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (nt ? nt : 1)));
	xyz = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nfree));
	sol = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)nfree));
	if (!rp || !col || !val || !xyz || !sol)
		goto oom;
	for (i = 0, j = 0; i < nt; i++)
	{
		if (j && T[i].i == T[i - 1].i && T[i].j == T[i - 1].j)
		{
			val[j - 1] += T[i].v;
			continue;
		}
		col[j] = T[i].j;
		val[j++] = T[i].v;
		rp[T[i].i + 1]++;
	}
	for (i = 0; i < nfree; i++)
		rp[i + 1] += rp[i];
	for (v = v0; v < nv; v++)
		if (RP_FREE(v))
			memcpy(xyz + 3 * (size_t)fi[v - v0], RP_X(m, v), sizeof(double) * 3);
	st = qaws_internal_ldlt_factor(&L, rp, col, val, xyz, nfree, 0);
	if (st == QAWS_STATUS_ALLOCATION_FAILURE)
		goto oom;
	if (st == QAWS_STATUS_OK)
		for (k = 0; k < 3; k++)
		{
			for (i = 0; i < nfree; i++)
				sol[i] = b[3 * (size_t)i + k];
			qaws_internal_ldlt_solve(sol, &L, sol);
			for (v = v0; v < nv; v++)
				if (RP_FREE(v))
					RP_X(m, v)[k] = sol[fi[v - v0]];
		}
	goto done;
oom:
	w->st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_ldlt_free(&L);
	qaws_internal_free(E); qaws_internal_free(deg); qaws_internal_free(off); qaws_internal_free(adj); qaws_internal_free(rp); qaws_internal_free(col); qaws_internal_free(fc); qaws_internal_free(val); qaws_internal_free(b); qaws_internal_free(xyz); qaws_internal_free(sol); qaws_internal_free(fw); qaws_internal_free(T); qaws_internal_free(fi);
}
#undef RP_FREE

/* boundary loops traced over the half-edges alone on their edges (a loop meeting itself at a vertex cut there) and filled */
static void rp_fill(rp_work* w, double perimeter, unsigned int flags, unsigned int const* shell, unsigned char const* bad)
{
	rp_mesh* m = &w->m;
	unsigned int nt0 = m->nt, nv0 = m->nv, ne, nb = 0, r, k, i, v;
	rp_fill_ctx c;
	rp_he* E = rp_edges(m, 0, nt0, &ne);
	unsigned int *be = NULL, *start = NULL, *cur = NULL, *ord = NULL, *inpath = NULL, *pv = NULL, *pr = NULL, *cnt = NULL;
	unsigned char* used = NULL;
	memset(&c, 0, sizeof c);
	c.sigma = (double*)qaws_internal_zalloc(nv0 ? nv0 : 1, sizeof(double));
	c.lock = (unsigned char*)qaws_internal_zalloc(nv0 ? nv0 : 1, 1);
	c.crease = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nv0 ? nv0 : 1)));
	c.crease_used = (unsigned char*)qaws_internal_scratch_zalloc(nv0 ? nv0 : 1, 1);
	c.ncrease = nv0;
	c.scap = nv0 ? nv0 : 1;
	cnt = (unsigned int*)qaws_internal_scratch_zalloc(nv0 ? nv0 : 1, sizeof(unsigned int));
	if (!E || !c.sigma || !cnt || !c.lock || !c.crease || !c.crease_used)
		goto oom;
	for (v = 0; v < nv0; v++)
		c.crease[v] = HM_NONE;
	/* creases: edges of two faces whose normals part by more than 40 degrees; a vertex on one keeps its far end, on several none */
	for (r = 0; r < ne; r = k)
	{
		k = rp_group_end(E, ne, r);
		if (k - r == 2)
		{
			unsigned int a = (unsigned int)(E[r].key >> 32), b = (unsigned int)E[r].key, f0 = E[r].f, f1 = E[r + 1].f, e;
			double n0[3], n1[3];
			rp_normal(RP_X(m, m->t[3 * (size_t)f0]), RP_X(m, m->t[3 * (size_t)f0 + 1]), RP_X(m, m->t[3 * (size_t)f0 + 2]), n0);
			rp_normal(RP_X(m, m->t[3 * (size_t)f1]), RP_X(m, m->t[3 * (size_t)f1 + 1]), RP_X(m, m->t[3 * (size_t)f1 + 2]), n1);
			if (rp_unit(n0) > 0 && rp_unit(n1) > 0 && mrp_dot(n0, n1) < 0.766)
				for (e = 0; e < 2; e++)
				{
					unsigned int p = e ? b : a, q = e ? a : b;
					if (!c.crease_used[p])
					{
						c.crease[p] = q;
						c.crease_used[p] = 1;
					}
					else
						c.crease[p] = HM_NONE;
				}
		}
	}
	memset(c.crease_used, 0, nv0 ? nv0 : 1);
	c.E = E;
	c.ne = ne;
	c.perimeter = perimeter;
	/* the lengths wanted: each vertex's mean edge */
	for (r = 0; r < ne; r = k)
	{
		unsigned int a = (unsigned int)(E[r].key >> 32), b = (unsigned int)E[r].key;
		double l = rp_dist(RP_X(m, a), RP_X(m, b));
		k = rp_group_end(E, ne, r);
		c.sigma[a] += l;
		c.sigma[b] += l;
		cnt[a]++;
		cnt[b]++;
		nb += k - r == 1;
	}
	for (v = 0; v < nv0; v++)
		c.sigma[v] = cnt[v] ? c.sigma[v] / cnt[v] : 0;
	be = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nb ? nb : 1)));
	start = (unsigned int*)qaws_internal_scratch_zalloc(nv0 + 1, sizeof(unsigned int));
	cur = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nv0 + 1)));
	ord = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nb ? nb : 1)));
	inpath = (unsigned int*)qaws_internal_scratch_zalloc(nv0 ? nv0 : 1, sizeof(unsigned int));
	pv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nb + 1)));
	pr = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nb + 1)));
	c.O = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nb + 1)));
	used = (unsigned char*)qaws_internal_zalloc(nb ? nb : 1, 1);
	if (!be || !start || !cur || !ord || !inpath || !pv || !pr || !c.O || !used)
		goto oom;
	for (r = 0, nb = 0; r < ne; r = k)
	{
		k = rp_group_end(E, ne, r);
		if (k - r == 1)
		{
			be[nb++] = r;
			start[m->t[3 * (size_t)E[r].f + E[r].k] + 1]++;
		}
	}
	for (v = 0; v < nv0; v++)
	{
		start[v + 1] += start[v];
		cur[v] = start[v];
	}
	for (i = 0; i < nb; i++)
		ord[cur[m->t[3 * (size_t)E[be[i]].f + E[be[i]].k]]++] = i;
	for (v = 0; v < nv0; v++)
	{
		cur[v] = start[v];
		cnt[v] = v;
	}
	/* a part with no consistent orientation keeps its border: one hole left per connected border */
	for (i = 0; i < nb; i++)
		if (bad[shell[E[be[i]].f]])
		{
			unsigned int a = (unsigned int)(E[be[i]].key >> 32), b = (unsigned int)E[be[i]].key;
			used[i] = 1;
			rp_union(cnt, a, b);
			inpath[a] = inpath[b] = 1;
		}
	for (v = 0; v < nv0; v++)
	{
		w->r->holes_left += inpath[v] && rp_find(cnt, v) == v;
		inpath[v] = 0;
	}
	for (i = 0; i < nb && w->st == QAWS_STATUS_OK; i++)
	{
		unsigned int len = 0, b0 = i, q;
		if (used[i])
			continue;
		for (;;)
		{
			rp_he const* e = E + be[b0];
			unsigned int a = m->t[3 * (size_t)e->f + e->k], b = m->t[3 * (size_t)e->f + (e->k + 1) % 3];
			used[b0] = 1;
			inpath[a] = len + 1;
			pv[len] = a;
			pr[len++] = b0;
			if (inpath[b])
			{
				unsigned int p = inpath[b] - 1;
				rp_fill_loop(w, &c, be, pv + p, pr + p, len - p);
				for (q = p; q < len; q++)
					inpath[pv[q]] = 0;
				len = p;
				if (!len || w->st != QAWS_STATUS_OK)
					break;
			}
			while (cur[b] < start[b + 1] && used[ord[cur[b]]])
				cur[b]++;
			if (cur[b] == start[b + 1])
			{
				/* a chain that does not close (an orientation the walk could not make consistent): left open */
				for (q = 0; q < len; q++)
					inpath[pv[q]] = 0;
				w->r->holes_left++;
				break;
			}
			b0 = ord[cur[b]];
		}
	}
	if (w->st == QAWS_STATUS_OK && m->nt > nt0 && !(flags & QAWS_MESH_REPAIR_COARSE_FILLS))
	{
		rp_refine(w, &c, nt0);
		if (w->st == QAWS_STATUS_OK && !(flags & QAWS_MESH_REPAIR_FLAT_FILLS))
			rp_fair(w, nv0, c.lock);
	}
	w->r->fill_faces += m->nt - nt0;
	goto done;
oom:
	w->st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(E); qaws_internal_free(be); qaws_internal_free(start); qaws_internal_free(cur); qaws_internal_free(ord); qaws_internal_free(inpath); qaws_internal_free(pv); qaws_internal_free(pr); qaws_internal_free(used); qaws_internal_free(cnt);
	qaws_internal_free(c.sigma); qaws_internal_free(c.O); qaws_internal_free(c.lock); qaws_internal_free(c.crease); qaws_internal_free(c.crease_used);
}

/* ---- closed parts by volume ---------------------------------------------- */

/* closed parts (every edge's directed uses cancel) holding less volume than min_volume removed */
static void rp_small_volumes(rp_work* w, double min_volume)
{
	rp_mesh* m = &w->m;
	unsigned int ne, r, k, i, v, f;
	rp_he* E = rp_edges(m, 0, m->nt, &ne);
	unsigned int* uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv ? m->nv : 1)));
	unsigned int* cnt = (unsigned int*)qaws_internal_scratch_zalloc(m->nv ? m->nv : 1, sizeof(unsigned int));
	double* vol = (double*)qaws_internal_scratch_zalloc(m->nv ? m->nv : 1, sizeof(double));
	unsigned char* open = (unsigned char*)qaws_internal_scratch_zalloc(m->nv ? m->nv : 1, 1);
	unsigned char* del = (unsigned char*)qaws_internal_scratch_zalloc(m->nt ? m->nt : 1, 1);
	if (!E || !uf || !cnt || !vol || !open || !del)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (v = 0; v < m->nv; v++)
		uf[v] = v;
	for (f = 0; f < m->nt; f++)
	{
		rp_union(uf, m->t[3 * (size_t)f], m->t[3 * (size_t)f + 1]);
		rp_union(uf, m->t[3 * (size_t)f], m->t[3 * (size_t)f + 2]);
	}
	for (f = 0; f < m->nt; f++)
	{
		double const *a = RP_X(m, m->t[3 * (size_t)f]), *b = RP_X(m, m->t[3 * (size_t)f + 1]), *c = RP_X(m, m->t[3 * (size_t)f + 2]);
		double n[3];
		unsigned int root = rp_find(uf, m->t[3 * (size_t)f]);
		rp_cross(b, c, n);
		vol[root] += mrp_dot(a, n) / 6;
		cnt[root]++;
	}
	for (r = 0; r < ne; r = k)
	{
		int net = 0;
		k = rp_group_end(E, ne, r);
		for (i = r; i < k; i++)
			net += rp_dir(m, E + i) ? 1 : -1;
		if (net)
			open[rp_find(uf, m->t[3 * (size_t)E[r].f])] = 1;
	}
	for (v = 0; v < m->nv; v++)
		if (cnt[v] && !open[v] && fabs(vol[v]) < min_volume)
		{
			w->r->removed_components++;
			w->r->removed_faces += cnt[v];
		}
	for (f = 0; f < m->nt; f++)
	{
		unsigned int root = rp_find(uf, m->t[3 * (size_t)f]);
		del[f] = !open[root] && fabs(vol[root]) < min_volume;
	}
	rp_compact(m, del, m->nt);
done:
	qaws_internal_free(E); qaws_internal_free(uf); qaws_internal_free(cnt); qaws_internal_free(vol); qaws_internal_free(open); qaws_internal_free(del);
}

/* ---- non-manifold edges and vertices split ------------------------------- */

typedef struct rp_round
{
	double angle;
	unsigned int r;             /* the half-edge record */
} rp_round;

static int rp_round_cmp(void const* p, void const* q)
{
	rp_round const *a = (rp_round const*)p, *b = (rp_round const*)q;
	if (a->angle != b->angle)
		return a->angle < b->angle ? -1 : 1;
	return a->r < b->r ? -1 : a->r > b->r;
}

/* the corner of triangle f at vertex v */
static unsigned int rp_corner_at(rp_mesh const* m, unsigned int f, unsigned int v)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
		if (m->t[3 * (size_t)f + k] == v)
			return 3 * f + k;
	return 3 * f;
}

/* Faces paired across each edge: two faces as they are; more, sorted by angle round the edge and paired with the next one
   round where the wedge between them is inside (the first face's normal pointing away from it), so solids touching along an
   edge come apart. Corners joined through the pairs; each vertex's corners in more than one set get a copy per extra set
   (non-manifold edges and vertices both). */
static void rp_split_non_manifold(rp_work* w)
{
	rp_mesh* m = &w->m;
	unsigned int ne, r, k, i, f, v, nv0 = m->nv;
	rp_he* E = rp_edges(m, 0, m->nt, &ne);
	unsigned int* uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (m->nt ? m->nt : 1)));
	unsigned int* first = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nv0 ? nv0 : 1)));
	unsigned int* copy = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (m->nt ? m->nt : 1)));
	rp_round* ring = NULL;
	unsigned int rcap = 0;
	if (!E || !uf || !first || !copy)
		goto oom;
	for (i = 0; i < 3 * m->nt; i++)
		uf[i] = i;
	for (r = 0; r < ne; r = k)
	{
		unsigned int a = (unsigned int)(E[r].key >> 32), b = (unsigned int)E[r].key, n, o, s, best = 0, score[2] = { 0, 0 };
		double u[3], ref[3], ln;
		k = rp_group_end(E, ne, r);
		n = k - r;
		if (n < 2)
			continue;
		if (n > rcap)
		{
			rp_round* g = (rp_round*)qaws_internal_realloc(ring, qaws_internal_mul2(sizeof(rp_round), n));
			if (!g)
				goto oom;
			ring = g;
			rcap = n;
		}
		mrp_sub(RP_X(m, b), RP_X(m, a), u);
		ln = rp_unit(u);
		for (i = 0; i < n; i++)
		{
			unsigned int g = E[r + i].f, c = m->t[3 * (size_t)g + (E[r + i].k + 2) % 3];
			double q[3], t[3], d;
			mrp_sub(RP_X(m, c), RP_X(m, a), q);
			d = mrp_dot(q, u);
			q[0] -= d * u[0];
			q[1] -= d * u[1];
			q[2] -= d * u[2];
			if (!i)
			{
				memcpy(ref, q, sizeof ref);
				rp_unit(ref);
			}
			rp_cross(u, ref, t);
			ring[i].angle = atan2(mrp_dot(q, t), mrp_dot(q, ref));
			ring[i].r = r + i;
		}
		qsort(ring, n, sizeof *ring, rp_round_cmp);
		/* the offset whose pairs open inward: face i's normal against the turn toward face i + 1 */
		for (o = 0; o < 2 && n > 2; o++)
			for (s = o; s + 1 < n + o; s += 2)
			{
				rp_he const* e = E + ring[s % n].r;
				unsigned int g = e->f, c = m->t[3 * (size_t)g + (e->k + 2) % 3];
				double nr[3], q[3], t[3];
				rp_normal(RP_X(m, m->t[3 * (size_t)g]), RP_X(m, m->t[3 * (size_t)g + 1]), RP_X(m, m->t[3 * (size_t)g + 2]), nr);
				mrp_sub(RP_X(m, c), RP_X(m, a), q);
				rp_cross(u, q, t);
				score[o] += mrp_dot(nr, t) < 0 && rp_dir(m, e) != rp_dir(m, E + ring[(s + 1) % n].r);
			}
		best = score[1] > score[0];
		(void)ln;
		for (s = n == 2 ? 0 : best; s + 1 < (n == 2 ? 2 : n + best); s += 2)
		{
			rp_he const *e1 = E + ring[s % n].r, *e2 = E + ring[(s + 1) % n].r;
			rp_union(uf, rp_corner_at(m, e1->f, a), rp_corner_at(m, e2->f, a));
			rp_union(uf, rp_corner_at(m, e1->f, b), rp_corner_at(m, e2->f, b));
		}
	}
	/* each vertex's first corner set keeps it; every further set gets a copy */
	for (v = 0; v < nv0; v++)
		first[v] = HM_NONE;
	for (i = 0; i < 3 * m->nt; i++)
		copy[i] = HM_NONE;
	for (f = 0; f < m->nt; f++)
		for (k = 0; k < 3; k++)
		{
			unsigned int c = 3 * f + k, root = rp_find(uf, c), vv = m->t[c];
			if (copy[root] == HM_NONE)
			{
				if (first[vv] == HM_NONE)
				{
					first[vv] = root;
					copy[root] = vv;
				}
				else
				{
					double p[3];
					memcpy(p, RP_X(m, vv), sizeof p);
					copy[root] = rp_add_v(m, p);
					if (copy[root] == HM_NONE)
						goto oom;
					w->r->split_vertices++;
				}
			}
		}
	for (i = 0; i < 3 * m->nt; i++)
		m->t[i] = copy[rp_find(uf, i)];
	goto done;
oom:
	w->st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(E); qaws_internal_free(uf); qaws_internal_free(first); qaws_internal_free(copy); qaws_internal_free(ring);
}

/* ---- self-intersections, by exact predicates ------------------------------ */

static int rp_o3(double const* a, double const* b, double const* c, double const* d)
{
	qaws_exact_sign s = QAWS_EXACT_ZERO;
	qaws_exact_orient3d(&s, NULL, a, b, c, d);
	return (int)s;
}

static int rp_o2(double const* a, double const* b, double const* c)
{
	qaws_exact_sign s = QAWS_EXACT_ZERO;
	qaws_exact_orient2d(&s, NULL, a, b, c);
	return (int)s;
}

/* segment ab against triangle pqr (closed), not coplanar */
static int rp_seg_tri(double const* a, double const* b, double const* p, double const* q, double const* r)
{
	int sa = rp_o3(p, q, r, a), sb = rp_o3(p, q, r, b), s1, s2, s3;
	if ((sa > 0 && sb > 0) || (sa < 0 && sb < 0) || (!sa && !sb))
		return 0;
	s1 = rp_o3(a, b, p, q);
	s2 = rp_o3(a, b, q, r);
	s3 = rp_o3(a, b, r, p);
	return (s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0);
}

/* 2D: segments ab and cd meet (closed, collinear overlaps included) */
static int rp_seg_seg2(double const* a, double const* b, double const* c, double const* d)
{
	int o1 = rp_o2(a, b, c), o2 = rp_o2(a, b, d), o3 = rp_o2(c, d, a), o4 = rp_o2(c, d, b);
	if (!o1 && !o2 && !o3 && !o4)
	{
		unsigned int k = fabs(b[0] - a[0]) + fabs(d[0] - c[0]) >= fabs(b[1] - a[1]) + fabs(d[1] - c[1]) ? 0 : 1;
		return fmax(fmin(a[k], b[k]), fmin(c[k], d[k])) <= fmin(fmax(a[k], b[k]), fmax(c[k], d[k]));
	}
	return o1 * o2 <= 0 && o3 * o4 <= 0;
}

/* 2D: p in triangle abc (closed, either turn) */
static int rp_in_tri2(double const* p, double const* a, double const* b, double const* c)
{
	int s1 = rp_o2(a, b, p), s2 = rp_o2(b, c, p), s3 = rp_o2(c, a, p);
	return (s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0);
}

/* two triangles in one plane, projected along its dominant axis; shared corners matched by index (s: shared count) */
static int rp_coplanar(double const (*P)[3], double const (*Q)[3], unsigned int const* ip, unsigned int const* iq)
{
	double n[3], p2[3][2], q2[3][2];
	unsigned int ax, i, j, k, x, y;
	rp_normal(P[0], P[1], P[2], n);
	ax = fabs(n[0]) > fabs(n[1]) ? (fabs(n[0]) > fabs(n[2]) ? 0 : 2) : (fabs(n[1]) > fabs(n[2]) ? 1 : 2);
	x = ax == 0 ? 1 : 0;
	y = ax == 2 ? 1 : 2;
	for (k = 0; k < 3; k++)
	{
		p2[k][0] = P[k][x];
		p2[k][1] = P[k][y];
		q2[k][0] = Q[k][x];
		q2[k][1] = Q[k][y];
	}
	/* a corner of one strictly inside the other (shared corners skipped), or edges crossing away from shared corners */
	for (i = 0; i < 3; i++)
	{
		int shared_p = ip[i] == iq[0] || ip[i] == iq[1] || ip[i] == iq[2];
		int shared_q = iq[i] == ip[0] || iq[i] == ip[1] || iq[i] == ip[2];
		if (!shared_p && rp_in_tri2(p2[i], q2[0], q2[1], q2[2]))
			return 1;
		if (!shared_q && rp_in_tri2(q2[i], p2[0], p2[1], p2[2]))
			return 1;
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			unsigned int a = ip[i], b = ip[(i + 1) % 3], c = iq[j], d = iq[(j + 1) % 3];
			if ((a == c || a == d) && (b == c || b == d))
				continue;   /* the shared edge itself */
			if (a == c || a == d || b == c || b == d)
			{
				/* edges from a shared corner: they overlap only when collinear and pointing the same way */
				unsigned int sh = (a == c || a == d) ? i : (i + 1) % 3, op = sh == i ? (i + 1) % 3 : i;
				unsigned int sq = (iq[j] == ip[sh]) ? j : (j + 1) % 3, oq = sq == j ? (j + 1) % 3 : j;
				if (!rp_o2(p2[sh], p2[op], q2[oq])
					&& (p2[op][0] - p2[sh][0]) * (q2[oq][0] - q2[sq][0]) + (p2[op][1] - p2[sh][1]) * (q2[oq][1] - q2[sq][1]) > 0)
					return 1;
				continue;
			}
			if (rp_seg_seg2(p2[i], p2[(i + 1) % 3], q2[j], q2[(j + 1) % 3]))
				return 1;
		}
	return 0;
}

/* triangles P and Q (vertex ids ip, iq) intersect beyond their shared corners and edge, by exact predicates */
int rp_tri_tri_at(double const (*P)[3], double const (*Q)[3], unsigned int const* ip, unsigned int const* iq)
{
	unsigned int shared = 0, i, j;
	int sq[3], cop;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			shared += ip[i] == iq[j];
	if (shared == 3)
		return 0;
	for (i = 0; i < 3; i++)   /* (a shared corner is on the plane: no exact zero to grind out) */
		sq[i] = iq[i] == ip[0] || iq[i] == ip[1] || iq[i] == ip[2] ? 0 : rp_o3(P[0], P[1], P[2], Q[i]);
	cop = !sq[0] && !sq[1] && !sq[2];
	if (cop)
		return rp_coplanar((double const(*)[3])P, (double const(*)[3])Q, ip, iq);
	if (shared == 2)
		return 0;   /* two faces on an edge meet only along it unless they lie in one plane */
	if ((sq[0] > 0 && sq[1] > 0 && sq[2] > 0) || (sq[0] < 0 && sq[1] < 0 && sq[2] < 0))
		return 0;
	/* the edges of each against the other; with a shared corner, only the edges away from it */
	for (i = 0; i < 3; i++)
	{
		unsigned int a = ip[i], b = ip[(i + 1) % 3];
		int sa = a == iq[0] || a == iq[1] || a == iq[2], sb = b == iq[0] || b == iq[1] || b == iq[2];
		if (!sa && !sb && rp_seg_tri(P[i], P[(i + 1) % 3], Q[0], Q[1], Q[2]))
			return 1;
		a = iq[i];
		b = iq[(i + 1) % 3];
		sa = a == ip[0] || a == ip[1] || a == ip[2];
		sb = b == ip[0] || b == ip[1] || b == ip[2];
		if (!sa && !sb && rp_seg_tri(Q[i], Q[(i + 1) % 3], P[0], P[1], P[2]))
			return 1;
	}
	return 0;
}

static int rp_tri_tri(rp_mesh const* m, unsigned int f, unsigned int g)
{
	unsigned int const *ip = m->t + 3 * (size_t)f, *iq = m->t + 3 * (size_t)g;
	double P[3][3], Q[3][3];
	unsigned int i;
	for (i = 0; i < 3; i++)
	{
		memcpy(P[i], RP_X(m, ip[i]), sizeof P[i]);
		memcpy(Q[i], RP_X(m, iq[i]), sizeof Q[i]);
	}
	return rp_tri_tri_at((double const(*)[3])P, (double const(*)[3])Q, ip, iq);
}

typedef struct rp_hit_ctx
{
	rp_mesh const* m;
	unsigned int* hits;         /* per triangle */
	unsigned int pairs;
} rp_hit_ctx;

static qaws_status rp_hit_visit(void* user, unsigned int i, unsigned int j)
{
	rp_hit_ctx* c = (rp_hit_ctx*)user;
	if (rp_tri_tri(c->m, i, j))
	{
		c->pairs++;
		c->hits[i]++;
		c->hits[j]++;
	}
	return QAWS_STATUS_OK;
}

static qaws_status rp_hit_visit_chunk(void* user, unsigned int chunk, unsigned int i, unsigned int j)
{
	(void)chunk;
	return rp_hit_visit(user, i, j);
}

/* every pair of triangles that cross (beyond shared corners and edges), counted per triangle into w's hits */
static void rp_intersections(rp_work* w, unsigned int** out_hits)
{
	rp_mesh* m = &w->m;
	qaws_bp_box* box = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (m->nt ? m->nt : 1)));
	rp_hit_ctx c;
	unsigned int f, k, j;
	c.m = m;
	c.pairs = 0;
	c.hits = (unsigned int*)qaws_internal_zalloc(m->nt ? m->nt : 1, sizeof(unsigned int));
	*out_hits = c.hits;
	if (!box || !c.hits)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		qaws_internal_free(box);
		return;
	}
	for (f = 0; f < m->nt; f++)
		for (k = 0; k < 3; k++)
		{
			double const* p = RP_X(m, m->t[3 * (size_t)f + k]);
			for (j = 0; j < 3; j++)
			{
				box[f].lo[j] = k ? fmin(box[f].lo[j], p[j]) : p[j];
				box[f].hi[j] = k ? fmax(box[f].hi[j], p[j]) : p[j];
			}
		}
	if (m->nt > 1)
	{
		/* cells fitted to the triangles: a surface (or a seam of pieces) fills little of its box */
		qaws_bp_grid* grid = NULL;
		w->st = qaws_internal_grid_create_fit(&grid, box, m->nt, 3, NULL);
		if (w->st == QAWS_STATUS_OK)
			w->st = qaws_internal_grid_pairs(NULL, grid, box, NULL, NULL, rp_hit_visit_chunk, &c, NULL);
		qaws_internal_grid_destroy(grid);
	}
	w->r->self_intersections = c.pairs;
	for (f = 0; f < m->nt; f++)
		w->r->intersecting_faces += c.hits[f] > 0;
	qaws_internal_free(box);
}

/* ---- the call ----------------------------------------------------------- */

/* what is left: boundary, non-manifold and open edges */
static void rp_measure(rp_work* w)
{
	rp_mesh* m = &w->m;
	unsigned int ne, r, k, i;
	rp_he* E = rp_edges(m, 0, m->nt, &ne);
	if (!E)
	{
		w->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	for (r = 0; r < ne; r = k)
	{
		int net = 0;
		k = rp_group_end(E, ne, r);
		for (i = r; i < k; i++)
			net += rp_dir(m, E + i) ? 1 : -1;
		w->r->boundary_edges += k - r == 1;
		w->r->non_manifold_edges += k - r > 2;
		w->r->open_edges += net != 0;
	}
	w->r->closed = m->nt > 0 && !w->r->open_edges;
	for (i = 0; i < m->nt; i++)
		w->r->flipped_faces += m->src[i] >= 0 && m->fl[i];
	qaws_internal_free(E);
}

/* the used vertices renumbered in order, the store built with the source_face attribute */
static qaws_status rp_output(rp_work* w, qaws_mesh_3d** out, unsigned int const* hits)
{
	rp_mesh* m = &w->m;
	unsigned int* map = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv ? m->nv : 1)));
	double* x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (m->nv ? m->nv : 1)));
	unsigned int v, n = 0, a;
	size_t i;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!map || !x)
		goto done;
	for (v = 0; v < m->nv; v++)
		map[v] = HM_NONE;
	for (i = 0; i < 3 * (size_t)m->nt; i++)
		map[m->t[i]] = 0;
	for (v = 0; v < m->nv; v++)
		if (map[v] == 0)
		{
			memcpy(x + 3 * (size_t)n, RP_X(m, v), sizeof(double) * 3);
			map[v] = n++;
		}
	for (i = 0; i < 3 * (size_t)m->nt; i++)
		m->t[i] = map[m->t[i]];
	st = hm_mesh_3d_create(x, n, m->t, NULL, m->nt, out);
	if (st != QAWS_STATUS_OK)
		goto done;
	a = hm_attr_add(&(*out)->hm, "source_face", QAWS_MESH_FACE, QAWS_ATTRIBUTE_NEAREST, 1);
	if (a == HM_NONE)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (v = 0; v < m->nt; v++)
		*hm_attr_at(&(*out)->hm, a, v) = m->src[v];
	if (hits)
	{
		a = hm_attr_add(&(*out)->hm, "intersecting", QAWS_MESH_FACE, QAWS_ATTRIBUTE_NEAREST, 1);
		if (a == HM_NONE)
		{
			qaws_mesh_3d_destroy(*out);
			*out = NULL;
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (v = 0; v < m->nt; v++)
			*hm_attr_at(&(*out)->hm, a, v) = hits[v];
	}
done:
	qaws_internal_free(map);
	qaws_internal_free(x);
	return st;
}

qaws_status qaws_mesh_repair(qaws_mesh_3d** out_mesh, qaws_mesh_repair_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_mesh_repair_desc const* desc)
{
	qaws_mesh_repair_desc d;
	qaws_mesh_repair_report local;
	rp_work w;
	rp_mesh* m = &w.m;
	unsigned int *uf = NULL, *poly = NULL, *stamp = NULL, *shell = NULL, *hits = NULL;
	unsigned char *used = NULL, *bad = NULL;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, D = 0, weld, degen, gap, perim;
	size_t off, i;
	unsigned int f, smax = 3, k, v;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	report->invalid_face = report->invalid_vertex = QAWS_MESH_NONE;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if ((vertex_count && !positions) || (face_count && !face_vertices))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc)
		d = *desc;
	else
		memset(&d, 0, sizeof d);
	if ((unsigned int)d.nested_shells > QAWS_NESTED_SHELLS_REMOVE || d.weld_tolerance != d.weld_tolerance
		|| d.degenerate_tolerance != d.degenerate_tolerance || d.gap_tolerance != d.gap_tolerance || d.hole_perimeter != d.hole_perimeter
		|| d.min_component_area != d.min_component_area || d.min_component_volume != d.min_component_volume)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * (size_t)vertex_count; i++)
		if (!rp_finite((double)positions[i]))
		{
			report->invalid_vertex = (unsigned int)(i / 3);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
	for (f = 0, off = 0; f < face_count; f++)
	{
		unsigned int s = face_sizes ? face_sizes[f] : 3;
		if (s < 3)
		{
			report->invalid_face = f;
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
		for (k = 0; k < s; k++)
			if (face_vertices[off + k] >= vertex_count)
			{
				report->invalid_face = f;
				return QAWS_STATUS_INVALID_ARGUMENT;
			}
		off += s;
		smax = s > smax ? s : smax;
	}
	memset(&w, 0, sizeof w);
	w.r = report;
	m->nv = m->capv = vertex_count;
	m->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (vertex_count ? vertex_count : 1)));
	used = (unsigned char*)qaws_internal_scratch_zalloc(vertex_count ? vertex_count : 1, 1);
	stamp = (unsigned int*)qaws_internal_scratch_zalloc(vertex_count ? vertex_count : 1, sizeof(unsigned int));
	poly = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), smax));
	if (!m->x || !used || !stamp || !poly)
		goto oom;
	for (i = 0; i < 3 * (size_t)vertex_count; i++)
		m->x[i] = (double)positions[i];
	for (i = 0; i < off; i++)
	{
		double const* p = RP_X(m, face_vertices[i]);
		used[face_vertices[i]] = 1;
		for (k = 0; k < 3; k++)
		{
			lo[k] = fmin(lo[k], p[k]);
			hi[k] = fmax(hi[k], p[k]);
		}
	}
	for (v = 0; v < vertex_count; v++)
		report->unused_vertices += !used[v];
	if (off)
		D = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
	weld = d.weld_tolerance > 0 ? (double)d.weld_tolerance : d.weld_tolerance == 0 ? 1e-6 * D : 0;
	degen = d.degenerate_tolerance > 0 ? (double)d.degenerate_tolerance : d.degenerate_tolerance < 0 ? -1 : d.weld_tolerance >= 0 ? weld : 1e-6 * D;
	gap = d.gap_tolerance > 0 ? (double)d.gap_tolerance : d.gap_tolerance == 0 ? 1e-3 * D : -1;
	perim = d.hole_perimeter > 0 ? (double)d.hole_perimeter : d.hole_perimeter == 0 ? HUGE_VAL : -1;
	uf = rp_weld(&w, used, weld);
	if (!uf)
		goto oom;
	/* the faces: a vertex twice as given or after welding dropped, the others cut into triangles */
	for (f = 0, off = 0; f < face_count && w.st == QAWS_STATUS_OK; f++)
	{
		unsigned int s = face_sizes ? face_sizes[f] : 3, n = 0, twice = 0;
		for (k = 0; k < s; k++)
		{
			v = face_vertices[off + k];
			twice |= stamp[v] == 2 * f + 1;
			stamp[v] = 2 * f + 1;
		}
		for (k = 0; k < s && !twice; k++)
		{
			v = rp_find(uf, face_vertices[off + k]);
			if (!n || poly[n - 1] != v)
				poly[n++] = v;
		}
		while (n > 1 && poly[n - 1] == poly[0])
			n--;
		for (k = 0; k < n && !twice; k++)
		{
			twice |= stamp[poly[k]] == 2 * f + 2;
			stamp[poly[k]] = 2 * f + 2;
		}
		off += s;
		if (twice || n < 3)
			report->dropped_faces++;
		else if (n == 3)
		{
			if (!rp_add_t(m, poly[0], poly[1], poly[2], (int)f, 0))
				goto oom;
		}
		else
			rp_ear_clip(&w, poly, n, (int)f);
	}
	if (w.st == QAWS_STATUS_OK)
		rp_clean(&w);
	if (w.st == QAWS_STATUS_OK && degen >= 0)
		rp_degenerate(&w, degen);
	if (w.st == QAWS_STATUS_OK && gap >= 0)
	{
		rp_stitch(&w, gap);
		if (w.st == QAWS_STATUS_OK && degen >= 0)
			rp_degenerate(&w, degen);
	}
	if (w.st == QAWS_STATUS_OK && (d.min_component_area > 0 || d.min_component_faces > 0))
		rp_components(&w, (double)d.min_component_area, d.min_component_faces);
	if (w.st == QAWS_STATUS_OK && perim >= 0)
	{
		shell = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nt ? m->nt : 1)));
		if (!shell)
			goto oom;
		rp_consistent(&w, shell, &bad);
		if (w.st == QAWS_STATUS_OK)
			rp_fill(&w, perim, d.flags, shell, bad);
	}
	if (w.st == QAWS_STATUS_OK)
		rp_outward(&w, d.nested_shells == QAWS_NESTED_SHELLS_REMOVE);
	if (w.st == QAWS_STATUS_OK && d.min_component_volume > 0)
		rp_small_volumes(&w, (double)d.min_component_volume);
	if (w.st == QAWS_STATUS_OK && (d.flags & QAWS_MESH_REPAIR_SPLIT_NON_MANIFOLD))
		rp_split_non_manifold(&w);
	if (w.st == QAWS_STATUS_OK && (d.flags & QAWS_MESH_REPAIR_FIND_INTERSECTIONS))
		rp_intersections(&w, &hits);
	if (w.st == QAWS_STATUS_OK)
		rp_measure(&w);
	if (w.st == QAWS_STATUS_OK)
		w.st = rp_output(&w, out_mesh, hits);
	goto done;
oom:
	w.st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(m->x); qaws_internal_free(m->t); qaws_internal_free(m->src); qaws_internal_free(m->fl);
	qaws_internal_free(uf); qaws_internal_free(poly); qaws_internal_free(stamp); qaws_internal_free(used); qaws_internal_free(shell); qaws_internal_free(bad); qaws_internal_free(hits);
	return w.st;
}
