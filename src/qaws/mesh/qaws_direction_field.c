/*
 * Direction and position fields after Instant Meshes (Jakob, Tarini,
 * Panozzo, Sorkine-Hornung 2015): a hierarchy of the vertices (pairs of
 * similar normals and areas merged level by level), a 4-symmetric
 * orientation field and a lattice position field smoothed coarse to fine by
 * Gauss-Seidel, feature edges and the border as constraints.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include "qaws_cross_field.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double mf_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void mf_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double mf_unit(double* v)
{
	double l = sqrt(mf_dot(v, v));
	if (l > 0)
	{
		v[0] /= l; v[1] /= l; v[2] /= l;
	}
	return l;
}

/* v projected on the plane of unit normal n, made unit (any tangent when it vanishes) */
static void mf_tangent(double* v, double const* n)
{
	double d = mf_dot(v, n);
	v[0] -= d * n[0]; v[1] -= d * n[1]; v[2] -= d * n[2];
	if (mf_unit(v) == 0)
	{
		double a[3] = { 1, 0, 0 };
		if (fabs(n[0]) > 0.9)
		{
			a[0] = 0;
			a[1] = 1;
		}
		mf_cross(n, a, v);
		mf_unit(v);
	}
}

void mf_compat_orientation(double const* q0, double const* n0, double const* q1, double const* n1, double* out0, double* out1)
{
	double A[2][3], B[2][3], best = -1, dp;
	int i, j, ba = 0, bb = 0, k;
	memcpy(A[0], q0, sizeof A[0]);
	mf_cross(n0, q0, A[1]);
	memcpy(B[0], q1, sizeof B[0]);
	mf_cross(n1, q1, B[1]);
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			double s = fabs(mf_dot(A[i], B[j]));
			if (s > best)
			{
				best = s;
				ba = i;
				bb = j;
			}
		}
	dp = mf_dot(A[ba], B[bb]);
	for (k = 0; k < 3; k++)
	{
		out0[k] = A[ba][k];
		out1[k] = dp < 0 ? -B[bb][k] : B[bb][k];
	}
}

int mf_rotation(double const* q0, double const* n0, double const* q1, double const* n1)
{
	double A[2][3], B[2][3], best = -1;
	int i, j, ba = 0, bb = 0;
	memcpy(A[0], q0, sizeof A[0]);
	mf_cross(n0, q0, A[1]);
	memcpy(B[0], q1, sizeof B[0]);
	mf_cross(n1, q1, B[1]);
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			double s = fabs(mf_dot(A[i], B[j]));
			if (s > best)
			{
				best = s;
				ba = i;
				bb = j;
			}
		}
	/* R^a q0 = R^(b + 2 (sign < 0)) q1: q1 turned by b' - a quarter turns is q0 */
	return ((bb + (mf_dot(A[ba], B[bb]) < 0 ? 2 : 0)) - ba + 4) % 4;
}

/* a point near both tangent planes (Instant Meshes' middle point) */
static void mf_middle(double const* p0, double const* n0, double const* p1, double const* n1, double* out)
{
	double n0p1 = mf_dot(n0, p1), n0p0 = mf_dot(n0, p0), n1p0 = mf_dot(n1, p0), n1p1 = mf_dot(n1, p1), n0n1 = mf_dot(n0, n1);
	double den = 1 / (1 - n0n1 * n0n1 + 1e-4);
	double l0 = 2 * (n0p1 - n0p0 - n0n1 * (n1p0 - n1p1)) * den, l1 = 2 * (n1p0 - n1p1 - n0n1 * (n0p1 - n0p0)) * den;
	int k;
	for (k = 0; k < 3; k++)
		out[k] = 0.5 * (p0[k] + p1[k]) - 0.25 * (n0[k] * l0 + n1[k] * l1);
}

void mf_unfold(double const* x0, double const* n0, double const* x1, double const* n1, double const* p, double* out)
{
	double mid[3], axis[3], v[3], w[3], c = mf_dot(n1, n0), s2, f;
	int k;
	mf_cross(n1, n0, axis);
	s2 = mf_dot(axis, axis);
	if (c > 0.9999 || s2 < 1e-12)
	{
		memcpy(out, p, sizeof(double) * 3);   /* (alike normals: nothing to unfold; opposite ones: no hinge) */
		return;
	}
	mf_middle(x0, n0, x1, n1, mid);
	for (k = 0; k < 3; k++)
		v[k] = p[k] - mid[k];
	/* Rodrigues with the unnormalized axis n1 x n0 (its length the sine) */
	f = (1 - c) / s2;
	mf_cross(axis, v, w);
	for (k = 0; k < 3; k++)
		out[k] = mid[k] + v[k] * c + w[k] + axis[k] * mf_dot(axis, v) * f;
}

/* the lattice point of (o, q, t) below p (floor) or nearest it (round) */
static void mf_lattice(double const* o, double const* q, double const* t, double const* p, double scale, int round_it, double* out)
{
	double d[3], a, b;
	int k;
	for (k = 0; k < 3; k++)
		d[k] = p[k] - o[k];
	a = mf_dot(q, d) / scale;
	b = mf_dot(t, d) / scale;
	a = round_it ? floor(a + 0.5) : floor(a);
	b = round_it ? floor(b + 0.5) : floor(b);
	for (k = 0; k < 3; k++)
		out[k] = o[k] + (q[k] * a + t[k] * b) * scale;
}

void mf_compat_position(double const* x0, double const* n0, double const* q0, double const* o0, double const* x1, double const* n1,
	double const* q1, double const* o1, double scale, double* out0, double* out1)
{
	double t0[3], t1[3], mid[3], a[3], b[3], best = HUGE_VAL;
	int i, j, k;
	mf_cross(n0, q0, t0);
	mf_cross(n1, q1, t1);
	mf_middle(x0, n0, x1, n1, mid);
	mf_lattice(o0, q0, t0, mid, scale, 0, a);
	mf_lattice(o1, q1, t1, mid, scale, 0, b);
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			double u[3], v[3], c = 0;
			for (k = 0; k < 3; k++)
			{
				u[k] = a[k] + (q0[k] * (i & 1) + t0[k] * ((i & 2) >> 1)) * scale;
				v[k] = b[k] + (q1[k] * (j & 1) + t1[k] * ((j & 2) >> 1)) * scale;
				c += (u[k] - v[k]) * (u[k] - v[k]);
			}
			if (c < best)
			{
				best = c;
				memcpy(out0, u, sizeof u);
				memcpy(out1, v, sizeof v);
			}
		}
}

/* ======================================================================== */
/*  The hierarchy                                                           */
/* ======================================================================== */

static void mf_level_free(mf_level* L)
{
	qaws_internal_free(L->x); qaws_internal_free(L->nrm); qaws_internal_free(L->area); qaws_internal_free(L->adj0); qaws_internal_free(L->adj); qaws_internal_free(L->parent);
	qaws_internal_free(L->child0); qaws_internal_free(L->child); qaws_internal_free(L->q); qaws_internal_free(L->o); qaws_internal_free(L->cq); qaws_internal_free(L->con); qaws_internal_free(L->sc); qaws_internal_free(L->wt);
	memset(L, 0, sizeof *L);
}

void mf_free(mf_hier* h)
{
	unsigned int i;
	for (i = 0; i < h->nl; i++)
		mf_level_free(&h->lv[i]);
	qaws_internal_free(h->lv);
	memset(h, 0, sizeof *h);
}

static int mf_level_alloc(mf_level* L, unsigned int n)
{
	memset(L, 0, sizeof *L);
	L->n = n;
	L->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (n + 1)));
	L->nrm = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (n + 1)));
	L->area = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (n + 1)));
	L->parent = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (n + 1)));
	L->q = (double*)qaws_internal_zalloc(3 * (size_t)n + 3, sizeof(double));
	L->o = (double*)qaws_internal_zalloc(3 * (size_t)n + 3, sizeof(double));
	L->cq = (double*)qaws_internal_zalloc(3 * (size_t)n + 3, sizeof(double));
	L->con = (unsigned char*)qaws_internal_zalloc((size_t)n + 1, 1);
	L->adj0 = (unsigned int*)qaws_internal_zalloc((size_t)n + 1, sizeof(unsigned int));
	return L->x && L->nrm && L->area && L->parent && L->q && L->o && L->cq && L->con && L->adj0;
}

typedef struct mf_pair
{
	double score;
	unsigned int i, j;
} mf_pair;

static int mf_pair_cmp(void const* a, void const* b)
{
	double sa = ((mf_pair const*)a)->score, sb = ((mf_pair const*)b)->score;
	return sa > sb ? -1 : sa < sb ? 1 : 0;
}

/* the next level: pairs of neighbours with like normals and areas merged */
static int mf_coarsen(mf_level* f, mf_level* c)
{
	unsigned int ne = f->adj0[f->n], i, k, np = 0, nc = 0, *stamp = NULL, total;
	mf_pair* pairs = (mf_pair*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(mf_pair), (ne / 2 + 1)));
	unsigned int* first = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (f->n + 1)));
	unsigned int* second = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (f->n + 1)));
	int ok = 0;
	if (!pairs || !first || !second)
		goto done;
	for (i = 0; i < f->n; i++)
		for (k = f->adj0[i]; k < f->adj0[i + 1]; k++)
		{
			unsigned int j = f->adj[k];
			double ai = f->area[i], aj = f->area[j];
			if (j <= i)
				continue;
			pairs[np].i = i;
			pairs[np].j = j;
			pairs[np].score = mf_dot(f->nrm + 3 * (size_t)i, f->nrm + 3 * (size_t)j) * (ai < aj ? ai / aj : aj / ai);
			np++;
		}
	qsort(pairs, np, sizeof(mf_pair), mf_pair_cmp);
	for (i = 0; i < f->n; i++)
		f->parent[i] = HM_NONE;
	for (k = 0; k < np; k++)
	{
		unsigned int a = pairs[k].i, b = pairs[k].j;
		if (f->parent[a] != HM_NONE || f->parent[b] != HM_NONE)
			continue;
		first[nc] = a;
		second[nc] = b;
		f->parent[a] = f->parent[b] = nc++;
	}
	for (i = 0; i < f->n; i++)
		if (f->parent[i] == HM_NONE)
		{
			first[nc] = i;
			second[nc] = HM_NONE;
			f->parent[i] = nc++;
		}
	if (!mf_level_alloc(c, nc))
		goto done;
	c->child0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nc + 1)));
	c->child = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (2 * (size_t)nc + 1)));
	stamp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nc + 1)));
	if (!c->child0 || !c->child || !stamp)
		goto done;
	for (i = 0, total = 0; i < nc; i++)
	{
		unsigned int a = first[i], b = second[i], q;
		double A = f->area[a] + (b != HM_NONE ? f->area[b] : 0);
		c->child0[i] = total;
		c->child[total++] = a;
		if (b != HM_NONE)
			c->child[total++] = b;
		c->area[i] = A;
		for (q = 0; q < 3; q++)
		{
			c->x[3 * i + q] = (f->area[a] * f->x[3 * a + q] + (b != HM_NONE ? f->area[b] * f->x[3 * b + q] : 0)) / (A > 0 ? A : 1);
			c->nrm[3 * i + q] = f->area[a] * f->nrm[3 * a + q] + (b != HM_NONE ? f->area[b] * f->nrm[3 * b + q] : 0);
		}
		if (A <= 0)
			memcpy(c->x + 3 * i, f->x + 3 * a, sizeof(double) * 3);
		if (mf_unit(c->nrm + 3 * i) == 0)
			memcpy(c->nrm + 3 * i, f->nrm + 3 * a, sizeof(double) * 3);
		/* an orientation constraint from a constrained child */
		for (q = 0; q < 2; q++)
		{
			unsigned int ch = q ? b : a;
			if (ch == HM_NONE || !(f->con[ch] & MF_ALONG))
				continue;
			memcpy(c->cq + 3 * i, f->cq + 3 * ch, sizeof(double) * 3);
			mf_tangent(c->cq + 3 * i, c->nrm + 3 * i);
			c->con[i] = MF_ALONG;
			break;
		}
	}
	c->child0[nc] = total;
	/* neighbours: the children's, as parents, once each */
	for (i = 0; i < nc; i++)
		stamp[i] = HM_NONE;
	{
		unsigned int cap = ne + 1, cnt = 0;
		c->adj = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), cap));
		if (!c->adj)
			goto done;
		for (i = 0; i < nc; i++)
		{
			unsigned int q;
			c->adj0[i] = cnt;
			stamp[i] = i;
			for (q = c->child0[i]; q < c->child0[i + 1]; q++)
			{
				unsigned int ch = c->child[q];
				for (k = f->adj0[ch]; k < f->adj0[ch + 1]; k++)
				{
					unsigned int p = f->parent[f->adj[k]];
					if (stamp[p] == i)
						continue;
					stamp[p] = i;
					c->adj[cnt++] = p;
				}
			}
		}
		c->adj0[nc] = cnt;
	}
	ok = 1;
done:
	qaws_internal_free(pairs); qaws_internal_free(first); qaws_internal_free(second); qaws_internal_free(stamp);
	return ok;
}

qaws_status mf_build(mf_hier* h, qaws_hm const* m, double const* normals)
{
	mf_level* L;
	unsigned int v, f, c, cnt = 0, cap = 8;
	memset(h, 0, sizeof *h);
	h->lv = (mf_level*)qaws_internal_zalloc(cap, sizeof(mf_level));
	if (!h->lv)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	h->nl = 1;
	L = &h->lv[0];
	if (!mf_level_alloc(L, m->nv))
		goto oom;
	memcpy(L->x, m->x, sizeof(double) * 3 * m->nv);
	/* normals and areas: angle-weighted face normals, a third of each face's area */
	for (v = 0; v < 3 * m->nv; v++)
		L->nrm[v] = 0;
	for (v = 0; v < m->nv; v++)
		L->area[v] = 0;
	for (f = 0; f < m->nf; f++)
	{
		double n[3], a;
		hm_face_normal(m, f, n);
		a = mf_unit(n);
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
		{
			unsigned int k, vv = m->cv[c];
			double const *p = m->x + 3 * (size_t)vv, *pa = m->x + 3 * (size_t)m->cv[hm_prev(m, c)], *pb = m->x + 3 * (size_t)hm_to(m, c);
			double u[3], t[3], ang, d;
			for (k = 0; k < 3; k++)
			{
				u[k] = pa[k] - p[k];
				t[k] = pb[k] - p[k];
			}
			if (mf_unit(u) == 0 || mf_unit(t) == 0)
				continue;
			d = mf_dot(u, t);
			ang = acos(d > 1 ? 1 : d < -1 ? -1 : d);
			for (k = 0; k < 3; k++)
				L->nrm[3 * (size_t)vv + k] += ang * n[k];
			L->area[vv] += a / m->fn[f];
		}
	}
	for (v = 0; v < m->nv; v++)
	{
		if (normals)
			memcpy(L->nrm + 3 * (size_t)v, normals + 3 * (size_t)v, sizeof(double) * 3);
		if (mf_unit(L->nrm + 3 * (size_t)v) == 0)
			L->nrm[3 * (size_t)v + 2] = 1;
	}
	/* neighbours */
	for (v = 0; v < m->nv; v++)
	{
		L->adj0[v] = cnt;
		cnt += hm_valence(m, v);
	}
	L->adj0[m->nv] = cnt;
	L->adj = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (cnt + 1)));
	if (!L->adj)
		goto oom;
	for (v = 0; v < m->nv; v++)
		hm_ring(m, v, L->adj + L->adj0[v], L->adj0[v + 1] - L->adj0[v]);
	/* constraints: along feature edges and the border; corners where lines meet, end or turn */
	for (v = 0; v < m->nv; v++)
	{
		unsigned int w[2], n = 0, last = HM_NONE;
		if (m->vc[v] == HM_NONE)
			continue;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			if (((m->cflag[c] & HM_FEATURE) || m->co[c] == HM_NONE) && n < 2)
				w[n] = hm_to(m, c);
			n += (m->cflag[c] & HM_FEATURE) || m->co[c] == HM_NONE;
			last = c;
		}
		if (last != HM_NONE && hm_is_border_vertex(m, v))
		{
			if (n < 2)
				w[n] = m->cv[hm_prev(m, last)];
			n++;
		}
		if (n == 0)
		{
			if (m->vflag[v] & HM_CORNER)
				L->con[v] = MF_PINNED;   /* (a corner off the lines: a pinned vertex) */
			continue;
		}
		{
			double* q = L->cq + 3 * (size_t)v;
			unsigned int k;
			for (k = 0; k < 3; k++)
				q[k] = n == 2 ? m->x[3 * (size_t)w[1] + k] - m->x[3 * (size_t)w[0] + k] : m->x[3 * (size_t)w[0] + k] - m->x[3 * (size_t)v + k];
			mf_tangent(q, L->nrm + 3 * (size_t)v);
			L->con[v] = MF_ALONG;
			if (n != 2 || (m->vflag[v] & HM_CORNER))
				L->con[v] |= MF_PINNED;
			else
			{
				/* a line turning by more than 60 degrees: a corner too */
				double u[3], s[3];
				for (k = 0; k < 3; k++)
				{
					u[k] = m->x[3 * (size_t)v + k] - m->x[3 * (size_t)w[0] + k];
					s[k] = m->x[3 * (size_t)w[1] + k] - m->x[3 * (size_t)v + k];
				}
				if (mf_dot(u, s) < 0.5 * sqrt(mf_dot(u, u) * mf_dot(s, s)))
					L->con[v] |= MF_PINNED;
				else
					L->con[v] |= MF_ON_LINE;
			}
		}
	}
	/* coarser levels down to a few dozen vertices */
	while (h->lv[h->nl - 1].n > 32)
	{
		mf_level *fine, *coarse;
		if (h->nl == cap)
		{
			unsigned int nc = qaws_internal_grow(cap, (size_t)h->nl + 1, 8);
			mf_level* g = nc ? (mf_level*)qaws_internal_realloc(h->lv, qaws_internal_mul2(sizeof(mf_level), nc)) : NULL;
			if (!g)
				goto oom;
			memset(g + cap, 0, sizeof(mf_level) * (nc - cap));
			h->lv = g;
			cap = nc;
		}
		fine = &h->lv[h->nl - 1];
		coarse = &h->lv[h->nl];
		if (!mf_coarsen(fine, coarse))
		{
			mf_level_free(coarse);
			goto oom;
		}
		h->nl++;
		if (coarse->n > 0.9 * fine->n)
			break;   /* (no longer shrinking) */
	}
	h->lv[h->nl - 1].parent[0] = HM_NONE;
	return QAWS_STATUS_OK;
oom:
	mf_free(h);
	return QAWS_STATUS_ALLOCATION_FAILURE;
}

/* ======================================================================== */
/*  Smoothing, coarse to fine                                               */
/* ======================================================================== */

static void mf_sweep_q(mf_level* L, unsigned int sweeps)
{
	unsigned int s, i, k;
	for (s = 0; s < sweeps; s++)
		for (i = 0; i < L->n; i++)
		{
			double const* n = L->nrm + 3 * (size_t)i;
			double qi[3], ws = 0;
			if (L->con[i] & MF_ALONG)
			{
				memcpy(L->q + 3 * (size_t)i, L->cq + 3 * (size_t)i, sizeof qi);
				continue;
			}
			memcpy(qi, L->q + 3 * (size_t)i, sizeof qi);
			for (k = L->adj0[i]; k < L->adj0[i + 1]; k++)
			{
				unsigned int j = L->adj[k], q;
				double a[3], b[3], wj = L->wt ? 0.5 * (L->wt[i] + L->wt[j]) : 1;
				mf_compat_orientation(qi, n, L->q + 3 * (size_t)j, L->nrm + 3 * (size_t)j, a, b);
				for (q = 0; q < 3; q++)
					qi[q] = a[q] * ws + b[q] * wj;
				ws += wj;
				mf_tangent(qi, n);
			}
			if (L->con[i] & MF_PREFER)
			{
				double a[3], b[3];
				unsigned int q;
				mf_compat_orientation(qi, n, L->cq + 3 * (size_t)i, n, a, b);
				for (q = 0; q < 3; q++)
					qi[q] = (a[q] + b[q]) * ws;   /* (b as long as cq: |cq| times the neighbours) */
				mf_tangent(qi, n);
			}
			memcpy(L->q + 3 * (size_t)i, qi, sizeof qi);
		}
}

void mf_orientations(mf_hier* h, unsigned int sweeps, double const* init, int pull)
{
	int l;
	unsigned int i, k, c;
	mf_level* top = &h->lv[h->nl - 1];
	/* a preferred start (level 0) carried up: each vertex its children's, matched to the first one's */
	if (init)
	{
		mf_level* L0 = &h->lv[0];
		memcpy(L0->q, init, sizeof(double) * 3 * L0->n);
		for (i = 0; i < L0->n; i++)
			if (pull && !(L0->con[i] & MF_ALONG) && (init[3 * (size_t)i] != 0 || init[3 * (size_t)i + 1] != 0 || init[3 * (size_t)i + 2] != 0))
			{
				memcpy(L0->cq + 3 * (size_t)i, init + 3 * (size_t)i, sizeof(double) * 3);
				L0->con[i] |= MF_PREFER;
			}
		for (l = 1; l < (int)h->nl; l++)
		{
			mf_level *L = &h->lv[l], *F = &h->lv[l - 1];
			for (i = 0; i < L->n; i++)
			{
				double* q = L->q + 3 * (size_t)i;
				double ref[3] = { 0, 0, 0 }, w = 0;
				int have = 0;
				q[0] = q[1] = q[2] = 0;
				for (c = L->child0[i]; c < L->child0[i + 1]; c++)
				{
					unsigned int ch = L->child[c];
					double const* cq = F->q + 3 * (size_t)ch;
					double a[3], b[3];
					if (F->con[ch] & MF_PREFER)
						w += sqrt(mf_dot(F->cq + 3 * (size_t)ch, F->cq + 3 * (size_t)ch)) / (L->child0[i + 1] - L->child0[i]);
					if (cq[0] == 0 && cq[1] == 0 && cq[2] == 0)
						continue;
					if (!have)
					{
						memcpy(ref, cq, sizeof ref);
						have = 1;
					}
					mf_compat_orientation(ref, L->nrm + 3 * (size_t)i, cq, F->nrm + 3 * (size_t)ch, a, b);
					for (k = 0; k < 3; k++)
						q[k] += b[k];
				}
				if (have)
				{
					mf_tangent(q, L->nrm + 3 * (size_t)i);
					if (pull && w > 0 && !(L->con[i] & MF_ALONG))
					{
						for (k = 0; k < 3; k++)
							L->cq[3 * (size_t)i + k] = q[k] * w;
						L->con[i] |= MF_PREFER;
					}
				}
			}
		}
	}
	/* the coarsest level: the start carried up, or any tangent */
	for (i = 0; i < top->n; i++)
	{
		double* q = top->q + 3 * (size_t)i;
		if (init && (q[0] != 0 || q[1] != 0 || q[2] != 0))
			continue;
		q[0] = 1 + 0.3 * sin(i * 12.9898);
		q[1] = 0.7 * cos(i * 78.233);
		q[2] = 0.5 * sin(i * 39.425);
		mf_tangent(q, top->nrm + 3 * (size_t)i);
	}
	mf_sweep_q(top, 4 * sweeps);
	for (l = (int)h->nl - 2; l >= 0; l--)
	{
		mf_level* L = &h->lv[l];
		for (i = 0; i < L->n; i++)
		{
			double* q = L->q + 3 * (size_t)i;
			memcpy(q, h->lv[l + 1].q + 3 * (size_t)L->parent[i], sizeof(double) * 3);
			mf_tangent(q, L->nrm + 3 * (size_t)i);
		}
		mf_sweep_q(L, sweeps);
	}
}

/* o moved onto its lattice point nearest x, kept on its constraint */
static void mf_settle(mf_level* L, unsigned int i, double* o, double scale)
{
	double const *x = L->x + 3 * (size_t)i, *n = L->nrm + 3 * (size_t)i, *q = L->q + 3 * (size_t)i;
	double t[3], d;
	unsigned int k;
	if (L->con[i] & MF_PINNED)
	{
		memcpy(o, x, sizeof(double) * 3);
		return;
	}
	d = (o[0] - x[0]) * n[0] + (o[1] - x[1]) * n[1] + (o[2] - x[2]) * n[2];
	for (k = 0; k < 3; k++)
		o[k] -= d * n[k];
	if (L->con[i] & MF_ON_LINE)
	{
		double const* c = L->cq + 3 * (size_t)i;
		d = (o[0] - x[0]) * c[0] + (o[1] - x[1]) * c[1] + (o[2] - x[2]) * c[2];
		for (k = 0; k < 3; k++)
			o[k] = x[k] + d * c[k];
	}
	mf_cross(n, q, t);
	mf_lattice(o, q, t, x, L->sc ? scale * L->sc[i] : scale, 1, o);
}

static void mf_sweep_o(mf_level* L, double scale, unsigned int sweeps)
{
	unsigned int s, i, k, q;
	for (s = 0; s < sweeps; s++)
		for (i = 0; i < L->n; i++)
		{
			double const *x = L->x + 3 * (size_t)i, *n = L->nrm + 3 * (size_t)i, *qi = L->q + 3 * (size_t)i;
			double sum[3], ws = 0;
			memcpy(sum, L->o + 3 * (size_t)i, sizeof sum);
			for (k = L->adj0[i]; k < L->adj0[i + 1]; k++)
			{
				unsigned int j = L->adj[k];
				double qa[3], qb[3], a[3], b[3], d;
				mf_compat_orientation(qi, n, L->q + 3 * (size_t)j, L->nrm + 3 * (size_t)j, qa, qb);
				mf_compat_position(x, n, qi, sum, L->x + 3 * (size_t)j, L->nrm + 3 * (size_t)j, qb, L->o + 3 * (size_t)j,
					mf_pair_scale(L, i, j, scale), a, b);
				for (q = 0; q < 3; q++)
					sum[q] = (a[q] * ws + b[q]) / (ws + 1);
				ws += 1;
				d = (sum[0] - x[0]) * n[0] + (sum[1] - x[1]) * n[1] + (sum[2] - x[2]) * n[2];
				for (q = 0; q < 3; q++)
					sum[q] -= d * n[q];
			}
			mf_settle(L, i, sum, scale);
			memcpy(L->o + 3 * (size_t)i, sum, sizeof sum);
		}
}

void mf_positions(mf_hier* h, double scale, unsigned int sweeps)
{
	int l;
	unsigned int i;
	mf_level* top = &h->lv[h->nl - 1];
	for (i = 0; i < top->n; i++)
		memcpy(top->o + 3 * (size_t)i, top->x + 3 * (size_t)i, sizeof(double) * 3);
	mf_sweep_o(top, scale, 4 * sweeps);
	for (l = (int)h->nl - 2; l >= 0; l--)
	{
		mf_level* L = &h->lv[l];
		for (i = 0; i < L->n; i++)
		{
			double* o = L->o + 3 * (size_t)i;
			memcpy(o, h->lv[l + 1].o + 3 * (size_t)L->parent[i], sizeof(double) * 3);
			mf_settle(L, i, o, scale);
		}
		mf_sweep_o(L, scale, sweeps);
	}
}

double mf_pair_scale(mf_level const* L, unsigned int i, unsigned int j, double scale)
{
	return L->sc ? 0.5 * scale * (L->sc[i] + L->sc[j]) : scale;
}

/* per level-0 vertex a value (sc or wt) carried up the levels, each coarse vertex its children's mean by area */
static int mf_carry(mf_hier* h, double const* v0, int weights)
{
	unsigned int l, i, c;
	for (l = 0; l < h->nl; l++)
	{
		mf_level *L = &h->lv[l], *F = l ? &h->lv[l - 1] : NULL;
		double **slot = weights ? &L->wt : &L->sc, *fv = F ? (weights ? F->wt : F->sc) : NULL;
		qaws_internal_free(*slot);
		*slot = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (L->n + 1)));
		if (!*slot)
			return 0;
		for (i = 0; i < L->n; i++)
		{
			double s = 0, a = 0;
			if (!F)
			{
				(*slot)[i] = v0[i];
				continue;
			}
			for (c = L->child0[i]; c < L->child0[i + 1]; c++)
			{
				unsigned int ch = L->child[c];
				double w = F->area[ch] > 0 ? F->area[ch] : 1e-300;
				s += w * fv[ch];
				a += w;
			}
			(*slot)[i] = s / a;
		}
	}
	return 1;
}

int mf_scales(mf_hier* h, double const* s0)
{
	return mf_carry(h, s0, 0);
}

int mf_weights(mf_hier* h, double const* w0)
{
	return mf_carry(h, w0, 1);
}

void mf_attract(double const* x, double const* area, unsigned int n, qaws_scalar const* points, qaws_scalar const* strengths, unsigned int count,
	double r, double* out)
{
	unsigned int i, a, k, pass;
	for (i = 0; i < 2 * n; i++)
		out[i] = i < n ? 1 : 0;
	for (a = 0; a < count; a++)
	{
		double sum = 0, sa = (double)strengths[a];
		for (pass = 0; pass < 2; pass++)
			for (i = 0; i < n; i++)
			{
				double d2 = 0, g;
				for (k = 0; k < 3; k++)
					d2 += (x[3 * (size_t)i + k] - (double)points[3 * (size_t)a + k]) * (x[3 * (size_t)i + k] - (double)points[3 * (size_t)a + k]);
				g = exp(-d2 / (r * r));
				if (!pass)
				{
					out[i] *= exp(-sa * g);
					sum += area[i] * g;
				}
				else if (sum > 0)
					out[n + i] += sa * area[i] * g / sum;
			}
	}
}

void mf_guide(mf_hier* h, double const* dir, double const* w)
{
	mf_level* L0 = &h->lv[0];
	unsigned int l, i, c, k;
	for (i = 0; i < L0->n; i++)
	{
		double* q = L0->cq + 3 * (size_t)i;
		double d = mf_dot(dir + 3 * (size_t)i, L0->nrm + 3 * (size_t)i), t[3];
		if (w[i] == 0 || w[i] != w[i] || (L0->con[i] & MF_ALONG))
			continue;
		for (k = 0; k < 3; k++)
			t[k] = dir[3 * (size_t)i + k] - d * L0->nrm[3 * (size_t)i + k];
		if (mf_unit(t) == 0)
			continue;   /* (along the normal: none) */
		for (k = 0; k < 3; k++)
			q[k] = w[i] > 0 && w[i] < HUGE_VAL ? t[k] * w[i] : t[k];
		L0->con[i] |= w[i] > 0 && w[i] < HUGE_VAL ? MF_PREFER : MF_ALONG;
	}
	/* up: held by a held child (as the lines are), else drawn to its drawn children's directions, as strongly as their mean */
	for (l = 1; l < h->nl; l++)
	{
		mf_level *L = &h->lv[l], *F = &h->lv[l - 1];
		for (i = 0; i < L->n; i++)
		{
			double sum[3] = { 0, 0, 0 }, ref[3], wsum = 0;
			unsigned int n = L->child0[i + 1] - L->child0[i], have = 0;
			if (L->con[i] & MF_ALONG)
				continue;
			for (c = L->child0[i]; c < L->child0[i + 1]; c++)
			{
				unsigned int ch = L->child[c];
				double const* cq = F->cq + 3 * (size_t)ch;
				double a[3], b[3], u[3];
				if (F->con[ch] & MF_ALONG)
				{
					memcpy(L->cq + 3 * (size_t)i, cq, sizeof(double) * 3);
					mf_tangent(L->cq + 3 * (size_t)i, L->nrm + 3 * (size_t)i);
					L->con[i] = (unsigned char)((L->con[i] & ~MF_PREFER) | MF_ALONG);
					break;
				}
				if (!(F->con[ch] & MF_PREFER))
					continue;
				wsum += sqrt(mf_dot(cq, cq)) / n;
				memcpy(u, cq, sizeof u);
				mf_unit(u);
				if (!have)
				{
					memcpy(ref, u, sizeof ref);
					have = 1;
				}
				mf_compat_orientation(ref, L->nrm + 3 * (size_t)i, u, F->nrm + 3 * (size_t)ch, a, b);
				for (k = 0; k < 3; k++)
					sum[k] += b[k];
			}
			if (c < L->child0[i + 1] || !have)
				continue;
			mf_tangent(sum, L->nrm + 3 * (size_t)i);
			for (k = 0; k < 3; k++)
				L->cq[3 * (size_t)i + k] = sum[k] * wsum;
			L->con[i] |= MF_PREFER;
		}
	}
}

/* ======================================================================== */
/*  The public cross field                                                  */
/* ======================================================================== */

qaws_status qaws_mesh_3d_direction_field(qaws_scalar* out_directions, int* out_face_index, qaws_mesh_3d const* mesh,
	qaws_direction_field_desc const* desc)
{
	mf_hier h;
	qaws_hm const* m;
	unsigned int v, f, c, k, n = desc && desc->symmetry ? desc->symmetry : 4;
	int smoothest = desc && desc->method == QAWS_DIRECTION_FIELD_SMOOTHEST;
	int* cone = NULL;
	double* vw = NULL;
	qaws_status st;
	if (!mesh || (desc && (unsigned int)desc->method > QAWS_DIRECTION_FIELD_SMOOTHEST))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (smoothest ? n != 2 && n != 4 && n != 6 : n != 4)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	m = &mesh->hm;
	if (desc && ((desc->cone_count && (!desc->cone_vertices || !desc->cone_indices))
		|| (desc->attractor_count && (!desc->attractors || !desc->attractor_strengths)) || !(desc->attractor_radius >= 0)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (v = 0; desc && v < desc->cone_count; v++)
		if (desc->cone_vertices[v] >= m->nv)
			return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc && desc->cone_count && !smoothest)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (!m->nv)
		return QAWS_STATUS_OK;
	st = mf_build(&h, m, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (desc && (desc->cone_count || desc->attractor_count))
	{
		/* per vertex its cones' indices summed; the attractors' weights */
		double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, r = (double)desc->attractor_radius;
		cone = desc->cone_count ? (int*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(int)) : NULL;
		vw = desc->attractor_count ? (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (m->nv + 1))) : NULL;
		if ((desc->cone_count && !cone) || (desc->attractor_count && !vw))
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (v = 0; v < desc->cone_count; v++)
			cone[desc->cone_vertices[v]] += desc->cone_indices[v];
		for (v = 0; v < m->nv; v++)
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], m->x[3 * (size_t)v + k]);
				hi[k] = fmax(hi[k], m->x[3 * (size_t)v + k]);
			}
		if (r == 0)
			r = 0.1 * sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
		if (vw)
		{
			mf_attract(m->x, h.lv[0].area, m->nv, desc->attractors, desc->attractor_strengths, desc->attractor_count, r > 0 ? r : 1, vw);
			if (!smoothest && !mf_weights(&h, vw))
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
		}
	}
	if (desc && desc->guide_directions)
	{
		/* the guidance on every level (SMOOTHEST takes level 0's: held values and penalties) */
		double *g = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (m->nv + 1))), *w = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (m->nv + 1)));
		if (!g || !w)
		{
			qaws_internal_free(g);
			qaws_internal_free(w);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (v = 0; v < m->nv; v++)
		{
			for (k = 0; k < 3; k++)
				g[3 * (size_t)v + k] = (double)desc->guide_directions[3 * (size_t)v + k];
			w[v] = desc->guide_weights ? (double)desc->guide_weights[v] : 1;
		}
		mf_guide(&h, g, w);
		qaws_internal_free(g);
		qaws_internal_free(w);
	}
	if (smoothest)
	{
		/* one global solve (qaws_cross_field.c), the hierarchy's level 0 for the lines and normals */
		st = cf_direction_field(mesh, &h, n, desc->curvature_aligned, vw ? vw + m->nv : NULL, cone, out_directions, out_face_index);
		goto done;
	}
	mf_orientations(&h, desc && desc->sweeps ? desc->sweeps : 8, NULL, 0);
	for (v = 0; out_directions && v < m->nv; v++)
		for (k = 0; k < 3; k++)
			out_directions[3 * (size_t)v + k] = (qaws_scalar)h.lv[0].q[3 * (size_t)v + k];
	for (f = 0; out_face_index && f < m->nf; f++)
	{
		int r = 0;
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
		{
			unsigned int a = m->cv[c], b = hm_to(m, c);
			/* turns taking b's onto a's, round the face */
			r += mf_rotation(h.lv[0].q + 3 * (size_t)a, h.lv[0].nrm + 3 * (size_t)a, h.lv[0].q + 3 * (size_t)b, h.lv[0].nrm + 3 * (size_t)b);
		}
		r %= 4;
		out_face_index[f] = r == 3 ? -1 : r;
	}
done:
	qaws_internal_free(cone);
	qaws_internal_free(vw);
	mf_free(&h);
	return st;
}
