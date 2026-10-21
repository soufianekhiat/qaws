/*
 * The mesh connector (qaws_mesh_connect.h): borders walked on the corner
 * store, matched by arc length, rings Hermite-blended between them along
 * the meshes' surfaces, zipped into faces, the result built anew with the
 * attributes and features carried; planar meshes through their surfaces;
 * volume patches matched through their Tutte embeddings and joined by
 * layers of prisms and hexahedra.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_connect.h"
#include "../qaws_mesh_smooth.h"
#include "../qaws_quality_report.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CN_NONE 0xFFFFFFFFu

/* a border: its vertices (in the border's direction, or a chain's), positions, the directions leaving the mesh across it, its
   arc length's share at each vertex (n + 1 for a loop: the last 1) */
typedef struct cn_loop
{
	unsigned int n;
	int closed;
	unsigned int* v;
	double *x, *o, *s, h;           /* h: the mean edge length */
} cn_loop;

static void cn_loop_free(cn_loop* L)
{
	qaws_internal_free(L->v);
	qaws_internal_free(L->x);
	qaws_internal_free(L->o);
	qaws_internal_free(L->s);
	memset(L, 0, sizeof *L);
}

static void cn_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double cn_norm(double* a)
{
	double l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	if (l > 0)
	{
		a[0] /= l;
		a[1] /= l;
		a[2] /= l;
	}
	return l;
}

static double cn_dist(double const* a, double const* b)
{
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

/* the border's shares by arc length and its mean edge (L->x set) */
static void cn_shares(cn_loop* L)
{
	unsigned int i, segs = L->closed ? L->n : L->n - 1;
	double len = 0;
	L->s[0] = 0;
	for (i = 0; i < segs; i++)
	{
		len += cn_dist(L->x + 3 * (size_t)i, L->x + 3 * (size_t)((i + 1) % L->n));
		L->s[i + 1] = len;
	}
	for (i = 0; i <= segs; i++)
		L->s[i] = len > 0 ? L->s[i] / len : 0;
	L->h = segs ? len / segs : 0;
}

/* the border's geometry from its vertices: positions, shares, the directions leaving the mesh (across the border in the vertex's
   tangent plane, away from the face at the vertex) */
static qaws_status cn_geom(cn_loop* L, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m = &mesh->hm;
	qaws_scalar* nrm = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)m->nv + 1));
	unsigned int i, k;
	L->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, L->n));
	L->o = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, L->n));
	L->s = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)L->n + 1));
	if (!nrm || !L->x || !L->o || !L->s || qaws_mesh_3d_get_vertex_normals(nrm, mesh) != QAWS_STATUS_OK)
	{
		qaws_internal_free(nrm);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < L->n; i++)
		memcpy(L->x + 3 * (size_t)i, m->x + 3 * (size_t)L->v[i], 3 * sizeof(double));
	cn_shares(L);
	for (i = 0; i < L->n; i++)
	{
		unsigned int ip = L->closed ? (i + L->n - 1) % L->n : (i ? i - 1 : 0), in = L->closed ? (i + 1) % L->n : (i + 1 < L->n ? i + 1 : i);
		double const *pp = L->x + 3 * (size_t)ip, *pn = L->x + 3 * (size_t)in, *x = L->x + 3 * (size_t)i;
		double t[3] = { pn[0] - pp[0], pn[1] - pp[1], pn[2] - pp[2] }, n[3], cen[3] = { 0, 0, 0 }, *o = L->o + 3 * (size_t)i;
		unsigned int f = m->cf[m->vc[L->v[i]]], c;
		for (k = 0; k < 3; k++)
			n[k] = (double)nrm[3 * (size_t)L->v[i] + k];
		cn_norm(t);
		cn_cross(t, n, o);
		cn_norm(o);
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			for (k = 0; k < 3; k++)
				cen[k] += m->x[3 * (size_t)m->cv[c] + k] / m->fn[f];
		if (o[0] * (cen[0] - x[0]) + o[1] * (cen[1] - x[1]) + o[2] * (cen[2] - x[2]) > 0)
			for (k = 0; k < 3; k++)
				o[k] = -o[k];
	}
	qaws_internal_free(nrm);
	return QAWS_STATUS_OK;
}

/* the border walk from v0 along the border's own corners: to stop (a chain, stop included) or back to v0 (a loop, stop CN_NONE);
   INVALID_ARGUMENT when v0 is not on a border or the walk misses its end */
static qaws_status cn_walk(cn_loop* L, qaws_mesh_3d const* mesh, unsigned int v0, unsigned int stop)
{
	qaws_hm const* m = &mesh->hm;
	unsigned int v;
	memset(L, 0, sizeof *L);
	L->closed = stop == CN_NONE;
	if (v0 >= m->nv || !hm_is_border_vertex(m, v0) || (!L->closed && (stop >= m->nv || stop == v0)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	L->v = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)m->nv + 1));
	if (!L->v)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	v = v0;
	for (;;)
	{
		unsigned int c = m->vc[v];
		if (L->n >= m->nv || c == HM_NONE || m->co[c] != HM_NONE)
		{
			cn_loop_free(L);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
		L->v[L->n++] = v;
		if (!L->closed && v == stop)
			break;
		v = hm_to(m, c);
		if (v == v0)
		{
			if (L->closed)
				break;
			cn_loop_free(L);   /* (around the loop without meeting the chain's end) */
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
	}
	if (L->n < (L->closed ? 3u : 2u))
	{
		cn_loop_free(L);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	return cn_geom(L, mesh);
}

/* the border reversed (a loop also started at its vertex j0), the shares recomputed */
static qaws_status cn_turn(cn_loop* L, unsigned int j0)
{
	unsigned int i, n = L->n;
	unsigned int* v = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * n);
	double* x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, n));
	double* o = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, n));
	if (!v || !x || !o)
	{
		qaws_internal_free(v); qaws_internal_free(x); qaws_internal_free(o);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < n; i++)
	{
		unsigned int j = L->closed ? (j0 + n - i) % n : n - 1 - i;
		v[i] = L->v[j];
		memcpy(x + 3 * (size_t)i, L->x + 3 * (size_t)j, 3 * sizeof(double));
		memcpy(o + 3 * (size_t)i, L->o + 3 * (size_t)j, 3 * sizeof(double));
	}
	qaws_internal_free(L->v); qaws_internal_free(L->x); qaws_internal_free(L->o);
	L->v = v;
	L->x = x;
	L->o = o;
	cn_shares(L);
	return QAWS_STATUS_OK;
}

/* the segment and fraction at share s */
static unsigned int cn_seg(cn_loop const* L, double s, double* t)
{
	unsigned int lo = 0, hi = L->closed ? L->n : L->n - 1;
	if (L->closed)
		s -= floor(s);
	else
		s = s < 0 ? 0 : s > 1 ? 1 : s;
	while (hi - lo > 1)
	{
		unsigned int mid = (lo + hi) / 2;
		if (L->s[mid] <= s)
			lo = mid;
		else
			hi = mid;
	}
	*t = L->s[lo + 1] > L->s[lo] ? (s - L->s[lo]) / (L->s[lo + 1] - L->s[lo]) : 0;
	*t = *t < 0 ? 0 : *t > 1 ? 1 : *t;
	return lo;
}

/* per-vertex values (dim each) along the border at share s, linear between its vertices */
static void cn_lerp(cn_loop const* L, double const* val, unsigned int dim, double s, double* out)
{
	double t;
	unsigned int lo = cn_seg(L, s, &t), hi = (lo + 1) % L->n, k;
	for (k = 0; k < dim; k++)
		out[k] = val[(size_t)dim * lo + k] + t * (val[(size_t)dim * hi + k] - val[(size_t)dim * lo + k]);
}

static void cn_at(cn_loop const* L, double s, double* x, double* o)
{
	cn_lerp(L, L->x, 3, s, x);
	if (o)
	{
		cn_lerp(L, L->o, 3, s, o);
		cn_norm(o);
	}
}

static void cn_hermite(double t, double* h)
{
	h[0] = 2 * t * t * t - 3 * t * t + 1;
	h[1] = t * t * t - 2 * t * t + t;
	h[2] = -2 * t * t * t + 3 * t * t;
	h[3] = t * t * t - t * t;
}

/* the position of corner c of hm in the faces' order (qaws_mesh_3d_get_faces): its face's first plus its rank */
static unsigned int cn_corner_rank(qaws_hm const* m, unsigned int const* first, unsigned int c)
{
	return first[m->cf[c]] + (c - m->fc[m->cf[c]]);
}

/* the attributes a and b share carried onto the result: a's and b's elements copied, the bridge's blended */
typedef struct cn_carry
{
	qaws_mesh_3d const *a, *b;
	int same;
	cn_loop const *A, *B;
	unsigned int nva, nvb, nfa, nfb, ca, cb, nbv, nf;
	double const *bs, *bt;          /* per bridge vertex: share and Hermite weight of b */
	unsigned int const *fv, *fs;    /* the result's faces */
	unsigned int const *frow;       /* per bridge face: nearer b (1) or a (0) */
	double const* fsh;              /* per bridge face: its share */
} cn_carry;

static qaws_status cn_carry_attributes(cn_carry const* C, qaws_mesh_3d* out)
{
	unsigned int na = qaws_mesh_3d_get_attribute_count(C->a), i, k, j;
	unsigned int *fa = NULL, *fb = NULL;
	qaws_status st = QAWS_STATUS_OK;
	fa = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)C->a->hm.nf + 1));
	fb = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)C->b->hm.nf + 1));
	if (!fa || !fb)
	{
		qaws_internal_free(fa);
		qaws_internal_free(fb);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* (each face's first corner in the faces' order) */
	for (i = 0, k = 0; i < C->a->hm.nf; i++)
	{
		fa[i] = k;
		k += C->a->hm.fn[i];
	}
	for (i = 0, k = 0; i < C->b->hm.nf; i++)
	{
		fb[i] = k;
		k += C->b->hm.fn[i];
	}
	for (i = 0; i < na && st == QAWS_STATUS_OK; i++)
	{
		qaws_attribute_desc da, db;
		unsigned int ib = i, dim, nel, ela, elb, at;
		qaws_scalar *va = NULL, *vb = NULL, *vo = NULL;
		double *la = NULL, *lb = NULL;
		if (qaws_mesh_3d_get_attribute_desc(&da, C->a, i) != QAWS_STATUS_OK)
			continue;
		if (!C->same)
		{
			if (!da.name || qaws_mesh_3d_find_attribute(&ib, C->b, da.name) != QAWS_STATUS_OK
				|| qaws_mesh_3d_get_attribute_desc(&db, C->b, ib) != QAWS_STATUS_OK || db.element != da.element || db.components != da.components)
				continue;
		}
		dim = da.components;
		ela = da.element == QAWS_MESH_VERTEX ? C->nva : da.element == QAWS_MESH_CORNER ? C->ca : C->nfa;
		elb = C->same ? 0 : da.element == QAWS_MESH_VERTEX ? C->nvb : da.element == QAWS_MESH_CORNER ? C->cb : C->nfb;
		nel = da.element == QAWS_MESH_VERTEX ? C->nva + C->nvb + C->nbv : da.element == QAWS_MESH_CORNER ? C->ca + C->cb : C->nf;
		if (da.element == QAWS_MESH_CORNER)
			for (j = C->nfa + C->nfb, nel = C->ca + C->cb; j < C->nf; j++)
				nel += C->fs[j];
		va = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)dim * ela + 1));
		vb = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)dim * elb + 1));
		vo = (qaws_scalar*)qaws_internal_zalloc((size_t)dim * nel + 1, sizeof(qaws_scalar));
		la = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), dim, (size_t)C->A->n + 1));
		lb = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), dim, (size_t)C->B->n + 1));
		if (!va || !vb || !vo || !la || !lb || qaws_mesh_3d_get_attribute(va, C->a, i) != QAWS_STATUS_OK
			|| (!C->same && qaws_mesh_3d_get_attribute(vb, C->b, ib) != QAWS_STATUS_OK))
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		if (st == QAWS_STATUS_OK)
		{
			qaws_scalar const* srcb = C->same ? va : vb;
			unsigned int const* fbb = C->same ? fa : fb;
			qaws_mesh_3d const* mb = C->same ? C->a : C->b;
			memcpy(vo, va, sizeof(qaws_scalar) * dim * ela);
			memcpy(vo + (size_t)dim * ela, vb, sizeof(qaws_scalar) * dim * elb);
			/* the borders' values per their vertices: a vertex's own, or its border corner's */
			for (j = 0; j < C->A->n; j++)
			{
				unsigned int v = C->A->v[j], e = da.element == QAWS_MESH_VERTEX ? v : da.element == QAWS_MESH_CORNER ? cn_corner_rank(&C->a->hm, fa, C->a->hm.vc[v])
					: C->a->hm.cf[C->a->hm.vc[v]];
				for (k = 0; k < dim; k++)
					la[(size_t)dim * j + k] = (double)va[(size_t)dim * e + k];
			}
			for (j = 0; j < C->B->n; j++)
			{
				unsigned int v = C->B->v[j], e = da.element == QAWS_MESH_VERTEX ? v : da.element == QAWS_MESH_CORNER ? cn_corner_rank(&mb->hm, fbb, mb->hm.vc[v])
					: mb->hm.cf[mb->hm.vc[v]];
				for (k = 0; k < dim; k++)
					lb[(size_t)dim * j + k] = (double)srcb[(size_t)dim * e + k];
			}
			if (da.element == QAWS_MESH_FACE)
				for (j = C->nfa + C->nfb; j < C->nf; j++)
				{
					/* a bridge face: the nearer border's face at its share */
					double t;
					double const* lv = C->frow[j - C->nfa - C->nfb] ? lb : la;
					cn_loop const* L = C->frow[j - C->nfa - C->nfb] ? C->B : C->A;
					unsigned int seg = cn_seg(L, C->fsh[j - C->nfa - C->nfb], &t), q = t < 0.5 ? seg : (seg + 1) % L->n;
					for (k = 0; k < dim; k++)
						vo[(size_t)dim * j + k] = (qaws_scalar)lv[(size_t)dim * q + k];
				}
			else
			{
				/* a bridge vertex (or a bridge face's corner at one): the borders' values at its share, blended by its weight */
				unsigned int nbase = da.element == QAWS_MESH_VERTEX ? C->nva + C->nvb : C->ca + C->cb, cnt = da.element == QAWS_MESH_VERTEX ? C->nbv : nel - nbase;
				for (j = 0; j < cnt; j++)
				{
					unsigned int v = da.element == QAWS_MESH_VERTEX ? nbase + j : C->fv[nbase + j];
					double xa[16], xb[16], w, s, l2 = 0;
					qaws_scalar* o = vo + (size_t)dim * (nbase + j);
					if (v >= C->nva + C->nvb)
					{
						s = C->bs[v - C->nva - C->nvb];
						w = C->bt[v - C->nva - C->nvb];
					}
					else
					{
						/* (a border vertex: its own border's value) */
						unsigned int q;
						int onb = 0;
						for (q = 0; q < C->B->n && !onb; q++)
							onb = (C->same ? C->B->v[q] : C->nva + C->B->v[q]) == v;
						s = 0;
						w = onb ? 1 : 0;
						for (q = 0; q < (onb ? C->B->n : C->A->n); q++)
							if ((onb ? (C->same ? C->B->v[q] : C->nva + C->B->v[q]) : C->A->v[q]) == v)
								s = onb ? C->B->s[q] : C->A->s[q];
					}
					cn_lerp(C->A, la, dim, s, xa);
					cn_lerp(C->B, lb, dim, s, xb);
					for (k = 0; k < dim; k++)
					{
						double y = da.rule == QAWS_ATTRIBUTE_NEAREST ? (w < 0.5 ? xa[k] : xb[k]) : (1 - w) * xa[k] + w * xb[k];
						xa[k] = y;
						l2 += y * y;
					}
					for (k = 0; k < dim; k++)
						o[k] = (qaws_scalar)(da.rule == QAWS_ATTRIBUTE_NORMALIZED && l2 > 0 ? xa[k] / sqrt(l2) : xa[k]);
				}
			}
			st = qaws_mesh_3d_add_attribute(&at, out, vo, &da);
		}
		qaws_internal_free(va); qaws_internal_free(vb); qaws_internal_free(vo); qaws_internal_free(la); qaws_internal_free(lb);
	}
	qaws_internal_free(fa);
	qaws_internal_free(fb);
	return st;
}

/* the bridge's faces measured (qaws_quality_report.h): smallest and largest angle, smallest scaled Jacobian */
static void cn_measure(qaws_connect_report* r, qaws_mesh_3d const* out, unsigned int f0, unsigned int nf)
{
	qaws_scalar* v = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * ((size_t)nf + 1));
	unsigned int f, k;
	if (!v)
		return;
	for (k = 0; k < 3; k++)
	{
		qaws_element_measure ms = k == 0 ? QAWS_MEASURE_MIN_ANGLE : k == 1 ? QAWS_MEASURE_MAX_ANGLE : QAWS_MEASURE_SCALED_JACOBIAN;
		double best = k == 1 ? -HUGE_VAL : HUGE_VAL;
		if (qaws_mesh_3d_measure_elements(v, out, ms) != QAWS_STATUS_OK)
			continue;
		for (f = f0; f < nf; f++)
			if (v[f] == v[f])
				best = k == 1 ? (v[f] > best ? v[f] : best) : (v[f] < best ? v[f] : best);
		if (k == 0)
			r->min_angle = (qaws_scalar)(best * 3.14159265358979323846 / 180);
		else if (k == 1)
			r->max_angle = (qaws_scalar)(best * 3.14159265358979323846 / 180);
		else
			r->min_scaled_jacobian = (qaws_scalar)best;
	}
	qaws_internal_free(v);
}

qaws_status qaws_mesh_3d_connect(qaws_mesh_3d** out, qaws_connect_report* report, qaws_mesh_3d const* a, unsigned int av,
	qaws_mesh_3d const* b, unsigned int bv, qaws_connect_desc const* desc)
{
	qaws_connect_report local;
	cn_loop A, B;
	unsigned int nva, nvb, R, r, i, j, k, nfa = 0, nfb = 0, ca = 0, cb = 0, nbv = 0, nbf = 0, nbc = 0, *ring_start = NULL, *ring_n = NULL,
		*fv = NULL, *fs = NULL, *afs = NULL, *bfs = NULL, *frow = NULL;
	double *bx = NULL, *bs = NULL, *bt = NULL, *fsh = NULL, gap = 0, *t = NULL;
	qaws_scalar* pos = NULL;
	int same = !b || b == a, chains = desc && (desc->flags & QAWS_CONNECT_CHAINS);
	qaws_status st;
	memset(&A, 0, sizeof A);
	memset(&B, 0, sizeof B);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!a)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (same)
		b = a;
	st = cn_walk(&A, a, av, chains ? desc->a_end : CN_NONE);
	if (st == QAWS_STATUS_OK)
	{
		/* b's chain runs against its border: walked from its end back to b_vertex, then turned */
		st = chains ? cn_walk(&B, b, desc->b_end, bv) : cn_walk(&B, b, bv, CN_NONE);
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	if (same)
		for (i = 0; i < A.n; i++)
			for (j = 0; j < B.n; j++)
				if (A.v[i] == B.v[j])
				{
					st = QAWS_STATUS_INVALID_ARGUMENT;   /* (one border named twice) */
					goto done;
				}
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (chains)
	{
		if (cn_turn(&B, 0) != QAWS_STATUS_OK)
			goto done;
	}
	else
	{
		/* b's loop turned the other way, started where the matched points are nearest */
		unsigned int best = 0;
		double bc = HUGE_VAL;
		for (j = 0; j < B.n; j++)
		{
			cn_loop T;
			double c = 0;
			memset(&T, 0, sizeof T);
			T.n = B.n;
			T.closed = 1;
			T.v = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * B.n);
			T.x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, B.n));
			T.o = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, B.n));
			T.s = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)B.n + 1));
			if (!T.v || !T.x || !T.o || !T.s)
			{
				cn_loop_free(&T);
				goto done;
			}
			memcpy(T.v, B.v, sizeof(unsigned int) * B.n);
			memcpy(T.x, B.x, sizeof(double) * 3 * B.n);
			memcpy(T.o, B.o, sizeof(double) * 3 * B.n);
			if (cn_turn(&T, j) != QAWS_STATUS_OK)
			{
				cn_loop_free(&T);
				goto done;
			}
			for (i = 0; i < A.n; i++)
			{
				double y[3];
				cn_at(&T, A.s[i], y, NULL);
				c += cn_dist(y, A.x + 3 * (size_t)i) * cn_dist(y, A.x + 3 * (size_t)i);
			}
			cn_loop_free(&T);
			if (c < bc)
			{
				bc = c;
				best = j;
			}
		}
		if (cn_turn(&B, best) != QAWS_STATUS_OK)
			goto done;
	}
	for (i = 0; i < A.n; i++)
	{
		double y[3];
		cn_at(&B, A.s[i], y, NULL);
		gap += cn_dist(y, A.x + 3 * (size_t)i);
	}
	gap /= A.n;
	R = desc && desc->rows ? desc->rows : (unsigned int)fmax(1, floor(gap / (0.5 * (A.h + B.h)) + 0.5));
	report->rows = R;
	report->loop_vertices[0] = A.n;
	report->loop_vertices[1] = B.n;
	report->gap = (qaws_scalar)gap;
	/* the rings: 0 is a's border, R is b's; between, counts from one to the other, spacing graded from a's edges to b's */
	ring_start = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)R + 2));
	ring_n = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)R + 2));
	t = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)R + 2));
	if (!ring_start || !ring_n || !t)
		goto done;
	{
		double tot = 0;
		t[0] = 0;
		for (r = 0; r < R; r++)
		{
			tot += A.h + (B.h - A.h) * (r + 0.5) / R;
			t[r + 1] = tot;
		}
		for (r = 0; r <= R; r++)
			t[r] /= tot;
	}
	nva = a->hm.nv;
	nvb = same ? 0 : b->hm.nv;
	for (r = 0; r <= R; r++)
	{
		ring_n[r] = r == 0 ? A.n : r == R ? B.n : (unsigned int)floor(A.n + ((double)B.n - A.n) * r / R + 0.5);
		if (r > 0 && r < R)
		{
			ring_start[r] = nva + nvb + nbv;
			nbv += ring_n[r];
		}
	}
	bx = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nbv + 1));
	bs = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nbv + 1));
	bt = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nbv + 1));
	if (!bx || !bs || !bt)
		goto done;
	for (r = 1, k = 0; r < R; r++)
		for (i = 0; i < ring_n[r]; i++, k++)
		{
			/* the Hermite blend at t[r] between the borders' points at share s, leaving a and arriving into b along their surfaces */
			double s = chains ? (double)i / (ring_n[r] - 1) : (double)i / ring_n[r], pa[3], oa[3], pb[3], ob[3], L, h[4];
			unsigned int q;
			cn_at(&A, s, pa, oa);
			cn_at(&B, s, pb, ob);
			L = cn_dist(pa, pb);
			cn_hermite(t[r], h);
			for (q = 0; q < 3; q++)
				bx[3 * (size_t)k + q] = h[0] * pa[q] + h[1] * L * oa[q] + h[2] * pb[q] - h[3] * L * ob[q];
			bs[k] = s;
			bt[k] = h[2];
		}
	/* faces: a's, b's, then the bridge's rows */
	nfa = qaws_mesh_3d_get_face_count(a);
	ca = qaws_mesh_3d_get_corner_count(a);
	if (!same)
	{
		nfb = qaws_mesh_3d_get_face_count(b);
		cb = qaws_mesh_3d_get_corner_count(b);
	}
	for (r = 0; r < R; r++)
	{
		unsigned int sl = chains ? ring_n[r] - 1 : ring_n[r], su = chains ? ring_n[r + 1] - 1 : ring_n[r + 1];
		int quads = desc && desc->quads && ring_n[r] == ring_n[r + 1];
		nbf += quads ? sl : sl + su;
		nbc += quads ? 4 * sl : 3 * (sl + su);
	}
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nva + nvb + nbv + 1));
	fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)ca + cb + nbc + 1));
	fs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nfa + nfb + nbf + 1));
	afs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nfa + 1));
	frow = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nbf + 1));
	fsh = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nbf + 1));
	if (!pos || !fv || !fs || !afs || !frow || !fsh)
		goto done;
	if (qaws_mesh_3d_get_positions(pos, a) != QAWS_STATUS_OK || qaws_mesh_3d_get_faces(fv, afs, a) != QAWS_STATUS_OK)
		goto done;
	memcpy(fs, afs, sizeof(unsigned int) * nfa);
	if (!same)
	{
		bfs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nfb + 1));
		if (!bfs || qaws_mesh_3d_get_positions(pos + 3 * (size_t)nva, b) != QAWS_STATUS_OK || qaws_mesh_3d_get_faces(fv + ca, bfs, b) != QAWS_STATUS_OK)
			goto done;
		for (i = 0; i < cb; i++)
			fv[ca + i] += nva;
		memcpy(fs + nfa, bfs, sizeof(unsigned int) * nfb);
	}
	for (i = 0; i < 3 * nbv; i++)
		pos[3 * ((size_t)nva + nvb) + i] = (qaws_scalar)bx[i];
	{
		unsigned int nf = nfa + nfb, nc = ca + cb;
		for (r = 0; r < R; r++)
		{
			/* ring r's vertices (lower, nearer a) and ring r + 1's (upper), with their shares; lower edges run as a's border does */
			unsigned int nl = ring_n[r], nu = ring_n[r + 1], il = 0, iu = 0, sl = chains ? nl - 1 : nl, su = chains ? nu - 1 : nu;
			unsigned int *lv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * nl), *uv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * nu);
			double *ls = (double*)qaws_internal_malloc(sizeof(double) * (nl + 1)), *us = (double*)qaws_internal_malloc(sizeof(double) * (nu + 1));
			int quads = desc && desc->quads && nl == nu, upper = 2 * r + 1 >= R;
			if (!lv || !uv || !ls || !us)
			{
				qaws_internal_free(lv); qaws_internal_free(uv); qaws_internal_free(ls); qaws_internal_free(us);
				goto done;
			}
			for (i = 0; i < nl; i++)
			{
				lv[i] = r == 0 ? A.v[i] : ring_start[r] + i;
				ls[i] = r == 0 ? A.s[i] : chains ? (double)i / (nl - 1) : (double)i / nl;
			}
			for (i = 0; i < nu; i++)
			{
				uv[i] = r + 1 == R ? (same ? B.v[i] : nva + B.v[i]) : ring_start[r + 1] + i;
				us[i] = r + 1 == R ? B.s[i] : chains ? (double)i / (nu - 1) : (double)i / nu;
			}
			ls[nl] = 1;
			us[nu] = 1;
			if (quads)
				for (i = 0; i < sl; i++)
				{
					fv[nc++] = lv[(i + 1) % nl];
					fv[nc++] = lv[i];
					fv[nc++] = uv[i];
					fv[nc++] = uv[(i + 1) % nu];
					frow[nf - nfa - nfb] = (unsigned int)upper;
					fsh[nf - nfa - nfb] = 0.5 * (ls[i] + ls[i + 1]);
					fs[nf++] = 4;
				}
			else
				while (il < sl || iu < su)
				{
					/* (the ring whose next vertex comes first by arc length advances) */
					double pl = il < sl ? ls[il + 1] : HUGE_VAL, pu = iu < su ? us[iu + 1] : HUGE_VAL;
					if (pl <= pu)
					{
						fv[nc++] = lv[(il + 1) % nl];
						fv[nc++] = lv[il % nl];
						fv[nc++] = uv[iu % nu];
						fsh[nf - nfa - nfb] = 0.5 * (ls[il] + ls[il + 1]);
						il++;
					}
					else
					{
						fv[nc++] = uv[iu % nu];
						fv[nc++] = uv[(iu + 1) % nu];
						fv[nc++] = lv[il % nl];
						fsh[nf - nfa - nfb] = 0.5 * (us[iu] + us[iu + 1]);
						iu++;
					}
					frow[nf - nfa - nfb] = (unsigned int)upper;
					fs[nf++] = 3;
				}
			qaws_internal_free(lv); qaws_internal_free(uv); qaws_internal_free(ls); qaws_internal_free(us);
		}
		report->bridge_faces = nf - nfa - nfb;
		{
			qaws_mesh_3d_desc d;
			memset(&d, 0, sizeof d);
			d.positions = pos;
			d.vertex_count = nva + nvb + nbv;
			d.face_vertices = fv;
			d.face_sizes = fs;
			d.face_count = nf;
			st = qaws_mesh_3d_create(out, &d);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
		{
			/* the features, the attributes, the faces' parts */
			unsigned int nea = qaws_mesh_3d_get_feature_edge_count(a), neb = same ? 0 : qaws_mesh_3d_get_feature_edge_count(b), at = 0, f;
			unsigned int* fe = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nea + neb + 1));
			qaws_scalar* part = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * ((size_t)nf + 1));
			qaws_attribute_desc ad;
			cn_carry C;
			if (!fe || !part)
			{
				qaws_internal_free(fe);
				qaws_internal_free(part);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			st = qaws_mesh_3d_get_feature_edges(fe, a);
			if (st == QAWS_STATUS_OK && neb)
			{
				st = qaws_mesh_3d_get_feature_edges(fe + 2 * (size_t)nea, b);
				for (i = 0; i < 2 * neb; i++)
					fe[2 * (size_t)nea + i] += nva;
			}
			if (st == QAWS_STATUS_OK && nea + neb)
				st = qaws_mesh_3d_set_feature_edges(*out, fe, nea + neb);
			qaws_internal_free(fe);
			memset(&C, 0, sizeof C);
			C.a = a;
			C.b = b;
			C.same = same;
			C.A = &A;
			C.B = &B;
			C.nva = nva;
			C.nvb = nvb;
			C.nfa = nfa;
			C.nfb = nfb;
			C.ca = ca;
			C.cb = cb;
			C.nbv = nbv;
			C.nf = nf;
			C.bs = bs;
			C.bt = bt;
			C.fv = fv;
			C.fs = fs;
			C.frow = frow;
			C.fsh = fsh;
			if (st == QAWS_STATUS_OK)
				st = cn_carry_attributes(&C, *out);
			for (f = 0; f < nf; f++)
				part[f] = (qaws_scalar)(f < nfa ? 0 : f < nfa + nfb ? 2 : 1);
			memset(&ad, 0, sizeof ad);
			ad.name = "part";
			ad.element = QAWS_MESH_FACE;
			ad.components = 1;
			ad.rule = QAWS_ATTRIBUTE_NEAREST;
			if (st == QAWS_STATUS_OK)
				st = qaws_mesh_3d_add_attribute(&at, *out, part, &ad);
			qaws_internal_free(part);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
		if (desc && (desc->flags & QAWS_CONNECT_FAIR) && nbv)
		{
			/* (a smoothing flow, not a thin plate: spanning two borders far apart, a thin plate pinches like a soap film; the flow
			   evens the bridge out and leaves its shape) */
			unsigned char* held = (unsigned char*)qaws_internal_malloc((size_t)nva + nvb + nbv);
			qaws_smooth_desc sd;
			if (!held)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			memset(held, 1, (size_t)nva + nvb);
			memset(held + nva + nvb, 0, nbv);
			memset(&sd, 0, sizeof sd);
			sd.method = QAWS_SMOOTH_BILAPLACIAN;
			sd.iterations = 3;
			sd.held = held;
			st = qaws_mesh_3d_smooth(*out, &sd);
			qaws_internal_free(held);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
		cn_measure(report, *out, nfa + nfb, nf);
	}
	st = QAWS_STATUS_OK;
done:
	if (st != QAWS_STATUS_OK && out && *out)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
	}
	cn_loop_free(&A);
	cn_loop_free(&B);
	qaws_internal_free(ring_start); qaws_internal_free(ring_n); qaws_internal_free(t); qaws_internal_free(bx); qaws_internal_free(bs);
	qaws_internal_free(bt); qaws_internal_free(pos); qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(afs);
	qaws_internal_free(bfs); qaws_internal_free(frow); qaws_internal_free(fsh);
	return st;
}

qaws_status qaws_mesh_2d_connect(qaws_mesh_2d** out, qaws_connect_report* report, qaws_mesh_2d const* a, unsigned int av,
	qaws_mesh_2d const* b, unsigned int bv, qaws_connect_desc const* desc)
{
	qaws_mesh_3d* m = NULL;
	qaws_scalar *p3 = NULL, *p2 = NULL;
	unsigned int *fv = NULL, *fs = NULL, nv, nf, i;
	qaws_mesh_2d_desc d;
	qaws_status st;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!a)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* (the planar meshes' surfaces at z = 0: the bridge stays in the plane) */
	st = qaws_mesh_3d_connect(&m, report, qaws_mesh_2d_get_mesh(a), av, b ? qaws_mesh_2d_get_mesh(b) : NULL, bv, desc);
	if (st != QAWS_STATUS_OK)
		return st;
	nv = qaws_mesh_3d_get_vertex_count(m);
	nf = qaws_mesh_3d_get_face_count(m);
	p3 = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	p2 = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 2, (size_t)nv + 1));
	fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)qaws_mesh_3d_get_corner_count(m) + 1));
	fs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (p3 && p2 && fv && fs && qaws_mesh_3d_get_positions(p3, m) == QAWS_STATUS_OK && qaws_mesh_3d_get_faces(fv, fs, m) == QAWS_STATUS_OK)
	{
		for (i = 0; i < nv; i++)
		{
			p2[2 * (size_t)i] = p3[3 * (size_t)i];
			p2[2 * (size_t)i + 1] = p3[3 * (size_t)i + 1];
		}
		memset(&d, 0, sizeof d);
		d.positions = p2;
		d.vertex_count = nv;
		d.face_vertices = fv;
		d.face_sizes = fs;
		d.face_count = nf;
		st = qaws_mesh_2d_create(out, &d);
	}
	qaws_mesh_3d_destroy(m);
	qaws_internal_free(p3); qaws_internal_free(p2); qaws_internal_free(fv); qaws_internal_free(fs);
	return st;
}

/* ---- volumes ------------------------------------------------------------ */

/* the local facets of each cell type (qaws_volume_mesh.h's numbering), at most 4 corners each, -1 past the end */
static int const cn_tet_f[4][4] = { { 1, 2, 3, -1 }, { 0, 3, 2, -1 }, { 0, 1, 3, -1 }, { 0, 2, 1, -1 } };
static int const cn_hex_f[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
static int const cn_pri_f[5][4] = { { 0, 2, 1, -1 }, { 3, 4, 5, -1 }, { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 2, 0, 3, 5 } };
static int const cn_pyr_f[5][4] = { { 0, 3, 2, 1 }, { 0, 1, 4, -1 }, { 1, 2, 4, -1 }, { 2, 3, 4, -1 }, { 3, 0, 4, -1 } };

/* a volume patch: its facets' vertices (volume indices, outward), its vertices, positions, normals, border loop, embedding */
typedef struct cn_patch
{
	unsigned int nf, nv, nb;
	unsigned int *fsz, *fvv;        /* per facet its size; its vertices (local indices), 4 per facet */
	unsigned int* vol;              /* per local vertex its volume vertex */
	double *x, *n, *e;              /* positions, outward normals, Tutte embedding (2 each) */
	unsigned int* loop;             /* the border loop (local indices), in the facets' direction */
	double h;                       /* mean edge */
} cn_patch;

static void cn_patch_free(cn_patch* P)
{
	qaws_internal_free(P->fsz); qaws_internal_free(P->fvv); qaws_internal_free(P->vol); qaws_internal_free(P->x); qaws_internal_free(P->n);
	qaws_internal_free(P->e); qaws_internal_free(P->loop);
	memset(P, 0, sizeof *P);
}

static int cn_cmp_u64(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

/* the patch from the given border facets; INVALID_ARGUMENT for a facet that is not on the border or a patch that is not a disk */
static qaws_status cn_patch_build(cn_patch* P, qaws_volume_mesh const* vm, unsigned int const* facets, unsigned int count)
{
	unsigned int nc = qaws_volume_mesh_get_cell_count(vm), nvv = qaws_volume_mesh_get_vertex_count(vm), nfac = qaws_volume_mesh_get_facet_count(vm);
	unsigned int *cv = NULL, *adj = NULL, *cfirst = NULL, *ffirst = NULL, *map = NULL, i, k, c, ne = 0;
	qaws_cell_type* ty = NULL;
	qaws_scalar* pos = NULL;
	unsigned long long* he = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(P, 0, sizeof *P);
	if (!facets || !count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	cv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)qaws_volume_mesh_get_corner_count(vm) + 1));
	ty = (qaws_cell_type*)qaws_internal_malloc(sizeof(qaws_cell_type) * ((size_t)nc + 1));
	adj = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nfac + 1));
	cfirst = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nc + 1));
	ffirst = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nc + 2));
	map = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nvv + 1));
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nvv + 1));
	P->fsz = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * count);
	P->fvv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), 4 * (size_t)count));
	P->vol = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (4 * (size_t)count + 1));
	if (!cv || !ty || !adj || !cfirst || !ffirst || !map || !pos || !P->fsz || !P->fvv || !P->vol
		|| qaws_volume_mesh_get_cells(cv, ty, vm) != QAWS_STATUS_OK || qaws_volume_mesh_get_cell_adjacency(adj, vm) != QAWS_STATUS_OK
		|| qaws_volume_mesh_get_positions(pos, vm) != QAWS_STATUS_OK)
		goto done;
	for (c = 0, k = 0, i = 0; c < nc; c++)
	{
		cfirst[c] = k;
		ffirst[c] = i;
		k += ty[c] == QAWS_CELL_TETRAHEDRON ? 4 : ty[c] == QAWS_CELL_HEXAHEDRON ? 8 : ty[c] == QAWS_CELL_PRISM ? 6 : 5;
		i += ty[c] == QAWS_CELL_TETRAHEDRON ? 4 : ty[c] == QAWS_CELL_HEXAHEDRON ? 6 : 5;
	}
	ffirst[nc] = i;
	for (i = 0; i < nvv; i++)
		map[i] = CN_NONE;
	for (i = 0; i < count; i++)
	{
		unsigned int f = facets[i], lo = 0, hi = nc, j;
		int const* lf;
		if (f >= nfac || adj[f] != QAWS_MESH_NONE)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		while (hi - lo > 1)
		{
			unsigned int mid = (lo + hi) / 2;
			if (ffirst[mid] <= f)
				lo = mid;
			else
				hi = mid;
		}
		c = lo;
		j = f - ffirst[c];
		lf = ty[c] == QAWS_CELL_TETRAHEDRON ? cn_tet_f[j] : ty[c] == QAWS_CELL_HEXAHEDRON ? cn_hex_f[j] : ty[c] == QAWS_CELL_PRISM ? cn_pri_f[j] : cn_pyr_f[j];
		P->fsz[i] = lf[3] < 0 ? 3 : 4;
		for (k = 0; k < P->fsz[i]; k++)
		{
			unsigned int v = cv[cfirst[c] + (unsigned int)lf[k]];
			if (map[v] == CN_NONE)
			{
				map[v] = P->nv;
				P->vol[P->nv++] = v;
			}
			P->fvv[4 * (size_t)i + k] = map[v];
		}
	}
	P->nf = count;
	P->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, P->nv));
	P->n = (double*)qaws_internal_zalloc(3 * (size_t)P->nv, sizeof(double));
	he = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(2 * sizeof(unsigned long long), 4, count));
	if (!P->x || !P->n || !he)
		goto done;
	for (i = 0; i < P->nv; i++)
		for (k = 0; k < 3; k++)
			P->x[3 * (size_t)i + k] = (double)pos[3 * (size_t)P->vol[i] + k];
	/* outward normals (area-weighted, Newell), the edges: their half-edges by key, a single use is the border */
	for (i = 0; i < count; i++)
	{
		double nn[3] = { 0, 0, 0 };
		for (k = 0; k < P->fsz[i]; k++)
		{
			double const *p = P->x + 3 * (size_t)P->fvv[4 * (size_t)i + k], *q = P->x + 3 * (size_t)P->fvv[4 * (size_t)i + (k + 1) % P->fsz[i]];
			unsigned long long a = P->fvv[4 * (size_t)i + k], b = P->fvv[4 * (size_t)i + (k + 1) % P->fsz[i]];
			nn[0] += (p[1] - q[1]) * (p[2] + q[2]);
			nn[1] += (p[2] - q[2]) * (p[0] + q[0]);
			nn[2] += (p[0] - q[0]) * (p[1] + q[1]);
			he[2 * (size_t)ne] = a < b ? (a << 32 | b) : (b << 32 | a);
			he[2 * (size_t)ne++ + 1] = a << 32 | b;
		}
		for (k = 0; k < P->fsz[i]; k++)
		{
			double* nv = P->n + 3 * (size_t)P->fvv[4 * (size_t)i + k];
			nv[0] += nn[0];
			nv[1] += nn[1];
			nv[2] += nn[2];
		}
	}
	for (i = 0; i < P->nv; i++)
		cn_norm(P->n + 3 * (size_t)i);
	qsort(he, ne, 2 * sizeof(unsigned long long), cn_cmp_u64);
	{
		/* the border's directed edges; a disk: one loop, V - E + F = 1 */
		unsigned int* nxt = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)P->nv + 1));
		unsigned int nb = 0, edges = 0, j, v;
		if (!nxt)
			goto done;
		for (i = 0; i < P->nv; i++)
			nxt[i] = CN_NONE;
		for (i = 0; i < ne; i = j)
		{
			for (j = i + 1; j < ne && he[2 * (size_t)j] == he[2 * (size_t)i]; j++)
				;
			edges++;
			if (j == i + 1)
			{
				unsigned int a = (unsigned int)(he[2 * (size_t)i + 1] >> 32), b = (unsigned int)he[2 * (size_t)i + 1];
				if (nxt[a] != CN_NONE)
				{
					qaws_internal_free(nxt);
					st = QAWS_STATUS_INVALID_ARGUMENT;
					goto done;
				}
				nxt[a] = b;
				nb++;
			}
		}
		if (nb < 3 || (int)P->nv - (int)edges + (int)count != 1)
		{
			qaws_internal_free(nxt);
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		P->loop = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * nb);
		if (!P->loop)
		{
			qaws_internal_free(nxt);
			goto done;
		}
		for (v = 0; v < P->nv && nxt[v] == CN_NONE; v++)
			;
		for (i = 0; i < nb && v != CN_NONE; i++)
		{
			P->loop[i] = v;
			v = nxt[v];
		}
		qaws_internal_free(nxt);
		if (i != nb || v != P->loop[0])
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		P->nb = nb;
		P->h = 0;
		for (i = 0; i < nb; i++)
			P->h += cn_dist(P->x + 3 * (size_t)P->loop[i], P->x + 3 * (size_t)P->loop[(i + 1) % nb]) / nb;
	}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(cv); qaws_internal_free(ty); qaws_internal_free(adj); qaws_internal_free(cfirst); qaws_internal_free(ffirst);
	qaws_internal_free(map); qaws_internal_free(pos); qaws_internal_free(he);
	if (st != QAWS_STATUS_OK)
		cn_patch_free(P);
	return st;
}

/* the patch's Tutte embedding: its border loop (in order[] from its start) on the unit circle by arc length, the inner vertices
   at their neighbours' mean (the graph's uniform Laplacian, solved) */
static qaws_status cn_tutte(cn_patch* P, unsigned int const* order)
{
	unsigned int i, k, n = P->nv, nt = 0;
	unsigned char* move = (unsigned char*)qaws_internal_malloc((size_t)n + 1);
	unsigned long long* key = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), 4, P->nf + 1));
	double *rhs = (double*)qaws_internal_zalloc(2 * (size_t)n + 1, sizeof(double)), len = 0, acc = 0;
	sm_csr K;
	sm_system S;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&K, 0, sizeof K);
	memset(&S, 0, sizeof S);
	P->e = (double*)qaws_internal_zalloc(2 * (size_t)n + 1, sizeof(double));
	if (!move || !key || !rhs || !P->e)
		goto done;
	for (i = 0; i < n; i++)
		move[i] = 1;
	for (i = 0; i < P->nb; i++)
		len += cn_dist(P->x + 3 * (size_t)order[i], P->x + 3 * (size_t)order[(i + 1) % P->nb]);
	for (i = 0; i < P->nb; i++)
	{
		double a = 2 * 3.14159265358979323846 * acc / len;
		P->e[2 * (size_t)order[i]] = cos(a);
		P->e[2 * (size_t)order[i] + 1] = sin(a);
		move[order[i]] = 0;
		acc += cn_dist(P->x + 3 * (size_t)order[i], P->x + 3 * (size_t)order[(i + 1) % P->nb]);
	}
	/* the uniform graph Laplacian over the facets' own edges */
	for (i = 0; i < P->nf; i++)
		for (k = 0; k < P->fsz[i]; k++)
		{
			unsigned long long a = P->fvv[4 * (size_t)i + k], b = P->fvv[4 * (size_t)i + (k + 1) % P->fsz[i]];
			key[nt++] = a < b ? (a << 32 | b) : (b << 32 | a);
		}
	qsort(key, nt, sizeof *key, cn_cmp_u64);
	{
		unsigned int ne = 0, *cnt;
		for (i = 0; i < nt; i++)
			if (!i || key[i] != key[i - 1])
				key[ne++] = key[i];
		K.n = n;
		K.rp = (unsigned int*)qaws_internal_zalloc((size_t)n + 2, sizeof(unsigned int));
		cnt = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
		if (!K.rp || !cnt)
		{
			qaws_internal_free(cnt);
			goto done;
		}
		for (i = 0; i < n; i++)
			K.rp[i + 1] = 1;
		for (i = 0; i < ne; i++)
		{
			K.rp[(unsigned int)(key[i] >> 32) + 1]++;
			K.rp[(unsigned int)key[i] + 1]++;
		}
		for (i = 0; i < n; i++)
			K.rp[i + 1] += K.rp[i];
		K.ci = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)K.rp[n] + 1));
		K.v = (double*)qaws_internal_zalloc((size_t)K.rp[n] + 1, sizeof(double));
		if (!K.ci || !K.v)
		{
			qaws_internal_free(cnt);
			goto done;
		}
		for (i = 0; i < n; i++)
		{
			cnt[i] = K.rp[i] + 1;
			K.ci[K.rp[i]] = i;
		}
		for (i = 0; i < ne; i++)
		{
			unsigned int a = (unsigned int)(key[i] >> 32), b = (unsigned int)key[i];
			K.ci[cnt[a]] = b;
			K.v[cnt[a]++] = -1;
			K.ci[cnt[b]] = a;
			K.v[cnt[b]++] = -1;
			K.v[K.rp[a]] += 1;
			K.v[K.rp[b]] += 1;
		}
		qaws_internal_free(cnt);
		/* (columns sorted within each row) */
		for (i = 0; i < n; i++)
		{
			unsigned int p, q;
			for (p = K.rp[i] + 1; p < K.rp[i + 1]; p++)
			{
				unsigned int c = K.ci[p];
				double w = K.v[p];
				for (q = p; q > K.rp[i] && K.ci[q - 1] > c; q--)
				{
					K.ci[q] = K.ci[q - 1];
					K.v[q] = K.v[q - 1];
				}
				K.ci[q] = c;
				K.v[q] = w;
			}
		}
	}
	st = sm_system_factor(&S, &K, P->x, move);
	if (st == QAWS_STATUS_OK)
		st = sm_system_solve(&S, P->e, rhs, 2);
done:
	sm_system_free(&S);
	sm_csr_free(&K);
	qaws_internal_free(move);
	qaws_internal_free(key);
	qaws_internal_free(rhs);
	return st;
}

qaws_status qaws_volume_mesh_connect(qaws_volume_mesh** out, qaws_connect_report* report, qaws_volume_mesh const* a, unsigned int const* af,
	qaws_volume_mesh const* b, unsigned int const* bf, unsigned int nafc, unsigned int nbfc, qaws_connect_desc const* desc)
{
	qaws_connect_report local;
	cn_patch PA, PB;
	unsigned int *orda = NULL, *ordb = NULL, *match = NULL, *used = NULL, *cells = NULL, i, k, r, R, nva, nvb, nca, ncb,
		ka, kb, nbc, nbv;
	qaws_cell_type* types = NULL;
	qaws_scalar* pos = NULL;
	int same = !b || b == a;
	double gap = 0;
	qaws_status st;
	memset(&PA, 0, sizeof PA);
	memset(&PB, 0, sizeof PB);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!a)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (same)
		b = a;
	st = cn_patch_build(&PA, a, af, nafc);
	if (st == QAWS_STATUS_OK)
		st = cn_patch_build(&PB, b, bf, nbfc);
	if (st != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_INVALID_ARGUMENT;
	if (PA.nv != PB.nv || PA.nf != PB.nf || PA.nb != PB.nb)
		goto done;   /* (not the same triangulation) */
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	orda = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * PA.nb);
	ordb = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * PB.nb);
	match = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * PA.nv);
	used = (unsigned int*)qaws_internal_zalloc(PB.nv, sizeof(unsigned int));
	if (!orda || !ordb || !match || !used)
		goto done;
	{
		/* b's border turned the other way (its facets face a), started where the vertices meet a's nearest */
		unsigned int best = 0, j;
		double bc = HUGE_VAL;
		for (j = 0; j < PB.nb; j++)
		{
			double c = 0;
			for (i = 0; i < PA.nb; i++)
				c += cn_dist(PA.x + 3 * (size_t)PA.loop[i], PB.x + 3 * (size_t)PB.loop[(j + PB.nb - i) % PB.nb]);
			if (c < bc)
			{
				bc = c;
				best = j;
			}
		}
		for (i = 0; i < PA.nb; i++)
		{
			orda[i] = PA.loop[i];
			ordb[i] = PB.loop[(best + PB.nb - i) % PB.nb];
		}
	}
	if (cn_tutte(&PA, orda) != QAWS_STATUS_OK || cn_tutte(&PB, ordb) != QAWS_STATUS_OK)
		goto done;
	/* vertex for vertex by the embeddings (the same where the triangulations are), then facet for facet */
	st = QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < PA.nv; i++)
	{
		unsigned int j, bj = CN_NONE;
		double bd = HUGE_VAL;
		for (j = 0; j < PB.nv; j++)
		{
			double du = PA.e[2 * (size_t)i] - PB.e[2 * (size_t)j], dv = PA.e[2 * (size_t)i + 1] - PB.e[2 * (size_t)j + 1], d = du * du + dv * dv;
			if (d < bd)
			{
				bd = d;
				bj = j;
			}
		}
		if (bj == CN_NONE || bd > 1e-12 || used[bj])
			goto done;
		used[bj] = 1;
		match[i] = bj;
		gap += cn_dist(PA.x + 3 * (size_t)i, PB.x + 3 * (size_t)bj) / PA.nv;
	}
	{
		/* each of a's facets, mapped, among b's (vertex sets compared sorted) */
		unsigned long long* keys = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), 2, PB.nf + 1));
		if (!keys)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (i = 0; i < PB.nf; i++)
		{
			unsigned int s[4] = { CN_NONE, CN_NONE, CN_NONE, CN_NONE }, p, q;
			for (p = 0; p < PB.fsz[i]; p++)
				s[p] = PB.fvv[4 * (size_t)i + p];
			for (p = 0; p < 4; p++)
				for (q = p + 1; q < 4; q++)
					if (s[q] < s[p])
					{
						unsigned int tmp = s[p];
						s[p] = s[q];
						s[q] = tmp;
					}
			keys[2 * (size_t)i] = (unsigned long long)s[0] << 32 | s[1];
			keys[2 * (size_t)i + 1] = (unsigned long long)s[2] << 32 | s[3];
		}
		for (i = 0; i < PA.nf; i++)
		{
			unsigned int s[4] = { CN_NONE, CN_NONE, CN_NONE, CN_NONE }, p, q, found = 0;
			unsigned long long k0, k1;
			for (p = 0; p < PA.fsz[i]; p++)
				s[p] = match[PA.fvv[4 * (size_t)i + p]];
			for (p = 0; p < 4; p++)
				for (q = p + 1; q < 4; q++)
					if (s[q] < s[p])
					{
						unsigned int tmp = s[p];
						s[p] = s[q];
						s[q] = tmp;
					}
			k0 = (unsigned long long)s[0] << 32 | s[1];
			k1 = (unsigned long long)s[2] << 32 | s[3];
			for (p = 0; p < PB.nf && !found; p++)
				found = keys[2 * (size_t)p] == k0 && keys[2 * (size_t)p + 1] == k1;
			if (!found)
			{
				qaws_internal_free(keys);
				goto done;
			}
		}
		qaws_internal_free(keys);
	}
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	R = desc && desc->rows ? desc->rows : (unsigned int)fmax(1, floor(gap / (0.5 * (PA.h + PB.h)) + 0.5));
	report->rows = R;
	report->loop_vertices[0] = PA.nv;
	report->loop_vertices[1] = PB.nv;
	report->gap = (qaws_scalar)gap;
	/* the result: a's vertices and cells, b's, the layers' vertices (1 to R - 1, the patch's in turn), the bridge's cells */
	nva = qaws_volume_mesh_get_vertex_count(a);
	nvb = same ? 0 : qaws_volume_mesh_get_vertex_count(b);
	nca = qaws_volume_mesh_get_cell_count(a);
	ncb = same ? 0 : qaws_volume_mesh_get_cell_count(b);
	ka = qaws_volume_mesh_get_corner_count(a);
	kb = same ? 0 : qaws_volume_mesh_get_corner_count(b);
	nbv = (R - 1) * PA.nv;
	nbc = R * PA.nf;
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nva + nvb + nbv + 1));
	cells = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)ka + kb + 8 * (size_t)nbc + 1));
	types = (qaws_cell_type*)qaws_internal_malloc(sizeof(qaws_cell_type) * ((size_t)nca + ncb + nbc + 1));
	if (!pos || !cells || !types || qaws_volume_mesh_get_positions(pos, a) != QAWS_STATUS_OK || qaws_volume_mesh_get_cells(cells, types, a) != QAWS_STATUS_OK)
		goto done;
	if (!same)
	{
		if (qaws_volume_mesh_get_positions(pos + 3 * (size_t)nva, b) != QAWS_STATUS_OK || qaws_volume_mesh_get_cells(cells + ka, types + nca, b) != QAWS_STATUS_OK)
			goto done;
		for (i = 0; i < kb; i++)
			cells[ka + i] += nva;
	}
	for (r = 1; r < R; r++)
	{
		double tt = (double)r / R, h[4];
		cn_hermite(tt, h);
		for (i = 0; i < PA.nv; i++)
		{
			double const *xa = PA.x + 3 * (size_t)i, *na = PA.n + 3 * (size_t)i, *xb = PB.x + 3 * (size_t)match[i], *nb = PB.n + 3 * (size_t)match[i];
			double L = cn_dist(xa, xb);
			for (k = 0; k < 3; k++)
				pos[3 * ((size_t)nva + nvb + (size_t)(r - 1) * PA.nv + i) + k] = (qaws_scalar)(h[0] * xa[k] + h[1] * L * na[k] + h[2] * xb[k] - h[3] * L * nb[k]);
		}
	}
	{
		/* each of a's facets extruded through the layers: its own order is the cells' base (its normal toward the bridge) */
		unsigned int nc = ka + kb, ncell = nca + ncb, f;
		for (r = 0; r < R; r++)
			for (f = 0; f < PA.nf; f++)
			{
				unsigned int s = PA.fsz[f], p;
				for (k = 0; k < 2; k++)
					for (p = 0; p < s; p++)
					{
						unsigned int lv = PA.fvv[4 * (size_t)f + p], layer = r + k, v;
						v = layer == 0 ? PA.vol[lv] : layer == R ? (same ? PB.vol[match[lv]] : nva + PB.vol[match[lv]]) : nva + nvb + (layer - 1) * PA.nv + lv;
						cells[nc++] = v;
					}
				types[ncell++] = s == 3 ? QAWS_CELL_PRISM : QAWS_CELL_HEXAHEDRON;
			}
		{
			qaws_volume_mesh_desc d;
			memset(&d, 0, sizeof d);
			d.positions = pos;
			d.vertex_count = nva + nvb + nbv;
			d.cell_vertices = cells;
			d.cell_types = types;
			d.cell_count = ncell;
			st = qaws_volume_mesh_create(out, &d);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
		report->bridge_faces = nbc;
		{
			/* the bridge's cells measured: dihedral angles, scaled Jacobians */
			qaws_scalar* v = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * ((size_t)ncell + 1));
			if (v)
			{
				unsigned int c;
				double lo = HUGE_VAL, hi = -HUGE_VAL, sj = HUGE_VAL;
				if (qaws_volume_mesh_measure_elements(v, *out, QAWS_MEASURE_MIN_ANGLE) == QAWS_STATUS_OK)
					for (c = nca + ncb; c < ncell; c++)
						lo = v[c] == v[c] && v[c] < lo ? v[c] : lo;
				if (qaws_volume_mesh_measure_elements(v, *out, QAWS_MEASURE_MAX_ANGLE) == QAWS_STATUS_OK)
					for (c = nca + ncb; c < ncell; c++)
						hi = v[c] == v[c] && v[c] > hi ? v[c] : hi;
				if (qaws_volume_mesh_measure_elements(v, *out, QAWS_MEASURE_SCALED_JACOBIAN) == QAWS_STATUS_OK)
					for (c = nca + ncb; c < ncell; c++)
						sj = v[c] == v[c] && v[c] < sj ? v[c] : sj;
				report->min_angle = (qaws_scalar)(lo * 3.14159265358979323846 / 180);
				report->max_angle = (qaws_scalar)(hi * 3.14159265358979323846 / 180);
				report->min_scaled_jacobian = (qaws_scalar)sj;
				qaws_internal_free(v);
			}
		}
	}
	st = QAWS_STATUS_OK;
done:
	cn_patch_free(&PA);
	cn_patch_free(&PB);
	qaws_internal_free(orda); qaws_internal_free(ordb); qaws_internal_free(match); qaws_internal_free(used);
	qaws_internal_free(cells); qaws_internal_free(types);
	qaws_internal_free(pos);
	return st;
}
