/*
 * Implicit sweeps (qaws_mesh_sweep.h): the swept solid as the level set of
 * the least over time of the moving shape's signed distance, sampled in time
 * (no vertex moving more than a cell between samples), each point's samples
 * pruned by their bounding spheres, a rigid motion's best time refined by
 * golden sections; contoured by the offsets' dual contouring
 * (wt_contour_field, blocks pruned: the field is 1-Lipschitz). Attributes
 * from each result corner's source: the shape's closest point at its best
 * time (swp_carry).
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_sweep.h"
#include "../qaws_mesh_bvh.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SWI_CHUNK 4096u   /* points per query batch */
#define SWI_GOLDEN 10u    /* golden section rounds (the bracket shrunk to 0.618^10, about 1 / 120) */

typedef struct swi_ctx
{
	qaws_motion const* m;
	int pv;                         /* a vertex motion */
	unsigned int nv, nt;
	unsigned int* t;
	qaws_scalar* x;                 /* the positions as given (rigid: the body frame) */
	qaws_scalar* buf;               /* nv points */
	qaws_mesh_bvh* b;               /* rigid: over x; vertex: refit to sample cur */
	unsigned int cur;
	unsigned int ns;
	double* ts;                     /* the samples' times */
	double* pose;                   /* rigid: per sample R (by rows) and T, 12 each */
	double* cen;                    /* per sample the shape's bounding sphere: centre (3) and radius */
	double band;                    /* points within it of the level have their time refined (rigid) */
	double cell;
	double tie;                     /* values within it equal: the earlier time kept (a point a face touches over a span of time:
	                                   the first touch stands) */
	qaws_batch_executor const* ex;
	/* per chunk */
	qaws_scalar *q, *d;
	qaws_mesh_hit* h;
	double *best, *bt, *hp, *hw;    /* the least value, its time, the closest point (world), its weights on its triangle */
	unsigned int *bj, *idx, *htri;
} swi_ctx;

static qaws_status swi_pose(qaws_motion const* m, double t, double* P)
{
	qaws_motion_state s;
	unsigned int k;
	qaws_status st = qaws_motion_evaluate(&s, m, (qaws_scalar)t);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < 9; k++)
		P[k] = (double)s.rotation[k];
	P[9] = s.position.x; P[10] = s.position.y; P[11] = s.position.z;
	return QAWS_STATUS_OK;
}

/* the vertex motion's positions at sample j into the BVH */
static qaws_status swi_refit(swi_ctx* c, unsigned int j)
{
	qaws_status st;
	if (c->cur == j)
		return QAWS_STATUS_OK;
	if ((st = qaws_motion_transform_points(c->buf, NULL, c->m, (qaws_scalar)c->ts[j], NULL, c->nv)) != QAWS_STATUS_OK)
		return st;
	c->cur = j;
	return qaws_mesh_bvh_refit(c->b, c->buf);
}

/* the signed distance at the k points idx[] of the chunk (positions X), each at its own pose (rigid: P, 12 per point) or all
   at the BVH's current sample (vertex: P NULL), kept where less than the best (times tt, samples jj) */
static qaws_status swi_query(swi_ctx* c, double const* X, unsigned int k, double const* P, double const* tt, unsigned int const* jj)
{
	unsigned int i, a;
	qaws_status st;
	for (i = 0; i < k; i++)
	{
		double const* x = X + 3 * (size_t)c->idx[i];
		if (P)
		{
			double const* R = P + 12 * (size_t)i;
			double y[3] = { x[0] - R[9], x[1] - R[10], x[2] - R[11] };
			for (a = 0; a < 3; a++)
				c->q[3 * (size_t)i + a] = (qaws_scalar)(R[a] * y[0] + R[3 + a] * y[1] + R[6 + a] * y[2]);
		}
		else
			for (a = 0; a < 3; a++)
				c->q[3 * (size_t)i + a] = (qaws_scalar)x[a];
	}
	if ((st = qaws_mesh_bvh_signed_distance(c->d, c->h, c->b, c->q, k, c->ex)) != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < k; i++)
	{
		unsigned int p = c->idx[i];
		if (!((double)c->d[i] < c->best[p] - c->tie || ((double)c->d[i] <= c->best[p] + c->tie && tt[i] < c->bt[p])))
			continue;
		c->best[p] = (double)c->d[i];
		c->bt[p] = tt[i];
		c->bj[p] = jj[i];
		c->htri[p] = c->h[i].element;
		for (a = 0; a < 3; a++)
		{
			double hq = (double)c->h[i].point[a];
			if (P)
			{
				double const* R = P + 12 * (size_t)i;
				hq = R[3 * a] * (double)c->h[i].point[0] + R[3 * a + 1] * (double)c->h[i].point[1] + R[3 * a + 2] * (double)c->h[i].point[2] + R[9 + a];
			}
			c->hp[3 * (size_t)p + a] = hq;
			c->hw[3 * (size_t)p + a] = (double)c->h[i].weights[a];
		}
	}
	return QAWS_STATUS_OK;
}

/* the samples' queries for the chunk's m points: those of each sample j whose sphere may hold a point nearer than the best
   (first the sample nearest each point, then the others), one sample at a time */
static qaws_status swi_samples(swi_ctx* c, double const* X, unsigned int m, double* P, double* tt, unsigned int* jj)
{
	unsigned int pass, j, i, k;
	qaws_status st = QAWS_STATUS_OK;
	for (i = 0; i < m; i++)
	{
		double bd = HUGE_VAL;
		c->best[i] = HUGE_VAL;
		c->bt[i] = HUGE_VAL;
		c->bj[i] = 0;
		for (j = 0; j < c->ns; j++)
		{
			double const* s = c->cen + 4 * (size_t)j;
			double e = sqrt((X[3 * i] - s[0]) * (X[3 * i] - s[0]) + (X[3 * i + 1] - s[1]) * (X[3 * i + 1] - s[1]) + (X[3 * i + 2] - s[2]) * (X[3 * i + 2] - s[2])) - s[3];
			if (e < bd)
			{
				bd = e;
				c->bj[i] = j;
			}
		}
	}
	for (pass = 0; pass < 2 && st == QAWS_STATUS_OK; pass++)
		for (j = 0; j < c->ns && st == QAWS_STATUS_OK; j++)
		{
			double const* s = c->cen + 4 * (size_t)j;
			for (i = k = 0; i < m; i++)
			{
				double e;
				if (pass == 0 ? c->bj[i] != j || c->best[i] < HUGE_VAL : c->bj[i] == j)
					continue;
				if (pass == 1)
				{
					e = sqrt((X[3 * i] - s[0]) * (X[3 * i] - s[0]) + (X[3 * i + 1] - s[1]) * (X[3 * i + 1] - s[1]) + (X[3 * i + 2] - s[2]) * (X[3 * i + 2] - s[2])) - s[3];
					if (e >= c->best[i])
						continue;
				}
				c->idx[k] = i;
				tt[k] = c->ts[j];
				jj[k] = j;
				if (!c->pv)
					memcpy(P + 12 * (size_t)k, c->pose + 12 * (size_t)j, sizeof(double) * 12);
				k++;
			}
			if (!k)
				continue;
			if (c->pv && (st = swi_refit(c, j)) != QAWS_STATUS_OK)
				break;
			st = swi_query(c, X, k, c->pv ? NULL : P, tt, jj);
		}
	return st;
}

/* a rigid motion's best times refined: golden sections over the samples on either side of each point's best one, for the
   points within the band of the level (G: 8 m scratch) */
static qaws_status swi_refine(swi_ctx* c, double const* X, unsigned int m, double* P, double* tt, unsigned int* jj, double* G)
{
	static double const gr = 0.61803398874989484820;
	unsigned int i, k, r, ns = 0;
	qaws_status st = QAWS_STATUS_OK;
	double *A = G, *B = G + m, *t1 = G + 2 * (size_t)m, *t2 = G + 3 * (size_t)m, *f1 = G + 4 * (size_t)m, *f2 = G + 5 * (size_t)m;
	unsigned int *sel = (unsigned int*)(G + 6 * (size_t)m), *which = (unsigned int*)(G + 7 * (size_t)m);
	for (i = 0; i < m; i++)
		if (fabs(c->best[i]) < c->band)
		{
			unsigned int j = c->bj[i];
			A[ns] = c->ts[j > 0 ? j - 1 : 0];
			B[ns] = c->ts[j + 1 < c->ns ? j + 1 : j];
			t1[ns] = B[ns] - gr * (B[ns] - A[ns]);
			t2[ns] = A[ns] + gr * (B[ns] - A[ns]);
			sel[ns++] = i;
		}
	/* round 0 evaluates t1, round 1 t2, then each round the bracket shrunk toward the lesser and its new point evaluated */
	for (r = 0; ns && r < SWI_GOLDEN + 2 && st == QAWS_STATUS_OK; r++)
	{
		for (k = 0; k < ns && st == QAWS_STATUS_OK; k++)
		{
			if (r < 2)
				which[k] = r + 1;
			else if (f1[k] < f2[k])
			{
				B[k] = t2[k]; t2[k] = t1[k]; f2[k] = f1[k];
				t1[k] = B[k] - gr * (B[k] - A[k]);
				which[k] = 1;
			}
			else
			{
				A[k] = t1[k]; t1[k] = t2[k]; f1[k] = f2[k];
				t2[k] = A[k] + gr * (B[k] - A[k]);
				which[k] = 2;
			}
			c->idx[k] = sel[k];
			tt[k] = which[k] == 1 ? t1[k] : t2[k];
			jj[k] = c->bj[sel[k]];
			st = swi_pose(c->m, tt[k], P + 12 * (size_t)k);
		}
		if (st == QAWS_STATUS_OK)
			st = swi_query(c, X, ns, P, tt, jj);
		for (k = 0; k < ns; k++)
			*(which[k] == 1 ? f1 + k : f2 + k) = (double)c->d[k];
	}
	return st;
}

static qaws_status swi_eval(void* vctx, double const* x, size_t n, double* phi, double* grad)
{
	swi_ctx* c = (swi_ctx*)vctx;
	size_t c0, nfd = 0, q;
	unsigned int i, a;
	qaws_status st = QAWS_STATUS_OK;
	double* P = (double*)qaws_internal_scratch_malloc(sizeof(double) * 12 * SWI_CHUNK);
	double* tt = (double*)qaws_internal_scratch_malloc(sizeof(double) * SWI_CHUNK);
	double* G = (double*)qaws_internal_scratch_malloc(sizeof(double) * 8 * SWI_CHUNK);
	unsigned int* jj = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * SWI_CHUNK);
	size_t* fd = grad ? (size_t*)qaws_internal_scratch_malloc(sizeof(size_t) * (n + 1)) : NULL;
	if (!P || !tt || !G || !jj || (grad && !fd))
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (c0 = 0; c0 < n && st == QAWS_STATUS_OK; c0 += SWI_CHUNK)
	{
		unsigned int m = (unsigned int)(n - c0 < SWI_CHUNK ? n - c0 : SWI_CHUNK);
		double const* X = x + 3 * c0;
		if ((st = swi_samples(c, X, m, P, tt, jj)) != QAWS_STATUS_OK)
			break;
		if (!c->pv && (st = swi_refine(c, X, m, P, tt, jj, G)) != QAWS_STATUS_OK)
			break;
		for (i = 0; i < m; i++)
		{
			double g[3], l = 0;
			phi[c0 + i] = c->best[i];
			if (!grad)
				continue;
			for (a = 0; a < 3; a++)
			{
				g[a] = X[3 * i + a] - c->hp[3 * (size_t)i + a];
				l += g[a] * g[a];
			}
			l = sqrt(l);
			for (a = 0; a < 3; a++)
				grad[3 * (c0 + i) + a] = l > 0 ? g[a] / (c->best[i] < 0 ? -l : l) : 0;
			if (l < 0.05 * c->cell)
				fd[nfd++] = c0 + i;
		}
	}
	if (st == QAWS_STATUS_OK && nfd)
	{
		/* on the surface the closest point may sit on an edge or corner of the shape, and the least time a hair off moves it
		   onto a face (an axis for a diagonal wall; which hit is right cannot be told there): the gradient from differences of
		   phi a fiftieth of a cell along each axis */
		double e = 0.02 * c->cell;
		double* Y = (double*)qaws_internal_scratch_malloc(sizeof(double) * 9 * nfd);
		double* v = (double*)qaws_internal_scratch_malloc(sizeof(double) * 3 * nfd);
		if (!Y || !v)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (q = 0; st == QAWS_STATUS_OK && q < nfd; q++)
			for (a = 0; a < 3; a++)
			{
				unsigned int b;
				for (b = 0; b < 3; b++)
					Y[9 * q + 3 * a + b] = x[3 * fd[q] + b] + (b == a ? e : 0);
			}
		if (st == QAWS_STATUS_OK)
			st = swi_eval(c, Y, 3 * nfd, v, NULL);
		for (q = 0; st == QAWS_STATUS_OK && q < nfd; q++)
		{
			double f0 = phi[fd[q]], g[3] = { v[3 * q] - f0, v[3 * q + 1] - f0, v[3 * q + 2] - f0 }, l = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
			if (l > 0)
				for (a = 0; a < 3; a++)
					grad[3 * fd[q] + a] = g[a] / l;
		}
		qaws_internal_free(Y); qaws_internal_free(v);
	}
	qaws_internal_free(P); qaws_internal_free(tt); qaws_internal_free(G); qaws_internal_free(jj); qaws_internal_free(fd);
	return st;
}

/* the weights on triangle f (rest positions) of the point there nearest q in its plane, clamped into it */
static void swi_weights(swi_ctx const* c, unsigned int f, double const* q, double* w)
{
	double V0[3], V1[3], V2[3], e0[3], e1[3], d[3], a00 = 0, a01 = 0, a11 = 0, b0 = 0, b1 = 0, den, s;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		V0[k] = (double)c->x[3 * (size_t)c->t[3 * (size_t)f] + k];
		V1[k] = (double)c->x[3 * (size_t)c->t[3 * (size_t)f + 1] + k];
		V2[k] = (double)c->x[3 * (size_t)c->t[3 * (size_t)f + 2] + k];
		e0[k] = V1[k] - V0[k];
		e1[k] = V2[k] - V0[k];
		d[k] = q[k] - V0[k];
		a00 += e0[k] * e0[k]; a01 += e0[k] * e1[k]; a11 += e1[k] * e1[k]; b0 += e0[k] * d[k]; b1 += e1[k] * d[k];
	}
	den = a00 * a11 - a01 * a01;
	w[1] = den > 0 ? (a11 * b0 - a01 * b1) / den : 0;
	w[2] = den > 0 ? (a00 * b1 - a01 * b0) / den : 0;
	w[0] = 1 - w[1] - w[2];
	for (k = 0, s = 0; k < 3; k++)
		s += (w[k] = fmax(0, w[k]));
	for (k = 0; k < 3; k++)
		w[k] = s > 0 ? w[k] / s : 1.0 / 3;
}

/* a rigid motion's ties pushed back to their first time: for the points within the band whose best is tied with an earlier
   time than the sample before it (a face touching them over a span), bisected between that sample (not tied) and the best
   time; the first touch stands */
static qaws_status swi_first(swi_ctx* c, double const* X, unsigned int m, double* P, double* tt, unsigned int* jj, double* G)
{
	double *lo = G, *hi = G + m, *ref = G + 2 * (size_t)m;
	unsigned int* sel = (unsigned int*)(G + 3 * (size_t)m);
	unsigned int i, k, r, ns = 0;
	qaws_status st = QAWS_STATUS_OK;
	for (i = 0; i < m; i++)
		if (fabs(c->best[i]) < c->band && c->bj[i] > 0 && c->bt[i] > c->ts[c->bj[i] - 1])
		{
			lo[ns] = c->ts[c->bj[i] - 1];
			hi[ns] = c->bt[i];
			ref[ns] = c->best[i];
			sel[ns++] = i;
		}
	for (r = 0; ns && r < SWI_GOLDEN && st == QAWS_STATUS_OK; r++)
	{
		for (k = 0; k < ns && st == QAWS_STATUS_OK; k++)
		{
			c->idx[k] = sel[k];
			tt[k] = 0.5 * (lo[k] + hi[k]);
			jj[k] = c->bj[sel[k]];
			st = swi_pose(c->m, tt[k], P + 12 * (size_t)k);
		}
		if (st == QAWS_STATUS_OK)
			st = swi_query(c, X, ns, P, tt, jj);
		for (k = 0; k < ns; k++)
		{
			if ((double)c->d[k] <= ref[k] + c->tie)
				hi[k] = tt[k];
			else
				lo[k] = tt[k];
		}
	}
	return st;
}

/* the sources of n points: per point its triangle, weights and best time (the field's evaluation, kept; ties at their
   first time) */
static qaws_status swi_sources(swi_ctx* c, double const* x, size_t n, unsigned int* tri, double* w, double* tb)
{
	size_t c0;
	unsigned int i, k;
	qaws_status st = QAWS_STATUS_OK;
	double* P = (double*)qaws_internal_scratch_malloc(sizeof(double) * 12 * SWI_CHUNK);
	double* tt = (double*)qaws_internal_scratch_malloc(sizeof(double) * SWI_CHUNK);
	double* G = (double*)qaws_internal_scratch_malloc(sizeof(double) * 8 * SWI_CHUNK);
	unsigned int* jj = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * SWI_CHUNK);
	if (!P || !tt || !G || !jj)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (c0 = 0; c0 < n && st == QAWS_STATUS_OK; c0 += SWI_CHUNK)
	{
		unsigned int m = (unsigned int)(n - c0 < SWI_CHUNK ? n - c0 : SWI_CHUNK);
		if ((st = swi_samples(c, x + 3 * c0, m, P, tt, jj)) != QAWS_STATUS_OK)
			break;
		if (!c->pv && ((st = swi_refine(c, x + 3 * c0, m, P, tt, jj, G)) != QAWS_STATUS_OK || (st = swi_first(c, x + 3 * c0, m, P, tt, jj, G)) != QAWS_STATUS_OK))
			break;
		for (i = 0; i < m; i++)
		{
			tri[c0 + i] = c->htri[i];
			tb[c0 + i] = c->bt[i];
			for (k = 0; k < 3; k++)
				w[3 * (c0 + i) + k] = c->hw[3 * (size_t)i + k];
		}
	}
	qaws_internal_free(P); qaws_internal_free(tt); qaws_internal_free(G); qaws_internal_free(jj);
	return st;
}

/* a sample at which the shape comes within reach of an output vertex */
typedef struct swi_near
{
	unsigned int v, j;
	double d;
} swi_near;

static int swi_cmp_near(void const* a, void const* b)
{
	swi_near const *x = (swi_near const*)a, *y = (swi_near const*)b;
	return x->v != y->v ? (x->v < y->v ? -1 : 1) : x->j < y->j ? -1 : x->j > y->j;
}

/* the shape's distance to x at time t (rigid) */
static qaws_status swi_point_distance(swi_ctx* c, double const* x, double t, double* d)
{
	double R[12], y[3];
	qaws_scalar q[3], s;
	unsigned int a;
	qaws_status st;
	if ((st = swi_pose(c->m, t, R)) != QAWS_STATUS_OK)
		return st;
	for (a = 0; a < 3; a++)
		y[a] = x[a] - R[9 + a];
	for (a = 0; a < 3; a++)
		q[a] = (qaws_scalar)(R[a] * y[0] + R[3 + a] * y[1] + R[6 + a] * y[2]);
	if ((st = qaws_mesh_bvh_signed_distance(&s, NULL, c->b, q, 1, c->ex)) != QAWS_STATUS_OK)
		return st;
	*d = fabs((double)s);
	return QAWS_STATUS_OK;
}

/* the shape's closest point to x at time t (rigid) or at sample j (a vertex motion): its triangle and weights */
static qaws_status swi_point_source(swi_ctx* c, double const* x, double t, unsigned int j, unsigned int* tri, double* w)
{
	qaws_mesh_hit h;
	qaws_scalar q[3];
	qaws_status st;
	unsigned int a;
	if (c->pv)
	{
		if ((st = swi_refit(c, j)) != QAWS_STATUS_OK)
			return st;
		for (a = 0; a < 3; a++)
			q[a] = (qaws_scalar)x[a];
	}
	else
	{
		double R[12], y[3];
		if ((st = swi_pose(c->m, t, R)) != QAWS_STATUS_OK)
			return st;
		for (a = 0; a < 3; a++)
			y[a] = x[a] - R[9 + a];
		for (a = 0; a < 3; a++)
			q[a] = (qaws_scalar)(R[a] * y[0] + R[3 + a] * y[1] + R[6 + a] * y[2]);
	}
	if ((st = qaws_mesh_bvh_closest(&h, c->b, q, 1, 0, c->ex)) != QAWS_STATUS_OK)
		return st;
	*tri = h.element;
	for (a = 0; a < 3; a++)
		w[a] = (double)h.weights[a];
	return QAWS_STATUS_OK;
}

/* each output vertex's passes (vpo: nvo + 1 offsets into the pass arrays): the runs of time samples at which the shape comes
   within 1.5 cells of it, each at its nearest sample (a parabola through it and its neighbours in the run), its dwell the
   run's length in time, its source the shape's closest point then; the run holding the first touch keeps the first touch's
   time and source (tb, src, w); a vertex without a run that alone */
static qaws_status swi_passes(swi_ctx* c, double const* X, unsigned int nvo, double const* tb, unsigned int const* src, double const* w,
	unsigned int** out_po, unsigned int** out_src, double** out_w, double** out_t, double** out_dw)
{
	double eps = 1.5 * c->cell, dt = c->ns > 1 ? 1.0 / (c->ns - 1) : 1;
	swi_near* nr = NULL;
	size_t ne = 0, cap = 0, i, r0;
	unsigned int j, k, np = 0, *po = NULL, *ps = NULL, *idx = NULL;
	double *pw = NULL, *pt = NULL, *pd = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	idx = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nvo + 1));
	if (!idx)
		goto done;
	/* the samples one by one: the vertices their sphere may reach, their distances, those within reach kept */
	for (j = 0; j < c->ns; j++)
	{
		double const* s = c->cen + 4 * (size_t)j;
		unsigned int m = 0, b0;
		for (k = 0; k < nvo; k++)
		{
			double const* x = X + 3 * (size_t)k;
			if (sqrt((x[0] - s[0]) * (x[0] - s[0]) + (x[1] - s[1]) * (x[1] - s[1]) + (x[2] - s[2]) * (x[2] - s[2])) - s[3] <= eps)
				idx[m++] = k;
		}
		if (!m)
			continue;
		if (c->pv && (st = swi_refit(c, j)) != QAWS_STATUS_OK)
			goto done;
		for (b0 = 0; b0 < m; b0 += SWI_CHUNK)
		{
			unsigned int bn = m - b0 < SWI_CHUNK ? m - b0 : SWI_CHUNK, q, a;
			double const* R = c->pv ? NULL : c->pose + 12 * (size_t)j;
			for (q = 0; q < bn; q++)
			{
				double const* x = X + 3 * (size_t)idx[b0 + q];
				if (R)
				{
					double y[3] = { x[0] - R[9], x[1] - R[10], x[2] - R[11] };
					for (a = 0; a < 3; a++)
						c->q[3 * (size_t)q + a] = (qaws_scalar)(R[a] * y[0] + R[3 + a] * y[1] + R[6 + a] * y[2]);
				}
				else
					for (a = 0; a < 3; a++)
						c->q[3 * (size_t)q + a] = (qaws_scalar)x[a];
			}
			if ((st = qaws_mesh_bvh_signed_distance(c->d, NULL, c->b, c->q, bn, c->ex)) != QAWS_STATUS_OK)
				goto done;
			for (q = 0; q < bn; q++)
				if (fabs((double)c->d[q]) <= eps)
				{
					if (ne + 1 > cap)
					{
						size_t c2 = cap ? 2 * cap : 4096;
						swi_near* g = (swi_near*)(nr ? qaws_internal_realloc(nr, qaws_internal_mul2(sizeof(swi_near), c2))
						                            : qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(swi_near), c2)));
						if (!g)
						{
							st = QAWS_STATUS_ALLOCATION_FAILURE;
							goto done;
						}
						nr = g;
						cap = c2;
					}
					nr[ne].v = idx[b0 + q];
					nr[ne].j = j;
					nr[ne++].d = fabs((double)c->d[q]);
				}
		}
	}
	if (nr)
		qsort(nr, ne, sizeof(swi_near), swi_cmp_near);
	/* the passes: per vertex its runs (at least one); first counted, then filled */
	po = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nvo + 1, sizeof(unsigned int));
	if (!po)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < ne; i = r0)
	{
		unsigned int runs = 0, v = nr[i].v;
		for (r0 = i; r0 < ne && nr[r0].v == v; r0++)
			runs += r0 == i || nr[r0].j != nr[r0 - 1].j + 1;
		po[v + 1] = runs;
	}
	for (k = 0; k < nvo; k++)
		po[k + 1] = (po[k + 1] ? po[k + 1] : 1) + po[k];
	np = po[nvo];
	ps = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)np + 1));
	pw = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)np + 1));
	pt = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)np + 1));
	pd = (double*)qaws_internal_scratch_zalloc((size_t)np + 1, sizeof(double));
	if (!ps || !pw || !pt || !pd)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (k = 0; k < nvo; k++)
	{
		/* (a vertex without a run: its first touch alone) */
		ps[po[k]] = src[k];
		memcpy(pw + 3 * (size_t)po[k], w + 3 * (size_t)k, sizeof(double) * 3);
		pt[po[k]] = tb[k];
	}
	for (i = 0; i < ne; i = r0)
	{
		unsigned int v = nr[i].v, p = po[v], first = 1;
		size_t a = i;
		for (r0 = i; r0 < ne && nr[r0].v == v; r0++)
			if (r0 + 1 == ne || nr[r0 + 1].v != v || nr[r0 + 1].j != nr[r0].j + 1)
			{
				/* the run [a, r0]: its nearest sample, refined by a parabola */
				size_t b = a, q;
				double t0, t1, tt;
				for (q = a; q <= r0; q++)
					if (nr[q].d < nr[b].d)
						b = q;
				tt = c->ts[nr[b].j];
				if (!c->pv)
				{
					/* (a rigid motion: golden sections over the samples beside it, as the first touch) */
					double lo = fmax(0, tt - dt), hi = fmin(1, tt + dt), g = 0.61803398874989484820, x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo), f1, f2;
					unsigned int it;
					if ((st = swi_point_distance(c, X + 3 * (size_t)v, x1, &f1)) != QAWS_STATUS_OK || (st = swi_point_distance(c, X + 3 * (size_t)v, x2, &f2)) != QAWS_STATUS_OK)
						goto done;
					for (it = 0; it < 20; it++)
					{
						if (f1 <= f2)
						{
							hi = x2; x2 = x1; f2 = f1;
							x1 = hi - g * (hi - lo);
							if ((st = swi_point_distance(c, X + 3 * (size_t)v, x1, &f1)) != QAWS_STATUS_OK)
								goto done;
						}
						else
						{
							lo = x1; x1 = x2; f1 = f2;
							x2 = lo + g * (hi - lo);
							if ((st = swi_point_distance(c, X + 3 * (size_t)v, x2, &f2)) != QAWS_STATUS_OK)
								goto done;
						}
					}
					tt = 0.5 * (lo + hi);
				}
				t0 = c->ts[nr[a].j] - dt;
				t1 = c->ts[nr[r0].j] + dt;
				pd[p] = (double)(r0 - a + 1) * dt;
				if (tb[v] >= t0 && tb[v] <= t1 && first)
				{
					/* (the first touch's run: its own time and source) */
					ps[p] = src[v];
					memcpy(pw + 3 * (size_t)p, w + 3 * (size_t)v, sizeof(double) * 3);
					pt[p] = tb[v];
					first = 0;
				}
				else if ((st = swi_point_source(c, X + 3 * (size_t)v, tt, nr[b].j, ps + p, pw + 3 * (size_t)p)) != QAWS_STATUS_OK)
					goto done;
				else
					pt[p] = tt;
				p++;
				a = r0 + 1;
			}
	}
	*out_po = po; *out_src = ps; *out_w = pw; *out_t = pt; *out_dw = pd;
	po = ps = NULL;
	pw = pt = pd = NULL;
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(idx); qaws_internal_free(nr); qaws_internal_free(po); qaws_internal_free(ps); qaws_internal_free(pw); qaws_internal_free(pt);
	qaws_internal_free(pd);
	return st;
}

/* a corner's source from its vertex's (triangle f, weights wv): on f, or, when the point lies on an edge or a corner of the
   shape, on the triangle there whose normal (turned by R) is nearest the corner's own triangle's (no) */
static void swi_corner_source(swi_ctx const* c, double const* tn, unsigned int const* vto, unsigned int const* vtl, unsigned int f, double const* wv,
	double const* R, double const* no, unsigned int* out_f, double* out_w)
{
	unsigned int best = f, on = 0, lead = 0, j, k;
	double bd = -2, u[3], q[3], l;
	for (k = 0; k < 3; k++)
		if (wv[k] < 1e-6)
			on = 1;
		else if (wv[k] > wv[lead])
			lead = k;
	if (on)
	{
		/* the triangles round the source's leading vertex holding every vertex the source weighs */
		unsigned int g = c->t[3 * (size_t)f + lead];
		for (j = vto[g]; j < vto[g + 1]; j++)
		{
			unsigned int t2 = vtl[j], m2, hold = 1;
			double d;
			for (k = 0; k < 3; k++)
				if (wv[k] >= 1e-6)
				{
					unsigned int vk = c->t[3 * (size_t)f + k];
					for (m2 = 0; m2 < 3 && c->t[3 * (size_t)t2 + m2] != vk; m2++)
						;
					hold &= m2 < 3;
				}
			if (!hold)
				continue;
			for (k = 0; k < 3; k++)
				u[k] = R[3 * k] * tn[3 * (size_t)t2] + R[3 * k + 1] * tn[3 * (size_t)t2 + 1] + R[3 * k + 2] * tn[3 * (size_t)t2 + 2];
			l = sqrt(no[0] * no[0] + no[1] * no[1] + no[2] * no[2]);
			d = l > 0 ? (u[0] * no[0] + u[1] * no[1] + u[2] * no[2]) / l : 0;
			if (d > bd)
			{
				bd = d;
				best = t2;
			}
		}
	}
	*out_f = best;
	if (best == f)
		memcpy(out_w, wv, sizeof(double) * 3);
	else
	{
		/* (the source point, at rest, onto that triangle) */
		q[0] = q[1] = q[2] = 0;
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				q[k] += wv[j] * (double)c->x[3 * (size_t)c->t[3 * (size_t)f + j] + k];
		swi_weights(c, best, q, out_w);
	}
}

/* the mesh's attributes onto the result: each vertex's source, the shape's closest point at its best time (and with passes
   wanted, at each later pass); each corner takes it through swi_corner_source: one face's values on each side of a seam,
   the same point (so the same values) where they are continuous */
static qaws_status swi_attributes(swi_ctx* c, qaws_mesh_3d* out, qaws_mesh_3d const* mesh, qaws_sweep_desc const* desc)
{
	unsigned int generate = desc ? desc->generate : 0;
	unsigned int nvo = qaws_mesh_3d_get_vertex_count(out), nto = qaws_mesh_3d_get_triangle_count(out), i, k;
	qaws_scalar* xo = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nvo + 1));
	unsigned int* tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nto + 1));
	double* X = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nvo + 1));
	unsigned int* src = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nvo + 1));
	double* w = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nvo + 1));
	double* tb = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)nvo + 1));
	double* tn = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)c->nt + 1));
	unsigned int *vto = (unsigned int*)qaws_internal_scratch_zalloc((size_t)c->nv + 2, sizeof(unsigned int)), *vtl = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)c->nt + 1));
	unsigned int *vpo = NULL, *vsrc = NULL, *cpo = NULL, *cs = NULL;
	double *vw = NULL, *vt = NULL, *vdw = NULL, *cw = NULL, *cR = NULL, *ct = NULL, *cdw = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!xo || !tri || !X || !src || !w || !tb || !tn || !vto || !vtl)
		goto done;
	qaws_mesh_3d_get_positions(xo, out);
	qaws_mesh_3d_get_triangles(tri, out);
	for (i = 0; i < 3 * nvo; i++)
		X[i] = (double)xo[i];
	/* the mesh's triangles round each vertex, and their normals at rest */
	for (i = 0; i < 3 * c->nt; i++)
		vto[c->t[i] + 1]++;
	for (i = 0; i < c->nv; i++)
		vto[i + 1] += vto[i];
	for (i = 0; i < 3 * c->nt; i++)
		vtl[vto[c->t[i]]++] = i / 3;
	for (i = c->nv; i > 0; i--)
		vto[i] = vto[i - 1];
	vto[0] = 0;
	for (i = 0; i < c->nt; i++)
	{
		double u[3], v[3], l;
		double* n = tn + 3 * (size_t)i;
		for (k = 0; k < 3; k++)
		{
			u[k] = (double)c->x[3 * (size_t)c->t[3 * (size_t)i + 1] + k] - (double)c->x[3 * (size_t)c->t[3 * (size_t)i] + k];
			v[k] = (double)c->x[3 * (size_t)c->t[3 * (size_t)i + 2] + k] - (double)c->x[3 * (size_t)c->t[3 * (size_t)i] + k];
		}
		n[0] = u[1] * v[2] - u[2] * v[1]; n[1] = u[2] * v[0] - u[0] * v[2]; n[2] = u[0] * v[1] - u[1] * v[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (k = 0; k < 3; k++)
			n[k] = l > 0 ? n[k] / l : 0;
	}
	if ((st = swi_sources(c, X, nvo, src, w, tb)) != QAWS_STATUS_OK)
		goto done;
	/* each vertex's passes: all of them when wanted, else its first touch */
	if (swp_passes_wanted(desc))
	{
		if ((st = swi_passes(c, X, nvo, tb, src, w, &vpo, &vsrc, &vw, &vt, &vdw)) != QAWS_STATUS_OK)
			goto done;
	}
	else
	{
		vpo = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nvo + 1));
		if (!vpo)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (i = 0; i <= nvo; i++)
			vpo[i] = i;
	}
	/* the corners' passes: their vertex's, each source taken to the corner's side */
	cpo = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nto + 1) + sizeof(unsigned int));
	if (!cpo)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0, k = 0; i < 3 * nto; i++)
	{
		cpo[i] = k;
		k += vpo[tri[i] + 1] - vpo[tri[i]];
	}
	cpo[3 * (size_t)nto] = k;
	cs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)k + 1));
	cw = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)k + 1));
	cR = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 9, (size_t)k + 1));
	ct = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)k + 1));
	cdw = (double*)qaws_internal_scratch_zalloc((size_t)k + 1, sizeof(double));
	if (!cs || !cw || !cR || !ct || !cdw)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < 3 * nto; i++)
	{
		unsigned int v = tri[i], pp;
		double no[3];
		{
			unsigned int const* tv = tri + 3 * (size_t)(i / 3);
			double a[3], b[3];
			for (k = 0; k < 3; k++)
			{
				a[k] = (double)xo[3 * (size_t)tv[1] + k] - (double)xo[3 * (size_t)tv[0] + k];
				b[k] = (double)xo[3 * (size_t)tv[2] + k] - (double)xo[3 * (size_t)tv[0] + k];
			}
			no[0] = a[1] * b[2] - a[2] * b[1]; no[1] = a[2] * b[0] - a[0] * b[2]; no[2] = a[0] * b[1] - a[1] * b[0];
		}
		for (pp = vpo[v]; pp < vpo[v + 1]; pp++)
		{
			unsigned int q = cpo[i] + (pp - vpo[v]), f = vsrc ? vsrc[pp] : src[v];
			double const* wv = vw ? vw + 3 * (size_t)pp : w + 3 * (size_t)v;
			double R[12], t = vt ? vt[pp] : tb[v];
			if (c->pv)
				for (k = 0; k < 9; k++)
					R[k] = k % 4 == 0;
			else if ((st = swi_pose(c->m, t, R)) != QAWS_STATUS_OK)
				goto done;
			memcpy(cR + 9 * (size_t)q, R, sizeof(double) * 9);
			swi_corner_source(c, tn, vto, vtl, f, wv, R, no, cs + q, cw + 3 * (size_t)q);
			ct[q] = t;
			cdw[q] = vdw ? vdw[pp] : 0;
		}
	}
	st = swp_carry(out, mesh, c->nt, cpo, cs, cw, c->pv ? NULL : cR, 0, ct, desc);
	if (st == QAWS_STATUS_OK && generate)
		st = swp_generate(out, mesh, c->m, generate, c->nt, cpo, cs, cw, ct, cdw, desc);
done:
	qaws_internal_free(xo); qaws_internal_free(tri); qaws_internal_free(X); qaws_internal_free(src); qaws_internal_free(w); qaws_internal_free(tb);
	qaws_internal_free(tn); qaws_internal_free(vto); qaws_internal_free(vtl); qaws_internal_free(vpo); qaws_internal_free(vsrc); qaws_internal_free(vw);
	qaws_internal_free(vt); qaws_internal_free(vdw); qaws_internal_free(cpo); qaws_internal_free(cs); qaws_internal_free(cw); qaws_internal_free(cR);
	qaws_internal_free(ct); qaws_internal_free(cdw);
	return st;
}

/* the context for a mesh under a motion: its triangles, positions and BVH buffers; the box of every position (65 times, grown by
   the most a vertex moves between them) and that speed (half again) */
static qaws_status swi_open(swi_ctx* c, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_sweep_desc const* desc, double* lo, double* hi, double* vmax)
{
	unsigned int i, j, k;
	double slack;
	qaws_scalar* vel = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(c, 0, sizeof *c);
	c->m = motion;
	c->pv = qaws_motion_is_per_vertex(motion);
	c->nv = qaws_mesh_3d_get_vertex_count(mesh);
	c->nt = qaws_mesh_3d_get_triangle_count(mesh);
	c->ex = desc ? desc->executor : NULL;
	c->cur = ~0u;
	if (!c->nt || (c->pv && qaws_motion_get_vertex_count(motion) != c->nv) || swp_tracks_check(desc, mesh) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_ARGUMENT;
	c->t = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), 3 * (size_t)c->nt + 1));
	c->x = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), 3 * (size_t)c->nv + 1));
	c->buf = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), 3 * (size_t)c->nv + 1));
	vel = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), 3 * (size_t)c->nv + 1));
	c->q = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 4 * SWI_CHUNK);
	c->h = (qaws_mesh_hit*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_hit) * SWI_CHUNK);
	c->best = (double*)qaws_internal_scratch_malloc(sizeof(double) * 8 * SWI_CHUNK);
	c->bj = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * 3 * SWI_CHUNK);
	if (!c->t || !c->x || !c->buf || !vel || !c->q || !c->h || !c->best || !c->bj)
		goto done;
	c->d = c->q + 3 * SWI_CHUNK;
	c->bt = c->best + SWI_CHUNK;
	c->hp = c->best + 2 * SWI_CHUNK;
	c->hw = c->best + 5 * SWI_CHUNK;
	c->idx = c->bj + SWI_CHUNK;
	c->htri = c->bj + 2 * SWI_CHUNK;
	qaws_mesh_3d_get_triangles(c->t, mesh);
	qaws_mesh_3d_get_positions(c->x, mesh);
	for (k = 0; k < 3; k++)
	{
		lo[k] = HUGE_VAL;
		hi[k] = -HUGE_VAL;
	}
	*vmax = 0;
	for (i = 0; i <= 64; i++)
	{
		if ((st = qaws_motion_transform_points(c->buf, vel, motion, (qaws_scalar)(i / 64.0), c->pv ? NULL : c->x, c->nv)) != QAWS_STATUS_OK)
			goto done;
		for (j = 0; j < c->nv; j++)
		{
			*vmax = fmax(*vmax, sqrt((double)vel[3 * j] * vel[3 * j] + (double)vel[3 * j + 1] * vel[3 * j + 1] + (double)vel[3 * j + 2] * vel[3 * j + 2]));
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], (double)c->buf[3 * (size_t)j + k]);
				hi[k] = fmax(hi[k], (double)c->buf[3 * (size_t)j + k]);
			}
		}
	}
	*vmax *= 1.5;
	slack = *vmax / 128;
	for (k = 0; k < 3; k++)
	{
		lo[k] -= slack;
		hi[k] += slack;
	}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(vel);
	return st;
}

/* the samples: no vertex moving past two cells between two (refined between them; a vertex motion: a quarter cell, its values
   not refined), at most maxs; their poses and bounding spheres; the BVH; ext the box's size */
static qaws_status swi_times(swi_ctx* c, double cell, double vmax, unsigned int maxs, double ext)
{
	unsigned int i, j, k;
	double step = c->pv ? cell / 4 : 2 * cell;
	qaws_status st;
	c->ns = vmax > 0 ? (unsigned int)fmin((double)maxs, ceil(vmax / step)) + 1 : 1;
	c->ts = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), c->ns));
	c->cen = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 4, c->ns));
	c->pose = c->pv ? NULL : (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 12, c->ns));
	if (!c->ts || !c->cen || (!c->pv && !c->pose))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (j = 0; j < c->ns; j++)
	{
		double* s = c->cen + 4 * (size_t)j;
		c->ts[j] = c->ns > 1 ? (double)j / (c->ns - 1) : 0;
		if ((st = qaws_motion_transform_points(c->buf, NULL, c->m, (qaws_scalar)c->ts[j], c->pv ? NULL : c->x, c->nv)) != QAWS_STATUS_OK ||
			(!c->pv && (st = swi_pose(c->m, c->ts[j], c->pose + 12 * (size_t)j)) != QAWS_STATUS_OK))
			return st;
		s[0] = s[1] = s[2] = s[3] = 0;
		for (i = 0; i < c->nv; i++)
			for (k = 0; k < 3; k++)
				s[k] += (double)c->buf[3 * (size_t)i + k] / c->nv;
		for (i = 0; i < c->nv; i++)
		{
			double e = 0;
			for (k = 0; k < 3; k++)
				e += ((double)c->buf[3 * (size_t)i + k] - s[k]) * ((double)c->buf[3 * (size_t)i + k] - s[k]);
			s[3] = fmax(s[3], sqrt(e));
		}
	}
	{
		qaws_mesh_bvh_desc bd;
		memset(&bd, 0, sizeof bd);
		bd.positions = c->x;
		bd.vertex_count = c->nv;
		bd.dimension = 3;
		bd.simplices = c->t;
		bd.simplex_size = 3;
		bd.simplex_count = c->nt;
		if ((st = qaws_mesh_bvh_create(&c->b, &bd)) != QAWS_STATUS_OK)
			return st;
	}
	c->band = 4 * cell;
	c->cell = cell;
	c->tie = 1e-9 * ext;
	return QAWS_STATUS_OK;
}

static void swi_close(swi_ctx* c)
{
	qaws_mesh_bvh_destroy(c->b);
	qaws_internal_free(c->t); qaws_internal_free(c->x); qaws_internal_free(c->buf); qaws_internal_free(c->q);
	qaws_internal_free(c->h); qaws_internal_free(c->best); qaws_internal_free(c->bj); qaws_internal_free(c->ts); qaws_internal_free(c->cen);
	qaws_internal_free(c->pose);
}

qaws_status qaws_mesh_3d_sweep_implicit(qaws_mesh_3d** out_mesh, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion,
	qaws_sweep_desc const* desc)
{
	qaws_sweep_report local;
	swi_ctx c;
	wt_contour_info info;
	double lo[3], hi[3], vmax = 0, ext = 0, cell;
	unsigned int k, res = desc && desc->resolution ? desc->resolution : 64, maxs = desc && desc->max_slabs ? desc->max_slabs : 4096;
	qaws_status st;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out_mesh || !mesh || !motion)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	memset(&info, 0, sizeof info);
	if ((st = swi_open(&c, mesh, motion, desc, lo, hi, &vmax)) != QAWS_STATUS_OK)
		goto done;
	for (k = 0; k < 3; k++)
		ext = fmax(ext, hi[k] - lo[k]);
	cell = ext / res;
	if ((st = swi_times(&c, cell, vmax, maxs, ext)) != QAWS_STATUS_OK)
		goto done;
	report->slabs = c.ns;
	report->deviation = (qaws_scalar)(c.ns > 1 ? 0.5 * vmax / (c.ns - 1) : 0);
	st = wt_contour_field(out_mesh, &info, swi_eval, &c, lo, hi, cell, res, 1.0);
	report->closed = info.closed;
	if (st == QAWS_STATUS_OK && *out_mesh && (qaws_mesh_3d_get_attribute_count(mesh) || (desc && desc->generate)))
		st = swi_attributes(&c, *out_mesh, mesh, desc);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(*out_mesh);
		*out_mesh = NULL;
	}
done:
	swi_close(&c);
	return st;
}

qaws_status swp_slab_passes(qaws_mesh_3d* out, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_sweep_desc const* desc, double reach)
{
	swi_ctx c;
	double lo[3], hi[3], vmax = 0, ext = 0;
	unsigned int k, na = qaws_mesh_3d_get_attribute_count(mesh), a;
	qaws_status st;
	if ((st = swi_open(&c, mesh, motion, desc, lo, hi, &vmax)) != QAWS_STATUS_OK)
		goto done;
	for (k = 0; k < 3; k++)
		ext = fmax(ext, hi[k] - lo[k]);
	if ((st = swi_times(&c, reach / 1.5, vmax, desc && desc->max_slabs ? desc->max_slabs : 4096, ext)) != QAWS_STATUS_OK)
		goto done;
	/* (the body's values carried by the slabs taken off: carried again pass by pass; edge values stay the slabs') */
	for (a = 0; a < na; a++)
	{
		qaws_attribute_desc d;
		unsigned int id;
		qaws_mesh_3d_get_attribute_desc(&d, mesh, a);
		if (d.name && strcmp(d.name, "source") && d.element != QAWS_MESH_EDGE && qaws_mesh_3d_find_attribute(&id, out, d.name) == QAWS_STATUS_OK)
			qaws_mesh_3d_remove_attribute(out, id);
	}
	st = swi_attributes(&c, out, mesh, desc);
done:
	swi_close(&c);
	return st;
}
