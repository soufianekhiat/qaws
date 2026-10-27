/*
 * Implicit sweeps (qaws_mesh_sweep.h): the swept solid as the level set of
 * the least over time of the moving shape's signed distance, sampled in time
 * (no vertex moving more than a cell between samples), each point's samples
 * pruned by their bounding spheres, a rigid motion's best time refined by
 * golden sections; contoured by the offsets' dual contouring
 * (wt_contour_field, blocks pruned: the field is 1-Lipschitz).
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
	qaws_batch_executor const* ex;
	/* per chunk */
	qaws_scalar *q, *d;
	qaws_mesh_hit* h;
	double *best, *bt, *hp;         /* the least value, its time, the closest point (world) */
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
		if (!((double)c->d[i] < c->best[p]))
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
	size_t c0;
	unsigned int i, a;
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
			if (l > 1e-12)
				for (a = 0; a < 3; a++)
					g[a] /= c->best[i] < 0 ? -l : l;
			else
			{
				/* on the surface: its triangle's normal at the best time */
				double V[3][3], u[3], w[3];
				unsigned int const* tv = c->t + 3 * (size_t)c->htri[i];
				if (c->pv)
				{
					if ((st = swi_refit(c, c->bj[i])) != QAWS_STATUS_OK)
						break;
					for (a = 0; a < 9; a++)
						V[a / 3][a % 3] = (double)c->buf[3 * (size_t)tv[a / 3] + a % 3];
				}
				else
				{
					double R[12];
					unsigned int b;
					if ((st = swi_pose(c->m, c->bt[i], R)) != QAWS_STATUS_OK)
						break;
					for (a = 0; a < 3; a++)
						for (b = 0; b < 3; b++)
							V[a][b] = R[3 * b] * (double)c->x[3 * (size_t)tv[a]] + R[3 * b + 1] * (double)c->x[3 * (size_t)tv[a] + 1] +
								R[3 * b + 2] * (double)c->x[3 * (size_t)tv[a] + 2];
				}
				for (a = 0; a < 3; a++)
				{
					u[a] = V[1][a] - V[0][a];
					w[a] = V[2][a] - V[0][a];
				}
				g[0] = u[1] * w[2] - u[2] * w[1]; g[1] = u[2] * w[0] - u[0] * w[2]; g[2] = u[0] * w[1] - u[1] * w[0];
				l = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
				for (a = 0; a < 3; a++)
					g[a] = l > 0 ? g[a] / l : 0;
			}
			for (a = 0; a < 3; a++)
				grad[3 * (c0 + i) + a] = g[a];
		}
	}
	qaws_internal_free(P); qaws_internal_free(tt); qaws_internal_free(G); qaws_internal_free(jj);
	return st;
}

qaws_status qaws_mesh_3d_sweep_implicit(qaws_mesh_3d** out_mesh, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion,
	qaws_sweep_desc const* desc)
{
	qaws_sweep_report local;
	swi_ctx c;
	wt_contour_info info;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, vmax = 0, ext = 0, cell, step, slack;
	unsigned int i, j, k, res = desc && desc->resolution ? desc->resolution : 64, maxs = desc && desc->max_slabs ? desc->max_slabs : 4096;
	qaws_scalar* vel = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out_mesh || !mesh || !motion)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	memset(&c, 0, sizeof c);
	memset(&info, 0, sizeof info);
	c.m = motion;
	c.pv = qaws_motion_is_per_vertex(motion);
	c.nv = qaws_mesh_3d_get_vertex_count(mesh);
	c.nt = qaws_mesh_3d_get_triangle_count(mesh);
	c.ex = desc ? desc->executor : NULL;
	c.cur = ~0u;
	if (!c.nt || (c.pv && qaws_motion_get_vertex_count(motion) != c.nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	c.t = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)c.nt + 1));
	c.x = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)c.nv + 1));
	c.buf = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)c.nv + 1));
	vel = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)c.nv + 1));
	c.q = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * 4 * SWI_CHUNK);
	c.h = (qaws_mesh_hit*)qaws_internal_malloc(sizeof(qaws_mesh_hit) * SWI_CHUNK);
	c.best = (double*)qaws_internal_malloc(sizeof(double) * 5 * SWI_CHUNK);
	c.bj = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * 3 * SWI_CHUNK);
	if (!c.t || !c.x || !c.buf || !vel || !c.q || !c.h || !c.best || !c.bj)
		goto done;
	c.d = c.q + 3 * SWI_CHUNK;
	c.bt = c.best + SWI_CHUNK;
	c.hp = c.best + 2 * SWI_CHUNK;
	c.idx = c.bj + SWI_CHUNK;
	c.htri = c.bj + 2 * SWI_CHUNK;
	qaws_mesh_3d_get_triangles(c.t, mesh);
	qaws_mesh_3d_get_positions(c.x, mesh);
	/* the fastest vertex and the box of every position (65 times, grown by the most a vertex moves between them) */
	for (i = 0; i <= 64; i++)
	{
		if ((st = qaws_motion_transform_points(c.buf, vel, motion, (qaws_scalar)(i / 64.0), c.pv ? NULL : c.x, c.nv)) != QAWS_STATUS_OK)
			goto done;
		for (j = 0; j < c.nv; j++)
		{
			vmax = fmax(vmax, sqrt((double)vel[3 * j] * vel[3 * j] + (double)vel[3 * j + 1] * vel[3 * j + 1] + (double)vel[3 * j + 2] * vel[3 * j + 2]));
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], (double)c.buf[3 * (size_t)j + k]);
				hi[k] = fmax(hi[k], (double)c.buf[3 * (size_t)j + k]);
			}
		}
	}
	vmax *= 1.5;
	slack = vmax / 128;
	for (k = 0; k < 3; k++)
	{
		lo[k] -= slack;
		hi[k] += slack;
		ext = fmax(ext, hi[k] - lo[k]);
	}
	cell = ext / res;
	/* the samples: no vertex moving past two cells between two (refined between them; a vertex motion: a quarter cell, its
	   values not refined) */
	step = c.pv ? cell / 4 : 2 * cell;
	c.ns = vmax > 0 ? (unsigned int)fmin((double)maxs, ceil(vmax / step)) + 1 : 1;
	report->slabs = c.ns;
	report->deviation = (qaws_scalar)(c.ns > 1 ? 0.5 * vmax / (c.ns - 1) : 0);
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	c.ts = (double*)qaws_internal_malloc(sizeof(double) * c.ns);
	c.cen = (double*)qaws_internal_malloc(sizeof(double) * 4 * (size_t)c.ns);
	c.pose = c.pv ? NULL : (double*)qaws_internal_malloc(sizeof(double) * 12 * (size_t)c.ns);
	if (!c.ts || !c.cen || (!c.pv && !c.pose))
		goto done;
	for (j = 0; j < c.ns; j++)
	{
		double* s = c.cen + 4 * (size_t)j;
		c.ts[j] = c.ns > 1 ? (double)j / (c.ns - 1) : 0;
		if ((st = qaws_motion_transform_points(c.buf, NULL, motion, (qaws_scalar)c.ts[j], c.pv ? NULL : c.x, c.nv)) != QAWS_STATUS_OK ||
			(!c.pv && (st = swi_pose(motion, c.ts[j], c.pose + 12 * (size_t)j)) != QAWS_STATUS_OK))
			goto done;
		s[0] = s[1] = s[2] = s[3] = 0;
		for (i = 0; i < c.nv; i++)
			for (k = 0; k < 3; k++)
				s[k] += (double)c.buf[3 * (size_t)i + k] / c.nv;
		for (i = 0; i < c.nv; i++)
		{
			double e = 0;
			for (k = 0; k < 3; k++)
				e += ((double)c.buf[3 * (size_t)i + k] - s[k]) * ((double)c.buf[3 * (size_t)i + k] - s[k]);
			s[3] = fmax(s[3], sqrt(e));
		}
	}
	{
		qaws_mesh_bvh_desc bd;
		memset(&bd, 0, sizeof bd);
		bd.positions = c.x;
		bd.vertex_count = c.nv;
		bd.dimension = 3;
		bd.simplices = c.t;
		bd.simplex_size = 3;
		bd.simplex_count = c.nt;
		if ((st = qaws_mesh_bvh_create(&c.b, &bd)) != QAWS_STATUS_OK)
			goto done;
	}
	c.band = 4 * cell;
	st = wt_contour_field(out_mesh, &info, swi_eval, &c, lo, hi, cell, res, 1.0);
	report->closed = info.closed;
done:
	qaws_mesh_bvh_destroy(c.b);
	qaws_internal_free(c.t); qaws_internal_free(c.x); qaws_internal_free(c.buf); qaws_internal_free(vel); qaws_internal_free(c.q);
	qaws_internal_free(c.h); qaws_internal_free(c.best); qaws_internal_free(c.bj); qaws_internal_free(c.ts); qaws_internal_free(c.cen);
	qaws_internal_free(c.pose);
	return st;
}
