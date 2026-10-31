#include "qaws_surface_batch.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_eval.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_broadphase.h"
#include "internal/qaws_internal_flatten.h"
#include "internal/qaws_internal_batch.h"
#include "internal/qaws_internal_parallel.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* split by tools/splitc.py: shared between qaws_surface_batch.c and its moved part */
#define sb_dot QAWS_SFX(qaws_internal_surface_batch_sb_dot)

/*
 * Batched curve / surface intersection.
 *
 *   1. flatten every curve into chord segments and every surface into
 *      (u, v) patches (qaws_internal_flatten_*), boxes inflated by twice the
 *      measured deviation;
 *   2. one grid over all pieces (qaws_internal_broadphase); only segment /
 *      patch pairs are kept;
 *   3. the segment against the patch's two triangles, with a slack scaled by
 *      both inflations, seeds (t, u, v);
 *   4. Newton on S(u, v) - C(t) = 0, polished to rounding level;
 *   5. coincident hits of one curve / surface pair are kept once.
 */

#define SB_NEWTON_ITERS 24

#if QAWS_SCALAR_IS_FLOAT
#define SB_POS_REL ((qaws_scalar)2e-5)
#define SB_PAR_REL ((qaws_scalar)1e-4)
#define SB_STEP_REL ((qaws_scalar)1e-7)   /* Newton steps below it: converged */
#else
#define SB_POS_REL ((qaws_scalar)1e-10)
#define SB_PAR_REL ((qaws_scalar)1e-7)
#define SB_STEP_REL ((qaws_scalar)1e-15)
#endif

qaws_scalar sb_dot(qaws_scalar const* a, qaws_scalar const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* ------------------------------------------------------------------ */
/*  Prepared surface sets                                              */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_set_create(qaws_surface_set** out_set, qaws_surface_batch_desc const* desc)
{
	qaws_surface_set* set;
	unsigned int i, n, cap = 0;
	qaws_status s = QAWS_STATUS_OK;
	if (!desc || !out_set || (desc->surface_count && !desc->surfaces))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_set = NULL;
	n = desc->surface_count;
	for (i = 0; i < n; i++)
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	set = (qaws_surface_set*)qaws_internal_zalloc(1, sizeof(qaws_surface_set));
	if (!set)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	set->desc = *desc;
	set->desc.surfaces = (qaws_surface const* const*)qaws_internal_malloc(qaws_internal_mul2((n ? n : 1), sizeof(qaws_surface*)));
	set->desc.families = desc->families ? (unsigned int const*)qaws_internal_malloc(qaws_internal_mul2((n ? n : 1), sizeof(unsigned int))) : NULL;
	if (!set->desc.surfaces || (desc->families && !set->desc.families))
	{
		qaws_surface_set_destroy(set);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy((void*)set->desc.surfaces, desc->surfaces, n * sizeof(qaws_surface*));
	if (desc->families)
		memcpy((void*)set->desc.families, desc->families, n * sizeof(unsigned int));
	set->ext = qaws_internal_flatten_extent(NULL, 0, 3, desc->surfaces, n);
	set->flat = desc->flatness > 0 ? desc->flatness : set->ext / 512;
	for (i = 0; i < n && s == QAWS_STATUS_OK; i++)
		s = qaws_internal_flatten_surface(&set->patches, &set->npatch, &cap, desc->surfaces[i], set->flat, i);
	if (s != QAWS_STATUS_OK)
	{
		qaws_surface_set_destroy(set);
		return s;
	}
	*out_set = set;
	return QAWS_STATUS_OK;
}

void qaws_surface_set_destroy(qaws_surface_set* set)
{
	if (!set)
		return;
	qaws_internal_free((void*)set->desc.surfaces);
	qaws_internal_free((void*)set->desc.families);
	qaws_internal_free(set->patches);
	qaws_internal_free(set);
}

unsigned int qaws_surface_set_get_patch_count(qaws_surface_set const* set)
{
	return set ? set->npatch : 0;
}

#define SS_MARCH_STEPS 96

/* ------------------------------------------------------------------ */
/*  Closest points on surfaces                                         */
/* ------------------------------------------------------------------ */

#define SC_ITERS 32

/* nearest point of triangle (a, b, c) to p: barycentrics w (Ericson's regions) */
static void sc_triangle(qaws_scalar const* p, qaws_scalar const* a, qaws_scalar const* b, qaws_scalar const* c, qaws_scalar* w)
{
	qaws_scalar ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, den, v, t;
	int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = b[k] - a[k];
		ac[k] = c[k] - a[k];
		ap[k] = p[k] - a[k];
		bp[k] = p[k] - b[k];
		cp[k] = p[k] - c[k];
	}
	d1 = sb_dot(ab, ap); d2 = sb_dot(ac, ap);
	if (d1 <= 0 && d2 <= 0) { w[0] = 1; w[1] = w[2] = 0; return; }
	d3 = sb_dot(ab, bp); d4 = sb_dot(ac, bp);
	if (d3 >= 0 && d4 <= d3) { w[1] = 1; w[0] = w[2] = 0; return; }
	vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0)
	{
		v = d1 / (d1 - d3);
		w[0] = 1 - v; w[1] = v; w[2] = 0;
		return;
	}
	d5 = sb_dot(ab, cp); d6 = sb_dot(ac, cp);
	if (d6 >= 0 && d5 <= d6) { w[2] = 1; w[0] = w[1] = 0; return; }
	vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0)
	{
		t = d2 / (d2 - d6);
		w[0] = 1 - t; w[2] = t; w[1] = 0;
		return;
	}
	va = d3 * d6 - d5 * d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
	{
		t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		w[1] = 1 - t; w[2] = t; w[0] = 0;
		return;
	}
	den = 1 / (va + vb + vc);
	w[1] = vb * den;
	w[2] = vc * den;
	w[0] = 1 - w[1] - w[2];
}

/* Newton on (S - p) . S_u = (S - p) . S_v = 0 from (u, v), clamped */
static qaws_scalar sc_refine(qaws_surface const* s, qaws_scalar const* p, qaws_scalar* u, qaws_scalar* v, qaws_scalar* pos)
{
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	unsigned int it, flags = QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV;
	qaws_surface_eval_result r;
	qaws_scalar w[3];
	for (it = 0; it < SC_ITERS; it++)
	{
		qaws_scalar su[3], sv[3], suu[3], suv[3], svv[3], g0, g1, h00, h01, h11, det, du, dv;
		if (qaws_surface_evaluate(&r, s, *u, *v, flags) != QAWS_STATUS_OK)
			break;
		w[0] = r.position.x - p[0]; w[1] = r.position.y - p[1]; w[2] = r.position.z - p[2];
		su[0] = r.du.x; su[1] = r.du.y; su[2] = r.du.z;
		sv[0] = r.dv.x; sv[1] = r.dv.y; sv[2] = r.dv.z;
		suu[0] = r.duu.x; suu[1] = r.duu.y; suu[2] = r.duu.z;
		suv[0] = r.duv.x; suv[1] = r.duv.y; suv[2] = r.duv.z;
		svv[0] = r.dvv.x; svv[1] = r.dvv.y; svv[2] = r.dvv.z;
		g0 = sb_dot(w, su);
		g1 = sb_dot(w, sv);
		h00 = sb_dot(su, su) + sb_dot(w, suu);
		h01 = sb_dot(su, sv) + sb_dot(w, suv);
		h11 = sb_dot(sv, sv) + sb_dot(w, svv);
		det = h00 * h11 - h01 * h01;
		if (!(det > 0 && h00 > 0))
		{
			/* not convex here: the Gauss-Newton part only */
			h00 = sb_dot(su, su);
			h01 = sb_dot(su, sv);
			h11 = sb_dot(sv, sv);
			det = h00 * h11 - h01 * h01;
			if (!(det > 0))
				break;
		}
		/* active set: a parameter on its bound with the descent pointing out
		   stays there and the other one takes the 1D step along that edge */
		{
			int hold_u = (*u <= ur.min_value && g0 > 0) || (*u >= ur.max_value && g0 < 0);
			int hold_v = (*v <= vr.min_value && g1 > 0) || (*v >= vr.max_value && g1 < 0);
			if (hold_u && hold_v)
				break;   /* a corner where the distance grows inward */
			if (hold_u || hold_v)
			{
				qaws_scalar e = hold_v ? sb_dot(su, su) + sb_dot(w, suu) : sb_dot(sv, sv) + sb_dot(w, svv);
				if (!(e > 0))
					e = hold_v ? sb_dot(su, su) : sb_dot(sv, sv);
				if (!(e > 0))
					break;
				du = hold_v ? g0 / e : 0;
				dv = hold_u ? g1 / e : 0;
				goto step;
			}
		}
		du = (h11 * g0 - h01 * g1) / det;
		dv = (h00 * g1 - h01 * g0) / det;
		/* a parameter the step would take out of the domain stops on its
		   bound and the other one takes the 1D step along that edge */
		if ((*u - du < ur.min_value || *u - du > ur.max_value) && h11 > 0)
		{
			qaws_scalar ub = *u - du < ur.min_value ? ur.min_value : ur.max_value;
			du = *u - ub;
			dv = g1 / h11;
		}
		else if ((*v - dv < vr.min_value || *v - dv > vr.max_value) && h00 > 0)
		{
			qaws_scalar vb = *v - dv < vr.min_value ? vr.min_value : vr.max_value;
			dv = *v - vb;
			du = g0 / h00;
		}
	step:
		*u -= du;
		*v -= dv;
		if (*u < ur.min_value) *u = ur.min_value;
		if (*u > ur.max_value) *u = ur.max_value;
		if (*v < vr.min_value) *v = vr.min_value;
		if (*v > vr.max_value) *v = vr.max_value;
		if (fabs(du) <= (ur.max_value - ur.min_value) * SB_STEP_REL && fabs(dv) <= (vr.max_value - vr.min_value) * SB_STEP_REL)
			break;
	}
	qaws_surface_evaluate(&r, s, *u, *v, QAWS_SURFACE_EVAL_POSITION);
	pos[0] = r.position.x; pos[1] = r.position.y; pos[2] = r.position.z;
	w[0] = pos[0] - p[0]; w[1] = pos[1] - p[1]; w[2] = pos[2] - p[2];
	return (qaws_scalar)sqrt(sb_dot(w, w));
}

typedef struct sc_cand
{
	unsigned int patch;
	qaws_scalar lb;
} sc_cand;

typedef struct sc_ctx
{
	qaws_flat_patch const* patches;
	unsigned int* stamp;
	unsigned int query;
	qaws_scalar const* p;
	qaws_scalar best_ub;
	sc_cand* cand;
	unsigned int ncand, capcand;
	int failed;
} sc_ctx;

static void sc_visit(void* user, unsigned int item)
{
	sc_ctx* x = (sc_ctx*)user;
	qaws_flat_patch const* P = &x->patches[item];
	qaws_scalar lb = 0, e;
	int k, i;
	if (x->stamp[item] == x->query)
		return;
	x->stamp[item] = x->query;
	for (k = 0; k < 3; k++)
	{
		e = x->p[k] < P->lo[k] ? P->lo[k] - x->p[k] : (x->p[k] > P->hi[k] ? x->p[k] - P->hi[k] : 0);
		lb += e * e;
	}
	lb = (qaws_scalar)sqrt(lb);
	if (lb > x->best_ub)
		return;
	/* the corners lie on the surface: an upper bound */
	for (i = 0; i < 4; i++)
	{
		qaws_scalar d[3], dd;
		for (k = 0; k < 3; k++)
			d[k] = P->p[i][k] - x->p[k];
		dd = (qaws_scalar)sqrt(sb_dot(d, d));
		if (dd < x->best_ub)
			x->best_ub = dd;
	}
	if (x->ncand == x->capcand)
	{
		unsigned int cap = qaws_internal_grow(x->capcand, (size_t)x->ncand + 1, 64);
		sc_cand* g = cap ? (sc_cand*)qaws_internal_realloc(x->cand, qaws_internal_mul2(cap, sizeof(sc_cand))) : NULL;
		if (!g)
		{
			x->failed = 1;
			return;
		}
		x->cand = g;
		x->capcand = cap;
	}
	x->cand[x->ncand].patch = item;
	x->cand[x->ncand].lb = lb;
	x->ncand++;
}

static int sc_cmp(void const* a, void const* b)
{
	qaws_scalar x = ((sc_cand const*)a)->lb, y = ((sc_cand const*)b)->lb;
	return x < y ? -1 : (x > y ? 1 : 0);
}

#define SC_GRAIN 32

typedef struct sc_job
{
	qaws_surface const* const* surfaces;
	qaws_flat_patch const* patches;
	unsigned int npatch;
	qaws_bp_grid* grid;
	qaws_scalar const* points;
	qaws_scalar max_distance;
	qaws_surface_batch_closest* out;
	qaws_surface_batch_stats* stats;   /* one per chunk */
} sc_job;

static qaws_status sc_chunk(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	static int const tri[2][3] = { { 0, 1, 3 }, { 0, 3, 2 } };
	sc_job const* J = (sc_job const*)ctx;
	qaws_surface_batch_stats* stats = &J->stats[chunk];
	sc_ctx x;
	unsigned int i, k;
	qaws_status s = QAWS_STATUS_OK;
	memset(&x, 0, sizeof(x));
	x.stamp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(J->npatch, sizeof(unsigned int)));
	if (!x.stamp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < J->npatch; i++)
		x.stamp[i] = ~0u;
	x.patches = J->patches;
	for (i = begin; i < end && s == QAWS_STATUS_OK; i++)
	{
		qaws_scalar const* p = J->points + 3 * i;
		qaws_scalar best = (qaws_scalar)HUGE_VAL, bu = 0, bv = 0, bp[3] = { 0, 0, 0 };
		unsigned int r, c, bsurf = QAWS_CURVE_BATCH_NONE;
		double pd[3];
		x.query = i;
		x.p = p;
		x.best_ub = J->max_distance > 0 ? J->max_distance : (qaws_scalar)HUGE_VAL;
		x.ncand = 0;
		pd[0] = (double)p[0]; pd[1] = (double)p[1]; pd[2] = (double)p[2];
		for (r = 0;; r++)
		{
			double b = qaws_internal_grid_ring(J->grid, pd, r, sc_visit, &x);
			if (x.failed)
			{
				s = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			if (b == HUGE_VAL || b > (double)x.best_ub)
				break;
		}
		if (s != QAWS_STATUS_OK)
			break;
		/* nearest box first, while a box can beat the best refined point */
		qsort(x.cand, x.ncand, sizeof(sc_cand), sc_cmp);
		for (c = 0; c < x.ncand; c++)
		{
			qaws_flat_patch const* P = &J->patches[x.cand[c].patch];
			qaws_scalar uc[4], vc[4], wb[3], u, v, d, pos[3], seed_d = (qaws_scalar)HUGE_VAL;
			int t, ts = 0;
			qaws_scalar wsel[3] = { 1, 0, 0 };
			if (x.cand[c].lb > best || x.cand[c].lb > x.best_ub)
				break;
			stats->candidate_count++;
			stats->newton_count++;
			uc[0] = P->u0; vc[0] = P->v0; uc[1] = P->u1; vc[1] = P->v0; uc[2] = P->u0; vc[2] = P->v1; uc[3] = P->u1; vc[3] = P->v1;
			/* the seed: the nearest point of the patch's two triangles */
			for (t = 0; t < 2; t++)
			{
				qaws_scalar q[3], dq[3], dd;
				sc_triangle(p, P->p[tri[t][0]], P->p[tri[t][1]], P->p[tri[t][2]], wb);
				for (k = 0; k < 3; k++)
				{
					q[k] = wb[0] * P->p[tri[t][0]][k] + wb[1] * P->p[tri[t][1]][k] + wb[2] * P->p[tri[t][2]][k];
					dq[k] = q[k] - p[k];
				}
				dd = sb_dot(dq, dq);
				if (dd < seed_d)
				{
					seed_d = dd;
					ts = t;
					wsel[0] = wb[0]; wsel[1] = wb[1]; wsel[2] = wb[2];
				}
			}
			u = wsel[0] * uc[tri[ts][0]] + wsel[1] * uc[tri[ts][1]] + wsel[2] * uc[tri[ts][2]];
			v = wsel[0] * vc[tri[ts][0]] + wsel[1] * vc[tri[ts][1]] + wsel[2] * vc[tri[ts][2]];
			d = sc_refine(J->surfaces[P->owner], p, &u, &v, pos);
			if (d < best || (d == best && P->owner < bsurf))
			{
				best = d;
				bsurf = P->owner;
				bu = u;
				bv = v;
				bp[0] = pos[0]; bp[1] = pos[1]; bp[2] = pos[2];
			}
		}
		if (bsurf != QAWS_CURVE_BATCH_NONE && (J->max_distance <= 0 || best <= J->max_distance))
		{
			J->out[i].surface = bsurf;
			J->out[i].u = bu;
			J->out[i].v = bv;
			J->out[i].distance = best;
			J->out[i].position.x = bp[0];
			J->out[i].position.y = bp[1];
			J->out[i].position.z = bp[2];
			stats->hit_count++;
		}
	}
	qaws_internal_free(x.stamp);
	qaws_internal_free(x.cand);
	return s;
}

static qaws_status sc_run(qaws_surface const* const* surfaces, qaws_flat_patch const* patches, unsigned int npatch, qaws_scalar const* points,
	unsigned int point_count, qaws_scalar max_distance, qaws_batch_executor const* executor, qaws_surface_batch_closest* out, qaws_surface_batch_stats* stats)
{
	sc_job J;
	qaws_bp_box* boxes;
	unsigned int i, k, nchunk;
	qaws_status s;
	memset(&J, 0, sizeof(J));
	stats->patch_count = npatch;
	for (i = 0; i < point_count; i++)
	{
		memset(&out[i], 0, sizeof(out[i]));
		out[i].surface = QAWS_CURVE_BATCH_NONE;
	}
	if (!npatch || !point_count)
		return QAWS_STATUS_OK;
	boxes = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(npatch, sizeof(qaws_bp_box)));
	nchunk = qaws_internal_chunk_count(point_count, SC_GRAIN);
	J.stats = (qaws_surface_batch_stats*)qaws_internal_zalloc(nchunk, sizeof(qaws_surface_batch_stats));
	if (!boxes || !J.stats)
	{
		qaws_internal_free(boxes);
		qaws_internal_free(J.stats);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < npatch; i++)
		for (k = 0; k < 3; k++)
		{
			boxes[i].lo[k] = (double)patches[i].lo[k];
			boxes[i].hi[k] = (double)patches[i].hi[k];
		}
	s = qaws_internal_grid_create(&J.grid, boxes, npatch, 3);
	qaws_internal_free(boxes);
	J.surfaces = surfaces;
	J.patches = patches;
	J.npatch = npatch;
	J.points = points;
	J.max_distance = max_distance;
	J.out = out;
	if (s == QAWS_STATUS_OK)
		s = qaws_internal_parallel(point_count, SC_GRAIN, sc_chunk, &J, executor);
	for (k = 0; k < nchunk; k++)
	{
		stats->candidate_count += J.stats[k].candidate_count;
		stats->newton_count += J.stats[k].newton_count;
		stats->hit_count += J.stats[k].hit_count;
	}
	qaws_internal_grid_destroy(J.grid);
	qaws_internal_free(J.stats);
	return s;
}

qaws_status qaws_surface_batch_find_closest(
	qaws_surface_batch_closest* out_points,
	qaws_surface_batch_stats* out_stats,
	qaws_surface_batch_closest_desc const* desc)
{
	qaws_surface_batch_stats st;
	qaws_flat_patch* patches = NULL;
	unsigned int npatch = 0, cap = 0, i;
	qaws_scalar ext, flat;
	qaws_status s = QAWS_STATUS_OK;
	memset(&st, 0, sizeof(st));
	if (out_stats)
		*out_stats = st;
	if (!desc || (desc->surface_count && !desc->surfaces) || (desc->point_count && (!desc->points || !out_points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->surface_count; i++)
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	ext = qaws_internal_flatten_extent(NULL, 0, 3, desc->surfaces, desc->surface_count);
	flat = desc->flatness > 0 ? desc->flatness : ext / 512;
	for (i = 0; i < desc->surface_count && s == QAWS_STATUS_OK; i++)
		s = qaws_internal_flatten_surface(&patches, &npatch, &cap, desc->surfaces[i], flat, i);
	if (s == QAWS_STATUS_OK)
		s = sc_run(desc->surfaces, patches, npatch, desc->points, desc->point_count, desc->max_distance, desc->executor, out_points, &st);
	if (out_stats)
		*out_stats = st;
	qaws_internal_free(patches);
	return s;
}

qaws_status qaws_surface_set_find_closest(
	qaws_surface_batch_closest* out_points,
	qaws_surface_batch_stats* out_stats,
	qaws_surface_set const* set,
	qaws_scalar const* points,
	unsigned int point_count,
	qaws_scalar max_distance)
{
	qaws_surface_batch_stats st;
	qaws_status s;
	memset(&st, 0, sizeof(st));
	if (out_stats)
		*out_stats = st;
	if (!set || (point_count && (!points || !out_points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	s = sc_run(set->desc.surfaces, set->patches, set->npatch, points, point_count, max_distance, set->desc.executor, out_points, &st);
	if (out_stats)
		*out_stats = st;
	return s;
}

#define SR_GRAIN 64
