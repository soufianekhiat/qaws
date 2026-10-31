#include "qaws_curve_batch.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_broadphase.h"
#include "internal/qaws_internal_flatten.h"
#include "internal/qaws_internal_batch.h"
#include "internal/qaws_internal_parallel.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* split by tools/splitc.py: shared between qaws_curve_batch.c and its moved part */
#define cb_dot QAWS_SFX(qaws_internal_curve_batch_cb_dot)
#define cb_flatten_all QAWS_SFX(qaws_internal_curve_batch_cb_flatten_all)
#define cb_check QAWS_SFX(qaws_internal_curve_batch_cb_check)

/*
 * Batched curve / curve intersection.
 *
 *   1. flatten: every curve, span by span, into chord segments; a piece is
 *      split while its midpoint or a quarter point lies farther than the
 *      flatness bound from the chord. A segment's box is its chord box
 *      inflated by twice the measured deviation.
 *   2. grid: one uniform grid over every segment of every curve, sized so a
 *      cell holds about one segment; filled by counting sort.
 *   3. pairs: inside each cell, segment pairs from curves that may meet and
 *      whose boxes overlap. A pair is handled only in the cell holding the
 *      low corner of the two boxes' overlap, so it is seen once.
 *   4. refine: chords closer than their two inflations seed Newton on the
 *      true curves (least squares in 3D); a bent piece running along the
 *      other within its deviation (a line through a shallow arc, two
 *      crossings in one pair) is seeded again from the segment ends.
 *   5. merge: hits of one curve pair that coincide are kept once; output is
 *      sorted by (curve_a, curve_b, parameter_a).
 */

#define CB_NEWTON_ITERS 24

#if QAWS_SCALAR_IS_FLOAT
#define CB_POS_REL   ((qaws_scalar)2e-5)
#define CB_PAR_REL   ((qaws_scalar)1e-4)
#define CB_STEP_REL  ((qaws_scalar)1e-7)   /* Newton steps below it: converged */
#else
#define CB_POS_REL   ((qaws_scalar)1e-10)
#define CB_PAR_REL   ((qaws_scalar)1e-7)
#define CB_STEP_REL  ((qaws_scalar)1e-15)
#endif

#define CB_PARALLEL_SIN ((qaws_scalar)0.25)  /* chords this parallel may overlap */
#define CB_TANGENT_SIN  ((qaws_scalar)1e-6)  /* hits this tangent are classified by their sides */
#if QAWS_SCALAR_IS_FLOAT
#define CB_MORE_SIN2    ((qaws_scalar)2.5e-3) /* squared sine: the extra seeds keep only crossings this clear (above a tangent contact's noise) */
#define CB_SIDE_MAX     ((qaws_scalar)1e-1)  /* farthest step off a contact to see its sides, and widest zone */
#else
#define CB_MORE_SIN2    ((qaws_scalar)1e-6)
#define CB_SIDE_MAX     ((qaws_scalar)1e-2)
#endif

typedef struct cb_hit
{
	unsigned int a, b, kind;
	qaws_scalar ta, tb;
	qaws_scalar ta1, tb1;       /* overlap ends; ta, tb for a point */
	qaws_scalar s2;             /* squared sine of the tangent angle at a point hit */
	qaws_scalar zlo, zhi;       /* contact zone on curve a around a point hit */
	qaws_scalar p[3];
} cb_hit;

typedef struct cb_ctx
{
	qaws_curve_batch_desc const* desc;
	unsigned int dim;
	qaws_flat_seg* segs;
	unsigned int nseg, capseg;
	unsigned int* seg_count;    /* per curve */
	cb_hit* hits;
	unsigned int nhit, caphit;
	qaws_scalar flat;
	qaws_scalar pos_tol;
	qaws_curve_batch_stats stats;
} cb_ctx;

qaws_scalar cb_dot(qaws_scalar const* a, qaws_scalar const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* 1. flatten every curve; segment counts per curve for the self-intersection rule */
qaws_status cb_flatten_all(cb_ctx* x, qaws_scalar* out_ext)
{
	qaws_curve_batch_desc const* d = x->desc;
	qaws_scalar ext;
	unsigned int i;
	qaws_status s = QAWS_STATUS_OK;
	/* scene extent from a cheap sampling of every curve, to scale the defaults */
	ext = qaws_internal_flatten_extent(d->curves, d->curve_count, x->dim, NULL, 0);
	x->flat = d->flatness > 0 ? d->flatness : ext / 1024;
	x->pos_tol = ext * CB_POS_REL;
	*out_ext = ext;
	x->seg_count = (unsigned int*)qaws_internal_zalloc(d->curve_count ? d->curve_count : 1, sizeof(unsigned int));
	if (!x->seg_count)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < d->curve_count && s == QAWS_STATUS_OK; i++)
	{
		unsigned int before = x->nseg;
		s = qaws_internal_flatten_curve(&x->segs, &x->nseg, &x->capseg, d->curves[i], x->dim, x->flat, i);
		x->seg_count[i] = x->nseg - before;
	}
	return s;
}

qaws_status cb_check(qaws_curve_batch_desc const* desc, unsigned int dim)
{
	unsigned int i;
	if (!desc || (desc->curve_count && !desc->curves))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->curve_count; i++)
	{
		if (!desc->curves[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (dim && desc->curves[i]->dimension != (dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D))
			return QAWS_STATUS_INVALID_DIMENSION;
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Prepared sets                                                      */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_set_create(qaws_curve_set** out_set, qaws_curve_batch_desc const* desc)
{
	qaws_curve_set* set;
	cb_ctx x;
	qaws_scalar ext;
	unsigned int n, dim;
	qaws_status s = cb_check(desc, 0);
	if (s != QAWS_STATUS_OK)
		return s;
	if (!out_set)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_set = NULL;
	n = desc->curve_count;
	dim = n ? (desc->curves[0]->dimension == QAWS_DIMENSION_2D ? 2 : 3) : 2;
	s = cb_check(desc, dim);
	if (s != QAWS_STATUS_OK)
		return s;
	set = (qaws_curve_set*)qaws_internal_zalloc(1, sizeof(qaws_curve_set));
	if (!set)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	set->desc = *desc;
	set->dim = dim;
	set->desc.curves = (qaws_curve const* const*)qaws_internal_malloc(qaws_internal_mul2((n ? n : 1), sizeof(qaws_curve*)));
	set->desc.families = desc->families ? (unsigned int const*)qaws_internal_malloc(qaws_internal_mul2((n ? n : 1), sizeof(unsigned int))) : NULL;
	if (!set->desc.curves || (desc->families && !set->desc.families))
	{
		qaws_curve_set_destroy(set);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy((void*)set->desc.curves, desc->curves, n * sizeof(qaws_curve*));
	if (desc->families)
		memcpy((void*)set->desc.families, desc->families, n * sizeof(unsigned int));
	memset(&x, 0, sizeof(x));
	x.desc = &set->desc;
	x.dim = dim;
	s = cb_flatten_all(&x, &ext);
	set->segs = x.segs;
	set->nseg = x.nseg;
	set->seg_count = x.seg_count;
	set->flat = x.flat;
	set->ext = ext;
	set->pos_tol = x.pos_tol;
	if (s != QAWS_STATUS_OK)
	{
		qaws_curve_set_destroy(set);
		return s;
	}
	*out_set = set;
	return QAWS_STATUS_OK;
}

void qaws_curve_set_destroy(qaws_curve_set* set)
{
	if (!set)
		return;
	qaws_internal_free((void*)set->desc.curves);
	qaws_internal_free((void*)set->desc.families);
	qaws_internal_free(set->segs);
	qaws_internal_free(set->seg_count);
	qaws_internal_free(set);
}

unsigned int qaws_curve_set_get_segment_count(qaws_curve_set const* set)
{
	return set ? set->nseg : 0;
}

/* ------------------------------------------------------------------ */
/*  Level crossings                                                    */
/* ------------------------------------------------------------------ */

#define LC_MAX_DEPTH 16
#define LC_ITERS 64

#define LC_GRAIN 64

/* ------------------------------------------------------------------ */
/*  Closest points                                                     */
/* ------------------------------------------------------------------ */

#define CP_ITERS 32

typedef struct cp_cand
{
	unsigned int seg;
	qaws_scalar lb;
} cp_cand;

typedef struct cp_ctx
{
	qaws_curve const* const* curves;
	qaws_flat_seg const* segs;
	unsigned int dim;
	unsigned int* stamp;        /* per segment: the last query that saw it */
	unsigned int query;
	qaws_scalar const* p;
	qaws_scalar best_ub;
	cp_cand* cand;
	unsigned int ncand, capcand;
	int failed;
} cp_ctx;

/* distance from p to the chord of s and the chord parameter of the nearest point */
static qaws_scalar cp_chord(qaws_flat_seg const* s, qaws_scalar const* p, qaws_scalar* u)
{
	qaws_scalar ab[3], ap[3], l2, w, d[3];
	int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = s->p1[k] - s->p0[k];
		ap[k] = p[k] - s->p0[k];
	}
	l2 = cb_dot(ab, ab);
	w = l2 > 0 ? cb_dot(ap, ab) / l2 : 0;
	if (w < 0) w = 0;
	if (w > 1) w = 1;
	for (k = 0; k < 3; k++)
		d[k] = ap[k] - w * ab[k];
	*u = w;
	return (qaws_scalar)sqrt(cb_dot(d, d));
}

static void cp_visit(void* user, unsigned int item)
{
	cp_ctx* x = (cp_ctx*)user;
	qaws_flat_seg const* s = &x->segs[item];
	qaws_scalar u, dc, lb, ub;
	if (x->stamp[item] == x->query)
		return;
	x->stamp[item] = x->query;
	dc = cp_chord(s, x->p, &u);
	lb = dc - s->r;
	ub = dc + s->r;
	if (lb > x->best_ub)
		return;
	if (ub < x->best_ub)
		x->best_ub = ub;
	if (x->ncand == x->capcand)
	{
		unsigned int cap = qaws_internal_grow(x->capcand, (size_t)x->ncand + 1, 64);
		cp_cand* g = cap ? (cp_cand*)qaws_internal_realloc(x->cand, qaws_internal_mul2(cap, sizeof(cp_cand))) : NULL;
		if (!g)
		{
			x->failed = 1;
			return;
		}
		x->cand = g;
		x->capcand = cap;
	}
	x->cand[x->ncand].seg = item;
	x->cand[x->ncand].lb = lb;
	x->ncand++;
}

static qaws_status cp_eval2(qaws_curve const* c, unsigned int dim, qaws_scalar t, qaws_scalar* P, qaws_scalar* D1, qaws_scalar* D2)
{
	unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;
	qaws_status s;
	if (dim == 2)
	{
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(&r, c, t, flags);
		P[0] = r.position.x; P[1] = r.position.y; P[2] = 0;
		D1[0] = r.d1.x; D1[1] = r.d1.y; D1[2] = 0;
		D2[0] = r.d2.x; D2[1] = r.d2.y; D2[2] = 0;
	}
	else
	{
		qaws_eval_result_3d r;
		s = qaws_curve_evaluate_3d(&r, c, t, flags);
		P[0] = r.position.x; P[1] = r.position.y; P[2] = r.position.z;
		D1[0] = r.d1.x; D1[1] = r.d1.y; D1[2] = r.d1.z;
		D2[0] = r.d2.x; D2[1] = r.d2.y; D2[2] = r.d2.z;
	}
	return s;
}

/* Newton on g(t) = (C(t) - p) . C'(t) from t, clamped to the domain. g also
   vanishes where C' = 0 (a cusp, a curve doubling back), which need not be
   nearest: a step that lands farther from p is halved back toward the
   nearest point so far, and that point is the answer (never worse than the
   start) */
static qaws_scalar cp_refine(qaws_curve const* c, unsigned int dim, qaws_scalar const* p, qaws_scalar* t, qaws_scalar* pos)
{
	qaws_scalar tmin = c->parameter_range.min_value, tmax = c->parameter_range.max_value, P[3], D1[3], D2[3], w[3];
	qaws_scalar best_t = *t, best_d = (qaws_scalar)HUGE_VAL;
	unsigned int it, k;
	for (it = 0; it < CP_ITERS; it++)
	{
		qaws_scalar g, gp, step, d;
		if (cp_eval2(c, dim, *t, P, D1, D2) != QAWS_STATUS_OK)
			break;
		for (k = 0; k < 3; k++)
			w[k] = P[k] - p[k];
		d = cb_dot(w, w);
		/* not clearly farther (beyond the rounding of the squared distance) */
		if (d <= best_d + (sizeof(qaws_scalar) == 4 ? 1e-5 : 1e-13) * (best_d + cb_dot(P, P) * (sizeof(qaws_scalar) == 4 ? 1e-6 : 1e-14)))
		{
			best_d = d;
			best_t = *t;
		}
		else
		{
			step = *t - best_t;
			*t = best_t + (qaws_scalar)0.5 * step;
			if (fabs(step) <= (tmax - tmin) * CB_STEP_REL)
				break;
			continue;
		}
		g = cb_dot(w, D1);
		gp = cb_dot(D1, D1) + cb_dot(w, D2);
		if (!(gp > 0))
			gp = cb_dot(D1, D1);
		if (!(gp > 0))
			break;
		step = g / gp;
		*t -= step;
		if (*t < tmin) *t = tmin;
		if (*t > tmax) *t = tmax;
		if (fabs(step) <= (tmax - tmin) * CB_STEP_REL)
			break;
	}
	*t = best_t;
	qaws_internal_curve_point(P, NULL, c, dim, *t, QAWS_EVAL_FLAG_POSITION);
	for (k = 0; k < 3; k++)
	{
		pos[k] = P[k];
		w[k] = P[k] - p[k];
	}
	return (qaws_scalar)sqrt(cb_dot(w, w));
}

#define CP_GRAIN 64

/* the shared state of a closest-point run; every chunk has its own scratch and statistics */
typedef struct cp_job
{
	qaws_curve const* const* curves;
	qaws_flat_seg const* segs;
	unsigned int nseg, dim;
	qaws_bp_grid* grid;
	qaws_scalar const* points;
	qaws_scalar max_distance;
	qaws_closest_point* out;
	qaws_curve_batch_stats* stats;   /* one per chunk */
} cp_job;

static qaws_status cp_chunk(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	cp_job const* J = (cp_job const*)ctx;
	qaws_curve_batch_stats* stats = &J->stats[chunk];
	unsigned int dim = J->dim, i;
	cp_ctx x;
	qaws_status s = QAWS_STATUS_OK;
	memset(&x, 0, sizeof(x));
	x.stamp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(J->nseg, sizeof(unsigned int)));
	if (!x.stamp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < J->nseg; i++)
		x.stamp[i] = ~0u;
	x.curves = J->curves;
	x.segs = J->segs;
	x.dim = dim;
	for (i = begin; i < end && s == QAWS_STATUS_OK; i++)
	{
		qaws_scalar p[3], best = (qaws_scalar)HUGE_VAL, pos[3];
		unsigned int r, c, best_curve = QAWS_CURVE_BATCH_NONE;
		qaws_scalar best_t = 0, best_pos[3] = { 0, 0, 0 };
		p[0] = J->points[i * dim];
		p[1] = J->points[i * dim + 1];
		p[2] = dim == 3 ? J->points[i * dim + 2] : 0;
		x.query = i;
		x.p = p;
		x.best_ub = J->max_distance > 0 ? J->max_distance : (qaws_scalar)HUGE_VAL;
		x.ncand = 0;
		/* rings outward until the next ring cannot beat the best bound */
		for (r = 0;; r++)
		{
			double pd[3], bound;
			pd[0] = (double)p[0]; pd[1] = (double)p[1]; pd[2] = (double)p[2];
			bound = qaws_internal_grid_ring(J->grid, pd, r, cp_visit, &x);
			if (x.failed)
			{
				s = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			if (bound == HUGE_VAL || bound > (double)x.best_ub)
				break;
		}
		/* refine the survivors */
		for (c = 0; c < x.ncand && s == QAWS_STATUS_OK; c++)
		{
			qaws_flat_seg const* sg = &J->segs[x.cand[c].seg];
			qaws_scalar u, t, d;
			if (x.cand[c].lb > x.best_ub)
				continue;
			stats->candidate_count++;
			stats->newton_count++;
			cp_chord(sg, p, &u);
			t = sg->t0 + u * (sg->t1 - sg->t0);
			d = cp_refine(J->curves[sg->owner], dim, p, &t, pos);
			if (d < best || (d == best && sg->owner < best_curve))
			{
				best = d;
				best_curve = sg->owner;
				best_t = t;
				best_pos[0] = pos[0]; best_pos[1] = pos[1]; best_pos[2] = pos[2];
			}
		}
		if (best_curve != QAWS_CURVE_BATCH_NONE && (J->max_distance <= 0 || best <= J->max_distance))
		{
			qaws_closest_point* o = &J->out[i];
			o->curve = best_curve;
			o->parameter = best_t;
			o->distance = best;
			o->position.x = best_pos[0];
			o->position.y = best_pos[1];
			o->position.z = best_pos[2];
			stats->hit_count++;
		}
	}
	qaws_internal_free(x.stamp);
	qaws_internal_free(x.cand);
	return s;
}

static qaws_status cp_run(qaws_curve const* const* curves, unsigned int dim, qaws_flat_seg const* segs, unsigned int nseg, qaws_scalar const* points,
	unsigned int point_count, qaws_scalar max_distance, qaws_batch_executor const* executor, qaws_closest_point* out, qaws_curve_batch_stats* stats)
{
	cp_job J;
	qaws_bp_box* boxes;
	unsigned int i, k, nchunk;
	qaws_status s;
	memset(&J, 0, sizeof(J));
	stats->segment_count = nseg;
	for (i = 0; i < point_count; i++)
	{
		out[i].curve = QAWS_CURVE_BATCH_NONE;
		out[i].parameter = out[i].distance = 0;
		out[i].position.x = out[i].position.y = out[i].position.z = 0;
	}
	if (!nseg || !point_count)
		return QAWS_STATUS_OK;
	boxes = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(nseg, sizeof(qaws_bp_box)));
	nchunk = qaws_internal_chunk_count(point_count, CP_GRAIN);
	J.stats = (qaws_curve_batch_stats*)qaws_internal_zalloc(nchunk, sizeof(qaws_curve_batch_stats));
	if (!boxes || !J.stats)
	{
		qaws_internal_free(boxes);
		qaws_internal_free(J.stats);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nseg; i++)
		for (k = 0; k < 3; k++)
		{
			boxes[i].lo[k] = (double)segs[i].lo[k];
			boxes[i].hi[k] = (double)segs[i].hi[k];
		}
	s = qaws_internal_grid_create(&J.grid, boxes, nseg, dim);
	qaws_internal_free(boxes);
	J.curves = curves;
	J.segs = segs;
	J.nseg = nseg;
	J.dim = dim;
	J.points = points;
	J.max_distance = max_distance;
	J.out = out;
	if (s == QAWS_STATUS_OK)
		s = qaws_internal_parallel(point_count, CP_GRAIN, cp_chunk, &J, executor);
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

qaws_status qaws_curve_batch_find_closest(
	qaws_closest_point* out_points,
	qaws_curve_batch_stats* out_stats,
	qaws_closest_desc const* desc)
{
	qaws_curve_batch_stats st;
	qaws_flat_seg* segs = NULL;
	unsigned int nseg = 0, capseg = 0, i, dim;
	qaws_scalar ext, flat;
	qaws_status s = QAWS_STATUS_OK;
	memset(&st, 0, sizeof(st));
	if (out_stats)
		*out_stats = st;
	if (!desc || (desc->curve_count && !desc->curves) || (desc->point_count && (!desc->points || !out_points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->curve_count; i++)
	{
		if (!desc->curves[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (desc->curves[i]->dimension != desc->curves[0]->dimension)
			return QAWS_STATUS_INVALID_DIMENSION;
	}
	dim = desc->curve_count && desc->curves[0]->dimension == QAWS_DIMENSION_3D ? 3 : 2;
	ext = qaws_internal_flatten_extent(desc->curves, desc->curve_count, dim, NULL, 0);
	flat = desc->flatness > 0 ? desc->flatness : ext / 1024;
	for (i = 0; i < desc->curve_count && s == QAWS_STATUS_OK; i++)
		s = qaws_internal_flatten_curve(&segs, &nseg, &capseg, desc->curves[i], dim, flat, i);
	if (s == QAWS_STATUS_OK)
		s = cp_run(desc->curves, dim, segs, nseg, desc->points, desc->point_count, desc->max_distance, desc->executor, out_points, &st);
	if (out_stats)
		*out_stats = st;
	qaws_internal_free(segs);
	return s;
}

qaws_status qaws_curve_set_find_closest(
	qaws_closest_point* out_points,
	qaws_curve_batch_stats* out_stats,
	qaws_curve_set const* set,
	qaws_scalar const* points,
	unsigned int point_count,
	qaws_scalar max_distance)
{
	qaws_curve_batch_stats st;
	qaws_status s;
	memset(&st, 0, sizeof(st));
	if (out_stats)
		*out_stats = st;
	if (!set || (point_count && (!points || !out_points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	s = cp_run(set->desc.curves, set->dim, set->segs, set->nseg, points, point_count, max_distance, set->desc.executor, out_points, &st);
	if (out_stats)
		*out_stats = st;
	return s;
}
