#include "qaws_curve_batch.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_broadphase.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

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
 *      true curves (least squares in 3D).
 *   5. merge: hits of one curve pair that coincide are kept once; output is
 *      sorted by (curve_a, curve_b, parameter_a).
 */

#define CB_MAX_DEPTH   22
#define CB_NEWTON_ITERS 24

#if QAWS_SCALAR_IS_FLOAT
#define CB_POS_REL   ((qaws_scalar)2e-5)
#define CB_PAR_REL   ((qaws_scalar)1e-4)
#else
#define CB_POS_REL   ((qaws_scalar)1e-10)
#define CB_PAR_REL   ((qaws_scalar)1e-7)
#endif

typedef struct cb_seg
{
	qaws_scalar t0, t1;
	qaws_scalar p0[3], p1[3];
	qaws_scalar lo[3], hi[3];   /* inflated box */
	qaws_scalar r;              /* inflation */
	unsigned int curve;
	unsigned int index;         /* position along its curve */
} cb_seg;

typedef struct cb_hit
{
	unsigned int a, b;
	qaws_scalar ta, tb;
	qaws_scalar p[3];
} cb_hit;

typedef struct cb_ctx
{
	qaws_curve_batch_desc const* desc;
	unsigned int dim;
	cb_seg* segs;
	unsigned int nseg, capseg;
	unsigned int* seg_count;    /* per curve */
	cb_hit* hits;
	unsigned int nhit, caphit;
	qaws_scalar flat;
	qaws_scalar pos_tol;
	qaws_curve_batch_stats stats;
} cb_ctx;

static qaws_status cb_eval(qaws_curve const* c, unsigned int dim, qaws_scalar t, unsigned int flags, qaws_scalar* p, qaws_scalar* d)
{
	qaws_status s;
	if (dim == 2)
	{
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(c, t, flags, &r);
		p[0] = r.position.x; p[1] = r.position.y; p[2] = 0;
		if (d) { d[0] = r.d1.x; d[1] = r.d1.y; d[2] = 0; }
	}
	else
	{
		qaws_eval_result_3d r;
		s = qaws_curve_evaluate_3d(c, t, flags, &r);
		p[0] = r.position.x; p[1] = r.position.y; p[2] = r.position.z;
		if (d) { d[0] = r.d1.x; d[1] = r.d1.y; d[2] = r.d1.z; }
	}
	return s;
}

static qaws_scalar cb_dot(qaws_scalar const* a, qaws_scalar const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* distance from m to the segment [a, b] */
static qaws_scalar cb_point_segment(qaws_scalar const* m, qaws_scalar const* a, qaws_scalar const* b)
{
	qaws_scalar ab[3], am[3], l2, u, d[3];
	int k;
	for (k = 0; k < 3; k++) { ab[k] = b[k] - a[k]; am[k] = m[k] - a[k]; }
	l2 = cb_dot(ab, ab);
	u = l2 > 0 ? cb_dot(am, ab) / l2 : 0;
	if (u < 0) u = 0;
	if (u > 1) u = 1;
	for (k = 0; k < 3; k++) d[k] = am[k] - u * ab[k];
	return (qaws_scalar)sqrt(cb_dot(d, d));
}

static qaws_status cb_push_seg(cb_ctx* x, unsigned int curve, qaws_scalar t0, qaws_scalar const* p0, qaws_scalar t1, qaws_scalar const* p1, qaws_scalar dev)
{
	cb_seg* s;
	int k;
	if (x->nseg == x->capseg)
	{
		unsigned int cap = x->capseg ? x->capseg * 2 : 1024;
		cb_seg* g = (cb_seg*)realloc(x->segs, cap * sizeof(cb_seg));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		x->segs = g;
		x->capseg = cap;
	}
	s = &x->segs[x->nseg];
	s->t0 = t0; s->t1 = t1;
	s->r = 2 * dev;
	for (k = 0; k < 3; k++)
	{
		s->p0[k] = p0[k];
		s->p1[k] = p1[k];
		s->lo[k] = (p0[k] < p1[k] ? p0[k] : p1[k]) - s->r;
		s->hi[k] = (p0[k] > p1[k] ? p0[k] : p1[k]) + s->r;
	}
	s->curve = curve;
	s->index = x->seg_count[curve]++;
	x->nseg++;
	return QAWS_STATUS_OK;
}

/* piece [t0, t1] with its known midpoint: the midpoint and both quarter
   points must lie within the flatness bound of the chord */
static qaws_status cb_flatten(cb_ctx* x, unsigned int curve, qaws_scalar t0, qaws_scalar const* p0, qaws_scalar tm, qaws_scalar const* pm,
	qaws_scalar t1, qaws_scalar const* p1, unsigned int depth)
{
	qaws_curve const* c = x->desc->curves[curve];
	qaws_scalar ta = (t0 + tm) / 2, tb = (tm + t1) / 2, pa[3], pb[3], dev, d;
	qaws_status s = cb_eval(c, x->dim, ta, QAWS_EVAL_FLAG_POSITION, pa, NULL);
	if (s == QAWS_STATUS_OK)
		s = cb_eval(c, x->dim, tb, QAWS_EVAL_FLAG_POSITION, pb, NULL);
	if (s != QAWS_STATUS_OK)
		return s;
	dev = cb_point_segment(pm, p0, p1);
	d = cb_point_segment(pa, p0, p1);
	if (d > dev) dev = d;
	d = cb_point_segment(pb, p0, p1);
	if (d > dev) dev = d;
	if (depth >= CB_MAX_DEPTH || dev <= x->flat)
		return cb_push_seg(x, curve, t0, p0, t1, p1, dev);
	s = cb_flatten(x, curve, t0, p0, ta, pa, tm, pm, depth + 1);
	if (s != QAWS_STATUS_OK)
		return s;
	return cb_flatten(x, curve, tm, pm, tb, pb, t1, p1, depth + 1);
}

static qaws_status cb_flatten_curve(cb_ctx* x, unsigned int curve)
{
	qaws_curve const* c = x->desc->curves[curve];
	qaws_scalar p0[3], p1[3], pm[3];
	unsigned int k, have = 0;
	for (k = 0; k < c->span_count; k++)
	{
		qaws_scalar t0 = c->span_boundaries[k], t1 = c->span_boundaries[k + 1], tm = (t0 + t1) / 2;
		qaws_status s = QAWS_STATUS_OK;
		if (!(t1 > t0))
			continue;
		/* a span starts where the previous one ended */
		if (!have)
			s = cb_eval(c, x->dim, t0, QAWS_EVAL_FLAG_POSITION, p0, NULL);
		if (s == QAWS_STATUS_OK)
			s = cb_eval(c, x->dim, t1, QAWS_EVAL_FLAG_POSITION, p1, NULL);
		if (s == QAWS_STATUS_OK)
			s = cb_eval(c, x->dim, tm, QAWS_EVAL_FLAG_POSITION, pm, NULL);
		if (s == QAWS_STATUS_OK)
			s = cb_flatten(x, curve, t0, p0, tm, pm, t1, p1, 0);
		if (s != QAWS_STATUS_OK)
			return s;
		memcpy(p0, p1, sizeof(p0));
		have = 1;
	}
	return QAWS_STATUS_OK;
}

/* closest points of the chords P(u) = a0 + u (a1 - a0), Q(v) = b0 + v (b1 - b0) */
static qaws_scalar cb_chords(cb_seg const* A, cb_seg const* B, qaws_scalar* u, qaws_scalar* v)
{
	qaws_scalar d1[3], d2[3], r[3], a, e, f, c, b, den, s, t, w[3];
	int k;
	for (k = 0; k < 3; k++)
	{
		d1[k] = A->p1[k] - A->p0[k];
		d2[k] = B->p1[k] - B->p0[k];
		r[k] = A->p0[k] - B->p0[k];
	}
	a = cb_dot(d1, d1); e = cb_dot(d2, d2); f = cb_dot(d2, r);
	if (a <= 0 && e <= 0)
		s = t = 0;
	else if (a <= 0)
	{
		s = 0;
		t = f / e;
		t = t < 0 ? 0 : (t > 1 ? 1 : t);
	}
	else
	{
		c = cb_dot(d1, r);
		if (e <= 0)
		{
			t = 0;
			s = -c / a;
			s = s < 0 ? 0 : (s > 1 ? 1 : s);
		}
		else
		{
			b = cb_dot(d1, d2);
			den = a * e - b * b;
			s = den > 0 ? (b * f - c * e) / den : 0;
			s = s < 0 ? 0 : (s > 1 ? 1 : s);
			t = (b * s + f) / e;
			if (t < 0)
			{
				t = 0;
				s = -c / a;
				s = s < 0 ? 0 : (s > 1 ? 1 : s);
			}
			else if (t > 1)
			{
				t = 1;
				s = (b - c) / a;
				s = s < 0 ? 0 : (s > 1 ? 1 : s);
			}
		}
	}
	for (k = 0; k < 3; k++)
		w[k] = r[k] + s * d1[k] - t * d2[k];
	*u = s;
	*v = t;
	return (qaws_scalar)sqrt(cb_dot(w, w));
}

/* Newton on C_a(ta) = C_b(tb); in 3D the 3x2 system by its normal equations */
static int cb_newton(cb_ctx* x, qaws_curve const* ca, qaws_curve const* cb, qaws_scalar* ta, qaws_scalar* tb, qaws_scalar* pos)
{
	qaws_scalar amin = ca->parameter_range.min_value, amax = ca->parameter_range.max_value;
	qaws_scalar bmin = cb->parameter_range.min_value, bmax = cb->parameter_range.max_value;
	unsigned int it, polish = 0;
	for (it = 0; it < CB_NEWTON_ITERS; it++)
	{
		qaws_scalar pa[3], da[3], pb[3], db[3], f[3], m00, m01, m11, g0, g1, det, dta, dtb;
		int k;
		if (cb_eval(ca, x->dim, *ta, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, pa, da) != QAWS_STATUS_OK
			|| cb_eval(cb, x->dim, *tb, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, pb, db) != QAWS_STATUS_OK)
			return 0;
		for (k = 0; k < 3; k++)
			f[k] = pa[k] - pb[k];
		/* converged: two more steps take a shallow crossing to rounding level */
		if (sqrt(cb_dot(f, f)) <= x->pos_tol && polish++ == 2)
		{
			for (k = 0; k < 3; k++)
				pos[k] = (pa[k] + pb[k]) / 2;
			return 1;
		}
		/* J = [da, -db]; solve (J^T J) d = J^T f */
		m00 = cb_dot(da, da);
		m01 = -cb_dot(da, db);
		m11 = cb_dot(db, db);
		g0 = cb_dot(da, f);
		g1 = -cb_dot(db, f);
		det = m00 * m11 - m01 * m01;
		if (!(det > (m00 * m11) * (qaws_scalar)1e-12))
		{
			if (!polish)
				return 0;
			for (k = 0; k < 3; k++)
				pos[k] = (pa[k] + pb[k]) / 2;
			return 1;
		}
		dta = (m11 * g0 - m01 * g1) / det;
		dtb = (m00 * g1 - m01 * g0) / det;
		*ta -= dta;
		*tb -= dtb;
		if (*ta < amin) *ta = amin;
		if (*ta > amax) *ta = amax;
		if (*tb < bmin) *tb = bmin;
		if (*tb > bmax) *tb = bmax;
	}
	return 0;
}

static qaws_status cb_push_hit(cb_ctx* x, unsigned int a, unsigned int b, qaws_scalar ta, qaws_scalar tb, qaws_scalar const* p)
{
	cb_hit* h;
	if (x->nhit == x->caphit)
	{
		unsigned int cap = x->caphit ? x->caphit * 2 : 256;
		cb_hit* g = (cb_hit*)realloc(x->hits, cap * sizeof(cb_hit));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		x->hits = g;
		x->caphit = cap;
	}
	h = &x->hits[x->nhit++];
	if (a > b || (a == b && ta > tb))
	{
		unsigned int ti = a; qaws_scalar tt = ta;
		a = b; b = ti;
		ta = tb; tb = tt;
	}
	h->a = a; h->b = b; h->ta = ta; h->tb = tb;
	h->p[0] = p[0]; h->p[1] = p[1]; h->p[2] = p[2];
	return QAWS_STATUS_OK;
}

static qaws_scalar cb_par_tol(qaws_curve const* c)
{
	return (c->parameter_range.max_value - c->parameter_range.min_value) * CB_PAR_REL;
}

static qaws_status cb_pair(cb_ctx* x, cb_seg const* A, cb_seg const* B)
{
	qaws_curve const* ca = x->desc->curves[A->curve];
	qaws_curve const* cb = x->desc->curves[B->curve];
	qaws_scalar u, v, ta, tb, pos[3];
	if (cb_chords(A, B, &u, &v) > A->r + B->r + x->pos_tol)
		return QAWS_STATUS_OK;
	x->stats.newton_count++;
	ta = A->t0 + u * (A->t1 - A->t0);
	tb = B->t0 + v * (B->t1 - B->t0);
	if (!cb_newton(x, ca, cb, &ta, &tb, pos))
		return QAWS_STATUS_OK;
	/* a curve meets itself trivially at ta == tb */
	if (A->curve == B->curve && fabs(ta - tb) <= 64 * cb_par_tol(ca))
		return QAWS_STATUS_OK;
	return cb_push_hit(x, A->curve, B->curve, ta, tb, pos);
}

/* may segments i and j meet as a pair to report? */
static int cb_allowed(cb_ctx const* x, cb_seg const* A, cb_seg const* B)
{
	qaws_curve_batch_desc const* d = x->desc;
	if (A->curve == B->curve)
	{
		unsigned int n, lo, hi;
		if (!(d->flags & QAWS_CURVE_BATCH_SELF))
			return 0;
		lo = A->index < B->index ? A->index : B->index;
		hi = A->index < B->index ? B->index : A->index;
		if (hi - lo <= 1)
			return 0;
		n = x->seg_count[A->curve];
		if (lo == 0 && hi == n - 1 && qaws_curve_is_closed(d->curves[A->curve]))
			return 0;
		return 1;
	}
	return !d->families || d->families[A->curve] != d->families[B->curve];
}

static int cb_accept(void* user, unsigned int i, unsigned int j)
{
	cb_ctx const* x = (cb_ctx const*)user;
	return cb_allowed(x, &x->segs[i], &x->segs[j]);
}

static qaws_status cb_visit(void* user, unsigned int i, unsigned int j)
{
	cb_ctx* x = (cb_ctx*)user;
	return cb_pair(x, &x->segs[i], &x->segs[j]);
}

static int cb_cmp_hit(void const* p, void const* q)
{
	cb_hit const* a = (cb_hit const*)p;
	cb_hit const* b = (cb_hit const*)q;
	if (a->a != b->a) return a->a < b->a ? -1 : 1;
	if (a->b != b->b) return a->b < b->b ? -1 : 1;
	if (a->ta != b->ta) return a->ta < b->ta ? -1 : 1;
	if (a->tb != b->tb) return a->tb < b->tb ? -1 : 1;
	return 0;
}

/* keep one hit per point of each curve pair */
static void cb_merge(cb_ctx* x)
{
	unsigned int i, n = 0, group = 0;
	qsort(x->hits, x->nhit, sizeof(cb_hit), cb_cmp_hit);
	for (i = 0; i < x->nhit; i++)
	{
		cb_hit const* h = &x->hits[i];
		qaws_scalar pa = cb_par_tol(x->desc->curves[h->a]), pb = cb_par_tol(x->desc->curves[h->b]);
		unsigned int j;
		int dup = 0;
		if (n == 0 || x->hits[n - 1].a != h->a || x->hits[n - 1].b != h->b)
			group = n;
		for (j = group; j < n && !dup; j++)
		{
			cb_hit const* k = &x->hits[j];
			qaws_scalar d[3];
			int c;
			for (c = 0; c < 3; c++)
				d[c] = k->p[c] - h->p[c];
			dup = (fabs(k->ta - h->ta) <= pa && fabs(k->tb - h->tb) <= pb) || sqrt(cb_dot(d, d)) <= 16 * x->pos_tol;
		}
		if (!dup)
			x->hits[n++] = *h;
	}
	x->nhit = n;
}

static qaws_status cb_run(cb_ctx* x)
{
	qaws_curve_batch_desc const* d = x->desc;
	qaws_scalar lo[3], hi[3], ext = 0;
	unsigned int i, k;
	qaws_status s = QAWS_STATUS_OK;

	/* scene extent from a cheap sampling of every curve, to scale the defaults */
	for (k = 0; k < 3; k++) { lo[k] = (qaws_scalar)HUGE_VAL; hi[k] = -(qaws_scalar)HUGE_VAL; }
	for (i = 0; i < d->curve_count; i++)
	{
		qaws_curve const* c = d->curves[i];
		unsigned int j;
		for (j = 0; j <= 16; j++)
		{
			qaws_scalar t = c->parameter_range.min_value + (c->parameter_range.max_value - c->parameter_range.min_value) * (qaws_scalar)j / 16, p[3];
			s = cb_eval(c, x->dim, t, QAWS_EVAL_FLAG_POSITION, p, NULL);
			if (s != QAWS_STATUS_OK)
				return s;
			for (k = 0; k < 3; k++)
			{
				if (p[k] < lo[k]) lo[k] = p[k];
				if (p[k] > hi[k]) hi[k] = p[k];
			}
		}
	}
	for (k = 0; k < 3; k++)
		if (hi[k] - lo[k] > ext)
			ext = hi[k] - lo[k];
	if (!(ext > 0))
		ext = 1;
	x->flat = d->flatness > 0 ? d->flatness : ext / 1024;
	x->pos_tol = ext * CB_POS_REL;

	/* 1. flatten */
	x->seg_count = (unsigned int*)calloc(d->curve_count, sizeof(unsigned int));
	if (!x->seg_count)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < d->curve_count && s == QAWS_STATUS_OK; i++)
		s = cb_flatten_curve(x, i);
	if (s != QAWS_STATUS_OK || x->nseg == 0)
		return s;
	x->stats.segment_count = x->nseg;

	/* 2-3. one grid over every segment; each overlapping allowed pair once */
	{
		qaws_bp_box* boxes = (qaws_bp_box*)malloc(x->nseg * sizeof(qaws_bp_box));
		qaws_bp_stats bs;
		if (!boxes)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < x->nseg; i++)
			for (k = 0; k < 3; k++)
			{
				boxes[i].lo[k] = (double)x->segs[i].lo[k];
				boxes[i].hi[k] = (double)x->segs[i].hi[k];
			}
		s = qaws_internal_broadphase(boxes, x->nseg, x->dim, cb_accept, cb_visit, x, &bs);
		free(boxes);
		x->stats.cell_count = bs.cell_count;
		x->stats.candidate_count = bs.candidate_count;
	}

	/* 5. merge */
	if (s == QAWS_STATUS_OK)
		cb_merge(x);
	return s;
}

static qaws_status cb_find(qaws_curve_batch_desc const* desc, unsigned int dim, void* out, unsigned int capacity,
	unsigned int* out_count, qaws_curve_batch_stats* out_stats)
{
	cb_ctx x;
	qaws_status s;
	unsigned int i;
	if (!desc || !out_count || (!out && capacity) || (desc->curve_count && !desc->curves))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->curve_count; i++)
	{
		if (!desc->curves[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (desc->curves[i]->dimension != (dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D))
			return QAWS_STATUS_INVALID_DIMENSION;
	}
	*out_count = 0;
	memset(&x, 0, sizeof(x));
	x.desc = desc;
	x.dim = dim;
	s = cb_run(&x);
	if (s == QAWS_STATUS_OK)
	{
		for (i = 0; i < x.nhit && i < capacity; i++)
		{
			cb_hit const* h = &x.hits[i];
			if (dim == 2)
			{
				qaws_curve_batch_hit_2d* o = (qaws_curve_batch_hit_2d*)out + i;
				o->curve_a = h->a; o->curve_b = h->b;
				o->parameter_a = h->ta; o->parameter_b = h->tb;
				o->position.x = h->p[0]; o->position.y = h->p[1];
			}
			else
			{
				qaws_curve_batch_hit_3d* o = (qaws_curve_batch_hit_3d*)out + i;
				o->curve_a = h->a; o->curve_b = h->b;
				o->parameter_a = h->ta; o->parameter_b = h->tb;
				o->position.x = h->p[0]; o->position.y = h->p[1]; o->position.z = h->p[2];
			}
		}
		*out_count = x.nhit;
		x.stats.hit_count = x.nhit;
	}
	if (out_stats)
		*out_stats = x.stats;
	free(x.segs);
	free(x.hits);
	free(x.seg_count);
	return s;
}

qaws_status qaws_curve_batch_find_intersections_2d(
	qaws_curve_batch_desc const* desc,
	qaws_curve_batch_hit_2d* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_curve_batch_stats* out_stats)
{
	return cb_find(desc, 2, out_hits, hit_capacity, out_count, out_stats);
}

qaws_status qaws_curve_batch_find_intersections_3d(
	qaws_curve_batch_desc const* desc,
	qaws_curve_batch_hit_3d* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_curve_batch_stats* out_stats)
{
	return cb_find(desc, 3, out_hits, hit_capacity, out_count, out_stats);
}
