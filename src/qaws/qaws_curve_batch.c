#include "qaws_curve_batch.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_broadphase.h"
#include "internal/qaws_internal_flatten.h"
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

#define CB_NEWTON_ITERS 24

#if QAWS_SCALAR_IS_FLOAT
#define CB_POS_REL   ((qaws_scalar)2e-5)
#define CB_PAR_REL   ((qaws_scalar)1e-4)
#else
#define CB_POS_REL   ((qaws_scalar)1e-10)
#define CB_PAR_REL   ((qaws_scalar)1e-7)
#endif

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
	qaws_flat_seg* segs;
	unsigned int nseg, capseg;
	unsigned int* seg_count;    /* per curve */
	cb_hit* hits;
	unsigned int nhit, caphit;
	qaws_scalar flat;
	qaws_scalar pos_tol;
	qaws_curve_batch_stats stats;
} cb_ctx;

static qaws_scalar cb_dot(qaws_scalar const* a, qaws_scalar const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* closest points of the chords P(u) = a0 + u (a1 - a0), Q(v) = b0 + v (b1 - b0) */
static qaws_scalar cb_chords(qaws_flat_seg const* A, qaws_flat_seg const* B, qaws_scalar* u, qaws_scalar* v)
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
		if (qaws_internal_curve_point(ca, x->dim, *ta, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, pa, da) != QAWS_STATUS_OK
			|| qaws_internal_curve_point(cb, x->dim, *tb, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, pb, db) != QAWS_STATUS_OK)
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

static qaws_status cb_pair(cb_ctx* x, qaws_flat_seg const* A, qaws_flat_seg const* B)
{
	qaws_curve const* ca = x->desc->curves[A->owner];
	qaws_curve const* cb = x->desc->curves[B->owner];
	qaws_scalar u, v, ta, tb, pos[3];
	if (cb_chords(A, B, &u, &v) > A->r + B->r + x->pos_tol)
		return QAWS_STATUS_OK;
	x->stats.newton_count++;
	ta = A->t0 + u * (A->t1 - A->t0);
	tb = B->t0 + v * (B->t1 - B->t0);
	if (!cb_newton(x, ca, cb, &ta, &tb, pos))
		return QAWS_STATUS_OK;
	/* a curve meets itself trivially at ta == tb */
	if (A->owner == B->owner && fabs(ta - tb) <= 64 * cb_par_tol(ca))
		return QAWS_STATUS_OK;
	return cb_push_hit(x, A->owner, B->owner, ta, tb, pos);
}

/* may segments i and j meet as a pair to report? */
static int cb_allowed(cb_ctx const* x, qaws_flat_seg const* A, qaws_flat_seg const* B)
{
	qaws_curve_batch_desc const* d = x->desc;
	if (A->owner == B->owner)
	{
		unsigned int n, lo, hi;
		if (!(d->flags & QAWS_CURVE_BATCH_SELF))
			return 0;
		lo = A->index < B->index ? A->index : B->index;
		hi = A->index < B->index ? B->index : A->index;
		if (hi - lo <= 1)
			return 0;
		n = x->seg_count[A->owner];
		if (lo == 0 && hi == n - 1 && qaws_curve_is_closed(d->curves[A->owner]))
			return 0;
		return 1;
	}
	return !d->families || d->families[A->owner] != d->families[B->owner];
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

/* 1. flatten every curve; segment counts per curve for the self-intersection rule */
static qaws_status cb_flatten_all(cb_ctx* x, qaws_scalar* out_ext)
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
	x->seg_count = (unsigned int*)calloc(d->curve_count ? d->curve_count : 1, sizeof(unsigned int));
	if (!x->seg_count)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < d->curve_count && s == QAWS_STATUS_OK; i++)
	{
		unsigned int before = x->nseg;
		s = qaws_internal_flatten_curve(d->curves[i], x->dim, x->flat, i, &x->segs, &x->nseg, &x->capseg);
		x->seg_count[i] = x->nseg - before;
	}
	return s;
}

/* 2-5. one grid over every segment; each overlapping allowed pair refined; merge */
static qaws_status cb_solve(cb_ctx* x)
{
	unsigned int i, k;
	qaws_status s;
	qaws_bp_box* boxes;
	qaws_bp_stats bs;
	x->stats.segment_count = x->nseg;
	if (x->nseg == 0)
		return QAWS_STATUS_OK;
	boxes = (qaws_bp_box*)malloc(x->nseg * sizeof(qaws_bp_box));
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
	if (s == QAWS_STATUS_OK)
		cb_merge(x);
	return s;
}

/* hits out; curve_b is moved back by b_offset (pairs across two sets) */
static void cb_emit(cb_ctx* x, unsigned int dim, void* out, unsigned int capacity, unsigned int* out_count, unsigned int b_offset)
{
	unsigned int i;
	for (i = 0; i < x->nhit && i < capacity; i++)
	{
		cb_hit const* h = &x->hits[i];
		if (dim == 2)
		{
			qaws_curve_batch_hit_2d* o = (qaws_curve_batch_hit_2d*)out + i;
			o->curve_a = h->a; o->curve_b = h->b - b_offset;
			o->parameter_a = h->ta; o->parameter_b = h->tb;
			o->position.x = h->p[0]; o->position.y = h->p[1];
		}
		else
		{
			qaws_curve_batch_hit_3d* o = (qaws_curve_batch_hit_3d*)out + i;
			o->curve_a = h->a; o->curve_b = h->b - b_offset;
			o->parameter_a = h->ta; o->parameter_b = h->tb;
			o->position.x = h->p[0]; o->position.y = h->p[1]; o->position.z = h->p[2];
		}
	}
	*out_count = x->nhit;
	x->stats.hit_count = x->nhit;
}

static qaws_status cb_check(qaws_curve_batch_desc const* desc, unsigned int dim)
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

static qaws_status cb_find(qaws_curve_batch_desc const* desc, unsigned int dim, void* out, unsigned int capacity,
	unsigned int* out_count, qaws_curve_batch_stats* out_stats)
{
	cb_ctx x;
	qaws_scalar ext;
	qaws_status s = cb_check(desc, dim);
	if (s != QAWS_STATUS_OK)
		return s;
	if (!out_count || (!out && capacity))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	memset(&x, 0, sizeof(x));
	x.desc = desc;
	x.dim = dim;
	s = cb_flatten_all(&x, &ext);
	if (s == QAWS_STATUS_OK)
		s = cb_solve(&x);
	if (s == QAWS_STATUS_OK)
		cb_emit(&x, dim, out, capacity, out_count, 0);
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

/* ------------------------------------------------------------------ */
/*  Prepared sets                                                      */
/* ------------------------------------------------------------------ */

struct qaws_curve_set
{
	qaws_curve_batch_desc desc;     /* curves and families owned by the set */
	unsigned int dim;
	qaws_flat_seg* segs;
	unsigned int nseg;
	unsigned int* seg_count;
	qaws_scalar flat, ext, pos_tol;
};

qaws_status qaws_curve_set_create(qaws_curve_batch_desc const* desc, qaws_curve_set** out_set)
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
	set = (qaws_curve_set*)calloc(1, sizeof(qaws_curve_set));
	if (!set)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	set->desc = *desc;
	set->dim = dim;
	set->desc.curves = (qaws_curve const* const*)malloc((n ? n : 1) * sizeof(qaws_curve*));
	set->desc.families = desc->families ? (unsigned int const*)malloc((n ? n : 1) * sizeof(unsigned int)) : NULL;
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
	free((void*)set->desc.curves);
	free((void*)set->desc.families);
	free(set->segs);
	free(set->seg_count);
	free(set);
}

unsigned int qaws_curve_set_get_segment_count(qaws_curve_set const* set)
{
	return set ? set->nseg : 0;
}

/* one set with itself (b NULL), or every pair across a and b */
static qaws_status cs_find(qaws_curve_set const* a, qaws_curve_set const* b, unsigned int dim, void* out, unsigned int capacity,
	unsigned int* out_count, qaws_curve_batch_stats* out_stats)
{
	cb_ctx x;
	qaws_curve_batch_desc d;
	qaws_curve const** curves = NULL;
	unsigned int* fam = NULL;
	qaws_flat_seg* segs = NULL;
	unsigned int* counts = NULL;
	unsigned int na, nb = 0, i;
	qaws_status s = QAWS_STATUS_OK;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (!a || !out_count || (!out && capacity))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if ((a->desc.curve_count && a->dim != dim) || (b && b->desc.curve_count && b->dim != dim))
		return QAWS_STATUS_INVALID_DIMENSION;
	memset(&x, 0, sizeof(x));
	x.dim = dim;
	na = a->desc.curve_count;
	if (!b)
	{
		/* the prepared segments as they are */
		x.desc = &a->desc;
		x.segs = a->segs;
		x.nseg = a->nseg;
		x.seg_count = a->seg_count;
		x.flat = a->flat;
		x.pos_tol = a->pos_tol;
	}
	else
	{
		/* both sets side by side: b's curves numbered after a's, in another family */
		nb = b->desc.curve_count;
		curves = (qaws_curve const**)malloc((na + nb ? na + nb : 1) * sizeof(qaws_curve*));
		fam = (unsigned int*)malloc((na + nb ? na + nb : 1) * sizeof(unsigned int));
		segs = (qaws_flat_seg*)malloc((a->nseg + b->nseg ? a->nseg + b->nseg : 1) * sizeof(qaws_flat_seg));
		counts = (unsigned int*)malloc((na + nb ? na + nb : 1) * sizeof(unsigned int));
		if (!curves || !fam || !segs || !counts)
		{
			s = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		memcpy((void*)curves, a->desc.curves, na * sizeof(qaws_curve*));
		memcpy((void*)(curves + na), b->desc.curves, nb * sizeof(qaws_curve*));
		for (i = 0; i < na + nb; i++)
		{
			fam[i] = i >= na;
			counts[i] = i < na ? a->seg_count[i] : b->seg_count[i - na];
		}
		memcpy(segs, a->segs, a->nseg * sizeof(qaws_flat_seg));
		memcpy(segs + a->nseg, b->segs, b->nseg * sizeof(qaws_flat_seg));
		for (i = a->nseg; i < a->nseg + b->nseg; i++)
			segs[i].owner += na;
		memset(&d, 0, sizeof(d));
		d.curves = curves;
		d.curve_count = na + nb;
		d.families = fam;
		x.desc = &d;
		x.segs = segs;
		x.nseg = a->nseg + b->nseg;
		x.seg_count = counts;
		x.flat = a->flat > b->flat ? a->flat : b->flat;
		x.pos_tol = (a->ext > b->ext ? a->ext : b->ext) * CB_POS_REL;
	}
	s = cb_solve(&x);
	if (s == QAWS_STATUS_OK)
		cb_emit(&x, dim, out, capacity, out_count, b ? na : 0);
done:
	if (out_stats)
		*out_stats = x.stats;
	free(x.hits);
	free((void*)curves);
	free(fam);
	free(segs);
	free(counts);
	return s;
}

qaws_status qaws_curve_set_find_intersections_2d(qaws_curve_set const* set, qaws_curve_set const* other, qaws_curve_batch_hit_2d* out_hits,
	unsigned int hit_capacity, unsigned int* out_count, qaws_curve_batch_stats* out_stats)
{
	return cs_find(set, other, 2, out_hits, hit_capacity, out_count, out_stats);
}

qaws_status qaws_curve_set_find_intersections_3d(qaws_curve_set const* set, qaws_curve_set const* other, qaws_curve_batch_hit_3d* out_hits,
	unsigned int hit_capacity, unsigned int* out_count, qaws_curve_batch_stats* out_stats)
{
	return cs_find(set, other, 3, out_hits, hit_capacity, out_count, out_stats);
}
