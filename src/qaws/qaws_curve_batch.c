#include "qaws_curve_batch.h"
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
#define CB_STEP_REL  ((qaws_scalar)1e-7)   /* Newton steps below it: converged */
#else
#define CB_POS_REL   ((qaws_scalar)1e-10)
#define CB_PAR_REL   ((qaws_scalar)1e-7)
#define CB_STEP_REL  ((qaws_scalar)1e-15)
#endif

#define CB_PARALLEL_SIN ((qaws_scalar)0.25)  /* chords this parallel may overlap */
#define CB_TANGENT_SIN  ((qaws_scalar)1e-6)  /* hits this tangent are classified by their sides */
#if QAWS_SCALAR_IS_FLOAT
#define CB_SIDE_MAX     ((qaws_scalar)1e-1)  /* farthest step off a contact to see its sides, and widest zone */
#else
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

/* Newton on C_a(ta) = C_b(tb); in 3D the 3x2 system by its normal equations.
   *sin2 receives the squared sine of the angle between the tangents. */
static int cb_newton(cb_ctx* x, qaws_curve const* ca, qaws_curve const* cb, qaws_scalar* ta, qaws_scalar* tb, qaws_scalar* pos,
	qaws_scalar* sin2)
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
		/* J = [da, -db]; solve (J^T J) d = J^T f */
		m00 = cb_dot(da, da);
		m01 = -cb_dot(da, db);
		m11 = cb_dot(db, db);
		det = m00 * m11 - m01 * m01;
		*sin2 = m00 * m11 > 0 ? det / (m00 * m11) : 0;
		/* converged: two more steps take a shallow crossing to rounding level */
		if (sqrt(cb_dot(f, f)) <= x->pos_tol && polish++ == 2)
		{
			for (k = 0; k < 3; k++)
				pos[k] = (pa[k] + pb[k]) / 2;
			return 1;
		}
		g0 = cb_dot(da, f);
		g1 = -cb_dot(db, f);
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

static qaws_status cb_push(cb_ctx* x, unsigned int a, unsigned int b, unsigned int kind,
	qaws_scalar ta, qaws_scalar tb, qaws_scalar ta1, qaws_scalar tb1, qaws_scalar const* p)
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
	if (a > b || (a == b && (ta < ta1 ? ta : ta1) > (tb < tb1 ? tb : tb1)))
	{
		unsigned int ti = a; qaws_scalar tt = ta, tt1 = ta1;
		a = b; b = ti;
		ta = tb; tb = tt;
		ta1 = tb1; tb1 = tt1;
	}
	if (ta1 < ta)
	{
		/* an overlap runs up curve a */
		qaws_scalar t = ta; ta = ta1; ta1 = t;
		t = tb; tb = tb1; tb1 = t;
	}
	h->a = a; h->b = b; h->kind = kind; h->s2 = 0;
	h->zlo = h->zhi = ta;
	h->ta = ta; h->tb = tb; h->ta1 = ta1; h->tb1 = tb1;
	h->p[0] = p[0]; h->p[1] = p[1]; h->p[2] = p[2];
	return QAWS_STATUS_OK;
}

static qaws_scalar cb_par_tol(qaws_curve const* c)
{
	return (c->parameter_range.max_value - c->parameter_range.min_value) * CB_PAR_REL;
}

/* closest point of c on [lo, hi] to q, by Gauss-Newton on (C - q).C' = 0 from
   *t; returns the distance, or -1 when the curve cannot be evaluated */
static qaws_scalar cb_project(cb_ctx const* x, qaws_curve const* c, qaws_scalar const* q,
	qaws_scalar lo, qaws_scalar hi, qaws_scalar* t)
{
	qaws_scalar p[3], d[3], f[3];
	unsigned int it, k;
	for (it = 0; it < CB_NEWTON_ITERS; it++)
	{
		qaws_scalar dd, step;
		if (qaws_internal_curve_point(c, x->dim, *t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, p, d) != QAWS_STATUS_OK)
			return -1;
		for (k = 0; k < 3; k++)
			f[k] = p[k] - q[k];
		dd = cb_dot(d, d);
		if (!(dd > 0))
			break;
		step = cb_dot(f, d) / dd;
		*t -= step;
		if (*t < lo) *t = lo;
		if (*t > hi) *t = hi;
		if (fabs(step) <= (hi - lo) * CB_STEP_REL)
			break;
	}
	if (qaws_internal_curve_point(c, x->dim, *t, QAWS_EVAL_FLAG_POSITION, p, NULL) != QAWS_STATUS_OK)
		return -1;
	for (k = 0; k < 3; k++)
		f[k] = p[k] - q[k];
	return (qaws_scalar)sqrt(cb_dot(f, f));
}

/* Do the curves of nearly parallel segments A and B share a stretch? Its
   ends are segment ends lying on the other curve; three points between them
   must lie on the other curve too. Returns 1 and pushes the overlap. */
static int cb_overlap(cb_ctx* x, qaws_flat_seg const* A, qaws_flat_seg const* B, qaws_status* st)
{
	qaws_curve const* ca = x->desc->curves[A->owner];
	qaws_curve const* cb = x->desc->curves[B->owner];
	qaws_scalar da[3], db[3], cr[3], na, nb, tol = 16 * x->pos_tol;
	qaws_scalar cta[4], ctb[4], sp[3], ep[3], gap[3];
	unsigned int nc = 0, lo = 0, hi = 0, i, k;
	for (k = 0; k < 3; k++)
	{
		da[k] = A->p1[k] - A->p0[k];
		db[k] = B->p1[k] - B->p0[k];
	}
	cr[0] = da[1] * db[2] - da[2] * db[1];
	cr[1] = da[2] * db[0] - da[0] * db[2];
	cr[2] = da[0] * db[1] - da[1] * db[0];
	na = (qaws_scalar)sqrt(cb_dot(da, da));
	nb = (qaws_scalar)sqrt(cb_dot(db, db));
	if (!(na > 0 && nb > 0) || (qaws_scalar)sqrt(cb_dot(cr, cr)) > CB_PARALLEL_SIN * na * nb)
		return 0;
	/* segment ends lying on the other curve */
	for (i = 0; i < 4; i++)
	{
		qaws_flat_seg const* S = i < 2 ? A : B;
		qaws_flat_seg const* O = i < 2 ? B : A;
		qaws_curve const* oc = i < 2 ? cb : ca;
		qaws_scalar const* q = (i & 1) ? S->p1 : S->p0;
		qaws_scalar w[3], u, t, dist;
		for (k = 0; k < 3; k++)
			w[k] = q[k] - O->p0[k];
		u = cb_dot(w, i < 2 ? db : da) / (i < 2 ? nb * nb : na * na);
		u = u < 0 ? 0 : (u > 1 ? 1 : u);
		t = O->t0 + u * (O->t1 - O->t0);
		dist = cb_project(x, oc, q, O->t0 < O->t1 ? O->t0 : O->t1, O->t0 < O->t1 ? O->t1 : O->t0, &t);
		if (dist < 0 || dist > tol)
			continue;
		cta[nc] = i < 2 ? ((i & 1) ? A->t1 : A->t0) : t;
		ctb[nc] = i < 2 ? t : ((i & 1) ? B->t1 : B->t0);
		nc++;
	}
	if (nc < 2)
		return 0;
	for (i = 1; i < nc; i++)
	{
		if (cta[i] < cta[lo]) lo = i;
		if (cta[i] > cta[hi]) hi = i;
	}
	if (qaws_internal_curve_point(ca, x->dim, cta[lo], QAWS_EVAL_FLAG_POSITION, sp, NULL) != QAWS_STATUS_OK ||
		qaws_internal_curve_point(ca, x->dim, cta[hi], QAWS_EVAL_FLAG_POSITION, ep, NULL) != QAWS_STATUS_OK)
		return 0;
	for (k = 0; k < 3; k++)
		gap[k] = ep[k] - sp[k];
	if ((qaws_scalar)sqrt(cb_dot(gap, gap)) <= 4 * tol)
		return 0;
	/* the stretch between must lie on curve b */
	for (i = 1; i <= 3; i++)
	{
		qaws_scalar f = (qaws_scalar)i / 4, q[3];
		qaws_scalar t = ctb[lo] + f * (ctb[hi] - ctb[lo]), dist;
		qaws_scalar blo = ctb[lo] < ctb[hi] ? ctb[lo] : ctb[hi], bhi = ctb[lo] < ctb[hi] ? ctb[hi] : ctb[lo];
		if (qaws_internal_curve_point(ca, x->dim, cta[lo] + f * (cta[hi] - cta[lo]), QAWS_EVAL_FLAG_POSITION, q, NULL) != QAWS_STATUS_OK)
			return 0;
		dist = cb_project(x, cb, q, blo, bhi, &t);
		if (dist < 0 || dist > tol)
			return 0;
	}
	*st = cb_push(x, A->owner, B->owner, QAWS_CURVE_HIT_OVERLAP, cta[lo], ctb[lo], cta[hi], ctb[hi], sp);
	return 1;
}

static qaws_status cb_pair(cb_ctx* x, qaws_flat_seg const* A, qaws_flat_seg const* B)
{
	qaws_curve const* ca = x->desc->curves[A->owner];
	qaws_curve const* cb = x->desc->curves[B->owner];
	qaws_scalar u, v, ta, tb, pos[3], s2 = 0;
	qaws_status st = QAWS_STATUS_OK;
	if (cb_chords(A, B, &u, &v) > A->r + B->r + x->pos_tol)
		return QAWS_STATUS_OK;
	if (cb_overlap(x, A, B, &st))
		return st;
	x->stats.newton_count++;
	ta = A->t0 + u * (A->t1 - A->t0);
	tb = B->t0 + v * (B->t1 - B->t0);
	if (!cb_newton(x, ca, cb, &ta, &tb, pos, &s2))
		return QAWS_STATUS_OK;
	/* a curve meets itself trivially at ta == tb */
	if (A->owner == B->owner && fabs(ta - tb) <= 64 * cb_par_tol(ca))
		return QAWS_STATUS_OK;
	st = cb_push(x, A->owner, B->owner, QAWS_CURVE_HIT_CROSSING, ta, tb, ta, tb, pos);
	if (st == QAWS_STATUS_OK)
		x->hits[x->nhit - 1].s2 = s2;
	return st;
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

/* signed distance of point q from curve c, searched from parameter t (left
   of the tangent is positive); 0 when it cannot be evaluated */
static qaws_scalar cb_side(cb_ctx const* x, qaws_curve const* c, qaws_scalar t, qaws_scalar const* q)
{
	qaws_scalar p[3], d[3], lo = c->parameter_range.min_value, hi = c->parameter_range.max_value, n;
	if (cb_project(x, c, q, lo, hi, &t) < 0 ||
		qaws_internal_curve_point(c, x->dim, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, p, d) != QAWS_STATUS_OK)
		return 0;
	n = (qaws_scalar)sqrt(d[0] * d[0] + d[1] * d[1]);
	return n > 0 ? (d[0] * (q[1] - p[1]) - d[1] * (q[0] - p[0])) / n : 0;
}

static int cb_at_open_end(qaws_curve const* c, qaws_scalar t, qaws_scalar tol)
{
	return (t <= c->parameter_range.min_value + tol || t >= c->parameter_range.max_value - tol) && !qaws_curve_is_closed(c);
}

/* 2D point hits. The end of an open curve lying on the other touches it.
   A transversal hit crosses. At a tangency, or across a stretch where the
   curves stay within tolerance (the zone [zlo, zhi] on curve a), step along
   curve a off both ends until it leaves the tolerance and compare the sides
   of curve b it lies on. */
static unsigned int cb_classify(cb_ctx const* x, cb_hit const* h)
{
	qaws_curve const* ca = x->desc->curves[h->a];
	qaws_curve const* cb = x->desc->curves[h->b];
	qaws_scalar da[3], q[3], la, ext = x->pos_tol / CB_POS_REL, side[2];
	qaws_scalar alo = ca->parameter_range.min_value, ahi = ca->parameter_range.max_value;
	int s;
	if (x->dim != 2)
		return QAWS_CURVE_HIT_CROSSING;
	if (cb_at_open_end(ca, h->ta, cb_par_tol(ca)) || cb_at_open_end(cb, h->tb, cb_par_tol(cb)))
		return QAWS_CURVE_HIT_TOUCH;
	if (h->zlo == h->zhi && h->s2 > CB_TANGENT_SIN * CB_TANGENT_SIN)
		return QAWS_CURVE_HIT_CROSSING;
	if (qaws_internal_curve_point(ca, 2, h->ta, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, q, da) != QAWS_STATUS_OK)
		return QAWS_CURVE_HIT_CROSSING;
	la = (qaws_scalar)sqrt(da[0] * da[0] + da[1] * da[1]);
	if (!(la > 0))
		return QAWS_CURVE_HIT_TOUCH;
	for (s = 0; s < 2; s++)
	{
		qaws_scalar dist = ext * (qaws_scalar)1e-6 > 64 * x->pos_tol ? ext * (qaws_scalar)1e-6 : 64 * x->pos_tol;
		side[s] = 0;
		for (; dist <= ext * CB_SIDE_MAX; dist *= 2)
		{
			qaws_scalar t = s ? h->zhi + dist / la : h->zlo - dist / la;
			if (t < alo || t > ahi)
			{
				/* closed curve a: wrap; open: one side only */
				if (!qaws_curve_is_closed(ca))
					return QAWS_CURVE_HIT_TOUCH;
				t += t < alo ? ahi - alo : alo - ahi;
			}
			if (qaws_internal_curve_point(ca, 2, t, QAWS_EVAL_FLAG_POSITION, q, NULL) != QAWS_STATUS_OK)
				return QAWS_CURVE_HIT_CROSSING;
			side[s] = cb_side(x, cb, h->tb, q);
			if (fabs(side[s]) > 4 * x->pos_tol)
				break;
		}
	}
	return fabs(side[0]) > 4 * x->pos_tol && fabs(side[1]) > 4 * x->pos_tol && (side[0] < 0) != (side[1] < 0)
		? QAWS_CURVE_HIT_CROSSING : QAWS_CURVE_HIT_TOUCH;
}

/* do the curves of hits g and h (one pair, g.ta < h.ta) stay within
   tolerance between them? */
static int cb_same_zone(cb_ctx* x, cb_hit const* g, cb_hit const* h)
{
	qaws_curve const* ca = x->desc->curves[g->a];
	qaws_curve const* cb = x->desc->curves[g->b];
	qaws_scalar d[3], ext = x->pos_tol / CB_POS_REL;
	unsigned int i, k;
	for (k = 0; k < 3; k++)
		d[k] = h->p[k] - g->p[k];
	if (sqrt(cb_dot(d, d)) > ext * CB_SIDE_MAX)
		return 0;
	for (i = 1; i <= 3; i++)
	{
		qaws_scalar f = (qaws_scalar)i / 4, q[3], t = g->tb + f * (h->tb - g->tb), dist;
		qaws_scalar blo = g->tb < h->tb ? g->tb : h->tb, bhi = g->tb < h->tb ? h->tb : g->tb;
		if (qaws_internal_curve_point(ca, x->dim, g->zhi + f * (h->ta - g->zhi), QAWS_EVAL_FLAG_POSITION, q, NULL) != QAWS_STATUS_OK)
			return 0;
		dist = cb_project(x, cb, q, blo, bhi, &t);
		if (dist < 0 || dist > 16 * x->pos_tol)
			return 0;
	}
	return 1;
}

/* is t within [lo, hi] (tolerance tol), also across the seam of a closed curve? */
static int cb_in(qaws_curve const* c, int closed, qaws_scalar t, qaws_scalar lo, qaws_scalar hi, qaws_scalar tol)
{
	qaws_scalar len = c->parameter_range.max_value - c->parameter_range.min_value;
	if (t >= lo - tol && t <= hi + tol)
		return 1;
	return closed && ((t + len >= lo - tol && t + len <= hi + tol) || (t - len >= lo - tol && t - len <= hi + tol));
}

/* overlaps of one curve pair joined into maximal stretches; point hits kept
   once per point and per contact zone, dropped on a stretch, and classified */
static qaws_status cb_merge(cb_ctx* x)
{
	unsigned int i, n = 0, start = 0;
	cb_hit* src;
	qsort(x->hits, x->nhit, sizeof(cb_hit), cb_cmp_hit);
	src = (cb_hit*)malloc((x->nhit ? x->nhit : 1) * sizeof(cb_hit));
	if (!src)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memcpy(src, x->hits, x->nhit * sizeof(cb_hit));
	while (start < x->nhit)
	{
		unsigned int end = start, j, novl = 0, first, group;
		qaws_curve const* ca = x->desc->curves[src[start].a];
		qaws_curve const* cb = x->desc->curves[src[start].b];
		qaws_scalar pa = cb_par_tol(ca), pb = cb_par_tol(cb);
		int closed_a = 0, closed_b = 0;
		while (end < x->nhit && src[end].a == src[start].a && src[end].b == src[start].b)
			end++;
		/* 1. overlaps, joined while they meet on curve a */
		first = n;
		for (j = start; j < end; j++)
		{
			cb_hit const* h = &src[j];
			if (h->kind != QAWS_CURVE_HIT_OVERLAP)
				continue;
			if (novl && h->ta <= x->hits[n - 1].ta1 + pa)
			{
				cb_hit* g = &x->hits[n - 1];
				if (h->ta1 > g->ta1)
				{
					g->ta1 = h->ta1;
					g->tb1 = h->tb1;
				}
				continue;
			}
			x->hits[n++] = *h;
			novl++;
		}
		if (novl)
		{
			closed_a = qaws_curve_is_closed(ca);
			closed_b = qaws_curve_is_closed(cb);
		}
		/* 2. point hits: none on a stretch, once per point, once per zone */
		group = n;
		for (j = start; j < end; j++)
		{
			cb_hit h = src[j];
			unsigned int k;
			int drop = 0;
			if (h.kind == QAWS_CURVE_HIT_OVERLAP)
				continue;
			h.zlo = h.zhi = h.ta;
			for (k = first; k < first + novl && !drop; k++)
			{
				cb_hit const* o = &x->hits[k];
				qaws_scalar blo = o->tb < o->tb1 ? o->tb : o->tb1, bhi = o->tb < o->tb1 ? o->tb1 : o->tb;
				drop = cb_in(ca, closed_a, h.ta, o->ta, o->ta1, pa) && cb_in(cb, closed_b, h.tb, blo, bhi, pb);
			}
			for (k = group; k < n && !drop; k++)
			{
				cb_hit const* g = &x->hits[k];
				qaws_scalar d[3];
				int c;
				for (c = 0; c < 3; c++)
					d[c] = g->p[c] - h.p[c];
				drop = (fabs(g->ta - h.ta) <= pa && fabs(g->tb - h.tb) <= pb) || sqrt(cb_dot(d, d)) <= 16 * x->pos_tol;
			}
			if (drop)
				continue;
			if (x->dim == 2 && n > group && cb_same_zone(x, &x->hits[n - 1], &h))
			{
				/* one contact: keep the hit nearest the middle of the zone */
				cb_hit* g = &x->hits[n - 1];
				qaws_scalar mid = (g->zlo + h.ta) / 2;
				if (fabs(h.ta - mid) < fabs(g->ta - mid))
				{
					qaws_scalar zlo = g->zlo;
					*g = h;
					g->zlo = zlo;
				}
				g->zhi = h.ta;
				continue;
			}
			x->hits[n++] = h;
		}
		start = end;
	}
	x->nhit = n;
	for (i = 0; i < n; i++)
		if (x->hits[i].kind != QAWS_CURVE_HIT_OVERLAP)
			x->hits[i].kind = cb_classify(x, &x->hits[i]);
	qsort(x->hits, x->nhit, sizeof(cb_hit), cb_cmp_hit);
	free(src);
	return QAWS_STATUS_OK;
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
static qaws_status cb_visit_chunk(void* user, unsigned int chunk, unsigned int i, unsigned int j)
{
	cb_ctx* ch = &((cb_ctx*)user)[chunk];
	return cb_pair(ch, &ch->segs[i], &ch->segs[j]);
}

/* 2-5. one grid over every segment; each overlapping allowed pair refined
   (cells in chunks, each chunk with its own hit list); merge */
static qaws_status cb_solve(cb_ctx* x)
{
	unsigned int i, k, n, cand = 0;
	qaws_status s;
	qaws_bp_box* boxes;
	qaws_bp_grid* grid = NULL;
	cb_ctx* ch = NULL;
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
	s = qaws_internal_grid_create(boxes, x->nseg, x->dim, &grid);
	n = s == QAWS_STATUS_OK ? qaws_internal_grid_pair_chunks(grid) : 0;
	if (s == QAWS_STATUS_OK && n)
	{
		ch = (cb_ctx*)malloc(n * sizeof(cb_ctx));
		if (!ch)
			s = QAWS_STATUS_ALLOCATION_FAILURE;
		for (k = 0; s == QAWS_STATUS_OK && k < n; k++)
		{
			ch[k] = *x;
			ch[k].hits = NULL;
			ch[k].nhit = ch[k].caphit = 0;
			memset(&ch[k].stats, 0, sizeof(ch[k].stats));
		}
		if (s == QAWS_STATUS_OK)
			s = qaws_internal_grid_pairs(grid, boxes, cb_accept, x, cb_visit_chunk, ch, x->desc->executor, &cand);
		for (k = 0; ch && k < n; k++)
		{
			if (s == QAWS_STATUS_OK && ch[k].nhit)
			{
				cb_hit* g = (cb_hit*)realloc(x->hits, (x->nhit + ch[k].nhit) * sizeof(cb_hit));
				if (!g)
					s = QAWS_STATUS_ALLOCATION_FAILURE;
				else
				{
					memcpy(g + x->nhit, ch[k].hits, ch[k].nhit * sizeof(cb_hit));
					x->hits = g;
					x->nhit += ch[k].nhit;
					x->caphit = x->nhit;
				}
			}
			x->stats.newton_count += ch[k].stats.newton_count;
			free(ch[k].hits);
		}
	}
	x->stats.cell_count = grid ? qaws_internal_grid_cells(grid) : 0;
	x->stats.candidate_count = cand;
	free(ch);
	free(boxes);
	qaws_internal_grid_destroy(grid);
	if (s == QAWS_STATUS_OK)
		s = cb_merge(x);
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
			o->kind = h->kind; o->parameter_a_end = h->ta1; o->parameter_b_end = h->tb1;
			o->position.x = h->p[0]; o->position.y = h->p[1];
		}
		else
		{
			qaws_curve_batch_hit_3d* o = (qaws_curve_batch_hit_3d*)out + i;
			o->curve_a = h->a; o->curve_b = h->b - b_offset;
			o->parameter_a = h->ta; o->parameter_b = h->tb;
			o->kind = h->kind; o->parameter_a_end = h->ta1; o->parameter_b_end = h->tb1;
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
		d.executor = a->desc.executor;
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

/* ------------------------------------------------------------------ */
/*  Level crossings                                                    */
/* ------------------------------------------------------------------ */

#define LC_MAX_DEPTH 16
#define LC_ITERS 64

typedef struct lc_ctx
{
	qaws_level_crossing_desc const* desc;
	unsigned int dim, curve;
	int has_gradient;
	qaws_scalar htol, ftol;
	qaws_level_crossing* out;
	unsigned int n, cap;
} lc_ctx;

/* h at the curve point of t, and dh/dt when wanted (the gradient dotted with dC/dt) */
static qaws_status lc_eval(lc_ctx* x, qaws_scalar t, qaws_scalar* p, qaws_scalar* h, qaws_scalar* dh)
{
	qaws_curve const* c = x->desc->curves[x->curve];
	qaws_scalar d[3], g[3];
	qaws_status s = qaws_internal_curve_point(c, x->dim, t, QAWS_EVAL_FLAG_POSITION | (dh ? QAWS_EVAL_FLAG_D1 : 0), p, dh ? d : NULL);
	if (s != QAWS_STATUS_OK)
		return s;
	g[0] = g[1] = g[2] = 0;
	*h = x->desc->field(x->desc->user, p, dh && x->has_gradient ? g : NULL);
	if (dh)
		*dh = g[0] * d[0] + g[1] * d[1] + g[2] * d[2];
	return QAWS_STATUS_OK;
}

static qaws_status lc_push(lc_ctx* x, unsigned int level, qaws_scalar t, qaws_scalar const* p)
{
	qaws_level_crossing* o;
	if (x->n == x->cap)
	{
		unsigned int cap = x->cap ? x->cap * 2 : 256;
		qaws_level_crossing* g = (qaws_level_crossing*)realloc(x->out, cap * sizeof(qaws_level_crossing));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		x->out = g;
		x->cap = cap;
	}
	o = &x->out[x->n++];
	o->curve = x->curve;
	o->level = level;
	o->parameter = t;
	o->position.x = p[0];
	o->position.y = p[1];
	o->position.z = p[2];
	return QAWS_STATUS_OK;
}

/* h(C(t)) = L on [ta, tb], fa = h - L at ta and fb at tb of opposite signs
   (or fa = 0): Newton (secant without a gradient) kept inside the bracket */
static qaws_status lc_solve(lc_ctx* x, unsigned int level, qaws_scalar ta, qaws_scalar fa, qaws_scalar tb, qaws_scalar fb)
{
	qaws_scalar L = x->desc->levels[level], t, p[3], h, dh = 0, f;
	unsigned int it;
	int side = 0;
	qaws_status s;
	if (fa == 0)
	{
		s = lc_eval(x, ta, p, &h, NULL);
		return s == QAWS_STATUS_OK ? lc_push(x, level, ta, p) : s;
	}
	t = ta + (tb - ta) * fa / (fa - fb);
	for (it = 0; it < LC_ITERS; it++)
	{
		qaws_scalar tn, lo, hi;
		s = lc_eval(x, t, p, &h, x->has_gradient ? &dh : NULL);
		if (s != QAWS_STATUS_OK)
			return s;
		f = h - L;
		if (fabs(f) <= x->ftol)
			break;
		/* shrink the bracket (Illinois weights for the secant) */
		if ((f < 0) == (fa < 0))
		{
			ta = t;
			fa = f;
			if (side == -1) fb /= 2;
			side = -1;
		}
		else
		{
			tb = t;
			fb = f;
			if (side == 1) fa /= 2;
			side = 1;
		}
		lo = ta < tb ? ta : tb;
		hi = ta < tb ? tb : ta;
		if (!(hi - lo > (fabs(lo) + fabs(hi)) * (qaws_scalar)1e-16))
			break;
		tn = x->has_gradient && dh != 0 ? t - f / dh : ta + (tb - ta) * fa / (fa - fb);
		if (!(tn > lo && tn < hi))
			tn = (lo + hi) / 2;
		t = tn;
	}
	return lc_push(x, level, t, p);
}

/* first level index with levels[i] >= v */
static unsigned int lc_lower(qaws_scalar const* lv, unsigned int n, qaws_scalar v)
{
	unsigned int lo = 0, hi = n;
	while (lo < hi)
	{
		unsigned int mid = (lo + hi) / 2;
		if (lv[mid] < v)
			lo = mid + 1;
		else
			hi = mid;
	}
	return lo;
}

/* piece [t0, t1] with field values h0, h1: split while h is not close to
   linear along it, then every level in [min, max) is crossed once */
static qaws_status lc_piece(lc_ctx* x, qaws_scalar t0, qaws_scalar h0, qaws_scalar t1, qaws_scalar h1, unsigned int depth)
{
	qaws_scalar const* lv = x->desc->levels;
	unsigned int n = x->desc->level_count, i, i0, i1;
	qaws_scalar lo = h0 < h1 ? h0 : h1, hi = h0 < h1 ? h1 : h0, tm = (t0 + t1) / 2, p[3], hm;
	qaws_status s;
	if (depth < LC_MAX_DEPTH)
	{
		s = lc_eval(x, tm, p, &hm, NULL);
		if (s != QAWS_STATUS_OK)
			return s;
		if (fabs(hm - (h0 + h1) / 2) > x->htol)
		{
			s = lc_piece(x, t0, h0, tm, hm, depth + 1);
			return s == QAWS_STATUS_OK ? lc_piece(x, tm, hm, t1, h1, depth + 1) : s;
		}
	}
	/* levels in [lo, hi): a value exactly on a level belongs to the piece
	   where it is the lower end, so a crossing through it is counted once */
	i0 = lc_lower(lv, n, lo);
	i1 = lc_lower(lv, n, hi);
	for (i = 0; i < i1 - i0; i++)
	{
		/* in the order of t */
		unsigned int k = h0 <= h1 ? i0 + i : i1 - 1 - i;
		s = lc_solve(x, k, t0, h0 - lv[k], t1, h1 - lv[k]);
		if (s != QAWS_STATUS_OK)
			return s;
	}
	return QAWS_STATUS_OK;
}

#define LC_GRAIN 64

typedef struct lc_job
{
	qaws_flat_seg const* segs;
	lc_ctx* ctx;                /* one per chunk */
	unsigned int nchunk;
} lc_job;

static qaws_status lc_chunk(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	lc_job const* J = (lc_job const*)ctx;
	lc_ctx* x = &J->ctx[chunk];
	unsigned int k;
	qaws_status s = QAWS_STATUS_OK;
	x->out = NULL;
	x->n = x->cap = 0;
	for (k = begin; k < end && s == QAWS_STATUS_OK; k++)
	{
		qaws_scalar p[3], h0 = 0, h1 = 0;
		x->curve = J->segs[k].owner;
		s = lc_eval(x, J->segs[k].t0, p, &h0, NULL);
		if (s == QAWS_STATUS_OK)
			s = lc_eval(x, J->segs[k].t1, p, &h1, NULL);
		if (s == QAWS_STATUS_OK)
			s = lc_piece(x, J->segs[k].t0, h0, J->segs[k].t1, h1, 0);
	}
	return s;
}

static int lc_cmp(void const* p, void const* q)
{
	qaws_level_crossing const* a = (qaws_level_crossing const*)p;
	qaws_level_crossing const* b = (qaws_level_crossing const*)q;
	if (a->curve != b->curve) return a->curve < b->curve ? -1 : 1;
	if (a->parameter != b->parameter) return a->parameter < b->parameter ? -1 : 1;
	return 0;
}

qaws_status qaws_curve_batch_find_level_crossings(
	qaws_level_crossing_desc const* desc,
	qaws_level_crossing* out_crossings,
	unsigned int capacity,
	unsigned int* out_count)
{
	lc_ctx x;
	qaws_flat_seg* segs = NULL;
	unsigned int nseg = 0, capseg = 0, i, k;
	qaws_scalar ext, flat, hmin = (qaws_scalar)HUGE_VAL, hmax = -(qaws_scalar)HUGE_VAL, dl = (qaws_scalar)HUGE_VAL, g[3];
	qaws_status s = QAWS_STATUS_OK;
	if (!desc || !desc->field || !out_count || (!out_crossings && capacity) || (desc->curve_count && !desc->curves) || (desc->level_count && !desc->levels))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	for (i = 1; i < desc->level_count; i++)
	{
		if (!(desc->levels[i] > desc->levels[i - 1]))
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (desc->levels[i] - desc->levels[i - 1] < dl)
			dl = desc->levels[i] - desc->levels[i - 1];
	}
	for (i = 0; i < desc->curve_count; i++)
	{
		if (!desc->curves[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (desc->curves[i]->dimension != desc->curves[0]->dimension)
			return QAWS_STATUS_INVALID_DIMENSION;
	}
	if (!desc->curve_count || !desc->level_count)
		return QAWS_STATUS_OK;
	memset(&x, 0, sizeof(x));
	x.desc = desc;
	x.dim = desc->curves[0]->dimension == QAWS_DIMENSION_2D ? 2 : 3;
	ext = qaws_internal_flatten_extent(desc->curves, desc->curve_count, x.dim, NULL, 0);
	flat = desc->flatness > 0 ? desc->flatness : ext / 1024;
	/* does the field fill the gradient? */
	{
		qaws_scalar p[3] = { 0, 0, 0 };
		qaws_internal_curve_point(desc->curves[0], x.dim, desc->curves[0]->parameter_range.min_value, QAWS_EVAL_FLAG_POSITION, p, NULL);
		g[0] = g[1] = g[2] = (qaws_scalar)HUGE_VAL;
		desc->field(desc->user, p, g);
		x.has_gradient = g[0] != (qaws_scalar)HUGE_VAL && g[1] != (qaws_scalar)HUGE_VAL;
	}
	for (i = 0; i < desc->curve_count && s == QAWS_STATUS_OK; i++)
	{
		unsigned int first = nseg;
		s = qaws_internal_flatten_curve(desc->curves[i], x.dim, flat, i, &segs, &nseg, &capseg);
		/* field range over the segment ends: scales the tolerances */
		for (k = first; k < nseg && s == QAWS_STATUS_OK; k++)
		{
			qaws_scalar h = desc->field(desc->user, segs[k].p0, NULL);
			if (h < hmin) hmin = h;
			if (h > hmax) hmax = h;
		}
	}
	if (s == QAWS_STATUS_OK && nseg)
	{
		qaws_scalar hr = hmax > hmin ? hmax - hmin : 1, scale = hr + (qaws_scalar)fabs(desc->levels[0]) + (qaws_scalar)fabs(desc->levels[desc->level_count - 1]);
		lc_job J;
		x.htol = (dl < hr ? dl : hr) / 8;
		x.ftol = scale * (QAWS_SCALAR_IS_FLOAT ? (qaws_scalar)1e-6 : (qaws_scalar)1e-14);
		/* segments in chunks, each with its own copy of the context and list */
		J.segs = segs;
		J.nchunk = qaws_internal_chunk_count(nseg, LC_GRAIN);
		J.ctx = (lc_ctx*)malloc(J.nchunk * sizeof(lc_ctx));
		if (!J.ctx)
			s = QAWS_STATUS_ALLOCATION_FAILURE;
		for (k = 0; s == QAWS_STATUS_OK && k < J.nchunk; k++)
			J.ctx[k] = x;
		if (s == QAWS_STATUS_OK)
			s = qaws_internal_parallel(desc->executor, nseg, LC_GRAIN, lc_chunk, &J);
		/* the chunks' lists together */
		for (k = 0; J.ctx && k < J.nchunk; k++)
		{
			if (s == QAWS_STATUS_OK && J.ctx[k].n)
			{
				unsigned int m;
				for (m = 0; m < J.ctx[k].n && s == QAWS_STATUS_OK; m++)
				{
					x.curve = J.ctx[k].out[m].curve;
					if (x.n == x.cap)
					{
						unsigned int cap = x.cap ? x.cap * 2 : 256;
						qaws_level_crossing* grown = (qaws_level_crossing*)realloc(x.out, cap * sizeof(qaws_level_crossing));
						if (!grown)
						{
							s = QAWS_STATUS_ALLOCATION_FAILURE;
							break;
						}
						x.out = grown;
						x.cap = cap;
					}
					x.out[x.n++] = J.ctx[k].out[m];
				}
			}
			free(J.ctx[k].out);
		}
		free(J.ctx);
	}
	if (s == QAWS_STATUS_OK)
	{
		qsort(x.out, x.n, sizeof(qaws_level_crossing), lc_cmp);
		for (i = 0; i < x.n && i < capacity; i++)
			out_crossings[i] = x.out[i];
		*out_count = x.n;
	}
	free(segs);
	free(x.out);
	return s;
}

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
		unsigned int cap = x->capcand ? x->capcand * 2 : 64;
		cp_cand* g = (cp_cand*)realloc(x->cand, cap * sizeof(cp_cand));
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
		s = qaws_curve_evaluate_2d(c, t, flags, &r);
		P[0] = r.position.x; P[1] = r.position.y; P[2] = 0;
		D1[0] = r.d1.x; D1[1] = r.d1.y; D1[2] = 0;
		D2[0] = r.d2.x; D2[1] = r.d2.y; D2[2] = 0;
	}
	else
	{
		qaws_eval_result_3d r;
		s = qaws_curve_evaluate_3d(c, t, flags, &r);
		P[0] = r.position.x; P[1] = r.position.y; P[2] = r.position.z;
		D1[0] = r.d1.x; D1[1] = r.d1.y; D1[2] = r.d1.z;
		D2[0] = r.d2.x; D2[1] = r.d2.y; D2[2] = r.d2.z;
	}
	return s;
}

/* Newton on g(t) = (C(t) - p) . C'(t) from t, clamped to the domain */
static qaws_scalar cp_refine(qaws_curve const* c, unsigned int dim, qaws_scalar const* p, qaws_scalar* t, qaws_scalar* pos)
{
	qaws_scalar tmin = c->parameter_range.min_value, tmax = c->parameter_range.max_value, P[3], D1[3], D2[3], w[3];
	unsigned int it, k;
	for (it = 0; it < CP_ITERS; it++)
	{
		qaws_scalar g, gp, step;
		if (cp_eval2(c, dim, *t, P, D1, D2) != QAWS_STATUS_OK)
			break;
		for (k = 0; k < 3; k++)
			w[k] = P[k] - p[k];
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
	qaws_internal_curve_point(c, dim, *t, QAWS_EVAL_FLAG_POSITION, P, NULL);
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
	x.stamp = (unsigned int*)malloc(J->nseg * sizeof(unsigned int));
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
	free(x.stamp);
	free(x.cand);
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
	boxes = (qaws_bp_box*)malloc(nseg * sizeof(qaws_bp_box));
	nchunk = qaws_internal_chunk_count(point_count, CP_GRAIN);
	J.stats = (qaws_curve_batch_stats*)calloc(nchunk, sizeof(qaws_curve_batch_stats));
	if (!boxes || !J.stats)
	{
		free(boxes);
		free(J.stats);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nseg; i++)
		for (k = 0; k < 3; k++)
		{
			boxes[i].lo[k] = (double)segs[i].lo[k];
			boxes[i].hi[k] = (double)segs[i].hi[k];
		}
	s = qaws_internal_grid_create(boxes, nseg, dim, &J.grid);
	free(boxes);
	J.curves = curves;
	J.segs = segs;
	J.nseg = nseg;
	J.dim = dim;
	J.points = points;
	J.max_distance = max_distance;
	J.out = out;
	if (s == QAWS_STATUS_OK)
		s = qaws_internal_parallel(executor, point_count, CP_GRAIN, cp_chunk, &J);
	for (k = 0; k < nchunk; k++)
	{
		stats->candidate_count += J.stats[k].candidate_count;
		stats->newton_count += J.stats[k].newton_count;
		stats->hit_count += J.stats[k].hit_count;
	}
	qaws_internal_grid_destroy(J.grid);
	free(J.stats);
	return s;
}

qaws_status qaws_curve_batch_find_closest(
	qaws_closest_desc const* desc,
	qaws_closest_point* out_points,
	qaws_curve_batch_stats* out_stats)
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
		s = qaws_internal_flatten_curve(desc->curves[i], dim, flat, i, &segs, &nseg, &capseg);
	if (s == QAWS_STATUS_OK)
		s = cp_run(desc->curves, dim, segs, nseg, desc->points, desc->point_count, desc->max_distance, desc->executor, out_points, &st);
	if (out_stats)
		*out_stats = st;
	free(segs);
	return s;
}

qaws_status qaws_curve_set_find_closest(
	qaws_curve_set const* set,
	qaws_scalar const* points,
	unsigned int point_count,
	qaws_scalar max_distance,
	qaws_closest_point* out_points,
	qaws_curve_batch_stats* out_stats)
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
