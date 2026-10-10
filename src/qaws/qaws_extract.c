/*
 * qaws_curve_extract: the piece of a curve between two parameters, as a new
 * standalone curve that lies exactly on the source.
 *
 *   Bezier, rational Bezier   same kind, de Casteljau (homogeneous when rational)
 *   B-spline, NURBS           same kind, knot insertion, parameters kept
 *   Polynomial                same kind, narrowed domain
 *   Arc, clothoid             same kind, new angles / origin and curvatures
 *   Composite                 one segment: that segment's piece; else a composite
 *   Hermite, Catmull-Rom,     cubic spans rebuilt as Bezier spans from four
 *   trajectory, subdivision   evaluations: one span gives a Bezier, more a
 *                             cubic B-spline whose knots are the source parameters
 *   Yuksel                    no polynomial form: cubic spans fitted to a relative
 *                             tolerance, the only kind that is not exact
 *   Reparameterized           the piece of its source curve
 *
 * t0 > t1 returns the piece reversed.
 */

#include "qaws_operations.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_bezier.h"
#include "qaws_rational_bezier.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_polynomial.h"
#include "qaws_arc.h"
#include "qaws_clothoid.h"
#include "qaws_composite.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_span.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if QAWS_SCALAR_IS_FLOAT
#define EX_REL 2e-6
#define EX_FIT_REL 1e-5
#else
#define EX_REL 1e-12
#define EX_FIT_REL 1e-9
#endif

#define EX_FIT_DEPTH 24

/* -------------------------------------------------------------------------- */
/*  Evaluation helpers                                                        */
/* -------------------------------------------------------------------------- */

static qaws_status ex_point(qaws_curve const* c, qaws_scalar t, qaws_scalar* out)
{
	qaws_status s;
	if (c->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION, &r);
		out[0] = r.position.x; out[1] = r.position.y; out[2] = QAWS_ZERO;
	}
	else
	{
		qaws_eval_result_3d r;
		s = qaws_curve_evaluate_3d(c, t, QAWS_EVAL_FLAG_POSITION, &r);
		out[0] = r.position.x; out[1] = r.position.y; out[2] = r.position.z;
	}
	return s;
}

static qaws_status ex_span_point(qaws_curve const* c, unsigned int span, qaws_scalar u, qaws_scalar* out)
{
	qaws_status s;
	if (c->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d r;
		s = c->vtable->eval_span_2d(c, span, u, QAWS_EVAL_FLAG_POSITION, &r);
		out[0] = r.position.x; out[1] = r.position.y; out[2] = QAWS_ZERO;
	}
	else
	{
		qaws_eval_result_3d r;
		s = c->vtable->eval_span_3d(c, span, u, QAWS_EVAL_FLAG_POSITION, &r);
		out[0] = r.position.x; out[1] = r.position.y; out[2] = r.position.z;
	}
	return s;
}

/* -------------------------------------------------------------------------- */
/*  Bezier and rational Bezier                                                */
/* -------------------------------------------------------------------------- */

/* Restricts the Bezier with control points cp (count points, stride scalars
   each) to [t0, t1] of its [0, 1] domain, in place. */
static void ex_bezier_restrict(qaws_scalar* cp, unsigned int count, unsigned int stride,
	qaws_scalar t0, qaws_scalar t1)
{
	unsigned int deg = count - 1, r, i, k;
	if (t1 < QAWS_ONE)
	{
		/* left part of the split at t1 */
		for (r = 1; r <= deg; r++)
			for (i = deg; i >= r; i--)
				for (k = 0; k < stride; k++)
					cp[i * stride + k] = (QAWS_ONE - t1) * cp[(i - 1) * stride + k] + t1 * cp[i * stride + k];
	}
	if (t0 > QAWS_ZERO)
	{
		/* right part of the split at t0 of the remaining piece */
		qaws_scalar s = t0 / t1;
		for (r = 1; r <= deg; r++)
			for (i = 0; i + r <= deg; i++)
				for (k = 0; k < stride; k++)
					cp[i * stride + k] = (QAWS_ONE - s) * cp[i * stride + k] + s * cp[(i + 1) * stride + k];
	}
}

static qaws_status ex_bezier(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, qaws_curve** out)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)c->impl;
	unsigned int dim = (unsigned int)c->dimension, n = impl->control_point_count;
	qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n * dim);
	qaws_bezier_desc d;
	qaws_status s;
	if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
	memcpy(cp, impl->control_points, sizeof(qaws_scalar) * n * dim);
	ex_bezier_restrict(cp, n, dim, t0, t1);
	memset(&d, 0, sizeof(d));
	d.dimension = c->dimension;
	d.degree = c->degree;
	d.control_points = cp;
	d.control_point_count = n;
	s = qaws_curve_create_bezier(&d, out);
	free(cp);
	return s;
}

static qaws_status ex_rational_bezier(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, qaws_curve** out)
{
	qaws_rational_bezier_impl const* impl = (qaws_rational_bezier_impl const*)c->impl;
	unsigned int dim = (unsigned int)c->dimension, n = impl->control_point_count, h = dim + 1, i, k;
	qaws_scalar* hp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n * (h + dim + 1));
	qaws_scalar *cp, *w;
	qaws_rational_bezier_desc d;
	qaws_status s;
	if (!hp) return QAWS_STATUS_ALLOCATION_FAILURE;
	cp = hp + n * h;
	w = cp + n * dim;
	for (i = 0; i < n; i++)
	{
		for (k = 0; k < dim; k++)
			hp[i * h + k] = impl->control_points[i * dim + k] * impl->weights[i];
		hp[i * h + dim] = impl->weights[i];
	}
	ex_bezier_restrict(hp, n, h, t0, t1);
	for (i = 0; i < n; i++)
	{
		w[i] = hp[i * h + dim];
		for (k = 0; k < dim; k++)
			cp[i * dim + k] = hp[i * h + k] / w[i];
	}
	memset(&d, 0, sizeof(d));
	d.dimension = c->dimension;
	d.degree = c->degree;
	d.control_points = cp;
	d.control_point_count = n;
	d.weights = w;
	d.weight_count = n;
	s = qaws_curve_create_rational_bezier(&d, out);
	free(hp);
	return s;
}

/* -------------------------------------------------------------------------- */
/*  B-spline and NURBS: knot insertion                                        */
/* -------------------------------------------------------------------------- */

typedef struct ex_spline
{
	qaws_scalar* knots;   /* knot_count */
	qaws_scalar* cp;      /* cp_count * stride */
	unsigned int knot_count, cp_count, stride, degree;
} ex_spline;

static unsigned int ex_multiplicity(ex_spline const* sp, qaws_scalar t)
{
	unsigned int i, m = 0;
	for (i = 0; i < sp->knot_count; i++)
		if (sp->knots[i] == t) m++;
	return m;
}

/* Inserts t once (Boehm). The arrays were allocated with room for it. */
static void ex_insert(ex_spline* sp, qaws_scalar t)
{
	unsigned int p = sp->degree, n = sp->cp_count, st = sp->stride, s = p, i, k;
	/* span s in [p, n - 1] with knots[s] <= t, the last such */
	while (s + 1 < n && sp->knots[s + 1] <= t) s++;
	for (k = 0; k < st; k++)
		sp->cp[n * st + k] = sp->cp[(n - 1) * st + k];
	for (i = n - 1; i > s; i--)
		for (k = 0; k < st; k++)
			sp->cp[i * st + k] = sp->cp[(i - 1) * st + k];
	for (i = s; i + p > s && i >= 1; i--)
	{
		qaws_scalar den = sp->knots[i + p] - sp->knots[i];
		qaws_scalar a = den > QAWS_ZERO ? (t - sp->knots[i]) / den : QAWS_ZERO;
		for (k = 0; k < st; k++)
			sp->cp[i * st + k] = a * sp->cp[i * st + k] + (QAWS_ONE - a) * sp->cp[(i - 1) * st + k];
	}
	for (i = sp->knot_count; i > s + 1; i--)
		sp->knots[i] = sp->knots[i - 1];
	sp->knots[s + 1] = t;
	sp->knot_count++;
	sp->cp_count++;
}

/* The insertion above walks i down from s while computing new points from
   old ones; the old P[i-1] it reads is still unmodified because the loop runs
   downward and each new Q[i] only needs P[i] and P[i-1]. */

/* Snaps t onto an existing knot when it is within rounding of it. */
static qaws_scalar ex_snap(ex_spline const* sp, qaws_scalar t, qaws_scalar tol)
{
	unsigned int i;
	for (i = 0; i < sp->knot_count; i++)
		if (QAWS_FABS(sp->knots[i] - t) <= tol)
			return sp->knots[i];
	return t;
}

/* Restricts sp to [t0, t1]; writes the clamped result into *res. */
static qaws_status ex_spline_restrict(ex_spline* sp, qaws_scalar t0, qaws_scalar t1, ex_spline* res)
{
	unsigned int p = sp->degree, st = sp->stride, L, F, i, k, nc;
	qaws_scalar tol = (sp->knots[sp->knot_count - 1] - sp->knots[0]) * (qaws_scalar)EX_REL;
	t0 = ex_snap(sp, t0, tol);
	t1 = ex_snap(sp, t1, tol);
	while (ex_multiplicity(sp, t0) < p) ex_insert(sp, t0);
	while (ex_multiplicity(sp, t1) < p) ex_insert(sp, t1);
	/* C(t0) = P[L - p] with L the last knot equal to t0; C(t1) = P[F - 1]
	   with F the first knot equal to t1. */
	L = 0;
	for (i = 0; i < sp->knot_count; i++)
		if (sp->knots[i] == t0) L = i;
	F = 0;
	for (i = sp->knot_count; i-- > 0;)
		if (sp->knots[i] == t1) F = i;
	if (L < p || F < 1 || F - 1 < L - p)
		return QAWS_STATUS_INTERNAL_ERROR;
	nc = F - L + p;
	res->degree = p;
	res->stride = st;
	res->cp_count = nc;
	res->knot_count = nc + p + 1;
	res->cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc * st);
	res->knots = (qaws_scalar*)malloc(sizeof(qaws_scalar) * res->knot_count);
	if (!res->cp || !res->knots)
	{
		free(res->cp); free(res->knots);
		res->cp = res->knots = NULL;
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nc; i++)
		for (k = 0; k < st; k++)
			res->cp[i * st + k] = sp->cp[(L - p + i) * st + k];
	for (i = 0; i <= p; i++)
	{
		res->knots[i] = t0;
		res->knots[res->knot_count - 1 - i] = t1;
	}
	for (i = L + 1; i < F; i++)
		res->knots[p + 1 + (i - L - 1)] = sp->knots[i];
	return QAWS_STATUS_OK;
}

static qaws_status ex_spline_kind(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, qaws_curve** out)
{
	int rational = c->kind == QAWS_CURVE_KIND_NURBS;
	unsigned int dim = (unsigned int)c->dimension, p = c->degree, i, k;
	qaws_scalar const* knots;
	qaws_scalar const* cps;
	qaws_scalar const* weights = NULL;
	unsigned int kc, n;
	ex_spline sp, res;
	qaws_status s;

	if (rational)
	{
		qaws_nurbs_impl const* impl = (qaws_nurbs_impl const*)c->impl;
		knots = impl->knots; kc = impl->knot_count;
		cps = impl->control_points; n = impl->control_point_count;
		weights = impl->weights;
	}
	else
	{
		qaws_bspline_impl const* impl = (qaws_bspline_impl const*)c->impl;
		knots = impl->knots; kc = impl->knot_count;
		cps = impl->control_points; n = impl->control_point_count;
	}

	sp.degree = p;
	sp.stride = dim + (rational ? 1u : 0u);
	sp.knot_count = kc;
	sp.cp_count = n;
	sp.knots = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (kc + 2 * p + 2));
	sp.cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (n + 2 * p + 2) * sp.stride);
	memset(&res, 0, sizeof(res));
	if (!sp.knots || !sp.cp)
	{
		free(sp.knots); free(sp.cp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(sp.knots, knots, sizeof(qaws_scalar) * kc);
	for (i = 0; i < n; i++)
	{
		qaws_scalar w = rational ? weights[i] : QAWS_ONE;
		for (k = 0; k < dim; k++)
			sp.cp[i * sp.stride + k] = cps[i * dim + k] * w;
		if (rational)
			sp.cp[i * sp.stride + dim] = w;
	}

	s = ex_spline_restrict(&sp, t0, t1, &res);
	free(sp.knots); free(sp.cp);
	if (s != QAWS_STATUS_OK)
		return s;

	if (rational)
	{
		qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * res.cp_count * (dim + 1));
		qaws_scalar* w;
		qaws_nurbs_desc d;
		if (!cp) { free(res.cp); free(res.knots); return QAWS_STATUS_ALLOCATION_FAILURE; }
		w = cp + res.cp_count * dim;
		for (i = 0; i < res.cp_count; i++)
		{
			w[i] = res.cp[i * res.stride + dim];
			for (k = 0; k < dim; k++)
				cp[i * dim + k] = res.cp[i * res.stride + k] / w[i];
		}
		memset(&d, 0, sizeof(d));
		d.dimension = c->dimension;
		d.degree = p;
		d.control_points = cp;
		d.control_point_count = res.cp_count;
		d.knots = res.knots;
		d.knot_count = res.knot_count;
		d.weights = w;
		d.weight_count = res.cp_count;
		s = qaws_curve_create_nurbs(&d, out);
		free(cp);
	}
	else
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = c->dimension;
		d.degree = p;
		d.control_points = res.cp;
		d.control_point_count = res.cp_count;
		d.knots = res.knots;
		d.knot_count = res.knot_count;
		s = qaws_curve_create_bspline(&d, out);
	}
	free(res.cp); free(res.knots);
	return s;
}

/* -------------------------------------------------------------------------- */
/*  Piecewise cubic kinds: Bezier spans from four evaluations                  */
/* -------------------------------------------------------------------------- */

typedef struct ex_cubics
{
	qaws_scalar* cp;      /* 3 * count + 1 points, 3 scalars each */
	qaws_scalar* breaks;  /* count + 1 parameters */
	unsigned int count, capacity;
	qaws_scalar tol;      /* fit tolerance, 0 = exact kinds */
	int status;
} ex_cubics;

static int ex_cubics_reserve(ex_cubics* q, unsigned int n)
{
	qaws_scalar *cp, *br;
	unsigned int cap = q->capacity ? q->capacity : 8;
	if (n <= q->capacity) return 1;
	while (cap < n) cap *= 2;
	cp = (qaws_scalar*)realloc(q->cp, sizeof(qaws_scalar) * 3 * (3 * cap + 1));
	if (!cp) return 0;
	q->cp = cp;
	br = (qaws_scalar*)realloc(q->breaks, sizeof(qaws_scalar) * (cap + 1));
	if (!br) return 0;
	q->breaks = br;
	q->capacity = cap;
	return 1;
}

static void ex_cubic_at(qaws_scalar const* b, qaws_scalar u, qaws_scalar* out)
{
	qaws_scalar v = QAWS_ONE - u;
	qaws_scalar w0 = v * v * v, w1 = 3 * u * v * v, w2 = 3 * u * u * v, w3 = u * u * u;
	unsigned int k;
	for (k = 0; k < 3; k++)
		out[k] = w0 * b[k] + w1 * b[3 + k] + w2 * b[6 + k] + w3 * b[9 + k];
}

/* Appends the cubic of span `span` on local [ua, ub] (global [ta, tb]). */
static void ex_cubic_piece(qaws_curve const* c, ex_cubics* q, unsigned int span,
	qaws_scalar ua, qaws_scalar ub, qaws_scalar ta, qaws_scalar tb, unsigned int depth)
{
	qaws_scalar p[4][3], b[12], h = ub - ua;
	unsigned int k, j;
	if (q->status != QAWS_STATUS_OK) return;
	for (j = 0; j < 4; j++)
	{
		qaws_status s = ex_span_point(c, span, j == 3 ? ub : ua + h * (qaws_scalar)j / 3, p[j]);
		if (s != QAWS_STATUS_OK) { q->status = s; return; }
	}
	/* interpolation at 0, 1/3, 2/3, 1:
	   27 P(1/3) = 8 b0 + 12 b1 + 6 b2 + b3, 27 P(2/3) = b0 + 6 b1 + 12 b2 + 8 b3 */
	for (k = 0; k < 3; k++)
	{
		qaws_scalar A = 27 * p[1][k] - 8 * p[0][k] - p[3][k];
		qaws_scalar B = 27 * p[2][k] - p[0][k] - 8 * p[3][k];
		b[k] = p[0][k];
		b[3 + k] = (2 * A - B) / 18;
		b[6 + k] = (2 * B - A) / 18;
		b[9 + k] = p[3][k];
	}
	if (q->tol > QAWS_ZERO && depth < EX_FIT_DEPTH)
	{
		static qaws_scalar const probe[3] = { (qaws_scalar)(1.0 / 6.0), (qaws_scalar)0.5, (qaws_scalar)(5.0 / 6.0) };
		qaws_scalar err = QAWS_ZERO, x[3], y[3];
		for (j = 0; j < 3; j++)
		{
			qaws_status s = ex_span_point(c, span, ua + h * probe[j], x);
			if (s != QAWS_STATUS_OK) { q->status = s; return; }
			ex_cubic_at(b, probe[j], y);
			for (k = 0; k < 3; k++)
				if (QAWS_FABS(x[k] - y[k]) > err) err = QAWS_FABS(x[k] - y[k]);
		}
		if (err > q->tol)
		{
			qaws_scalar um = (ua + ub) * (qaws_scalar)0.5, tm = (ta + tb) * (qaws_scalar)0.5;
			ex_cubic_piece(c, q, span, ua, um, ta, tm, depth + 1);
			ex_cubic_piece(c, q, span, um, ub, tm, tb, depth + 1);
			return;
		}
	}
	if (!ex_cubics_reserve(q, q->count + 1)) { q->status = QAWS_STATUS_ALLOCATION_FAILURE; return; }
	if (q->count == 0)
	{
		memcpy(q->cp, b, sizeof(qaws_scalar) * 3);
		q->breaks[0] = ta;
	}
	/* the shared end point stays the previous piece's */
	memcpy(&q->cp[(3 * q->count + 1) * 3], &b[3], sizeof(qaws_scalar) * 9);
	q->count++;
	q->breaks[q->count] = tb;
}

static qaws_status ex_cubics_curve(qaws_curve const* c, ex_cubics* q, qaws_curve** out)
{
	unsigned int dim = (unsigned int)c->dimension, n = 3 * q->count + 1, i, k;
	qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (n * dim + n + 4));
	qaws_scalar* knots;
	qaws_status s;
	if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
	knots = cp + n * dim;
	for (i = 0; i < n; i++)
		for (k = 0; k < dim; k++)
			cp[i * dim + k] = q->cp[i * 3 + k];
	if (q->count == 1)
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = c->dimension;
		d.degree = 3;
		d.control_points = cp;
		d.control_point_count = 4;
		s = qaws_curve_create_bezier(&d, out);
	}
	else
	{
		/* piecewise Bezier as a cubic B-spline: interior knots of multiplicity 3 */
		qaws_bspline_desc d;
		unsigned int kc = 0;
		for (i = 0; i < 4; i++) knots[kc++] = q->breaks[0];
		for (i = 1; i < q->count; i++)
			for (k = 0; k < 3; k++) knots[kc++] = q->breaks[i];
		for (i = 0; i < 4; i++) knots[kc++] = q->breaks[q->count];
		memset(&d, 0, sizeof(d));
		d.dimension = c->dimension;
		d.degree = 3;
		d.control_points = cp;
		d.control_point_count = n;
		d.knots = knots;
		d.knot_count = kc;
		s = qaws_curve_create_bspline(&d, out);
	}
	free(cp);
	return s;
}

static qaws_scalar ex_extent(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1)
{
	qaws_scalar lo[3], hi[3], p[3], e = QAWS_ZERO;
	unsigned int i, k;
	for (i = 0; i <= 32; i++)
	{
		ex_point(c, t0 + (t1 - t0) * (qaws_scalar)i / 32, p);
		for (k = 0; k < 3; k++)
		{
			if (i == 0 || p[k] < lo[k]) lo[k] = p[k];
			if (i == 0 || p[k] > hi[k]) hi[k] = p[k];
		}
	}
	for (k = 0; k < 3; k++)
		e += (hi[k] - lo[k]) * (hi[k] - lo[k]);
	return (qaws_scalar)sqrt((double)e);
}

static qaws_status ex_cubic_kind(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, int fit, qaws_curve** out)
{
	ex_cubics q;
	qaws_scalar u0, u1;
	unsigned int s0 = qaws_internal_find_span(c, t0, &u0);
	unsigned int s1 = qaws_internal_find_span(c, t1, &u1);
	unsigned int s;
	qaws_status st;

	if (s1 > s0 && u1 <= QAWS_ZERO) { s1--; u1 = QAWS_ONE; }
	memset(&q, 0, sizeof(q));
	q.status = QAWS_STATUS_OK;
	if (fit)
	{
		q.tol = ex_extent(c, t0, t1) * (qaws_scalar)EX_FIT_REL;
		if (q.tol <= QAWS_ZERO) q.tol = (qaws_scalar)EX_FIT_REL;
	}
	for (s = s0; s <= s1; s++)
	{
		qaws_scalar ua = s == s0 ? u0 : QAWS_ZERO, ub = s == s1 ? u1 : QAWS_ONE;
		qaws_scalar ta = s == s0 ? t0 : c->span_boundaries[s];
		qaws_scalar tb = s == s1 ? t1 : c->span_boundaries[s + 1];
		if (ub > ua)
			ex_cubic_piece(c, &q, s, ua, ub, ta, tb, 0);
	}
	st = q.status;
	if (st == QAWS_STATUS_OK)
		st = q.count ? ex_cubics_curve(c, &q, out) : QAWS_STATUS_INVALID_ARGUMENT;
	free(q.cp); free(q.breaks);
	return st;
}

/* -------------------------------------------------------------------------- */
/*  Polynomial, arc, clothoid                                                 */
/* -------------------------------------------------------------------------- */

static qaws_status ex_polynomial(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, qaws_curve** out)
{
	qaws_polynomial_impl const* impl = (qaws_polynomial_impl const*)c->impl;
	qaws_polynomial_desc d;
	memset(&d, 0, sizeof(d));
	d.dimension = c->dimension;
	d.degree = c->degree;
	d.coefficients = impl->coefficients;
	d.coefficient_count = impl->coefficient_count;
	d.t_min = t0;
	d.t_max = t1;
	return qaws_curve_create_polynomial(&d, out);
}

static qaws_status ex_arc(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, qaws_curve** out)
{
	qaws_arc_impl const* impl = (qaws_arc_impl const*)c->impl;
	qaws_scalar u0, u1;
	unsigned int s0 = qaws_internal_find_span(c, t0, &u0);
	unsigned int s1 = qaws_internal_find_span(c, t1, &u1);
	unsigned int s, n = 0;
	qaws_arc_segment* seg;
	qaws_arc_desc d;
	qaws_status st;

	if (s1 > s0 && u1 <= QAWS_ZERO) { s1--; u1 = QAWS_ONE; }
	seg = (qaws_arc_segment*)malloc(sizeof(qaws_arc_segment) * (s1 - s0 + 1));
	if (!seg) return QAWS_STATUS_ALLOCATION_FAILURE;
	for (s = s0; s <= s1; s++)
	{
		qaws_arc_segment g = impl->segments[s];
		qaws_scalar sweep = g.angle_end - g.angle_start;
		qaws_scalar ua = s == s0 ? u0 : QAWS_ZERO, ub = s == s1 ? u1 : QAWS_ONE;
		if (ub <= ua) continue;
		g.angle_end = g.angle_start + ub * sweep;
		g.angle_start = g.angle_start + ua * sweep;
		seg[n++] = g;
	}
	memset(&d, 0, sizeof(d));
	d.dimension = c->dimension;
	d.segments = seg;
	d.segment_count = n;
	st = n ? qaws_curve_create_arc(&d, out) : QAWS_STATUS_INVALID_ARGUMENT;
	free(seg);
	return st;
}

static qaws_status ex_clothoid(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, qaws_curve** out)
{
	qaws_clothoid_impl const* impl = (qaws_clothoid_impl const*)c->impl;
	qaws_clothoid_desc d;
	qaws_scalar p[3];
	qaws_status st = ex_point(c, t0, p);
	if (st != QAWS_STATUS_OK) return st;
	memset(&d, 0, sizeof(d));
	d.origin_x = p[0];
	d.origin_y = p[1];
	d.start_angle = impl->start_angle + impl->kappa_0 * t0 + impl->rate * t0 * t0 * (qaws_scalar)0.5;
	d.start_curvature = impl->kappa_0 + impl->rate * t0;
	d.end_curvature = impl->kappa_0 + impl->rate * t1;
	d.length = t1 - t0;
	return qaws_curve_create_clothoid(&d, out);
}

/* -------------------------------------------------------------------------- */
/*  Composite                                                                 */
/* -------------------------------------------------------------------------- */

static qaws_status ex_composite(qaws_curve const* c, qaws_scalar t0, qaws_scalar t1, qaws_curve** out)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)c->impl;
	qaws_scalar u0, u1;
	unsigned int s0 = qaws_internal_find_span(c, t0, &u0);
	unsigned int s1 = qaws_internal_find_span(c, t1, &u1);
	unsigned int s, n = 0, i;
	qaws_curve** parts;
	qaws_composite_desc d;
	qaws_status st = QAWS_STATUS_OK;

	if (s1 > s0 && u1 <= QAWS_ZERO) { s1--; u1 = QAWS_ONE; }
	parts = (qaws_curve**)malloc(sizeof(qaws_curve*) * (s1 - s0 + 1));
	if (!parts) return QAWS_STATUS_ALLOCATION_FAILURE;
	for (s = s0; s <= s1 && st == QAWS_STATUS_OK; s++)
	{
		qaws_curve const* g = impl->segments[s];
		qaws_scalar a = g->parameter_range.min_value, b = g->parameter_range.max_value;
		qaws_scalar ua = s == s0 ? u0 : QAWS_ZERO, ub = s == s1 ? u1 : QAWS_ONE;
		if (ub <= ua) continue;
		st = qaws_curve_extract(g, a + (b - a) * ua, ub >= QAWS_ONE ? b : a + (b - a) * ub, &parts[n]);
		if (st == QAWS_STATUS_OK) n++;
	}
	if (st == QAWS_STATUS_OK && n == 0)
		st = QAWS_STATUS_INVALID_ARGUMENT;
	if (st != QAWS_STATUS_OK)
	{
		for (i = 0; i < n; i++) qaws_curve_destroy(parts[i]);
		free(parts);
		return st;
	}
	if (n == 1)
	{
		*out = parts[0];
		free(parts);
		return QAWS_STATUS_OK;
	}
	memset(&d, 0, sizeof(d));
	d.dimension = c->dimension;
	d.segments = parts;
	d.segment_count = n;
	st = qaws_curve_create_composite(&d, out);
	if (st != QAWS_STATUS_OK)
		for (i = 0; i < n; i++) qaws_curve_destroy(parts[i]);
	free(parts);
	return st;
}

/* -------------------------------------------------------------------------- */
/*  Public API                                                                */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_extract(
	qaws_curve const* curve,
	qaws_scalar t0,
	qaws_scalar t1,
	qaws_curve** out_curve)
{
	qaws_scalar a, b, tol;
	qaws_status st;

	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;

	a = curve->parameter_range.min_value;
	b = curve->parameter_range.max_value;
	tol = (b - a) * (qaws_scalar)EX_REL;
	if (t0 < a - tol || t0 > b + tol || t1 < a - tol || t1 > b + tol)
		return QAWS_STATUS_OUT_OF_RANGE;
	if (t0 < a) t0 = a;
	if (t0 > b) t0 = b;
	if (t1 < a) t1 = a;
	if (t1 > b) t1 = b;
	if (QAWS_FABS(t1 - t0) <= tol)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (t0 > t1)
	{
		qaws_curve* fwd = NULL;
		st = qaws_curve_extract(curve, t1, t0, &fwd);
		if (st != QAWS_STATUS_OK)
			return st;
		st = qaws_curve_reverse(fwd, out_curve);
		qaws_curve_destroy(fwd);
		return st;
	}

	switch (curve->kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
		return ex_bezier(curve, (t0 - a) / (b - a), (t1 - a) / (b - a), out_curve);
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
		return ex_rational_bezier(curve, (t0 - a) / (b - a), (t1 - a) / (b - a), out_curve);
	case QAWS_CURVE_KIND_BSPLINE:
	case QAWS_CURVE_KIND_NURBS:
		return ex_spline_kind(curve, t0, t1, out_curve);
	case QAWS_CURVE_KIND_HERMITE:
	case QAWS_CURVE_KIND_CATMULL_ROM:
	case QAWS_CURVE_KIND_TRAJECTORY:
	case QAWS_CURVE_KIND_SUBDIVISION:
		return ex_cubic_kind(curve, t0, t1, 0, out_curve);
	case QAWS_CURVE_KIND_YUKSEL:
		return ex_cubic_kind(curve, t0, t1, 1, out_curve);
	case QAWS_CURVE_KIND_POLYNOMIAL:
		return ex_polynomial(curve, t0, t1, out_curve);
	case QAWS_CURVE_KIND_ARC:
		return ex_arc(curve, t0, t1, out_curve);
	case QAWS_CURVE_KIND_CLOTHOID:
		return ex_clothoid(curve, t0, t1, out_curve);
	case QAWS_CURVE_KIND_COMPOSITE:
		return ex_composite(curve, t0, t1, out_curve);
	case QAWS_CURVE_KIND_REPARAMETERIZED:
	{
		qaws_curve const* src0 = NULL;
		qaws_curve const* src1 = NULL;
		qaws_scalar s0, s1;
		if (!qaws_internal_reparam_source(curve, t0, &src0, &s0) ||
		    !qaws_internal_reparam_source(curve, t1, &src1, &s1))
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
		return qaws_curve_extract(src0, s0, s1, out_curve);
	}
	default:
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	}
}
