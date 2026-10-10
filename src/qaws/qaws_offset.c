/*
 * Offsetting paths of curves: raw offset loops (pieces, joins, caps), then
 * one positive union through qaws_clip.
 */

#include "qaws_offset.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "qaws_bezier.h"
#include "qaws_bspline.h"
#include "qaws_arc.h"
#include "qaws_composite.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_flatten.h"
#include "internal/qaws_internal_span.h"
#include "internal/qaws_internal_offset.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define OS_PI 3.14159265358979323846
#define OS_FIT_DEPTH 20

/* ======================================================================== */
/*  Elements: the pieces of a path, in travel order                         */
/* ======================================================================== */

typedef struct os_elem
{
	qaws_curve const* c;
	double t0, t1;          /* travelled from t0 to t1 (t1 < t0: backwards) */
	int line, arc;          /* a straight line; one arc segment */
	double p0[2], p1[2];    /* end points */
	double d0[2], d1[2];    /* unit tangents at the ends, in travel direction */
} os_elem;

typedef struct os_elems
{
	os_elem* e;
	unsigned int n, cap;
} os_elems;

typedef struct os_ctx
{
	qaws_offset_desc const* d;
	double tol;
	double miter_lim;              /* 2 / ml^2 */
	unsigned int group, path;
	qaws_curve** curves;           /* every curve made (owned) */
	unsigned int ncurves, capcurves;
	unsigned int* loop_first;      /* loop i: curves [loop_first[i], loop_first[i + 1]) */
	unsigned int nloops, caploops;
} os_ctx;

static qaws_status os_eval(qaws_curve const* c, double t, double* p, double* d1, double* d2)
{
	qaws_eval_result_2d r;
	unsigned int fl = QAWS_EVAL_FLAG_POSITION | (d1 ? QAWS_EVAL_FLAG_D1 : 0) | (d2 ? QAWS_EVAL_FLAG_D2 : 0);
	qaws_status s = qaws_curve_evaluate_2d(c, (qaws_scalar)t, fl, &r);
	p[0] = r.position.x; p[1] = r.position.y;
	if (d1) { d1[0] = r.d1.x; d1[1] = r.d1.y; }
	if (d2) { d2[0] = r.d2.x; d2[1] = r.d2.y; }
	return s;
}

/* position and derivatives at t; at a span boundary, side < 0 takes the span
   ending there (the end of a piece), side > 0 the one starting there */
static qaws_status os_eval_side(qaws_curve const* c, double t, int side, double* p, double* d1, double* d2)
{
	qaws_eval_result_2d r;
	unsigned int fl = QAWS_EVAL_FLAG_POSITION | (d1 ? QAWS_EVAL_FLAG_D1 : 0) | (d2 ? QAWS_EVAL_FLAG_D2 : 0);
	qaws_scalar local = 0;
	unsigned int s = qaws_internal_find_span(c, (qaws_scalar)t, &local);
	qaws_status st;
	if (side < 0 && s > 0 && local <= 0 && (qaws_scalar)t > c->parameter_range.min_value)
	{
		s--;
		local = 1;
	}
	st = c->vtable->eval_span_2d(c, s, local, fl, &r);
	p[0] = r.position.x; p[1] = r.position.y;
	if (d1) { d1[0] = r.d1.x; d1[1] = r.d1.y; }
	if (d2) { d2[0] = r.d2.x; d2[1] = r.d2.y; }
	return st;
}

static void os_unit(double* v)
{
	double n = hypot(v[0], v[1]);
	if (n > 0) { v[0] /= n; v[1] /= n; }
}

static int os_push_elem(os_elems* q, os_elem const* e)
{
	if (q->n == q->cap)
	{
		unsigned int nc = q->cap ? 2 * q->cap : 32;
		os_elem* g = (os_elem*)realloc(q->e, sizeof(os_elem) * nc);
		if (!g) return 0;
		q->e = g;
		q->cap = nc;
	}
	q->e[q->n++] = *e;
	return 1;
}

/* a piece from t0 to t1 of c, with its ends and tangents */
static qaws_status os_elem_make(qaws_curve const* c, double t0, double t1, int line, int arc, os_elems* q)
{
	os_elem e;
	double dd[2], step = (t1 - t0) * 1e-7;
	memset(&e, 0, sizeof(e));
	e.c = c; e.t0 = t0; e.t1 = t1; e.line = line; e.arc = arc;
	os_eval_side(c, t0, t1 >= t0 ? 1 : -1, e.p0, e.d0, NULL);
	os_eval_side(c, t1, t1 >= t0 ? -1 : 1, e.p1, e.d1, NULL);
	if (t1 < t0)
	{
		e.d0[0] = -e.d0[0]; e.d0[1] = -e.d0[1];
		e.d1[0] = -e.d1[0]; e.d1[1] = -e.d1[1];
	}
	/* a vanishing derivative at an end (a cusp): the direction of a nearby point */
	if (!(hypot(e.d0[0], e.d0[1]) > 0))
	{
		os_eval(c, t0 + step, dd, NULL, NULL);
		e.d0[0] = dd[0] - e.p0[0]; e.d0[1] = dd[1] - e.p0[1];
	}
	if (!(hypot(e.d1[0], e.d1[1]) > 0))
	{
		os_eval(c, t1 - step, dd, NULL, NULL);
		e.d1[0] = e.p1[0] - dd[0]; e.d1[1] = e.p1[1] - dd[1];
	}
	if (line)
	{
		/* a line's direction is its chord (an evaluator at a knot may hand
		   back the next span's derivative) */
		if (e.p0[0] == e.p1[0] && e.p0[1] == e.p1[1])
			return QAWS_STATUS_OK;
		e.d0[0] = e.d1[0] = e.p1[0] - e.p0[0];
		e.d0[1] = e.d1[1] = e.p1[1] - e.p0[1];
	}
	os_unit(e.d0);
	os_unit(e.d1);
	return os_push_elem(q, &e) ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
}

/* a curve's elements: polylines by segment, composites by segment, arcs by
   arc segment, anything else whole */
static qaws_status os_curve_elems(qaws_curve const* c, os_elems* q)
{
	qaws_status s = QAWS_STATUS_OK;
	unsigned int k;
	if (c->kind == QAWS_CURVE_KIND_BSPLINE && c->degree == 1)
	{
		qaws_bspline_impl const* impl = (qaws_bspline_impl const*)c->impl;
		for (k = 0; k + 1 < impl->control_point_count && s == QAWS_STATUS_OK; k++)
			if (impl->knots[k + 2] > impl->knots[k + 1])
				s = os_elem_make(c, impl->knots[k + 1], impl->knots[k + 2], 1, 0, q);
		return s;
	}
	if (c->kind == QAWS_CURVE_KIND_BEZIER && c->degree == 1)
		return os_elem_make(c, c->parameter_range.min_value, c->parameter_range.max_value, 1, 0, q);
	if (c->kind == QAWS_CURVE_KIND_COMPOSITE)
	{
		qaws_composite_impl const* impl = (qaws_composite_impl const*)c->impl;
		for (k = 0; k < impl->segment_count && s == QAWS_STATUS_OK; k++)
			s = os_curve_elems(impl->segments[k], q);
		return s;
	}
	if (c->kind == QAWS_CURVE_KIND_ARC)
	{
		for (k = 0; k < c->span_count && s == QAWS_STATUS_OK; k++)
			if (c->span_boundaries[k + 1] > c->span_boundaries[k])
				s = os_elem_make(c, c->span_boundaries[k], c->span_boundaries[k + 1], 0, 1, q);
		return s;
	}
	return os_elem_make(c, c->parameter_range.min_value, c->parameter_range.max_value, 0, 0, q);
}

/* ======================================================================== */
/*  Output curves and loops                                                 */
/* ======================================================================== */

static qaws_status os_keep(os_ctx* x, qaws_curve* c)
{
	if (x->ncurves == x->capcurves)
	{
		unsigned int nc = x->capcurves ? 2 * x->capcurves : 256;
		qaws_curve** g = (qaws_curve**)realloc(x->curves, sizeof(qaws_curve*) * nc);
		if (!g) { qaws_curve_destroy(c); return QAWS_STATUS_ALLOCATION_FAILURE; }
		x->curves = g;
		x->capcurves = nc;
	}
	x->curves[x->ncurves++] = c;
	return QAWS_STATUS_OK;
}

static qaws_status os_line(os_ctx* x, double const* a, double const* b)
{
	qaws_scalar cp[4];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	qaws_status s;
	if (a[0] == b[0] && a[1] == b[1])
		return QAWS_STATUS_OK;
	cp[0] = (qaws_scalar)a[0]; cp[1] = (qaws_scalar)a[1];
	cp[2] = (qaws_scalar)b[0]; cp[3] = (qaws_scalar)b[1];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
	s = qaws_curve_create_bezier(&d, &c);
	return s == QAWS_STATUS_OK ? os_keep(x, c) : s;
}

/* an arc round center from angle a0 sweeping by sw */
static qaws_status os_arc(os_ctx* x, double const* center, double r, double a0, double sw)
{
	qaws_arc_segment g;
	qaws_arc_desc d;
	qaws_curve* c = NULL;
	qaws_status s;
	if (!(r > 0) || fabs(sw) * r <= 1e-12 * (r + fabs(center[0]) + fabs(center[1])))
		return QAWS_STATUS_OK;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)center[0]; g.center[1] = (qaws_scalar)center[1];
	g.radius = (qaws_scalar)r;
	g.angle_start = (qaws_scalar)a0;
	g.angle_end = (qaws_scalar)(a0 + sw);
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
	s = qaws_curve_create_arc(&d, &c);
	return s == QAWS_STATUS_OK ? os_keep(x, c) : s;
}

static qaws_status os_loop_end(os_ctx* x, unsigned int first)
{
	if (x->nloops + 2 > x->caploops)
	{
		unsigned int nc = x->caploops ? 2 * x->caploops : 64;
		unsigned int* g = (unsigned int*)realloc(x->loop_first, sizeof(unsigned int) * nc);
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->loop_first = g;
		x->caploops = nc;
	}
	if (x->ncurves == first)
		return QAWS_STATUS_OK;
	x->loop_first[x->nloops] = first;
	x->nloops++;
	x->loop_first[x->nloops] = x->ncurves;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Offsets of elements                                                     */
/* ======================================================================== */

/* delta at a point with the offset side's unit normal */
static double os_delta(os_ctx const* x, double delta, double const* p, double const* n)
{
	qaws_vec2 q, nn;
	if (!x->d->delta_fn)
		return delta;
	q.x = (qaws_scalar)p[0]; q.y = (qaws_scalar)p[1];
	nn.x = (qaws_scalar)n[0]; nn.y = (qaws_scalar)n[1];
	return (double)x->d->delta_fn(x->d->delta_user, x->group, x->path, q, nn) * (delta < 0 ? -1.0 : 1.0);
}

/* right normal of a unit tangent */
static void os_right(double const* t, double* n)
{
	n[0] = t[1];
	n[1] = -t[0];
}

/* the offset point at parameter t in span `span`: C + delta N, N the right
   normal of the travel direction (only the direction of C' is used, so the
   derivative's scale does not matter) */
static void os_offset_at(os_ctx const* x, os_elem const* e, double delta, double t, unsigned int span, double* o)
{
	qaws_curve const* c = e->c;
	qaws_eval_result_2d r;
	double a = c->span_boundaries[span], b = c->span_boundaries[span + 1], local = b > a ? (t - a) / (b - a) : 0, d[2], n[2], dl;
	if (local < 0) local = 0;
	if (local > 1) local = 1;
	c->vtable->eval_span_2d(c, span, (qaws_scalar)local, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
	d[0] = r.d1.x; d[1] = r.d1.y;
	if (e->t1 < e->t0) { d[0] = -d[0]; d[1] = -d[1]; }
	if (!(hypot(d[0], d[1]) > 0))
	{
		/* a cusp: the chord to a nearby point */
		qaws_eval_result_2d q;
		double l2 = local + (local < 0.5 ? 1e-7 : -1e-7);
		c->vtable->eval_span_2d(c, span, (qaws_scalar)l2, QAWS_EVAL_FLAG_POSITION, &q);
		d[0] = (q.position.x - r.position.x) * (l2 > local ? 1 : -1);
		d[1] = (q.position.y - r.position.y) * (l2 > local ? 1 : -1);
		if (e->t1 < e->t0) { d[0] = -d[0]; d[1] = -d[1]; }
	}
	os_unit(d);
	n[0] = d[1]; n[1] = -d[0];
	{
		double p[2] = { r.position.x, r.position.y };
		dl = os_delta(x, delta, p, n);
	}
	o[0] = r.position.x + dl * n[0];
	o[1] = r.position.y + dl * n[1];
}

typedef struct os_fit
{
	qaws_scalar* cp;      /* 3 count + 1 points */
	unsigned int count, cap;
	int status;
	int pend;             /* a link from the evolute back to the offset, held at a span end */
	double pe[2], po[2];
} os_fit;

static void os_fit_push(os_fit* f, double const* b);
static void os_link_flush(os_fit* f);
static qaws_status os_join(os_ctx* x, qaws_join_type jt, double const* v, double const* tin, double const* tout,
	double const* P, double const* Q, double delta);

/* the cubic through the offset at ta, ta + h/3, ta + 2h/3, tb, split while a
   probe between misses the true offset by more than the tolerance */
static void os_fit_piece(os_ctx const* x, os_elem const* e, double delta, double ta, double tb, unsigned int span,
	os_fit* f, unsigned int depth)
{
	double o[4][2], b[8], err = 0.0;
	unsigned int j, k;
	if (f->status != QAWS_STATUS_OK) return;
	for (j = 0; j < 4; j++)
		os_offset_at(x, e, delta, j == 3 ? tb : ta + (tb - ta) * j / 3.0, span, o[j]);
	/* 27 O(1/3) = 8 b0 + 12 b1 + 6 b2 + b3, 27 O(2/3) = b0 + 6 b1 + 12 b2 + 8 b3 */
	for (k = 0; k < 2; k++)
	{
		double A = 27 * o[1][k] - 8 * o[0][k] - o[3][k], B = 27 * o[2][k] - o[0][k] - 8 * o[3][k];
		b[k] = o[0][k];
		b[2 + k] = (2 * A - B) / 18;
		b[4 + k] = (2 * B - A) / 18;
		b[6 + k] = o[3][k];
	}
	if (depth < OS_FIT_DEPTH)
	{
		static double const probe[3] = { 1.0 / 6.0, 0.5, 5.0 / 6.0 };
		for (j = 0; j < 3; j++)
		{
			double u = probe[j], v = 1 - u, q[2], y[2];
			os_offset_at(x, e, delta, ta + (tb - ta) * u, span, q);
			for (k = 0; k < 2; k++)
				y[k] = v * v * v * b[k] + 3 * u * v * v * b[2 + k] + 3 * u * u * v * b[4 + k] + u * u * u * b[6 + k];
			if (hypot(q[0] - y[0], q[1] - y[1]) > err) err = hypot(q[0] - y[0], q[1] - y[1]);
		}
		if (err > x->tol)
		{
			double tm = 0.5 * (ta + tb);
			os_fit_piece(x, e, delta, ta, tm, span, f, depth + 1);
			os_fit_piece(x, e, delta, tm, tb, span, f, depth + 1);
			return;
		}
	}
	os_fit_push(f, b);
}

/* the chain of cubics as one cubic B-spline (knots 0..count, multiplicity
   3), kept; the chain is emptied */
static qaws_status os_fit_emit(os_ctx* x, os_fit* f)
{
	qaws_status st;
	os_link_flush(f);
	st = f->status;
	if (st == QAWS_STATUS_OK && f->count)
	{
		unsigned int np = 3 * f->count + 1, kc = 0, i;
		qaws_scalar* kn = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (np + 4));
		qaws_bspline_desc d;
		qaws_curve* c = NULL;
		if (!kn) st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			for (i = 0; i < 4; i++) kn[kc++] = 0;
			for (i = 1; i < f->count; i++) { kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; }
			for (i = 0; i < 4; i++) kn[kc++] = (qaws_scalar)f->count;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = f->cp; d.control_point_count = np;
			d.knots = kn; d.knot_count = kc;
			st = qaws_curve_create_bspline(&d, &c);
			if (st == QAWS_STATUS_OK)
				st = os_keep(x, c);
			free(kn);
		}
	}
	free(f->cp);
	memset(f, 0, sizeof(*f));
	f->status = QAWS_STATUS_OK;
	return st;
}

/* the offset of a curve element under a variable delta: cubics through the
   offset points, span by span */
static qaws_status os_fit_elem(os_ctx* x, os_elem const* e, double delta)
{
	os_fit f;
	double lo = e->t0 < e->t1 ? e->t0 : e->t1, hi = e->t0 < e->t1 ? e->t1 : e->t0;
	unsigned int n;
	memset(&f, 0, sizeof(f));
	f.status = QAWS_STATUS_OK;
	/* span by span, in travel order */
	for (n = 0; n < e->c->span_count; n++)
	{
		unsigned int k = e->t1 >= e->t0 ? n : e->c->span_count - 1 - n;
		double a = e->c->span_boundaries[k], b = e->c->span_boundaries[k + 1];
		if (a < lo) a = lo;
		if (b > hi) b = hi;
		if (!(b > a)) continue;
		if (e->t1 >= e->t0) os_fit_piece(x, e, delta, a, b, k, &f, 0);
		else os_fit_piece(x, e, delta, b, a, k, &f, 0);
	}
	return os_fit_emit(x, &f);
}

/* ======================================================================== */
/*  Offsets of curves at a constant delta                                   */
/* ======================================================================== */

/*
 * After "Fast GPU stroke expansion" (Levien and Uguray, HPG 2024), with
 * curves out instead of line segments.
 *
 * Each smooth span is cut, in travel order, wherever it stops looking like
 * an Euler spiral: the paper's estimate on the end points and derivatives
 * (as it converts cubics to Euler spirals), plus the gap at the middle to
 * the Hermite cubic of those ends, for kinds that are not cubics. Between
 * two cuts the curvature k is then close to linear.
 *
 * The offset O = C + delta n (n the right normal) has O' = |C'| g T with
 * g = 1 + delta k: it runs forward where g > 0, backward where g < 0, and
 * has a cusp at each root of g. With k linear between cuts, each root is
 * bracketed by two cuts of opposite sign, then polished on the true curve.
 *
 * Forward runs become cubics fitted to the true offset: the ends and their
 * tangents (T) held, so neighbouring cubics meet G1. The two arm lengths
 * start from least squares on samples or from the Hermite arms (the true
 * speed |C'| |g| per unit of the curve's parameter, which vanishes at a
 * cusp: the cubic's end is then a cusp too, where free arms along unit
 * tangents would need many short cubics), whichever is closer, then
 * Gauss-Newton on the samples' distances. A cubic is split while it misses
 * the offset by more than the tolerance.
 *
 * A backward run is replaced by the evolute E = C - n / k, the centres of
 * curvature, which the offset meets at its cusps (g = 0 there). The chain
 * runs O up to the cusp, along E, then on along O from the next cusp: the
 * smooth form of Clipper's link through an inner corner, and that link in
 * the limit of a sharp corner. A backward run alone already winds with the
 * sign of delta, so the positive union of one side would be right without
 * the evolute; the evolute is smooth where the backward run ends in two
 * cusps, so it takes fewer cubics (a parabola toward its focus at 1e-9:
 * 30 cubics, 54 with the backward run). E' = n k' / k^2 runs along the
 * normal and has its own cusps where k peaks; its cubics are cut there.
 */

#define OS_ES_DEPTH 22          /* finest Euler cut: 2^-22 of a span */
#define OS_FIT_SAMPLES 10       /* least-squares samples per cubic */

typedef struct os_run
{
	qaws_curve const* c;
	unsigned int span;
	double a, b;            /* span-local parameters, travelled from a to b */
	double scale;           /* the evaluator's derivatives per unit of the local parameter */
	double delta;
	double* cusps;          /* optional: the offset's cusps, as curve parameters */
	unsigned int cusp_cap, *cusp_count;
} os_run;

typedef struct os_geo
{
	double p[2], t[2], n[2];    /* point, unit travel tangent, right normal */
	double k;                   /* curvature, positive turning left along travel */
	double q[2];                /* derivative along travel per unit of tau */
} os_geo;

static void os_geo_raw(os_run const* r, double tau, os_geo* g)
{
	qaws_eval_result_2d e;
	double s = r->b - r->a, local = r->a + tau * s, q2[2], v;
	if (local < 0) local = 0;
	if (local > 1) local = 1;
	r->c->vtable->eval_span_2d(r->c, r->span, (qaws_scalar)local,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &e);
	g->p[0] = e.position.x; g->p[1] = e.position.y;
	g->q[0] = e.d1.x * r->scale * s; g->q[1] = e.d1.y * r->scale * s;
	q2[0] = e.d2.x * r->scale * r->scale * s * s; q2[1] = e.d2.y * r->scale * r->scale * s * s;
	v = hypot(g->q[0], g->q[1]);
	g->t[0] = v > 0 ? g->q[0] / v : 0;
	g->t[1] = v > 0 ? g->q[1] / v : 0;
	g->n[0] = g->t[1]; g->n[1] = -g->t[0];
	g->k = v > 0 ? (g->q[0] * q2[1] - g->q[1] * q2[0]) / (v * v * v) : 0;
}

/* the geometry at tau in [0, 1]; where the derivative vanishes (a cusp of
   the curve) the tangent and curvature of a point just inside */
static void os_geo_at(os_run const* r, double tau, os_geo* g)
{
	os_geo_raw(r, tau, g);
	if (!(hypot(g->q[0], g->q[1]) > 1e-13 * (fabs(g->p[0]) + fabs(g->p[1]) + 1)))
	{
		os_geo h;
		double p[2] = { g->p[0], g->p[1] };
		os_geo_raw(r, tau + (tau < 0.5 ? 1e-7 : -1e-7), &h);
		*g = h;
		g->p[0] = p[0]; g->p[1] = p[1];
	}
}

/* the paper's estimate (CubicParams::from_points_derivs, est_euler_err) of
   how far the cubic with these ends and end derivatives (per unit of the
   piece) is from an Euler spiral, times the chord */
static double os_euler_err(double const* p0, double const* p1, double const* q0, double const* q1)
{
	double ch[2] = { p1[0] - p0[0], p1[1] - p0[1] }, ch2 = ch[0] * ch[0] + ch[1] * ch[1], sc;
	double h0x, h0y, h1x, h1y, th0, th1, d0, d1, e0, e1, s0, s1, s01, amin, a, sym, asym, dist, s2;
	if (!(ch2 > 0))
		return hypot(q0[0], q0[1]) + hypot(q1[0], q1[1]);
	/* the cubic's arms are the derivatives / 3 */
	sc = 1.0 / (3.0 * ch2);
	h0x = q0[0] * ch[0] + q0[1] * ch[1]; h0y = q0[1] * ch[0] - q0[0] * ch[1];
	h1x = q1[0] * ch[0] + q1[1] * ch[1]; h1y = q1[0] * ch[1] - q1[1] * ch[0];
	th0 = atan2(h0y, h0x); d0 = hypot(h0x, h0y) * sc;
	th1 = atan2(h1y, h1x); d1 = hypot(h1x, h1y) * sc;
	e0 = (2.0 / 3.0) / fmax(1e-6, 1 + cos(th0));
	e1 = (2.0 / 3.0) / fmax(1e-6, 1 + cos(th1));
	s0 = sin(th0); s1 = sin(th1); s01 = sin(th0 + th1);
	amin = 0.15 * (2 * e0 * s0 + 2 * e1 * s1 - e0 * e1 * s01);
	a = 0.15 * (2 * d0 * s0 + 2 * d1 * s1 - d0 * d1 * s01);
	sym = fabs(th0 + th1); asym = fabs(th0 - th1);
	dist = hypot(d0 - e0, d1 - e1);
	s2 = sym * sym;
	return sqrt(ch2) * (1.25 * (3.7e-6 * s2 * s2 * sym + 6e-3 * asym * s2) + 1.55 * fabs(a - amin) + 5e-3 * sym * dist + 7e-2 * asym * dist);
}

typedef struct os_cuts
{
	double* t;
	os_geo* g;
	unsigned int n, cap;
} os_cuts;

static int os_cut_push(os_cuts* c, double t, os_geo const* g)
{
	if (c->n == c->cap)
	{
		unsigned int nc = c->cap ? 2 * c->cap : 32;
		double* a = (double*)realloc(c->t, sizeof(double) * nc);
		os_geo* b;
		if (!a) return 0;
		c->t = a;
		b = (os_geo*)realloc(c->g, sizeof(os_geo) * nc);
		if (!b) return 0;
		c->g = b;
		c->cap = nc;
	}
	c->t[c->n] = t;
	c->g[c->n] = *g;
	c->n++;
	return 1;
}

/* the Euler cuts of a run (the paper's dyadic walk: halve until the piece
   passes, then step on, doubling back up when aligned) */
static int os_euler_cuts(os_run const* r, double es_tol, os_cuts* out)
{
	os_geo L, G, M;
	double lt = 0, dt = 1;
	unsigned long long t0u = 0;
	unsigned int depth = 0;
	os_geo_at(r, 0, &L);
	if (!os_cut_push(out, 0, &L)) return 0;
	while ((double)t0u * dt < 1)
	{
		double t0 = (double)t0u * dt, t1 = t0 + dt, w, err, h[2];
		if (t1 > 1) t1 = 1;
		os_geo_at(r, t1, &G);
		w = t1 - lt;
		{
			double q0[2] = { L.q[0] * w, L.q[1] * w }, q1[2] = { G.q[0] * w, G.q[1] * w };
			err = os_euler_err(L.p, G.p, q0, q1);
			/* the Hermite cubic of these ends at the middle, against the curve */
			os_geo_at(r, 0.5 * (lt + t1), &M);
			h[0] = 0.5 * (L.p[0] + G.p[0]) + (q0[0] - q1[0]) / 8;
			h[1] = 0.5 * (L.p[1] + G.p[1]) + (q0[1] - q1[1]) / 8;
			err += hypot(M.p[0] - h[0], M.p[1] - h[1]);
		}
		if (err <= es_tol || depth >= OS_ES_DEPTH)
		{
			if (!os_cut_push(out, t1, &G)) return 0;
			L = G;
			lt = t1;
			t0u++;
			while (depth > 0 && (t0u & 1) == 0) { t0u >>= 1; dt *= 2; depth--; }
		}
		else
		{
			t0u *= 2;
			dt *= 0.5;
			depth++;
		}
	}
	return 1;
}

static double os_g(os_run const* r, double tau)
{
	os_geo g;
	os_geo_at(r, tau, &g);
	return 1 + r->delta * g.k;
}

/* the root of g in [ta, tb] (g(ta) ga, g(tb) gb of opposite signs):
   Illinois steps from the linear guess (the Euler spiral's exact root) */
static double os_g_root(os_run const* r, double ta, double ga, double tb, double gb)
{
	unsigned int i;
	int side = 0;
	for (i = 0; i < 100 && tb - ta > 1e-15; i++)
	{
		double t = (ta * gb - tb * ga) / (gb - ga), gt;
		if (!(t > ta && t < tb)) t = 0.5 * (ta + tb);
		gt = os_g(r, t);
		if (gt == 0) return t;
		if ((gt > 0) == (gb > 0))
		{
			tb = t; gb = gt;
			if (side == -1) ga *= 0.5;
			side = -1;
		}
		else
		{
			ta = t; ga = gt;
			if (side == 1) gb *= 0.5;
			side = 1;
		}
	}
	return 0.5 * (ta + tb);
}

/* the extremum of k in [ta, tb] (golden section; maximum or minimum) */
static double os_k_peak(os_run const* r, double ta, double tb, int maximum)
{
	double const ig = 0.6180339887498949;
	double c = tb - ig * (tb - ta), d = ta + ig * (tb - ta), kc, kd;
	os_geo g;
	unsigned int i;
	os_geo_at(r, c, &g); kc = maximum ? g.k : -g.k;
	os_geo_at(r, d, &g); kd = maximum ? g.k : -g.k;
	for (i = 0; i < 80 && tb - ta > 1e-14; i++)
	{
		if (kc > kd)
		{
			tb = d; d = c; kd = kc;
			c = tb - ig * (tb - ta);
			os_geo_at(r, c, &g); kc = maximum ? g.k : -g.k;
		}
		else
		{
			ta = c; c = d; kc = kd;
			d = ta + ig * (tb - ta);
			os_geo_at(r, d, &g); kd = maximum ? g.k : -g.k;
		}
	}
	return 0.5 * (ta + tb);
}

/* a point and the unit direction of travel of the offset (mode 0) or of the
   evolute (mode 1, sgn the sign of k' on the piece) */
static void os_track(os_run const* r, int mode, double sgn, double tau, double* p, double* d)
{
	os_geo g;
	os_geo_at(r, tau, &g);
	if (mode == 0)
	{
		p[0] = g.p[0] + r->delta * g.n[0];
		p[1] = g.p[1] + r->delta * g.n[1];
		d[0] = g.t[0]; d[1] = g.t[1];
	}
	else
	{
		double rho = g.k != 0 ? -1.0 / g.k : -r->delta;
		p[0] = g.p[0] + rho * g.n[0];
		p[1] = g.p[1] + rho * g.n[1];
		d[0] = sgn * g.n[0]; d[1] = sgn * g.n[1];
	}
}

/* the speed of the offset (|C'| |g|) or of the evolute (|k'| / k^2) per unit
   of tau */
static double os_speed(os_run const* r, int mode, double tau)
{
	os_geo g;
	os_geo_at(r, tau, &g);
	if (mode == 0)
		return hypot(g.q[0], g.q[1]) * fabs(1 + r->delta * g.k);
	{
		double h = 1e-6, a = tau - h < 0 ? 0 : tau - h, b = tau + h > 1 ? 1 : tau + h;
		os_geo ga, gb;
		os_geo_at(r, a, &ga);
		os_geo_at(r, b, &gb);
		return g.k != 0 ? fabs((gb.k - ga.k) / (b - a)) / (g.k * g.k) : 0;
	}
}

static void os_fit_push(os_fit* f, double const* b)
{
	unsigned int k;
	if (f->status != QAWS_STATUS_OK) return;
	if (3 * (f->count + 1) + 1 > f->cap)
	{
		unsigned int nc = f->cap ? 2 * f->cap : 64;
		qaws_scalar* g;
		while (nc < 3 * (f->count + 1) + 1) nc *= 2;
		g = (qaws_scalar*)realloc(f->cp, sizeof(qaws_scalar) * 2 * nc);
		if (!g) { f->status = QAWS_STATUS_ALLOCATION_FAILURE; return; }
		f->cp = g;
		f->cap = nc;
	}
	if (f->count == 0)
	{
		f->cp[0] = (qaws_scalar)b[0];
		f->cp[1] = (qaws_scalar)b[1];
	}
	for (k = 0; k < 6; k++)
		f->cp[2 * (3 * f->count + 1) + k] = (qaws_scalar)b[2 + k];
	f->count++;
}

/* a straight piece of the chain (a degenerate cubic) */
static void os_push_line(os_fit* f, double const* a, double const* b)
{
	double c[8];
	c[0] = a[0]; c[1] = a[1];
	c[2] = (2 * a[0] + b[0]) / 3; c[3] = (2 * a[1] + b[1]) / 3;
	c[4] = (a[0] + 2 * b[0]) / 3; c[5] = (a[1] + 2 * b[1]) / 3;
	c[6] = b[0]; c[7] = b[1];
	os_fit_push(f, c);
}

/* the held link, if any, into the chain */
static void os_link_flush(os_fit* f)
{
	if (f->pend)
	{
		f->pend = 0;
		os_push_line(f, f->pe, f->po);
	}
}

static void os_bez(double const* b, double u, double* p, double* d)
{
	double v = 1 - u;
	unsigned int i;
	for (i = 0; i < 2; i++)
	{
		p[i] = v * v * v * b[i] + 3 * u * v * v * b[2 + i] + 3 * u * u * v * b[4 + i] + u * u * u * b[6 + i];
		if (d) d[i] = 3 * (v * v * (b[2 + i] - b[i]) + 2 * u * v * (b[4 + i] - b[2 + i]) + u * u * (b[6 + i] - b[4 + i]));
	}
}

/* u on cubic b nearest q, from u (Newton on (B - q) . B' = 0) */
static double os_bez_project(double const* b, double const* q, double u, unsigned int steps)
{
	unsigned int i;
	for (i = 0; i < steps; i++)
	{
		double p[2], d[2], v = 1 - u, dd[2], num, den;
		unsigned int k;
		os_bez(b, u, p, d);
		for (k = 0; k < 2; k++)
			dd[k] = 6 * (v * (b[4 + k] - 2 * b[2 + k] + b[k]) + u * (b[6 + k] - 2 * b[4 + k] + b[2 + k]));
		num = (p[0] - q[0]) * d[0] + (p[1] - q[1]) * d[1];
		den = d[0] * d[0] + d[1] * d[1] + (p[0] - q[0]) * dd[0] + (p[1] - q[1]) * dd[1];
		if (!(den > 0)) break;
		u -= num / den;
		if (u < 0) u = 0;
		if (u > 1) u = 1;
	}
	return u;
}

/* the signed distances of points X to cubic b (their parameters u refined
   in place), and their sum of squares */
static double os_bez_resid(double const* b, double const (*X)[2], double* u, unsigned int n, double* res)
{
	unsigned int i;
	double s = 0;
	for (i = 0; i < n; i++)
	{
		double p[2], d[2], l;
		u[i] = os_bez_project(b, X[i], u[i], 4);
		os_bez(b, u[i], p, d);
		l = hypot(d[0], d[1]);
		res[i] = l > 0 ? ((X[i][0] - p[0]) * -d[1] + (X[i][1] - p[1]) * d[0]) / l : hypot(X[i][0] - p[0], X[i][1] - p[1]);
		/* off the ends, the distance itself */
		if (u[i] <= 0 || u[i] >= 1)
		{
			double e = hypot(X[i][0] - p[0], X[i][1] - p[1]);
			res[i] = res[i] < 0 ? -e : e;
		}
		s += res[i] * res[i];
	}
	return s;
}

/* the cubic from ta to tb along the offset or the evolute: ends and end
   directions held, arms by least squares, split while off by more than tol */
static void os_g1(os_ctx const* x, os_run const* r, int mode, double sgn, double ta, double tb, os_fit* f, unsigned int depth)
{
	double P0[2], D0[2], P3[2], D3[2], Q[OS_FIT_SAMPLES + 1][2], U[OS_FIT_SAMPLES + 1], M[OS_FIT_SAMPLES][2];
	double b[8], L, err = 0, al = 0, be = 0, acc = 0, flat = 0;
	unsigned int j, it, m = OS_FIT_SAMPLES;
	if (f->status != QAWS_STATUS_OK) return;
	os_track(r, mode, sgn, ta, P0, D0);
	os_track(r, mode, sgn, tb, P3, D3);
	L = hypot(P3[0] - P0[0], P3[1] - P0[1]);
	/* samples Q_1..Q_m-1, and probes M_j halfway between */
	for (j = 1; j < m; j++)
	{
		double d[2];
		os_track(r, mode, sgn, ta + (tb - ta) * j / m, Q[j], d);
	}
	for (j = 0; j < m; j++)
	{
		double d[2];
		os_track(r, mode, sgn, ta + (tb - ta) * (j + 0.5) / m, M[j], d);
	}
	/* straight within the tolerance: a line */
	for (j = 1; j < m; j++)
	{
		double e = L > 0 ? fabs((Q[j][0] - P0[0]) * (P3[1] - P0[1]) - (Q[j][1] - P0[1]) * (P3[0] - P0[0])) / L
			: hypot(Q[j][0] - P0[0], Q[j][1] - P0[1]);
		if (e > flat) flat = e;
	}
	for (j = 0; j < m; j++)
	{
		double e = L > 0 ? fabs((M[j][0] - P0[0]) * (P3[1] - P0[1]) - (M[j][1] - P0[1]) * (P3[0] - P0[0])) / L
			: hypot(M[j][0] - P0[0], M[j][1] - P0[1]);
		if (e > flat) flat = e;
	}
	if (flat <= x->tol)
	{
		/* and monotone along the chord (no doubling back) */
		int mono = 1;
		double last = 0;
		for (j = 1; j < m && mono; j++)
		{
			double s = (Q[j][0] - P0[0]) * (P3[0] - P0[0]) + (Q[j][1] - P0[1]) * (P3[1] - P0[1]);
			if (s < last - x->tol * L) mono = 0;
			last = s;
		}
		if (mono && last <= L * L + x->tol * L)
		{
			os_push_line(f, P0, P3);
			return;
		}
	}
	/* chord-length parameters */
	U[0] = 0;
	{
		double prev[2] = { P0[0], P0[1] };
		for (j = 1; j < m; j++)
		{
			acc += hypot(Q[j][0] - prev[0], Q[j][1] - prev[1]);
			U[j] = acc;
			prev[0] = Q[j][0]; prev[1] = Q[j][1];
		}
		acc += hypot(P3[0] - prev[0], P3[1] - prev[1]);
		for (j = 1; j < m; j++) U[j] = acc > 0 ? U[j] / acc : (double)j / m;
	}
	b[0] = P0[0]; b[1] = P0[1]; b[6] = P3[0]; b[7] = P3[1];
	for (it = 0; it < 4; it++)
	{
		double s11 = 0, s12 = 0, s22 = 0, r1 = 0, r2 = 0, c = D0[0] * D3[0] + D0[1] * D3[1], det;
		for (j = 1; j < m; j++)
		{
			double u = U[j], v = 1 - u, b0 = v * v * v, b1 = 3 * u * v * v, b2 = 3 * u * u * v, b3 = u * u * u, A[2], R[2];
			A[0] = (b0 + b1) * P0[0] + (b2 + b3) * P3[0];
			A[1] = (b0 + b1) * P0[1] + (b2 + b3) * P3[1];
			R[0] = Q[j][0] - A[0]; R[1] = Q[j][1] - A[1];
			s11 += b1 * b1; s22 += b2 * b2; s12 += b1 * b2;
			r1 += b1 * (D0[0] * R[0] + D0[1] * R[1]);
			r2 += b2 * (D3[0] * R[0] + D3[1] * R[1]);
		}
		/* al s11 - be c s12 = r1, al c s12 - be s22 = r2 */
		det = -s11 * s22 + c * c * s12 * s12;
		if (fabs(det) > 1e-12 * s11 * s22)
		{
			al = (r1 * -s22 - (-c * s12) * r2) / det;
			be = (s11 * r2 - c * s12 * r1) / det;
		}
		else
			al = be = -1;
		if (!(al > 0) || !(be > 0) || al > 4 * L + acc || be > 4 * L + acc)
			al = be = (acc > 0 ? acc : L) / 3;
		b[2] = P0[0] + al * D0[0]; b[3] = P0[1] + al * D0[1];
		b[4] = P3[0] - be * D3[0]; b[5] = P3[1] - be * D3[1];
		for (j = 1; j < m; j++)
			U[j] = os_bez_project(b, Q[j], U[j], 2);
	}
	/* Gauss-Newton on the two arms, on the distances of all the points
	   (samples and probes, in order) to the cubic */
	{
		double X[2 * OS_FIT_SAMPLES][2], ux[2 * OS_FIT_SAMPLES], uh[2 * OS_FIT_SAMPLES];
		double r0[2 * OS_FIT_SAMPLES], ra[2 * OS_FIT_SAMPLES], rb[2 * OS_FIT_SAMPLES], cost, h = 1e-7 * (L + acc) + 1e-300;
		unsigned int nx = 0, k;
		U[m] = 1;
		for (j = 0; j < m; j++)
		{
			X[nx][0] = M[j][0]; X[nx][1] = M[j][1]; ux[nx++] = 0.5 * (U[j] + U[j + 1]);
			if (j + 1 < m) { X[nx][0] = Q[j + 1][0]; X[nx][1] = Q[j + 1][1]; ux[nx++] = U[j + 1]; }
		}
		cost = os_bez_resid(b, X, ux, nx, r0);
		/* the Hermite start: arms from the true speeds per unit of tau, which
		   vanish at a cusp (there the cubic's end is a cusp too) */
		{
			double ha = os_speed(r, mode, ta) * (tb - ta) / 3, hb = os_speed(r, mode, tb) * (tb - ta) / 3, c[8], hc;
			double uh2[2 * OS_FIT_SAMPLES], rh[2 * OS_FIT_SAMPLES];
			memcpy(c, b, sizeof(c));
			c[2] = P0[0] + ha * D0[0]; c[3] = P0[1] + ha * D0[1];
			c[4] = P3[0] - hb * D3[0]; c[5] = P3[1] - hb * D3[1];
			for (k = 0; k < nx; k++) uh2[k] = (k + 1) / (2.0 * m);
			hc = os_bez_resid(c, X, uh2, nx, rh);
			if (hc < cost)
			{
				memcpy(b, c, sizeof(c));
				memcpy(ux, uh2, sizeof(double) * nx);
				memcpy(r0, rh, sizeof(double) * nx);
				al = ha; be = hb; cost = hc;
			}
		}
		for (it = 0; it < 12; it++)
		{
			double c[8], jaa = 0, jab = 0, jbb = 0, ga = 0, gb = 0, det, da, db, step = 1;
			memcpy(c, b, sizeof(c));
			c[2] = P0[0] + (al + h) * D0[0]; c[3] = P0[1] + (al + h) * D0[1];
			memcpy(uh, ux, sizeof(double) * nx);
			os_bez_resid(c, X, uh, nx, ra);
			memcpy(c, b, sizeof(c));
			c[4] = P3[0] - (be + h) * D3[0]; c[5] = P3[1] - (be + h) * D3[1];
			memcpy(uh, ux, sizeof(double) * nx);
			os_bez_resid(c, X, uh, nx, rb);
			for (k = 0; k < nx; k++)
			{
				double a1 = (ra[k] - r0[k]) / h, b1 = (rb[k] - r0[k]) / h;
				jaa += a1 * a1; jab += a1 * b1; jbb += b1 * b1;
				ga += a1 * r0[k]; gb += b1 * r0[k];
			}
			det = jaa * jbb - jab * jab;
			if (!(fabs(det) > 1e-30 * (jaa * jbb + 1e-300))) break;
			da = -(jbb * ga - jab * gb) / det;
			db = -(jaa * gb - jab * ga) / det;
			for (k = 0; k < 6; k++, step *= 0.5)
			{
				double na = al + step * da, nb = be + step * db, nc;
				if (!(na >= 0) || !(nb >= 0)) continue;
				memcpy(c, b, sizeof(c));
				c[2] = P0[0] + na * D0[0]; c[3] = P0[1] + na * D0[1];
				c[4] = P3[0] - nb * D3[0]; c[5] = P3[1] - nb * D3[1];
				memcpy(uh, ux, sizeof(double) * nx);
				nc = os_bez_resid(c, X, uh, nx, ra);
				if (nc < cost)
				{
					memcpy(b, c, sizeof(c));
					memcpy(ux, uh, sizeof(double) * nx);
					memcpy(r0, ra, sizeof(double) * nx);
					al = na; be = nb; cost = nc;
					break;
				}
			}
			if (k == 6 || fabs(step * da) + fabs(step * db) < 1e-13 * (L + acc)) break;
		}
		for (k = 0; k < nx; k++)
			if (fabs(r0[k]) > err) err = fabs(r0[k]);
	}
	if (err > x->tol && depth < OS_FIT_DEPTH && tb - ta > 1e-13)
	{
		double tm = 0.5 * (ta + tb);
		os_g1(x, r, mode, sgn, ta, tm, f, depth + 1);
		os_g1(x, r, mode, sgn, tm, tb, f, depth + 1);
		return;
	}
	os_fit_push(f, b);
}

static void os_record_cusp(os_run const* r, double tau)
{
	if (r->cusps && *r->cusp_count < r->cusp_cap)
	{
		double a = r->c->span_boundaries[r->span], b = r->c->span_boundaries[r->span + 1];
		r->cusps[*r->cusp_count] = a + (r->a + tau * (r->b - r->a)) * (b - a);
	}
	if (r->cusp_count) (*r->cusp_count)++;
}

/* a backward run [ta, tb]: O(ta) to E(ta), E to tb cut at the peaks of k,
   E(tb) to O(tb) (the links vanish at cusps) */
static void os_evolute_run(os_ctx const* x, os_run const* r, os_cuts const* cuts, double ta, double tb, os_fit* f)
{
	double cut[64], o[2], e[2], d[2];
	unsigned int nc = 0, i, j;
	os_geo ga, gb;
	os_geo_at(r, ta, &ga);
	os_geo_at(r, tb, &gb);
	/* the peaks of k between cuts of the run (k' changes sign) */
	cut[nc++] = ta;
	{
		double kt[66], tt[66];
		unsigned int n = 0;
		tt[n] = ta; kt[n++] = ga.k;
		for (i = 0; i < cuts->n && n < 65; i++)
			if (cuts->t[i] > ta && cuts->t[i] < tb) { tt[n] = cuts->t[i]; kt[n++] = cuts->g[i].k; }
		tt[n] = tb; kt[n++] = gb.k;
		for (j = 1; j + 1 < n && nc < 62; j++)
		{
			double l = kt[j] - kt[j - 1], rr = kt[j + 1] - kt[j];
			if (l * rr < 0)
				cut[nc++] = os_k_peak(r, tt[j - 1], tt[j + 1], l > 0);
		}
	}
	cut[nc++] = tb;
	os_track(r, 0, 1, ta, o, d);
	os_track(r, 1, 1, ta, e, d);
	if (f->pend && ta == 0 && hypot(e[0] - f->pe[0], e[1] - f->pe[1]) <= x->tol)
		f->pend = 0;    /* the run goes on across a span boundary: no link out and back */
	else
	{
		os_link_flush(f);
		if (hypot(o[0] - e[0], o[1] - e[1]) > x->tol)
			os_push_line(f, o, e);
	}
	for (i = 0; i + 1 < nc; i++)
	{
		os_geo g0, g1;
		double sgn;
		if (!(cut[i + 1] > cut[i])) continue;
		os_geo_at(r, cut[i], &g0);
		os_geo_at(r, cut[i + 1], &g1);
		sgn = g1.k > g0.k ? 1.0 : -1.0;
		os_g1(x, r, 1, sgn, cut[i], cut[i + 1], f, 0);
	}
	os_track(r, 1, 1, tb, e, d);
	os_track(r, 0, 1, tb, o, d);
	if (hypot(o[0] - e[0], o[1] - e[1]) > x->tol)
	{
		if (tb >= 1)
		{
			f->pend = 1;
			f->pe[0] = e[0]; f->pe[1] = e[1];
			f->po[0] = o[0]; f->po[1] = o[1];
		}
		else
			os_push_line(f, e, o);
	}
}

/* the chain of one run: forward runs fitted, backward runs through the evolute */
static qaws_status os_run_chain(os_ctx const* x, os_run const* r, os_fit* f)
{
	os_cuts cuts;
	double es_tol = x->tol > 1e-4 * fabs(r->delta) ? x->tol : 1e-4 * fabs(r->delta);
	double bound[256];
	unsigned int nb = 0, i;
	memset(&cuts, 0, sizeof(cuts));
	if (!os_euler_cuts(r, es_tol, &cuts))
	{
		free(cuts.t); free(cuts.g);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	bound[nb++] = 0;
	for (i = 0; i + 1 < cuts.n && nb < 254; i++)
	{
		double g0 = 1 + r->delta * cuts.g[i].k, g1 = 1 + r->delta * cuts.g[i + 1].k;
		if (g0 * g1 < 0)
		{
			double t = os_g_root(r, cuts.t[i], g0, cuts.t[i + 1], g1);
			if (t > bound[nb - 1])
			{
				bound[nb++] = t;
				os_record_cusp(r, t);
			}
		}
	}
	bound[nb++] = 1;
	for (i = 0; i + 1 < nb; i++)
	{
		double ta = bound[i], tb = bound[i + 1];
		if (!(tb > ta)) continue;
		if (os_g(r, 0.5 * (ta + tb)) >= 0)
		{
			os_link_flush(f);
			os_g1(x, r, 0, 1, ta, tb, f, 0);
		}
		else
			os_evolute_run(x, r, &cuts, ta, tb, f);
	}
	free(cuts.t);
	free(cuts.g);
	return f->status;
}

/* the evaluator's derivative scale on a span: a central difference against
   the first derivative at the middle */
static double os_span_scale(qaws_curve const* c, unsigned int span)
{
	qaws_eval_result_2d e0, e1, em;
	double h = 1e-4, d, v;
	c->vtable->eval_span_2d(c, span, (qaws_scalar)(0.5 - h), QAWS_EVAL_FLAG_POSITION, &e0);
	c->vtable->eval_span_2d(c, span, (qaws_scalar)(0.5 + h), QAWS_EVAL_FLAG_POSITION, &e1);
	c->vtable->eval_span_2d(c, span, (qaws_scalar)0.5, QAWS_EVAL_FLAG_D1, &em);
	d = hypot(e1.position.x - e0.position.x, e1.position.y - e0.position.y) / (2 * h);
	v = hypot(em.d1.x, em.d1.y);
	return v > 0 && d > 0 ? d / v : 1.0;
}

/* the offset of a curve element at a constant delta: one chain per run of
   smooth spans; a corner between spans gets the join */
static qaws_status os_curve_elem(os_ctx* x, os_elem const* e, double delta, qaws_join_type jt,
	double* cusps, unsigned int cusp_cap, unsigned int* cusp_count)
{
	os_fit f;
	double lo = e->t0 < e->t1 ? e->t0 : e->t1, hi = e->t0 < e->t1 ? e->t1 : e->t0;
	double last_o[2] = { 0, 0 }, last_t[2] = { 0, 0 }, last_p[2] = { 0, 0 };
	unsigned int n, started = 0;
	qaws_status st = QAWS_STATUS_OK;
	memset(&f, 0, sizeof(f));
	f.status = QAWS_STATUS_OK;
	for (n = 0; n < e->c->span_count && st == QAWS_STATUS_OK; n++)
	{
		unsigned int k = e->t1 >= e->t0 ? n : e->c->span_count - 1 - n;
		double a = e->c->span_boundaries[k], b = e->c->span_boundaries[k + 1], w = b - a, la, lb;
		os_run r;
		os_geo g0, g1;
		double o0[2];
		if (!(w > 0)) continue;
		la = ((a < lo ? lo : a) - a) / w;
		lb = ((b > hi ? hi : b) - a) / w;
		if (!(lb > la)) continue;
		memset(&r, 0, sizeof(r));
		r.c = e->c; r.span = k; r.delta = delta;
		r.a = e->t1 >= e->t0 ? la : lb;
		r.b = e->t1 >= e->t0 ? lb : la;
		r.scale = os_span_scale(e->c, k);
		r.cusps = cusps; r.cusp_cap = cusp_cap; r.cusp_count = cusp_count;
		os_geo_at(&r, 0, &g0);
		o0[0] = g0.p[0] + delta * g0.n[0];
		o0[1] = g0.p[1] + delta * g0.n[1];
		if (started && hypot(o0[0] - last_o[0], o0[1] - last_o[1]) > x->tol)
		{
			/* a corner inside the element (a gap within the tolerance is rounding
			   between spans: the chain runs on) */
			st = os_fit_emit(x, &f);
			if (st == QAWS_STATUS_OK)
				st = os_join(x, jt, last_p, last_t, g0.t, last_o, o0, delta);
			if (st != QAWS_STATUS_OK) break;
		}
		st = os_run_chain(x, &r, &f);
		os_geo_at(&r, 1, &g1);
		last_o[0] = g1.p[0] + delta * g1.n[0];
		last_o[1] = g1.p[1] + delta * g1.n[1];
		last_t[0] = g1.t[0]; last_t[1] = g1.t[1];
		last_p[0] = g1.p[0]; last_p[1] = g1.p[1];
		started = 1;
	}
	if (st != QAWS_STATUS_OK)
	{
		free(f.cp);
		return st;
	}
	return os_fit_emit(x, &f);
}

/* the offset of an element; *s, *e receive its start and end points */
static qaws_status os_offset_elem(os_ctx* x, os_elem const* el, double delta, qaws_join_type jt, double* sp, double* ep)
{
	double n0[2], n1[2], d0, d1;
	os_right(el->d0, n0);
	os_right(el->d1, n1);
	d0 = os_delta(x, delta, el->p0, n0);
	d1 = os_delta(x, delta, el->p1, n1);
	sp[0] = el->p0[0] + d0 * n0[0]; sp[1] = el->p0[1] + d0 * n0[1];
	ep[0] = el->p1[0] + d1 * n1[0]; ep[1] = el->p1[1] + d1 * n1[1];
	if (el->line)
		return os_line(x, sp, ep);
	if (el->arc && !x->d->delta_fn)
	{
		/* the arc piece's own segment: concentric, radius moved by delta */
		qaws_curve* piece = NULL;
		qaws_status s = qaws_curve_extract(el->c, (qaws_scalar)el->t0, (qaws_scalar)el->t1, &piece);
		if (s == QAWS_STATUS_OK && piece->kind == QAWS_CURVE_KIND_ARC)
		{
			qaws_arc_impl const* impl = (qaws_arc_impl const*)piece->impl;
			qaws_arc_segment const* g = &impl->segments[0];
			double sw = (double)(g->angle_end - g->angle_start), c[2] = { g->center[0], g->center[1] };
			/* travelling counter-clockwise, the right is outward */
			double r = (double)g->radius + (sw > 0 ? delta : -delta);
			if (impl->segment_count == 1 && r > 0)
			{
				s = os_arc(x, c, r, (double)g->angle_start, sw);
				qaws_curve_destroy(piece);
				return s;
			}
		}
		qaws_curve_destroy(piece);
	}
	if (x->d->delta_fn)
		return os_fit_elem(x, el, delta);
	return os_curve_elem(x, el, delta, jt, NULL, 0, NULL);
}

/* ======================================================================== */
/*  Joins and caps                                                          */
/* ======================================================================== */

static double os_cross(double const* a, double const* b) { return a[0] * b[1] - a[1] * b[0]; }
static double os_dot(double const* a, double const* b) { return a[0] * b[0] + a[1] * b[1]; }

/* from P (end of the incoming offset) to Q (start of the outgoing one) round
   corner v, tangents tin / tout, offset delta to the right */
static qaws_status os_join(os_ctx* x, qaws_join_type jt, double const* v, double const* tin, double const* tout,
	double const* P, double const* Q, double delta)
{
	double turn = os_cross(tin, tout), u0[2] = { P[0] - v[0], P[1] - v[1] }, u1[2] = { Q[0] - v[0], Q[1] - v[1] };
	double r0 = hypot(u0[0], u0[1]), r1 = hypot(u1[0], u1[1]), r = 0.5 * (r0 + r1);
	double eps = 1e-12 * (fabs(v[0]) + fabs(v[1]) + r + 1);
	qaws_status s;
	if (hypot(P[0] - Q[0], P[1] - Q[1]) <= eps)
		return QAWS_STATUS_OK;
	/* inner side: through the corner (the union removes the overlap) */
	if (turn * delta < 0 || (fabs(turn) < 1e-12 && os_dot(tin, tout) > 0))
	{
		if (fabs(turn) < 1e-12 && os_dot(tin, tout) > 0)
			return os_line(x, P, Q);
		s = os_line(x, P, v);
		return s == QAWS_STATUS_OK ? os_line(x, v, Q) : s;
	}
	if (r0 > 0) { u0[0] /= r0; u0[1] /= r0; }
	if (r1 > 0) { u1[0] /= r1; u1[1] /= r1; }
	switch (jt)
	{
	case QAWS_JOIN_ROUND:
	{
		double a0 = atan2(u0[1], u0[0]), sw = atan2(os_cross(u0, u1), os_dot(u0, u1));
		/* a half turn: round the side the path came from (ahead of P) */
		if (fabs(os_cross(u0, u1)) < 1e-12 && os_dot(u0, u1) < 0)
			sw = os_cross(u0, tin) > 0 ? OS_PI : -OS_PI;
		return os_arc(x, v, r, a0, sw);
	}
	case QAWS_JOIN_MITER:
	{
		double cosa = os_dot(u0, u1);
		if (1 + cosa > x->miter_lim)
		{
			double q = r / (1 + cosa), m[2];
			m[0] = v[0] + (u0[0] + u1[0]) * q;
			m[1] = v[1] + (u0[1] + u1[1]) * q;
			s = os_line(x, P, m);
			return s == QAWS_STATUS_OK ? os_line(x, m, Q) : s;
		}
	}
	/* beyond the miter limit: square */
	/* fall through */
	case QAWS_JOIN_SQUARE:
	{
		double b[2] = { u0[0] + u1[0], u0[1] + u1[1] }, c[2], la, lb, A[2], B[2], nb;
		nb = hypot(b[0], b[1]);
		if (nb < 1e-12) { b[0] = tin[0]; b[1] = tin[1]; nb = 1; }
		b[0] /= nb; b[1] /= nb;
		c[0] = v[0] + r * b[0]; c[1] = v[1] + r * b[1];
		if (fabs(os_dot(tin, b)) < 1e-12 || fabs(os_dot(tout, b)) < 1e-12)
			return os_line(x, P, Q);
		la = ((c[0] - P[0]) * b[0] + (c[1] - P[1]) * b[1]) / os_dot(tin, b);
		lb = ((c[0] - Q[0]) * b[0] + (c[1] - Q[1]) * b[1]) / os_dot(tout, b);
		A[0] = P[0] + la * tin[0]; A[1] = P[1] + la * tin[1];
		B[0] = Q[0] + lb * tout[0]; B[1] = Q[1] + lb * tout[1];
		s = os_line(x, P, A);
		if (s == QAWS_STATUS_OK) s = os_line(x, A, B);
		if (s == QAWS_STATUS_OK) s = os_line(x, B, Q);
		return s;
	}
	default:
		return os_line(x, P, Q);
	}
}

/* the cap at an open end E (tangent t pointing out of the path) from R (the
   right offset arriving) to L (the left one leaving) */
static qaws_status os_cap(os_ctx* x, qaws_end_type et, double const* E, double const* t, double const* R, double const* L, double delta)
{
	double d = fabs(delta);
	qaws_status s;
	switch (et)
	{
	case QAWS_END_SQUARE:
	{
		double a[2] = { R[0] + d * t[0], R[1] + d * t[1] }, b[2] = { L[0] + d * t[0], L[1] + d * t[1] };
		s = os_line(x, R, a);
		if (s == QAWS_STATUS_OK) s = os_line(x, a, b);
		if (s == QAWS_STATUS_OK) s = os_line(x, b, L);
		return s;
	}
	case QAWS_END_ROUND:
	{
		double u[2] = { R[0] - E[0], R[1] - E[1] }, r = hypot(u[0], u[1]);
		return os_arc(x, E, r, atan2(u[1], u[0]), OS_PI);
	}
	default:
		return os_line(x, R, L);
	}
}

/* ======================================================================== */
/*  Paths                                                                   */
/* ======================================================================== */

/* the start of element e's offset, from its geometry (what os_offset_elem makes) */
static void os_start(os_ctx const* x, os_elem const* e, double delta, double* q)
{
	double n0[2], dl;
	os_right(e->d0, n0);
	dl = os_delta(x, delta, e->p0, n0);
	q[0] = e->p0[0] + dl * n0[0];
	q[1] = e->p0[1] + dl * n0[1];
}

/* one loop round the elements e[0..n) (cyclic) offset by delta: offset i,
   then the join from it to offset i + 1 */
static qaws_status os_closed_loop(os_ctx* x, os_elem const* e, unsigned int n, double delta, qaws_join_type jt)
{
	unsigned int first = x->ncurves, i;
	qaws_status s = QAWS_STATUS_OK;
	for (i = 0; i < n && s == QAWS_STATUS_OK; i++)
	{
		double sp[2], ep[2], q[2];
		s = os_offset_elem(x, &e[i], delta, jt, sp, ep);
		if (s != QAWS_STATUS_OK) break;
		os_start(x, &e[(i + 1) % n], delta, q);
		s = os_join(x, jt, e[i].p1, e[i].d1, e[(i + 1) % n].d0, ep, q, delta);
	}
	return s == QAWS_STATUS_OK ? os_loop_end(x, first) : s;
}

/* a cap needs the left point too: computed from the end's geometry */
static qaws_status os_cap_points(os_ctx* x, qaws_end_type et, double const* E, double const* t, double const* R, double delta)
{
	double n[2], L[2];
	os_right(t, n);
	/* t points along travel into the end; the left offset is -n */
	L[0] = E[0] - fabs(delta) * n[0];
	L[1] = E[1] - fabs(delta) * n[1];
	return os_cap(x, et, E, t, R, L, delta);
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

static qaws_status os_path(os_ctx* x, qaws_path_2d const* p, double delta, qaws_join_type jt, qaws_end_type et)
{
	os_elems q;
	unsigned int i;
	qaws_status s = QAWS_STATUS_OK;
	memset(&q, 0, sizeof(q));
	for (i = 0; i < p->curve_count && s == QAWS_STATUS_OK; i++)
	{
		if (!p->curves[i] || p->curves[i]->dimension != QAWS_DIMENSION_2D)
			s = QAWS_STATUS_INVALID_ARGUMENT;
		else
			s = os_curve_elems(p->curves[i], &q);
	}
	if (s == QAWS_STATUS_OK && q.n == 0 && p->curve_count)
	{
		/* a point: a circle (round) or a square of half-side |delta| */
		double c[2], d = fabs(delta);
		unsigned int first = x->ncurves;
		os_eval(p->curves[0], p->curves[0]->parameter_range.min_value, c, NULL, NULL);
		if (jt == QAWS_JOIN_ROUND || et == QAWS_END_ROUND)
			s = os_arc(x, c, d, 0, 2 * OS_PI);
		else
		{
			double a[2] = { c[0] - d, c[1] - d }, b[2] = { c[0] + d, c[1] - d }, e2[2] = { c[0] + d, c[1] + d }, f[2] = { c[0] - d, c[1] + d };
			s = os_line(x, a, b);
			if (s == QAWS_STATUS_OK) s = os_line(x, b, e2);
			if (s == QAWS_STATUS_OK) s = os_line(x, e2, f);
			if (s == QAWS_STATUS_OK) s = os_line(x, f, a);
		}
		if (s == QAWS_STATUS_OK) s = os_loop_end(x, first);
	}
	else if (s == QAWS_STATUS_OK && q.n)
	{
		if (et == QAWS_END_POLYGON)
			s = os_closed_loop(x, q.e, q.n, delta, jt);
		else if (et == QAWS_END_JOINED)
		{
			/* both sides joined: the path there and back as one closed loop */
			os_elems both;
			memset(&both, 0, sizeof(both));
			for (i = 0; i < q.n && s == QAWS_STATUS_OK; i++)
				if (!os_push_elem(&both, &q.e[i])) s = QAWS_STATUS_ALLOCATION_FAILURE;
			for (i = q.n; i-- > 0 && s == QAWS_STATUS_OK;)
			{
				os_elem r = q.e[i];
				double t = r.t0;
				r.t0 = r.t1; r.t1 = t;
				memcpy(r.p0, q.e[i].p1, sizeof(r.p0));
				memcpy(r.p1, q.e[i].p0, sizeof(r.p1));
				r.d0[0] = -q.e[i].d1[0]; r.d0[1] = -q.e[i].d1[1];
				r.d1[0] = -q.e[i].d0[0]; r.d1[1] = -q.e[i].d0[1];
				if (!os_push_elem(&both, &r)) s = QAWS_STATUS_ALLOCATION_FAILURE;
			}
			if (s == QAWS_STATUS_OK)
				s = os_closed_loop(x, both.e, both.n, fabs(delta), jt);
			free(both.e);
		}
		else
		{
			unsigned int first = x->ncurves, n = q.n;
			double* se = (double*)malloc(sizeof(double) * 4 * (2 * n + 1));
			os_elem* both = (os_elem*)malloc(sizeof(os_elem) * (2 * n + 1));
			if (!se || !both) s = QAWS_STATUS_ALLOCATION_FAILURE;
			for (i = 0; s == QAWS_STATUS_OK && i < n; i++)
			{
				os_elem r = q.e[n - 1 - i];
				double t = r.t0;
				both[i] = q.e[i];
				r.t0 = r.t1; r.t1 = t;
				memcpy(r.p0, q.e[n - 1 - i].p1, sizeof(r.p0));
				memcpy(r.p1, q.e[n - 1 - i].p0, sizeof(r.p1));
				r.d0[0] = -q.e[n - 1 - i].d1[0]; r.d0[1] = -q.e[n - 1 - i].d1[1];
				r.d1[0] = -q.e[n - 1 - i].d0[0]; r.d1[1] = -q.e[n - 1 - i].d0[1];
				both[n + i] = r;
			}
			/* offsets, joins between them, caps at the two turnarounds */
			for (i = 0; s == QAWS_STATUS_OK && i < 2 * n; i++)
			{
				s = os_offset_elem(x, &both[i], fabs(delta), jt, &se[4 * i], &se[4 * i + 2]);
				if (s != QAWS_STATUS_OK) break;
				if (i + 1 == n || i + 1 == 2 * n)
					s = os_cap_points(x, et, both[i].p1, both[i].d1, &se[4 * i + 2], delta);
				else
				{
					double n0[2], qq[2], dl;
					os_right(both[i + 1].d0, n0);
					dl = os_delta(x, fabs(delta), both[i + 1].p0, n0);
					qq[0] = both[i + 1].p0[0] + dl * n0[0];
					qq[1] = both[i + 1].p0[1] + dl * n0[1];
					s = os_join(x, jt, both[i].p1, both[i].d1, both[i + 1].d0, &se[4 * i + 2], qq, fabs(delta));
				}
			}
			free(se);
			free(both);
			if (s == QAWS_STATUS_OK) s = os_loop_end(x, first);
		}
	}
	free(q.e);
	return s;
}

qaws_status qaws_offset_execute(qaws_offset_desc const* desc, qaws_clip_result** out_result)
{
	os_ctx x;
	qaws_status s = QAWS_STATUS_OK;
	unsigned int g, i, first_reversed = 0;
	qaws_path_2d* loops = NULL;
	qaws_curve const** views = NULL;
	double ml, ext = 1.0;
	if (!desc || !out_result || (desc->group_count && !desc->groups))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_result = NULL;
	memset(&x, 0, sizeof(x));
	x.d = desc;
	ml = desc->miter_limit > 0 ? (double)desc->miter_limit : 2.0;
	x.miter_lim = ml <= 1 ? 2.0 : 2.0 / (ml * ml);
	/* the fit tolerance from the inputs' extent */
	{
		unsigned int n = 0, k, j;
		qaws_curve const** all;
		for (g = 0; g < desc->group_count; g++)
			for (k = 0; k < desc->groups[g].path_count; k++)
				n += desc->groups[g].paths[k].curve_count;
		all = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (n + 1));
		if (!all) return QAWS_STATUS_ALLOCATION_FAILURE;
		n = 0;
		for (g = 0; g < desc->group_count; g++)
			for (k = 0; k < desc->groups[g].path_count; k++)
				for (j = 0; j < desc->groups[g].paths[k].curve_count; j++)
					if (desc->groups[g].paths[k].curves[j])
						all[n++] = desc->groups[g].paths[k].curves[j];
		if (n) ext = (double)qaws_internal_flatten_extent(all, n, 2, NULL, 0);
		free((void*)all);
	}
	x.tol = desc->tolerance > 0 ? (double)desc->tolerance : 1e-6 * (fabs((double)desc->delta) > ext ? fabs((double)desc->delta) : ext);
	for (g = 0; g < desc->group_count && s == QAWS_STATUS_OK; g++)
	{
		qaws_offset_group const* gr = &desc->groups[g];
		double delta = (double)desc->delta;
		unsigned int loops_before = x.nloops;
		int reversed = 0;
		/* a closed group whose outermost path runs clockwise is reversed */
		if (gr->end_type == QAWS_END_POLYGON && gr->path_count)
		{
			double lowest = 1e300;
			unsigned int best = 0;
			for (i = 0; i < gr->path_count; i++)
			{
				qaws_vec2 lo, hi;
				if (gr->paths[i].curve_count && qaws_path_compute_bounds_2d(&gr->paths[i], &lo, &hi) == QAWS_STATUS_OK && lo.y < lowest)
				{
					lowest = lo.y;
					best = i;
				}
			}
			reversed = !qaws_path_is_positive_2d(&gr->paths[best]);
			if (reversed)
				delta = -delta;
		}
		if (g == 0)
			first_reversed = reversed;
		x.group = g;
		for (i = 0; i < gr->path_count && s == QAWS_STATUS_OK; i++)
		{
			x.path = i;
			s = os_path(&x, &gr->paths[i], delta, gr->join_type, gr->end_type);
		}
		/* a reversed group's loops run clockwise: turn them round for the union */
		if (reversed && s == QAWS_STATUS_OK)
		{
			unsigned int l;
			for (l = loops_before; l < x.nloops && s == QAWS_STATUS_OK; l++)
			{
				unsigned int a = x.loop_first[l], b = x.loop_first[l + 1], k;
				for (k = a; k < b && s == QAWS_STATUS_OK; k++)
				{
					qaws_curve* rv = NULL;
					s = qaws_curve_reverse(x.curves[k], &rv);
					if (s == QAWS_STATUS_OK)
					{
						qaws_curve_destroy(x.curves[k]);
						x.curves[k] = rv;
					}
				}
				for (k = 0; k < (b - a) / 2; k++)
				{
					qaws_curve* t = x.curves[a + k];
					x.curves[a + k] = x.curves[b - 1 - k];
					x.curves[b - 1 - k] = t;
				}
			}
		}
	}
	if (s == QAWS_STATUS_OK)
	{
		qaws_clip_desc cd;
		loops = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (x.nloops + 1));
		views = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (x.ncurves + 1));
		if (!loops || !views)
			s = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			for (i = 0; i < x.ncurves; i++) views[i] = x.curves[i];
			for (i = 0; i < x.nloops; i++)
			{
				loops[i].curves = views + x.loop_first[i];
				loops[i].curve_count = x.loop_first[i + 1] - x.loop_first[i];
				loops[i].closed = 1;
			}
			memset(&cd, 0, sizeof(cd));
			cd.subjects = loops;
			cd.subject_count = x.nloops;
			cd.clip_type = QAWS_CLIP_UNION;
			cd.fill_rule = QAWS_FILL_POSITIVE;
			cd.flags = desc->flags ^ (first_reversed ? QAWS_CLIP_REVERSE_SOLUTION : 0u);
			cd.executor = desc->executor;
			s = qaws_clip_execute(&cd, out_result);
		}
	}
	for (i = 0; i < x.ncurves; i++)
		qaws_curve_destroy(x.curves[i]);
	free(x.curves);
	free(x.loop_first);
	free(loops);
	free((void*)views);
	return s;
}

qaws_status qaws_offset_paths(qaws_path_2d const* paths, unsigned int path_count, qaws_scalar delta,
	qaws_join_type join_type, qaws_end_type end_type, qaws_clip_result** out_result)
{
	qaws_offset_group g;
	qaws_offset_desc d;
	g.paths = paths;
	g.path_count = path_count;
	g.join_type = join_type;
	g.end_type = end_type;
	memset(&d, 0, sizeof(d));
	d.groups = &g;
	d.group_count = 1;
	d.delta = delta;
	return qaws_offset_execute(&d, out_result);
}

qaws_status qaws_internal_offset_curve(qaws_curve const* curve, double delta, double tolerance,
	qaws_curve** out_chain, unsigned int* out_curve_count,
	double* cusps, unsigned int cusp_capacity, unsigned int* cusp_count)
{
	os_ctx x;
	qaws_offset_desc d;
	os_elems q;
	unsigned int i, nc = 0;
	qaws_status s;
	if (!curve || !out_chain || curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_chain = NULL;
	if (cusp_count) *cusp_count = 0;
	memset(&x, 0, sizeof(x));
	memset(&d, 0, sizeof(d));
	memset(&q, 0, sizeof(q));
	x.d = &d;
	x.tol = tolerance;
	s = os_elem_make(curve, curve->parameter_range.min_value, curve->parameter_range.max_value, 0, 0, &q);
	if (s == QAWS_STATUS_OK && q.n)
		s = os_curve_elem(&x, &q.e[0], delta, QAWS_JOIN_ROUND, cusps, cusp_capacity, cusp_count ? cusp_count : &nc);
	if (s == QAWS_STATUS_OK && x.ncurves)
	{
		*out_chain = x.curves[0];
		for (i = 1; i < x.ncurves; i++)
			qaws_curve_destroy(x.curves[i]);
	}
	else
		for (i = 0; i < x.ncurves; i++)
			qaws_curve_destroy(x.curves[i]);
	if (out_curve_count) *out_curve_count = x.ncurves;
	free(x.curves);
	free(x.loop_first);
	free(q.e);
	return s;
}
