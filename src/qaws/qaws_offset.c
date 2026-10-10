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
} os_fit;

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

/* the offset of a curve element: a cubic B-spline (knots 0..count, multiplicity 3) */
static qaws_status os_fit_elem(os_ctx* x, os_elem const* e, double delta)
{
	os_fit f;
	double lo = e->t0 < e->t1 ? e->t0 : e->t1, hi = e->t0 < e->t1 ? e->t1 : e->t0;
	unsigned int n;
	qaws_status st = QAWS_STATUS_OK;
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
	if (f.status == QAWS_STATUS_OK && f.count)
	{
		unsigned int np = 3 * f.count + 1, kc = 0, i;
		qaws_scalar* kn = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (np + 4));
		qaws_bspline_desc d;
		qaws_curve* c = NULL;
		if (!kn) { free(f.cp); return QAWS_STATUS_ALLOCATION_FAILURE; }
		for (i = 0; i < 4; i++) kn[kc++] = 0;
		for (i = 1; i < f.count; i++) { kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; }
		for (i = 0; i < 4; i++) kn[kc++] = (qaws_scalar)f.count;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = f.cp; d.control_point_count = np;
		d.knots = kn; d.knot_count = kc;
		st = qaws_curve_create_bspline(&d, &c);
		if (st == QAWS_STATUS_OK)
			st = os_keep(x, c);
		free(kn);
	}
	else
		st = f.status;
	free(f.cp);
	return st;
}

/* the offset of an element; *s, *e receive its start and end points */
static qaws_status os_offset_elem(os_ctx* x, os_elem const* el, double delta, double* sp, double* ep)
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
	return os_fit_elem(x, el, delta);
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
		s = os_offset_elem(x, &e[i], delta, sp, ep);
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
				s = os_offset_elem(x, &both[i], fabs(delta), &se[4 * i], &se[4 * i + 2]);
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
