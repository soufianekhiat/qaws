/*
 * Borrowed curves (qaws_view.h): array views over host memory, prepared
 * evaluators that own only derived data (the knots in the computing
 * precision, the span bounds), and qaws_curve objects over them. The points
 * and weights are read where the host keeps them, one value at a time,
 * through the same basis kernels as the owned B-spline and NURBS curves.
 */

#include "qaws_view.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_basis.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_platform.h"
#include "qaws_surface_bspline.h"
#include "qaws_surface_nurbs.h"
#include "internal/qaws_internal_surface.h"
#include "qaws_alloc.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_polynomial.h"
#include "qaws_prepare.h"
#include "core/qaws_cubic_poly_core.h"
#include "core/qaws_horner_core.h"
#include "internal/qaws_internal_diff.h"
#include "internal/qaws_internal_closest.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- array views -------------------------------------------------------- */

static size_t vw_size(qaws_element_type t)
{
	return t == QAWS_ELEMENT_F32 ? 4 : 8;
}

static size_t vw_step(void const* data, size_t stride, qaws_element_type t, unsigned int components)
{
	return stride ? stride : data ? vw_size(t) * components : vw_size(t);
}

/* element i's component c */
static double vw_get(qaws_array_view const* v, unsigned int i, unsigned int c)
{
	size_t s = vw_step(v->data, v->stride, v->type, v->components);
	unsigned char const* p = v->data ? (unsigned char const*)v->data + (size_t)i * s + c * vw_size(v->type)
		: (unsigned char const*)v->planes[c] + (size_t)i * s;
	return v->type == QAWS_ELEMENT_F32 ? (double)*(float const*)p : *(double const*)p;
}

static void vw_put(qaws_array_view_mut const* v, unsigned int i, unsigned int c, double x)
{
	size_t s = vw_step(v->data, v->stride, v->type, v->components);
	unsigned char* p = v->data ? (unsigned char*)v->data + (size_t)i * s + c * vw_size(v->type) : (unsigned char*)v->planes[c] + (size_t)i * s;
	if (v->type == QAWS_ELEMENT_F32)
		*(float*)p = (float)x;
	else
		*(double*)p = x;
}

static qaws_array_view vw_const(qaws_array_view_mut const* m)
{
	qaws_array_view v;
	unsigned int c;
	v.data = m->data;
	for (c = 0; c < 4; c++)
		v.planes[c] = m->planes[c];
	v.type = m->type;
	v.components = m->components;
	v.count = m->count;
	v.stride = m->stride;
	return v;
}

qaws_status qaws_array_view_validate(qaws_array_view const* view)
{
	size_t es;
	unsigned int c;
	if (!view || (view->type != QAWS_ELEMENT_F32 && view->type != QAWS_ELEMENT_F64) || view->components < 1 || view->components > 4)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (view->count == 0)
		return QAWS_STATUS_OK;
	es = vw_size(view->type);
	if (view->stride % es)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (view->data)
		return (uintptr_t)view->data % es || (view->stride && view->stride < es * view->components) ? QAWS_STATUS_INVALID_ARGUMENT
			: QAWS_STATUS_OK;
	for (c = 0; c < view->components; c++)
		if (!view->planes[c] || (uintptr_t)view->planes[c] % es)
			return QAWS_STATUS_INVALID_ARGUMENT;
	return QAWS_STATUS_OK;
}

/* [first, first + span) inside [base, base + bytes)? */
static int vw_inside(void const* first, size_t span, void const* base, size_t bytes)
{
	uintptr_t f = (uintptr_t)first, b = (uintptr_t)base;
	return f >= b && f - b <= bytes && span <= bytes - (f - b);
}

qaws_status qaws_array_view_check_extent(qaws_array_view const* view, void const* base, size_t bytes)
{
	qaws_status st = qaws_array_view_validate(view);
	size_t es, s, span;
	unsigned int c;
	if (st != QAWS_STATUS_OK)
		return st;
	if (view->count == 0)
		return QAWS_STATUS_OK;
	if (!base)
		return QAWS_STATUS_INVALID_ARGUMENT;
	es = vw_size(view->type);
	s = vw_step(view->data, view->stride, view->type, view->components);
	if ((size_t)(view->count - 1) > ((size_t)-1 - es * view->components) / s)
		return QAWS_STATUS_OUT_OF_RANGE;
	if (view->data)
	{
		span = (size_t)(view->count - 1) * s + es * view->components;
		return vw_inside(view->data, span, base, bytes) ? QAWS_STATUS_OK : QAWS_STATUS_OUT_OF_RANGE;
	}
	span = (size_t)(view->count - 1) * s + es;
	for (c = 0; c < view->components; c++)
		if (!vw_inside(view->planes[c], span, base, bytes))
			return QAWS_STATUS_OUT_OF_RANGE;
	return QAWS_STATUS_OK;
}

/* ---- curve views -------------------------------------------------------- */

/* n points of dim components (or dim + 1 homogeneous) with optional weights: layouts, counts, finite values, positive weights */
static qaws_status vw_check_points(qaws_array_view const* points, qaws_array_view const* weights, unsigned int dim, unsigned int n)
{
	unsigned int i, c;
	qaws_status st;
	if ((st = qaws_array_view_validate(points)) != QAWS_STATUS_OK)
		return st;
	if (points->components != dim && points->components != dim + 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (points->count != n)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (weights->count)
	{
		if (points->components != dim || weights->count != n || weights->components != 1)
			return QAWS_STATUS_INVALID_WEIGHT_COUNT;
		if ((st = qaws_array_view_validate(weights)) != QAWS_STATUS_OK)
			return st;
	}
	for (i = 0; i < n; i++)
	{
		double w = weights->count ? vw_get(weights, i, 0) : points->components == dim + 1 ? vw_get(points, i, dim) : 1;
		if (!(w > 0) || !isfinite(w))
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (c = 0; c < points->components; c++)
			if (!isfinite(vw_get(points, i, c)))
				return QAWS_STATUS_INVALID_ARGUMENT;
	}
	return QAWS_STATUS_OK;
}

/* knots for n points of degree d: none (count 0), or n + d + 1 finite non-decreasing values over a non-empty domain */
static qaws_status vw_check_knots(qaws_array_view const* knots, unsigned int n, unsigned int d)
{
	unsigned int i;
	qaws_status st;
	if (!knots->count)
		return QAWS_STATUS_OK;
	if (knots->count != n + d + 1 || knots->components != 1)
		return QAWS_STATUS_INVALID_KNOT_COUNT;
	if ((st = qaws_array_view_validate(knots)) != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < knots->count; i++)
	{
		double k = vw_get(knots, i, 0);
		if (!isfinite(k) || (i && k < vw_get(knots, i - 1, 0)))
			return QAWS_STATUS_INVALID_KNOT_VECTOR;
	}
	return vw_get(knots, d, 0) < vw_get(knots, n, 0) ? QAWS_STATUS_OK : QAWS_STATUS_INVALID_KNOT_VECTOR;
}

/* the n + d + 1 knots in the computing precision: the view's, or clamped uniform ones as the owned constructors make them
   (curves 0, 1, ..., n - d; surfaces over [0, 1]) */
static void vw_knots(qaws_array_view const* knots, unsigned int n, unsigned int d, int unit, qaws_scalar* out)
{
	unsigned int k, kc = n + d + 1;
	if (knots->count)
		for (k = 0; k < kc; k++)
			out[k] = (qaws_scalar)vw_get(knots, k, 0);
	else if (unit)
		for (k = 0; k < kc; k++)
			out[k] = k <= d ? QAWS_ZERO : k >= n ? QAWS_ONE : (qaws_scalar)(k - d) / (qaws_scalar)(n - d);
	else
		for (k = 0; k < kc; k++)
			out[k] = (qaws_scalar)(k <= d ? 0 : k >= n ? n - d : k - d);
}

qaws_status qaws_curve_view_validate(qaws_curve_view const* view)
{
	unsigned int dim, n, d;
	qaws_status st;
	if (!view)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (view->dimension != QAWS_DIMENSION_2D && view->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	dim = (unsigned int)view->dimension;
	n = view->points.count;
	d = view->degree;
	if (d < 1 || d > QAWS_VIEW_MAX_DEGREE)
		return QAWS_STATUS_INVALID_DEGREE;
	switch (view->family)
	{
	case QAWS_VIEW_SPLINE:
		if (n < d + 1)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		if ((st = vw_check_points(&view->points, &view->weights, dim, n)) != QAWS_STATUS_OK)
			return st;
		return vw_check_knots(&view->knots, n, d);
	case QAWS_VIEW_HERMITE:
	case QAWS_VIEW_CATMULL_ROM:
		if (d != 3)
			return QAWS_STATUS_INVALID_DEGREE;
		if (n < (view->family == QAWS_VIEW_HERMITE ? 2u : view->closed ? 3u : 4u))
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		if (view->weights.count || view->knots.count || view->points.components != dim)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if ((st = vw_check_points(&view->points, &view->weights, dim, n)) != QAWS_STATUS_OK)
			return st;
		if (view->family == QAWS_VIEW_CATMULL_ROM)
			return (unsigned int)view->parameterization <= QAWS_PARAMETERIZATION_CENTRIPETAL ? QAWS_STATUS_OK : QAWS_STATUS_INVALID_ARGUMENT;
		if (view->tangents.count != n)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		if (view->tangents.components != dim)
			return QAWS_STATUS_INVALID_ARGUMENT;
		return vw_check_points(&view->tangents, &view->weights, dim, n);
	case QAWS_VIEW_POLYNOMIAL:
		if (n != d + 1)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		if (view->weights.count || view->points.components != dim)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if ((st = vw_check_points(&view->points, &view->weights, dim, n)) != QAWS_STATUS_OK)
			return st;
		if (!view->knots.count)
			return QAWS_STATUS_OK;
		if (view->knots.count != 2 || view->knots.components != 1)
			return QAWS_STATUS_INVALID_KNOT_COUNT;
		if ((st = qaws_array_view_validate(&view->knots)) != QAWS_STATUS_OK)
			return st;
		return isfinite(vw_get(&view->knots, 0, 0)) && vw_get(&view->knots, 0, 0) < vw_get(&view->knots, 1, 0)
			&& isfinite(vw_get(&view->knots, 1, 0)) ? QAWS_STATUS_OK : QAWS_STATUS_INVALID_PARAMETER_RANGE;
	default:
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
}

/* ---- prepared curves ---------------------------------------------------- */

struct qaws_curve_prepared
{
	qaws_curve_view view;               /* pointers into host memory */
	qaws_allocator const* allocator;
	int owned;                          /* the block is ours to free */
	int rational;                       /* splines: 0 polynomial, 1 weights, 2 homogeneous points */
	unsigned int n, kc, span_count;
	qaws_scalar t0, t1;                 /* the domain */
	qaws_scalar* knots;                 /* splines: kc values, the computing precision */
	qaws_scalar* bounds;                /* span_count + 1 span bounds over the domain */
};

#define PR_ALIGN(x) (((x) + 15) & ~(size_t)15)

/* 16-byte aligned blocks from an allocator: one extra block, the object at the first aligned address, the raw pointer before it */
static void* vw_alloc16(qaws_allocator const* allocator, size_t bytes)
{
	void* raw = qaws_internal_alloc(allocator, (unsigned long)(bytes + 16 + sizeof(void*)));
	uintptr_t at;
	if (!raw)
		return NULL;
	at = ((uintptr_t)raw + sizeof(void*) + 15) & ~(uintptr_t)15;
	((void**)at)[-1] = raw;
	return (void*)at;
}

static void vw_free16(qaws_allocator const* allocator, void* p)
{
	qaws_internal_dealloc(allocator, ((void**)p)[-1]);
}

/* the prepared block's knots (kc, splines only) and span bounds (nb) for a valid view */
static void pr_counts(qaws_curve_view const* v, unsigned int* kc, unsigned int* nb)
{
	unsigned int n = v->points.count, d = v->degree;
	*kc = v->family == QAWS_VIEW_SPLINE ? n + d + 1 : 0;
	*nb = v->family == QAWS_VIEW_SPLINE ? n - d + 1 : v->family == QAWS_VIEW_HERMITE ? n : v->family == QAWS_VIEW_CATMULL_ROM
		? (v->closed ? n + 1 : n - 2) : 2;
}

static size_t pr_bytes(qaws_curve_view const* v)
{
	unsigned int kc, nb;
	pr_counts(v, &kc, &nb);
	return PR_ALIGN(sizeof(struct qaws_curve_prepared)) + PR_ALIGN(sizeof(qaws_scalar) * kc) + PR_ALIGN(sizeof(qaws_scalar) * nb);
}

size_t qaws_curve_prepared_size(qaws_curve_view const* view)
{
	return qaws_curve_view_validate(view) == QAWS_STATUS_OK ? pr_bytes(view) : 0;
}

/* the object laid out in memory (pr_bytes of it), from a valid view */
static qaws_curve_prepared* pr_fill(qaws_curve_view const* view, void* memory)
{
	qaws_curve_prepared* p = (qaws_curve_prepared*)memory;
	unsigned int n = view->points.count, d = view->degree, k, nb;
	memset(p, 0, sizeof *p);
	p->view = *view;
	p->n = n;
	pr_counts(view, &p->kc, &nb);
	p->knots = (qaws_scalar*)((unsigned char*)memory + PR_ALIGN(sizeof *p));
	p->bounds = (qaws_scalar*)((unsigned char*)p->knots + PR_ALIGN(sizeof(qaws_scalar) * p->kc));
	if (view->family == QAWS_VIEW_SPLINE)
	{
		p->rational = view->weights.count ? 1 : view->points.components == (unsigned int)view->dimension + 1 ? 2 : 0;
		vw_knots(&view->knots, n, d, 0, p->knots);
		p->bounds[0] = p->knots[d];
		for (k = d; k < n; k++)
			if (p->knots[k + 1] > p->knots[k])
				p->bounds[++p->span_count] = p->knots[k + 1];
	}
	else if (view->family == QAWS_VIEW_POLYNOMIAL)
	{
		p->span_count = 1;
		p->bounds[0] = view->knots.count ? (qaws_scalar)vw_get(&view->knots, 0, 0) : QAWS_ZERO;
		p->bounds[1] = view->knots.count ? (qaws_scalar)vw_get(&view->knots, 1, 0) : QAWS_ONE;
	}
	else
	{
		/* Hermite and Catmull-Rom: unit spans 0, 1, ..., as the owned curves */
		p->span_count = nb - 1;
		for (k = 0; k < nb; k++)
			p->bounds[k] = (qaws_scalar)k;
	}
	p->t0 = p->bounds[0];
	p->t1 = p->bounds[p->span_count];
	return p;
}

qaws_status qaws_curve_prepare_into(qaws_curve_view const* view, void* memory, size_t bytes, qaws_curve_prepared** out_prepared)
{
	qaws_status st;
	if (!out_prepared)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_prepared = NULL;
	if ((st = qaws_curve_view_validate(view)) != QAWS_STATUS_OK)
		return st;
	if (!memory || (uintptr_t)memory % 16)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (bytes < pr_bytes(view))
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	*out_prepared = pr_fill(view, memory);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_prepare(qaws_curve_view const* view, qaws_allocator const* allocator, qaws_curve_prepared** out_prepared)
{
	qaws_status st;
	void* memory;
	if (!out_prepared)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_prepared = NULL;
	if ((st = qaws_curve_view_validate(view)) != QAWS_STATUS_OK)
		return st;
	if (!(memory = vw_alloc16(allocator, pr_bytes(view))))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	*out_prepared = pr_fill(view, memory);
	(*out_prepared)->allocator = allocator;
	(*out_prepared)->owned = 1;
	return QAWS_STATUS_OK;
}

void qaws_curve_prepared_destroy(qaws_curve_prepared* prepared)
{
	if (prepared && prepared->owned)
		vw_free16(prepared->allocator, prepared);
}

unsigned int qaws_curve_prepared_dependencies(qaws_curve_prepared const* prepared)
{
	return prepared ? QAWS_DEPENDS_LAYOUT | (prepared->view.knots.count ? QAWS_DEPENDS_KNOTS : 0u) : 0u;
}

qaws_range qaws_curve_prepared_get_range(qaws_curve_prepared const* prepared)
{
	qaws_range r;
	r.min_value = prepared ? prepared->t0 : QAWS_ZERO;
	r.max_value = prepared ? prepared->t1 : QAWS_ZERO;
	return r;
}

#define PR_ORDERS 4

static qaws_vec3 pr_vec3(qaws_array_view const* v, unsigned int i, unsigned int dim)
{
	qaws_vec3 r;
	r.x = (qaws_scalar)vw_get(v, i, 0);
	r.y = (qaws_scalar)vw_get(v, i, 1);
	r.z = dim > 2 ? (qaws_scalar)vw_get(v, i, 2) : QAWS_ZERO;
	return r;
}

static void pr_put(qaws_scalar out[3], qaws_vec3 v)
{
	out[0] = v.x;
	out[1] = v.y;
	out[2] = v.z;
}

/* Hermite and Catmull-Rom span `span` at local u in [0, 1]: its coefficients from the borrowed points (and tangents) by the owned
   curves' own routines (qaws_hermite_prepare_3d, qaws_catmull_rom_prepare_3d on the span's points: the same operations, a zero
   third coordinate in 2D), then the cubic core; D3 = 6 a */
static void pr_eval_cubic(qaws_curve_prepared const* p, unsigned int span, qaws_scalar u, qaws_scalar out[PR_ORDERS][3])
{
	unsigned int dim = (unsigned int)p->view.dimension, n = p->n, idx[4], k, cnt = 4, sc;
	qaws_vec3 P[4], T[2], a, b, c, d;
	qaws_eval_3d e;
	if (p->view.family == QAWS_VIEW_HERMITE)
	{
		idx[0] = span;
		idx[1] = span + 1;
		cnt = 2;
	}
	else if (p->view.closed)
	{
		idx[0] = span ? span - 1 : n - 1;
		idx[1] = span;
		idx[2] = (span + 1) % n;
		idx[3] = (span + 2) % n;
	}
	else
		for (k = 0; k < 4; k++)
			idx[k] = span + k;
	for (k = 0; k < cnt; k++)
		P[k] = pr_vec3(&p->view.points, idx[k], dim);
	if (p->view.family == QAWS_VIEW_HERMITE)
	{
		T[0] = pr_vec3(&p->view.tangents, idx[0], dim);
		T[1] = pr_vec3(&p->view.tangents, idx[1], dim);
		qaws_hermite_prepare_3d(P, T, 2, &a, &b, &c, &d);
	}
	else
		qaws_catmull_rom_prepare_3d(P, 4, p->view.parameterization, 0, &a, &b, &c, &d, &sc);
	e = qaws_cubic_eval_3d(a, b, c, d, u, 7);
	pr_put(out[0], e.position);
	pr_put(out[1], e.d1);
	pr_put(out[2], e.d2);
	out[3][0] = QAWS_LITERAL(6.0) * a.x;
	out[3][1] = QAWS_LITERAL(6.0) * a.y;
	out[3][2] = QAWS_LITERAL(6.0) * a.z;
}

/* sum c_i t^i: the Horner core as the owned polynomial, D3 by its own Horner pass */
static void pr_eval_poly(qaws_curve_prepared const* p, qaws_scalar t, qaws_scalar out[PR_ORDERS][3])
{
	qaws_scalar cf[(QAWS_VIEW_MAX_DEGREE + 1) * 3];
	unsigned int dim = (unsigned int)p->view.dimension, d = p->view.degree, i, c;
	int k;
	qaws_eval_3d e;
	memset(cf, 0, sizeof cf);
	for (i = 0; i <= d; i++)
		for (c = 0; c < dim; c++)
			cf[i * 3 + c] = (qaws_scalar)vw_get(&p->view.points, i, c);
	e = qaws_horner_eval_3d(cf, (int)d, t, 7);
	pr_put(out[0], e.position);
	pr_put(out[1], e.d1);
	pr_put(out[2], e.d2);
	for (c = 0; c < 3; c++)
	{
		qaws_scalar s = QAWS_ZERO;
		for (k = (int)d; k >= 3; k--)
			s = s * t + QAWS_ITOF(k * (k - 1) * (k - 2)) * cf[(unsigned int)k * 3 + c];
		out[3][c] = s;
	}
}

/* position and derivatives 1 to kmax at t (clamped to the domain), rows of out. Splines: the polynomial part as the owned B-spline
   (same span search, same basis, same sums), the rational part by the quotient rule. The other families as theirs. */
static void pr_eval(qaws_curve_prepared const* p, qaws_scalar t, unsigned int kmax, qaws_scalar out[PR_ORDERS][3])
{
	unsigned int d = p->view.degree, dim = (unsigned int)p->view.dimension, span, j, r, c, kb = kmax < d ? kmax : d;
	qaws_scalar ders[PR_ORDERS * (QAWS_VIEW_MAX_DEGREE + 1)], A[PR_ORDERS][3], W[PR_ORDERS];
	if (t < p->t0)
		t = p->t0;
	if (t > p->t1)
		t = p->t1;
	if (p->view.family == QAWS_VIEW_HERMITE || p->view.family == QAWS_VIEW_CATMULL_ROM)
	{
		/* the owned curves' span lookup on unit spans */
		span = t >= p->t1 ? p->span_count - 1 : (unsigned int)floor((double)t);
		pr_eval_cubic(p, span, (t - p->bounds[span]) / (p->bounds[span + 1] - p->bounds[span]), out);
		return;
	}
	if (p->view.family == QAWS_VIEW_POLYNOMIAL)
	{
		pr_eval_poly(p, t, out);
		return;
	}
	span = qaws_internal_find_knot_span(p->knots, p->kc, d, p->n, t);
	qaws_internal_bspline_basis_derivs(p->knots, p->kc, d, span, t, kb, ders);
	memset(A, 0, sizeof A);
	memset(W, 0, sizeof W);
	for (j = 0; j <= d; j++)
	{
		unsigned int i = span - d + j;
		qaws_scalar x[3], w = QAWS_ONE;
		for (c = 0; c < dim; c++)
			x[c] = (qaws_scalar)vw_get(&p->view.points, i, c);
		if (p->rational == 1)
			w = (qaws_scalar)vw_get(&p->view.weights, i, 0);
		else if (p->rational == 2)
			w = (qaws_scalar)vw_get(&p->view.points, i, dim);
		for (r = 0; r <= kb; r++)
		{
			/* weights: (N w) P as the owned NURBS; homogeneous points are already w P */
			qaws_scalar b = ders[r * (d + 1) + j], bw = p->rational == 1 ? b * w : b;
			for (c = 0; c < dim; c++)
				A[r][c] += bw * x[c];
			W[r] += b * w;
		}
	}
	if (!p->rational)
	{
		memcpy(out, A, sizeof A);
		return;
	}
	{
		qaws_scalar inv = QAWS_ONE / W[0];
		for (c = 0; c < dim; c++)
		{
			out[0][c] = A[0][c] * inv;
			out[1][c] = (A[1][c] - out[0][c] * W[1]) * inv;
			out[2][c] = (A[2][c] - 2 * out[1][c] * W[1] - out[0][c] * W[2]) * inv;
			out[3][c] = (A[3][c] - 3 * out[2][c] * W[1] - 3 * out[1][c] * W[2] - out[0][c] * W[3]) * inv;
		}
	}
}

/* span `span` at local u: the cubic families directly (as the owned ones, no round trip through the global parameter), the others at
   the global parameter as the owned B-spline and polynomial map it */
static void pr_eval_span(qaws_curve_prepared const* p, unsigned int span, qaws_scalar u, unsigned int kmax, qaws_scalar out[PR_ORDERS][3])
{
	if (p->view.family == QAWS_VIEW_HERMITE || p->view.family == QAWS_VIEW_CATMULL_ROM)
		pr_eval_cubic(p, span, u, out);
	else
		pr_eval(p, p->bounds[span] + u * (p->bounds[span + 1] - p->bounds[span]), kmax, out);
}

static unsigned int pr_order(unsigned int flags)
{
	return flags & QAWS_EVAL_FLAG_D3 ? 3 : flags & QAWS_EVAL_FLAG_D2 ? 2 : flags & QAWS_EVAL_FLAG_D1 ? 1 : 0;
}

static void pr_store_3d(qaws_scalar const v[PR_ORDERS][3], unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_vec3* dst[PR_ORDERS];
	unsigned int r;
	memset(out_result, 0, sizeof *out_result);
	dst[0] = &out_result->position;
	dst[1] = &out_result->d1;
	dst[2] = &out_result->d2;
	dst[3] = &out_result->d3;
	for (r = 0; r < PR_ORDERS; r++)
		if (eval_flags & (1u << r))
		{
			dst[r]->x = v[r][0];
			dst[r]->y = v[r][1];
			dst[r]->z = v[r][2];
			out_result->valid_flags |= 1u << r;
		}
}

static void pr_store_2d(qaws_scalar const v[PR_ORDERS][3], unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_vec2* dst[PR_ORDERS];
	unsigned int r;
	memset(out_result, 0, sizeof *out_result);
	dst[0] = &out_result->position;
	dst[1] = &out_result->d1;
	dst[2] = &out_result->d2;
	dst[3] = &out_result->d3;
	for (r = 0; r < PR_ORDERS; r++)
		if (eval_flags & (1u << r))
		{
			dst[r]->x = v[r][0];
			dst[r]->y = v[r][1];
			out_result->valid_flags |= 1u << r;
		}
}

qaws_status qaws_curve_prepared_evaluate_3d(qaws_curve_prepared const* prepared, qaws_scalar parameter, unsigned int eval_flags,
	qaws_eval_result_3d* out_result)
{
	qaws_scalar v[PR_ORDERS][3];
	if (!prepared || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (prepared->view.dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	pr_eval(prepared, parameter, pr_order(eval_flags), v);
	pr_store_3d(v, eval_flags, out_result);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_prepared_evaluate_2d(qaws_curve_prepared const* prepared, qaws_scalar parameter, unsigned int eval_flags,
	qaws_eval_result_2d* out_result)
{
	qaws_scalar v[PR_ORDERS][3];
	if (!prepared || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (prepared->view.dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	pr_eval(prepared, parameter, pr_order(eval_flags), v);
	pr_store_2d(v, eval_flags, out_result);
	return QAWS_STATUS_OK;
}

static int pr_present(qaws_array_view_mut const* o)
{
	return o->data || o->planes[0];
}

qaws_status qaws_curve_prepared_evaluate_range(qaws_curve_prepared const* prepared, qaws_array_view const* parameters, unsigned int begin,
	unsigned int end, qaws_curve_outputs const* outputs)
{
	qaws_array_view_mut const* o[PR_ORDERS];
	unsigned int r, i, c, dim, kmax = 0, any = 0;
	qaws_status st;
	if (!prepared || !parameters || !outputs || begin > end)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = qaws_array_view_validate(parameters)) != QAWS_STATUS_OK)
		return st;
	if (parameters->components != 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (parameters->count < end)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	dim = (unsigned int)prepared->view.dimension;
	o[0] = &outputs->position;
	o[1] = &outputs->d1;
	o[2] = &outputs->d2;
	o[3] = &outputs->d3;
	for (r = 0; r < PR_ORDERS; r++)
		if (pr_present(o[r]))
		{
			qaws_array_view v = vw_const(o[r]);
			if ((st = qaws_array_view_validate(&v)) != QAWS_STATUS_OK)
				return st;
			if (o[r]->components != dim)
				return QAWS_STATUS_INVALID_DIMENSION;
			if (o[r]->count < end)
				return QAWS_STATUS_BUFFER_TOO_SMALL;
			kmax = r;
			any = 1;
		}
	if (!any)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = begin; i < end; i++)
	{
		qaws_scalar v[PR_ORDERS][3];
		pr_eval(prepared, (qaws_scalar)vw_get(parameters, i, 0), kmax, v);
		for (r = 0; r <= kmax; r++)
			if (pr_present(o[r]))
				for (c = 0; c < dim; c++)
					vw_put(o[r], i, c, (double)v[r][c]);
	}
	return QAWS_STATUS_OK;
}

/* ---- borrowed qaws_curve ------------------------------------------------ */

static qaws_status bw_eval_span_2d(qaws_curve const* curve, unsigned int span, qaws_scalar local_t, unsigned int flags, qaws_eval_result_2d* out)
{
	qaws_scalar v[PR_ORDERS][3];
	pr_eval_span((qaws_curve_prepared const*)curve->impl, span, local_t, pr_order(flags), v);
	pr_store_2d(v, flags, out);
	return QAWS_STATUS_OK;
}

static qaws_status bw_eval_span_3d(qaws_curve const* curve, unsigned int span, qaws_scalar local_t, unsigned int flags, qaws_eval_result_3d* out)
{
	qaws_scalar v[PR_ORDERS][3];
	pr_eval_span((qaws_curve_prepared const*)curve->impl, span, local_t, pr_order(flags), v);
	pr_store_3d(v, flags, out);
	return QAWS_STATUS_OK;
}

static void bw_destroy(void* impl, qaws_allocator const* allocator)
{
	(void)allocator;
	qaws_curve_prepared_destroy((qaws_curve_prepared*)impl);
}

static int bw_is_closed(qaws_curve const* curve)
{
	qaws_curve_prepared const* p = (qaws_curve_prepared const*)curve->impl;
	return (p->view.family == QAWS_VIEW_CATMULL_ROM && p->view.closed) || qaws_internal_curve_ends_meet(curve);
}

static int bw_is_periodic(qaws_curve const* curve)
{
	qaws_curve_prepared const* p = (qaws_curve_prepared const*)curve->impl;
	return p->view.family == QAWS_VIEW_CATMULL_ROM && p->view.closed;
}

static int bw_is_rational(qaws_curve const* curve)
{
	return ((qaws_curve_prepared const*)curve->impl)->rational != 0;
}

/* splines: the degree less the largest multiplicity of a knot inside the domain; the others as the owned families */
static qaws_continuity bw_get_continuity(qaws_curve const* curve)
{
	qaws_curve_prepared const* p = (qaws_curve_prepared const*)curve->impl;
	unsigned int d = p->view.degree, k, run = 1, worst = 0;
	int smooth;
	if (p->view.family == QAWS_VIEW_POLYNOMIAL)
		return QAWS_CONTINUITY_C3;
	if (p->view.family != QAWS_VIEW_SPLINE)
		return QAWS_CONTINUITY_C1;
	for (k = d + 1; k < p->n; k++)
		if (p->knots[k] > p->knots[d] && p->knots[k] < p->knots[p->n])
		{
			run = p->knots[k] == p->knots[k - 1] ? run + 1 : 1;
			worst = run > worst ? run : worst;
		}
	smooth = (int)d - (int)worst;
	return smooth <= 0 ? QAWS_CONTINUITY_C0 : smooth == 1 ? QAWS_CONTINUITY_C1 : smooth == 2 ? QAWS_CONTINUITY_C2 : QAWS_CONTINUITY_C3;
}

static qaws_curve_vtable const bw_vtable = {
	bw_eval_span_2d, bw_eval_span_3d, bw_destroy, bw_is_closed, bw_is_periodic, bw_is_rational, bw_get_continuity, NULL
};

qaws_status qaws_curve_create_borrowed(qaws_curve_view const* view, qaws_allocator const* allocator, qaws_curve** out_curve)
{
	qaws_curve_prepared* p;
	qaws_curve* curve;
	qaws_range range;
	qaws_status st;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if ((st = qaws_curve_prepare(view, allocator, &p)) != QAWS_STATUS_OK)
		return st;
	range = qaws_curve_prepared_get_range(p);
	curve = qaws_internal_curve_alloc_ex(QAWS_CURVE_KIND_VIEW, view->dimension, view->degree, p->span_count, range, &bw_vtable, allocator);
	if (!curve)
	{
		qaws_curve_prepared_destroy(p);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(curve->span_boundaries, p->bounds, sizeof(qaws_scalar) * ((size_t)p->span_count + 1));
	curve->impl = p;
	*out_curve = curve;
	return QAWS_STATUS_OK;
}

int qaws_internal_view_span_window(qaws_curve const* curve, unsigned int s, unsigned int* degree, double* window, double* knots)
{
	qaws_curve_prepared const* p = (qaws_curve_prepared const*)curve->impl;
	unsigned int d, dim, k, j, c;
	if (curve->kind != QAWS_CURVE_KIND_VIEW || p->view.family != QAWS_VIEW_SPLINE)
		return 0;
	d = p->view.degree;
	dim = (unsigned int)p->view.dimension;
	k = qaws_internal_find_knot_span(p->knots, p->kc, d, p->n, (p->bounds[s] + p->bounds[s + 1]) / 2);
	for (j = 0; j <= d; j++)
	{
		unsigned int i = k - d + j;
		double w = p->rational == 1 ? vw_get(&p->view.weights, i, 0) : p->rational == 2 ? vw_get(&p->view.points, i, dim) : 1.0;
		/* homogeneous points are already w P */
		for (c = 0; c < dim; c++)
			window[j * (dim + 1) + c] = (p->rational == 1 ? w : 1.0) * vw_get(&p->view.points, i, c);
		window[j * (dim + 1) + dim] = w;
	}
	for (j = 0; j < 2 * d; j++)
		knots[j] = p->knots[k - d + 1 + j];
	*degree = d;
	return 1;
}

/* ---- owned copies ------------------------------------------------------- */

/* the owned curve of the view's family: splines a B-spline (polynomial) or NURBS curve with the prepared knots (homogeneous points
   divided by w); Hermite with its tangents; Catmull-Rom; polynomial over the prepared range */
static qaws_status pr_copy(qaws_curve_prepared const* p, qaws_allocator const* allocator, qaws_curve** out_curve)
{
	unsigned int dim = (unsigned int)p->view.dimension, i, c;
	qaws_scalar* pts = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)(sizeof(qaws_scalar) * (size_t)p->n * (2 * dim + 1)));
	qaws_scalar *wts, *tans;
	qaws_status st;
	if (!pts)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	wts = pts + (size_t)p->n * dim;
	tans = wts + p->n;
	for (i = 0; i < p->n; i++)
	{
		wts[i] = p->rational == 1 ? (qaws_scalar)vw_get(&p->view.weights, i, 0)
			: p->rational == 2 ? (qaws_scalar)vw_get(&p->view.points, i, dim) : QAWS_ONE;
		for (c = 0; c < dim; c++)
		{
			pts[i * dim + c] = (qaws_scalar)vw_get(&p->view.points, i, c) / (p->rational == 2 ? wts[i] : QAWS_ONE);
			tans[i * dim + c] = p->view.family == QAWS_VIEW_HERMITE ? (qaws_scalar)vw_get(&p->view.tangents, i, c) : QAWS_ZERO;
		}
	}
	if (p->view.family == QAWS_VIEW_HERMITE)
	{
		qaws_hermite_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = p->view.dimension;
		d.degree = 3;
		d.points = pts;
		d.derivatives = tans;
		d.point_count = p->n;
		d.derivative_count = p->n;
		st = qaws_curve_create_hermite_ex(&d, allocator, out_curve);
	}
	else if (p->view.family == QAWS_VIEW_CATMULL_ROM)
	{
		qaws_catmull_rom_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = p->view.dimension;
		d.control_points = pts;
		d.control_point_count = p->n;
		d.parameterization = p->view.parameterization;
		d.closed = p->view.closed;
		st = qaws_curve_create_catmull_rom_ex(&d, allocator, out_curve);
	}
	else if (p->view.family == QAWS_VIEW_POLYNOMIAL)
	{
		qaws_polynomial_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = p->view.dimension;
		d.degree = p->view.degree;
		d.coefficients = pts;
		d.coefficient_count = p->n;
		d.t_min = p->t0;
		d.t_max = p->t1;
		st = qaws_curve_create_polynomial_ex(&d, allocator, out_curve);
	}
	else if (p->rational)
	{
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = p->view.dimension;
		d.degree = p->view.degree;
		d.control_points = pts;
		d.control_point_count = p->n;
		d.knots = p->knots;
		d.knot_count = p->kc;
		d.weights = wts;
		d.weight_count = p->n;
		st = qaws_curve_create_nurbs_ex(&d, allocator, out_curve);
	}
	else
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = p->view.dimension;
		d.degree = p->view.degree;
		d.control_points = pts;
		d.control_point_count = p->n;
		d.knots = p->knots;
		d.knot_count = p->kc;
		st = qaws_curve_create_bspline_ex(&d, allocator, out_curve);
	}
	qaws_internal_dealloc(allocator, pts);
	return st;
}

qaws_status qaws_curve_create_from_view(qaws_curve_view const* view, qaws_allocator const* allocator, qaws_curve** out_curve)
{
	qaws_curve_prepared* p;
	qaws_status st;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if ((st = qaws_curve_prepare(view, allocator, &p)) != QAWS_STATUS_OK)
		return st;
	st = pr_copy(p, allocator, out_curve);
	qaws_curve_prepared_destroy(p);
	return st;
}

qaws_status qaws_internal_curve_view_copy(qaws_curve const* curve, qaws_curve** out_curve)
{
	return pr_copy((qaws_curve_prepared const*)curve->impl, curve->allocator, out_curve);
}

/* ---- surfaces ----------------------------------------------------------- */

qaws_status qaws_surface_view_validate(qaws_surface_view const* view)
{
	qaws_status st;
	if (!view)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (view->u_degree < 1 || view->u_degree > QAWS_VIEW_MAX_DEGREE || view->v_degree < 1 || view->v_degree > QAWS_VIEW_MAX_DEGREE)
		return QAWS_STATUS_INVALID_DEGREE;
	if (view->u_count < view->u_degree + 1 || view->v_count < view->v_degree + 1 || view->u_count > 0xffffffffu / view->v_count)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if ((st = vw_check_points(&view->points, &view->weights, 3, view->u_count * view->v_count)) != QAWS_STATUS_OK)
		return st;
	if ((st = vw_check_knots(&view->u_knots, view->u_count, view->u_degree)) != QAWS_STATUS_OK)
		return st;
	return vw_check_knots(&view->v_knots, view->v_count, view->v_degree);
}

struct qaws_surface_prepared
{
	qaws_surface_view view;             /* pointers into host memory */
	qaws_allocator const* allocator;
	int owned;
	int rational;                       /* 0: polynomial, 1: weights, 2: homogeneous points */
	qaws_scalar* u_knots;               /* u_count + u_degree + 1, the computing precision */
	qaws_scalar* v_knots;
};

static size_t sp_bytes(qaws_surface_view const* v)
{
	return PR_ALIGN(sizeof(struct qaws_surface_prepared)) + PR_ALIGN(sizeof(qaws_scalar) * ((size_t)v->u_count + v->u_degree + 1))
		+ PR_ALIGN(sizeof(qaws_scalar) * ((size_t)v->v_count + v->v_degree + 1));
}

size_t qaws_surface_prepared_size(qaws_surface_view const* view)
{
	return qaws_surface_view_validate(view) == QAWS_STATUS_OK ? sp_bytes(view) : 0;
}

static qaws_surface_prepared* sp_fill(qaws_surface_view const* view, void* memory)
{
	qaws_surface_prepared* p = (qaws_surface_prepared*)memory;
	memset(p, 0, sizeof *p);
	p->view = *view;
	p->rational = view->weights.count ? 1 : view->points.components == 4 ? 2 : 0;
	p->u_knots = (qaws_scalar*)((unsigned char*)memory + PR_ALIGN(sizeof *p));
	p->v_knots = (qaws_scalar*)((unsigned char*)p->u_knots + PR_ALIGN(sizeof(qaws_scalar) * ((size_t)view->u_count + view->u_degree + 1)));
	vw_knots(&view->u_knots, view->u_count, view->u_degree, 1, p->u_knots);
	vw_knots(&view->v_knots, view->v_count, view->v_degree, 1, p->v_knots);
	return p;
}

qaws_status qaws_surface_prepare_into(qaws_surface_view const* view, void* memory, size_t bytes, qaws_surface_prepared** out_prepared)
{
	qaws_status st;
	if (!out_prepared)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_prepared = NULL;
	if ((st = qaws_surface_view_validate(view)) != QAWS_STATUS_OK)
		return st;
	if (!memory || (uintptr_t)memory % 16)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (bytes < sp_bytes(view))
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	*out_prepared = sp_fill(view, memory);
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_prepare(qaws_surface_view const* view, qaws_allocator const* allocator, qaws_surface_prepared** out_prepared)
{
	qaws_status st;
	void* memory;
	if (!out_prepared)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_prepared = NULL;
	if ((st = qaws_surface_view_validate(view)) != QAWS_STATUS_OK)
		return st;
	if (!(memory = vw_alloc16(allocator, sp_bytes(view))))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	*out_prepared = sp_fill(view, memory);
	(*out_prepared)->allocator = allocator;
	(*out_prepared)->owned = 1;
	return QAWS_STATUS_OK;
}

void qaws_surface_prepared_destroy(qaws_surface_prepared* prepared)
{
	if (prepared && prepared->owned)
		vw_free16(prepared->allocator, prepared);
}

unsigned int qaws_surface_prepared_dependencies(qaws_surface_prepared const* prepared)
{
	return prepared ? QAWS_DEPENDS_LAYOUT | (prepared->view.u_knots.count || prepared->view.v_knots.count ? QAWS_DEPENDS_KNOTS : 0u) : 0u;
}

qaws_range qaws_surface_prepared_get_u_range(qaws_surface_prepared const* prepared)
{
	qaws_range r;
	r.min_value = prepared ? prepared->u_knots[prepared->view.u_degree] : QAWS_ZERO;
	r.max_value = prepared ? prepared->u_knots[prepared->view.u_count] : QAWS_ZERO;
	return r;
}

qaws_range qaws_surface_prepared_get_v_range(qaws_surface_prepared const* prepared)
{
	qaws_range r;
	r.min_value = prepared ? prepared->v_knots[prepared->view.v_degree] : QAWS_ZERO;
	r.max_value = prepared ? prepared->v_knots[prepared->view.v_count] : QAWS_ZERO;
	return r;
}

#define SP_ORDERS 4

/* S[k][l]: the partial derivative of orders k in u and l in v, for k <= mu, l <= mv, k + l <= top <= 3 (others zero); the sums as
   the owned B-spline and NURBS surfaces (same spans, basis, products and order of terms), the rational part by the quotient rule
   (Piegl and Tiller A4.4, the same operations as the owned NURBS surface up to second order) */
static void sp_eval(qaws_surface_prepared const* p, qaws_scalar u, qaws_scalar v, unsigned int mu, unsigned int mv, unsigned int top,
	qaws_scalar S[SP_ORDERS][SP_ORDERS][3])
{
	static qaws_scalar const bin[SP_ORDERS][SP_ORDERS] = { { 1, 0, 0, 0 }, { 1, 1, 0, 0 }, { 1, 2, 1, 0 }, { 1, 3, 3, 1 } };
	unsigned int ud = p->view.u_degree, vd = p->view.v_degree, uc = p->view.u_count, vc = p->view.v_count;
	unsigned int bu = mu < ud ? mu : ud, bv = mv < vd ? mv : vd, us, vs, i, j, ku, kv, c;
	qaws_scalar uders[SP_ORDERS * (QAWS_VIEW_MAX_DEGREE + 1)], vders[SP_ORDERS * (QAWS_VIEW_MAX_DEGREE + 1)], W[SP_ORDERS][SP_ORDERS];
	u = u < p->u_knots[ud] ? p->u_knots[ud] : u > p->u_knots[uc] ? p->u_knots[uc] : u;
	v = v < p->v_knots[vd] ? p->v_knots[vd] : v > p->v_knots[vc] ? p->v_knots[vc] : v;
	us = qaws_internal_find_knot_span(p->u_knots, uc + ud + 1, ud, uc, u);
	vs = qaws_internal_find_knot_span(p->v_knots, vc + vd + 1, vd, vc, v);
	qaws_internal_bspline_basis_derivs(p->u_knots, uc + ud + 1, ud, us, u, bu, uders);
	qaws_internal_bspline_basis_derivs(p->v_knots, vc + vd + 1, vd, vs, v, bv, vders);
	memset(S, 0, sizeof(qaws_scalar) * SP_ORDERS * SP_ORDERS * 3);
	memset(W, 0, sizeof W);
	for (i = 0; i <= ud; i++)
		for (j = 0; j <= vd; j++)
		{
			unsigned int e = (us - ud + i) * vc + (vs - vd + j);
			qaws_scalar x[3], w = QAWS_ONE;
			for (c = 0; c < 3; c++)
				x[c] = (qaws_scalar)vw_get(&p->view.points, e, c);
			if (p->rational == 1)
				w = (qaws_scalar)vw_get(&p->view.weights, e, 0);
			else if (p->rational == 2)
				w = (qaws_scalar)vw_get(&p->view.points, e, 3);
			for (ku = 0; ku <= bu; ku++)
				for (kv = 0; kv <= bv && ku + kv <= top; kv++)
				{
					/* weights: (Nu Nv w) P as the owned NURBS surface; homogeneous points are already w P */
					qaws_scalar b = uders[ku * (ud + 1) + i] * vders[kv * (vd + 1) + j], bw = p->rational == 1 ? b * w : b;
					W[ku][kv] += b * w;
					for (c = 0; c < 3; c++)
						S[ku][kv][c] += bw * x[c];
				}
		}
	if (p->rational)
	{
		qaws_scalar inv = QAWS_ONE / W[0][0];
		unsigned int k, l;
		/* in order of k then l: every term read is already a quotient */
		for (k = 0; k <= mu; k++)
			for (l = 0; l <= mv && k + l <= top; l++)
				for (c = 0; c < 3; c++)
				{
					qaws_scalar x = S[k][l][c];
					for (j = 1; j <= l; j++)
						x -= bin[l][j] * W[0][j] * S[k][l - j][c];
					for (i = 1; i <= k; i++)
					{
						qaws_scalar y = QAWS_ZERO;
						x -= bin[k][i] * W[i][0] * S[k - i][l][c];
						for (j = 1; j <= l; j++)
							y += bin[l][j] * W[i][j] * S[k - i][l - j][c];
						x -= bin[k][i] * y;
					}
					S[k][l][c] = x * inv;
				}
	}
}

/* the orders sp_eval needs for the partials in `bits` of a table of (k, l) orders */
static void sp_orders(unsigned int bits, unsigned char const (*orders)[2], unsigned int n, unsigned int* mu, unsigned int* mv, unsigned int* top)
{
	unsigned int b;
	*mu = *mv = *top = 0;
	for (b = 0; b < n; b++)
		if (bits & (1u << b))
		{
			unsigned int sum = (unsigned int)orders[b][0] + orders[b][1];
			*mu = orders[b][0] > *mu ? orders[b][0] : *mu;
			*mv = orders[b][1] > *mv ? orders[b][1] : *mv;
			*top = sum > *top ? sum : *top;
		}
}

/* QAWS_SURFACE_EVAL_* bit b: the (k, l) orders of its partial derivative (bit 3, the normal, is made from du and dv) */
static unsigned char const sp_order[7][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 0, 0 }, { 2, 0 }, { 1, 1 }, { 0, 2 } };

/* qaws_surface_jet channel b (QAWS_SJET_*): its orders */
static unsigned char const sp_jet[QAWS_SURFACE_JET_COUNT][2] = {
	{ 0, 0 }, { 1, 0 }, { 0, 1 }, { 2, 0 }, { 1, 1 }, { 0, 2 }, { 3, 0 }, { 2, 1 }, { 1, 2 }, { 0, 3 }
};

qaws_status qaws_surface_prepared_evaluate(qaws_surface_prepared const* prepared, qaws_scalar u, qaws_scalar v, unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_scalar S[SP_ORDERS][SP_ORDERS][3];
	qaws_vec3* dst[7];
	unsigned int mu, mv, top, b;
	if (!prepared || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(out_result, 0, sizeof *out_result);
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
		eval_flags |= QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV;
	sp_orders(eval_flags, sp_order, 7, &mu, &mv, &top);
	sp_eval(prepared, u, v, mu, mv, top, S);
	dst[0] = &out_result->position;
	dst[1] = &out_result->du;
	dst[2] = &out_result->dv;
	dst[3] = &out_result->normal;
	dst[4] = &out_result->duu;
	dst[5] = &out_result->duv;
	dst[6] = &out_result->dvv;
	for (b = 0; b < 7; b++)
		if (b != 3 && (eval_flags & (1u << b)))
		{
			qaws_scalar const* s = S[sp_order[b][0]][sp_order[b][1]];
			dst[b]->x = s[0];
			dst[b]->y = s[1];
			dst[b]->z = s[2];
			out_result->valid_flags |= 1u << b;
		}
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_prepared_evaluate_range(qaws_surface_prepared const* prepared, qaws_array_view const* parameters,
	unsigned int begin, unsigned int end, qaws_surface_outputs const* outputs)
{
	qaws_array_view_mut const* o[7];
	unsigned int b, i, flags = 0;
	qaws_status st;
	if (!prepared || !parameters || !outputs || begin > end)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = qaws_array_view_validate(parameters)) != QAWS_STATUS_OK)
		return st;
	if (parameters->components != 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (parameters->count < end)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	o[0] = &outputs->position;
	o[1] = &outputs->du;
	o[2] = &outputs->dv;
	o[3] = &outputs->normal;
	o[4] = &outputs->duu;
	o[5] = &outputs->duv;
	o[6] = &outputs->dvv;
	for (b = 0; b < 7; b++)
		if (pr_present(o[b]))
		{
			qaws_array_view v = vw_const(o[b]);
			if ((st = qaws_array_view_validate(&v)) != QAWS_STATUS_OK)
				return st;
			if (o[b]->components != 3)
				return QAWS_STATUS_INVALID_DIMENSION;
			if (o[b]->count < end)
				return QAWS_STATUS_BUFFER_TOO_SMALL;
			flags |= 1u << b;
		}
	if (!flags)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = begin; i < end; i++)
	{
		qaws_surface_eval_result r;
		qaws_vec3 const* src[7];
		qaws_surface_prepared_evaluate(prepared, (qaws_scalar)vw_get(parameters, i, 0), (qaws_scalar)vw_get(parameters, i, 1), flags, &r);
		src[0] = &r.position;
		src[1] = &r.du;
		src[2] = &r.dv;
		src[3] = &r.normal;
		src[4] = &r.duu;
		src[5] = &r.duv;
		src[6] = &r.dvv;
		for (b = 0; b < 7; b++)
			if (flags & (1u << b))
			{
				vw_put(o[b], i, 0, (double)src[b]->x);
				vw_put(o[b], i, 1, (double)src[b]->y);
				vw_put(o[b], i, 2, (double)src[b]->z);
			}
	}
	return QAWS_STATUS_OK;
}

static qaws_status bs_eval(qaws_surface const* surface, qaws_scalar u, qaws_scalar v, unsigned int flags, qaws_surface_eval_result* out)
{
	return qaws_surface_prepared_evaluate((qaws_surface_prepared const*)surface->impl, u, v, flags, out);
}

static void bs_destroy(void* impl, qaws_allocator const* allocator)
{
	(void)allocator;
	qaws_surface_prepared_destroy((qaws_surface_prepared*)impl);
}

static int bs_is_rational(qaws_surface const* surface)
{
	return ((qaws_surface_prepared const*)surface->impl)->rational != 0;
}

/* the primal jet (closest points, curvatures, offsets of it ...); no derivatives with respect to the host's arrays */
static qaws_status bs_eval_jet(qaws_surface const* surface, qaws_scalar u, qaws_scalar v, unsigned int channels, qaws_surface_jet* out)
{
	qaws_scalar S[SP_ORDERS][SP_ORDERS][3];
	unsigned int mu, mv, top, b;
	sp_orders(channels, sp_jet, QAWS_SURFACE_JET_COUNT, &mu, &mv, &top);
	sp_eval((qaws_surface_prepared const*)surface->impl, u, v, mu, mv, top, S);
	memset(out, 0, sizeof *out);
	for (b = 0; b < QAWS_SURFACE_JET_COUNT; b++)
		if (channels & (1u << b))
		{
			qaws_scalar const* s = S[sp_jet[b][0]][sp_jet[b][1]];
			out->d[b].x = s[0];
			out->d[b].y = s[1];
			out->d[b].z = s[2];
		}
	out->channels = channels;
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const bs_diff = { 0, QAWS_DIFF_PIECEWISE_SMOOTH, NULL, NULL, NULL, bs_eval_jet, NULL, NULL, NULL };

static qaws_surface_vtable const bs_vtable = { bs_eval, bs_destroy, bs_is_rational, &bs_diff };

qaws_status qaws_surface_create_borrowed(qaws_surface_view const* view, qaws_allocator const* allocator, qaws_surface** out_surface)
{
	qaws_surface_prepared* p;
	qaws_surface* s;
	qaws_status st;
	if (!out_surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_surface = NULL;
	if ((st = qaws_surface_prepare(view, allocator, &p)) != QAWS_STATUS_OK)
		return st;
	s = qaws_internal_surface_alloc_ex(QAWS_SURFACE_KIND_VIEW, view->u_degree, view->v_degree, qaws_surface_prepared_get_u_range(p),
		qaws_surface_prepared_get_v_range(p), &bs_vtable, allocator);
	if (!s)
	{
		qaws_surface_prepared_destroy(p);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	s->impl = p;
	*out_surface = s;
	return QAWS_STATUS_OK;
}

int qaws_internal_surface_view_net(qaws_surface const* surface, unsigned int* degrees, unsigned int* counts, qaws_scalar const** u_knots,
	qaws_scalar const** v_knots)
{
	qaws_surface_prepared const* p = (qaws_surface_prepared const*)surface->impl;
	if (surface->kind != QAWS_SURFACE_KIND_VIEW)
		return 0;
	degrees[0] = p->view.u_degree;
	degrees[1] = p->view.v_degree;
	counts[0] = p->view.u_count;
	counts[1] = p->view.v_count;
	*u_knots = p->u_knots;
	*v_knots = p->v_knots;
	return 1;
}

void qaws_internal_surface_view_point(qaws_surface const* surface, unsigned int i, unsigned int j, double* h)
{
	qaws_surface_prepared const* p = (qaws_surface_prepared const*)surface->impl;
	unsigned int e = i * p->view.v_count + j, c;
	double w = p->rational == 1 ? vw_get(&p->view.weights, e, 0) : p->rational == 2 ? vw_get(&p->view.points, e, 3) : 1.0;
	/* homogeneous points are already w P */
	for (c = 0; c < 3; c++)
		h[c] = (p->rational == 1 ? w : 1.0) * vw_get(&p->view.points, e, c);
	h[3] = w;
}

qaws_status qaws_surface_create_from_view(qaws_surface_view const* view, qaws_allocator const* allocator, qaws_surface** out_surface)
{
	qaws_surface_prepared* p;
	qaws_vec3* pts;
	qaws_scalar* wts;
	unsigned int n, i;
	qaws_status st;
	if (!out_surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_surface = NULL;
	if ((st = qaws_surface_prepare(view, allocator, &p)) != QAWS_STATUS_OK)
		return st;
	n = view->u_count * view->v_count;
	pts = (qaws_vec3*)qaws_internal_alloc(allocator, (unsigned long)((sizeof(qaws_vec3) + sizeof(qaws_scalar)) * (size_t)n));
	if (!pts)
	{
		qaws_surface_prepared_destroy(p);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	wts = (qaws_scalar*)(pts + n);
	for (i = 0; i < n; i++)
	{
		qaws_scalar dw;
		wts[i] = p->rational == 1 ? (qaws_scalar)vw_get(&view->weights, i, 0) : p->rational == 2 ? (qaws_scalar)vw_get(&view->points, i, 3) : QAWS_ONE;
		dw = p->rational == 2 ? wts[i] : QAWS_ONE;
		pts[i].x = (qaws_scalar)vw_get(&view->points, i, 0) / dw;
		pts[i].y = (qaws_scalar)vw_get(&view->points, i, 1) / dw;
		pts[i].z = (qaws_scalar)vw_get(&view->points, i, 2) / dw;
	}
	if (p->rational)
	{
		qaws_surface_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = view->u_degree;
		d.v_degree = view->v_degree;
		d.control_points = pts;
		d.u_point_count = view->u_count;
		d.v_point_count = view->v_count;
		d.weights = wts;
		d.u_knots = p->u_knots;
		d.u_knot_count = view->u_count + view->u_degree + 1;
		d.v_knots = p->v_knots;
		d.v_knot_count = view->v_count + view->v_degree + 1;
		st = qaws_surface_create_nurbs_ex(&d, allocator, out_surface);
	}
	else
	{
		qaws_surface_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = view->u_degree;
		d.v_degree = view->v_degree;
		d.control_points = pts;
		d.u_point_count = view->u_count;
		d.v_point_count = view->v_count;
		d.u_knots = p->u_knots;
		d.u_knot_count = view->u_count + view->u_degree + 1;
		d.v_knots = p->v_knots;
		d.v_knot_count = view->v_count + view->v_degree + 1;
		st = qaws_surface_create_bspline_ex(&d, allocator, out_surface);
	}
	qaws_internal_dealloc(allocator, pts);
	qaws_surface_prepared_destroy(p);
	return st;
}
