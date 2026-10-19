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

qaws_status qaws_curve_view_validate(qaws_curve_view const* view)
{
	unsigned int n, d, dim, i, c;
	qaws_status st;
	if (!view)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (view->dimension != QAWS_DIMENSION_2D && view->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	dim = (unsigned int)view->dimension;
	d = view->degree;
	if (d < 1 || d > QAWS_VIEW_MAX_DEGREE)
		return QAWS_STATUS_INVALID_DEGREE;
	if ((st = qaws_array_view_validate(&view->points)) != QAWS_STATUS_OK)
		return st;
	if (view->points.components != dim && view->points.components != dim + 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = view->points.count;
	if (n < d + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (view->weights.count)
	{
		if (view->points.components != dim || view->weights.count != n || view->weights.components != 1)
			return QAWS_STATUS_INVALID_WEIGHT_COUNT;
		if ((st = qaws_array_view_validate(&view->weights)) != QAWS_STATUS_OK)
			return st;
	}
	for (i = 0; i < n; i++)
	{
		double w = view->weights.count ? vw_get(&view->weights, i, 0) : view->points.components == dim + 1 ? vw_get(&view->points, i, dim) : 1;
		if (!(w > 0) || !isfinite(w))
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (c = 0; c < view->points.components; c++)
			if (!isfinite(vw_get(&view->points, i, c)))
				return QAWS_STATUS_INVALID_ARGUMENT;
	}
	if (view->knots.count)
	{
		if (view->knots.count != n + d + 1 || view->knots.components != 1)
			return QAWS_STATUS_INVALID_KNOT_COUNT;
		if ((st = qaws_array_view_validate(&view->knots)) != QAWS_STATUS_OK)
			return st;
		for (i = 0; i < view->knots.count; i++)
		{
			double k = vw_get(&view->knots, i, 0);
			if (!isfinite(k) || (i && k < vw_get(&view->knots, i - 1, 0)))
				return QAWS_STATUS_INVALID_KNOT_VECTOR;
		}
		if (!(vw_get(&view->knots, d, 0) < vw_get(&view->knots, n, 0)))
			return QAWS_STATUS_INVALID_KNOT_VECTOR;
	}
	return QAWS_STATUS_OK;
}

/* ---- prepared curves ---------------------------------------------------- */

struct qaws_curve_prepared
{
	qaws_curve_view view;               /* pointers into host memory */
	qaws_allocator const* allocator;
	int owned;                          /* the block is ours to free */
	int rational;                       /* 0: polynomial, 1: weights, 2: homogeneous points */
	unsigned int n, kc, span_count;
	qaws_scalar* knots;                 /* kc values, the computing precision */
	qaws_scalar* bounds;                /* span_count + 1 distinct knot values over the domain */
};

#define PR_ALIGN(x) (((x) + 15) & ~(size_t)15)

static size_t pr_bytes(unsigned int n, unsigned int d)
{
	return PR_ALIGN(sizeof(struct qaws_curve_prepared)) + PR_ALIGN(sizeof(qaws_scalar) * ((size_t)n + d + 1))
		+ PR_ALIGN(sizeof(qaws_scalar) * ((size_t)n - d + 1));
}

size_t qaws_curve_prepared_size(qaws_curve_view const* view)
{
	return qaws_curve_view_validate(view) == QAWS_STATUS_OK ? pr_bytes(view->points.count, view->degree) : 0;
}

/* the object laid out in memory (pr_bytes of it), from a valid view */
static qaws_curve_prepared* pr_fill(qaws_curve_view const* view, void* memory)
{
	qaws_curve_prepared* p = (qaws_curve_prepared*)memory;
	unsigned int n = view->points.count, d = view->degree, k;
	memset(p, 0, sizeof *p);
	p->view = *view;
	p->n = n;
	p->kc = n + d + 1;
	p->rational = view->weights.count ? 1 : view->points.components == (unsigned int)view->dimension + 1 ? 2 : 0;
	p->knots = (qaws_scalar*)((unsigned char*)memory + PR_ALIGN(sizeof *p));
	p->bounds = (qaws_scalar*)((unsigned char*)p->knots + PR_ALIGN(sizeof(qaws_scalar) * p->kc));
	if (view->knots.count)
		for (k = 0; k < p->kc; k++)
			p->knots[k] = (qaws_scalar)vw_get(&view->knots, k, 0);
	else
		/* clamped uniform 0, 1, ..., n - d, as qaws_curve_create_bspline */
		for (k = 0; k < p->kc; k++)
			p->knots[k] = (qaws_scalar)(k <= d ? 0 : k >= n ? n - d : k - d);
	p->bounds[0] = p->knots[d];
	for (k = d; k < n; k++)
		if (p->knots[k + 1] > p->knots[k])
			p->bounds[++p->span_count] = p->knots[k + 1];
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
	if (bytes < pr_bytes(view->points.count, view->degree))
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	*out_prepared = pr_fill(view, memory);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_prepare(qaws_curve_view const* view, qaws_allocator const* allocator, qaws_curve_prepared** out_prepared)
{
	qaws_status st;
	size_t bytes;
	void* memory;
	if (!out_prepared)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_prepared = NULL;
	if ((st = qaws_curve_view_validate(view)) != QAWS_STATUS_OK)
		return st;
	bytes = pr_bytes(view->points.count, view->degree);
	/* (16-byte aligned: one extra block, the object at the first aligned address, the raw pointer just before it) */
	memory = qaws_internal_alloc(allocator, (unsigned long)(bytes + 16 + sizeof(void*)));
	if (!memory)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	{
		uintptr_t at = ((uintptr_t)memory + sizeof(void*) + 15) & ~(uintptr_t)15;
		qaws_curve_prepared* p = pr_fill(view, (void*)at);
		((void**)at)[-1] = memory;
		p->allocator = allocator;
		p->owned = 1;
		*out_prepared = p;
	}
	return QAWS_STATUS_OK;
}

void qaws_curve_prepared_destroy(qaws_curve_prepared* prepared)
{
	if (prepared && prepared->owned)
		qaws_internal_dealloc(prepared->allocator, ((void**)prepared)[-1]);
}

unsigned int qaws_curve_prepared_dependencies(qaws_curve_prepared const* prepared)
{
	return prepared ? QAWS_DEPENDS_LAYOUT | (prepared->view.knots.count ? QAWS_DEPENDS_KNOTS : 0u) : 0u;
}

qaws_range qaws_curve_prepared_get_range(qaws_curve_prepared const* prepared)
{
	qaws_range r;
	r.min_value = prepared ? prepared->knots[prepared->view.degree] : QAWS_ZERO;
	r.max_value = prepared ? prepared->knots[prepared->n] : QAWS_ZERO;
	return r;
}

#define PR_ORDERS 4

/* position and derivatives 1 to kmax at t (clamped to the domain), rows of out; the polynomial part as the owned B-spline
   (same span search, same basis, same sums), the rational part by the quotient rule */
static void pr_eval(qaws_curve_prepared const* p, qaws_scalar t, unsigned int kmax, qaws_scalar out[PR_ORDERS][3])
{
	unsigned int d = p->view.degree, dim = (unsigned int)p->view.dimension, span, j, r, c, kb = kmax < d ? kmax : d;
	qaws_scalar ders[PR_ORDERS * (QAWS_VIEW_MAX_DEGREE + 1)], A[PR_ORDERS][3], W[PR_ORDERS];
	if (t < p->knots[d])
		t = p->knots[d];
	if (t > p->knots[p->n])
		t = p->knots[p->n];
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

static unsigned int pr_order(unsigned int flags)
{
	return flags & QAWS_EVAL_FLAG_D3 ? 3 : flags & QAWS_EVAL_FLAG_D2 ? 2 : flags & QAWS_EVAL_FLAG_D1 ? 1 : 0;
}

qaws_status qaws_curve_prepared_evaluate_3d(qaws_curve_prepared const* prepared, qaws_scalar parameter, unsigned int eval_flags,
	qaws_eval_result_3d* out_result)
{
	qaws_scalar v[PR_ORDERS][3];
	qaws_vec3* dst[PR_ORDERS];
	unsigned int r;
	if (!prepared || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (prepared->view.dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	memset(out_result, 0, sizeof *out_result);
	pr_eval(prepared, parameter, pr_order(eval_flags), v);
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
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_prepared_evaluate_2d(qaws_curve_prepared const* prepared, qaws_scalar parameter, unsigned int eval_flags,
	qaws_eval_result_2d* out_result)
{
	qaws_scalar v[PR_ORDERS][3];
	qaws_vec2* dst[PR_ORDERS];
	unsigned int r;
	if (!prepared || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (prepared->view.dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	memset(out_result, 0, sizeof *out_result);
	pr_eval(prepared, parameter, pr_order(eval_flags), v);
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

static qaws_scalar bw_global(qaws_curve const* curve, unsigned int span, qaws_scalar local_t)
{
	return curve->span_boundaries[span] + local_t * (curve->span_boundaries[span + 1] - curve->span_boundaries[span]);
}

static qaws_status bw_eval_span_2d(qaws_curve const* curve, unsigned int span, qaws_scalar local_t, unsigned int flags, qaws_eval_result_2d* out)
{
	return qaws_curve_prepared_evaluate_2d((qaws_curve_prepared const*)curve->impl, bw_global(curve, span, local_t), flags, out);
}

static qaws_status bw_eval_span_3d(qaws_curve const* curve, unsigned int span, qaws_scalar local_t, unsigned int flags, qaws_eval_result_3d* out)
{
	return qaws_curve_prepared_evaluate_3d((qaws_curve_prepared const*)curve->impl, bw_global(curve, span, local_t), flags, out);
}

static void bw_destroy(void* impl, qaws_allocator const* allocator)
{
	(void)allocator;
	qaws_curve_prepared_destroy((qaws_curve_prepared*)impl);
}

static int bw_is_closed(qaws_curve const* curve)
{
	return qaws_internal_curve_ends_meet(curve);
}

static int bw_is_periodic(qaws_curve const* curve)
{
	(void)curve;
	return 0;
}

static int bw_is_rational(qaws_curve const* curve)
{
	return ((qaws_curve_prepared const*)curve->impl)->rational != 0;
}

/* the degree less the largest multiplicity of a knot inside the domain */
static qaws_continuity bw_get_continuity(qaws_curve const* curve)
{
	qaws_curve_prepared const* p = (qaws_curve_prepared const*)curve->impl;
	unsigned int d = p->view.degree, k, run = 1, worst = 0;
	int smooth;
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

/* ---- owned copies ------------------------------------------------------- */

/* an owned B-spline (polynomial) or NURBS curve with the prepared knots and the view's points (homogeneous ones divided by w) */
static qaws_status pr_copy(qaws_curve_prepared const* p, qaws_allocator const* allocator, qaws_curve** out_curve)
{
	unsigned int dim = (unsigned int)p->view.dimension, i, c;
	qaws_scalar* pts = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)(sizeof(qaws_scalar) * (size_t)p->n * (dim + 1)));
	qaws_scalar* wts;
	qaws_status st;
	if (!pts)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	wts = pts + (size_t)p->n * dim;
	for (i = 0; i < p->n; i++)
	{
		wts[i] = p->rational == 1 ? (qaws_scalar)vw_get(&p->view.weights, i, 0)
			: p->rational == 2 ? (qaws_scalar)vw_get(&p->view.points, i, dim) : QAWS_ONE;
		for (c = 0; c < dim; c++)
			pts[i * dim + c] = (qaws_scalar)vw_get(&p->view.points, i, c) / (p->rational == 2 ? wts[i] : QAWS_ONE);
	}
	if (p->rational)
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
