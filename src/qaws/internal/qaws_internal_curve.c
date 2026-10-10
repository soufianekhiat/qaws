#include "qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>

void* qaws_internal_alloc(qaws_allocator const* allocator, unsigned long size)
{
	if (allocator && allocator->alloc)
		return allocator->alloc(size, allocator->user_data);
	return malloc((size_t)size);
}

void qaws_internal_dealloc(qaws_allocator const* allocator, void* ptr)
{
	if (allocator && allocator->dealloc)
	{
		allocator->dealloc(ptr, allocator->user_data);
		return;
	}
	free(ptr);
}

qaws_curve* qaws_internal_curve_alloc_ex(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable,
	qaws_allocator const* allocator)
{
	qaws_curve* curve = (qaws_curve*)qaws_internal_alloc(allocator, sizeof(qaws_curve));
	if (!curve) return NULL;

	curve->span_boundaries = (qaws_scalar*)qaws_internal_alloc(
		allocator, (unsigned long)(span_count + 1) * sizeof(qaws_scalar));
	if (!curve->span_boundaries) {
		qaws_internal_dealloc(allocator, curve);
		return NULL;
	}

	memset(curve->span_boundaries, 0, (size_t)(span_count + 1) * sizeof(qaws_scalar));

	curve->kind = kind;
	curve->dimension = dimension;
	curve->degree = degree;
	curve->span_count = span_count;
	curve->parameter_range = parameter_range;
	curve->vtable = vtable;
	curve->impl = NULL;
	curve->allocator = allocator;

	return curve;
}

qaws_curve* qaws_internal_curve_alloc(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable)
{
	return qaws_internal_curve_alloc_ex(
		kind, dimension, degree, span_count,
		parameter_range, vtable, NULL);
}

void qaws_internal_curve_free(qaws_curve* curve)
{
	if (!curve) return;
	qaws_internal_dealloc(curve->allocator, curve->span_boundaries);
	qaws_internal_dealloc(curve->allocator, curve);
}

#include "qaws_internal_kinds.h"
#include "qaws_internal_span.h"
#include <math.h>

static int ends_point(qaws_curve const* curve, qaws_scalar t, double out[3])
{
	qaws_scalar local_t;
	unsigned int span = qaws_internal_find_span(curve, t, &local_t);
	out[0] = out[1] = out[2] = 0.0;
	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d r;
		if (curve->vtable->eval_span_2d(curve, span, local_t, QAWS_EVAL_FLAG_POSITION, &r) != QAWS_STATUS_OK)
			return 0;
		out[0] = (double)r.position.x; out[1] = (double)r.position.y;
	}
	else
	{
		qaws_eval_result_3d r;
		if (curve->vtable->eval_span_3d(curve, span, local_t, QAWS_EVAL_FLAG_POSITION, &r) != QAWS_STATUS_OK)
			return 0;
		out[0] = (double)r.position.x; out[1] = (double)r.position.y; out[2] = (double)r.position.z;
	}
	return 1;
}

int qaws_internal_curve_ends_meet(qaws_curve const* curve)
{
	/* The ends meet within a few rounding steps of the curve's extent,
	   measured on 16 samples. */
	double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
	double p0[3], p1[3], p[3], ext = 0.0, d = 0.0;
	qaws_scalar a = curve->parameter_range.min_value, b = curve->parameter_range.max_value;
	unsigned int i, k;
#if QAWS_SCALAR_IS_FLOAT
	double const rel = 1e-5;
#else
	double const rel = 1e-10;
#endif
	if (!ends_point(curve, a, p0) || !ends_point(curve, b, p1))
		return 0;
	for (i = 0; i <= 16; i++)
	{
		if (!ends_point(curve, a + (b - a) * (qaws_scalar)i / (qaws_scalar)16, p))
			return 0;
		for (k = 0; k < 3; k++)
		{
			if (p[k] < lo[k]) lo[k] = p[k];
			if (p[k] > hi[k]) hi[k] = p[k];
		}
	}
	for (k = 0; k < 3; k++)
	{
		ext += (hi[k] - lo[k]) * (hi[k] - lo[k]);
		d += (p1[k] - p0[k]) * (p1[k] - p0[k]);
	}
	return sqrt(d) <= rel * sqrt(ext) ? 1 : 0;
}
