#include "qaws_surface_trim.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_trim_loop_internal
{
	qaws_curve const** curves;
	unsigned int curve_count;
	int is_outer;
} qaws_trim_loop_internal;

typedef struct qaws_surface_trim_impl
{
	qaws_surface const* base;
	qaws_trim_loop_internal* loops;
	unsigned int loop_count;
} qaws_surface_trim_impl;

/* The winding contribution of the polygon edge (x1, y1) -> (x2, y2) around (px, py): +1 an upward crossing with the point left of
   the edge, -1 a downward one with it right, else 0 */
static int winding_edge(qaws_scalar x1, qaws_scalar y1, qaws_scalar x2, qaws_scalar y2, qaws_scalar px, qaws_scalar py)
{
	qaws_scalar cross = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1);
	if (y1 <= py)
		return y2 > py && cross > QAWS_ZERO ? 1 : 0;
	return y2 <= py && cross < QAWS_ZERO ? -1 : 0;
}

/* Sample a trim loop into a polygon and test containment: the winding number summed edge by edge as the samples come (closed by the
   edge from the last sample back to the first), so a query allocates nothing */
static int loop_contains_point(
	qaws_trim_loop_internal const* loop,
	qaws_scalar u, qaws_scalar v)
{
	unsigned int samples_per_curve = 100;
	unsigned int ci, count = 0;
	qaws_scalar x0 = QAWS_ZERO, y0 = QAWS_ZERO, xp = QAWS_ZERO, yp = QAWS_ZERO;
	int winding = 0;

	for (ci = 0; ci < loop->curve_count; ++ci)
	{
		qaws_curve const* crv = loop->curves[ci];
		qaws_range range = qaws_curve_get_parameter_range(crv);
		qaws_scalar span = range.max_value - range.min_value;
		unsigned int si;

		for (si = 0; si < samples_per_curve; ++si)
		{
			qaws_scalar t_frac = (qaws_scalar)si / (qaws_scalar)samples_per_curve;
			qaws_scalar t = range.min_value + t_frac * span;
			qaws_eval_result_2d r2d;
			qaws_scalar x = QAWS_ZERO, y = QAWS_ZERO;

			memset(&r2d, 0, sizeof(r2d));
			if (qaws_curve_evaluate_2d(&r2d, crv, t, QAWS_EVAL_FLAG_POSITION) == QAWS_STATUS_OK)
			{
				x = r2d.position.x;
				y = r2d.position.y;
			}
			if (count++ == 0)
			{
				x0 = x;
				y0 = y;
			}
			else
				winding += winding_edge(xp, yp, x, y, u, v);
			xp = x;
			yp = y;
		}
	}
	if (count)
		winding += winding_edge(xp, yp, x0, y0, u, v);

	return (winding != 0) ? 1 : 0;
}

/* Trimmed surface evaluation: delegates entirely to base surface. */
static qaws_status trimmed_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_trim_impl const* impl =
		(qaws_surface_trim_impl const*)surface->impl;
	return qaws_surface_evaluate(out_result, impl->base, u, v, eval_flags);
}

static void trimmed_surface_destroy(void* impl_ptr, qaws_allocator const* allocator)
{
	qaws_surface_trim_impl* impl = (qaws_surface_trim_impl*)impl_ptr;
	if (impl)
	{
		unsigned int i;
		if (impl->loops)
		{
			for (i = 0; i < impl->loop_count; ++i)
			{
				if (impl->loops[i].curves)
					qaws_internal_dealloc(impl->loops[i].curves, allocator);
			}
			qaws_internal_dealloc(impl->loops, allocator);
		}
	}
	qaws_internal_dealloc(impl, allocator);
}

static int trimmed_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const trimmed_surface_vtable = {
	trimmed_surface_eval,
	trimmed_surface_destroy,
	trimmed_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_trimmed_ex(
	qaws_surface_trim_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_trim_impl* impl;
	qaws_range u_range, v_range;
	unsigned int i;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->base) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->loop_count > 0 && !desc->loops)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate trim loops */
	for (i = 0; i < desc->loop_count; ++i)
	{
		unsigned int ci;
		if (desc->loops[i].curve_count == 0 || !desc->loops[i].curves)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (ci = 0; ci < desc->loops[i].curve_count; ++ci)
		{
			if (!desc->loops[i].curves[ci])
				return QAWS_STATUS_INVALID_ARGUMENT;
			if (qaws_curve_get_dimension(desc->loops[i].curves[ci]) != QAWS_DIMENSION_2D)
				return QAWS_STATUS_INVALID_DIMENSION;
		}
	}

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc_ex(QAWS_SURFACE_KIND_TRIMMED,
		0, 0, u_range, v_range,
		&trimmed_surface_vtable, allocator);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_trim_impl*)qaws_internal_alloc(sizeof(qaws_surface_trim_impl), allocator);
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->base = desc->base;
	impl->loop_count = desc->loop_count;
	impl->loops = NULL;

	if (desc->loop_count > 0)
	{
		impl->loops = (qaws_trim_loop_internal*)qaws_internal_alloc(qaws_internal_mul2(desc->loop_count, sizeof(qaws_trim_loop_internal)), allocator);
		if (!impl->loops)
		{
			qaws_internal_dealloc(impl, allocator);
			qaws_internal_surface_free(surface);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		for (i = 0; i < desc->loop_count; ++i)
		{
			unsigned int ci;
			qaws_trim_loop const* src = &desc->loops[i];

			impl->loops[i].curve_count = src->curve_count;
			impl->loops[i].is_outer = src->is_outer;
			impl->loops[i].curves = (qaws_curve const**)qaws_internal_alloc(qaws_internal_mul2(src->curve_count, sizeof(qaws_curve const*)), allocator);
			if (!impl->loops[i].curves)
			{
				/* Clean up previously allocated loops */
				unsigned int j;
				for (j = 0; j < i; ++j)
					qaws_internal_dealloc(impl->loops[j].curves, allocator);
				qaws_internal_dealloc(impl->loops, allocator);
				qaws_internal_dealloc(impl, allocator);
				qaws_internal_surface_free(surface);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			for (ci = 0; ci < src->curve_count; ++ci)
				impl->loops[i].curves[ci] = src->curves[ci];
		}
	}

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_create_trimmed(
	qaws_surface** out_surface,
	qaws_surface_trim_desc const* desc)
{
	return qaws_surface_create_trimmed_ex(desc, NULL, out_surface);
}

qaws_status qaws_surface_trim_contains(
	int* out_inside,
	qaws_surface const* trimmed_surface, qaws_scalar u,
	qaws_scalar v)
{
	qaws_surface_trim_impl const* impl;
	int inside = 1;
	int has_outer = 0;
	unsigned int i;

	if (!trimmed_surface || !out_inside)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_surface_get_kind(trimmed_surface) != QAWS_SURFACE_KIND_TRIMMED)
		return QAWS_STATUS_INVALID_ARGUMENT;

	impl = (qaws_surface_trim_impl const*)trimmed_surface->impl;

	/* Check outer loops first */
	for (i = 0; i < impl->loop_count; ++i)
	{
		if (impl->loops[i].is_outer)
		{
			has_outer = 1;
			if (!loop_contains_point(&impl->loops[i], u, v))
			{
				inside = 0;
				break;
			}
		}
	}

	/* If inside outer (or no outer loop), check hole loops */
	if (inside)
	{
		for (i = 0; i < impl->loop_count; ++i)
		{
			if (!impl->loops[i].is_outer)
			{
				if (loop_contains_point(&impl->loops[i], u, v))
				{
					inside = 0;
					break;
				}
			}
		}
	}

	(void)has_outer;

	*out_inside = inside;
	return QAWS_STATUS_OK;
}
