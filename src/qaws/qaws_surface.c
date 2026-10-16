#include "qaws_surface.h"
#include "internal/qaws_internal_surface.h"
#include <math.h>
#include <string.h>

qaws_surface_kind qaws_surface_get_kind(qaws_surface const* surface)
{
	return surface ? surface->kind : QAWS_SURFACE_KIND_INVALID;
}

unsigned int qaws_surface_get_u_degree(qaws_surface const* surface)
{
	return surface ? surface->u_degree : 0;
}

unsigned int qaws_surface_get_v_degree(qaws_surface const* surface)
{
	return surface ? surface->v_degree : 0;
}

qaws_range qaws_surface_get_u_range(qaws_surface const* surface)
{
	qaws_range r = {0, 0};
	if (surface) r = surface->u_range;
	return r;
}

qaws_range qaws_surface_get_v_range(qaws_surface const* surface)
{
	qaws_range r = {0, 0};
	if (surface) r = surface->v_range;
	return r;
}

int qaws_surface_is_rational(qaws_surface const* surface)
{
	if (!surface || !surface->vtable || !surface->vtable->is_rational) return 0;
	return surface->vtable->is_rational(surface);
}

qaws_status qaws_surface_evaluate(
	qaws_surface_eval_result* out_result,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags)
{
	if (!surface || !out_result) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!surface->vtable || !surface->vtable->eval_3d)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	/* Clamp to domain */
	if (u < surface->u_range.min_value) u = surface->u_range.min_value;
	if (u > surface->u_range.max_value) u = surface->u_range.max_value;
	if (v < surface->v_range.min_value) v = surface->v_range.min_value;
	if (v > surface->v_range.max_value) v = surface->v_range.max_value;

	memset(out_result, 0, sizeof(*out_result));

	/* If normal is requested, we also need du and dv */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
		eval_flags |= QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV;

	return surface->vtable->eval_3d(surface, u, v, eval_flags, out_result);
}

qaws_status qaws_surface_get_normal(
	qaws_vec3* out_normal,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v)
{
	qaws_surface_eval_result r;
	qaws_status s;
	if (!surface || !out_normal) return QAWS_STATUS_INVALID_ARGUMENT;
	s = qaws_surface_evaluate(&r, surface, u,
		v, QAWS_SURFACE_EVAL_NORMAL);
	if (s != QAWS_STATUS_OK) return s;
	*out_normal = r.normal;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_sample(
	qaws_vec3* out_positions,
	unsigned int* out_position_count,
	qaws_surface const* surface,
	unsigned int u_count,
	unsigned int v_count,
	unsigned int position_capacity)
{
	unsigned int total, ui, vi, idx;
	qaws_scalar u_min, u_len, v_min, v_len;
	if (!surface || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (u_count < 2 || v_count < 2) return QAWS_STATUS_INVALID_ARGUMENT;

	total = u_count * v_count;
	if (total > position_capacity) return QAWS_STATUS_BUFFER_TOO_SMALL;

	u_min = surface->u_range.min_value;
	u_len = surface->u_range.max_value - u_min;
	v_min = surface->v_range.min_value;
	v_len = surface->v_range.max_value - v_min;

	idx = 0;
	for (ui = 0; ui < u_count; ui++)
	{
		qaws_scalar u = u_min + u_len * (qaws_scalar)ui / (qaws_scalar)(u_count - 1);
		for (vi = 0; vi < v_count; vi++)
		{
			qaws_scalar v_val = v_min + v_len * (qaws_scalar)vi / (qaws_scalar)(v_count - 1);
			qaws_surface_eval_result r;
			qaws_status s = qaws_surface_evaluate(&r, surface, u,
				v_val, QAWS_SURFACE_EVAL_POSITION);
			if (s != QAWS_STATUS_OK) return s;
			out_positions[idx++] = r.position;
		}
	}
	*out_position_count = total;
	return QAWS_STATUS_OK;
}

void qaws_surface_destroy(qaws_surface* surface)
{
	qaws_internal_surface_free(surface);
}
