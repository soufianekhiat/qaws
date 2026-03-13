#include "qaws_surface_offset.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_offset_impl
{
	qaws_surface const* base;
	qaws_scalar distance;
} qaws_surface_offset_impl;

static void compute_normal(qaws_vec3 du, qaws_vec3 dv, qaws_vec3* out)
{
	qaws_scalar nx = du.y * dv.z - du.z * dv.y;
	qaws_scalar ny = du.z * dv.x - du.x * dv.z;
	qaws_scalar nz = du.x * dv.y - du.y * dv.x;
	qaws_scalar len = QAWS_SQRT(nx * nx + ny * ny + nz * nz);
	if (len > QAWS_LITERAL(1e-12))
	{
		out->x = nx / len; out->y = ny / len; out->z = nz / len;
	}
	else
	{
		out->x = 0; out->y = 0; out->z = 1;
	}
}

/* Helper: evaluate offset position at (u,v).
   S_off(u,v) = S_base(u,v) + d * N_base(u,v) */
static qaws_status offset_eval_position(
	qaws_surface const* base,
	qaws_scalar distance,
	qaws_scalar u,
	qaws_scalar v,
	qaws_vec3* out_pos)
{
	qaws_surface_eval_result base_r;
	qaws_status status;
	unsigned int flags = QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU
		| QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL;

	memset(&base_r, 0, sizeof(base_r));
	status = qaws_surface_evaluate(base, u, v, flags, &base_r);
	if (status != QAWS_STATUS_OK) return status;

	out_pos->x = base_r.position.x + distance * base_r.normal.x;
	out_pos->y = base_r.position.y + distance * base_r.normal.y;
	out_pos->z = base_r.position.z + distance * base_r.normal.z;
	return QAWS_STATUS_OK;
}

/* Offset surface evaluation.
   Position: S_off(u,v) = S_base(u,v) + d * N_base(u,v)
   Derivatives: finite differences of offset position to avoid Weingarten map. */
static qaws_status offset_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_offset_impl const* impl =
		(qaws_surface_offset_impl const*)surface->impl;
	qaws_vec3 pos;
	qaws_status status;

	/* Always compute position (needed by derivatives too) */
	status = offset_eval_position(impl->base, impl->distance, u, v, &pos);
	if (status != QAWS_STATUS_OK) return status;

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* First derivatives via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_scalar u_lo = u - h;
		qaws_scalar u_hi = u + h;
		qaws_scalar hu;
		qaws_vec3 p_lo, p_hi;

		if (u_lo < QAWS_ZERO) u_lo = QAWS_ZERO;
		if (u_hi > QAWS_ONE) u_hi = QAWS_ONE;
		hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

		status = offset_eval_position(impl->base, impl->distance, u_lo, v, &p_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_hi, v, &p_hi);
		if (status != QAWS_STATUS_OK) return status;

		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
			out_result->du.x = (p_hi.x - p_lo.x) * inv2h;
			out_result->du.y = (p_hi.y - p_lo.y) * inv2h;
			out_result->du.z = (p_hi.z - p_lo.z) * inv2h;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}

		/* duu via central 2nd-order finite difference */
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar inv_h2 = QAWS_ONE / (hu * hu);
			out_result->duu.x = (p_hi.x - QAWS_LITERAL(2.0) * pos.x + p_lo.x) * inv_h2;
			out_result->duu.y = (p_hi.y - QAWS_LITERAL(2.0) * pos.y + p_lo.y) * inv_h2;
			out_result->duu.z = (p_hi.z - QAWS_LITERAL(2.0) * pos.z + p_lo.z) * inv_h2;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}
	}

	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_scalar v_lo = v - h;
		qaws_scalar v_hi = v + h;
		qaws_scalar hv;
		qaws_vec3 p_lo, p_hi;

		if (v_lo < QAWS_ZERO) v_lo = QAWS_ZERO;
		if (v_hi > QAWS_ONE) v_hi = QAWS_ONE;
		hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

		status = offset_eval_position(impl->base, impl->distance, u, v_lo, &p_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u, v_hi, &p_hi);
		if (status != QAWS_STATUS_OK) return status;

		if (eval_flags & QAWS_SURFACE_EVAL_DV)
		{
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
			out_result->dv.x = (p_hi.x - p_lo.x) * inv2h;
			out_result->dv.y = (p_hi.y - p_lo.y) * inv2h;
			out_result->dv.z = (p_hi.z - p_lo.z) * inv2h;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
		}

		/* dvv via central 2nd-order finite difference */
		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		{
			qaws_scalar inv_h2 = QAWS_ONE / (hv * hv);
			out_result->dvv.x = (p_hi.x - QAWS_LITERAL(2.0) * pos.x + p_lo.x) * inv_h2;
			out_result->dvv.y = (p_hi.y - QAWS_LITERAL(2.0) * pos.y + p_lo.y) * inv_h2;
			out_result->dvv.z = (p_hi.z - QAWS_LITERAL(2.0) * pos.z + p_lo.z) * inv_h2;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}
	}

	/* duv via finite difference of dv w.r.t. u */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_scalar u_lo = u - h;
		qaws_scalar u_hi = u + h;
		qaws_scalar hu;
		qaws_scalar hv_val = QAWS_LITERAL(1e-5);
		qaws_scalar v_lo, v_hi, hv;
		qaws_vec3 pv_lo_lo, pv_lo_hi, pv_hi_lo, pv_hi_hi;
		qaws_vec3 dv_lo, dv_hi;
		qaws_scalar inv2h;

		if (u_lo < QAWS_ZERO) u_lo = QAWS_ZERO;
		if (u_hi > QAWS_ONE) u_hi = QAWS_ONE;
		hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

		v_lo = v - hv_val;
		v_hi = v + hv_val;
		if (v_lo < QAWS_ZERO) v_lo = QAWS_ZERO;
		if (v_hi > QAWS_ONE) v_hi = QAWS_ONE;
		hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

		/* dv at u_lo and u_hi via finite differences in v */
		status = offset_eval_position(impl->base, impl->distance, u_lo, v_lo, &pv_lo_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_lo, v_hi, &pv_lo_hi);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_hi, v_lo, &pv_hi_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_hi, v_hi, &pv_hi_hi);
		if (status != QAWS_STATUS_OK) return status;

		{
			qaws_scalar inv2hv = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
			dv_lo.x = (pv_lo_hi.x - pv_lo_lo.x) * inv2hv;
			dv_lo.y = (pv_lo_hi.y - pv_lo_lo.y) * inv2hv;
			dv_lo.z = (pv_lo_hi.z - pv_lo_lo.z) * inv2hv;
			dv_hi.x = (pv_hi_hi.x - pv_hi_lo.x) * inv2hv;
			dv_hi.y = (pv_hi_hi.y - pv_hi_lo.y) * inv2hv;
			dv_hi.z = (pv_hi_hi.z - pv_hi_lo.z) * inv2hv;
		}

		inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
		out_result->duv.x = (dv_hi.x - dv_lo.x) * inv2h;
		out_result->duv.y = (dv_hi.y - dv_lo.y) * inv2h;
		out_result->duv.z = (dv_hi.z - dv_lo.z) * inv2h;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
	}

	/* Normal */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		compute_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

static void offset_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int offset_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const offset_surface_vtable = {
	offset_surface_eval,
	offset_surface_destroy,
	offset_surface_is_rational
};

qaws_status qaws_surface_create_offset(
	qaws_surface_offset_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_offset_impl* impl;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->base) return QAWS_STATUS_INVALID_ARGUMENT;
	if (QAWS_FABS(desc->distance) < QAWS_LITERAL(1e-12))
		return QAWS_STATUS_INVALID_ARGUMENT;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_OFFSET,
		0, 0, u_range, v_range,
		&offset_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_offset_impl*)malloc(sizeof(qaws_surface_offset_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->base = desc->base;
	impl->distance = desc->distance;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
