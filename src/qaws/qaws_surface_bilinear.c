#include "qaws_surface_bilinear.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_bilinear_impl
{
	qaws_vec3 p00;
	qaws_vec3 p10;
	qaws_vec3 p01;
	qaws_vec3 p11;
} qaws_surface_bilinear_impl;

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

/* Bilinear patch:
   S(u,v) = (1-u)(1-v)P00 + u(1-v)P10 + (1-u)vP01 + uvP11
   dS/du  = (1-v)(P10-P00) + v(P11-P01)
   dS/dv  = (1-u)(P01-P00) + u(P11-P10)
   d2S/du2  = 0
   d2S/dv2  = 0
   d2S/dudv = P00 - P10 - P01 + P11 */
static qaws_status bilinear_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_bilinear_impl const* impl =
		(qaws_surface_bilinear_impl const*)surface->impl;
	qaws_scalar one_minus_u = QAWS_ONE - u;
	qaws_scalar one_minus_v = QAWS_ONE - v;

	/* Position */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = one_minus_u * one_minus_v * impl->p00.x
			+ u * one_minus_v * impl->p10.x
			+ one_minus_u * v * impl->p01.x
			+ u * v * impl->p11.x;
		out_result->position.y = one_minus_u * one_minus_v * impl->p00.y
			+ u * one_minus_v * impl->p10.y
			+ one_minus_u * v * impl->p01.y
			+ u * v * impl->p11.y;
		out_result->position.z = one_minus_u * one_minus_v * impl->p00.z
			+ u * one_minus_v * impl->p10.z
			+ one_minus_u * v * impl->p01.z
			+ u * v * impl->p11.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = (1-v)(P10-P00) + v(P11-P01) */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL))
	{
		out_result->du.x = one_minus_v * (impl->p10.x - impl->p00.x)
			+ v * (impl->p11.x - impl->p01.x);
		out_result->du.y = one_minus_v * (impl->p10.y - impl->p00.y)
			+ v * (impl->p11.y - impl->p01.y);
		out_result->du.z = one_minus_v * (impl->p10.z - impl->p00.z)
			+ v * (impl->p11.z - impl->p01.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = (1-u)(P01-P00) + u(P11-P10) */
	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL))
	{
		out_result->dv.x = one_minus_u * (impl->p01.x - impl->p00.x)
			+ u * (impl->p11.x - impl->p10.x);
		out_result->dv.y = one_minus_u * (impl->p01.y - impl->p00.y)
			+ u * (impl->p11.y - impl->p10.y);
		out_result->dv.z = one_minus_u * (impl->p01.z - impl->p00.z)
			+ u * (impl->p11.z - impl->p10.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* d2S/du2 = 0 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		out_result->duu.x = 0; out_result->duu.y = 0; out_result->duu.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* d2S/dv2 = 0 */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv.x = 0; out_result->dvv.y = 0; out_result->dvv.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* d2S/dudv = P00 - P10 - P01 + P11 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		out_result->duv.x = impl->p00.x - impl->p10.x - impl->p01.x + impl->p11.x;
		out_result->duv.y = impl->p00.y - impl->p10.y - impl->p01.y + impl->p11.y;
		out_result->duv.z = impl->p00.z - impl->p10.z - impl->p01.z + impl->p11.z;
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

static void bilinear_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int bilinear_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const bilinear_surface_vtable = {
	bilinear_surface_eval,
	bilinear_surface_destroy,
	bilinear_surface_is_rational
};

qaws_status qaws_surface_create_bilinear(
	qaws_surface_bilinear_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_bilinear_impl* impl;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_BILINEAR,
		1, 1, u_range, v_range,
		&bilinear_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_bilinear_impl*)malloc(sizeof(qaws_surface_bilinear_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->p00 = desc->p00;
	impl->p10 = desc->p10;
	impl->p01 = desc->p01;
	impl->p11 = desc->p11;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
