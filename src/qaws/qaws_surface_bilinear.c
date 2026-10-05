#include "qaws_surface_bilinear.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_diff.h"
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
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
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

/* -------------------------------------------------------------------------- */
/*  Differential rules: linear in the four corners                            */
/* -------------------------------------------------------------------------- */

#define BILINEAR_SURFACE_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT | QAWS_CAP_LINEAR)

static unsigned int bilinear_surface_describe_fields(qaws_surface const* surface,
	qaws_field_desc* out, unsigned int capacity)
{
	(void)surface;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CONTROL_POINTS, QAWS_VALUE_VEC3,
			4, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, BILINEAR_SURFACE_DIFF_CAPS);
	return 1;
}

/* Corners are stored contiguously as P00, P10, P01, P11:
   element = v_index * 2 + u_index. */
static qaws_status bilinear_surface_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_surface_bilinear_impl const* impl = (qaws_surface_bilinear_impl const*)surface->impl;
	if (field != QAWS_FIELD_CONTROL_POINTS)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_data = &impl->p00.x;
	*out_count = 4;
	*out_components = 3;
	return QAWS_STATUS_OK;
}

static qaws_status bilinear_surface_linear_support(qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int order, qaws_surface_support* out)
{
	unsigned int r;
	(void)surface;
	if (order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	out->kind = QAWS_SUPPORT_GLOBAL;
	out->field = QAWS_FIELD_CONTROL_POINTS;
	out->u_first = 0;
	out->u_count = 2;
	out->v_first = 0;
	out->v_count = 2;
	out->u_stride = 1;
	out->v_stride = 2;
	out->on_boundary = 0;
	out->has_weights = 1;
	out->order = order;
	for (r = 0; r <= order; r++)
	{
		out->u_weights[r][0] = (r == 0) ? QAWS_ONE - u : (r == 1) ? -QAWS_ONE : QAWS_ZERO;
		out->u_weights[r][1] = (r == 0) ? u : (r == 1) ? QAWS_ONE : QAWS_ZERO;
		out->v_weights[r][0] = (r == 0) ? QAWS_ONE - v : (r == 1) ? -QAWS_ONE : QAWS_ZERO;
		out->v_weights[r][1] = (r == 0) ? v : (r == 1) ? QAWS_ONE : QAWS_ZERO;
	}
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const bilinear_surface_diff_vtable = {
	BILINEAR_SURFACE_DIFF_CAPS,
	QAWS_DIFF_SMOOTH,
	bilinear_surface_describe_fields,
	bilinear_surface_primal_field,
	bilinear_surface_linear_support,
	NULL
};

static qaws_surface_vtable const bilinear_surface_vtable = {
	bilinear_surface_eval,
	bilinear_surface_destroy,
	bilinear_surface_is_rational,
	&bilinear_surface_diff_vtable
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
