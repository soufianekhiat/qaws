#include "qaws_surface_ruled.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_ruled_impl
{
	qaws_curve const* curve_a;
	qaws_curve const* curve_b;
	qaws_range range_a;
	qaws_range range_b;
} qaws_surface_ruled_impl;

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

/* S(u,v) = (1-v) * A(u_a) + v * B(u_b)
   dS/du  = (1-v) * A'(u_a)*s_a + v * B'(u_b)*s_b
   dS/dv  = B(u_b) - A(u_a)
   dS/duu = (1-v) * A''(u_a)*s_a^2 + v * B''(u_b)*s_b^2
   dS/dvv = 0
   dS/duv = B'(u_b)*s_b - A'(u_a)*s_a
   where s_a = range_a length, s_b = range_b length */
static qaws_status ruled_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_ruled_impl const* impl =
		(qaws_surface_ruled_impl const*)surface->impl;
	qaws_scalar one_minus_v = QAWS_ONE - v;
	qaws_scalar s_a = impl->range_a.max_value - impl->range_a.min_value;
	qaws_scalar s_b = impl->range_b.max_value - impl->range_b.min_value;
	qaws_scalar u_a = impl->range_a.min_value + u * s_a;
	qaws_scalar u_b = impl->range_b.min_value + u * s_b;

	unsigned int curve_flags = 0;
	qaws_eval_result_3d ra, rb;
	qaws_status sa_status, sb_status;

	/* Determine what we need from the curves */
	if (eval_flags & (QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DV))
		curve_flags |= QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUV))
		curve_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		curve_flags |= QAWS_EVAL_FLAG_D2;

	memset(&ra, 0, sizeof(ra));
	memset(&rb, 0, sizeof(rb));

	sa_status = qaws_curve_evaluate_3d(impl->curve_a, u_a, curve_flags, &ra);
	sb_status = qaws_curve_evaluate_3d(impl->curve_b, u_b, curve_flags, &rb);
	if (sa_status != QAWS_STATUS_OK) return sa_status;
	if (sb_status != QAWS_STATUS_OK) return sb_status;

	/* Position: (1-v)*A + v*B */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = one_minus_v * ra.position.x + v * rb.position.x;
		out_result->position.y = one_minus_v * ra.position.y + v * rb.position.y;
		out_result->position.z = one_minus_v * ra.position.z + v * rb.position.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* du: (1-v)*A'*s_a + v*B'*s_b */
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = one_minus_v * ra.d1.x * s_a + v * rb.d1.x * s_b;
		out_result->du.y = one_minus_v * ra.d1.y * s_a + v * rb.d1.y * s_b;
		out_result->du.z = one_minus_v * ra.d1.z * s_a + v * rb.d1.z * s_b;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dv: B - A */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = rb.position.x - ra.position.x;
		out_result->dv.y = rb.position.y - ra.position.y;
		out_result->dv.z = rb.position.z - ra.position.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* duu: (1-v)*A''*s_a^2 + v*B''*s_b^2 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		qaws_scalar sa2 = s_a * s_a, sb2 = s_b * s_b;
		out_result->duu.x = one_minus_v * ra.d2.x * sa2 + v * rb.d2.x * sb2;
		out_result->duu.y = one_minus_v * ra.d2.y * sa2 + v * rb.d2.y * sb2;
		out_result->duu.z = one_minus_v * ra.d2.z * sa2 + v * rb.d2.z * sb2;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* dvv: 0 (linear in v) */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv.x = 0; out_result->dvv.y = 0; out_result->dvv.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* duv: B'*s_b - A'*s_a */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		out_result->duv.x = rb.d1.x * s_b - ra.d1.x * s_a;
		out_result->duv.y = rb.d1.y * s_b - ra.d1.y * s_a;
		out_result->duv.z = rb.d1.z * s_b - ra.d1.z * s_a;
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

static void ruled_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	/* We don't own the curves - just free the impl struct */
	qaws_internal_dealloc(allocator, impl);
}

static int ruled_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const ruled_surface_vtable = {
	ruled_surface_eval,
	ruled_surface_destroy,
	ruled_surface_is_rational
};

qaws_status qaws_surface_create_ruled(
	qaws_surface_ruled_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_ruled_impl* impl;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->curve_a || !desc->curve_b) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(desc->curve_a) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->curve_b) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_RULED,
		0, 0, u_range, v_range,
		&ruled_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_ruled_impl*)malloc(sizeof(qaws_surface_ruled_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->curve_a = desc->curve_a;
	impl->curve_b = desc->curve_b;
	impl->range_a = qaws_curve_get_parameter_range(desc->curve_a);
	impl->range_b = qaws_curve_get_parameter_range(desc->curve_b);

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
