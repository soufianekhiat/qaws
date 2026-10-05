#include "qaws_surface_bezier.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_diff.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* Evaluate a row of Bezier CPs at parameter t using de Casteljau.
   cp: flat array of 3 scalars per point, stride = v_count * 3.
   degree: polynomial degree in this direction.
   dim: 3 (always 3D for surfaces).
   Returns interpolated 3-vector in out. */
static void decasteljau_row(
	qaws_scalar const* cp,
	unsigned int degree,
	unsigned int stride,
	qaws_scalar t,
	qaws_scalar* out)
{
	qaws_scalar buf[64 * 3]; /* max degree 63 */
	unsigned int i, k;
	for (i = 0; i <= degree; i++)
	{
		buf[i * 3 + 0] = cp[i * stride + 0];
		buf[i * 3 + 1] = cp[i * stride + 1];
		buf[i * 3 + 2] = cp[i * stride + 2];
	}
	for (k = 1; k <= degree; k++)
		for (i = 0; i <= degree - k; i++)
		{
			buf[i * 3 + 0] = (QAWS_ONE - t) * buf[i * 3 + 0] + t * buf[(i + 1) * 3 + 0];
			buf[i * 3 + 1] = (QAWS_ONE - t) * buf[i * 3 + 1] + t * buf[(i + 1) * 3 + 1];
			buf[i * 3 + 2] = (QAWS_ONE - t) * buf[i * 3 + 2] + t * buf[(i + 1) * 3 + 2];
		}
	out[0] = buf[0]; out[1] = buf[1]; out[2] = buf[2];
}

/* Compute derivative control points: D[i] = degree * (P[i+1] - P[i]) */
static void deriv_cps(
	qaws_scalar const* cp,
	unsigned int degree,
	unsigned int stride,
	unsigned int count,
	qaws_scalar* out,
	unsigned int out_stride)
{
	unsigned int i;
	qaws_scalar d = (qaws_scalar)degree;
	for (i = 0; i < count; i++)
	{
		out[i * out_stride + 0] = d * (cp[(i + 1) * stride + 0] - cp[i * stride + 0]);
		out[i * out_stride + 1] = d * (cp[(i + 1) * stride + 1] - cp[i * stride + 1]);
		out[i * out_stride + 2] = d * (cp[(i + 1) * stride + 2] - cp[i * stride + 2]);
	}
}

static qaws_status bezier_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_bezier_impl const* impl =
		(qaws_surface_bezier_impl const*)surface->impl;
	unsigned int ud = surface->u_degree;
	unsigned int vd = surface->v_degree;
	unsigned int uc = impl->u_count;
	unsigned int vc = impl->v_count;
	qaws_scalar const* cp = impl->control_points;
	unsigned int row_stride = vc * 3; /* scalars per row */

	/* Step 1: For each u-row, evaluate at v to get a column of (ud+1) points */
	qaws_scalar col[64 * 3]; /* evaluated column */
	unsigned int i;

	for (i = 0; i < uc; i++)
		decasteljau_row(cp + i * row_stride, vd, 3, v, col + i * 3);

	/* Position: evaluate column at u */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		qaws_scalar pt[3];
		qaws_internal_decasteljau(col, ud, 3, u, pt);
		out_result->position.x = pt[0];
		out_result->position.y = pt[1];
		out_result->position.z = pt[2];
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* Partial derivative du: differentiate in u direction */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU))
	{
		qaws_scalar du_col[64 * 3]; /* derivative column */
		qaws_scalar pt[3];
		deriv_cps(col, ud, 3, ud, du_col, 3);
		qaws_internal_decasteljau(du_col, ud - 1, 3, u, pt);
		out_result->du.x = pt[0]; out_result->du.y = pt[1]; out_result->du.z = pt[2];
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;

		if ((eval_flags & QAWS_SURFACE_EVAL_DUU) && ud >= 2)
		{
			qaws_scalar duu_col[64 * 3];
			deriv_cps(du_col, ud - 1, 3, ud - 1, duu_col, 3);
			qaws_internal_decasteljau(duu_col, ud - 2, 3, u, pt);
			out_result->duu.x = pt[0]; out_result->duu.y = pt[1]; out_result->duu.z = pt[2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}
	}

	/* Partial derivative dv: differentiate in v direction, then evaluate in u */
	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DVV))
	{
		qaws_scalar dv_col[64 * 3];
		qaws_scalar pt[3];
		/* For each u-row, evaluate dv of the row at v */
		for (i = 0; i < uc; i++)
		{
			/* Compute v-derivative CPs for this row, then evaluate at v */
			qaws_scalar dv_row[64 * 3];
			deriv_cps(cp + i * row_stride, vd, 3, vd, dv_row, 3);
			qaws_internal_decasteljau(dv_row, vd - 1, 3, v, dv_col + i * 3);
		}
		/* Evaluate the dv column at u */
		qaws_internal_decasteljau(dv_col, ud, 3, u, pt);
		out_result->dv.x = pt[0]; out_result->dv.y = pt[1]; out_result->dv.z = pt[2];
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;

		if ((eval_flags & QAWS_SURFACE_EVAL_DVV) && vd >= 2)
		{
			qaws_scalar dvv_col[64 * 3];
			for (i = 0; i < uc; i++)
			{
				qaws_scalar dv_row[64 * 3], dvv_row[64 * 3];
				deriv_cps(cp + i * row_stride, vd, 3, vd, dv_row, 3);
				deriv_cps(dv_row, vd - 1, 3, vd - 1, dvv_row, 3);
				qaws_internal_decasteljau(dvv_row, vd - 2, 3, v, dvv_col + i * 3);
			}
			qaws_internal_decasteljau(dvv_col, ud, 3, u, pt);
			out_result->dvv.x = pt[0]; out_result->dvv.y = pt[1]; out_result->dvv.z = pt[2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}
	}

	/* Mixed partial duv */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		/* Differentiate in v first, then in u */
		qaws_scalar dv_col[64 * 3], duv_col[64 * 3];
		qaws_scalar pt[3];
		for (i = 0; i < uc; i++)
		{
			qaws_scalar dv_row[64 * 3];
			deriv_cps(cp + i * row_stride, vd, 3, vd, dv_row, 3);
			qaws_internal_decasteljau(dv_row, vd - 1, 3, v, dv_col + i * 3);
		}
		deriv_cps(dv_col, ud, 3, ud, duv_col, 3);
		qaws_internal_decasteljau(duv_col, ud - 1, 3, u, pt);
		out_result->duv.x = pt[0]; out_result->duv.y = pt[1]; out_result->duv.z = pt[2];
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

static void bezier_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_bezier_impl* bi = (qaws_surface_bezier_impl*)impl;
	if (bi)
	{
		qaws_internal_dealloc(allocator, bi->control_points);
		qaws_internal_dealloc(allocator, bi);
	}
}

static int bezier_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

/* -------------------------------------------------------------------------- */
/*  Differential rules: linear in the control net                             */
/* -------------------------------------------------------------------------- */

#define BEZIER_SURFACE_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT | QAWS_CAP_LINEAR)

static unsigned int bezier_surface_describe_fields(qaws_surface const* surface,
	qaws_field_desc* out, unsigned int capacity)
{
	qaws_surface_bezier_impl const* impl = (qaws_surface_bezier_impl const*)surface->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CONTROL_POINTS, QAWS_VALUE_VEC3,
			impl->u_count * impl->v_count, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, BEZIER_SURFACE_DIFF_CAPS);
	return 1;
}

static qaws_status bezier_surface_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_surface_bezier_impl const* impl = (qaws_surface_bezier_impl const*)surface->impl;
	if (field != QAWS_FIELD_CONTROL_POINTS)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_data = impl->control_points;
	*out_count = impl->u_count * impl->v_count;
	*out_components = 3;
	return QAWS_STATUS_OK;
}

static qaws_status bezier_surface_linear_support(qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int order, qaws_surface_support* out)
{
	qaws_surface_bezier_impl const* impl = (qaws_surface_bezier_impl const*)surface->impl;
	qaws_scalar bu[(QAWS_DIFF_MAX_ORDER + 1) * QAWS_DIFF_MAX_SUPPORT];
	qaws_scalar bv[(QAWS_DIFF_MAX_ORDER + 1) * QAWS_DIFF_MAX_SUPPORT];
	unsigned int ud = surface->u_degree, vd = surface->v_degree, r, i;

	if (ud + 1 > QAWS_DIFF_MAX_SUPPORT || vd + 1 > QAWS_DIFF_MAX_SUPPORT || order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	qaws_internal_bernstein_derivs(ud, u, order, bu);
	qaws_internal_bernstein_derivs(vd, v, order, bv);

	out->kind = QAWS_SUPPORT_GLOBAL;
	out->field = QAWS_FIELD_CONTROL_POINTS;
	out->u_first = 0;
	out->u_count = ud + 1;
	out->v_first = 0;
	out->v_count = vd + 1;
	out->u_stride = impl->v_count;
	out->v_stride = 1;
	out->on_boundary = 0;
	out->has_weights = 1;
	out->order = order;
	for (r = 0; r <= order; r++)
	{
		for (i = 0; i <= ud; i++)
			out->u_weights[r][i] = bu[r * (ud + 1) + i];
		for (i = 0; i <= vd; i++)
			out->v_weights[r][i] = bv[r * (vd + 1) + i];
	}
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const bezier_surface_diff_vtable = {
	BEZIER_SURFACE_DIFF_CAPS,
	QAWS_DIFF_SMOOTH,
	bezier_surface_describe_fields,
	bezier_surface_primal_field,
	bezier_surface_linear_support,
	NULL
};

static qaws_surface_vtable const bezier_surface_vtable = {
	bezier_surface_eval,
	bezier_surface_destroy,
	bezier_surface_is_rational,
	&bezier_surface_diff_vtable
};

qaws_status qaws_surface_create_bezier_ex(
	qaws_surface_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_bezier_impl* impl;
	size_t cp_size;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->u_degree < 1 || desc->v_degree < 1) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->u_point_count != desc->u_degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (desc->v_point_count != desc->v_degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points) return QAWS_STATUS_INVALID_ARGUMENT;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc_ex(
		QAWS_SURFACE_KIND_BEZIER,
		desc->u_degree, desc->v_degree,
		u_range, v_range,
		&bezier_surface_vtable,
		allocator);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_bezier_impl*)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_surface_bezier_impl));
	if (!impl) { qaws_internal_surface_free(surface); return QAWS_STATUS_ALLOCATION_FAILURE; }

	cp_size = (size_t)desc->u_point_count * desc->v_point_count * 3 * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)cp_size);
	if (!impl->control_points)
	{
		qaws_internal_dealloc(allocator, impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points, cp_size);
	impl->u_count = desc->u_point_count;
	impl->v_count = desc->v_point_count;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_create_bezier(
	qaws_surface_bezier_desc const* desc,
	qaws_surface** out_surface)
{
	return qaws_surface_create_bezier_ex(desc, NULL, out_surface);
}
