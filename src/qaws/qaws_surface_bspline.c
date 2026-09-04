#include "qaws_surface_bspline.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* Evaluate a tensor-product B-spline surface at (u,v).
   Uses basis function derivatives from qaws_internal_bspline_basis_derivs. */
static qaws_status bspline_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_bspline_impl const* impl =
		(qaws_surface_bspline_impl const*)surface->impl;
	unsigned int ud = surface->u_degree;
	unsigned int vd = surface->v_degree;
	unsigned int uc = impl->u_count;
	unsigned int vc = impl->v_count;
	qaws_scalar const* cp = impl->control_points;
	unsigned int u_order = ud + 1;
	unsigned int v_order = vd + 1;

	unsigned int u_span, v_span;
	unsigned int max_u_deriv = 0, max_v_deriv = 0;

	/* Determine derivative orders needed */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU))
		max_u_deriv = (eval_flags & QAWS_SURFACE_EVAL_DUU) ? 2 : 1;
	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DVV))
		max_v_deriv = (eval_flags & QAWS_SURFACE_EVAL_DVV) ? 2 : 1;
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		if (max_u_deriv < 1) max_u_deriv = 1;
		if (max_v_deriv < 1) max_v_deriv = 1;
	}
	if (max_u_deriv > ud) max_u_deriv = ud;
	if (max_v_deriv > vd) max_v_deriv = vd;

	/* Find knot spans */
	u_span = qaws_internal_find_knot_span(
		impl->u_knots, impl->u_knot_count, ud, uc, u);
	v_span = qaws_internal_find_knot_span(
		impl->v_knots, impl->v_knot_count, vd, vc, v);

	/* Compute basis function derivatives */
	{
		qaws_scalar u_ders[3 * 64]; /* (max_u_deriv+1) * u_order */
		qaws_scalar v_ders[3 * 64]; /* (max_v_deriv+1) * v_order */
		unsigned int ku, kv, i, j;
		qaws_scalar S[3][3][3]; /* S[ku][kv] = d^(ku+kv) S / du^ku dv^kv */

		qaws_internal_bspline_basis_derivs(
			impl->u_knots, impl->u_knot_count,
			ud, u_span, u, max_u_deriv, u_ders);
		qaws_internal_bspline_basis_derivs(
			impl->v_knots, impl->v_knot_count,
			vd, v_span, v, max_v_deriv, v_ders);

		/* Compute S[ku][kv] = sum_i sum_j N'_i^(ku)(u) * N'_j^(kv)(v) * P[i][j] */
		memset(S, 0, sizeof(S));
		for (ku = 0; ku <= max_u_deriv; ku++)
		{
			for (kv = 0; kv <= max_v_deriv; kv++)
			{
				for (i = 0; i < u_order; i++)
				{
					unsigned int ci = u_span - ud + i;
					qaws_scalar Nu = u_ders[ku * u_order + i];
					for (j = 0; j < v_order; j++)
					{
						unsigned int cj = v_span - vd + j;
						qaws_scalar Nv = v_ders[kv * v_order + j];
						qaws_scalar w = Nu * Nv;
						unsigned int cp_idx = (ci * vc + cj) * 3;
						S[ku][kv][0] += w * cp[cp_idx + 0];
						S[ku][kv][1] += w * cp[cp_idx + 1];
						S[ku][kv][2] += w * cp[cp_idx + 2];
					}
				}
			}
		}

		/* Position */
		if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
		{
			out_result->position.x = S[0][0][0];
			out_result->position.y = S[0][0][1];
			out_result->position.z = S[0][0][2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
		}

		/* du */
		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			out_result->du.x = S[1][0][0];
			out_result->du.y = S[1][0][1];
			out_result->du.z = S[1][0][2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}

		/* dv */
		if (eval_flags & QAWS_SURFACE_EVAL_DV)
		{
			out_result->dv.x = S[0][1][0];
			out_result->dv.y = S[0][1][1];
			out_result->dv.z = S[0][1][2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
		}

		/* duu */
		if ((eval_flags & QAWS_SURFACE_EVAL_DUU) && max_u_deriv >= 2)
		{
			out_result->duu.x = S[2][0][0];
			out_result->duu.y = S[2][0][1];
			out_result->duu.z = S[2][0][2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		/* dvv */
		if ((eval_flags & QAWS_SURFACE_EVAL_DVV) && max_v_deriv >= 2)
		{
			out_result->dvv.x = S[0][2][0];
			out_result->dvv.y = S[0][2][1];
			out_result->dvv.z = S[0][2][2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}

		/* duv */
		if ((eval_flags & QAWS_SURFACE_EVAL_DUV) && max_u_deriv >= 1 && max_v_deriv >= 1)
		{
			out_result->duv.x = S[1][1][0];
			out_result->duv.y = S[1][1][1];
			out_result->duv.z = S[1][1][2];
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
		}

		/* Normal */
		if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
		{
			qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
		}
	}

	return QAWS_STATUS_OK;
}

static void bspline_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_bspline_impl* bi = (qaws_surface_bspline_impl*)impl;
	if (bi)
	{
		qaws_internal_dealloc(allocator, bi->control_points);
		qaws_internal_dealloc(allocator, bi->u_knots);
		qaws_internal_dealloc(allocator, bi->v_knots);
		qaws_internal_dealloc(allocator, bi);
	}
}

static int bspline_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const bspline_surface_vtable = {
	bspline_surface_eval,
	bspline_surface_destroy,
	bspline_surface_is_rational
};

qaws_status qaws_surface_create_bspline_ex(
	qaws_surface_bspline_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_bspline_impl* impl;
	size_t cp_size;
	qaws_range u_range, v_range;
	qaws_scalar* u_knots = NULL;
	qaws_scalar* v_knots = NULL;
	unsigned int u_knot_count, v_knot_count;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->u_degree < 1 || desc->v_degree < 1) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->u_point_count < desc->u_degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (desc->v_point_count < desc->v_degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points) return QAWS_STATUS_INVALID_ARGUMENT;

	/* Generate or copy U knots */
	if (desc->u_knot_count == 0 || !desc->u_knots)
	{
		u_knot_count = desc->u_point_count + desc->u_degree + 1;
		u_knots = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)(sizeof(qaws_scalar) * u_knot_count));
		if (!u_knots) return QAWS_STATUS_ALLOCATION_FAILURE;
		qaws_internal_surface_uniform_knots(
			desc->u_degree, desc->u_point_count, u_knots, u_knot_count);
	}
	else
	{
		u_knot_count = desc->u_knot_count;
		if (u_knot_count != desc->u_point_count + desc->u_degree + 1)
		{
			return QAWS_STATUS_INVALID_KNOT_VECTOR;
		}
		u_knots = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)(sizeof(qaws_scalar) * u_knot_count));
		if (!u_knots) return QAWS_STATUS_ALLOCATION_FAILURE;
		memcpy(u_knots, desc->u_knots, sizeof(qaws_scalar) * u_knot_count);
	}

	/* Generate or copy V knots */
	if (desc->v_knot_count == 0 || !desc->v_knots)
	{
		v_knot_count = desc->v_point_count + desc->v_degree + 1;
		v_knots = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)(sizeof(qaws_scalar) * v_knot_count));
		if (!v_knots) { qaws_internal_dealloc(allocator, u_knots); return QAWS_STATUS_ALLOCATION_FAILURE; }
		qaws_internal_surface_uniform_knots(
			desc->v_degree, desc->v_point_count, v_knots, v_knot_count);
	}
	else
	{
		v_knot_count = desc->v_knot_count;
		if (v_knot_count != desc->v_point_count + desc->v_degree + 1)
		{
			qaws_internal_dealloc(allocator, u_knots);
			return QAWS_STATUS_INVALID_KNOT_VECTOR;
		}
		v_knots = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)(sizeof(qaws_scalar) * v_knot_count));
		if (!v_knots) { qaws_internal_dealloc(allocator, u_knots); return QAWS_STATUS_ALLOCATION_FAILURE; }
		memcpy(v_knots, desc->v_knots, sizeof(qaws_scalar) * v_knot_count);
	}

	/* Domain range */
	u_range.min_value = u_knots[desc->u_degree];
	u_range.max_value = u_knots[desc->u_point_count];
	v_range.min_value = v_knots[desc->v_degree];
	v_range.max_value = v_knots[desc->v_point_count];

	surface = qaws_internal_surface_alloc_ex(
		QAWS_SURFACE_KIND_BSPLINE,
		desc->u_degree, desc->v_degree,
		u_range, v_range,
		&bspline_surface_vtable,
		allocator);
	if (!surface)
	{
		qaws_internal_dealloc(allocator, u_knots); qaws_internal_dealloc(allocator, v_knots);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl = (qaws_surface_bspline_impl*)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_surface_bspline_impl));
	if (!impl)
	{
		qaws_internal_dealloc(allocator, u_knots); qaws_internal_dealloc(allocator, v_knots);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	cp_size = (size_t)desc->u_point_count * desc->v_point_count * 3 * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)cp_size);
	if (!impl->control_points)
	{
		qaws_internal_dealloc(allocator, u_knots); qaws_internal_dealloc(allocator, v_knots); qaws_internal_dealloc(allocator, impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy control points from vec3 array to flat scalar array */
	{
		unsigned int total = desc->u_point_count * desc->v_point_count;
		unsigned int k;
		for (k = 0; k < total; k++)
		{
			impl->control_points[k * 3 + 0] = desc->control_points[k].x;
			impl->control_points[k * 3 + 1] = desc->control_points[k].y;
			impl->control_points[k * 3 + 2] = desc->control_points[k].z;
		}
	}

	impl->u_count = desc->u_point_count;
	impl->v_count = desc->v_point_count;
	impl->u_knots = u_knots;
	impl->u_knot_count = u_knot_count;
	impl->v_knots = v_knots;
	impl->v_knot_count = v_knot_count;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_create_bspline(
	qaws_surface_bspline_desc const* desc,
	qaws_surface** out_surface)
{
	return qaws_surface_create_bspline_ex(desc, NULL, out_surface);
}
