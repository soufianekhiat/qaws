#include "qaws_surface_nurbs.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* Evaluate a tensor-product NURBS surface at (u,v) using the quotient rule.
   S(u,v) = A(u,v) / W(u,v), where:
     A(u,v) = sum_ij N_i(u) N_j(v) w_ij P_ij
     W(u,v) = sum_ij N_i(u) N_j(v) w_ij
   Derivatives use the quotient rule:
     dS/du = (dA/du - S * dW/du) / W  etc. */
static qaws_status nurbs_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_nurbs_impl const* impl =
		(qaws_surface_nurbs_impl const*)surface->impl;
	unsigned int ud = surface->u_degree;
	unsigned int vd = surface->v_degree;
	unsigned int uc = impl->u_count;
	unsigned int vc = impl->v_count;
	qaws_scalar const* cp = impl->control_points;
	qaws_scalar const* wt = impl->weights;
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

	{
		qaws_scalar u_ders[3 * 64];
		qaws_scalar v_ders[3 * 64];
		unsigned int ku, kv, i, j;
		/* A[ku][kv] = weighted point accumulator, W[ku][kv] = weight accumulator */
		qaws_scalar A[3][3][3];
		qaws_scalar W[3][3];
		qaws_scalar inv_w0;
		qaws_vec3 S00; /* position */
		qaws_vec3 Su, Sv; /* first partials */
		qaws_vec3 Suu, Svv, Suv; /* second partials */

		qaws_internal_bspline_basis_derivs(
			impl->u_knots, impl->u_knot_count,
			ud, u_span, u, max_u_deriv, u_ders);
		qaws_internal_bspline_basis_derivs(
			impl->v_knots, impl->v_knot_count,
			vd, v_span, v, max_v_deriv, v_ders);

		memset(A, 0, sizeof(A));
		memset(W, 0, sizeof(W));

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
						qaws_scalar w = wt[ci * vc + cj];
						qaws_scalar Nw = Nu * Nv * w;
						unsigned int cp_idx = (ci * vc + cj) * 3;
						W[ku][kv] += Nw;
						A[ku][kv][0] += Nw * cp[cp_idx + 0];
						A[ku][kv][1] += Nw * cp[cp_idx + 1];
						A[ku][kv][2] += Nw * cp[cp_idx + 2];
					}
				}
			}
		}

		/* Check for degenerate weight */
		if (W[0][0] < QAWS_LITERAL(1e-15) && W[0][0] > -QAWS_LITERAL(1e-15))
			return QAWS_STATUS_DEGENERATE_CURVE;

		inv_w0 = QAWS_ONE / W[0][0];

		/* Position: S = A[0][0] / W[0][0] */
		S00.x = A[0][0][0] * inv_w0;
		S00.y = A[0][0][1] * inv_w0;
		S00.z = A[0][0][2] * inv_w0;

		if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
		{
			out_result->position = S00;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
		}

		/* First partial du: Su = (A[1][0] - S * W[1][0]) / W[0][0] */
		Su.x = Su.y = Su.z = 0;
		if (max_u_deriv >= 1)
		{
			Su.x = (A[1][0][0] - S00.x * W[1][0]) * inv_w0;
			Su.y = (A[1][0][1] - S00.y * W[1][0]) * inv_w0;
			Su.z = (A[1][0][2] - S00.z * W[1][0]) * inv_w0;
		}
		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			out_result->du = Su;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}

		/* First partial dv: Sv = (A[0][1] - S * W[0][1]) / W[0][0] */
		Sv.x = Sv.y = Sv.z = 0;
		if (max_v_deriv >= 1)
		{
			Sv.x = (A[0][1][0] - S00.x * W[0][1]) * inv_w0;
			Sv.y = (A[0][1][1] - S00.y * W[0][1]) * inv_w0;
			Sv.z = (A[0][1][2] - S00.z * W[0][1]) * inv_w0;
		}
		if (eval_flags & QAWS_SURFACE_EVAL_DV)
		{
			out_result->dv = Sv;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
		}

		/* Second partial duu: Suu = (A[2][0] - 2*Su*W[1][0] - S*W[2][0]) / W[0][0] */
		Suu.x = Suu.y = Suu.z = 0;
		if ((eval_flags & QAWS_SURFACE_EVAL_DUU) && max_u_deriv >= 2)
		{
			Suu.x = (A[2][0][0] - QAWS_LITERAL(2.0) * Su.x * W[1][0] - S00.x * W[2][0]) * inv_w0;
			Suu.y = (A[2][0][1] - QAWS_LITERAL(2.0) * Su.y * W[1][0] - S00.y * W[2][0]) * inv_w0;
			Suu.z = (A[2][0][2] - QAWS_LITERAL(2.0) * Su.z * W[1][0] - S00.z * W[2][0]) * inv_w0;
			out_result->duu = Suu;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		/* Second partial dvv: Svv = (A[0][2] - 2*Sv*W[0][1] - S*W[0][2]) / W[0][0] */
		Svv.x = Svv.y = Svv.z = 0;
		if ((eval_flags & QAWS_SURFACE_EVAL_DVV) && max_v_deriv >= 2)
		{
			Svv.x = (A[0][2][0] - QAWS_LITERAL(2.0) * Sv.x * W[0][1] - S00.x * W[0][2]) * inv_w0;
			Svv.y = (A[0][2][1] - QAWS_LITERAL(2.0) * Sv.y * W[0][1] - S00.y * W[0][2]) * inv_w0;
			Svv.z = (A[0][2][2] - QAWS_LITERAL(2.0) * Sv.z * W[0][1] - S00.z * W[0][2]) * inv_w0;
			out_result->dvv = Svv;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}

		/* Mixed partial duv: Suv = (A[1][1] - Su*W[0][1] - Sv*W[1][0] - S*W[1][1]) / W[0][0] */
		Suv.x = Suv.y = Suv.z = 0;
		if ((eval_flags & QAWS_SURFACE_EVAL_DUV) && max_u_deriv >= 1 && max_v_deriv >= 1)
		{
			Suv.x = (A[1][1][0] - Su.x * W[0][1] - Sv.x * W[1][0] - S00.x * W[1][1]) * inv_w0;
			Suv.y = (A[1][1][1] - Su.y * W[0][1] - Sv.y * W[1][0] - S00.y * W[1][1]) * inv_w0;
			Suv.z = (A[1][1][2] - Su.z * W[0][1] - Sv.z * W[1][0] - S00.z * W[1][1]) * inv_w0;
			out_result->duv = Suv;
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

static void nurbs_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_nurbs_impl* ni = (qaws_surface_nurbs_impl*)impl;
	if (ni)
	{
		qaws_internal_dealloc(allocator, ni->control_points);
		qaws_internal_dealloc(allocator, ni->weights);
		qaws_internal_dealloc(allocator, ni->u_knots);
		qaws_internal_dealloc(allocator, ni->v_knots);
		qaws_internal_dealloc(allocator, ni);
	}
}

static int nurbs_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 1;
}

static qaws_surface_vtable const nurbs_surface_vtable = {
	nurbs_surface_eval,
	nurbs_surface_destroy,
	nurbs_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_nurbs_ex(
	qaws_surface_nurbs_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_nurbs_impl* impl;
	size_t cp_size;
	qaws_range u_range, v_range;
	qaws_scalar* u_knots = NULL;
	qaws_scalar* v_knots = NULL;
	unsigned int u_knot_count, v_knot_count;
	unsigned int total;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->u_degree < 1 || desc->v_degree < 1) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->u_point_count < desc->u_degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (desc->v_point_count < desc->v_degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points || !desc->weights)
		return QAWS_STATUS_INVALID_ARGUMENT;

	total = desc->u_point_count * desc->v_point_count;

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
			return QAWS_STATUS_INVALID_KNOT_VECTOR;
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
		QAWS_SURFACE_KIND_NURBS,
		desc->u_degree, desc->v_degree,
		u_range, v_range,
		&nurbs_surface_vtable,
		allocator);
	if (!surface)
	{
		qaws_internal_dealloc(allocator, u_knots); qaws_internal_dealloc(allocator, v_knots);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl = (qaws_surface_nurbs_impl*)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_surface_nurbs_impl));
	if (!impl)
	{
		qaws_internal_dealloc(allocator, u_knots); qaws_internal_dealloc(allocator, v_knots);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	cp_size = (size_t)total * 3 * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)cp_size);
	if (!impl->control_points)
	{
		qaws_internal_dealloc(allocator, u_knots); qaws_internal_dealloc(allocator, v_knots); qaws_internal_dealloc(allocator, impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy control points from vec3 to flat scalar array */
	{
		unsigned int k;
		for (k = 0; k < total; k++)
		{
			impl->control_points[k * 3 + 0] = desc->control_points[k].x;
			impl->control_points[k * 3 + 1] = desc->control_points[k].y;
			impl->control_points[k * 3 + 2] = desc->control_points[k].z;
		}
	}

	/* Copy weights */
	impl->weights = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)(sizeof(qaws_scalar) * total));
	if (!impl->weights)
	{
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, u_knots); qaws_internal_dealloc(allocator, v_knots); qaws_internal_dealloc(allocator, impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->weights, desc->weights, sizeof(qaws_scalar) * total);

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

qaws_status qaws_surface_create_nurbs(
	qaws_surface_nurbs_desc const* desc,
	qaws_surface** out_surface)
{
	return qaws_surface_create_nurbs_ex(desc, NULL, out_surface);
}
