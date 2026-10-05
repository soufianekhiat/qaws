#include "qaws_surface_tspline.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_tspline_impl
{
	qaws_vec3* positions;        /* cp_count positions */
	qaws_scalar* weights;        /* cp_count weights */
	qaws_scalar* u_knots_all;    /* cp_count * u_knot_len knots (contiguous) */
	qaws_scalar* v_knots_all;    /* cp_count * v_knot_len knots (contiguous) */
	unsigned int cp_count;
	unsigned int u_degree;
	unsigned int v_degree;
	unsigned int u_knot_len;     /* 2*u_degree+2 */
	unsigned int v_knot_len;     /* 2*v_degree+2 */
} qaws_surface_tspline_impl;

/* Cox-de Boor recursion for a single B-spline basis function.
   knots: local knot vector of length (degree+2).
   Returns N_{0,degree}(t) defined over knots[0..degree+1]. */
static qaws_scalar eval_basis(
	qaws_scalar const* knots,
	unsigned int degree,
	qaws_scalar t)
{
	qaws_scalar N[16]; /* max degree 15 */
	unsigned int d, j;
	unsigned int n = degree + 1; /* number of intervals */

	/* Degree 0: piecewise constant */
	for (j = 0; j < n; j++)
	{
		if (knots[j] < knots[j + 1])
			N[j] = (t >= knots[j] && t < knots[j + 1]) ? QAWS_ONE : QAWS_ZERO;
		else
			N[j] = QAWS_ZERO;
	}
	/* Special case: include right endpoint */
	if (t >= knots[n] && QAWS_FABS(t - knots[n]) < QAWS_LITERAL(1e-10))
	{
		for (j = n; j > 0; j--)
		{
			if (knots[j - 1] < knots[j])
			{
				N[j - 1] = QAWS_ONE;
				break;
			}
		}
	}

	/* Build up degrees */
	for (d = 1; d <= degree; d++)
	{
		for (j = 0; j + d < n; j++)
		{
			qaws_scalar left = QAWS_ZERO, right = QAWS_ZERO;
			qaws_scalar denom;

			denom = knots[j + d] - knots[j];
			if (denom > QAWS_LITERAL(1e-14))
				left = (t - knots[j]) / denom * N[j];

			denom = knots[j + d + 1] - knots[j + 1];
			if (denom > QAWS_LITERAL(1e-14))
				right = (knots[j + d + 1] - t) / denom * N[j + 1];

			N[j] = left + right;
		}
	}

	return N[0];
}

/* Derivative of B-spline basis function.
   d/dt N_{0,p}(t) = p * [N_{0,p-1}(t)/(knots[p]-knots[0])
                         - N_{1,p-1}(t)/(knots[p+1]-knots[1])]
   where N_{0,p-1} uses knots[0..p] and N_{1,p-1} uses knots[1..p+1]. */
static qaws_scalar eval_basis_deriv(
	qaws_scalar const* knots,
	unsigned int degree,
	qaws_scalar t)
{
	qaws_scalar left = QAWS_ZERO, right = QAWS_ZERO;
	qaws_scalar denom;

	if (degree == 0) return QAWS_ZERO;

	denom = knots[degree] - knots[0];
	if (denom > QAWS_LITERAL(1e-14))
		left = eval_basis(knots, degree - 1, t) / denom;

	denom = knots[degree + 1] - knots[1];
	if (denom > QAWS_LITERAL(1e-14))
		right = eval_basis(knots + 1, degree - 1, t) / denom;

	return (qaws_scalar)degree * (left - right);
}

static qaws_status tspline_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_tspline_impl const* impl =
		(qaws_surface_tspline_impl const*)surface->impl;
	unsigned int i;

	qaws_scalar sum_wB = QAWS_ZERO;
	qaws_vec3 sum_wBP = {0, 0, 0};

	qaws_scalar sum_wB_du = QAWS_ZERO, sum_wB_dv = QAWS_ZERO;
	qaws_vec3 sum_wBP_du = {0, 0, 0}, sum_wBP_dv = {0, 0, 0};

	int need_deriv = (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV |
	                                QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DUU |
	                                QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV));

	qaws_scalar inv_wB;

	for (i = 0; i < impl->cp_count; i++)
	{
		qaws_scalar const* uk = impl->u_knots_all + i * impl->u_knot_len;
		qaws_scalar const* vk = impl->v_knots_all + i * impl->v_knot_len;

		qaws_scalar Nu = eval_basis(uk, impl->u_degree, u);
		qaws_scalar Nv = eval_basis(vk, impl->v_degree, v);
		qaws_scalar B = Nu * Nv;
		qaws_scalar wB = impl->weights[i] * B;

		sum_wB += wB;
		sum_wBP.x += wB * impl->positions[i].x;
		sum_wBP.y += wB * impl->positions[i].y;
		sum_wBP.z += wB * impl->positions[i].z;

		if (need_deriv)
		{
			qaws_scalar dNu = eval_basis_deriv(uk, impl->u_degree, u);
			qaws_scalar dNv = eval_basis_deriv(vk, impl->v_degree, v);
			qaws_scalar dB_du = dNu * Nv;
			qaws_scalar dB_dv = Nu * dNv;
			qaws_scalar wdB_du = impl->weights[i] * dB_du;
			qaws_scalar wdB_dv = impl->weights[i] * dB_dv;

			sum_wB_du += wdB_du;
			sum_wB_dv += wdB_dv;
			sum_wBP_du.x += wdB_du * impl->positions[i].x;
			sum_wBP_du.y += wdB_du * impl->positions[i].y;
			sum_wBP_du.z += wdB_du * impl->positions[i].z;
			sum_wBP_dv.x += wdB_dv * impl->positions[i].x;
			sum_wBP_dv.y += wdB_dv * impl->positions[i].y;
			sum_wBP_dv.z += wdB_dv * impl->positions[i].z;
		}
	}

	/* Avoid division by zero */
	if (QAWS_FABS(sum_wB) < QAWS_LITERAL(1e-14))
		sum_wB = QAWS_LITERAL(1e-14);

	inv_wB = QAWS_ONE / sum_wB;

	/* Position: S = sum(wBP) / sum(wB) */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = sum_wBP.x * inv_wB;
		out_result->position.y = sum_wBP.y * inv_wB;
		out_result->position.z = sum_wBP.z * inv_wB;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* Derivatives: quotient rule d/du [f/g] = (f'g - fg') / g^2 */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar inv_wB2 = inv_wB * inv_wB;
		out_result->du.x = (sum_wBP_du.x * sum_wB - sum_wBP.x * sum_wB_du) * inv_wB2;
		out_result->du.y = (sum_wBP_du.y * sum_wB - sum_wBP.y * sum_wB_du) * inv_wB2;
		out_result->du.z = (sum_wBP_du.z * sum_wB - sum_wBP.z * sum_wB_du) * inv_wB2;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar inv_wB2 = inv_wB * inv_wB;
		out_result->dv.x = (sum_wBP_dv.x * sum_wB - sum_wBP.x * sum_wB_dv) * inv_wB2;
		out_result->dv.y = (sum_wBP_dv.y * sum_wB - sum_wBP.y * sum_wB_dv) * inv_wB2;
		out_result->dv.z = (sum_wBP_dv.z * sum_wB - sum_wBP.z * sum_wB_dv) * inv_wB2;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* Second derivatives via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_surface_eval_result rp, rm;

		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			memset(&rp, 0, sizeof(rp));
			memset(&rm, 0, sizeof(rm));
			tspline_surface_eval(surface, u + h, v, QAWS_SURFACE_EVAL_DU, &rp);
			tspline_surface_eval(surface, u - h, v, QAWS_SURFACE_EVAL_DU, &rm);
			out_result->duu.x = (rp.du.x - rm.du.x) / (QAWS_LITERAL(2.0) * h);
			out_result->duu.y = (rp.du.y - rm.du.y) / (QAWS_LITERAL(2.0) * h);
			out_result->duu.z = (rp.du.z - rm.du.z) / (QAWS_LITERAL(2.0) * h);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}
		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		{
			memset(&rp, 0, sizeof(rp));
			memset(&rm, 0, sizeof(rm));
			tspline_surface_eval(surface, u, v + h, QAWS_SURFACE_EVAL_DV, &rp);
			tspline_surface_eval(surface, u, v - h, QAWS_SURFACE_EVAL_DV, &rm);
			out_result->dvv.x = (rp.dv.x - rm.dv.x) / (QAWS_LITERAL(2.0) * h);
			out_result->dvv.y = (rp.dv.y - rm.dv.y) / (QAWS_LITERAL(2.0) * h);
			out_result->dvv.z = (rp.dv.z - rm.dv.z) / (QAWS_LITERAL(2.0) * h);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}
		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			memset(&rp, 0, sizeof(rp));
			memset(&rm, 0, sizeof(rm));
			tspline_surface_eval(surface, u + h, v, QAWS_SURFACE_EVAL_DV, &rp);
			tspline_surface_eval(surface, u - h, v, QAWS_SURFACE_EVAL_DV, &rm);
			out_result->duv.x = (rp.dv.x - rm.dv.x) / (QAWS_LITERAL(2.0) * h);
			out_result->duv.y = (rp.dv.y - rm.dv.y) / (QAWS_LITERAL(2.0) * h);
			out_result->duv.z = (rp.dv.z - rm.dv.z) / (QAWS_LITERAL(2.0) * h);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
		}
	}

	/* Normal from cross product du x dv */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

static void tspline_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_tspline_impl* ti = (qaws_surface_tspline_impl*)impl;
	if (ti)
	{
		qaws_internal_dealloc(allocator, ti->positions);
		qaws_internal_dealloc(allocator, ti->weights);
		qaws_internal_dealloc(allocator, ti->u_knots_all);
		qaws_internal_dealloc(allocator, ti->v_knots_all);
		qaws_internal_dealloc(allocator, ti);
	}
}

static int tspline_surface_is_rational(qaws_surface const* s)
{
	qaws_surface_tspline_impl const* impl =
		(qaws_surface_tspline_impl const*)s->impl;
	unsigned int i;
	for (i = 0; i < impl->cp_count; i++)
	{
		if (QAWS_FABS(impl->weights[i] - QAWS_ONE) > QAWS_LITERAL(1e-12))
			return 1;
	}
	return 0;
}

static qaws_surface_vtable const tspline_surface_vtable = {
	tspline_surface_eval,
	tspline_surface_destroy,
	tspline_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_tspline(
	qaws_surface_tspline_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_tspline_impl* impl;
	unsigned int u_knot_len, v_knot_len;
	unsigned int i;
	size_t pos_size, wt_size, uk_size, vk_size;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->control_point_count == 0) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->u_degree < 1) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->v_degree < 1) return QAWS_STATUS_INVALID_DEGREE;

	u_knot_len = 2 * desc->u_degree + 2;
	v_knot_len = 2 * desc->v_degree + 2;

	/* Validate that all control points have non-NULL knot vectors */
	for (i = 0; i < desc->control_point_count; i++)
	{
		if (!desc->control_points[i].u_knots || !desc->control_points[i].v_knots)
			return QAWS_STATUS_INVALID_ARGUMENT;
	}

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_TSPLINE,
		desc->u_degree, desc->v_degree,
		desc->u_range, desc->v_range,
		&tspline_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_tspline_impl*)malloc(sizeof(qaws_surface_tspline_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->cp_count = desc->control_point_count;
	impl->u_degree = desc->u_degree;
	impl->v_degree = desc->v_degree;
	impl->u_knot_len = u_knot_len;
	impl->v_knot_len = v_knot_len;
	impl->positions = NULL;
	impl->weights = NULL;
	impl->u_knots_all = NULL;
	impl->v_knots_all = NULL;

	pos_size = sizeof(qaws_vec3) * impl->cp_count;
	wt_size = sizeof(qaws_scalar) * impl->cp_count;
	uk_size = sizeof(qaws_scalar) * impl->cp_count * u_knot_len;
	vk_size = sizeof(qaws_scalar) * impl->cp_count * v_knot_len;

	impl->positions = (qaws_vec3*)malloc(pos_size);
	impl->weights = (qaws_scalar*)malloc(wt_size);
	impl->u_knots_all = (qaws_scalar*)malloc(uk_size);
	impl->v_knots_all = (qaws_scalar*)malloc(vk_size);

	if (!impl->positions || !impl->weights || !impl->u_knots_all || !impl->v_knots_all)
	{
		free(impl->positions);
		free(impl->weights);
		free(impl->u_knots_all);
		free(impl->v_knots_all);
		free(impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy data from descriptor */
	for (i = 0; i < impl->cp_count; i++)
	{
		impl->positions[i] = desc->control_points[i].position;
		impl->weights[i] = desc->control_points[i].weight;
		memcpy(impl->u_knots_all + i * u_knot_len,
		       desc->control_points[i].u_knots,
		       sizeof(qaws_scalar) * u_knot_len);
		memcpy(impl->v_knots_all + i * v_knot_len,
		       desc->control_points[i].v_knots,
		       sizeof(qaws_scalar) * v_knot_len);
	}

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
