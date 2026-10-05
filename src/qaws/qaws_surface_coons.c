#include "qaws_surface_coons.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_coons_impl
{
	qaws_curve const* c0;
	qaws_curve const* c1;
	qaws_curve const* d0;
	qaws_curve const* d1;
	qaws_range range_c0;
	qaws_range range_c1;
	qaws_range range_d0;
	qaws_range range_d1;
	qaws_vec3 P00;  /* c0(0) = d0(0) */
	qaws_vec3 P10;  /* c0(1) = d1(0) */
	qaws_vec3 P01;  /* c1(0) = d0(1) */
	qaws_vec3 P11;  /* c1(1) = d1(1) */
} qaws_surface_coons_impl;

/* Coons patch: S(u,v) = Lc(u,v) + Ld(u,v) - B(u,v)
   Lc(u,v) = (1-v)*c0(u) + v*c1(u)
   Ld(u,v) = (1-u)*d0(v) + u*d1(v)
   B(u,v)  = (1-u)*(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11

   dS/du = (1-v)*c0'(u)*s_c0 + v*c1'(u)*s_c1
         - d0(v) + d1(v)
         - [-(1-v)*P00 + (1-v)*P10 - v*P01 + v*P11]

   dS/dv = -c0(u) + c1(u)
         + (1-u)*d0'(v)*s_d0 + u*d1'(v)*s_d1
         - [(1-u)*(-P00+P01) + u*(-P10+P11)]

   Second derivatives use central finite differences. */
static qaws_status coons_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_coons_impl const* impl =
		(qaws_surface_coons_impl const*)surface->impl;
	qaws_scalar one_minus_u = QAWS_ONE - u;
	qaws_scalar one_minus_v = QAWS_ONE - v;
	qaws_scalar s_c0 = impl->range_c0.max_value - impl->range_c0.min_value;
	qaws_scalar s_c1 = impl->range_c1.max_value - impl->range_c1.min_value;
	qaws_scalar s_d0 = impl->range_d0.max_value - impl->range_d0.min_value;
	qaws_scalar s_d1 = impl->range_d1.max_value - impl->range_d1.min_value;
	qaws_scalar t_c0 = impl->range_c0.min_value + u * s_c0;
	qaws_scalar t_c1 = impl->range_c1.min_value + u * s_c1;
	qaws_scalar t_d0 = impl->range_d0.min_value + v * s_d0;
	qaws_scalar t_d1 = impl->range_d1.min_value + v * s_d1;

	unsigned int c_flags = 0;
	unsigned int d_flags = 0;
	qaws_eval_result_3d rc0, rc1, rd0, rd1;
	qaws_status status;
	qaws_vec3 pos;

	/* Determine what we need from the curves */
	if (eval_flags & (QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DV))
		c_flags |= QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & (QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU))
		d_flags |= QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
		c_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
		d_flags |= QAWS_EVAL_FLAG_D1;

	/* For second derivatives and normal, we need position + D1 from all */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		c_flags |= QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
		d_flags |= QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
	}

	memset(&rc0, 0, sizeof(rc0));
	memset(&rc1, 0, sizeof(rc1));
	memset(&rd0, 0, sizeof(rd0));
	memset(&rd1, 0, sizeof(rd1));

	status = qaws_curve_evaluate_3d(impl->c0, t_c0, c_flags, &rc0);
	if (status != QAWS_STATUS_OK) return status;
	status = qaws_curve_evaluate_3d(impl->c1, t_c1, c_flags, &rc1);
	if (status != QAWS_STATUS_OK) return status;
	status = qaws_curve_evaluate_3d(impl->d0, t_d0, d_flags, &rd0);
	if (status != QAWS_STATUS_OK) return status;
	status = qaws_curve_evaluate_3d(impl->d1, t_d1, d_flags, &rd1);
	if (status != QAWS_STATUS_OK) return status;

	/* Position: (1-v)*c0 + v*c1 + (1-u)*d0 + u*d1
	   - [(1-u)*(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11] */
	pos.x = one_minus_v * rc0.position.x + v * rc1.position.x
		+ one_minus_u * rd0.position.x + u * rd1.position.x
		- (one_minus_u * one_minus_v * impl->P00.x
		   + u * one_minus_v * impl->P10.x
		   + one_minus_u * v * impl->P01.x
		   + u * v * impl->P11.x);
	pos.y = one_minus_v * rc0.position.y + v * rc1.position.y
		+ one_minus_u * rd0.position.y + u * rd1.position.y
		- (one_minus_u * one_minus_v * impl->P00.y
		   + u * one_minus_v * impl->P10.y
		   + one_minus_u * v * impl->P01.y
		   + u * v * impl->P11.y);
	pos.z = one_minus_v * rc0.position.z + v * rc1.position.z
		+ one_minus_u * rd0.position.z + u * rd1.position.z
		- (one_minus_u * one_minus_v * impl->P00.z
		   + u * one_minus_v * impl->P10.z
		   + one_minus_u * v * impl->P01.z
		   + u * v * impl->P11.z);

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = (1-v)*c0'(u)*s_c0 + v*c1'(u)*s_c1
	         - d0(v) + d1(v)
	         - [-(1-v)*P00 + (1-v)*P10 - v*P01 + v*P11] */
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = one_minus_v * rc0.d1.x * s_c0 + v * rc1.d1.x * s_c1
			- rd0.position.x + rd1.position.x
			- (-(one_minus_v) * impl->P00.x + one_minus_v * impl->P10.x
			   - v * impl->P01.x + v * impl->P11.x);
		out_result->du.y = one_minus_v * rc0.d1.y * s_c0 + v * rc1.d1.y * s_c1
			- rd0.position.y + rd1.position.y
			- (-(one_minus_v) * impl->P00.y + one_minus_v * impl->P10.y
			   - v * impl->P01.y + v * impl->P11.y);
		out_result->du.z = one_minus_v * rc0.d1.z * s_c0 + v * rc1.d1.z * s_c1
			- rd0.position.z + rd1.position.z
			- (-(one_minus_v) * impl->P00.z + one_minus_v * impl->P10.z
			   - v * impl->P01.z + v * impl->P11.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = -c0(u) + c1(u)
	         + (1-u)*d0'(v)*s_d0 + u*d1'(v)*s_d1
	         - [(1-u)*(-P00+P01) + u*(-P10+P11)] */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = -rc0.position.x + rc1.position.x
			+ one_minus_u * rd0.d1.x * s_d0 + u * rd1.d1.x * s_d1
			- (one_minus_u * (-impl->P00.x + impl->P01.x)
			   + u * (-impl->P10.x + impl->P11.x));
		out_result->dv.y = -rc0.position.y + rc1.position.y
			+ one_minus_u * rd0.d1.y * s_d0 + u * rd1.d1.y * s_d1
			- (one_minus_u * (-impl->P00.y + impl->P01.y)
			   + u * (-impl->P10.y + impl->P11.y));
		out_result->dv.z = -rc0.position.z + rc1.position.z
			+ one_minus_u * rd0.d1.z * s_d0 + u * rd1.d1.z * s_d1
			- (one_minus_u * (-impl->P00.z + impl->P01.z)
			   + u * (-impl->P10.z + impl->P11.z));
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* Second derivatives via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_surface_eval_result r_lo, r_hi;

		/* duu via finite difference of du w.r.t. u */
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			coons_surface_eval(surface, u_lo, v, QAWS_SURFACE_EVAL_DU, &r_lo);
			coons_surface_eval(surface, u_hi, v, QAWS_SURFACE_EVAL_DU, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duu.x = (r_hi.du.x - r_lo.du.x) * inv2h;
				out_result->duu.y = (r_hi.du.y - r_lo.du.y) * inv2h;
				out_result->duu.z = (r_hi.du.z - r_lo.du.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		/* dvv via finite difference of dv w.r.t. v */
		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		{
			qaws_scalar v_lo = v - h, v_hi = v + h;
			qaws_scalar hv;
			if (v_lo < 0) v_lo = 0;
			if (v_hi > 1) v_hi = 1;
			hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			coons_surface_eval(surface, u, v_lo, QAWS_SURFACE_EVAL_DV, &r_lo);
			coons_surface_eval(surface, u, v_hi, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
				out_result->dvv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->dvv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->dvv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}

		/* duv via finite difference of dv w.r.t. u */
		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			coons_surface_eval(surface, u_lo, v, QAWS_SURFACE_EVAL_DV, &r_lo);
			coons_surface_eval(surface, u_hi, v, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->duv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->duv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
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

static void coons_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int coons_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const coons_surface_vtable = {
	coons_surface_eval,
	coons_surface_destroy,
	coons_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_coons(
	qaws_surface_coons_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_coons_impl* impl;
	qaws_range u_range, v_range;
	qaws_eval_result_3d corner_r;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->c0 || !desc->c1 || !desc->d0 || !desc->d1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(desc->c0) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->c1) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->d0) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->d1) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_COONS,
		0, 0, u_range, v_range,
		&coons_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_coons_impl*)malloc(sizeof(qaws_surface_coons_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->c0 = desc->c0;
	impl->c1 = desc->c1;
	impl->d0 = desc->d0;
	impl->d1 = desc->d1;
	impl->range_c0 = qaws_curve_get_parameter_range(desc->c0);
	impl->range_c1 = qaws_curve_get_parameter_range(desc->c1);
	impl->range_d0 = qaws_curve_get_parameter_range(desc->d0);
	impl->range_d1 = qaws_curve_get_parameter_range(desc->d1);

	/* Cache corner points: P00=c0(start), P10=c0(end), P01=c1(start), P11=c1(end) */
	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c0, impl->range_c0.min_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P00 = corner_r.position;

	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c0, impl->range_c0.max_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P10 = corner_r.position;

	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c1, impl->range_c1.min_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P01 = corner_r.position;

	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c1, impl->range_c1.max_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P11 = corner_r.position;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
