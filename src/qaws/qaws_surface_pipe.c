#include "qaws_surface_pipe.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_pipe_impl
{
	qaws_curve const* path;
	qaws_range path_range;
	qaws_scalar radius_x;
	qaws_scalar radius_y;
} qaws_surface_pipe_impl;

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

/* Pipe surface: S(u,v) = P(u) + rx * cos(theta) * N(u) + ry * sin(theta) * B(u)
   where theta = 2 * pi * v.
   u follows the path, v parameterizes the circular/elliptical cross-section.

   dS/dv is computed analytically.
   dS/du uses central finite differences (frame derivatives are complex). */
static qaws_status pipe_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_pipe_impl const* impl =
		(qaws_surface_pipe_impl const*)surface->impl;
	qaws_scalar s_path = impl->path_range.max_value - impl->path_range.min_value;
	qaws_scalar t_path = impl->path_range.min_value + u * s_path;
	qaws_scalar rx = impl->radius_x;
	qaws_scalar ry = impl->radius_y;
	qaws_scalar pi = QAWS_LITERAL(3.14159265358979323846);
	qaws_scalar two_pi = QAWS_LITERAL(2.0) * pi;
	qaws_scalar theta = two_pi * v;
	qaws_scalar cos_theta = QAWS_COS(theta);
	qaws_scalar sin_theta = QAWS_SIN(theta);

	qaws_eval_result_3d path_r;
	qaws_vec3 T, N, B, pos;

	/* Evaluate path position */
	memset(&path_r, 0, sizeof(path_r));
	{
		qaws_status s = qaws_curve_evaluate_3d(impl->path, t_path,
			QAWS_EVAL_FLAG_POSITION, &path_r);
		if (s != QAWS_STATUS_OK) return s;
	}

	/* Get Frenet frame at path point */
	{
		qaws_status s = qaws_curve_compute_frenet_frame_3d(
			impl->path, t_path, &T, &N, &B);
		if (s != QAWS_STATUS_OK) return s;
	}

	/* Position: P(u) + rx * cos(theta) * N(u) + ry * sin(theta) * B(u) */
	pos.x = path_r.position.x + rx * cos_theta * N.x + ry * sin_theta * B.x;
	pos.y = path_r.position.y + rx * cos_theta * N.y + ry * sin_theta * B.y;
	pos.z = path_r.position.z + rx * cos_theta * N.z + ry * sin_theta * B.z;

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/dv = 2*pi * (-rx * sin(theta) * N(u) + ry * cos(theta) * B(u)) */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		qaws_scalar dnx = -rx * sin_theta;
		qaws_scalar dny = ry * cos_theta;
		out_result->dv.x = two_pi * (dnx * N.x + dny * B.x);
		out_result->dv.y = two_pi * (dnx * N.y + dny * B.y);
		out_result->dv.z = two_pi * (dnx * N.z + dny * B.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* d2S/dv2 = (2*pi)^2 * (-rx * cos(theta) * N(u) - ry * sin(theta) * B(u)) */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		qaws_scalar tp2 = two_pi * two_pi;
		qaws_scalar ddnx = -rx * cos_theta;
		qaws_scalar ddny = -ry * sin_theta;
		out_result->dvv.x = tp2 * (ddnx * N.x + ddny * B.x);
		out_result->dvv.y = tp2 * (ddnx * N.y + ddny * B.y);
		out_result->dvv.z = tp2 * (ddnx * N.z + ddny * B.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* dS/du via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_scalar u_lo = u - h, u_hi = u + h;
		qaws_vec3 p_lo, p_hi;

		if (u_lo < 0) u_lo = 0;
		if (u_hi > 1) u_hi = 1;
		h = (u_hi - u_lo) * QAWS_LITERAL(0.5);

		/* Evaluate position at u-h and u+h */
		{
			qaws_scalar t_lo = impl->path_range.min_value + u_lo * s_path;
			qaws_scalar t_hi = impl->path_range.min_value + u_hi * s_path;
			qaws_eval_result_3d pr;
			qaws_vec3 T_lo, N_lo, B_lo, T_hi, N_hi, B_hi;

			memset(&pr, 0, sizeof(pr));
			qaws_curve_evaluate_3d(impl->path, t_lo, QAWS_EVAL_FLAG_POSITION, &pr);
			qaws_curve_compute_frenet_frame_3d(impl->path, t_lo, &T_lo, &N_lo, &B_lo);
			p_lo.x = pr.position.x + rx * cos_theta * N_lo.x + ry * sin_theta * B_lo.x;
			p_lo.y = pr.position.y + rx * cos_theta * N_lo.y + ry * sin_theta * B_lo.y;
			p_lo.z = pr.position.z + rx * cos_theta * N_lo.z + ry * sin_theta * B_lo.z;

			memset(&pr, 0, sizeof(pr));
			qaws_curve_evaluate_3d(impl->path, t_hi, QAWS_EVAL_FLAG_POSITION, &pr);
			qaws_curve_compute_frenet_frame_3d(impl->path, t_hi, &T_hi, &N_hi, &B_hi);
			p_hi.x = pr.position.x + rx * cos_theta * N_hi.x + ry * sin_theta * B_hi.x;
			p_hi.y = pr.position.y + rx * cos_theta * N_hi.y + ry * sin_theta * B_hi.y;
			p_hi.z = pr.position.z + rx * cos_theta * N_hi.z + ry * sin_theta * B_hi.z;
		}

		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * h);
			out_result->du.x = (p_hi.x - p_lo.x) * inv2h;
			out_result->du.y = (p_hi.y - p_lo.y) * inv2h;
			out_result->du.z = (p_hi.z - p_lo.z) * inv2h;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}

		/* duu via central 2nd-order finite difference: (p_hi - 2*p + p_lo) / h^2 */
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar inv_h2 = QAWS_ONE / (h * h);
			out_result->duu.x = (p_hi.x - 2 * pos.x + p_lo.x) * inv_h2;
			out_result->duu.y = (p_hi.y - 2 * pos.y + p_lo.y) * inv_h2;
			out_result->duu.z = (p_hi.z - 2 * pos.z + p_lo.z) * inv_h2;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		/* duv via finite difference of dv w.r.t. u */
		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			qaws_scalar t_lo = impl->path_range.min_value + u_lo * s_path;
			qaws_scalar t_hi = impl->path_range.min_value + u_hi * s_path;
			qaws_vec3 T_tmp, N_lo, B_lo, N_hi, B_hi;
			qaws_scalar dnx = -rx * sin_theta;
			qaws_scalar dny = ry * cos_theta;
			qaws_vec3 dv_lo, dv_hi;
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * h);

			qaws_curve_compute_frenet_frame_3d(impl->path, t_lo, &T_tmp, &N_lo, &B_lo);
			qaws_curve_compute_frenet_frame_3d(impl->path, t_hi, &T_tmp, &N_hi, &B_hi);

			dv_lo.x = two_pi * (dnx * N_lo.x + dny * B_lo.x);
			dv_lo.y = two_pi * (dnx * N_lo.y + dny * B_lo.y);
			dv_lo.z = two_pi * (dnx * N_lo.z + dny * B_lo.z);

			dv_hi.x = two_pi * (dnx * N_hi.x + dny * B_hi.x);
			dv_hi.y = two_pi * (dnx * N_hi.y + dny * B_hi.y);
			dv_hi.z = two_pi * (dnx * N_hi.z + dny * B_hi.z);

			out_result->duv.x = (dv_hi.x - dv_lo.x) * inv2h;
			out_result->duv.y = (dv_hi.y - dv_lo.y) * inv2h;
			out_result->duv.z = (dv_hi.z - dv_lo.z) * inv2h;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
		}
	}

	/* Normal */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		compute_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

static void pipe_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int pipe_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const pipe_surface_vtable = {
	pipe_surface_eval,
	pipe_surface_destroy,
	pipe_surface_is_rational
};

qaws_status qaws_surface_create_pipe(
	qaws_surface_pipe_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_pipe_impl* impl;
	qaws_range u_range, v_range;
	qaws_scalar ry;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->path) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(desc->path) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (desc->radius_x <= 0) return QAWS_STATUS_INVALID_ARGUMENT;

	ry = (desc->radius_y > 0) ? desc->radius_y : desc->radius_x;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_PIPE,
		0, 0, u_range, v_range,
		&pipe_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_pipe_impl*)malloc(sizeof(qaws_surface_pipe_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->path = desc->path;
	impl->path_range = qaws_curve_get_parameter_range(desc->path);
	impl->radius_x = desc->radius_x;
	impl->radius_y = ry;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
