#include "qaws_surface_swept.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_swept_impl
{
	qaws_curve const* path;
	qaws_curve const* profile;
	qaws_scalar scale;
	qaws_range path_range;
	qaws_range prof_range;
} qaws_surface_swept_impl;

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

/* Swept surface: S(u,v) = P(u) + scale * [profile_x(v) * N(u) + profile_y(v) * B(u)]
   u follows the path, v follows the profile cross-section.
   N, B are Frenet normal and binormal of the path curve.

   dS/dv is computed analytically.
   dS/du uses central finite differences (frame derivatives are complex). */
static qaws_status swept_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_swept_impl const* impl =
		(qaws_surface_swept_impl const*)surface->impl;
	qaws_scalar s_path = impl->path_range.max_value - impl->path_range.min_value;
	qaws_scalar s_prof = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_scalar t_path = impl->path_range.min_value + u * s_path;
	qaws_scalar t_prof = impl->prof_range.min_value + v * s_prof;
	qaws_scalar scale = impl->scale;

	qaws_eval_result_3d path_r;
	qaws_eval_result_2d prof_r;
	qaws_vec3 T, N, B, pos;
	qaws_scalar px, py;

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

	/* Evaluate profile (2D) - need position and optionally D1 */
	{
		unsigned int pf = QAWS_EVAL_FLAG_POSITION;
		qaws_status s;
		if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DVV))
			pf |= QAWS_EVAL_FLAG_D1;
		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
			pf |= QAWS_EVAL_FLAG_D2;

		memset(&prof_r, 0, sizeof(prof_r));
		s = qaws_curve_evaluate_2d(impl->profile, t_prof, pf, &prof_r);
		if (s != QAWS_STATUS_OK) return s;
	}

	px = prof_r.position.x * scale;
	py = prof_r.position.y * scale;

	/* Position: P(u) + px * N(u) + py * B(u) */
	pos.x = path_r.position.x + px * N.x + py * B.x;
	pos.y = path_r.position.y + px * N.y + py * B.y;
	pos.z = path_r.position.z + px * N.z + py * B.z;

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/dv = s_prof * scale * [profile_x'(v) * N(u) + profile_y'(v) * B(u)] */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		qaws_scalar dpx = prof_r.d1.x * scale * s_prof;
		qaws_scalar dpy = prof_r.d1.y * scale * s_prof;
		out_result->dv.x = dpx * N.x + dpy * B.x;
		out_result->dv.y = dpx * N.y + dpy * B.y;
		out_result->dv.z = dpx * N.z + dpy * B.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* dvv = s_prof^2 * scale * [profile_x''(v) * N(u) + profile_y''(v) * B(u)] */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		qaws_scalar sp2 = s_prof * s_prof;
		qaws_scalar ddpx = prof_r.d2.x * scale * sp2;
		qaws_scalar ddpy = prof_r.d2.y * scale * sp2;
		out_result->dvv.x = ddpx * N.x + ddpy * B.x;
		out_result->dvv.y = ddpx * N.y + ddpy * B.y;
		out_result->dvv.z = ddpx * N.z + ddpy * B.z;
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
			p_lo.x = pr.position.x + px * N_lo.x + py * B_lo.x;
			p_lo.y = pr.position.y + px * N_lo.y + py * B_lo.y;
			p_lo.z = pr.position.z + px * N_lo.z + py * B_lo.z;

			memset(&pr, 0, sizeof(pr));
			qaws_curve_evaluate_3d(impl->path, t_hi, QAWS_EVAL_FLAG_POSITION, &pr);
			qaws_curve_compute_frenet_frame_3d(impl->path, t_hi, &T_hi, &N_hi, &B_hi);
			p_hi.x = pr.position.x + px * N_hi.x + py * B_hi.x;
			p_hi.y = pr.position.y + px * N_hi.y + py * B_hi.y;
			p_hi.z = pr.position.z + px * N_hi.z + py * B_hi.z;
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
			qaws_scalar dpx = prof_r.d1.x * scale * s_prof;
			qaws_scalar dpy = prof_r.d1.y * scale * s_prof;
			qaws_vec3 dv_lo, dv_hi;
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * h);

			qaws_curve_compute_frenet_frame_3d(impl->path, t_lo, &T_tmp, &N_lo, &B_lo);
			qaws_curve_compute_frenet_frame_3d(impl->path, t_hi, &T_tmp, &N_hi, &B_hi);

			dv_lo.x = dpx * N_lo.x + dpy * B_lo.x;
			dv_lo.y = dpx * N_lo.y + dpy * B_lo.y;
			dv_lo.z = dpx * N_lo.z + dpy * B_lo.z;

			dv_hi.x = dpx * N_hi.x + dpy * B_hi.x;
			dv_hi.y = dpx * N_hi.y + dpy * B_hi.y;
			dv_hi.z = dpx * N_hi.z + dpy * B_hi.z;

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

static void swept_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int swept_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const swept_surface_vtable = {
	swept_surface_eval,
	swept_surface_destroy,
	swept_surface_is_rational
};

qaws_status qaws_surface_create_swept(
	qaws_surface_swept_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_swept_impl* impl;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->path || !desc->profile) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(desc->path) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->profile) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_SWEPT,
		0, 0, u_range, v_range,
		&swept_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_swept_impl*)malloc(sizeof(qaws_surface_swept_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->path = desc->path;
	impl->profile = desc->profile;
	impl->scale = (desc->scale > 0) ? desc->scale : QAWS_ONE;
	impl->path_range = qaws_curve_get_parameter_range(desc->path);
	impl->prof_range = qaws_curve_get_parameter_range(desc->profile);

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
