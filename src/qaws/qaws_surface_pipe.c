#include "qaws_surface_pipe.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

#define PIPE_FRAME_SAMPLES 256

typedef struct qaws_surface_pipe_impl
{
	qaws_curve const* path;
	qaws_range path_range;
	qaws_scalar radius_x;
	qaws_scalar radius_y;
	/* Precomputed parallel transport (Bishop) frames at uniform u samples */
	qaws_vec3 frame_N[PIPE_FRAME_SAMPLES + 1];
	qaws_vec3 frame_B[PIPE_FRAME_SAMPLES + 1];
} qaws_surface_pipe_impl;

/* ------------------------------------------------------------------ */
/*  Vector helpers (local to this file)                                */
/* ------------------------------------------------------------------ */

static qaws_scalar v3_len(qaws_vec3 a)
{
	return QAWS_SQRT(a.x * a.x + a.y * a.y + a.z * a.z);
}

static void v3_norm(qaws_vec3* a)
{
	qaws_scalar len = v3_len(*a);
	if (len > QAWS_LITERAL(1e-15))
	{
		a->x /= len; a->y /= len; a->z /= len;
	}
}

static qaws_scalar v3_dot_l(qaws_vec3 a, qaws_vec3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static qaws_vec3 v3_cross_l(qaws_vec3 a, qaws_vec3 b)
{
	qaws_vec3 r;
	r.x = a.y * b.z - a.z * b.y;
	r.y = a.z * b.x - a.x * b.z;
	r.z = a.x * b.y - a.y * b.x;
	return r;
}

/* Rotate vector v around unit axis by angle (cos_a, sin_a) via Rodrigues */
static qaws_vec3 v3_rotate_l(qaws_vec3 v, qaws_vec3 axis,
	qaws_scalar cos_a, qaws_scalar sin_a)
{
	qaws_scalar d = v3_dot_l(axis, v);
	qaws_vec3 cr = v3_cross_l(axis, v);
	qaws_vec3 r;
	r.x = v.x * cos_a + cr.x * sin_a + axis.x * d * (QAWS_ONE - cos_a);
	r.y = v.y * cos_a + cr.y * sin_a + axis.y * d * (QAWS_ONE - cos_a);
	r.z = v.z * cos_a + cr.z * sin_a + axis.z * d * (QAWS_ONE - cos_a);
	return r;
}

/* Build an initial perpendicular frame for a given tangent */
static void build_frame(qaws_vec3 T, qaws_vec3* N, qaws_vec3* B)
{
	qaws_vec3 up;
	if (QAWS_FABS(T.x) < QAWS_LITERAL(0.9))
	{
		up.x = QAWS_ONE; up.y = QAWS_ZERO; up.z = QAWS_ZERO;
	}
	else
	{
		up.x = QAWS_ZERO; up.y = QAWS_ONE; up.z = QAWS_ZERO;
	}
	*B = v3_cross_l(T, up);
	v3_norm(B);
	*N = v3_cross_l(*B, T);
	v3_norm(N);
}

/* ------------------------------------------------------------------ */
/*  Precompute parallel transport frames along the path                */
/* ------------------------------------------------------------------ */

static void precompute_bishop_frames(qaws_surface_pipe_impl* impl)
{
	unsigned int i;
	qaws_scalar s_path = impl->path_range.max_value - impl->path_range.min_value;
	qaws_vec3 prev_T = {0,0,0}, prev_N = {0,0,0}, prev_B = {0,0,0};
	int initialized = 0;

	for (i = 0; i <= PIPE_FRAME_SAMPLES; i++)
	{
		qaws_scalar u = (qaws_scalar)i / (qaws_scalar)PIPE_FRAME_SAMPLES;
		qaws_scalar t_path = impl->path_range.min_value + u * s_path;
		qaws_vec3 T;

		if (qaws_curve_compute_tangent_3d(impl->path, t_path, &T) != QAWS_STATUS_OK)
		{
			/* Fallback: use arbitrary frame */
			T.x = QAWS_ZERO; T.y = QAWS_ZERO; T.z = QAWS_ONE;
		}
		v3_norm(&T);

		if (!initialized)
		{
			build_frame(T, &impl->frame_N[i], &impl->frame_B[i]);
			prev_T = T;
			prev_N = impl->frame_N[i];
			prev_B = impl->frame_B[i];
			initialized = 1;
		}
		else
		{
			/* Parallel transport: rotate previous N,B to align with new T */
			qaws_scalar d = v3_dot_l(prev_T, T);
			if (d > QAWS_ONE) d = QAWS_ONE;
			if (d < -QAWS_ONE) d = -QAWS_ONE;

			if (d > QAWS_LITERAL(0.9999999))
			{
				/* Nearly parallel: keep previous frame */
				impl->frame_N[i] = prev_N;
				impl->frame_B[i] = prev_B;
			}
			else
			{
				qaws_vec3 rot_axis = v3_cross_l(prev_T, T);
				qaws_scalar sin_a, cos_a;
				v3_norm(&rot_axis);
				cos_a = d;
				sin_a = QAWS_SQRT(QAWS_ONE - d * d);
				impl->frame_N[i] = v3_rotate_l(prev_N, rot_axis, cos_a, sin_a);
				impl->frame_B[i] = v3_rotate_l(prev_B, rot_axis, cos_a, sin_a);
				v3_norm(&impl->frame_N[i]);
				v3_norm(&impl->frame_B[i]);
			}
			prev_T = T;
			prev_N = impl->frame_N[i];
			prev_B = impl->frame_B[i];
		}
	}
}

/* Interpolate the precomputed Bishop frame at a given u in [0,1] */
static void lookup_bishop_frame(
	qaws_surface_pipe_impl const* impl,
	qaws_scalar u,
	qaws_vec3* out_N,
	qaws_vec3* out_B)
{
	qaws_scalar u_idx = u * (qaws_scalar)PIPE_FRAME_SAMPLES;
	unsigned int lo, hi;
	qaws_scalar frac;

	if (u_idx < QAWS_ZERO) u_idx = QAWS_ZERO;
	if (u_idx > (qaws_scalar)PIPE_FRAME_SAMPLES) u_idx = (qaws_scalar)PIPE_FRAME_SAMPLES;

	lo = (unsigned int)u_idx;
	if (lo >= PIPE_FRAME_SAMPLES) lo = PIPE_FRAME_SAMPLES;
	hi = lo + 1;
	if (hi > PIPE_FRAME_SAMPLES) hi = PIPE_FRAME_SAMPLES;
	frac = u_idx - (qaws_scalar)lo;

	out_N->x = (QAWS_ONE - frac) * impl->frame_N[lo].x + frac * impl->frame_N[hi].x;
	out_N->y = (QAWS_ONE - frac) * impl->frame_N[lo].y + frac * impl->frame_N[hi].y;
	out_N->z = (QAWS_ONE - frac) * impl->frame_N[lo].z + frac * impl->frame_N[hi].z;
	out_B->x = (QAWS_ONE - frac) * impl->frame_B[lo].x + frac * impl->frame_B[hi].x;
	out_B->y = (QAWS_ONE - frac) * impl->frame_B[lo].y + frac * impl->frame_B[hi].y;
	out_B->z = (QAWS_ONE - frac) * impl->frame_B[lo].z + frac * impl->frame_B[hi].z;

	v3_norm(out_N);
	v3_norm(out_B);
}

/* ------------------------------------------------------------------ */
/*  Surface evaluation                                                 */
/* ------------------------------------------------------------------ */

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
   N(u), B(u) are from the precomputed Bishop (parallel transport) frame.
   u follows the path, v parameterizes the circular/elliptical cross-section.

   dS/dv is computed analytically.
   dS/du uses central finite differences. */
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
	qaws_vec3 N, B, pos;

	/* Evaluate path position */
	memset(&path_r, 0, sizeof(path_r));
	{
		qaws_status s = qaws_curve_evaluate_3d(impl->path, t_path,
			QAWS_EVAL_FLAG_POSITION, &path_r);
		if (s != QAWS_STATUS_OK) return s;
	}

	/* Get Bishop frame at this u */
	lookup_bishop_frame(impl, u, &N, &B);

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

		/* Evaluate position at u-h and u+h using Bishop frame */
		{
			qaws_scalar t_lo = impl->path_range.min_value + u_lo * s_path;
			qaws_scalar t_hi = impl->path_range.min_value + u_hi * s_path;
			qaws_eval_result_3d pr;
			qaws_vec3 N_lo, B_lo, N_hi, B_hi;

			memset(&pr, 0, sizeof(pr));
			qaws_curve_evaluate_3d(impl->path, t_lo, QAWS_EVAL_FLAG_POSITION, &pr);
			lookup_bishop_frame(impl, u_lo, &N_lo, &B_lo);
			p_lo.x = pr.position.x + rx * cos_theta * N_lo.x + ry * sin_theta * B_lo.x;
			p_lo.y = pr.position.y + rx * cos_theta * N_lo.y + ry * sin_theta * B_lo.y;
			p_lo.z = pr.position.z + rx * cos_theta * N_lo.z + ry * sin_theta * B_lo.z;

			memset(&pr, 0, sizeof(pr));
			qaws_curve_evaluate_3d(impl->path, t_hi, QAWS_EVAL_FLAG_POSITION, &pr);
			lookup_bishop_frame(impl, u_hi, &N_hi, &B_hi);
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
			qaws_vec3 N_lo, B_lo, N_hi, B_hi;
			qaws_scalar dnx = -rx * sin_theta;
			qaws_scalar dny = ry * cos_theta;
			qaws_vec3 dv_lo, dv_hi;
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * h);

			lookup_bishop_frame(impl, u_lo, &N_lo, &B_lo);
			lookup_bishop_frame(impl, u_hi, &N_hi, &B_hi);

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

	/* Precompute Bishop frames along the path */
	precompute_bishop_frames(impl);

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
