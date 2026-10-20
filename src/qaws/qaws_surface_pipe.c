#include "qaws_surface_pipe.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include "core/qaws_frame_core.h"
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
	/* rotation-minimizing frames at uniform u samples, by double reflection (Wang et al. 2008): the path's point and unit
	   tangent there, the frame's normal and binormal */
	qaws_vec3 frame_P[PIPE_FRAME_SAMPLES + 1];
	qaws_vec3 frame_T[PIPE_FRAME_SAMPLES + 1];
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

/* the path at u: D1 and D2 in its parameter, its point and unit tangent (fallback where D1 vanishes) */
static void pipe_path(qaws_surface_pipe_impl const* impl, qaws_scalar u, qaws_vec3 fallback, qaws_eval_result_3d* c, qaws_vec3* T)
{
	qaws_scalar s_path = impl->path_range.max_value - impl->path_range.min_value;
	memset(c, 0, sizeof *c);
	qaws_curve_evaluate_3d(c, impl->path,
		impl->path_range.min_value + u * s_path, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
	*T = c->d1;
	if (v3_len(*T) > QAWS_LITERAL(1e-15))
		v3_norm(T);
	else
		*T = fallback;
}

/* N from the double reflection of the frame (P0, T0, N0) to (P1, T1), made unit and perpendicular to T1; B = T1 x N */
static void pipe_reflect(qaws_vec3 p0, qaws_vec3 t0, qaws_vec3 n0, qaws_vec3 p1, qaws_vec3 t1, qaws_vec3* N, qaws_vec3* B)
{
	qaws_vec3 r = qaws_double_reflection(p0, t0, n0, p1, t1);
	qaws_scalar d = v3_dot_l(r, t1);
	r.x -= d * t1.x;
	r.y -= d * t1.y;
	r.z -= d * t1.z;
	v3_norm(&r);
	*N = r;
	*B = v3_cross_l(t1, r);
}

static void precompute_frames(qaws_surface_pipe_impl* impl)
{
	unsigned int i;
	qaws_vec3 z = { 0, 0, 1 };
	for (i = 0; i <= PIPE_FRAME_SAMPLES; i++)
	{
		qaws_eval_result_3d c;
		qaws_vec3 T;
		pipe_path(impl, (qaws_scalar)i / (qaws_scalar)PIPE_FRAME_SAMPLES, i ? impl->frame_T[i - 1] : z, &c, &T);
		impl->frame_P[i] = c.position;
		impl->frame_T[i] = T;
		if (i == 0)
			build_frame(T, &impl->frame_N[0], &impl->frame_B[0]);
		else
			pipe_reflect(impl->frame_P[i - 1], impl->frame_T[i - 1], impl->frame_N[i - 1], c.position, T, &impl->frame_N[i], &impl->frame_B[i]);
	}
}

/* the frame at u: one double reflection from the sample at or below u, so it meets the samples and is smooth between */
static void pipe_frame(qaws_surface_pipe_impl const* impl, qaws_scalar u, qaws_eval_result_3d* c, qaws_vec3* T, qaws_vec3* N, qaws_vec3* B)
{
	qaws_scalar ui = u * (qaws_scalar)PIPE_FRAME_SAMPLES;
	unsigned int i = ui <= QAWS_ZERO ? 0 : ui >= (qaws_scalar)PIPE_FRAME_SAMPLES ? PIPE_FRAME_SAMPLES - 1 : (unsigned int)ui;
	pipe_path(impl, u, impl->frame_T[i], c, T);
	pipe_reflect(impl->frame_P[i], impl->frame_T[i], impl->frame_N[i], c->position, *T, N, B);
}

/* dS/du at (u, cos and sin of the angle), and the frame's rates in u: a rotation-minimizing frame turns with the tangent
   only, N' = -(T' . N) T, B' = -(T' . B) T, with T' = (C'' - (C'' . T) T) / |C'| */
static void pipe_du(qaws_surface_pipe_impl const* impl, qaws_scalar u, qaws_scalar ca, qaws_scalar sa, qaws_vec3* out_du,
	qaws_vec3* out_dn, qaws_vec3* out_db)
{
	qaws_scalar s_path = impl->path_range.max_value - impl->path_range.min_value;
	qaws_eval_result_3d c;
	qaws_vec3 T, N, B, dt = { 0, 0, 0 }, dn, db;
	qaws_scalar sp, kn, kb;
	pipe_frame(impl, u, &c, &T, &N, &B);
	sp = v3_len(c.d1);
	if (sp > QAWS_LITERAL(1e-15))
	{
		qaws_scalar a = v3_dot_l(c.d2, T);
		dt.x = (c.d2.x - a * T.x) / sp;
		dt.y = (c.d2.y - a * T.y) / sp;
		dt.z = (c.d2.z - a * T.z) / sp;
	}
	kn = v3_dot_l(dt, N);
	kb = v3_dot_l(dt, B);
	dn.x = -kn * T.x * s_path; dn.y = -kn * T.y * s_path; dn.z = -kn * T.z * s_path;
	db.x = -kb * T.x * s_path; db.y = -kb * T.y * s_path; db.z = -kb * T.z * s_path;
	out_du->x = s_path * c.d1.x + impl->radius_x * ca * dn.x + impl->radius_y * sa * db.x;
	out_du->y = s_path * c.d1.y + impl->radius_x * ca * dn.y + impl->radius_y * sa * db.y;
	out_du->z = s_path * c.d1.z + impl->radius_x * ca * dn.z + impl->radius_y * sa * db.z;
	if (out_dn)
	{
		*out_dn = dn;
		*out_db = db;
	}
}

/* ------------------------------------------------------------------ */
/*  Surface evaluation                                                 */
/* ------------------------------------------------------------------ */

/* Pipe surface: S(u,v) = P(u) + rx * cos(theta) * N(u) + ry * sin(theta) * B(u)
   where theta = 2 * pi * v.
   N(u), B(u) are a rotation-minimizing frame (pipe_frame).
   u follows the path, v parameterizes the circular/elliptical cross-section.

   dS/dv, dS/du and the mixed derivative are analytic (pipe_du); duu is the central difference of the analytic dS/du. */
static qaws_status pipe_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_pipe_impl const* impl =
		(qaws_surface_pipe_impl const*)surface->impl;
	qaws_scalar rx = impl->radius_x;
	qaws_scalar ry = impl->radius_y;
	qaws_scalar pi = QAWS_LITERAL(3.14159265358979323846);
	qaws_scalar two_pi = QAWS_LITERAL(2.0) * pi;
	qaws_scalar theta = two_pi * v;
	qaws_scalar cos_theta = QAWS_COS(theta);
	qaws_scalar sin_theta = QAWS_SIN(theta);

	qaws_eval_result_3d path_r;
	qaws_vec3 T, N, B, pos;

	pipe_frame(impl, u, &path_r, &T, &N, &B);

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

	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV))
	{
		qaws_vec3 du, dn, db;
		pipe_du(impl, u, cos_theta, sin_theta, &du, &dn, &db);

		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			out_result->du = du;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}

		/* duu: the central difference of the analytic dS/du (one-sided at the ends) */
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar h = QAWS_LITERAL(1e-4);
			qaws_scalar u_lo = u - h < 0 ? 0 : u - h, u_hi = u + h > 1 ? 1 : u + h;
			qaws_vec3 d_lo, d_hi;
			pipe_du(impl, u_lo, cos_theta, sin_theta, &d_lo, NULL, NULL);
			pipe_du(impl, u_hi, cos_theta, sin_theta, &d_hi, NULL, NULL);
			out_result->duu.x = (d_hi.x - d_lo.x) / (u_hi - u_lo);
			out_result->duu.y = (d_hi.y - d_lo.y) / (u_hi - u_lo);
			out_result->duu.z = (d_hi.z - d_lo.z) / (u_hi - u_lo);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		/* duv = 2*pi * (-rx * sin(theta) * N'(u) + ry * cos(theta) * B'(u)) */
		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			out_result->duv.x = two_pi * (-rx * sin_theta * dn.x + ry * cos_theta * db.x);
			out_result->duv.y = two_pi * (-rx * sin_theta * dn.y + ry * cos_theta * db.y);
			out_result->duv.z = two_pi * (-rx * sin_theta * dn.z + ry * cos_theta * db.z);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
		}
	}

	/* Normal */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(&out_result->normal, out_result->du, out_result->dv);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

static void pipe_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(impl, allocator);
}

static int pipe_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const pipe_surface_vtable = {
	pipe_surface_eval,
	pipe_surface_destroy,
	pipe_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_pipe_ex(
	qaws_surface_pipe_desc const* desc,
	qaws_allocator const* allocator,
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

	surface = qaws_internal_surface_alloc_ex(QAWS_SURFACE_KIND_PIPE,
		0, 0, u_range, v_range,
		&pipe_surface_vtable, allocator);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_pipe_impl*)qaws_internal_alloc(sizeof(qaws_surface_pipe_impl), allocator);
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->path = desc->path;
	impl->path_range = qaws_curve_get_parameter_range(desc->path);
	impl->radius_x = desc->radius_x;
	impl->radius_y = ry;

	/* rotation-minimizing frames along the path */
	precompute_frames(impl);

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_create_pipe(
	qaws_surface** out_surface,
	qaws_surface_pipe_desc const* desc)
{
	return qaws_surface_create_pipe_ex(desc, NULL, out_surface);
}
