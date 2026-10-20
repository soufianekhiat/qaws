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

/* dS/du at u for the profile point (px, py) (scaled), and the frame's rates N', B' in u: the Frenet-Serret formulas in the
   path parameter t, N' = |C'| (-k T + t B), B' = -|C'| t N (k, t: curvature and torsion, from D1 to D3); where D1 x D2
   vanishes the frame is the constant one qaws_curve_compute_frenet_frame_3d picks, so it does not turn */
static qaws_status swept_du(qaws_surface_swept_impl const* impl, qaws_scalar u, qaws_scalar px, qaws_scalar py, qaws_vec3* out_du,
	qaws_vec3* out_dn, qaws_vec3* out_db)
{
	qaws_scalar s_path = impl->path_range.max_value - impl->path_range.min_value;
	qaws_scalar t_path = impl->path_range.min_value + u * s_path;
	qaws_eval_result_3d c;
	qaws_vec3 T, N, B, dn = { 0, 0, 0 }, db = { 0, 0, 0 };
	double x[3], xl2, sp, kappa_s, tau;
	qaws_status st;
	memset(&c, 0, sizeof c);
	if ((st = qaws_curve_evaluate_3d(&c, impl->path, t_path, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) != QAWS_STATUS_OK)
		return st;
	if ((st = qaws_curve_compute_frenet_frame_3d(&T, &N, &B, impl->path, t_path)) != QAWS_STATUS_OK)
		return st;
	x[0] = (double)c.d1.y * c.d2.z - (double)c.d1.z * c.d2.y;
	x[1] = (double)c.d1.z * c.d2.x - (double)c.d1.x * c.d2.z;
	x[2] = (double)c.d1.x * c.d2.y - (double)c.d1.y * c.d2.x;
	xl2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	sp = sqrt((double)c.d1.x * c.d1.x + (double)c.d1.y * c.d1.y + (double)c.d1.z * c.d1.z);
	/* the threshold of qaws_curve_compute_frenet_frame_3d */
	if (sqrt(xl2) >= 1.0e-12 && sp > 0)
	{
		kappa_s = sqrt(xl2) / (sp * sp);      /* |C'| k */
		tau = (x[0] * c.d3.x + x[1] * c.d3.y + x[2] * c.d3.z) / xl2;
		dn.x = (qaws_scalar)(-kappa_s * T.x + sp * tau * B.x);
		dn.y = (qaws_scalar)(-kappa_s * T.y + sp * tau * B.y);
		dn.z = (qaws_scalar)(-kappa_s * T.z + sp * tau * B.z);
		db.x = (qaws_scalar)(-sp * tau * N.x);
		db.y = (qaws_scalar)(-sp * tau * N.y);
		db.z = (qaws_scalar)(-sp * tau * N.z);
	}
	out_du->x = s_path * (c.d1.x + px * dn.x + py * db.x);
	out_du->y = s_path * (c.d1.y + px * dn.y + py * db.y);
	out_du->z = s_path * (c.d1.z + px * dn.z + py * db.z);
	if (out_dn)
	{
		out_dn->x = s_path * dn.x; out_dn->y = s_path * dn.y; out_dn->z = s_path * dn.z;
		out_db->x = s_path * db.x; out_db->y = s_path * db.y; out_db->z = s_path * db.z;
	}
	return QAWS_STATUS_OK;
}

/* Swept surface: S(u,v) = P(u) + scale * [profile_x(v) * N(u) + profile_y(v) * B(u)]
   u follows the path, v follows the profile cross-section.
   N, B are Frenet normal and binormal of the path curve.

   dS/dv, dS/du and the mixed derivative are analytic (swept_du); duu is the central difference of the analytic dS/du
   (the torsion's rate would need the path's fourth derivative). */
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
		qaws_status s = qaws_curve_evaluate_3d(&path_r, impl->path,
			t_path, QAWS_EVAL_FLAG_POSITION);
		if (s != QAWS_STATUS_OK) return s;
	}

	/* Get Frenet frame at path point */
	{
		qaws_status s = qaws_curve_compute_frenet_frame_3d(
			&T, &N, &B, impl->path, t_path);
		if (s != QAWS_STATUS_OK) return s;
	}

	/* Evaluate profile (2D) - need position and optionally D1 */
	{
		unsigned int pf = QAWS_EVAL_FLAG_POSITION;
		qaws_status s;
		if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DVV | QAWS_SURFACE_EVAL_DUV))
			pf |= QAWS_EVAL_FLAG_D1;
		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
			pf |= QAWS_EVAL_FLAG_D2;

		memset(&prof_r, 0, sizeof(prof_r));
		s = qaws_curve_evaluate_2d(&prof_r, impl->profile, t_prof, pf);
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

	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_vec3 du, dn, db;
		qaws_status s = swept_du(impl, u, px, py, &du, &dn, &db);
		if (s != QAWS_STATUS_OK) return s;

		if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL))
		{
			out_result->du = du;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}

		/* duu: the central difference of the analytic dS/du (one-sided at the ends) */
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar h = QAWS_LITERAL(1e-4);
			qaws_scalar u_lo = u - h < 0 ? 0 : u - h, u_hi = u + h > 1 ? 1 : u + h;
			qaws_vec3 d_lo = { 0, 0, 0 }, d_hi = { 0, 0, 0 };
			s = swept_du(impl, u_lo, px, py, &d_lo, NULL, NULL);
			if (s == QAWS_STATUS_OK)
				s = swept_du(impl, u_hi, px, py, &d_hi, NULL, NULL);
			if (s != QAWS_STATUS_OK) return s;
			out_result->duu.x = (d_hi.x - d_lo.x) / (u_hi - u_lo);
			out_result->duu.y = (d_hi.y - d_lo.y) / (u_hi - u_lo);
			out_result->duu.z = (d_hi.z - d_lo.z) / (u_hi - u_lo);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		/* duv = s_prof * scale * [profile_x'(v) * N'(u) + profile_y'(v) * B'(u)] */
		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			qaws_scalar dpx = prof_r.d1.x * scale * s_prof;
			qaws_scalar dpy = prof_r.d1.y * scale * s_prof;
			out_result->duv.x = dpx * dn.x + dpy * db.x;
			out_result->duv.y = dpx * dn.y + dpy * db.y;
			out_result->duv.z = dpx * dn.z + dpy * db.z;
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

static void swept_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(impl, allocator);
}

static int swept_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const swept_surface_vtable = {
	swept_surface_eval,
	swept_surface_destroy,
	swept_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_swept_ex(
	qaws_surface_swept_desc const* desc,
	qaws_allocator const* allocator,
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

	surface = qaws_internal_surface_alloc_ex(QAWS_SURFACE_KIND_SWEPT,
		0, 0, u_range, v_range,
		&swept_surface_vtable, allocator);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_swept_impl*)qaws_internal_alloc(sizeof(qaws_surface_swept_impl), allocator);
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

qaws_status qaws_surface_create_swept(
	qaws_surface** out_surface,
	qaws_surface_swept_desc const* desc)
{
	return qaws_surface_create_swept_ex(desc, NULL, out_surface);
}
