#include "qaws_surface_extrusion.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_extrusion_impl
{
	qaws_curve const* profile;
	qaws_range prof_range;
	qaws_vec3 dir;          /* normalized direction */
	qaws_scalar length;
	int profile_is_2d;
	qaws_vec3 extent;       /* dir * length: the differentiable extrusion vector */
} qaws_surface_extrusion_impl;

/* S(u,v) = profile_3d(u) + v * dir * length
   dS/du = profile'(u) * s_prof  (lifted to 3D if 2D)
   dS/dv = dir * length
   dS/duu = profile''(u) * s_prof^2
   dS/dvv = 0
   dS/duv = 0 */
static qaws_status extrusion_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_extrusion_impl const* impl =
		(qaws_surface_extrusion_impl const*)surface->impl;
	qaws_scalar s_prof = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_scalar t_prof = impl->prof_range.min_value + u * s_prof;
	qaws_vec3 prof_pos;
	qaws_vec3 prof_d1;
	qaws_vec3 prof_d2;

	prof_pos.x = 0; prof_pos.y = 0; prof_pos.z = 0;
	prof_d1.x = 0; prof_d1.y = 0; prof_d1.z = 0;
	prof_d2.x = 0; prof_d2.y = 0; prof_d2.z = 0;

	if (impl->profile_is_2d)
	{
		unsigned int pf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_2d pr;
		qaws_status s;

		if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUV
			| QAWS_SURFACE_EVAL_NORMAL))
			pf |= QAWS_EVAL_FLAG_D1;
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
			pf |= QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;

		memset(&pr, 0, sizeof(pr));
		s = qaws_curve_evaluate_2d(impl->profile, t_prof, pf, &pr);
		if (s != QAWS_STATUS_OK) return s;

		/* Lift 2D to 3D: (x,y) -> (x, y, 0) */
		prof_pos.x = pr.position.x;
		prof_pos.y = pr.position.y;
		prof_d1.x = pr.d1.x;
		prof_d1.y = pr.d1.y;
		prof_d2.x = pr.d2.x;
		prof_d2.y = pr.d2.y;
	}
	else
	{
		unsigned int pf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d pr;
		qaws_status s;

		if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUV
			| QAWS_SURFACE_EVAL_NORMAL))
			pf |= QAWS_EVAL_FLAG_D1;
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
			pf |= QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;

		memset(&pr, 0, sizeof(pr));
		s = qaws_curve_evaluate_3d(impl->profile, t_prof, pf, &pr);
		if (s != QAWS_STATUS_OK) return s;

		prof_pos = pr.position;
		prof_d1 = pr.d1;
		prof_d2 = pr.d2;
	}

	/* Position: profile_3d(u) + v * dir * length */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = prof_pos.x + v * impl->dir.x * impl->length;
		out_result->position.y = prof_pos.y + v * impl->dir.y * impl->length;
		out_result->position.z = prof_pos.z + v * impl->dir.z * impl->length;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = profile'(u) * s_prof */
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = prof_d1.x * s_prof;
		out_result->du.y = prof_d1.y * s_prof;
		out_result->du.z = prof_d1.z * s_prof;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = dir * length */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = impl->dir.x * impl->length;
		out_result->dv.y = impl->dir.y * impl->length;
		out_result->dv.z = impl->dir.z * impl->length;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* dS/duu = profile''(u) * s_prof^2 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		qaws_scalar sp2 = s_prof * s_prof;
		out_result->duu.x = prof_d2.x * sp2;
		out_result->duu.y = prof_d2.y * sp2;
		out_result->duu.z = prof_d2.z * sp2;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* dS/dvv = 0 */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv.x = 0; out_result->dvv.y = 0; out_result->dvv.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* dS/duv = 0 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		out_result->duv.x = 0; out_result->duv.y = 0; out_result->duv.z = 0;
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

static void extrusion_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int extrusion_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

/* -------------------------------------------------------------------------- */
/*  Differential rules                                                        */
/*                                                                            */
/*  S(u,v) = P(t) + v E,  t = t0 + u s,  E = direction * length.              */
/*  Own field: direction (the extrusion vector E). Child 0: the profile.      */
/* -------------------------------------------------------------------------- */

#define EXTRUSION_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned char const g_extrusion_jet_a[QAWS_SURFACE_JET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_extrusion_jet_b[QAWS_SURFACE_JET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

static unsigned int extrusion_describe_fields(qaws_surface const* surface,
	qaws_field_desc* out, unsigned int capacity)
{
	(void)surface;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_DIRECTION, QAWS_VALUE_VEC3, 1,
			QAWS_DOMAIN_DIRECTION, QAWS_CONSTRAINT_NONE, QAWS_DIFF_SMOOTH, EXTRUSION_DIFF_CAPS);
	return 1;
}

static qaws_status extrusion_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_surface_extrusion_impl const* impl = (qaws_surface_extrusion_impl const*)surface->impl;
	if (field != QAWS_FIELD_DIRECTION)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_data = &impl->extent.x;
	*out_count = 1;
	*out_components = 3;
	return QAWS_STATUS_OK;
}

static unsigned int extrusion_children(qaws_surface const* surface, qaws_diff_child* out, unsigned int capacity)
{
	qaws_surface_extrusion_impl const* impl = (qaws_surface_extrusion_impl const*)surface->impl;
	if (capacity >= 1)
	{
		out[0].curve = impl->profile;
		out[0].surface = NULL;
	}
	return 1;
}

static qaws_scalar extrusion_power(qaws_scalar s, unsigned int a)
{
	qaws_scalar r = QAWS_ONE;
	while (a--)
		r *= s;
	return r;
}

static qaws_status extrusion_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_surface_jet* primal, qaws_surface_jet* tangent, qaws_surface_jet* tangent2)
{
	qaws_surface_extrusion_impl const* impl = (qaws_surface_extrusion_impl const*)surface->impl;
	qaws_scalar s = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_scalar t = impl->prof_range.min_value + u * s;
	qaws_field_view const* ev = views ? qaws_diff_views_find(views, QAWS_FIELD_DIRECTION) : NULL;
	qaws_curve_jet_3d cp, ct, ctt;
	qaws_vec3 e = impl->extent, e_dot = qaws_v3_zero();
	unsigned int ch;
	qaws_status st;

	if (ev)
	{
		qaws_scalar tmp[3];
		qaws_internal_view_read(ev, 0, 3, tmp);
		e_dot = qaws_v3(tmp[0], tmp[1], tmp[2]);
	}

	st = qaws_internal_curve_tangent_any(ctx, impl->profile, t, s * u_dot, 0xFu,
		qaws_internal_child_views(views, 0), &cp, &ct, tangent2 ? &ctt : NULL);
	if (st != QAWS_STATUS_OK)
		return st;

	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));

	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		unsigned int a = g_extrusion_jet_a[ch], b = g_extrusion_jet_b[ch];
		if (!(channels & (1u << ch)))
			continue;
		if (b == 0)
		{
			qaws_scalar sa = extrusion_power(s, a);
			primal->d[ch] = qaws_v3_scale(cp.d[a], sa);
			tangent->d[ch] = qaws_v3_scale(ct.d[a], sa);
			if (tangent2)
				tangent2->d[ch] = qaws_v3_scale(ctt.d[a], sa);
			if (a == 0)
			{
				primal->d[ch] = qaws_v3_axpy(primal->d[ch], e, v);
				tangent->d[ch] = qaws_v3_axpy(qaws_v3_axpy(tangent->d[ch], e, v_dot), e_dot, v);
				if (tangent2)
					tangent2->d[ch] = qaws_v3_axpy(tangent2->d[ch], e_dot, QAWS_LITERAL(2.0) * v_dot);
			}
		}
		else if (a == 0 && b == 1)
		{
			primal->d[ch] = e;
			tangent->d[ch] = e_dot;
		}
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

static qaws_status extrusion_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* ybar, qaws_diff_views* views,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	qaws_surface_extrusion_impl const* impl = (qaws_surface_extrusion_impl const*)surface->impl;
	qaws_scalar s = impl->prof_range.max_value - impl->prof_range.min_value;
	qaws_scalar t = impl->prof_range.min_value + u * s;
	qaws_field_view* ev = views ? qaws_diff_views_find(views, QAWS_FIELD_DIRECTION) : NULL;
	qaws_curve_jet_3d cbar;
	qaws_vec3 ebar = qaws_v3_zero();
	qaws_scalar tbar = QAWS_ZERO;
	unsigned int ch;
	qaws_status st;

	memset(&cbar, 0, sizeof(cbar));
	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		unsigned int a = g_extrusion_jet_a[ch], b = g_extrusion_jet_b[ch];
		if (!(channels & (1u << ch)))
			continue;
		if (b == 0)
		{
			cbar.d[a] = qaws_v3_add(cbar.d[a], qaws_v3_scale(ybar->d[ch], extrusion_power(s, a)));
			cbar.channels |= 1u << a;
			if (a == 0)
			{
				ebar = qaws_v3_axpy(ebar, ybar->d[ch], v);
				if (v_adjoint)
					*v_adjoint += qaws_v3_dot(ybar->d[ch], impl->extent);
			}
		}
		else if (a == 0 && b == 1)
			ebar = qaws_v3_add(ebar, ybar->d[ch]);
	}

	st = qaws_internal_curve_adjoint_any(ctx, impl->profile, t, cbar.channels, &cbar,
		(qaws_diff_views*)qaws_internal_child_views(views, 0), u_adjoint ? &tbar : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (u_adjoint)
		*u_adjoint += s * tbar;
	if (ev)
	{
		qaws_scalar g[3];
		g[0] = ebar.x; g[1] = ebar.y; g[2] = ebar.z;
		qaws_internal_view_add(ev, 0, 3, g);
	}
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const extrusion_surface_diff_vtable = {
	EXTRUSION_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	extrusion_describe_fields,
	extrusion_primal_field,
	NULL,
	NULL,
	extrusion_tangent,
	extrusion_adjoint,
	extrusion_children
};

static qaws_surface_vtable const extrusion_surface_vtable = {
	extrusion_surface_eval,
	extrusion_surface_destroy,
	extrusion_surface_is_rational,
	&extrusion_surface_diff_vtable
};

qaws_status qaws_surface_create_extrusion(
	qaws_surface_extrusion_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_extrusion_impl* impl;
	qaws_range u_range, v_range;
	qaws_scalar dir_len;
	qaws_dimension dim;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->profile) return QAWS_STATUS_INVALID_ARGUMENT;

	dim = qaws_curve_get_dimension(desc->profile);
	if (dim != QAWS_DIMENSION_2D && dim != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	dir_len = QAWS_SQRT(desc->direction.x * desc->direction.x
		+ desc->direction.y * desc->direction.y
		+ desc->direction.z * desc->direction.z);
	if (dir_len < QAWS_LITERAL(1e-12))
		return QAWS_STATUS_INVALID_ARGUMENT;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_EXTRUSION,
		0, 0, u_range, v_range,
		&extrusion_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_extrusion_impl*)malloc(sizeof(qaws_surface_extrusion_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->profile = desc->profile;
	impl->prof_range = qaws_curve_get_parameter_range(desc->profile);
	impl->profile_is_2d = (dim == QAWS_DIMENSION_2D) ? 1 : 0;

	/* Normalize direction */
	impl->dir.x = desc->direction.x / dir_len;
	impl->dir.y = desc->direction.y / dir_len;
	impl->dir.z = desc->direction.z / dir_len;

	/* Length: if desc->length > 0 use it, otherwise use the magnitude of direction */
	impl->length = (desc->length > QAWS_ZERO) ? desc->length : dir_len;
	impl->extent.x = impl->dir.x * impl->length;
	impl->extent.y = impl->dir.y * impl->length;
	impl->extent.z = impl->dir.z * impl->length;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
