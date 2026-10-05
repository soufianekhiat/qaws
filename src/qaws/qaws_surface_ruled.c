#include "qaws_surface_ruled.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_ruled_impl
{
	qaws_curve const* curve_a;
	qaws_curve const* curve_b;
	qaws_range range_a;
	qaws_range range_b;
} qaws_surface_ruled_impl;

/* S(u,v) = (1-v) * A(u_a) + v * B(u_b)
   dS/du  = (1-v) * A'(u_a)*s_a + v * B'(u_b)*s_b
   dS/dv  = B(u_b) - A(u_a)
   dS/duu = (1-v) * A''(u_a)*s_a^2 + v * B''(u_b)*s_b^2
   dS/dvv = 0
   dS/duv = B'(u_b)*s_b - A'(u_a)*s_a
   where s_a = range_a length, s_b = range_b length */
static qaws_status ruled_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_ruled_impl const* impl =
		(qaws_surface_ruled_impl const*)surface->impl;
	qaws_scalar one_minus_v = QAWS_ONE - v;
	qaws_scalar s_a = impl->range_a.max_value - impl->range_a.min_value;
	qaws_scalar s_b = impl->range_b.max_value - impl->range_b.min_value;
	qaws_scalar u_a = impl->range_a.min_value + u * s_a;
	qaws_scalar u_b = impl->range_b.min_value + u * s_b;

	unsigned int curve_flags = 0;
	qaws_eval_result_3d ra, rb;
	qaws_status sa_status, sb_status;

	/* Determine what we need from the curves */
	if (eval_flags & (QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DV))
		curve_flags |= QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUV))
		curve_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		curve_flags |= QAWS_EVAL_FLAG_D2;

	memset(&ra, 0, sizeof(ra));
	memset(&rb, 0, sizeof(rb));

	sa_status = qaws_curve_evaluate_3d(impl->curve_a, u_a, curve_flags, &ra);
	sb_status = qaws_curve_evaluate_3d(impl->curve_b, u_b, curve_flags, &rb);
	if (sa_status != QAWS_STATUS_OK) return sa_status;
	if (sb_status != QAWS_STATUS_OK) return sb_status;

	/* Position: (1-v)*A + v*B */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = one_minus_v * ra.position.x + v * rb.position.x;
		out_result->position.y = one_minus_v * ra.position.y + v * rb.position.y;
		out_result->position.z = one_minus_v * ra.position.z + v * rb.position.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* du: (1-v)*A'*s_a + v*B'*s_b */
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = one_minus_v * ra.d1.x * s_a + v * rb.d1.x * s_b;
		out_result->du.y = one_minus_v * ra.d1.y * s_a + v * rb.d1.y * s_b;
		out_result->du.z = one_minus_v * ra.d1.z * s_a + v * rb.d1.z * s_b;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dv: B - A */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = rb.position.x - ra.position.x;
		out_result->dv.y = rb.position.y - ra.position.y;
		out_result->dv.z = rb.position.z - ra.position.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* duu: (1-v)*A''*s_a^2 + v*B''*s_b^2 */
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		qaws_scalar sa2 = s_a * s_a, sb2 = s_b * s_b;
		out_result->duu.x = one_minus_v * ra.d2.x * sa2 + v * rb.d2.x * sb2;
		out_result->duu.y = one_minus_v * ra.d2.y * sa2 + v * rb.d2.y * sb2;
		out_result->duu.z = one_minus_v * ra.d2.z * sa2 + v * rb.d2.z * sb2;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* dvv: 0 (linear in v) */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv.x = 0; out_result->dvv.y = 0; out_result->dvv.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* duv: B'*s_b - A'*s_a */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		out_result->duv.x = rb.d1.x * s_b - ra.d1.x * s_a;
		out_result->duv.y = rb.d1.y * s_b - ra.d1.y * s_a;
		out_result->duv.z = rb.d1.z * s_b - ra.d1.z * s_a;
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

static void ruled_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	/* We don't own the curves - just free the impl struct */
	qaws_internal_dealloc(allocator, impl);
}

static int ruled_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

/* -------------------------------------------------------------------------- */
/*  Differential rules                                                        */
/*                                                                            */
/*  S(u,v) = (1 - v) A(ta) + v B(tb),  ta = a0 + u sa,  tb = b0 + u sb.       */
/*  No own fields. Child 0: curve A, child 1: curve B.                        */
/* -------------------------------------------------------------------------- */

#define RULED_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned char const g_ruled_jet_a[QAWS_SURFACE_JET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_ruled_jet_b[QAWS_SURFACE_JET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

static unsigned int ruled_describe_fields(qaws_surface const* surface,
	qaws_field_desc* out, unsigned int capacity)
{
	(void)surface;
	(void)out;
	(void)capacity;
	return 0;
}

static qaws_status ruled_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	(void)surface;
	(void)field;
	(void)out_data;
	(void)out_count;
	(void)out_components;
	return QAWS_STATUS_INVALID_ARGUMENT;
}

static unsigned int ruled_children(qaws_surface const* surface, qaws_diff_child* out, unsigned int capacity)
{
	qaws_surface_ruled_impl const* impl = (qaws_surface_ruled_impl const*)surface->impl;
	if (capacity >= 2)
	{
		out[0].curve = impl->curve_a;
		out[0].surface = NULL;
		out[1].curve = impl->curve_b;
		out[1].surface = NULL;
	}
	return 2;
}

static qaws_scalar ruled_power(qaws_scalar s, unsigned int a)
{
	qaws_scalar r = QAWS_ONE;
	while (a--)
		r *= s;
	return r;
}

static qaws_status ruled_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_surface_jet* primal, qaws_surface_jet* tangent, qaws_surface_jet* tangent2)
{
	qaws_surface_ruled_impl const* impl = (qaws_surface_ruled_impl const*)surface->impl;
	qaws_scalar sa = impl->range_a.max_value - impl->range_a.min_value;
	qaws_scalar sb = impl->range_b.max_value - impl->range_b.min_value;
	qaws_curve_jet_3d ap, at, att, bp, bt, btt;
	unsigned int ch;
	qaws_status st;

	st = qaws_internal_curve_tangent_any(ctx, impl->curve_a, impl->range_a.min_value + u * sa, sa * u_dot,
		0xFu, qaws_internal_child_views(views, 0), &ap, &at, tangent2 ? &att : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	st = qaws_internal_curve_tangent_any(ctx, impl->curve_b, impl->range_b.min_value + u * sb, sb * u_dot,
		0xFu, qaws_internal_child_views(views, 1), &bp, &bt, tangent2 ? &btt : NULL);
	if (st != QAWS_STATUS_OK)
		return st;

	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));

	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		unsigned int a = g_ruled_jet_a[ch], b = g_ruled_jet_b[ch];
		qaws_scalar pa = ruled_power(sa, a), pb = ruled_power(sb, a);
		if (!(channels & (1u << ch)) || b > 1)
			continue;
		if (b == 0)
		{
			/* (1-v) A_a + v B_a; v' couples with the v-derivative B_a - A_a. */
			qaws_vec3 diff_p = qaws_v3_sub(qaws_v3_scale(bp.d[a], pb), qaws_v3_scale(ap.d[a], pa));
			qaws_vec3 diff_t = qaws_v3_sub(qaws_v3_scale(bt.d[a], pb), qaws_v3_scale(at.d[a], pa));
			primal->d[ch] = qaws_v3_add(qaws_v3_scale(ap.d[a], (QAWS_ONE - v) * pa), qaws_v3_scale(bp.d[a], v * pb));
			tangent->d[ch] = qaws_v3_axpy(qaws_v3_add(qaws_v3_scale(at.d[a], (QAWS_ONE - v) * pa),
				qaws_v3_scale(bt.d[a], v * pb)), diff_p, v_dot);
			if (tangent2)
				tangent2->d[ch] = qaws_v3_axpy(qaws_v3_add(qaws_v3_scale(att.d[a], (QAWS_ONE - v) * pa),
					qaws_v3_scale(btt.d[a], v * pb)), diff_t, QAWS_LITERAL(2.0) * v_dot);
		}
		else
		{
			primal->d[ch] = qaws_v3_sub(qaws_v3_scale(bp.d[a], pb), qaws_v3_scale(ap.d[a], pa));
			tangent->d[ch] = qaws_v3_sub(qaws_v3_scale(bt.d[a], pb), qaws_v3_scale(at.d[a], pa));
			if (tangent2)
				tangent2->d[ch] = qaws_v3_sub(qaws_v3_scale(btt.d[a], pb), qaws_v3_scale(att.d[a], pa));
		}
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

static qaws_status ruled_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* ybar, qaws_diff_views* views,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	qaws_surface_ruled_impl const* impl = (qaws_surface_ruled_impl const*)surface->impl;
	qaws_scalar sa = impl->range_a.max_value - impl->range_a.min_value;
	qaws_scalar sb = impl->range_b.max_value - impl->range_b.min_value;
	qaws_scalar ta = impl->range_a.min_value + u * sa;
	qaws_scalar tb = impl->range_b.min_value + u * sb;
	qaws_curve_jet_3d abar, bbar;
	qaws_scalar tabar = QAWS_ZERO, tbbar = QAWS_ZERO;
	unsigned int ch;
	qaws_status st;

	memset(&abar, 0, sizeof(abar));
	memset(&bbar, 0, sizeof(bbar));

	if (v_adjoint)
	{
		/* dS_{a,0}/dv = S_{a,1} = B_a sb^a - A_a sa^a */
		qaws_curve_jet_3d ap, at, bp, bt;
		st = qaws_internal_curve_tangent_any(ctx, impl->curve_a, ta, QAWS_ZERO, 0xFu, NULL, &ap, &at, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		st = qaws_internal_curve_tangent_any(ctx, impl->curve_b, tb, QAWS_ZERO, 0xFu, NULL, &bp, &bt, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			unsigned int a = g_ruled_jet_a[ch];
			if (!(channels & (1u << ch)) || g_ruled_jet_b[ch] != 0)
				continue;
			*v_adjoint += qaws_v3_dot(ybar->d[ch],
				qaws_v3_sub(qaws_v3_scale(bp.d[a], ruled_power(sb, a)), qaws_v3_scale(ap.d[a], ruled_power(sa, a))));
		}
	}

	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		unsigned int a = g_ruled_jet_a[ch], b = g_ruled_jet_b[ch];
		qaws_scalar pa = ruled_power(sa, a), pb = ruled_power(sb, a);
		if (!(channels & (1u << ch)) || b > 1)
			continue;
		if (b == 0)
		{
			abar.d[a] = qaws_v3_axpy(abar.d[a], ybar->d[ch], (QAWS_ONE - v) * pa);
			bbar.d[a] = qaws_v3_axpy(bbar.d[a], ybar->d[ch], v * pb);
		}
		else
		{
			abar.d[a] = qaws_v3_axpy(abar.d[a], ybar->d[ch], -pa);
			bbar.d[a] = qaws_v3_axpy(bbar.d[a], ybar->d[ch], pb);
		}
		abar.channels |= 1u << a;
		bbar.channels |= 1u << a;
	}

	st = qaws_internal_curve_adjoint_any(ctx, impl->curve_a, ta, abar.channels, &abar,
		(qaws_diff_views*)qaws_internal_child_views(views, 0), u_adjoint ? &tabar : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	st = qaws_internal_curve_adjoint_any(ctx, impl->curve_b, tb, bbar.channels, &bbar,
		(qaws_diff_views*)qaws_internal_child_views(views, 1), u_adjoint ? &tbbar : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (u_adjoint)
		*u_adjoint += sa * tabar + sb * tbbar;
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const ruled_surface_diff_vtable = {
	RULED_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	ruled_describe_fields,
	ruled_primal_field,
	NULL,
	NULL,
	ruled_tangent,
	ruled_adjoint,
	ruled_children
};

static qaws_surface_vtable const ruled_surface_vtable = {
	ruled_surface_eval,
	ruled_surface_destroy,
	ruled_surface_is_rational,
	&ruled_surface_diff_vtable
};

qaws_status qaws_surface_create_ruled(
	qaws_surface_ruled_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_ruled_impl* impl;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->curve_a || !desc->curve_b) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(desc->curve_a) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->curve_b) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_RULED,
		0, 0, u_range, v_range,
		&ruled_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_ruled_impl*)malloc(sizeof(qaws_surface_ruled_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->curve_a = desc->curve_a;
	impl->curve_b = desc->curve_b;
	impl->range_a = qaws_curve_get_parameter_range(desc->curve_a);
	impl->range_b = qaws_curve_get_parameter_range(desc->curve_b);

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
