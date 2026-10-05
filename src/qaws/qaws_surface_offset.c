#include "qaws_surface_offset.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "internal/qaws_internal_bijet.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_offset_impl
{
	qaws_surface const* base;
	qaws_scalar distance;
} qaws_surface_offset_impl;

/* Helper: evaluate offset position at (u,v).
   S_off(u,v) = S_base(u,v) + d * N_base(u,v) */
static qaws_status offset_eval_position(
	qaws_surface const* base,
	qaws_scalar distance,
	qaws_scalar u,
	qaws_scalar v,
	qaws_vec3* out_pos)
{
	qaws_surface_eval_result base_r;
	qaws_status status;
	unsigned int flags = QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU
		| QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL;

	memset(&base_r, 0, sizeof(base_r));
	status = qaws_surface_evaluate(base, u, v, flags, &base_r);
	if (status != QAWS_STATUS_OK) return status;

	out_pos->x = base_r.position.x + distance * base_r.normal.x;
	out_pos->y = base_r.position.y + distance * base_r.normal.y;
	out_pos->z = base_r.position.z + distance * base_r.normal.z;
	return QAWS_STATUS_OK;
}

/* Sampled offset evaluation, kept for bases without analytic jets.
   Position: S_off(u,v) = S_base(u,v) + d * N_base(u,v)
   Derivatives: central differences of the offset position. */
static qaws_status offset_surface_eval_sampled(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_offset_impl const* impl =
		(qaws_surface_offset_impl const*)surface->impl;
	qaws_vec3 pos;
	qaws_status status;

	/* Always compute position (needed by derivatives too) */
	status = offset_eval_position(impl->base, impl->distance, u, v, &pos);
	if (status != QAWS_STATUS_OK) return status;

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* First derivatives via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_scalar u_lo = u - h;
		qaws_scalar u_hi = u + h;
		qaws_scalar hu;
		qaws_vec3 p_lo, p_hi;

		if (u_lo < QAWS_ZERO) u_lo = QAWS_ZERO;
		if (u_hi > QAWS_ONE) u_hi = QAWS_ONE;
		hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

		status = offset_eval_position(impl->base, impl->distance, u_lo, v, &p_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_hi, v, &p_hi);
		if (status != QAWS_STATUS_OK) return status;

		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
			out_result->du.x = (p_hi.x - p_lo.x) * inv2h;
			out_result->du.y = (p_hi.y - p_lo.y) * inv2h;
			out_result->du.z = (p_hi.z - p_lo.z) * inv2h;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}

		/* duu via central 2nd-order finite difference */
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar inv_h2 = QAWS_ONE / (hu * hu);
			out_result->duu.x = (p_hi.x - QAWS_LITERAL(2.0) * pos.x + p_lo.x) * inv_h2;
			out_result->duu.y = (p_hi.y - QAWS_LITERAL(2.0) * pos.y + p_lo.y) * inv_h2;
			out_result->duu.z = (p_hi.z - QAWS_LITERAL(2.0) * pos.z + p_lo.z) * inv_h2;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}
	}

	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_scalar v_lo = v - h;
		qaws_scalar v_hi = v + h;
		qaws_scalar hv;
		qaws_vec3 p_lo, p_hi;

		if (v_lo < QAWS_ZERO) v_lo = QAWS_ZERO;
		if (v_hi > QAWS_ONE) v_hi = QAWS_ONE;
		hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

		status = offset_eval_position(impl->base, impl->distance, u, v_lo, &p_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u, v_hi, &p_hi);
		if (status != QAWS_STATUS_OK) return status;

		if (eval_flags & QAWS_SURFACE_EVAL_DV)
		{
			qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
			out_result->dv.x = (p_hi.x - p_lo.x) * inv2h;
			out_result->dv.y = (p_hi.y - p_lo.y) * inv2h;
			out_result->dv.z = (p_hi.z - p_lo.z) * inv2h;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
		}

		/* dvv via central 2nd-order finite difference */
		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		{
			qaws_scalar inv_h2 = QAWS_ONE / (hv * hv);
			out_result->dvv.x = (p_hi.x - QAWS_LITERAL(2.0) * pos.x + p_lo.x) * inv_h2;
			out_result->dvv.y = (p_hi.y - QAWS_LITERAL(2.0) * pos.y + p_lo.y) * inv_h2;
			out_result->dvv.z = (p_hi.z - QAWS_LITERAL(2.0) * pos.z + p_lo.z) * inv_h2;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}
	}

	/* duv via finite difference of dv w.r.t. u */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_scalar u_lo = u - h;
		qaws_scalar u_hi = u + h;
		qaws_scalar hu;
		qaws_scalar hv_val = QAWS_LITERAL(1e-5);
		qaws_scalar v_lo, v_hi, hv;
		qaws_vec3 pv_lo_lo, pv_lo_hi, pv_hi_lo, pv_hi_hi;
		qaws_vec3 dv_lo, dv_hi;
		qaws_scalar inv2h;

		if (u_lo < QAWS_ZERO) u_lo = QAWS_ZERO;
		if (u_hi > QAWS_ONE) u_hi = QAWS_ONE;
		hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

		v_lo = v - hv_val;
		v_hi = v + hv_val;
		if (v_lo < QAWS_ZERO) v_lo = QAWS_ZERO;
		if (v_hi > QAWS_ONE) v_hi = QAWS_ONE;
		hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

		/* dv at u_lo and u_hi via finite differences in v */
		status = offset_eval_position(impl->base, impl->distance, u_lo, v_lo, &pv_lo_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_lo, v_hi, &pv_lo_hi);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_hi, v_lo, &pv_hi_lo);
		if (status != QAWS_STATUS_OK) return status;
		status = offset_eval_position(impl->base, impl->distance, u_hi, v_hi, &pv_hi_hi);
		if (status != QAWS_STATUS_OK) return status;

		{
			qaws_scalar inv2hv = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
			dv_lo.x = (pv_lo_hi.x - pv_lo_lo.x) * inv2hv;
			dv_lo.y = (pv_lo_hi.y - pv_lo_lo.y) * inv2hv;
			dv_lo.z = (pv_lo_hi.z - pv_lo_lo.z) * inv2hv;
			dv_hi.x = (pv_hi_hi.x - pv_hi_lo.x) * inv2hv;
			dv_hi.y = (pv_hi_hi.y - pv_hi_lo.y) * inv2hv;
			dv_hi.z = (pv_hi_hi.z - pv_hi_lo.z) * inv2hv;
		}

		inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
		out_result->duv.x = (dv_hi.x - dv_lo.x) * inv2h;
		out_result->duv.y = (dv_hi.y - dv_lo.y) * inv2h;
		out_result->duv.z = (dv_hi.z - dv_lo.z) * inv2h;
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

/* -------------------------------------------------------------------------- */
/*  Analytic offset jets                                                      */
/*                                                                            */
/*  O = S + d N,  N = n / |n|,  n = Su x Sv.                                  */
/*  Partials of N up to order 2 need base partials up to order 3; they are    */
/*  formed with bivariate Leibniz / quotient recurrences on dual numbers, so  */
/*  parameter and coordinate tangents come with the same evaluation.         */
/* -------------------------------------------------------------------------- */

#define OFFSET_JET_ORDER 2
#define OFFSET_JET_COUNT 6

static qaws_status offset_bijet(qaws_dual3 const* S, qaws_dual1 d, qaws_dual3* O)
{
	qaws_dual3 su[OFFSET_JET_COUNT], sv[OFFSET_JET_COUNT], n[OFFSET_JET_COUNT], N[OFFSET_JET_COUNT];
	qaws_dual1 q[OFFSET_JET_COUNT], l[OFFSET_JET_COUNT];
	qaws_scalar scale;
	unsigned int k;

	qaws_bijet_shift_u(S, OFFSET_JET_ORDER, su);
	qaws_bijet_shift_v(S, OFFSET_JET_ORDER, sv);
	qaws_bijet_cross(su, sv, OFFSET_JET_ORDER, n);
	qaws_bijet_dot(n, n, OFFSET_JET_ORDER, q);

	scale = qaws_v3_dot(su[0].v, su[0].v) * qaws_v3_dot(sv[0].v, sv[0].v);
	if (!(q[0].v > QAWS_EPSILON * QAWS_EPSILON * scale) || !(scale > QAWS_ZERO))
		return QAWS_STATUS_DEGENERATE_CURVE;

	qaws_bijet_sqrt(q, OFFSET_JET_ORDER, l);
	qaws_bijet_div(n, l, OFFSET_JET_ORDER, N);
	for (k = 0; k < OFFSET_JET_COUNT; k++)
		O[k] = qaws_dual3_add(S[k], qaws_dual3_scale(N[k], d));
	return QAWS_STATUS_OK;
}

/* Base jet up to third order; fails when the base has no analytic jets. */
static qaws_status offset_base_jet(qaws_surface const* base, qaws_scalar u, qaws_scalar v, qaws_dual3* S)
{
	qaws_surface_jet j;
	unsigned int k;
	qaws_status st = qaws_surface_eval_jet(base, u, v, QAWS_SJET_ORDER3, &j);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
		S[k] = qaws_dual3_const(j.d[k]);
	return QAWS_STATUS_OK;
}

static qaws_status offset_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_offset_impl const* impl =
		(qaws_surface_offset_impl const*)surface->impl;
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT], O[OFFSET_JET_COUNT];
	qaws_status st;

	/* Bases without analytic jets keep the sampled derivatives. */
	if (offset_base_jet(impl->base, u, v, S) != QAWS_STATUS_OK)
		return offset_surface_eval_sampled(surface, u, v, eval_flags, out_result);

	st = offset_bijet(S, qaws_dual1_const(impl->distance), O);
	if (st != QAWS_STATUS_OK)
		return offset_surface_eval_sampled(surface, u, v, eval_flags, out_result);

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = O[0].v;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL))
	{
		out_result->du = O[1].v;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}
	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL))
	{
		out_result->dv = O[2].v;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		out_result->duu = O[3].v;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		out_result->duv = O[4].v;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
	}
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv = O[5].v;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}
	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Differential rules                                                        */
/*                                                                            */
/*  Own field: offset_distance. Child 0: the base surface. Jets are limited   */
/*  to second order (third-order offset partials would need fourth-order      */
/*  base partials).                                                           */
/* -------------------------------------------------------------------------- */

#define OFFSET_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned int offset_describe_fields(qaws_surface const* surface,
	qaws_field_desc* out, unsigned int capacity)
{
	(void)surface;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_OFFSET_DISTANCE, QAWS_VALUE_SCALAR, 1,
			QAWS_DOMAIN_LENGTH, QAWS_CONSTRAINT_NONE, QAWS_DIFF_SMOOTH, OFFSET_DIFF_CAPS);
	return 1;
}

static qaws_status offset_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_surface_offset_impl const* impl = (qaws_surface_offset_impl const*)surface->impl;
	if (field != QAWS_FIELD_OFFSET_DISTANCE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_data = &impl->distance;
	*out_count = 1;
	*out_components = 1;
	return QAWS_STATUS_OK;
}

static unsigned int offset_children(qaws_surface const* surface, qaws_diff_child* out, unsigned int capacity)
{
	qaws_surface_offset_impl const* impl = (qaws_surface_offset_impl const*)surface->impl;
	if (capacity >= 1)
	{
		out[0].curve = NULL;
		out[0].surface = impl->base;
	}
	return 1;
}

static qaws_status offset_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_surface_jet* primal, qaws_surface_jet* tangent, qaws_surface_jet* tangent2)
{
	qaws_surface_offset_impl const* impl = (qaws_surface_offset_impl const*)surface->impl;
	qaws_field_view const* dv = views ? qaws_diff_views_find(views, QAWS_FIELD_OFFSET_DISTANCE) : NULL;
	qaws_surface_jet bp, bt, btt;
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT], O[OFFSET_JET_COUNT];
	qaws_scalar d_dot = QAWS_ZERO;
	unsigned int k;
	qaws_status st;

	if (channels & ~(unsigned int)QAWS_SJET_ORDER2)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (dv)
		qaws_internal_view_read(dv, 0, 1, &d_dot);

	if (tangent2)
		st = qaws_surface_eval_batch_tangent2(ctx, impl->base, &u, &v, &u_dot, &v_dot, 1, QAWS_SJET_ORDER3,
			qaws_internal_child_views(views, 0), &bp, &bt, &btt);
	else
		st = qaws_surface_eval_tangent(ctx, impl->base, u, v, u_dot, v_dot, QAWS_SJET_ORDER3,
			qaws_internal_child_views(views, 0), &bp, &bt);
	if (st != QAWS_STATUS_OK)
		return st;

	for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
		S[k] = qaws_dual3_make(bp.d[k], bt.d[k], tangent2 ? btt.d[k] : qaws_v3_zero());
	st = offset_bijet(S, qaws_dual1_make(impl->distance, d_dot, QAWS_ZERO), O);
	if (st != QAWS_STATUS_OK)
		return st;

	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));
	for (k = 0; k < OFFSET_JET_COUNT; k++)
	{
		if (!(channels & (1u << k)))
			continue;
		primal->d[k] = O[k].v;
		tangent->d[k] = O[k].t;
		if (tangent2)
			tangent2->d[k] = O[k].tt;
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

/* Exact pullback: one forward seed per base jet component and one for the
   distance, contracted with the output adjoint. */
static qaws_status offset_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* ybar, qaws_diff_views* views,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	qaws_surface_offset_impl const* impl = (qaws_surface_offset_impl const*)surface->impl;
	qaws_field_view* dv = views ? qaws_diff_views_find(views, QAWS_FIELD_OFFSET_DISTANCE) : NULL;
	qaws_surface_jet bj, bbar;
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT], O[OFFSET_JET_COUNT];
	qaws_scalar dbar = QAWS_ZERO;
	unsigned int k, c, q;
	qaws_status st;

	if (channels & ~(unsigned int)QAWS_SJET_ORDER2)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	st = qaws_surface_eval_jet(impl->base, u, v, QAWS_SJET_ORDER3, &bj);
	if (st != QAWS_STATUS_OK)
		return st;

	memset(&bbar, 0, sizeof(bbar));
	for (k = 0; k <= QAWS_SURFACE_JET_COUNT; k++)
	{
		for (c = 0; c < (k < QAWS_SURFACE_JET_COUNT ? 3u : 1u); c++)
		{
			qaws_dual1 d = qaws_dual1_make(impl->distance, k == QAWS_SURFACE_JET_COUNT ? QAWS_ONE : QAWS_ZERO, QAWS_ZERO);
			qaws_scalar acc = QAWS_ZERO;
			for (q = 0; q < QAWS_SURFACE_JET_COUNT; q++)
			{
				qaws_vec3 seed = qaws_v3_zero();
				if (q == k)
					seed = qaws_v3(c == 0 ? QAWS_ONE : QAWS_ZERO, c == 1 ? QAWS_ONE : QAWS_ZERO, c == 2 ? QAWS_ONE : QAWS_ZERO);
				S[q] = qaws_dual3_make(bj.d[q], seed, qaws_v3_zero());
			}
			st = offset_bijet(S, d, O);
			if (st != QAWS_STATUS_OK)
				return st;
			for (q = 0; q < OFFSET_JET_COUNT; q++)
				if (channels & (1u << q))
					acc += qaws_v3_dot(ybar->d[q], O[q].t);
			if (k == QAWS_SURFACE_JET_COUNT)
				dbar += acc;
			else if (c == 0)
				bbar.d[k].x += acc;
			else if (c == 1)
				bbar.d[k].y += acc;
			else
				bbar.d[k].z += acc;
		}
	}
	bbar.channels = QAWS_SJET_ORDER3;

	st = qaws_surface_eval_adjoint(ctx, impl->base, u, v, QAWS_SJET_ORDER3, &bbar,
		(qaws_diff_views*)qaws_internal_child_views(views, 0), u_adjoint, v_adjoint);
	if (st != QAWS_STATUS_OK)
		return st;
	if (dv)
		qaws_internal_view_add(dv, 0, 1, &dbar);
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const offset_surface_diff_vtable = {
	OFFSET_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	offset_describe_fields,
	offset_primal_field,
	NULL,
	NULL,
	offset_tangent,
	offset_adjoint,
	offset_children
};

static void offset_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int offset_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const offset_surface_vtable = {
	offset_surface_eval,
	offset_surface_destroy,
	offset_surface_is_rational,
	&offset_surface_diff_vtable
};

qaws_status qaws_surface_create_offset(
	qaws_surface_offset_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_offset_impl* impl;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->base) return QAWS_STATUS_INVALID_ARGUMENT;
	if (QAWS_FABS(desc->distance) < QAWS_LITERAL(1e-12))
		return QAWS_STATUS_INVALID_ARGUMENT;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_OFFSET,
		0, 0, u_range, v_range,
		&offset_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_offset_impl*)malloc(sizeof(qaws_surface_offset_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->base = desc->base;
	impl->distance = desc->distance;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
