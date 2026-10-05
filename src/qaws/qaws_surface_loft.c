#include "qaws_surface_loft.h"
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

typedef struct qaws_surface_loft_impl
{
	qaws_curve const** sections;   /* copied array of N section curve pointers */
	unsigned int section_count;
	qaws_scalar* v_params;         /* copied array of N v-parameter values */
	qaws_range* section_ranges;    /* copied array of N parameter ranges */
} qaws_surface_loft_impl;

/* Catmull-Rom interpolation through points at non-uniform parameters.
   pts[n_pts], params[n_pts], evaluate at parameter t.
   Returns blended position. If out_deriv is non-NULL, also returns derivative w.r.t. t. */
/* Evaluate the loft at (u,v) by evaluating each section at u
   then blending in v with Catmull-Rom. */
static qaws_status loft_surface_eval_sampled(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_loft_impl const* impl =
		(qaws_surface_loft_impl const*)surface->impl;
	unsigned int n = impl->section_count;
	unsigned int i;
	qaws_vec3* sec_pos = NULL;
	qaws_vec3* sec_du = NULL;
	qaws_vec3 pos, dv_vec, du_vec;
	qaws_status status;
	int need_pos, need_du, need_dv;

	memset(&pos, 0, sizeof(pos));
	memset(&dv_vec, 0, sizeof(dv_vec));
	memset(&du_vec, 0, sizeof(du_vec));

	need_pos = (eval_flags & QAWS_SURFACE_EVAL_POSITION) != 0;
	need_du = (eval_flags & QAWS_SURFACE_EVAL_DU) != 0;
	need_dv = (eval_flags & QAWS_SURFACE_EVAL_DV) != 0;

	/* For second derivatives and normal, we need du and dv */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		need_pos = 1;
		need_du = 1;
		need_dv = 1;
	}

	/* Allocate temporary arrays for section evaluations */
	sec_pos = (qaws_vec3*)malloc(n * sizeof(qaws_vec3));
	if (!sec_pos) return QAWS_STATUS_ALLOCATION_FAILURE;

	if (need_du)
	{
		sec_du = (qaws_vec3*)malloc(n * sizeof(qaws_vec3));
		if (!sec_du)
		{
			free(sec_pos);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}

	/* Evaluate each section curve at parameter u */
	for (i = 0; i < n; i++)
	{
		qaws_scalar s_range = impl->section_ranges[i].max_value
			- impl->section_ranges[i].min_value;
		qaws_scalar t_curve = impl->section_ranges[i].min_value + u * s_range;
		unsigned int cf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d cr;

		if (need_du)
			cf |= QAWS_EVAL_FLAG_D1;

		memset(&cr, 0, sizeof(cr));
		status = qaws_curve_evaluate_3d(impl->sections[i], t_curve, cf, &cr);
		if (status != QAWS_STATUS_OK)
		{
			free(sec_pos);
			if (sec_du) free(sec_du);
			return status;
		}

		sec_pos[i] = cr.position;
		if (need_du)
		{
			sec_du[i].x = cr.d1.x * s_range;
			sec_du[i].y = cr.d1.y * s_range;
			sec_du[i].z = cr.d1.z * s_range;
		}
	}

	/* Blend section positions in v via Catmull-Rom */
	if (need_pos || need_dv)
	{
		qaws_internal_surface_catmull_rom_blend(sec_pos, impl->v_params, n, v,
			&pos, need_dv ? &dv_vec : NULL);
	}

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv = dv_vec;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* dS/du: blend section derivatives in v with same weights */
	if (need_du)
	{
		qaws_internal_surface_catmull_rom_blend(sec_du, impl->v_params, n, v,
			&du_vec, NULL);
		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			out_result->du = du_vec;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}
	}

	/* Second derivatives via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_surface_eval_result r_lo, r_hi;

		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			loft_surface_eval_sampled(surface, u_lo, v, QAWS_SURFACE_EVAL_DU, &r_lo);
			loft_surface_eval_sampled(surface, u_hi, v, QAWS_SURFACE_EVAL_DU, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duu.x = (r_hi.du.x - r_lo.du.x) * inv2h;
				out_result->duu.y = (r_hi.du.y - r_lo.du.y) * inv2h;
				out_result->duu.z = (r_hi.du.z - r_lo.du.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		{
			qaws_scalar v_lo = v - h, v_hi = v + h;
			qaws_scalar hv;
			if (v_lo < 0) v_lo = 0;
			if (v_hi > 1) v_hi = 1;
			hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			loft_surface_eval_sampled(surface, u, v_lo, QAWS_SURFACE_EVAL_DV, &r_lo);
			loft_surface_eval_sampled(surface, u, v_hi, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
				out_result->dvv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->dvv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->dvv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}

		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			loft_surface_eval_sampled(surface, u_lo, v, QAWS_SURFACE_EVAL_DV, &r_lo);
			loft_surface_eval_sampled(surface, u_hi, v, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->duv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->duv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
		}

		if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
		{
			qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
		}
	}

	free(sec_pos);
	if (sec_du) free(sec_du);
	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Analytic jets and differential rules                                      */
/*                                                                            */
/*  S(u,v) = sum_k W_k(v) Q_k(t_k(u)): sections blended by Catmull-Rom        */
/*  weights at fixed v-parameters. Child i is section i; only four sections */
/*  act at a given v.                                                         */
/* -------------------------------------------------------------------------- */

#define LOFT_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned char const g_loft_jet_a[QAWS_SURFACE_JET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_loft_jet_b[QAWS_SURFACE_JET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

static qaws_scalar loft_power(qaws_scalar s, unsigned int a)
{
	qaws_scalar r = QAWS_ONE;
	while (a--)
		r *= s;
	return r;
}

static qaws_status loft_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_surface_jet* primal, qaws_surface_jet* tangent, qaws_surface_jet* tangent2)
{
	qaws_surface_loft_impl const* impl = (qaws_surface_loft_impl const*)surface->impl;
	unsigned int idx[4], k, ch;
	qaws_scalar w[QAWS_INTERNAL_CR_ROWS][4];
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT];

	qaws_internal_catmull_rom_weights(impl->v_params, impl->section_count, v, idx, w);
	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		S[ch] = qaws_dual3_const(qaws_v3_zero());

	for (k = 0; k < 4; k++)
	{
		qaws_range r = impl->section_ranges[idx[k]];
		qaws_scalar s = r.max_value - r.min_value;
		qaws_curve_jet_3d p, t, tt;
		qaws_status st;
		if (w[0][k] == QAWS_ZERO && w[1][k] == QAWS_ZERO && w[2][k] == QAWS_ZERO && w[3][k] == QAWS_ZERO)
			continue;
		st = qaws_internal_curve_tangent_any(ctx, impl->sections[idx[k]], r.min_value + u * s, s * u_dot, 0xFu,
			qaws_internal_child_views(views, idx[k]), &p, &t, tangent2 ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			unsigned int a = g_loft_jet_a[ch], b = g_loft_jet_b[ch];
			qaws_scalar sa = loft_power(s, a);
			qaws_dual3 q = qaws_dual3_make(qaws_v3_scale(p.d[a], sa), qaws_v3_scale(t.d[a], sa),
				tangent2 ? qaws_v3_scale(tt.d[a], sa) : qaws_v3_zero());
			qaws_dual1 wd = qaws_dual1_make(w[b][k], w[b + 1][k] * v_dot, w[b + 2][k] * v_dot * v_dot);
			S[ch] = qaws_dual3_add(S[ch], qaws_dual3_scale(q, wd));
		}
	}

	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));
	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		if (!(channels & (1u << ch)))
			continue;
		primal->d[ch] = S[ch].v;
		tangent->d[ch] = S[ch].t;
		if (tangent2)
			tangent2->d[ch] = S[ch].tt;
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

static qaws_status loft_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* ybar, qaws_diff_views* views,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	qaws_surface_loft_impl const* impl = (qaws_surface_loft_impl const*)surface->impl;
	unsigned int idx[4], k, ch;
	qaws_scalar w[QAWS_INTERNAL_CR_ROWS][4];

	qaws_internal_catmull_rom_weights(impl->v_params, impl->section_count, v, idx, w);
	for (k = 0; k < 4; k++)
	{
		qaws_range r = impl->section_ranges[idx[k]];
		qaws_scalar s = r.max_value - r.min_value, t = r.min_value + u * s, tbar = QAWS_ZERO;
		qaws_curve_jet_3d jbar;
		qaws_status st;
		if (w[0][k] == QAWS_ZERO && w[1][k] == QAWS_ZERO && w[2][k] == QAWS_ZERO && w[3][k] == QAWS_ZERO)
			continue;
		memset(&jbar, 0, sizeof(jbar));
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			unsigned int a = g_loft_jet_a[ch], b = g_loft_jet_b[ch];
			if (!(channels & (1u << ch)))
				continue;
			jbar.d[a] = qaws_v3_axpy(jbar.d[a], ybar->d[ch], w[b][k] * loft_power(s, a));
		}
		jbar.channels = 0xFu;
		if (v_adjoint)
		{
			/* d/dv of the weight: v enters only through W_k. */
			qaws_curve_jet_3d p, tg;
			st = qaws_internal_curve_tangent_any(ctx, impl->sections[idx[k]], t, QAWS_ZERO, 0xFu, NULL, &p, &tg, NULL);
			if (st != QAWS_STATUS_OK)
				return st;
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				unsigned int a = g_loft_jet_a[ch], b = g_loft_jet_b[ch];
				if (!(channels & (1u << ch)))
					continue;
				*v_adjoint += w[b + 1][k] * loft_power(s, a) * qaws_v3_dot(ybar->d[ch], p.d[a]);
			}
		}
		st = qaws_internal_curve_adjoint_any(ctx, impl->sections[idx[k]], t, 0xFu, &jbar,
			(qaws_diff_views*)qaws_internal_child_views(views, idx[k]), u_adjoint ? &tbar : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		if (u_adjoint)
			*u_adjoint += s * tbar;
	}
	return QAWS_STATUS_OK;
}

static unsigned int loft_describe_fields(qaws_surface const* surface, qaws_field_desc* out, unsigned int capacity)
{
	(void)surface; (void)out; (void)capacity;
	return 0;
}

static qaws_status loft_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	(void)surface; (void)field; (void)out_data; (void)out_count; (void)out_components;
	return QAWS_STATUS_INVALID_ARGUMENT;
}

static unsigned int loft_children(qaws_surface const* surface, qaws_diff_child* out, unsigned int capacity)
{
	qaws_surface_loft_impl const* impl = (qaws_surface_loft_impl const*)surface->impl;
	unsigned int i;
	for (i = 0; i < impl->section_count && i < capacity; i++)
	{
		out[i].curve = impl->sections[i];
		out[i].surface = NULL;
	}
	return impl->section_count;
}

static qaws_surface_diff_vtable const loft_surface_diff_vtable = {
	LOFT_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	loft_describe_fields,
	loft_primal_field,
	NULL,
	NULL,
	loft_tangent,
	loft_adjoint,
	loft_children
};

static qaws_status loft_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_jet p, t;

	/* Sections without analytic jets keep the sampled derivatives. */
	if (loft_tangent(NULL, surface, u, v, QAWS_ZERO, QAWS_ZERO, QAWS_SJET_ORDER2, NULL, &p, &t, NULL) != QAWS_STATUS_OK)
		return loft_surface_eval_sampled(surface, u, v, eval_flags, out_result);

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION) { out_result->position = p.d[0]; out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION; }
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL)) { out_result->du = p.d[1]; out_result->valid_flags |= QAWS_SURFACE_EVAL_DU; }
	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL)) { out_result->dv = p.d[2]; out_result->valid_flags |= QAWS_SURFACE_EVAL_DV; }
	if (eval_flags & QAWS_SURFACE_EVAL_DUU) { out_result->duu = p.d[3]; out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU; }
	if (eval_flags & QAWS_SURFACE_EVAL_DUV) { out_result->duv = p.d[4]; out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV; }
	if (eval_flags & QAWS_SURFACE_EVAL_DVV) { out_result->dvv = p.d[5]; out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV; }
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}
	return QAWS_STATUS_OK;
}

static void loft_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_loft_impl* li = (qaws_surface_loft_impl*)impl;
	if (li)
	{
		if (li->sections) free(li->sections);
		if (li->v_params) free(li->v_params);
		if (li->section_ranges) free(li->section_ranges);
	}
	qaws_internal_dealloc(allocator, impl);
}

static int loft_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const loft_surface_vtable = {
	loft_surface_eval,
	loft_surface_destroy,
	loft_surface_is_rational,
	&loft_surface_diff_vtable
};

qaws_status qaws_surface_create_loft(
	qaws_surface_loft_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_loft_impl* impl;
	qaws_range u_range, v_range;
	unsigned int i;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->sections) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->section_count < 2) return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate all sections are 3D */
	for (i = 0; i < desc->section_count; i++)
	{
		if (!desc->sections[i]) return QAWS_STATUS_INVALID_ARGUMENT;
		if (qaws_curve_get_dimension(desc->sections[i]) != QAWS_DIMENSION_3D)
			return QAWS_STATUS_INVALID_DIMENSION;
	}

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_LOFT,
		0, 0, u_range, v_range,
		&loft_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_loft_impl*)malloc(sizeof(qaws_surface_loft_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(impl, 0, sizeof(qaws_surface_loft_impl));

	impl->section_count = desc->section_count;

	/* Copy sections array */
	impl->sections = (qaws_curve const**)malloc(
		desc->section_count * sizeof(qaws_curve const*));
	if (!impl->sections)
	{
		free(impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->sections, desc->sections,
		desc->section_count * sizeof(qaws_curve const*));

	/* Copy or generate v_parameters */
	impl->v_params = (qaws_scalar*)malloc(
		desc->section_count * sizeof(qaws_scalar));
	if (!impl->v_params)
	{
		free(impl->sections);
		free(impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	if (desc->v_parameters)
	{
		memcpy(impl->v_params, desc->v_parameters,
			desc->section_count * sizeof(qaws_scalar));
	}
	else
	{
		/* Uniform spacing */
		for (i = 0; i < desc->section_count; i++)
		{
			impl->v_params[i] = (qaws_scalar)i
				/ (qaws_scalar)(desc->section_count - 1);
		}
	}

	/* Cache section parameter ranges */
	impl->section_ranges = (qaws_range*)malloc(
		desc->section_count * sizeof(qaws_range));
	if (!impl->section_ranges)
	{
		free(impl->v_params);
		free(impl->sections);
		free(impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < desc->section_count; i++)
	{
		impl->section_ranges[i] = qaws_curve_get_parameter_range(desc->sections[i]);
	}

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
