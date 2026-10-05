#include "qaws_surface_gordon.h"
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

typedef struct qaws_surface_gordon_impl
{
	qaws_curve const** u_curves;    /* copied array of M u-curve pointers */
	unsigned int u_curve_count;      /* M */
	qaws_scalar* v_params;           /* copied array of M v-parameter values */
	qaws_range* u_curve_ranges;      /* copied array of M parameter ranges */

	qaws_curve const** v_curves;    /* copied array of N v-curve pointers */
	unsigned int v_curve_count;      /* N */
	qaws_scalar* u_params;           /* copied array of N u-parameter values */
	qaws_range* v_curve_ranges;      /* copied array of N parameter ranges */

	qaws_vec3* Q;                    /* precomputed intersection grid Q[M*N], row-major */
} qaws_surface_gordon_impl;

/* Catmull-Rom interpolation through points at non-uniform parameters.
   pts[n_pts], params[n_pts], evaluate at parameter t.
   Returns blended position. If out_deriv is non-NULL, also returns derivative w.r.t. t. */
/* Evaluate the tensor product term T(u,v) by blending Q[M*N] grid.
   First blend each row (fixed i) in u, then blend those results in v. */
static void tensor_product_eval(
	qaws_surface_gordon_impl const* impl,
	qaws_scalar u, qaws_scalar v,
	qaws_vec3* out_pos, qaws_vec3* out_du, qaws_vec3* out_dv)
{
	unsigned int M = impl->u_curve_count;
	unsigned int N = impl->v_curve_count;
	unsigned int i;
	qaws_vec3* row_blended;
	qaws_vec3* row_du;

	row_blended = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
	if (!row_blended) return;

	row_du = NULL;
	if (out_du)
	{
		row_du = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
		if (!row_du)
		{
			free(row_blended);
			return;
		}
	}

	/* For each u-curve index i, blend the N intersection points Q[i][0..N-1] at u */
	for (i = 0; i < M; i++)
	{
		qaws_vec3 const* row = &impl->Q[i * N];
		qaws_internal_surface_catmull_rom_blend(row, impl->u_params, N, u,
			&row_blended[i], row_du ? &row_du[i] : NULL);
	}

	/* Blend the M results in v */
	qaws_internal_surface_catmull_rom_blend(row_blended, impl->v_params, M, v,
		out_pos, out_dv);

	/* For dT/du: blend the M du-results in v */
	if (out_du)
	{
		qaws_internal_surface_catmull_rom_blend(row_du, impl->v_params, M, v,
			out_du, NULL);
	}

	free(row_blended);
	if (row_du) free(row_du);
}

/* Gordon surface: S(u,v) = L_u(u,v) + L_v(u,v) - T(u,v)
   L_u: evaluate each u-curve at u, blend in v
   L_v: evaluate each v-curve at v, blend in u
   T: tensor product of intersection grid */
static qaws_status gordon_surface_eval_sampled(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_gordon_impl const* impl =
		(qaws_surface_gordon_impl const*)surface->impl;
	unsigned int M = impl->u_curve_count;
	unsigned int N = impl->v_curve_count;
	unsigned int i;
	qaws_vec3* u_pts = NULL;
	qaws_vec3* u_du = NULL;
	qaws_vec3* v_pts = NULL;
	qaws_vec3* v_dv = NULL;
	qaws_vec3 Lu_pos, Lu_dv, Lu_du;
	qaws_vec3 Lv_pos, Lv_du, Lv_dv;
	qaws_vec3 T_pos, T_du, T_dv;
	qaws_status status;
	int need_pos, need_du, need_dv;

	need_pos = (eval_flags & QAWS_SURFACE_EVAL_POSITION) != 0;
	need_du = (eval_flags & QAWS_SURFACE_EVAL_DU) != 0;
	need_dv = (eval_flags & QAWS_SURFACE_EVAL_DV) != 0;

	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		need_pos = 1;
		need_du = 1;
		need_dv = 1;
	}

	/* Initialize vectors */
	memset(&Lu_pos, 0, sizeof(Lu_pos));
	memset(&Lu_dv, 0, sizeof(Lu_dv));
	memset(&Lu_du, 0, sizeof(Lu_du));
	memset(&Lv_pos, 0, sizeof(Lv_pos));
	memset(&Lv_du, 0, sizeof(Lv_du));
	memset(&Lv_dv, 0, sizeof(Lv_dv));
	memset(&T_pos, 0, sizeof(T_pos));
	memset(&T_du, 0, sizeof(T_du));
	memset(&T_dv, 0, sizeof(T_dv));

	/* === L_u(u,v): evaluate each u-curve at u, blend in v === */
	u_pts = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
	if (!u_pts) return QAWS_STATUS_ALLOCATION_FAILURE;

	u_du = NULL;
	if (need_du)
	{
		u_du = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
		if (!u_du)
		{
			free(u_pts);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}

	for (i = 0; i < M; i++)
	{
		qaws_scalar s_range = impl->u_curve_ranges[i].max_value
			- impl->u_curve_ranges[i].min_value;
		qaws_scalar t_curve = impl->u_curve_ranges[i].min_value + u * s_range;
		unsigned int cf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d cr;

		if (need_du)
			cf |= QAWS_EVAL_FLAG_D1;

		memset(&cr, 0, sizeof(cr));
		status = qaws_curve_evaluate_3d(impl->u_curves[i], t_curve, cf, &cr);
		if (status != QAWS_STATUS_OK)
		{
			free(u_pts);
			if (u_du) free(u_du);
			return status;
		}

		u_pts[i] = cr.position;
		if (need_du)
		{
			u_du[i].x = cr.d1.x * s_range;
			u_du[i].y = cr.d1.y * s_range;
			u_du[i].z = cr.d1.z * s_range;
		}
	}

	qaws_internal_surface_catmull_rom_blend(u_pts, impl->v_params, M, v,
		&Lu_pos, need_dv ? &Lu_dv : NULL);

	if (need_du)
	{
		qaws_internal_surface_catmull_rom_blend(u_du, impl->v_params, M, v,
			&Lu_du, NULL);
	}

	free(u_pts);
	if (u_du) free(u_du);

	/* === L_v(u,v): evaluate each v-curve at v, blend in u === */
	v_pts = (qaws_vec3*)malloc(N * sizeof(qaws_vec3));
	if (!v_pts) return QAWS_STATUS_ALLOCATION_FAILURE;

	v_dv = NULL;
	if (need_dv)
	{
		v_dv = (qaws_vec3*)malloc(N * sizeof(qaws_vec3));
		if (!v_dv)
		{
			free(v_pts);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}

	for (i = 0; i < N; i++)
	{
		qaws_scalar s_range = impl->v_curve_ranges[i].max_value
			- impl->v_curve_ranges[i].min_value;
		qaws_scalar t_curve = impl->v_curve_ranges[i].min_value + v * s_range;
		unsigned int cf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d cr;

		if (need_dv)
			cf |= QAWS_EVAL_FLAG_D1;

		memset(&cr, 0, sizeof(cr));
		status = qaws_curve_evaluate_3d(impl->v_curves[i], t_curve, cf, &cr);
		if (status != QAWS_STATUS_OK)
		{
			free(v_pts);
			if (v_dv) free(v_dv);
			return status;
		}

		v_pts[i] = cr.position;
		if (need_dv)
		{
			v_dv[i].x = cr.d1.x * s_range;
			v_dv[i].y = cr.d1.y * s_range;
			v_dv[i].z = cr.d1.z * s_range;
		}
	}

	qaws_internal_surface_catmull_rom_blend(v_pts, impl->u_params, N, u,
		&Lv_pos, need_du ? &Lv_du : NULL);

	if (need_dv)
	{
		qaws_internal_surface_catmull_rom_blend(v_dv, impl->u_params, N, u,
			&Lv_dv, NULL);
	}

	free(v_pts);
	if (v_dv) free(v_dv);

	/* === T(u,v): tensor product of intersection grid === */
	tensor_product_eval(impl, u, v, &T_pos,
		need_du ? &T_du : NULL,
		need_dv ? &T_dv : NULL);

	/* === Combine: S = L_u + L_v - T === */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = Lu_pos.x + Lv_pos.x - T_pos.x;
		out_result->position.y = Lu_pos.y + Lv_pos.y - T_pos.y;
		out_result->position.z = Lu_pos.z + Lv_pos.z - T_pos.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = Lu_du.x + Lv_du.x - T_du.x;
		out_result->du.y = Lu_du.y + Lv_du.y - T_du.y;
		out_result->du.z = Lu_du.z + Lv_du.z - T_du.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = Lu_dv.x + Lv_dv.x - T_dv.x;
		out_result->dv.y = Lu_dv.y + Lv_dv.y - T_dv.y;
		out_result->dv.z = Lu_dv.z + Lv_dv.z - T_dv.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
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
			gordon_surface_eval_sampled(surface, u_lo, v, QAWS_SURFACE_EVAL_DU, &r_lo);
			gordon_surface_eval_sampled(surface, u_hi, v, QAWS_SURFACE_EVAL_DU, &r_hi);

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
			gordon_surface_eval_sampled(surface, u, v_lo, QAWS_SURFACE_EVAL_DV, &r_lo);
			gordon_surface_eval_sampled(surface, u, v_hi, QAWS_SURFACE_EVAL_DV, &r_hi);

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
			gordon_surface_eval_sampled(surface, u_lo, v, QAWS_SURFACE_EVAL_DV, &r_lo);
			gordon_surface_eval_sampled(surface, u_hi, v, QAWS_SURFACE_EVAL_DV, &r_hi);

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

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Analytic jets and differential rules                                      */
/*                                                                            */
/*  S = Lu + Lv - T with Catmull-Rom weights Wv (over v_params) and Wu (over  */
/*  u_params):                                                                */
/*    Lu = sum_i Wv_i(v) U_i(u),  Lv = sum_j Wu_j(u) V_j(v),                  */
/*    T  = sum_ij Wv_i(v) Wu_j(u) Q_ij,  Q_ij = U_i(u_params[j]) (live).      */
/*  Children 0..M-1 are the u-curves, M..M+N-1 the v-curves.                  */
/* -------------------------------------------------------------------------- */

#define GORDON_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned char const g_gordon_jet_a[QAWS_SURFACE_JET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_gordon_jet_b[QAWS_SURFACE_JET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

static qaws_scalar gordon_power(qaws_scalar s, unsigned int a)
{
	qaws_scalar r = QAWS_ONE;
	while (a--)
		r *= s;
	return r;
}

static int gordon_slot_active(qaws_scalar (*w)[4], unsigned int k)
{
	return w[0][k] != QAWS_ZERO || w[1][k] != QAWS_ZERO || w[2][k] != QAWS_ZERO || w[3][k] != QAWS_ZERO;
}

/* Dual weight of derivative order d along a coordinate moving at rate c. */
static qaws_dual1 gordon_weight(qaws_scalar (*w)[4], unsigned int d, unsigned int k, qaws_scalar c)
{
	return qaws_dual1_make(w[d][k], w[d + 1][k] * c, w[d + 2][k] * c * c);
}

static qaws_status gordon_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_surface_jet* primal, qaws_surface_jet* tangent, qaws_surface_jet* tangent2)
{
	qaws_surface_gordon_impl const* impl = (qaws_surface_gordon_impl const*)surface->impl;
	unsigned int M = impl->u_curve_count, N = impl->v_curve_count;
	unsigned int iv[4], ju[4], k, l, ch;
	qaws_scalar wv[QAWS_INTERNAL_CR_ROWS][4], wu[QAWS_INTERNAL_CR_ROWS][4];
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT];
	qaws_curve_jet_3d p, t, tt;
	qaws_status st;

	qaws_internal_catmull_rom_weights(impl->v_params, M, v, iv, wv);
	qaws_internal_catmull_rom_weights(impl->u_params, N, u, ju, wu);
	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		S[ch] = qaws_dual3_const(qaws_v3_zero());

	/* Lu and T (both driven by the u-curves) */
	for (k = 0; k < 4; k++)
	{
		unsigned int i = iv[k];
		qaws_range r = impl->u_curve_ranges[i];
		qaws_scalar s = r.max_value - r.min_value;
		if (!gordon_slot_active(wv, k))
			continue;
		st = qaws_internal_curve_tangent_any(ctx, impl->u_curves[i], r.min_value + u * s, s * u_dot, 0xFu,
			qaws_internal_child_views(views, i), &p, &t, tangent2 ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
			qaws_scalar sa = gordon_power(s, a);
			qaws_dual3 q = qaws_dual3_make(qaws_v3_scale(p.d[a], sa), qaws_v3_scale(t.d[a], sa),
				tangent2 ? qaws_v3_scale(tt.d[a], sa) : qaws_v3_zero());
			S[ch] = qaws_dual3_add(S[ch], qaws_dual3_scale(q, gordon_weight(wv, b, k, v_dot)));
		}
		for (l = 0; l < 4; l++)
		{
			unsigned int j = ju[l];
			qaws_dual3 Q;
			if (!gordon_slot_active(wu, l))
				continue;
			st = qaws_internal_curve_tangent_any(ctx, impl->u_curves[i], r.min_value + impl->u_params[j] * s,
				QAWS_ZERO, QAWS_EVAL_FLAG_POSITION, qaws_internal_child_views(views, i), &p, &t, tangent2 ? &tt : NULL);
			if (st != QAWS_STATUS_OK)
				return st;
			Q = qaws_dual3_make(p.d[0], t.d[0], tangent2 ? tt.d[0] : qaws_v3_zero());
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
				qaws_dual1 wgt = qaws_dual1_mul(gordon_weight(wv, b, k, v_dot), gordon_weight(wu, a, l, u_dot));
				S[ch] = qaws_dual3_sub(S[ch], qaws_dual3_scale(Q, wgt));
			}
		}
	}

	/* Lv (driven by the v-curves) */
	for (l = 0; l < 4; l++)
	{
		unsigned int j = ju[l];
		qaws_range r = impl->v_curve_ranges[j];
		qaws_scalar s = r.max_value - r.min_value;
		if (!gordon_slot_active(wu, l))
			continue;
		st = qaws_internal_curve_tangent_any(ctx, impl->v_curves[j], r.min_value + v * s, s * v_dot, 0xFu,
			qaws_internal_child_views(views, M + j), &p, &t, tangent2 ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
			qaws_scalar sb = gordon_power(s, b);
			qaws_dual3 q = qaws_dual3_make(qaws_v3_scale(p.d[b], sb), qaws_v3_scale(t.d[b], sb),
				tangent2 ? qaws_v3_scale(tt.d[b], sb) : qaws_v3_zero());
			S[ch] = qaws_dual3_add(S[ch], qaws_dual3_scale(q, gordon_weight(wu, a, l, u_dot)));
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

static qaws_status gordon_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* ybar, qaws_diff_views* views,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	qaws_surface_gordon_impl const* impl = (qaws_surface_gordon_impl const*)surface->impl;
	unsigned int M = impl->u_curve_count, N = impl->v_curve_count;
	unsigned int iv[4], ju[4], k, l, ch;
	qaws_scalar wv[QAWS_INTERNAL_CR_ROWS][4], wu[QAWS_INTERNAL_CR_ROWS][4];
	qaws_curve_jet_3d p, t;
	qaws_status st;
	(void)N;

	qaws_internal_catmull_rom_weights(impl->v_params, M, v, iv, wv);
	qaws_internal_catmull_rom_weights(impl->u_params, impl->v_curve_count, u, ju, wu);

	for (k = 0; k < 4; k++)
	{
		unsigned int i = iv[k];
		qaws_range r = impl->u_curve_ranges[i];
		qaws_scalar s = r.max_value - r.min_value, tbar = QAWS_ZERO;
		qaws_curve_jet_3d jbar;
		qaws_diff_views* cv = (qaws_diff_views*)qaws_internal_child_views(views, i);
		if (!gordon_slot_active(wv, k))
			continue;

		/* Lu: weights Wv_i^(b) on U_i^(a) s^a */
		memset(&jbar, 0, sizeof(jbar));
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
			if (channels & (1u << ch))
				jbar.d[a] = qaws_v3_axpy(jbar.d[a], ybar->d[ch], wv[b][k] * gordon_power(s, a));
		}
		jbar.channels = 0xFu;
		if (v_adjoint)
		{
			st = qaws_internal_curve_tangent_any(ctx, impl->u_curves[i], r.min_value + u * s, QAWS_ZERO, 0xFu, NULL, &p, &t, NULL);
			if (st != QAWS_STATUS_OK)
				return st;
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
				if (channels & (1u << ch))
					*v_adjoint += wv[b + 1][k] * gordon_power(s, a) * qaws_v3_dot(ybar->d[ch], p.d[a]);
			}
		}
		st = qaws_internal_curve_adjoint_any(ctx, impl->u_curves[i], r.min_value + u * s, 0xFu, &jbar, cv,
			u_adjoint ? &tbar : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		if (u_adjoint)
			*u_adjoint += s * tbar;

		/* T: -Wv_i^(b) Wu_j^(a) on Q_ij */
		for (l = 0; l < 4; l++)
		{
			unsigned int j = ju[l];
			qaws_scalar tq = r.min_value + impl->u_params[j] * s;
			qaws_curve_jet_3d qbar;
			if (!gordon_slot_active(wu, l))
				continue;
			memset(&qbar, 0, sizeof(qbar));
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
				if (channels & (1u << ch))
					qbar.d[0] = qaws_v3_axpy(qbar.d[0], ybar->d[ch], -wv[b][k] * wu[a][l]);
			}
			qbar.channels = QAWS_EVAL_FLAG_POSITION;
			if (u_adjoint || v_adjoint)
			{
				st = qaws_internal_curve_tangent_any(ctx, impl->u_curves[i], tq, QAWS_ZERO, QAWS_EVAL_FLAG_POSITION, NULL, &p, &t, NULL);
				if (st != QAWS_STATUS_OK)
					return st;
				for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
				{
					unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
					qaws_scalar d = qaws_v3_dot(ybar->d[ch], p.d[0]);
					if (!(channels & (1u << ch)))
						continue;
					if (u_adjoint)
						*u_adjoint -= wv[b][k] * wu[a + 1][l] * d;
					if (v_adjoint)
						*v_adjoint -= wv[b + 1][k] * wu[a][l] * d;
				}
			}
			st = qaws_internal_curve_adjoint_any(ctx, impl->u_curves[i], tq, QAWS_EVAL_FLAG_POSITION, &qbar, cv, NULL);
			if (st != QAWS_STATUS_OK)
				return st;
		}
	}

	for (l = 0; l < 4; l++)
	{
		unsigned int j = ju[l];
		qaws_range r = impl->v_curve_ranges[j];
		qaws_scalar s = r.max_value - r.min_value, tbar = QAWS_ZERO;
		qaws_curve_jet_3d jbar;
		if (!gordon_slot_active(wu, l))
			continue;
		memset(&jbar, 0, sizeof(jbar));
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
			if (channels & (1u << ch))
				jbar.d[b] = qaws_v3_axpy(jbar.d[b], ybar->d[ch], wu[a][l] * gordon_power(s, b));
		}
		jbar.channels = 0xFu;
		if (u_adjoint)
		{
			st = qaws_internal_curve_tangent_any(ctx, impl->v_curves[j], r.min_value + v * s, QAWS_ZERO, 0xFu, NULL, &p, &t, NULL);
			if (st != QAWS_STATUS_OK)
				return st;
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				unsigned int a = g_gordon_jet_a[ch], b = g_gordon_jet_b[ch];
				if (channels & (1u << ch))
					*u_adjoint += wu[a + 1][l] * gordon_power(s, b) * qaws_v3_dot(ybar->d[ch], p.d[b]);
			}
		}
		st = qaws_internal_curve_adjoint_any(ctx, impl->v_curves[j], r.min_value + v * s, 0xFu, &jbar,
			(qaws_diff_views*)qaws_internal_child_views(views, M + j), v_adjoint ? &tbar : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		if (v_adjoint)
			*v_adjoint += s * tbar;
	}
	return QAWS_STATUS_OK;
}

static unsigned int gordon_describe_fields(qaws_surface const* surface, qaws_field_desc* out, unsigned int capacity)
{
	(void)surface; (void)out; (void)capacity;
	return 0;
}

static qaws_status gordon_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	(void)surface; (void)field; (void)out_data; (void)out_count; (void)out_components;
	return QAWS_STATUS_INVALID_ARGUMENT;
}

static unsigned int gordon_children(qaws_surface const* surface, qaws_diff_child* out, unsigned int capacity)
{
	qaws_surface_gordon_impl const* impl = (qaws_surface_gordon_impl const*)surface->impl;
	unsigned int i, M = impl->u_curve_count, N = impl->v_curve_count;
	for (i = 0; i < M + N && i < capacity; i++)
	{
		out[i].curve = i < M ? impl->u_curves[i] : impl->v_curves[i - M];
		out[i].surface = NULL;
	}
	return M + N;
}

static qaws_surface_diff_vtable const gordon_surface_diff_vtable = {
	GORDON_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	gordon_describe_fields,
	gordon_primal_field,
	NULL,
	NULL,
	gordon_tangent,
	gordon_adjoint,
	gordon_children
};

static qaws_status gordon_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_jet p, t;

	/* Curves without analytic jets keep the sampled derivatives. */
	if (gordon_tangent(NULL, surface, u, v, QAWS_ZERO, QAWS_ZERO, QAWS_SJET_ORDER2, NULL, &p, &t, NULL) != QAWS_STATUS_OK)
		return gordon_surface_eval_sampled(surface, u, v, eval_flags, out_result);

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

static void gordon_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_gordon_impl* gi = (qaws_surface_gordon_impl*)impl;
	if (gi)
	{
		if (gi->u_curves) free(gi->u_curves);
		if (gi->v_params) free(gi->v_params);
		if (gi->u_curve_ranges) free(gi->u_curve_ranges);
		if (gi->v_curves) free(gi->v_curves);
		if (gi->u_params) free(gi->u_params);
		if (gi->v_curve_ranges) free(gi->v_curve_ranges);
		if (gi->Q) free(gi->Q);
	}
	qaws_internal_dealloc(allocator, impl);
}

static int gordon_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const gordon_surface_vtable = {
	gordon_surface_eval,
	gordon_surface_destroy,
	gordon_surface_is_rational,
	&gordon_surface_diff_vtable
};

qaws_status qaws_surface_create_gordon(
	qaws_surface_gordon_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_gordon_impl* impl;
	qaws_range u_range, v_range;
	unsigned int M, N, i, j;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->u_curves || !desc->v_curves) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->v_params || !desc->u_params) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->u_curve_count < 2) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->v_curve_count < 2) return QAWS_STATUS_INVALID_ARGUMENT;

	M = desc->u_curve_count;
	N = desc->v_curve_count;

	/* Validate all curves are 3D */
	for (i = 0; i < M; i++)
	{
		if (!desc->u_curves[i]) return QAWS_STATUS_INVALID_ARGUMENT;
		if (qaws_curve_get_dimension(desc->u_curves[i]) != QAWS_DIMENSION_3D)
			return QAWS_STATUS_INVALID_DIMENSION;
	}
	for (i = 0; i < N; i++)
	{
		if (!desc->v_curves[i]) return QAWS_STATUS_INVALID_ARGUMENT;
		if (qaws_curve_get_dimension(desc->v_curves[i]) != QAWS_DIMENSION_3D)
			return QAWS_STATUS_INVALID_DIMENSION;
	}

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_GORDON,
		0, 0, u_range, v_range,
		&gordon_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_gordon_impl*)malloc(sizeof(qaws_surface_gordon_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(impl, 0, sizeof(qaws_surface_gordon_impl));

	impl->u_curve_count = M;
	impl->v_curve_count = N;

	/* Copy u_curves array */
	impl->u_curves = (qaws_curve const**)malloc(M * sizeof(qaws_curve const*));
	if (!impl->u_curves) goto fail;
	memcpy(impl->u_curves, desc->u_curves, M * sizeof(qaws_curve const*));

	/* Copy v_params */
	impl->v_params = (qaws_scalar*)malloc(M * sizeof(qaws_scalar));
	if (!impl->v_params) goto fail;
	memcpy(impl->v_params, desc->v_params, M * sizeof(qaws_scalar));

	/* Cache u-curve parameter ranges */
	impl->u_curve_ranges = (qaws_range*)malloc(M * sizeof(qaws_range));
	if (!impl->u_curve_ranges) goto fail;
	for (i = 0; i < M; i++)
		impl->u_curve_ranges[i] = qaws_curve_get_parameter_range(desc->u_curves[i]);

	/* Copy v_curves array */
	impl->v_curves = (qaws_curve const**)malloc(N * sizeof(qaws_curve const*));
	if (!impl->v_curves) goto fail;
	memcpy(impl->v_curves, desc->v_curves, N * sizeof(qaws_curve const*));

	/* Copy u_params */
	impl->u_params = (qaws_scalar*)malloc(N * sizeof(qaws_scalar));
	if (!impl->u_params) goto fail;
	memcpy(impl->u_params, desc->u_params, N * sizeof(qaws_scalar));

	/* Cache v-curve parameter ranges */
	impl->v_curve_ranges = (qaws_range*)malloc(N * sizeof(qaws_range));
	if (!impl->v_curve_ranges) goto fail;
	for (i = 0; i < N; i++)
		impl->v_curve_ranges[i] = qaws_curve_get_parameter_range(desc->v_curves[i]);

	/* Precompute intersection grid Q[M*N].
	   Q[i][j] = u_curve[i] evaluated at the u-parameter corresponding to u_params[j].
	   u_curve[i] runs along u, so evaluate at:
	     t = range.min + u_params[j] * (range.max - range.min) */
	impl->Q = (qaws_vec3*)malloc(M * N * sizeof(qaws_vec3));
	if (!impl->Q) goto fail;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			qaws_scalar s_range = impl->u_curve_ranges[i].max_value
				- impl->u_curve_ranges[i].min_value;
			qaws_scalar t_curve = impl->u_curve_ranges[i].min_value
				+ desc->u_params[j] * s_range;
			qaws_eval_result_3d cr;
			qaws_status status;

			memset(&cr, 0, sizeof(cr));
			status = qaws_curve_evaluate_3d(desc->u_curves[i], t_curve,
				QAWS_EVAL_FLAG_POSITION, &cr);
			if (status != QAWS_STATUS_OK)
			{
				gordon_surface_destroy(impl, NULL);
				qaws_internal_surface_free(surface);
				return status;
			}
			impl->Q[i * N + j] = cr.position;
		}
	}

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;

fail:
	gordon_surface_destroy(impl, NULL);
	qaws_internal_surface_free(surface);
	return QAWS_STATUS_ALLOCATION_FAILURE;
}
