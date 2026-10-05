#include "qaws_surface_coons.h"
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

typedef struct qaws_surface_coons_impl
{
	qaws_curve const* c0;
	qaws_curve const* c1;
	qaws_curve const* d0;
	qaws_curve const* d1;
	qaws_range range_c0;
	qaws_range range_c1;
	qaws_range range_d0;
	qaws_range range_d1;
	qaws_vec3 P00;  /* c0(0) = d0(0) */
	qaws_vec3 P10;  /* c0(1) = d1(0) */
	qaws_vec3 P01;  /* c1(0) = d0(1) */
	qaws_vec3 P11;  /* c1(1) = d1(1) */
} qaws_surface_coons_impl;

/* Coons patch: S(u,v) = Lc(u,v) + Ld(u,v) - B(u,v)
   Lc(u,v) = (1-v)*c0(u) + v*c1(u)
   Ld(u,v) = (1-u)*d0(v) + u*d1(v)
   B(u,v)  = (1-u)*(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11

   dS/du = (1-v)*c0'(u)*s_c0 + v*c1'(u)*s_c1
         - d0(v) + d1(v)
         - [-(1-v)*P00 + (1-v)*P10 - v*P01 + v*P11]

   dS/dv = -c0(u) + c1(u)
         + (1-u)*d0'(v)*s_d0 + u*d1'(v)*s_d1
         - [(1-u)*(-P00+P01) + u*(-P10+P11)]

   Second derivatives use central finite differences. */
static qaws_status coons_surface_eval_sampled(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_coons_impl const* impl =
		(qaws_surface_coons_impl const*)surface->impl;
	qaws_scalar one_minus_u = QAWS_ONE - u;
	qaws_scalar one_minus_v = QAWS_ONE - v;
	qaws_scalar s_c0 = impl->range_c0.max_value - impl->range_c0.min_value;
	qaws_scalar s_c1 = impl->range_c1.max_value - impl->range_c1.min_value;
	qaws_scalar s_d0 = impl->range_d0.max_value - impl->range_d0.min_value;
	qaws_scalar s_d1 = impl->range_d1.max_value - impl->range_d1.min_value;
	qaws_scalar t_c0 = impl->range_c0.min_value + u * s_c0;
	qaws_scalar t_c1 = impl->range_c1.min_value + u * s_c1;
	qaws_scalar t_d0 = impl->range_d0.min_value + v * s_d0;
	qaws_scalar t_d1 = impl->range_d1.min_value + v * s_d1;

	unsigned int c_flags = 0;
	unsigned int d_flags = 0;
	qaws_eval_result_3d rc0, rc1, rd0, rd1;
	qaws_status status;
	qaws_vec3 pos;

	/* Determine what we need from the curves */
	if (eval_flags & (QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DV))
		c_flags |= QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & (QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU))
		d_flags |= QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
		c_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
		d_flags |= QAWS_EVAL_FLAG_D1;

	/* For second derivatives and normal, we need position + D1 from all */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		c_flags |= QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
		d_flags |= QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
	}

	memset(&rc0, 0, sizeof(rc0));
	memset(&rc1, 0, sizeof(rc1));
	memset(&rd0, 0, sizeof(rd0));
	memset(&rd1, 0, sizeof(rd1));

	status = qaws_curve_evaluate_3d(impl->c0, t_c0, c_flags, &rc0);
	if (status != QAWS_STATUS_OK) return status;
	status = qaws_curve_evaluate_3d(impl->c1, t_c1, c_flags, &rc1);
	if (status != QAWS_STATUS_OK) return status;
	status = qaws_curve_evaluate_3d(impl->d0, t_d0, d_flags, &rd0);
	if (status != QAWS_STATUS_OK) return status;
	status = qaws_curve_evaluate_3d(impl->d1, t_d1, d_flags, &rd1);
	if (status != QAWS_STATUS_OK) return status;

	/* Position: (1-v)*c0 + v*c1 + (1-u)*d0 + u*d1
	   - [(1-u)*(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11] */
	pos.x = one_minus_v * rc0.position.x + v * rc1.position.x
		+ one_minus_u * rd0.position.x + u * rd1.position.x
		- (one_minus_u * one_minus_v * impl->P00.x
		   + u * one_minus_v * impl->P10.x
		   + one_minus_u * v * impl->P01.x
		   + u * v * impl->P11.x);
	pos.y = one_minus_v * rc0.position.y + v * rc1.position.y
		+ one_minus_u * rd0.position.y + u * rd1.position.y
		- (one_minus_u * one_minus_v * impl->P00.y
		   + u * one_minus_v * impl->P10.y
		   + one_minus_u * v * impl->P01.y
		   + u * v * impl->P11.y);
	pos.z = one_minus_v * rc0.position.z + v * rc1.position.z
		+ one_minus_u * rd0.position.z + u * rd1.position.z
		- (one_minus_u * one_minus_v * impl->P00.z
		   + u * one_minus_v * impl->P10.z
		   + one_minus_u * v * impl->P01.z
		   + u * v * impl->P11.z);

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = (1-v)*c0'(u)*s_c0 + v*c1'(u)*s_c1
	         - d0(v) + d1(v)
	         - [-(1-v)*P00 + (1-v)*P10 - v*P01 + v*P11] */
	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = one_minus_v * rc0.d1.x * s_c0 + v * rc1.d1.x * s_c1
			- rd0.position.x + rd1.position.x
			- (-(one_minus_v) * impl->P00.x + one_minus_v * impl->P10.x
			   - v * impl->P01.x + v * impl->P11.x);
		out_result->du.y = one_minus_v * rc0.d1.y * s_c0 + v * rc1.d1.y * s_c1
			- rd0.position.y + rd1.position.y
			- (-(one_minus_v) * impl->P00.y + one_minus_v * impl->P10.y
			   - v * impl->P01.y + v * impl->P11.y);
		out_result->du.z = one_minus_v * rc0.d1.z * s_c0 + v * rc1.d1.z * s_c1
			- rd0.position.z + rd1.position.z
			- (-(one_minus_v) * impl->P00.z + one_minus_v * impl->P10.z
			   - v * impl->P01.z + v * impl->P11.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = -c0(u) + c1(u)
	         + (1-u)*d0'(v)*s_d0 + u*d1'(v)*s_d1
	         - [(1-u)*(-P00+P01) + u*(-P10+P11)] */
	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = -rc0.position.x + rc1.position.x
			+ one_minus_u * rd0.d1.x * s_d0 + u * rd1.d1.x * s_d1
			- (one_minus_u * (-impl->P00.x + impl->P01.x)
			   + u * (-impl->P10.x + impl->P11.x));
		out_result->dv.y = -rc0.position.y + rc1.position.y
			+ one_minus_u * rd0.d1.y * s_d0 + u * rd1.d1.y * s_d1
			- (one_minus_u * (-impl->P00.y + impl->P01.y)
			   + u * (-impl->P10.y + impl->P11.y));
		out_result->dv.z = -rc0.position.z + rc1.position.z
			+ one_minus_u * rd0.d1.z * s_d0 + u * rd1.d1.z * s_d1
			- (one_minus_u * (-impl->P00.z + impl->P01.z)
			   + u * (-impl->P10.z + impl->P11.z));
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* Second derivatives via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_surface_eval_result r_lo, r_hi;

		/* duu via finite difference of du w.r.t. u */
		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			coons_surface_eval_sampled(surface, u_lo, v, QAWS_SURFACE_EVAL_DU, &r_lo);
			coons_surface_eval_sampled(surface, u_hi, v, QAWS_SURFACE_EVAL_DU, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duu.x = (r_hi.du.x - r_lo.du.x) * inv2h;
				out_result->duu.y = (r_hi.du.y - r_lo.du.y) * inv2h;
				out_result->duu.z = (r_hi.du.z - r_lo.du.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		/* dvv via finite difference of dv w.r.t. v */
		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		{
			qaws_scalar v_lo = v - h, v_hi = v + h;
			qaws_scalar hv;
			if (v_lo < 0) v_lo = 0;
			if (v_hi > 1) v_hi = 1;
			hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			coons_surface_eval_sampled(surface, u, v_lo, QAWS_SURFACE_EVAL_DV, &r_lo);
			coons_surface_eval_sampled(surface, u, v_hi, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
				out_result->dvv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->dvv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->dvv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}

		/* duv via finite difference of dv w.r.t. u */
		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			coons_surface_eval_sampled(surface, u_lo, v, QAWS_SURFACE_EVAL_DV, &r_lo);
			coons_surface_eval_sampled(surface, u_hi, v, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->duv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->duv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
		}

		/* Normal */
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
/*  S = (1-v) c0(u) + v c1(u) + (1-u) d0(v) + u d1(v) - B(u,v), where the     */
/*  bilinear corners come live from c0 and c1 so their derivatives flow into  */
/*  those curves. No own fields; children 0..3 are c0, c1, d0, d1.           */
/* -------------------------------------------------------------------------- */

#define COONS_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned char const g_coons_jet_a[QAWS_SURFACE_JET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_coons_jet_b[QAWS_SURFACE_JET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

/* Dual inputs of the Coons formula. */
typedef struct coons_inputs
{
	qaws_dual3 c0[4], c1[4], d0[4], d1[4];  /* jets in each curve's own parameter */
	qaws_dual3 P00, P10, P01, P11;
	qaws_dual1 u, v;
} coons_inputs;

static qaws_dual3 coons_scaled(qaws_dual3 x, qaws_dual1 k, qaws_scalar c)
{
	return qaws_dual3_mul_const(qaws_dual3_scale(x, k), c);
}

static void coons_jets(qaws_surface_coons_impl const* impl, coons_inputs const* in, qaws_dual3* S)
{
	qaws_scalar sc0 = impl->range_c0.max_value - impl->range_c0.min_value;
	qaws_scalar sc1 = impl->range_c1.max_value - impl->range_c1.min_value;
	qaws_scalar sd0 = impl->range_d0.max_value - impl->range_d0.min_value;
	qaws_scalar sd1 = impl->range_d1.max_value - impl->range_d1.min_value;
	qaws_dual1 one = qaws_dual1_const(QAWS_ONE), mone = qaws_dual1_const(-QAWS_ONE);
	qaws_dual1 omu = qaws_dual1_sub(one, in->u), omv = qaws_dual1_sub(one, in->v);
	unsigned int k;

	for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
	{
		unsigned int a = g_coons_jet_a[k], b = g_coons_jet_b[k], i;
		qaws_scalar pc0 = QAWS_ONE, pc1 = QAWS_ONE, pd0 = QAWS_ONE, pd1 = QAWS_ONE;
		qaws_dual3 acc = qaws_dual3_const(qaws_v3_zero());
		for (i = 0; i < a; i++) { pc0 *= sc0; pc1 *= sc1; }
		for (i = 0; i < b; i++) { pd0 *= sd0; pd1 *= sd1; }

		if (b == 0)
			acc = qaws_dual3_add(coons_scaled(in->c0[a], omv, pc0), coons_scaled(in->c1[a], in->v, pc1));
		else if (b == 1)
			acc = qaws_dual3_add(coons_scaled(in->c0[a], mone, pc0), coons_scaled(in->c1[a], one, pc1));

		if (a == 0)
			acc = qaws_dual3_add(acc, qaws_dual3_add(coons_scaled(in->d0[b], omu, pd0), coons_scaled(in->d1[b], in->u, pd1)));
		else if (a == 1)
			acc = qaws_dual3_add(acc, qaws_dual3_add(coons_scaled(in->d0[b], mone, pd0), coons_scaled(in->d1[b], one, pd1)));

		if (a <= 1 && b <= 1)
		{
			qaws_dual1 lu0 = a == 0 ? omu : mone, lu1 = a == 0 ? in->u : one;
			qaws_dual1 lv0 = b == 0 ? omv : mone, lv1 = b == 0 ? in->v : one;
			qaws_dual3 bl = qaws_dual3_scale(in->P00, qaws_dual1_mul(lu0, lv0));
			bl = qaws_dual3_add(bl, qaws_dual3_scale(in->P10, qaws_dual1_mul(lu1, lv0)));
			bl = qaws_dual3_add(bl, qaws_dual3_scale(in->P01, qaws_dual1_mul(lu0, lv1)));
			bl = qaws_dual3_add(bl, qaws_dual3_scale(in->P11, qaws_dual1_mul(lu1, lv1)));
			acc = qaws_dual3_sub(acc, bl);
		}
		S[k] = acc;
	}
}

static void coons_fill(qaws_dual3* dst, qaws_curve_jet_3d const* p, qaws_curve_jet_3d const* t, qaws_curve_jet_3d const* tt)
{
	unsigned int k;
	for (k = 0; k < 4; k++)
		dst[k] = qaws_dual3_make(p->d[k], t ? t->d[k] : qaws_v3_zero(), tt ? tt->d[k] : qaws_v3_zero());
}

/* Evaluates every child (with optional tangents) into dual inputs. */
static qaws_status coons_gather(
	qaws_diff_context const* ctx, qaws_surface_coons_impl const* impl,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	qaws_diff_views const* views, int with_tangent, int with_second, coons_inputs* in)
{
	qaws_curve const* curves[4];
	qaws_range ranges[4];
	qaws_scalar coord[4], coord_dot[4];
	qaws_dual3* slots[4];
	qaws_curve_jet_3d p, t, tt;
	unsigned int i;
	qaws_status st;

	curves[0] = impl->c0; curves[1] = impl->c1; curves[2] = impl->d0; curves[3] = impl->d1;
	ranges[0] = impl->range_c0; ranges[1] = impl->range_c1; ranges[2] = impl->range_d0; ranges[3] = impl->range_d1;
	coord[0] = coord[1] = u; coord[2] = coord[3] = v;
	coord_dot[0] = coord_dot[1] = u_dot; coord_dot[2] = coord_dot[3] = v_dot;
	slots[0] = in->c0; slots[1] = in->c1; slots[2] = in->d0; slots[3] = in->d1;

	for (i = 0; i < 4; i++)
	{
		qaws_scalar s = ranges[i].max_value - ranges[i].min_value;
		st = qaws_internal_curve_tangent_any(ctx, curves[i], ranges[i].min_value + coord[i] * s,
			with_tangent ? s * coord_dot[i] : QAWS_ZERO, 0xFu,
			with_tangent ? qaws_internal_child_views(views, i) : NULL, &p, &t, with_second ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		coons_fill(slots[i], &p, with_tangent ? &t : NULL, with_second ? &tt : NULL);
	}

	/* Corners at the ends of c0 and c1. */
	for (i = 0; i < 4; i++)
	{
		unsigned int ci = i / 2;
		qaws_scalar tc = (i % 2) ? ranges[ci].max_value : ranges[ci].min_value;
		qaws_dual3* dst = i == 0 ? &in->P00 : i == 1 ? &in->P10 : i == 2 ? &in->P01 : &in->P11;
		st = qaws_internal_curve_tangent_any(ctx, curves[ci], tc, QAWS_ZERO, QAWS_EVAL_FLAG_POSITION,
			with_tangent ? qaws_internal_child_views(views, ci) : NULL, &p, &t, with_second ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		*dst = qaws_dual3_make(p.d[0], with_tangent ? t.d[0] : qaws_v3_zero(), with_second ? tt.d[0] : qaws_v3_zero());
	}
	in->u = qaws_dual1_make(u, with_tangent ? u_dot : QAWS_ZERO, QAWS_ZERO);
	in->v = qaws_dual1_make(v, with_tangent ? v_dot : QAWS_ZERO, QAWS_ZERO);
	return QAWS_STATUS_OK;
}

static qaws_status coons_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_coons_impl const* impl = (qaws_surface_coons_impl const*)surface->impl;
	coons_inputs in;
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT];

	/* Boundary curves without analytic jets keep the sampled derivatives. */
	if (coons_gather(NULL, impl, u, v, 0, 0, NULL, 0, 0, &in) != QAWS_STATUS_OK)
		return coons_surface_eval_sampled(surface, u, v, eval_flags, out_result);
	coons_jets(impl, &in, S);

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION) { out_result->position = S[0].v; out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION; }
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL)) { out_result->du = S[1].v; out_result->valid_flags |= QAWS_SURFACE_EVAL_DU; }
	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL)) { out_result->dv = S[2].v; out_result->valid_flags |= QAWS_SURFACE_EVAL_DV; }
	if (eval_flags & QAWS_SURFACE_EVAL_DUU) { out_result->duu = S[3].v; out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU; }
	if (eval_flags & QAWS_SURFACE_EVAL_DUV) { out_result->duv = S[4].v; out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV; }
	if (eval_flags & QAWS_SURFACE_EVAL_DVV) { out_result->dvv = S[5].v; out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV; }
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}
	return QAWS_STATUS_OK;
}

static unsigned int coons_describe_fields(qaws_surface const* surface, qaws_field_desc* out, unsigned int capacity)
{
	(void)surface;
	(void)out;
	(void)capacity;
	return 0;
}

static qaws_status coons_primal_field(qaws_surface const* surface, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	(void)surface; (void)field; (void)out_data; (void)out_count; (void)out_components;
	return QAWS_STATUS_INVALID_ARGUMENT;
}

static unsigned int coons_children(qaws_surface const* surface, qaws_diff_child* out, unsigned int capacity)
{
	qaws_surface_coons_impl const* impl = (qaws_surface_coons_impl const*)surface->impl;
	if (capacity >= 4)
	{
		out[0].curve = impl->c0; out[0].surface = NULL;
		out[1].curve = impl->c1; out[1].surface = NULL;
		out[2].curve = impl->d0; out[2].surface = NULL;
		out[3].curve = impl->d1; out[3].surface = NULL;
	}
	return 4;
}

static qaws_status coons_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_dot, qaws_scalar v_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_surface_jet* primal, qaws_surface_jet* tangent, qaws_surface_jet* tangent2)
{
	qaws_surface_coons_impl const* impl = (qaws_surface_coons_impl const*)surface->impl;
	coons_inputs in;
	qaws_dual3 S[QAWS_SURFACE_JET_COUNT];
	unsigned int k;
	qaws_status st = coons_gather(ctx, impl, u, v, u_dot, v_dot, views, 1, tangent2 != NULL, &in);
	if (st != QAWS_STATUS_OK)
		return st;
	coons_jets(impl, &in, S);

	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));
	for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
	{
		if (!(channels & (1u << k)))
			continue;
		primal->d[k] = S[k].v;
		tangent->d[k] = S[k].t;
		if (tangent2)
			tangent2->d[k] = S[k].tt;
	}
	primal->channels = tangent->channels = channels;
	if (tangent2)
		tangent2->channels = channels;
	return QAWS_STATUS_OK;
}

/* S is linear in the curve jets and corners; their adjoints are formed by
   seeding each input component, then pulled back through the curves. */
static qaws_status coons_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* ybar, qaws_diff_views* views,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	qaws_surface_coons_impl const* impl = (qaws_surface_coons_impl const*)surface->impl;
	coons_inputs base, seeded;
	qaws_dual3* seed_slots[20];
	qaws_vec3 bars[20];
	qaws_scalar coord_bar[2] = { 0, 0 };
	qaws_curve const* curves[4];
	qaws_range ranges[4];
	unsigned int i, c, k;
	qaws_status st = coons_gather(ctx, impl, u, v, 0, 0, NULL, 0, 0, &base);
	if (st != QAWS_STATUS_OK)
		return st;

	for (i = 0; i < 22; i++)
	{
		for (c = 0; c < (i < 20 ? 3u : 1u); c++)
		{
			qaws_dual3 S[QAWS_SURFACE_JET_COUNT];
			qaws_scalar acc = QAWS_ZERO;
			seeded = base;
			for (k = 0; k < 4; k++)
			{
				seed_slots[k] = &seeded.c0[k]; seed_slots[4 + k] = &seeded.c1[k];
				seed_slots[8 + k] = &seeded.d0[k]; seed_slots[12 + k] = &seeded.d1[k];
			}
			seed_slots[16] = &seeded.P00; seed_slots[17] = &seeded.P10; seed_slots[18] = &seeded.P01; seed_slots[19] = &seeded.P11;
			if (i < 20)
				seed_slots[i]->t = qaws_v3(c == 0 ? QAWS_ONE : QAWS_ZERO, c == 1 ? QAWS_ONE : QAWS_ZERO, c == 2 ? QAWS_ONE : QAWS_ZERO);
			else if (i == 20)
				seeded.u.t = QAWS_ONE;
			else
				seeded.v.t = QAWS_ONE;
			coons_jets(impl, &seeded, S);
			for (k = 0; k < QAWS_SURFACE_JET_COUNT; k++)
				if (channels & (1u << k))
					acc += qaws_v3_dot(ybar->d[k], S[k].t);
			if (i >= 20)
				coord_bar[i - 20] = acc;
			else if (c == 0)
				bars[i].x = acc;
			else if (c == 1)
				bars[i].y = acc;
			else
				bars[i].z = acc;
		}
	}

	curves[0] = impl->c0; curves[1] = impl->c1; curves[2] = impl->d0; curves[3] = impl->d1;
	ranges[0] = impl->range_c0; ranges[1] = impl->range_c1; ranges[2] = impl->range_d0; ranges[3] = impl->range_d1;
	for (i = 0; i < 4; i++)
	{
		qaws_curve_jet_3d jbar;
		qaws_scalar s = ranges[i].max_value - ranges[i].min_value;
		qaws_scalar coord = i < 2 ? u : v;
		qaws_scalar tbar = QAWS_ZERO;
		int want = i < 2 ? u_adjoint != NULL : v_adjoint != NULL;
		qaws_diff_views* cv = (qaws_diff_views*)qaws_internal_child_views(views, i);
		for (k = 0; k < 4; k++)
			jbar.d[k] = bars[4 * i + k];
		jbar.channels = 0xFu;
		st = qaws_internal_curve_adjoint_any(ctx, curves[i], ranges[i].min_value + coord * s, 0xFu, &jbar, cv,
			want ? &tbar : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		coord_bar[i < 2 ? 0 : 1] += s * tbar;
		if (i < 2)
		{
			/* corner adjoints: start and end of c0 (i = 0) and c1 (i = 1) */
			unsigned int e;
			for (e = 0; e < 2; e++)
			{
				qaws_curve_jet_3d cbar;
				memset(&cbar, 0, sizeof(cbar));
				cbar.d[0] = bars[16 + 2 * i + e];
				cbar.channels = QAWS_EVAL_FLAG_POSITION;
				st = qaws_internal_curve_adjoint_any(ctx, curves[i], e ? ranges[i].max_value : ranges[i].min_value,
					QAWS_EVAL_FLAG_POSITION, &cbar, cv, NULL);
				if (st != QAWS_STATUS_OK)
					return st;
			}
		}
	}
	if (u_adjoint)
		*u_adjoint += coord_bar[0];
	if (v_adjoint)
		*v_adjoint += coord_bar[1];
	return QAWS_STATUS_OK;
}

static qaws_surface_diff_vtable const coons_surface_diff_vtable = {
	COONS_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	coons_describe_fields,
	coons_primal_field,
	NULL,
	NULL,
	coons_tangent,
	coons_adjoint,
	coons_children
};

static void coons_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int coons_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const coons_surface_vtable = {
	coons_surface_eval,
	coons_surface_destroy,
	coons_surface_is_rational,
	&coons_surface_diff_vtable
};

qaws_status qaws_surface_create_coons(
	qaws_surface_coons_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_coons_impl* impl;
	qaws_range u_range, v_range;
	qaws_eval_result_3d corner_r;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->c0 || !desc->c1 || !desc->d0 || !desc->d1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(desc->c0) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->c1) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->d0) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (qaws_curve_get_dimension(desc->d1) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_COONS,
		0, 0, u_range, v_range,
		&coons_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_coons_impl*)malloc(sizeof(qaws_surface_coons_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->c0 = desc->c0;
	impl->c1 = desc->c1;
	impl->d0 = desc->d0;
	impl->d1 = desc->d1;
	impl->range_c0 = qaws_curve_get_parameter_range(desc->c0);
	impl->range_c1 = qaws_curve_get_parameter_range(desc->c1);
	impl->range_d0 = qaws_curve_get_parameter_range(desc->d0);
	impl->range_d1 = qaws_curve_get_parameter_range(desc->d1);

	/* Cache corner points: P00=c0(start), P10=c0(end), P01=c1(start), P11=c1(end) */
	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c0, impl->range_c0.min_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P00 = corner_r.position;

	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c0, impl->range_c0.max_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P10 = corner_r.position;

	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c1, impl->range_c1.min_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P01 = corner_r.position;

	memset(&corner_r, 0, sizeof(corner_r));
	qaws_curve_evaluate_3d(desc->c1, impl->range_c1.max_value,
		QAWS_EVAL_FLAG_POSITION, &corner_r);
	impl->P11 = corner_r.position;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
