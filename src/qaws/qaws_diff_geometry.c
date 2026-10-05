#include "qaws_diff_geometry.h"
#include "core/qaws_dual_core.h"
#include <string.h>

/*
 * All quantities are evaluated once with second-order dual numbers. The
 * adjoint of these small fixed-size maps is formed exactly by seeding one
 * jet component at a time (forward-mode Jacobian) and contracting with the
 * quantity adjoint; no finite differences are involved.
 */

static qaws_dual1 dual1_zero(void)
{
	return qaws_dual1_const(QAWS_ZERO);
}

static qaws_dual3 dual3_zero(void)
{
	return qaws_dual3_const(qaws_v3_zero());
}

static qaws_dual1 dual1_scale(qaws_dual1 a, qaws_scalar s)
{
	return qaws_dual1_make(a.v * s, a.t * s, a.tt * s);
}

/* Keep the value, drop the derivatives (quantity not differentiable here). */
static qaws_dual1 dual1_freeze(qaws_dual1 a)
{
	return qaws_dual1_const(a.v);
}

/* ================================================================== */
/*  Curves                                                            */
/* ================================================================== */

typedef struct curve_geometry_dual
{
	qaws_dual3 T, N, B;
	qaws_dual1 speed, kappa, tau;
	qaws_diff_validity validity;
} curve_geometry_dual;

static void curve_geometry_dual_3d(qaws_dual3 d1, qaws_dual3 d2, qaws_dual3 d3, curve_geometry_dual* g)
{
	qaws_dual3 b;
	qaws_dual1 bl, speed3;

	g->validity = QAWS_DIFF_VALID;
	g->speed = qaws_dual3_length(d1);
	if (!(g->speed.v > QAWS_EPSILON))
	{
		memset(g, 0, sizeof(*g));
		g->validity = QAWS_DIFF_INVALID;
		return;
	}
	g->T = qaws_dual3_div(d1, g->speed);

	b = qaws_dual3_cross(d1, d2);
	bl = qaws_dual3_length(b);
	speed3 = qaws_dual1_mul(qaws_dual1_mul(g->speed, g->speed), g->speed);
	g->kappa = qaws_dual1_div(bl, speed3);

	if (!(bl.v > QAWS_EPSILON * g->speed.v * g->speed.v))
	{
		/* Straight point: |C' x C''| is not differentiable at zero and the
		   frame and torsion are undefined. */
		g->kappa = dual1_freeze(g->kappa);
		g->B = dual3_zero();
		g->N = dual3_zero();
		g->tau = dual1_zero();
		g->validity = QAWS_DIFF_ILL_CONDITIONED;
		return;
	}

	g->B = qaws_dual3_div(b, bl);
	g->N = qaws_dual3_cross(g->B, g->T);
	g->tau = qaws_dual1_div(qaws_dual3_dot(b, d3), qaws_dual1_mul(bl, bl));
}

/* 2D: inputs embedded with z = 0, normal is +90 degree rotation. */
static void curve_geometry_dual_2d(qaws_dual3 d1, qaws_dual3 d2, curve_geometry_dual* g)
{
	qaws_dual1 speed3, num;

	memset(g, 0, sizeof(*g));
	g->speed = qaws_dual3_length(d1);
	if (!(g->speed.v > QAWS_EPSILON))
	{
		g->validity = QAWS_DIFF_INVALID;
		return;
	}
	g->validity = QAWS_DIFF_VALID;
	g->T = qaws_dual3_div(d1, g->speed);
	g->N = qaws_dual3_make(qaws_v3(-g->T.v.y, g->T.v.x, 0),
	                       qaws_v3(-g->T.t.y, g->T.t.x, 0),
	                       qaws_v3(-g->T.tt.y, g->T.tt.x, 0));
	num = qaws_dual1_make(qaws_v3_cross(d1.v, d2.v).z,
		qaws_v3_cross(d1.t, d2.v).z + qaws_v3_cross(d1.v, d2.t).z,
		qaws_v3_cross(d1.tt, d2.v).z + QAWS_LITERAL(2.0) * qaws_v3_cross(d1.t, d2.t).z + qaws_v3_cross(d1.v, d2.tt).z);
	speed3 = qaws_dual1_mul(qaws_dual1_mul(g->speed, g->speed), g->speed);
	g->kappa = qaws_dual1_div(num, speed3);
}

static qaws_vec3 vec3_of2(qaws_vec2 a)
{
	return qaws_v3(a.x, a.y, QAWS_ZERO);
}

static qaws_vec2 vec2_of3(qaws_vec3 a)
{
	qaws_vec2 r;
	r.x = a.x;
	r.y = a.y;
	return r;
}

static void store_curve_3d(curve_geometry_dual const* g, int which, qaws_curve_geometry_3d* out)
{
	if (!out)
		return;
	if (which == 0)
	{
		out->tangent = g->T.v; out->normal = g->N.v; out->binormal = g->B.v;
		out->speed = g->speed.v; out->curvature = g->kappa.v; out->torsion = g->tau.v;
	}
	else if (which == 1)
	{
		out->tangent = g->T.t; out->normal = g->N.t; out->binormal = g->B.t;
		out->speed = g->speed.t; out->curvature = g->kappa.t; out->torsion = g->tau.t;
	}
	else
	{
		out->tangent = g->T.tt; out->normal = g->N.tt; out->binormal = g->B.tt;
		out->speed = g->speed.tt; out->curvature = g->kappa.tt; out->torsion = g->tau.tt;
	}
}

static void store_curve_2d(curve_geometry_dual const* g, int which, qaws_curve_geometry_2d* out)
{
	if (!out)
		return;
	if (which == 0)
	{
		out->tangent = vec2_of3(g->T.v); out->normal = vec2_of3(g->N.v);
		out->speed = g->speed.v; out->curvature = g->kappa.v;
	}
	else if (which == 1)
	{
		out->tangent = vec2_of3(g->T.t); out->normal = vec2_of3(g->N.t);
		out->speed = g->speed.t; out->curvature = g->kappa.t;
	}
	else
	{
		out->tangent = vec2_of3(g->T.tt); out->normal = vec2_of3(g->N.tt);
		out->speed = g->speed.tt; out->curvature = g->kappa.tt;
	}
}

qaws_status qaws_curve_geometry_eval_3d(
	qaws_curve_jet_3d const* primal,
	qaws_curve_jet_3d const* tangent,
	qaws_curve_jet_3d const* tangent2,
	qaws_curve_geometry_3d* out_value,
	qaws_curve_geometry_3d* out_tangent,
	qaws_curve_geometry_3d* out_tangent2,
	qaws_diff_validity* out_validity)
{
	qaws_dual3 d[3];
	curve_geometry_dual g;
	unsigned int k;

	if (!primal || !out_value || (out_tangent && !tangent) || (out_tangent2 && (!tangent || !tangent2)))
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (k = 0; k < 3; k++)
		d[k] = qaws_dual3_make(primal->d[k + 1],
			tangent ? tangent->d[k + 1] : qaws_v3_zero(),
			tangent2 ? tangent2->d[k + 1] : qaws_v3_zero());
	curve_geometry_dual_3d(d[0], d[1], d[2], &g);
	store_curve_3d(&g, 0, out_value);
	store_curve_3d(&g, 1, out_tangent);
	store_curve_3d(&g, 2, out_tangent2);
	if (out_validity)
		*out_validity = g.validity;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_geometry_eval_2d(
	qaws_curve_jet_2d const* primal,
	qaws_curve_jet_2d const* tangent,
	qaws_curve_jet_2d const* tangent2,
	qaws_curve_geometry_2d* out_value,
	qaws_curve_geometry_2d* out_tangent,
	qaws_curve_geometry_2d* out_tangent2,
	qaws_diff_validity* out_validity)
{
	qaws_dual3 d[2];
	curve_geometry_dual g;
	unsigned int k;

	if (!primal || !out_value || (out_tangent && !tangent) || (out_tangent2 && (!tangent || !tangent2)))
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (k = 0; k < 2; k++)
		d[k] = qaws_dual3_make(vec3_of2(primal->d[k + 1]),
			tangent ? vec3_of2(tangent->d[k + 1]) : qaws_v3_zero(),
			tangent2 ? vec3_of2(tangent2->d[k + 1]) : qaws_v3_zero());
	curve_geometry_dual_2d(d[0], d[1], &g);
	store_curve_2d(&g, 0, out_value);
	store_curve_2d(&g, 1, out_tangent);
	store_curve_2d(&g, 2, out_tangent2);
	if (out_validity)
		*out_validity = g.validity;
	return QAWS_STATUS_OK;
}

static double curve_adjoint_dot_3d(qaws_curve_geometry_3d const* a, curve_geometry_dual const* g)
{
	return (double)qaws_v3_dot(a->tangent, g->T.t) + qaws_v3_dot(a->normal, g->N.t) +
	       qaws_v3_dot(a->binormal, g->B.t) +
	       (double)a->speed * g->speed.t + (double)a->curvature * g->kappa.t + (double)a->torsion * g->tau.t;
}

static double curve_adjoint_dot_2d(qaws_curve_geometry_2d const* a, curve_geometry_dual const* g)
{
	return (double)a->tangent.x * g->T.t.x + (double)a->tangent.y * g->T.t.y +
	       (double)a->normal.x * g->N.t.x + (double)a->normal.y * g->N.t.y +
	       (double)a->speed * g->speed.t + (double)a->curvature * g->kappa.t;
}

static void set_component(qaws_vec3* v, unsigned int c, qaws_scalar x)
{
	if (c == 0) v->x = x;
	else if (c == 1) v->y = x;
	else v->z = x;
}

qaws_status qaws_curve_geometry_adjoint_3d(
	qaws_curve_jet_3d const* primal,
	qaws_curve_geometry_3d const* adjoint,
	qaws_curve_jet_3d* inout_jet_adjoint,
	qaws_diff_validity* out_validity)
{
	unsigned int k, c, q;
	qaws_diff_validity validity = QAWS_DIFF_VALID;

	if (!primal || !adjoint || !inout_jet_adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (k = 0; k < 3; k++)
	{
		for (c = 0; c < 3; c++)
		{
			qaws_dual3 d[3];
			curve_geometry_dual g;
			qaws_vec3 seed = qaws_v3_zero();
			set_component(&seed, c, QAWS_ONE);
			for (q = 0; q < 3; q++)
				d[q] = qaws_dual3_make(primal->d[q + 1], q == k ? seed : qaws_v3_zero(), qaws_v3_zero());
			curve_geometry_dual_3d(d[0], d[1], d[2], &g);
			validity = g.validity;
			if (c == 0) inout_jet_adjoint->d[k + 1].x += (qaws_scalar)curve_adjoint_dot_3d(adjoint, &g);
			else if (c == 1) inout_jet_adjoint->d[k + 1].y += (qaws_scalar)curve_adjoint_dot_3d(adjoint, &g);
			else inout_jet_adjoint->d[k + 1].z += (qaws_scalar)curve_adjoint_dot_3d(adjoint, &g);
		}
	}
	inout_jet_adjoint->channels |= QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;
	if (out_validity)
		*out_validity = validity;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_geometry_adjoint_2d(
	qaws_curve_jet_2d const* primal,
	qaws_curve_geometry_2d const* adjoint,
	qaws_curve_jet_2d* inout_jet_adjoint,
	qaws_diff_validity* out_validity)
{
	unsigned int k, c, q;
	qaws_diff_validity validity = QAWS_DIFF_VALID;

	if (!primal || !adjoint || !inout_jet_adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (k = 0; k < 2; k++)
	{
		for (c = 0; c < 2; c++)
		{
			qaws_dual3 d[2];
			curve_geometry_dual g;
			qaws_vec3 seed = qaws_v3_zero();
			double dot;
			set_component(&seed, c, QAWS_ONE);
			for (q = 0; q < 2; q++)
				d[q] = qaws_dual3_make(vec3_of2(primal->d[q + 1]), q == k ? seed : qaws_v3_zero(), qaws_v3_zero());
			curve_geometry_dual_2d(d[0], d[1], &g);
			validity = g.validity;
			dot = curve_adjoint_dot_2d(adjoint, &g);
			if (c == 0) inout_jet_adjoint->d[k + 1].x += (qaws_scalar)dot;
			else inout_jet_adjoint->d[k + 1].y += (qaws_scalar)dot;
		}
	}
	inout_jet_adjoint->channels |= QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;
	if (out_validity)
		*out_validity = validity;
	return QAWS_STATUS_OK;
}

/* ================================================================== */
/*  Surfaces                                                          */
/* ================================================================== */

typedef struct surface_geometry_dual
{
	qaws_dual3 normal;
	qaws_dual1 E, F, G, L, M, N;
	qaws_dual1 K, H, k1, k2;
	qaws_diff_validity validity;
} surface_geometry_dual;

/* Jet slots used: 1 = Su, 2 = Sv, 3 = Suu, 4 = Suv, 5 = Svv. */
static void surface_geometry_dual_eval(qaws_dual3 const* d, surface_geometry_dual* g)
{
	qaws_dual3 su = d[1], sv = d[2];
	qaws_dual3 n_raw = qaws_dual3_cross(su, sv);
	qaws_dual1 nl = qaws_dual3_length(n_raw);
	qaws_dual1 det, two_det, disc, num_h;
	qaws_scalar scale;

	memset(g, 0, sizeof(*g));
	scale = QAWS_SQRT(qaws_v3_dot(su.v, su.v) * qaws_v3_dot(sv.v, sv.v));
	if (!(nl.v > QAWS_EPSILON * scale) || !(scale > QAWS_ZERO))
	{
		g->validity = QAWS_DIFF_INVALID;
		return;
	}
	g->validity = QAWS_DIFF_VALID;
	g->normal = qaws_dual3_div(n_raw, nl);

	g->E = qaws_dual3_dot(su, su);
	g->F = qaws_dual3_dot(su, sv);
	g->G = qaws_dual3_dot(sv, sv);
	g->L = qaws_dual3_dot(d[3], g->normal);
	g->M = qaws_dual3_dot(d[4], g->normal);
	g->N = qaws_dual3_dot(d[5], g->normal);

	det = qaws_dual1_sub(qaws_dual1_mul(g->E, g->G), qaws_dual1_mul(g->F, g->F));
	g->K = qaws_dual1_div(qaws_dual1_sub(qaws_dual1_mul(g->L, g->N), qaws_dual1_mul(g->M, g->M)), det);
	num_h = qaws_dual1_add(qaws_dual1_sub(qaws_dual1_mul(g->E, g->N),
		dual1_scale(qaws_dual1_mul(g->F, g->M), QAWS_LITERAL(2.0))), qaws_dual1_mul(g->G, g->L));
	two_det = dual1_scale(det, QAWS_LITERAL(2.0));
	g->H = qaws_dual1_div(num_h, two_det);

	disc = qaws_dual1_sub(qaws_dual1_mul(g->H, g->H), g->K);
	if (!(disc.v > QAWS_LITERAL(100.0) * QAWS_EPSILON * (g->H.v * g->H.v + QAWS_FABS(g->K.v) + QAWS_EPSILON)))
	{
		/* Umbilic: both principal curvatures equal H, their square-root
		   branch is not differentiable. */
		g->k1 = dual1_freeze(g->H);
		g->k2 = dual1_freeze(g->H);
		g->validity = QAWS_DIFF_ILL_CONDITIONED;
		return;
	}
	{
		qaws_dual1 s = qaws_dual1_sqrt(disc);
		g->k1 = qaws_dual1_add(g->H, s);
		g->k2 = qaws_dual1_sub(g->H, s);
	}
}

static void store_surface(surface_geometry_dual const* g, int which, qaws_surface_geometry* out)
{
#define SG_PICK(x) (which == 0 ? (x).v : which == 1 ? (x).t : (x).tt)
	if (!out)
		return;
	out->normal = SG_PICK(g->normal);
	out->E = SG_PICK(g->E);
	out->F = SG_PICK(g->F);
	out->G = SG_PICK(g->G);
	out->L = SG_PICK(g->L);
	out->M = SG_PICK(g->M);
	out->N = SG_PICK(g->N);
	out->gaussian = SG_PICK(g->K);
	out->mean = SG_PICK(g->H);
	out->kappa1 = SG_PICK(g->k1);
	out->kappa2 = SG_PICK(g->k2);
#undef SG_PICK
}

qaws_status qaws_surface_geometry_eval(
	qaws_surface_jet const* primal,
	qaws_surface_jet const* tangent,
	qaws_surface_jet const* tangent2,
	qaws_surface_geometry* out_value,
	qaws_surface_geometry* out_tangent,
	qaws_surface_geometry* out_tangent2,
	qaws_diff_validity* out_validity)
{
	qaws_dual3 d[6];
	surface_geometry_dual g;
	unsigned int k;

	if (!primal || !out_value || (out_tangent && !tangent) || (out_tangent2 && (!tangent || !tangent2)))
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (k = 0; k < 6; k++)
		d[k] = qaws_dual3_make(primal->d[k],
			tangent ? tangent->d[k] : qaws_v3_zero(),
			tangent2 ? tangent2->d[k] : qaws_v3_zero());
	surface_geometry_dual_eval(d, &g);
	store_surface(&g, 0, out_value);
	store_surface(&g, 1, out_tangent);
	store_surface(&g, 2, out_tangent2);
	if (out_validity)
		*out_validity = g.validity;
	return QAWS_STATUS_OK;
}

static double surface_adjoint_dot(qaws_surface_geometry const* a, surface_geometry_dual const* g)
{
	return (double)qaws_v3_dot(a->normal, g->normal.t) +
	       (double)a->E * g->E.t + (double)a->F * g->F.t + (double)a->G * g->G.t +
	       (double)a->L * g->L.t + (double)a->M * g->M.t + (double)a->N * g->N.t +
	       (double)a->gaussian * g->K.t + (double)a->mean * g->H.t +
	       (double)a->kappa1 * g->k1.t + (double)a->kappa2 * g->k2.t;
}

qaws_status qaws_surface_geometry_adjoint(
	qaws_surface_jet const* primal,
	qaws_surface_geometry const* adjoint,
	qaws_surface_jet* inout_jet_adjoint,
	qaws_diff_validity* out_validity)
{
	unsigned int k, c, q;
	qaws_diff_validity validity = QAWS_DIFF_VALID;

	if (!primal || !adjoint || !inout_jet_adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (k = 1; k < 6; k++)
	{
		for (c = 0; c < 3; c++)
		{
			qaws_dual3 d[6];
			surface_geometry_dual g;
			qaws_vec3 seed = qaws_v3_zero();
			double dot;
			set_component(&seed, c, QAWS_ONE);
			for (q = 0; q < 6; q++)
				d[q] = qaws_dual3_make(primal->d[q], q == k ? seed : qaws_v3_zero(), qaws_v3_zero());
			surface_geometry_dual_eval(d, &g);
			validity = g.validity;
			dot = surface_adjoint_dot(adjoint, &g);
			if (c == 0) inout_jet_adjoint->d[k].x += (qaws_scalar)dot;
			else if (c == 1) inout_jet_adjoint->d[k].y += (qaws_scalar)dot;
			else inout_jet_adjoint->d[k].z += (qaws_scalar)dot;
		}
	}
	inout_jet_adjoint->channels |= QAWS_SJET_U | QAWS_SJET_V | QAWS_SJET_UU | QAWS_SJET_UV | QAWS_SJET_VV;
	if (out_validity)
		*out_validity = validity;
	return QAWS_STATUS_OK;
}
