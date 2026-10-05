#include "qaws_diff_ops.h"
#include "qaws_inspect.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include <string.h>
#include <math.h>

#define OPS_GAP_SAMPLES 128
#define OPS_GAP_GRID 24
#define OPS_NEWTON_STEPS 6

/* ================================================================== */
/*  Report                                                            */
/* ================================================================== */

static void report_implicit(
	qaws_diff_context const* ctx,
	qaws_diff_validity validity,
	unsigned int frozen,
	qaws_scalar condition,
	qaws_scalar residual,
	qaws_scalar gap)
{
	qaws_diff_report* r;
	qaws_internal_diff_report_note(ctx,
		frozen ? QAWS_DIFF_ACTIVE_SET : QAWS_DIFF_SMOOTH, validity, frozen, 0);
	if (!ctx || !ctx->report)
		return;
	r = ctx->report;
	if (condition > r->condition_number) r->condition_number = condition;
	if (residual > r->residual) r->residual = residual;
	if (gap < r->branch_gap) r->branch_gap = gap;
}

/* Validity from conditioning, competing solutions and active bounds. */
static qaws_diff_validity implicit_validity(int singular, qaws_scalar gap, qaws_scalar distance, int frozen)
{
	qaws_scalar scale = distance > QAWS_LITERAL(1e-6) ? distance : QAWS_LITERAL(1e-6);
	if (singular)
		return QAWS_DIFF_ILL_CONDITIONED;
	if (gap < QAWS_LITERAL(0.005) * scale)
		return QAWS_DIFF_AMBIGUOUS;
	if (frozen)
		return QAWS_DIFF_VALID_LOCALLY;
	return QAWS_DIFF_VALID;
}

/* ================================================================== */
/*  Curves                                                            */
/* ================================================================== */

typedef struct curve_solution
{
	qaws_scalar t;
	qaws_curve_jet_3d jet;   /* P, D1, D2 at t */
	qaws_vec3 r;             /* P - q */
	qaws_scalar g, gt, distance;
	int at_bound, singular;
	qaws_scalar gap;
} curve_solution;

static qaws_status curve_jet_at(qaws_curve const* curve, qaws_scalar t, qaws_curve_jet_3d* jet)
{
	qaws_curve_jet_3d unused;
	return qaws_internal_curve_tangent_any(NULL, curve, t, QAWS_ZERO,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, NULL, jet, &unused, NULL);
}

/* Newton refinement on g(t) = C' . (C - q) starting at t0, so the implicit
   derivative is taken at an accurate solution. */
static qaws_status curve_refine(qaws_curve const* curve, qaws_vec3 q, qaws_scalar t0, curve_solution* s)
{
	qaws_range range = curve->parameter_range;
	unsigned int it;
	qaws_status st;

	s->t = t0;
	for (it = 0; it <= OPS_NEWTON_STEPS; it++)
	{
		st = curve_jet_at(curve, s->t, &s->jet);
		if (st != QAWS_STATUS_OK)
			return st;
		s->r = qaws_v3_sub(s->jet.d[0], q);
		s->g = qaws_v3_dot(s->jet.d[1], s->r);
		s->gt = qaws_v3_dot(s->jet.d[2], s->r) + qaws_v3_dot(s->jet.d[1], s->jet.d[1]);
		if (it == OPS_NEWTON_STEPS || !(s->gt > QAWS_ZERO))
			break;
		s->t -= s->g / s->gt;
		if (s->t < range.min_value) s->t = range.min_value;
		if (s->t > range.max_value) s->t = range.max_value;
	}
	s->distance = QAWS_SQRT(qaws_v3_dot(s->r, s->r));
	return QAWS_STATUS_OK;
}

static qaws_status curve_solve(qaws_curve const* curve, qaws_vec3 q, curve_solution* s)
{
	qaws_range range = curve->parameter_range;
	qaws_scalar eps = (range.max_value - range.min_value) * QAWS_LITERAL(1e-7);
	qaws_scalar step = (range.max_value - range.min_value) / (qaws_scalar)(OPS_GAP_SAMPLES - 1);
	qaws_scalar d[OPS_GAP_SAMPLES], t0;
	unsigned int i, best = 0;
	qaws_status st;

	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		qaws_vec2 q2;
		q2.x = q.x;
		q2.y = q.y;
		q.z = QAWS_ZERO;
		st = qaws_curve_find_closest_parameter_2d(curve, q2, &t0);
	}
	else
		st = qaws_curve_find_closest_parameter_3d(curve, q, &t0);
	if (st != QAWS_STATUS_OK)
		return st;
	st = curve_refine(curve, q, t0, s);
	if (st != QAWS_STATUS_OK)
		return st;

	/* Global check on a sampled distance profile: when a sampled basin is
	   closer than the solver's answer, refine from it instead. */
	for (i = 0; i < OPS_GAP_SAMPLES; i++)
	{
		qaws_curve_jet_3d j;
		qaws_vec3 r;
		if (curve_jet_at(curve, range.min_value + step * (qaws_scalar)i, &j) != QAWS_STATUS_OK)
			return QAWS_STATUS_INTERNAL_ERROR;
		r = qaws_v3_sub(j.d[0], q);
		d[i] = QAWS_SQRT(qaws_v3_dot(r, r));
		if (d[i] < d[best])
			best = i;
	}
	if (d[best] < s->distance)
	{
		curve_solution alt;
		st = curve_refine(curve, q, range.min_value + step * (qaws_scalar)best, &alt);
		if (st != QAWS_STATUS_OK)
			return st;
		if (alt.distance < s->distance)
			*s = alt;
	}

	s->at_bound = (s->t <= range.min_value + eps) || (s->t >= range.max_value - eps);
	s->singular = !s->at_bound && !(s->gt > QAWS_LITERAL(1e-8) * qaws_v3_dot(s->jet.d[1], s->jet.d[1]));

	/* Best competing local minimum: separated from the solution by a ridge. */
	s->gap = QAWS_DIFF_UNBOUNDED;
	for (i = 0; i < OPS_GAP_SAMPLES; i++)
	{
		qaws_scalar ti = range.min_value + step * (qaws_scalar)i;
		int local_min = (i == 0 || d[i] <= d[i - 1]) && (i + 1 == OPS_GAP_SAMPLES || d[i] <= d[i + 1]);
		unsigned int home = (unsigned int)((s->t - range.min_value) / step + QAWS_LITERAL(0.5));
		unsigned int lo, hi, k;
		qaws_scalar ridge = QAWS_ZERO;
		if (!local_min || QAWS_FABS(ti - s->t) <= QAWS_LITERAL(2.5) * step)
			continue;
		if (home >= OPS_GAP_SAMPLES)
			home = OPS_GAP_SAMPLES - 1;
		lo = i < home ? i : home;
		hi = i < home ? home : i;
		for (k = lo; k <= hi; k++)
			if (d[k] > ridge)
				ridge = d[k];
		if (!(ridge > (d[i] > s->distance ? d[i] : s->distance) * (QAWS_ONE + QAWS_LITERAL(1e-6))))
			continue;
		if (d[i] - s->distance < s->gap)
			s->gap = d[i] - s->distance;
	}
	if (s->gap < QAWS_ZERO)
		s->gap = QAWS_ZERO;
	return QAWS_STATUS_OK;
}

static void curve_report(qaws_diff_context const* ctx, curve_solution const* s)
{
	qaws_scalar speed2 = qaws_v3_dot(s->jet.d[1], s->jet.d[1]);
	qaws_scalar cond = (s->gt > QAWS_ZERO) ? speed2 / s->gt : QAWS_DIFF_UNBOUNDED;
	qaws_scalar res = QAWS_FABS(s->g) / (QAWS_SQRT(speed2) * s->distance + QAWS_LITERAL(1e-30));
	report_implicit(ctx, implicit_validity(s->singular, s->gap, s->distance, s->at_bound),
		s->at_bound ? (unsigned int)QAWS_FREEZE_ACTIVE_SET : 0u,
		s->at_bound ? QAWS_ONE : cond, s->at_bound ? QAWS_ZERO : res, s->gap);
}

qaws_status qaws_curve_closest_point_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 query,
	qaws_vec3 const* query_tangent,
	qaws_diff_views const* param_tangent,
	qaws_curve_closest_point* out_value,
	qaws_curve_closest_point* out_tangent)
{
	curve_solution s;
	qaws_curve_jet_3d p, tg;
	qaws_vec3 qd = query_tangent ? *query_tangent : qaws_v3_zero();
	qaws_scalar t_dot = QAWS_ZERO;
	qaws_vec3 p_dot;
	qaws_status st;

	if (!curve || !out_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		query.z = QAWS_ZERO;
		qd.z = QAWS_ZERO;
	}
	st = curve_solve(curve, query, &s);
	if (st != QAWS_STATUS_OK)
		return st;

	out_value->t = s.t;
	out_value->position = s.jet.d[0];
	out_value->distance = s.distance;
	curve_report(ctx, &s);
	if (!out_tangent)
		return QAWS_STATUS_OK;

	/* Parameter tangents of C and C' at the fixed solution parameter. */
	st = qaws_internal_curve_tangent_any(NULL, curve, s.t, QAWS_ZERO,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, param_tangent, &p, &tg, NULL);
	if (st != QAWS_STATUS_OK)
		return st;

	if (!s.at_bound && !s.singular)
	{
		/* dg = C'_dot . r + C' . (C_dot - q_dot);  t_dot = -dg / g_t */
		qaws_scalar dg = qaws_v3_dot(tg.d[1], s.r) + qaws_v3_dot(s.jet.d[1], qaws_v3_sub(tg.d[0], qd));
		t_dot = -dg / s.gt;
	}
	p_dot = qaws_v3_axpy(tg.d[0], s.jet.d[1], t_dot);
	out_tangent->t = t_dot;
	out_tangent->position = p_dot;
	out_tangent->distance = s.distance > QAWS_ZERO
		? qaws_v3_dot(s.r, qaws_v3_sub(p_dot, qd)) / s.distance : QAWS_ZERO;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_closest_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 query,
	qaws_curve_closest_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* query_adjoint)
{
	curve_solution s;
	qaws_curve_jet_3d jbar;
	qaws_vec3 n = qaws_v3_zero();
	qaws_scalar t_total, mu = QAWS_ZERO;
	qaws_status st;

	if (!curve || !adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension == QAWS_DIMENSION_2D)
		query.z = QAWS_ZERO;
	st = curve_solve(curve, query, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	curve_report(ctx, &s);

	if (s.distance > QAWS_ZERO)
		n = qaws_v3_scale(s.r, QAWS_ONE / s.distance);
	t_total = adjoint->t + qaws_v3_dot(adjoint->position, s.jet.d[1]) + adjoint->distance * qaws_v3_dot(n, s.jet.d[1]);
	if (!s.at_bound && !s.singular)
		mu = -t_total / s.gt;

	memset(&jbar, 0, sizeof(jbar));
	jbar.d[0] = qaws_v3_axpy(qaws_v3_axpy(adjoint->position, n, adjoint->distance), s.jet.d[1], mu);
	jbar.d[1] = qaws_v3_scale(s.r, mu);
	jbar.channels = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
	st = qaws_internal_curve_adjoint_any(NULL, curve, s.t, jbar.channels, &jbar, param_adjoint, NULL);
	if (st != QAWS_STATUS_OK)
		return st;

	if (query_adjoint)
	{
		qaws_vec3 qb = qaws_v3_add(qaws_v3_scale(n, -adjoint->distance), qaws_v3_scale(s.jet.d[1], -mu));
		query_adjoint->x += qb.x;
		query_adjoint->y += qb.y;
		if (curve->dimension == QAWS_DIMENSION_3D)
			query_adjoint->z += qb.z;
	}
	return QAWS_STATUS_OK;
}

/* ================================================================== */
/*  Surfaces                                                          */
/* ================================================================== */

typedef struct surface_solution
{
	qaws_scalar u, v;
	qaws_surface_jet jet;     /* S, Su, Sv, Suu, Suv, Svv */
	qaws_vec3 r;
	qaws_scalar g[2], J[2][2], distance, gap;
	int fixed_u, fixed_v, singular;
} surface_solution;

static void surface_system(surface_solution* s, qaws_vec3 q)
{
	qaws_vec3 const* d = s->jet.d;
	s->r = qaws_v3_sub(d[0], q);
	s->g[0] = qaws_v3_dot(d[1], s->r);
	s->g[1] = qaws_v3_dot(d[2], s->r);
	s->J[0][0] = qaws_v3_dot(d[3], s->r) + qaws_v3_dot(d[1], d[1]);
	s->J[0][1] = s->J[1][0] = qaws_v3_dot(d[4], s->r) + qaws_v3_dot(d[1], d[2]);
	s->J[1][1] = qaws_v3_dot(d[5], s->r) + qaws_v3_dot(d[2], d[2]);
}

/* Solves J x = -b on the free coordinates (fixed ones stay zero).
   Returns 0 when the reduced system is singular. */
static int surface_solve_reduced(surface_solution const* s, qaws_scalar const* b, qaws_scalar* x)
{
	x[0] = x[1] = QAWS_ZERO;
	if (s->fixed_u && s->fixed_v)
		return 1;
	if (s->fixed_u)
	{
		if (!(s->J[1][1] > QAWS_ZERO)) return 0;
		x[1] = -b[1] / s->J[1][1];
		return 1;
	}
	if (s->fixed_v)
	{
		if (!(s->J[0][0] > QAWS_ZERO)) return 0;
		x[0] = -b[0] / s->J[0][0];
		return 1;
	}
	{
		qaws_scalar det = s->J[0][0] * s->J[1][1] - s->J[0][1] * s->J[1][0];
		if (!(det > QAWS_ZERO) || !(s->J[0][0] > QAWS_ZERO))
			return 0;
		x[0] = -(s->J[1][1] * b[0] - s->J[0][1] * b[1]) / det;
		x[1] = -(-s->J[1][0] * b[0] + s->J[0][0] * b[1]) / det;
		return 1;
	}
}

/* Newton refinement of the 2x2 optimality system from (u0, v0), holding
   coordinates that sit on an active boundary. */
static qaws_status surface_refine(qaws_surface const* surface, qaws_vec3 q, qaws_scalar u0, qaws_scalar v0, surface_solution* s)
{
	qaws_range ur = surface->u_range, vr = surface->v_range;
	qaws_scalar eu = (ur.max_value - ur.min_value) * QAWS_LITERAL(1e-7);
	qaws_scalar ev = (vr.max_value - vr.min_value) * QAWS_LITERAL(1e-7);
	unsigned int it;
	qaws_status st;

	s->u = u0;
	s->v = v0;
	for (it = 0; it <= OPS_NEWTON_STEPS; it++)
	{
		qaws_scalar x[2];
		st = qaws_surface_eval_jet(surface, s->u, s->v, QAWS_SJET_ORDER2, &s->jet);
		if (st != QAWS_STATUS_OK)
			return st;
		surface_system(s, q);
		s->fixed_u = (s->u <= ur.min_value + eu && s->g[0] > QAWS_ZERO) || (s->u >= ur.max_value - eu && s->g[0] < QAWS_ZERO);
		s->fixed_v = (s->v <= vr.min_value + ev && s->g[1] > QAWS_ZERO) || (s->v >= vr.max_value - ev && s->g[1] < QAWS_ZERO);
		if (it == OPS_NEWTON_STEPS || !surface_solve_reduced(s, s->g, x))
			break;
		s->u += x[0];
		s->v += x[1];
		if (s->u < ur.min_value) s->u = ur.min_value;
		if (s->u > ur.max_value) s->u = ur.max_value;
		if (s->v < vr.min_value) s->v = vr.min_value;
		if (s->v > vr.max_value) s->v = vr.max_value;
	}
	s->distance = QAWS_SQRT(qaws_v3_dot(s->r, s->r));
	return QAWS_STATUS_OK;
}

static qaws_status surface_solve(qaws_surface const* surface, qaws_vec3 q, surface_solution* s)
{
	static qaws_scalar d[OPS_GAP_GRID][OPS_GAP_GRID];
	qaws_range ur = surface->u_range, vr = surface->v_range;
	qaws_scalar du = (ur.max_value - ur.min_value) / (OPS_GAP_GRID - 1);
	qaws_scalar dv = (vr.max_value - vr.min_value) / (OPS_GAP_GRID - 1);
	qaws_scalar u0, v0;
	qaws_vec3 pt;
	qaws_status st;
	int i, j, bi = 0, bj = 0;

	st = qaws_surface_find_closest_point(surface, q, &u0, &v0, &pt);
	if (st != QAWS_STATUS_OK)
		return st;
	st = surface_refine(surface, q, u0, v0, s);
	if (st != QAWS_STATUS_OK)
		return st;

	/* Global check on a sample grid: refine from a closer basin if any. */
	for (i = 0; i < OPS_GAP_GRID; i++)
		for (j = 0; j < OPS_GAP_GRID; j++)
		{
			qaws_surface_eval_result e;
			qaws_vec3 r;
			qaws_surface_evaluate(surface, ur.min_value + du * (qaws_scalar)i, vr.min_value + dv * (qaws_scalar)j,
				QAWS_SURFACE_EVAL_POSITION, &e);
			r = qaws_v3_sub(e.position, q);
			d[i][j] = QAWS_SQRT(qaws_v3_dot(r, r));
			if (d[i][j] < d[bi][bj])
			{
				bi = i;
				bj = j;
			}
		}
	if (d[bi][bj] < s->distance)
	{
		surface_solution alt;
		st = surface_refine(surface, q, ur.min_value + du * (qaws_scalar)bi, vr.min_value + dv * (qaws_scalar)bj, &alt);
		if (st != QAWS_STATUS_OK)
			return st;
		if (alt.distance < s->distance)
			*s = alt;
	}

	{
		qaws_scalar dummy[2] = { 0, 0 }, x[2];
		qaws_scalar metric = qaws_v3_dot(s->jet.d[1], s->jet.d[1]) + qaws_v3_dot(s->jet.d[2], s->jet.d[2]);
		qaws_scalar det = s->J[0][0] * s->J[1][1] - s->J[0][1] * s->J[1][0];
		s->singular = !surface_solve_reduced(s, dummy, x) ||
			(!s->fixed_u && !s->fixed_v && !(det > QAWS_LITERAL(1e-10) * metric * metric));
	}

	/* Best competing local minimum on the grid, away from the solution. */
	s->gap = QAWS_DIFF_UNBOUNDED;
	for (i = 0; i < OPS_GAP_GRID; i++)
		for (j = 0; j < OPS_GAP_GRID; j++)
		{
			int a, b, local_min = 1;
			qaws_scalar ui = ur.min_value + du * (qaws_scalar)i, vj = vr.min_value + dv * (qaws_scalar)j;
			for (a = -1; a <= 1 && local_min; a++)
				for (b = -1; b <= 1; b++)
				{
					int ii = i + a, jj = j + b;
					if ((a || b) && ii >= 0 && jj >= 0 && ii < OPS_GAP_GRID && jj < OPS_GAP_GRID && d[ii][jj] < d[i][j])
					{
						local_min = 0;
						break;
					}
				}
			if (!local_min || (QAWS_FABS(ui - s->u) <= QAWS_LITERAL(2.5) * du && QAWS_FABS(vj - s->v) <= QAWS_LITERAL(2.5) * dv))
				continue;
			if (d[i][j] - s->distance < s->gap)
				s->gap = d[i][j] - s->distance;
		}
	if (s->gap < QAWS_ZERO)
		s->gap = QAWS_ZERO;
	return QAWS_STATUS_OK;
}

static void surface_report(qaws_diff_context const* ctx, surface_solution const* s)
{
	qaws_scalar tr = s->J[0][0] + s->J[1][1];
	qaws_scalar det = s->J[0][0] * s->J[1][1] - s->J[0][1] * s->J[1][0];
	qaws_scalar disc = QAWS_SQRT(qaws_max(tr * tr - 4 * det, QAWS_ZERO));
	qaws_scalar lmax = (tr + disc) / 2, lmin = (tr - disc) / 2;
	qaws_scalar cond = (lmin > QAWS_ZERO) ? lmax / lmin : QAWS_DIFF_UNBOUNDED;
	qaws_scalar scale = QAWS_SQRT(qaws_v3_dot(s->jet.d[1], s->jet.d[1]) + qaws_v3_dot(s->jet.d[2], s->jet.d[2]));
	qaws_scalar res = QAWS_SQRT(s->g[0] * s->g[0] + s->g[1] * s->g[1]) / (scale * s->distance + QAWS_LITERAL(1e-30));
	int frozen = s->fixed_u || s->fixed_v;
	report_implicit(ctx, implicit_validity(s->singular, s->gap, s->distance, frozen),
		frozen ? (unsigned int)QAWS_FREEZE_ACTIVE_SET : 0u, cond, frozen ? QAWS_ZERO : res, s->gap);
}

qaws_status qaws_surface_closest_point_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_vec3 query,
	qaws_vec3 const* query_tangent,
	qaws_diff_views const* param_tangent,
	qaws_surface_closest_point* out_value,
	qaws_surface_closest_point* out_tangent)
{
	surface_solution s;
	qaws_surface_jet p, tg;
	qaws_vec3 qd = query_tangent ? *query_tangent : qaws_v3_zero();
	qaws_scalar dg[2], x[2] = { 0, 0 };
	qaws_vec3 p_dot;
	qaws_status st;

	if (!surface || !out_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = surface_solve(surface, query, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	out_value->u = s.u;
	out_value->v = s.v;
	out_value->position = s.jet.d[0];
	out_value->distance = s.distance;
	surface_report(ctx, &s);
	if (!out_tangent)
		return QAWS_STATUS_OK;

	st = qaws_surface_eval_tangent(NULL, surface, s.u, s.v, QAWS_ZERO, QAWS_ZERO, QAWS_SJET_ORDER1, param_tangent, &p, &tg);
	if (st != QAWS_STATUS_OK)
		return st;
	dg[0] = qaws_v3_dot(tg.d[1], s.r) + qaws_v3_dot(s.jet.d[1], qaws_v3_sub(tg.d[0], qd));
	dg[1] = qaws_v3_dot(tg.d[2], s.r) + qaws_v3_dot(s.jet.d[2], qaws_v3_sub(tg.d[0], qd));
	if (!s.singular)
		surface_solve_reduced(&s, dg, x);
	p_dot = qaws_v3_axpy(qaws_v3_axpy(tg.d[0], s.jet.d[1], x[0]), s.jet.d[2], x[1]);
	out_tangent->u = x[0];
	out_tangent->v = x[1];
	out_tangent->position = p_dot;
	out_tangent->distance = s.distance > QAWS_ZERO
		? qaws_v3_dot(s.r, qaws_v3_sub(p_dot, qd)) / s.distance : QAWS_ZERO;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_closest_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_vec3 query,
	qaws_surface_closest_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* query_adjoint)
{
	surface_solution s;
	qaws_surface_jet jbar;
	qaws_vec3 n = qaws_v3_zero();
	qaws_scalar total[2], lambda[2] = { 0, 0 };
	qaws_status st;

	if (!surface || !adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = surface_solve(surface, query, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	surface_report(ctx, &s);

	if (s.distance > QAWS_ZERO)
		n = qaws_v3_scale(s.r, QAWS_ONE / s.distance);
	total[0] = adjoint->u + qaws_v3_dot(adjoint->position, s.jet.d[1]) + adjoint->distance * qaws_v3_dot(n, s.jet.d[1]);
	total[1] = adjoint->v + qaws_v3_dot(adjoint->position, s.jet.d[2]) + adjoint->distance * qaws_v3_dot(n, s.jet.d[2]);
	/* J is symmetric: lambda = -J^-1 total on the free coordinates. */
	if (!s.singular)
		surface_solve_reduced(&s, total, lambda);

	memset(&jbar, 0, sizeof(jbar));
	jbar.d[0] = qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(adjoint->position, n, adjoint->distance),
		s.jet.d[1], lambda[0]), s.jet.d[2], lambda[1]);
	jbar.d[1] = qaws_v3_scale(s.r, lambda[0]);
	jbar.d[2] = qaws_v3_scale(s.r, lambda[1]);
	jbar.channels = QAWS_SJET_ORDER1;
	st = qaws_surface_eval_adjoint(NULL, surface, s.u, s.v, jbar.channels, &jbar, param_adjoint, NULL, NULL);
	if (st != QAWS_STATUS_OK)
		return st;

	if (query_adjoint)
	{
		qaws_vec3 qb = qaws_v3_scale(n, -adjoint->distance);
		qb = qaws_v3_axpy(qaws_v3_axpy(qb, s.jet.d[1], -lambda[0]), s.jet.d[2], -lambda[1]);
		*query_adjoint = qaws_v3_add(*query_adjoint, qb);
	}
	return QAWS_STATUS_OK;
}

/* ================================================================== */
/*  Roots seeded by the discrete finders                              */
/* ================================================================== */

#define ROOT_CHANNELS (QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)

static qaws_vec3 v3_make(qaws_scalar x, qaws_scalar y, qaws_scalar z)
{
	qaws_vec3 r;
	r.x = x;
	r.y = y;
	r.z = z;
	return r;
}

static qaws_status root_jet(qaws_curve const* curve, qaws_scalar t, qaws_diff_views const* views,
	qaws_curve_jet_3d* primal, qaws_curve_jet_3d* tangent)
{
	qaws_curve_jet_3d unused;
	return qaws_internal_curve_tangent_any(NULL, curve, t, QAWS_ZERO, ROOT_CHANNELS, views, primal,
		tangent ? tangent : &unused, NULL);
}

static qaws_scalar clamp_range(qaws_range r, qaws_scalar t)
{
	if (t < r.min_value) return r.min_value;
	if (t > r.max_value) return r.max_value;
	return t;
}

/* Solves A x = b for n <= 3 (Gaussian elimination, partial pivoting) on a
   copy; transpose selects A^T. Returns 0 when singular relative to scale. */
static int small_solve(qaws_scalar const A[3][3], unsigned int n, int transpose, qaws_scalar const* b, qaws_scalar* x)
{
	double M[3][4], scale = 0;
	unsigned int i, j, k;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			M[i][j] = transpose ? A[j][i] : A[i][j];
			if (fabs(M[i][j]) > scale) scale = fabs(M[i][j]);
		}
		M[i][n] = b[i];
	}
	if (!(scale > 0))
		return 0;
	for (i = 0; i < n; i++)
	{
		unsigned int piv = i;
		for (k = i + 1; k < n; k++)
			if (fabs(M[k][i]) > fabs(M[piv][i]))
				piv = k;
		if (!(fabs(M[piv][i]) > 1e-10 * scale))
			return 0;
		if (piv != i)
			for (j = 0; j <= n; j++)
			{
				double tmp = M[i][j];
				M[i][j] = M[piv][j];
				M[piv][j] = tmp;
			}
		for (k = i + 1; k < n; k++)
		{
			double f = M[k][i] / M[i][i];
			for (j = i; j <= n; j++)
				M[k][j] -= f * M[i][j];
		}
	}
	for (i = n; i-- > 0;)
	{
		double v = M[i][n];
		for (j = i + 1; j < n; j++)
			v -= M[i][j] * (double)x[j];
		x[i] = (qaws_scalar)(v / M[i][i]);
	}
	return 1;
}

/* Frobenius estimate of |A| |A^-1| from n solves with unit vectors. */
static qaws_scalar small_condition(qaws_scalar const A[3][3], unsigned int n)
{
	double na = 0, ni = 0;
	unsigned int i, j;
	for (i = 0; i < n; i++)
	{
		qaws_scalar e[3] = { 0, 0, 0 }, x[3] = { 0, 0, 0 };
		e[i] = QAWS_ONE;
		if (!small_solve(A, n, 0, e, x))
			return QAWS_DIFF_UNBOUNDED;
		for (j = 0; j < n; j++)
		{
			na += (double)A[i][j] * A[i][j];
			ni += (double)x[j] * x[j];
		}
	}
	return (qaws_scalar)sqrt(na * ni);
}

static void report_root(qaws_diff_context const* ctx, int singular, int frozen, qaws_scalar condition, qaws_scalar residual)
{
	report_implicit(ctx, singular ? QAWS_DIFF_ILL_CONDITIONED : (frozen ? QAWS_DIFF_VALID_LOCALLY : QAWS_DIFF_VALID),
		(unsigned int)QAWS_FREEZE_ACTIVE_SET, singular ? QAWS_DIFF_UNBOUNDED : condition, residual, QAWS_DIFF_UNBOUNDED);
}

/* ------------------------------------------------------------------ */
/*  Curve pair                                                        */
/* ------------------------------------------------------------------ */

typedef struct pair_solution
{
	qaws_scalar t[2];
	qaws_curve_jet_3d ja, jb;
	qaws_vec3 r;
	qaws_scalar g[2], H[3][3], distance;
	int fixed[2], singular;
} pair_solution;

static void pair_system(pair_solution* s)
{
	qaws_vec3 const* a = s->ja.d;
	qaws_vec3 const* b = s->jb.d;
	memset(s->H, 0, sizeof(s->H));
	s->r = qaws_v3_sub(a[0], b[0]);
	s->g[0] = qaws_v3_dot(a[1], s->r);
	s->g[1] = -qaws_v3_dot(b[1], s->r);
	s->H[0][0] = qaws_v3_dot(a[1], a[1]) + qaws_v3_dot(a[2], s->r);
	s->H[0][1] = s->H[1][0] = -qaws_v3_dot(a[1], b[1]);
	s->H[1][1] = qaws_v3_dot(b[1], b[1]) - qaws_v3_dot(b[2], s->r);
}

/* x = -H^-1 b on the free coordinates; fixed ones stay zero. */
static int pair_solve(pair_solution const* s, qaws_scalar const* b, qaws_scalar* x)
{
	x[0] = x[1] = QAWS_ZERO;
	if (s->fixed[0] && s->fixed[1])
		return 1;
	if (s->fixed[0] || s->fixed[1])
	{
		unsigned int i = s->fixed[0] ? 1u : 0u;
		if (!(s->H[i][i] > QAWS_LITERAL(1e-12) * (QAWS_ONE + QAWS_FABS(b[i]))))
			return 0;
		x[i] = -b[i] / s->H[i][i];
		return 1;
	}
	{
		qaws_scalar nb[3] = { -b[0], -b[1], 0 };
		return small_solve(s->H, 2, 0, nb, x);
	}
}

static qaws_status pair_refine(qaws_curve const* ca, qaws_curve const* cb, qaws_scalar ta, qaws_scalar tb, pair_solution* s)
{
	qaws_range ra = ca->parameter_range, rb = cb->parameter_range;
	qaws_scalar ea = (ra.max_value - ra.min_value) * QAWS_LITERAL(1e-7);
	qaws_scalar eb = (rb.max_value - rb.min_value) * QAWS_LITERAL(1e-7);
	unsigned int it;
	qaws_status st;
	s->t[0] = clamp_range(ra, ta);
	s->t[1] = clamp_range(rb, tb);
	for (it = 0; it <= OPS_NEWTON_STEPS + 2; it++)
	{
		qaws_scalar x[2];
		st = root_jet(ca, s->t[0], NULL, &s->ja, NULL);
		if (st == QAWS_STATUS_OK)
			st = root_jet(cb, s->t[1], NULL, &s->jb, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		pair_system(s);
		s->fixed[0] = (s->t[0] <= ra.min_value + ea && s->g[0] > QAWS_ZERO) || (s->t[0] >= ra.max_value - ea && s->g[0] < QAWS_ZERO);
		s->fixed[1] = (s->t[1] <= rb.min_value + eb && s->g[1] > QAWS_ZERO) || (s->t[1] >= rb.max_value - eb && s->g[1] < QAWS_ZERO);
		s->singular = !pair_solve(s, s->g, x);
		if (it == OPS_NEWTON_STEPS + 2 || s->singular)
			break;
		s->t[0] = clamp_range(ra, s->t[0] + x[0]);
		s->t[1] = clamp_range(rb, s->t[1] + x[1]);
	}
	s->distance = QAWS_SQRT(qaws_v3_dot(s->r, s->r));
	return QAWS_STATUS_OK;
}

static void pair_report(qaws_diff_context const* ctx, pair_solution const* s)
{
	qaws_scalar speed = QAWS_SQRT(qaws_v3_dot(s->ja.d[1], s->ja.d[1]) + qaws_v3_dot(s->jb.d[1], s->jb.d[1]));
	qaws_scalar res = (QAWS_FABS(s->g[0]) + QAWS_FABS(s->g[1])) / (speed * (s->distance + QAWS_LITERAL(1e-30)) + QAWS_LITERAL(1e-30));
	report_root(ctx, s->singular, s->fixed[0] || s->fixed[1], small_condition(s->H, 2), s->distance > QAWS_ZERO ? res : QAWS_ZERO);
}

static qaws_vec3 unit_or_zero(qaws_vec3 r, qaws_scalar len)
{
	return len > QAWS_ZERO ? qaws_v3_scale(r, QAWS_ONE / len) : qaws_v3_zero();
}

qaws_status qaws_curve_pair_point_tangent(
	qaws_diff_context const* ctx, qaws_curve const* curve_a, qaws_curve const* curve_b,
	qaws_scalar t_a_seed, qaws_scalar t_b_seed,
	qaws_diff_views const* tangent_a, qaws_diff_views const* tangent_b,
	qaws_curve_pair_point* out_value, qaws_curve_pair_point* out_tangent)
{
	pair_solution s;
	qaws_curve_jet_3d pa, ta, pb, tb;
	qaws_scalar dg[2], x[2] = { 0, 0 };
	qaws_vec3 dr, n;
	qaws_status st;

	if (!curve_a || !curve_b || !out_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = pair_refine(curve_a, curve_b, t_a_seed, t_b_seed, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	out_value->t_a = s.t[0];
	out_value->t_b = s.t[1];
	out_value->position_a = s.ja.d[0];
	out_value->position_b = s.jb.d[0];
	out_value->distance = s.distance;
	pair_report(ctx, &s);
	if (!out_tangent)
		return QAWS_STATUS_OK;

	st = root_jet(curve_a, s.t[0], tangent_a, &pa, &ta);
	if (st == QAWS_STATUS_OK)
		st = root_jet(curve_b, s.t[1], tangent_b, &pb, &tb);
	if (st != QAWS_STATUS_OK)
		return st;
	dr = qaws_v3_sub(ta.d[0], tb.d[0]);
	dg[0] = qaws_v3_dot(ta.d[1], s.r) + qaws_v3_dot(s.ja.d[1], dr);
	dg[1] = -qaws_v3_dot(tb.d[1], s.r) - qaws_v3_dot(s.jb.d[1], dr);
	if (!s.singular)
		pair_solve(&s, dg, x);
	out_tangent->t_a = x[0];
	out_tangent->t_b = x[1];
	out_tangent->position_a = qaws_v3_axpy(ta.d[0], s.ja.d[1], x[0]);
	out_tangent->position_b = qaws_v3_axpy(tb.d[0], s.jb.d[1], x[1]);
	n = unit_or_zero(s.r, s.distance);
	out_tangent->distance = qaws_v3_dot(n, qaws_v3_sub(out_tangent->position_a, out_tangent->position_b));
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_pair_point_adjoint(
	qaws_diff_context const* ctx, qaws_curve const* curve_a, qaws_curve const* curve_b,
	qaws_scalar t_a_seed, qaws_scalar t_b_seed, qaws_curve_pair_point const* adjoint,
	qaws_diff_views* adjoint_a, qaws_diff_views* adjoint_b)
{
	pair_solution s;
	qaws_curve_jet_3d ja, jb;
	qaws_scalar tau[2], mu[2] = { 0, 0 };
	qaws_vec3 n;
	unsigned int ch = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
	qaws_status st;

	if (!curve_a || !curve_b || !adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = pair_refine(curve_a, curve_b, t_a_seed, t_b_seed, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	pair_report(ctx, &s);
	n = unit_or_zero(s.r, s.distance);
	tau[0] = adjoint->t_a + qaws_v3_dot(adjoint->position_a, s.ja.d[1]) + adjoint->distance * qaws_v3_dot(n, s.ja.d[1]);
	tau[1] = adjoint->t_b + qaws_v3_dot(adjoint->position_b, s.jb.d[1]) - adjoint->distance * qaws_v3_dot(n, s.jb.d[1]);
	if (!s.singular)
		pair_solve(&s, tau, mu);   /* mu = -H^-1 tau (H symmetric) */

	memset(&ja, 0, sizeof(ja));
	memset(&jb, 0, sizeof(jb));
	ja.d[0] = qaws_v3_add(qaws_v3_axpy(adjoint->position_a, n, adjoint->distance),
		qaws_v3_sub(qaws_v3_scale(s.ja.d[1], mu[0]), qaws_v3_scale(s.jb.d[1], mu[1])));
	ja.d[1] = qaws_v3_scale(s.r, mu[0]);
	ja.channels = ch;
	jb.d[0] = qaws_v3_sub(qaws_v3_axpy(adjoint->position_b, n, -adjoint->distance),
		qaws_v3_sub(qaws_v3_scale(s.ja.d[1], mu[0]), qaws_v3_scale(s.jb.d[1], mu[1])));
	jb.d[1] = qaws_v3_scale(s.r, -mu[1]);
	jb.channels = ch;
	if (adjoint_a)
	{
		st = qaws_internal_curve_adjoint_any(NULL, curve_a, s.t[0], ch, &ja, adjoint_a, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	if (adjoint_b)
		st = qaws_internal_curve_adjoint_any(NULL, curve_b, s.t[1], ch, &jb, adjoint_b, NULL);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Scalar roots in t: plane crossing, extremum, inflection           */
/* ------------------------------------------------------------------ */

typedef qaws_scalar (*root_fn)(void const* user, qaws_curve_jet_3d const* j, qaws_scalar* g_t);

static qaws_status root_refine(qaws_curve const* curve, qaws_scalar t0, root_fn fn, void const* user,
	qaws_scalar* t, qaws_curve_jet_3d* jet, qaws_scalar* g, qaws_scalar* g_t)
{
	qaws_range range = curve->parameter_range;
	unsigned int it;
	*t = clamp_range(range, t0);
	for (it = 0; it <= OPS_NEWTON_STEPS + 2; it++)
	{
		qaws_status st = root_jet(curve, *t, NULL, jet, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		*g = fn(user, jet, g_t);
		if (it == OPS_NEWTON_STEPS + 2 || !(QAWS_FABS(*g_t) > QAWS_ZERO))
			break;
		*t = clamp_range(range, *t - *g / *g_t);
	}
	return QAWS_STATUS_OK;
}

static int root_singular(qaws_scalar g_t, qaws_scalar scale)
{
	return !(QAWS_FABS(g_t) > QAWS_LITERAL(1e-8) * scale);
}

/* plane */
static qaws_scalar plane_fn(void const* user, qaws_curve_jet_3d const* j, qaws_scalar* g_t)
{
	qaws_plane const* p = (qaws_plane const*)user;
	*g_t = qaws_v3_dot(p->normal, j->d[1]);
	return qaws_v3_dot(p->normal, qaws_v3_sub(j->d[0], p->point));
}

qaws_status qaws_curve_plane_point_tangent(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_plane const* plane, qaws_scalar t_seed,
	qaws_plane const* plane_tangent, qaws_diff_views const* param_tangent,
	qaws_curve_plane_point* out_value, qaws_curve_plane_point* out_tangent)
{
	qaws_curve_jet_3d j, p, tg;
	qaws_scalar t, g, g_t, t_dot = QAWS_ZERO, scale;
	int singular;
	qaws_status st;

	if (!curve || !plane || !out_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = root_refine(curve, t_seed, plane_fn, plane, &t, &j, &g, &g_t);
	if (st != QAWS_STATUS_OK)
		return st;
	scale = QAWS_SQRT(qaws_v3_dot(plane->normal, plane->normal) * qaws_v3_dot(j.d[1], j.d[1]));
	singular = root_singular(g_t, scale);
	report_root(ctx, singular, 0, singular ? QAWS_DIFF_UNBOUNDED : scale / QAWS_FABS(g_t),
		QAWS_FABS(g) / (QAWS_SQRT(qaws_v3_dot(plane->normal, plane->normal)) + QAWS_LITERAL(1e-30)));
	out_value->t = t;
	out_value->position = j.d[0];
	if (!out_tangent)
		return QAWS_STATUS_OK;
	st = root_jet(curve, t, param_tangent, &p, &tg);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!singular)
	{
		qaws_scalar dg = qaws_v3_dot(plane->normal, tg.d[0]);
		if (plane_tangent)
			dg += qaws_v3_dot(plane_tangent->normal, qaws_v3_sub(j.d[0], plane->point)) -
				qaws_v3_dot(plane->normal, plane_tangent->point);
		t_dot = -dg / g_t;
	}
	out_tangent->t = t_dot;
	out_tangent->position = qaws_v3_axpy(tg.d[0], j.d[1], t_dot);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_plane_point_adjoint(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_plane const* plane, qaws_scalar t_seed,
	qaws_curve_plane_point const* adjoint, qaws_diff_views* param_adjoint, qaws_plane* plane_adjoint)
{
	qaws_curve_jet_3d j, jbar;
	qaws_scalar t, g, g_t, mu = QAWS_ZERO, scale;
	int singular;
	qaws_status st;

	if (!curve || !plane || !adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = root_refine(curve, t_seed, plane_fn, plane, &t, &j, &g, &g_t);
	if (st != QAWS_STATUS_OK)
		return st;
	scale = QAWS_SQRT(qaws_v3_dot(plane->normal, plane->normal) * qaws_v3_dot(j.d[1], j.d[1]));
	singular = root_singular(g_t, scale);
	report_root(ctx, singular, 0, singular ? QAWS_DIFF_UNBOUNDED : scale / QAWS_FABS(g_t), QAWS_ZERO);
	if (!singular)
		mu = -(adjoint->t + qaws_v3_dot(adjoint->position, j.d[1])) / g_t;
	memset(&jbar, 0, sizeof(jbar));
	jbar.d[0] = qaws_v3_axpy(adjoint->position, plane->normal, mu);
	jbar.channels = QAWS_EVAL_FLAG_POSITION;
	if (plane_adjoint)
	{
		plane_adjoint->normal = qaws_v3_axpy(plane_adjoint->normal, qaws_v3_sub(j.d[0], plane->point), mu);
		plane_adjoint->point = qaws_v3_axpy(plane_adjoint->point, plane->normal, -mu);
	}
	if (!param_adjoint)
		return QAWS_STATUS_OK;
	return qaws_internal_curve_adjoint_any(NULL, curve, t, jbar.channels, &jbar, param_adjoint, NULL);
}

/* extremum */
static qaws_scalar extremum_fn(void const* user, qaws_curve_jet_3d const* j, qaws_scalar* g_t)
{
	qaws_vec3 const* e = (qaws_vec3 const*)user;
	*g_t = qaws_v3_dot(*e, j->d[2]);
	return qaws_v3_dot(*e, j->d[1]);
}

static qaws_vec3 flat_if_2d(qaws_curve const* curve, qaws_vec3 v)
{
	if (curve->dimension == QAWS_DIMENSION_2D)
		v.z = QAWS_ZERO;
	return v;
}

qaws_status qaws_curve_extremum_tangent(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_vec3 direction, qaws_scalar t_seed,
	qaws_vec3 const* direction_tangent, qaws_diff_views const* param_tangent,
	qaws_curve_extremum* out_value, qaws_curve_extremum* out_tangent)
{
	qaws_curve_jet_3d j, p, tg;
	qaws_scalar t, g, g_t, t_dot = QAWS_ZERO, scale;
	qaws_vec3 de = direction_tangent ? flat_if_2d(curve, *direction_tangent) : qaws_v3_zero();
	int singular;
	qaws_status st;

	if (!curve || !out_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	direction = flat_if_2d(curve, direction);
	st = root_refine(curve, t_seed, extremum_fn, &direction, &t, &j, &g, &g_t);
	if (st != QAWS_STATUS_OK)
		return st;
	scale = QAWS_SQRT(qaws_v3_dot(direction, direction) * qaws_v3_dot(j.d[2], j.d[2]));
	singular = root_singular(g_t, scale);
	report_root(ctx, singular, 0, singular ? QAWS_DIFF_UNBOUNDED : scale / QAWS_FABS(g_t),
		QAWS_FABS(g) / (QAWS_SQRT(qaws_v3_dot(direction, direction) * qaws_v3_dot(j.d[1], j.d[1])) + QAWS_LITERAL(1e-30)));
	out_value->t = t;
	out_value->position = j.d[0];
	out_value->value = qaws_v3_dot(direction, j.d[0]);
	if (!out_tangent)
		return QAWS_STATUS_OK;
	st = root_jet(curve, t, param_tangent, &p, &tg);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!singular)
		t_dot = -(qaws_v3_dot(direction, tg.d[1]) + qaws_v3_dot(de, j.d[1])) / g_t;
	out_tangent->t = t_dot;
	out_tangent->position = qaws_v3_axpy(tg.d[0], j.d[1], t_dot);
	out_tangent->value = qaws_v3_dot(de, j.d[0]) + qaws_v3_dot(direction, out_tangent->position);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_extremum_adjoint(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_vec3 direction, qaws_scalar t_seed,
	qaws_curve_extremum const* adjoint, qaws_diff_views* param_adjoint, qaws_vec3* direction_adjoint)
{
	qaws_curve_jet_3d j, jbar;
	qaws_scalar t, g, g_t, mu = QAWS_ZERO, scale, tau;
	int singular;
	qaws_status st;

	if (!curve || !adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	direction = flat_if_2d(curve, direction);
	st = root_refine(curve, t_seed, extremum_fn, &direction, &t, &j, &g, &g_t);
	if (st != QAWS_STATUS_OK)
		return st;
	scale = QAWS_SQRT(qaws_v3_dot(direction, direction) * qaws_v3_dot(j.d[2], j.d[2]));
	singular = root_singular(g_t, scale);
	report_root(ctx, singular, 0, singular ? QAWS_DIFF_UNBOUNDED : scale / QAWS_FABS(g_t), QAWS_ZERO);
	/* value = e . C(t): its t-derivative e . C' is the residual g */
	tau = adjoint->t + qaws_v3_dot(adjoint->position, j.d[1]) + adjoint->value * g;
	if (!singular)
		mu = -tau / g_t;
	memset(&jbar, 0, sizeof(jbar));
	jbar.d[0] = qaws_v3_axpy(adjoint->position, direction, adjoint->value);
	jbar.d[1] = qaws_v3_scale(direction, mu);
	jbar.channels = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
	if (direction_adjoint)
	{
		qaws_vec3 eb = qaws_v3_axpy(qaws_v3_scale(j.d[1], mu), j.d[0], adjoint->value);
		eb = flat_if_2d(curve, eb);
		*direction_adjoint = qaws_v3_add(*direction_adjoint, eb);
	}
	if (!param_adjoint)
		return QAWS_STATUS_OK;
	return qaws_internal_curve_adjoint_any(NULL, curve, t, jbar.channels, &jbar, param_adjoint, NULL);
}

/* inflection */
static qaws_scalar inflection_fn(void const* user, qaws_curve_jet_3d const* j, qaws_scalar* g_t)
{
	(void)user;
	*g_t = j->d[1].x * j->d[3].y - j->d[1].y * j->d[3].x;
	return j->d[1].x * j->d[2].y - j->d[1].y * j->d[2].x;
}

static qaws_scalar inflection_scale(qaws_curve_jet_3d const* j)
{
	return QAWS_SQRT((j->d[1].x * j->d[1].x + j->d[1].y * j->d[1].y) * (j->d[3].x * j->d[3].x + j->d[3].y * j->d[3].y));
}

qaws_status qaws_curve_inflection_tangent(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_scalar t_seed,
	qaws_diff_views const* param_tangent, qaws_curve_inflection* out_value, qaws_curve_inflection* out_tangent)
{
	qaws_curve_jet_3d j, p, tg;
	qaws_scalar t, g, g_t, t_dot = QAWS_ZERO, scale;
	int singular;
	qaws_status st;

	if (!curve || !out_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = root_refine(curve, t_seed, inflection_fn, NULL, &t, &j, &g, &g_t);
	if (st != QAWS_STATUS_OK)
		return st;
	scale = inflection_scale(&j);
	singular = root_singular(g_t, scale);
	report_root(ctx, singular, 0, singular ? QAWS_DIFF_UNBOUNDED : scale / QAWS_FABS(g_t), QAWS_ZERO);
	out_value->t = t;
	out_value->position = j.d[0];
	if (!out_tangent)
		return QAWS_STATUS_OK;
	st = root_jet(curve, t, param_tangent, &p, &tg);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!singular)
	{
		qaws_scalar dg = tg.d[1].x * j.d[2].y + j.d[1].x * tg.d[2].y - tg.d[1].y * j.d[2].x - j.d[1].y * tg.d[2].x;
		t_dot = -dg / g_t;
	}
	out_tangent->t = t_dot;
	out_tangent->position = qaws_v3_axpy(tg.d[0], j.d[1], t_dot);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_inflection_adjoint(
	qaws_diff_context const* ctx, qaws_curve const* curve, qaws_scalar t_seed,
	qaws_curve_inflection const* adjoint, qaws_diff_views* param_adjoint)
{
	qaws_curve_jet_3d j, jbar;
	qaws_scalar t, g, g_t, mu = QAWS_ZERO, scale;
	int singular;
	qaws_status st;

	if (!curve || !adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = root_refine(curve, t_seed, inflection_fn, NULL, &t, &j, &g, &g_t);
	if (st != QAWS_STATUS_OK)
		return st;
	scale = inflection_scale(&j);
	singular = root_singular(g_t, scale);
	report_root(ctx, singular, 0, singular ? QAWS_DIFF_UNBOUNDED : scale / QAWS_FABS(g_t), QAWS_ZERO);
	if (!singular)
		mu = -(adjoint->t + qaws_v3_dot(adjoint->position, j.d[1])) / g_t;
	memset(&jbar, 0, sizeof(jbar));
	jbar.d[0] = adjoint->position;
	jbar.d[1] = v3_make(mu * j.d[2].y, -mu * j.d[2].x, QAWS_ZERO);
	jbar.d[2] = v3_make(-mu * j.d[1].y, mu * j.d[1].x, QAWS_ZERO);
	jbar.channels = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;
	if (!param_adjoint)
		return QAWS_STATUS_OK;
	return qaws_internal_curve_adjoint_any(NULL, curve, t, jbar.channels, &jbar, param_adjoint, NULL);
}

/* ------------------------------------------------------------------ */
/*  Surface-curve piercing                                            */
/* ------------------------------------------------------------------ */

typedef struct pierce_solution
{
	qaws_scalar u, v, t;
	qaws_surface_jet sj;
	qaws_curve_jet_3d cj;
	qaws_scalar J[3][3], g[3];
	int singular;
} pierce_solution;

static void pierce_system(pierce_solution* s)
{
	qaws_vec3 r = qaws_v3_sub(s->sj.d[0], s->cj.d[0]);
	qaws_vec3 const* su = &s->sj.d[1];
	qaws_vec3 const* sv = &s->sj.d[2];
	qaws_vec3 const* ct = &s->cj.d[1];
	s->g[0] = r.x; s->g[1] = r.y; s->g[2] = r.z;
	s->J[0][0] = su->x; s->J[0][1] = sv->x; s->J[0][2] = -ct->x;
	s->J[1][0] = su->y; s->J[1][1] = sv->y; s->J[1][2] = -ct->y;
	s->J[2][0] = su->z; s->J[2][1] = sv->z; s->J[2][2] = -ct->z;
}

static qaws_status pierce_refine(qaws_surface const* surface, qaws_curve const* curve,
	qaws_scalar u0, qaws_scalar v0, qaws_scalar t0, pierce_solution* s)
{
	unsigned int it;
	qaws_status st;
	s->u = clamp_range(surface->u_range, u0);
	s->v = clamp_range(surface->v_range, v0);
	s->t = clamp_range(curve->parameter_range, t0);
	for (it = 0; it <= OPS_NEWTON_STEPS + 2; it++)
	{
		qaws_scalar x[3] = { 0, 0, 0 }, nb[3];
		st = qaws_surface_eval_jet(surface, s->u, s->v, QAWS_SJET_ORDER1, &s->sj);
		if (st == QAWS_STATUS_OK)
			st = root_jet(curve, s->t, NULL, &s->cj, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		pierce_system(s);
		nb[0] = -s->g[0]; nb[1] = -s->g[1]; nb[2] = -s->g[2];
		s->singular = !small_solve(s->J, 3, 0, nb, x);
		if (it == OPS_NEWTON_STEPS + 2 || s->singular)
			break;
		s->u = clamp_range(surface->u_range, s->u + x[0]);
		s->v = clamp_range(surface->v_range, s->v + x[1]);
		s->t = clamp_range(curve->parameter_range, s->t + x[2]);
	}
	return QAWS_STATUS_OK;
}

static void pierce_report(qaws_diff_context const* ctx, pierce_solution const* s)
{
	qaws_scalar res = QAWS_SQRT(s->g[0] * s->g[0] + s->g[1] * s->g[1] + s->g[2] * s->g[2]);
	report_root(ctx, s->singular, 0, s->singular ? QAWS_DIFF_UNBOUNDED : small_condition(s->J, 3), res);
}

qaws_status qaws_surface_curve_point_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_curve const* curve,
	qaws_scalar u_seed, qaws_scalar v_seed, qaws_scalar t_seed,
	qaws_diff_views const* surface_tangent, qaws_diff_views const* curve_tangent,
	qaws_surface_curve_point* out_value, qaws_surface_curve_point* out_tangent)
{
	pierce_solution s;
	qaws_surface_jet sp, stg;
	qaws_curve_jet_3d cp, ctg;
	qaws_scalar b[3], x[3] = { 0, 0, 0 };
	qaws_status st;

	if (!surface || !curve || !out_value)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	st = pierce_refine(surface, curve, u_seed, v_seed, t_seed, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	pierce_report(ctx, &s);
	out_value->u = s.u;
	out_value->v = s.v;
	out_value->t = s.t;
	out_value->position = s.cj.d[0];
	if (!out_tangent)
		return QAWS_STATUS_OK;
	st = qaws_surface_eval_tangent(NULL, surface, s.u, s.v, QAWS_ZERO, QAWS_ZERO, QAWS_SJET_P, surface_tangent, &sp, &stg);
	if (st == QAWS_STATUS_OK)
		st = root_jet(curve, s.t, curve_tangent, &cp, &ctg);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!surface_tangent)
		stg.d[0] = qaws_v3_zero();
	b[0] = -(stg.d[0].x - ctg.d[0].x);
	b[1] = -(stg.d[0].y - ctg.d[0].y);
	b[2] = -(stg.d[0].z - ctg.d[0].z);
	if (!s.singular)
		small_solve(s.J, 3, 0, b, x);
	out_tangent->u = x[0];
	out_tangent->v = x[1];
	out_tangent->t = x[2];
	out_tangent->position = qaws_v3_axpy(ctg.d[0], s.cj.d[1], x[2]);
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_curve_point_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_curve const* curve,
	qaws_scalar u_seed, qaws_scalar v_seed, qaws_scalar t_seed, qaws_surface_curve_point const* adjoint,
	qaws_diff_views* surface_adjoint, qaws_diff_views* curve_adjoint)
{
	pierce_solution s;
	qaws_scalar tau[3], lam[3] = { 0, 0, 0 };
	qaws_vec3 l;
	qaws_status st;

	if (!surface || !curve || !adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	st = pierce_refine(surface, curve, u_seed, v_seed, t_seed, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	pierce_report(ctx, &s);
	/* lambda = -J^-T tau */
	tau[0] = -adjoint->u;
	tau[1] = -adjoint->v;
	tau[2] = -(adjoint->t + qaws_v3_dot(adjoint->position, s.cj.d[1]));
	if (!s.singular)
		small_solve(s.J, 3, 1, tau, lam);
	l = v3_make(lam[0], lam[1], lam[2]);
	if (surface_adjoint)
	{
		qaws_surface_jet sb;
		memset(&sb, 0, sizeof(sb));
		sb.d[0] = l;
		sb.channels = QAWS_SJET_P;
		st = qaws_surface_eval_adjoint(NULL, surface, s.u, s.v, QAWS_SJET_P, &sb, surface_adjoint, NULL, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	if (curve_adjoint)
	{
		qaws_curve_jet_3d cb;
		memset(&cb, 0, sizeof(cb));
		cb.d[0] = qaws_v3_sub(adjoint->position, l);
		cb.channels = QAWS_EVAL_FLAG_POSITION;
		st = qaws_internal_curve_adjoint_any(NULL, curve, s.t, cb.channels, &cb, curve_adjoint, NULL);
	}
	return st;
}
