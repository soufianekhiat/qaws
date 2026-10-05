#include "qaws_diff_ops.h"
#include "qaws_inspect.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include <string.h>

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
