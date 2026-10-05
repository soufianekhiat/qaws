#ifndef QAWS_DIFF_OPS_H
#define QAWS_DIFF_OPS_H

#include "qaws_diff_types.h"
#include "qaws_surface_types.h"

/*
 * Differentiable geometric operations defined by a solve.
 *
 * The solution x* of g(x, theta, q) = 0 is differentiated with the
 * implicit function theorem, never by differentiating solver iterations:
 *
 *     d(x*) / d(theta) = -(dg / dx)^-1 (dg / dtheta)
 *
 * Every call reports through ctx->report (when set):
 *   condition_number  conditioning of dg/dx relative to the metric
 *   residual          relative |g| at the returned solution
 *   branch_gap        distance gap to the best competing local solution
 *   validity          ILL_CONDITIONED near singular dg/dx,
 *                     AMBIGUOUS when a competing solution is within 0.5%
 *                     of the distance, VALID_LOCALLY (+ FREEZE_ACTIVE_SET) when
 *                     the solution sits on the domain boundary and the
 *                     boundary coordinate is held fixed.
 */

/* ===================================================================
 * Closest point on a curve
 *
 * Optimality: C'(t) . (C(t) - q) = 0. Works for 2D curves too (z of the
 * query and of all outputs is then zero).
 * =================================================================== */

typedef struct qaws_curve_closest_point
{
	qaws_scalar t;
	qaws_vec3 position;
	qaws_scalar distance;
} qaws_curve_closest_point;

/* Solves the closest point to `query` and returns its tangent for a query
   tangent and/or a parameter tangent (either may be NULL / zero). */
qaws_status qaws_curve_closest_point_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 query,
	qaws_vec3 const* query_tangent,
	qaws_diff_views const* param_tangent,
	qaws_curve_closest_point* out_value,
	qaws_curve_closest_point* out_tangent);

/* Solves the closest point and accumulates (+=) the pullback of an output
   adjoint (t, position, distance) into parameter and query adjoints. */
qaws_status qaws_curve_closest_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 query,
	qaws_curve_closest_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* query_adjoint);

/* ===================================================================
 * Closest point on a surface
 *
 * Optimality: Su . (S - q) = 0 and Sv . (S - q) = 0.
 * =================================================================== */

typedef struct qaws_surface_closest_point
{
	qaws_scalar u, v;
	qaws_vec3 position;
	qaws_scalar distance;
} qaws_surface_closest_point;

qaws_status qaws_surface_closest_point_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_vec3 query,
	qaws_vec3 const* query_tangent,
	qaws_diff_views const* param_tangent,
	qaws_surface_closest_point* out_value,
	qaws_surface_closest_point* out_tangent);

qaws_status qaws_surface_closest_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_vec3 query,
	qaws_surface_closest_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* query_adjoint);

#endif /* QAWS_DIFF_OPS_H */
