#ifndef QAWS_DIFF_OPS_H
#define QAWS_DIFF_OPS_H

#include "qaws_diff_types.h"
#include "qaws_surface_types.h"
#include "qaws_inspect.h"

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

/* ===================================================================
 * Roots seeded by the discrete finders
 *
 * The operations below differentiate one solution: the caller passes a
 * seed (typically from qaws_curve_find_intersections_*,
 * qaws_curve_find_plane_intersections, qaws_surface_find_curve_
 * intersections, qaws_curve_find_extrema or qaws_curve_find_inflection_
 * points), the solution is refined by Newton and differentiated by the
 * implicit function theorem. Which solution exists, and how many, is
 * frozen (QAWS_FREEZE_ACTIVE_SET); tangential configurations report
 * QAWS_DIFF_ILL_CONDITIONED and return zero coordinate tangents.
 * =================================================================== */

/* Closest approach of two curves (an intersection when distance = 0).
   Stationarity of |A(ta) - B(tb)|^2 / 2:
     A' . (A - B) = 0,  B' . (A - B) = 0.
   A parameter on its range end with the gradient pointing outward is
   held fixed. Two views select the parameters of A and of B. */
typedef struct qaws_curve_pair_point
{
	qaws_scalar t_a, t_b;
	qaws_vec3 position_a, position_b;
	qaws_scalar distance;
} qaws_curve_pair_point;

qaws_status qaws_curve_pair_point_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_scalar t_a_seed,
	qaws_scalar t_b_seed,
	qaws_diff_views const* tangent_a,
	qaws_diff_views const* tangent_b,
	qaws_curve_pair_point* out_value,
	qaws_curve_pair_point* out_tangent);

qaws_status qaws_curve_pair_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_scalar t_a_seed,
	qaws_scalar t_b_seed,
	qaws_curve_pair_point const* adjoint,
	qaws_diff_views* adjoint_a,
	qaws_diff_views* adjoint_b);

/* Curve crossing a plane: n . (C(t) - p) = 0. The plane is an input with
   its own tangent / adjoint (point and normal). */
typedef struct qaws_curve_plane_point
{
	qaws_scalar t;
	qaws_vec3 position;
} qaws_curve_plane_point;

qaws_status qaws_curve_plane_point_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_plane const* plane,
	qaws_scalar t_seed,
	qaws_plane const* plane_tangent,
	qaws_diff_views const* param_tangent,
	qaws_curve_plane_point* out_value,
	qaws_curve_plane_point* out_tangent);

qaws_status qaws_curve_plane_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_plane const* plane,
	qaws_scalar t_seed,
	qaws_curve_plane_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_plane* plane_adjoint);

/* Curve piercing a surface: S(u, v) - C(t) = 0. */
typedef struct qaws_surface_curve_point
{
	qaws_scalar u, v, t;
	qaws_vec3 position;
} qaws_surface_curve_point;

qaws_status qaws_surface_curve_point_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_curve const* curve,
	qaws_scalar u_seed,
	qaws_scalar v_seed,
	qaws_scalar t_seed,
	qaws_diff_views const* surface_tangent,
	qaws_diff_views const* curve_tangent,
	qaws_surface_curve_point* out_value,
	qaws_surface_curve_point* out_tangent);

qaws_status qaws_surface_curve_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_curve const* curve,
	qaws_scalar u_seed,
	qaws_scalar v_seed,
	qaws_scalar t_seed,
	qaws_surface_curve_point const* adjoint,
	qaws_diff_views* surface_adjoint,
	qaws_diff_views* curve_adjoint);

/* Extremum of the height e . C(t) along a direction: e . C'(t) = 0.
   value = e . C(t). The direction is an input (tangent / adjoint). */
typedef struct qaws_curve_extremum
{
	qaws_scalar t;
	qaws_vec3 position;
	qaws_scalar value;
} qaws_curve_extremum;

qaws_status qaws_curve_extremum_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 direction,
	qaws_scalar t_seed,
	qaws_vec3 const* direction_tangent,
	qaws_diff_views const* param_tangent,
	qaws_curve_extremum* out_value,
	qaws_curve_extremum* out_tangent);

qaws_status qaws_curve_extremum_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 direction,
	qaws_scalar t_seed,
	qaws_curve_extremum const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* direction_adjoint);

/* Inflection of the xy projection: x'y'' - y'x'' = 0. */
typedef struct qaws_curve_inflection
{
	qaws_scalar t;
	qaws_vec3 position;
} qaws_curve_inflection;

qaws_status qaws_curve_inflection_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t_seed,
	qaws_diff_views const* param_tangent,
	qaws_curve_inflection* out_value,
	qaws_curve_inflection* out_tangent);

qaws_status qaws_curve_inflection_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t_seed,
	qaws_curve_inflection const* adjoint,
	qaws_diff_views* param_adjoint);

#endif /* QAWS_DIFF_OPS_H */
