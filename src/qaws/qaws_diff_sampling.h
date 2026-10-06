#ifndef QAWS_DIFF_SAMPLING_H
#define QAWS_DIFF_SAMPLING_H

#include "qaws_diff_types.h"
#include "qaws_surface_types.h"

/*
 * Inverse-CDF sampling of curves with exact first and second order
 * derivatives, forward and backward.
 *
 * A measure along the curve, M(t) = integral from t_min to t of m(tau) dtau
 * with m = rho |C'|, plays the role of a CDF:
 *
 *   QAWS_MEASURE_ARC_LENGTH  rho = 1: constant-speed (arc-length) sampling
 *   QAWS_MEASURE_CURVATURE   rho = sqrt(floor^2 + kappa^2): samples gather
 *                            where the curve bends, floor sets the share of
 *                            plain arc length (smooth, unlike kappa itself)
 *   QAWS_MEASURE_DENSITY     rho = density(C(t)), a positive user field in
 *                            space with its gradient and Hessian (importance
 *                            maps, distance fields, ...)
 *
 * A target asks for the point at measure sigma = distance + fraction * M_total:
 * fraction 0 is an absolute measure, distance 0 a normalized one
 * (i / (n - 1) gives n samples equidistributed in the measure), distance -d
 * with fraction 1 a measure from the end. The sample parameter t solves
 * M(t) = sigma. M is integrated span by span with a composite Gauss-Legendre
 * rule (every span and the last partial span split into 8 pieces of
 * `quadrature` points, 0 = 8); the forward value and every derivative use
 * that same rule. Implicit differentiation of M(t) = sigma gives, along a
 * parameter direction (and a distance direction d'):
 *
 *     t'  = (sigma' - M'(t)) / m(t)
 *     t'' = (sigma'' - M''(t) - 2 m_e(t) t' - m_t(t) t'^2) / m(t)
 *
 * with M', M'' the directional derivatives of the integral, m_e the rate of
 * the integrand along the direction at fixed t and m_t its t-derivative; the
 * sample point follows by the chain rule. A sample adjoint (t_bar, p_bar)
 * pulls back as J^T p_bar plus lambda (fraction grad M_total - grad M(t)),
 * lambda = (t_bar + p_bar . C') / m, accumulated for the whole batch in one
 * quadrature pass.
 *
 *   tangent  : values, first and (optionally) second directional derivatives
 *   adjoint  : gradient (+=) into parameter views and distance adjoints
 *   hvp      : Hessian-vector product (+=) of sum(t_bar t + p_bar . p);
 *              direct for families linear in their fields (QAWS_CAP_LINEAR),
 *              QAWS_STATUS_UNSUPPORTED_OPERATION otherwise
 *
 * measure may be NULL (arc length). Knots move the span boundaries of the
 * quadrature and are not supported as parameters here
 * (QAWS_STATUS_UNSUPPORTED_OPERATION). 2D curves work too (z of the outputs
 * is zero; a density then reads z = 0). Samples where m vanishes are
 * reported as QAWS_DIFF_ILL_CONDITIONED through ctx->report.
 */

typedef enum qaws_sample_measure
{
	QAWS_MEASURE_ARC_LENGTH = 0,
	QAWS_MEASURE_CURVATURE,
	QAWS_MEASURE_DENSITY
} qaws_sample_measure;

/* Positive density at a point, with its gradient and Hessian (xx, xy, xz,
   yy, yz, zz); the outputs are always given. */
typedef qaws_scalar (*qaws_density_fn)(qaws_vec3 position, void* user_data, qaws_vec3* out_gradient,
	qaws_scalar* out_hessian);

typedef struct qaws_sample_measure_desc
{
	qaws_sample_measure kind;
	qaws_scalar curvature_floor;     /* QAWS_MEASURE_CURVATURE: rho = sqrt(floor^2 + kappa^2) */
	qaws_density_fn density;         /* QAWS_MEASURE_DENSITY */
	void* density_user_data;
} qaws_sample_measure_desc;

typedef struct qaws_cdf_target
{
	qaws_scalar distance;    /* measure from the start */
	qaws_scalar fraction;    /* plus this fraction of the total measure */
} qaws_cdf_target;

typedef struct qaws_cdf_sample
{
	qaws_scalar t;           /* curve parameter */
	qaws_vec3 position;
} qaws_cdf_sample;

/* Solves every target and returns the samples with their tangents along
   param_tangent (may be NULL) and the per-target distance rates
   distance_tangent / distance_tangent2 (first and second derivative of the
   target distance along the same path, may be NULL; time traversal gives
   them, see qaws_traversal_cdf_targets). out_tangent and out_tangent2 may be
   NULL; out_total (may be NULL) receives M_total. */
qaws_status qaws_curve_cdf_sample_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	qaws_scalar const* distance_tangent,
	qaws_scalar const* distance_tangent2,
	unsigned int count,
	unsigned int quadrature,
	qaws_diff_views const* param_tangent,
	qaws_cdf_sample* out_value,
	qaws_cdf_sample* out_tangent,
	qaws_cdf_sample* out_tangent2,
	qaws_scalar* out_total);

/* Accumulates (+=) the pullback of per-sample adjoints into param_adjoint
   (may be NULL) and distance_adjoint (per target, may be NULL). */
qaws_status qaws_curve_cdf_sample_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_cdf_sample const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* distance_adjoint);

/* Accumulates (+=) H direction, H the Hessian with respect to the
   parameters of sum_i (adjoint_i.t t_i + adjoint_i.position . p_i). */
qaws_status qaws_curve_cdf_sample_hvp(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_cdf_sample const* adjoint,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);

/* Traversal inputs as arc-length targets: the distance of each input (time
   through easing and motion profile, or arc length) after the wrap mode,
   written as distance + fraction * L_total so that the derivatives of the
   total length are carried (loop: fraction -k, ping-pong: (s, -2k) going
   out and (-s, 2k + 2) coming back, clamp: (0, 1) past the end).
   out_rate / out_rate2 (may be NULL) receive the first and second
   derivatives of the target distance with respect to the input (speed and
   acceleration of the profile through the easing); pass them as
   distance_tangent / distance_tangent2 scaled by the input rate. Sampling
   then goes through qaws_curve_cdf_sample_* with a NULL measure (exact
   arc length, where the traversal itself reads a table). Parameter-mode and
   multi-curve traversals are refused (QAWS_STATUS_UNSUPPORTED_OPERATION). */
qaws_status qaws_traversal_cdf_targets(
	qaws_traversal const* traversal,
	qaws_scalar const* inputs,
	unsigned int count,
	qaws_cdf_target* out_targets,
	qaws_scalar* out_rate,
	qaws_scalar* out_rate2);

/* ===================================================================
 * Surfaces: inverse-CDF warping of the unit square
 *
 * The measure w = rho(S) |S_u x S_v| (QAWS_MEASURE_AREA: rho = 1, or
 * QAWS_MEASURE_DENSITY) over the parameter domain gives a marginal CDF in u
 * and a conditional CDF in v. A point xi = (xi_u, xi_v) of the unit square
 * (stratified, blue noise, low discrepancy ...) maps to the sample (u, v):
 *
 *   A(u)    = xi_u A_total,    A(u) = int_{u0}^{u} int_{v0}^{v1} w
 *   B(v; u) = xi_v B(v1; u),   B(v; u) = int_{v0}^{v} w(u, .)
 *
 * so that samples are distributed by the measure and stratification of xi
 * carries over. Both integrals use a grid of `cells` x `cells` cells
 * (0 = 6) with `quadrature` Gauss points per cell and direction (0 = 8).
 * Derivatives are exact for these discrete equations: the forward pass
 * solves them in dual numbers (the unknown and the quadrature nodes that
 * move with it carry first and second order rates), the adjoint pulls back
 * through the conditional then the marginal equation, and the HVP polarizes
 * the second order forward pass over the parameters of out_hv (about one
 * pass per parameter; works for rational families too).
 * Requires QAWS_CAP_TANGENT2. QAWS_MEASURE_CURVATURE is not defined on
 * surfaces (QAWS_STATUS_INVALID_ARGUMENT); knots are refused as parameters.
 * =================================================================== */

#define QAWS_MEASURE_AREA QAWS_MEASURE_ARC_LENGTH

typedef struct qaws_surface_cdf_sample
{
	qaws_scalar u;
	qaws_scalar v;
	qaws_vec3 position;
} qaws_surface_cdf_sample;

/* Warps count points xi (2 scalars each) and returns the samples with their
   tangents along param_tangent (may be NULL) and xi_tangent (2 per point,
   may be NULL). out_tangent / out_tangent2 / out_total may be NULL. */
qaws_status qaws_surface_cdf_sample_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi,
	qaws_scalar const* xi_tangent,
	unsigned int count,
	unsigned int cells,
	unsigned int quadrature,
	qaws_diff_views const* param_tangent,
	qaws_surface_cdf_sample* out_value,
	qaws_surface_cdf_sample* out_tangent,
	qaws_surface_cdf_sample* out_tangent2,
	qaws_scalar* out_total);

/* Accumulates (+=) the pullback of per-sample adjoints into param_adjoint
   (may be NULL) and xi_adjoint (2 per point, may be NULL). */
qaws_status qaws_surface_cdf_sample_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi,
	unsigned int count,
	unsigned int cells,
	unsigned int quadrature,
	qaws_surface_cdf_sample const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* xi_adjoint);

/* Accumulates (+=) H direction for the parameters present in out_hv (top
   level fields; views with children are refused), H the Hessian of
   sum_i (adjoint_i.u u_i + adjoint_i.v v_i + adjoint_i.position . p_i). */
qaws_status qaws_surface_cdf_sample_hvp(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi,
	unsigned int count,
	unsigned int cells,
	unsigned int quadrature,
	qaws_surface_cdf_sample const* adjoint,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);

#endif /* QAWS_DIFF_SAMPLING_H */
