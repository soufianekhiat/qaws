#ifndef QAWS_DIFF_SAMPLING_H
#define QAWS_DIFF_SAMPLING_H

#include "qaws_diff_types.h"

/*
 * Constant-speed (arc-length) sampling with exact first and second order
 * derivatives, forward and backward.
 *
 * A target asks for the point at arc length
 *
 *     sigma = distance + fraction * L_total
 *
 * from the start of the curve: fraction 0 is an absolute distance,
 * distance 0 a normalized one (i / (n - 1) gives n constant-speed samples),
 * distance -d with fraction 1 a distance from the end. The parameter t of the
 * sample solves L(t) = sigma, where L(t) is the length of [t_min, t]
 * integrated span by span with a composite Gauss-Legendre rule (every span,
 * and the last partial span, split into 4 pieces of `quadrature` points,
 * 0 = 8: near machine precision for the square root of the speed); the
 * forward value and every derivative use the same rule. Implicit differentiation of L(t) = sigma gives, along a parameter
 * direction (and a distance direction d'):
 *
 *     t'  = (sigma' - L'(t)) / |C'(t)|
 *     t'' = (sigma'' - L''(t) - 2 g(t) t' - (C'.C'' / |C'|) t'^2) / |C'(t)|
 *
 * with L', L'' the directional derivatives of the length integral and
 * g = C'.dC'/|C'| its integrand's rate at t; the sample point follows by
 * the chain rule. Adjoints use the same relation: a sample adjoint
 * (t_bar, p_bar) pulls back as J^T p_bar plus
 * lambda (fraction grad L_total - grad L(t)), lambda = (t_bar + p_bar.C') / |C'|,
 * accumulated for the whole batch in one quadrature pass.
 *
 *   tangent  : value, first and (optionally) second directional derivative
 *   adjoint  : gradient (+=) into parameter views and distance adjoints
 *   hvp      : Hessian-vector product (+=) of sum(t_bar t + p_bar . p);
 *              direct for families linear in their fields (QAWS_CAP_LINEAR),
 *              QAWS_STATUS_UNSUPPORTED_OPERATION otherwise
 *
 * Knots move the span boundaries of the quadrature and are not supported
 * as parameters here (QAWS_STATUS_UNSUPPORTED_OPERATION). 2D curves work
 * too (z of the outputs is zero). Samples where |C'| vanishes are reported
 * as QAWS_DIFF_ILL_CONDITIONED through ctx->report.
 */

typedef struct qaws_arc_length_target
{
	qaws_scalar distance;    /* arc length from the start */
	qaws_scalar fraction;    /* plus this fraction of the total length */
} qaws_arc_length_target;

typedef struct qaws_arc_length_sample
{
	qaws_scalar t;           /* curve parameter */
	qaws_vec3 position;
} qaws_arc_length_sample;

/* Solves every target and returns the samples with their tangents along
   param_tangent (may be NULL) and distance_tangent (per target, may be
   NULL). out_tangent and out_tangent2 may be NULL; out_total_length (may be
   NULL) receives L_total. */
qaws_status qaws_curve_arc_length_sample_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_arc_length_target const* targets,
	qaws_scalar const* distance_tangent,
	unsigned int count,
	unsigned int quadrature,
	qaws_diff_views const* param_tangent,
	qaws_arc_length_sample* out_value,
	qaws_arc_length_sample* out_tangent,
	qaws_arc_length_sample* out_tangent2,
	qaws_scalar* out_total_length);

/* Accumulates (+=) the pullback of per-sample adjoints into param_adjoint
   (may be NULL) and distance_adjoint (per target, may be NULL). */
qaws_status qaws_curve_arc_length_sample_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_arc_length_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_arc_length_sample const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* distance_adjoint);

/* Accumulates (+=) H direction, H the Hessian with respect to the
   parameters of sum_i (adjoint_i.t t_i + adjoint_i.position . p_i). */
qaws_status qaws_curve_arc_length_sample_hvp(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_arc_length_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_arc_length_sample const* adjoint,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);

#endif /* QAWS_DIFF_SAMPLING_H */
