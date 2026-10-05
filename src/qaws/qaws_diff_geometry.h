#ifndef QAWS_DIFF_GEOMETRY_H
#define QAWS_DIFF_GEOMETRY_H

#include "qaws_diff_types.h"

/*
 * Differentiable local differential geometry.
 *
 * Every quantity is a closed-form function of a spatial jet, so the
 * rules compose with the curve and surface rules:
 *
 *   parameters --eval_tangent--> jet tangent --geometry_eval--> quantity tangent
 *   quantity adjoint --geometry_adjoint--> jet adjoint --eval_adjoint--> parameters
 *
 * Validity:
 *   QAWS_DIFF_VALID            all quantities are differentiable here
 *   QAWS_DIFF_ILL_CONDITIONED  a quantity is undefined or not differentiable
 *                              (curve: zero curvature leaves normal, binormal
 *                              and torsion undefined; surface: umbilic point
 *                              for principal curvatures). Those quantities
 *                              and their derivatives are returned as zero.
 *   QAWS_DIFF_INVALID          singular jet (zero speed, degenerate normal):
 *                              everything is returned as zero.
 */

typedef struct qaws_curve_geometry_2d
{
	qaws_vec2 tangent;       /* unit tangent */
	qaws_vec2 normal;        /* tangent rotated by +90 degrees */
	qaws_scalar speed;       /* |C'| */
	qaws_scalar curvature;   /* signed: (x'y'' - y'x'') / |C'|^3 */
} qaws_curve_geometry_2d;

typedef struct qaws_curve_geometry_3d
{
	qaws_vec3 tangent;       /* T */
	qaws_vec3 normal;        /* N = B x T */
	qaws_vec3 binormal;      /* B = (C' x C'') / |C' x C''| */
	qaws_scalar speed;
	qaws_scalar curvature;   /* |C' x C''| / |C'|^3 */
	qaws_scalar torsion;     /* (C' x C'') . C''' / |C' x C''|^2 */
} qaws_curve_geometry_3d;

typedef struct qaws_surface_geometry
{
	qaws_vec3 normal;        /* (Su x Sv) / |Su x Sv| */
	qaws_scalar E, F, G;     /* first fundamental form */
	qaws_scalar L, M, N;     /* second fundamental form */
	qaws_scalar gaussian;    /* (LN - M^2) / (EG - F^2) */
	qaws_scalar mean;        /* (EN + GL - 2FM) / (2 (EG - F^2)) */
	qaws_scalar kappa1;      /* H + sqrt(H^2 - K) */
	qaws_scalar kappa2;      /* H - sqrt(H^2 - K) */
} qaws_surface_geometry;

/*
 * Forward evaluation. Required primal channels: curves D1, D2 (and D3 for
 * torsion in 3D); surfaces U, V, UU, UV, VV. tangent / tangent2 inputs and
 * outputs are optional (NULL): a tangent output needs a tangent input, a
 * tangent2 output needs both.
 */
qaws_status qaws_curve_geometry_eval_2d(
	qaws_curve_jet_2d const* primal,
	qaws_curve_jet_2d const* tangent,
	qaws_curve_jet_2d const* tangent2,
	qaws_curve_geometry_2d* out_value,
	qaws_curve_geometry_2d* out_tangent,
	qaws_curve_geometry_2d* out_tangent2,
	qaws_diff_validity* out_validity);

qaws_status qaws_curve_geometry_eval_3d(
	qaws_curve_jet_3d const* primal,
	qaws_curve_jet_3d const* tangent,
	qaws_curve_jet_3d const* tangent2,
	qaws_curve_geometry_3d* out_value,
	qaws_curve_geometry_3d* out_tangent,
	qaws_curve_geometry_3d* out_tangent2,
	qaws_diff_validity* out_validity);

qaws_status qaws_surface_geometry_eval(
	qaws_surface_jet const* primal,
	qaws_surface_jet const* tangent,
	qaws_surface_jet const* tangent2,
	qaws_surface_geometry* out_value,
	qaws_surface_geometry* out_tangent,
	qaws_surface_geometry* out_tangent2,
	qaws_diff_validity* out_validity);

/*
 * Reverse propagation: accumulates (+=) into the jet adjoint the pullback
 * of the quantity adjoint. Only the jet channels the quantities depend on
 * are touched; the jet adjoint channels mask is extended accordingly.
 */
qaws_status qaws_curve_geometry_adjoint_2d(
	qaws_curve_jet_2d const* primal,
	qaws_curve_geometry_2d const* adjoint,
	qaws_curve_jet_2d* inout_jet_adjoint,
	qaws_diff_validity* out_validity);

qaws_status qaws_curve_geometry_adjoint_3d(
	qaws_curve_jet_3d const* primal,
	qaws_curve_geometry_3d const* adjoint,
	qaws_curve_jet_3d* inout_jet_adjoint,
	qaws_diff_validity* out_validity);

qaws_status qaws_surface_geometry_adjoint(
	qaws_surface_jet const* primal,
	qaws_surface_geometry const* adjoint,
	qaws_surface_jet* inout_jet_adjoint,
	qaws_diff_validity* out_validity);

#endif /* QAWS_DIFF_GEOMETRY_H */
