#ifndef QAWS_INTERNAL_BASIS_H
#define QAWS_INTERNAL_BASIS_H

#include "../qaws_types.h"

/* The degrees these kernels hold on the stack: B-spline basis functions and their derivatives (left, right and ndu arrays
   of order 32), de Casteljau (192 scalars: 48 homogeneous 3D points). Constructors refuse more. */
#define QAWS_INTERNAL_BASIS_MAX_DEGREE 31
#define QAWS_INTERNAL_DECASTELJAU_MAX_DEGREE 47

/* De Casteljau evaluation for Bezier curves.
   control_points: flat array of dim_count scalars per point.
   degree: polynomial degree (point_count - 1).
   t: parameter in [0,1].
   out_point: receives dim_count scalars. */
void qaws_internal_decasteljau(
	qaws_scalar* out_point,
	qaws_scalar const* control_points,
	unsigned int dim_count,
	unsigned int degree,
	qaws_scalar t);

/* Bezier derivative control points.
   Given degree-N control points, writes degree-(N-1) differenced control points.
   out_deriv_points must hold degree * dim_count scalars. */
void qaws_internal_bezier_derivative_points(
	qaws_scalar* out_deriv_points,
	qaws_scalar const* control_points,
	unsigned int dim_count,
	unsigned int degree);

/* Cubic Hermite basis values and derivatives.
   t: local parameter in [0,1].
   out_h: receives 4 basis values  [h00, h10, h01, h11].
   out_dh: if non-NULL, receives 4 first derivative values.
   out_d2h: if non-NULL, receives 4 second derivative values.
   out_d3h: if non-NULL, receives 4 third derivative values. */
void qaws_internal_hermite_basis(
	qaws_scalar* out_h,
	qaws_scalar* out_dh,
	qaws_scalar* out_d2h,
	qaws_scalar* out_d3h,
	qaws_scalar t);

/* B-spline basis functions N_{i,p}(t).
   knots: full knot vector.
   knot_count: number of knots.
   degree: spline degree.
   span: knot span index such that knots[span] <= t < knots[span+1].
   t: parameter value.
   out_basis: receives (degree+1) basis values for indices [span-degree .. span]. */
void qaws_internal_bspline_basis(
	qaws_scalar* out_basis,
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int span,
	qaws_scalar t);

/* B-spline basis functions and their derivatives up to order k.
   out_ders: (k+1) x (degree+1) row-major array.
   out_ders[0] = basis values, out_ders[1] = first derivs, etc. */
void qaws_internal_bspline_basis_derivs(
	qaws_scalar* out_ders,
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int span,
	qaws_scalar t,
	unsigned int k);

/* Find knot span index such that knots[span] <= t < knots[span+1].
   Returns the span index. Clamps to valid range. */
unsigned int qaws_internal_find_knot_span(
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int num_control_points,
	qaws_scalar t);

/* The basis functions and derivatives of qaws_internal_bspline_basis_derivs for several parameters at once, one per
   SIMD lane: lane i at t[i] in knot span span[i], up to order k (at most 3), with the same operations in the same order
   (each lane equal to the one-parameter call bit for bit). qaws_internal_bspline_basis_lanes(degree): the lanes filled,
   0 where this does not apply (a degree past 5, no SIMD for the precision). t holds that many parameters, 32-byte
   aligned; out_ders[(d * (degree + 1) + j) * QAWS_INTERNAL_BASIS_LANES_MAX + i], 32-byte aligned. */
#define QAWS_INTERNAL_BASIS_LANES_MAX 8
unsigned int qaws_internal_bspline_basis_lanes(unsigned int degree);
void qaws_internal_bspline_basis_derivs_lanes(
	qaws_scalar* out_ders,
	qaws_scalar const* knots,
	unsigned int const* span,
	qaws_scalar const* t,
	unsigned int degree,
	unsigned int k);

#endif /* QAWS_INTERNAL_BASIS_H */
