#ifndef QAWS_INTERNAL_CLOSEST_H
#define QAWS_INTERNAL_CLOSEST_H

#include "qaws_internal_types.h"

/* Span s of a borrowed spline curve as the B-spline points of its window (degree + 1 homogeneous points, dim + 1 doubles
   each, the weight last) and the 2 * degree knots around it; 0 for the other borrowed families. */
int qaws_internal_view_span_window(qaws_curve const* curve, unsigned int s, unsigned int* degree, double* window, double* knots);

/* Certifies a closest point on a curve with Bezier spans (Bezier, rational Bezier, B-spline, NURBS, borrowed splines):
   branch and bound from the parameter *io_t, which it replaces when it finds a closer point. The global minimum distance
   lies in [*out_distance - *out_gap, *out_distance]. QAWS_STATUS_UNSUPPORTED_OPERATION for the other families (nothing
   written). q: dim coordinates. */
qaws_status qaws_internal_closest_certify(qaws_curve const* curve, double const* q, qaws_scalar* io_t, double* out_distance,
	double* out_gap);

#endif /* QAWS_INTERNAL_CLOSEST_H */
