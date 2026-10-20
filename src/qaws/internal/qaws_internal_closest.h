#ifndef QAWS_INTERNAL_CLOSEST_H
#define QAWS_INTERNAL_CLOSEST_H

#include "qaws_internal_types.h"
#include "../qaws_surface_types.h"

/* Span s of a borrowed spline curve as the B-spline points of its window (degree + 1 homogeneous points, dim + 1 doubles
   each, the weight last) and the 2 * degree knots around it; 0 for the other borrowed families. */
int qaws_internal_view_span_window(unsigned int* degree, double* window, double* knots, qaws_curve const* curve, unsigned int s);

/* Certifies a closest point on a curve with Bezier spans (Bezier, rational Bezier, B-spline, NURBS, borrowed splines):
   branch and bound from the parameter *io_t, which it replaces when it finds a closer point. The global minimum distance
   lies in [*out_distance - *out_gap, *out_distance]. QAWS_STATUS_UNSUPPORTED_OPERATION for the other families (nothing
   written). q: dim coordinates. */
qaws_status qaws_internal_closest_certify(double* out_distance, double* out_gap, qaws_scalar* io_t, qaws_curve const* curve,
	double const* q);

/* A borrowed surface's net: degrees and point counts (u, v), its knots; 0 for other surfaces. Point (i, j) homogeneous
   (w x, w y, w z, w). */
int qaws_internal_surface_view_net(qaws_scalar const** u_knots, qaws_scalar const** v_knots, unsigned int* degrees, unsigned int* counts,
	qaws_surface const* surface);
void qaws_internal_surface_view_point(double* h, qaws_surface const* surface, unsigned int i, unsigned int j);

/* The same for surfaces with Bezier patches (Bezier, B-spline, NURBS, borrowed; degrees up to 15) from (*io_u, *io_v). */
qaws_status qaws_internal_surface_closest_certify(double* out_distance, double* out_gap, qaws_scalar* io_u, qaws_scalar* io_v,
	qaws_surface const* surface, double const* q);

#endif /* QAWS_INTERNAL_CLOSEST_H */
