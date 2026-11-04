#ifndef QAWS_INTERNAL_CLOSEST_H
#define QAWS_INTERNAL_CLOSEST_H

#include "qaws_internal_types.h"
#include "../qaws_surface_types.h"

/* Span s of a borrowed spline curve as the B-spline points of its window (degree + 1 homogeneous points, dim + 1 doubles
   each, the weight last) and the 2 * degree knots around it; 0 for the other borrowed families. */
int qaws_internal_view_span_window(unsigned int* degree, double* window, double* knots, qaws_curve const* curve, unsigned int s);

/* The same for a curve with Bezier spans (Bezier, rational Bezier, B-spline, NURBS, borrowed splines): span s's window,
   its degree and knots; 0 for the other families. */
int qaws_internal_curve_window(qaws_curve const* c, unsigned int s, unsigned int* p, double* win, double* uk);

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

/* A surface's control net for its Bezier patches (Bezier, B-spline, NURBS, borrowed; degrees under
   QAWS_INTERNAL_NET_MAX): degrees, counts, knots (u, v; a Bezier surface's the clamped knots of [0, 1]); its points
   through qaws_internal_surface_net_point, homogeneous (w x, w y, w z, w), u index major. */
#define QAWS_INTERNAL_NET_MAX 16
typedef struct qaws_internal_surface_net
{
	qaws_surface const* s;
	unsigned int p[2], n[2];
	qaws_scalar const* k[2];
	qaws_scalar const* cp;              /* 3 scalars per point, u-major (NULL: a borrowed net) */
	qaws_scalar const* w;
	qaws_scalar bez[2][2 * QAWS_INTERNAL_NET_MAX];
} qaws_internal_surface_net;
/* 0 for other surfaces or degrees too high */
int qaws_internal_surface_net_init(qaws_surface const* s, qaws_internal_surface_net* m);
void qaws_internal_surface_net_point(qaws_internal_surface_net const* m, unsigned int i, unsigned int j, double* h);

/* The same for surfaces with Bezier patches (Bezier, B-spline, NURBS, borrowed; degrees up to 15) from (*io_u, *io_v). */
qaws_status qaws_internal_surface_closest_certify(double* out_distance, double* out_gap, qaws_scalar* io_u, qaws_scalar* io_v,
	qaws_surface const* surface, double const* q);

#endif /* QAWS_INTERNAL_CLOSEST_H */
