#ifndef QAWS_INTERNAL_BEZIER_CHAIN_H
#define QAWS_INTERNAL_BEZIER_CHAIN_H

#include "../qaws_types.h"
#include "../qaws_status.h"
#include "../qaws_allocator.h"

/* A composite of `segments` 2D Bezier curves of one degree, each starting where the last ended: points holds
   2 (degree segments + 1) doubles, x then y. */
qaws_status qaws_internal_bezier_chain(qaws_curve** out_curve, double const* points, unsigned int degree, unsigned int segments,
	qaws_allocator const* allocator);

/* Curvature-matching spline solvers. points: 2 n doubles, x then y, no two consecutive equal (nor the last and first when
   closed), scaled into about [-1, 1], n >= 3. On success *out_points (2 (3 segments + 1) doubles, a chain of cubics, from
   qaws_internal_malloc, the caller frees it) starts at the first point and, closed, ends at it. */
qaws_status qaws_internal_interpolate_cubic_g2(double** out_points, unsigned int* out_segments, double const* points, unsigned int n, int closed);
qaws_status qaws_internal_interpolate_s9(double** out_points, unsigned int* out_segments, double const* points, unsigned int n, int closed);

#endif /* QAWS_INTERNAL_BEZIER_CHAIN_H */
