#ifndef QAWS_INTERNAL_OFFSET_H
#define QAWS_INTERNAL_OFFSET_H

#include "../qaws_offset.h"

/* The offset chain of one curve at a constant delta (to the right of
   travel), before any union: G1 cubics along the offset, the evolute through
   its backward runs, joins at corners (round). The chain is one cubic
   B-spline when the curve is smooth (straight pieces are lines); out_curves
   receives its first curve_capacity curves, *out_curve_count their number.
   cusps (optional) receives up to cusp_capacity parameters of the offset's
   cusps; *cusp_count their number. */
qaws_status qaws_internal_offset_curve(qaws_curve** out_curves, unsigned int* out_curve_count, double* cusps,
	unsigned int* cusp_count, qaws_curve const* curve, unsigned int curve_capacity,
	unsigned int cusp_capacity, double delta, double tolerance);

#endif /* QAWS_INTERNAL_OFFSET_H */
