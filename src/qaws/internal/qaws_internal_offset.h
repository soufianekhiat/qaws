#ifndef QAWS_INTERNAL_OFFSET_H
#define QAWS_INTERNAL_OFFSET_H

#include "../qaws_offset.h"

/* The offset chain of one curve at a constant delta (to the right of
   travel), before any union: G1 cubics along the offset, the evolute through
   its backward runs, joins at corners (round). The chain is one cubic
   B-spline when the curve is smooth; *out_chain receives the first curve of
   the chain. cusps (optional) receives up to cusp_capacity parameters of the
   offset's cusps; *cusp_count their number. */
qaws_status qaws_internal_offset_curve(qaws_curve const* curve, double delta, double tolerance,
	qaws_curve** out_chain, unsigned int* out_curve_count,
	double* cusps, unsigned int cusp_capacity, unsigned int* cusp_count);

#endif /* QAWS_INTERNAL_OFFSET_H */
