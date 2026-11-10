#ifndef QAWS_INTERNAL_BEZIER_CHAIN_H
#define QAWS_INTERNAL_BEZIER_CHAIN_H

#include "../qaws_types.h"
#include "../qaws_status.h"
#include "../qaws_allocator.h"

/* A composite of `segments` 2D Bezier curves of one degree, each starting where the last ended: points holds
   2 (degree segments + 1) doubles, x then y. */
qaws_status qaws_internal_bezier_chain(qaws_curve** out_curve, double const* points, unsigned int degree, unsigned int segments,
	qaws_allocator const* allocator);

#endif /* QAWS_INTERNAL_BEZIER_CHAIN_H */
