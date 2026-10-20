#ifndef QAWS_INTERNAL_SPAN_H
#define QAWS_INTERNAL_SPAN_H

#include "qaws_internal_types.h"

/* Find the span index for a global parameter value.
   Also writes the local parameter within that span to out_local_t.
   Returns the span index (clamped to [0, span_count-1]). */
unsigned int qaws_internal_find_span(
	qaws_scalar* out_local_t,
	qaws_curve const* curve,
	qaws_scalar parameter);

#endif /* QAWS_INTERNAL_SPAN_H */
