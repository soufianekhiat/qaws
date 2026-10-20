#ifndef QAWS_CURVE_H
#define QAWS_CURVE_H

#include "qaws_types.h"
#include "qaws_status.h"

/* Thread-safe: curve is not accessed after destruction. */
void qaws_curve_destroy(qaws_curve* curve);

qaws_status qaws_curve_reverse(
	qaws_curve** out_reversed,
	qaws_curve const* curve);

#endif /* QAWS_CURVE_H */
