#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_CURVE_H_F32
#define QAWS_CURVE_H_F32
#define QAWS_CURVE_H_BODY
#endif
#else
#ifndef QAWS_CURVE_H_F64
#define QAWS_CURVE_H_F64
#define QAWS_CURVE_H_BODY
#endif
#endif
#ifdef QAWS_CURVE_H_BODY
#undef QAWS_CURVE_H_BODY
#define QAWS_CURVE_H

#include "qaws_types.h"
#include "qaws_status.h"

/* Thread-safe: curve is not accessed after destruction. */
void qaws_curve_destroy(qaws_curve* curve);

qaws_status qaws_curve_reverse(
	qaws_curve** out_reversed,
	qaws_curve const* curve);

#endif /* QAWS_CURVE_H */
