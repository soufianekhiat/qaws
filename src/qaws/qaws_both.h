#ifndef QAWS_BOTH_H
#define QAWS_BOTH_H

/*
 * Both precisions in one translation unit: every public function and type
 * declared as name_f32 and name_f64 (qaws.h, once per precision). The
 * unsuffixed names, qaws_scalar and the scalar macros keep meaning the
 * caller's precision (QAWS_SCALAR_IS_FLOAT, default 1).
 */

#include "qaws_platform.h"

#if QAWS_SCALAR_IS_FLOAT
#define QAWS_BOTH_CALLER_IS_FLOAT 1
#else
#define QAWS_BOTH_CALLER_IS_FLOAT 0
#endif

#undef QAWS_SCALAR_IS_FLOAT
#define QAWS_SCALAR_IS_FLOAT 1
#include "qaws_precision.h"
#include "qaws.h"

#undef QAWS_SCALAR_IS_FLOAT
#define QAWS_SCALAR_IS_FLOAT 0
#include "qaws_precision.h"
#include "qaws.h"

#undef QAWS_SCALAR_IS_FLOAT
#if QAWS_BOTH_CALLER_IS_FLOAT
#define QAWS_SCALAR_IS_FLOAT 1
#else
#define QAWS_SCALAR_IS_FLOAT 0
#endif
#include "qaws_precision.h"

#endif /* QAWS_BOTH_H */
