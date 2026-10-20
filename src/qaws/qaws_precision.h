/*
 * qaws_precision.h - the names of the current precision (C backend).
 *
 * No include guard: included again after QAWS_SCALAR_IS_FLOAT changes, it
 * switches every precision-bearing name to the other precision's form.
 * Qaws's library holds both: every scalar-bearing function and type exists
 * as name_f32 and name_f64 (qaws_names.h maps name to QAWS_SFX(name)).
 * QAWS_SCALAR_IS_FLOAT (default 1) picks which of the two the unsuffixed
 * names, qaws_scalar and the scalar macros below mean. qaws_both.h
 * declares both precisions in one translation unit.
 */

#if QAWS_BACKEND == QAWS_BACKEND_C

#ifndef QAWS_SCALAR_IS_FLOAT
#define QAWS_SCALAR_IS_FLOAT 1
#endif

#undef QAWS_SFX
#undef qaws_scalar
#undef QAWS_LITERAL
#undef QAWS_EPSILON
#undef QAWS_SQRT
#undef QAWS_FABS
#undef QAWS_FLOOR
#undef QAWS_CEIL
#undef QAWS_SIN
#undef QAWS_COS
#undef QAWS_ATAN2
#undef QAWS_POW
#undef QAWS_FMOD

#if QAWS_SCALAR_IS_FLOAT
#define QAWS_SFX(name) name##_f32
#define qaws_scalar qaws_f32
#define QAWS_LITERAL(x) x##f
#define QAWS_EPSILON 1e-6f
#define QAWS_SQRT(x) sqrtf(x)
#define QAWS_FABS(x) fabsf(x)
#define QAWS_FLOOR(x) floorf(x)
#define QAWS_CEIL(x) ceilf(x)
#define QAWS_SIN(x) sinf(x)
#define QAWS_COS(x) cosf(x)
#define QAWS_ATAN2(y, x) atan2f(y, x)
#define QAWS_POW(x, y) powf(x, y)
#define QAWS_FMOD(x, y) fmodf(x, y)
#else
#define QAWS_SFX(name) name##_f64
#define qaws_scalar qaws_f64
#define QAWS_LITERAL(x) x
#define QAWS_EPSILON 1e-12
#define QAWS_SQRT(x) sqrt(x)
#define QAWS_FABS(x) fabs(x)
#define QAWS_FLOOR(x) floor(x)
#define QAWS_CEIL(x) ceil(x)
#define QAWS_SIN(x) sin(x)
#define QAWS_COS(x) cos(x)
#define QAWS_ATAN2(y, x) atan2(y, x)
#define QAWS_POW(x, y) pow(x, y)
#define QAWS_FMOD(x, y) fmod(x, y)
#endif

#include "qaws_names.h"

/* the scalar helpers, once per precision */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_PRECISION_HELPERS_F32
#define QAWS_PRECISION_HELPERS_F32
#define QAWS_PRECISION_HELPERS_BODY
#endif
#else
#ifndef QAWS_PRECISION_HELPERS_F64
#define QAWS_PRECISION_HELPERS_F64
#define QAWS_PRECISION_HELPERS_BODY
#endif
#endif
#ifdef QAWS_PRECISION_HELPERS_BODY
#undef QAWS_PRECISION_HELPERS_BODY
QAWS_INLINE qaws_scalar qaws_min(qaws_scalar a, qaws_scalar b) { return (a < b) ? a : b; }
QAWS_INLINE qaws_scalar qaws_max(qaws_scalar a, qaws_scalar b) { return (a > b) ? a : b; }
QAWS_INLINE qaws_scalar qaws_clamp(qaws_scalar x, qaws_scalar lo, qaws_scalar hi) { return (x < lo) ? lo : (x > hi) ? hi : x; }
QAWS_INLINE qaws_scalar qaws_lerp(qaws_scalar a, qaws_scalar b, qaws_scalar t) { return (QAWS_ONE - t) * a + t * b; }
#endif

#endif /* QAWS_BACKEND == QAWS_BACKEND_C */
