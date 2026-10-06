#ifndef QAWS_EXACT_H
#define QAWS_EXACT_H

/*
 * Exact / certified geometry (C only, built with QAWS_ENABLE_EXACT).
 *
 * Predicates on raw doubles: a finite double is an exact dyadic rational
 * m 2^e, so the sign of a polynomial in double inputs is a well defined
 * exact quantity. Each predicate evaluates in f64 first with a proven error
 * bound; when the bound does not prove the sign, it decodes the inputs
 * exactly and evaluates in multi-limb integers. The f64 filter never
 * decides a zero. Results are certified: never a wrong sign.
 *
 * QAWS_STATUS_EXACT_RANGE_EXCEEDED: the inputs span more binary exponents
 * than the integer budget covers (about 2^900 between the smallest and the
 * largest magnitude); QAWS_STATUS_INVALID_ARGUMENT for non-finite inputs.
 */

#include "qaws_status.h"

typedef enum qaws_exact_sign
{
	QAWS_EXACT_NEGATIVE = -1,
	QAWS_EXACT_ZERO = 0,
	QAWS_EXACT_POSITIVE = 1
} qaws_exact_sign;

/* How a predicate decided (optional output). */
typedef enum qaws_exact_path
{
	QAWS_EXACT_PATH_FILTER = 0,   /* the f64 evaluation and its error bound */
	QAWS_EXACT_PATH_EXACT         /* exact integer evaluation */
} qaws_exact_path;

/* Sign of (a - c) x (b - c): positive when a, b, c turn counterclockwise. */
qaws_status qaws_exact_orient2d(double const a[2], double const b[2], double const c[2], qaws_exact_sign* out_sign,
	qaws_exact_path* out_path);

/* Sign of the volume (a - d) . ((b - d) x (c - d)). */
qaws_status qaws_exact_orient3d(double const a[3], double const b[3], double const c[3], double const d[3],
	qaws_exact_sign* out_sign, qaws_exact_path* out_path);

/* Sign of a / b - c / d (b, d non-zero). */
qaws_status qaws_exact_compare_ratio(double a, double b, double c, double d, qaws_exact_sign* out_sign, qaws_exact_path* out_path);

#endif /* QAWS_EXACT_H */
