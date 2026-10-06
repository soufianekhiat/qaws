#ifndef QAWS_EXACT_CURVE_H
#define QAWS_EXACT_CURVE_H

/* Internal: exact curve layout and the exact rational evaluation. */

#include "../qaws_exact.h"
#include "qaws_exact_int.h"

#define QAWS_EXACT_MAX_DEGREE 16

struct qaws_exact_curve
{
	int dimension;                 /* 2 or 3 */
	unsigned int degree;
	unsigned int param_bits;
	int space_exp2;
	qaws_exact_desc desc;
	/* (degree + 1) homogeneous points (w x, w y, [w z], w): (dimension + 1) each */
	qaws_exact_int* h;
};

/*
 * C^(k)(T 2^-param_bits) exactly: num[c] / den per component c (den > 0 is
 * shared), k <= 3, 0 <= T <= 2^param_bits. Values are in lattice units
 * (multiply by 2^space_exp2 for world units).
 */
qaws_status qaws_exact_curve_eval_rational(qaws_exact_curve const* curve, int64_t T, unsigned int k, qaws_exact_int* num,
	qaws_exact_int* den);

#endif /* QAWS_EXACT_CURVE_H */
