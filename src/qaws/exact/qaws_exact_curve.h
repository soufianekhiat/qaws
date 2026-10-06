#ifndef QAWS_EXACT_CURVE_H
#define QAWS_EXACT_CURVE_H

/* Internal: exact curve layout and the exact rational evaluation. */

#include "../qaws_exact.h"
#include "qaws_exact_int.h"

#define QAWS_EXACT_MAX_DEGREE 16

/* One span: an integer homogeneous Bezier over the parameter lattice
   interval [a, b] (t = T 2^-param_shift). */
typedef struct qaws_exact_span
{
	unsigned int degree;
	int64_t a, b;
	qaws_exact_int* h;    /* (degree + 1) points of (w x, w y, [w z], w) */
} qaws_exact_span;

struct qaws_exact_curve
{
	int dimension;                 /* 2 or 3 */
	unsigned int span_count;
	qaws_exact_span* spans;
	int param_shift;               /* t = T 2^-param_shift */
	int space_exp2;
	qaws_exact_desc desc;
};

/*
 * C^(k)(T 2^-param_shift) exactly: num[c] / den per component c (den > 0 is
 * shared), k <= 3, T clamped to the domain. The span holding T is the one
 * with a <= T < b, the last one at the end of the domain. Values are in
 * lattice units (multiply by 2^space_exp2 for world units); derivatives
 * are with respect to the public parameter t.
 */
qaws_status qaws_exact_curve_eval_rational(qaws_exact_curve const* curve, int64_t T, unsigned int k, qaws_exact_int* num,
	qaws_exact_int* den);

#endif /* QAWS_EXACT_CURVE_H */
