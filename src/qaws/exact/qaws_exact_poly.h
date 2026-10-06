#ifndef QAWS_EXACT_POLY_H
#define QAWS_EXACT_POLY_H

/* Internal: integer polynomials (power basis) for the exact intersections. */

#include "qaws_exact_int.h"
#include "qaws_exact_roots.h"

#define QAWS_EXACT_POLY_COEF (QAWS_EXACT_ROOTS_MAX_DEGREE + 1)

typedef struct qaws_exact_poly
{
	unsigned int deg;
	qaws_exact_int c[QAWS_EXACT_POLY_COEF];
} qaws_exact_poly;

void qaws_exact_poly_zero(qaws_exact_poly* p, unsigned int deg);
void qaws_exact_poly_const(qaws_exact_poly* p, int64_t v);
int qaws_exact_poly_is_zero(qaws_exact_poly const* p);
/* Drops leading zero coefficients (the zero polynomial keeps degree 0). */
void qaws_exact_poly_trim(qaws_exact_poly* p);
/* Raises the stored degree to deg with zero coefficients. */
qaws_status qaws_exact_poly_pad(qaws_exact_poly* p, unsigned int deg);

qaws_status qaws_exact_poly_mul(qaws_exact_poly* r, qaws_exact_poly const* a, qaws_exact_poly const* b);
/* r += sign a */
qaws_status qaws_exact_poly_acc(qaws_exact_poly* r, qaws_exact_poly const* a, int sign);
/* r = k a (k an exact integer) */
qaws_status qaws_exact_poly_scale(qaws_exact_poly* r, qaws_exact_poly const* a, qaws_exact_int const* k);

int64_t qaws_exact_binom64(unsigned int n, unsigned int k);

/* Bernstein h_i (stride D, component c) on [0, 1] -> power basis. */
qaws_status qaws_exact_poly_from_bernstein(qaws_exact_int const* h, unsigned int n, unsigned int D, unsigned int c, qaws_exact_poly* out);
/* Power basis -> Bernstein of degree p->deg (times a positive integer). */
qaws_status qaws_exact_poly_to_bernstein(qaws_exact_poly const* p, qaws_exact_int* b);

/* Determinant of the n x n polynomial matrix M (row-major) over rows row..n-1 and columns outside mask. */
qaws_status qaws_exact_poly_det(qaws_exact_poly const* M, unsigned int n, unsigned int row, unsigned int mask, qaws_exact_poly* out);

/*
 * Primitive gcd of a and b (non-zero): Brown's modular algorithm (gcds
 * modulo 31-bit primes, Chinese remaindering, exact trial division as the
 * proof). A constant result proves a and b have no common root.
 * QAWS_STATUS_EXACT_RANGE_EXCEEDED past the integer budget.
 */
qaws_status qaws_exact_poly_gcd(qaws_exact_poly const* a, qaws_exact_poly const* b, qaws_exact_poly* out);

#endif /* QAWS_EXACT_POLY_H */
