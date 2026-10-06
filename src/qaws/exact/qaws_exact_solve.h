#ifndef QAWS_EXACT_SOLVE_H
#define QAWS_EXACT_SOLVE_H

/*
 * Internal: certified solutions of two polynomial equations P(s, r) =
 * Q(s, r) = 0 eliminated through a Bezout matrix in s (entries polynomials
 * in r): R(r) = det, s = C01 / C00 at a root (the kernel (1, s, s^2, ...)).
 */

#include "qaws_exact_poly.h"
#include "qaws_exact_roots.h"

#define QAWS_EXACT_SOLVE_MAX_ROOTS (QAWS_EXACT_ROOTS_MAX_DEGREE + 2)

/* A solution in local parameters: r (a root of R) and s = C01 / C00 there. */
typedef struct qaws_exact_local_hit
{
	qaws_exact_root r;
	int s_exact;
	qaws_exact_int s_num, s_den;   /* s exactly (s_den > 0), when s_exact */
	double s_lo, s_hi;             /* otherwise its enclosure */
} qaws_exact_local_hit;

enum
{
	QAWS_EXACT_CON_SIGN_OF_C00 = 1,   /* H~(r) has the sign of C00(r): (s - r) > 0 for H~ = C01 - r C00 */
	QAWS_EXACT_CON_ZERO = 2           /* H~(r) == 0 (an exact zero test) */
};

typedef struct qaws_exact_constraint
{
	int kind;
	qaws_exact_poly const* H;
} qaws_exact_constraint;

/* Bezout matrix of P, Q (degree n, integer coefficients), n x n, symmetric. */
qaws_status qaws_exact_bezout(qaws_exact_poly const* P, qaws_exact_poly const* Q, unsigned int n, qaws_exact_int* B);

/* The same with polynomial coefficients: P[i], Q[i] the coefficients of s^i (degree k in s), E k x k. */
qaws_status qaws_exact_bezout_poly(qaws_exact_poly const* P, qaws_exact_poly const* Q, unsigned int k, qaws_exact_poly* E);

/* det M and the row-0 cofactors C00, C01 (n >= 2). */
qaws_status qaws_exact_det_cofactors(qaws_exact_poly const* M, unsigned int n, qaws_exact_poly* det, qaws_exact_poly* C00, qaws_exact_poly* C01);

/* out = sum_i H[i] p^i q^(k - i): the substitution s = p / q, times q^k. */
qaws_status qaws_exact_substitute(qaws_exact_poly const* H, unsigned int k, qaws_exact_poly const* p, qaws_exact_poly const* q, qaws_exact_poly* out);

/* Enclosure of num / den (den != 0) by doubles, rounded outward; exact when dyadic and representable. */
void qaws_exact_ratio_enclose(qaws_exact_int const* num, qaws_exact_int const* den, double* lo, double* hi, int* exact);

/*
 * Roots r in [0, 1] of R with s = C01 / C00 proven in [0, 1] and every
 * constraint proven, refined to the narrowest interval the integer budget
 * allows. *common is set (with QAWS_STATUS_CERTIFICATION_FAILED) when R
 * vanishes identically.
 */
qaws_status qaws_exact_solve_system(qaws_exact_poly const* R, qaws_exact_poly const* C00, qaws_exact_poly const* C01,
	qaws_exact_constraint const* cons, unsigned int ncons, qaws_exact_local_hit* out, unsigned int capacity, unsigned int* count, int* common);

#endif /* QAWS_EXACT_SOLVE_H */
