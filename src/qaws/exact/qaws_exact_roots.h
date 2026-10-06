#ifndef QAWS_EXACT_ROOTS_H
#define QAWS_EXACT_ROOTS_H

/* Internal: certified real root isolation of integer Bernstein polynomials. */

#include "qaws_exact_int.h"
#include <stdint.h>

#define QAWS_EXACT_ROOTS_MAX_DEGREE 64
#define QAWS_EXACT_ROOTS_MAX_DEPTH 60

/*
 * A root of p(s) = sum b_i B_i^n(s) on [0, 1], located on the dyadic grid:
 *   exact == 1: the root is s = index / 2^depth exactly;
 *   exact == 0: the open interval (index, index + 1) / 2^depth holds exactly
 *               one root, a simple one (p changes sign across it).
 */
typedef struct qaws_exact_root
{
	uint64_t index;
	int depth;
	int exact;
} qaws_exact_root;

/*
 * Every root of p in [0, 1], sorted. Isolation subdivides at midpoints and
 * counts the sign variations of the Bernstein coefficients (Descartes' rule
 * in Bernstein form: zero variations, no root in the open interval; one,
 * exactly one simple root). A cluster that does not separate within
 * QAWS_EXACT_ROOTS_MAX_DEPTH halvings (a multiple root: a tangency) gives
 * QAWS_STATUS_CERTIFICATION_FAILED; the zero polynomial too (every s is a
 * root). QAWS_STATUS_BUFFER_TOO_SMALL when more than `capacity` roots.
 */
qaws_status qaws_exact_bernstein_isolate(qaws_exact_int const* b, unsigned int n, qaws_exact_root* out, unsigned int capacity,
	unsigned int* out_count);

/* Coefficients of p on [index, index + 1] / 2^depth (scaled by a positive factor). */
qaws_status qaws_exact_bernstein_restrict(qaws_exact_int const* b, unsigned int n, uint64_t index, int depth, qaws_exact_int* out);

/* Two polynomials of degree n restricted together: one common positive factor (their ratio is kept). */
qaws_status qaws_exact_bernstein_restrict_pair(qaws_exact_int const* b1, qaws_exact_int const* b2, unsigned int n, uint64_t index, int depth,
	qaws_exact_int* out1, qaws_exact_int* out2);

/* Shrinks an isolating interval to `depth` (by sign bisection); it may land on the exact root. */
qaws_status qaws_exact_bernstein_refine(qaws_exact_int const* b, unsigned int n, qaws_exact_root* root, int depth);

#endif /* QAWS_EXACT_ROOTS_H */
