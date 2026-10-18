#ifndef QAWS_INTERNAL_LDLT_H
#define QAWS_INTERNAL_LDLT_H

/*
 * Sparse symmetric LDL^T: the rows ordered by geometric nested dissection
 * (George 1973; each set split at the median of its box's longest axis,
 * the left rows touching the right ones last), the elimination tree and
 * column counts as Davis's LDL (2005), the columns grouped into
 * fundamental supernodes (a chain of columns sharing their rows below) and
 * factored left-looking (Ng, Peyton 1993): each supernode a dense panel,
 * the updates of the supernodes below it dense products. Double precision
 * whatever qaws_scalar is.
 */

#include <stddef.h>
#include "../qaws_types.h"
#include "../qaws_status.h"

typedef struct qaws_internal_ldlt
{
	unsigned int n, ns;             /* rows, supernodes */
	unsigned int* perm;             /* the factor's row k is the matrix's row perm[k] */
	unsigned int *sf, *rp, *ri;     /* per supernode its first column, its rows' start (ns + 1 each); the rows */
	size_t* lp;                     /* per supernode its panel's start (ns + 1): column-major, its rows by its width */
	double *lx, *d;                 /* the panels (L below the diagonal, its unit diagonal implicit), D */
	double* w;                      /* the solve's workspace (2 n) */
	double nnz;                     /* L's entries below the diagonal */
} qaws_internal_ldlt;

/* A: a full symmetric pattern in CSR (only the lower triangle's values are read); xyz: 3 per row for the ordering (two rows
   of one 2 x 2 block given the same point stay together), NULL for the natural order; shift added to the diagonal; on a
   pivot that is not positive the factor is redone with 10 times the shift (a 1e-12 relative one from 0), INVALID_ARGUMENT
   after 3 retries */
qaws_status qaws_internal_ldlt_factor(qaws_internal_ldlt* f, unsigned int n, unsigned int const* row_ptr, unsigned int const* col,
	double const* val, double const* xyz, double shift);
/* A x = b; b and x may alias */
void qaws_internal_ldlt_solve(qaws_internal_ldlt const* f, double const* b, double* x);
void qaws_internal_ldlt_free(qaws_internal_ldlt* f);

#endif /* QAWS_INTERNAL_LDLT_H */
