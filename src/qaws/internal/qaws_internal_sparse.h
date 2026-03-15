#ifndef QAWS_INTERNAL_SPARSE_H
#define QAWS_INTERNAL_SPARSE_H

#include "qaws_internal_types.h"
#include "qaws_internal_mesh.h"

/*
 * CSR sparse matrix and conjugate gradient solver.
 * Used internally for the heat method geodesic distance computation.
 */

typedef struct qaws_internal_csr {
	unsigned int rows, cols, nnz;
	unsigned int* row_ptr;  /* rows+1 */
	unsigned int* col_idx;  /* nnz */
	qaws_scalar* values;    /* nnz */
} qaws_internal_csr;

/* Allocate CSR matrix with given dimensions and nnz. Returns NULL fields on failure. */
qaws_status qaws_internal_csr_alloc(
	qaws_internal_csr* mat,
	unsigned int rows,
	unsigned int cols,
	unsigned int nnz);

void qaws_internal_csr_free(qaws_internal_csr* mat);

/* Conjugate gradient solve: A * x = b.
   x must be pre-allocated [rows]. b is [rows].
   Returns QAWS_STATUS_OK on convergence. */
qaws_status qaws_internal_cg_solve(
	qaws_internal_csr const* A,
	qaws_scalar const* b,
	qaws_scalar* x,
	unsigned int max_iterations,
	qaws_scalar tolerance);

/* Build cotangent Laplacian L and lumped mass diagonal A from mesh.
   L is [n x n] symmetric negative semi-definite.
   mass_diag is [n] (Voronoi area per vertex). */
qaws_status qaws_internal_build_laplacian(
	qaws_internal_mesh const* mesh,
	qaws_internal_csr* out_L,
	qaws_scalar* out_mass_diag);

#endif /* QAWS_INTERNAL_SPARSE_H */
