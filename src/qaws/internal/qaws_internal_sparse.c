#include "qaws_internal_sparse.h"
#include "../qaws_platform.h"
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  CSR allocation                                                     */
/* ------------------------------------------------------------------ */

qaws_status qaws_internal_csr_alloc(
	qaws_internal_csr* mat,
	unsigned int rows,
	unsigned int cols,
	unsigned int nnz)
{
	if (!mat) return QAWS_STATUS_INVALID_ARGUMENT;

	memset(mat, 0, sizeof(*mat));
	mat->rows = rows;
	mat->cols = cols;
	mat->nnz = nnz;

	mat->row_ptr = (unsigned int*)calloc(rows + 1, sizeof(unsigned int));
	mat->col_idx = (unsigned int*)malloc(nnz * sizeof(unsigned int));
	mat->values = (qaws_scalar*)calloc(nnz, sizeof(qaws_scalar));

	if (!mat->row_ptr || !mat->col_idx || !mat->values)
	{
		qaws_internal_csr_free(mat);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

void qaws_internal_csr_free(qaws_internal_csr* mat)
{
	if (!mat) return;
	free(mat->row_ptr);
	free(mat->col_idx);
	free(mat->values);
	memset(mat, 0, sizeof(*mat));
}

/* ------------------------------------------------------------------ */
/*  CSR matrix-vector multiply: y = A * x                             */
/* ------------------------------------------------------------------ */

static void csr_matvec(
	qaws_internal_csr const* A,
	qaws_scalar const* x,
	qaws_scalar* y)
{
	unsigned int i, j;
	for (i = 0; i < A->rows; i++)
	{
		qaws_scalar sum = QAWS_ZERO;
		for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
			sum += A->values[j] * x[A->col_idx[j]];
		y[i] = sum;
	}
}

/* ------------------------------------------------------------------ */
/*  Conjugate gradient with Jacobi preconditioner                      */
/* ------------------------------------------------------------------ */

qaws_status qaws_internal_cg_solve(
	qaws_internal_csr const* A,
	qaws_scalar const* b,
	qaws_scalar* x,
	unsigned int max_iterations,
	qaws_scalar tolerance)
{
	unsigned int n = A->rows;
	qaws_scalar* r = NULL;
	qaws_scalar* p = NULL;
	qaws_scalar* Ap = NULL;
	qaws_scalar* M_inv = NULL; /* Jacobi preconditioner (diagonal) */
	qaws_scalar* z = NULL;
	qaws_scalar rz, rz_new, alpha, beta, pAp;
	unsigned int iter, i, j;
	qaws_scalar tol_sq;

	if (!A || !b || !x)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (max_iterations == 0)
		max_iterations = 2 * n;

#if QAWS_SCALAR_IS_FLOAT
	if (tolerance <= QAWS_ZERO)
		tolerance = QAWS_LITERAL(1e-4);
#else
	if (tolerance <= QAWS_ZERO)
		tolerance = QAWS_LITERAL(1e-6);
#endif

	tol_sq = tolerance * tolerance;

	r = (qaws_scalar*)malloc(n * sizeof(qaws_scalar));
	p = (qaws_scalar*)malloc(n * sizeof(qaws_scalar));
	Ap = (qaws_scalar*)malloc(n * sizeof(qaws_scalar));
	M_inv = (qaws_scalar*)malloc(n * sizeof(qaws_scalar));
	z = (qaws_scalar*)malloc(n * sizeof(qaws_scalar));

	if (!r || !p || !Ap || !M_inv || !z)
	{
		free(r); free(p); free(Ap); free(M_inv); free(z);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Extract diagonal for Jacobi preconditioner */
	for (i = 0; i < n; i++)
	{
		qaws_scalar diag = QAWS_ZERO;
		for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
		{
			if (A->col_idx[j] == i)
			{
				diag = A->values[j];
				break;
			}
		}
		M_inv[i] = (QAWS_FABS(diag) > QAWS_LITERAL(1e-15)) ?
			(QAWS_ONE / diag) : QAWS_ONE;
	}

	/* Initial: x = 0, r = b - A*x = b */
	memset(x, 0, n * sizeof(qaws_scalar));
	memcpy(r, b, n * sizeof(qaws_scalar));

	/* z = M^{-1} r */
	rz = QAWS_ZERO;
	for (i = 0; i < n; i++)
	{
		z[i] = M_inv[i] * r[i];
		rz += r[i] * z[i];
	}
	memcpy(p, z, n * sizeof(qaws_scalar));

	for (iter = 0; iter < max_iterations; iter++)
	{
		/* Check convergence: ||r||^2 */
		{
			qaws_scalar r_norm_sq = QAWS_ZERO;
			for (i = 0; i < n; i++)
				r_norm_sq += r[i] * r[i];
			if (r_norm_sq < tol_sq)
				break;
		}

		/* Ap = A * p */
		csr_matvec(A, p, Ap);

		/* alpha = rz / (p . Ap) */
		pAp = QAWS_ZERO;
		for (i = 0; i < n; i++)
			pAp += p[i] * Ap[i];

		if (QAWS_FABS(pAp) < QAWS_LITERAL(1e-30))
			break;

		alpha = rz / pAp;

		/* x = x + alpha * p, r = r - alpha * Ap */
		for (i = 0; i < n; i++)
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		/* z = M^{-1} r, rz_new = r . z */
		rz_new = QAWS_ZERO;
		for (i = 0; i < n; i++)
		{
			z[i] = M_inv[i] * r[i];
			rz_new += r[i] * z[i];
		}

		if (QAWS_FABS(rz) < QAWS_LITERAL(1e-30))
			break;

		beta = rz_new / rz;
		rz = rz_new;

		/* p = z + beta * p */
		for (i = 0; i < n; i++)
			p[i] = z[i] + beta * p[i];
	}

	free(r); free(p); free(Ap); free(M_inv); free(z);
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Build cotangent Laplacian + mass diagonal                          */
/* ------------------------------------------------------------------ */

qaws_status qaws_internal_build_laplacian(
	qaws_internal_mesh const* mesh,
	qaws_internal_csr* out_L,
	qaws_scalar* out_mass_diag)
{
	unsigned int n = mesh->vertex_count;
	unsigned int ei;
	unsigned int* degree = NULL;
	unsigned int* row_fill = NULL;
	unsigned int nnz;
	unsigned int vi, ti;
	qaws_status status;

	if (!mesh || !out_L || !out_mass_diag)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Count non-zeros: each edge contributes 2 off-diagonal entries, plus n diagonal */
	nnz = n + 2 * mesh->edge_count;

	status = qaws_internal_csr_alloc(out_L, n, n, nnz);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Count entries per row: 1 (diagonal) + degree */
	degree = (unsigned int*)calloc(n, sizeof(unsigned int));
	row_fill = (unsigned int*)calloc(n, sizeof(unsigned int));
	if (!degree || !row_fill)
	{
		free(degree); free(row_fill);
		qaws_internal_csr_free(out_L);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (ei = 0; ei < mesh->edge_count; ei++)
	{
		degree[mesh->edges[ei].v[0]]++;
		degree[mesh->edges[ei].v[1]]++;
	}

	/* Build row_ptr */
	out_L->row_ptr[0] = 0;
	for (vi = 0; vi < n; vi++)
		out_L->row_ptr[vi + 1] = out_L->row_ptr[vi] + 1 + degree[vi];

	/* Fill diagonal entries first (at position 0 of each row) */
	for (vi = 0; vi < n; vi++)
	{
		unsigned int pos = out_L->row_ptr[vi];
		out_L->col_idx[pos] = vi;
		out_L->values[pos] = QAWS_ZERO;
		row_fill[vi] = 1;
	}

	/* Fill off-diagonal entries from cotangent weights */
	for (ei = 0; ei < mesh->edge_count; ei++)
	{
		unsigned int va = mesh->edges[ei].v[0];
		unsigned int vb = mesh->edges[ei].v[1];
		qaws_scalar w = qaws_internal_mesh_cotan_weight(mesh, ei);
		unsigned int pos_a, pos_b;

		/* Entry (va, vb) = w */
		pos_a = out_L->row_ptr[va] + row_fill[va]++;
		out_L->col_idx[pos_a] = vb;
		out_L->values[pos_a] = w;

		/* Entry (vb, va) = w */
		pos_b = out_L->row_ptr[vb] + row_fill[vb]++;
		out_L->col_idx[pos_b] = va;
		out_L->values[pos_b] = w;

		/* Diagonal: L[va,va] -= w, L[vb,vb] -= w */
		out_L->values[out_L->row_ptr[va]] -= w;
		out_L->values[out_L->row_ptr[vb]] -= w;
	}

	free(degree);
	free(row_fill);

	/* Compute lumped mass diagonal (1/3 of area of adjacent triangles) */
	memset(out_mass_diag, 0, n * sizeof(qaws_scalar));

	for (ti = 0; ti < mesh->tri_count; ti++)
	{
		unsigned int va = mesh->tris[ti].v[0];
		unsigned int vb = mesh->tris[ti].v[1];
		unsigned int vc = mesh->tris[ti].v[2];
		qaws_scalar const* pa = mesh->positions + va * 3;
		qaws_scalar const* pb = mesh->positions + vb * 3;
		qaws_scalar const* pc = mesh->positions + vc * 3;

		/* Cross product for triangle area */
		qaws_scalar abx = pb[0] - pa[0], aby = pb[1] - pa[1], abz = pb[2] - pa[2];
		qaws_scalar acx = pc[0] - pa[0], acy = pc[1] - pa[1], acz = pc[2] - pa[2];
		qaws_scalar cx = aby * acz - abz * acy;
		qaws_scalar cy = abz * acx - abx * acz;
		qaws_scalar cz = abx * acy - aby * acx;
		qaws_scalar area = QAWS_LITERAL(0.5) * QAWS_SQRT(cx * cx + cy * cy + cz * cz);
		qaws_scalar third_area = area / QAWS_LITERAL(3.0);

		out_mass_diag[va] += third_area;
		out_mass_diag[vb] += third_area;
		out_mass_diag[vc] += third_area;
	}

	/* Ensure no zero mass (degenerate triangles) */
	for (vi = 0; vi < n; vi++)
	{
		if (out_mass_diag[vi] < QAWS_LITERAL(1e-15))
			out_mass_diag[vi] = QAWS_LITERAL(1e-15);
	}

	return QAWS_STATUS_OK;
}
