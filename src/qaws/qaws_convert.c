#include "qaws_alloc.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_convert.h"
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "internal/qaws_internal_types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* -------------------------------------------------------------------------- */
/*  Hermite to Bezier (per span)                                              */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_hermite_to_bezier(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_curve **out_bezier)
{
	qaws_hermite_impl const *impl;
	unsigned int dim_count;
	qaws_scalar const *P0;
	qaws_scalar const *P1;
	qaws_scalar const *M0;
	qaws_scalar const *M1;
	qaws_scalar *cp;
	size_t cp_size;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int d;

	if (!curve || !out_bezier)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_bezier = NULL;

	if (curve->kind != QAWS_CURVE_KIND_HERMITE)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	if (span_index >= curve->span_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	impl = (qaws_hermite_impl const *)curve->impl;
	dim_count = (unsigned int)curve->dimension;

	P0 = impl->points + span_index * dim_count;
	P1 = impl->points + (span_index + 1) * dim_count;
	M0 = impl->tangents + span_index * dim_count;
	M1 = impl->tangents + (span_index + 1) * dim_count;

	cp_size = qaws_internal_mul3(4, dim_count, sizeof(qaws_scalar));
	cp = (qaws_scalar *)qaws_internal_malloc(cp_size);
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (d = 0; d < dim_count; d++) {
		cp[0 * dim_count + d] = P0[d];
		cp[1 * dim_count + d] = P0[d] + M0[d] / (qaws_scalar)3.0;
		cp[2 * dim_count + d] = P1[d] - M1[d] / (qaws_scalar)3.0;
		cp[3 * dim_count + d] = P1[d];
	}

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	status = qaws_curve_create_bezier_ex(&desc, curve->allocator, out_bezier);
	qaws_internal_free(cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Catmull-Rom to Bezier (per span)                                          */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_catmull_rom_to_bezier(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_curve **out_bezier)
{
	qaws_catmull_rom_impl const *impl;
	unsigned int dim_count;
	qaws_scalar const *coeffs;
	qaws_scalar *cp;
	size_t cp_size;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int d;

	if (!curve || !out_bezier)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_bezier = NULL;

	if (curve->kind != QAWS_CURVE_KIND_CATMULL_ROM)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	if (span_index >= curve->span_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	impl = (qaws_catmull_rom_impl const *)curve->impl;
	dim_count = (unsigned int)curve->dimension;

	cp_size = qaws_internal_mul3(4, dim_count, sizeof(qaws_scalar));
	cp = (qaws_scalar *)qaws_internal_malloc(cp_size);
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/*
	 * Coefficients are stored per span as dim_count groups of 4 scalars.
	 * For dimension d: coeffs[d*4+0]=a, [d*4+1]=b, [d*4+2]=c, [d*4+3]=d_coeff
	 *
	 * Polynomial: P(t) = a*t^3 + b*t^2 + c*t + d_coeff  over [0,1]
	 *
	 * Hermite-to-Bezier mapping:
	 *   B0 = d_coeff
	 *   B1 = d_coeff + c/3
	 *   B2 = d_coeff + 2c/3 + b/3
	 *   B3 = a + b + c + d_coeff
	 */
	coeffs = impl->segment_coeffs + span_index * dim_count * 4;

	for (d = 0; d < dim_count; d++) {
		qaws_scalar a = coeffs[d * 4 + 0];
		qaws_scalar b = coeffs[d * 4 + 1];
		qaws_scalar c = coeffs[d * 4 + 2];
		qaws_scalar d_coeff = coeffs[d * 4 + 3];

		cp[0 * dim_count + d] = d_coeff;
		cp[1 * dim_count + d] = d_coeff + c / (qaws_scalar)3.0;
		cp[2 * dim_count + d] = d_coeff + (qaws_scalar)2.0 * c / (qaws_scalar)3.0
		                      + b / (qaws_scalar)3.0;
		cp[3 * dim_count + d] = a + b + c + d_coeff;
	}

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	status = qaws_curve_create_bezier_ex(&desc, curve->allocator, out_bezier);
	qaws_internal_free(cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Bezier to B-Spline                                                        */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_bezier_to_bspline(
	qaws_curve const *curve,
	qaws_curve **out_bspline)
{
	qaws_bezier_impl const *impl;
	unsigned int degree;
	unsigned int knot_count;
	qaws_scalar *knots;
	qaws_bspline_desc desc;
	qaws_status status;
	unsigned int i;

	if (!curve || !out_bspline)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_bspline = NULL;

	if (curve->kind != QAWS_CURVE_KIND_BEZIER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bezier_impl const *)curve->impl;
	degree = curve->degree;

	/* Clamped knot vector: [0,...,0, 1,...,1] with (degree+1) repeats each */
	knot_count = 2 * (degree + 1);
	knots = (qaws_scalar *)qaws_internal_malloc(qaws_internal_mul2((size_t)knot_count, sizeof(qaws_scalar)));
	if (!knots)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i <= degree; i++)
		knots[i] = (qaws_scalar)0.0;
	for (i = degree + 1; i < knot_count; i++)
		knots[i] = (qaws_scalar)1.0;

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = degree;
	desc.control_points = impl->control_points;
	desc.control_point_count = impl->control_point_count;
	desc.knots = knots;
	desc.knot_count = knot_count;
	desc.is_uniform = 0;
	desc.is_closed = 0;

	status = qaws_curve_create_bspline_ex(&desc, curve->allocator, out_bspline);
	qaws_internal_free(knots);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  B-Spline to NURBS                                                         */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_bspline_to_nurbs(
	qaws_curve const *curve,
	qaws_curve **out_nurbs)
{
	qaws_bspline_impl const *impl;
	unsigned int cp_count;
	qaws_scalar *weights;
	qaws_nurbs_desc desc;
	qaws_status status;
	unsigned int i;

	if (!curve || !out_nurbs)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_nurbs = NULL;

	if (curve->kind != QAWS_CURVE_KIND_BSPLINE)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bspline_impl const *)curve->impl;
	cp_count = impl->control_point_count;

	/* All weights = 1.0 for a non-rational B-Spline */
	weights = (qaws_scalar *)qaws_internal_malloc(qaws_internal_mul2((size_t)cp_count, sizeof(qaws_scalar)));
	if (!weights)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i < cp_count; i++)
		weights[i] = (qaws_scalar)1.0;

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = curve->degree;
	desc.control_points = impl->control_points;
	desc.control_point_count = cp_count;
	desc.knots = impl->knots;
	desc.knot_count = impl->knot_count;
	desc.weights = weights;
	desc.weight_count = cp_count;
	desc.is_closed = 0;

	status = qaws_curve_create_nurbs_ex(&desc, curve->allocator, out_nurbs);
	qaws_internal_free(weights);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Degree elevation (Bezier only)                                            */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_elevate_degree(
	qaws_curve const *curve,
	qaws_curve **out_elevated)
{
	qaws_bezier_impl const *impl;
	unsigned int n;
	unsigned int dim_count;
	unsigned int new_count;
	qaws_scalar *new_cp;
	size_t new_cp_size;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int i, d;

	if (!curve || !out_elevated)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_elevated = NULL;

	if (curve->kind != QAWS_CURVE_KIND_BEZIER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bezier_impl const *)curve->impl;
	n = curve->degree;
	dim_count = (unsigned int)curve->dimension;
	new_count = n + 2;

	new_cp_size = qaws_internal_mul3(new_count, dim_count, sizeof(qaws_scalar));
	new_cp = (qaws_scalar *)qaws_internal_malloc(new_cp_size);
	if (!new_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Q[0] = P[0] */
	for (d = 0; d < dim_count; d++)
		new_cp[0 * dim_count + d] = impl->control_points[0 * dim_count + d];

	/* Q[i] = (i/(n+1)) * P[i-1] + (1 - i/(n+1)) * P[i]  for i = 1..n */
	for (i = 1; i <= n; i++) {
		qaws_scalar alpha = (qaws_scalar)i / (qaws_scalar)(n + 1);
		qaws_scalar one_minus_alpha = (qaws_scalar)1.0 - alpha;
		for (d = 0; d < dim_count; d++) {
			new_cp[i * dim_count + d] =
				alpha * impl->control_points[(i - 1) * dim_count + d]
				+ one_minus_alpha * impl->control_points[i * dim_count + d];
		}
	}

	/* Q[n+1] = P[n] */
	for (d = 0; d < dim_count; d++)
		new_cp[(n + 1) * dim_count + d] = impl->control_points[n * dim_count + d];

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = n + 1;
	desc.control_points = new_cp;
	desc.control_point_count = new_count;

	status = qaws_curve_create_bezier_ex(&desc, curve->allocator, out_elevated);
	qaws_internal_free(new_cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Helper: binomial coefficient C(n, k) using multiplicative formula         */
/* -------------------------------------------------------------------------- */

static double deg_reduce_binom(unsigned int n, unsigned int k)
{
	double result = 1.0;
	unsigned int i;

	if (k > n)
		return 0.0;
	if (k > n - k)
		k = n - k;

	for (i = 0; i < k; i++) {
		result *= (double)(n - i);
		result /= (double)(i + 1);
	}
	return result;
}

/* -------------------------------------------------------------------------- */
/*  Helper: Gaussian elimination with partial pivoting                        */
/*  Solves A * x = b in-place. A is sys_n x sys_n, b is sys_n x 1.           */
/*  On return, b contains the solution.                                       */
/*  Returns 0 on success, -1 on singular matrix.                              */
/* -------------------------------------------------------------------------- */

static int deg_reduce_gauss_solve(double *A, double *b, unsigned int sys_n)
{
	unsigned int col, row, i;

	for (col = 0; col < sys_n; col++) {
		/* Partial pivoting: find row with largest absolute value in column */
		unsigned int pivot_row = col;
		double pivot_val = fabs(A[col * sys_n + col]);

		for (row = col + 1; row < sys_n; row++) {
			double v = fabs(A[row * sys_n + col]);
			if (v > pivot_val) {
				pivot_val = v;
				pivot_row = row;
			}
		}

		if (pivot_val < 1e-30)
			return -1; /* singular */

		/* Swap rows if needed */
		if (pivot_row != col) {
			for (i = 0; i < sys_n; i++) {
				double tmp = A[col * sys_n + i];
				A[col * sys_n + i] = A[pivot_row * sys_n + i];
				A[pivot_row * sys_n + i] = tmp;
			}
			{
				double tmp = b[col];
				b[col] = b[pivot_row];
				b[pivot_row] = tmp;
			}
		}

		/* Eliminate below */
		for (row = col + 1; row < sys_n; row++) {
			double factor = A[row * sys_n + col] / A[col * sys_n + col];
			for (i = col; i < sys_n; i++)
				A[row * sys_n + i] -= factor * A[col * sys_n + i];
			b[row] -= factor * b[col];
		}
	}

	/* Back substitution */
	for (i = sys_n; i > 0; i--) {
		unsigned int idx = i - 1;
		unsigned int j;
		double sum = b[idx];

		for (j = idx + 1; j < sys_n; j++)
			sum -= A[idx * sys_n + j] * b[j];
		b[idx] = sum / A[idx * sys_n + idx];
	}

	return 0;
}

/* -------------------------------------------------------------------------- */
/*  Degree reduction (Bezier only, n -> n-1)                                  */
/*  Least-squares with endpoint interpolation.                                */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_reduce_degree(
	qaws_curve const *curve,
	qaws_curve **out_reduced)
{
	qaws_bezier_impl const *impl;
	unsigned int n;       /* original degree */
	unsigned int m;       /* reduced degree = n - 1 */
	unsigned int dim_count;
	unsigned int new_count; /* m + 1 */
	qaws_scalar *new_cp;
	size_t new_cp_size;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int i, j, k, d;
	unsigned int sys_n;   /* interior unknowns count = m - 1 */
	double *gram;         /* Gram matrix for interior points */
	double *cross_mat;    /* Cross-term matrix (m+1) x (n+1) */
	double *rhs;          /* RHS vector per dimension solve */
	double *A_copy;       /* working copy of Gram for each solve */

	if (!curve || !out_reduced)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_reduced = NULL;

	if (curve->kind != QAWS_CURVE_KIND_BEZIER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bezier_impl const *)curve->impl;
	n = curve->degree;
	dim_count = (unsigned int)curve->dimension;

	if (n < 2)
		return QAWS_STATUS_INVALID_DEGREE;

	m = n - 1;        /* reduced degree */
	new_count = m + 1; /* number of reduced control points */
	sys_n = m - 1;     /* number of interior unknowns (indices 1..m-1) */

	new_cp_size = qaws_internal_mul3(new_count, dim_count, sizeof(qaws_scalar));
	new_cp = (qaws_scalar *)qaws_internal_malloc(new_cp_size);
	if (!new_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Endpoint interpolation: Q[0] = P[0], Q[m] = P[n] */
	for (d = 0; d < dim_count; d++) {
		new_cp[0 * dim_count + d] = impl->control_points[0 * dim_count + d];
		new_cp[m * dim_count + d] = impl->control_points[n * dim_count + d];
	}

	/* Special case: degree 2 -> 1. No interior points to solve. */
	if (sys_n == 0) {
		memset(&desc, 0, sizeof(desc));
		desc.dimension = curve->dimension;
		desc.degree = m;
		desc.control_points = new_cp;
		desc.control_point_count = new_count;
		status = qaws_curve_create_bezier_ex(&desc, curve->allocator, out_reduced);
		qaws_internal_free(new_cp);
		return status;
	}

	/*
	 * Build the full cross-term matrix Cross[i][k] for i in [0..m], k in [0..n].
	 * Cross[i][k] = <B_i^m, B_k^n> = C(m,i)*C(n,k) / ((m+n+1)*C(m+n, i+k))
	 *
	 * Note: m + n = (n-1) + n = 2n - 1, so denominator = 2n * C(2n-1, i+k).
	 */
	cross_mat = (double *)qaws_internal_malloc(qaws_internal_mul3(((size_t)(new_count)), ((size_t)(n + 1)), sizeof(double)));
	if (!cross_mat) {
		qaws_internal_free(new_cp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i <= m; i++) {
		for (k = 0; k <= n; k++) {
			cross_mat[i * (n + 1) + k] =
				deg_reduce_binom(m, i) * deg_reduce_binom(n, k)
				/ ((double)(m + n + 1) * deg_reduce_binom(m + n, i + k));
		}
	}

	/*
	 * Build the interior Gram matrix: sys_n x sys_n.
	 * Gram[a][b] = <B_{a+1}^m, B_{b+1}^m>
	 *            = C(m, a+1)*C(m, b+1) / ((2m+1)*C(2m, a+1+b+1))
	 */
	gram = (double *)qaws_internal_malloc(qaws_internal_mul3(((size_t)(sys_n)), ((size_t)(sys_n)), sizeof(double)));
	if (!gram) {
		qaws_internal_free(cross_mat);
		qaws_internal_free(new_cp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < sys_n; i++) {
		unsigned int gi = i + 1; /* actual index in reduced basis */
		for (j = 0; j < sys_n; j++) {
			unsigned int gj = j + 1;
			gram[i * sys_n + j] =
				deg_reduce_binom(m, gi) * deg_reduce_binom(m, gj)
				/ ((double)(2 * m + 1) * deg_reduce_binom(2 * m, gi + gj));
		}
	}

	/* Allocate workspace for RHS and Gram copy */
	rhs = (double *)qaws_internal_malloc(qaws_internal_mul2((size_t)sys_n, sizeof(double)));
	A_copy = (double *)qaws_internal_malloc(qaws_internal_mul3(((size_t)(sys_n)), ((size_t)(sys_n)), sizeof(double)));
	if (!rhs || !A_copy) {
		qaws_internal_free(A_copy);
		qaws_internal_free(rhs);
		qaws_internal_free(gram);
		qaws_internal_free(cross_mat);
		qaws_internal_free(new_cp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Solve per dimension */
	for (d = 0; d < dim_count; d++) {
		/*
		 * RHS[a] = sum_{k=0}^{n} Cross[a+1][k] * P_k[d]
		 *        - Gram[a][0-1=skip] * Q_0[d]      (from fixing Q_0)
		 *        - Gram[a][sys_n-1+1-1=skip] * Q_m[d] (from fixing Q_m)
		 *
		 * More precisely:
		 * Full system: M * Q_interior = R_interior
		 * where R[a] = sum_k Cross[a+1][k]*P_k[d]
		 *            - M_full[a+1][0] * Q_0[d]
		 *            - M_full[a+1][m] * Q_m[d]
		 *
		 * M_full[i][j] = <B_i^m, B_j^m>
		 * M_full[a+1][0] = <B_{a+1}^m, B_0^m>
		 *   = C(m, a+1)*C(m, 0) / ((2m+1)*C(2m, a+1))
		 *   = C(m, a+1) / ((2m+1)*C(2m, a+1))
		 * M_full[a+1][m] = <B_{a+1}^m, B_m^m>
		 *   = C(m, a+1)*C(m, m) / ((2m+1)*C(2m, a+1+m))
		 *   = C(m, a+1) / ((2m+1)*C(2m, a+1+m))
		 */
		double Q0d = (double)impl->control_points[0 * dim_count + d];
		double Qmd = (double)impl->control_points[n * dim_count + d];

		for (i = 0; i < sys_n; i++) {
			unsigned int gi = i + 1;
			double sum = 0.0;
			double gram_i0, gram_im;

			/* Cross-term contribution */
			for (k = 0; k <= n; k++)
				sum += cross_mat[gi * (n + 1) + k]
				       * (double)impl->control_points[k * dim_count + d];

			/* Subtract endpoint contributions from Gram */
			gram_i0 = deg_reduce_binom(m, gi) * deg_reduce_binom(m, 0)
			          / ((double)(2 * m + 1) * deg_reduce_binom(2 * m, gi + 0));
			gram_im = deg_reduce_binom(m, gi) * deg_reduce_binom(m, m)
			          / ((double)(2 * m + 1) * deg_reduce_binom(2 * m, gi + m));

			rhs[i] = sum - gram_i0 * Q0d - gram_im * Qmd;
		}

		/* Copy Gram matrix (solver modifies it in-place) */
		memcpy(A_copy, gram, (size_t)(sys_n) * (size_t)(sys_n) * sizeof(double));

		/* Solve */
		if (deg_reduce_gauss_solve(A_copy, rhs, sys_n) != 0) {
			qaws_internal_free(A_copy);
			qaws_internal_free(rhs);
			qaws_internal_free(gram);
			qaws_internal_free(cross_mat);
			qaws_internal_free(new_cp);
			return QAWS_STATUS_NUMERICAL_FAILURE;
		}

		/* Store interior control points */
		for (i = 0; i < sys_n; i++)
			new_cp[(i + 1) * dim_count + d] = (qaws_scalar)rhs[i];
	}

	qaws_internal_free(A_copy);
	qaws_internal_free(rhs);
	qaws_internal_free(gram);
	qaws_internal_free(cross_mat);

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = m;
	desc.control_points = new_cp;
	desc.control_point_count = new_count;

	status = qaws_curve_create_bezier_ex(&desc, curve->allocator, out_reduced);
	qaws_internal_free(new_cp);
	return status;
}
