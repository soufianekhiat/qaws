#include "qaws_internal_fit.h"
#include "../qaws_bspline.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../qaws_platform.h"

/* ========================================================================== */
/*  Internal: band-diagonal solver for B-spline fitting                       */
/* ========================================================================== */

/* Evaluate B-spline basis function N_{i,p}(t) using de Boor recursion.
   knots[] has (n + p + 1) entries where n = control_point_count. */
static qaws_scalar basis_function(
	unsigned int i, unsigned int p,
	qaws_scalar t,
	qaws_scalar const* knots,
	unsigned int knot_count)
{
	qaws_scalar* left = NULL;
	qaws_scalar* right = NULL;
	qaws_scalar result;
	unsigned int j;

	(void)knot_count;

	if (p == 0)
		return (t >= knots[i] && t < knots[i + 1]) ? QAWS_ONE : QAWS_ZERO;

	left = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(p + 1));
	right = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(p + 1));
	if (!left || !right)
	{
		free(left);
		free(right);
		return QAWS_ZERO;
	}

	/* Triangular table algorithm */
	{
		qaws_scalar* N = (qaws_scalar*)calloc((size_t)(p + 1), sizeof(qaws_scalar));
		unsigned int k;
		if (!N)
		{
			free(left);
			free(right);
			return QAWS_ZERO;
		}

		/* Degree 0 */
		for (j = 0; j <= p; j++)
		{
			unsigned int idx = i + j;
			N[j] = (t >= knots[idx] && t < knots[idx + 1]) ? QAWS_ONE : QAWS_ZERO;
		}

		/* Build up to degree p */
		for (k = 1; k <= p; k++)
		{
			for (j = 0; j <= p - k; j++)
			{
				unsigned int idx = i + j;
				qaws_scalar denom_l = knots[idx + k] - knots[idx];
				qaws_scalar denom_r = knots[idx + k + 1] - knots[idx + 1];
				qaws_scalar val = QAWS_ZERO;
				if (QAWS_FABS(denom_l) > QAWS_LITERAL(1e-30))
					val += (t - knots[idx]) / denom_l * N[j];
				if (QAWS_FABS(denom_r) > QAWS_LITERAL(1e-30))
					val += (knots[idx + k + 1] - t) / denom_r * N[j + 1];
				N[j] = val;
			}
		}

		result = N[0];
		free(N);
	}

	free(left);
	free(right);
	return result;
}

/* Build clamped uniform knot vector for degree p, n control points.
   knot_count = n + p + 1.  Knots are clamped: first p+1 = t_min, last p+1 = t_max,
   interior knots uniformly spaced. */
static void build_clamped_knots(
	unsigned int n, unsigned int p,
	qaws_scalar t_min, qaws_scalar t_max,
	qaws_scalar* knots)
{
	unsigned int knot_count = n + p + 1;
	unsigned int interior = n - p; /* number of interior knot spans */
	unsigned int i;

	for (i = 0; i <= p; i++)
		knots[i] = t_min;

	for (i = 1; i < interior; i++)
	{
		qaws_scalar frac = (qaws_scalar)i / (qaws_scalar)interior;
		knots[p + i] = t_min + frac * (t_max - t_min);
	}

	for (i = 0; i <= p; i++)
		knots[knot_count - 1 - i] = t_max;
}

/* Solve a symmetric positive definite banded system A * x = b.
   A is n×n with half-bandwidth bw (A[i][j] = 0 for |i-j| > bw).
   A is stored as a[i * (bw+1) + j - i] for j >= i, j - i <= bw.
   Uses Cholesky decomposition for the banded case. */
static int solve_banded_spd(
	qaws_scalar* a_band,   /* n * (bw+1) upper triangular band storage */
	qaws_scalar* b,         /* n rhs values (overwritten with solution) */
	unsigned int n,
	unsigned int bw)
{
	unsigned int stride = bw + 1;
	unsigned int i, j, k;
	qaws_scalar sum;

	/* Cholesky factorization: L * L^T = A, stored in band form */
	for (i = 0; i < n; i++)
	{
		unsigned int jmax = (i + bw < n) ? i + bw : n - 1;

		/* Diagonal element */
		sum = a_band[i * stride];
		for (k = (i > bw ? i - bw : 0); k < i; k++)
		{
			qaws_scalar lik = a_band[k * stride + (i - k)];
			sum -= lik * lik;
		}
		if (sum <= QAWS_ZERO)
			return 0; /* not positive definite */
		a_band[i * stride] = QAWS_SQRT(sum);

		/* Off-diagonal elements in column i */
		for (j = i + 1; j <= jmax; j++)
		{
			sum = a_band[i * stride + (j - i)];
			for (k = (j > bw ? j - bw : 0); k < i; k++)
			{
				if (i >= k && (i - k) <= bw && j >= k && (j - k) <= bw)
				{
					qaws_scalar lik = a_band[k * stride + (i - k)];
					qaws_scalar ljk = a_band[k * stride + (j - k)];
					sum -= lik * ljk;
				}
			}
			a_band[i * stride + (j - i)] = sum / a_band[i * stride];
		}
	}

	/* Forward substitution: L * y = b */
	for (i = 0; i < n; i++)
	{
		sum = b[i];
		for (k = (i > bw ? i - bw : 0); k < i; k++)
			sum -= a_band[k * stride + (i - k)] * b[k];
		b[i] = sum / a_band[i * stride];
	}

	/* Back substitution: L^T * x = y */
	for (i = n; i > 0; i--)
	{
		unsigned int ii = i - 1;
		unsigned int jmax = (ii + bw < n) ? ii + bw : n - 1;
		sum = b[ii];
		for (j = ii + 1; j <= jmax; j++)
			sum -= a_band[ii * stride + (j - ii)] * b[j];
		b[ii] = sum / a_band[ii * stride];
	}

	return 1;
}

/* ========================================================================== */
/*  Public: fit B-spline through sample points                                */
/* ========================================================================== */

qaws_status qaws_internal_fit_bspline_range(
	qaws_dimension dimension,
	unsigned int degree,
	qaws_scalar const* sample_params,
	qaws_scalar const* sample_coords,
	unsigned int sample_count,
	unsigned int control_point_count,
	qaws_scalar t_min,
	qaws_scalar t_max,
	qaws_curve** out_curve)
{
	unsigned int n = control_point_count;
	unsigned int p = degree;
	unsigned int dim = (unsigned int)dimension;
	unsigned int knot_count = n + p + 1;
	unsigned int bw = p; /* bandwidth of the normal equations */
	unsigned int stride = bw + 1;
	qaws_scalar* knots = NULL;
	qaws_scalar* ata_band = NULL; /* n * stride */
	qaws_scalar* atb = NULL;      /* n * dim */
	qaws_scalar* basis_row = NULL; /* n basis values for one sample */
	qaws_scalar* cps = NULL;
	qaws_bspline_desc desc;
	unsigned int si, ci, d, di;
	qaws_status status;

	if (!sample_params || !sample_coords || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (sample_count < 2 || control_point_count < degree + 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (dimension != QAWS_DIMENSION_2D && dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_curve = NULL;

	/* Allocate working buffers */
	knots = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)knot_count);
	ata_band = (qaws_scalar*)calloc((size_t)(n * stride), sizeof(qaws_scalar));
	atb = (qaws_scalar*)calloc((size_t)(n * dim), sizeof(qaws_scalar));
	basis_row = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n);
	cps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));

	if (!knots || !ata_band || !atb || !basis_row || !cps)
	{
		free(knots); free(ata_band); free(atb); free(basis_row); free(cps);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Build clamped knot vector */
	build_clamped_knots(n, p, t_min, t_max, knots);

	/* Assemble normal equations A^T A and A^T b.
	   For each sample point, evaluate all basis functions and accumulate. */
	for (si = 0; si < sample_count; si++)
	{
		qaws_scalar t = sample_params[si];

		/* Clamp the last sample to be slightly inside the domain
		   so the last basis function gets a non-zero value */
		if (t >= t_max)
			t = t_max - QAWS_LITERAL(1e-10);
		if (t < t_min)
			t = t_min;

		/* Evaluate all basis functions at t */
		for (ci = 0; ci < n; ci++)
			basis_row[ci] = basis_function(ci, p, t, knots, knot_count);

		/* Accumulate A^T A (banded upper triangle) */
		for (ci = 0; ci < n; ci++)
		{
			unsigned int cj;
			unsigned int jmax = (ci + bw < n) ? ci + bw : n - 1;
			for (cj = ci; cj <= jmax; cj++)
				ata_band[ci * stride + (cj - ci)] += basis_row[ci] * basis_row[cj];

			/* Accumulate A^T b for each dimension */
			for (di = 0; di < dim; di++)
				atb[ci * dim + di] += basis_row[ci] * sample_coords[si * dim + di];
		}
	}

	/* Solve per-dimension: A^T A * cp_d = atb_d */
	for (d = 0; d < dim; d++)
	{
		/* Copy the band matrix (it gets overwritten by Cholesky) */
		qaws_scalar* ata_copy = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(n * stride));
		qaws_scalar* rhs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n);
		if (!ata_copy || !rhs)
		{
			free(ata_copy); free(rhs);
			free(knots); free(ata_band); free(atb); free(basis_row); free(cps);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		memcpy(ata_copy, ata_band, sizeof(qaws_scalar) * (size_t)(n * stride));

		for (ci = 0; ci < n; ci++)
			rhs[ci] = atb[ci * dim + d];

		if (!solve_banded_spd(ata_copy, rhs, n, bw))
		{
			/* Fallback: add regularization and retry */
			memcpy(ata_copy, ata_band, sizeof(qaws_scalar) * (size_t)(n * stride));
			for (ci = 0; ci < n; ci++)
				ata_copy[ci * stride] += QAWS_LITERAL(1e-6);
			for (ci = 0; ci < n; ci++)
				rhs[ci] = atb[ci * dim + d];
			if (!solve_banded_spd(ata_copy, rhs, n, bw))
			{
				free(ata_copy); free(rhs);
				free(knots); free(ata_band); free(atb); free(basis_row); free(cps);
				return QAWS_STATUS_DEGENERATE_CURVE;
			}
		}

		for (ci = 0; ci < n; ci++)
			cps[ci * dim + d] = rhs[ci];

		free(ata_copy);
		free(rhs);
	}

	/* Create B-spline curve from the computed control points */
	memset(&desc, 0, sizeof(desc));
	desc.dimension = dimension;
	desc.degree = degree;
	desc.control_points = cps;
	desc.control_point_count = n;
	desc.knots = knots;
	desc.knot_count = knot_count;
	desc.is_uniform = 0;
	desc.is_closed = 0;

	status = qaws_curve_create_bspline(&desc, out_curve);

	free(knots);
	free(ata_band);
	free(atb);
	free(basis_row);
	free(cps);

	return status;
}

qaws_status qaws_internal_fit_bspline(
	qaws_dimension dimension,
	unsigned int degree,
	qaws_scalar const* sample_params,
	qaws_scalar const* sample_coords,
	unsigned int sample_count,
	unsigned int control_point_count,
	qaws_curve** out_curve)
{
	qaws_scalar t_min, t_max;
	if (!sample_params || sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	t_min = sample_params[0];
	t_max = sample_params[sample_count - 1];

	return qaws_internal_fit_bspline_range(
		dimension, degree,
		sample_params, sample_coords,
		sample_count, control_point_count,
		t_min, t_max, out_curve);
}
