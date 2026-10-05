#include "qaws_internal_diff.h"
#include "qaws_internal_basis.h"
#include <string.h>

void qaws_internal_bspline_basis_derivs_any(
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int span,
	qaws_scalar t,
	unsigned int k,
	qaws_scalar* out_ders)
{
	unsigned int stride = degree + 1;
	unsigned int kk = (k < degree) ? k : degree;
	unsigned int r, j;

	qaws_internal_bspline_basis_derivs(knots, knot_count, degree, span, t, kk, out_ders);

	for (r = kk + 1; r <= k; r++)
		for (j = 0; j < stride; j++)
			out_ders[r * stride + j] = QAWS_ZERO;
}

void qaws_internal_bernstein_derivs(
	unsigned int n,
	qaws_scalar t,
	unsigned int k,
	qaws_scalar* out)
{
	/* levels[m][j] = B_{j,m}(t) for m = 0..n */
	qaws_scalar levels[QAWS_DIFF_MAX_SUPPORT][QAWS_DIFF_MAX_SUPPORT];
	qaws_scalar u = QAWS_ONE - t;
	unsigned int m, j, r, i;

	levels[0][0] = QAWS_ONE;
	for (m = 1; m <= n; m++)
	{
		levels[m][0] = u * levels[m - 1][0];
		for (j = 1; j < m; j++)
			levels[m][j] = u * levels[m - 1][j] + t * levels[m - 1][j - 1];
		levels[m][m] = t * levels[m - 1][m - 1];
	}

	for (r = 0; r <= k; r++)
	{
		for (i = 0; i <= n; i++)
			out[r * (n + 1) + i] = QAWS_ZERO;
		if (r > n)
			continue;

		{
			/* d^r B_{i,n} = n!/(n-r)! * sum_c (-1)^(r-c) C(r,c) B_{i-c, n-r} */
			qaws_scalar falling = QAWS_ONE;
			unsigned int mm = n - r;
			for (j = 0; j < r; j++)
				falling *= (qaws_scalar)(n - j);

			for (i = 0; i <= n; i++)
			{
				qaws_scalar sum = QAWS_ZERO;
				qaws_scalar binom = QAWS_ONE;
				unsigned int c;
				for (c = 0; c <= r; c++)
				{
					if (c > 0)
						binom = binom * (qaws_scalar)(r - c + 1) / (qaws_scalar)c;
					if (i >= c && i - c <= mm)
					{
						qaws_scalar b = levels[mm][i - c];
						sum += (((r - c) & 1u) ? -binom : binom) * b;
					}
				}
				out[r * (n + 1) + i] = falling * sum;
			}
		}
	}
}
