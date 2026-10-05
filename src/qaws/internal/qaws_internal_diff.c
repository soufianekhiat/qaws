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

/*
 * Basis derivatives in dual numbers along a knot direction (Piegl A2.3
 * with every knot carrying v + e kdot). kdot[i] seeds knot span-p+1+i,
 * the 2p knots the span's basis depends on.
 */
void qaws_internal_bspline_basis_derivs_knot(
	qaws_scalar const* knots,
	unsigned int degree,
	unsigned int span,
	qaws_scalar t,
	qaws_scalar const* kdot,
	unsigned int k,
	qaws_dual1* out_ders)
{
	qaws_dual1 ndu[QAWS_INTERNAL_KNOT_MAX_DEGREE + 1][QAWS_INTERNAL_KNOT_MAX_DEGREE + 1];
	qaws_dual1 a[2][QAWS_INTERNAL_KNOT_MAX_DEGREE + 1];
	qaws_dual1 left[QAWS_INTERNAL_KNOT_MAX_DEGREE + 1], right[QAWS_INTERNAL_KNOT_MAX_DEGREE + 1];
	qaws_dual1 const tt = qaws_dual1_const(t);
	unsigned int p = degree, stride = degree + 1, n = (k < degree) ? k : degree;
	unsigned int j, r, kk;
	int s1, s2;

#define KNOT_DUAL(idx) qaws_dual1_make(knots[idx], kdot[(idx) + p - 1 - span], QAWS_ZERO)

	ndu[0][0] = qaws_dual1_const(QAWS_ONE);
	for (j = 1; j <= p; j++)
	{
		qaws_dual1 saved = qaws_dual1_const(QAWS_ZERO);
		left[j] = qaws_dual1_sub(tt, KNOT_DUAL(span + 1 - j));
		right[j] = qaws_dual1_sub(KNOT_DUAL(span + j), tt);
		for (r = 0; r < j; r++)
		{
			qaws_dual1 temp;
			ndu[j][r] = qaws_dual1_add(right[r + 1], left[j - r]);
			temp = qaws_dual1_div(ndu[r][j - 1], ndu[j][r]);
			ndu[r][j] = qaws_dual1_add(saved, qaws_dual1_mul(right[r + 1], temp));
			saved = qaws_dual1_mul(left[j - r], temp);
		}
		ndu[j][j] = saved;
	}
#undef KNOT_DUAL

	for (j = 0; j <= p; j++)
		out_ders[j] = ndu[j][p];

	for (r = 0; r <= p; r++)
	{
		s1 = 0;
		s2 = 1;
		a[0][0] = qaws_dual1_const(QAWS_ONE);
		for (kk = 1; kk <= n; kk++)
		{
			qaws_dual1 d = qaws_dual1_const(QAWS_ZERO);
			int rk = (int)r - (int)kk, pk = (int)p - (int)kk;
			int j1, j2, jj;
			if (rk >= 0)
			{
				a[s2][0] = qaws_dual1_div(a[s1][0], ndu[pk + 1][rk]);
				d = qaws_dual1_mul(a[s2][0], ndu[rk][pk]);
			}
			j1 = (rk >= -1) ? 1 : -rk;
			j2 = ((int)r - 1 <= pk) ? (int)kk - 1 : (int)p - (int)r;
			for (jj = j1; jj <= j2; jj++)
			{
				a[s2][jj] = qaws_dual1_div(qaws_dual1_sub(a[s1][jj], a[s1][jj - 1]), ndu[pk + 1][rk + jj]);
				d = qaws_dual1_add(d, qaws_dual1_mul(a[s2][jj], ndu[rk + jj][pk]));
			}
			if ((int)r <= pk)
			{
				a[s2][kk] = qaws_dual1_div(qaws_dual1_make(-a[s1][kk - 1].v, -a[s1][kk - 1].t, -a[s1][kk - 1].tt),
					ndu[pk + 1][r]);
				d = qaws_dual1_add(d, qaws_dual1_mul(a[s2][kk], ndu[r][pk]));
			}
			out_ders[kk * stride + r] = d;
			jj = s1;
			s1 = s2;
			s2 = jj;
		}
	}

	{
		qaws_scalar f = (qaws_scalar)p;
		for (kk = 1; kk <= n; kk++)
		{
			for (j = 0; j <= p; j++)
			{
				out_ders[kk * stride + j].v *= f;
				out_ders[kk * stride + j].t *= f;
				out_ders[kk * stride + j].tt *= f;
			}
			f *= (qaws_scalar)(p - kk);
		}
	}
	for (kk = n + 1; kk <= k; kk++)
		for (j = 0; j <= p; j++)
			out_ders[kk * stride + j] = qaws_dual1_const(QAWS_ZERO);
}
