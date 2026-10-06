#include "qaws_internal_basis.h"
#include "../qaws_platform.h"
#include <string.h>
#include <math.h>

void qaws_internal_decasteljau(
	qaws_scalar const* control_points,
	unsigned int degree,
	unsigned int dim_count,
	qaws_scalar t,
	qaws_scalar* out_point)
{
	qaws_scalar work[192]; /* max degree 64, max dim 3 */
	unsigned int num_points = degree + 1;
	unsigned int total = num_points * dim_count;

	memcpy(work, control_points, total * sizeof(qaws_scalar));

	qaws_scalar one_minus_t = QAWS_ONE - t;

	for (unsigned int r = 1; r <= degree; r++) {
		for (unsigned int i = 0; i <= degree - r; i++) {
			for (unsigned int d = 0; d < dim_count; d++) {
				work[i * dim_count + d] =
					one_minus_t * work[i * dim_count + d] +
					t * work[(i + 1) * dim_count + d];
			}
		}
	}

	memcpy(out_point, work, dim_count * sizeof(qaws_scalar));
}

void qaws_internal_bezier_derivative_points(
	qaws_scalar const* control_points,
	unsigned int degree,
	unsigned int dim_count,
	qaws_scalar* out_deriv_points)
{
	qaws_scalar n = (qaws_scalar)degree;

	for (unsigned int i = 0; i < degree; i++) {
		for (unsigned int d = 0; d < dim_count; d++) {
			out_deriv_points[i * dim_count + d] =
				n * (control_points[(i + 1) * dim_count + d] -
				     control_points[i * dim_count + d]);
		}
	}
}

void qaws_internal_hermite_basis(
	qaws_scalar t,
	qaws_scalar* out_h,
	qaws_scalar* out_dh,
	qaws_scalar* out_d2h,
	qaws_scalar* out_d3h)
{
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;

	out_h[0] = QAWS_LITERAL(2.0) * t3 - QAWS_LITERAL(3.0) * t2 + QAWS_ONE;
	out_h[1] = t3 - QAWS_LITERAL(2.0) * t2 + t;
	out_h[2] = -QAWS_LITERAL(2.0) * t3 + QAWS_LITERAL(3.0) * t2;
	out_h[3] = t3 - t2;

	if (out_dh) {
		out_dh[0] = QAWS_LITERAL(6.0) * t2 - QAWS_LITERAL(6.0) * t;
		out_dh[1] = QAWS_LITERAL(3.0) * t2 - QAWS_LITERAL(4.0) * t + QAWS_ONE;
		out_dh[2] = -QAWS_LITERAL(6.0) * t2 + QAWS_LITERAL(6.0) * t;
		out_dh[3] = QAWS_LITERAL(3.0) * t2 - QAWS_LITERAL(2.0) * t;
	}

	if (out_d2h) {
		out_d2h[0] = QAWS_LITERAL(12.0) * t - QAWS_LITERAL(6.0);
		out_d2h[1] = QAWS_LITERAL(6.0) * t - QAWS_LITERAL(4.0);
		out_d2h[2] = -QAWS_LITERAL(12.0) * t + QAWS_LITERAL(6.0);
		out_d2h[3] = QAWS_LITERAL(6.0) * t - QAWS_LITERAL(2.0);
	}

	if (out_d3h) {
		out_d3h[0] = QAWS_LITERAL(12.0);
		out_d3h[1] = QAWS_LITERAL(6.0);
		out_d3h[2] = -QAWS_LITERAL(12.0);
		out_d3h[3] = QAWS_LITERAL(6.0);
	}
}

unsigned int qaws_internal_find_knot_span(
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int num_control_points,
	qaws_scalar t)
{
	unsigned int n = num_control_points - 1;

	(void)knot_count;

	/* at the domain ends, the nearest non-empty span (a knot repeated at the
	   end of an unclamped vector leaves empty spans there) */
	if (t >= knots[n + 1]) {
		unsigned int s = n;
		while (s > degree && knots[s] >= knots[n + 1])
			s--;
		return s;
	}
	if (t <= knots[degree]) {
		unsigned int s = degree;
		while (s < n && knots[s + 1] <= knots[degree])
			s++;
		return s;
	}

	unsigned int low = degree;
	unsigned int high = n;

	while (low < high) {
		unsigned int mid = (low + high) / 2;
		if (t < knots[mid]) {
			high = mid;
		} else if (t >= knots[mid + 1]) {
			low = mid + 1;
		} else {
			return mid;
		}
	}

	return low;
}

void qaws_internal_bspline_basis(
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int span,
	qaws_scalar t,
	qaws_scalar* out_basis)
{
	qaws_scalar left[32];
	qaws_scalar right[32];

	(void)knot_count;

	out_basis[0] = QAWS_ONE;

	for (unsigned int j = 1; j <= degree; j++) {
		left[j] = t - knots[span + 1 - j];
		right[j] = knots[span + j] - t;
		qaws_scalar saved = QAWS_ZERO;

		for (unsigned int r = 0; r < j; r++) {
			qaws_scalar denom = right[r + 1] + left[j - r];
			qaws_scalar temp = (denom != QAWS_ZERO) ? out_basis[r] / denom : QAWS_ZERO;
			out_basis[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}

		out_basis[j] = saved;
	}
}

void qaws_internal_bspline_basis_derivs(
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int span,
	qaws_scalar t,
	unsigned int k,
	qaws_scalar* out_ders)
{
	unsigned int p = degree;
	unsigned int order = p + 1;

	qaws_scalar ndu[1024]; /* max (32*32) */
	qaws_scalar left_arr[32];
	qaws_scalar right_arr[32];
	qaws_scalar a[2][32];

	(void)knot_count;

	ndu[0] = QAWS_ONE;

	for (unsigned int j = 1; j <= p; j++) {
		left_arr[j] = t - knots[span + 1 - j];
		right_arr[j] = knots[span + j] - t;
		qaws_scalar saved = QAWS_ZERO;

		for (unsigned int r = 0; r < j; r++) {
			ndu[j * order + r] = right_arr[r + 1] + left_arr[j - r];
			qaws_scalar temp = ndu[r * order + (j - 1)] / ndu[j * order + r];

			ndu[r * order + j] = saved + right_arr[r + 1] * temp;
			saved = left_arr[j - r] * temp;
		}

		ndu[j * order + j] = saved;
	}

	for (unsigned int j = 0; j <= p; j++) {
		out_ders[0 * order + j] = ndu[j * order + p];
	}

	for (unsigned int r = 0; r <= p; r++) {
		int s1 = 0;
		int s2 = 1;
		a[0][0] = QAWS_ONE;

		for (unsigned int d = 1; d <= k; d++) {
			qaws_scalar dd = QAWS_ZERO;
			int rj = (int)r - (int)d;
			int pk = (int)p - (int)d;

			if (rj >= 0) {
				a[s2][0] = a[s1][0] / ndu[(pk + 1) * order + rj];
				dd = a[s2][0] * ndu[rj * order + pk];
			}

			int j1 = (rj >= -1) ? 1 : (-(int)rj);
			int j2 = ((int)r - 1 <= pk) ? ((int)d - 1) : ((int)p - (int)r);

			for (int j = j1; j <= j2; j++) {
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[(pk + 1) * order + (rj + j)];
				dd += a[s2][j] * ndu[(rj + j) * order + pk];
			}

			if ((int)r <= pk) {
				a[s2][d] = -a[s1][d - 1] / ndu[(pk + 1) * order + r];
				dd += a[s2][d] * ndu[r * order + pk];
			}

			out_ders[d * order + r] = dd;

			{ int tmp = s1; s1 = s2; s2 = tmp; }
		}
	}

	qaws_scalar factor = (qaws_scalar)p;
	for (unsigned int d = 1; d <= k; d++) {
		for (unsigned int j = 0; j <= p; j++) {
			out_ders[d * order + j] *= factor;
		}
		factor *= (qaws_scalar)(p - d);
	}
}
