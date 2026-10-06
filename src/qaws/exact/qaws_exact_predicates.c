#include "../qaws_exact.h"
#include "qaws_exact_int.h"
#include <math.h>

/*
 * Filters follow Shewchuk, "Adaptive Precision Floating-Point Arithmetic
 * and Fast Robust Geometric Predicates" (1997): with eps = 2^-53 the f64
 * determinant differs from the exact one by at most err_A * (sum of the
 * magnitudes of its terms), err_A = (3 + 16 eps) eps for orient2d and
 * (7 + 56 eps) eps for orient3d, as long as nothing underflows or
 * overflows. Below QAWS_EXACT_FILTER_TINY, or when a term overflowed (not
 * finite), the filter is skipped. The exact path scales every input to the smallest
 * binary exponent among them: a positive common factor that leaves every
 * sign unchanged.
 */

#define QAWS_EXACT_EPS (1.0 / 9007199254740992.0)   /* 2^-53 */
#define QAWS_EXACT_FILTER_TINY 1e-280

static int finite_all(double const* v, int n)
{
	int i;
	for (i = 0; i < n; i++)
		if (!(v[i] - v[i] == 0.0))
			return 0;
	return 1;
}

/* Exact integers v_i 2^-emin for the n doubles. */
static qaws_status to_common(double const* v, int n, qaws_exact_int* out)
{
	int64_t m[16];
	int e[16], i, emin = 0, have = 0;
	for (i = 0; i < n; i++)
	{
		qaws_exact_split_double(v[i], &m[i], &e[i]);
		if (m[i] != 0 && (!have || e[i] < emin))
		{
			emin = e[i];
			have = 1;
		}
	}
	for (i = 0; i < n; i++)
	{
		qaws_status st;
		qaws_exact_int_from_i64(&out[i], m[i]);
		if (m[i] == 0)
			continue;
		st = qaws_exact_int_shl(&out[i], &out[i], (unsigned int)(e[i] - emin));
		if (st != QAWS_STATUS_OK)
			return st;
	}
	return QAWS_STATUS_OK;
}

static qaws_exact_sign sign_of(double v)
{
	return v > 0 ? QAWS_EXACT_POSITIVE : (v < 0 ? QAWS_EXACT_NEGATIVE : QAWS_EXACT_ZERO);
}

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

qaws_status qaws_exact_orient2d(double const a[2], double const b[2], double const c[2], qaws_exact_sign* out_sign,
	qaws_exact_path* out_path)
{
	double v[6], detleft, detright, det, detsum;
	qaws_exact_int x[6], acx, bcy, acy, bcx, l, r;
	qaws_status st;
	if (!a || !b || !c || !out_sign)
		return QAWS_STATUS_INVALID_ARGUMENT;
	v[0] = a[0]; v[1] = a[1]; v[2] = b[0]; v[3] = b[1]; v[4] = c[0]; v[5] = c[1];
	if (!finite_all(v, 6))
		return QAWS_STATUS_INVALID_ARGUMENT;
	detleft = (a[0] - c[0]) * (b[1] - c[1]);
	detright = (a[1] - c[1]) * (b[0] - c[0]);
	det = detleft - detright;
	detsum = fabs(detleft) + fabs(detright);
	if (detsum > QAWS_EXACT_FILTER_TINY && detsum - detsum == 0.0 &&
	    fabs(det) > (3.0 + 16.0 * QAWS_EXACT_EPS) * QAWS_EXACT_EPS * detsum)
	{
		*out_sign = sign_of(det);
		if (out_path) *out_path = QAWS_EXACT_PATH_FILTER;
		return QAWS_STATUS_OK;
	}
	TRY(to_common(v, 6, x));
	TRY(qaws_exact_int_sub(&acx, &x[0], &x[4]));
	TRY(qaws_exact_int_sub(&bcy, &x[3], &x[5]));
	TRY(qaws_exact_int_sub(&acy, &x[1], &x[5]));
	TRY(qaws_exact_int_sub(&bcx, &x[2], &x[4]));
	TRY(qaws_exact_int_mul(&l, &acx, &bcy));
	TRY(qaws_exact_int_mul(&r, &acy, &bcx));
	*out_sign = (qaws_exact_sign)qaws_exact_int_cmp(&l, &r);
	if (out_path) *out_path = QAWS_EXACT_PATH_EXACT;
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_orient3d(double const a[3], double const b[3], double const c[3], double const d[3],
	qaws_exact_sign* out_sign, qaws_exact_path* out_path)
{
	double v[12], adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz, det, permanent;
	qaws_exact_int x[12], D[9], t1, t2, m1, m2, m3, acc;
	qaws_status st;
	int i;
	if (!a || !b || !c || !d || !out_sign)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3; i++)
	{
		v[i] = a[i];
		v[3 + i] = b[i];
		v[6 + i] = c[i];
		v[9 + i] = d[i];
	}
	if (!finite_all(v, 12))
		return QAWS_STATUS_INVALID_ARGUMENT;
	adx = a[0] - d[0]; bdx = b[0] - d[0]; cdx = c[0] - d[0];
	ady = a[1] - d[1]; bdy = b[1] - d[1]; cdy = c[1] - d[1];
	adz = a[2] - d[2]; bdz = b[2] - d[2]; cdz = c[2] - d[2];
	det = adz * (bdx * cdy - cdx * bdy) + bdz * (cdx * ady - adx * cdy) + cdz * (adx * bdy - bdx * ady);
	permanent = (fabs(bdx * cdy) + fabs(cdx * bdy)) * fabs(adz) + (fabs(cdx * ady) + fabs(adx * cdy)) * fabs(bdz) +
		(fabs(adx * bdy) + fabs(bdx * ady)) * fabs(cdz);
	if (permanent > QAWS_EXACT_FILTER_TINY && permanent - permanent == 0.0 &&
	    fabs(det) > (7.0 + 56.0 * QAWS_EXACT_EPS) * QAWS_EXACT_EPS * permanent)
	{
		*out_sign = sign_of(det);
		if (out_path) *out_path = QAWS_EXACT_PATH_FILTER;
		return QAWS_STATUS_OK;
	}
	TRY(to_common(v, 12, x));
	/* D = (a - d, b - d, c - d) by components */
	for (i = 0; i < 3; i++)
	{
		TRY(qaws_exact_int_sub(&D[i], &x[i], &x[9 + i]));
		TRY(qaws_exact_int_sub(&D[3 + i], &x[3 + i], &x[9 + i]));
		TRY(qaws_exact_int_sub(&D[6 + i], &x[6 + i], &x[9 + i]));
	}
	/* adz (bdx cdy - cdx bdy) + bdz (cdx ady - adx cdy) + cdz (adx bdy - bdx ady) */
	TRY(qaws_exact_int_mul(&t1, &D[3], &D[7]));
	TRY(qaws_exact_int_mul(&t2, &D[6], &D[4]));
	TRY(qaws_exact_int_sub(&m1, &t1, &t2));
	TRY(qaws_exact_int_mul(&m1, &m1, &D[2]));
	TRY(qaws_exact_int_mul(&t1, &D[6], &D[1]));
	TRY(qaws_exact_int_mul(&t2, &D[0], &D[7]));
	TRY(qaws_exact_int_sub(&m2, &t1, &t2));
	TRY(qaws_exact_int_mul(&m2, &m2, &D[5]));
	TRY(qaws_exact_int_mul(&t1, &D[0], &D[4]));
	TRY(qaws_exact_int_mul(&t2, &D[3], &D[1]));
	TRY(qaws_exact_int_sub(&m3, &t1, &t2));
	TRY(qaws_exact_int_mul(&m3, &m3, &D[8]));
	TRY(qaws_exact_int_add(&acc, &m1, &m2));
	TRY(qaws_exact_int_add(&acc, &acc, &m3));
	*out_sign = (qaws_exact_sign)qaws_exact_int_sign(&acc);
	if (out_path) *out_path = QAWS_EXACT_PATH_EXACT;
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_compare_ratio(double a, double b, double c, double d, qaws_exact_sign* out_sign, qaws_exact_path* out_path)
{
	double v[4], ad, cb, det, sum;
	qaws_exact_int x[4], l, r;
	qaws_exact_sign den;
	qaws_status st;
	if (!out_sign)
		return QAWS_STATUS_INVALID_ARGUMENT;
	v[0] = a; v[1] = b; v[2] = c; v[3] = d;
	if (!finite_all(v, 4) || b == 0.0 || d == 0.0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* sign(a/b - c/d) = sign(a d - c b) sign(b d) */
	den = (qaws_exact_sign)(sign_of(b) * sign_of(d));
	ad = a * d;
	cb = c * b;
	det = ad - cb;
	sum = fabs(ad) + fabs(cb);
	if (sum > QAWS_EXACT_FILTER_TINY && sum - sum == 0.0 && fabs(det) > 3.0 * QAWS_EXACT_EPS * sum)
	{
		*out_sign = (qaws_exact_sign)(sign_of(det) * den);
		if (out_path) *out_path = QAWS_EXACT_PATH_FILTER;
		return QAWS_STATUS_OK;
	}
	TRY(to_common(v, 4, x));
	TRY(qaws_exact_int_mul(&l, &x[0], &x[3]));
	TRY(qaws_exact_int_mul(&r, &x[2], &x[1]));
	*out_sign = (qaws_exact_sign)(qaws_exact_int_cmp(&l, &r) * den);
	if (out_path) *out_path = QAWS_EXACT_PATH_EXACT;
	return QAWS_STATUS_OK;
}
