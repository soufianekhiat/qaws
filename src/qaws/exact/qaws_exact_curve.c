#include "qaws_exact_curve.h"
#include "../qaws_inspect.h"
#include "../qaws_diff.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

void qaws_exact_desc_default(qaws_exact_desc* desc)
{
	desc->space_exp2 = -20;
	desc->coord_bits = 26;
	desc->param_bits = 24;
	desc->weight_bits = 24;
}

/* x -> nearest lattice integer (ties to even) and its rounding error. */
static qaws_status quantize(double x, int exp2, unsigned int bits, int64_t* out, double* err)
{
	double v = ldexp(x, -exp2), r = nearbyint(v);
	if (!(fabs(r) < ldexp(1.0, (int)bits)))
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	*out = (int64_t)r;
	*err = fabs(x - ldexp(r, exp2));
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_curve_prepare(qaws_exact_desc const* desc, qaws_curve const* curve, qaws_exact_curve** out_curve,
	qaws_exact_report* out_report)
{
	qaws_exact_desc d;
	qaws_curve_kind kind;
	qaws_field_desc fields[8];
	qaws_scalar cps[(QAWS_EXACT_MAX_DEGREE + 1) * 3], ws[QAWS_EXACT_MAX_DEGREE + 1];
	unsigned int nf = 0, f, n = 0, got = 0, i, c, dim = 0, storage = 0;
	int rational, wexp = 0;
	double max_pos = 0, max_w = 0, wmax = 0;
	qaws_exact_curve* ec;
	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if (desc)
		d = *desc;
	else
		qaws_exact_desc_default(&d);
	if (d.param_bits == 0 || d.param_bits > 60 || d.coord_bits == 0 || d.coord_bits > 32 || d.weight_bits == 0 || d.weight_bits > 30)
		return QAWS_STATUS_INVALID_ARGUMENT;
	kind = qaws_curve_get_kind(curve);
	if (kind != QAWS_CURVE_KIND_BEZIER && kind != QAWS_CURVE_KIND_RATIONAL_BEZIER)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	rational = kind == QAWS_CURVE_KIND_RATIONAL_BEZIER;
	if (qaws_curve_describe_fields(curve, fields, 8, &nf) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (f = 0; f < nf && f < 8; f++)
		if (fields[f].field == QAWS_FIELD_CONTROL_POINTS)
			n = fields[f].count;
	dim = (unsigned int)curve->dimension;
	if (n < 1 || n > QAWS_EXACT_MAX_DEGREE + 1 || (dim != 2 && dim != 3))
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	if (qaws_curve_read_field(curve, QAWS_FIELD_CONTROL_POINTS, cps, n * dim, &got) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (i = 0; i < n; i++)
		ws[i] = 1;
	if (rational && qaws_curve_read_field(curve, QAWS_FIELD_WEIGHTS, ws, n, &got) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (i = 0; i < n; i++)
	{
		if (!(ws[i] > 0))
			return QAWS_STATUS_EXACT_UNSUPPORTED;
		if (ws[i] > wmax)
			wmax = ws[i];
	}
	/* a common power of two brings the largest weight below 2^weight_bits */
	frexp(wmax, &wexp);
	wexp = (int)d.weight_bits - wexp;

	ec = (qaws_exact_curve*)qaws_internal_alloc(NULL, (unsigned long)sizeof(qaws_exact_curve));
	if (!ec)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(ec, 0, sizeof(*ec));
	ec->h = (qaws_exact_int*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_int) * n * (dim + 1)));
	if (!ec->h)
	{
		qaws_internal_dealloc(NULL, ec);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	ec->dimension = (int)dim;
	ec->degree = n - 1;
	ec->param_bits = d.param_bits;
	ec->space_exp2 = d.space_exp2;
	ec->desc = d;
	for (i = 0; i < n; i++)
	{
		int64_t W = 1, X;
		double err;
		qaws_status st;
		if (rational)
		{
			double r = nearbyint(ldexp((double)ws[i], wexp));
			if (!(r >= 1))
			{
				qaws_exact_curve_destroy(ec);
				return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
			}
			W = (int64_t)r;
			err = fabs(ldexp(r, -wexp) - (double)ws[i]) / (double)ws[i];
			if (err > max_w) max_w = err;
		}
		for (c = 0; c < dim; c++)
		{
			st = quantize((double)cps[i * dim + c], d.space_exp2, d.coord_bits, &X, &err);
			if (st != QAWS_STATUS_OK)
			{
				qaws_exact_curve_destroy(ec);
				return st;
			}
			if (err > max_pos) max_pos = err;
			qaws_exact_int_from_i64(&ec->h[i * (dim + 1) + c], W * X);
			if (qaws_exact_int_bits(&ec->h[i * (dim + 1) + c]) > storage)
				storage = qaws_exact_int_bits(&ec->h[i * (dim + 1) + c]);
		}
		qaws_exact_int_from_i64(&ec->h[i * (dim + 1) + dim], W);
	}
	if (out_report)
	{
		memset(out_report, 0, sizeof(*out_report));
		out_report->quality = QAWS_NUMERIC_EXACT_RATIONAL;
		out_report->flags = (max_pos > 0 || max_w > 0) ? QAWS_EXACT_FLAG_INPUT_QUANTIZED : QAWS_EXACT_FLAG_NONE;
		out_report->storage_bits = storage;
		out_report->max_position_quantization_error = max_pos;
		out_report->max_weight_quantization_error = max_w;
	}
	*out_curve = ec;
	return QAWS_STATUS_OK;
}

void qaws_exact_curve_destroy(qaws_exact_curve* curve)
{
	if (!curve)
		return;
	qaws_internal_dealloc(NULL, curve->h);
	qaws_internal_dealloc(NULL, curve);
}

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

/*
 * d^j H / ds^j at s = a / b, b = 2^P, scaled by b^degree: the j-th
 * derivative polygon n (n-1) ... (n-j+1) Delta^j H, then division-free De
 * Casteljau Q_i <- (b - a) Q_i + a Q_(i+1), which leaves the common factor
 * b^(n-j); a shift by j P bits brings it to b^n.
 */
static qaws_status homogeneous_derivative(qaws_exact_curve const* ec, int64_t a, unsigned int j, qaws_exact_int* out)
{
	qaws_exact_int pts[(QAWS_EXACT_MAX_DEGREE + 1) * 4], t1, t2;
	unsigned int n = ec->degree, D = (unsigned int)ec->dimension + 1, i, r, c, m;
	int64_t b = (int64_t)1 << ec->param_bits;
	qaws_status st;
	if (j > n)
	{
		for (c = 0; c < D; c++)
			qaws_exact_int_zero(&out[c]);
		return QAWS_STATUS_OK;
	}
	for (i = 0; i < (n + 1) * D; i++)
		pts[i] = ec->h[i];
	/* derivative polygon */
	for (r = 0; r < j; r++)
		for (i = 0; i + r < n; i++)
			for (c = 0; c < D; c++)
			{
				TRY(qaws_exact_int_sub(&pts[i * D + c], &pts[(i + 1) * D + c], &pts[i * D + c]));
				TRY(qaws_exact_int_mul_i64(&pts[i * D + c], &pts[i * D + c], (int64_t)(n - r)));
			}
	m = n - j;
	/* scaled De Casteljau */
	for (r = 1; r <= m; r++)
		for (i = 0; i + r <= m; i++)
			for (c = 0; c < D; c++)
			{
				TRY(qaws_exact_int_mul_i64(&t1, &pts[i * D + c], b - a));
				TRY(qaws_exact_int_mul_i64(&t2, &pts[(i + 1) * D + c], a));
				TRY(qaws_exact_int_add(&pts[i * D + c], &t1, &t2));
			}
	for (c = 0; c < D; c++)
		TRY(qaws_exact_int_shl(&out[c], &pts[c], j * ec->param_bits));
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_curve_eval_rational(qaws_exact_curve const* curve, int64_t T, unsigned int k, qaws_exact_int* num,
	qaws_exact_int* den)
{
	/* H[j] = d^j H / ds^j (common scale), N[j][c] = numerators of C^(j) over W^(j+1) */
	qaws_exact_int H[4][4], N[4][3], Wp[5], t1, t2;
	static int const binom[4][4] = { { 1, 0, 0, 0 }, { 1, 1, 0, 0 }, { 1, 2, 1, 0 }, { 1, 3, 3, 1 } };
	unsigned int dim, j, p, c, i;
	qaws_status st;
	if (!curve || !num || !den || k > 3)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (T < 0 || T > ((int64_t)1 << curve->param_bits))
		return QAWS_STATUS_OUT_OF_RANGE;
	dim = (unsigned int)curve->dimension;
	for (j = 0; j <= k; j++)
		TRY(homogeneous_derivative(curve, T, j, H[j]));
	/* powers of W = H[0][dim] */
	qaws_exact_int_from_i64(&Wp[0], 1);
	for (p = 1; p <= k + 1; p++)
		TRY(qaws_exact_int_mul(&Wp[p], &Wp[p - 1], &H[0][dim]));
	/* N_j = X^(j) W^j - sum_{i<j} C(j, i) N_i W^(j-i) W^(j-i-1) */
	for (j = 0; j <= k; j++)
		for (c = 0; c < dim; c++)
		{
			TRY(qaws_exact_int_mul(&N[j][c], &H[j][c], &Wp[j]));
			for (i = 0; i < j; i++)
			{
				TRY(qaws_exact_int_mul(&t1, &N[i][c], &H[j - i][dim]));
				TRY(qaws_exact_int_mul(&t1, &t1, &Wp[j - i - 1]));
				TRY(qaws_exact_int_mul_i64(&t2, &t1, binom[j][i]));
				TRY(qaws_exact_int_sub(&N[j][c], &N[j][c], &t2));
			}
		}
	for (c = 0; c < dim; c++)
		num[c] = N[k][c];
	*den = Wp[k + 1];
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_curve_evaluate(qaws_exact_curve const* curve, double t, unsigned int order, double* out,
	qaws_exact_report* out_report)
{
	double scale, r;
	int64_t T, top;
	unsigned int j, c;
	qaws_status st;
	if (!curve || !out || order > 3 || !(t - t == 0.0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	top = (int64_t)1 << curve->param_bits;
	r = nearbyint(ldexp(t, (int)curve->param_bits));
	T = r < 0 ? 0 : (r > (double)top ? top : (int64_t)r);
	if (out_report)
	{
		out_report->quality = QAWS_NUMERIC_EXACT_RATIONAL;
		out_report->parameter_quantization_error = fabs(t - ldexp((double)T, -(int)curve->param_bits));
	}
	scale = ldexp(1.0, curve->space_exp2);
	for (j = 0; j <= order; j++)
	{
		qaws_exact_int num[3], den;
		TRY(qaws_exact_curve_eval_rational(curve, T, j, num, &den));
		for (c = 0; c < (unsigned int)curve->dimension; c++)
			out[j * (unsigned int)curve->dimension + c] = qaws_exact_ratio_to_double(&num[c], &den) * scale;
	}
	return QAWS_STATUS_OK;
}
