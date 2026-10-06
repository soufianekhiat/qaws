#include "qaws_exact_curve.h"
#include "../qaws_inspect.h"
#include "../qaws_diff.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define QAWS_EXACT_MAX_KNOTS 256

void qaws_exact_desc_default(qaws_exact_desc* desc)
{
	desc->space_exp2 = -20;
	desc->coord_bits = 26;
	desc->param_bits = 24;
	desc->weight_bits = 24;
}

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

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

/* ------------------------------------------------------------------ */
/*  Fractions of homogeneous points: num[0..D) / den, den > 0          */
/* ------------------------------------------------------------------ */

typedef struct hfrac
{
	qaws_exact_int num[4];
	qaws_exact_int den;
} hfrac;

/* Divides every numerator and the denominator by their gcd. */
static qaws_status hfrac_reduce(hfrac* f, unsigned int D)
{
	qaws_exact_int g;
	unsigned int c;
	qaws_status st;
	g = f->den;
	for (c = 0; c < D; c++)
		qaws_exact_int_gcd(&g, &g, &f->num[c]);
	if (qaws_exact_int_bits(&g) <= 1)
		return QAWS_STATUS_OK;
	for (c = 0; c < D; c++)
		TRY(qaws_exact_int_divmod(&f->num[c], NULL, &f->num[c], &g));
	TRY(qaws_exact_int_divmod(&f->den, NULL, &f->den, &g));
	return QAWS_STATUS_OK;
}

/* r = ((dd - nu) x + nu y) / dd for integers nu, dd > 0 */
static qaws_status hfrac_lerp(hfrac* r, hfrac const* x, hfrac const* y, int64_t nu, int64_t dd, unsigned int D)
{
	qaws_exact_int t1, t2, den;
	hfrac out;
	unsigned int c;
	qaws_status st;
	for (c = 0; c < D; c++)
	{
		TRY(qaws_exact_int_mul(&t1, &x->num[c], &y->den));
		TRY(qaws_exact_int_mul_i64(&t1, &t1, dd - nu));
		TRY(qaws_exact_int_mul(&t2, &y->num[c], &x->den));
		TRY(qaws_exact_int_mul_i64(&t2, &t2, nu));
		TRY(qaws_exact_int_add(&out.num[c], &t1, &t2));
	}
	TRY(qaws_exact_int_mul(&den, &x->den, &y->den));
	TRY(qaws_exact_int_mul_i64(&out.den, &den, dd));
	TRY(hfrac_reduce(&out, D));
	*r = out;
	return QAWS_STATUS_OK;
}

/*
 * Blossom of a degree-p B-spline piece on span s:
 * f(a^(p - j), b^j) by the de Boor triangle with the argument of level r
 * equal to a for the first p - j levels and b after them:
 *   d_i^r = (1 - alpha) d_(i-1)^(r-1) + alpha d_i^(r-1),
 *   alpha = (t_r - u_i) / (u_(i+p+1-r) - u_i),  i = s - p + r .. s.
 * These are the Bezier control points of the span (j = 0..p).
 */
static qaws_status blossom(hfrac const* local, int64_t const* K, unsigned int s, unsigned int p, unsigned int j, unsigned int D,
	hfrac* out)
{
	hfrac d[QAWS_EXACT_MAX_DEGREE + 1];
	unsigned int r, i;
	qaws_status st;
	for (i = 0; i <= p; i++)
		d[i] = local[i];
	for (r = 1; r <= p; r++)
	{
		int64_t t = r <= p - j ? K[s] : K[s + 1];
		for (i = p; i >= r; i--)
		{
			int64_t lo = K[s - p + i], hi = K[s - p + i + p + 1 - r];
			TRY(hfrac_lerp(&d[i], &d[i - 1], &d[i], t - lo, hi - lo, D));
		}
	}
	*out = d[p];
	return QAWS_STATUS_OK;
}

/* Clears the denominators of the span's Bezier points with their lcm. */
static qaws_status clear_denominators(hfrac* b, unsigned int p, unsigned int D, qaws_exact_int* h)
{
	qaws_exact_int lcm, g, q;
	unsigned int j, c;
	qaws_status st;
	lcm = b[0].den;
	for (j = 1; j <= p; j++)
	{
		qaws_exact_int_gcd(&g, &lcm, &b[j].den);
		TRY(qaws_exact_int_divmod(&q, NULL, &b[j].den, &g));
		TRY(qaws_exact_int_mul(&lcm, &lcm, &q));
	}
	for (j = 0; j <= p; j++)
	{
		TRY(qaws_exact_int_divmod(&q, NULL, &lcm, &b[j].den));
		for (c = 0; c < D; c++)
			TRY(qaws_exact_int_mul(&h[j * D + c], &b[j].num[c], &q));
	}
	/* a common factor of every entry is free in homogeneous coordinates */
	qaws_exact_int_zero(&g);
	for (j = 0; j < (p + 1) * D; j++)
		qaws_exact_int_gcd(&g, &g, &h[j]);
	if (qaws_exact_int_bits(&g) > 1)
		for (j = 0; j < (p + 1) * D; j++)
			TRY(qaws_exact_int_divmod(&h[j], NULL, &h[j], &g));
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Preparation                                                        */
/* ------------------------------------------------------------------ */

static void exact_curve_free(qaws_exact_curve* ec)
{
	unsigned int s;
	if (!ec)
		return;
	for (s = 0; s < ec->span_count && ec->spans; s++)
		qaws_internal_dealloc(NULL, ec->spans[s].h);
	qaws_internal_dealloc(NULL, ec->spans);
	qaws_internal_dealloc(NULL, ec);
}

void qaws_exact_curve_destroy(qaws_exact_curve* curve)
{
	exact_curve_free(curve);
}

qaws_status qaws_exact_curve_prepare(qaws_exact_desc const* desc, qaws_curve const* curve, qaws_exact_curve** out_curve,
	qaws_exact_report* out_report)
{
	qaws_exact_desc d;
	qaws_curve_kind kind;
	qaws_field_desc fields[8];
	qaws_scalar cps[QAWS_EXACT_MAX_KNOTS * 3], ws[QAWS_EXACT_MAX_KNOTS], kn[QAWS_EXACT_MAX_KNOTS];
	int64_t K[QAWS_EXACT_MAX_KNOTS];
	unsigned int nf = 0, f, n = 0, nk = 0, got = 0, i, c, dim, D, p, storage = 0, s;
	int rational, spline, wexp = 0;
	double max_pos = 0, max_w = 0, max_k = 0, wmax = 0;
	hfrac* local = NULL;
	qaws_exact_curve* ec;
	qaws_status st = QAWS_STATUS_OK;
	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if (desc)
		d = *desc;
	else
		qaws_exact_desc_default(&d);
	if (d.param_bits == 0 || d.param_bits > 56 || d.coord_bits == 0 || d.coord_bits > 32 || d.weight_bits == 0 || d.weight_bits > 30)
		return QAWS_STATUS_INVALID_ARGUMENT;
	kind = qaws_curve_get_kind(curve);
	rational = kind == QAWS_CURVE_KIND_RATIONAL_BEZIER || kind == QAWS_CURVE_KIND_NURBS;
	spline = kind == QAWS_CURVE_KIND_BSPLINE || kind == QAWS_CURVE_KIND_NURBS;
	if (kind != QAWS_CURVE_KIND_BEZIER && kind != QAWS_CURVE_KIND_RATIONAL_BEZIER && !spline)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	if (qaws_curve_describe_fields(curve, fields, 8, &nf) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (f = 0; f < nf && f < 8; f++)
	{
		if (fields[f].field == QAWS_FIELD_CONTROL_POINTS) n = fields[f].count;
		if (fields[f].field == QAWS_FIELD_KNOTS) nk = fields[f].count;
	}
	dim = (unsigned int)curve->dimension;
	D = dim + 1;
	p = spline ? curve->degree : n - 1;
	if (n < 1 || n > QAWS_EXACT_MAX_KNOTS || p > QAWS_EXACT_MAX_DEGREE || (dim != 2 && dim != 3) ||
	    (spline && (nk != n + p + 1 || nk > QAWS_EXACT_MAX_KNOTS)))
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
	frexp(wmax, &wexp);
	wexp = (int)d.weight_bits - wexp;

	ec = (qaws_exact_curve*)qaws_internal_alloc(NULL, (unsigned long)sizeof(qaws_exact_curve));
	if (!ec)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(ec, 0, sizeof(*ec));
	ec->dimension = (int)dim;
	ec->space_exp2 = d.space_exp2;
	ec->desc = d;

	/* parameter lattice: t = T 2^-shift with |T| <= 2^param_bits */
	{
		double tmax = 1;
		int e = 0;
		if (spline)
		{
			if (qaws_curve_read_field(curve, QAWS_FIELD_KNOTS, kn, nk, &got) != QAWS_STATUS_OK)
			{
				exact_curve_free(ec);
				return QAWS_STATUS_EXACT_UNSUPPORTED;
			}
			tmax = 0;
			for (i = 0; i < nk; i++)
				if (fabs((double)kn[i]) > tmax)
					tmax = fabs((double)kn[i]);
			if (tmax == 0)
				tmax = 1;
		}
		frexp(tmax, &e);
		ec->param_shift = (int)d.param_bits + 1 - e;
		if (spline)
			for (i = 0; i < nk; i++)
			{
				double err;
				if (quantize((double)kn[i], -ec->param_shift, d.param_bits + 1, &K[i], &err) != QAWS_STATUS_OK)
				{
					exact_curve_free(ec);
					return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
				}
				if (err > max_k) max_k = err;
			}
	}

	/* homogeneous integer control points */
	local = (hfrac*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(hfrac) * n));
	if (!local)
	{
		exact_curve_free(ec);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
	{
		int64_t W = 1, X;
		double err;
		if (rational)
		{
			double r = nearbyint(ldexp((double)ws[i], wexp));
			if (!(r >= 1))
				st = QAWS_STATUS_EXACT_RANGE_EXCEEDED;
			W = (int64_t)r;
			err = fabs(ldexp(r, -wexp) - (double)ws[i]) / (double)ws[i];
			if (err > max_w) max_w = err;
		}
		for (c = 0; c < dim && st == QAWS_STATUS_OK; c++)
		{
			st = quantize((double)cps[i * dim + c], d.space_exp2, d.coord_bits, &X, &err);
			if (err > max_pos) max_pos = err;
			qaws_exact_int_from_i64(&local[i].num[c], W * X);
		}
		qaws_exact_int_from_i64(&local[i].num[dim], W);
		qaws_exact_int_from_i64(&local[i].den, 1);
	}

	/* spans */
	if (st == QAWS_STATUS_OK)
	{
		unsigned int count = 0;
		if (spline)
		{
			for (s = p; s < n; s++)
				if (K[s] < K[s + 1])
					count++;
		}
		else
			count = 1;
		ec->spans = (qaws_exact_span*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_span) * (count ? count : 1)));
		if (!ec->spans || count == 0)
			st = ec->spans ? QAWS_STATUS_EXACT_UNSUPPORTED : QAWS_STATUS_ALLOCATION_FAILURE;
		else
			memset(ec->spans, 0, sizeof(qaws_exact_span) * count);
		ec->span_count = 0;
		for (s = spline ? p : 0; st == QAWS_STATUS_OK && (spline ? s < n : s < 1); s++)
		{
			qaws_exact_span* sp;
			if (spline && !(K[s] < K[s + 1]))
				continue;
			sp = &ec->spans[ec->span_count++];
			sp->degree = p;
			sp->h = (qaws_exact_int*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_int) * (p + 1) * D));
			if (!sp->h)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			if (!spline)
			{
				sp->a = 0;
				sp->b = (int64_t)1 << ec->param_shift;
				for (i = 0; i <= p; i++)
					for (c = 0; c < D; c++)
						sp->h[i * D + c] = local[i].num[c];
			}
			else
			{
				hfrac b[QAWS_EXACT_MAX_DEGREE + 1];
				unsigned int j;
				sp->a = K[s];
				sp->b = K[s + 1];
				for (j = 0; j <= p && st == QAWS_STATUS_OK; j++)
					st = blossom(&local[s - p], K, s, p, j, D, &b[j]);
				if (st == QAWS_STATUS_OK)
					st = clear_denominators(b, p, D, sp->h);
			}
			for (i = 0; i < (p + 1) * D && st == QAWS_STATUS_OK; i++)
				if (qaws_exact_int_bits(&sp->h[i]) > storage)
					storage = qaws_exact_int_bits(&sp->h[i]);
		}
	}
	qaws_internal_dealloc(NULL, local);
	if (st != QAWS_STATUS_OK)
	{
		exact_curve_free(ec);
		return st;
	}
	if (out_report)
	{
		memset(out_report, 0, sizeof(*out_report));
		out_report->quality = QAWS_NUMERIC_EXACT_RATIONAL;
		out_report->flags = (max_pos > 0 || max_w > 0 || max_k > 0) ? QAWS_EXACT_FLAG_INPUT_QUANTIZED : QAWS_EXACT_FLAG_NONE;
		out_report->storage_bits = storage;
		out_report->max_position_quantization_error = max_pos;
		out_report->max_weight_quantization_error = max_w;
		out_report->parameter_quantization_error = max_k;
	}
	*out_curve = ec;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Evaluation                                                         */
/* ------------------------------------------------------------------ */

/* The span holding T: a <= T < b, the last span at the domain end. */
static qaws_exact_span const* find_span(qaws_exact_curve const* ec, int64_t* T)
{
	unsigned int s;
	if (*T < ec->spans[0].a) *T = ec->spans[0].a;
	if (*T > ec->spans[ec->span_count - 1].b) *T = ec->spans[ec->span_count - 1].b;
	for (s = 0; s + 1 < ec->span_count; s++)
		if (*T < ec->spans[s].b)
			return &ec->spans[s];
	return &ec->spans[ec->span_count - 1];
}

/*
 * d^j H / ds^j at s = x / L (x = T - a, L = b - a), scaled by L^degree:
 * the j-th derivative polygon n (n-1) ... (n-j+1) Delta^j H, then
 * division-free De Casteljau Q_i <- (L - x) Q_i + x Q_(i+1), which leaves
 * L^(n-j); a factor L^j brings it to L^n.
 */
static qaws_status homogeneous_derivative(qaws_exact_span const* sp, unsigned int D, int64_t x, unsigned int j, qaws_exact_int* out)
{
	qaws_exact_int pts[(QAWS_EXACT_MAX_DEGREE + 1) * 4], t1, t2;
	unsigned int n = sp->degree, i, r, c, m;
	int64_t L = sp->b - sp->a;
	qaws_status st;
	if (j > n)
	{
		for (c = 0; c < D; c++)
			qaws_exact_int_zero(&out[c]);
		return QAWS_STATUS_OK;
	}
	for (i = 0; i < (n + 1) * D; i++)
		pts[i] = sp->h[i];
	for (r = 0; r < j; r++)
		for (i = 0; i + r < n; i++)
			for (c = 0; c < D; c++)
			{
				TRY(qaws_exact_int_sub(&pts[i * D + c], &pts[(i + 1) * D + c], &pts[i * D + c]));
				TRY(qaws_exact_int_mul_i64(&pts[i * D + c], &pts[i * D + c], (int64_t)(n - r)));
			}
	m = n - j;
	for (r = 1; r <= m; r++)
		for (i = 0; i + r <= m; i++)
			for (c = 0; c < D; c++)
			{
				TRY(qaws_exact_int_mul_i64(&t1, &pts[i * D + c], L - x));
				TRY(qaws_exact_int_mul_i64(&t2, &pts[(i + 1) * D + c], x));
				TRY(qaws_exact_int_add(&pts[i * D + c], &t1, &t2));
			}
	for (c = 0; c < D; c++)
	{
		out[c] = pts[c];
		for (r = 0; r < j; r++)
			TRY(qaws_exact_int_mul_i64(&out[c], &out[c], L));
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_curve_eval_rational(qaws_exact_curve const* curve, int64_t T, unsigned int k, qaws_exact_int* num,
	qaws_exact_int* den)
{
	/* H[j] = d^j H / ds^j (common scale), N[j][c] = numerators of C^(j) over W^(j+1) */
	qaws_exact_int H[4][4], N[4][3], Wp[5], t1, t2;
	static int const binom[4][4] = { { 1, 0, 0, 0 }, { 1, 1, 0, 0 }, { 1, 2, 1, 0 }, { 1, 3, 3, 1 } };
	qaws_exact_span const* sp;
	unsigned int dim, D, j, p, c, i;
	qaws_status st;
	if (!curve || !num || !den || k > 3)
		return QAWS_STATUS_INVALID_ARGUMENT;
	dim = (unsigned int)curve->dimension;
	D = dim + 1;
	sp = find_span(curve, &T);
	for (j = 0; j <= k; j++)
		TRY(homogeneous_derivative(sp, D, T - sp->a, j, H[j]));
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
	*den = Wp[k + 1];
	for (c = 0; c < dim; c++)
		num[c] = N[k][c];
	/* d/dt = (2^shift / L) d/ds */
	for (j = 0; j < k; j++)
	{
		for (c = 0; c < dim; c++)
			TRY(qaws_exact_int_shl(&num[c], &num[c], (unsigned int)curve->param_shift));
		TRY(qaws_exact_int_mul_i64(den, den, sp->b - sp->a));
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_curve_evaluate(qaws_exact_curve const* curve, double t, unsigned int order, double* out,
	qaws_exact_report* out_report)
{
	double scale, r;
	int64_t T;
	unsigned int j, c;
	qaws_status st;
	if (!curve || !out || order > 3 || !(t - t == 0.0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	r = nearbyint(ldexp(t, curve->param_shift));
	if (r < -9.0e18) r = -9.0e18;
	if (r > 9.0e18) r = 9.0e18;
	T = (int64_t)r;
	(void)find_span(curve, &T);
	if (out_report)
	{
		out_report->quality = QAWS_NUMERIC_EXACT_RATIONAL;
		out_report->parameter_quantization_error = fabs(t - ldexp((double)T, -curve->param_shift));
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
