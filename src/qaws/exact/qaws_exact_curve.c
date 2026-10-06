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

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

/* x -> nearest lattice integer (ties to even) and its rounding error. */
qaws_status qaws_exact_quantize(double x, int exp2, unsigned int bits, int64_t* out, double* err)
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

/* Divides every numerator and the denominator by their gcd. */
static qaws_status hfrac_reduce(qaws_exact_hfrac* f, unsigned int D)
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
static qaws_status hfrac_lerp(qaws_exact_hfrac* r, qaws_exact_hfrac const* x, qaws_exact_hfrac const* y, int64_t nu, int64_t dd, unsigned int D)
{
	qaws_exact_int t1, t2, den;
	qaws_exact_hfrac out;
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
qaws_status qaws_exact_blossom(qaws_exact_hfrac const* local, int64_t const* K, unsigned int s, unsigned int p, unsigned int j, unsigned int D,
	qaws_exact_hfrac* out)
{
	qaws_exact_hfrac d[QAWS_EXACT_MAX_DEGREE + 1];
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
qaws_status qaws_exact_clear_denominators(qaws_exact_hfrac* b, unsigned int p, unsigned int D, qaws_exact_int* h)
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

/* Allocates an exact curve of `count` spans of degree p (dimension dim). */
static qaws_exact_curve* exact_curve_alloc(qaws_exact_desc const* d, unsigned int dim, unsigned int count, unsigned int p)
{
	qaws_exact_curve* ec = (qaws_exact_curve*)qaws_internal_alloc(NULL, (unsigned long)sizeof(qaws_exact_curve));
	unsigned int s;
	if (!ec)
		return NULL;
	memset(ec, 0, sizeof(*ec));
	ec->dimension = (int)dim;
	ec->space_exp2 = d->space_exp2;
	ec->desc = *d;
	ec->spans = (qaws_exact_span*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_span) * count));
	if (!ec->spans)
	{
		qaws_internal_dealloc(NULL, ec);
		return NULL;
	}
	memset(ec->spans, 0, sizeof(qaws_exact_span) * count);
	ec->span_count = count;
	for (s = 0; s < count; s++)
	{
		ec->spans[s].degree = p;
		ec->spans[s].h = (qaws_exact_int*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_int) * (p + 1) * (dim + 1)));
		if (!ec->spans[s].h)
		{
			exact_curve_free(ec);
			return NULL;
		}
	}
	return ec;
}

/* Parameter shift of a domain whose largest magnitude is tmax. */
int qaws_exact_param_shift_for(double tmax, unsigned int param_bits)
{
	int e = 0;
	frexp(tmax > 0 ? tmax : 1.0, &e);
	return (int)param_bits + 1 - e;
}

/*
 * Hermite (cubic, unit spans [i, i + 1]) and uniform Catmull-Rom (unit
 * spans) become cubic Beziers with one integer common factor:
 *   Hermite:      3 (P0, P0 + M0/3, P1 - M1/3, P1) = (3 P0, 3 P0 + M0, 3 P1 - M1, 3 P1), W = 3
 *   Catmull-Rom:  6 (P1, P1 + (P2 - P0)/6, P2 - (P3 - P1)/6, P2)
 *               = (6 P1, 6 P1 + P2 - P0, 6 P2 - P3 + P1, 6 P2), W = 6
 * Points and tangents are quantized onto the space lattice (tangents per
 * unit parameter).
 */
static qaws_status prepare_cubic_family(qaws_exact_desc const* d, qaws_curve const* curve, qaws_curve_kind kind, qaws_exact_curve** out,
	double* max_pos)
{
	qaws_field_desc fields[8];
	qaws_scalar pts[QAWS_EXACT_MAX_KNOTS * 3], tan[QAWS_EXACT_MAX_KNOTS * 3];
	int64_t P[QAWS_EXACT_MAX_KNOTS * 3], M[QAWS_EXACT_MAX_KNOTS * 3];
	unsigned int nf = 0, f, n = 0, got = 0, dim = (unsigned int)curve->dimension, D = dim + 1, spans, s, i, c;
	int closed = 0;
	qaws_exact_curve* ec;
	qaws_status st;
	if (kind == QAWS_CURVE_KIND_CATMULL_ROM)
	{
		qaws_catmull_rom_impl const* impl = (qaws_catmull_rom_impl const*)curve->impl;
		if (impl->parameterization != QAWS_PARAMETERIZATION_UNIFORM)
			return QAWS_STATUS_EXACT_UNSUPPORTED;   /* distances: sqrt */
		closed = impl->closed;
	}
	if (qaws_curve_describe_fields(curve, fields, 8, &nf) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (f = 0; f < nf && f < 8; f++)
		if (fields[f].field == QAWS_FIELD_POINTS)
			n = fields[f].count;
	if (n < 2 || n > QAWS_EXACT_MAX_KNOTS || (dim != 2 && dim != 3))
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	if (qaws_curve_read_field(curve, QAWS_FIELD_POINTS, pts, n * dim, &got) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	if (kind == QAWS_CURVE_KIND_HERMITE &&
	    qaws_curve_read_field(curve, QAWS_FIELD_DERIVATIVES, tan, n * dim, &got) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (i = 0; i < n * dim; i++)
	{
		double err;
		TRY(qaws_exact_quantize((double)pts[i], d->space_exp2, d->coord_bits, &P[i], &err));
		if (err > *max_pos) *max_pos = err;
		if (kind == QAWS_CURVE_KIND_HERMITE)
		{
			TRY(qaws_exact_quantize((double)tan[i], d->space_exp2, d->coord_bits, &M[i], &err));
			if (err > *max_pos) *max_pos = err;
		}
	}
	spans = kind == QAWS_CURVE_KIND_HERMITE ? n - 1 : (closed ? n : (n >= 4 ? n - 3 : 0));
	if (spans == 0)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	ec = exact_curve_alloc(d, dim, spans, 3);
	if (!ec)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	ec->param_shift = qaws_exact_param_shift_for((double)spans, d->param_bits);
	for (s = 0; s < spans; s++)
	{
		qaws_exact_span* sp = &ec->spans[s];
		int64_t W = kind == QAWS_CURVE_KIND_HERMITE ? 3 : 6;
		sp->a = (int64_t)s << ec->param_shift;
		sp->b = (int64_t)(s + 1) << ec->param_shift;
		for (c = 0; c < dim; c++)
		{
			int64_t v[4];
			if (kind == QAWS_CURVE_KIND_HERMITE)
			{
				int64_t p0 = P[s * dim + c], p1 = P[(s + 1) * dim + c], m0 = M[s * dim + c], m1 = M[(s + 1) * dim + c];
				v[0] = 3 * p0;
				v[1] = 3 * p0 + m0;
				v[2] = 3 * p1 - m1;
				v[3] = 3 * p1;
			}
			else
			{
				/* open: span s runs P(s+1) -> P(s+2); closed: P(s) -> P(s+1), indices mod n */
				unsigned int b0 = closed ? s + n - 1 : s;
				int64_t q0 = P[((b0) % n) * dim + c], q1 = P[((b0 + 1) % n) * dim + c];
				int64_t q2 = P[((b0 + 2) % n) * dim + c], q3 = P[((b0 + 3) % n) * dim + c];
				v[0] = 6 * q1;
				v[1] = 6 * q1 + q2 - q0;
				v[2] = 6 * q2 - q3 + q1;
				v[3] = 6 * q2;
			}
			for (i = 0; i < 4; i++)
				qaws_exact_int_from_i64(&sp->h[i * D + c], v[i]);
		}
		for (i = 0; i < 4; i++)
			qaws_exact_int_from_i64(&sp->h[i * D + dim], W);
	}
	*out = ec;
	return QAWS_STATUS_OK;
}

/*
 * Polynomial C(t) = sum a_k t^k on [t0, t1], taken exactly: each double
 * coefficient is the dyadic m 2^e, so nothing is quantized but t0, t1 onto
 * the parameter lattice. With a_k = A_k / 2^F, t0 = T0 / 2^S, h = L / 2^S:
 *   C(t0 + h s) = sum_j b_j s^j,  b_j = B_j / 2^(F + S n),
 *   B_j = sum_{k >= j} A_k C(k, j) T0^(k-j) L^j 2^(S (n - k)),
 * and in Bernstein form c_i = sum_{j <= i} C(i, j) / C(n, j) b_j, cleared by
 * M = lcm_j C(n, j).
 */
static qaws_status prepare_polynomial(qaws_exact_desc const* d, qaws_curve const* curve, qaws_exact_curve** out, double* max_param)
{
	qaws_field_desc fields[8];
	qaws_scalar co[(QAWS_EXACT_MAX_DEGREE + 1) * 3];
	unsigned int nf = 0, f, ncoef = 0, got = 0, dim = (unsigned int)curve->dimension, D = dim + 1, n, i, j, k, c;
	int64_t m[(QAWS_EXACT_MAX_DEGREE + 1) * 3], T0, T1, L, M = 1;
	int e[(QAWS_EXACT_MAX_DEGREE + 1) * 3], F = 0, S, E;
	int64_t binom[QAWS_EXACT_MAX_DEGREE + 1][QAWS_EXACT_MAX_DEGREE + 1];
	double t0 = (double)curve->parameter_range.min_value, t1 = (double)curve->parameter_range.max_value, err;
	qaws_exact_int A[(QAWS_EXACT_MAX_DEGREE + 1) * 3], B[(QAWS_EXACT_MAX_DEGREE + 1) * 3], t, pw;
	qaws_exact_curve* ec;
	qaws_status st;
	if (qaws_curve_describe_fields(curve, fields, 8, &nf) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (f = 0; f < nf && f < 8; f++)
		if (fields[f].field == QAWS_FIELD_COEFFICIENTS)
			ncoef = fields[f].count;
	if (ncoef < 1 || ncoef > QAWS_EXACT_MAX_DEGREE + 1 || (dim != 2 && dim != 3))
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	n = ncoef - 1;
	if (qaws_curve_read_field(curve, QAWS_FIELD_COEFFICIENTS, co, ncoef * dim, &got) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (i = 0; i <= QAWS_EXACT_MAX_DEGREE; i++)
		for (j = 0; j <= QAWS_EXACT_MAX_DEGREE; j++)
			binom[i][j] = j > i ? 0 : (j == 0 || j == i ? 1 : binom[i - 1][j - 1] + binom[i - 1][j]);
	/* exact dyadic coefficients over a common 2^F */
	for (i = 0; i < ncoef * dim; i++)
	{
		qaws_exact_split_double((double)co[i], &m[i], &e[i]);
		if (m[i] != 0 && -e[i] > F)
			F = -e[i];
	}
	for (i = 0; i < ncoef * dim; i++)
	{
		qaws_exact_int_from_i64(&A[i], m[i]);
		if (m[i] != 0)
			TRY(qaws_exact_int_shl(&A[i], &A[i], (unsigned int)(e[i] + F)));
	}
	S = qaws_exact_param_shift_for(fabs(t0) > fabs(t1) ? fabs(t0) : fabs(t1), d->param_bits);
	TRY(qaws_exact_quantize(t0, -S, d->param_bits + 1, &T0, &err));
	*max_param = err;
	TRY(qaws_exact_quantize(t1, -S, d->param_bits + 1, &T1, &err));
	if (err > *max_param) *max_param = err;
	L = T1 - T0;
	if (L <= 0)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	/* B_j per component */
	for (j = 0; j <= n; j++)
		for (c = 0; c < dim; c++)
		{
			qaws_exact_int_zero(&B[j * dim + c]);
			for (k = j; k <= n; k++)
			{
				unsigned int r;
				TRY(qaws_exact_int_mul_i64(&t, &A[k * dim + c], binom[k][j]));
				for (r = 0; r < k - j; r++)
					TRY(qaws_exact_int_mul_i64(&t, &t, T0));
				for (r = 0; r < j; r++)
					TRY(qaws_exact_int_mul_i64(&t, &t, L));
				TRY(qaws_exact_int_shl(&t, &t, (unsigned int)(S * (int)(n - k))));
				TRY(qaws_exact_int_add(&B[j * dim + c], &B[j * dim + c], &t));
			}
		}
	for (j = 0; j <= n; j++)
	{
		int64_t a = M, b = binom[n][j];
		while (b) { int64_t r = a % b; a = b; b = r; }
		M = M / a * binom[n][j];
	}
	ec = exact_curve_alloc(d, dim, 1, n);
	if (!ec)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	ec->param_shift = S;
	ec->spans[0].a = T0;
	ec->spans[0].b = T1;
	/* x_lattice = c_i / (M 2^E), E = F + S n + space_exp2 */
	E = F + S * (int)n + d->space_exp2;
	for (i = 0; i <= n && st == QAWS_STATUS_OK; i++)
	{
		for (c = 0; c < dim; c++)
		{
			qaws_exact_int acc;
			qaws_exact_int_zero(&acc);
			for (j = 0; j <= i; j++)
			{
				st = qaws_exact_int_mul_i64(&t, &B[j * dim + c], binom[i][j] * (M / binom[n][j]));
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&acc, &acc, &t);
			}
			if (st == QAWS_STATUS_OK && E < 0)
				st = qaws_exact_int_shl(&acc, &acc, (unsigned int)-E);
			ec->spans[0].h[i * D + c] = acc;
		}
		qaws_exact_int_from_i64(&pw, M);
		if (st == QAWS_STATUS_OK && E > 0)
			st = qaws_exact_int_shl(&pw, &pw, (unsigned int)E);
		ec->spans[0].h[i * D + dim] = pw;
	}
	if (st != QAWS_STATUS_OK)
	{
		exact_curve_free(ec);
		return st;
	}
	*out = ec;
	return QAWS_STATUS_OK;
}

/*
 * Chordal / centripetal Catmull-Rom: the knot spacing needs square roots,
 * so the runtime's own preparation (per segment cubics a s^3 + b s^2 + c s
 * + d, unit spans) is frozen and taken exactly: in Bernstein form times 3,
 * (3d, 3d + c, 3d + 2c + b, 3 (a + b + c + d)), W = 3, on a common power of
 * two. QAWS_EXACT_FLAG_PREP_QUANTIZED: exact relative to that preparation.
 */
static qaws_status prepare_cr_frozen(qaws_exact_desc const* d, qaws_curve const* curve, qaws_exact_curve** out)
{
	qaws_catmull_rom_impl const* impl = (qaws_catmull_rom_impl const*)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension, D = dim + 1, n = impl->control_point_count, spans, s, c, i;
	int kpow = 0, first = 1;
	qaws_exact_curve* ec;
	qaws_status st = QAWS_STATUS_OK;
	spans = impl->closed ? n : (n >= 4 ? n - 3 : 0);
	if (spans == 0 || !impl->segment_coeffs || (dim != 2 && dim != 3))
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	/* the common scale: every coefficient m 2^e as an integer in lattice units times 2^kpow */
	for (i = 0; i < spans * dim * 4; i++)
	{
		int64_t m;
		int e;
		if (!qaws_exact_split_double((double)impl->segment_coeffs[i], &m, &e))
			return QAWS_STATUS_EXACT_UNSUPPORTED;
		if (m != 0 && (first || -(e - d->space_exp2) > kpow))
		{
			if (-(e - d->space_exp2) > kpow)
				kpow = -(e - d->space_exp2);
			first = 0;
		}
	}
	if (kpow > 1500)
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	ec = exact_curve_alloc(d, dim, spans, 3);
	if (!ec)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	ec->param_shift = qaws_exact_param_shift_for((double)spans, d->param_bits);
	for (s = 0; s < spans && st == QAWS_STATUS_OK; s++)
	{
		qaws_exact_span* sp = &ec->spans[s];
		sp->a = (int64_t)s << ec->param_shift;
		sp->b = (int64_t)(s + 1) << ec->param_shift;
		for (c = 0; c < dim && st == QAWS_STATUS_OK; c++)
		{
			qaws_exact_int k[4], t;   /* a, b, c, d */
			for (i = 0; i < 4 && st == QAWS_STATUS_OK; i++)
			{
				int64_t m;
				int e;
				qaws_exact_split_double((double)impl->segment_coeffs[(s * dim + c) * 4 + i], &m, &e);
				qaws_exact_int_from_i64(&k[i], m);
				if (m != 0)
					st = qaws_exact_int_shl(&k[i], &k[i], (unsigned int)(e - d->space_exp2 + kpow));
			}
			if (st != QAWS_STATUS_OK)
				break;
			/* H0 = 3d, H1 = 3d + c, H2 = 3d + 2c + b, H3 = 3 (a + b + c + d) */
			st = qaws_exact_int_mul_i64(&sp->h[c], &k[3], 3);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&sp->h[D + c], &sp->h[c], &k[2]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&sp->h[2 * D + c], &sp->h[D + c], &k[2]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&sp->h[2 * D + c], &sp->h[2 * D + c], &k[1]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&t, &k[0], &k[1]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&t, &t, &k[2]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&t, &t, &k[3]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul_i64(&sp->h[3 * D + c], &t, 3);
		}
		for (i = 0; i < 4 && st == QAWS_STATUS_OK; i++)
		{
			qaws_exact_int_from_i64(&sp->h[i * D + dim], 3);
			st = qaws_exact_int_shl(&sp->h[i * D + dim], &sp->h[i * D + dim], (unsigned int)kpow);
		}
	}
	if (st != QAWS_STATUS_OK)
	{
		exact_curve_free(ec);
		return st;
	}
	*out = ec;
	return QAWS_STATUS_OK;
}

/*
 * Composite: segment i occupies [i, i + 1]. Each segment is prepared
 * exactly; a span's integer Bezier does not depend on its parameter
 * interval, so the spans are kept as they are and only their bounds move:
 * x on the segment's lattice maps to T = i 2^S + (x - a0) 2^S / (bN - a0).
 * Bounds that are not on the composite lattice are rounded (the geometry
 * stays exact; the parameterization error is reported).
 */
static qaws_status prepare_composite(qaws_exact_desc const* d, qaws_curve const* curve, qaws_exact_curve** out, qaws_exact_report* rep)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	qaws_exact_curve** seg = NULL;
	qaws_exact_curve* ec = NULL;
	unsigned int i, k, total = 0, idx = 0;
	int S;
	double perr = 0;
	qaws_status st = QAWS_STATUS_OK;
	if (!impl || impl->segment_count == 0)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	seg = (qaws_exact_curve**)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_curve*) * impl->segment_count));
	if (!seg)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(seg, 0, sizeof(qaws_exact_curve*) * impl->segment_count);
	if (rep)
		memset(rep, 0, sizeof(*rep));
	for (i = 0; i < impl->segment_count && st == QAWS_STATUS_OK; i++)
	{
		qaws_exact_report r;
		st = qaws_exact_curve_prepare(d, impl->segments[i], &seg[i], &r);
		if (st == QAWS_STATUS_OK && seg[i]->dimension != curve->dimension)
			st = QAWS_STATUS_INVALID_ARGUMENT;
		if (st == QAWS_STATUS_OK)
		{
			total += seg[i]->span_count;
			if (rep)
			{
				rep->flags |= r.flags;
				if (r.storage_bits > rep->storage_bits) rep->storage_bits = r.storage_bits;
				if (r.max_position_quantization_error > rep->max_position_quantization_error)
					rep->max_position_quantization_error = r.max_position_quantization_error;
				if (r.max_weight_quantization_error > rep->max_weight_quantization_error)
					rep->max_weight_quantization_error = r.max_weight_quantization_error;
				if (r.parameter_quantization_error > rep->parameter_quantization_error)
					rep->parameter_quantization_error = r.parameter_quantization_error;
			}
		}
	}
	if (st == QAWS_STATUS_OK)
	{
		ec = (qaws_exact_curve*)qaws_internal_alloc(NULL, (unsigned long)sizeof(qaws_exact_curve));
		if (ec)
		{
			memset(ec, 0, sizeof(*ec));
			ec->spans = (qaws_exact_span*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_span) * total));
		}
		if (!ec || !ec->spans)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (st == QAWS_STATUS_OK)
	{
		ec->dimension = curve->dimension;
		ec->space_exp2 = d->space_exp2;
		ec->desc = *d;
		S = qaws_exact_param_shift_for((double)impl->segment_count, d->param_bits);
		ec->param_shift = S;
		for (i = 0; i < impl->segment_count && st == QAWS_STATUS_OK; i++)
		{
			qaws_exact_curve* s = seg[i];
			int64_t a0 = s->spans[0].a, len = s->spans[s->span_count - 1].b - a0;
			for (k = 0; k < s->span_count && st == QAWS_STATUS_OK; k++)
			{
				qaws_exact_span* sp = &ec->spans[idx];
				int64_t bound[2];
				unsigned int e;
				for (e = 0; e < 2; e++)
				{
					/* (x - a0) 2^S / len, rounded to nearest; exact when it divides */
					int64_t x = (e == 0 ? s->spans[k].a : s->spans[k].b) - a0;
					double q = ldexp((double)x, S) / (double)len, r = nearbyint(q);
					qaws_exact_int num, den, quo, rem;
					qaws_exact_int_from_i64(&num, x);
					qaws_exact_int_shl(&num, &num, (unsigned int)S);
					qaws_exact_int_from_i64(&den, len);
					qaws_exact_int_divmod(&quo, &rem, &num, &den);
					if (!qaws_exact_int_is_zero(&rem))
					{
						double err = fabs(q - r) * ldexp(1.0, -S);
						if (err > perr) perr = err;
					}
					else
						r = qaws_exact_int_to_double(&quo);
					bound[e] = ((int64_t)i << S) + (int64_t)r;
				}
				if (bound[1] <= bound[0] || (idx > 0 && bound[0] != ec->spans[idx - 1].b))
					st = QAWS_STATUS_EXACT_RANGE_EXCEEDED;   /* a span collapsed on the composite lattice */
				sp->degree = s->spans[k].degree;
				sp->a = bound[0];
				sp->b = bound[1];
				sp->h = s->spans[k].h;   /* moved */
				s->spans[k].h = NULL;
				idx++;
				ec->span_count = idx;
			}
		}
	}
	for (i = 0; i < impl->segment_count; i++)
		exact_curve_free(seg[i]);
	qaws_internal_dealloc(NULL, seg);
	if (st != QAWS_STATUS_OK)
	{
		exact_curve_free(ec);
		return st;
	}
	if (rep)
	{
		rep->quality = QAWS_NUMERIC_EXACT_RATIONAL;
		if (perr > 0)
		{
			rep->flags |= QAWS_EXACT_FLAG_INPUT_QUANTIZED;
			if (perr > rep->parameter_quantization_error)
				rep->parameter_quantization_error = perr;
		}
	}
	*out = ec;
	return QAWS_STATUS_OK;
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
	qaws_exact_hfrac* local = NULL;
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
	if (kind == QAWS_CURVE_KIND_COMPOSITE)
		return prepare_composite(&d, curve, out_curve, out_report);
	if (kind == QAWS_CURVE_KIND_CATMULL_ROM && ((qaws_catmull_rom_impl const*)curve->impl)->parameterization != QAWS_PARAMETERIZATION_UNIFORM)
	{
		qaws_exact_curve* ec2 = NULL;
		qaws_status st2 = prepare_cr_frozen(&d, curve, &ec2);
		unsigned int s2, i2;
		if (st2 != QAWS_STATUS_OK)
			return st2;
		if (out_report)
		{
			memset(out_report, 0, sizeof(*out_report));
			out_report->quality = QAWS_NUMERIC_EXACT_RATIONAL;
			out_report->flags = QAWS_EXACT_FLAG_PREP_QUANTIZED;
			for (s2 = 0; s2 < ec2->span_count; s2++)
				for (i2 = 0; i2 < 4 * (unsigned int)(ec2->dimension + 1); i2++)
					if (qaws_exact_int_bits(&ec2->spans[s2].h[i2]) > out_report->storage_bits)
						out_report->storage_bits = qaws_exact_int_bits(&ec2->spans[s2].h[i2]);
		}
		*out_curve = ec2;
		return QAWS_STATUS_OK;
	}
	if (kind == QAWS_CURVE_KIND_HERMITE || kind == QAWS_CURVE_KIND_CATMULL_ROM || kind == QAWS_CURVE_KIND_POLYNOMIAL)
	{
		double qpos = 0, qpar = 0;
		qaws_exact_curve* ec2 = NULL;
		qaws_status st2 = kind == QAWS_CURVE_KIND_POLYNOMIAL ? prepare_polynomial(&d, curve, &ec2, &qpar)
		                                                     : prepare_cubic_family(&d, curve, kind, &ec2, &qpos);
		if (st2 != QAWS_STATUS_OK)
			return st2;
		if (out_report)
		{
			unsigned int s2, i2;
			memset(out_report, 0, sizeof(*out_report));
			out_report->quality = QAWS_NUMERIC_EXACT_RATIONAL;
			out_report->flags = (qpos > 0 || qpar > 0) ? QAWS_EXACT_FLAG_INPUT_QUANTIZED : QAWS_EXACT_FLAG_NONE;
			out_report->max_position_quantization_error = qpos;
			out_report->parameter_quantization_error = qpar;
			for (s2 = 0; s2 < ec2->span_count; s2++)
				for (i2 = 0; i2 < (ec2->spans[s2].degree + 1) * (unsigned int)(ec2->dimension + 1); i2++)
					if (qaws_exact_int_bits(&ec2->spans[s2].h[i2]) > out_report->storage_bits)
						out_report->storage_bits = qaws_exact_int_bits(&ec2->spans[s2].h[i2]);
		}
		*out_curve = ec2;
		return QAWS_STATUS_OK;
	}
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
				if (qaws_exact_quantize((double)kn[i], -ec->param_shift, d.param_bits + 1, &K[i], &err) != QAWS_STATUS_OK)
				{
					exact_curve_free(ec);
					return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
				}
				if (err > max_k) max_k = err;
			}
	}

	/* homogeneous integer control points */
	local = (qaws_exact_hfrac*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_hfrac) * n));
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
			st = qaws_exact_quantize((double)cps[i * dim + c], d.space_exp2, d.coord_bits, &X, &err);
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
				qaws_exact_hfrac b[QAWS_EXACT_MAX_DEGREE + 1];
				unsigned int j;
				sp->a = K[s];
				sp->b = K[s + 1];
				for (j = 0; j <= p && st == QAWS_STATUS_OK; j++)
					st = qaws_exact_blossom(&local[s - p], K, s, p, j, D, &b[j]);
				if (st == QAWS_STATUS_OK)
					st = qaws_exact_clear_denominators(b, p, D, sp->h);
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
qaws_status qaws_exact_homogeneous_derivative(qaws_exact_span const* sp, unsigned int D, int64_t x, unsigned int j, qaws_exact_int* out)
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
		TRY(qaws_exact_homogeneous_derivative(sp, D, T - sp->a, j, H[j]));
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

unsigned int qaws_exact_curve_span_count(qaws_exact_curve const* curve)
{
	return curve ? curve->span_count : 0;
}

qaws_status qaws_exact_curve_span_bezier(qaws_exact_curve const* curve, unsigned int s, unsigned int* out_degree, double* out_t0,
	double* out_t1, double* out_points, double* out_weights)
{
	qaws_exact_span const* sp;
	unsigned int D, i, c, dim;
	double scale;
	qaws_exact_int wmax;
	if (!curve || s >= curve->span_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	sp = &curve->spans[s];
	dim = (unsigned int)curve->dimension;
	D = dim + 1;
	scale = ldexp(1.0, curve->space_exp2);
	if (out_degree) *out_degree = sp->degree;
	if (out_t0) *out_t0 = ldexp((double)sp->a, -curve->param_shift);
	if (out_t1) *out_t1 = ldexp((double)sp->b, -curve->param_shift);
	wmax = sp->h[dim];
	for (i = 1; i <= sp->degree; i++)
		if (qaws_exact_int_cmp(&sp->h[i * D + dim], &wmax) > 0)
			wmax = sp->h[i * D + dim];
	for (i = 0; i <= sp->degree; i++)
	{
		for (c = 0; c < dim && out_points; c++)
			out_points[i * dim + c] = qaws_exact_ratio_to_double(&sp->h[i * D + c], &sp->h[i * D + dim]) * scale;
		if (out_weights)
			out_weights[i] = qaws_exact_ratio_to_double(&sp->h[i * D + dim], &wmax);
	}
	return QAWS_STATUS_OK;
}
