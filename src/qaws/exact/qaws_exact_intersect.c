#include "qaws_exact_curve.h"
#include "qaws_exact_roots.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

/* ------------------------------------------------------------------ */
/*  Dyadic rationals m 2^e                                             */
/* ------------------------------------------------------------------ */

typedef struct dyadic
{
	qaws_exact_int m;
	int e;
} dyadic;

static qaws_status dy_from_double(dyadic* r, double x)
{
	int64_t m;
	if (!(x - x == 0.0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_exact_split_double(x, &m, &r->e);
	qaws_exact_int_from_i64(&r->m, m);
	if (m == 0)
		r->e = 0;
	return QAWS_STATUS_OK;
}

/* r = a + sign b (exponents aligned to the smaller one) */
static qaws_status dy_add(dyadic* r, dyadic const* a, dyadic const* b, int sign)
{
	qaws_exact_int x = a->m, y = b->m;
	int e = a->e < b->e ? a->e : b->e;
	qaws_status st;
	if (qaws_exact_int_is_zero(&a->m))
		e = b->e;
	else if (qaws_exact_int_is_zero(&b->m))
		e = a->e;
	if (!qaws_exact_int_is_zero(&x) && a->e > e)
		TRY(qaws_exact_int_shl(&x, &x, (unsigned int)(a->e - e)));
	if (!qaws_exact_int_is_zero(&y) && b->e > e)
		TRY(qaws_exact_int_shl(&y, &y, (unsigned int)(b->e - e)));
	r->e = e;
	return sign > 0 ? qaws_exact_int_add(&r->m, &x, &y) : qaws_exact_int_sub(&r->m, &x, &y);
}

static qaws_status dy_mul_int(dyadic* r, dyadic const* a, qaws_exact_int const* k, int e_extra)
{
	r->e = a->e + e_extra;
	return qaws_exact_int_mul(&r->m, &a->m, k);
}

/* ------------------------------------------------------------------ */
/*  Hits                                                               */
/* ------------------------------------------------------------------ */

/* Trailing zero bits of a non-zero magnitude. */
static unsigned int ctz_int(qaws_exact_int const* x)
{
	unsigned int i, b;
	for (i = 0; i < (unsigned int)x->size; i++)
		if (x->limb[i])
		{
			for (b = 0; !((x->limb[i] >> b) & 1u); b++)
				;
			return 32 * i + b;
		}
	return 0;
}

/* t = (a 2^d + (b - a) index) 2^-(d + shift): nearest double, and whether it is exact. */
qaws_status qaws_exact_span_param_to_double(qaws_exact_span const* sp, int shift, uint64_t index, int depth, double* out, int* exact)
{
	qaws_exact_int num, t;
	qaws_status st;
	unsigned int bits;
	qaws_exact_int_from_i64(&num, sp->a);
	TRY(qaws_exact_int_shl(&num, &num, (unsigned int)depth));
	/* index < 2^62: split to stay in int64 */
	qaws_exact_int_from_i64(&t, (int64_t)(index >> 31));
	TRY(qaws_exact_int_mul_i64(&t, &t, sp->b - sp->a));
	TRY(qaws_exact_int_shl(&t, &t, 31));
	TRY(qaws_exact_int_add(&num, &num, &t));
	qaws_exact_int_from_i64(&t, (int64_t)(index & 0x7FFFFFFFu));
	TRY(qaws_exact_int_mul_i64(&t, &t, sp->b - sp->a));
	TRY(qaws_exact_int_add(&num, &num, &t));
	*out = ldexp(qaws_exact_int_to_double(&num), -(depth + shift));
	bits = qaws_exact_int_bits(&num);
	*exact = qaws_exact_int_is_zero(&num) || (bits - ctz_int(&num) <= 53 && fabs(*out) >= 2.2250738585072014e-308);
	return QAWS_STATUS_OK;
}

static qaws_status emit(qaws_exact_hit* out, unsigned int capacity, unsigned int* count, qaws_exact_hit_kind kind, double lo, double hi)
{
	if (*count >= capacity)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	out[*count].kind = kind;
	out[*count].t_lo = lo;
	out[*count].t_hi = hi;
	(*count)++;
	return QAWS_STATUS_OK;
}

/*
 * Hits of the hyperplane N . (x - P) = 0 (N, P exact dyadics in world
 * units). On a span, W (N . (C - P)) has the integer Bernstein coefficients
 *   g_i = sum_c N_c (H_ic 2^space_exp2 - P_c W_i)
 * (W > 0), brought to one exponent.
 */
static qaws_status hyperplane_hits(qaws_exact_curve const* curve, dyadic const* N, dyadic const* P, double width, qaws_exact_hit* out,
	unsigned int capacity, unsigned int* out_count)
{
	unsigned int dim = (unsigned int)curve->dimension, D = dim + 1, s, i, c, count = 0;
	qaws_status st = QAWS_STATUS_OK;
	*out_count = 0;
	for (s = 0; s < curve->span_count && st == QAWS_STATUS_OK; s++)
	{
		qaws_exact_span const* sp = &curve->spans[s];
		unsigned int n = sp->degree;
		dyadic g[QAWS_EXACT_MAX_DEGREE + 1];
		qaws_exact_int b[QAWS_EXACT_MAX_DEGREE + 1];
		qaws_exact_root roots[QAWS_EXACT_MAX_DEGREE + 2];
		unsigned int nr = 0, r;
		int emin = 0, any = 0, zero = 1;
		double t_a, t_b;
		int ex;
		for (i = 0; i <= n; i++)
		{
			qaws_exact_int_zero(&g[i].m);
			g[i].e = 0;
			for (c = 0; c < dim; c++)
			{
				dyadic t1, t2, NP;
				qaws_exact_int pm;
				TRY(dy_mul_int(&t1, &N[c], &sp->h[i * D + c], curve->space_exp2));
				pm = P[c].m;
				NP.e = N[c].e + P[c].e;
				TRY(qaws_exact_int_mul(&NP.m, &N[c].m, &pm));
				TRY(dy_mul_int(&t2, &NP, &sp->h[i * D + dim], 0));
				TRY(dy_add(&t1, &t1, &t2, -1));
				TRY(dy_add(&g[i], &g[i], &t1, 1));
			}
			if (!qaws_exact_int_is_zero(&g[i].m))
			{
				if (!any || g[i].e < emin)
					emin = g[i].e;
				any = 1;
				zero = 0;
			}
		}
		TRY(qaws_exact_span_param_to_double(sp, curve->param_shift, 0, 0, &t_a, &ex));
		TRY(qaws_exact_span_param_to_double(sp, curve->param_shift, 1, 0, &t_b, &ex));
		if (zero)
		{
			/* the whole span lies on the hyperplane: merge with an overlap or a point that ends at t_a */
			if (count > 0 && out[count - 1].t_hi == t_a && out[count - 1].kind != QAWS_EXACT_HIT_CROSSING)
			{
				out[count - 1].kind = QAWS_EXACT_HIT_OVERLAP;
				out[count - 1].t_hi = t_b;
			}
			else
				TRY(emit(out, capacity, &count, QAWS_EXACT_HIT_OVERLAP, t_a, t_b));
			*out_count = count;
			continue;
		}
		for (i = 0; i <= n; i++)
		{
			b[i] = g[i].m;
			if (!qaws_exact_int_is_zero(&b[i]) && g[i].e > emin)
				TRY(qaws_exact_int_shl(&b[i], &b[i], (unsigned int)(g[i].e - emin)));
		}
		st = qaws_exact_bernstein_isolate(b, n, roots, QAWS_EXACT_MAX_DEGREE + 2, &nr);
		if (st != QAWS_STATUS_OK)
			break;
		for (r = 0; r < nr && st == QAWS_STATUS_OK; r++)
		{
			qaws_exact_root rt = roots[r];
			double lo, hi;
			int ex_lo, ex_hi;
			if (rt.exact && rt.index == 0)
			{
				/* s = 0: already reported as the end of the previous span (or inside its overlap) */
				if (count > 0 && out[count - 1].t_hi == t_a && out[count - 1].kind != QAWS_EXACT_HIT_CROSSING)
					continue;
			}
			if (!rt.exact)
			{
				int depth = 60;
				if (width > 0)
				{
					double len = t_b - t_a;
					depth = (int)ceil(log2(len / width));
					if (depth < rt.depth) depth = rt.depth;
					if (depth > 60) depth = 60;
				}
				TRY(qaws_exact_bernstein_refine(b, n, &rt, depth));
			}
			TRY(qaws_exact_span_param_to_double(sp, curve->param_shift, rt.index, rt.depth, &lo, &ex_lo));
			if (rt.exact)
			{
				hi = lo;
				if (!ex_lo)
				{
					lo = nextafter(lo, -HUGE_VAL);
					hi = nextafter(hi, HUGE_VAL);
				}
				st = emit(out, capacity, &count, QAWS_EXACT_HIT_POINT, lo, hi);
			}
			else
			{
				TRY(qaws_exact_span_param_to_double(sp, curve->param_shift, rt.index + 1, rt.depth, &hi, &ex_hi));
				if (!ex_lo) lo = nextafter(lo, -HUGE_VAL);
				if (!ex_hi) hi = nextafter(hi, HUGE_VAL);
				st = emit(out, capacity, &count, QAWS_EXACT_HIT_CROSSING, lo, hi);
			}
		}
		*out_count = count;
	}
	*out_count = count;
	return st;
}

qaws_status qaws_exact_curve_line_hits(qaws_exact_curve const* curve, double const p0[2], double const p1[2], double width,
	qaws_exact_hit* out_hits, unsigned int capacity, unsigned int* out_count)
{
	dyadic N[2], P[2], a, b;
	qaws_status st;
	if (!curve || !p0 || !p1 || !out_count || (!out_hits && capacity) || curve->dimension != 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	/* N = (-(p1y - p0y), p1x - p0x), exactly */
	TRY(dy_from_double(&P[0], p0[0]));
	TRY(dy_from_double(&P[1], p0[1]));
	TRY(dy_from_double(&a, p1[1]));
	TRY(dy_add(&b, &P[1], &a, -1));   /* p0y - p1y */
	N[0] = b;
	TRY(dy_from_double(&a, p1[0]));
	TRY(dy_add(&N[1], &a, &P[0], -1));
	if (qaws_exact_int_is_zero(&N[0].m) && qaws_exact_int_is_zero(&N[1].m))
		return QAWS_STATUS_INVALID_ARGUMENT;
	return hyperplane_hits(curve, N, P, width, out_hits, capacity, out_count);
}

qaws_status qaws_exact_curve_plane_hits(qaws_exact_curve const* curve, double const point[3], double const normal[3], double width,
	qaws_exact_hit* out_hits, unsigned int capacity, unsigned int* out_count)
{
	dyadic N[3], P[3];
	unsigned int c;
	qaws_status st;
	if (!curve || !point || !normal || !out_count || (!out_hits && capacity) || curve->dimension != 3)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	for (c = 0; c < 3; c++)
	{
		TRY(dy_from_double(&P[c], point[c]));
		TRY(dy_from_double(&N[c], normal[c]));
	}
	if (qaws_exact_int_is_zero(&N[0].m) && qaws_exact_int_is_zero(&N[1].m) && qaws_exact_int_is_zero(&N[2].m))
		return QAWS_STATUS_INVALID_ARGUMENT;
	return hyperplane_hits(curve, N, P, width, out_hits, capacity, out_count);
}
