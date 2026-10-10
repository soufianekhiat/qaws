#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define CC_COEF QAWS_EXACT_POLY_COEF
#define CC_MAX_ROOTS QAWS_EXACT_SOLVE_MAX_ROOTS

typedef qaws_exact_poly poly;
static void* sv_alloc(size_t bytes)
{
	return qaws_internal_alloc(NULL, (unsigned long)bytes);
}

static void sv_free(void* p)
{
	qaws_internal_dealloc(NULL, p);
}

/* ------------------------------------------------------------------ */
/*  Bezout matrices                                                    */
/* ------------------------------------------------------------------ */

/*
 * Bezout matrix of P, Q (degree n, integer coefficients):
 * (P(t) Q(s) - P(s) Q(t)) / (t - s) = sum B_ij t^i s^j. A term
 * c_ab (t^a s^b - t^b s^a), a > b, c_ab = p_a q_b - p_b q_a, divides into
 * t^b s^b sum_{k < a - b} t^(a-b-1-k) s^k. Symmetric.
 */
qaws_status qaws_exact_bezout(poly const* P, poly const* Q, unsigned int n, qaws_exact_int* B)
{
	qaws_exact_int c, t;
	unsigned int a, b, k;
	qaws_status st;
	for (a = 0; a < n * n; a++)
		qaws_exact_int_zero(&B[a]);
	for (a = 1; a <= n; a++)
		for (b = 0; b < a; b++)
		{
			TRY(qaws_exact_int_mul(&c, &P->c[a], &Q->c[b]));
			TRY(qaws_exact_int_mul(&t, &P->c[b], &Q->c[a]));
			TRY(qaws_exact_int_sub(&c, &c, &t));
			if (qaws_exact_int_is_zero(&c))
				continue;
			for (k = 0; k + b < a; k++)
				TRY(qaws_exact_int_add(&B[(a - 1 - k) * n + b + k], &B[(a - 1 - k) * n + b + k], &c));
		}
	return QAWS_STATUS_OK;
}

/* The same with polynomial coefficients: P, Q of degree k in s, P[i] the coefficient of s^i. */
qaws_status qaws_exact_bezout_poly(poly const* P, poly const* Q, unsigned int k, poly* E)
{
	poly* t;
	unsigned int a, b, j;
	qaws_status st = QAWS_STATUS_OK;
	t = (poly*)sv_alloc(sizeof(poly) * 2);
	if (!t)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (a = 0; a < k * k; a++)
		qaws_exact_poly_zero(&E[a], 0);
	for (a = 1; a <= k && st == QAWS_STATUS_OK; a++)
		for (b = 0; b < a && st == QAWS_STATUS_OK; b++)
		{
			st = qaws_exact_poly_mul(&t[0], &P[a], &Q[b]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_mul(&t[1], &P[b], &Q[a]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(&t[0], &t[1], -1);
			for (j = 0; j + b < a && st == QAWS_STATUS_OK; j++)
				st = qaws_exact_poly_acc(&E[(a - 1 - j) * k + b + j], &t[0], 1);
		}
	sv_free(t);
	return st;
}

/* det, and the row-0 cofactors C00 and C01 (C01 = -minor_01). */
qaws_status qaws_exact_det_cofactors(poly const* M, unsigned int n, poly* det, poly* C00, poly* C01)
{
	poly* minor;
	unsigned int j, k;
	qaws_status st = QAWS_STATUS_OK;
	minor = (poly*)sv_alloc(sizeof(poly) * 2);
	if (!minor)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	qaws_exact_poly_zero(det, 0);
	for (j = 0; j < n && st == QAWS_STATUS_OK; j++)
	{
		st = qaws_exact_poly_det(M, n, 1, 1u << j, &minor[0]);
		if (st == QAWS_STATUS_OK && j == 0)
			*C00 = minor[0];
		if (st == QAWS_STATUS_OK && j == 1)
		{
			*C01 = minor[0];
			for (k = 0; k <= C01->deg; k++)
				qaws_exact_int_neg(&C01->c[k], &C01->c[k]);
		}
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_mul(&minor[1], &M[j], &minor[0]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(det, &minor[1], (j & 1) ? -1 : 1);
	}
	sv_free(minor);
	return st;
}

/* out = sum_i H[i] p^i q^(k - i) (the substitution s = p / q, times q^k) */
qaws_status qaws_exact_substitute(poly const* H, unsigned int k, poly const* p, poly const* q, poly* out)
{
	poly* pw;   /* p^0..p^k, q^0..q^k, scratch */
	unsigned int i;
	qaws_status st = QAWS_STATUS_OK;
	pw = (poly*)sv_alloc(sizeof(poly) * (2 * (k + 1) + 1));
	if (!pw)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	qaws_exact_poly_const(&pw[0], 1);
	qaws_exact_poly_const(&pw[k + 1], 1);
	for (i = 1; i <= k && st == QAWS_STATUS_OK; i++)
	{
		st = qaws_exact_poly_mul(&pw[i], &pw[i - 1], p);
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_mul(&pw[k + 1 + i], &pw[k + i], q);
	}
	qaws_exact_poly_zero(out, 0);
	for (i = 0; i <= k && st == QAWS_STATUS_OK; i++)
	{
		poly* t = &pw[2 * (k + 1)];
		st = qaws_exact_poly_mul(t, &pw[i], &pw[k + 1 + (k - i)]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_mul(t, t, &H[i]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(out, t, 1);
	}
	sv_free(pw);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Root processing                                                    */
/* ------------------------------------------------------------------ */

/* Sign on the root interval (+1 / -1 when every coefficient there has it, else 0); at an exact root, the exact sign. */
static qaws_status sign_on(qaws_exact_int const* b, unsigned int n, qaws_exact_root const* r, int* out)
{
	qaws_exact_int c[CC_COEF];
	unsigned int i;
	int s;
	qaws_status st;
	if (r->exact)
	{
		int at_end = r->index == ((uint64_t)1 << r->depth);
		TRY(qaws_exact_bernstein_restrict(b, n, r->index - (uint64_t)at_end, r->depth, c));
		*out = qaws_exact_int_sign(&c[at_end ? n : 0]);
		return QAWS_STATUS_OK;
	}
	TRY(qaws_exact_bernstein_restrict(b, n, r->index, r->depth, c));
	s = qaws_exact_int_sign(&c[0]);
	for (i = 1; i <= n && s != 0; i++)
		if (qaws_exact_int_sign(&c[i]) != s)
			s = 0;
	*out = s;
	return QAWS_STATUS_OK;
}

/* Sign variations on the (open) root interval. */
static qaws_status variations_on(qaws_exact_int const* b, unsigned int n, qaws_exact_root const* r, unsigned int* out)
{
	qaws_exact_int c[CC_COEF];
	unsigned int i, v = 0;
	int last = 0;
	qaws_status st;
	TRY(qaws_exact_bernstein_restrict(b, n, r->index, r->depth, c));
	for (i = 0; i <= n; i++)
	{
		int s = qaws_exact_int_sign(&c[i]);
		if (s == 0)
			continue;
		if (last != 0 && s != last)
			v++;
		last = s;
	}
	*out = v;
	return QAWS_STATUS_OK;
}

/* Enclosure of num / den (den != 0) by doubles; exact when dyadic and representable. */
void qaws_exact_ratio_enclose(qaws_exact_int const* num, qaws_exact_int const* den, double* lo, double* hi, int* exact)
{
	qaws_exact_int g, n2, d2, p2;
	double v;
	unsigned int db;
	if (qaws_exact_int_is_zero(num))
	{
		*lo = *hi = 0;
		*exact = 1;
		return;
	}
	qaws_exact_int_gcd(&g, num, den);
	qaws_exact_int_divmod(&n2, NULL, num, &g);
	qaws_exact_int_divmod(&d2, NULL, den, &g);
	if (qaws_exact_int_sign(&d2) < 0)
	{
		qaws_exact_int_neg(&n2, &n2);
		qaws_exact_int_neg(&d2, &d2);
	}
	v = qaws_exact_ratio_to_double(&n2, &d2);
	db = qaws_exact_int_bits(&d2);
	qaws_exact_int_from_i64(&p2, 1);
	qaws_exact_int_shl(&p2, &p2, db - 1);
	*exact = qaws_exact_int_cmp(&p2, &d2) == 0 && qaws_exact_int_bits(&n2) <= 53 && fabs(v) >= 2.2250738585072014e-308;
	*lo = *exact ? v : nextafter(v, -HUGE_VAL);
	*hi = *exact ? v : nextafter(v, HUGE_VAL);
}

typedef struct sys_work
{
	poly R, C0, C1, D, G[2], Z[2];
	qaws_exact_int bR[CC_COEF], b0[CC_COEF], b1[CC_COEF], bD[CC_COEF], bH[2][CC_COEF], bG[2][CC_COEF], bZ[2][CC_COEF];
	unsigned int dH[2], dG[2];
	int hzero[2], gconst[2], zstate[2];
} sys_work;

/*
 * Is C01 (which 0: s = 0) or C00 - C01 (which 1: s = 1) exactly zero at the
 * irrational root isolated by rt? Z = gcd(R, .) divides R, which has one root
 * there: an odd number of sign variations of Z on the interval means the
 * root is Z's (*out 1), an even one that it is not (*out 0, the refinement
 * then decides the sign); -1 while undecided.
 */
static qaws_status zero_at(sys_work* w, int which, qaws_exact_root const* rt, int* out)
{
	unsigned int v = 0;
	qaws_status st = QAWS_STATUS_OK;
	*out = -1;
	if (w->zstate[which] < 0)
	{
		poly h = which ? w->D : w->C1;
		qaws_exact_poly_trim(&h);
		st = qaws_exact_poly_gcd(&w->R, &h, &w->Z[which]);
		if (st != QAWS_STATUS_OK)
			return st;
		w->zstate[which] = w->Z[which].deg > 0;
		if (w->zstate[which])
			TRY(qaws_exact_poly_to_bernstein(&w->Z[which], w->bZ[which]));
	}
	if (!w->zstate[which])
	{
		*out = 0;
		return QAWS_STATUS_OK;
	}
	TRY(variations_on(w->bZ[which], w->Z[which].deg, rt, &v));
	if (v == 1)
		*out = 1;
	else if (v == 0)
		*out = 0;
	return QAWS_STATUS_OK;
}

/*
 * Roots r in [0, 1] of R with s = C01 / C00 proven in [0, 1] and every
 * constraint proven, refined to the narrowest interval. *common is set
 * (with QAWS_STATUS_CERTIFICATION_FAILED) when R vanishes identically.
 */
qaws_status qaws_exact_solve_system(poly const* R0, poly const* C00, poly const* C01, qaws_exact_constraint const* cons, unsigned int ncons, qaws_exact_local_hit* out,
	unsigned int capacity, unsigned int* count, int* common)
{
	sys_work* w;
	qaws_exact_root roots[CC_MAX_ROOTS];
	unsigned int N, N1, nr = 0, k, c;
	qaws_status st = QAWS_STATUS_OK;
	*count = 0;
	*common = 0;
	w = (sys_work*)sv_alloc(sizeof(sys_work));
	if (!w)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	w->R = *R0;
	w->zstate[0] = w->zstate[1] = -1;
	qaws_exact_poly_trim(&w->R);
	if (qaws_exact_poly_is_zero(&w->R))
	{
		*common = 1;
		sv_free(w);
		return QAWS_STATUS_CERTIFICATION_FAILED;
	}
	N = w->R.deg;
	if (N == 0 || N > QAWS_EXACT_ROOTS_MAX_DEGREE)
	{
		sv_free(w);
		return N == 0 ? QAWS_STATUS_OK : QAWS_STATUS_EXACT_UNSUPPORTED;
	}
	w->C0 = *C00;
	w->C1 = *C01;
	N1 = w->C0.deg > w->C1.deg ? w->C0.deg : w->C1.deg;
	st = qaws_exact_poly_pad(&w->C0, N1);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_pad(&w->C1, N1);
	if (st == QAWS_STATUS_OK)
	{
		w->D = w->C0;
		st = qaws_exact_poly_acc(&w->D, &w->C1, -1);
	}
	if (st == QAWS_STATUS_OK && N1 > QAWS_EXACT_ROOTS_MAX_DEGREE)
		st = QAWS_STATUS_EXACT_UNSUPPORTED;
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_to_bernstein(&w->R, w->bR);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_to_bernstein(&w->C0, w->b0);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_to_bernstein(&w->C1, w->b1);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_to_bernstein(&w->D, w->bD);
	for (c = 0; c < ncons && st == QAWS_STATUS_OK; c++)
	{
		poly h = *cons[c].H;
		qaws_exact_poly_trim(&h);
		w->hzero[c] = qaws_exact_poly_is_zero(&h);
		w->gconst[c] = 0;
		w->dH[c] = h.deg;
		if (w->hzero[c])
			continue;
		if (h.deg > QAWS_EXACT_ROOTS_MAX_DEGREE)
			st = QAWS_STATUS_EXACT_UNSUPPORTED;
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_to_bernstein(&h, w->bH[c]);
		if (st == QAWS_STATUS_OK && cons[c].kind == QAWS_EXACT_CON_ZERO)
		{
			st = qaws_exact_poly_gcd(&w->R, &h, &w->G[c]);
			if (st == QAWS_STATUS_OK)
			{
				w->dG[c] = w->G[c].deg;
				w->gconst[c] = w->G[c].deg == 0;
				if (!w->gconst[c])
					st = qaws_exact_poly_to_bernstein(&w->G[c], w->bG[c]);
			}
		}
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_exact_bernstein_isolate(w->bR, N, roots, CC_MAX_ROOTS, &nr);

	for (k = 0; k < nr && st == QAWS_STATUS_OK; k++)
	{
		qaws_exact_root rt = roots[k];
		int s00 = 0, s01 = 0, sd = 0, inside = -1, cst[2] = { -1, -1 }, accept = -1, s_end = -1;
		while (st == QAWS_STATUS_OK && accept < 0)
		{
			int undecided = 0;
			if (inside < 0)
			{
				st = sign_on(w->b0, N1, &rt, &s00);
				if (st == QAWS_STATUS_OK) st = sign_on(w->b1, N1, &rt, &s01);
				if (st == QAWS_STATUS_OK) st = sign_on(w->bD, N1, &rt, &sd);
				if (st != QAWS_STATUS_OK)
					break;
				if (s00 != 0 && ((s01 != 0 && s01 != s00) || (sd != 0 && sd != s00)))
					inside = 0;
				else if (s00 != 0 && (s01 == s00 || (rt.exact && s01 == 0)) && (sd == s00 || (rt.exact && sd == 0)))
					inside = 1;
				else if (s00 != 0 && !rt.exact && (s01 == 0) != (sd == 0))
				{
					/* s on a span end at an irrational root: an exact zero test */
					int z = -1;
					st = zero_at(w, s01 == 0 ? 0 : 1, &rt, &z);
					if (st != QAWS_STATUS_OK)
						break;
					if (z == 1 && (s01 == 0 ? sd == s00 : s01 == s00))
					{
						inside = 1;
						s_end = s01 == 0 ? 0 : 1;
					}
				}
			}
			if (inside == 0)
			{
				accept = 0;
				break;
			}
			for (c = 0; c < ncons && st == QAWS_STATUS_OK; c++)
			{
				if (cst[c] >= 0)
					continue;
				if (cons[c].kind == QAWS_EXACT_CON_SIGN_OF_C00)
				{
					int sh = 0;
					if (w->hzero[c])
						cst[c] = 0;
					else if (s00 != 0)
					{
						st = sign_on(w->bH[c], w->dH[c], &rt, &sh);
						if (rt.exact)
							cst[c] = sh == s00;
						else if (sh != 0)
							cst[c] = sh == s00;
					}
				}
				else
				{
					if (w->hzero[c])
						cst[c] = 1;
					else if (w->gconst[c])
						cst[c] = 0;
					else if (rt.exact)
					{
						int sg = 0;
						st = sign_on(w->bG[c], w->dG[c], &rt, &sg);
						cst[c] = sg == 0;
					}
					else
					{
						unsigned int v = 0;
						st = variations_on(w->bG[c], w->dG[c], &rt, &v);
						if (v == 0)
							cst[c] = 0;
						else if (v == 1)
							cst[c] = 1;   /* G | R and R has one root here: it is G's */
					}
				}
			}
			if (st != QAWS_STATUS_OK)
				break;
			for (c = 0; c < ncons; c++)
			{
				if (cst[c] == 0)
					accept = 0;
				if (cst[c] < 0)
					undecided = 1;
			}
			if (accept == 0)
				break;
			if (inside == 1 && !undecided)
			{
				accept = 1;
				break;
			}
			if (rt.exact || rt.depth >= QAWS_EXACT_ROOTS_MAX_DEPTH)
				st = QAWS_STATUS_CERTIFICATION_FAILED;   /* singular point, tangency, span end at an irrational root, cusp */
			else
				st = qaws_exact_bernstein_refine(w->bR, N, &rt, rt.depth + 1);
		}
		if (st != QAWS_STATUS_OK || accept != 1)
			continue;
		if (!rt.exact)
		{
			/* tighten; past the integer budget the certified interval stays wider */
			qaws_exact_root keep = rt;
			st = qaws_exact_bernstein_refine(w->bR, N, &rt, QAWS_EXACT_ROOTS_MAX_DEPTH);
			if (st == QAWS_STATUS_EXACT_RANGE_EXCEEDED)
			{
				st = QAWS_STATUS_OK;
				if (rt.depth < keep.depth)
					rt = keep;
			}
		}
		if (st == QAWS_STATUS_OK && *count >= capacity)
			st = QAWS_STATUS_BUFFER_TOO_SMALL;
		if (st == QAWS_STATUS_OK)
		{
			qaws_exact_local_hit* h = &out[*count];
			qaws_exact_int c00[CC_COEF], c01[CC_COEF];
			int at_end = rt.exact && rt.index == ((uint64_t)1 << rt.depth);
			unsigned int i;
			h->r = rt;
			st = qaws_exact_bernstein_restrict_pair(w->b0, w->b1, N1, rt.index - (uint64_t)at_end, rt.depth, c00, c01);
			if (st != QAWS_STATUS_OK)
				break;
			if (s_end >= 0)
			{
				h->s_exact = 1;
				qaws_exact_int_from_i64(&h->s_num, s_end);
				qaws_exact_int_from_i64(&h->s_den, 1);
			}
			else if (rt.exact)
			{
				unsigned int e = at_end ? N1 : 0;
				h->s_exact = 1;
				h->s_num = c01[e];
				h->s_den = c00[e];
				if (qaws_exact_int_sign(&h->s_den) < 0)
				{
					qaws_exact_int_neg(&h->s_num, &h->s_num);
					qaws_exact_int_neg(&h->s_den, &h->s_den);
				}
			}
			else
			{
				/* s = sum c01_i B_i / sum c00_i B_i, c00 of one sign: within the coefficient ratios */
				h->s_exact = 0;
				h->s_lo = HUGE_VAL;
				h->s_hi = -HUGE_VAL;
				for (i = 0; i <= N1; i++)
				{
					double lo, hi;
					int ex;
					qaws_exact_ratio_enclose(&c01[i], &c00[i], &lo, &hi, &ex);
					if (lo < h->s_lo) h->s_lo = lo;
					if (hi > h->s_hi) h->s_hi = hi;
				}
				if (h->s_lo < 0) h->s_lo = 0;
				if (h->s_hi > 1) h->s_hi = 1;
			}
			(*count)++;
		}
	}
	sv_free(w);
	return st;
}
