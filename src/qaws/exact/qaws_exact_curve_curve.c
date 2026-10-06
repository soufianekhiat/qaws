#include "qaws_exact_curve.h"
#include "qaws_exact_roots.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define CC_MAX_IMPLICIT 6
#define CC_MAX_COEF (QAWS_EXACT_ROOTS_MAX_DEGREE + 1)
#define CC_MAX_ROOTS (QAWS_EXACT_ROOTS_MAX_DEGREE + 2)

/* ------------------------------------------------------------------ */
/*  Integer polynomials in the power basis                             */
/* ------------------------------------------------------------------ */

typedef struct poly
{
	unsigned int deg;
	qaws_exact_int c[CC_MAX_COEF];
} poly;

static void poly_zero(poly* p, unsigned int deg)
{
	unsigned int i;
	p->deg = deg;
	for (i = 0; i <= deg; i++)
		qaws_exact_int_zero(&p->c[i]);
}

static qaws_status poly_mul(poly* r, poly const* a, poly const* b)
{
	qaws_exact_int t;
	unsigned int i, j;
	qaws_status st;
	if (a->deg + b->deg >= CC_MAX_COEF)
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	poly_zero(r, a->deg + b->deg);
	for (i = 0; i <= a->deg; i++)
	{
		if (qaws_exact_int_is_zero(&a->c[i]))
			continue;
		for (j = 0; j <= b->deg; j++)
		{
			TRY(qaws_exact_int_mul(&t, &a->c[i], &b->c[j]));
			TRY(qaws_exact_int_add(&r->c[i + j], &r->c[i + j], &t));
		}
	}
	return QAWS_STATUS_OK;
}

/* r += sign a (r's degree grows to a's) */
static qaws_status poly_acc(poly* r, poly const* a, int sign)
{
	unsigned int i;
	qaws_status st;
	while (r->deg < a->deg)
		qaws_exact_int_zero(&r->c[++r->deg]);
	for (i = 0; i <= a->deg; i++)
		TRY(sign > 0 ? qaws_exact_int_add(&r->c[i], &r->c[i], &a->c[i]) : qaws_exact_int_sub(&r->c[i], &r->c[i], &a->c[i]));
	return QAWS_STATUS_OK;
}

static int poly_is_zero(poly const* p)
{
	unsigned int i;
	for (i = 0; i <= p->deg; i++)
		if (!qaws_exact_int_is_zero(&p->c[i]))
			return 0;
	return 1;
}

static int64_t binom64(unsigned int n, unsigned int k)
{
	int64_t r = 1;
	unsigned int i;
	if (k > n)
		return 0;
	if (k > n - k)
		k = n - k;
	for (i = 1; i <= k; i++)
		r = r * (int64_t)(n - k + i) / (int64_t)i;
	return r;
}

/* Bernstein h_i (stride D, component c) on [0, 1] -> power basis. */
static qaws_status bernstein_to_power(qaws_exact_int const* h, unsigned int n, unsigned int D, unsigned int c, poly* out)
{
	qaws_exact_int t;
	unsigned int i, k;
	qaws_status st;
	poly_zero(out, n);
	for (k = 0; k <= n; k++)
		for (i = 0; i <= k; i++)
		{
			int64_t f = binom64(n, i) * binom64(n - i, k - i) * (((k - i) & 1) ? -1 : 1);
			TRY(qaws_exact_int_mul_i64(&t, &h[i * D + c], f));
			TRY(qaws_exact_int_add(&out->c[k], &out->c[k], &t));
		}
	return QAWS_STATUS_OK;
}

/* Power basis -> Bernstein of degree N = p->deg, times L = lcm_j C(N, j): b_i = sum_{j <= i} C(i, j) (L / C(N, j)) p_j. */
static qaws_status power_to_bernstein(poly const* p, qaws_exact_int* b)
{
	qaws_exact_int L, g, q, t, cj;
	unsigned int N = p->deg, i, j;
	qaws_status st;
	qaws_exact_int_from_i64(&L, 1);
	for (j = 0; j <= N; j++)
	{
		qaws_exact_int_from_i64(&cj, binom64(N, j));
		qaws_exact_int_gcd(&g, &L, &cj);
		TRY(qaws_exact_int_divmod(&q, NULL, &cj, &g));
		TRY(qaws_exact_int_mul(&L, &L, &q));
	}
	for (i = 0; i <= N; i++)
		qaws_exact_int_zero(&b[i]);
	for (j = 0; j <= N; j++)
	{
		if (qaws_exact_int_is_zero(&p->c[j]))
			continue;
		qaws_exact_int_from_i64(&cj, binom64(N, j));
		TRY(qaws_exact_int_divmod(&q, NULL, &L, &cj));
		TRY(qaws_exact_int_mul(&q, &q, &p->c[j]));
		for (i = j; i <= N; i++)
		{
			TRY(qaws_exact_int_mul_i64(&t, &q, binom64(i, j)));
			TRY(qaws_exact_int_add(&b[i], &b[i], &t));
		}
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Implicitization                                                    */
/* ------------------------------------------------------------------ */

/*
 * Bezout matrix of P, Q (degree n, power basis):
 * (P(t) Q(s) - P(s) Q(t)) / (t - s) = sum B_ij t^i s^j. A term
 * c_ab (t^a s^b - t^b s^a), a > b, c_ab = p_a q_b - p_b q_a, divides into
 * t^b s^b sum_{k < a - b} t^(a-b-1-k) s^k.
 */
static qaws_status bezout(poly const* P, poly const* Q, unsigned int n, qaws_exact_int* B)
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
			{
				unsigned int i = a - 1 - k, j = b + k;
				TRY(qaws_exact_int_add(&B[i * n + j], &B[i * n + j], &c));
			}
		}
	return QAWS_STATUS_OK;
}

/* Laplace expansion over rows row..n-1 and the columns outside mask. */
static qaws_status det_rec(poly const* M, unsigned int n, unsigned int row, unsigned int mask, poly* out)
{
	poly* sub;
	poly* t;
	unsigned int j;
	int sign = 1;
	qaws_status st = QAWS_STATUS_OK;
	if (row == n)
	{
		poly_zero(out, 0);
		qaws_exact_int_from_i64(&out->c[0], 1);
		return QAWS_STATUS_OK;
	}
	sub = (poly*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(poly) * 2));
	if (!sub)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	t = sub + 1;
	poly_zero(out, 0);
	for (j = 0; j < n && st == QAWS_STATUS_OK; j++)
	{
		if (mask & (1u << j))
			continue;
		st = det_rec(M, n, row + 1, mask | (1u << j), sub);
		if (st == QAWS_STATUS_OK)
			st = poly_mul(t, &M[row * n + j], sub);
		if (st == QAWS_STATUS_OK)
			st = poly_acc(out, t, sign);
		sign = -sign;
	}
	qaws_internal_dealloc(NULL, sub);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Span against span                                                  */
/* ------------------------------------------------------------------ */

/* Signs on the root interval: +1 / -1 when every coefficient there has that strict sign, 0 otherwise. */
static qaws_status sign_on(qaws_exact_int const* b, unsigned int n, qaws_exact_root const* r, int* out)
{
	qaws_exact_int c[CC_MAX_COEF];
	unsigned int i;
	int s;
	qaws_status st;
	if (r->exact)
	{
		/* the value at the point (a positive multiple): the first coefficient of
		   the interval starting there, or the last of the one ending at s = 1 */
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

/* Enclosure of the rational num / den (den != 0) by doubles; exact when dyadic and representable. */
static void ratio_enclose(qaws_exact_int const* num, qaws_exact_int const* den, double* lo, double* hi, int* exact)
{
	qaws_exact_int g, n2, d2;
	double v;
	unsigned int db;
	qaws_exact_int_gcd(&g, num, den);
	if (qaws_exact_int_is_zero(num))
	{
		*lo = *hi = 0;
		*exact = 1;
		return;
	}
	qaws_exact_int_divmod(&n2, NULL, num, &g);
	qaws_exact_int_divmod(&d2, NULL, den, &g);
	if (qaws_exact_int_sign(&d2) < 0)
	{
		qaws_exact_int_neg(&n2, &n2);
		qaws_exact_int_neg(&d2, &d2);
	}
	v = qaws_exact_ratio_to_double(&n2, &d2);
	db = qaws_exact_int_bits(&d2);
	{
		/* den a power of two and the numerator within 53 bits */
		qaws_exact_int p2;
		qaws_exact_int_from_i64(&p2, 1);
		qaws_exact_int_shl(&p2, &p2, db - 1);
		*exact = qaws_exact_int_cmp_abs(&p2, &d2) == 0 && qaws_exact_int_bits(&n2) <= 53 && fabs(v) >= 2.2250738585072014e-308;
	}
	*lo = *exact ? v : nextafter(v, -HUGE_VAL);
	*hi = *exact ? v : nextafter(v, HUGE_VAL);
}

/* Local s in [s_lo, s_hi] of a span -> the curve parameter, rounded outward. */
static void span_param_enclose(qaws_exact_span const* sp, int shift, double s_lo, double s_hi, int exact, double* t_lo, double* t_hi)
{
	double a = ldexp((double)sp->a, -shift), L = ldexp((double)(sp->b - sp->a), -shift);
	*t_lo = a + L * s_lo;
	*t_hi = a + L * s_hi;
	if (!exact)
	{
		*t_lo = nextafter(nextafter(*t_lo, -HUGE_VAL), -HUGE_VAL);
		*t_hi = nextafter(nextafter(*t_hi, HUGE_VAL), HUGE_VAL);
	}
}

typedef struct span_hit
{
	double imp_lo, imp_hi;   /* implicitized span's parameter */
	double sub_lo, sub_hi;   /* substituted span's parameter */
	int exact;
} span_hit;

/*
 * Intersections of span I (implicitized, degree n >= 2) with span S
 * (substituted, degree m): g(r) = det M(S(r)) and the parameter of I at a
 * root, s = C01 / C00 (the adjugate's first column is the kernel
 * (1, s, s^2, ...) of the symmetric Bezout matrix).
 */
static qaws_status span_pair(qaws_exact_curve const* ci, unsigned int si, qaws_exact_curve const* cs, unsigned int ss, int skip_sub_start,
	int skip_imp_start, span_hit* out, unsigned int capacity, unsigned int* count)
{
	qaws_exact_span const* I = &ci->spans[si];
	qaws_exact_span const* S = &cs->spans[ss];
	unsigned int n = I->degree, m = S->degree, i, j, k, nr = 0, N, N1;
	poly XI, YI, WI, XS, YS, WS;
	qaws_exact_int BYW[CC_MAX_IMPLICIT * CC_MAX_IMPLICIT], BWX[CC_MAX_IMPLICIT * CC_MAX_IMPLICIT], BXY[CC_MAX_IMPLICIT * CC_MAX_IMPLICIT];
	poly* M = NULL;
	poly* g = NULL;
	poly* C = NULL;   /* C[0] = C00, C[1] = C01, C[2] = C00 - C01 */
	qaws_exact_int* bg = NULL;
	qaws_exact_int* bc = NULL;   /* 3 rows of CC_MAX_COEF */
	qaws_exact_root roots[CC_MAX_ROOTS];
	qaws_status st = QAWS_STATUS_OK;
	if (n < 2 || n > CC_MAX_IMPLICIT || n * m > QAWS_EXACT_ROOTS_MAX_DEGREE)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	TRY(bernstein_to_power(I->h, n, 3, 0, &XI));
	TRY(bernstein_to_power(I->h, n, 3, 1, &YI));
	TRY(bernstein_to_power(I->h, n, 3, 2, &WI));
	TRY(bernstein_to_power(S->h, m, 3, 0, &XS));
	TRY(bernstein_to_power(S->h, m, 3, 1, &YS));
	TRY(bernstein_to_power(S->h, m, 3, 2, &WS));
	TRY(bezout(&YI, &WI, n, BYW));
	TRY(bezout(&WI, &XI, n, BWX));
	TRY(bezout(&XI, &YI, n, BXY));
	M = (poly*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(poly) * (n * n + 1 + 3)));
	bg = (qaws_exact_int*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_int) * CC_MAX_COEF * 4));
	if (!M || !bg)
	{
		qaws_internal_dealloc(NULL, M);
		qaws_internal_dealloc(NULL, bg);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	g = M + n * n;
	C = g + 1;
	bc = bg + CC_MAX_COEF;
	/* M_ij(r) = X(r) BYW_ij + Y(r) BWX_ij + W(r) BXY_ij */
	for (i = 0; i < n * n && st == QAWS_STATUS_OK; i++)
	{
		poly_zero(&M[i], m);
		for (k = 0; k <= m && st == QAWS_STATUS_OK; k++)
		{
			qaws_exact_int t;
			st = qaws_exact_int_mul(&t, &XS.c[k], &BYW[i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&M[i].c[k], &M[i].c[k], &t);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &YS.c[k], &BWX[i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&M[i].c[k], &M[i].c[k], &t);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &WS.c[k], &BXY[i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&M[i].c[k], &M[i].c[k], &t);
		}
	}
	/* cofactors of row 0 and the determinant */
	if (st == QAWS_STATUS_OK)
	{
		poly* minor = (poly*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(poly) * 2));
		if (!minor)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			poly_zero(g, 0);
			for (j = 0; j < n && st == QAWS_STATUS_OK; j++)
			{
				st = det_rec(M, n, 1, 1u << j, &minor[0]);
				if (st == QAWS_STATUS_OK && j < 2)
				{
					C[j] = minor[0];
					if (j == 1)
						for (k = 0; k <= C[1].deg; k++)
							qaws_exact_int_neg(&C[1].c[k], &C[1].c[k]);
				}
				if (st == QAWS_STATUS_OK) st = poly_mul(&minor[1], &M[j], &minor[0]);
				if (st == QAWS_STATUS_OK) st = poly_acc(g, &minor[1], (j & 1) ? -1 : 1);
			}
			qaws_internal_dealloc(NULL, minor);
		}
	}
	if (st == QAWS_STATUS_OK)
	{
		C[2] = C[0];
		st = poly_acc(&C[2], &C[1], -1);
	}
	if (st == QAWS_STATUS_OK && poly_is_zero(g))
		st = QAWS_STATUS_CERTIFICATION_FAILED;   /* a common component */
	/* fixed degrees: g to n m, the cofactors to (n - 1) m */
	N = n * m;
	N1 = (n - 1) * m;
	for (k = 0; k < 3 && st == QAWS_STATUS_OK; k++)
		while (C[k].deg < N1)
			qaws_exact_int_zero(&C[k].c[++C[k].deg]);
	while (st == QAWS_STATUS_OK && g->deg < N)
		qaws_exact_int_zero(&g->c[++g->deg]);
	if (st == QAWS_STATUS_OK && (g->deg != N || C[0].deg != N1))
		st = QAWS_STATUS_INTERNAL_ERROR;
	if (st == QAWS_STATUS_OK) st = power_to_bernstein(g, bg);
	for (k = 0; k < 3 && st == QAWS_STATUS_OK; k++)
		st = power_to_bernstein(&C[k], &bc[k * CC_MAX_COEF]);
	if (st == QAWS_STATUS_OK)
		st = qaws_exact_bernstein_isolate(bg, N, roots, CC_MAX_ROOTS, &nr);

	for (k = 0; k < nr && st == QAWS_STATUS_OK; k++)
	{
		qaws_exact_root rt = roots[k];
		int s00 = 0, s01 = 0, sd = 0, inside = -1;
		if (rt.exact && rt.index == 0 && skip_sub_start)
			continue;
		/* decide 0 <= C01 / C00 <= 1, refining until the signs settle */
		while (st == QAWS_STATUS_OK && inside < 0)
		{
			st = sign_on(&bc[0], N1, &rt, &s00);
			if (st == QAWS_STATUS_OK) st = sign_on(&bc[CC_MAX_COEF], N1, &rt, &s01);
			if (st == QAWS_STATUS_OK) st = sign_on(&bc[2 * CC_MAX_COEF], N1, &rt, &sd);
			if (st != QAWS_STATUS_OK)
				break;
			if (s00 != 0 && ((s01 != 0 && s01 != s00) || (sd != 0 && sd != s00)))
				inside = 0;
			else if (s00 != 0 && (s01 == s00 || (rt.exact && s01 == 0)) && (sd == s00 || (rt.exact && sd == 0)))
				inside = 1;
			else if (rt.exact || rt.depth >= QAWS_EXACT_ROOTS_MAX_DEPTH)
				st = QAWS_STATUS_CERTIFICATION_FAILED;   /* singular point, or a span end at an irrational root */
			else
				st = qaws_exact_bernstein_refine(bg, N, &rt, rt.depth + 1);
		}
		if (st != QAWS_STATUS_OK || inside != 1)
			continue;
		if (!rt.exact)
			st = qaws_exact_bernstein_refine(bg, N, &rt, QAWS_EXACT_ROOTS_MAX_DEPTH);
		if (st != QAWS_STATUS_OK)
			break;
		{
			span_hit h;
			qaws_exact_int c00[CC_MAX_COEF], c01[CC_MAX_COEF];
			double s_lo = HUGE_VAL, s_hi = -HUGE_VAL, r_lo, r_hi = 0;
			int ex_lo, ex_hi = 1, s_exact = 0, at_end = rt.exact && rt.index == ((uint64_t)1 << rt.depth);
			/* at an exact root: the interval starting there (ending there at r = 1) */
			st = qaws_exact_bernstein_restrict_pair(&bc[0], &bc[CC_MAX_COEF], N1, rt.index - (uint64_t)at_end, rt.depth, c00, c01);
			if (st != QAWS_STATUS_OK)
				break;
			if (rt.exact)
			{
				unsigned int e = at_end ? N1 : 0;
				ratio_enclose(&c01[e], &c00[e], &s_lo, &s_hi, &s_exact);
				if (skip_imp_start && qaws_exact_int_is_zero(&c01[e]))
					continue;   /* s = 0: the end of the previous implicitized span */
			}
			else
				/* s = sum c01_i B_i / sum c00_i B_i, c00 of one sign: within the coefficient ratios */
				for (i = 0; i <= N1; i++)
				{
					double lo, hi;
					int ex;
					ratio_enclose(&c01[i], &c00[i], &lo, &hi, &ex);
					if (lo < s_lo) s_lo = lo;
					if (hi > s_hi) s_hi = hi;
				}
			if (s_lo < 0) s_lo = 0;
			if (s_hi > 1) s_hi = 1;
			span_param_enclose(I, ci->param_shift, s_lo, s_hi, s_exact, &h.imp_lo, &h.imp_hi);
			st = qaws_exact_span_param_to_double(S, cs->param_shift, rt.index, rt.depth, &r_lo, &ex_lo);
			if (st == QAWS_STATUS_OK && !rt.exact)
				st = qaws_exact_span_param_to_double(S, cs->param_shift, rt.index + 1, rt.depth, &r_hi, &ex_hi);
			if (st != QAWS_STATUS_OK)
				break;
			if (rt.exact)
			{
				r_hi = r_lo;
				if (!ex_lo)
				{
					r_lo = nextafter(r_lo, -HUGE_VAL);
					r_hi = nextafter(r_hi, HUGE_VAL);
				}
			}
			else
			{
				if (!ex_lo) r_lo = nextafter(r_lo, -HUGE_VAL);
				if (!ex_hi) r_hi = nextafter(r_hi, HUGE_VAL);
			}
			h.sub_lo = r_lo;
			h.sub_hi = r_hi;
			h.exact = rt.exact && ex_lo && s_exact;
			if (*count >= capacity)
				st = QAWS_STATUS_BUFFER_TOO_SMALL;
			else
				out[(*count)++] = h;
		}
	}
	qaws_internal_dealloc(NULL, M);
	qaws_internal_dealloc(NULL, bg);
	return st;
}

/* Two segments (degree-1 spans, homogeneous): the exact linear roots. */
static qaws_status segment_pair(qaws_exact_curve const* ca, unsigned int sa, qaws_exact_curve const* cb, unsigned int sb, int skip_a, int skip_b,
	span_hit* out, unsigned int capacity, unsigned int* count)
{
	qaws_exact_span const* A = &ca->spans[sa];
	qaws_exact_span const* B = &cb->spans[sb];
	qaws_exact_int fa[2], fb[2], num, den, t1, t2;
	unsigned int e, c;
	int sa_ex, sb_ex;
	double a_lo, a_hi, b_lo, b_hi, s_lo, s_hi, r_lo, r_hi;
	qaws_status st;
	/* f(Q) = det(A0, A1, Q) (homogeneous): zero on A's line. fb_e = f(B_e), fa_e likewise. */
	for (e = 0; e < 2; e++)
	{
		qaws_exact_int const* P0 = &A->h[0];
		qaws_exact_int const* P1 = &A->h[3];
		qaws_exact_int const* Q = &B->h[e * 3];
		qaws_exact_int_zero(&fb[e]);
		for (c = 0; c < 3; c++)
		{
			unsigned int c1 = (c + 1) % 3, c2 = (c + 2) % 3;
			TRY(qaws_exact_int_mul(&t1, &P0[c1], &P1[c2]));
			TRY(qaws_exact_int_mul(&t2, &P0[c2], &P1[c1]));
			TRY(qaws_exact_int_sub(&t1, &t1, &t2));
			TRY(qaws_exact_int_mul(&t1, &t1, &Q[c]));
			TRY(qaws_exact_int_add(&fb[e], &fb[e], &t1));
		}
		P0 = &B->h[0];
		P1 = &B->h[3];
		Q = &A->h[e * 3];
		qaws_exact_int_zero(&fa[e]);
		for (c = 0; c < 3; c++)
		{
			unsigned int c1 = (c + 1) % 3, c2 = (c + 2) % 3;
			TRY(qaws_exact_int_mul(&t1, &P0[c1], &P1[c2]));
			TRY(qaws_exact_int_mul(&t2, &P0[c2], &P1[c1]));
			TRY(qaws_exact_int_sub(&t1, &t1, &t2));
			TRY(qaws_exact_int_mul(&t1, &t1, &Q[c]));
			TRY(qaws_exact_int_add(&fa[e], &fa[e], &t1));
		}
	}
	/* f along a span is (1 - u) f_0 + u f_1 (homogeneous Bernstein): root u = f0 / (f0 - f1) */
	if (qaws_exact_int_is_zero(&fb[0]) && qaws_exact_int_is_zero(&fb[1]))
		return QAWS_STATUS_CERTIFICATION_FAILED;   /* collinear */
	if (qaws_exact_int_sign(&fb[0]) * qaws_exact_int_sign(&fb[1]) > 0 || qaws_exact_int_sign(&fa[0]) * qaws_exact_int_sign(&fa[1]) > 0)
		return QAWS_STATUS_OK;
	TRY(qaws_exact_int_sub(&den, &fb[0], &fb[1]));
	ratio_enclose(&fb[0], &den, &r_lo, &r_hi, &sb_ex);
	if (skip_b && qaws_exact_int_is_zero(&fb[0]))
		return QAWS_STATUS_OK;
	num = fa[0];
	TRY(qaws_exact_int_sub(&den, &fa[0], &fa[1]));
	ratio_enclose(&num, &den, &s_lo, &s_hi, &sa_ex);
	if (skip_a && qaws_exact_int_is_zero(&fa[0]))
		return QAWS_STATUS_OK;
	span_param_enclose(A, ca->param_shift, s_lo, s_hi, sa_ex, &a_lo, &a_hi);
	span_param_enclose(B, cb->param_shift, r_lo, r_hi, sb_ex, &b_lo, &b_hi);
	if (*count >= capacity)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	out[*count].imp_lo = a_lo;
	out[*count].imp_hi = a_hi;
	out[*count].sub_lo = b_lo;
	out[*count].sub_hi = b_hi;
	out[*count].exact = sa_ex && sb_ex;
	(*count)++;
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_curve_curve_hits(qaws_exact_curve const* a, qaws_exact_curve const* b, qaws_exact_pair* out_pairs, unsigned int capacity,
	unsigned int* out_count)
{
	unsigned int ia, ib, k, count = 0;
	qaws_status st = QAWS_STATUS_OK;
	span_hit hits[CC_MAX_ROOTS];
	if (!a || !b || !out_count || (!out_pairs && capacity) || a->dimension != 2 || b->dimension != 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (a->space_exp2 != b->space_exp2)
		return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
	for (ia = 0; ia < a->span_count && st == QAWS_STATUS_OK; ia++)
		for (ib = 0; ib < b->span_count && st == QAWS_STATUS_OK; ib++)
		{
			unsigned int na = a->spans[ia].degree, nb = b->spans[ib].degree, nh = 0;
			int swap;
			if (na == 1 && nb == 1)
			{
				st = segment_pair(a, ia, b, ib, ia > 0, ib > 0, hits, CC_MAX_ROOTS, &nh);
				swap = 0;
			}
			else
			{
				/* implicitize the lower degree (at least 2); when that span is
				   degenerate (its image of lower degree, e.g. a line traced by a
				   quadratic: a zero or squared implicit form), the other one */
				swap = !(na >= 2 && (na <= nb || nb < 2));
				st = swap ? span_pair(b, ib, a, ia, ia > 0, ib > 0, hits, CC_MAX_ROOTS, &nh)
				          : span_pair(a, ia, b, ib, ib > 0, ia > 0, hits, CC_MAX_ROOTS, &nh);
				if ((st == QAWS_STATUS_CERTIFICATION_FAILED || st == QAWS_STATUS_EXACT_UNSUPPORTED) && (swap ? na : nb) >= 2)
				{
					qaws_status st2;
					nh = 0;
					swap = !swap;
					st2 = swap ? span_pair(b, ib, a, ia, ia > 0, ib > 0, hits, CC_MAX_ROOTS, &nh)
					           : span_pair(a, ia, b, ib, ib > 0, ia > 0, hits, CC_MAX_ROOTS, &nh);
					if (st2 == QAWS_STATUS_OK)
						st = st2;
				}
			}
			for (k = 0; k < nh && st == QAWS_STATUS_OK; k++)
			{
				qaws_exact_pair p;
				p.kind = hits[k].exact ? QAWS_EXACT_HIT_POINT : QAWS_EXACT_HIT_CROSSING;
				p.a_lo = swap ? hits[k].sub_lo : hits[k].imp_lo;
				p.a_hi = swap ? hits[k].sub_hi : hits[k].imp_hi;
				p.b_lo = swap ? hits[k].imp_lo : hits[k].sub_lo;
				p.b_hi = swap ? hits[k].imp_hi : hits[k].sub_hi;
				if (count >= capacity)
					st = QAWS_STATUS_BUFFER_TOO_SMALL;
				else
					out_pairs[count++] = p;
			}
		}
	/* sort by a's parameter */
	for (ia = 1; ia < count; ia++)
	{
		qaws_exact_pair key = out_pairs[ia];
		int j = (int)ia - 1;
		while (j >= 0 && out_pairs[j].a_lo > key.a_lo)
		{
			out_pairs[j + 1] = out_pairs[j];
			j--;
		}
		out_pairs[j + 1] = key;
	}
	*out_count = count;
	return st;
}
