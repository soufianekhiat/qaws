#include "qaws_exact_curve.h"
#include "qaws_exact_roots.h"
#include "qaws_exact_poly.h"
#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define CC_MAX_IMPLICIT 6
#define CC_COEF QAWS_EXACT_POLY_COEF
#define CC_MAX_ROOTS QAWS_EXACT_SOLVE_MAX_ROOTS

typedef qaws_exact_poly poly;

static void* cc_alloc(size_t bytes)
{
	return qaws_internal_alloc(NULL, (unsigned long)bytes);
}

static void cc_free(void* p)
{
	qaws_internal_dealloc(NULL, p);
}


/* ------------------------------------------------------------------ */
/*  Parameters                                                         */
/* ------------------------------------------------------------------ */

static qaws_status map_r(qaws_exact_span const* sp, int shift, qaws_exact_root const* rt, double* lo, double* hi, int* exact)
{
	int ex_lo, ex_hi;
	qaws_status st;
	TRY(qaws_exact_span_param_to_double(sp, shift, rt->index, rt->depth, lo, &ex_lo));
	if (rt->exact)
	{
		*hi = *lo;
		if (!ex_lo)
		{
			*lo = nextafter(*lo, -HUGE_VAL);
			*hi = nextafter(*hi, HUGE_VAL);
		}
		*exact = ex_lo;
		return QAWS_STATUS_OK;
	}
	TRY(qaws_exact_span_param_to_double(sp, shift, rt->index + 1, rt->depth, hi, &ex_hi));
	if (!ex_lo) *lo = nextafter(*lo, -HUGE_VAL);
	if (!ex_hi) *hi = nextafter(*hi, HUGE_VAL);
	*exact = 0;
	return QAWS_STATUS_OK;
}

/* t = (a + (b - a) s) 2^-shift for the local s of a hit, rounded outward (exact when it is a double). */
static qaws_status map_s(qaws_exact_span const* sp, int shift, qaws_exact_local_hit const* h, double* lo, double* hi, int* exact)
{
	if (h->s_exact)
	{
		qaws_exact_int num, den, t;
		qaws_status st;
		TRY(qaws_exact_int_mul_i64(&num, &h->s_den, sp->a));
		TRY(qaws_exact_int_mul_i64(&t, &h->s_num, sp->b - sp->a));
		TRY(qaws_exact_int_add(&num, &num, &t));
		TRY(qaws_exact_int_shl(&den, &h->s_den, (unsigned int)shift));
		qaws_exact_ratio_enclose(&num, &den, lo, hi, exact);
		return QAWS_STATUS_OK;
	}
	{
		double a = ldexp((double)sp->a, -shift), L = ldexp((double)(sp->b - sp->a), -shift);
		*lo = nextafter(nextafter(a + L * h->s_lo, -HUGE_VAL), -HUGE_VAL);
		*hi = nextafter(nextafter(a + L * h->s_hi, HUGE_VAL), HUGE_VAL);
		*exact = 0;
	}
	return QAWS_STATUS_OK;
}

static int s_is(qaws_exact_local_hit const* h, int v)
{
	if (!h->s_exact)
		return 0;
	return v == 0 ? qaws_exact_int_is_zero(&h->s_num) : qaws_exact_int_cmp(&h->s_num, &h->s_den) == 0;
}

/* ------------------------------------------------------------------ */
/*  Systems                                                            */
/* ------------------------------------------------------------------ */

/* Components of a projection: the plane (c0, c1), the weight w, and the checked coordinate cz (-1: none). */
typedef struct proj
{
	unsigned int c0, c1, w;
	int cz;
} proj;

static unsigned int projections(unsigned int dim, proj* out)
{
	if (dim == 2)
	{
		out[0].c0 = 0; out[0].c1 = 1; out[0].w = 2; out[0].cz = -1;
		return 1;
	}
	out[0].c0 = 0; out[0].c1 = 1; out[0].w = 3; out[0].cz = 2;
	out[1].c0 = 0; out[1].c1 = 2; out[1].w = 3; out[1].cz = 1;
	out[2].c0 = 1; out[2].c1 = 2; out[2].w = 3; out[2].cz = 0;
	return 3;
}

/*
 * Span I implicitized in the projection, span S substituted:
 * M(x, y, w) = x Bez(Y, W) + y Bez(W, X) + w Bez(X, Y) is singular on I's
 * curve; R(r) = det M(S(r)), s = C01 / C00 (the kernel (1, s, s^2, ...) of
 * the symmetric M). H (3D): Z_I(s) W_S(r) - Z_S(r) W_I(s) times C00^n at s.
 */
static qaws_status build_cross(qaws_exact_span const* I, qaws_exact_span const* S, unsigned int D, proj const* pj, poly* R, poly* C00, poly* C01,
	poly* H)
{
	unsigned int n = I->degree, m = S->degree, i, k;
	poly* P;   /* XI YI WI XS YS WS ZI ZS, then M (n n), then Hc (n + 1) */
	qaws_exact_int* B;
	qaws_status st = QAWS_STATUS_OK;
	P = (poly*)cc_alloc(sizeof(poly) * (8 + n * n + n + 1));
	B = (qaws_exact_int*)cc_alloc(sizeof(qaws_exact_int) * 3 * n * n);
	if (!P || !B)
	{
		cc_free(P);
		cc_free(B);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	st = qaws_exact_poly_from_bernstein(I->h, n, D, pj->c0, &P[0]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(I->h, n, D, pj->c1, &P[1]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(I->h, n, D, pj->w, &P[2]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(S->h, m, D, pj->c0, &P[3]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(S->h, m, D, pj->c1, &P[4]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(S->h, m, D, pj->w, &P[5]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&P[1], &P[2], n, &B[0]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&P[2], &P[0], n, &B[n * n]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&P[0], &P[1], n, &B[2 * n * n]);
	for (i = 0; i < n * n && st == QAWS_STATUS_OK; i++)
	{
		poly* Mi = &P[8 + i];
		qaws_exact_poly_zero(Mi, m);
		for (k = 0; k <= m && st == QAWS_STATUS_OK; k++)
		{
			qaws_exact_int t;
			st = qaws_exact_int_mul(&t, &P[3].c[k], &B[i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&Mi->c[k], &Mi->c[k], &t);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &P[4].c[k], &B[n * n + i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&Mi->c[k], &Mi->c[k], &t);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &P[5].c[k], &B[2 * n * n + i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&Mi->c[k], &Mi->c[k], &t);
		}
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_exact_det_cofactors(&P[8], n, R, C00, C01);
	if (st == QAWS_STATUS_OK && pj->cz >= 0 && H)
	{
		/* H_i(r) = zI_i WS(r) - wI_i ZS(r) */
		poly* Hc = &P[8 + n * n];
		st = qaws_exact_poly_from_bernstein(I->h, n, D, (unsigned int)pj->cz, &P[6]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(S->h, m, D, (unsigned int)pj->cz, &P[7]);
		for (i = 0; i <= n && st == QAWS_STATUS_OK; i++)
		{
			poly t;
			st = qaws_exact_poly_scale(&Hc[i], &P[5], &P[6].c[i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_scale(&t, &P[7], &P[2].c[i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(&Hc[i], &t, -1);
		}
		if (st == QAWS_STATUS_OK)
			st = qaws_exact_substitute(Hc, n, C01, C00, H);
	}
	cc_free(P);
	cc_free(B);
	return st;
}

/*
 * Self-intersections inside one span: s != t with C(s) = C(t). The divided
 * differences P(s, t) = (X(s) W(t) - X(t) W(s)) / (s - t) = sum Bez(X, W)_ij
 * s^i t^j (and Q from Y) vanish there; their Bezout matrix in s (entries
 * polynomials in t) gives R(t) and s = C01 / C00. Hgt = C01 - t C00 has the
 * sign of C00 exactly when s > t; Hz the third divided difference (3D).
 */
static qaws_status build_self(qaws_exact_span const* sp, unsigned int D, proj const* pj, poly* R, poly* C00, poly* C01, poly* Hgt, poly* Hz)
{
	unsigned int n = sp->degree, k = n - 1, i, j;
	poly* P;   /* X Y W Z, then Ps (k + 1), Qs (k + 1), Zs (k + 1), E (k k) */
	qaws_exact_int* B;
	qaws_status st = QAWS_STATUS_OK;
	P = (poly*)cc_alloc(sizeof(poly) * (4 + 3 * (k + 1) + k * k));
	B = (qaws_exact_int*)cc_alloc(sizeof(qaws_exact_int) * 3 * n * n);
	if (!P || !B)
	{
		cc_free(P);
		cc_free(B);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	st = qaws_exact_poly_from_bernstein(sp->h, n, D, pj->c0, &P[0]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(sp->h, n, D, pj->c1, &P[1]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(sp->h, n, D, pj->w, &P[2]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&P[0], &P[2], n, &B[0]);
	if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&P[1], &P[2], n, &B[n * n]);
	if (st == QAWS_STATUS_OK && pj->cz >= 0)
	{
		st = qaws_exact_poly_from_bernstein(sp->h, n, D, (unsigned int)pj->cz, &P[3]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&P[3], &P[2], n, &B[2 * n * n]);
	}
	/* coefficient of s^i: polynomials in t of degree k */
	for (i = 0; i <= k && st == QAWS_STATUS_OK; i++)
	{
		poly* Ps = &P[4 + i];
		poly* Qs = &P[4 + (k + 1) + i];
		poly* Zs = &P[4 + 2 * (k + 1) + i];
		qaws_exact_poly_zero(Ps, k);
		qaws_exact_poly_zero(Qs, k);
		qaws_exact_poly_zero(Zs, k);
		for (j = 0; j <= k; j++)
		{
			Ps->c[j] = B[i * n + j];
			Qs->c[j] = B[n * n + i * n + j];
			if (pj->cz >= 0)
				Zs->c[j] = B[2 * n * n + i * n + j];
		}
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_exact_bezout_poly(&P[4], &P[4 + (k + 1)], k, &P[4 + 3 * (k + 1)]);
	if (st == QAWS_STATUS_OK)
		st = qaws_exact_det_cofactors(&P[4 + 3 * (k + 1)], k, R, C00, C01);
	if (st == QAWS_STATUS_OK)
	{
		/* Hgt = C01 - t C00 */
		poly tq;
		qaws_exact_poly_zero(&tq, 0);
		st = qaws_exact_poly_pad(&tq, C00->deg + 1);
		for (i = 0; i <= C00->deg && st == QAWS_STATUS_OK; i++)
			tq.c[i + 1] = C00->c[i];
		*Hgt = *C01;
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(Hgt, &tq, -1);
	}
	if (st == QAWS_STATUS_OK && pj->cz >= 0 && Hz)
		st = qaws_exact_substitute(&P[4 + 2 * (k + 1)], k, C01, C00, Hz);
	cc_free(P);
	cc_free(B);
	return st;
}

/* Is some coordinate (or the sum / difference of two) strictly monotone on the span? Then it cannot self-intersect. */
static qaws_status span_monotone(qaws_exact_span const* sp, unsigned int dim, int* out)
{
	unsigned int n = sp->degree, D = dim + 1, c, d, i;
	poly* P;   /* X_0.. X_dim-1, W, dX_c (dim), dW, N_c (dim), comb, scratch */
	qaws_exact_int b[CC_COEF];
	qaws_status st = QAWS_STATUS_OK;
	*out = 0;
	if (n <= 1)
	{
		*out = 1;
		return QAWS_STATUS_OK;
	}
	if (2 * n - 1 > QAWS_EXACT_ROOTS_MAX_DEGREE)
		return QAWS_STATUS_OK;
	P = (poly*)cc_alloc(sizeof(poly) * (3 * dim + 4));
	if (!P)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (c = 0; c <= dim && st == QAWS_STATUS_OK; c++)
		st = qaws_exact_poly_from_bernstein(sp->h, n, D, c, &P[c]);   /* P[dim] = W */
	/* N_c = X_c' W - X_c W' */
	{
		poly* dW = &P[2 * dim + 1];
		qaws_exact_poly_zero(dW, n - 1);
		for (i = 1; i <= n && st == QAWS_STATUS_OK; i++)
			st = qaws_exact_int_mul_i64(&dW->c[i - 1], &P[dim].c[i], (int64_t)i);
		for (c = 0; c < dim && st == QAWS_STATUS_OK; c++)
		{
			poly* dX = &P[dim + 1 + c];
			poly* Nc = &P[2 * dim + 2 + c];
			poly t;
			qaws_exact_poly_zero(dX, n - 1);
			for (i = 1; i <= n && st == QAWS_STATUS_OK; i++)
				st = qaws_exact_int_mul_i64(&dX->c[i - 1], &P[c].c[i], (int64_t)i);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_mul(Nc, dX, &P[dim]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_mul(&t, &P[c], dW);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(Nc, &t, -1);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_pad(Nc, 2 * n - 1);
		}
	}
	/* candidate directions: each coordinate (d == c), then sums and differences of two */
	for (c = 0; c < dim && st == QAWS_STATUS_OK && !*out; c++)
		for (d = c; d < dim && st == QAWS_STATUS_OK && !*out; d++)
		{
			int sgn;
			for (sgn = (d == c ? 1 : -1); sgn <= 1 && st == QAWS_STATUS_OK && !*out; sgn += 2)
			{
				poly* comb = &P[3 * dim + 2];
				int neg = 0, pos = 0;
				*comb = P[2 * dim + 2 + c];
				if (d != c)
					st = qaws_exact_poly_acc(comb, &P[2 * dim + 2 + d], sgn);
				if (st == QAWS_STATUS_OK) st = qaws_exact_poly_to_bernstein(comb, b);
				for (i = 0; i <= comb->deg && st == QAWS_STATUS_OK; i++)
				{
					int s = qaws_exact_int_sign(&b[i]);
					if (s < 0) neg = 1;
					if (s > 0) pos = 1;
				}
				/* one strict sign, the rest zero: a non-zero polynomial of one sign, the coordinate strictly monotone */
				if (neg != pos)
					*out = 1;
			}
		}
	cc_free(P);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Span against span                                                  */
/* ------------------------------------------------------------------ */

typedef struct span_hit
{
	double a_lo, a_hi, b_lo, b_hi;   /* parameters on the first and second span */
	int exact;
	int a_exact_end, b_exact_start;  /* a exactly at its span's end / b exactly at its start, a exactly at start, b at end */
	int a_exact_start, b_exact_end;
} span_hit;

static qaws_status emit_cross(qaws_exact_span const* I, int ishift, qaws_exact_span const* S, int sshift, qaws_exact_local_hit const* lh, int i_is_a,
	span_hit* out, unsigned int capacity, unsigned int* count)
{
	double s_lo, s_hi, r_lo, r_hi;
	int s_ex, r_ex;
	span_hit h;
	qaws_status st;
	TRY(map_s(I, ishift, lh, &s_lo, &s_hi, &s_ex));
	TRY(map_r(S, sshift, &lh->r, &r_lo, &r_hi, &r_ex));
	h.exact = s_ex && r_ex && lh->s_exact && lh->r.exact;
	if (i_is_a)
	{
		h.a_lo = s_lo; h.a_hi = s_hi; h.b_lo = r_lo; h.b_hi = r_hi;
		h.a_exact_start = s_is(lh, 0);
		h.a_exact_end = s_is(lh, 1);
		h.b_exact_start = lh->r.exact && lh->r.index == 0;
		h.b_exact_end = lh->r.exact && lh->r.index == ((uint64_t)1 << lh->r.depth);
	}
	else
	{
		h.a_lo = r_lo; h.a_hi = r_hi; h.b_lo = s_lo; h.b_hi = s_hi;
		h.a_exact_start = lh->r.exact && lh->r.index == 0;
		h.a_exact_end = lh->r.exact && lh->r.index == ((uint64_t)1 << lh->r.depth);
		h.b_exact_start = s_is(lh, 0);
		h.b_exact_end = s_is(lh, 1);
	}
	if (*count >= capacity)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	out[(*count)++] = h;
	return QAWS_STATUS_OK;
}

/* Implicitize I, substitute S, over the projections until one is not degenerate. */
static qaws_status cross_spans(qaws_exact_span const* I, int ishift, qaws_exact_span const* S, int sshift, unsigned int dim, int i_is_a, span_hit* out,
	unsigned int capacity, unsigned int* count, int* common)
{
	proj pj[3];
	unsigned int np = projections(dim, pj), p, nl = 0, k;
	poly* sys;   /* R C00 C01 H */
	qaws_exact_local_hit* lh;
	qaws_status st = QAWS_STATUS_CERTIFICATION_FAILED;
	*common = 1;
	if (I->degree < 2 || I->degree > CC_MAX_IMPLICIT || I->degree * S->degree > QAWS_EXACT_ROOTS_MAX_DEGREE)
	{
		*common = 0;
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	}
	sys = (poly*)cc_alloc(sizeof(poly) * 4);
	lh = (qaws_exact_local_hit*)cc_alloc(sizeof(qaws_exact_local_hit) * CC_MAX_ROOTS);
	if (!sys || !lh)
	{
		cc_free(sys);
		cc_free(lh);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (p = 0; p < np; p++)
	{
		qaws_exact_constraint con;
		int com = 0;
		st = build_cross(I, S, dim + 1, &pj[p], &sys[0], &sys[1], &sys[2], &sys[3]);
		if (st != QAWS_STATUS_OK)
		{
			*common = 0;
			break;
		}
		con.kind = QAWS_EXACT_CON_ZERO;
		con.H = &sys[3];
		st = qaws_exact_solve_system(&sys[0], &sys[1], &sys[2], &con, pj[p].cz >= 0 ? 1u : 0u, lh, CC_MAX_ROOTS, &nl, &com);
		if (com)
			continue;   /* this projection is degenerate (or shared): the next one */
		*common = 0;
		for (k = 0; k < nl && st == QAWS_STATUS_OK; k++)
			st = emit_cross(I, ishift, S, sshift, &lh[k], i_is_a, out, capacity, count);
		break;
	}
	cc_free(sys);
	cc_free(lh);
	return st;
}

/* Two segments: the exact linear roots in a projection where they are not parallel, the remaining coordinate checked exactly. */
static qaws_status cross_segments(qaws_exact_span const* A, int ashift, qaws_exact_span const* B, int bshift, unsigned int dim, span_hit* out,
	unsigned int capacity, unsigned int* count, int* common)
{
	proj pj[3];
	unsigned int np = projections(dim, pj), p, D = dim + 1, e, c;
	qaws_status st;
	*common = 1;
	for (p = 0; p < np; p++)
	{
		unsigned int idx[3];
		qaws_exact_int fa[2], fb[2], den, t1, t2;
		qaws_exact_local_hit ha, hb;
		idx[0] = pj[p].c0;
		idx[1] = pj[p].c1;
		idx[2] = pj[p].w;
		/* f(Q) = det(P0, P1, Q) in the projection: zero on the other segment's line */
		for (e = 0; e < 2; e++)
		{
			unsigned int side;
			for (side = 0; side < 2; side++)
			{
				qaws_exact_span const* L = side == 0 ? A : B;
				qaws_exact_span const* T = side == 0 ? B : A;
				qaws_exact_int* f = side == 0 ? &fb[e] : &fa[e];
				qaws_exact_int_zero(f);
				for (c = 0; c < 3; c++)
				{
					unsigned int c1 = idx[(c + 1) % 3], c2 = idx[(c + 2) % 3];
					TRY(qaws_exact_int_mul(&t1, &L->h[c1], &L->h[D + c2]));
					TRY(qaws_exact_int_mul(&t2, &L->h[c2], &L->h[D + c1]));
					TRY(qaws_exact_int_sub(&t1, &t1, &t2));
					TRY(qaws_exact_int_mul(&t1, &t1, &T->h[e * D + idx[c]]));
					TRY(qaws_exact_int_add(f, f, &t1));
				}
			}
		}
		if (qaws_exact_int_is_zero(&fb[0]) && qaws_exact_int_is_zero(&fb[1]))
			continue;   /* collinear in this projection */
		*common = 0;
		if (qaws_exact_int_sign(&fb[0]) * qaws_exact_int_sign(&fb[1]) > 0 || qaws_exact_int_sign(&fa[0]) * qaws_exact_int_sign(&fa[1]) > 0)
			return QAWS_STATUS_OK;
		/* r on B: f0 / (f0 - f1); s on A likewise */
		memset(&ha, 0, sizeof(ha));
		memset(&hb, 0, sizeof(hb));
		hb.s_exact = 1;
		hb.s_num = fb[0];
		TRY(qaws_exact_int_sub(&den, &fb[0], &fb[1]));
		hb.s_den = den;
		ha.s_exact = 1;
		ha.s_num = fa[0];
		TRY(qaws_exact_int_sub(&den, &fa[0], &fa[1]));
		ha.s_den = den;
		if (qaws_exact_int_sign(&ha.s_den) < 0) { qaws_exact_int_neg(&ha.s_num, &ha.s_num); qaws_exact_int_neg(&ha.s_den, &ha.s_den); }
		if (qaws_exact_int_sign(&hb.s_den) < 0) { qaws_exact_int_neg(&hb.s_num, &hb.s_num); qaws_exact_int_neg(&hb.s_den, &hb.s_den); }
		if (pj[p].cz >= 0)
		{
			/* the remaining coordinate: Z_A(s) W_B(r) == Z_B(r) W_A(s) with homogeneous s = sn : (sd - sn) */
			qaws_exact_int za, wa, zb, wb, l, r, u;
			unsigned int z = (unsigned int)pj[p].cz, w = pj[p].w;
			TRY(qaws_exact_int_sub(&u, &ha.s_den, &ha.s_num));
			TRY(qaws_exact_int_mul(&za, &A->h[z], &u));
			TRY(qaws_exact_int_mul(&t1, &A->h[D + z], &ha.s_num));
			TRY(qaws_exact_int_add(&za, &za, &t1));
			TRY(qaws_exact_int_mul(&wa, &A->h[w], &u));
			TRY(qaws_exact_int_mul(&t1, &A->h[D + w], &ha.s_num));
			TRY(qaws_exact_int_add(&wa, &wa, &t1));
			TRY(qaws_exact_int_sub(&u, &hb.s_den, &hb.s_num));
			TRY(qaws_exact_int_mul(&zb, &B->h[z], &u));
			TRY(qaws_exact_int_mul(&t1, &B->h[D + z], &hb.s_num));
			TRY(qaws_exact_int_add(&zb, &zb, &t1));
			TRY(qaws_exact_int_mul(&wb, &B->h[w], &u));
			TRY(qaws_exact_int_mul(&t1, &B->h[D + w], &hb.s_num));
			TRY(qaws_exact_int_add(&wb, &wb, &t1));
			TRY(qaws_exact_int_mul(&l, &za, &wb));
			TRY(qaws_exact_int_mul(&r, &zb, &wa));
			if (qaws_exact_int_cmp(&l, &r) != 0)
				return QAWS_STATUS_OK;
		}
		{
			span_hit h;
			int ea = 0, eb = 0;
			TRY(map_s(A, ashift, &ha, &h.a_lo, &h.a_hi, &ea));
			TRY(map_s(B, bshift, &hb, &h.b_lo, &h.b_hi, &eb));
			h.exact = ea && eb;
			h.a_exact_start = s_is(&ha, 0);
			h.a_exact_end = s_is(&ha, 1);
			h.b_exact_start = s_is(&hb, 0);
			h.b_exact_end = s_is(&hb, 1);
			if (*count >= capacity)
				return QAWS_STATUS_BUFFER_TOO_SMALL;
			out[(*count)++] = h;
		}
		return QAWS_STATUS_OK;
	}
	return QAWS_STATUS_CERTIFICATION_FAILED;   /* collinear in every projection */
}

/* Span A against span B (any degrees): implicitize the lower degree, the other on failure. */
static qaws_status span_pair(qaws_exact_span const* A, int ashift, qaws_exact_span const* B, int bshift, unsigned int dim, span_hit* out,
	unsigned int capacity, unsigned int* count, int* common)
{
	unsigned int na = A->degree, nb = B->degree, start = *count;
	int a_first, com2 = 0;
	qaws_status st;
	if (na == 1 && nb == 1)
		return cross_segments(A, ashift, B, bshift, dim, out, capacity, count, common);
	a_first = na >= 2 && (na <= nb || nb < 2);
	st = a_first ? cross_spans(A, ashift, B, bshift, dim, 1, out, capacity, count, common)
	             : cross_spans(B, bshift, A, ashift, dim, 0, out, capacity, count, common);
	if ((st == QAWS_STATUS_CERTIFICATION_FAILED || st == QAWS_STATUS_EXACT_UNSUPPORTED) && (a_first ? nb : na) >= 2)
	{
		qaws_status st2;
		*count = start;
		st2 = a_first ? cross_spans(B, bshift, A, ashift, dim, 0, out, capacity, count, &com2)
		              : cross_spans(A, ashift, B, bshift, dim, 1, out, capacity, count, &com2);
		if (st2 == QAWS_STATUS_OK)
		{
			*common = 0;
			return st2;
		}
		*common = *common && com2;
		*count = start;
	}
	return st;
}

/* ------------------------------------------------------------------ */
/*  Point inversion (for spans on one algebraic curve)                 */
/* ------------------------------------------------------------------ */

/* Homogeneous point of a span at local s = 0, 1 or 1/2 (times a positive factor). */
static qaws_status span_point(qaws_exact_span const* sp, unsigned int D, int which, qaws_exact_int* P)
{
	unsigned int n = sp->degree, i, c;
	qaws_exact_int t;
	qaws_status st;
	for (c = 0; c < D; c++)
	{
		if (which == 0)
			P[c] = sp->h[c];
		else if (which == 1)
			P[c] = sp->h[n * D + c];
		else
		{
			qaws_exact_int_zero(&P[c]);
			for (i = 0; i <= n; i++)
			{
				TRY(qaws_exact_int_mul_i64(&t, &sp->h[i * D + c], qaws_exact_binom64(n, i)));
				TRY(qaws_exact_int_add(&P[c], &P[c], &t));
			}
		}
	}
	return QAWS_STATUS_OK;
}

/*
 * The local parameter of the point E on span I's curve: s = num / den.
 * *on = 0 when E is not on the curve. Segments: E = l I0 + m I1, s = m / (l + m);
 * higher degrees: the cofactors of M(E). QAWS_STATUS_CERTIFICATION_FAILED at
 * a singular point.
 */
static qaws_status invert_point(qaws_exact_span const* I, unsigned int dim, qaws_exact_int const* E, qaws_exact_int* num, qaws_exact_int* den, int* on)
{
	unsigned int D = dim + 1, n = I->degree, i, c;
	proj pj[3];
	unsigned int np = projections(dim, pj), p;
	qaws_status st = QAWS_STATUS_OK;
	*on = 0;
	if (n == 1)
	{
		/* in a (c, w) plane where I0, I1 differ: l = det(E, I1), m = det(I0, E) */
		for (c = 0; c < dim; c++)
		{
			qaws_exact_int d, l, m, t1, t2;
			TRY(qaws_exact_int_mul(&t1, &I->h[c], &I->h[D + dim]));
			TRY(qaws_exact_int_mul(&t2, &I->h[dim], &I->h[D + c]));
			TRY(qaws_exact_int_sub(&d, &t1, &t2));
			if (qaws_exact_int_is_zero(&d))
				continue;
			TRY(qaws_exact_int_mul(&t1, &E[c], &I->h[D + dim]));
			TRY(qaws_exact_int_mul(&t2, &E[dim], &I->h[D + c]));
			TRY(qaws_exact_int_sub(&l, &t1, &t2));
			TRY(qaws_exact_int_mul(&t1, &I->h[c], &E[dim]));
			TRY(qaws_exact_int_mul(&t2, &I->h[dim], &E[c]));
			TRY(qaws_exact_int_sub(&m, &t1, &t2));
			/* every component: d E == l I0 + m I1 */
			for (i = 0; i < D; i++)
			{
				qaws_exact_int lhs, rhs;
				TRY(qaws_exact_int_mul(&lhs, &d, &E[i]));
				TRY(qaws_exact_int_mul(&rhs, &l, &I->h[i]));
				TRY(qaws_exact_int_mul(&t1, &m, &I->h[D + i]));
				TRY(qaws_exact_int_add(&rhs, &rhs, &t1));
				if (qaws_exact_int_cmp(&lhs, &rhs) != 0)
					return QAWS_STATUS_OK;
			}
			TRY(qaws_exact_int_add(den, &l, &m));
			*num = m;
			if (qaws_exact_int_is_zero(den))
				return QAWS_STATUS_OK;   /* the point at infinity of the line */
			if (qaws_exact_int_sign(den) < 0)
			{
				qaws_exact_int_neg(num, num);
				qaws_exact_int_neg(den, den);
			}
			*on = 1;
			return QAWS_STATUS_OK;
		}
		return QAWS_STATUS_CERTIFICATION_FAILED;
	}
	for (p = 0; p < np; p++)
	{
		poly* M;
		qaws_exact_int* B;
		poly* sys;
		M = (poly*)cc_alloc(sizeof(poly) * (n * n + 3 + 3));
		B = (qaws_exact_int*)cc_alloc(sizeof(qaws_exact_int) * 3 * n * n);
		if (!M || !B)
		{
			cc_free(M);
			cc_free(B);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		sys = &M[n * n];
		st = qaws_exact_poly_from_bernstein(I->h, n, D, pj[p].c0, &sys[3]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(I->h, n, D, pj[p].c1, &sys[4]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_from_bernstein(I->h, n, D, pj[p].w, &sys[5]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&sys[4], &sys[5], n, &B[0]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&sys[5], &sys[3], n, &B[n * n]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_bezout(&sys[3], &sys[4], n, &B[2 * n * n]);
		for (i = 0; i < n * n && st == QAWS_STATUS_OK; i++)
		{
			qaws_exact_int t;
			qaws_exact_poly_zero(&M[i], 0);
			st = qaws_exact_int_mul(&t, &E[pj[p].c0], &B[i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&M[i].c[0], &M[i].c[0], &t);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &E[pj[p].c1], &B[n * n + i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&M[i].c[0], &M[i].c[0], &t);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &E[pj[p].w], &B[2 * n * n + i]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&M[i].c[0], &M[i].c[0], &t);
		}
		if (st == QAWS_STATUS_OK)
			st = qaws_exact_det_cofactors(M, n, &sys[0], &sys[1], &sys[2]);
		if (st == QAWS_STATUS_OK)
		{
			int degenerate = 1;
			for (i = 0; i < n * n && degenerate; i++)
				if (!qaws_exact_int_is_zero(&B[i]) || !qaws_exact_int_is_zero(&B[n * n + i]) || !qaws_exact_int_is_zero(&B[2 * n * n + i]))
					degenerate = 0;
			if (!qaws_exact_int_is_zero(&sys[0].c[0]))
			{
				/* not on the curve's projection */
				cc_free(M);
				cc_free(B);
				return QAWS_STATUS_OK;
			}
			if (qaws_exact_int_is_zero(&sys[1].c[0]))
			{
				/* kernel (0, ..., 0, 1): the parameter's point at infinity, on the curve but off
				   the span. M's last row and column vanish, adj(M) is a multiple of e e^T and its
				   one non-zero entry is the leading (n - 1) x (n - 1) minor. */
				int at_inf = 0;
				if (!degenerate && n >= 2)
				{
					poly* sub = (poly*)cc_alloc(sizeof(poly) * ((n - 1) * (n - 1) + 1));
					unsigned int r2, c2;
					if (sub)
					{
						for (r2 = 0; r2 + 1 < n; r2++)
							for (c2 = 0; c2 + 1 < n; c2++)
								sub[r2 * (n - 1) + c2] = M[r2 * n + c2];
						if (qaws_exact_poly_det(sub, n - 1, 0, 0, &sub[(n - 1) * (n - 1)]) == QAWS_STATUS_OK &&
						    !qaws_exact_int_is_zero(&sub[(n - 1) * (n - 1)].c[0]))
						{
							/* and the last column of M is zero: e_(n-1) is the kernel */
							at_inf = 1;
							for (r2 = 0; r2 < n; r2++)
								if (!qaws_exact_int_is_zero(&M[r2 * n + n - 1].c[0]))
									at_inf = 0;
						}
						cc_free(sub);
					}
				}
				cc_free(M);
				cc_free(B);
				if (at_inf)
				{
					qaws_exact_int_from_i64(num, 1);
					qaws_exact_int_zero(den);
					*on = 1;
					return QAWS_STATUS_OK;
				}
				if (degenerate && p + 1 < np)
					continue;
				return QAWS_STATUS_CERTIFICATION_FAILED;   /* singular point, or a degenerate span */
			}
			*num = sys[2].c[0];
			*den = sys[1].c[0];
			if (qaws_exact_int_sign(den) < 0)
			{
				qaws_exact_int_neg(num, num);
				qaws_exact_int_neg(den, den);
			}
			/* 3D: the remaining coordinate at s, exactly: Z(s) W_E == Z_E W(s) */
			if (pj[p].cz >= 0)
			{
				qaws_exact_int z, w, pw, t, l, r;
				poly Zp;
				unsigned int k;
				st = qaws_exact_poly_from_bernstein(I->h, n, D, (unsigned int)pj[p].cz, &Zp);
				qaws_exact_int_zero(&z);
				qaws_exact_int_zero(&w);
				for (k = 0; k <= n && st == QAWS_STATUS_OK; k++)
				{
					/* num^k den^(n - k) */
					unsigned int j;
					qaws_exact_int_from_i64(&pw, 1);
					for (j = 0; j < k && st == QAWS_STATUS_OK; j++) st = qaws_exact_int_mul(&pw, &pw, num);
					for (j = k; j < n && st == QAWS_STATUS_OK; j++) st = qaws_exact_int_mul(&pw, &pw, den);
					if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &Zp.c[k], &pw);
					if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&z, &z, &t);
					if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &sys[5].c[k], &pw);
					if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&w, &w, &t);
				}
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&l, &z, &E[pj[p].w]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&r, &E[pj[p].cz], &w);
				if (st == QAWS_STATUS_OK && qaws_exact_int_cmp(&l, &r) != 0)
				{
					cc_free(M);
					cc_free(B);
					return QAWS_STATUS_OK;
				}
			}
			*on = st == QAWS_STATUS_OK;
		}
		cc_free(M);
		cc_free(B);
		return st;
	}
	return QAWS_STATUS_CERTIFICATION_FAILED;
}

/* where in [0, 1]: -1 outside, 0 at 0, 1 inside, 2 at 1 */
static int classify01(qaws_exact_int const* num, qaws_exact_int const* den)
{
	int s = qaws_exact_int_sign(num), c = qaws_exact_int_cmp(num, den);
	if (qaws_exact_int_is_zero(den))
		return -1;   /* the parameter at infinity */
	if (s < 0 || c > 0)
		return -1;
	if (s == 0)
		return 0;
	return c == 0 ? 2 : 1;
}

/*
 * Two spans on one algebraic curve: they overlap, or touch at endpoints.
 * Each endpoint of one is inverted on the other; an inversion strictly
 * inside is an overlap (QAWS_STATUS_CERTIFICATION_FAILED), endpoint-to-
 * endpoint contacts are hits. When every endpoint lands on an endpoint the
 * arcs coincide or complete each other: B's midpoint decides.
 */
static qaws_status common_spans(qaws_exact_span const* A, int ashift, qaws_exact_span const* B, int bshift, unsigned int dim, span_hit* out,
	unsigned int capacity, unsigned int* count)
{
	unsigned int D = dim + 1, e;
	qaws_exact_int P[4], num, den;
	int on, cls, all_ends = 1;
	qaws_status st;
	/* endpoints of B on A */
	for (e = 0; e < 2; e++)
	{
		TRY(span_point(B, D, (int)e, P));
		st = invert_point(A, dim, P, &num, &den, &on);
		if (st != QAWS_STATUS_OK) return st;
		if (!on)
		{
			all_ends = 0;
			continue;
		}
		cls = classify01(&num, &den);
		if (cls == 1)
			return QAWS_STATUS_CERTIFICATION_FAILED;   /* overlap */
		if (cls < 0)
		{
			all_ends = 0;
			continue;
		}
		{
			span_hit h;
			memset(&h, 0, sizeof(h));
			h.a_lo = h.a_hi = ldexp((double)(cls == 0 ? A->a : A->b), -ashift);
			h.b_lo = h.b_hi = ldexp((double)(e == 0 ? B->a : B->b), -bshift);
			h.exact = 1;
			h.a_exact_start = cls == 0;
			h.a_exact_end = cls == 2;
			h.b_exact_start = e == 0;
			h.b_exact_end = e == 1;
			if (*count >= capacity)
				return QAWS_STATUS_BUFFER_TOO_SMALL;
			out[(*count)++] = h;
		}
	}
	/* endpoints of A strictly inside B: overlap */
	for (e = 0; e < 2; e++)
	{
		TRY(span_point(A, D, (int)e, P));
		st = invert_point(B, dim, P, &num, &den, &on);
		if (st != QAWS_STATUS_OK) return st;
		if (on && classify01(&num, &den) == 1)
			return QAWS_STATUS_CERTIFICATION_FAILED;
	}
	if (all_ends)
	{
		TRY(span_point(B, D, 2, P));
		st = invert_point(A, dim, P, &num, &den, &on);
		if (st != QAWS_STATUS_OK) return st;
		if (on && classify01(&num, &den) == 1)
			return QAWS_STATUS_CERTIFICATION_FAILED;   /* the same arc */
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Public                                                             */
/* ------------------------------------------------------------------ */

static void sort_pairs(qaws_exact_pair* p, unsigned int n)
{
	unsigned int i;
	for (i = 1; i < n; i++)
	{
		qaws_exact_pair key = p[i];
		int j = (int)i - 1;
		while (j >= 0 && (p[j].a_lo > key.a_lo || (p[j].a_lo == key.a_lo && p[j].b_lo > key.b_lo)))
		{
			p[j + 1] = p[j];
			j--;
		}
		p[j + 1] = key;
	}
}

static qaws_status push_pair(qaws_exact_pair* out, unsigned int capacity, unsigned int* count, span_hit const* h, int swap)
{
	qaws_exact_pair p;
	if (*count >= capacity)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	p.kind = h->exact ? QAWS_EXACT_HIT_POINT : QAWS_EXACT_HIT_CROSSING;
	p.a_lo = swap ? h->b_lo : h->a_lo;
	p.a_hi = swap ? h->b_hi : h->a_hi;
	p.b_lo = swap ? h->a_lo : h->b_lo;
	p.b_hi = swap ? h->a_hi : h->b_hi;
	out[(*count)++] = p;
	return QAWS_STATUS_OK;
}

static int curve_closed(qaws_exact_curve const* c);

void qaws_exact_span_box(qaws_exact_span const* sp, unsigned int dim, double lo[3], double hi[3])
{
	unsigned int D = dim + 1, j, k;
	for (k = 0; k < 3; k++)
	{
		lo[k] = k < dim ? HUGE_VAL : 0;
		hi[k] = k < dim ? -HUGE_VAL : 0;
	}
	for (j = 0; j <= sp->degree; j++)
		if (qaws_exact_int_sign(&sp->h[j * D + dim]) <= 0)
		{
			for (k = 0; k < dim; k++)
			{
				lo[k] = HUGE_VAL;
				hi[k] = -HUGE_VAL;
			}
			return;
		}
	for (j = 0; j <= sp->degree; j++)
		for (k = 0; k < dim; k++)
		{
			/* correctly rounded (twice at worst, for subnormals): two ulps out are sound */
			double v = qaws_exact_ratio_to_double(&sp->h[j * D + k], &sp->h[j * D + dim]);
			double l = nextafter(nextafter(v, -HUGE_VAL), -HUGE_VAL), h = nextafter(nextafter(v, HUGE_VAL), HUGE_VAL);
			if (l < lo[k]) lo[k] = l;
			if (h > hi[k]) hi[k] = h;
		}
}

static int box_disjoint(double const* alo, double const* ahi, double const* blo, double const* bhi, unsigned int dim)
{
	unsigned int k;
	for (k = 0; k < dim; k++)
		if (alo[k] <= ahi[k] && blo[k] <= bhi[k] && (ahi[k] < blo[k] || bhi[k] < alo[k]))
			return 1;
	return 0;
}

qaws_status qaws_exact_span_pair_hits(qaws_exact_curve const* a, unsigned int ia, qaws_exact_curve const* b, unsigned int ib,
	qaws_exact_pair* out, unsigned int capacity, unsigned int* count, int* common)
{
	unsigned int nh = 0, k, dim = (unsigned int)a->dimension;
	span_hit* hits = (span_hit*)cc_alloc(sizeof(span_hit) * CC_MAX_ROOTS);
	qaws_status st;
	if (!hits)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	*common = 0;
	st = span_pair(&a->spans[ia], a->param_shift, &b->spans[ib], b->param_shift, dim, hits, CC_MAX_ROOTS, &nh, common);
	for (k = 0; k < nh && st == QAWS_STATUS_OK; k++)
	{
		/* a knot point is reported by the span it ends, the start of a
		   closed curve by its last span */
		if ((hits[k].b_exact_start && (ib > 0 || curve_closed(b))) || (hits[k].a_exact_start && (ia > 0 || curve_closed(a))))
			continue;
		st = push_pair(out, capacity, count, &hits[k], 0);
	}
	cc_free(hits);
	return st;
}

qaws_status qaws_exact_curve_curve_hits(qaws_exact_curve const* a, qaws_exact_curve const* b, qaws_exact_pair* out_pairs, unsigned int capacity,
	unsigned int* out_count)
{
	unsigned int ia, ib, count = 0, dim;
	double (*box)[2][3];
	qaws_status st = QAWS_STATUS_OK;
	if (!a || !b || !out_count || (!out_pairs && capacity) || a->dimension != b->dimension || (a->dimension != 2 && a->dimension != 3))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (a->space_exp2 != b->space_exp2)
		return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
	dim = (unsigned int)a->dimension;
	/* control boxes of b's spans: span pairs with disjoint boxes cannot meet */
	box = (double (*)[2][3])cc_alloc(sizeof(double) * 6 * (b->span_count + 1));
	if (!box)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (ib = 0; ib < b->span_count; ib++)
		qaws_exact_span_box(&b->spans[ib], dim, box[ib][0], box[ib][1]);
	for (ia = 0; ia < a->span_count && st == QAWS_STATUS_OK; ia++)
	{
		double alo[3], ahi[3];
		qaws_exact_span_box(&a->spans[ia], dim, alo, ahi);
		for (ib = 0; ib < b->span_count && st == QAWS_STATUS_OK; ib++)
		{
			int common = 0;
			if (box_disjoint(alo, ahi, box[ib][0], box[ib][1], dim))
				continue;
			st = qaws_exact_span_pair_hits(a, ia, b, ib, out_pairs, capacity, &count, &common);
			if (st == QAWS_STATUS_CERTIFICATION_FAILED && common)
				st = QAWS_STATUS_CERTIFICATION_FAILED;   /* a common component: the curves overlap or share their support */
		}
	}
	cc_free(box);
	sort_pairs(out_pairs, count);
	*out_count = count;
	return st;
}

/* Is the curve closed (its start and end the same point, exactly)? */
static int curve_closed(qaws_exact_curve const* c)
{
	unsigned int D = (unsigned int)c->dimension + 1, i;
	qaws_exact_span const* f = &c->spans[0];
	qaws_exact_span const* l = &c->spans[c->span_count - 1];
	qaws_exact_int x, y;
	for (i = 0; i < D - 1; i++)
	{
		/* f0_i l_w == l_i f0_w */
		if (qaws_exact_int_mul(&x, &f->h[i], &l->h[l->degree * D + D - 1]) != QAWS_STATUS_OK ||
		    qaws_exact_int_mul(&y, &l->h[l->degree * D + i], &f->h[D - 1]) != QAWS_STATUS_OK)
			return 0;
		if (qaws_exact_int_cmp(&x, &y) != 0)
			return 0;
	}
	return 1;
}

qaws_status qaws_exact_curve_self_hits(qaws_exact_curve const* curve, qaws_exact_pair* out_pairs, unsigned int capacity, unsigned int* out_count)
{
	unsigned int i, j, k, count = 0, dim;
	span_hit* hits;
	poly* sys;
	qaws_exact_local_hit* lh;
	int closed;
	double t_start, t_end;
	qaws_status st = QAWS_STATUS_OK;
	if (!curve || !out_count || (!out_pairs && capacity) || (curve->dimension != 2 && curve->dimension != 3))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	dim = (unsigned int)curve->dimension;
	closed = curve_closed(curve);
	t_start = ldexp((double)curve->spans[0].a, -curve->param_shift);
	t_end = ldexp((double)curve->spans[curve->span_count - 1].b, -curve->param_shift);
	hits = (span_hit*)cc_alloc(sizeof(span_hit) * CC_MAX_ROOTS);
	sys = (poly*)cc_alloc(sizeof(poly) * 5);
	lh = (qaws_exact_local_hit*)cc_alloc(sizeof(qaws_exact_local_hit) * CC_MAX_ROOTS);
	if (!hits || !sys || !lh)
	{
		cc_free(hits);
		cc_free(sys);
		cc_free(lh);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < curve->span_count && st == QAWS_STATUS_OK; i++)
	{
		qaws_exact_span const* sp = &curve->spans[i];
		int mono = 0;
		/* inside the span */
		st = span_monotone(sp, dim, &mono);
		if (st == QAWS_STATUS_OK && !mono)
		{
			if (sp->degree == 2)
			{
				/* a conic arc is injective unless it is a line folding back */
				int on_line = 0;
				poly R0;
				qaws_exact_poly_zero(&R0, 0);
				{
					qaws_exact_int P[4], num, den;
					int on;
					qaws_exact_span seg;
					qaws_exact_int segh[8];
					unsigned int c, D = dim + 1;
					for (c = 0; c < D; c++)
					{
						segh[c] = sp->h[c];
						segh[D + c] = sp->h[2 * D + c];
					}
					seg.degree = 1;
					seg.a = 0;
					seg.b = 1;
					seg.h = segh;
					for (c = 0; c < D; c++)
						P[c] = sp->h[D + c];
					st = invert_point(&seg, dim, P, &num, &den, &on);
					/* a failed inversion: equal ends, the arc goes out and back */
					on_line = st == QAWS_STATUS_CERTIFICATION_FAILED || (st == QAWS_STATUS_OK && on);
					if (st == QAWS_STATUS_CERTIFICATION_FAILED)
						st = QAWS_STATUS_OK;
				}
				if (on_line)
					st = QAWS_STATUS_CERTIFICATION_FAILED;   /* a degenerate arc overlapping itself */
			}
			else if (sp->degree >= 3)
			{
				proj pj[3];
				unsigned int np = projections(dim, pj), p, nl = 0;
				int com = 1;
				if (sp->degree - 1 > CC_MAX_IMPLICIT)
					st = QAWS_STATUS_EXACT_UNSUPPORTED;
				for (p = 0; p < np && com && st == QAWS_STATUS_OK; p++)
				{
					qaws_exact_constraint con[2];
					st = build_self(sp, dim + 1, &pj[p], &sys[0], &sys[1], &sys[2], &sys[3], &sys[4]);
					if (st != QAWS_STATUS_OK)
						break;
					con[0].kind = QAWS_EXACT_CON_SIGN_OF_C00;
					con[0].H = &sys[3];
					con[1].kind = QAWS_EXACT_CON_ZERO;
					con[1].H = &sys[4];
					st = qaws_exact_solve_system(&sys[0], &sys[1], &sys[2], con, pj[p].cz >= 0 ? 2u : 1u, lh, CC_MAX_ROOTS, &nl, &com);
					if (com)
						st = QAWS_STATUS_OK;
				}
				if (st == QAWS_STATUS_OK && com)
					st = QAWS_STATUS_CERTIFICATION_FAILED;
				for (k = 0; k < nl && st == QAWS_STATUS_OK; k++)
				{
					/* r = t (the smaller parameter), s the larger */
					span_hit h;
					int ea = 0, eb = 0;
					st = map_r(sp, curve->param_shift, &lh[k].r, &h.a_lo, &h.a_hi, &ea);
					if (st == QAWS_STATUS_OK) st = map_s(sp, curve->param_shift, &lh[k], &h.b_lo, &h.b_hi, &eb);
					h.exact = ea && eb && lh[k].r.exact && lh[k].s_exact;
					if (st == QAWS_STATUS_OK)
					{
						if (closed && h.exact && h.a_lo == t_start && h.b_lo == t_end)
							continue;
						if (i > 0 && lh[k].r.exact && lh[k].r.index == 0)
							continue;   /* a knot point: reported by the previous span */
						st = push_pair(out_pairs, capacity, &count, &h, 0);
					}
				}
			}
		}
		/* against the later spans */
		for (j = i + 1; j < curve->span_count && st == QAWS_STATUS_OK; j++)
		{
			unsigned int nh = 0;
			int common = 0;
			st = span_pair(sp, curve->param_shift, &curve->spans[j], curve->param_shift, dim, hits, CC_MAX_ROOTS, &nh, &common);
			if (st == QAWS_STATUS_CERTIFICATION_FAILED && common)
			{
				nh = 0;
				st = common_spans(sp, curve->param_shift, &curve->spans[j], curve->param_shift, dim, hits, CC_MAX_ROOTS, &nh);
			}
			for (k = 0; k < nh && st == QAWS_STATUS_OK; k++)
			{
				span_hit const* h = &hits[k];
				if (h->exact && h->a_lo == h->b_lo)
					continue;   /* the shared knot point */
				if (closed && h->exact && h->a_lo == t_start && h->b_lo == t_end)
					continue;   /* the closing point */
				if ((h->a_exact_start && i > 0) || h->b_exact_start)
					continue;   /* a knot point: reported by the span ending there */
				st = push_pair(out_pairs, capacity, &count, h, 0);
			}
		}
	}
	cc_free(hits);
	cc_free(sys);
	cc_free(lh);
	sort_pairs(out_pairs, count);
	*out_count = count;
	return st;
}
