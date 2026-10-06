#include "qaws_exact_surface.h"
#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

typedef qaws_exact_poly poly;

static void* sh_alloc(size_t bytes)
{
	return qaws_internal_alloc(NULL, (unsigned long)bytes);
}

static void sh_free(void* p)
{
	qaws_internal_dealloc(NULL, p);
}

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

static qaws_status dy_mul(dyadic* r, dyadic const* a, dyadic const* b)
{
	r->e = a->e + b->e;
	return qaws_exact_int_mul(&r->m, &a->m, &b->m);
}

/* ------------------------------------------------------------------ */
/*  Patch polynomials                                                  */
/* ------------------------------------------------------------------ */

/*
 * The plane N . (x - P) over a patch: W (N . (S - P)) has the integer
 * Bernstein coefficients F_ab = sum_c N_c (H_abc 2^E - P_c W_ab), brought to
 * one exponent.
 */
static qaws_status plane_net(qaws_exact_surface const* s, qaws_exact_int const* h, dyadic const* N, dyadic const* P, qaws_exact_int* F)
{
	unsigned int n = (s->p + 1) * (s->q + 1), i, c;
	dyadic* g;
	int emin = 0, any = 0;
	qaws_status st = QAWS_STATUS_OK;
	g = (dyadic*)sh_alloc(sizeof(dyadic) * n);
	if (!g)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
	{
		qaws_exact_int_zero(&g[i].m);
		g[i].e = 0;
		for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
		{
			dyadic hx, w, t1, t2;
			hx.m = h[i * 4 + c];
			hx.e = s->space_exp2;
			w.m = h[i * 4 + 3];
			w.e = 0;
			st = dy_mul(&t1, &N[c], &hx);
			if (st == QAWS_STATUS_OK) st = dy_mul(&t2, &N[c], &P[c]);
			if (st == QAWS_STATUS_OK) st = dy_mul(&t2, &t2, &w);
			if (st == QAWS_STATUS_OK) st = dy_add(&t1, &t1, &t2, -1);
			if (st == QAWS_STATUS_OK) st = dy_add(&g[i], &g[i], &t1, 1);
		}
		if (st == QAWS_STATUS_OK && !qaws_exact_int_is_zero(&g[i].m))
		{
			if (!any || g[i].e < emin)
				emin = g[i].e;
			any = 1;
		}
	}
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
	{
		F[i] = g[i].m;
		if (!qaws_exact_int_is_zero(&F[i]) && g[i].e > emin)
			st = qaws_exact_int_shl(&F[i], &F[i], (unsigned int)(g[i].e - emin));
	}
	sh_free(g);
	return st;
}

/*
 * Power basis in the eliminated variable s with coefficients polynomials
 * (power basis) in the other one r: out[i] is the coefficient of s^i.
 * F is indexed [a (u) * (q + 1) + b (v)].
 */
static qaws_status to_power2(qaws_exact_int const* F, unsigned int p, unsigned int q, int elim_u, poly* out)
{
	unsigned int ns = elim_u ? p : q, nr = elim_u ? q : p, a, i, j;
	poly* rows;
	qaws_exact_int* line;
	qaws_status st = QAWS_STATUS_OK;
	rows = (poly*)sh_alloc(sizeof(poly) * (ns + 1));
	line = (qaws_exact_int*)sh_alloc(sizeof(qaws_exact_int) * (nr + 1));
	if (!rows || !line)
	{
		sh_free(rows);
		sh_free(line);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* along r first: rows[a] = the r-power polynomial of Bernstein row a (in s) */
	for (a = 0; a <= ns && st == QAWS_STATUS_OK; a++)
	{
		for (j = 0; j <= nr; j++)
			line[j] = elim_u ? F[a * (q + 1) + j] : F[j * (q + 1) + a];
		st = qaws_exact_poly_from_bernstein(line, nr, 1, 0, &rows[a]);
	}
	/* then along s: out_i = sum_a C(ns, a) C(ns - a, i - a) (-1)^(i - a) rows[a] */
	for (i = 0; i <= ns && st == QAWS_STATUS_OK; i++)
	{
		qaws_exact_poly_zero(&out[i], nr);
		for (a = 0; a <= i && st == QAWS_STATUS_OK; a++)
		{
			qaws_exact_int f;
			poly t;
			qaws_exact_int_from_i64(&f, qaws_exact_binom64(ns, a) * qaws_exact_binom64(ns - a, i - a) * (((i - a) & 1) ? -1 : 1));
			st = qaws_exact_poly_scale(&t, &rows[a], &f);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(&out[i], &t, 1);
		}
	}
	sh_free(rows);
	sh_free(line);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Exact sub-patches                                                  */
/* ------------------------------------------------------------------ */

/*
 * Split the control points Q[0..n] (stride, 4 components) at x = a / b,
 * integer De Casteljau with the weights (b - a, a): both halves come out
 * times b^n. keep_right selects the half kept in place.
 */
static qaws_status split_at(qaws_exact_int* Q, unsigned int n, unsigned int stride, qaws_exact_int const* a, qaws_exact_int const* b, int keep_right)
{
	qaws_exact_int d[17 * 4], out[17 * 4], ba, t1, t2;
	unsigned int r, i, c;
	qaws_status st;
	TRY(qaws_exact_int_sub(&ba, b, a));
	for (i = 0; i <= n; i++)
		for (c = 0; c < 4; c++)
			d[i * 4 + c] = Q[i * stride + c];
	/* level 0 */
	for (c = 0; c < 4; c++)
	{
		if (!keep_right)
			out[c] = d[c];
		else
			out[n * 4 + c] = d[n * 4 + c];
	}
	for (r = 1; r <= n; r++)
	{
		for (i = 0; i + r <= n; i++)
			for (c = 0; c < 4; c++)
			{
				TRY(qaws_exact_int_mul(&t1, &d[i * 4 + c], &ba));
				TRY(qaws_exact_int_mul(&t2, &d[(i + 1) * 4 + c], a));
				TRY(qaws_exact_int_add(&d[i * 4 + c], &t1, &t2));
			}
		for (c = 0; c < 4; c++)
		{
			if (!keep_right)
				out[r * 4 + c] = d[c];               /* left_r = d_0^r b^(n-r) */
			else
				out[(n - r) * 4 + c] = d[(n - r) * 4 + c];   /* right_(n-r) = d_(n-r)^r b^(n-r) */
		}
	}
	/* bring every point to the scale b^n */
	for (i = 0; i <= n; i++)
	{
		unsigned int lev = keep_right ? n - i : i, k;   /* the level that produced point i */
		for (c = 0; c < 4; c++)
		{
			t1 = out[i * 4 + c];
			for (k = lev; k < n; k++)
				TRY(qaws_exact_int_mul(&t1, &t1, b));
			Q[i * stride + c] = t1;
		}
	}
	return QAWS_STATUS_OK;
}

/* Restrict a curve of control points to [lo, hi] (rationals ln / ld, hn / hd, 0 <= lo <= hi <= 1). */
static qaws_status restrict_curve(qaws_exact_int* Q, unsigned int n, unsigned int stride, qaws_exact_int const* ln, qaws_exact_int const* ld,
	qaws_exact_int const* hn, qaws_exact_int const* hd)
{
	qaws_exact_int bn, bd, t1, t2;
	qaws_status st;
	if (!qaws_exact_int_is_zero(ln))
		TRY(split_at(Q, n, stride, ln, ld, 1));   /* keep [lo, 1] */
	/* hi relative to [lo, 1]: (hn ld - ln hd) / (hd (ld - ln)) */
	TRY(qaws_exact_int_mul(&t1, hn, ld));
	TRY(qaws_exact_int_mul(&t2, ln, hd));
	TRY(qaws_exact_int_sub(&bn, &t1, &t2));
	TRY(qaws_exact_int_sub(&t1, ld, ln));
	if (qaws_exact_int_is_zero(&t1))
	{
		/* lo = 1: the end point */
		unsigned int i, c;
		for (i = 0; i < n; i++)
			for (c = 0; c < 4; c++)
				Q[i * stride + c] = Q[n * stride + c];
		return QAWS_STATUS_OK;
	}
	TRY(qaws_exact_int_mul(&bd, hd, &t1));
	if (qaws_exact_int_cmp(&bn, &bd) != 0)
		TRY(split_at(Q, n, stride, &bn, &bd, 0));   /* keep [0, beta] */
	return QAWS_STATUS_OK;
}

/* double (a dyadic) -> num / den with den a power of two */
static qaws_status double_ratio(double x, qaws_exact_int* num, qaws_exact_int* den)
{
	int64_t m;
	int e;
	qaws_status st;
	qaws_exact_split_double(x, &m, &e);
	qaws_exact_int_from_i64(num, m);
	qaws_exact_int_from_i64(den, 1);
	if (m == 0)
		return QAWS_STATUS_OK;
	if (e >= 0)
		return qaws_exact_int_shl(num, num, (unsigned int)e);
	TRY(qaws_exact_int_shl(den, den, (unsigned int)-e));
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Line / surface                                                     */
/* ------------------------------------------------------------------ */

typedef struct line_ctx
{
	dyadic P[3], d[3], n1[3], n2[3], dd;   /* dd = |d|^2 */
} line_ctx;

/* t enclosure over the sub-patch [u_lo, u_hi] x [v_lo, v_hi] (local, rationals): t = d . (S - P) / |d|^2 */
static qaws_status t_enclose(qaws_exact_surface const* s, qaws_exact_int const* h, line_ctx const* L, qaws_exact_int const* ul, qaws_exact_int const* ud,
	qaws_exact_int const* uh, qaws_exact_int const* uhd, qaws_exact_int const* vl, qaws_exact_int const* vd, qaws_exact_int const* vh,
	qaws_exact_int const* vhd, double* t_lo, double* t_hi, int* exact)
{
	unsigned int p = s->p, q = s->q, n = (p + 1) * (q + 1), a, b, i, c;
	qaws_exact_int* net;
	qaws_status st = QAWS_STATUS_OK;
	net = (qaws_exact_int*)sh_alloc(sizeof(qaws_exact_int) * n * 4);
	if (!net)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n * 4; i++)
		net[i] = h[i];
	/* rows (fixed u index a): the v curves; then columns: the u curves */
	for (a = 0; a <= p && st == QAWS_STATUS_OK; a++)
		st = restrict_curve(&net[a * (q + 1) * 4], q, 4, vl, vd, vh, vhd);
	for (b = 0; b <= q && st == QAWS_STATUS_OK; b++)
		st = restrict_curve(&net[b * 4], p, (q + 1) * 4, ul, ud, uh, uhd);
	*t_lo = HUGE_VAL;
	*t_hi = -HUGE_VAL;
	*exact = 1;
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
	{
		dyadic num, den = { { { 0 }, 0, 0 }, 0 }, w;
		qaws_exact_int nn, dn;
		double lo, hi;
		int ex;
		qaws_exact_int_zero(&num.m);
		num.e = 0;
		w.m = net[i * 4 + 3];
		w.e = 0;
		for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
		{
			dyadic hx, t1, t2;
			hx.m = net[i * 4 + c];
			hx.e = s->space_exp2;
			st = dy_mul(&t1, &L->d[c], &hx);
			if (st == QAWS_STATUS_OK) st = dy_mul(&t2, &L->d[c], &L->P[c]);
			if (st == QAWS_STATUS_OK) st = dy_mul(&t2, &t2, &w);
			if (st == QAWS_STATUS_OK) st = dy_add(&t1, &t1, &t2, -1);
			if (st == QAWS_STATUS_OK) st = dy_add(&num, &num, &t1, 1);
		}
		if (st == QAWS_STATUS_OK) st = dy_mul(&den, &w, &L->dd);
		if (st != QAWS_STATUS_OK)
			break;
		/* num / den = (nm / dm) 2^(ne - de) */
		nn = num.m;
		dn = den.m;
		if (num.e > den.e)
			st = qaws_exact_int_shl(&nn, &nn, (unsigned int)(num.e - den.e));
		else if (den.e > num.e)
			st = qaws_exact_int_shl(&dn, &dn, (unsigned int)(den.e - num.e));
		if (st != QAWS_STATUS_OK)
			break;
		qaws_exact_ratio_enclose(&nn, &dn, &lo, &hi, &ex);
		if (lo < *t_lo) *t_lo = lo;
		if (hi > *t_hi) *t_hi = hi;
		if (!ex || (i > 0 && lo != *t_lo))
			*exact = 0;
	}
	if (*t_lo != *t_hi)
		*exact = 0;
	sh_free(net);
	return st;
}

/* Local parameter of a root interval / exact root, as rationals [lo, hi]. */
static qaws_status root_bounds(qaws_exact_root const* r, qaws_exact_int* ln, qaws_exact_int* ld, qaws_exact_int* hn, qaws_exact_int* hd)
{
	qaws_status st;
	qaws_exact_int_from_i64(ln, (int64_t)r->index);   /* index <= 2^62 */
	qaws_exact_int_from_i64(ld, 1);
	TRY(qaws_exact_int_shl(ld, ld, (unsigned int)r->depth));
	*hd = *ld;
	qaws_exact_int_from_i64(hn, (int64_t)(r->exact ? r->index : r->index + 1));
	return QAWS_STATUS_OK;
}

static qaws_status patch_param(int64_t a, int64_t b, int shift, int exact, qaws_exact_int const* ln, qaws_exact_int const* ld, qaws_exact_int const* hn,
	qaws_exact_int const* hd, double* lo, double* hi, int* ex)
{
	/* t = (a + (b - a) x) 2^-shift at x = n / d */
	qaws_exact_int num, den, t;
	double l2, h2;
	int e1, e2;
	qaws_status st;
	TRY(qaws_exact_int_mul_i64(&num, ld, a));
	TRY(qaws_exact_int_mul_i64(&t, ln, b - a));
	TRY(qaws_exact_int_add(&num, &num, &t));
	TRY(qaws_exact_int_shl(&den, ld, (unsigned int)shift));
	qaws_exact_ratio_enclose(&num, &den, lo, &l2, &e1);
	if (exact)
	{
		*hi = l2;
		*ex = e1;
		return QAWS_STATUS_OK;
	}
	TRY(qaws_exact_int_mul_i64(&num, hd, a));
	TRY(qaws_exact_int_mul_i64(&t, hn, b - a));
	TRY(qaws_exact_int_add(&num, &num, &t));
	TRY(qaws_exact_int_shl(&den, hd, (unsigned int)shift));
	qaws_exact_ratio_enclose(&num, &den, &h2, hi, &e2);
	*ex = 0;
	return QAWS_STATUS_OK;
}

/* One patch: eliminate one parameter, solve, map. */
static qaws_status patch_hits(qaws_exact_surface const* s, unsigned int iu, unsigned int iv, line_ctx const* L, qaws_exact_surface_hit* out,
	unsigned int capacity, unsigned int* count)
{
	qaws_exact_int const* h = s->patch[iu * s->nv + iv];
	unsigned int p = s->p, q = s->q, n = (p + 1) * (q + 1), k, attempt;
	qaws_exact_int* F;
	poly* sys;   /* P (k + 1), Q (k + 1), E (k k), R, C00, C01 */
	qaws_exact_local_hit* lh;
	qaws_status st = QAWS_STATUS_OK;
	int elim_u = p <= q, common = 1;
	unsigned int maxk = p > q ? p : q, nl = 0;
	F = (qaws_exact_int*)sh_alloc(sizeof(qaws_exact_int) * n * 2);
	sys = (poly*)sh_alloc(sizeof(poly) * (2 * (maxk + 1) + maxk * maxk + 3));
	lh = (qaws_exact_local_hit*)sh_alloc(sizeof(qaws_exact_local_hit) * QAWS_EXACT_SOLVE_MAX_ROOTS);
	if (!F || !sys || !lh)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK) st = plane_net(s, h, L->n1, L->P, &F[0]);
	if (st == QAWS_STATUS_OK) st = plane_net(s, h, L->n2, L->P, &F[n]);
	if (st == QAWS_STATUS_OK)
	{
		/* a plane with the whole control net strictly on one side: the line misses the patch */
		unsigned int e, k2;
		for (e = 0; e < 2; e++)
		{
			int s0 = qaws_exact_int_sign(&F[e * n]);
			for (k2 = 1; k2 < n && s0 != 0; k2++)
				if (qaws_exact_int_sign(&F[e * n + k2]) != s0)
					s0 = 0;
			if (s0 != 0)
			{
				sh_free(F);
				sh_free(sys);
				sh_free(lh);
				return QAWS_STATUS_OK;
			}
		}
	}
	for (attempt = 0; attempt < 2 && st == QAWS_STATUS_OK && common; attempt++, elim_u = !elim_u)
	{
		unsigned int kk = elim_u ? p : q;
		poly* P = &sys[0];
		poly* Q = &sys[kk + 1];
		poly* E = &sys[2 * (kk + 1)];
		poly* R = &sys[2 * (kk + 1) + kk * kk];
		st = to_power2(&F[0], p, q, elim_u, P);
		if (st == QAWS_STATUS_OK) st = to_power2(&F[n], p, q, elim_u, Q);
		if (st == QAWS_STATUS_OK && kk == 1)
		{
			/* R = P1 Q0 - P0 Q1, s = -P0 / P1 (or -Q0 / Q1 when P does not involve s) */
			poly t;
			poly const* S1 = qaws_exact_poly_is_zero(&P[1]) ? Q : P;
			st = qaws_exact_poly_mul(&R[0], &P[1], &Q[0]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_mul(&t, &P[0], &Q[1]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_acc(&R[0], &t, -1);
			R[1] = S1[1];
			R[2] = S1[0];
			for (k = 0; k <= R[2].deg; k++)
				qaws_exact_int_neg(&R[2].c[k], &R[2].c[k]);
		}
		else if (st == QAWS_STATUS_OK)
		{
			st = qaws_exact_bezout_poly(P, Q, kk, E);
			if (st == QAWS_STATUS_OK) st = qaws_exact_det_cofactors(E, kk, &R[0], &R[1], &R[2]);
		}
		if (st == QAWS_STATUS_OK)
		{
			st = qaws_exact_solve_system(&R[0], &R[1], &R[2], NULL, 0, lh, QAWS_EXACT_SOLVE_MAX_ROOTS, &nl, &common);
			if (common)
				st = QAWS_STATUS_OK;
		}
		if (st == QAWS_STATUS_OK && !common)
			break;
	}
	if (st == QAWS_STATUS_OK && common)
		st = QAWS_STATUS_CERTIFICATION_FAILED;   /* the line lies on the surface */
	for (k = 0; k < nl && st == QAWS_STATUS_OK; k++)
	{
		qaws_exact_local_hit const* x = &lh[k];
		qaws_exact_int rl, rd, rh, rhd, sl, sd, sh, shd;
		qaws_exact_surface_hit hit;
		double tl = 0, th = 0;
		int r_start = x->r.exact && x->r.index == 0;
		int s_start = x->s_exact && qaws_exact_int_is_zero(&x->s_num);
		int eu = 0, ev = 0, et = 0;
		/* a hit on a patch edge belongs to the patch before it */
		if (elim_u ? ((r_start && iv > 0) || (s_start && iu > 0)) : ((r_start && iu > 0) || (s_start && iv > 0)))
			continue;
		st = root_bounds(&x->r, &rl, &rd, &rh, &rhd);
		if (st != QAWS_STATUS_OK)
			break;
		if (x->s_exact)
		{
			sl = x->s_num;
			sd = x->s_den;
			sh = sl;
			shd = sd;
		}
		else
		{
			st = double_ratio(x->s_lo, &sl, &sd);
			if (st == QAWS_STATUS_OK) st = double_ratio(x->s_hi, &sh, &shd);
			if (st != QAWS_STATUS_OK)
				break;
		}
		if (elim_u)
		{
			st = patch_param(s->ub[iu], s->ub[iu + 1], s->u_shift, x->s_exact, &sl, &sd, &sh, &shd, &hit.u_lo, &hit.u_hi, &eu);
			if (st == QAWS_STATUS_OK)
				st = patch_param(s->vb[iv], s->vb[iv + 1], s->v_shift, x->r.exact, &rl, &rd, &rh, &rhd, &hit.v_lo, &hit.v_hi, &ev);
			if (st == QAWS_STATUS_OK)
				st = t_enclose(s, h, L, &sl, &sd, &sh, &shd, &rl, &rd, &rh, &rhd, &tl, &th, &et);
		}
		else
		{
			st = patch_param(s->ub[iu], s->ub[iu + 1], s->u_shift, x->r.exact, &rl, &rd, &rh, &rhd, &hit.u_lo, &hit.u_hi, &eu);
			if (st == QAWS_STATUS_OK)
				st = patch_param(s->vb[iv], s->vb[iv + 1], s->v_shift, x->s_exact, &sl, &sd, &sh, &shd, &hit.v_lo, &hit.v_hi, &ev);
			if (st == QAWS_STATUS_OK)
				st = t_enclose(s, h, L, &rl, &rd, &rh, &rhd, &sl, &sd, &sh, &shd, &tl, &th, &et);
		}
		if (st != QAWS_STATUS_OK)
			break;
		hit.t_lo = tl;
		hit.t_hi = th;
		hit.kind = (eu && ev && et && x->r.exact && x->s_exact) ? QAWS_EXACT_HIT_POINT : QAWS_EXACT_HIT_CROSSING;
		if (*count >= capacity)
			st = QAWS_STATUS_BUFFER_TOO_SMALL;
		else
			out[(*count)++] = hit;
	}
	sh_free(F);
	sh_free(sys);
	sh_free(lh);
	return st;
}

qaws_status qaws_exact_surface_line_hits(qaws_exact_surface const* surface, double const p0[3], double const p1[3], qaws_exact_surface_hit* out_hits,
	unsigned int capacity, unsigned int* out_count)
{
	line_ctx L;
	dyadic q1[3], e[3];
	unsigned int c, k, iu, iv, count = 0;
	double ad[3];
	qaws_status st = QAWS_STATUS_OK;
	if (!surface || !p0 || !p1 || !out_count || (!out_hits && capacity))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	for (c = 0; c < 3; c++)
	{
		TRY(dy_from_double(&L.P[c], p0[c]));
		TRY(dy_from_double(&q1[c], p1[c]));
		TRY(dy_add(&L.d[c], &q1[c], &L.P[c], -1));
		ad[c] = fabs(p1[c] - p0[c]);
	}
	if (qaws_exact_int_is_zero(&L.d[0].m) && qaws_exact_int_is_zero(&L.d[1].m) && qaws_exact_int_is_zero(&L.d[2].m))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* n1 = d x e_k (k: the smallest |d_k|), n2 = d x n1 */
	k = ad[0] <= ad[1] && ad[0] <= ad[2] ? 0 : (ad[1] <= ad[2] ? 1 : 2);
	for (c = 0; c < 3; c++)
	{
		qaws_exact_int_from_i64(&e[c].m, c == k ? 1 : 0);
		e[c].e = 0;
	}
	for (c = 0; c < 3; c++)
	{
		unsigned int c1 = (c + 1) % 3, c2 = (c + 2) % 3;
		dyadic t1, t2;
		TRY(dy_mul(&t1, &L.d[c1], &e[c2]));
		TRY(dy_mul(&t2, &L.d[c2], &e[c1]));
		TRY(dy_add(&L.n1[c], &t1, &t2, -1));
	}
	for (c = 0; c < 3; c++)
	{
		unsigned int c1 = (c + 1) % 3, c2 = (c + 2) % 3;
		dyadic t1, t2;
		TRY(dy_mul(&t1, &L.d[c1], &L.n1[c2]));
		TRY(dy_mul(&t2, &L.d[c2], &L.n1[c1]));
		TRY(dy_add(&L.n2[c], &t1, &t2, -1));
	}
	qaws_exact_int_zero(&L.dd.m);
	L.dd.e = 0;
	for (c = 0; c < 3; c++)
	{
		dyadic t;
		TRY(dy_mul(&t, &L.d[c], &L.d[c]));
		TRY(dy_add(&L.dd, &L.dd, &t, 1));
	}
	for (iu = 0; iu < surface->nu && st == QAWS_STATUS_OK; iu++)
		for (iv = 0; iv < surface->nv && st == QAWS_STATUS_OK; iv++)
			st = patch_hits(surface, iu, iv, &L, out_hits, capacity, &count);
	/* sort by t */
	for (iu = 1; iu < count; iu++)
	{
		qaws_exact_surface_hit key = out_hits[iu];
		int j = (int)iu - 1;
		while (j >= 0 && out_hits[j].t_lo > key.t_lo)
		{
			out_hits[j + 1] = out_hits[j];
			j--;
		}
		out_hits[j + 1] = key;
	}
	*out_count = count;
	return st;
}
