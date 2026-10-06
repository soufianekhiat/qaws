#include "qaws_exact_poly.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

void qaws_exact_poly_zero(qaws_exact_poly* p, unsigned int deg)
{
	unsigned int i;
	p->deg = deg;
	for (i = 0; i <= deg; i++)
		qaws_exact_int_zero(&p->c[i]);
}

void qaws_exact_poly_const(qaws_exact_poly* p, int64_t v)
{
	p->deg = 0;
	qaws_exact_int_from_i64(&p->c[0], v);
}

int qaws_exact_poly_is_zero(qaws_exact_poly const* p)
{
	unsigned int i;
	for (i = 0; i <= p->deg; i++)
		if (!qaws_exact_int_is_zero(&p->c[i]))
			return 0;
	return 1;
}

void qaws_exact_poly_trim(qaws_exact_poly* p)
{
	while (p->deg > 0 && qaws_exact_int_is_zero(&p->c[p->deg]))
		p->deg--;
}

qaws_status qaws_exact_poly_pad(qaws_exact_poly* p, unsigned int deg)
{
	if (deg >= QAWS_EXACT_POLY_COEF)
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	while (p->deg < deg)
		qaws_exact_int_zero(&p->c[++p->deg]);
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_poly_mul(qaws_exact_poly* r, qaws_exact_poly const* a, qaws_exact_poly const* b)
{
	qaws_exact_int t;
	qaws_exact_poly out;
	unsigned int i, j;
	qaws_status st;
	if (a->deg + b->deg >= QAWS_EXACT_POLY_COEF)
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	qaws_exact_poly_zero(&out, a->deg + b->deg);
	for (i = 0; i <= a->deg; i++)
	{
		if (qaws_exact_int_is_zero(&a->c[i]))
			continue;
		for (j = 0; j <= b->deg; j++)
		{
			if (qaws_exact_int_is_zero(&b->c[j]))
				continue;
			TRY(qaws_exact_int_mul(&t, &a->c[i], &b->c[j]));
			TRY(qaws_exact_int_add(&out.c[i + j], &out.c[i + j], &t));
		}
	}
	*r = out;
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_poly_acc(qaws_exact_poly* r, qaws_exact_poly const* a, int sign)
{
	unsigned int i;
	qaws_status st;
	TRY(qaws_exact_poly_pad(r, a->deg));
	for (i = 0; i <= a->deg; i++)
		TRY(sign > 0 ? qaws_exact_int_add(&r->c[i], &r->c[i], &a->c[i]) : qaws_exact_int_sub(&r->c[i], &r->c[i], &a->c[i]));
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_poly_scale(qaws_exact_poly* r, qaws_exact_poly const* a, qaws_exact_int const* k)
{
	unsigned int i;
	qaws_status st;
	r->deg = a->deg;
	for (i = 0; i <= a->deg; i++)
		TRY(qaws_exact_int_mul(&r->c[i], &a->c[i], k));
	return QAWS_STATUS_OK;
}

int64_t qaws_exact_binom64(unsigned int n, unsigned int k)
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

qaws_status qaws_exact_poly_from_bernstein(qaws_exact_int const* h, unsigned int n, unsigned int D, unsigned int c, qaws_exact_poly* out)
{
	qaws_exact_int t;
	unsigned int i, k;
	qaws_status st;
	qaws_exact_poly_zero(out, n);
	for (k = 0; k <= n; k++)
		for (i = 0; i <= k; i++)
		{
			int64_t f = qaws_exact_binom64(n, i) * qaws_exact_binom64(n - i, k - i) * (((k - i) & 1) ? -1 : 1);
			TRY(qaws_exact_int_mul_i64(&t, &h[i * D + c], f));
			TRY(qaws_exact_int_add(&out->c[k], &out->c[k], &t));
		}
	return QAWS_STATUS_OK;
}

/* b_i = sum_{j <= i} C(i, j) (L / C(N, j)) p_j with L = lcm_j C(N, j) */
qaws_status qaws_exact_poly_to_bernstein(qaws_exact_poly const* p, qaws_exact_int* b)
{
	qaws_exact_int L, g, q, t, cj;
	unsigned int N = p->deg, i, j;
	qaws_status st;
	qaws_exact_int_from_i64(&L, 1);
	for (j = 0; j <= N; j++)
	{
		qaws_exact_int_from_i64(&cj, qaws_exact_binom64(N, j));
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
		qaws_exact_int_from_i64(&cj, qaws_exact_binom64(N, j));
		TRY(qaws_exact_int_divmod(&q, NULL, &L, &cj));
		TRY(qaws_exact_int_mul(&q, &q, &p->c[j]));
		for (i = j; i <= N; i++)
		{
			TRY(qaws_exact_int_mul_i64(&t, &q, qaws_exact_binom64(i, j)));
			TRY(qaws_exact_int_add(&b[i], &b[i], &t));
		}
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_poly_det(qaws_exact_poly const* M, unsigned int n, unsigned int row, unsigned int mask, qaws_exact_poly* out)
{
	qaws_exact_poly* sub;
	qaws_exact_poly* t;
	unsigned int j;
	int sign = 1;
	qaws_status st = QAWS_STATUS_OK;
	if (row == n)
	{
		qaws_exact_poly_const(out, 1);
		return QAWS_STATUS_OK;
	}
	sub = (qaws_exact_poly*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_poly) * 2));
	if (!sub)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	t = sub + 1;
	qaws_exact_poly_zero(out, 0);
	for (j = 0; j < n && st == QAWS_STATUS_OK; j++)
	{
		if (mask & (1u << j))
			continue;
		st = qaws_exact_poly_det(M, n, row + 1, mask | (1u << j), sub);
		if (st == QAWS_STATUS_OK)
			st = qaws_exact_poly_mul(t, &M[row * n + j], sub);
		if (st == QAWS_STATUS_OK)
			st = qaws_exact_poly_acc(out, t, sign);
		sign = -sign;
	}
	qaws_internal_dealloc(NULL, sub);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Modular gcd                                                        */
/* ------------------------------------------------------------------ */

static uint32_t mod_int(qaws_exact_int const* x, uint32_t p)
{
	uint64_t r = 0;
	int i;
	for (i = x->size - 1; i >= 0; i--)
		r = ((r << 32) | x->limb[i]) % p;
	if (x->sign < 0 && r)
		r = p - r;
	return (uint32_t)r;
}

static uint32_t pow_mod(uint64_t b, uint64_t e, uint32_t p)
{
	uint64_t r = 1;
	b %= p;
	while (e)
	{
		if (e & 1)
			r = r * b % p;
		b = b * b % p;
		e >>= 1;
	}
	return (uint32_t)r;
}

/* Deterministic Miller-Rabin for n < 2^32 (bases 2, 7, 61). */
static int is_prime32(uint32_t n)
{
	static uint32_t const bases[3] = { 2, 7, 61 };
	uint32_t d = n - 1;
	int s = 0, i, k;
	if (n < 2 || (n % 2) == 0)
		return n == 2;
	while ((d & 1) == 0)
	{
		d >>= 1;
		s++;
	}
	for (i = 0; i < 3; i++)
	{
		uint64_t x;
		if (bases[i] % n == 0)
			continue;
		x = pow_mod(bases[i], d, n);
		if (x == 1 || x == n - 1)
			continue;
		for (k = 1; k < s; k++)
		{
			x = x * x % n;
			if (x == n - 1)
				break;
		}
		if (k == s)
			return 0;
	}
	return 1;
}

static void trim_mod(uint32_t const* x, unsigned int* d)
{
	while (*d > 0 && x[*d] == 0)
		(*d)--;
}

/* a <- a mod b (b != 0) */
static void rem_mod(uint32_t* a, unsigned int* da, uint32_t const* b, unsigned int db, uint32_t p)
{
	uint64_t inv = pow_mod(b[db], p - 2, p);
	unsigned int i;
	trim_mod(a, da);
	while (!(*da == 0 && a[0] == 0) && *da >= db)
	{
		uint64_t f = (uint64_t)a[*da] * inv % p;
		unsigned int sh = *da - db;
		for (i = 0; i <= db; i++)
			a[i + sh] = (uint32_t)((a[i + sh] + (uint64_t)p - f * b[i] % p) % p);
		if (*da == 0)
			break;
		trim_mod(a, da);
	}
}

/* Monic gcd of a and b mod p into out (a, b destroyed); returns its degree (0: coprime). */
static unsigned int gcd_mod(uint32_t* a, unsigned int da, uint32_t* b, unsigned int db, uint32_t p, uint32_t* out)
{
	uint64_t inv;
	unsigned int i;
	trim_mod(a, &da);
	trim_mod(b, &db);
	while (!(db == 0 && b[0] == 0))
	{
		uint32_t* t = a;
		unsigned int dt;
		rem_mod(a, &da, b, db, p);   /* a: the remainder, degree da */
		a = b;
		b = t;
		dt = da;
		da = db;
		db = dt;
	}
	if (da == 0)
	{
		out[0] = 1;
		return 0;
	}
	inv = pow_mod(a[da], p - 2, p);
	for (i = 0; i <= da; i++)
		out[i] = (uint32_t)((uint64_t)a[i] * inv % p);
	return da;
}

static qaws_status primitive(qaws_exact_poly* p)
{
	qaws_exact_int g;
	unsigned int i;
	qaws_status st;
	qaws_exact_poly_trim(p);
	qaws_exact_int_zero(&g);
	for (i = 0; i <= p->deg; i++)
		qaws_exact_int_gcd(&g, &g, &p->c[i]);
	if (qaws_exact_int_is_zero(&g))
		return QAWS_STATUS_OK;
	if (qaws_exact_int_sign(&p->c[p->deg]) < 0)
		qaws_exact_int_neg(&g, &g);
	for (i = 0; i <= p->deg; i++)
		TRY(qaws_exact_int_divmod(&p->c[i], NULL, &p->c[i], &g));
	return QAWS_STATUS_OK;
}

/* Does g divide a exactly over the integers (g primitive)? */
static qaws_status divides(qaws_exact_poly const* a, qaws_exact_poly const* g, int* out)
{
	qaws_exact_poly r = *a;
	qaws_exact_int q, rem, t;
	int k;
	unsigned int i;
	qaws_status st;
	*out = 0;
	if (g->deg > r.deg)
	{
		*out = qaws_exact_poly_is_zero(&r);
		return QAWS_STATUS_OK;
	}
	for (k = (int)(r.deg - g->deg); k >= 0; k--)
	{
		qaws_exact_int* lead = &r.c[(unsigned int)k + g->deg];
		if (qaws_exact_int_is_zero(lead))
			continue;
		TRY(qaws_exact_int_divmod(&q, &rem, lead, &g->c[g->deg]));
		if (!qaws_exact_int_is_zero(&rem))
			return QAWS_STATUS_OK;
		for (i = 0; i <= g->deg; i++)
		{
			TRY(qaws_exact_int_mul(&t, &q, &g->c[i]));
			TRY(qaws_exact_int_sub(&r.c[(unsigned int)k + i], &r.c[(unsigned int)k + i], &t));
		}
	}
	*out = qaws_exact_poly_is_zero(&r);
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_poly_gcd(qaws_exact_poly const* a, qaws_exact_poly const* b, qaws_exact_poly* out)
{
	qaws_exact_poly* A = NULL;
	qaws_exact_poly* B = NULL;
	qaws_exact_poly* G = NULL;
	qaws_exact_poly* cand = NULL;
	qaws_exact_int gamma, M, t, half;
	uint32_t ap[QAWS_EXACT_POLY_COEF], bp[QAWS_EXACT_POLY_COEF], gp[QAWS_EXACT_POLY_COEF];
	uint32_t p = 0x7FFFFFFFu;
	unsigned int dbest = 0, i, used = 0;
	int have = 0;
	qaws_status st = QAWS_STATUS_OK;
	A = (qaws_exact_poly*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_poly) * 4));
	if (!A)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	B = A + 1;
	G = A + 2;
	cand = A + 3;
	*A = *a;
	*B = *b;
	st = primitive(A);
	if (st == QAWS_STATUS_OK) st = primitive(B);
	if (st == QAWS_STATUS_OK && (qaws_exact_poly_is_zero(A) || qaws_exact_poly_is_zero(B)))
		st = QAWS_STATUS_INVALID_ARGUMENT;
	if (st == QAWS_STATUS_OK && (A->deg == 0 || B->deg == 0))
	{
		qaws_exact_poly_const(out, 1);
		qaws_internal_dealloc(NULL, A);
		return QAWS_STATUS_OK;
	}
	if (st == QAWS_STATUS_OK)
		qaws_exact_int_gcd(&gamma, &A->c[A->deg], &B->c[B->deg]);
	while (st == QAWS_STATUS_OK)
	{
		unsigned int d;
		uint32_t gm;
		/* the next prime not dividing the leading coefficients */
		do
			p -= 2;
		while (p > 3 && (!is_prime32(p) || mod_int(&A->c[A->deg], p) == 0 || mod_int(&B->c[B->deg], p) == 0));
		if (p <= 3 || ++used > 64)
		{
			st = QAWS_STATUS_EXACT_RANGE_EXCEEDED;
			break;
		}
		for (i = 0; i <= A->deg; i++) ap[i] = mod_int(&A->c[i], p);
		for (i = 0; i <= B->deg; i++) bp[i] = mod_int(&B->c[i], p);
		d = gcd_mod(ap, A->deg, bp, B->deg, p, gp);
		if (d == 0)
		{
			qaws_exact_poly_const(out, 1);   /* coprime modulo p (lc not 0): coprime over Z */
			break;
		}
		gm = mod_int(&gamma, p);
		for (i = 0; i <= d; i++)
			gp[i] = (uint32_t)((uint64_t)gp[i] * gm % p);
		if (have && d > dbest)
			continue;   /* an unlucky prime */
		if (!have || d < dbest)
		{
			dbest = d;
			have = 1;
			G->deg = d;
			for (i = 0; i <= d; i++)
				qaws_exact_int_from_i64(&G->c[i], (int64_t)gp[i]);
			qaws_exact_int_from_i64(&M, (int64_t)p);
		}
		else
		{
			/* Chinese remaindering: G_i += M ((gp_i - G_i) / M mod p) */
			uint32_t minv = pow_mod(mod_int(&M, p), p - 2, p);
			for (i = 0; i <= d && st == QAWS_STATUS_OK; i++)
			{
				uint32_t gi = mod_int(&G->c[i], p);
				uint64_t k = (uint64_t)((gp[i] + (uint64_t)p - gi) % p) * minv % p;
				st = qaws_exact_int_mul_i64(&t, &M, (int64_t)k);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&G->c[i], &G->c[i], &t);
			}
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul_i64(&M, &M, (int64_t)p);
			if (st != QAWS_STATUS_OK)
				break;
		}
		/* candidate: symmetric residues, primitive part, then the proof by division */
		qaws_exact_int_shr(&half, &M, 1);
		cand->deg = G->deg;
		for (i = 0; i <= G->deg && st == QAWS_STATUS_OK; i++)
		{
			cand->c[i] = G->c[i];
			if (qaws_exact_int_cmp(&cand->c[i], &half) > 0)
				st = qaws_exact_int_sub(&cand->c[i], &cand->c[i], &M);
		}
		if (st == QAWS_STATUS_OK) st = primitive(cand);
		if (st == QAWS_STATUS_OK && cand->deg == dbest)
		{
			int da = 0, db = 0;
			st = divides(A, cand, &da);
			if (st == QAWS_STATUS_OK && da) st = divides(B, cand, &db);
			if (st == QAWS_STATUS_OK && da && db)
			{
				*out = *cand;
				break;
			}
		}
		if (st == QAWS_STATUS_OK && qaws_exact_int_bits(&M) > QAWS_EXACT_MAX_BITS - 96)
			st = QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	}
	qaws_internal_dealloc(NULL, A);
	return st;
}
