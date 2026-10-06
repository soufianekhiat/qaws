#include "qaws_exact_surface.h"
#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define CS_MAX_SIZE 216          /* (p + 1)(q + 1)(m + 1) */
#define CS_MAX_SPLITS 96         /* subdivision depth, all directions */
#define CS_MAX_BOXES 400000
#define CS_TARGET_DEPTH 44       /* refinement: about 2^-44 of a span / patch */

/*
 * Three tensors F_c (c = x, y, z) of Bernstein coefficients on a box of
 * the local (u, v, t) cube, degrees n[0..2]; index ((a (n1 + 1) + b) (n2 + 1) + k).
 */
typedef struct cs_box
{
	uint64_t lo[3], hi[3];   /* the box: [lo, hi] 2^-dep per direction */
	int dep[3];
	qaws_exact_int* F;   /* 3 * size */
} cs_box;

typedef struct cs_ctx
{
	unsigned int n[3], size, stride[3];
} cs_ctx;

static void* cs_alloc(size_t b)
{
	return qaws_internal_alloc(NULL, (unsigned long)b);
}

static void cs_free(void* p)
{
	qaws_internal_dealloc(NULL, p);
}

/* Divides one tensor by the gcd of its coefficients (a positive factor). */
/*
 * Divides a tensor by the common power of two of its coefficients (a
 * positive factor: signs and roots are kept). Subdivision only ever
 * introduces powers of two (2^n halving, 16^n the 7/16 cut), so this keeps
 * the integers as small as a full gcd would, at the cost of a shift.
 */
static qaws_status normalize(qaws_exact_int* T, unsigned int size)
{
	unsigned int i, k = ~0u;
	for (i = 0; i < size; i++)
		if (!qaws_exact_int_is_zero(&T[i]))
		{
			unsigned int z = qaws_exact_int_ctz(&T[i]);
			if (z < k)
				k = z;
			if (k == 0)
				return QAWS_STATUS_OK;
		}
	if (k == ~0u)
		return QAWS_STATUS_OK;
	for (i = 0; i < size; i++)
		qaws_exact_int_shr(&T[i], &T[i], k);
	return QAWS_STATUS_OK;
}

/* Halves the tensor T along direction d at 1/2 into L and R (both scaled by 2^n_d). */
static qaws_status split_dir(cs_ctx const* cx, qaws_exact_int const* T, unsigned int d, qaws_exact_int* L, qaws_exact_int* R)
{
	unsigned int n = cx->n[d], st_d = cx->stride[d], i, r, o;
	qaws_exact_int line[17];
	qaws_status st;
	for (o = 0; o < cx->size; o++)
	{
		/* o is the start of a line along d when its d-index is 0 */
		if ((o / st_d) % (n + 1) != 0)
			continue;
		for (i = 0; i <= n; i++)
			line[i] = T[o + i * st_d];
		L[o] = line[0];
		R[o + n * st_d] = line[n];
		for (r = 1; r <= n; r++)
		{
			for (i = 0; i + r <= n; i++)
				TRY(qaws_exact_int_add(&line[i], &line[i], &line[i + 1]));
			TRY(qaws_exact_int_shl(&L[o + r * st_d], &line[0], n - r));
			TRY(qaws_exact_int_shl(&R[o + (n - r) * st_d], &line[n - r], n - r));
		}
		TRY(qaws_exact_int_shl(&L[o], &L[o], n));
		TRY(qaws_exact_int_shl(&R[o + n * st_d], &R[o + n * st_d], n));
	}
	return QAWS_STATUS_OK;
}

/*
 * T along direction d restricted to [x, 1] (keep_right) or [0, x], x = a / b:
 * integer De Casteljau with the weights (b - a, a); every point times b^n.
 */
static qaws_status split_dir_at(cs_ctx const* cx, qaws_exact_int* T, unsigned int d, int64_t a, int64_t b, int keep_right)
{
	unsigned int n = cx->n[d], st_d = cx->stride[d], i, r, o, k;
	qaws_exact_int line[17], out[17], t1, t2;
	qaws_status st;
	for (o = 0; o < cx->size; o++)
	{
		if ((o / st_d) % (n + 1) != 0)
			continue;
		for (i = 0; i <= n; i++)
			line[i] = T[o + i * st_d];
		if (keep_right)
			out[n] = line[n];
		else
			out[0] = line[0];
		for (r = 1; r <= n; r++)
		{
			for (i = 0; i + r <= n; i++)
			{
				TRY(qaws_exact_int_mul_i64(&t1, &line[i], b - a));
				TRY(qaws_exact_int_mul_i64(&t2, &line[i + 1], a));
				TRY(qaws_exact_int_add(&line[i], &t1, &t2));
			}
			if (keep_right)
				out[n - r] = line[n - r];
			else
				out[r] = line[0];
		}
		/* point i came from level (keep_right ? n - i : i): bring all to b^n */
		for (i = 0; i <= n; i++)
		{
			unsigned int lev = keep_right ? n - i : i;
			for (k = lev; k < n; k++)
				TRY(qaws_exact_int_mul_i64(&out[i], &out[i], b));
			T[o + i * st_d] = out[i];
		}
	}
	return QAWS_STATUS_OK;
}

/* Some equation of one strict sign over the box: no root. */
static int excluded(cs_ctx const* cx, qaws_exact_int const* F)
{
	unsigned int c, i;
	for (c = 0; c < 3; c++)
	{
		int s = qaws_exact_int_sign(&F[c * cx->size]);
		if (s == 0)
			continue;
		for (i = 1; i < cx->size; i++)
			if (qaws_exact_int_sign(&F[c * cx->size + i]) != s)
				break;
		if (i == cx->size)
			return 1;
	}
	return 0;
}

/* Value of a double tensor at the box center (de Casteljau at 1/2 in every direction). */
static double center_value(double* T, unsigned int const* n)
{
	unsigned int a, b, k, r;
	unsigned int n0 = n[0], n1 = n[1], n2 = n[2];
	/* along t */
	for (a = 0; a <= n0; a++)
		for (b = 0; b <= n1; b++)
		{
			double* L = &T[(a * (n1 + 1) + b) * (n2 + 1)];
			for (r = 1; r <= n2; r++)
				for (k = 0; k + r <= n2; k++)
					L[k] = 0.5 * (L[k] + L[k + 1]);
		}
	/* along v */
	for (a = 0; a <= n0; a++)
		for (r = 1; r <= n1; r++)
			for (b = 0; b + r <= n1; b++)
				T[(a * (n1 + 1) + b) * (n2 + 1)] = 0.5 * (T[(a * (n1 + 1) + b) * (n2 + 1)] + T[(a * (n1 + 1) + b + 1) * (n2 + 1)]);
	/* along u */
	for (r = 1; r <= n0; r++)
		for (a = 0; a + r <= n0; a++)
			T[a * (n1 + 1) * (n2 + 1)] = 0.5 * (T[a * (n1 + 1) * (n2 + 1)] + T[(a + 1) * (n1 + 1) * (n2 + 1)]);
	return T[0];
}

/* Jacobian of the (scaled) equations at the box center, in doubles: J[c][j]. Each row scaled by 2^-shift[c]. */
static void center_jacobian(cs_ctx const* cx, qaws_exact_int const* F, double J[3][3], unsigned int shift[3])
{
	double* tmp = (double*)cs_alloc(sizeof(double) * CS_MAX_SIZE);
	unsigned int c, j, i;
	for (c = 0; c < 3; c++)
	{
		unsigned int maxb = 0;
		for (i = 0; i < cx->size; i++)
			if (qaws_exact_int_bits(&F[c * cx->size + i]) > maxb)
				maxb = qaws_exact_int_bits(&F[c * cx->size + i]);
		shift[c] = maxb > 900 ? maxb - 900 : 0;
		for (j = 0; j < 3; j++)
		{
			/* derivative tensor along j: n_j (T[k + 1] - T[k]), degree n_j - 1 in j */
			unsigned int nd[3] = { cx->n[0], cx->n[1], cx->n[2] }, a, b, k, o = 0;
			nd[j] = cx->n[j] - 1;
			for (a = 0; a <= nd[0]; a++)
				for (b = 0; b <= nd[1]; b++)
					for (k = 0; k <= nd[2]; k++)
					{
						unsigned int p0 = (a * (cx->n[1] + 1) + b) * (cx->n[2] + 1) + k;
						qaws_exact_int x0 = F[c * cx->size + p0], x1 = F[c * cx->size + p0 + cx->stride[j]];
						qaws_exact_int_shr(&x0, &x0, shift[c]);
						qaws_exact_int_shr(&x1, &x1, shift[c]);
						tmp[o++] = (double)cx->n[j] * (qaws_exact_int_to_double(&x1) - qaws_exact_int_to_double(&x0));
					}
			J[c][j] = center_value(tmp, nd);
		}
	}
	cs_free(tmp);
}

static int invert3(double const A[3][3], double Y[3][3])
{
	double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
	             A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
	unsigned int i, j;
	if (!(fabs(det) > 0) || !(det - det == 0))
		return 0;
	Y[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
	Y[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / det;
	Y[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;
	Y[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / det;
	Y[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
	Y[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / det;
	Y[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
	Y[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / det;
	Y[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			if (!(Y[i][j] - Y[i][j] == 0))
				return 0;
	return 1;
}

/*
 * G = Y F with integer Y (any Y keeps the test exact): G_i = sum_c Yint_ic F_c 2^(smax - s_c).
 * Then: existence (Miranda: G_i <= 0 on one face of direction i, >= 0 on the
 * other) and, if want_unique, uniqueness (the interval Jacobian of G,
 * bounded by Bernstein coefficients, strictly diagonally dominant).
 */
static qaws_status precondition_test(cs_ctx const* cx, qaws_exact_int const* F, int want_unique, int* out)
{
	double J[3][3], Y[3][3];
	int64_t Yi[3][3];
	unsigned int shift[3], smax = 0, i, c, j, o;
	qaws_exact_int* G;
	qaws_status st = QAWS_STATUS_OK;
	*out = 0;
	center_jacobian(cx, F, J, shift);
	if (!invert3(J, Y))
		return QAWS_STATUS_OK;
	/* integer Y, each row on its own scale (a row's small entries must not all round away) */
	for (i = 0; i < 3; i++)
	{
		double ymax = 0;
		int ys;
		if (shift[i] > smax) smax = shift[i];
		for (j = 0; j < 3; j++)
			if (fabs(Y[i][j]) > ymax) ymax = fabs(Y[i][j]);
		frexp(ymax, &ys);
		for (j = 0; j < 3; j++)
			Yi[i][j] = (int64_t)nearbyint(ldexp(Y[i][j], 30 - ys));
	}
	/* the existence and uniqueness tests hold for G = Y F only with Y invertible: check det Y != 0 exactly */
	{
		qaws_exact_int det, t, u;
		static unsigned int const perm[6][3] = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 0, 2, 1 }, { 2, 1, 0 }, { 1, 0, 2 } };
		unsigned int p;
		qaws_exact_int_zero(&det);
		for (p = 0; p < 6; p++)
		{
			qaws_exact_int_from_i64(&t, Yi[0][perm[p][0]]);
			qaws_exact_int_mul_i64(&t, &t, Yi[1][perm[p][1]]);
			qaws_exact_int_mul_i64(&u, &t, Yi[2][perm[p][2]]);
			if (p < 3)
				qaws_exact_int_add(&det, &det, &u);
			else
				qaws_exact_int_sub(&det, &det, &u);
		}
		if (qaws_exact_int_is_zero(&det))
			want_unique = -1;   /* still usable for exclusion only */
	}
	G = (qaws_exact_int*)cs_alloc(sizeof(qaws_exact_int) * 3 * cx->size);
	if (!G)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 && st == QAWS_STATUS_OK; i++)
		for (o = 0; o < cx->size && st == QAWS_STATUS_OK; o++)
		{
			qaws_exact_int_zero(&G[i * cx->size + o]);
			for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
			{
				qaws_exact_int t;
				int64_t y = Yi[i][c];
				if (y == 0)
					continue;
				st = qaws_exact_int_mul_i64(&t, &F[c * cx->size + o], y);
				if (st == QAWS_STATUS_OK && smax > shift[c]) st = qaws_exact_int_shl(&t, &t, smax - shift[c]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&G[i * cx->size + o], &G[i * cx->size + o], &t);
			}
		}
	if (st == QAWS_STATUS_OK)
	{
		int ok = 1;
		/* preconditioned exclusion: a root of F is a root of G, so some G_i of one strict sign proves the box empty */
		for (i = 0; i < 3; i++)
		{
			int s0 = qaws_exact_int_sign(&G[i * cx->size]);
			for (o = 1; o < cx->size && s0 != 0; o++)
				if (qaws_exact_int_sign(&G[i * cx->size + o]) != s0)
					s0 = 0;
			if (s0 != 0)
			{
				cs_free(G);
				*out = -1;
				return QAWS_STATUS_OK;
			}
		}
		/* a singular integer Y proves nothing about F's roots */
		if (want_unique < 0)
			ok = 0;
		/* Miranda: G_i against the faces of direction i */
		for (i = 0; i < 3 && ok; i++)
		{
			int lo_le = 1, lo_ge = 1, hi_le = 1, hi_ge = 1;
			for (o = 0; o < cx->size; o++)
			{
				unsigned int di = (o / cx->stride[i]) % (cx->n[i] + 1);
				int s = qaws_exact_int_sign(&G[i * cx->size + o]);
				if (di == 0)
				{
					if (s > 0) lo_le = 0;
					if (s < 0) lo_ge = 0;
				}
				if (di == cx->n[i])
				{
					if (s > 0) hi_le = 0;
					if (s < 0) hi_ge = 0;
				}
			}
			ok = (lo_le && hi_ge) || (lo_ge && hi_le);
		}
		/* uniqueness: strict diagonal dominance of the interval Jacobian (rows of G) */
		for (i = 0; i < 3 && ok && want_unique; i++)
		{
			qaws_exact_int lower, sum, mn[3], mx[3];
			for (j = 0; j < 3 && st == QAWS_STATUS_OK; j++)
			{
				int first = 1;
				for (o = 0; o < cx->size && st == QAWS_STATUS_OK; o++)
				{
					qaws_exact_int dlt;
					if ((o / cx->stride[j]) % (cx->n[j] + 1) == cx->n[j])
						continue;
					st = qaws_exact_int_sub(&dlt, &G[i * cx->size + o + cx->stride[j]], &G[i * cx->size + o]);
					if (st != QAWS_STATUS_OK)
						break;
					if (first || qaws_exact_int_cmp(&dlt, &mn[j]) < 0) mn[j] = dlt;
					if (first || qaws_exact_int_cmp(&dlt, &mx[j]) > 0) mx[j] = dlt;
					first = 0;
				}
				/* the derivative is n_j times the differences */
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul_i64(&mn[j], &mn[j], (int64_t)cx->n[j]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul_i64(&mx[j], &mx[j], (int64_t)cx->n[j]);
			}
			if (st != QAWS_STATUS_OK)
				break;
			if (qaws_exact_int_sign(&mn[i]) > 0)
				lower = mn[i];
			else if (qaws_exact_int_sign(&mx[i]) < 0)
				qaws_exact_int_neg(&lower, &mx[i]);
			else
			{
				ok = 0;
				break;
			}
			qaws_exact_int_zero(&sum);
			for (j = 0; j < 3 && st == QAWS_STATUS_OK; j++)
			{
				qaws_exact_int a, b;
				if (j == i)
					continue;
				a = mn[j];
				b = mx[j];
				if (qaws_exact_int_sign(&a) < 0) qaws_exact_int_neg(&a, &a);
				if (qaws_exact_int_sign(&b) < 0) qaws_exact_int_neg(&b, &b);
				st = qaws_exact_int_add(&sum, &sum, qaws_exact_int_cmp(&a, &b) > 0 ? &a : &b);
			}
			if (st == QAWS_STATUS_OK && qaws_exact_int_cmp(&lower, &sum) <= 0)
				ok = 0;
		}
		*out = ok && st == QAWS_STATUS_OK;
	}
	cs_free(G);
	return st;
}

/* The halves of box B along d: their tensors, normalized. */
static qaws_status split_box(cs_ctx const* cx, cs_box const* B, unsigned int d, cs_box* L, cs_box* R)
{
	unsigned int c, k;
	qaws_status st;
	for (c = 0; c < 3; c++)
	{
		TRY(split_dir(cx, &B->F[c * cx->size], d, &L->F[c * cx->size], &R->F[c * cx->size]));
		TRY(normalize(&L->F[c * cx->size], cx->size));
		TRY(normalize(&R->F[c * cx->size], cx->size));
	}
	for (k = 0; k < 3; k++)
	{
		L->lo[k] = R->lo[k] = B->lo[k];
		L->hi[k] = R->hi[k] = B->hi[k];
		L->dep[k] = R->dep[k] = B->dep[k];
	}
	L->lo[d] = 2 * B->lo[d];
	L->hi[d] = B->lo[d] + B->hi[d];
	R->lo[d] = B->lo[d] + B->hi[d];
	R->hi[d] = 2 * B->hi[d];
	L->dep[d] = R->dep[d] = B->dep[d] + 1;
	return QAWS_STATUS_OK;
}

static unsigned int shallowest(cs_box const* B)
{
	unsigned int d = 0, k;
	for (k = 1; k < 3; k++)
		if (B->dep[k] < B->dep[d])
			d = k;
	return d;
}

/* Shrink a box holding exactly one root: keep the half proven to hold it (exclusion of the other, or Miranda). */
static qaws_status refine(cs_ctx const* cx, cs_box* B, cs_box* L, cs_box* R)
{
	qaws_status st;
	int stuck[3] = { 0, 0, 0 };
	for (;;)
	{
		unsigned int d = 3, k;
		int inl = 0, inr = 0;
		cs_box* keep;
		/* the shallowest direction that is not stuck */
		for (k = 0; k < 3; k++)
			if (!stuck[k] && B->dep[k] < CS_TARGET_DEPTH && (d == 3 || B->dep[k] < B->dep[d]))
				d = k;
		if (d == 3)
			return QAWS_STATUS_OK;
		st = split_box(cx, B, d, L, R);
		if (st == QAWS_STATUS_EXACT_RANGE_EXCEEDED)
			return QAWS_STATUS_OK;   /* the budget: keep the certified box */
		if (st != QAWS_STATUS_OK)
			return st;
		if (excluded(cx, L->F))
			keep = R;
		else if (excluded(cx, R->F))
			keep = L;
		else
		{
			/* 1: the root is there (Miranda); -1: proven empty, so the root is in the other half */
			TRY(precondition_test(cx, L->F, 0, &inl));
			if (inl == -1)
			{
				inl = 0;
				inr = 1;
			}
			else if (inl == 0)
			{
				TRY(precondition_test(cx, R->F, 0, &inr));
				if (inr == -1)
				{
					inr = 0;
					inl = 1;
				}
			}
			if (!inl && !inr)
			{
				/* the root sits near the cut: the middle half [1/4, 3/4] holds it well inside */
				unsigned int c;
				int inm = 0;
				memcpy(L->F, B->F, sizeof(qaws_exact_int) * 3 * cx->size);
				for (c = 0; c < 3; c++)
				{
					st = split_dir_at(cx, &L->F[c * cx->size], d, 1, 4, 1);   /* [1/4, 1] */
					if (st == QAWS_STATUS_OK) st = split_dir_at(cx, &L->F[c * cx->size], d, 2, 3, 0);   /* its [0, 2/3]: [1/4, 3/4] */
					if (st == QAWS_STATUS_OK) st = normalize(&L->F[c * cx->size], cx->size);
					if (st != QAWS_STATUS_OK)
						return st == QAWS_STATUS_EXACT_RANGE_EXCEEDED ? QAWS_STATUS_OK : st;
				}
				TRY(precondition_test(cx, L->F, 0, &inm));
				if (inm != 1)
				{
					/* this direction is stuck for now: shrink the others first */
					stuck[d] = 1;
					continue;
				}
				memcpy(L->lo, B->lo, sizeof(B->lo));
				memcpy(L->hi, B->hi, sizeof(B->hi));
				memcpy(L->dep, B->dep, sizeof(B->dep));
				L->lo[d] = 3 * B->lo[d] + B->hi[d];
				L->hi[d] = B->lo[d] + 3 * B->hi[d];
				L->dep[d] = B->dep[d] + 2;
				inl = 1;
			}
			keep = inl ? L : R;
		}
		memcpy(B->lo, keep->lo, sizeof(B->lo));
		memcpy(B->hi, keep->hi, sizeof(B->hi));
		memcpy(B->dep, keep->dep, sizeof(B->dep));
		memcpy(B->F, keep->F, sizeof(qaws_exact_int) * 3 * cx->size);
		/* a smaller box linearizes better: stuck directions get another try */
		for (k = 0; k < 3; k++)
			if (k != d)
				stuck[k] = 0;
	}
}

static qaws_status local_to_param(int64_t a, int64_t b, int shift, uint64_t ilo, uint64_t ihi, int dep, double* lo, double* hi)
{
	qaws_exact_span sp;
	int e1, e2;
	qaws_status st;
	sp.degree = 1;
	sp.a = a;
	sp.b = b;
	sp.h = NULL;
	TRY(qaws_exact_span_param_to_double(&sp, shift, ilo, dep, lo, &e1));
	TRY(qaws_exact_span_param_to_double(&sp, shift, ihi, dep, hi, &e2));
	if (!e1) *lo = nextafter(*lo, -HUGE_VAL);
	if (!e2) *hi = nextafter(*hi, HUGE_VAL);
	return QAWS_STATUS_OK;
}

static int overlap(double a0, double a1, double b0, double b1)
{
	return a0 <= b1 && b0 <= a1;
}


/* Is the box B inside the box O (both dyadic [lo, hi] 2^-dep per direction)? Exact. */
static int box_inside(cs_box const* B, qaws_exact_box3 const* O)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		/* O.lo / 2^od <= B.lo / 2^bd and B.hi / 2^bd <= O.hi / 2^od */
		qaws_exact_int x, y;
		qaws_exact_int_from_i64(&x, (int64_t)O->lo[k]);
		qaws_exact_int_from_i64(&y, (int64_t)B->lo[k]);
		qaws_exact_int_shl(&x, &x, (unsigned int)B->dep[k]);
		qaws_exact_int_shl(&y, &y, (unsigned int)O->dep[k]);
		if (qaws_exact_int_cmp(&x, &y) > 0)
			return 0;
		qaws_exact_int_from_i64(&x, (int64_t)B->hi[k]);
		qaws_exact_int_from_i64(&y, (int64_t)O->hi[k]);
		qaws_exact_int_shl(&x, &x, (unsigned int)O->dep[k]);
		qaws_exact_int_shl(&y, &y, (unsigned int)B->dep[k]);
		if (qaws_exact_int_cmp(&x, &y) > 0)
			return 0;
	}
	return 1;
}

/* Do two local boxes (dyadic [lo, hi] 2^-dep per direction) intersect? Exact. */
static int box_overlap(qaws_exact_box3 const* a, qaws_exact_box3 const* b)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		/* a.lo / 2^da <= b.hi / 2^db and b.lo / 2^db <= a.hi / 2^da */
		qaws_exact_int x, y;
		qaws_exact_int_from_i64(&x, (int64_t)a->lo[k]);
		qaws_exact_int_from_i64(&y, (int64_t)b->hi[k]);
		qaws_exact_int_shl(&x, &x, (unsigned int)b->dep[k]);
		qaws_exact_int_shl(&y, &y, (unsigned int)a->dep[k]);
		if (qaws_exact_int_cmp(&x, &y) > 0)
			return 0;
		qaws_exact_int_from_i64(&x, (int64_t)b->lo[k]);
		qaws_exact_int_from_i64(&y, (int64_t)a->hi[k]);
		qaws_exact_int_shl(&x, &x, (unsigned int)a->dep[k]);
		qaws_exact_int_shl(&y, &y, (unsigned int)b->dep[k]);
		if (qaws_exact_int_cmp(&x, &y) > 0)
			return 0;
	}
	return 1;
}

qaws_status qaws_exact_solve3(unsigned int const n[3], qaws_exact_int const* F, qaws_exact_box3* out, unsigned int capacity, unsigned int* out_count,
	unsigned int max_boxes)
{
	cs_ctx cx;
	unsigned int nstack, boxes = 0, depth_cap = CS_MAX_SPLITS, count = 0, c, ncov = 0;
	cs_box* stack = NULL;
	qaws_exact_int* store = NULL;
	qaws_exact_int* tstore = NULL;
	qaws_exact_int* mstore = NULL;
	qaws_exact_box3 cov[256];
	cs_box tmpL, tmpR;
	qaws_status st = QAWS_STATUS_OK;
	*out_count = 0;
	cx.n[0] = n[0];
	cx.n[1] = n[1];
	cx.n[2] = n[2];
	cx.size = (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
	cx.stride[2] = 1;
	cx.stride[1] = n[2] + 1;
	cx.stride[0] = (n[1] + 1) * (n[2] + 1);
	if (cx.size > CS_MAX_SIZE || n[0] > 16 || n[1] > 16 || n[2] > 16)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	stack = (cs_box*)cs_alloc(sizeof(cs_box) * (2 * depth_cap + 4));
	store = (qaws_exact_int*)cs_alloc(sizeof(qaws_exact_int) * 3 * cx.size * (2 * depth_cap + 4));
	tstore = (qaws_exact_int*)cs_alloc(sizeof(qaws_exact_int) * 3 * cx.size * 2);
	mstore = (qaws_exact_int*)cs_alloc(sizeof(qaws_exact_int) * 3 * cx.size);
	if (!stack || !store || !tstore || !mstore)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK)
	{
		tmpL.F = tstore;
		tmpR.F = tstore + 3 * cx.size;
		memcpy(store, F, sizeof(qaws_exact_int) * 3 * cx.size);
		for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
			st = normalize(&store[c * cx.size], cx.size);
		nstack = 1;
		stack[0].F = store;
		stack[0].lo[0] = stack[0].lo[1] = stack[0].lo[2] = 0;
		stack[0].hi[0] = stack[0].hi[1] = stack[0].hi[2] = 1;
		stack[0].dep[0] = stack[0].dep[1] = stack[0].dep[2] = 0;
		while (st == QAWS_STATUS_OK && nstack > 0)
		{
			cs_box B = stack[--nstack];
			int ok = 0;
			if (++boxes > max_boxes)
			{
				st = QAWS_STATUS_CERTIFICATION_FAILED;
				break;
			}
			if (excluded(&cx, B.F))
				continue;
			{
				/* inside a region already holding exactly one (recorded) root */
				unsigned int q;
				int inside = 0;
				for (q = 0; q < ncov && !inside; q++)
					inside = box_inside(&B, &cov[q]);
				if (inside)
					continue;
			}
			st = precondition_test(&cx, B.F, 1, &ok);
			if (st != QAWS_STATUS_OK)
				break;
			if (ok == -1)
				continue;   /* proven empty after preconditioning */
			if (ok == 1)
			{
				/* one root: shrink it in place, then record once (a root on a cut is found twice) */
				qaws_exact_box3 b;
				unsigned int q, dup = 0, k;
				st = refine(&cx, &B, &tmpL, &tmpR);
				if (st != QAWS_STATUS_OK)
					break;
				for (k = 0; k < 3; k++)
				{
					b.lo[k] = B.lo[k];
					b.hi[k] = B.hi[k];
					b.dep[k] = B.dep[k];
				}
				for (q = 0; q < count && !dup; q++)
					dup = box_overlap(&out[q], &b);
				if (dup)
					continue;
				if (count >= capacity)
				{
					st = QAWS_STATUS_BUFFER_TOO_SMALL;
					break;
				}
				out[count++] = b;
				continue;
			}
			if ((unsigned int)(B.dep[0] + B.dep[1] + B.dep[2]) >= depth_cap)
			{
				st = QAWS_STATUS_CERTIFICATION_FAILED;   /* a tangency or an overlap */
				break;
			}
			if (B.dep[0] + B.dep[1] + B.dep[2] >= 3 && ncov < 256)
			{
				/* a root sitting on one of B's cuts (a symmetric or exactly dyadic position)
				   never certifies in the halves: try B's middle [1/4, 3/4]^3, where it is well inside */
				cs_box M;
				unsigned int k;
				int okm = 0;
				M.F = mstore;
				memcpy(M.F, B.F, sizeof(qaws_exact_int) * 3 * cx.size);
				for (k = 0; k < 3 && st == QAWS_STATUS_OK; k++)
					for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
					{
						st = split_dir_at(&cx, &M.F[c * cx.size], k, 1, 4, 1);
						if (st == QAWS_STATUS_OK) st = split_dir_at(&cx, &M.F[c * cx.size], k, 2, 3, 0);
						if (st == QAWS_STATUS_OK) st = normalize(&M.F[c * cx.size], cx.size);
					}
				if (st == QAWS_STATUS_EXACT_RANGE_EXCEEDED)
					st = QAWS_STATUS_OK;
				else if (st == QAWS_STATUS_OK)
				{
					for (k = 0; k < 3; k++)
					{
						M.lo[k] = 3 * B.lo[k] + B.hi[k];
						M.hi[k] = B.lo[k] + 3 * B.hi[k];
						M.dep[k] = B.dep[k] + 2;
					}
					st = precondition_test(&cx, M.F, 1, &okm);
					if (st == QAWS_STATUS_OK && okm == 1)
					{
						/* exactly one root in M: record it; descendants inside M are covered */
						qaws_exact_box3 b;
						unsigned int q, dup = 0;
						for (k = 0; k < 3; k++)
						{
							cov[ncov].lo[k] = M.lo[k];
							cov[ncov].hi[k] = M.hi[k];
							cov[ncov].dep[k] = M.dep[k];
						}
						ncov++;
						st = refine(&cx, &M, &tmpL, &tmpR);
						for (k = 0; k < 3; k++)
						{
							b.lo[k] = M.lo[k];
							b.hi[k] = M.hi[k];
							b.dep[k] = M.dep[k];
						}
						for (q = 0; q < count && !dup; q++)
							dup = box_overlap(&out[q], &b);
						if (st == QAWS_STATUS_OK && !dup)
						{
							if (count >= capacity)
								st = QAWS_STATUS_BUFFER_TOO_SMALL;
							else
								out[count++] = b;
						}
					}
				}
				if (st != QAWS_STATUS_OK)
					break;
			}
			{
				/* the children go back on the stack, each in its own slot */
				unsigned int d = shallowest(&B), slot = nstack;
				cs_box L, R;
				L.F = store + 3 * cx.size * (slot + 1);
				R.F = store + 3 * cx.size * (slot + 2);
				st = split_box(&cx, &B, d, &L, &R);
				if (st != QAWS_STATUS_OK)
					break;
				memmove(store + 3 * cx.size * slot, R.F, sizeof(qaws_exact_int) * 3 * cx.size);
				R.F = store + 3 * cx.size * slot;
				stack[nstack++] = R;
				stack[nstack++] = L;
			}
		}
	}
	cs_free(stack);
	cs_free(store);
	cs_free(tstore);
	cs_free(mstore);
	*out_count = count;
	return st;
}

qaws_status qaws_exact_curve_surface_hits(qaws_exact_curve const* curve, qaws_exact_surface const* surface, qaws_exact_curve_surface_hit* out_hits,
	unsigned int capacity, unsigned int* out_count)
{
	unsigned int ks, iu, iv, count = 0;
	qaws_exact_int* F = NULL;
	qaws_exact_box3* boxes = NULL;
	qaws_status st = QAWS_STATUS_OK;
	if (!curve || !surface || !out_count || (!out_hits && capacity) || curve->dimension != 3)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (curve->space_exp2 != surface->space_exp2)
		return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
	F = (qaws_exact_int*)cs_alloc(sizeof(qaws_exact_int) * 3 * CS_MAX_SIZE);
	boxes = (qaws_exact_box3*)cs_alloc(sizeof(qaws_exact_box3) * 64);
	if (!F || !boxes)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (ks = 0; ks < curve->span_count && st == QAWS_STATUS_OK; ks++)
	{
		qaws_exact_span const* cs = &curve->spans[ks];
		unsigned int n[3], size;
		n[0] = surface->p;
		n[1] = surface->q;
		n[2] = cs->degree;
		size = (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
		if (size > CS_MAX_SIZE)
		{
			st = QAWS_STATUS_EXACT_UNSUPPORTED;
			break;
		}
		for (iu = 0; iu < surface->nu && st == QAWS_STATUS_OK; iu++)
			for (iv = 0; iv < surface->nv && st == QAWS_STATUS_OK; iv++)
			{
				qaws_exact_int const* h = surface->patch[iu * surface->nv + iv];
				unsigned int a, b, k, c, nb = 0, q;
				/* F_c[a][b][k] = X_c[a][b] w_k - C_c[k] W[a][b] */
				for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
					for (a = 0; a <= n[0] && st == QAWS_STATUS_OK; a++)
						for (b = 0; b <= n[1] && st == QAWS_STATUS_OK; b++)
							for (k = 0; k <= n[2] && st == QAWS_STATUS_OK; k++)
							{
								qaws_exact_int t1, t2;
								unsigned int ab = a * (n[1] + 1) + b;
								st = qaws_exact_int_mul(&t1, &h[ab * 4 + c], &cs->h[k * 4 + 3]);
								if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t2, &cs->h[k * 4 + c], &h[ab * 4 + 3]);
								if (st == QAWS_STATUS_OK) st = qaws_exact_int_sub(&F[c * size + ab * (n[2] + 1) + k], &t1, &t2);
							}
				if (st == QAWS_STATUS_OK)
					st = qaws_exact_solve3(n, F, boxes, 64, &nb, CS_MAX_BOXES);
				for (q = 0; q < nb && st == QAWS_STATUS_OK; q++)
				{
					qaws_exact_curve_surface_hit hit = { 0, 0, 0, 0, 0, 0 };
					unsigned int r, dup = 0;
					st = local_to_param(surface->ub[iu], surface->ub[iu + 1], surface->u_shift, boxes[q].lo[0], boxes[q].hi[0], boxes[q].dep[0], &hit.u_lo,
						&hit.u_hi);
					if (st == QAWS_STATUS_OK)
						st = local_to_param(surface->vb[iv], surface->vb[iv + 1], surface->v_shift, boxes[q].lo[1], boxes[q].hi[1], boxes[q].dep[1],
							&hit.v_lo, &hit.v_hi);
					if (st == QAWS_STATUS_OK)
						st = local_to_param(cs->a, cs->b, curve->param_shift, boxes[q].lo[2], boxes[q].hi[2], boxes[q].dep[2], &hit.t_lo, &hit.t_hi);
					if (st != QAWS_STATUS_OK)
						break;
					/* a root on a patch or span edge is found by both: once */
					for (r = 0; r < count && !dup; r++)
						dup = overlap(out_hits[r].t_lo, out_hits[r].t_hi, hit.t_lo, hit.t_hi) && overlap(out_hits[r].u_lo, out_hits[r].u_hi, hit.u_lo, hit.u_hi) &&
						      overlap(out_hits[r].v_lo, out_hits[r].v_hi, hit.v_lo, hit.v_hi);
					if (dup)
						continue;
					if (count >= capacity)
					{
						st = QAWS_STATUS_BUFFER_TOO_SMALL;
						break;
					}
					out_hits[count++] = hit;
				}
			}
	}
	cs_free(F);
	cs_free(boxes);
	*out_count = count;
	return st;
}
