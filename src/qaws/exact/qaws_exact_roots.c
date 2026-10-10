#include "qaws_exact_roots.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

/* Divides every coefficient by their common gcd (keeps every sign). */
static qaws_status rt_normalize(qaws_exact_int* c, unsigned int n)
{
	/* halving introduces only powers of two: strip the common one (a shift, no gcd) */
	unsigned int i, k = ~0u;
	for (i = 0; i <= n; i++)
		if (!qaws_exact_int_is_zero(&c[i]))
		{
			unsigned int z = qaws_exact_int_ctz(&c[i]);
			if (z < k)
				k = z;
			if (k == 0)
				return QAWS_STATUS_OK;
		}
	if (k == ~0u)
		return QAWS_STATUS_OK;
	for (i = 0; i <= n; i++)
		qaws_exact_int_shr(&c[i], &c[i], k);
	return QAWS_STATUS_OK;
}

/*
 * Halves at s = 1/2 without division: pairwise sums give 2^r times the
 * De Casteljau points, so left_i = d_0^i 2^(n-i) and right_j = d_j^(n-j) 2^j
 * share the scale 2^n (a positive factor: signs and roots are kept).
 */
static qaws_status split_half_raw(qaws_exact_int const* c, unsigned int n, qaws_exact_int* left, qaws_exact_int* right)
{
	qaws_exact_int d[QAWS_EXACT_ROOTS_MAX_DEGREE + 1];
	unsigned int r, i;
	qaws_status st;
	for (i = 0; i <= n; i++)
		d[i] = c[i];
	left[0] = d[0];
	right[n] = d[n];
	for (r = 1; r <= n; r++)
	{
		for (i = 0; i + r <= n; i++)
			TRY(qaws_exact_int_add(&d[i], &d[i], &d[i + 1]));
		TRY(qaws_exact_int_shl(&left[r], &d[0], n - r));
		TRY(qaws_exact_int_shl(&right[n - r], &d[n - r], n - r));
	}
	TRY(qaws_exact_int_shl(&left[0], &left[0], n));
	TRY(qaws_exact_int_shl(&right[n], &right[n], n));
	return QAWS_STATUS_OK;
}

/* The halves, each divided by its own gcd. */
static qaws_status rt_split_half(qaws_exact_int const* c, unsigned int n, qaws_exact_int* left, qaws_exact_int* right)
{
	qaws_status st;
	TRY(split_half_raw(c, n, left, right));
	TRY(rt_normalize(left, n));
	return rt_normalize(right, n);
}

/* Sign variations, zeros ignored. */
static unsigned int variations(qaws_exact_int const* c, unsigned int n)
{
	unsigned int i, v = 0;
	int last = 0;
	for (i = 0; i <= n; i++)
	{
		int s = qaws_exact_int_sign(&c[i]);
		if (s == 0)
			continue;
		if (last != 0 && s != last)
			v++;
		last = s;
	}
	return v;
}

/*
 * c restricted to [x, 1] (keep_right) or [0, x], x = a / b (0 < a < b, b <= 2^62):
 * integer De Casteljau with the weights (b - a, a), in place; every point
 * comes out times b^n.
 */
static qaws_status split_ratio(qaws_exact_int* c, unsigned int n, int64_t a, int64_t b, int keep_right)
{
	qaws_exact_int d[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], out[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], t1, t2;
	unsigned int r, i, k;
	qaws_status st;
	for (i = 0; i <= n; i++)
		d[i] = c[i];
	if (keep_right)
		out[n] = d[n];
	else
		out[0] = d[0];
	for (r = 1; r <= n; r++)
	{
		for (i = 0; i + r <= n; i++)
		{
			TRY(qaws_exact_int_mul_i64(&t1, &d[i], b - a));
			TRY(qaws_exact_int_mul_i64(&t2, &d[i + 1], a));
			TRY(qaws_exact_int_add(&d[i], &t1, &t2));
		}
		if (keep_right)
			out[n - r] = d[n - r];
		else
			out[r] = d[0];
	}
	for (i = 0; i <= n; i++)
	{
		unsigned int lev = keep_right ? n - i : i;
		for (k = lev; k < n; k++)
			TRY(qaws_exact_int_mul_i64(&out[i], &out[i], b));
		c[i] = out[i];
	}
	return QAWS_STATUS_OK;
}

/* Divides two arrays by the gcd of all their entries (a positive factor). */
static qaws_status normalize_two(qaws_exact_int* c1, qaws_exact_int* c2, unsigned int n)
{
	qaws_exact_int g;
	unsigned int i;
	qaws_status st;
	qaws_exact_int_zero(&g);
	for (i = 0; i <= n; i++)
	{
		qaws_exact_int_gcd(&g, &g, &c1[i]);
		if (c2)
			qaws_exact_int_gcd(&g, &g, &c2[i]);
	}
	if (qaws_exact_int_bits(&g) <= 1)
		return QAWS_STATUS_OK;
	for (i = 0; i <= n; i++)
	{
		TRY(qaws_exact_int_divmod(&c1[i], NULL, &c1[i], &g));
		if (c2)
			TRY(qaws_exact_int_divmod(&c2[i], NULL, &c2[i], &g));
	}
	return QAWS_STATUS_OK;
}

/*
 * Direct restriction to [index, index + 1] / 2^depth: two rational cuts
 * instead of depth halvings (when the intermediate integers stay well in
 * the budget), then one gcd.
 */
static qaws_status restrict_direct(qaws_exact_int* c1, qaws_exact_int* c2, unsigned int n, uint64_t index, int depth)
{
	int64_t den = (int64_t)1 << depth, rest = den - (int64_t)index;
	qaws_status st;
	if (index > 0)
	{
		TRY(split_ratio(c1, n, (int64_t)index, den, 1));   /* [index / 2^depth, 1] */
		if (c2) TRY(split_ratio(c2, n, (int64_t)index, den, 1));
	}
	if (rest > 1)
	{
		TRY(split_ratio(c1, n, 1, rest, 0));   /* its first 1 / rest */
		if (c2) TRY(split_ratio(c2, n, 1, rest, 0));
	}
	return normalize_two(c1, c2, n);
}

/* Coefficients of p on (index, index + 1) / 2^depth. */
qaws_status qaws_exact_bernstein_restrict(qaws_exact_int const* b, unsigned int n, uint64_t index, int depth, qaws_exact_int* out)
{
	qaws_exact_int left[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], right[QAWS_EXACT_ROOTS_MAX_DEGREE + 1];
	unsigned int i, maxb = 0;
	int k;
	qaws_status st;
	for (i = 0; i <= n; i++)
	{
		out[i] = b[i];
		if (qaws_exact_int_bits(&b[i]) > maxb)
			maxb = qaws_exact_int_bits(&b[i]);
	}
	if (depth > 0 && depth <= 62 && maxb + 2 * n * (unsigned int)depth + 64 < QAWS_EXACT_MAX_BITS)
		return restrict_direct(out, NULL, n, index, depth);
	/* by halvings along the bits of index (each normalized: the integers stay small) */
	for (k = depth - 1; k >= 0; k--)
	{
		TRY(rt_split_half(out, n, left, right));
		for (i = 0; i <= n; i++)
			out[i] = ((index >> k) & 1) ? right[i] : left[i];
	}
	return QAWS_STATUS_OK;
}


/* Both halves of b1 and b2 at once, divided by one common gcd. */
qaws_status qaws_exact_bernstein_restrict_pair(qaws_exact_int const* b1, qaws_exact_int const* b2, unsigned int n, uint64_t index, int depth,
	qaws_exact_int* out1, qaws_exact_int* out2)
{
	qaws_exact_int left[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], right[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], g;
	unsigned int i, maxb = 0;
	int k;
	qaws_status st;
	for (i = 0; i <= n; i++)
	{
		out1[i] = b1[i];
		out2[i] = b2[i];
		if (qaws_exact_int_bits(&b1[i]) > maxb) maxb = qaws_exact_int_bits(&b1[i]);
		if (qaws_exact_int_bits(&b2[i]) > maxb) maxb = qaws_exact_int_bits(&b2[i]);
	}
	if (depth > 0 && depth <= 62 && maxb + 2 * n * (unsigned int)depth + 64 < QAWS_EXACT_MAX_BITS)
		return restrict_direct(out1, out2, n, index, depth);
	for (k = depth - 1; k >= 0; k--)
	{
		int right_half = (int)((index >> k) & 1);
		TRY(split_half_raw(out1, n, left, right));
		for (i = 0; i <= n; i++)
			out1[i] = right_half ? right[i] : left[i];
		TRY(split_half_raw(out2, n, left, right));
		for (i = 0; i <= n; i++)
			out2[i] = right_half ? right[i] : left[i];
		qaws_exact_int_zero(&g);
		for (i = 0; i <= n; i++)
		{
			qaws_exact_int_gcd(&g, &g, &out1[i]);
			qaws_exact_int_gcd(&g, &g, &out2[i]);
		}
		if (qaws_exact_int_bits(&g) > 1)
			for (i = 0; i <= n; i++)
			{
				TRY(qaws_exact_int_divmod(&out1[i], NULL, &out1[i], &g));
				TRY(qaws_exact_int_divmod(&out2[i], NULL, &out2[i], &g));
			}
	}
	return QAWS_STATUS_OK;
}
typedef struct root_item
{
	uint64_t index;
	int depth;
} root_item;

static int root_less(qaws_exact_root const* a, qaws_exact_root const* b)
{
	int d = a->depth > b->depth ? a->depth : b->depth;
	uint64_t x = a->index << (d - a->depth), y = b->index << (d - b->depth);
	if (x != y)
		return x < y;
	return a->exact && !b->exact;   /* a point before the interval starting at it */
}

static qaws_status push_root(qaws_exact_root* out, unsigned int capacity, unsigned int* count, uint64_t index, int depth, int exact)
{
	if (*count >= capacity)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	out[*count].index = index;
	out[*count].depth = depth;
	out[*count].exact = exact;
	(*count)++;
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_bernstein_isolate(qaws_exact_int const* b, unsigned int n, qaws_exact_root* out, unsigned int capacity,
	unsigned int* out_count)
{
	qaws_exact_int* stack = NULL;
	root_item* items = NULL;
	qaws_exact_int left[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], right[QAWS_EXACT_ROOTS_MAX_DEGREE + 1];
	unsigned int top = 0, i, count = 0, max_items = 2 * QAWS_EXACT_ROOTS_MAX_DEPTH + 4;
	qaws_status st = QAWS_STATUS_OK;
	if (!b || !out_count || n > QAWS_EXACT_ROOTS_MAX_DEGREE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	for (i = 0; i <= n && qaws_exact_int_is_zero(&b[i]); i++)
		;
	if (i > n)
		return QAWS_STATUS_CERTIFICATION_FAILED;   /* p == 0 */
	/* the domain ends */
	if (qaws_exact_int_is_zero(&b[0]))
		TRY(push_root(out, capacity, &count, 0, 0, 1));
	if (qaws_exact_int_is_zero(&b[n]))
		TRY(push_root(out, capacity, &count, 1, 0, 1));
	stack = (qaws_exact_int*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_int) * (n + 1) * max_items));
	items = (root_item*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(root_item) * max_items));
	if (!stack || !items)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK)
	{
		for (i = 0; i <= n; i++)
			stack[i] = b[i];
		items[0].index = 0;
		items[0].depth = 0;
		top = 1;
	}
	while (st == QAWS_STATUS_OK && top > 0)
	{
		qaws_exact_int* c;
		root_item it;
		unsigned int v;
		top--;
		c = &stack[top * (n + 1)];
		it = items[top];
		v = variations(c, n);
		if (v == 0)
			continue;
		if (v == 1)
		{
			st = push_root(out, capacity, &count, it.index, it.depth, 0);
			continue;
		}
		if (it.depth >= QAWS_EXACT_ROOTS_MAX_DEPTH)
		{
			st = QAWS_STATUS_CERTIFICATION_FAILED;
			break;
		}
		st = rt_split_half(c, n, left, right);
		if (st != QAWS_STATUS_OK)
			break;
		/* the midpoint value is right[0] (= left[n] up to the common scale) */
		if (qaws_exact_int_is_zero(&right[0]))
			st = push_root(out, capacity, &count, 2 * it.index + 1, it.depth + 1, 1);
		if (st != QAWS_STATUS_OK || top + 2 > max_items)
		{
			if (st == QAWS_STATUS_OK)
				st = QAWS_STATUS_INTERNAL_ERROR;
			break;
		}
		/* depth-first, left half on top */
		for (i = 0; i <= n; i++)
			stack[top * (n + 1) + i] = right[i];
		items[top].index = 2 * it.index + 1;
		items[top].depth = it.depth + 1;
		top++;
		for (i = 0; i <= n; i++)
			stack[top * (n + 1) + i] = left[i];
		items[top].index = 2 * it.index;
		items[top].depth = it.depth + 1;
		top++;
	}
	qaws_internal_dealloc(NULL, stack);
	qaws_internal_dealloc(NULL, items);
	if (st != QAWS_STATUS_OK)
		return st;
	/* sort (insertion: few roots) */
	for (i = 1; i < count; i++)
	{
		qaws_exact_root key = out[i];
		int j = (int)i - 1;
		while (j >= 0 && root_less(&key, &out[j]))
		{
			out[j + 1] = out[j];
			j--;
		}
		out[j + 1] = key;
	}
	*out_count = count;
	return QAWS_STATUS_OK;
}

/*
 * The sign of c (Bernstein on [0, 1]) at x = j / 2^k, exactly: De Casteljau
 * with the weights (2^k - j, j) (the value times 2^(k n), a positive factor).
 */
static qaws_status sign_at(qaws_exact_int const* c, unsigned int n, uint64_t j, int k, int* out)
{
	qaws_exact_int d[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], t1, t2;
	int64_t den = (int64_t)1 << k, a = (int64_t)j;
	unsigned int r, i;
	qaws_status st;
	for (i = 0; i <= n; i++)
		d[i] = c[i];
	if (a == 0)
	{
		*out = qaws_exact_int_sign(&d[0]);
		return QAWS_STATUS_OK;
	}
	if (a == den)
	{
		*out = qaws_exact_int_sign(&d[n]);
		return QAWS_STATUS_OK;
	}
	for (r = 1; r <= n; r++)
		for (i = 0; i + r <= n; i++)
		{
			TRY(qaws_exact_int_mul_i64(&t1, &d[i], den - a));
			TRY(qaws_exact_int_mul_i64(&t2, &d[i + 1], a));
			TRY(qaws_exact_int_add(&d[i], &t1, &t2));
		}
	*out = qaws_exact_int_sign(&d[0]);
	return QAWS_STATUS_OK;
}

/* Value and derivative of a double Bernstein polynomial at x. */
static void eval_double(double const* c, unsigned int n, double x, double* v, double* dv)
{
	double d[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], e[QAWS_EXACT_ROOTS_MAX_DEGREE + 1];
	unsigned int r, i;
	for (i = 0; i <= n; i++)
		d[i] = c[i];
	for (i = 0; i < n; i++)
		e[i] = (double)n * (c[i + 1] - c[i]);
	for (r = 1; r <= n; r++)
		for (i = 0; i + r <= n; i++)
			d[i] = (1 - x) * d[i] + x * d[i + 1];
	for (r = 1; r + 1 <= n; r++)
		for (i = 0; i + r < n; i++)
			e[i] = (1 - x) * e[i] + x * e[i + 1];
	*v = d[0];
	*dv = n > 0 ? e[0] : 0;
}

/*
 * The root of c in (0, 1) is unique and simple (an isolating interval):
 * sign s0 left of it, -s0 right of it. Newton in doubles gives x; the cell
 * [j, j + 1] / 2^k around it holds the root when the exact signs at its
 * ends are s0 and -s0 (a zero there: the exact root). A few neighbouring
 * cells are tried; *found = 0 leaves the caller's bisection to it.
 */
static qaws_status newton_bracket(qaws_exact_int const* c, unsigned int n, int s0, qaws_exact_root* root, int depth, int* found)
{
	double cd[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], x = 0.5;
	unsigned int i, maxb = 0, it;
	int k = depth - root->depth, tries;
	uint64_t cells, j;
	qaws_status st;
	*found = 0;
	if (k <= 0 || k > 60 || s0 == 0)
		return QAWS_STATUS_OK;
	for (i = 0; i <= n; i++)
		if (qaws_exact_int_bits(&c[i]) > maxb)
			maxb = qaws_exact_int_bits(&c[i]);
	if (maxb + n * (unsigned int)k + 64 >= QAWS_EXACT_MAX_BITS)
		return QAWS_STATUS_OK;
	for (i = 0; i <= n; i++)
	{
		qaws_exact_int t;
		qaws_exact_int_shr(&t, &c[i], maxb > 900 ? maxb - 900 : 0);
		cd[i] = qaws_exact_int_to_double(&t);
	}
	/* Newton, kept inside (0, 1) (bisection steps on the double signs when it leaves) */
	{
		double lo = 0, hi = 1;
		for (it = 0; it < 80; it++)
		{
			double v, dv, nx;
			eval_double(cd, n, x, &v, &dv);
			if (v == 0)
				break;
			if ((v > 0) == (s0 > 0))
				lo = x;
			else
				hi = x;
			nx = dv != 0 ? x - v / dv : 0.5 * (lo + hi);
			if (!(nx > lo && nx < hi))
				nx = 0.5 * (lo + hi);
			if (nx == x)
				break;
			x = nx;
		}
	}
	cells = (uint64_t)1 << k;
	j = (uint64_t)(x * (double)cells);
	if (j >= cells)
		j = cells - 1;
	for (tries = 0; tries < 6; tries++)
	{
		int sl = 0, sr = 0;
		TRY(sign_at(c, n, j, k, &sl));
		TRY(sign_at(c, n, j + 1, k, &sr));
		/* the ends of [0, 1] carry the signs just inside them */
		if (j == 0 && sl == 0) sl = s0;
		if (j + 1 == cells && sr == 0) sr = -s0;
		if (sl == 0 || sr == 0)
		{
			root->index = (root->index << k) + (sl == 0 ? j : j + 1);
			root->depth = depth;
			root->exact = 1;
			/* in lowest terms, as the bisection would have found it */
			while (root->depth > 0 && !(root->index & 1))
			{
				root->index >>= 1;
				root->depth--;
			}
			*found = 1;
			return QAWS_STATUS_OK;
		}
		if (sl == s0 && sr == -s0)
		{
			root->index = (root->index << k) + j;
			root->depth = depth;
			*found = 1;
			return QAWS_STATUS_OK;
		}
		if (sl == -s0 && j > 0)
			j--;   /* the root lies left of the cell */
		else if (sr == s0 && j + 1 < cells)
			j++;
		else
			break;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_bernstein_refine(qaws_exact_int const* b, unsigned int n, qaws_exact_root* root, int depth)
{
	qaws_exact_int c[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], left[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], right[QAWS_EXACT_ROOTS_MAX_DEGREE + 1];
	unsigned int i;
	int s0;
	qaws_status st;
	if (!b || !root || n > QAWS_EXACT_ROOTS_MAX_DEGREE || depth > 62)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (root->exact || root->depth >= depth)
		return QAWS_STATUS_OK;
	TRY(qaws_exact_bernstein_restrict(b, n, root->index, root->depth, c));
	/* the sign just right of the left end: its first non-zero coefficient
	   (the left end may itself be an exact root) */
	s0 = 0;
	for (i = 0; i <= n && s0 == 0; i++)
		s0 = qaws_exact_int_sign(&c[i]);
	{
		/* fast path: a double Newton estimate, then a target-width bracket certified by two exact signs */
		int found = 0;
		TRY(newton_bracket(c, n, s0, root, depth, &found));
		if (found)
			return QAWS_STATUS_OK;
	}
	while (root->depth < depth)
	{
		int sm;
		TRY(rt_split_half(c, n, left, right));
		sm = qaws_exact_int_sign(&right[0]);
		if (sm == 0)
		{
			root->index = 2 * root->index + 1;
			root->depth++;
			root->exact = 1;
			return QAWS_STATUS_OK;
		}
		root->depth++;
		if (sm != s0)
		{
			root->index = 2 * root->index;
			for (i = 0; i <= n; i++)
				c[i] = left[i];
		}
		else
		{
			root->index = 2 * root->index + 1;
			for (i = 0; i <= n; i++)
				c[i] = right[i];
		}
	}
	return QAWS_STATUS_OK;
}
