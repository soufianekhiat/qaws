#include "qaws_exact_roots.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

/* Divides every coefficient by their common gcd (keeps every sign). */
static qaws_status normalize(qaws_exact_int* c, unsigned int n)
{
	qaws_exact_int g;
	unsigned int i;
	qaws_status st;
	qaws_exact_int_zero(&g);
	for (i = 0; i <= n; i++)
		qaws_exact_int_gcd(&g, &g, &c[i]);
	if (qaws_exact_int_bits(&g) <= 1)
		return QAWS_STATUS_OK;
	for (i = 0; i <= n; i++)
		TRY(qaws_exact_int_divmod(&c[i], NULL, &c[i], &g));
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
static qaws_status split_half(qaws_exact_int const* c, unsigned int n, qaws_exact_int* left, qaws_exact_int* right)
{
	qaws_status st;
	TRY(split_half_raw(c, n, left, right));
	TRY(normalize(left, n));
	return normalize(right, n);
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

/* Coefficients of p on (index, index + 1) / 2^depth, by halvings along the bits of index. */
qaws_status qaws_exact_bernstein_restrict(qaws_exact_int const* b, unsigned int n, uint64_t index, int depth, qaws_exact_int* out)
{
	qaws_exact_int left[QAWS_EXACT_ROOTS_MAX_DEGREE + 1], right[QAWS_EXACT_ROOTS_MAX_DEGREE + 1];
	unsigned int i;
	int k;
	qaws_status st;
	for (i = 0; i <= n; i++)
		out[i] = b[i];
	for (k = depth - 1; k >= 0; k--)
	{
		TRY(split_half(out, n, left, right));
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
	unsigned int i;
	int k;
	qaws_status st;
	for (i = 0; i <= n; i++)
	{
		out1[i] = b1[i];
		out2[i] = b2[i];
	}
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
		st = split_half(c, n, left, right);
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
	while (root->depth < depth)
	{
		int sm;
		TRY(split_half(c, n, left, right));
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
