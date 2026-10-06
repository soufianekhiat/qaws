#include "qaws_exact_curve.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <string.h>

/*
 * Winding number of a closed loop around p by the crossings of the ray
 * x > p_x. For a piece with homogeneous Bernstein coefficients (X_i, Y_i,
 * W_i), W > 0, the polynomials
 *   f(t) = Y(t) - p_y W(t),  g(t) = X(t) - p_x W(t)
 * have exact integer Bernstein coefficients (p is scaled onto the integer
 * lattice by a common power of two). On an interval:
 *   - every f coefficient of one strict sign: f has no root, no crossing;
 *   - every g coefficient > 0: the curve stays right of p, so the signed
 *     crossings of y = p_y are up(f(end)) - up(f(start)), up(v) = v > 0
 *     (zero counts as below: a consistent half-open convention);
 *   - every g coefficient <= 0 and f(end), f(start) decided: the piece
 *     stays left of or on x = p_x... handled by subdividing until one of
 *     the two cases above, or until f and g vanish together (p on the
 *     curve) - reported, never guessed.
 * Subdivision at 1/2 is exact: de Casteljau with every level doubled, then
 * the common power of two is removed.
 */

#define QAWS_WIND_MAX_DEPTH 200

typedef struct wind_poly
{
	unsigned int n;
	qaws_exact_int c[2][QAWS_EXACT_MAX_DEGREE + 1];   /* f, g */
} wind_poly;

static int all_sign(qaws_exact_int const* c, unsigned int n, int s)
{
	unsigned int i;
	for (i = 0; i <= n; i++)
		if (qaws_exact_int_sign(&c[i]) != s)
			return 0;
	return 1;
}

/* Halves at 1/2 into left and right, both scaled by 2^n, then reduced. */
static qaws_status split_half(wind_poly const* p, wind_poly* left, wind_poly* right)
{
	qaws_exact_int tri[QAWS_EXACT_MAX_DEGREE + 1];
	unsigned int n = p->n, k, r, i;
	qaws_status st;
	left->n = right->n = n;
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i <= n; i++)
			tri[i] = p->c[k][i];
		/* level r: tri[i] = tri[i] + tri[i + 1] (each level doubles the scale) */
		left->c[k][0] = tri[0];
		right->c[k][n] = tri[n];
		for (r = 1; r <= n; r++)
		{
			for (i = 0; i + r <= n; i++)
			{
				st = qaws_exact_int_add(&tri[i], &tri[i], &tri[i + 1]);
				if (st != QAWS_STATUS_OK)
					return st;
			}
			/* left_r = sum at level r of the first entry, scaled 2^r; bring to 2^n */
			st = qaws_exact_int_shl(&left->c[k][r], &tri[0], n - r);
			if (st != QAWS_STATUS_OK)
				return st;
			st = qaws_exact_int_shl(&right->c[k][n - r], &tri[n - r], n - r);
			if (st != QAWS_STATUS_OK)
				return st;
		}
		st = qaws_exact_int_shl(&left->c[k][0], &left->c[k][0], n);
		if (st != QAWS_STATUS_OK)
			return st;
		st = qaws_exact_int_shl(&right->c[k][n], &right->c[k][n], n);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	/* remove the common power of two of each half (a positive scale) */
	for (k = 0; k < 2; k++)
	{
		wind_poly* h = k ? right : left;
		unsigned int shift = 1000, j;
		for (j = 0; j < 2; j++)
			for (i = 0; i <= n; i++)
			{
				qaws_exact_int const* x = &h->c[j][i];
				unsigned int z = 0;
				int l;
				if (x->sign == 0)
					continue;
				for (l = 0; l < x->size && x->limb[l] == 0; l++)
					z += 32;
				{
					uint32_t v = x->limb[l];
					while (!(v & 1u))
					{
						v >>= 1;
						z++;
					}
				}
				if (z < shift)
					shift = z;
			}
		if (shift > 0 && shift < 1000)
			for (j = 0; j < 2; j++)
				for (i = 0; i <= n; i++)
					qaws_exact_int_shr(&h->c[j][i], &h->c[j][i], shift);
	}
	return QAWS_STATUS_OK;
}

static int up(qaws_exact_int const* v)
{
	return qaws_exact_int_sign(v) > 0;
}

typedef struct wind_item
{
	wind_poly p;
	unsigned int depth;
} wind_item;

/* Signed crossings of the ray on this piece; CERTIFICATION_FAILED past the
   depth. Depth-first over a heap stack (a branch holds at most one pending
   sibling per level). */
static qaws_status crossings(wind_poly const* p0, int* out)
{
	wind_item* stack = (wind_item*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(wind_item) * (QAWS_WIND_MAX_DEPTH + 2)));
	unsigned int top = 0;
	qaws_status st = QAWS_STATUS_OK;
	*out = 0;
	if (!stack)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	stack[top].p = *p0;
	stack[top++].depth = 0;
	while (top > 0 && st == QAWS_STATUS_OK)
	{
		wind_item* it = &stack[--top];
		wind_poly const* p = &it->p;
		unsigned int n = p->n, depth = it->depth;
		if (all_sign(p->c[0], n, 1) || all_sign(p->c[0], n, -1))
			continue;
		if (all_sign(p->c[1], n, 1))
		{
			*out += up(&p->c[0][n]) - up(&p->c[0][0]);
			continue;
		}
		/* entirely left of p (g < 0 everywhere): no crossing of the ray */
		if (all_sign(p->c[1], n, -1))
			continue;
		if (depth >= QAWS_WIND_MAX_DEPTH)
		{
			st = QAWS_STATUS_CERTIFICATION_FAILED;
			break;
		}
		{
			/* the halves replace the item: right below, left on top */
			wind_poly whole = *p;
			st = split_half(&whole, &stack[top + 1].p, &stack[top].p);
			stack[top].depth = depth + 1;
			stack[top + 1].depth = depth + 1;
			top += 2;
		}
	}
	qaws_internal_dealloc(NULL, stack);
	return st;
}

/* Homogeneous endpoints equal: X1 W2 = X2 W1, Y1 W2 = Y2 W1. */
static int same_point(qaws_exact_int const* a, qaws_exact_int const* b)
{
	qaws_exact_int l, r;
	int c;
	for (c = 0; c < 2; c++)
	{
		if (qaws_exact_int_mul(&l, &a[c], &b[2]) != QAWS_STATUS_OK || qaws_exact_int_mul(&r, &b[c], &a[2]) != QAWS_STATUS_OK)
			return 0;
		if (qaws_exact_int_cmp(&l, &r) != 0)
			return 0;
	}
	return 1;
}

qaws_status qaws_exact_winding_2d(qaws_exact_curve const* const* pieces, unsigned int count, double const p[2], int* out_winding)
{
	int64_t m[2];
	int e[2], c, kscale = 0, total = 0;
	unsigned int k, i;
	qaws_exact_int P[2];
	qaws_status st;
	if (!pieces || count == 0 || !p || !out_winding)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (k = 0; k < count; k++)
	{
		qaws_exact_curve const* a;
		qaws_exact_curve const* b;
		qaws_exact_span const* last;
		if (!pieces[k] || pieces[k]->dimension != 2 || pieces[k]->span_count == 0)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (pieces[k]->space_exp2 != pieces[0]->space_exp2)
			return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
		/* closed: the end of piece k is the start of piece k + 1 */
		a = pieces[k];
		b = pieces[(k + 1) % count];
		last = &a->spans[a->span_count - 1];
		if (!b->spans || b->span_count == 0 || !same_point(&last->h[last->degree * 3], &b->spans[0].h[0]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	/* p in lattice units m 2^(e - s); a common 2^kscale makes it integral */
	for (c = 0; c < 2; c++)
	{
		if (!qaws_exact_split_double(p[c], &m[c], &e[c]))
			return QAWS_STATUS_INVALID_ARGUMENT;
		e[c] -= pieces[0]->space_exp2;
		if (m[c] != 0 && -e[c] > kscale)
			kscale = -e[c];
	}
	if (kscale > 1500)
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	{
		qaws_exact_int H[3];
		for (c = 0; c < 2; c++)
		{
			qaws_exact_int_from_i64(&P[c], m[c]);
			if (m[c] != 0)
			{
				st = qaws_exact_int_shl(&P[c], &P[c], (unsigned int)(e[c] + kscale));
				if (st != QAWS_STATUS_OK)
					return st;
			}
			H[c] = P[c];
		}
		qaws_exact_int_from_i64(&H[2], 1);
		st = qaws_exact_int_shl(&H[2], &H[2], (unsigned int)kscale);
		if (st != QAWS_STATUS_OK)
			return st;
		(void)total;
		(void)i;
		return qaws_exact_winding_2d_hom(pieces, count, H, out_winding);
	}
}

qaws_status qaws_exact_winding_2d_hom(qaws_exact_curve const* const* pieces, unsigned int count, qaws_exact_int const* P, int* out_winding)
{
	unsigned int k, i;
	int total = 0;
	qaws_status st;
	if (qaws_exact_int_sign(&P[2]) <= 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* every span of every piece is one integer homogeneous Bezier */
	for (k = 0; k < count; k++)
	{
		unsigned int s;
		for (s = 0; s < pieces[k]->span_count; s++)
		{
			qaws_exact_span const* sp = &pieces[k]->spans[s];
			wind_poly w;
			int cr;
			w.n = sp->degree;
			for (i = 0; i <= sp->degree; i++)
			{
				qaws_exact_int x, y, t;
				/* f_i = Y_i P_w - P_y W_i, g_i = X_i P_w - P_x W_i (P_w > 0) */
				st = qaws_exact_int_mul(&y, &sp->h[i * 3 + 1], &P[2]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &P[1], &sp->h[i * 3 + 2]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_sub(&w.c[0][i], &y, &t);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&x, &sp->h[i * 3], &P[2]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &P[0], &sp->h[i * 3 + 2]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_sub(&w.c[1][i], &x, &t);
				if (st != QAWS_STATUS_OK)
					return st;
			}
			st = crossings(&w, &cr);
			if (st != QAWS_STATUS_OK)
				return st;
			total += cr;
		}
	}
	*out_winding = total;
	return QAWS_STATUS_OK;
}
