#ifndef QAWS_INTERNAL_WIDE_H
#define QAWS_INTERNAL_WIDE_H

/*
 * Fixed-width signed integers for exact predicates on int64 coordinates:
 * qaws_wide holds 384 bits in six 64-bit limbs (two's complement), enough
 * for products of a 64-bit coordinate, two 128-bit cross products and their
 * comparisons. Header-only; the 64 x 64 -> 128 multiply uses the compiler's
 * intrinsic where there is one.
 */

#include "../qaws_platform.h"
#include <stdint.h>
#include <string.h>
#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
#endif

#define QAWS_WIDE_LIMBS 6

typedef struct qaws_wide
{
	uint64_t l[QAWS_WIDE_LIMBS];   /* little endian, two's complement */
} qaws_wide;

QAWS_INLINE void qaws_mul_u64(uint64_t* hi, uint64_t* lo, uint64_t a, uint64_t b)
{
#if defined(_MSC_VER) && defined(_M_X64)
	*lo = _umul128(a, b, hi);
#elif defined(__SIZEOF_INT128__)
	__extension__ unsigned __int128 p = (unsigned __int128)a * b;
	*lo = (uint64_t)p;
	*hi = (uint64_t)(p >> 64);
#else
	uint64_t a0 = a & 0xFFFFFFFFu, a1 = a >> 32, b0 = b & 0xFFFFFFFFu, b1 = b >> 32;
	uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
	uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFFu) + (p10 & 0xFFFFFFFFu);
	*lo = (p00 & 0xFFFFFFFFu) | (mid << 32);
	*hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
#endif
}

QAWS_INLINE qaws_wide qaws_wide_from_i64(int64_t v)
{
	qaws_wide w;
	int i;
	w.l[0] = (uint64_t)v;
	for (i = 1; i < QAWS_WIDE_LIMBS; i++)
		w.l[i] = v < 0 ? ~(uint64_t)0 : 0;
	return w;
}

QAWS_INLINE int qaws_wide_sign(qaws_wide const* a)
{
	int i;
	if ((int64_t)a->l[QAWS_WIDE_LIMBS - 1] < 0)
		return -1;
	for (i = 0; i < QAWS_WIDE_LIMBS; i++)
		if (a->l[i])
			return 1;
	return 0;
}

QAWS_INLINE qaws_wide qaws_wide_add(qaws_wide a, qaws_wide const* b)
{
	uint64_t carry = 0;
	int i;
	for (i = 0; i < QAWS_WIDE_LIMBS; i++)
	{
		uint64_t s = a.l[i] + carry;
		uint64_t c1 = s < carry;
		uint64_t t = s + b->l[i];
		carry = c1 + (t < s);
		a.l[i] = t;
	}
	return a;
}

QAWS_INLINE qaws_wide qaws_wide_neg(qaws_wide a)
{
	uint64_t carry = 1;
	int i;
	for (i = 0; i < QAWS_WIDE_LIMBS; i++)
	{
		uint64_t t = ~a.l[i] + carry;
		carry = (carry && t == 0) ? 1 : 0;
		a.l[i] = t;
	}
	return a;
}

QAWS_INLINE qaws_wide qaws_wide_sub(qaws_wide a, qaws_wide const* b)
{
	uint64_t borrow = 0;
	int i;
	for (i = 0; i < QAWS_WIDE_LIMBS; i++)
	{
		uint64_t x = a.l[i], y = b->l[i];
		uint64_t d = x - y;
		uint64_t b1 = x < y;
		uint64_t d2 = d - borrow;
		borrow = b1 | (d < borrow);
		a.l[i] = d2;
	}
	return a;
}

/* a * b, truncated to 384 bits (callers stay inside: the largest product
   in the int64 engine is about 2^375) */
QAWS_INLINE qaws_wide qaws_wide_mul(qaws_wide const* a, qaws_wide const* b)
{
	qaws_wide x = *a, y = *b, r;
	int neg = 0, i, j, nx, ny;
	if ((int64_t)x.l[QAWS_WIDE_LIMBS - 1] < 0) { x = qaws_wide_neg(x); neg = !neg; }
	if ((int64_t)y.l[QAWS_WIDE_LIMBS - 1] < 0) { y = qaws_wide_neg(y); neg = !neg; }
	memset(&r, 0, sizeof(r));
	/* only the limbs up to the highest non-zero one of each (most operands
	   are one or two limbs), then the carry */
	nx = QAWS_WIDE_LIMBS;
	while (nx > 0 && !x.l[nx - 1]) nx--;
	ny = QAWS_WIDE_LIMBS;
	while (ny > 0 && !y.l[ny - 1]) ny--;
	for (i = 0; i < nx; i++)
	{
		uint64_t carry = 0;
		if (!x.l[i]) continue;
		for (j = 0; j < ny && i + j < QAWS_WIDE_LIMBS; j++)
		{
			uint64_t hi, lo, s;
			qaws_mul_u64(&hi, &lo, x.l[i], y.l[j]);
			s = r.l[i + j] + lo;
			hi += s < lo;
			s += carry;
			hi += s < carry;
			r.l[i + j] = s;
			carry = hi;
		}
		for (j = i + ny; carry && j < QAWS_WIDE_LIMBS; j++)
		{
			uint64_t s = r.l[j] + carry;
			carry = s < carry;
			r.l[j] = s;
		}
	}
	return neg ? qaws_wide_neg(r) : r;
}

QAWS_INLINE int qaws_wide_cmp(qaws_wide const* a, qaws_wide const* b)
{
	int i;
	int64_t ha = (int64_t)a->l[QAWS_WIDE_LIMBS - 1], hb = (int64_t)b->l[QAWS_WIDE_LIMBS - 1];
	if (ha != hb)
		return ha < hb ? -1 : 1;
	for (i = QAWS_WIDE_LIMBS - 2; i >= 0; i--)
		if (a->l[i] != b->l[i])
			return a->l[i] < b->l[i] ? -1 : 1;
	return 0;
}

/* the signed 128-bit product of two int64, sign extended to the wide */
QAWS_INLINE void qaws_mul_i64(uint64_t* hi, uint64_t* lo, int64_t a, int64_t b)
{
#if defined(_MSC_VER) && defined(_M_X64)
	int64_t h;
	*lo = (uint64_t)_mul128(a, b, &h);
	*hi = (uint64_t)h;
#elif defined(__SIZEOF_INT128__)
	__extension__ __int128 p = (__int128)a * b;
	*lo = (uint64_t)p;
	*hi = (uint64_t)(p >> 64);
#else
	/* unsigned product, then the two's complement corrections */
	uint64_t h, l;
	qaws_mul_u64(&h, &l, (uint64_t)a, (uint64_t)b);
	if (a < 0) h -= (uint64_t)b;
	if (b < 0) h -= (uint64_t)a;
	*hi = h;
	*lo = l;
#endif
}

/* p + s q for 128-bit signed (hi, lo) pairs, s = +1 or -1, into a wide (129 bits at most) */
QAWS_INLINE qaws_wide qaws_wide_from_i128_pair(uint64_t ph, uint64_t pl, uint64_t qh, uint64_t ql, int s)
{
	qaws_wide w;
	uint64_t lo, hi, c, ext_p = (int64_t)ph < 0 ? ~(uint64_t)0 : 0, ext_q = (int64_t)qh < 0 ? ~(uint64_t)0 : 0, top;
	int i;
	if (s > 0)
	{
		lo = pl + ql;
		c = lo < pl;
		hi = ph + qh + c;
		/* carry into the third limb: from the sign extensions */
		top = ext_p + ext_q + ((hi < ph) || (c && hi == ph) ? 1 : 0);
	}
	else
	{
		uint64_t b;
		lo = pl - ql;
		b = pl < ql;
		hi = ph - qh - b;
		top = ext_p - ext_q - ((ph < qh) || (b && ph == qh) ? 1 : 0);
	}
	w.l[0] = lo;
	w.l[1] = hi;
	w.l[2] = top;
	for (i = 3; i < QAWS_WIDE_LIMBS; i++)
		w.l[i] = (int64_t)top < 0 ? ~(uint64_t)0 : 0;
	return w;
}

/* a x b - c x d style products of int64 values, exactly */
QAWS_INLINE qaws_wide qaws_wide_cross(int64_t ax, int64_t ay, int64_t bx, int64_t by)
{
	uint64_t ph, pl, qh, ql;
	qaws_mul_i64(&ph, &pl, ax, by);
	qaws_mul_i64(&qh, &ql, ay, bx);
	return qaws_wide_from_i128_pair(ph, pl, qh, ql, -1);
}

QAWS_INLINE qaws_wide qaws_wide_dot(int64_t ax, int64_t ay, int64_t bx, int64_t by)
{
	uint64_t ph, pl, qh, ql;
	qaws_mul_i64(&ph, &pl, ax, bx);
	qaws_mul_i64(&qh, &ql, ay, by);
	return qaws_wide_from_i128_pair(ph, pl, qh, ql, 1);
}

QAWS_INLINE double qaws_wide_to_double(qaws_wide const* a)
{
	qaws_wide x = *a;
	double r = 0.0, scale = 1.0;
	int neg = qaws_wide_sign(&x) < 0, i;
	if (neg) x = qaws_wide_neg(x);
	for (i = 0; i < QAWS_WIDE_LIMBS; i++)
	{
		r += (double)x.l[i] * scale;
		scale *= 18446744073709551616.0;
	}
	return neg ? -r : r;
}


#endif /* QAWS_INTERNAL_WIDE_H */
