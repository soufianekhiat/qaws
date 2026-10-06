#include "qaws_exact_int.h"
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Magnitudes                                                        */
/* ------------------------------------------------------------------ */

/* x = y, copying only the used limbs (an assignment would copy all of them). */
static void copy_int(qaws_exact_int* x, qaws_exact_int const* y)
{
	if (x == y)
		return;
	memcpy(x->limb, y->limb, sizeof(uint32_t) * (size_t)y->size);
	x->size = y->size;
	x->sign = y->sign;
}

static void normalize(qaws_exact_int* x)
{
	while (x->size > 0 && x->limb[x->size - 1] == 0)
		x->size--;
	if (x->size == 0)
		x->sign = 0;
}

static int mag_cmp(qaws_exact_int const* a, qaws_exact_int const* b)
{
	int i;
	if (a->size != b->size)
		return a->size < b->size ? -1 : 1;
	for (i = a->size - 1; i >= 0; i--)
		if (a->limb[i] != b->limb[i])
			return a->limb[i] < b->limb[i] ? -1 : 1;
	return 0;
}

/* |r| = |a| + |b|; r may alias a or b. */
static qaws_status mag_add(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b)
{
	int n = a->size > b->size ? a->size : b->size, i;
	uint64_t carry = 0;
	for (i = 0; i < n; i++)
	{
		uint64_t s = carry + (i < a->size ? a->limb[i] : 0u) + (i < b->size ? b->limb[i] : 0u);
		r->limb[i] = (uint32_t)s;
		carry = s >> 32;
	}
	if (carry)
	{
		if (n >= QAWS_EXACT_LIMBS)
			return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
		r->limb[n++] = (uint32_t)carry;
	}
	r->size = n;
	return QAWS_STATUS_OK;
}

/* |r| = |a| - |b| with |a| >= |b|; r may alias a or b. */
static void mag_sub(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b)
{
	int i, n = a->size;
	int64_t borrow = 0;
	for (i = 0; i < n; i++)
	{
		int64_t d = (int64_t)a->limb[i] - (int64_t)(i < b->size ? b->limb[i] : 0u) - borrow;
		borrow = d < 0;
		r->limb[i] = (uint32_t)(d + (borrow ? ((int64_t)1 << 32) : 0));
	}
	r->size = n;
}

/* |a| / d in place, returns the remainder. */
static uint32_t mag_divmod_u32(qaws_exact_int* a, uint32_t d)
{
	uint64_t rem = 0;
	int i;
	for (i = a->size - 1; i >= 0; i--)
	{
		uint64_t cur = (rem << 32) | a->limb[i];
		a->limb[i] = (uint32_t)(cur / d);
		rem = cur % d;
	}
	normalize(a);
	return (uint32_t)rem;
}

/* ------------------------------------------------------------------ */
/*  Construction and inspection                                       */
/* ------------------------------------------------------------------ */

void qaws_exact_int_zero(qaws_exact_int* x)
{
	x->size = 0;
	x->sign = 0;
}

void qaws_exact_int_from_i64(qaws_exact_int* x, int64_t v)
{
	uint64_t m = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
	x->limb[0] = (uint32_t)m;
	x->limb[1] = (uint32_t)(m >> 32);
	x->size = 2;
	x->sign = v < 0 ? -1 : 1;
	normalize(x);
}

int qaws_exact_int_sign(qaws_exact_int const* x)
{
	return x->sign;
}

int qaws_exact_int_is_zero(qaws_exact_int const* x)
{
	return x->sign == 0;
}

int qaws_exact_int_cmp_abs(qaws_exact_int const* a, qaws_exact_int const* b)
{
	return mag_cmp(a, b);
}

int qaws_exact_int_cmp(qaws_exact_int const* a, qaws_exact_int const* b)
{
	if (a->sign != b->sign)
		return a->sign < b->sign ? -1 : 1;
	return a->sign >= 0 ? mag_cmp(a, b) : -mag_cmp(a, b);
}

unsigned int qaws_exact_int_bits(qaws_exact_int const* x)
{
	uint32_t top;
	unsigned int n;
	if (x->size == 0)
		return 0;
	top = x->limb[x->size - 1];
	n = 32u * (unsigned int)(x->size - 1);
	while (top)
	{
		n++;
		top >>= 1;
	}
	return n;
}

/* ------------------------------------------------------------------ */
/*  Arithmetic                                                        */
/* ------------------------------------------------------------------ */

void qaws_exact_int_neg(qaws_exact_int* r, qaws_exact_int const* a)
{
	if (r != a)
		copy_int(r, a);
	r->sign = -r->sign;
}

/* r = a + s b, s = +1 or -1 */
static qaws_status add_signed(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b, int s)
{
	int bs = b->sign * s;
	qaws_status st;
	if (b->sign == 0)
	{
		if (r != a)
			copy_int(r, a);
		return QAWS_STATUS_OK;
	}
	if (a->sign == 0)
	{
		if (r != b)
			copy_int(r, b);
		r->sign = bs;
		return QAWS_STATUS_OK;
	}
	if (a->sign == bs)
	{
		st = mag_add(r, a, b);
		r->sign = bs;
		return st;
	}
	if (mag_cmp(a, b) >= 0)
	{
		int as = a->sign;
		mag_sub(r, a, b);
		r->sign = as;
	}
	else
	{
		mag_sub(r, b, a);
		r->sign = bs;
	}
	normalize(r);
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_int_add(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b)
{
	return add_signed(r, a, b, 1);
}

qaws_status qaws_exact_int_sub(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b)
{
	return add_signed(r, a, b, -1);
}

qaws_status qaws_exact_int_mul(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b)
{
	qaws_exact_int t;
	int i, j, n;
	if (a->sign == 0 || b->sign == 0)
	{
		qaws_exact_int_zero(r);
		return QAWS_STATUS_OK;
	}
	n = a->size + b->size;
	if (n > QAWS_EXACT_LIMBS)
	{
		/* the product may still fit one limb short of the sum of sizes */
		if (qaws_exact_int_bits(a) + qaws_exact_int_bits(b) - 1 > QAWS_EXACT_MAX_BITS)
			return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
		n = QAWS_EXACT_LIMBS;
	}
	memset(t.limb, 0, sizeof(uint32_t) * (size_t)n);
	for (i = 0; i < a->size; i++)
	{
		uint64_t carry = 0, ai = a->limb[i];
		for (j = 0; j < b->size; j++)
		{
			uint64_t cur;
			if (i + j >= n)
			{
				if ((ai && b->limb[j]) || carry)
					return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
				continue;
			}
			cur = ai * b->limb[j] + t.limb[i + j] + carry;
			t.limb[i + j] = (uint32_t)cur;
			carry = cur >> 32;
		}
		for (j = i + b->size; carry; j++)
		{
			uint64_t cur;
			if (j >= n)
				return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
			cur = (uint64_t)t.limb[j] + carry;
			t.limb[j] = (uint32_t)cur;
			carry = cur >> 32;
		}
	}
	t.size = n;
	t.sign = a->sign * b->sign;
	normalize(&t);
	copy_int(r, &t);
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_int_mul_i64(qaws_exact_int* r, qaws_exact_int const* a, int64_t s)
{
	uint64_t m = s < 0 ? (uint64_t)0 - (uint64_t)s : (uint64_t)s;
	int sign = s < 0 ? -a->sign : a->sign;
	if (a->sign == 0 || s == 0)
	{
		qaws_exact_int_zero(r);
		return QAWS_STATUS_OK;
	}
	if (m <= 0xFFFFFFFFu)
	{
		/* one limb: a single pass, in place */
		uint64_t carry = 0;
		int i, n = a->size;
		for (i = 0; i < n; i++)
		{
			uint64_t cur = (uint64_t)a->limb[i] * m + carry;
			r->limb[i] = (uint32_t)cur;
			carry = cur >> 32;
		}
		if (carry)
		{
			if (n >= QAWS_EXACT_LIMBS)
				return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
			r->limb[n++] = (uint32_t)carry;
		}
		r->size = n;
		r->sign = sign;
		return QAWS_STATUS_OK;
	}
	{
		/* two limbs: two passes into a temporary */
		uint32_t t[QAWS_EXACT_LIMBS + 2];
		uint64_t m0 = m & 0xFFFFFFFFu, m1 = m >> 32, carry;
		int i, n = a->size;
		if (n + 2 > QAWS_EXACT_LIMBS + 1)
		{
			qaws_exact_int b;
			qaws_exact_int_from_i64(&b, s);
			return qaws_exact_int_mul(r, a, &b);
		}
		carry = 0;
		for (i = 0; i < n; i++)
		{
			uint64_t cur = (uint64_t)a->limb[i] * m0 + carry;
			t[i] = (uint32_t)cur;
			carry = cur >> 32;
		}
		t[n] = (uint32_t)carry;
		t[n + 1] = 0;
		carry = 0;
		for (i = 0; i < n; i++)
		{
			uint64_t cur = (uint64_t)a->limb[i] * m1 + t[i + 1] + carry;
			t[i + 1] = (uint32_t)cur;
			carry = cur >> 32;
		}
		t[n + 1] = (uint32_t)((uint64_t)t[n + 1] + carry);
		n += 2;
		while (n > 0 && t[n - 1] == 0)
			n--;
		if (n > QAWS_EXACT_LIMBS)
			return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
		memcpy(r->limb, t, sizeof(uint32_t) * (size_t)n);
		r->size = n;
		r->sign = sign;
		return QAWS_STATUS_OK;
	}
}

qaws_status qaws_exact_int_shl(qaws_exact_int* r, qaws_exact_int const* a, unsigned int k)
{
	unsigned int limbs = k / 32, bits = k % 32;
	int i, n;
	if (a->sign == 0)
	{
		qaws_exact_int_zero(r);
		return QAWS_STATUS_OK;
	}
	if (qaws_exact_int_bits(a) + k > QAWS_EXACT_MAX_BITS)
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;
	n = a->size + (int)limbs + 1;
	if (n > QAWS_EXACT_LIMBS)
		n = QAWS_EXACT_LIMBS;
	{
		qaws_exact_int t;
		memset(t.limb, 0, sizeof(uint32_t) * (size_t)n);
		for (i = 0; i < a->size; i++)
		{
			uint64_t v = (uint64_t)a->limb[i] << bits;
			int at = i + (int)limbs;
			t.limb[at] |= (uint32_t)v;
			if (at + 1 < n)
				t.limb[at + 1] |= (uint32_t)(v >> 32);
		}
		t.size = n;
		t.sign = a->sign;
		normalize(&t);
		copy_int(r, &t);
	}
	return QAWS_STATUS_OK;
}

void qaws_exact_int_shr(qaws_exact_int* r, qaws_exact_int const* a, unsigned int k)
{
	unsigned int limbs = k / 32, bits = k % 32;
	qaws_exact_int t;
	int i, n = a->size - (int)limbs;
	if (n <= 0)
	{
		qaws_exact_int_zero(r);
		return;
	}
	for (i = 0; i < n; i++)
	{
		uint64_t lo = a->limb[i + (int)limbs];
		uint64_t hi = i + (int)limbs + 1 < a->size ? a->limb[i + (int)limbs + 1] : 0u;
		t.limb[i] = (uint32_t)(((hi << 32) | lo) >> bits);
	}
	t.size = n;
	t.sign = a->sign;
	normalize(&t);
	copy_int(r, &t);
}

qaws_status qaws_exact_int_divexact_u32(qaws_exact_int* r, qaws_exact_int const* a, uint32_t d)
{
	qaws_exact_int t = *a;
	int sign = a->sign;
	if (d == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (mag_divmod_u32(&t, d) != 0)
		return QAWS_STATUS_INTERNAL_ERROR;
	if (t.size)
		t.sign = sign;
	copy_int(r, &t);
	return QAWS_STATUS_OK;
}

static unsigned int nlz32(uint32_t x)
{
	unsigned int n = 0;
	if (x == 0)
		return 32;
	while (!(x & 0x80000000u))
	{
		x <<= 1;
		n++;
	}
	return n;
}

/*
 * Magnitudes: quo = |a| / |b|, rem = |a| mod |b| (Knuth, TAOCP 4.3.1,
 * algorithm D, on 32-bit limbs). |b| != 0.
 */
static void mag_divmod(qaws_exact_int* quo, qaws_exact_int* rem, qaws_exact_int const* a, qaws_exact_int const* b)
{
	uint32_t un[QAWS_EXACT_LIMBS + 1], vn[QAWS_EXACT_LIMBS];
	int m = a->size, n = b->size, i, j;
	unsigned int s;
	if (mag_cmp(a, b) < 0)
	{
		qaws_exact_int_zero(quo);
		*rem = *a;
		rem->sign = rem->size ? 1 : 0;
		return;
	}
	if (n == 1)
	{
		/* one limb: short division */
		uint64_t r = 0, d = b->limb[0];
		for (i = m - 1; i >= 0; i--)
		{
			uint64_t cur = (r << 32) | a->limb[i];
			quo->limb[i] = (uint32_t)(cur / d);
			r = cur % d;
		}
		quo->size = m;
		quo->sign = 1;
		normalize(quo);
		rem->limb[0] = (uint32_t)r;
		rem->size = 1;
		rem->sign = 1;
		normalize(rem);
		return;
	}
	/* normalize: the divisor's top bit set */
	s = nlz32(b->limb[n - 1]);
	for (i = n - 1; i > 0; i--)
		vn[i] = (b->limb[i] << s) | (s ? (uint32_t)((uint64_t)b->limb[i - 1] >> (32 - s)) : 0u);
	vn[0] = b->limb[0] << s;
	un[m] = s ? (uint32_t)((uint64_t)a->limb[m - 1] >> (32 - s)) : 0u;
	for (i = m - 1; i > 0; i--)
		un[i] = (a->limb[i] << s) | (s ? (uint32_t)((uint64_t)a->limb[i - 1] >> (32 - s)) : 0u);
	un[0] = a->limb[0] << s;
	for (j = m - n; j >= 0; j--)
	{
		uint64_t num = ((uint64_t)un[j + n] << 32) | un[j + n - 1];
		uint64_t qhat = num / vn[n - 1], rhat = num % vn[n - 1];
		int64_t t, k;
		while (qhat >= ((uint64_t)1 << 32) || qhat * vn[n - 2] > ((rhat << 32) | un[j + n - 2]))
		{
			qhat--;
			rhat += vn[n - 1];
			if (rhat >= ((uint64_t)1 << 32))
				break;
		}
		/* multiply and subtract */
		k = 0;
		for (i = 0; i < n; i++)
		{
			uint64_t p = qhat * vn[i];
			t = (int64_t)un[i + j] - k - (int64_t)(p & 0xFFFFFFFFu);
			un[i + j] = (uint32_t)t;
			k = (int64_t)(p >> 32) - (t >> 32);
		}
		t = (int64_t)un[j + n] - k;
		un[j + n] = (uint32_t)t;
		quo->limb[j] = (uint32_t)qhat;
		if (t < 0)
		{
			/* qhat was one too large: add back */
			uint64_t c = 0;
			quo->limb[j]--;
			for (i = 0; i < n; i++)
			{
				uint64_t sum = (uint64_t)un[i + j] + vn[i] + c;
				un[i + j] = (uint32_t)sum;
				c = sum >> 32;
			}
			un[j + n] = (uint32_t)((uint64_t)un[j + n] + c);
		}
	}
	quo->size = m - n + 1;
	quo->sign = 1;
	normalize(quo);
	/* unnormalize the remainder */
	for (i = 0; i < n; i++)
		rem->limb[i] = (un[i] >> s) | (s ? (uint32_t)((uint64_t)un[i + 1] << (32 - s)) : 0u);
	rem->size = n;
	rem->sign = 1;
	normalize(rem);
}

qaws_status qaws_exact_int_divmod(qaws_exact_int* q, qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b)
{
	qaws_exact_int quo, rem;
	int as = a->sign, bs = b->sign;
	if (b->sign == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (a->sign == 0)
	{
		qaws_exact_int_zero(&quo);
		qaws_exact_int_zero(&rem);
	}
	else
		mag_divmod(&quo, &rem, a, b);
	if (quo.size)
		quo.sign = as * bs;
	if (rem.size)
		rem.sign = as;
	if (q) copy_int(q, &quo);
	if (r) copy_int(r, &rem);
	return QAWS_STATUS_OK;
}

/* Trailing zero bits of a non-zero magnitude. */
unsigned int qaws_exact_int_ctz(qaws_exact_int const* x)
{
	unsigned int i, b;
	for (i = 0; i < (unsigned int)x->size; i++)
		if (x->limb[i])
		{
			uint32_t v = x->limb[i];
			for (b = 0; !(v & 1u); b++)
				v >>= 1;
			return 32 * i + b;
		}
	return 0;
}

void qaws_exact_int_gcd(qaws_exact_int* g, qaws_exact_int const* a, qaws_exact_int const* b)
{
	/* binary gcd (Stein) on magnitudes: shifts and subtractions only */
	qaws_exact_int x, y;
	unsigned int k, kx, ky;
	if (a->sign == 0 || b->sign == 0)
	{
		*g = a->sign == 0 ? *b : *a;
		g->sign = g->size ? 1 : 0;
		return;
	}
	x = *a;
	y = *b;
	x.sign = 1;
	y.sign = 1;
	/* a one-limb operand: Euclid in machine words after one reduction */
	kx = qaws_exact_int_ctz(&x);
	ky = qaws_exact_int_ctz(&y);
	k = kx < ky ? kx : ky;
	qaws_exact_int_shr(&x, &x, kx);
	for (;;)
	{
		qaws_exact_int_shr(&y, &y, qaws_exact_int_ctz(&y));
		if (x.size == 1 && y.size == 1)
		{
			uint32_t u = x.limb[0], v = y.limb[0];
			while (v)
			{
				uint32_t t = u % v;
				u = v;
				v = t;
			}
			qaws_exact_int_from_i64(&x, (int64_t)u);
			break;
		}
		if (mag_cmp(&x, &y) > 0)
		{
			qaws_exact_int t = x;
			x = y;
			y = t;
		}
		mag_sub(&y, &y, &x);
		normalize(&y);
		if (y.size == 0)
			break;
	}
	qaws_exact_int_shl(g, &x, k);
}

/* ------------------------------------------------------------------ */
/*  Conversions                                                       */
/* ------------------------------------------------------------------ */

double qaws_exact_int_to_double(qaws_exact_int const* x)
{
	unsigned int n = qaws_exact_int_bits(x);
	uint64_t m;
	int sticky = 0, i;
	double v;
	if (n == 0)
		return 0.0;
	if (n <= 64)
	{
		m = x->limb[0] | (x->size > 1 ? (uint64_t)x->limb[1] << 32 : 0u);
		v = (double)m;
	}
	else
	{
		/* the top 64 bits, with every lower bit folded into a sticky bit
		   below the rounding position: the uint64 -> double conversion then
		   rounds to nearest even exactly as the full value would */
		qaws_exact_int t;
		qaws_exact_int_shr(&t, x, n - 64);
		m = t.limb[0] | ((uint64_t)t.limb[1] << 32);
		for (i = 0; i < (int)((n - 64) / 32) && !sticky; i++)
			sticky = x->limb[i] != 0;
		if (!sticky && (n - 64) % 32)
			sticky = (x->limb[(n - 64) / 32] & ((1u << ((n - 64) % 32)) - 1u)) != 0;
		if (sticky)
			m |= 1u;
		v = ldexp((double)m, (int)(n - 64));
	}
	return x->sign < 0 ? -v : v;
}

double qaws_exact_ratio_to_double(qaws_exact_int const* num, qaws_exact_int const* den)
{
	qaws_exact_int n = *num, d = *den, step;
	int nb, db, k, i, sign;
	uint64_t q = 0;
	if (den->sign == 0)
		return NAN;
	if (num->sign == 0)
		return 0.0;
	sign = num->sign * den->sign;
	n.sign = 1;
	d.sign = 1;
	/* scale so that the quotient has 63 or 64 bits: q in [2^62, 2^64) */
	nb = (int)qaws_exact_int_bits(&n);
	db = (int)qaws_exact_int_bits(&d);
	k = 63 + db - nb;
	if (k > 0)
	{
		if (qaws_exact_int_shl(&n, &n, (unsigned int)k) != QAWS_STATUS_OK)
			return (double)sign * (qaws_exact_int_to_double(num) / qaws_exact_int_to_double(den));
	}
	else if (k < 0 && qaws_exact_int_shl(&d, &d, (unsigned int)-k) != QAWS_STATUS_OK)
		return (double)sign * (qaws_exact_int_to_double(num) / qaws_exact_int_to_double(den));
	/* binary long division over the 64 possible quotient bits */
	if (qaws_exact_int_shl(&step, &d, 63) != QAWS_STATUS_OK)
		return (double)sign * (qaws_exact_int_to_double(num) / qaws_exact_int_to_double(den));
	for (i = 63; i >= 0; i--)
	{
		if (mag_cmp(&n, &step) >= 0)
		{
			mag_sub(&n, &n, &step);
			normalize(&n);
			q |= (uint64_t)1 << i;
		}
		qaws_exact_int_shr(&step, &step, 1);
	}
	if (n.sign != 0)
		q |= 1u;   /* sticky: far below the rounding position */
	return (double)sign * ldexp((double)q, -k);
}

int qaws_exact_split_double(double d, int64_t* m, int* e)
{
	int ex;
	double f;
	if (!(d == d) || d - d != 0.0)
		return 0;
	if (d == 0.0)
	{
		*m = 0;
		*e = 0;
		return 1;
	}
	f = frexp(d, &ex);
	*m = (int64_t)ldexp(f, 53);
	*e = ex - 53;
	while ((*m & 1) == 0)
	{
		*m /= 2;
		(*e)++;
	}
	return 1;
}

qaws_status qaws_exact_int_from_text(qaws_exact_int* x, char const* text)
{
	int neg = 0;
	qaws_exact_int ten9;
	qaws_status st;
	if (!text)
		return QAWS_STATUS_INVALID_ARGUMENT;
	while (*text == ' ')
		text++;
	if (*text == '-' || *text == '+')
		neg = *text++ == '-';
	if (*text < '0' || *text > '9')
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_exact_int_zero(x);
	qaws_exact_int_from_i64(&ten9, 1000000000);
	while (*text >= '0' && *text <= '9')
	{
		/* nine digits at a time */
		int64_t chunk = 0, scale = 1;
		int k;
		qaws_exact_int c;
		for (k = 0; k < 9 && *text >= '0' && *text <= '9'; k++)
		{
			chunk = chunk * 10 + (*text++ - '0');
			scale *= 10;
		}
		st = qaws_exact_int_mul_i64(x, x, scale);
		if (st != QAWS_STATUS_OK)
			return st;
		qaws_exact_int_from_i64(&c, chunk);
		st = qaws_exact_int_add(x, x, &c);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	(void)ten9;
	if (*text != '\0')
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (neg)
		x->sign = -x->sign;
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_int_to_text(qaws_exact_int const* x, char* buf, size_t cap)
{
	/* base-1e9 digits, least significant first */
	uint32_t chunks[QAWS_EXACT_LIMBS * 2];
	qaws_exact_int t = *x;
	int nc = 0, i;
	size_t len = 0;
	char tmp[16];
	if (!buf || cap == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (x->sign == 0)
	{
		if (cap < 2)
			return QAWS_STATUS_BUFFER_TOO_SMALL;
		buf[0] = '0';
		buf[1] = '\0';
		return QAWS_STATUS_OK;
	}
	while (t.size)
		chunks[nc++] = mag_divmod_u32(&t, 1000000000u);
	if (x->sign < 0)
	{
		if (len + 1 >= cap)
			return QAWS_STATUS_BUFFER_TOO_SMALL;
		buf[len++] = '-';
	}
	for (i = nc - 1; i >= 0; i--)
	{
		int k, w;
		uint32_t c = chunks[i];
		/* the leading chunk without padding, the others on nine digits */
		w = 0;
		do
		{
			tmp[w++] = (char)('0' + c % 10);
			c /= 10;
		} while (c);
		if (i != nc - 1)
			while (w < 9)
				tmp[w++] = '0';
		if (len + (size_t)w >= cap)
			return QAWS_STATUS_BUFFER_TOO_SMALL;
		for (k = w - 1; k >= 0; k--)
			buf[len++] = tmp[k];
	}
	buf[len] = '\0';
	return QAWS_STATUS_OK;
}
