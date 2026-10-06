/*
 * Test 63: Exact multi-limb integers
 *
 *   - add, sub, mul, compare, shifts and the nearest double against
 *     Mathematica exact integers (tests/reference/63_exact_int.wls), values
 *     up to about 1900 bits, carries and borrows across every limb
 *   - randomized identities: (a + b) - b = a, (a * b) / b = a, distributivity,
 *     text round trips
 *   - the capacity edge: results past 2048 bits are reported, never wrapped
 *   - exact decomposition of doubles into m 2^e
 */

#include "test_common.h"
#include "exact/qaws_exact_int.h"
#include "reference/63_exact_int.h"
#include <stdlib.h>
#include <string.h>

#define EI_TEXT 1024

static int ei_equal_text(qaws_exact_int const* x, char const* expect)
{
	char buf[EI_TEXT];
	if (qaws_exact_int_to_text(x, buf, sizeof(buf)) != QAWS_STATUS_OK)
		return 0;
	return strcmp(buf, expect) == 0;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_int) / sizeof(g_ref_exact_int[0])), i;
	int ok_parse = 1, ok_add = 1, ok_sub = 1, ok_mul = 1, ok_cmp = 1, ok_shl = 1, ok_shr = 1, ok_dbl = 1, ok_div = 1, ok_gcd = 1;
	char msg[160];
	for (i = 0; i < n; i++)
	{
		ref_exact_int_row const* r = &g_ref_exact_int[i];
		qaws_exact_int a, b, x;
		ok_parse &= qaws_exact_int_from_text(&a, r->a) == QAWS_STATUS_OK && qaws_exact_int_from_text(&b, r->b) == QAWS_STATUS_OK;
		ok_parse &= ei_equal_text(&a, r->a) && ei_equal_text(&b, r->b);
		ok_add &= qaws_exact_int_add(&x, &a, &b) == QAWS_STATUS_OK && ei_equal_text(&x, r->sum);
		ok_sub &= qaws_exact_int_sub(&x, &a, &b) == QAWS_STATUS_OK && ei_equal_text(&x, r->diff);
		ok_mul &= qaws_exact_int_mul(&x, &a, &b) == QAWS_STATUS_OK && ei_equal_text(&x, r->prod);
		ok_cmp &= qaws_exact_int_cmp(&a, &b) == r->cmp && qaws_exact_int_cmp(&b, &a) == -r->cmp;
		ok_shl &= qaws_exact_int_shl(&x, &a, r->k) == QAWS_STATUS_OK && ei_equal_text(&x, r->shl);
		qaws_exact_int_shr(&x, &a, r->k);
		ok_shr &= ei_equal_text(&x, r->shr);
		/* nearest double, bit for bit (strtod rounds the 40-digit text) */
		ok_dbl &= qaws_exact_int_to_double(&a) == strtod(r->nearest, NULL);
		/* truncated division and gcd */
		if (b.sign != 0)
		{
			qaws_exact_int q, rem;
			ok_div &= qaws_exact_int_divmod(&q, &rem, &a, &b) == QAWS_STATUS_OK && ei_equal_text(&q, r->quo) && ei_equal_text(&rem, r->rem);
		}
		qaws_exact_int_gcd(&x, &a, &b);
		ok_gcd &= ei_equal_text(&x, r->gcd);
	}
	sprintf(msg, "Mathematica reference: %u rows parsed and printed back", n);
	TEST_ASSERT(ok_parse, msg);
	TEST_ASSERT(ok_add, "Mathematica reference: a + b");
	TEST_ASSERT(ok_sub, "Mathematica reference: a - b");
	TEST_ASSERT(ok_mul, "Mathematica reference: a * b");
	TEST_ASSERT(ok_cmp, "Mathematica reference: compare");
	TEST_ASSERT(ok_shl, "Mathematica reference: a * 2^k");
	TEST_ASSERT(ok_shr, "Mathematica reference: sign(a) floor(|a| / 2^k)");
	TEST_ASSERT(ok_dbl, "Mathematica reference: nearest double, bit for bit");
	TEST_ASSERT(ok_div, "Mathematica reference: truncated quotient and remainder");
	TEST_ASSERT(ok_gcd, "Mathematica reference: gcd");
}

/* xorshift64 */
static uint64_t g_ei_state = 0x9E3779B97F4A7C15ull;
static uint64_t ei_next(void)
{
	g_ei_state ^= g_ei_state << 13;
	g_ei_state ^= g_ei_state >> 7;
	g_ei_state ^= g_ei_state << 17;
	return g_ei_state;
}

static void ei_random(qaws_exact_int* x, unsigned int bits)
{
	int i, n = (int)((bits + 31) / 32);
	qaws_exact_int_zero(x);
	for (i = 0; i < n; i++)
		x->limb[i] = (uint32_t)ei_next();
	if (bits % 32)
		x->limb[n - 1] &= (1u << (bits % 32)) - 1u;
	x->size = n;
	x->sign = (ei_next() & 1) ? -1 : 1;
	while (x->size > 0 && x->limb[x->size - 1] == 0)
		x->size--;
	if (x->size == 0)
		x->sign = 0;
}

static void test_identities(void)
{
	int ok_addsub = 1, ok_div = 1, ok_dist = 1, ok_text = 1, i;
	for (i = 0; i < 2000; i++)
	{
		qaws_exact_int a, b, c, s, t, u, v;
		char buf[EI_TEXT];
		unsigned int ba = 1 + (unsigned int)(ei_next() % 1000), bb = 1 + (unsigned int)(ei_next() % 1000);
		uint32_t d = (uint32_t)(ei_next() | 1u);
		ei_random(&a, ba);
		ei_random(&b, bb);
		ei_random(&c, 1 + (unsigned int)(ei_next() % 1000));
		/* (a + b) - b = a */
		qaws_exact_int_add(&s, &a, &b);
		qaws_exact_int_sub(&s, &s, &b);
		ok_addsub &= qaws_exact_int_cmp(&s, &a) == 0;
		/* (a * d) / d = a */
		qaws_exact_int_mul_i64(&t, &a, (int64_t)d);
		ok_div &= qaws_exact_int_divexact_u32(&t, &t, d) == QAWS_STATUS_OK && qaws_exact_int_cmp(&t, &a) == 0;
		/* a (b + c) = a b + a c */
		qaws_exact_int_add(&u, &b, &c);
		qaws_exact_int_mul(&u, &a, &u);
		qaws_exact_int_mul(&v, &a, &b);
		qaws_exact_int_mul(&t, &a, &c);
		qaws_exact_int_add(&v, &v, &t);
		ok_dist &= qaws_exact_int_cmp(&u, &v) == 0;
		/* text round trip */
		qaws_exact_int_to_text(&u, buf, sizeof(buf));
		qaws_exact_int_from_text(&t, buf);
		ok_text &= qaws_exact_int_cmp(&t, &u) == 0;
	}
	TEST_ASSERT(ok_addsub, "2000 random pairs: (a + b) - b = a");
	TEST_ASSERT(ok_div, "2000 random pairs: (a * d) / d = a");
	TEST_ASSERT(ok_dist, "2000 random triples: a (b + c) = a b + a c");
	TEST_ASSERT(ok_text, "2000 random values: decimal text round trip");
}

static void test_capacity(void)
{
	qaws_exact_int a, b, x;
	/* 2^1023 * 2^1024 = 2^2047 fits; 2^1024 * 2^1024 = 2^2048 does not */
	qaws_exact_int_from_i64(&a, 1);
	qaws_exact_int_shl(&a, &a, 1023);
	qaws_exact_int_from_i64(&b, 1);
	qaws_exact_int_shl(&b, &b, 1024);
	TEST_ASSERT(qaws_exact_int_mul(&x, &a, &b) == QAWS_STATUS_OK && qaws_exact_int_bits(&x) == 2048, "2^1023 * 2^1024 fits 2048 bits");
	TEST_ASSERT(qaws_exact_int_mul(&x, &b, &b) == QAWS_STATUS_EXACT_RANGE_EXCEEDED, "2^1024 * 2^1024 is reported");
	TEST_ASSERT(qaws_exact_int_shl(&x, &b, 1024) == QAWS_STATUS_EXACT_RANGE_EXCEEDED, "2^1024 * 2^1024 by shift is reported");
	/* 2^2048 - 1 + 1 overflows */
	qaws_exact_int_from_i64(&a, 1);
	qaws_exact_int_shl(&a, &a, 2047);
	qaws_exact_int_sub(&x, &a, &b);   /* any 2048-bit value */
	qaws_exact_int_add(&x, &a, &a);
	TEST_ASSERT(qaws_exact_int_add(&x, &a, &a) == QAWS_STATUS_EXACT_RANGE_EXCEEDED, "2^2047 + 2^2047 is reported");
}

static void test_split_double(void)
{
	static double const values[] = { 1.0, -0.75, 3.0e-310, 1.0e300, 0.1, -1234567.125 };
	int ok = 1;
	unsigned int i;
	for (i = 0; i < sizeof(values) / sizeof(values[0]); i++)
	{
		int64_t m;
		int e;
		ok &= qaws_exact_split_double(values[i], &m, &e) && (m & 1) && ldexp((double)m, e) == values[i];
	}
	TEST_ASSERT(ok, "doubles split exactly into odd m * 2^e (normal, subnormal, huge)");
}

int test_63_exact_int_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 63: Exact multi-limb integers\n");
	test_reference();
	test_identities();
	test_capacity();
	test_split_double();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
