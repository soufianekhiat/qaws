/*
 * Test 64: Certified predicates on doubles
 *
 *   - orient2d, orient3d and ratio comparison against Mathematica exact
 *     signs (tests/reference/64_exact_predicates.wls): the Kettner et al.
 *     near-degenerate grid, points rounded onto random lines and planes,
 *     exact zeros, nearly equal ratios
 *   - the naive f64 determinant is shown to be wrong on the same sets
 *   - 200000 random near-collinear triples: the predicate (filter or exact
 *     path) against an independent exact integer determinant
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "exact/qaws_exact_int.h"
#include "reference/64_exact_predicates.h"
#include <math.h>

static double ep_value(ref_me v)
{
	return ldexp((double)v.m, v.e);
}

static int ep_naive2(double const* a, double const* b, double const* c)
{
	double det = (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]);
	return det > 0 ? 1 : (det < 0 ? -1 : 0);
}

static void test_orient2d_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_o2) / sizeof(g_ref_o2[0])), i, naive_wrong = 0, exact_path = 0;
	int ok = 1;
	char msg[200];
	for (i = 0; i < n; i++)
	{
		double a[2], b[2], c[2];
		qaws_exact_sign s;
		qaws_exact_path path;
		a[0] = ep_value(g_ref_o2[i].p[0]); a[1] = ep_value(g_ref_o2[i].p[1]);
		b[0] = ep_value(g_ref_o2[i].p[2]); b[1] = ep_value(g_ref_o2[i].p[3]);
		c[0] = ep_value(g_ref_o2[i].p[4]); c[1] = ep_value(g_ref_o2[i].p[5]);
		ok &= qaws_exact_orient2d(a, b, c, &s, &path) == QAWS_STATUS_OK && (int)s == g_ref_o2[i].sign;
		naive_wrong += ep_naive2(a, b, c) != g_ref_o2[i].sign;
		exact_path += path == QAWS_EXACT_PATH_EXACT;
	}
	printf("    orient2d: %u cases, naive f64 wrong on %u, exact fallback on %u\n", n, naive_wrong, exact_path);
	sprintf(msg, "orient2d matches the Mathematica exact sign on all %u cases", n);
	TEST_ASSERT(ok, msg);
	TEST_ASSERT(naive_wrong > 0, "the naive f64 orientation is wrong on part of the reference set");
}

static void test_orient3d_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_o3) / sizeof(g_ref_o3[0])), i, k, exact_path = 0, naive_wrong = 0;
	int ok = 1;
	char msg[200];
	for (i = 0; i < n; i++)
	{
		double p[12], det;
		qaws_exact_sign s;
		qaws_exact_path path;
		for (k = 0; k < 12; k++)
			p[k] = ep_value(g_ref_o3[i].p[k]);
		ok &= qaws_exact_orient3d(p, p + 3, p + 6, p + 9, &s, &path) == QAWS_STATUS_OK && (int)s == g_ref_o3[i].sign;
		exact_path += path == QAWS_EXACT_PATH_EXACT;
		{
			double adx = p[0] - p[9], ady = p[1] - p[10], adz = p[2] - p[11];
			double bdx = p[3] - p[9], bdy = p[4] - p[10], bdz = p[5] - p[11];
			double cdx = p[6] - p[9], cdy = p[7] - p[10], cdz = p[8] - p[11];
			det = adz * (bdx * cdy - cdx * bdy) + bdz * (cdx * ady - adx * cdy) + cdz * (adx * bdy - bdx * ady);
			naive_wrong += (det > 0 ? 1 : (det < 0 ? -1 : 0)) != g_ref_o3[i].sign;
		}
	}
	printf("    orient3d: %u near-coplanar cases, naive f64 wrong on %u, exact fallback on %u\n", n, naive_wrong, exact_path);
	sprintf(msg, "orient3d matches the Mathematica exact sign on all %u cases", n);
	TEST_ASSERT(ok, msg);
}

static void test_ratio_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_ratio) / sizeof(g_ref_ratio[0])), i, naive_wrong = 0;
	int ok = 1;
	char msg[200];
	for (i = 0; i < n; i++)
	{
		double a = ep_value(g_ref_ratio[i].p[0]), b = ep_value(g_ref_ratio[i].p[1]), c = ep_value(g_ref_ratio[i].p[2]),
			d = ep_value(g_ref_ratio[i].p[3]), q;
		qaws_exact_sign s;
		ok &= qaws_exact_compare_ratio(a, b, c, d, &s, NULL) == QAWS_STATUS_OK && (int)s == g_ref_ratio[i].sign;
		q = a / b - c / d;
		naive_wrong += (q > 0 ? 1 : (q < 0 ? -1 : 0)) != g_ref_ratio[i].sign;
	}
	printf("    ratios: %u nearly equal pairs, naive f64 wrong on %u\n", n, naive_wrong);
	sprintf(msg, "ratio comparison matches the Mathematica exact sign on all %u cases", n);
	TEST_ASSERT(ok, msg);
}

/* xorshift64 */
static uint64_t g_ep_state = 0xD1B54A32D192ED03ull;
static double ep_uniform(void)
{
	g_ep_state ^= g_ep_state << 13;
	g_ep_state ^= g_ep_state >> 7;
	g_ep_state ^= g_ep_state << 17;
	return (double)(g_ep_state >> 11) / 9007199254740992.0;
}

/* Independent exact sign: every coordinate scaled by 2^400 (all values
   here are above 2^-400 or zero) and the determinant in exact integers. */
static int ep_exact2(double const* a, double const* b, double const* c)
{
	double v[6];
	qaws_exact_int x[6], p, q, r, s, l, rr;
	int i;
	v[0] = a[0]; v[1] = a[1]; v[2] = b[0]; v[3] = b[1]; v[4] = c[0]; v[5] = c[1];
	for (i = 0; i < 6; i++)
	{
		int64_t m;
		int e;
		qaws_exact_split_double(v[i], &m, &e);
		qaws_exact_int_from_i64(&x[i], m);
		qaws_exact_int_shl(&x[i], &x[i], (unsigned int)(e + 400));
	}
	qaws_exact_int_sub(&p, &x[0], &x[4]);
	qaws_exact_int_sub(&q, &x[3], &x[5]);
	qaws_exact_int_sub(&r, &x[1], &x[5]);
	qaws_exact_int_sub(&s, &x[2], &x[4]);
	qaws_exact_int_mul(&l, &p, &q);
	qaws_exact_int_mul(&rr, &r, &s);
	return qaws_exact_int_cmp(&l, &rr);
}

static void test_fuzz(void)
{
	unsigned int i, n = 200000, exact_path = 0, naive_wrong = 0;
	int ok = 1;
	char msg[200];
	for (i = 0; i < n; i++)
	{
		double a[2], b[2], c[2], t;
		qaws_exact_sign s;
		qaws_exact_path path;
		int truth;
		a[0] = ep_uniform() * 200 - 100;
		a[1] = ep_uniform() * 200 - 100;
		b[0] = ep_uniform() * 200 - 100;
		b[1] = ep_uniform() * 200 - 100;
		t = ep_uniform() * 4 - 1.5;
		/* c on the line ab up to rounding, sometimes nudged by an ulp */
		c[0] = a[0] + t * (b[0] - a[0]);
		c[1] = a[1] + t * (b[1] - a[1]);
		if (i % 3 == 1) c[0] = nextafter(c[0], 1e300);
		if (i % 3 == 2) c[1] = nextafter(c[1], -1e300);
		truth = ep_exact2(a, b, c);
		ok &= qaws_exact_orient2d(a, b, c, &s, &path) == QAWS_STATUS_OK && (int)s == truth;
		exact_path += path == QAWS_EXACT_PATH_EXACT;
		naive_wrong += ep_naive2(a, b, c) != truth;
	}
	printf("    fuzz: %u near-collinear triples, naive f64 wrong on %u, exact fallback on %u\n", n, naive_wrong, exact_path);
	sprintf(msg, "%u random near-collinear triples: no wrong sign", n);
	TEST_ASSERT(ok, msg);
}

static void test_arguments(void)
{
	double a[2] = { 0, 0 }, b[2] = { 1, 0 }, c[2] = { HUGE_VAL, 0 };
	qaws_exact_sign s;
	TEST_ASSERT(qaws_exact_orient2d(a, b, c, &s, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "non-finite input is refused");
	TEST_ASSERT(qaws_exact_compare_ratio(1, 0, 1, 2, &s, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "zero denominator is refused");
	c[0] = 1e-300;
	c[1] = 1e300;
	TEST_ASSERT(qaws_exact_orient2d(a, b, c, &s, NULL) == QAWS_STATUS_OK && s == QAWS_EXACT_POSITIVE,
		"orient2d across 600 decades of magnitude");
}

int test_64_exact_predicates_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 64: Certified predicates on doubles\n");
	test_orient2d_reference();
	test_orient3d_reference();
	test_ratio_reference();
	test_fuzz();
	test_arguments();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
