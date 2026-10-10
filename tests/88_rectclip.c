/*
 * Test 88: Rectangle clipping and Minkowski sums
 *
 *   - a polygon and a circle clipped to a rectangle (areas in closed form),
 *     paths clipped on their own (two overlapping squares stay two)
 *   - an open zig-zag clipped to a rectangle: the pieces inside
 *   - int64: the same rectangles exactly
 *   - Minkowski: a square swept along a square path (band area 12^2 - 8^2),
 *     an open segment (a stadium of squares), sum and difference, int64
 */

#include "test_common.h"
#include <stdint.h>
#include <math.h>
#include <string.h>

static double rct_area(qaws_clip_result const* r)
{
	unsigned int i, n = qaws_clip_result_get_path_count(r);
	double s = 0.0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d p;
		qaws_scalar a = 0;
		qaws_clip_result_get_path(r, i, &p);
		qaws_path_compute_area_2d(&p, &a);
		s += a;
	}
	return s;
}

static double rct_area64(qaws_clip64_result const* r)
{
	unsigned int i, n = qaws_clip64_result_get_path_count(r);
	double s = 0.0;
	for (i = 0; i < n; i++)
	{
		int64_t const* p;
		unsigned int m;
		qaws_clip64_result_get_path(r, i, &p, &m);
		s += 0.5 * qaws_path64_area2(p, m);
	}
	return s;
}

static void test_rect(void)
{
	qaws_scalar const rect[4] = { 0, 0, 10, 5 };
	qaws_scalar sa[8] = { -5, -5, 5, -5, 5, 5, -5, 5 }, sb[8] = { 2, 2, 8, 2, 8, 8, 2, 8 };
	qaws_curve *a = NULL, *b = NULL, *circle = NULL;
	qaws_curve const* cs[3];
	qaws_path_2d p[3];
	qaws_clip_result* r = NULL;
	char msg[256];
	unsigned int i;
	qaws_curve_create_polyline_2d(sa, 4, 1, &a);
	qaws_curve_create_polyline_2d(sb, 4, 1, &b);
	qaws_curve_create_ellipse_2d(10, 5, 2, 2, 0, &circle);
	cs[0] = a; cs[1] = b; cs[2] = circle;
	for (i = 0; i < 3; i++) { p[i].curves = &cs[i]; p[i].curve_count = 1; p[i].closed = 1; }
	qaws_rect_clip_2d(rect, p, 3, &r);
	/* a: [0,5]x[0,5] = 25; b: [2,8]x[2,5] = 18; circle at the corner: a quarter, pi */
	sprintf(msg, "rect clip: 3 paths on their own (%u), areas 25 + 18 + pi (%.12g)", r ? qaws_clip_result_get_path_count(r) : 0, r ? rct_area(r) : 0);
	TEST_ASSERT(r && qaws_clip_result_get_path_count(r) == 3 && fabs(rct_area(r) - (43 + 3.14159265358979323846)) < 1e-6, msg);
	qaws_clip_result_destroy(r);
	{
		qaws_scalar z[12] = { -2, 1, 3, 4, 6, 1, 9, 4, 12, 1, 14, 3 };
		qaws_curve* zz = NULL;
		qaws_curve const* zc;
		qaws_path_2d zp;
		unsigned int n;
		qaws_curve_create_polyline_2d(z, 6, 0, &zz);
		zc = zz;
		zp.curves = &zc; zp.curve_count = 1; zp.closed = 0;
		r = NULL;
		qaws_rect_clip_lines_2d(rect, &zp, 1, &r);
		n = r ? qaws_clip_result_get_open_path_count(r) : 0;
		sprintf(msg, "rect clip lines: a zig-zag leaves one piece inside (%u)", n);
		TEST_ASSERT(n == 1, msg);
		qaws_clip_result_destroy(r);
		qaws_curve_destroy(zz);
	}
	{
		int64_t const r64[4] = { 0, 0, 10, 5 };
		int64_t pa[8] = { -5, -5, 5, -5, 5, 5, -5, 5 }, pb[8] = { 2, 2, 8, 2, 8, 8, 2, 8 };
		qaws_path64 q[2];
		qaws_clip64_result* x = NULL;
		q[0].points = pa; q[0].point_count = 4;
		q[1].points = pb; q[1].point_count = 4;
		qaws_rect_clip64(r64, q, 2, &x);
		sprintf(msg, "int64 rect clip: 2 paths, area 43 (%u, %g)", x ? qaws_clip64_result_get_path_count(x) : 0, x ? rct_area64(x) : 0);
		TEST_ASSERT(x && qaws_clip64_result_get_path_count(x) == 2 && rct_area64(x) == 43, msg);
		qaws_clip64_result_destroy(x);
	}
	qaws_curve_destroy(a); qaws_curve_destroy(b); qaws_curve_destroy(circle);
}

static void test_minkowski(void)
{
	qaws_scalar pat[8] = { -1, -1, 1, -1, 1, 1, -1, 1 }, sq[8] = { 0, 0, 10, 0, 10, 10, 0, 10 }, seg[4] = { 0, 0, 10, 0 };
	qaws_clip_result* r = NULL;
	char msg[256];
	qaws_minkowski_sum_2d(pat, 4, sq, 4, 1, &r);
	sprintf(msg, "Minkowski sum of a square along a square path: the band 144 - 64 (%.12g), 2 paths (%u)",
		r ? rct_area(r) : 0, r ? qaws_clip_result_get_path_count(r) : 0);
	TEST_ASSERT(r && fabs(rct_area(r) - 80) < 1e-9 && qaws_clip_result_get_path_count(r) == 2, msg);
	qaws_clip_result_destroy(r);
	r = NULL;
	qaws_minkowski_sum_2d(pat, 4, seg, 2, 0, &r);
	sprintf(msg, "Minkowski sum of a square along a segment: 12 x 2 (%.12g)", r ? rct_area(r) : 0);
	TEST_ASSERT(r && fabs(rct_area(r) - 24) < 1e-9, msg);
	qaws_clip_result_destroy(r);
	{
		/* a triangle pattern: sum and difference are point reflections of each other */
		qaws_scalar tri[6] = { 0, 0, 2, 0, 0, 1 };
		qaws_clip_result *s1 = NULL, *s2 = NULL;
		qaws_minkowski_sum_2d(tri, 3, seg, 2, 0, &s1);
		qaws_minkowski_diff_2d(tri, 3, seg, 2, 0, &s2);
		sprintf(msg, "Minkowski sum and difference of a triangle along a segment: equal areas (%.12g, %.12g)",
			s1 ? rct_area(s1) : 0, s2 ? rct_area(s2) : 0);
		TEST_ASSERT(s1 && s2 && fabs(rct_area(s1) - rct_area(s2)) < 1e-9 && rct_area(s1) > 0, msg);
		qaws_clip_result_destroy(s1);
		qaws_clip_result_destroy(s2);
	}
	{
		int64_t p64[8] = { -1, -1, 1, -1, 1, 1, -1, 1 }, s64[8] = { 0, 0, 10, 0, 10, 10, 0, 10 };
		qaws_clip64_result* x = NULL;
		qaws_minkowski_sum64(p64, 4, s64, 4, 1, &x);
		sprintf(msg, "int64 Minkowski band: area 80 exactly (%g)", x ? rct_area64(x) : 0);
		TEST_ASSERT(x && rct_area64(x) == 80, msg);
		qaws_clip64_result_destroy(x);
	}
}

int test_88_rectclip_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 88: Rectangle clipping and Minkowski sums\n");
	test_rect();
	test_minkowski();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
