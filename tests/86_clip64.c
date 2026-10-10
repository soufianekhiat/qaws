/*
 * Test 86: Exact Boolean operations on int64 polygons (qaws_clip64)
 *
 *   - two squares: every clip type, exact areas and counts
 *   - squares sharing an edge, touching at a corner, identical, collinear
 *     overlapping edges at 2^60
 *   - a pentagram under every fill rule (crossings rounded)
 *   - nesting depths, reverse solution, open subjects
 *   - Clipper2's TestIsCollinear #831: a bow tie at 0x4000000000000
 */

#include "test_common.h"
#include <stdint.h>
#include <math.h>
#include <string.h>

static qaws_clip64_result* c6_run(qaws_clip_type ct, qaws_fill_rule fr, int64_t const* const* s, unsigned int const* sn,
	unsigned int ns, int64_t const* const* c, unsigned int const* cn, unsigned int nc, unsigned int flags)
{
	qaws_path64 sp[8], cp[8];
	qaws_clip64_desc d;
	qaws_clip64_result* r = NULL;
	unsigned int i;
	qaws_status st;
	for (i = 0; i < ns; i++) { sp[i].points = s[i]; sp[i].point_count = sn[i]; }
	for (i = 0; i < nc; i++) { cp[i].points = c[i]; cp[i].point_count = cn[i]; }
	memset(&d, 0, sizeof(d));
	d.subjects = sp; d.subject_count = ns;
	d.clips = cp; d.clip_count = nc;
	d.clip_type = ct; d.fill_rule = fr; d.flags = flags;
	st = qaws_clip64_execute(&d, &r);
	if (st != QAWS_STATUS_OK)
		printf("    clip64 status %d\n", (int)st);
	return r;
}

static double c6_area(qaws_clip64_result const* r)
{
	unsigned int i, n = qaws_clip64_result_get_path_count(r);
	double sum = 0.0;
	for (i = 0; i < n; i++)
	{
		int64_t const* p;
		unsigned int m;
		qaws_clip64_result_get_path(r, i, &p, &m);
		sum += 0.5 * qaws_path64_area2(p, m);
	}
	return sum;
}

static unsigned int c6_points(qaws_clip64_result const* r, unsigned int i)
{
	int64_t const* p;
	unsigned int m = 0;
	qaws_clip64_result_get_path(r, i, &p, &m);
	return m;
}

static void test_squares(void)
{
	static int64_t const a[8] = { 0, 0, 20, 0, 20, 20, 0, 20 }, b[8] = { 10, 10, 30, 10, 30, 30, 10, 30 };
	static double const areas[5] = { 0, 100, 700, 300, 600 };
	static unsigned int const counts[5] = { 0, 1, 1, 1, 2 };
	int64_t const* s = a;
	int64_t const* c = b;
	unsigned int n4 = 4, ct;
	char msg[200];
	for (ct = 1; ct <= 4; ct++)
	{
		qaws_clip64_result* r = c6_run((qaws_clip_type)ct, QAWS_FILL_NON_ZERO, &s, &n4, 1, &c, &n4, 1, 0);
		sprintf(msg, "int64 squares, clip type %u: area %g (%g), %u paths (%u)", ct, areas[ct], r ? c6_area(r) : -1.0,
			counts[ct], r ? qaws_clip64_result_get_path_count(r) : 99);
		TEST_ASSERT(r && c6_area(r) == areas[ct] && qaws_clip64_result_get_path_count(r) == counts[ct], msg);
		if (r && ct == QAWS_CLIP_UNION)
			TEST_ASSERT(c6_points(r, 0) == 8 && !qaws_clip64_result_is_hole(r, 0), "union: one outer path of 8 corners");
		qaws_clip64_result_destroy(r);
	}
	{
		qaws_clip64_result* r = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n4, 1, &c, &n4, 1, QAWS_CLIP_REVERSE_SOLUTION);
		TEST_ASSERT(r && c6_area(r) == -700, "reverse solution: area -700");
		qaws_clip64_result_destroy(r);
	}
}

static void test_degenerate(void)
{
	char msg[256];
	unsigned int n4 = 4;
	/* sharing an edge, collinear corners dropped or kept */
	{
		static int64_t const a[8] = { 0, 0, 10, 0, 10, 10, 0, 10 }, b[8] = { 10, 0, 20, 0, 20, 10, 10, 10 };
		int64_t const* s = a;
		int64_t const* c = b;
		qaws_clip64_result* r = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n4, 1, &c, &n4, 1, 0);
		qaws_clip64_result* k = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n4, 1, &c, &n4, 1, QAWS_CLIP_PRESERVE_COLLINEAR);
		qaws_clip64_result* i = c6_run(QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO, &s, &n4, 1, &c, &n4, 1, 0);
		sprintf(msg, "int64 shared edge: one rectangle of area 200 with 4 corners (%u, %g, %u), 6 kept (%u); intersection empty (%u)",
			r ? qaws_clip64_result_get_path_count(r) : 99, r ? c6_area(r) : 0, r ? c6_points(r, 0) : 0, k ? c6_points(k, 0) : 0,
			i ? qaws_clip64_result_get_path_count(i) : 99);
		TEST_ASSERT(r && qaws_clip64_result_get_path_count(r) == 1 && c6_area(r) == 200 && c6_points(r, 0) == 4 &&
			k && c6_points(k, 0) == 6 && i && qaws_clip64_result_get_path_count(i) == 0, msg);
		qaws_clip64_result_destroy(r); qaws_clip64_result_destroy(k); qaws_clip64_result_destroy(i);
	}
	/* touching at a corner; identical */
	{
		static int64_t const a[8] = { 0, 0, 10, 0, 10, 10, 0, 10 }, b[8] = { 10, 10, 20, 10, 20, 20, 10, 20 };
		int64_t const* s = a;
		int64_t const* c = b;
		int64_t const* same = a;
		qaws_clip64_result* r = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n4, 1, &c, &n4, 1, 0);
		qaws_clip64_result* u = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n4, 1, &same, &n4, 1, 0);
		qaws_clip64_result* x = c6_run(QAWS_CLIP_XOR, QAWS_FILL_NON_ZERO, &s, &n4, 1, &same, &n4, 1, 0);
		sprintf(msg, "int64 corner touch: 2 paths (%u); identical: union 1 of area 100 (%u, %g), xor empty (%u)",
			r ? qaws_clip64_result_get_path_count(r) : 99, u ? qaws_clip64_result_get_path_count(u) : 99, u ? c6_area(u) : 0,
			x ? qaws_clip64_result_get_path_count(x) : 99);
		TEST_ASSERT(r && qaws_clip64_result_get_path_count(r) == 2 && u && qaws_clip64_result_get_path_count(u) == 1 &&
			c6_area(u) == 100 && x && qaws_clip64_result_get_path_count(x) == 0, msg);
		qaws_clip64_result_destroy(r); qaws_clip64_result_destroy(u); qaws_clip64_result_destroy(x);
	}
	/* overlapping collinear edges at 2^60: [0, 2 L] x [0, L] and [L, 3 L] x [0, L] share [L, 2 L] on both edges */
	{
		int64_t L = (int64_t)1 << 60;
		int64_t a[8], b[8];
		int64_t const* s = a;
		int64_t const* c = b;
		qaws_clip64_result* r;
		a[0] = 0; a[1] = 0; a[2] = 2 * L / 2 * 1; a[3] = 0; a[4] = 2 * L / 2; a[5] = L / 2; a[6] = 0; a[7] = L / 2;
		b[0] = L / 2; b[1] = 0; b[2] = 3 * (L / 2); b[3] = 0; b[4] = 3 * (L / 2); b[5] = L / 2; b[6] = L / 2; b[7] = L / 2;
		r = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n4, 1, &c, &n4, 1, 0);
		sprintf(msg, "rectangles at 2^60 overlapping on half their edges: one rectangle of 4 corners (%u paths, %u corners), area exact (%.17g)",
			r ? qaws_clip64_result_get_path_count(r) : 99, r ? c6_points(r, 0) : 0, r ? c6_area(r) : 0);
		TEST_ASSERT(r && qaws_clip64_result_get_path_count(r) == 1 && c6_points(r, 0) == 4 &&
			c6_area(r) == (double)(3 * (L / 2)) * (double)(L / 2), msg);
		qaws_clip64_result_destroy(r);
	}
	/* Clipper2 #831: a bow tie at 0x4000000000000 is two triangles */
	{
		int64_t const v = 0x4000000000000LL;
		int64_t bt[8];
		int64_t const* s = bt;
		qaws_clip64_result* r;
		bt[0] = 0; bt[1] = 0; bt[2] = v; bt[3] = v; bt[4] = v; bt[5] = 0; bt[6] = 0; bt[7] = v;
		r = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n4, 1, NULL, NULL, 0, 0);
		sprintf(msg, "Clipper2 #831: a bow tie at 0x4000000000000 unites into 2 triangles (%u)", r ? qaws_clip64_result_get_path_count(r) : 99);
		TEST_ASSERT(r && qaws_clip64_result_get_path_count(r) == 2, msg);
		qaws_clip64_result_destroy(r);
	}
}

static void test_fill_and_nesting(void)
{
	char msg[256];
	/* pentagram, radius 1000 */
	{
		int64_t st[10];
		int64_t const* s = st;
		unsigned int i, n5 = 5;
		qaws_clip64_result *nz, *eo, *po, *ne;
		for (i = 0; i < 5; i++)
		{
			double a = 3.14159265358979323846 / 2 + 2 * 3.14159265358979323846 * (2 * i % 5) / 5;
			st[2 * i] = (int64_t)floor(1000 * cos(a) + 0.5);
			st[2 * i + 1] = (int64_t)floor(1000 * sin(a) + 0.5);
		}
		nz = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &s, &n5, 1, NULL, NULL, 0, 0);
		eo = c6_run(QAWS_CLIP_UNION, QAWS_FILL_EVEN_ODD, &s, &n5, 1, NULL, NULL, 0, 0);
		po = c6_run(QAWS_CLIP_UNION, QAWS_FILL_POSITIVE, &s, &n5, 1, NULL, NULL, 0, 0);
		ne = c6_run(QAWS_CLIP_UNION, QAWS_FILL_NEGATIVE, &s, &n5, 1, NULL, NULL, 0, 0);
		sprintf(msg, "int64 pentagram: non-zero 1 path of 10 corners (%u, %u), even-odd 5 (%u), positive 1 (%u), negative 0 (%u)",
			nz ? qaws_clip64_result_get_path_count(nz) : 99, nz ? c6_points(nz, 0) : 0, eo ? qaws_clip64_result_get_path_count(eo) : 99,
			po ? qaws_clip64_result_get_path_count(po) : 99, ne ? qaws_clip64_result_get_path_count(ne) : 99);
		TEST_ASSERT(nz && qaws_clip64_result_get_path_count(nz) == 1 && c6_points(nz, 0) == 10 && eo &&
			qaws_clip64_result_get_path_count(eo) == 5 && po && qaws_clip64_result_get_path_count(po) == 1 &&
			ne && qaws_clip64_result_get_path_count(ne) == 0, msg);
		qaws_clip64_result_destroy(nz); qaws_clip64_result_destroy(eo); qaws_clip64_result_destroy(po); qaws_clip64_result_destroy(ne);
	}
	/* square, hole, island (even-odd) */
	{
		static int64_t const o[8] = { 0, 0, 100, 0, 100, 100, 0, 100 }, h[8] = { 20, 20, 80, 20, 80, 80, 20, 80 };
		static int64_t const i[8] = { 40, 40, 60, 40, 60, 60, 40, 60 };
		int64_t const* s[3];
		unsigned int sn[3] = { 4, 4, 4 }, k, n, d[3] = { 0, 0, 0 }, ok = 1;
		qaws_clip64_result* r;
		s[0] = o; s[1] = h; s[2] = i;
		r = c6_run(QAWS_CLIP_UNION, QAWS_FILL_EVEN_ODD, s, sn, 3, NULL, NULL, 0, 0);
		n = r ? qaws_clip64_result_get_path_count(r) : 0;
		for (k = 0; k < n; k++)
		{
			unsigned int dep = qaws_clip64_result_get_depth(r, k), p = qaws_clip64_result_get_parent(r, k);
			if (dep < 3) d[dep]++;
			if (dep > 0 && (p == QAWS_CLIP_NONE_INDEX || qaws_clip64_result_get_depth(r, p) != dep - 1)) ok = 0;
		}
		sprintf(msg, "int64 nesting: 3 paths at depths 0, 1, 2 (%u: %u %u %u), area 10000 - 3600 + 400 (%g)", n, d[0], d[1], d[2], r ? c6_area(r) : 0);
		TEST_ASSERT(n == 3 && d[0] == 1 && d[1] == 1 && d[2] == 1 && ok && c6_area(r) == 6800, msg);
		qaws_clip64_result_destroy(r);
	}
	/* open subject: a line through a square */
	{
		static int64_t const sq[8] = { 0, 0, 20, 0, 20, 20, 0, 20 }, ln[4] = { -10, 10, 30, 10 };
		qaws_path64 clip, open;
		qaws_clip64_desc d;
		qaws_clip64_result* r = NULL;
		unsigned int ni = 0, no = 0;
		clip.points = sq; clip.point_count = 4;
		open.points = ln; open.point_count = 2;
		memset(&d, 0, sizeof(d));
		d.clips = &clip; d.clip_count = 1;
		d.open_subjects = &open; d.open_subject_count = 1;
		d.clip_type = QAWS_CLIP_INTERSECTION; d.fill_rule = QAWS_FILL_NON_ZERO;
		qaws_clip64_execute(&d, &r);
		if (r)
		{
			int64_t const* p;
			unsigned int m = 0;
			ni = qaws_clip64_result_get_open_path_count(r);
			if (ni) qaws_clip64_result_get_open_path(r, 0, &p, &m);
			ni = ni == 1 && m == 2 && p[0] == 0 && p[2] == 20 ? 1 : 0;
		}
		qaws_clip64_result_destroy(r);
		r = NULL;
		d.clip_type = QAWS_CLIP_DIFFERENCE;
		qaws_clip64_execute(&d, &r);
		no = r ? qaws_clip64_result_get_open_path_count(r) : 0;
		qaws_clip64_result_destroy(r);
		sprintf(msg, "int64 open line: intersection keeps 0..20 (%u), difference the 2 outside pieces (%u)", ni, no);
		TEST_ASSERT(ni == 1 && no == 2, msg);
	}
}

int test_86_clip64_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 86: Exact Boolean operations on int64 polygons\n");
	test_squares();
	test_degenerate();
	test_fill_and_nesting();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
