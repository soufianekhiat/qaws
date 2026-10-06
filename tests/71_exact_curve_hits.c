/*
 * Test 71: Certified curve / curve intersections
 *
 *   - against Mathematica (tests/reference/71_exact_curve_hits.wls): every
 *     intersection of two Bezier / rational Bezier curves (degrees 1..4)
 *     is found, both parameters enclosed, dyadic ones exact
 *   - B-spline against B-spline: a crossing on a knot is reported once
 *   - overlapping curves are refused (a common component)
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "reference/71_exact_curve_hits.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static qaws_curve* ec_curve(int degree, int rational, long long const* cp, long long const* w)
{
	qaws_scalar cps[5 * 2], ws[5];
	qaws_curve* c = NULL;
	int i, n = degree + 1;
	for (i = 0; i < n * 2; i++)
		cps[i] = (qaws_scalar)ldexp((double)cp[i], -10);
	for (i = 0; i < n; i++)
		ws[i] = (qaws_scalar)w[i];
	if (rational)
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = (unsigned int)degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)n;
		d.weights = ws;
		d.weight_count = (unsigned int)n;
		qaws_curve_create_rational_bezier(&d, &c);
	}
	else
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = (unsigned int)degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)n;
		qaws_curve_create_bezier(&d, &c);
	}
	return c;
}

static int ec_inside(double lo, double hi, char const* text)
{
	double v = strtod(text, NULL);
	return lo <= v && v <= hi;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_curve_hits) / sizeof(g_ref_exact_curve_hits[0])), i, total = 0, exact = 0, bad = 0;
	double widest = 0;
	char msg[200];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	for (i = 0; i < n; i++)
	{
		ref_exact_curve_hits const* r = &g_ref_exact_curve_hits[i];
		qaws_curve* ca = ec_curve(r->da, r->ra, r->acp, r->aw);
		qaws_curve* cb = ec_curve(r->db, r->rb, r->bcp, r->bw);
		qaws_exact_curve* ea = NULL;
		qaws_exact_curve* eb = NULL;
		qaws_exact_pair hits[32];
		unsigned int count = 0, k;
		qaws_status st;
		int ok;
		qaws_exact_curve_prepare(&desc, ca, &ea, NULL);
		qaws_exact_curve_prepare(&desc, cb, &eb, NULL);
		st = qaws_exact_curve_curve_hits(ea, eb, hits, 32, &count);
		ok = st == QAWS_STATUS_OK && count == (unsigned int)r->n;
		for (k = 0; ok && k < count; k++)
		{
			ok &= ec_inside(hits[k].a_lo, hits[k].a_hi, r->s[k]) && ec_inside(hits[k].b_lo, hits[k].b_hi, r->r[k]);
			if (r->sdy[k] && r->rdy[k])
				ok &= hits[k].kind == QAWS_EXACT_HIT_POINT && hits[k].a_lo == hits[k].a_hi && hits[k].b_lo == hits[k].b_hi;
			exact += hits[k].kind == QAWS_EXACT_HIT_POINT;
			if (hits[k].a_hi - hits[k].a_lo > widest) widest = hits[k].a_hi - hits[k].a_lo;
			if (hits[k].b_hi - hits[k].b_lo > widest) widest = hits[k].b_hi - hits[k].b_lo;
		}
		total += (unsigned int)r->n;
		if (!ok)
		{
			bad++;
			printf("      case %u (deg %d%s x %d%s): status %d, %u hits, expected %d\n", i, r->da, r->ra ? "r" : "", r->db, r->rb ? "r" : "", (int)st,
				count, r->n);
			for (k = 0; k < count && k < 4; k++)
				printf("        a [%.17g, %.17g] b [%.17g, %.17g]  vs  %s, %s\n", hits[k].a_lo, hits[k].a_hi, hits[k].b_lo, hits[k].b_hi,
					k < (unsigned int)r->n ? r->s[k] : "-", k < (unsigned int)r->n ? r->r[k] : "-");
		}
		qaws_exact_curve_destroy(ea);
		qaws_exact_curve_destroy(eb);
		qaws_curve_destroy(ca);
		qaws_curve_destroy(cb);
	}
	printf("    %u curve pairs: %u intersections (%u exact points), widest enclosure %.1e\n", n, total, exact, widest);
	sprintf(msg, "every intersection found and enclosed on all %u Mathematica pairs (degrees 1..4, rational)", n);
	TEST_ASSERT(bad == 0, msg);
	TEST_ASSERT(widest < 1e-12, "both parameters enclosed tightly");
}

static void test_splines(void)
{
	/* two degree-2 B-splines meeting at a knot point of each, (2, 2), and crossing
	   once more at t = 0.6 on both (on [0, 1] both have x = 2t, y_b - y_a = (5t - 3)(t - 1)) */
	qaws_scalar ca_pts[8] = { 0, 0, 1, 2, 3, 2, 4, 0 }, cb_pts[8] = { 0, 3, 1, 1, 3, 3, 4, 1 }, knots[7] = { 0, 0, 0, 1, 2, 2, 2 };
	qaws_scalar la[4] = { 0, 0, 4, 4 }, lb[4] = { 0, 4, 2, 2 };
	qaws_bspline_desc d;
	qaws_bezier_desc bd;
	qaws_curve* ca = NULL;
	qaws_curve* cb = NULL;
	qaws_curve* sa = NULL;
	qaws_curve* sb = NULL;
	qaws_exact_curve* ea = NULL;
	qaws_exact_curve* eb = NULL;
	qaws_exact_curve* esa = NULL;
	qaws_exact_curve* esb = NULL;
	qaws_exact_pair hits[8];
	unsigned int count = 0;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = ca_pts;
	d.control_point_count = 4;
	d.knots = knots;
	d.knot_count = 7;
	qaws_curve_create_bspline(&d, &ca);
	d.control_points = cb_pts;
	qaws_curve_create_bspline(&d, &cb);
	qaws_exact_curve_prepare(NULL, ca, &ea, NULL);
	qaws_exact_curve_prepare(NULL, cb, &eb, NULL);
	TEST_ASSERT(qaws_exact_curve_curve_hits(ea, eb, hits, 8, &count) == QAWS_STATUS_OK && count == 2 && hits[0].kind == QAWS_EXACT_HIT_CROSSING &&
		hits[0].a_lo <= 0.6 && 0.6 <= hits[0].a_hi && hits[0].b_lo <= 0.6 && 0.6 <= hits[0].b_hi && hits[1].kind == QAWS_EXACT_HIT_POINT &&
		hits[1].a_lo == 1 && hits[1].b_lo == 1, "B-spline / B-spline: the crossing at 0.6, and the knot point (1, 1) reported once");
	/* segments: exact rational intersection, and the shared-end rule */
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D;
	bd.degree = 1;
	bd.control_points = la;
	bd.control_point_count = 2;
	qaws_curve_create_bezier(&bd, &sa);
	bd.control_points = lb;
	qaws_curve_create_bezier(&bd, &sb);
	qaws_exact_curve_prepare(NULL, sa, &esa, NULL);
	qaws_exact_curve_prepare(NULL, sb, &esb, NULL);
	TEST_ASSERT(qaws_exact_curve_curve_hits(esa, esb, hits, 8, &count) == QAWS_STATUS_OK && count == 1 && hits[0].kind == QAWS_EXACT_HIT_POINT &&
		hits[0].a_lo == 0.5 && hits[0].b_lo == 1, "segment / segment: the exact point (2, 2) at s = 1/2, r = 1");
	/* a curve against itself: a common component */
	TEST_ASSERT(qaws_exact_curve_curve_hits(ea, ea, hits, 8, &count) == QAWS_STATUS_CERTIFICATION_FAILED,
		"overlapping curves (a common component) are refused");
	qaws_exact_curve_destroy(ea);
	qaws_exact_curve_destroy(eb);
	qaws_exact_curve_destroy(esa);
	qaws_exact_curve_destroy(esb);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
	qaws_curve_destroy(sa);
	qaws_curve_destroy(sb);
}

int test_71_exact_curve_hits_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 71: Certified curve / curve intersections\n");
	test_reference();
	test_splines();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
