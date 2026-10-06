/*
 * Test 70: Certified curve / line and curve / plane intersections
 *
 *   - against Mathematica (tests/reference/70_exact_hits.wls): every root
 *     is found and enclosed, dyadic roots come back as exact points,
 *     others as intervals holding exactly one crossing; a non-dyadic
 *     tangency is reported as uncertifiable; an overlap as one interval
 *   - B-spline polylines: a hit on a knot is reported once, an overlapping
 *     segment merges with the hits at its ends
 *   - root isolation of integer Bernstein polynomials: separated clusters
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "exact/qaws_exact_curve.h"
#include "exact/qaws_exact_roots.h"
#include "reference/70_exact_hits.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static qaws_curve* eh_curve(ref_exact_hits const* r)
{
	qaws_scalar cps[9 * 3], ws[9];
	qaws_curve* c = NULL;
	int i, n = r->degree + 1;
	for (i = 0; i < n * r->dim; i++)
		cps[i] = (qaws_scalar)ldexp((double)r->cp[i], -10);
	for (i = 0; i < n; i++)
		ws[i] = (qaws_scalar)r->w[i];
	if (r->rational)
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
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
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)n;
		qaws_curve_create_bezier(&d, &c);
	}
	return c;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_hits) / sizeof(g_ref_exact_hits[0])), i;
	unsigned int roots = 0, points = 0, crossings = 0, bad_rows = 0, fails = 0, overlaps = 0;
	double widest = 0;
	char msg[200];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	for (i = 0; i < n; i++)
	{
		ref_exact_hits const* r = &g_ref_exact_hits[i];
		qaws_curve* c = eh_curve(r);
		qaws_exact_curve* e = NULL;
		qaws_exact_hit hits[32];
		unsigned int count = 0, k;
		double a[3], b[3];
		qaws_status st;
		int ok = 1;
		for (k = 0; k < 3; k++)
		{
			a[k] = ldexp((double)r->a[k], -10);
			b[k] = r->dim == 2 ? ldexp((double)r->b[k], -10) : (double)r->b[k];
		}
		if (qaws_exact_curve_prepare(&desc, c, &e, NULL) != QAWS_STATUS_OK)
		{
			bad_rows++;
			qaws_curve_destroy(c);
			continue;
		}
		st = r->dim == 2 ? qaws_exact_curve_line_hits(e, a, b, 0, hits, 32, &count) : qaws_exact_curve_plane_hits(e, a, b, 0, hits, 32, &count);
		if (r->expect == 1)
		{
			ok = st == QAWS_STATUS_CERTIFICATION_FAILED;
			fails += ok;
		}
		else if (r->expect == 2)
		{
			ok = st == QAWS_STATUS_OK && count == 1 && hits[0].kind == QAWS_EXACT_HIT_OVERLAP && hits[0].t_lo == 0 && hits[0].t_hi == 1;
			overlaps += ok;
		}
		else
		{
			ok = st == QAWS_STATUS_OK && count == (unsigned int)r->nroots;
			for (k = 0; ok && k < count; k++)
			{
				double rv = strtod(r->root[k], NULL);
				ok &= hits[k].t_lo <= rv && rv <= hits[k].t_hi;
				if (r->dyadic[k])
					ok &= hits[k].kind == QAWS_EXACT_HIT_POINT && hits[k].t_lo == rv && hits[k].t_hi == rv;
				else
				{
					ok &= hits[k].kind == QAWS_EXACT_HIT_CROSSING && hits[k].t_lo < hits[k].t_hi;
					if (hits[k].t_hi - hits[k].t_lo > widest) widest = hits[k].t_hi - hits[k].t_lo;
				}
				points += hits[k].kind == QAWS_EXACT_HIT_POINT;
				crossings += hits[k].kind == QAWS_EXACT_HIT_CROSSING;
			}
			roots += (unsigned int)r->nroots;
		}
		if (!ok)
		{
			bad_rows++;
			printf("      row %u: status %d, %u hits (expected %d roots, case %d)\n", i, (int)st, count, r->nroots, r->expect);
			for (k = 0; k < count && k < (unsigned int)r->nroots; k++)
				printf("        kind %d [%.17g, %.17g] vs %s (dyadic %d, mult %d)\n", (int)hits[k].kind, hits[k].t_lo, hits[k].t_hi, r->root[k],
					r->dyadic[k], r->mult[k]);
		}
		qaws_exact_curve_destroy(e);
		qaws_curve_destroy(c);
	}
	printf("    %u cases: %u roots (%u exact points, %u crossing intervals, widest %.1e), %u tangency refused, %u overlap\n", n, roots, points,
		crossings, widest, fails, overlaps);
	sprintf(msg, "every root found and enclosed on all %u Mathematica cases (degrees 2..8, 2D lines, 3D planes)", n);
	TEST_ASSERT(bad_rows == 0, msg);
	TEST_ASSERT(fails == 1 && overlaps == 1, "a non-dyadic tangency is refused, an overlap reported whole");
	TEST_ASSERT(widest < 1e-15, "crossing intervals shrink below 1e-15");
}

static void test_polyline(void)
{
	/* degree 1 B-spline through (0,0) (1,1) (2,1) (3,0) (4,2), t = 0..4 */
	qaws_scalar cps[10] = { 0, 0, 1, 1, 2, 1, 3, 0, 4, 2 }, knots[7] = { 0, 0, 1, 2, 3, 4, 4 };
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_hit hits[8];
	unsigned int count = 0;
	double y1a[2] = { -5, 1 }, y1b[2] = { 5, 1 }, y0a[2] = { -5, 0 }, y0b[2] = { 5, 0 };
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 1;
	d.control_points = cps;
	d.control_point_count = 5;
	d.knots = knots;
	d.knot_count = 7;
	qaws_curve_create_bspline(&d, &c);
	qaws_exact_curve_prepare(NULL, c, &e, NULL);
	TEST_ASSERT(qaws_exact_curve_line_hits(e, y1a, y1b, 0, hits, 8, &count) == QAWS_STATUS_OK && count == 2 &&
		hits[0].kind == QAWS_EXACT_HIT_OVERLAP && hits[0].t_lo == 1 && hits[0].t_hi == 2 && hits[1].kind == QAWS_EXACT_HIT_POINT &&
		hits[1].t_lo == 3.5, "y = 1: the segment on the line is one overlap [1, 2] (its end hits merged), then t = 3.5");
	TEST_ASSERT(qaws_exact_curve_line_hits(e, y0a, y0b, 0, hits, 8, &count) == QAWS_STATUS_OK && count == 2 &&
		hits[0].kind == QAWS_EXACT_HIT_POINT && hits[0].t_lo == 0 && hits[1].kind == QAWS_EXACT_HIT_POINT && hits[1].t_lo == 3,
		"y = 0: the hit on the knot t = 3 is reported once");
	TEST_ASSERT(qaws_exact_curve_line_hits(e, y0a, y0b, 0, hits, 1, &count) == QAWS_STATUS_BUFFER_TOO_SMALL && count == 1,
		"capacity overflow is reported");
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
}

static void test_isolation(void)
{
	/* (s - 1/3)(s - 1/3 - 2^-40) (s - 3/4): two roots 2^-40 apart and one more */
	qaws_exact_int b[4];
	qaws_exact_root roots[8];
	unsigned int count = 0, i;
	qaws_exact_int p, t;
	/* power form with denominators cleared: 3 * 2^40 * 4 * [(s - 1/3)(s - (2^40 + 3)/(3 * 2^40))(s - 3/4)]
	   = (3 s - 1) (3 * 2^40 s - 2^40 - 3) (4 s - 3) */
	int64_t A = (int64_t)1 << 40;
	int64_t c1[2] = { -1, 3 }, c2[2] = { -(A + 3), 3 * A }, c3[2] = { -3, 4 };
	int64_t bin[4] = { 1, 3, 3, 1 };
	int k, j;
	/* power form of c1 c2 c3 (beyond int64: exact ints) */
	qaws_exact_int q[4], r2[4];
	for (k = 0; k < 4; k++)
		qaws_exact_int_zero(&q[k]);
	for (k = 0; k < 2; k++)
		for (j = 0; j < 2; j++)
		{
			qaws_exact_int_from_i64(&t, c1[k]);
			qaws_exact_int_mul_i64(&t, &t, c2[j]);
			qaws_exact_int_add(&q[k + j], &q[k + j], &t);
		}
	for (k = 0; k < 4; k++)
		qaws_exact_int_zero(&r2[k]);
	for (k = 0; k < 3; k++)
		for (j = 0; j < 2; j++)
		{
			qaws_exact_int_mul_i64(&t, &q[k], c3[j]);
			qaws_exact_int_add(&r2[k + j], &r2[k + j], &t);
		}
	/* power -> Bernstein (cubic): b_i = sum_{j <= i} C(i, j) / C(3, j) a_j, cleared by 3 */
	for (i = 0; i < 4; i++)
	{
		qaws_exact_int_zero(&b[i]);
		for (k = 0; k <= (int)i; k++)
		{
			int64_t cij = (k == 0 || k == (int)i) ? 1 : (i == 2 ? 2 : 3);
			qaws_exact_int_mul_i64(&p, &r2[k], cij * (3 / bin[k]));
			qaws_exact_int_add(&b[i], &b[i], &p);
		}
	}
	TEST_ASSERT(qaws_exact_bernstein_isolate(b, 3, roots, 8, &count) == QAWS_STATUS_OK && count == 3 && !roots[0].exact &&
		!roots[1].exact && roots[0].index + 1 <= roots[1].index && roots[0].depth == roots[1].depth,
		"roots 2^-40 apart are separated into disjoint intervals");
	printf("    close pair isolated at depth %d (2^-40 apart)\n", roots[0].depth);
	TEST_ASSERT(count == 3 && qaws_exact_bernstein_refine(b, 3, &roots[2], 20) == QAWS_STATUS_OK && roots[2].exact && roots[2].index == 3 &&
		roots[2].depth == 2, "refinement lands on the dyadic root 3/4 exactly");
	/* the pair: refining both keeps them apart and each encloses its root */
	TEST_ASSERT(qaws_exact_bernstein_refine(b, 3, &roots[0], 60) == QAWS_STATUS_OK &&
		qaws_exact_bernstein_refine(b, 3, &roots[1], 60) == QAWS_STATUS_OK && roots[0].index < roots[1].index &&
		ldexp((double)roots[0].index, -60) <= 1.0 / 3 && 1.0 / 3 <= ldexp((double)(roots[0].index + 1), -60),
		"refined to 2^-60, the close pair stays separated and 1/3 is enclosed");
}

int test_70_exact_hits_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 70: Certified curve / line and curve / plane intersections\n");
	test_reference();
	test_polyline();
	test_isolation();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
