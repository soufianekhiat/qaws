/*
 * Test 73: Certified self-intersections
 *
 *   - against Mathematica (tests/reference/73_exact_self_hits.wls): cubic
 *     to quintic Bezier and rational Bezier loops in 2D, planar loops and
 *     space curves in 3D, B-splines and NURBS of degree 1..3 (crossings
 *     inside spans, across spans, on knots), a closed B-spline
 *   - a closed NURBS ellipse of three arcs on one conic: no self-intersection
 *     (the arcs touch at their ends); a doubled arc is refused (overlap)
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "reference/73_exact_self_hits.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static qaws_curve* es_curve(ref_exact_self_hits const* r)
{
	qaws_scalar cps[8 * 3], ws[8], kn[16];
	qaws_curve* c = NULL;
	int i;
	for (i = 0; i < r->n * r->dim; i++)
		cps[i] = (qaws_scalar)ldexp((double)r->cp[i], -10);
	for (i = 0; i < r->n; i++)
		ws[i] = (qaws_scalar)r->w[i];
	for (i = 0; i < r->nk; i++)
		kn[i] = (qaws_scalar)r->knots[i];
	if (r->kind == 0)
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)r->n;
		d.weights = ws;
		d.weight_count = (unsigned int)r->n;
		qaws_curve_create_rational_bezier(&d, &c);
	}
	else
	{
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)r->n;
		d.weights = ws;
		d.weight_count = (unsigned int)r->n;
		d.knots = kn;
		d.knot_count = (unsigned int)r->nk;
		qaws_curve_create_nurbs(&d, &c);
	}
	return c;
}

static int es_inside(double lo, double hi, char const* text)
{
	double v = strtod(text, NULL);
	return lo <= v && v <= hi;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_self_hits) / sizeof(g_ref_exact_self_hits[0])), i, total = 0, exact = 0, bad = 0;
	double widest = 0;
	char msg[200];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	for (i = 0; i < n; i++)
	{
		ref_exact_self_hits const* r = &g_ref_exact_self_hits[i];
		qaws_curve* c = es_curve(r);
		qaws_exact_curve* e = NULL;
		qaws_exact_pair hits[32];
		unsigned int count = 0, k;
		qaws_status st;
		int ok;
		qaws_exact_curve_prepare(&desc, c, &e, NULL);
		st = qaws_exact_curve_self_hits(e, hits, 32, &count);
		ok = st == QAWS_STATUS_OK && count == (unsigned int)r->count;
		for (k = 0; ok && k < count; k++)
		{
			ok &= hits[k].a_hi < hits[k].b_lo || (hits[k].a_lo < hits[k].b_lo);
			ok &= es_inside(hits[k].a_lo, hits[k].a_hi, r->a[k]) && es_inside(hits[k].b_lo, hits[k].b_hi, r->b[k]);
			if (r->ady[k] && r->bdy[k])
				ok &= hits[k].kind == QAWS_EXACT_HIT_POINT;
			exact += hits[k].kind == QAWS_EXACT_HIT_POINT;
			if (hits[k].a_hi - hits[k].a_lo > widest) widest = hits[k].a_hi - hits[k].a_lo;
			if (hits[k].b_hi - hits[k].b_lo > widest) widest = hits[k].b_hi - hits[k].b_lo;
		}
		total += (unsigned int)r->count;
		if (!ok)
		{
			bad++;
			printf("      case %u (kind %d, %dD, degree %d%s, %d points): status %d, %u hits, expected %d\n", i, r->kind, r->dim, r->degree,
				r->rational ? " rational" : "", r->n, (int)st, count, r->count);
			for (k = 0; k < count && k < 6; k++)
				printf("        a [%.17g, %.17g] b [%.17g, %.17g]  vs  %s, %s\n", hits[k].a_lo, hits[k].a_hi, hits[k].b_lo, hits[k].b_hi,
					k < (unsigned int)r->count ? r->a[k] : "-", k < (unsigned int)r->count ? r->b[k] : "-");
		}
		qaws_exact_curve_destroy(e);
		qaws_curve_destroy(c);
	}
	printf("    %u curves: %u self-intersections (%u exact points), widest enclosure %.1e\n", n, total, exact, widest);
	sprintf(msg, "every self-intersection found once and enclosed on all %u Mathematica curves", n);
	TEST_ASSERT(bad == 0, msg);
	TEST_ASSERT(widest < 1e-11, "both parameters enclosed tightly");
}

static void test_conic(void)
{
	/* three 120-degree arcs of one conic (an affine image of a circle, weights 2 1 2), closed */
	qaws_scalar cps[14] = { 2, 0, 2, 4, -1, 2, -4, 0, -1, -2, 2, -4, 2, 0 };
	qaws_scalar ws[7] = { 2, 1, 2, 1, 2, 1, 2 }, kn[10] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 };
	qaws_scalar cps2[18] = { 2, 0, 2, 4, -1, 2, -4, 0, -1, -2, 2, -4, 2, 0, 2, 4, -1, 2 };
	qaws_scalar ws2[9] = { 2, 1, 2, 1, 2, 1, 2, 1, 2 }, kn2[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_pair hits[8];
	unsigned int count = 99;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = cps;
	d.control_point_count = 7;
	d.weights = ws;
	d.weight_count = 7;
	d.knots = kn;
	d.knot_count = 10;
	qaws_curve_create_nurbs(&d, &c);
	qaws_exact_curve_prepare(NULL, c, &e, NULL);
	TEST_ASSERT(qaws_exact_curve_self_hits(e, hits, 8, &count) == QAWS_STATUS_OK && count == 0,
		"a closed NURBS conic of three arcs: no self-intersection (arcs on one conic touch at their ends)");
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
	/* four arcs: the first one traced twice */
	d.control_points = cps2;
	d.control_point_count = 9;
	d.weights = ws2;
	d.weight_count = 9;
	d.knots = kn2;
	d.knot_count = 12;
	c = NULL;
	e = NULL;
	qaws_curve_create_nurbs(&d, &c);
	qaws_exact_curve_prepare(NULL, c, &e, NULL);
	TEST_ASSERT(qaws_exact_curve_self_hits(e, hits, 8, &count) == QAWS_STATUS_CERTIFICATION_FAILED,
		"an arc traced twice overlaps itself: refused, not guessed");
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
}

int test_73_exact_self_hits_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 73: Certified self-intersections\n");
	test_reference();
	test_conic();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
