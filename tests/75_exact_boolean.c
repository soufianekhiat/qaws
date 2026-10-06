/*
 * Test 75: Certified 2D Boolean operations
 *
 *   - two overlapping squares (crossings at exact points): union,
 *     intersection, difference; areas 7, 1, 3 exactly (shoelace over the
 *     pieces' vertices)
 *   - two overlapping conics (irrational crossings): the area of every
 *     result, sampled along its pieces, against an independent grid of
 *     certified winding numbers
 *   - disjoint and nested regions: whole boundaries kept or dropped
 *   - a non-simple boundary is refused
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "qaws_boolean_2d.h"
#include <math.h>
#include <string.h>

static qaws_curve* bo_polygon(double const* xy, unsigned int n)
{
	/* closed degree-1 B-spline through n points (the first repeated at the end) */
	qaws_scalar cps[16 * 2], kn[20];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		cps[2 * i] = (qaws_scalar)xy[2 * i];
		cps[2 * i + 1] = (qaws_scalar)xy[2 * i + 1];
	}
	cps[2 * n] = cps[0];
	cps[2 * n + 1] = cps[1];
	kn[0] = 0;
	for (i = 1; i <= n; i++)
		kn[i] = (qaws_scalar)(i - 1);
	kn[n + 1] = (qaws_scalar)n;
	kn[n + 2] = (qaws_scalar)n;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 1;
	d.control_points = cps;
	d.control_point_count = n + 1;
	d.knots = kn;
	d.knot_count = n + 3;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static qaws_curve* bo_conic(double dx, double dy)
{
	/* three arcs of one conic (an affine image of a circle), closed */
	qaws_scalar cps[14] = { 2, 0, 2, 4, -1, 2, -4, 0, -1, -2, 2, -4, 2, 0 };
	qaws_scalar ws[7] = { 2, 1, 2, 1, 2, 1, 2 }, kn[10] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 7; i++)
	{
		cps[2 * i] += (qaws_scalar)dx;
		cps[2 * i + 1] += (qaws_scalar)dy;
	}
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
	return c;
}

/* Area of the result loops, sampled along each piece (polygons: their vertices, exactly). */
static double bo_area(qaws_exact_curve const* const* curves, qaws_exact_piece const* pieces, qaws_exact_loop const* loops, unsigned int nloops,
	unsigned int samples)
{
	double area = 0;
	unsigned int l, k, i;
	for (l = 0; l < nloops; l++)
	{
		double first[2] = { 0, 0 }, prev[2] = { 0, 0 };
		int have = 0;
		for (k = loops[l].first; k < loops[l].first + loops[l].count; k++)
		{
			qaws_exact_piece const* p = &pieces[k];
			double t0 = 0.5 * (p->t0_lo + p->t0_hi), t1 = 0.5 * (p->t1_lo + p->t1_hi);
			qaws_exact_curve const* c = curves[p->region];
			double span;
			/* a closed curve: a piece may wrap through the closing point */
			{
				double ts, te;
				unsigned int sc = qaws_exact_curve_span_count(c);
				qaws_exact_curve_span_bezier(c, 0, NULL, &ts, NULL, NULL, NULL);
				qaws_exact_curve_span_bezier(c, sc - 1, NULL, NULL, &te, NULL, NULL);
				span = te - ts;
				if (!p->reversed && t1 <= t0) t1 += span;
				if (p->reversed && t1 >= t0) t1 -= span;
				/* polygons (samples == 0): the piece start, then every knot (integer) strictly inside */
				unsigned int nsamp = samples ? samples : (unsigned int)fabs(ceil(fabs(t1 - t0))) + 2;
				for (i = 0; i < nsamp; i++)
				{
					double t, q[2];
					if (samples)
						t = t0 + (t1 - t0) * i / samples;
					else
					{
						/* knots strictly between t0 and t1 along the traversal */
						t = i == 0 ? t0 : (t1 > t0 ? floor(t0) + i : ceil(t0) - i);
						if (i > 0 && (t1 > t0 ? (t <= t0 || t >= t1) : (t >= t0 || t <= t1)))
							continue;
					}
					while (t >= te) t -= span;
					while (t < ts) t += span;
					qaws_exact_curve_evaluate(c, t, 0, q, NULL);
					if (have)
						area += prev[0] * q[1] - q[0] * prev[1];
					else
					{
						first[0] = q[0];
						first[1] = q[1];
						have = 1;
					}
					prev[0] = q[0];
					prev[1] = q[1];
				}
			}
		}
		area += prev[0] * first[1] - first[0] * prev[1];
	}
	return 0.5 * fabs(area);
}

static void test_squares(void)
{
	double sa[8] = { 0, 0, 2, 0, 2, 2, 0, 2 }, sb[8] = { 1, 1, 3, 1, 3, 3, 1, 3 };
	static double const expected[3] = { 7, 1, 3 };
	static char const* const names[3] = { "union", "intersection", "difference" };
	qaws_curve* ca = bo_polygon(sa, 4);
	qaws_curve* cb = bo_polygon(sb, 4);
	qaws_exact_curve* ea = NULL;
	qaws_exact_curve* eb = NULL;
	qaws_exact_curve const* curves[2];
	qaws_exact_piece pieces[32];
	qaws_exact_loop loops[8];
	unsigned int np = 0, nl = 0, op;
	char msg[160];
	qaws_exact_curve_prepare(NULL, ca, &ea, NULL);
	qaws_exact_curve_prepare(NULL, cb, &eb, NULL);
	curves[0] = ea;
	curves[1] = eb;
	for (op = 0; op < 3; op++)
	{
		qaws_status st = qaws_exact_boolean_2d(ea, eb, op, pieces, 32, &np, loops, 8, &nl);
		/* polygons: the area from the pieces' vertices, exactly */
		double area = bo_area(curves, pieces, loops, nl, 0);
		printf("    squares %-12s: status %d, %u loops, %u pieces, area %.12f\n", names[op], (int)st, nl, np, area);
		sprintf(msg, "squares %s: one loop of exact pieces, area %g", names[op], expected[op]);
		TEST_ASSERT(st == QAWS_STATUS_OK && nl == 1 && fabs(area - expected[op]) < 1e-9, msg);
	}
	qaws_exact_curve_destroy(ea);
	qaws_exact_curve_destroy(eb);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
}

static void test_conics(void)
{
	static char const* const names[3] = { "union", "intersection", "difference" };
	qaws_curve* ca = bo_conic(0, 0);
	qaws_curve* cb = bo_conic(1.5, 0.75);
	qaws_exact_curve* ea = NULL;
	qaws_exact_curve* eb = NULL;
	qaws_exact_curve const* curves[2];
	qaws_exact_piece pieces[32];
	qaws_exact_loop loops[8];
	unsigned int np = 0, nl = 0, op, i, j, cnt[3] = { 0, 0, 0 }, G = 240;
	double x0 = -4.5, x1 = 4.5, y0 = -4.5, y1 = 5.0, cell;
	char msg[160];
	qaws_exact_curve_prepare(NULL, ca, &ea, NULL);
	qaws_exact_curve_prepare(NULL, cb, &eb, NULL);
	curves[0] = ea;
	curves[1] = eb;
	/* the grid estimate: certified windings at cell centers */
	cell = (x1 - x0) / G;
	for (j = 0; j < G; j++)
		for (i = 0; i < G; i++)
		{
			double p[2] = { x0 + (i + 0.5) * cell, y0 + (j + 0.5) * (y1 - y0) / G };
			int wa = 0, wb = 0, ia, ib;
			qaws_exact_winding_2d(&curves[0], 1, p, &wa);
			qaws_exact_winding_2d(&curves[1], 1, p, &wb);
			ia = wa != 0;
			ib = wb != 0;
			cnt[0] += ia || ib;
			cnt[1] += ia && ib;
			cnt[2] += ia && !ib;
		}
	for (op = 0; op < 3; op++)
	{
		qaws_status st = qaws_exact_boolean_2d(ea, eb, op, pieces, 32, &np, loops, 8, &nl);
		double area = bo_area(curves, pieces, loops, nl, 2000), grid = cnt[op] * cell * (y1 - y0) / G;
		printf("    conics  %-12s: status %d, %u loops, %u pieces, area %.5f, certified-winding grid %.5f\n", names[op], (int)st, nl, np, area, grid);
		sprintf(msg, "conics %s: one loop, area agrees with the winding grid", names[op]);
		TEST_ASSERT(st == QAWS_STATUS_OK && nl == 1 && fabs(area - grid) < 0.01 * grid + 0.05, msg);
	}
	qaws_exact_curve_destroy(ea);
	qaws_exact_curve_destroy(eb);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
}

static void test_disjoint_nested(void)
{
	double big[8] = { -4, -4, 4, -4, 4, 4, -4, 4 }, small[8] = { -1, -1, 1, -1, 1, 1, -1, 1 }, far[8] = { 10, 10, 12, 10, 12, 12, 10, 12 };
	double bowtie[8] = { 0, 0, 2, 2, 2, 0, 0, 2 };
	qaws_curve* cbig = bo_polygon(big, 4);
	qaws_curve* csmall = bo_polygon(small, 4);
	qaws_curve* cfar = bo_polygon(far, 4);
	qaws_curve* cbow = bo_polygon(bowtie, 4);
	qaws_exact_curve* ebig = NULL;
	qaws_exact_curve* esmall = NULL;
	qaws_exact_curve* efar = NULL;
	qaws_exact_curve* ebow = NULL;
	qaws_exact_piece pieces[16];
	qaws_exact_loop loops[8];
	unsigned int np = 0, nl = 0, n0, n1, n2;
	qaws_exact_curve_prepare(NULL, cbig, &ebig, NULL);
	qaws_exact_curve_prepare(NULL, csmall, &esmall, NULL);
	qaws_exact_curve_prepare(NULL, cfar, &efar, NULL);
	qaws_exact_curve_prepare(NULL, cbow, &ebow, NULL);
	qaws_exact_boolean_2d(ebig, efar, QAWS_BOOLEAN_UNION, pieces, 16, &np, loops, 8, &n0);
	qaws_exact_boolean_2d(ebig, efar, QAWS_BOOLEAN_INTERSECTION, pieces, 16, &np, loops, 8, &n1);
	qaws_exact_boolean_2d(ebig, efar, QAWS_BOOLEAN_DIFFERENCE, pieces, 16, &np, loops, 8, &n2);
	TEST_ASSERT(n0 == 2 && n1 == 0 && n2 == 1, "disjoint regions: union both, intersection none, difference the first");
	qaws_exact_boolean_2d(ebig, esmall, QAWS_BOOLEAN_UNION, pieces, 16, &np, loops, 8, &n0);
	qaws_exact_boolean_2d(ebig, esmall, QAWS_BOOLEAN_INTERSECTION, pieces, 16, &np, loops, 8, &n1);
	qaws_exact_boolean_2d(ebig, esmall, QAWS_BOOLEAN_DIFFERENCE, pieces, 16, &np, loops, 8, &n2);
	TEST_ASSERT(n0 == 1 && n1 == 1 && n2 == 2 && pieces[0].region == 0 && pieces[1].region == 1,
		"nested regions: union the outer, intersection the inner, difference both (a hole)");
	TEST_ASSERT(qaws_exact_boolean_2d(ebow, esmall, QAWS_BOOLEAN_UNION, pieces, 16, &np, loops, 8, &nl) == QAWS_STATUS_CERTIFICATION_FAILED,
		"a self-intersecting boundary is refused");
	qaws_exact_curve_destroy(ebig);
	qaws_exact_curve_destroy(esmall);
	qaws_exact_curve_destroy(efar);
	qaws_exact_curve_destroy(ebow);
	qaws_curve_destroy(cbig);
	qaws_curve_destroy(csmall);
	qaws_curve_destroy(cfar);
	qaws_curve_destroy(cbow);
}

int test_75_exact_boolean_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 75: Certified 2D Boolean operations\n");
	test_squares();
	test_conics();
	test_disjoint_nested();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
