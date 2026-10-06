/*
 * Test 72: Certified curve / curve intersections in 3D
 *
 *   - against Mathematica (tests/reference/72_exact_curve_hits_3d.wls):
 *     pairs in one plane (irrational intersections), pairs in space whose
 *     projections cross without the curves meeting (the exact zero test of
 *     the third coordinate rejects them), pairs meeting at an exact point
 *   - a vertical pair: its xy projections overlap, another plane is used
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "reference/72_exact_curve_hits_3d.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static qaws_curve* e3_curve(int degree, int rational, long long const* cp, long long const* w)
{
	qaws_scalar cps[4 * 3], ws[4];
	qaws_curve* c = NULL;
	int i, n = degree + 1;
	for (i = 0; i < n * 3; i++)
		cps[i] = (qaws_scalar)ldexp((double)cp[i], -10);
	for (i = 0; i < n; i++)
		ws[i] = (qaws_scalar)w[i];
	if (rational)
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
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
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = (unsigned int)degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)n;
		qaws_curve_create_bezier(&d, &c);
	}
	return c;
}

static int e3_inside(double lo, double hi, char const* text)
{
	double v = strtod(text, NULL);
	return lo <= v && v <= hi;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_curve_hits_3d) / sizeof(g_ref_exact_curve_hits_3d[0])), i, total = 0, exact = 0, bad = 0;
	double widest = 0;
	char msg[200];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	for (i = 0; i < n; i++)
	{
		ref_exact_curve_hits_3d const* r = &g_ref_exact_curve_hits_3d[i];
		qaws_curve* ca = e3_curve(r->da, r->ra, r->acp, r->aw);
		qaws_curve* cb = e3_curve(r->db, r->rb, r->bcp, r->bw);
		qaws_exact_curve* ea = NULL;
		qaws_exact_curve* eb = NULL;
		qaws_exact_pair hits[16];
		unsigned int count = 0, k;
		qaws_status st;
		int ok;
		qaws_exact_curve_prepare(&desc, ca, &ea, NULL);
		qaws_exact_curve_prepare(&desc, cb, &eb, NULL);
		st = qaws_exact_curve_curve_hits(ea, eb, hits, 16, &count);
		ok = st == QAWS_STATUS_OK && count == (unsigned int)r->n;
		for (k = 0; ok && k < count; k++)
		{
			ok &= e3_inside(hits[k].a_lo, hits[k].a_hi, r->s[k]) && e3_inside(hits[k].b_lo, hits[k].b_hi, r->r[k]);
			if (r->sdy[k] && r->rdy[k])
				ok &= hits[k].kind == QAWS_EXACT_HIT_POINT;
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
		}
		qaws_exact_curve_destroy(ea);
		qaws_exact_curve_destroy(eb);
		qaws_curve_destroy(ca);
		qaws_curve_destroy(cb);
	}
	printf("    %u 3D pairs: %u intersections (%u exact points), widest enclosure %.1e\n", n, total, exact, widest);
	sprintf(msg, "every 3D intersection found and enclosed, crossings of projections rejected, on all %u Mathematica pairs", n);
	TEST_ASSERT(bad == 0, msg);
	TEST_ASSERT(widest < 1e-12, "both parameters enclosed tightly");
}

static void test_vertical(void)
{
	/* two curves in the vertical plane x = y: their xy projections lie on one line */
	qaws_scalar a3[9] = { 0, 0, 0, 1, 1, 4, 2, 2, 0 }, b3[9] = { 0, 0, 2, 1, 1, -2, 2, 2, 2 };
	qaws_bezier_desc d;
	qaws_curve* ca = NULL;
	qaws_curve* cb = NULL;
	qaws_exact_curve* ea = NULL;
	qaws_exact_curve* eb = NULL;
	qaws_exact_pair hits[8];
	unsigned int count = 0;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 2;
	d.control_points = a3;
	d.control_point_count = 3;
	qaws_curve_create_bezier(&d, &ca);
	d.control_points = b3;
	qaws_curve_create_bezier(&d, &cb);
	qaws_exact_curve_prepare(NULL, ca, &ea, NULL);
	qaws_exact_curve_prepare(NULL, cb, &eb, NULL);
	/* z_a = 4 t (1 - t) * 2, z_b = 2 (1 - t)^2 - 4 t (1 - t) + 2 t^2, same x = y = 2t: equal at 8 t (1 - t) = 2 - 8 t (1 - t): t (1 - t) = 1/8 */
	TEST_ASSERT(qaws_exact_curve_curve_hits(ea, eb, hits, 8, &count) == QAWS_STATUS_OK && count == 2 && hits[0].a_lo <= 0.5 - sqrt(0.125) + 1e-12 &&
		hits[0].a_hi >= 0.5 - sqrt(0.125) - 1e-12 && hits[0].b_lo <= hits[0].a_hi && hits[0].b_hi >= hits[0].a_lo,
		"curves in a vertical plane: the degenerate xy projection is skipped, both crossings found");
	qaws_exact_curve_destroy(ea);
	qaws_exact_curve_destroy(eb);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
}

int test_72_exact_curve_hits_3d_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 72: Certified curve / curve intersections in 3D\n");
	test_reference();
	test_vertical();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
