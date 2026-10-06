/*
 * Test 66: Certified winding numbers of exact curve loops
 *
 *   - a loop of four rational quadratic conics around points 2^-28 lattice
 *     units inside and outside the boundary, against Mathematica (exact
 *     algebraic ray crossings, tests/reference/66_exact_winding.wls)
 *   - random polygons (degree-1 Beziers on the integer lattice) against an
 *     independent crossing count built on qaws_exact_orient2d
 *   - points exactly on the boundary are reported (CERTIFICATION_FAILED),
 *     loops that do not close are refused
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "reference/66_exact_winding.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

static qaws_exact_desc ew_desc(void)
{
	qaws_exact_desc d;
	qaws_exact_desc_default(&d);
	d.space_exp2 = 0;
	return d;
}

static qaws_exact_curve* ew_piece(qaws_scalar const* cps, qaws_scalar const* ws, unsigned int degree)
{
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_desc desc = ew_desc();
	if (ws)
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = degree;
		d.control_points = cps;
		d.control_point_count = degree + 1;
		d.weights = ws;
		d.weight_count = degree + 1;
		qaws_curve_create_rational_bezier(&d, &c);
	}
	else
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = degree;
		d.control_points = cps;
		d.control_point_count = degree + 1;
		qaws_curve_create_bezier(&d, &c);
	}
	qaws_exact_curve_prepare(&desc, c, &e, NULL);
	qaws_curve_destroy(c);
	return e;
}

static void test_conics(void)
{
	qaws_exact_curve* pieces[4];
	unsigned int n = (unsigned int)(sizeof(g_ref_wind) / sizeof(g_ref_wind[0])), i, k;
	int ok = 1;
	char msg[200];
	for (k = 0; k < 4; k++)
	{
		qaws_scalar cps[6], ws[3];
		for (i = 0; i < 3; i++)
		{
			cps[2 * i] = (qaws_scalar)g_ref_wind_pieces[k][i][0];
			cps[2 * i + 1] = (qaws_scalar)g_ref_wind_pieces[k][i][1];
			ws[i] = (qaws_scalar)g_ref_wind_weights[i];
		}
		pieces[k] = ew_piece(cps, ws, 2);
	}
	for (i = 0; i < n; i++)
	{
		double p[2];
		int w = -99;
		p[0] = ldexp((double)g_ref_wind[i].mx, g_ref_wind[i].ex);
		p[1] = ldexp((double)g_ref_wind[i].my, g_ref_wind[i].ey);
		ok &= qaws_exact_winding_2d((qaws_exact_curve const* const*)pieces, 4, p, &w) == QAWS_STATUS_OK && w == g_ref_wind[i].winding;
	}
	sprintf(msg, "conic loop: winding equals Mathematica for all %u points (2^-28 lattice units from the boundary)", n);
	TEST_ASSERT(ok, msg);
	{
		/* a point exactly on the boundary: the start of a piece */
		double p[2];
		int w;
		p[0] = (double)g_ref_wind_pieces[1][0][0];
		p[1] = (double)g_ref_wind_pieces[1][0][1];
		TEST_ASSERT(qaws_exact_winding_2d((qaws_exact_curve const* const*)pieces, 4, p, &w) == QAWS_STATUS_CERTIFICATION_FAILED,
			"a point on the conic loop is reported, not guessed");
		/* an open chain is refused */
		TEST_ASSERT(qaws_exact_winding_2d((qaws_exact_curve const* const*)pieces, 3, p, &w) == QAWS_STATUS_INVALID_ARGUMENT,
			"a loop that does not close is refused");
	}
	for (k = 0; k < 4; k++)
		qaws_exact_curve_destroy(pieces[k]);
}

static uint64_t g_ew_state = 0x2545F4914F6CDD1Dull;
static int ew_rand(int lo, int hi)
{
	g_ew_state ^= g_ew_state << 13;
	g_ew_state ^= g_ew_state >> 7;
	g_ew_state ^= g_ew_state << 17;
	return lo + (int)(g_ew_state % (uint64_t)(hi - lo + 1));
}

/* Crossing count of the ray x > p_x with the polygon by exact orientation:
   +1 for an edge going up through (y_a <= p_y < y_b) with p on its left,
   -1 going down with p on its right; 2 when p lies on an edge. */
static int ew_polygon(double const* v, int n, double const* p)
{
	int i, w = 0;
	for (i = 0; i < n; i++)
	{
		double const* a = &v[2 * i];
		double const* b = &v[2 * ((i + 1) % n)];
		qaws_exact_sign s;
		qaws_exact_orient2d(a, b, p, &s, NULL);
		if (s == QAWS_EXACT_ZERO && fmin(a[0], b[0]) <= p[0] && p[0] <= fmax(a[0], b[0]) && fmin(a[1], b[1]) <= p[1] &&
		    p[1] <= fmax(a[1], b[1]))
			return 1000;
		if (a[1] <= p[1] && p[1] < b[1] && s == QAWS_EXACT_POSITIVE)
			w++;
		else if (b[1] <= p[1] && p[1] < a[1] && s == QAWS_EXACT_NEGATIVE)
			w--;
	}
	return w;
}

static void test_polygons(void)
{
	int poly, ok = 1, on_edge = 0, cases = 0;
	char msg[200];
	for (poly = 0; poly < 40; poly++)
	{
		int n = ew_rand(3, 12), i, q;
		double v[2 * 12];
		qaws_exact_curve* pieces[12];
		/* a star-shaped polygon around the origin (may self-touch rarely) */
		for (i = 0; i < n; i++)
		{
			double a = 2 * 3.14159265358979 * (i + 0.4 * ew_rand(0, 100) / 100.0) / n, r = ew_rand(50, 1000);
			v[2 * i] = floor(r * cos(a));
			v[2 * i + 1] = floor(r * sin(a));
		}
		for (i = 0; i < n; i++)
		{
			qaws_scalar cps[4];
			cps[0] = (qaws_scalar)v[2 * i];
			cps[1] = (qaws_scalar)v[2 * i + 1];
			cps[2] = (qaws_scalar)v[2 * ((i + 1) % n)];
			cps[3] = (qaws_scalar)v[2 * ((i + 1) % n) + 1];
			pieces[i] = ew_piece(cps, NULL, 1);
		}
		for (q = 0; q < 200; q++)
		{
			double p[2];
			int w = -99, expect;
			qaws_status st;
			if (q % 10 == 0)
			{
				/* a lattice point exactly on an edge midpoint (when integral) or a vertex */
				int e = ew_rand(0, n - 1);
				p[0] = v[2 * e];
				p[1] = v[2 * e + 1];
			}
			else
			{
				/* lattice points and half-lattice points */
				p[0] = ew_rand(-1100, 1100) + 0.5 * ew_rand(0, 1);
				p[1] = ew_rand(-1100, 1100) + 0.5 * ew_rand(0, 1);
			}
			expect = ew_polygon(v, n, p);
			st = qaws_exact_winding_2d((qaws_exact_curve const* const*)pieces, (unsigned int)n, p, &w);
			cases++;
			if (expect == 1000)
			{
				on_edge++;
				ok &= st == QAWS_STATUS_CERTIFICATION_FAILED;
			}
			else
				ok &= st == QAWS_STATUS_OK && w == expect;
		}
		for (i = 0; i < n; i++)
			qaws_exact_curve_destroy(pieces[i]);
	}
	sprintf(msg, "%d points against 40 random polygons agree with exact orientation (%d on the boundary reported)", cases, on_edge);
	TEST_ASSERT(ok, msg);
}

int test_66_exact_winding_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 66: Certified winding numbers of exact curve loops\n");
	test_conics();
	test_polygons();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
