/*
 * Test 96: exact contacts of curves (tangencies and overlaps)
 *
 *   - tangency at a dyadic parameter: a parabola touching a line, one exact
 *     point of multiplicity 2
 *   - tangencies at irrational parameters: a quartic touching a line at
 *     x = +-1/sqrt(2), two enclosed points of multiplicity 2 (the resultant's
 *     square-free factor isolates them)
 *   - an inflection on a line: one crossing of multiplicity 3
 *   - overlaps: collinear segments (same and opposite sense), arcs of one
 *     parabola (inside, and overhanging), with exact ends
 */

#include "test_common.h"
#include "qaws_exact.h"
#include <math.h>
#include <string.h>

static qaws_curve* ct_bezier(unsigned int degree, double const* p)
{
	qaws_scalar cps[10];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * (degree + 1); i++)
		cps[i] = (qaws_scalar)p[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = degree;
	d.control_points = cps;
	d.control_point_count = degree + 1;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

typedef struct ct_case
{
	qaws_curve *a, *b;
	qaws_exact_curve *ea, *eb;
	qaws_exact_pair hits[16];
	unsigned int count;
	qaws_status status;
} ct_case;

static void ct_run(ct_case* k, unsigned int da, double const* pa, unsigned int db, double const* pb)
{
	k->a = ct_bezier(da, pa);
	k->b = ct_bezier(db, pb);
	k->ea = k->eb = NULL;
	qaws_exact_curve_prepare(NULL, k->a, &k->ea, NULL);
	qaws_exact_curve_prepare(NULL, k->b, &k->eb, NULL);
	k->count = 0;
	k->status = qaws_exact_curve_curve_hits(k->ea, k->eb, k->hits, 16, &k->count);
}

static void ct_free(ct_case* k)
{
	qaws_exact_curve_destroy(k->ea);
	qaws_exact_curve_destroy(k->eb);
	qaws_curve_destroy(k->a);
	qaws_curve_destroy(k->b);
}

static void ct_print(char* msg, char const* name, ct_case const* k)
{
	unsigned int i;
	size_t len;
	sprintf(msg, "%s: status %d, %u pairs", name, (int)k->status, k->count);
	for (i = 0; i < k->count && i < 3; i++)
	{
		len = strlen(msg);
		sprintf(msg + len, " [kind %d mult %u a %.9g..%.9g b %.9g..%.9g%s]", (int)k->hits[i].kind, k->hits[i].multiplicity, k->hits[i].a_lo,
			k->hits[i].a_hi, k->hits[i].b_lo, k->hits[i].b_hi, k->hits[i].opposite ? " opposite" : "");
	}
}

static int ct_in(double lo, double hi, double v, double slack)
{
	return lo - slack <= v && v <= hi + slack && hi - lo < 1e-9;
}

int test_96_exact_contacts_main(void)
{
	char msg[600];
	g_pass = 0;
	g_fail = 0;
	printf("Test 96: exact contacts (tangencies, overlaps)\n");
	{
		/* y = 2x - x^2 / ... : the parabola (0,0) (1,2) (2,0) has its top (1, 1) at t = 1/2; the line y = 1 touches it there */
		static double const par[6] = { 0, 0, 1, 2, 2, 0 }, line[4] = { 0, 1, 2, 1 };
		ct_case k;
		ct_run(&k, 2, par, 1, line);
		ct_print(msg, "parabola touching a line at its top", &k);
		TEST_ASSERT(k.status == QAWS_STATUS_OK && k.count == 1 && k.hits[0].kind == QAWS_EXACT_HIT_POINT && k.hits[0].multiplicity == 2 &&
			k.hits[0].a_lo == 0.5 && k.hits[0].b_lo == 0.5, msg);
		ct_free(&k);
	}
	{
		/* y = 3 (x^2 - 1/2)^2 on x = 2t - 1 (dyadic Bernstein form), touching y = 0 at x = +-1/sqrt(2) */
		static double const quartic[10] = { -1, 0.75, -0.5, -2.25, 0, 4.75, 0.5, -2.25, 1, 0.75 }, line[4] = { -2, 0, 2, 0 };
		double r = 1 / sqrt(2.0);
		ct_case k;
		ct_run(&k, 4, quartic, 1, line);
		ct_print(msg, "quartic touching a line twice at irrational parameters", &k);
		TEST_ASSERT(k.status == QAWS_STATUS_OK && k.count == 2 && k.hits[0].kind == QAWS_EXACT_HIT_CROSSING && k.hits[0].multiplicity == 2 &&
			k.hits[1].multiplicity == 2 && ct_in(k.hits[0].a_lo, k.hits[0].a_hi, (1 - r) / 2, 1e-15) && ct_in(k.hits[1].a_lo, k.hits[1].a_hi, (1 + r) / 2, 1e-15) &&
			ct_in(k.hits[0].b_lo, k.hits[0].b_hi, (2 - r) / 4, 1e-15) && ct_in(k.hits[1].b_lo, k.hits[1].b_hi, (2 + r) / 4, 1e-15), msg);
		ct_free(&k);
	}
	{
		/* (0,-1) (0,1) (1,-1) (1,1): y = (2t - 1)^3, an inflection on y = 0 at (1/2, 0) */
		static double const cubic[8] = { 0, -1, 0, 1, 1, -1, 1, 1 }, line[4] = { -1, 0, 2, 0 };
		ct_case k;
		ct_run(&k, 3, cubic, 1, line);
		ct_print(msg, "inflection on a line", &k);
		TEST_ASSERT(k.status == QAWS_STATUS_OK && k.count == 1 && k.hits[0].kind == QAWS_EXACT_HIT_POINT && k.hits[0].multiplicity == 3 &&
			k.hits[0].a_lo == 0.5 && k.hits[0].b_lo == 0.5, msg);
		ct_free(&k);
	}
	{
		static double const s1[4] = { 0, 0, 2, 0 }, s2[4] = { 1, 0, 3, 0 }, s3[4] = { 3, 0, 1, 0 };
		ct_case k;
		ct_run(&k, 1, s1, 1, s2);
		ct_print(msg, "collinear segments overlapping", &k);
		TEST_ASSERT(k.status == QAWS_STATUS_OK && k.count == 1 && k.hits[0].kind == QAWS_EXACT_HIT_OVERLAP && !k.hits[0].opposite &&
			k.hits[0].a_lo == 0.5 && k.hits[0].a_hi == 1 && k.hits[0].b_lo == 0 && k.hits[0].b_hi == 0.5, msg);
		ct_free(&k);
		ct_run(&k, 1, s1, 1, s3);
		ct_print(msg, "collinear segments overlapping in opposite senses", &k);
		TEST_ASSERT(k.status == QAWS_STATUS_OK && k.count == 1 && k.hits[0].kind == QAWS_EXACT_HIT_OVERLAP && k.hits[0].opposite &&
			k.hits[0].a_lo == 0.5 && k.hits[0].a_hi == 1 && k.hits[0].b_lo == 0.5 && k.hits[0].b_hi == 1, msg);
		ct_free(&k);
	}
	{
		/* the parabola (2t, 4t (1 - t)): its arc t in [1/4, 3/4], and its arc t in [1/2, 3/2] */
		static double const par[6] = { 0, 0, 1, 2, 2, 0 }, inner[6] = { 0.5, 0.75, 1, 1.25, 1.5, 0.75 }, over[6] = { 1, 1, 2, 1, 3, -3 };
		ct_case k;
		ct_run(&k, 2, par, 2, inner);
		ct_print(msg, "an arc inside the parabola", &k);
		TEST_ASSERT(k.status == QAWS_STATUS_OK && k.count == 1 && k.hits[0].kind == QAWS_EXACT_HIT_OVERLAP && !k.hits[0].opposite &&
			k.hits[0].a_lo == 0.25 && k.hits[0].a_hi == 0.75 && k.hits[0].b_lo == 0 && k.hits[0].b_hi == 1, msg);
		ct_free(&k);
		ct_run(&k, 2, par, 2, over);
		ct_print(msg, "an arc overhanging the parabola's end", &k);
		TEST_ASSERT(k.status == QAWS_STATUS_OK && k.count == 1 && k.hits[0].kind == QAWS_EXACT_HIT_OVERLAP && !k.hits[0].opposite &&
			k.hits[0].a_lo == 0.5 && k.hits[0].a_hi == 1 && k.hits[0].b_lo == 0 && k.hits[0].b_hi == 0.5, msg);
		ct_free(&k);
	}
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
