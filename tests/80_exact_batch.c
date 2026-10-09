/*
 * Test 80: Certified batched curve / curve intersection
 *
 *   - heightfield x^2 + 2 y^2: contour ellipses (NURBS conics) x gradient
 *     parabolas (polynomials) with two families; every crossing found by the
 *     float batch lies in exactly one certified enclosure, and the counts agree
 *   - the same pairs as qaws_exact_curve_curve_hits over every contour /
 *     gradient pair, with timing
 *   - a duplicated (overlapping) curve in another family: that pair is
 *     uncertified and counted, every other pair still reported
 *   - self-intersections of a looped cubic: (1/2, 3/7) at t = (7 -+ sqrt 21) / 14
 *   - span boxes are sound: random points of every span inside its box
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "qaws_curve_batch.h"
#include "exact/qaws_exact_curve.h"
#include <math.h>
#include <string.h>
#include <time.h>

#define EBT_LEVELS 6
#define EBT_LINES 7
#define EBT_N (EBT_LEVELS + EBT_LINES)

static qaws_curve* ebt_ellipse(double a, double b)
{
	static double const ux[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, uy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	qaws_scalar cps[18], ws[9], kn[12] = { 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 9; i++)
	{
		cps[2 * i] = (qaws_scalar)(a * ux[i]);
		cps[2 * i + 1] = (qaws_scalar)(b * uy[i]);
		ws[i] = (qaws_scalar)(i % 2 ? sqrt(0.5) : 1);
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = cps;
	d.control_point_count = 9;
	d.knots = kn;
	d.knot_count = 12;
	d.weights = ws;
	d.weight_count = 9;
	qaws_curve_create_nurbs(&d, &c);
	return c;
}

static qaws_curve* ebt_parabola(double c, int vertical, double R)
{
	qaws_scalar co[6];
	qaws_polynomial_desc d;
	qaws_curve* out = NULL;
	memset(co, 0, sizeof(co));
	if (vertical)
		co[3] = 1;
	else
	{
		co[2] = 1;
		co[5] = (qaws_scalar)c;
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.coefficients = co;
	d.coefficient_count = 3;
	d.t_min = (qaws_scalar)-R;
	d.t_max = (qaws_scalar)R;
	qaws_curve_create_polynomial(&d, &out);
	return out;
}

static double ebt_now(void)
{
	return (double)clock() / CLOCKS_PER_SEC;
}

static void ebt_scene(qaws_curve** cs, unsigned int* fam)
{
	unsigned int i;
	double R = sqrt(0.5 + 0.75 * (EBT_LEVELS - 1)) * 1.0625;
	for (i = 0; i < EBT_LEVELS; i++)
	{
		double L = 0.5 + 0.75 * i;
		cs[i] = ebt_ellipse(sqrt(L), sqrt(L / 2));
		fam[i] = 0;
	}
	for (i = 0; i < EBT_LINES; i++)
	{
		/* dyadic slopes */
		double c = (double)((int)i - 3) / 4;
		cs[EBT_LEVELS + i] = ebt_parabola(c, i == EBT_LINES - 1, R);
		fam[EBT_LEVELS + i] = 1;
	}
}

static void test_against_float(void)
{
	qaws_curve* cs[EBT_N];
	qaws_exact_curve* ec[EBT_N];
	unsigned int fam[EBT_N], i, j, ne = 0, nf = 0, inside = 0, pairwise = 0, same = 1;
	static qaws_exact_batch_hit he[512];
	static qaws_curve_batch_hit_2d hf[512];
	static qaws_exact_pair hp[64];
	qaws_exact_batch_desc d;
	qaws_curve_batch_desc fd;
	qaws_exact_batch_stats st;
	qaws_exact_desc ed;
	qaws_status s;
	double t0, tb, tp;
	char msg[200];
	ebt_scene(cs, fam);
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -20;
	for (i = 0; i < EBT_N; i++)
		qaws_exact_curve_prepare(&ed, cs[i], &ec[i], NULL);
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_exact_curve const* const*)ec;
	d.curve_count = EBT_N;
	d.families = fam;
	t0 = ebt_now();
	s = qaws_exact_curve_batch_hits(&d, he, 512, &ne, &st);
	tb = ebt_now() - t0;
	memset(&fd, 0, sizeof(fd));
	fd.curves = (qaws_curve const* const*)cs;
	fd.curve_count = EBT_N;
	fd.families = fam;
	qaws_curve_batch_find_intersections_2d(&fd, hf, 512, &nf, NULL);
	/* each float hit in exactly one enclosure of its curve pair */
	for (i = 0; i < nf; i++)
	{
		unsigned int m = 0;
		for (j = 0; j < ne; j++)
		{
			qaws_exact_pair const* p = &he[j].pair;
			/* the exact curves are the quantized ones (2^-20 lattice, 24-bit weights) */
			double ea = 1e-5 * (1 + fabs(hf[i].parameter_a)), eb = 1e-5 * (1 + fabs(hf[i].parameter_b));
			/* a contour (curve_a) is closed on [0, 1]: its start is reported as t = 1 */
			double ta = hf[i].parameter_a < ea ? hf[i].parameter_a + 1 : hf[i].parameter_a;
			if (he[j].curve_a == hf[i].curve_a && he[j].curve_b == hf[i].curve_b
				&& ((hf[i].parameter_a >= p->a_lo - ea && hf[i].parameter_a <= p->a_hi + ea) || (ta >= p->a_lo - ea && ta <= p->a_hi + ea))
				&& hf[i].parameter_b >= p->b_lo - eb && hf[i].parameter_b <= p->b_hi + eb)
				m++;
		}
		inside += m == 1;
	}
	/* the pairwise exact call over every contour / gradient pair */
	t0 = ebt_now();
	for (i = 0; i < EBT_LEVELS; i++)
		for (j = EBT_LEVELS; j < EBT_N; j++)
		{
			unsigned int np = 0, k, m;
			qaws_exact_curve_curve_hits(ec[i], ec[j], hp, 64, &np);
			for (k = 0; k < np; k++)
			{
				int found = 0;
				for (m = 0; m < ne && !found; m++)
					found = he[m].curve_a == i && he[m].curve_b == j && he[m].pair.a_lo == hp[k].a_lo && he[m].pair.a_hi == hp[k].a_hi
						&& he[m].pair.b_lo == hp[k].b_lo && he[m].pair.b_hi == hp[k].b_hi;
				same &= found;
			}
			pairwise += np;
		}
	tp = ebt_now() - t0;
	printf("    %u exact curves: %u spans, %u cells, %u candidate span pairs, %u certified hits in %.3f s; float batch %u hits, %u inside one enclosure; pairwise exact %u hits in %.3f s\n",
		EBT_N, st.span_count, st.cell_count, st.candidate_count, ne, tb, nf, inside, pairwise, tp);
	sprintf(msg, "certified batch: all %u contour / gradient crossings, each float hit in exactly one enclosure", 2 * EBT_LEVELS * EBT_LINES);
	TEST_ASSERT(s == QAWS_STATUS_OK && ne == 2 * EBT_LEVELS * EBT_LINES && nf == ne && inside == nf, msg);
	TEST_ASSERT(pairwise == ne && same, "the same enclosures as the pairwise exact call");
	for (i = 0; i < EBT_N; i++)
	{
		qaws_exact_curve_destroy(ec[i]);
		qaws_curve_destroy(cs[i]);
	}
}

static void test_uncertified(void)
{
	qaws_curve* cs[EBT_N + 1];
	qaws_exact_curve* ec[EBT_N + 1];
	unsigned int fam[EBT_N + 1], i, ne = 0, dup = 0;
	static qaws_exact_batch_hit he[512];
	qaws_exact_batch_desc d;
	qaws_exact_batch_stats st;
	qaws_exact_desc ed;
	qaws_status s;
	ebt_scene(cs, fam);
	/* contour 2 again, as a gradient-family curve: it overlaps contour 2 */
	cs[EBT_N] = ebt_ellipse(sqrt(0.5 + 0.75 * 2), sqrt((0.5 + 0.75 * 2) / 2));
	fam[EBT_N] = 1;
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -20;
	for (i = 0; i <= EBT_N; i++)
		qaws_exact_curve_prepare(&ed, cs[i], &ec[i], NULL);
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_exact_curve const* const*)ec;
	d.curve_count = EBT_N + 1;
	d.families = fam;
	s = qaws_exact_curve_batch_hits(&d, he, 512, &ne, &st);
	for (i = 0; i < ne; i++)
		dup += he[i].curve_b == EBT_N;
	printf("    with an overlapping duplicate: status %d, %u uncertified curve pair(s), %u hits (%u on the duplicate)\n", (int)s, st.uncertified_count, ne, dup);
	/* the duplicate is nested with the other contours and in the gradient family: no hits of its own */
	TEST_ASSERT(s == QAWS_STATUS_CERTIFICATION_FAILED && st.uncertified_count == 1 && ne == 2 * EBT_LEVELS * EBT_LINES
		&& dup == 0, "an overlapping pair is refused and counted; the other pairs are still certified");
	for (i = 0; i <= EBT_N; i++)
	{
		qaws_exact_curve_destroy(ec[i]);
		qaws_curve_destroy(cs[i]);
	}
}

static void test_self(void)
{
	qaws_scalar cp[8] = { 0, 0, 1.5, 1, -0.5, 1, 1, 0 };
	qaws_bezier_desc bd;
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_batch_desc d;
	qaws_exact_batch_hit h[4];
	unsigned int n = 0;
	double t1 = (7 - sqrt(21.0)) / 14, t2 = (7 + sqrt(21.0)) / 14;
	qaws_status s;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D;
	bd.degree = 3;
	bd.control_points = cp;
	bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &c);
	qaws_exact_curve_prepare(NULL, c, &e, NULL);
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_exact_curve const* const*)&e;
	d.curve_count = 1;
	d.flags = QAWS_EXACT_BATCH_SELF;
	s = qaws_exact_curve_batch_hits(&d, h, 4, &n, NULL);
	printf("    looped cubic: %u certified self-intersection(s), t in [%.12f, %.12f] / [%.12f, %.12f]\n", n,
		n ? h[0].pair.a_lo : 0, n ? h[0].pair.a_hi : 0, n ? h[0].pair.b_lo : 0, n ? h[0].pair.b_hi : 0);
	TEST_ASSERT(s == QAWS_STATUS_OK && n == 1 && h[0].curve_a == 0 && h[0].curve_b == 0 && h[0].pair.a_lo <= t1 && t1 <= h[0].pair.a_hi
		&& h[0].pair.b_lo <= t2 && t2 <= h[0].pair.b_hi, "certified self-intersection at t = (7 -+ sqrt 21) / 14");
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
}

/* every span's box holds sampled points of the span */
static void test_boxes(void)
{
	qaws_curve* cs[EBT_N];
	qaws_exact_curve* ec[EBT_N];
	unsigned int fam[EBT_N], i, k, j, out = 0, total = 0;
	qaws_exact_desc ed;
	ebt_scene(cs, fam);
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -20;
	for (i = 0; i < EBT_N; i++)
	{
		qaws_exact_curve_prepare(&ed, cs[i], &ec[i], NULL);
		for (k = 0; k < ec[i]->span_count; k++)
		{
			qaws_exact_span const* sp = &ec[i]->spans[k];
			double lo[3], hi[3];
			qaws_exact_span_box(sp, 2, lo, hi);
			for (j = 0; j <= 32; j++)
			{
				/* T on the parameter lattice inside the span */
				int64_t T = sp->a + (sp->b - sp->a) * (int64_t)j / 32;
				qaws_exact_int num[3], den;
				if (qaws_exact_curve_eval_rational(ec[i], T, 0, num, &den) != QAWS_STATUS_OK)
					continue;
				total++;
				if (qaws_exact_ratio_to_double(&num[0], &den) < lo[0] || qaws_exact_ratio_to_double(&num[0], &den) > hi[0]
					|| qaws_exact_ratio_to_double(&num[1], &den) < lo[1] || qaws_exact_ratio_to_double(&num[1], &den) > hi[1])
					out++;
			}
		}
		qaws_exact_curve_destroy(ec[i]);
		qaws_curve_destroy(cs[i]);
	}
	printf("    span boxes: %u sampled points, %u outside their box\n", total, out);
	TEST_ASSERT(total > 0 && out == 0, "every sampled point of a span lies in its control box");
}

int test_80_exact_batch_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 80: Certified batched curve / curve intersection\n");
	test_against_float();
	test_uncertified();
	test_self();
	test_boxes();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
