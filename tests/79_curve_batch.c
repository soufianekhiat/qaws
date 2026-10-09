/*
 * Test 79: Batched curve / curve intersection
 *
 *   - heightfield h = x^2 + 2 y^2: contours h = L are ellipses (NURBS conics),
 *     gradient lines are the parabolas y = c x^2 and the y axis. With the two
 *     families, every contour / gradient pair crosses twice at
 *     x^2 = (sqrt(1 + 8 c^2 L) - 1) / (4 c^2), y = c x^2 (closed form); the
 *     gradient lines' common point at the critical point is skipped
 *   - the same curves lifted onto the surface z = h(x, y) in 3D: the same
 *     crossings at height L
 *   - random cubic Beziers without families: every pair as found by the
 *     pairwise qaws_curve_find_intersections_2d
 *   - self-intersections: a looped cubic against qaws_curve_find_self_intersections_2d
 *   - timing: the batch against all pairs of the pairwise call
 */

#include "test_common.h"
#include "qaws_curve_batch.h"
#include <math.h>
#include <string.h>
#include <time.h>

#define CBT_PI 3.14159265358979323846

static qaws_curve* cbt_ellipse(double a, double b, int lift, double z)
{
	static double const ux[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, uy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	qaws_scalar cps[27], ws[9], kn[12] = { 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i, dim = lift ? 3 : 2;
	for (i = 0; i < 9; i++)
	{
		cps[i * dim] = (qaws_scalar)(a * ux[i]);
		cps[i * dim + 1] = (qaws_scalar)(b * uy[i]);
		if (lift)
			cps[i * dim + 2] = (qaws_scalar)z;
		ws[i] = (qaws_scalar)(i % 2 ? sqrt(0.5) : 1);
	}
	memset(&d, 0, sizeof(d));
	d.dimension = lift ? QAWS_DIMENSION_3D : QAWS_DIMENSION_2D;
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

/* gradient line y = c x^2 for x in [-R, R], or the y axis when vertical;
   lifted: z = x^2 + 2 y^2 */
static qaws_curve* cbt_gradient(double c, int vertical, double R, int lift)
{
	qaws_scalar co[15];
	qaws_polynomial_desc d;
	qaws_curve* out = NULL;
	unsigned int dim = lift ? 3 : 2, deg = lift ? 4 : 2;
	memset(co, 0, sizeof(co));
	if (vertical)
	{
		co[1 * dim + 1] = 1;                       /* y = t */
		if (lift)
			co[2 * dim + 2] = 2;                   /* z = 2 t^2 */
	}
	else
	{
		co[1 * dim + 0] = 1;                       /* x = t */
		co[2 * dim + 1] = (qaws_scalar)c;          /* y = c t^2 */
		if (lift)
		{
			co[2 * dim + 2] = 1;                   /* z = t^2 + 2 c^2 t^4 */
			co[4 * dim + 2] = (qaws_scalar)(2 * c * c);
		}
	}
	memset(&d, 0, sizeof(d));
	d.dimension = lift ? QAWS_DIMENSION_3D : QAWS_DIMENSION_2D;
	d.degree = deg;
	d.coefficients = co;
	d.coefficient_count = deg + 1;
	d.t_min = (qaws_scalar)-R;
	d.t_max = (qaws_scalar)R;
	qaws_curve_create_polynomial(&d, &out);
	return out;
}

#define CBT_LEVELS 12
#define CBT_LINES 13

static double cbt_level(unsigned int i) { return 0.25 + 0.5 * i; }
static double cbt_slope(unsigned int j) { return j == CBT_LINES - 1 ? 0 : tan(CBT_PI * ((double)j / (CBT_LINES - 1) - 0.5) * 0.9) * 0.6; }

/* expected crossing of contour L and gradient line j: (x >= 0, y) and (-x, y) */
static void cbt_expect(double L, unsigned int j, double* x, double* y)
{
	if (j == CBT_LINES - 1)
	{
		*x = 0;
		*y = sqrt(L / 2);
	}
	else
	{
		double c = cbt_slope(j);
		double x2 = fabs(c) < 1e-12 ? L : (sqrt(1 + 8 * c * c * L) - 1) / (4 * c * c);
		*x = sqrt(x2);
		*y = c * x2;
	}
}

static void test_heightfield(int lift)
{
	qaws_curve* cs[CBT_LEVELS + CBT_LINES];
	unsigned int fam[CBT_LEVELS + CBT_LINES], i, j, n = 0, expected = 0, matched = 0;
	static qaws_curve_batch_hit_2d h2[1024];
	static qaws_curve_batch_hit_3d h3[1024];
	qaws_curve_batch_desc d;
	qaws_curve_batch_stats st;
	qaws_status s;
	double worst = 0, R = sqrt(cbt_level(CBT_LEVELS - 1)) * 1.05;
	char msg[200];
	for (i = 0; i < CBT_LEVELS; i++)
	{
		double L = cbt_level(i);
		cs[i] = cbt_ellipse(sqrt(L), sqrt(L / 2), lift, L);
		fam[i] = 0;
	}
	for (j = 0; j < CBT_LINES; j++)
	{
		cs[CBT_LEVELS + j] = cbt_gradient(cbt_slope(j), j == CBT_LINES - 1, j == CBT_LINES - 1 ? R : R, lift);
		fam[CBT_LEVELS + j] = 1;
	}
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_curve const* const*)cs;
	d.curve_count = CBT_LEVELS + CBT_LINES;
	d.families = fam;
	s = lift ? qaws_curve_batch_find_intersections_3d(&d, h3, 1024, &n, &st)
		: qaws_curve_batch_find_intersections_2d(&d, h2, 1024, &n, &st);
	/* every expected crossing must be found once */
	for (i = 0; i < CBT_LEVELS; i++)
		for (j = 0; j < CBT_LINES; j++)
		{
			double L = cbt_level(i), x, y;
			int side;
			cbt_expect(L, j, &x, &y);
			for (side = 0; side < 2; side++)
			{
				double ex = j == CBT_LINES - 1 ? 0 : (side ? -x : x), ey = j == CBT_LINES - 1 ? (side ? -y : y) : y, best = 1e30;
				unsigned int k, hits = 0;
				expected++;
				for (k = 0; k < n && k < 1024; k++)
				{
					unsigned int a = lift ? h3[k].curve_a : h2[k].curve_a, b = lift ? h3[k].curve_b : h2[k].curve_b;
					double px = lift ? h3[k].position.x : h2[k].position.x, py = lift ? h3[k].position.y : h2[k].position.y;
					double dd;
					if (a != i || b != CBT_LEVELS + j)
						continue;
					dd = hypot(px - ex, py - ey);
					if (lift)
						dd = fmax(dd, fabs(h3[k].position.z - L));
					if (dd < 1e-3)
						hits++;
					if (dd < best)
						best = dd;
				}
				if (hits == 1)
					matched++;
				if (best > worst)
					worst = best;
			}
		}
	printf("    %s: %u curves, %u segments, %u cells, %u candidate segment pairs, %u refined, %u hits; worst error %.2e\n",
		lift ? "3D lifted" : "2D", d.curve_count, st.segment_count, st.cell_count, st.candidate_count, st.newton_count, n, worst);
	sprintf(msg, "%s heightfield: every one of the %u contour / gradient crossings found once", lift ? "3D" : "2D", expected);
	TEST_ASSERT(s == QAWS_STATUS_OK && matched == expected && n == expected, msg);
	TEST_ASSERT(worst < (QAWS_SCALAR_IS_FLOAT ? 1e-3 : 1e-9), "crossings agree with the closed form");
	for (i = 0; i < CBT_LEVELS + CBT_LINES; i++)
		qaws_curve_destroy(cs[i]);
}

static unsigned long long g_cbt_state = 0x9E3779B97F4A7C15ull;
static double cbt_rand(void)
{
	g_cbt_state ^= g_cbt_state << 13; g_cbt_state ^= g_cbt_state >> 7; g_cbt_state ^= g_cbt_state << 17;
	return (double)(g_cbt_state >> 11) / 9007199254740992.0;
}

static qaws_curve* cbt_cubic(void)
{
	qaws_scalar cp[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int k;
	for (k = 0; k < 8; k++)
		cp[k] = (qaws_scalar)cbt_rand();
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = cp;
	d.control_point_count = 4;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

#define CBT_RANDOM 24

static void test_against_pairwise(void)
{
	qaws_curve* cs[CBT_RANDOM];
	static qaws_curve_batch_hit_2d hb[4096];
	qaws_intersection_2d hp[64];
	qaws_curve_batch_desc d;
	unsigned int i, j, nb = 0, total = 0, found = 0;
	char msg[200];
	for (i = 0; i < CBT_RANDOM; i++)
		cs[i] = cbt_cubic();
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_curve const* const*)cs;
	d.curve_count = CBT_RANDOM;
	qaws_curve_batch_find_intersections_2d(&d, hb, 4096, &nb, NULL);
	for (i = 0; i < CBT_RANDOM; i++)
		for (j = i + 1; j < CBT_RANDOM; j++)
		{
			unsigned int np = 0, k, m;
			qaws_curve_find_intersections_2d(cs[i], cs[j], hp, 64, &np);
			for (k = 0; k < np && k < 64; k++)
			{
				total++;
				for (m = 0; m < nb; m++)
					if (hb[m].curve_a == i && hb[m].curve_b == j
						&& hypot(hb[m].position.x - hp[k].position.x, hb[m].position.y - hp[k].position.y) < 1e-3)
					{
						found++;
						break;
					}
			}
		}
	printf("    %u random cubics: batch %u hits, pairwise %u, pairwise hits found by the batch %u\n", CBT_RANDOM, nb, total, found);
	sprintf(msg, "batch finds every one of the %u pairwise intersections", total);
	TEST_ASSERT(found == total && nb >= total, msg);
	for (i = 0; i < CBT_RANDOM; i++)
		qaws_curve_destroy(cs[i]);
}

static void test_self(void)
{
	qaws_scalar cp[8] = { 0, 0, 1.5, 1, -0.5, 1, 1, 0 };
	qaws_bezier_desc bd;
	qaws_curve* c = NULL;
	qaws_curve* cs[2];
	qaws_curve_batch_desc d;
	qaws_curve_batch_hit_2d hb[8];
	qaws_intersection_2d hp[8];
	unsigned int nb = 0, np = 0, nn = 0;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D;
	bd.degree = 3;
	bd.control_points = cp;
	bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &c);
	cs[0] = c;
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_curve const* const*)cs;
	d.curve_count = 1;
	qaws_curve_batch_find_intersections_2d(&d, hb, 8, &nn, NULL);
	d.flags = QAWS_CURVE_BATCH_SELF;
	qaws_curve_batch_find_intersections_2d(&d, hb, 8, &nb, NULL);
	qaws_curve_find_self_intersections_2d(c, hp, 8, &np);
	printf("    looped cubic: batch self %u at (%.6f, %.6f), t %.6f / %.6f; pairwise self %u\n", nb,
		nb ? hb[0].position.x : 0, nb ? hb[0].position.y : 0, nb ? hb[0].parameter_a : 0, nb ? hb[0].parameter_b : 0, np);
	TEST_ASSERT(nn == 0, "without the self flag a single curve has no hits");
	TEST_ASSERT(nb == 1 && np == 1 && hypot(hb[0].position.x - hp[0].position.x, hb[0].position.y - hp[0].position.y) < 1e-3
		&& fabs(hb[0].position.x - 0.5) < 1e-6 && fabs(hb[0].position.y - 3.0 / 7) < 1e-6
		&& fabs(hb[0].parameter_a - (7 - sqrt(21.0)) / 14) < 1e-6 && fabs(hb[0].parameter_b - (7 + sqrt(21.0)) / 14) < 1e-6,
		"self-intersection of the looped cubic, once, at (1/2, 3/7), t = (7 -+ sqrt 21) / 14");
	qaws_curve_destroy(c);
}

static double cbt_now(void)
{
	return (double)clock() / CLOCKS_PER_SEC;
}

/* the heightfield with n contours and n gradient lines: batch against all pairs */
static void test_timing(void)
{
	unsigned int sizes[3] = { 16, 48, 128 }, z;
	for (z = 0; z < 3; z++)
	{
		unsigned int n = sizes[z], i, j, nb = 0, np = 0;
		qaws_curve** cs = (qaws_curve**)malloc(2 * n * sizeof(qaws_curve*));
		unsigned int* fam = (unsigned int*)malloc(2 * n * sizeof(unsigned int));
		qaws_curve_batch_hit_2d* hb = (qaws_curve_batch_hit_2d*)malloc(4 * n * n * sizeof(qaws_curve_batch_hit_2d));
		qaws_intersection_2d hp[16];
		qaws_curve_batch_desc d;
		qaws_curve_batch_stats st;
		double t0, tb, tp, R = sqrt(0.25 + 0.5 * (n - 1)) * 1.05;
		char msg[200];
		for (i = 0; i < n; i++)
		{
			double L = 0.25 + 0.5 * i, c = tan(CBT_PI * ((double)i / n - 0.5) * 0.9) * 0.6;
			cs[i] = cbt_ellipse(sqrt(L), sqrt(L / 2), 0, 0);
			cs[n + i] = cbt_gradient(c, 0, R, 0);
			fam[i] = 0;
			fam[n + i] = 1;
		}
		memset(&d, 0, sizeof(d));
		d.curves = (qaws_curve const* const*)cs;
		d.curve_count = 2 * n;
		d.families = fam;
		t0 = cbt_now();
		qaws_curve_batch_find_intersections_2d(&d, hb, 4 * n * n, &nb, &st);
		tb = cbt_now() - t0;
		t0 = cbt_now();
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				unsigned int k = 0;
				qaws_curve_find_intersections_2d(cs[i], cs[n + j], hp, 16, &k);
				np += k;
			}
		tp = cbt_now() - t0;
		printf("    %u contours x %u gradient lines: batch %u hits in %.3f s (%u segments, %u candidates), pairwise %u hits in %.3f s: x%.0f\n",
			n, n, nb, tb, st.segment_count, st.candidate_count, np, tp, tp / (tb > 1e-4 ? tb : 1e-4));
		sprintf(msg, "%u x %u: batch finds the %u crossings", n, n, 2 * n * n);
		TEST_ASSERT(nb == 2 * n * n, msg);
		for (i = 0; i < 2 * n; i++)
			qaws_curve_destroy(cs[i]);
		free(cs);
		free(fam);
		free(hb);
	}
}

#define CBT_SET_N 64
#define CBT_FRAMES 6
#define CBT_SET_FLAT 0.004   /* one flatness for both paths */

/* fixed contours prepared once, queried against gradient lines that change */
static void test_sets(void)
{
	qaws_curve* contours[CBT_SET_N];
	qaws_curve* grads[CBT_SET_N];
	qaws_curve* all[2 * CBT_SET_N];
	unsigned int fam[2 * CBT_SET_N], i, f, same = 1, total = 0, self_same;
	qaws_curve_set* cset = NULL;
	qaws_curve_batch_desc d;
	qaws_curve_batch_hit_2d* hs = (qaws_curve_batch_hit_2d*)malloc(8192 * sizeof(qaws_curve_batch_hit_2d));
	qaws_curve_batch_hit_2d* ho = (qaws_curve_batch_hit_2d*)malloc(8192 * sizeof(qaws_curve_batch_hit_2d));
	double R = sqrt(0.25 + 0.5 * (CBT_SET_N - 1)) * 1.05, t0, t_set = 0, t_once = 0;
	char msg[200];
	for (i = 0; i < CBT_SET_N; i++)
	{
		double L = 0.25 + 0.5 * i;
		contours[i] = cbt_ellipse(sqrt(L), sqrt(L / 2), 0, 0);
	}
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_curve const* const*)contours;
	d.curve_count = CBT_SET_N;
	d.flatness = (qaws_scalar)CBT_SET_FLAT;
	t0 = cbt_now();
	qaws_curve_set_create(&d, &cset);
	t_set += cbt_now() - t0;
	for (f = 0; f < CBT_FRAMES; f++)
	{
		qaws_curve_set* gset = NULL;
		qaws_curve_batch_desc gd;
		unsigned int ns = 0, no = 0, k;
		for (i = 0; i < CBT_SET_N; i++)
			grads[i] = cbt_gradient(tan(CBT_PI * ((i + 0.37 * f) / CBT_SET_N - 0.5) * 0.9) * 0.6, 0, R, 0);
		memset(&gd, 0, sizeof(gd));
		gd.curves = (qaws_curve const* const*)grads;
		gd.curve_count = CBT_SET_N;
		gd.flatness = (qaws_scalar)CBT_SET_FLAT;
		/* the prepared contours against this frame's gradient lines */
		t0 = cbt_now();
		qaws_curve_set_create(&gd, &gset);
		qaws_curve_set_find_intersections_2d(cset, gset, hs, 8192, &ns, NULL);
		t_set += cbt_now() - t0;
		/* the one-shot call on everything, two families */
		for (i = 0; i < CBT_SET_N; i++)
		{
			all[i] = contours[i];
			all[CBT_SET_N + i] = grads[i];
			fam[i] = 0;
			fam[CBT_SET_N + i] = 1;
		}
		memset(&gd, 0, sizeof(gd));
		gd.curves = (qaws_curve const* const*)all;
		gd.curve_count = 2 * CBT_SET_N;
		gd.families = fam;
		gd.flatness = (qaws_scalar)CBT_SET_FLAT;
		t0 = cbt_now();
		qaws_curve_batch_find_intersections_2d(&gd, ho, 8192, &no, NULL);
		t_once += cbt_now() - t0;
		same &= ns == no;
		for (k = 0; k < ns && k < no && same; k++)
			same &= hs[k].curve_a == ho[k].curve_a && hs[k].curve_b + CBT_SET_N == ho[k].curve_b
				&& fabs(hs[k].parameter_a - ho[k].parameter_a) < 1e-9 && fabs(hs[k].parameter_b - ho[k].parameter_b) < 1e-9;
		total += ns;
		qaws_curve_set_destroy(gset);
		for (i = 0; i < CBT_SET_N; i++)
			qaws_curve_destroy(grads[i]);
	}
	/* a set with itself: the one-shot call on the same curves */
	{
		unsigned int ns = 0, no = 0, k;
		qaws_curve_set_find_intersections_2d(cset, NULL, hs, 8192, &ns, NULL);
		qaws_curve_batch_find_intersections_2d(&d, ho, 8192, &no, NULL);
		self_same = ns == no;
		for (k = 0; k < ns && k < no; k++)
			self_same &= hs[k].curve_a == ho[k].curve_a && hs[k].curve_b == ho[k].curve_b;
	}
	printf("    prepared contours (%u segments) x %u frames of %u gradient lines: %u hits, %.3f s; one-shot calls %.3f s\n",
		qaws_curve_set_get_segment_count(cset), CBT_FRAMES, CBT_SET_N, total, t_set, t_once);
	sprintf(msg, "prepared set queries equal the one-shot batch on all %u frames (%u hits)", CBT_FRAMES, total);
	TEST_ASSERT(same && total == CBT_FRAMES * 2 * CBT_SET_N * CBT_SET_N, msg);
	TEST_ASSERT(self_same, "a set with itself equals the one-shot call");
	qaws_curve_set_destroy(cset);
	for (i = 0; i < CBT_SET_N; i++)
		qaws_curve_destroy(contours[i]);
	free(hs);
	free(ho);
}

int test_79_curve_batch_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 79: Batched curve / curve intersection\n");
	test_heightfield(0);
	test_heightfield(1);
	test_against_pairwise();
	test_self();
	test_timing();
	test_sets();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
