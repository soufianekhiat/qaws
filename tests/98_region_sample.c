/*
 * Test 98: uniform sampling of regions bounded by curves (qaws_region_sampler)
 *
 *   Three regions against Mathematica (tests/reference/98_region_sample.wls):
 *   a crescent (a disc minus an overlapping disc: concave, two chords at
 *   some x), a disc with a hole, and a bean of four cubic Beziers with a
 *   dent. For each:
 *     - the area and the area left of x, A(x)
 *     - the map keeps area: A(map(u, v).x) = u * area
 *     - every point lies in the region
 *     - Sobol points: the centroid and second moments
 *     - random points: the distribution of x against A(x) / area
 */

#include "test_common.h"
#include "qaws_region_sample.h"
#include "reference/98_region_sample.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RST_PI 3.14159265358979323846
#define RST_F32 (sizeof(qaws_scalar) == 4)

static qaws_curve* rst_circle(double cx, double cy, double r, int ccw)
{
	qaws_arc_segment g;
	qaws_arc_desc d;
	qaws_curve* c = NULL;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)cx; g.center[1] = (qaws_scalar)cy; g.radius = (qaws_scalar)r;
	g.angle_start = ccw ? 0 : (qaws_scalar)(2 * RST_PI);
	g.angle_end = ccw ? (qaws_scalar)(2 * RST_PI) : 0;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
	qaws_curve_create_arc(&d, &c);
	return c;
}

static qaws_curve* rst_cubic(double const* xy)
{
	qaws_scalar p[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++) p[i] = (qaws_scalar)xy[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p; d.control_point_count = 4;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

/* inside test for the checks: discs (a minus b), or the bean by crossings of a flattened copy */
typedef struct rst_region
{
	char const* name;
	double const* moments;
	double const* left;
	unsigned int nleft;
	double c0[3], c1[3];        /* discs: centre, radius */
	double* poly;               /* bean: flattened */
	unsigned int npoly;
} rst_region;

static int rst_inside(rst_region const* r, double x, double y, double margin)
{
	if (!r->poly)
	{
		double d0 = hypot(x - r->c0[0], y - r->c0[1]), d1 = hypot(x - r->c1[0], y - r->c1[1]);
		return d0 <= r->c0[2] + margin && d1 >= r->c1[2] - margin;
	}
	{
		unsigned int i;
		int in = 0;
		for (i = 0; i < r->npoly; i++)
		{
			double const* a = &r->poly[2 * i];
			double const* b = &r->poly[2 * ((i + 1) % r->npoly)];
			if ((a[1] > y) != (b[1] > y) && x < a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]))
				in = !in;
		}
		return in;
	}
}

static void rst_check(qaws_region_sampler const* s, rst_region const* rg)
{
	double area = (double)qaws_region_sampler_get_area(s), worst = 0, scale = sqrt(rg->moments[0]);
	double atol = RST_F32 ? 2e-4 : 1e-10;
	unsigned int i, outside = 0, n = 1u << 16, nr = 20000;
	qaws_vec2* pts = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	double m[6] = { 0, 0, 0, 0, 0, 0 }, ks = 0;
	char msg[300];
	/* the area and A(x) */
	sprintf(msg, "%s: area %.15f (Mathematica %.15f)", rg->name, area, rg->moments[0]);
	TEST_ASSERT(fabs(area - rg->moments[0]) <= atol * rg->moments[0], msg);
	for (i = 0; i < rg->nleft; i++)
	{
		qaws_scalar a = 0;
		qaws_region_sampler_area_left(s, (qaws_scalar)rg->left[2 * i], &a);
		if (fabs((double)a - rg->left[2 * i + 1]) > worst) worst = fabs((double)a - rg->left[2 * i + 1]);
	}
	sprintf(msg, "%s: area left of x at %u abscissas, worst error %.3g", rg->name, rg->nleft, worst);
	TEST_ASSERT(worst <= atol * rg->moments[0], msg);
	/* the map keeps area */
	worst = 0;
	for (i = 0; i <= 64; i++)
	{
		qaws_vec2 p;
		qaws_scalar a = 0;
		double u = i / 64.0;
		qaws_region_sampler_map(s, (qaws_scalar)u, (qaws_scalar)0.37, &p);
		qaws_region_sampler_area_left(s, p.x, &a);
		if (fabs((double)a - u * area) > worst) worst = fabs((double)a - u * area);
	}
	sprintf(msg, "%s: A(map(u, v).x) = u area at 65 u, worst error %.3g", rg->name, worst);
	TEST_ASSERT(worst <= (RST_F32 ? 1e-3 : 1e-9) * area, msg);
	/* Sobol points: inside, and their moments */
	qaws_region_sampler_generate(s, n, QAWS_REGION_SAMPLE_SOBOL, 7, pts);
	for (i = 0; i < n; i++)
	{
		double x = pts[i].x, y = pts[i].y;
		if (!rst_inside(rg, x, y, RST_F32 ? 1e-4 : 1e-9)) outside++;
		m[0] += 1; m[1] += x; m[2] += y; m[3] += x * x; m[4] += y * y; m[5] += x * y;
	}
	worst = 0;
	for (i = 1; i < 6; i++)
	{
		double e = fabs(m[i] / n - rg->moments[i] / rg->moments[0]) / (i < 3 ? scale : scale * scale);
		if (e > worst) worst = e;
	}
	sprintf(msg, "%s: %u Sobol points, %u outside; centroid (%.5f, %.5f) (Mathematica (%.5f, %.5f)), worst moment error %.2g of the size",
		rg->name, n, outside, m[1] / n, m[2] / n, rg->moments[1] / rg->moments[0], rg->moments[2] / rg->moments[0], worst);
	TEST_ASSERT(outside == 0 && worst < 2e-3, msg);
	/* random points: x against A(x) / area (the largest gap of the empirical distribution) */
	qaws_region_sampler_generate(s, nr, QAWS_REGION_SAMPLE_RANDOM, 12345, pts);
	for (i = 0; i < rg->nleft; i++)
	{
		unsigned int k, below = 0;
		for (k = 0; k < nr; k++) if (pts[k].x <= rg->left[2 * i]) below++;
		if (fabs((double)below / nr - rg->left[2 * i + 1] / rg->moments[0]) > ks) ks = fabs((double)below / nr - rg->left[2 * i + 1] / rg->moments[0]);
	}
	sprintf(msg, "%s: %u random points, x against A(x) / area: largest gap %.4f (KS bound at 1%%: %.4f)", rg->name, nr, ks, 1.63 / sqrt((double)nr));
	TEST_ASSERT(ks < 1.63 / sqrt((double)nr), msg);
	free(pts);
}

static void test_discs(void)
{
	static double const cres[2][3] = { { 0, 0, 1 }, { 0.45, 0.1, 0.7 } }, holed[2][3] = { { 0, 0, 1 }, { 0.3, 0.2, 0.35 } };
	unsigned int k;
	for (k = 0; k < 2; k++)
	{
		double const (*c)[3] = k ? holed : cres;
		qaws_curve* cs[2];
		qaws_curve const* cv[2];
		qaws_path_2d p[2];
		qaws_region_sampler_desc d;
		qaws_region_sampler* s = NULL;
		qaws_status st;
		rst_region rg;
		char msg[200];
		cs[0] = rst_circle(c[0][0], c[0][1], c[0][2], 1);
		cs[1] = rst_circle(c[1][0], c[1][1], c[1][2], 0);
		cv[0] = cs[0]; cv[1] = cs[1];
		p[0].curves = &cv[0]; p[0].curve_count = 1; p[0].closed = 1;
		p[1].curves = &cv[1]; p[1].curve_count = 1; p[1].closed = 1;
		memset(&d, 0, sizeof(d));
		d.paths = p; d.path_count = 2; d.fill_rule = QAWS_FILL_POSITIVE;
		st = qaws_region_sampler_create(&d, &s);
		memset(&rg, 0, sizeof(rg));
		rg.name = k ? "a disc with a hole" : "a crescent";
		rg.moments = k ? ref_holed_moments : ref_crescent_moments;
		rg.left = k ? ref_holed_left : ref_crescent_left;
		rg.nleft = (unsigned int)((k ? sizeof(ref_holed_left) : sizeof(ref_crescent_left)) / sizeof(double) / 2);
		memcpy(rg.c0, c[0], sizeof(rg.c0));
		memcpy(rg.c1, c[1], sizeof(rg.c1));
		sprintf(msg, "%s: sampler created (status %d)", rg.name, (int)st);
		TEST_ASSERT(st == QAWS_STATUS_OK && s, msg);
		if (s) rst_check(s, &rg);
		qaws_region_sampler_destroy(s);
		qaws_curve_destroy(cs[0]);
		qaws_curve_destroy(cs[1]);
	}
}

static void test_bean(void)
{
	static double const segs[4][8] = {
		{ 0, 0, 1, -0.9, 2.4, -0.5, 2.2, 0.6 },
		{ 2.2, 0.6, 2, 1.6, 1.4, 1.1, 1, 1 },
		{ 1, 1, 0.6, 0.9, 0.2, 1.8, -0.5, 1.2 },
		{ -0.5, 1.2, -1.2, 0.6, -0.6, 0.2, 0, 0 },
	};
	qaws_curve* cs[4];
	qaws_curve const* cv[4];
	qaws_path_2d p;
	qaws_region_sampler_desc d;
	qaws_region_sampler* s = NULL;
	qaws_status st;
	rst_region rg;
	unsigned int i, k;
	char msg[200];
	for (i = 0; i < 4; i++) { cs[i] = rst_cubic(segs[i]); cv[i] = cs[i]; }
	p.curves = cv; p.curve_count = 4; p.closed = 1;
	memset(&d, 0, sizeof(d));
	d.paths = &p; d.path_count = 1; d.fill_rule = QAWS_FILL_NON_ZERO;
	st = qaws_region_sampler_create(&d, &s);
	memset(&rg, 0, sizeof(rg));
	rg.name = "a bean of four cubics";
	rg.moments = ref_bean_moments;
	rg.left = ref_bean_left;
	rg.nleft = (unsigned int)(sizeof(ref_bean_left) / sizeof(double) / 2);
	/* the inside test on a fine flattening */
	rg.npoly = 4 * 4096;
	rg.poly = (double*)malloc(sizeof(double) * 2 * rg.npoly);
	for (i = 0; i < 4; i++)
		for (k = 0; k < 4096; k++)
		{
			double t = k / 4096.0, b0 = (1 - t) * (1 - t) * (1 - t), b1 = 3 * t * (1 - t) * (1 - t), b2 = 3 * t * t * (1 - t), b3 = t * t * t;
			rg.poly[2 * (4096 * i + k)] = b0 * segs[i][0] + b1 * segs[i][2] + b2 * segs[i][4] + b3 * segs[i][6];
			rg.poly[2 * (4096 * i + k) + 1] = b0 * segs[i][1] + b1 * segs[i][3] + b2 * segs[i][5] + b3 * segs[i][7];
		}
	sprintf(msg, "%s: sampler created (status %d)", rg.name, (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && s, msg);
	if (s) rst_check(s, &rg);
	qaws_region_sampler_destroy(s);
	free(rg.poly);
	for (i = 0; i < 4; i++) qaws_curve_destroy(cs[i]);
}

int test_98_region_sample_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 98: uniform sampling of regions bounded by curves\n");
	test_discs();
	test_bean();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
