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
 *     - every pattern (random, Sobol, stratified, Halton, R2): inside, the
 *       centroid and second moments; blue noise: inside and spaced
 *     - the two discs under even-odd (xor): the sliver gets its share
 *     - random points: the distribution of x against A(x) / area
 *   - a sweep, nine shapes x six patterns: the area against the united
 *     region, closed forms and Monte Carlo by winding numbers; every
 *     pattern's cell areas against the exact cell areas (the region
 *     clipped by each cell of a 6 x 6 grid)
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
		/* outside: still in when within the margin of the boundary */
		for (i = 0; !in && i < r->npoly; i++)
		{
			double const* a = &r->poly[2 * i];
			double const* b = &r->poly[2 * ((i + 1) % r->npoly)];
			double dx = b[0] - a[0], dy = b[1] - a[1], t = ((x - a[0]) * dx + (y - a[1]) * dy) / (dx * dx + dy * dy + 1e-300);
			if (t < 0) t = 0;
			if (t > 1) t = 1;
			if (hypot(x - a[0] - t * dx, y - a[1] - t * dy) <= margin) in = 1;
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
	/* every pattern: inside, and its moments (random to its 1 / sqrt(n),
	   the others closer); the same seed gives the same points */
	{
		static unsigned int const kinds[5] = { QAWS_REGION_SAMPLE_RANDOM, QAWS_REGION_SAMPLE_SOBOL, QAWS_REGION_SAMPLE_STRATIFIED,
			QAWS_REGION_SAMPLE_HALTON, QAWS_REGION_SAMPLE_R2 };
		static char const* const names[5] = { "random", "Sobol", "stratified", "Halton", "R2" };
		unsigned int kk;
		for (kk = 0; kk < 5; kk++)
		{
			qaws_vec2 again;
			memset(m, 0, sizeof(m));
			outside = 0;
			qaws_region_sampler_generate(s, n, kinds[kk], 7, pts);
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
			qaws_region_sampler_generate(s, 1, kinds[kk], 7, &again);
			sprintf(msg, "%s, %s: %u points, %u outside; centroid (%.5f, %.5f) (Mathematica (%.5f, %.5f)), worst moment error %.2g of the size",
				rg->name, names[kk], n, outside, m[1] / n, m[2] / n, rg->moments[1] / rg->moments[0], rg->moments[2] / rg->moments[0], worst);
			TEST_ASSERT(outside == 0 && worst < (kk == 0 ? 1.5e-2 : 2e-3) && (kk == 2 || (again.x == pts[0].x && again.y == pts[0].y)), msg);
		}
	}
	/* blue noise: inside, uniform to the boundary's bias, and spaced: the
	   nearest neighbour never far below the packing radius */
	{
		unsigned int nb = 4096, j, closest = 0;
		double rmax = sqrt(rg->moments[0] / (2 * sqrt(3.0) * nb)), dmin = 1e300, dmin_random = 1e300;
		memset(m, 0, sizeof(m));
		outside = 0;
		qaws_region_sampler_generate(s, nb, QAWS_REGION_SAMPLE_BLUE_NOISE, 11, pts);
		for (i = 0; i < nb; i++)
		{
			double x = pts[i].x, y = pts[i].y;
			if (!rst_inside(rg, x, y, RST_F32 ? 1e-4 : 1e-9)) outside++;
			m[1] += x; m[2] += y;
			for (j = i + 1; j < nb; j++)
			{
				double d = hypot(pts[j].x - x, pts[j].y - y);
				if (d < dmin) { dmin = d; closest = i; }
			}
		}
		(void)closest;
		worst = fabs(m[1] / nb - rg->moments[1] / rg->moments[0]) / scale;
		if (fabs(m[2] / nb - rg->moments[2] / rg->moments[0]) / scale > worst) worst = fabs(m[2] / nb - rg->moments[2] / rg->moments[0]) / scale;
		qaws_region_sampler_generate(s, nb, QAWS_REGION_SAMPLE_RANDOM, 11, pts);
		for (i = 0; i < nb; i++)
			for (j = i + 1; j < nb; j++)
			{
				double d = hypot(pts[j].x - pts[i].x, pts[j].y - pts[i].y);
				if (d < dmin_random) dmin_random = d;
			}
		sprintf(msg, "%s, blue noise: %u points, %u outside, centroid error %.2g of the size, closest pair %.3f of the packing radius (random: %.3f)",
			rg->name, nb, outside, worst, dmin / rmax, dmin_random / rmax);
		TEST_ASSERT(outside == 0 && worst < 1e-2 && dmin > 0.5 * rmax && dmin > 10 * dmin_random, msg);
	}
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

/* the two discs of the crescent under even-odd: their xor, the sliver of
   the small disc outside the big one included, in its share */
static void test_xor(void)
{
	qaws_curve* cs[2];
	qaws_curve const* cv[2];
	qaws_path_2d p[2];
	qaws_region_sampler_desc d;
	qaws_region_sampler* s = NULL;
	qaws_vec2* pts;
	double lens = RST_PI - ref_crescent_moments[0], sliver = 0.49 * RST_PI - lens;
	double want = RST_PI + 0.49 * RST_PI - 2 * lens, area;
	unsigned int i, n = 1u << 16, in_sliver = 0, in_lens = 0;
	char msg[300];
	cs[0] = rst_circle(0, 0, 1, 1);
	cs[1] = rst_circle(0.45, 0.1, 0.7, 1);
	for (i = 0; i < 2; i++) { cv[i] = cs[i]; p[i].curves = &cv[i]; p[i].curve_count = 1; p[i].closed = 1; }
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = 2; d.fill_rule = QAWS_FILL_EVEN_ODD;
	qaws_region_sampler_create(&d, &s);
	area = s ? (double)qaws_region_sampler_get_area(s) : 0;
	pts = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	if (s) qaws_region_sampler_generate(s, n, QAWS_REGION_SAMPLE_SOBOL, 5, pts);
	for (i = 0; s && i < n; i++)
	{
		double d0 = hypot(pts[i].x, pts[i].y), d1 = hypot(pts[i].x - 0.45, pts[i].y - 0.1);
		if (d0 > 1 && d1 < 0.7) in_sliver++;
		if (d0 < 1 - 1e-6 && d1 < 0.7 - 1e-6) in_lens++;
	}
	sprintf(msg, "two discs under even-odd (xor): area %.12f (want %.12f), %.4f of the points in the sliver outside the big disc (want %.4f), %u in their overlap",
		area, want, (double)in_sliver / n, sliver / want, in_lens);
	TEST_ASSERT(s && fabs(area - want) < (RST_F32 ? 2e-4 : 1e-10) * want && fabs((double)in_sliver / n - sliver / want) < 2e-3 && in_lens == 0, msg);
	free(pts);
	qaws_region_sampler_destroy(s);
	qaws_curve_destroy(cs[0]);
	qaws_curve_destroy(cs[1]);
}

/* ------------------------------------------------------------------------ */
/*  The sweep: shapes x patterns, areas from the samples against qaws's     */
/*  exact areas                                                              */
/* ------------------------------------------------------------------------ */

static qaws_curve* clt_rect_curve(double x0, double y0, double x1, double y1)
{
	qaws_scalar p[8];
	qaws_curve* c = NULL;
	p[0] = (qaws_scalar)x0; p[1] = (qaws_scalar)y0; p[2] = (qaws_scalar)x1; p[3] = (qaws_scalar)y0;
	p[4] = (qaws_scalar)x1; p[5] = (qaws_scalar)y1; p[6] = (qaws_scalar)x0; p[7] = (qaws_scalar)y1;
	qaws_curve_create_polyline_2d(p, 4, 1, &c);
	return c;
}

/* the test's own uniform numbers (splitmix64), apart from the library's */
static double rs_test_unit(unsigned long long* st)
{
	unsigned long long z = (*st += 0x9E3779B97F4A7C15ull);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	z ^= z >> 31;
	return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}

typedef struct rst_shape
{
	char const* name;
	qaws_curve* curves[12];
	unsigned int ncurves;
	qaws_curve const* views[12];
	qaws_path_2d paths[12];
	unsigned int npaths;
	qaws_fill_rule rule;
	double closed_form;          /* the area in closed form, or 0 */
} rst_shape;

/* one closed path of one curve */
static void rst_add(rst_shape* sh, qaws_curve* c)
{
	sh->curves[sh->ncurves] = c;
	sh->views[sh->ncurves] = c;
	sh->paths[sh->npaths].curves = &sh->views[sh->ncurves];
	sh->paths[sh->npaths].curve_count = 1;
	sh->paths[sh->npaths].closed = 1;
	sh->ncurves++;
	sh->npaths++;
}

static qaws_curve* rst_polar(unsigned int n, double r0, double a1, unsigned int k1, double a2, unsigned int k2, double rot)
{
	qaws_scalar p[2 * 256];
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * RST_PI * i / n, r = r0 * (1 + a1 * cos(k1 * t + rot) + a2 * sin(k2 * t));
		p[2 * i] = (qaws_scalar)(r * cos(t));
		p[2 * i + 1] = (qaws_scalar)(r * sin(t));
	}
	qaws_curve_create_polyline_2d(p, n, 1, &c);
	return c;
}

static unsigned int rst_shapes(rst_shape* sh)
{
	static double const segs[4][8] = {
		{ 0, 0, 1, -0.9, 2.4, -0.5, 2.2, 0.6 }, { 2.2, 0.6, 2, 1.6, 1.4, 1.1, 1, 1 },
		{ 1, 1, 0.6, 0.9, 0.2, 1.8, -0.5, 1.2 }, { -0.5, 1.2, -1.2, 0.6, -0.6, 0.2, 0, 0 } };
	unsigned int n = 0, i;
	double lens = RST_PI - ref_crescent_moments[0];
	memset(sh, 0, sizeof(rst_shape) * 9);
	/* an ellipse, NURBS, turned */
	sh[n].name = "ellipse (NURBS)"; sh[n].rule = QAWS_FILL_NON_ZERO; sh[n].closed_form = RST_PI * 1.3 * 0.6;
	{ qaws_curve* c = NULL; qaws_curve_create_ellipse_2d((qaws_scalar)0.2, (qaws_scalar)-0.1, (qaws_scalar)1.3, (qaws_scalar)0.6, (qaws_scalar)0.5, &c); rst_add(&sh[n], c); }
	n++;
	sh[n].name = "disc minus disc"; sh[n].rule = QAWS_FILL_POSITIVE; sh[n].closed_form = ref_crescent_moments[0];
	rst_add(&sh[n], rst_circle(0, 0, 1, 1)); rst_add(&sh[n], rst_circle(0.45, 0.1, 0.7, 0));
	n++;
	sh[n].name = "disc xor disc"; sh[n].rule = QAWS_FILL_EVEN_ODD; sh[n].closed_form = 1.49 * RST_PI - 2 * lens;
	rst_add(&sh[n], rst_circle(0, 0, 1, 1)); rst_add(&sh[n], rst_circle(0.45, 0.1, 0.7, 1));
	n++;
	sh[n].name = "disc with a hole"; sh[n].rule = QAWS_FILL_POSITIVE; sh[n].closed_form = RST_PI * (1 - 0.35 * 0.35);
	rst_add(&sh[n], rst_circle(0, 0, 1, 1)); rst_add(&sh[n], rst_circle(0.3, 0.2, 0.35, 0));
	n++;
	sh[n].name = "bean (cubics)"; sh[n].rule = QAWS_FILL_NON_ZERO; sh[n].closed_form = ref_bean_moments[0];
	for (i = 0; i < 4; i++) { sh[n].curves[i] = rst_cubic(segs[i]); sh[n].views[i] = sh[n].curves[i]; }
	sh[n].ncurves = 4;
	sh[n].paths[0].curves = sh[n].views; sh[n].paths[0].curve_count = 4; sh[n].paths[0].closed = 1;
	sh[n].npaths = 1;
	n++;
	/* a star of 10 corners: concave polygon, area 10 r R sin(pi / 5) / 2 */
	sh[n].name = "star (polygon)"; sh[n].rule = QAWS_FILL_NON_ZERO; sh[n].closed_form = 5 * 1.0 * 0.45 * sin(RST_PI / 5);
	{
		qaws_scalar p[20];
		qaws_curve* c = NULL;
		for (i = 0; i < 10; i++)
		{
			double r = (i % 2) ? 0.45 : 1.0, t = RST_PI / 2 + RST_PI * i / 5;
			p[2 * i] = (qaws_scalar)(r * cos(t)); p[2 * i + 1] = (qaws_scalar)(r * sin(t));
		}
		qaws_curve_create_polyline_2d(p, 10, 1, &c);
		rst_add(&sh[n], c);
	}
	n++;
	/* five overlapping circles under non-zero: their union */
	sh[n].name = "five circles' union"; sh[n].rule = QAWS_FILL_NON_ZERO;
	for (i = 0; i < 5; i++)
		rst_add(&sh[n], rst_circle(0.7 * cos(1.3 * i), 0.5 * sin(2.1 * i + 0.4), 0.35 + 0.1 * (i % 3), 1));
	n++;
	sh[n].name = "wavy blob (polygon)"; sh[n].rule = QAWS_FILL_NON_ZERO;
	rst_add(&sh[n], rst_polar(240, 1.0, 0.25, 3, 0.15, 7, 0.3));
	n++;
	/* an ellipse with three holes */
	sh[n].name = "ellipse with 3 holes"; sh[n].rule = QAWS_FILL_POSITIVE; sh[n].closed_form = RST_PI * (1.5 * 1.0 - 0.25 * 0.25 - 0.2 * 0.2 - 0.3 * 0.3);
	{ qaws_curve* c = NULL; qaws_curve_create_ellipse_2d(0, 0, (qaws_scalar)1.5, (qaws_scalar)1.0, 0, &c); rst_add(&sh[n], c); }
	rst_add(&sh[n], rst_circle(-0.7, 0.1, 0.25, 0));
	rst_add(&sh[n], rst_circle(0.1, 0.45, 0.2, 0));
	rst_add(&sh[n], rst_circle(0.6, -0.2, 0.3, 0));
	n++;
	return n;
}

static void test_sweep(void)
{
	enum { G = 6 };
	static unsigned int const kinds[6] = { QAWS_REGION_SAMPLE_RANDOM, QAWS_REGION_SAMPLE_SOBOL, QAWS_REGION_SAMPLE_STRATIFIED,
		QAWS_REGION_SAMPLE_HALTON, QAWS_REGION_SAMPLE_R2, QAWS_REGION_SAMPLE_BLUE_NOISE };
	static char const* const knames[6] = { "random", "Sobol", "stratified", "Halton", "R2", "blue" };
	rst_shape sh[9];
	unsigned int ns = rst_shapes(sh), si, k, i;
	/* the full sweep under QAWS_RST_FULL; the suite runs a fifth of it */
	int full = getenv("QAWS_RST_FULL") != NULL;
	unsigned int m = full ? 16384 : 8192, nmc = full ? 200000 : 40000, nblue = full ? 4096 : 2048;
	qaws_vec2* pts = (qaws_vec2*)malloc(sizeof(qaws_vec2) * m);
	char msg[400];
	printf("    shape                 area (sampler)    area (clip region)  closed form       Monte Carlo (+- 4 sigma)    total variation of the cell areas from the samples:\n");
	printf("                                                                                                             random  Sobol   strat.  Halton  R2      blue\n");
	for (si = 0; si < ns; si++)
	{
		rst_shape* s = &sh[si];
		qaws_region_sampler_desc d;
		qaws_region_sampler* smp = NULL;
		qaws_clip_result* region = NULL;
		qaws_vec2 lo, hi;
		double area, carea = 0, cells[G * G], mc, sigma, tv[6], box;
		unsigned int inside = 0, outside[6];
		unsigned long long st = 99;
		memset(&d, 0, sizeof(d));
		d.paths = s->paths; d.path_count = s->npaths; d.fill_rule = s->rule;
		qaws_region_sampler_create(&d, &smp);
		if (!smp) { TEST_ASSERT(0, s->name); continue; }
		area = (double)qaws_region_sampler_get_area(smp);
		/* the region's area from the united paths (qaws_path_compute_area_2d) */
		qaws_clip_boolean(QAWS_CLIP_UNION, s->rule, s->paths, s->npaths, NULL, 0, &region);
		for (i = 0; region && i < qaws_clip_result_get_path_count(region); i++)
		{
			qaws_path_2d p;
			qaws_scalar a = 0;
			qaws_clip_result_get_path(region, i, &p);
			qaws_path_compute_area_2d(&p, &a);
			carea += a;
		}
		/* the box, and the exact area of the region in each cell of a G x G grid */
		lo.x = lo.y = (qaws_scalar)1e30; hi.x = hi.y = (qaws_scalar)-1e30;
		for (i = 0; region && i < qaws_clip_result_get_path_count(region); i++)
		{
			qaws_path_2d p;
			qaws_vec2 a, b;
			qaws_clip_result_get_path(region, i, &p);
			qaws_path_compute_bounds_2d(&p, &a, &b);
			if (a.x < lo.x) lo.x = a.x;
			if (a.y < lo.y) lo.y = a.y;
			if (b.x > hi.x) hi.x = b.x;
			if (b.y > hi.y) hi.y = b.y;
		}
		box = ((double)hi.x - lo.x) * ((double)hi.y - lo.y);
		for (k = 0; k < G * G; k++)
		{
			double x0 = lo.x + (hi.x - lo.x) * (double)(k % G) / G, x1 = lo.x + (hi.x - lo.x) * (double)(k % G + 1) / G;
			double y0 = lo.y + (hi.y - lo.y) * (double)(k / G) / G, y1 = lo.y + (hi.y - lo.y) * (double)(k / G + 1) / G;
			qaws_curve* rc = clt_rect_curve(x0, y0, x1, y1);
			qaws_curve const* rv = rc;
			qaws_path_2d rp;
			qaws_clip_result* cut = NULL;
			rp.curves = &rv; rp.curve_count = 1; rp.closed = 1;
			qaws_clip_boolean(QAWS_CLIP_INTERSECTION, s->rule, s->paths, s->npaths, &rp, 1, &cut);
			cells[k] = 0;
			for (i = 0; cut && i < qaws_clip_result_get_path_count(cut); i++)
			{
				qaws_path_2d p;
				qaws_scalar a = 0;
				qaws_clip_result_get_path(cut, i, &p);
				qaws_path_compute_area_2d(&p, &a);
				cells[k] += a;
			}
			qaws_clip_result_destroy(cut);
			qaws_curve_destroy(rc);
		}
		/* the area by Monte Carlo over the box, located by winding numbers on the input paths */
		for (i = 0; i < nmc; i++)
		{
			qaws_vec2 q;
			qaws_point_location loc = QAWS_POINT_OUTSIDE;
			q.x = (qaws_scalar)(lo.x + (hi.x - lo.x) * rs_test_unit(&st));
			q.y = (qaws_scalar)(lo.y + (hi.y - lo.y) * rs_test_unit(&st));
			qaws_region_locate_point_2d(s->paths, s->npaths, s->rule, q, 0, &loc);
			if (loc != QAWS_POINT_OUTSIDE) inside++;
		}
		mc = box * inside / nmc;
		sigma = box * sqrt((double)inside / nmc * (1 - (double)inside / nmc) / nmc);
		/* every pattern: the cell areas from its samples */
		for (k = 0; k < 6; k++)
		{
			unsigned int n = kinds[k] == QAWS_REGION_SAMPLE_BLUE_NOISE ? nblue : m, c;
			unsigned int counts[G * G];
			memset(counts, 0, sizeof(counts));
			outside[k] = 0;
			qaws_region_sampler_generate(smp, n, kinds[k], 21, pts);
			for (i = 0; i < n; i++)
			{
				int cx = (int)((pts[i].x - lo.x) / (hi.x - lo.x) * G), cy = (int)((pts[i].y - lo.y) / (hi.y - lo.y) * G);
				qaws_point_location loc = QAWS_POINT_OUTSIDE;
				qaws_region_locate_point_2d(s->paths, s->npaths, s->rule, pts[i], 0, &loc);
				if (loc == QAWS_POINT_OUTSIDE) outside[k]++;
				if (cx < 0) cx = 0;
				if (cy < 0) cy = 0;
				if (cx >= G) cx = G - 1;
				if (cy >= G) cy = G - 1;
				counts[cy * G + cx]++;
			}
			tv[k] = 0;
			for (c = 0; c < G * G; c++)
				tv[k] += fabs(area * counts[c] / n - cells[c]);
			tv[k] *= 0.5 / area;
		}
		printf("    %-20s  %.12f  %.12f    %-16s  %.5f +- %.5f        %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n", s->name, area, carea,
			s->closed_form > 0 ? "" : "-", mc, 4 * sigma, tv[0], tv[1], tv[2], tv[3], tv[4], tv[5]);
		if (s->closed_form > 0) printf("    %-20s  %38s%.12f\n", "", "", s->closed_form);
		sprintf(msg, "%s: area %.12f, the united region's %.12f, closed form %.12f, Monte Carlo %.5f +- %.5f (4 sigma)", s->name, area, carea, s->closed_form, mc, 4 * sigma);
		TEST_ASSERT(fabs(area - carea) <= (RST_F32 ? 1e-4 : 1e-9) * area && (s->closed_form == 0 || fabs(area - s->closed_form) <= (RST_F32 ? 2e-4 : 1e-9) * area)
			&& fabs(area - mc) <= 4 * sigma, msg);
		sprintf(msg, "%s: points outside %u %u %u %u %u %u; total variation of the cell areas from the samples %.4f %.4f %.4f %.4f %.4f %.4f (random, Sobol, stratified, Halton, R2, blue noise)",
			s->name, outside[0], outside[1], outside[2], outside[3], outside[4], outside[5], tv[0], tv[1], tv[2], tv[3], tv[4], tv[5]);
		TEST_ASSERT(outside[0] + outside[1] + outside[2] + outside[3] + outside[4] + outside[5] == 0
			&& tv[0] < (full ? 0.03 : 0.045) && tv[1] < (full ? 0.005 : 0.008) && tv[2] < (full ? 0.005 : 0.008) && tv[3] < (full ? 0.005 : 0.008) && tv[4] < (full ? 0.005 : 0.008) && tv[5] < (full ? 0.03 : 0.045), msg);
		qaws_clip_result_destroy(region);
		qaws_region_sampler_destroy(smp);
	}
	for (si = 0; si < ns; si++)
		for (i = 0; i < sh[si].ncurves; i++)
			qaws_curve_destroy(sh[si].curves[i]);
	free(pts);
}

int test_98_region_sample_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 98: uniform sampling of regions bounded by curves\n");
	test_discs();
	test_bean();
	test_xor();
	test_sweep();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
