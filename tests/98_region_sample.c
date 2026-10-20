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
 *     clipped by each cell of a 6 x 6 grid); with the density 2 + x^2 + y,
 *     the mass against Green's theorem and three patterns' cell masses
 *   - the density 1 + x^2 + y over the three regions, against Mathematica:
 *     the mass, the mass left of x, the map keeping mass, every pattern's
 *     weighted centroid
 *   - executors: chunks backwards and on threads give the serial points bit
 *     for bit
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
	qaws_curve_create_arc(&c, &d);
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
	qaws_curve_create_bezier(&c, &d);
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
		qaws_region_sampler_area_left(&a, s, (qaws_scalar)rg->left[2 * i]);
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
		qaws_region_sampler_map(&p, s, (qaws_scalar)u, (qaws_scalar)0.37);
		qaws_region_sampler_area_left(&a, s, p.x);
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
			qaws_region_sampler_generate(pts, s, n, kinds[kk], 7);
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
			qaws_region_sampler_generate(&again, s, 1, kinds[kk], 7);
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
		qaws_region_sampler_generate(pts, s, nb, QAWS_REGION_SAMPLE_BLUE_NOISE, 11);
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
		qaws_region_sampler_generate(pts, s, nb, QAWS_REGION_SAMPLE_RANDOM, 11);
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
	qaws_region_sampler_generate(pts, s, nr, QAWS_REGION_SAMPLE_RANDOM, 12345);
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
		st = qaws_region_sampler_create(&s, &d);
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
	st = qaws_region_sampler_create(&s, &d);
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
	qaws_region_sampler_create(&s, &d);
	area = s ? (double)qaws_region_sampler_get_area(s) : 0;
	pts = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	if (s) qaws_region_sampler_generate(pts, s, n, QAWS_REGION_SAMPLE_SOBOL, 5);
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
	qaws_curve_create_polyline_2d(&c, p, 4, 1);
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
	qaws_curve_create_polyline_2d(&c, p, n, 1);
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
	{ qaws_curve* c = NULL; qaws_curve_create_ellipse_2d(&c, (qaws_scalar)0.2, (qaws_scalar)-0.1, (qaws_scalar)1.3, (qaws_scalar)0.6, (qaws_scalar)0.5); rst_add(&sh[n], c); }
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
		qaws_curve_create_polyline_2d(&c, p, 10, 1);
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
	{ qaws_curve* c = NULL; qaws_curve_create_ellipse_2d(&c, 0, 0, (qaws_scalar)1.5, (qaws_scalar)1.0, 0); rst_add(&sh[n], c); }
	rst_add(&sh[n], rst_circle(-0.7, 0.1, 0.25, 0));
	rst_add(&sh[n], rst_circle(0.1, 0.45, 0.2, 0));
	rst_add(&sh[n], rst_circle(0.6, -0.2, 0.3, 0));
	n++;
	return n;
}

/* the sweep's density 2 + x^2 + y (positive on every shape), and the mass of
   a region bounded by paths by Green: the integral of P dy, P = 2x + x^3/3 + x y */
static qaws_scalar rst_rho2(void* user, qaws_scalar x, qaws_scalar y)
{
	(void)user;
	return (qaws_scalar)(2 + (double)x * x + y);
}

static double rst_green_mass(qaws_clip_result const* r)
{
	static double const gx[8] = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
		0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
	static double const gw[8] = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
		0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
	double sum = 0;
	unsigned int i, k, j, sp, g;
	for (i = 0; r && i < qaws_clip_result_get_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(&p, r, i);
		for (k = 0; k < p.curve_count; k++)
		{
			qaws_curve const* c = p.curves[k];
			qaws_range rg = qaws_curve_get_parameter_range(c);
			unsigned int ns = qaws_curve_get_span_count(c);
			if (ns == 0) ns = 1;
			for (sp = 0; sp < ns; sp++)
				for (j = 0; j < 16; j++)
				{
					double t0 = rg.min_value + (rg.max_value - rg.min_value) * (sp * 16.0 + j) / (16.0 * ns);
					double t1 = rg.min_value + (rg.max_value - rg.min_value) * (sp * 16.0 + j + 1) / (16.0 * ns);
					for (g = 0; g < 8; g++)
					{
						qaws_eval_result_2d e;
						double x, y;
						qaws_curve_evaluate_2d(&e, c, (qaws_scalar)(0.5 * (t0 + t1) + 0.5 * (t1 - t0) * gx[g]), QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
						x = e.position.x; y = e.position.y;
						sum += gw[g] * 0.5 * (t1 - t0) * (2 * x + x * x * x / 3 + x * y) * e.d1.y;
					}
				}
		}
	}
	return sum;
}

static void test_sweep(void)
{
	enum { G = 6 };
	enum { K = 11 };
	/* the unit-square patterns through the map, then over the box, then blue noise */
	static unsigned int const kinds[K] = { QAWS_REGION_SAMPLE_RANDOM, QAWS_REGION_SAMPLE_SOBOL, QAWS_REGION_SAMPLE_STRATIFIED,
		QAWS_REGION_SAMPLE_HALTON, QAWS_REGION_SAMPLE_R2,
		QAWS_REGION_SAMPLE_RANDOM | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_SOBOL | QAWS_REGION_SAMPLE_BOX,
		QAWS_REGION_SAMPLE_STRATIFIED | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_HALTON | QAWS_REGION_SAMPLE_BOX,
		QAWS_REGION_SAMPLE_R2 | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_BLUE_NOISE };
	rst_shape sh[9];
	unsigned int ns = rst_shapes(sh), si, k, i;
	/* the full sweep under QAWS_RST_FULL; the suite runs a fifth of it */
	int full = getenv("QAWS_RST_FULL") != NULL;
	unsigned int m = full ? 16384 : 8192, nmc = full ? 200000 : 40000, nblue = full ? 4096 : 2048;
	qaws_vec2* pts = (qaws_vec2*)malloc(sizeof(qaws_vec2) * m);
	char msg[400];
	printf("    shape                 area (sampler)    area (clip region)  closed form       Monte Carlo (+- 4 sigma)    total variation of the cell areas from the samples:\n");
	printf("                                                                                                             map: random Sobol  strat. Halton R2     box: random Sobol  strat. Halton R2     blue\n");
	for (si = 0; si < ns; si++)
	{
		rst_shape* s = &sh[si];
		qaws_region_sampler_desc d;
		qaws_region_sampler* smp = NULL;
		qaws_clip_result* region = NULL;
		qaws_vec2 lo, hi;
		double area, carea = 0, cells[G * G], mc, sigma, tv[K], box;
		double cmass[G * G], gmass;
		unsigned int inside = 0, outside[K], disagree = 0, out_total = 0;
		int ok_tv = 1;
		unsigned long long st = 99;
		memset(&d, 0, sizeof(d));
		d.paths = s->paths; d.path_count = s->npaths; d.fill_rule = s->rule;
		qaws_region_sampler_create(&smp, &d);
		if (!smp) { TEST_ASSERT(0, s->name); continue; }
		area = (double)qaws_region_sampler_get_area(smp);
		/* the region's area from the united paths (qaws_path_compute_area_2d) */
		qaws_clip_boolean(&region, s->paths, NULL, s->npaths, 0, QAWS_CLIP_UNION, s->rule);
		for (i = 0; region && i < qaws_clip_result_get_path_count(region); i++)
		{
			qaws_path_2d p;
			qaws_scalar a = 0;
			qaws_clip_result_get_path(&p, region, i);
			qaws_path_compute_area_2d(&a, &p);
			carea += a;
		}
		/* the box, and the exact area of the region in each cell of a G x G grid */
		lo.x = lo.y = (qaws_scalar)1e30; hi.x = hi.y = (qaws_scalar)-1e30;
		for (i = 0; region && i < qaws_clip_result_get_path_count(region); i++)
		{
			qaws_path_2d p;
			qaws_vec2 a, b;
			qaws_clip_result_get_path(&p, region, i);
			qaws_path_compute_bounds_2d(&a, &b, &p);
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
			qaws_clip_boolean(&cut, s->paths, &rp, s->npaths, 1, QAWS_CLIP_INTERSECTION, s->rule);
			cells[k] = 0;
			for (i = 0; cut && i < qaws_clip_result_get_path_count(cut); i++)
			{
				qaws_path_2d p;
				qaws_scalar a = 0;
				qaws_clip_result_get_path(&p, cut, i);
				qaws_path_compute_area_2d(&a, &p);
				cells[k] += a;
			}
			cmass[k] = rst_green_mass(cut);
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
			qaws_region_locate_point_2d(&loc, s->paths, s->npaths, s->rule, q, 0);
			if (loc != QAWS_POINT_OUTSIDE) inside++;
			/* the sampler's own inside test agrees, but within rounding of the boundary */
			if (loc != QAWS_POINT_ON)
			{
				int in = 0;
				qaws_region_sampler_contains(&in, smp, q);
				if (in != (loc == QAWS_POINT_INSIDE)) disagree++;
			}
		}
		mc = box * inside / nmc;
		sigma = box * sqrt((double)inside / nmc * (1 - (double)inside / nmc) / nmc);
		/* every pattern: the cell areas from its samples */
		for (k = 0; k < K; k++)
		{
			unsigned int n = kinds[k] == QAWS_REGION_SAMPLE_BLUE_NOISE ? nblue : m, c;
			unsigned int counts[G * G];
			memset(counts, 0, sizeof(counts));
			outside[k] = 0;
			qaws_region_sampler_generate(pts, smp, n, kinds[k], 21);
			for (i = 0; i < n; i++)
			{
				int cx = (int)((pts[i].x - lo.x) / (hi.x - lo.x) * G), cy = (int)((pts[i].y - lo.y) / (hi.y - lo.y) * G);
				qaws_point_location loc = QAWS_POINT_OUTSIDE;
				qaws_region_locate_point_2d(&loc, s->paths, s->npaths, s->rule, pts[i], 0);
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
		printf("    %-20s  %.12f  %.12f    %-16s  %.5f +- %.5f       ", s->name, area, carea, s->closed_form > 0 ? "" : "-", mc, 4 * sigma);
		for (k = 0; k < K; k++) printf(" %.4f", tv[k]);
		printf("\n");
		if (s->closed_form > 0) printf("    %-20s  %38s%.12f\n", "", "", s->closed_form);
		sprintf(msg, "%s: area %.12f, the united region's %.12f, closed form %.12f, Monte Carlo %.5f +- %.5f (4 sigma), the sampler's inside test against winding numbers: %u disagree",
			s->name, area, carea, s->closed_form, mc, 4 * sigma, disagree);
		TEST_ASSERT(fabs(area - carea) <= (RST_F32 ? 1e-4 : 1e-9) * area && (s->closed_form == 0 || fabs(area - s->closed_form) <= (RST_F32 ? 2e-4 : 1e-9) * area)
			&& fabs(area - mc) <= 4 * sigma && disagree == 0, msg);
		/* random (map and box) and blue noise to their noise, the others closer */
		for (k = 0; k < K; k++)
		{
			unsigned int base = kinds[k] & ~QAWS_REGION_SAMPLE_BOX;
			double lim = (base == QAWS_REGION_SAMPLE_RANDOM || base == QAWS_REGION_SAMPLE_BLUE_NOISE) ? (full ? 0.03 : 0.045) : (full ? 0.005 : 0.008);
			if (tv[k] >= lim) ok_tv = 0;
			out_total += outside[k];
		}
		sprintf(msg, "%s: %u points outside; total variation of the cell areas from the samples: map %.4f %.4f %.4f %.4f %.4f, box %.4f %.4f %.4f %.4f %.4f, blue noise %.4f (random, Sobol, stratified, Halton, R2)",
			s->name, out_total, tv[0], tv[1], tv[2], tv[3], tv[4], tv[5], tv[6], tv[7], tv[8], tv[9], tv[10]);
		TEST_ASSERT(out_total == 0 && ok_tv, msg);
		/* the density 2 + x^2 + y: the mass against Green's theorem on the united
		   region, and three patterns' cell masses against the exact ones */
		{
			static unsigned int const dk[3] = { QAWS_REGION_SAMPLE_SOBOL, QAWS_REGION_SAMPLE_STRATIFIED | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_BLUE_NOISE };
			qaws_region_sampler* ds = NULL;
			double dmass, dtv[3];
			unsigned int dout = 0;
			d.density = rst_rho2;
			qaws_region_sampler_create(&ds, &d);
			d.density = NULL;
			gmass = rst_green_mass(region);
			dmass = ds ? (double)qaws_region_sampler_get_mass(ds) : 0;
			for (k = 0; ds && k < 3; k++)
			{
				unsigned int nn = dk[k] == QAWS_REGION_SAMPLE_BLUE_NOISE ? nblue : m / 2, c;
				unsigned int counts[G * G];
				memset(counts, 0, sizeof(counts));
				qaws_region_sampler_generate(pts, ds, nn, dk[k], 23);
				for (i = 0; i < nn; i++)
				{
					int cx = (int)((pts[i].x - lo.x) / (hi.x - lo.x) * G), cy = (int)((pts[i].y - lo.y) / (hi.y - lo.y) * G);
					qaws_point_location loc = QAWS_POINT_OUTSIDE;
					qaws_region_locate_point_2d(&loc, s->paths, s->npaths, s->rule, pts[i], 0);
					if (loc == QAWS_POINT_OUTSIDE) dout++;
					if (cx < 0) cx = 0;
					if (cy < 0) cy = 0;
					if (cx >= G) cx = G - 1;
					if (cy >= G) cy = G - 1;
					counts[cy * G + cx]++;
				}
				dtv[k] = 0;
				for (c = 0; c < G * G; c++)
					dtv[k] += fabs(gmass * counts[c] / nn - cmass[c]);
				dtv[k] *= 0.5 / gmass;
			}
			printf("    %-20s  density 2 + x^2 + y: mass %.12f (Green %.12f); cell masses' total variation: Sobol %.4f, stratified over the box %.4f, blue noise %.4f\n",
				"", dmass, gmass, ds ? dtv[0] : 0, ds ? dtv[1] : 0, ds ? dtv[2] : 0);
			sprintf(msg, "%s, density 2 + x^2 + y: mass %.12f (Green on the region %.12f), %u points outside, total variation of the cell masses from the samples %.4f %.4f %.4f (Sobol, stratified over the box, blue noise)",
				s->name, dmass, gmass, dout, ds ? dtv[0] : 0, ds ? dtv[1] : 0, ds ? dtv[2] : 0);
			TEST_ASSERT(ds && fabs(dmass - gmass) <= (RST_F32 ? 2e-4 : 1e-8) * gmass && dout == 0
				&& dtv[0] < (full ? 0.006 : 0.01) && dtv[1] < (full ? 0.03 : 0.045) && dtv[2] < (full ? 0.04 : 0.06), msg);
			qaws_region_sampler_destroy(ds);
		}
		qaws_clip_result_destroy(region);
		qaws_region_sampler_destroy(smp);
	}
	for (si = 0; si < ns; si++)
		for (i = 0; i < sh[si].ncurves; i++)
			qaws_curve_destroy(sh[si].curves[i]);
	free(pts);
}

/* ------------------------------------------------------------------------ */
/*  Density: rho = 1 + x^2 + y over the three regions                       */
/* ------------------------------------------------------------------------ */

static qaws_scalar rst_rho(void* user, qaws_scalar x, qaws_scalar y)
{
	(void)user;
	return (qaws_scalar)(1 + (double)x * x + y);
}

static void test_density(void)
{
	static unsigned int const kinds[9] = { QAWS_REGION_SAMPLE_RANDOM, QAWS_REGION_SAMPLE_SOBOL, QAWS_REGION_SAMPLE_STRATIFIED,
		QAWS_REGION_SAMPLE_HALTON, QAWS_REGION_SAMPLE_R2, QAWS_REGION_SAMPLE_RANDOM | QAWS_REGION_SAMPLE_BOX,
		QAWS_REGION_SAMPLE_SOBOL | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_STRATIFIED | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_BLUE_NOISE };
	static char const* const knames[9] = { "random", "Sobol", "stratified", "Halton", "R2", "random over the box",
		"Sobol over the box", "stratified over the box", "blue noise" };
	static double const segs[4][8] = {
		{ 0, 0, 1, -0.9, 2.4, -0.5, 2.2, 0.6 }, { 2.2, 0.6, 2, 1.6, 1.4, 1.1, 1, 1 },
		{ 1, 1, 0.6, 0.9, 0.2, 1.8, -0.5, 1.2 }, { -0.5, 1.2, -1.2, 0.6, -0.6, 0.2, 0, 0 } };
	unsigned int r, k, i, n = 16384;
	qaws_vec2* pts = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	char msg[400];
	for (r = 0; r < 3; r++)
	{
		qaws_curve* cs[4];
		qaws_curve const* cv[4];
		qaws_path_2d p[2];
		qaws_region_sampler_desc d;
		qaws_region_sampler* s = NULL;
		double const* m = r == 0 ? ref_crescent_moments : r == 1 ? ref_holed_moments : ref_bean_moments;
		double const* m3 = r == 0 ? ref_crescent_moments3 : r == 1 ? ref_holed_moments3 : ref_bean_moments3;
		double const* ml = r == 0 ? ref_crescent_mass_left : r == 1 ? ref_holed_mass_left : ref_bean_mass_left;
		unsigned int nml = (unsigned int)((r == 2 ? sizeof(ref_bean_mass_left) : sizeof(ref_crescent_mass_left)) / sizeof(double) / 2);
		char const* name = r == 0 ? "crescent" : r == 1 ? "disc with a hole" : "bean";
		/* mass, and the density's centroid: integrals of rho, x rho, y rho */
		double mass = m[0] + m[3] + m[2], cx = (m[1] + m3[0] + m[5]) / mass, cy = (m[2] + m3[1] + m[4]) / mass, size = sqrt(m[0]);
		double got_mass, worst = 0;
		rst_region rg;
		memset(&rg, 0, sizeof(rg));
		if (r < 2)
		{
			double cc[2][3] = { { 0, 0, 1 }, { 0.45, 0.1, 0.7 } };
			if (r == 1) { cc[1][0] = 0.3; cc[1][1] = 0.2; cc[1][2] = 0.35; }
			cs[0] = rst_circle(cc[0][0], cc[0][1], cc[0][2], 1);
			cs[1] = rst_circle(cc[1][0], cc[1][1], cc[1][2], 0);
			for (i = 0; i < 2; i++) { cv[i] = cs[i]; p[i].curves = &cv[i]; p[i].curve_count = 1; p[i].closed = 1; }
			memcpy(rg.c0, cc[0], sizeof(rg.c0));
			memcpy(rg.c1, cc[1], sizeof(rg.c1));
		}
		else
		{
			for (i = 0; i < 4; i++) { cs[i] = rst_cubic(segs[i]); cv[i] = cs[i]; }
			p[0].curves = cv; p[0].curve_count = 4; p[0].closed = 1;
			rg.npoly = 4 * 4096;
			rg.poly = (double*)malloc(sizeof(double) * 2 * rg.npoly);
			for (i = 0; i < 4 * 4096; i++)
			{
				unsigned int sg = i / 4096;
				double t = (i % 4096) / 4096.0, b0 = (1 - t) * (1 - t) * (1 - t), b1 = 3 * t * (1 - t) * (1 - t), b2 = 3 * t * t * (1 - t), b3 = t * t * t;
				rg.poly[2 * i] = b0 * segs[sg][0] + b1 * segs[sg][2] + b2 * segs[sg][4] + b3 * segs[sg][6];
				rg.poly[2 * i + 1] = b0 * segs[sg][1] + b1 * segs[sg][3] + b2 * segs[sg][5] + b3 * segs[sg][7];
			}
		}
		memset(&d, 0, sizeof(d));
		d.paths = p; d.path_count = r < 2 ? 2 : 1; d.fill_rule = r < 2 ? QAWS_FILL_POSITIVE : QAWS_FILL_NON_ZERO;
		d.density = rst_rho;
		qaws_region_sampler_create(&s, &d);
		got_mass = s ? (double)qaws_region_sampler_get_mass(s) : 0;
		for (i = 0; s && i < nml; i++)
		{
			qaws_scalar a = 0;
			qaws_region_sampler_mass_left(&a, s, (qaws_scalar)ml[2 * i]);
			if (fabs((double)a - ml[2 * i + 1]) > worst) worst = fabs((double)a - ml[2 * i + 1]);
		}
		sprintf(msg, "density 1 + x^2 + y over the %s: mass %.12f (Mathematica %.12f), mass left of x at %u abscissas, worst error %.3g",
			name, got_mass, mass, nml, worst);
		TEST_ASSERT(s && fabs(got_mass - mass) <= (RST_F32 ? 2e-4 : 1e-9) * mass && worst <= (RST_F32 ? 2e-4 : 1e-9) * mass, msg);
		/* the map keeps mass */
		worst = 0;
		for (i = 0; s && i <= 32; i++)
		{
			qaws_vec2 q;
			qaws_scalar a = 0;
			qaws_region_sampler_map(&q, s, (qaws_scalar)(i / 32.0), (qaws_scalar)0.61);
			qaws_region_sampler_mass_left(&a, s, q.x);
			if (fabs((double)a - i / 32.0 * mass) > worst) worst = fabs((double)a - i / 32.0 * mass);
		}
		sprintf(msg, "density over the %s: mass left of map(u, v).x = u mass at 33 u, worst error %.3g", name, worst);
		TEST_ASSERT(worst <= (RST_F32 ? 1e-3 : 1e-8) * mass, msg);
		/* every pattern: inside, and the density's centroid */
		for (k = 0; s && k < 9; k++)
		{
			unsigned int base = kinds[k] & ~QAWS_REGION_SAMPLE_BOX, outside = 0, nn = base == QAWS_REGION_SAMPLE_BLUE_NOISE ? 4096 : n;
			double sx = 0, sy = 0, e, lim;
			int exact = !(kinds[k] & QAWS_REGION_SAMPLE_BOX) && base != QAWS_REGION_SAMPLE_RANDOM && base != QAWS_REGION_SAMPLE_BLUE_NOISE;
			qaws_region_sampler_generate(pts, s, nn, kinds[k], 5);
			for (i = 0; i < nn; i++)
			{
				if (!rst_inside(&rg, pts[i].x, pts[i].y, RST_F32 ? 1e-4 : 1e-9)) outside++;
				sx += pts[i].x; sy += pts[i].y;
			}
			e = hypot(sx / nn - cx, sy / nn - cy) / size;
			/* through the map the low-discrepancy patterns keep their order; random,
			   rejection and blue noise (adaptive radii) to their noise */
			lim = exact ? 2e-3 : (base == QAWS_REGION_SAMPLE_BLUE_NOISE ? 3e-2 : 1.5e-2);
			sprintf(msg, "density over the %s, %s: %u points, %u outside, centroid (%.5f, %.5f) (Mathematica (%.5f, %.5f)), error %.2g of the size",
				name, knames[k], nn, outside, sx / nn, sy / nn, cx, cy, e);
			TEST_ASSERT(outside == 0 && e < lim, msg);
		}
		qaws_region_sampler_destroy(s);
		free(rg.poly);
		for (i = 0; i < (r < 2 ? 2u : 4u); i++) qaws_curve_destroy(cs[i]);
	}
	free(pts);
}

/* ------------------------------------------------------------------------ */
/*  Executors: the points do not depend on how the chunks run               */
/* ------------------------------------------------------------------------ */

/* one thread, the chunks backwards in uneven ranges */
static void rst_backwards_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int e = count;
	(void)user;
	while (e > 0)
	{
		unsigned int b = e > 3 ? e - 3 : 0;
		task(ctx, b, e);
		e = b;
	}
}

#define RST_THREADS 4

typedef struct rst_thr
{
	qaws_batch_task_fn task;
	void* ctx;
	unsigned int begin, end;
} rst_thr;

#ifdef _WIN32
#include <windows.h>
#include <process.h>
static unsigned __stdcall rst_thr_main(void* p)
{
	rst_thr* a = (rst_thr*)p;
	a->task(a->ctx, a->begin, a->end);
	return 0;
}
#else
#include <pthread.h>
static void* rst_thr_main(void* p)
{
	rst_thr* a = (rst_thr*)p;
	a->task(a->ctx, a->begin, a->end);
	return NULL;
}
#endif

/* RST_THREADS threads, contiguous ranges */
static void rst_threads_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	rst_thr arg[RST_THREADS];
	unsigned int t;
#ifdef _WIN32
	HANDLE h[RST_THREADS];
#else
	pthread_t h[RST_THREADS];
#endif
	(void)user;
	for (t = 0; t < RST_THREADS; t++)
	{
		arg[t].task = task;
		arg[t].ctx = ctx;
		arg[t].begin = count * t / RST_THREADS;
		arg[t].end = count * (t + 1) / RST_THREADS;
	}
	for (t = 1; t < RST_THREADS; t++)
#ifdef _WIN32
		h[t] = (HANDLE)_beginthreadex(NULL, 0, rst_thr_main, &arg[t], 0, NULL);
#else
		pthread_create(&h[t], NULL, rst_thr_main, &arg[t]);
#endif
	task(ctx, arg[0].begin, arg[0].end);
	for (t = 1; t < RST_THREADS; t++)
	{
#ifdef _WIN32
		WaitForSingleObject(h[t], INFINITE);
		CloseHandle(h[t]);
#else
		pthread_join(h[t], NULL);
#endif
	}
}

static void test_executor(void)
{
	static unsigned int const kinds[9] = { QAWS_REGION_SAMPLE_RANDOM, QAWS_REGION_SAMPLE_SOBOL, QAWS_REGION_SAMPLE_STRATIFIED,
		QAWS_REGION_SAMPLE_HALTON, QAWS_REGION_SAMPLE_R2, QAWS_REGION_SAMPLE_RANDOM | QAWS_REGION_SAMPLE_BOX,
		QAWS_REGION_SAMPLE_SOBOL | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_STRATIFIED | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_BLUE_NOISE };
	unsigned int n = 3000, e, d, k, i, differ = 0, failed = 0;
	qaws_batch_executor ex[2];
	qaws_curve* cs[2];
	qaws_curve const* cv[2];
	qaws_path_2d p[2];
	qaws_vec2* a = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	qaws_vec2* b = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	qaws_vec2* uv = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	char msg[200];
	ex[0].parallel_for = rst_backwards_for; ex[0].user = NULL;
	ex[1].parallel_for = rst_threads_for; ex[1].user = NULL;
	/* the crescent */
	cs[0] = rst_circle(0, 0, 1, 1);
	cs[1] = rst_circle(0.45, 0.1, 0.7, 0);
	for (i = 0; i < 2; i++) { cv[i] = cs[i]; p[i].curves = &cv[i]; p[i].curve_count = 1; p[i].closed = 1; }
	for (i = 0; i < n; i++) { uv[i].x = (qaws_scalar)((i * 0.6180339887498949) - floor(i * 0.6180339887498949)); uv[i].y = (qaws_scalar)((i + 0.5) / n); }
	for (d = 0; d < 2; d++)
	{
		qaws_region_sampler_desc desc;
		qaws_region_sampler* s0 = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.paths = p; desc.path_count = 2; desc.fill_rule = QAWS_FILL_POSITIVE;
		desc.density = d ? rst_rho : NULL;
		if (qaws_region_sampler_create(&s0, &desc) != QAWS_STATUS_OK) { failed++; continue; }
		for (e = 0; e < 2; e++)
		{
			qaws_region_sampler* s1 = NULL;
			desc.executor = &ex[e];
			if (qaws_region_sampler_create(&s1, &desc) != QAWS_STATUS_OK) { failed++; continue; }
			if (qaws_region_sampler_get_mass(s1) != qaws_region_sampler_get_mass(s0)) differ++;
			for (k = 0; k < 9; k++)
			{
				if (qaws_region_sampler_generate(a, s0, n, kinds[k], 11) != QAWS_STATUS_OK ||
					qaws_region_sampler_generate(b, s1, n, kinds[k], 11) != QAWS_STATUS_OK) { failed++; continue; }
				if (memcmp(a, b, sizeof(qaws_vec2) * n)) differ++;
			}
			if (qaws_region_sampler_map_many(a, s0, uv, n) != QAWS_STATUS_OK || qaws_region_sampler_map_many(b, s1, uv, n) != QAWS_STATUS_OK)
				failed++;
			else if (memcmp(a, b, sizeof(qaws_vec2) * n))
				differ++;
			qaws_region_sampler_destroy(s1);
		}
		qaws_region_sampler_destroy(s0);
	}
	sprintf(msg, "executors (chunks backwards, %d threads), uniform and density, 9 kinds and map_many: %u results differ from serial, %u failed",
		RST_THREADS, differ, failed);
	TEST_ASSERT(differ == 0 && failed == 0, msg);
	for (i = 0; i < 2; i++) qaws_curve_destroy(cs[i]);
	free(a); free(b); free(uv);
}

int test_98_region_sample_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 98: uniform sampling of regions bounded by curves\n");
	test_discs();
	test_bean();
	test_xor();
	test_density();
	test_executor();
	test_sweep();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
