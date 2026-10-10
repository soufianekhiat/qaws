/*
 * Test 90: offsets of curves at a constant delta (qaws_offset), after
 * "Fast GPU stroke expansion" (Levien and Uguray, HPG 2024)
 *
 *   - the offset's cusps (1 + delta k = 0) against Mathematica
 *     (tests/reference/90_offset_curves.wls): a parabola toward its focus,
 *     a cubic S on both sides; each cusp is a joint of the chain
 *   - the chain lies on the true offset or on the evolute (closed forms for
 *     the parabola), G1 everywhere but at the cusps of the offset and of
 *     the evolute
 *   - regions grown and shrunk past their radii of curvature, and an open
 *     curve with round ends, against a grid oracle: a point is in the
 *     offset iff its distance to the region (or curve) is below delta
 *   - few cubics: a NURBS circle grown by 1
 */

#include "test_common.h"
#include "internal/qaws_internal_offset.h"
#include "reference/90_offset_curves.h"
#include <math.h>
#include <string.h>

#define OCT_PI 3.14159265358979323846
#define OCT_F32 (sizeof(qaws_scalar) == 4)

static qaws_curve* oct_bezier(double const* xy, unsigned int n)
{
	qaws_scalar p[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = n - 1; d.control_points = p; d.control_point_count = n;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

/* the joints of a cubic chain (control points 0, 3, 6, ...) and how many of
   the inner ones are corners */
static unsigned int oct_chain(qaws_curve const* c, double* cp, unsigned int cap, unsigned int* corners)
{
	unsigned int n = 0, i;
	qaws_scalar* tmp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * cap);
	qaws_curve_get_control_points(c, tmp, cap, &n);
	for (i = 0; i < 2 * n; i++) cp[i] = tmp[i];
	free(tmp);
	*corners = 0;
	for (i = 3; i + 3 < n; i += 3)
	{
		double a[2] = { cp[2 * i] - cp[2 * i - 2], cp[2 * i + 1] - cp[2 * i - 1] };
		double b[2] = { cp[2 * i + 2] - cp[2 * i], cp[2 * i + 3] - cp[2 * i + 1] };
		double la = hypot(a[0], a[1]), lb = hypot(b[0], b[1]);
		if (la > 0 && lb > 0 && (fabs(a[0] * b[1] - a[1] * b[0]) > (OCT_F32 ? 1e-3 : 1e-6) * la * lb || a[0] * b[0] + a[1] * b[1] < 0))
			(*corners)++;
	}
	return n;
}

static double oct_nearest_joint(double const* cp, unsigned int n, double x, double y)
{
	unsigned int i;
	double best = 1e300;
	for (i = 0; i < n; i += 3)
	{
		double d = hypot(cp[2 * i] - x, cp[2 * i + 1] - y);
		if (d < best) best = d;
	}
	return best;
}

/* the parabola (x, x^2): its offset 1 toward the focus, its evolute */
static void oct_par_offset(double x, double* p)
{
	double s = sqrt(1 + 4 * x * x);
	p[0] = x - 2 * x / s;
	p[1] = x * x + 1 / s;
}

static void oct_par_evolute(double x, double* p)
{
	p[0] = -4 * x * x * x;
	p[1] = 0.5 + 3 * x * x;
}

/* the distance to a parametric curve on x in [-1, 1]: every local minimum of
   a dense sampling refined (near a cusp two branches run side by side) */
static double oct_dist_to(void (*f)(double, double*), double px, double py)
{
	enum { N = 4000 };
	static double d[N + 1];
	unsigned int i, k;
	double best = 1e300, q[2];
	for (i = 0; i <= N; i++)
	{
		f(-1 + 2.0 * i / N, q);
		d[i] = hypot(q[0] - px, q[1] - py);
	}
	for (i = 0; i <= N; i++)
	{
		double a, b;
		if ((i > 0 && d[i - 1] < d[i]) || (i < N && d[i + 1] < d[i])) continue;
		a = -1 + 2.0 * (i ? i - 1 : 0) / N;
		b = -1 + 2.0 * (i < N ? i + 1 : N) / N;
		for (k = 0; k < 80; k++)
		{
			double m1 = a + (b - a) / 3, m2 = b - (b - a) / 3, q1[2], q2[2];
			f(m1, q1); f(m2, q2);
			if (hypot(q1[0] - px, q1[1] - py) < hypot(q2[0] - px, q2[1] - py)) b = m2; else a = m1;
		}
		f(0.5 * (a + b), q);
		if (hypot(q[0] - px, q[1] - py) < best) best = hypot(q[0] - px, q[1] - py);
	}
	return best;
}

static void test_parabola(void)
{
	static double const xy[6] = { -1, 1, 0, -1, 1, 1 };
	qaws_curve* c = oct_bezier(xy, 3);
	qaws_curve* chain = NULL;
	double cusps[8], cp[2 * 4096], worst = 0, tol = OCT_F32 ? 1e-5 : 1e-9;
	unsigned int nc = 0, ncurves = 0, n, corners, i;
	char msg[300];
	qaws_internal_offset_curve(c, -1, tol, &chain, &ncurves, cusps, 8, &nc);
	sprintf(msg, "parabola offset 1 toward its focus: 2 cusps (%u) at the reference parameters (%.15g, %.15g vs %.15g, %.15g)",
		nc, nc > 0 ? cusps[0] : 0, nc > 1 ? cusps[1] : 0, ref_parabola_t[0], ref_parabola_t[1]);
	TEST_ASSERT(nc == 2 && fabs(cusps[0] - ref_parabola_t[0]) < (OCT_F32 ? 1e-4 : 1e-10)
		&& fabs(cusps[1] - ref_parabola_t[1]) < (OCT_F32 ? 1e-4 : 1e-10), msg);
	TEST_ASSERT(chain && ncurves == 1, "the parabola's offset is one chain");
	if (!chain) { qaws_curve_destroy(c); return; }
	n = oct_chain(chain, cp, 4096, &corners);
	{
		double d0 = oct_nearest_joint(cp, n, ref_parabola_p[0], ref_parabola_p[1]);
		double d1 = oct_nearest_joint(cp, n, ref_parabola_p[2], ref_parabola_p[3]);
		double dv = oct_nearest_joint(cp, n, ref_parabola_evolute_vertex[0], ref_parabola_evolute_vertex[1]);
		sprintf(msg, "the cusps and the evolute's vertex (0, 1/2) are joints of the chain (%.2e, %.2e, %.2e)", d0, d1, dv);
		TEST_ASSERT(d0 < 1e3 * tol && d1 < 1e3 * tol && dv < 1e3 * tol, msg);
	}
	sprintf(msg, "G1 everywhere but at the 2 cusps and the evolute's vertex (%u corners, %u cubics)", corners, (n - 1) / 3);
	TEST_ASSERT(corners == 3, msg);
	/* every point of the chain on the offset or the evolute */
	{
		qaws_range rg = qaws_curve_get_parameter_range(chain);
		for (i = 0; i <= 2000; i++)
		{
			qaws_eval_result_2d e;
			double d, d2;
			qaws_curve_evaluate_2d(chain, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * i / 2000), QAWS_EVAL_FLAG_POSITION, &e);
			d = oct_dist_to(oct_par_offset, e.position.x, e.position.y);
			d2 = oct_dist_to(oct_par_evolute, e.position.x, e.position.y);
			if (d2 < d) d = d2;
			if (d > worst) worst = d;
		}
		sprintf(msg, "the chain lies on the true offset or the evolute (worst %.2e, tolerance %.0e)", worst, tol);
		TEST_ASSERT(worst < (OCT_F32 ? 1e-4 : 1e-8), msg);
	}
	qaws_curve_destroy(chain);
	qaws_curve_destroy(c);
}

static void test_s_curve(void)
{
	static double const xy[8] = { 0, 0, 1, 2, 2, -2, 3, 0 };
	qaws_curve* c = oct_bezier(xy, 4);
	unsigned int side;
	char msg[300];
	for (side = 0; side < 2; side++)
	{
		double const* rt = side ? ref_s_left_t : ref_s_right_t;
		double const* rp = side ? ref_s_left_p : ref_s_right_p;
		double cusps[8], cp[2 * 4096], tol = OCT_F32 ? 1e-5 : 1e-9, d0 = 1e300, d1 = 1e300;
		unsigned int nc = 0, ncurves = 0, n, corners;
		qaws_curve* chain = NULL;
		qaws_internal_offset_curve(c, side ? -0.6 : 0.6, tol, &chain, &ncurves, cusps, 8, &nc);
		if (chain)
		{
			n = oct_chain(chain, cp, 4096, &corners);
			d0 = oct_nearest_joint(cp, n, rp[0], rp[1]);
			d1 = oct_nearest_joint(cp, n, rp[2], rp[3]);
		}
		sprintf(msg, "cubic S offset 0.6 to the %s: 2 cusps (%u) at the reference parameters (%.15g, %.15g vs %.15g, %.15g), joints at the reference points (%.2e, %.2e)",
			side ? "left" : "right", nc, nc > 0 ? cusps[0] : 0, nc > 1 ? cusps[1] : 0, rt[0], rt[1], d0, d1);
		TEST_ASSERT(nc == 2 && fabs(cusps[0] - rt[0]) < (OCT_F32 ? 1e-4 : 1e-10) && fabs(cusps[1] - rt[1]) < (OCT_F32 ? 1e-4 : 1e-10)
			&& d0 < 1e3 * tol && d1 < 1e3 * tol, msg);
		qaws_curve_destroy(chain);
	}
	qaws_curve_destroy(c);
}

/* a closed C2 star: the uniform periodic cubic B-spline of m points at
   alternating radii, written as m Bezier spans in a clamped B-spline */
static qaws_curve* oct_star(unsigned int m, double r0, double r1)
{
	double Q[64][2];
	qaws_scalar cp[2 * (3 * 32 + 1)], kn[3 * 32 + 5];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i, k, n = 0, kc = 0;
	for (i = 0; i < m; i++)
	{
		double a = 2 * OCT_PI * i / m, r = (i & 1) ? r1 : r0;
		Q[i][0] = r * cos(a); Q[i][1] = r * sin(a);
	}
	for (i = 0; i < m; i++)
	{
		double const* qm = Q[(i + m - 1) % m];
		double const* q0 = Q[i];
		double const* q1 = Q[(i + 1) % m];
		if (i == 0)
			for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((qm[k] + 4 * q0[k] + q1[k]) / 6);
		if (i == 0) n++;
		for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((2 * q0[k] + q1[k]) / 3);
		n++;
		for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((q0[k] + 2 * q1[k]) / 3);
		n++;
		{
			double const* q2 = Q[(i + 2) % m];
			for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((q0[k] + 4 * q1[k] + q2[k]) / 6);
		}
		n++;
	}
	/* the seam: the last point is the first */
	cp[2 * (n - 1)] = cp[0];
	cp[2 * (n - 1) + 1] = cp[1];
	for (i = 0; i < 4; i++) kn[kc++] = 0;
	for (i = 1; i < m; i++) { kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; }
	for (i = 0; i < 4; i++) kn[kc++] = (qaws_scalar)m;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = n;
	d.knots = kn; d.knot_count = kc;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

/* closed paths as rings of points, 256 per span (the band the oracle skips
   is far wider than the chords' sag) */
typedef struct oct_ring
{
	double* xy;
	unsigned int* start;    /* ring i: points [start[i], start[i + 1]) */
	unsigned int rings, count;
} oct_ring;

static void oct_rings(qaws_path_2d const* paths, unsigned int n, oct_ring* o)
{
	unsigned int i, k, j, total = 0, cap;
	memset(o, 0, sizeof(*o));
	for (i = 0; i < n; i++)
		for (k = 0; k < paths[i].curve_count; k++)
			total += 256 * qaws_curve_get_span_count(paths[i].curves[k]) + 1;
	cap = total + 1;
	o->xy = (double*)malloc(sizeof(double) * 2 * cap);
	o->start = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	for (i = 0; i < n; i++)
	{
		o->start[i] = o->count;
		for (k = 0; k < paths[i].curve_count; k++)
		{
			qaws_curve const* c = paths[i].curves[k];
			qaws_range rg = qaws_curve_get_parameter_range(c);
			unsigned int m = 256 * qaws_curve_get_span_count(c);
			for (j = 0; j < m; j++)
			{
				qaws_eval_result_2d e;
				qaws_curve_evaluate_2d(c, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * j / m), QAWS_EVAL_FLAG_POSITION, &e);
				o->xy[2 * o->count] = e.position.x;
				o->xy[2 * o->count + 1] = e.position.y;
				o->count++;
			}
		}
	}
	o->start[n] = o->count;
	o->rings = n;
}

/* even-odd crossings of a rightward ray */
static int oct_inside(oct_ring const* o, double x, double y)
{
	unsigned int i, a;
	int in = 0;
	for (i = 0; i < o->rings; i++)
	{
		unsigned int s = o->start[i], e = o->start[i + 1];
		for (a = s; a < e; a++)
		{
			unsigned int b = a + 1 < e ? a + 1 : s;
			double ax = o->xy[2 * a], ay = o->xy[2 * a + 1], bx = o->xy[2 * b], by = o->xy[2 * b + 1];
			if ((ay > y) != (by > y) && x < ax + (y - ay) * (bx - ax) / (by - ay))
				in = !in;
		}
	}
	return in;
}

/* grid points where the offset result disagrees with the distance oracle:
   in iff (closed: inside the region grown, or inside and far from its
   boundary shrunk; open: near the curve) */
static unsigned int oct_oracle(qaws_curve const* c, int closed, double delta, qaws_clip_result const* r,
	unsigned int grid, unsigned int* checked)
{
	qaws_path_2d src, *res;
	oct_ring rr, sr;
	qaws_vec2 lo, hi;
	unsigned int i, j, np = qaws_clip_result_get_path_count(r), bad = 0;
	double ext, band;
	qaws_scalar* pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * grid * grid);
	qaws_closest_point* cl = (qaws_closest_point*)malloc(sizeof(qaws_closest_point) * grid * grid);
	qaws_closest_desc cd;
	src.curves = &c; src.curve_count = 1; src.closed = closed;
	qaws_path_compute_bounds_2d(&src, &lo, &hi);
	ext = fabs(delta) + 0.2;
	lo.x -= (qaws_scalar)ext; lo.y -= (qaws_scalar)ext; hi.x += (qaws_scalar)ext; hi.y += (qaws_scalar)ext;
	band = OCT_F32 ? 2e-3 : 1e-4;
	for (j = 0; j < grid; j++)
		for (i = 0; i < grid; i++)
		{
			pts[2 * (j * grid + i)] = (qaws_scalar)(lo.x + (hi.x - lo.x) * (i + 0.5) / grid);
			pts[2 * (j * grid + i) + 1] = (qaws_scalar)(lo.y + (hi.y - lo.y) * (j + 0.5) / grid);
		}
	memset(&cd, 0, sizeof(cd));
	cd.curves = &c; cd.curve_count = 1; cd.points = pts; cd.point_count = grid * grid;
	qaws_curve_batch_find_closest(&cd, cl, NULL);
	res = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (np + 1));
	for (i = 0; i < np; i++) qaws_clip_result_get_path(r, i, &res[i]);
	memset(&sr, 0, sizeof(sr));
	oct_rings(res, np, &rr);
	if (closed) oct_rings(&src, 1, &sr);
	*checked = 0;
	for (i = 0; i < grid * grid; i++)
	{
		double px = pts[2 * i], py = pts[2 * i + 1], dist = (double)cl[i].distance;
		int want, got;
		if (fabs(dist - fabs(delta)) < band) continue;
		if (closed)
		{
			int in;
			if (dist < band) continue;
			in = oct_inside(&sr, px, py);
			want = delta > 0 ? (in || dist < delta) : (in && dist > -delta);
		}
		else
			want = dist < fabs(delta);
		got = oct_inside(&rr, px, py);
		(*checked)++;
		if (want != got) bad++;
	}
	free(rr.xy); free(rr.start);
	free(sr.xy); free(sr.start);
	free(res);
	free(pts);
	free(cl);
	return bad;
}

static void test_regions(void)
{
	qaws_curve* star = oct_star(10, 3, 0.6);
	qaws_curve const* cc = star;
	qaws_path_2d path;
	static double const deltas[2] = { 1.0, -0.45 };
	unsigned int k;
	char msg[300];
	path.curves = &cc; path.curve_count = 1; path.closed = 1;
	for (k = 0; k < 2; k++)
	{
		qaws_clip_result* r = NULL;
		qaws_curve* chain = NULL;
		double cusps[64];
		unsigned int nc = 0, nk = 0, bad = 0, checked = 0;
		qaws_internal_offset_curve(star, deltas[k], 1e-6, &chain, &nk, cusps, 64, &nc);
		qaws_curve_destroy(chain);
		qaws_offset_paths(&path, 1, (qaws_scalar)deltas[k], QAWS_JOIN_ROUND, QAWS_END_POLYGON, &r);
		if (r) bad = oct_oracle(star, 1, deltas[k], r, 180, &checked);
		sprintf(msg, "star %s by %g past its radii of curvature (%u cusps): %u of %u grid points against the distance oracle",
			deltas[k] > 0 ? "grown" : "shrunk", fabs(deltas[k]), nc, bad, checked);
		TEST_ASSERT(r && nc > 0 && bad == 0 && checked > 10000, msg);
		qaws_clip_result_destroy(r);
	}
	qaws_curve_destroy(star);
}

static void test_open(void)
{
	static double const xy[8] = { 0, 0, 1, 2, 2, -2, 3, 0 };
	qaws_curve* c = oct_bezier(xy, 4);
	qaws_curve const* cc = c;
	qaws_path_2d path;
	qaws_clip_result* r = NULL;
	unsigned int bad = 0, checked = 0;
	char msg[300];
	path.curves = &cc; path.curve_count = 1; path.closed = 0;
	qaws_offset_paths(&path, 1, (qaws_scalar)0.6, QAWS_JOIN_ROUND, QAWS_END_ROUND, &r);
	if (r) bad = oct_oracle(c, 0, 0.6, r, 180, &checked);
	sprintf(msg, "cubic S stroked 0.6 with round ends (cusps on both sides): %u of %u grid points against the distance oracle", bad, checked);
	TEST_ASSERT(r && bad == 0 && checked > 10000, msg);
	qaws_clip_result_destroy(r);
	qaws_curve_destroy(c);
}

static void test_simple(void)
{
	qaws_curve* c = NULL;
	qaws_curve const* cc;
	qaws_path_2d path;
	qaws_clip_result* r = NULL;
	unsigned int i, spans = 0;
	double a = 0;
	char msg[300];
	qaws_curve_create_ellipse_2d(0, 0, 5, 5, 0, &c);
	cc = c;
	path.curves = &cc; path.curve_count = 1; path.closed = 1;
	qaws_offset_paths(&path, 1, 1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &r);
	for (i = 0; r && i < qaws_clip_result_get_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_scalar pa = 0;
		unsigned int k;
		qaws_clip_result_get_path(r, i, &p);
		for (k = 0; k < p.curve_count; k++) spans += qaws_curve_get_span_count(p.curves[k]);
		qaws_path_compute_area_2d(&p, &pa);
		a += pa;
	}
	sprintf(msg, "NURBS circle r 5 grown by 1 at the default tolerance (1e-6 of the extent): %u cubics, area %.12g vs %.12g", spans, a, 36 * OCT_PI);
	TEST_ASSERT(r && spans <= 24 && fabs(a - 36 * OCT_PI) < (OCT_F32 ? 1e-3 : 1e-4), msg);
	qaws_clip_result_destroy(r);
	qaws_curve_destroy(c);
}

int test_90_offset_curves_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 90: offsets of curves (cusps, evolutes, G1 cubics)\n");
	test_parabola();
	test_s_curve();
	test_regions();
	test_open();
	test_simple();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
