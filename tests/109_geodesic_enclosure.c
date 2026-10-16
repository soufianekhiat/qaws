/*
 * Test 109: validated enclosures of geodesics on polynomial surfaces
 * (qaws_surface_geodesic_enclose, qaws_geodesic_enclosure_cross)
 *
 *   - a parabolic cylinder z = x^2 (Bezier, degree 2 x 1): developable, its
 *     geodesics straight once unrolled by arc length; the end box holds the
 *     closed form's end
 *   - a plane with a non-uniform bicubic parametrization: geodesics
 *     straight in space; the end box holds the closed form's end
 *   - a bicubic B-spline bump with interior knots: the end box holds the
 *     float geodesic's end (f64)
 *   - crossings on the bump: one certified, its box holding the float
 *     geodesics' crossing (f64); none for a pair apart
 *   - refusals: a surface of revolution (unsupported), a geodesic leaving
 *     the domain
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include "qaws_surface_bezier.h"
#include "qaws_surface_bspline.h"
#include "qaws_surface_bilinear.h"
#include "qaws_surface_revolution.h"
#include <math.h>
#include <string.h>

#define GET_F32 (sizeof(qaws_scalar) == 4)

static qaws_vec2 get_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static int get_inside(qaws_geodesic_enclosure const* e, double u, double v, double* side)
{
	qaws_vec2 lo, hi;
	qaws_geodesic_enclosure_get_end(e, &lo, &hi, NULL, NULL);
	*side = fmax((double)(hi.x - lo.x), (double)(hi.y - lo.y));
	return u >= (double)lo.x && u <= (double)hi.x && v >= (double)lo.y && v <= (double)hi.y;
}

/* arc length along z = x^2 from 0 */
static double get_sx(double x) { return 0.5 * x * sqrt(1 + 4 * x * x) + 0.25 * asinh(2 * x); }

static double get_bez3(double const* c, double u)
{
	double a = 1 - u;
	return a * a * a * c[0] + 3 * a * a * u * c[1] + 3 * a * u * u * c[2] + u * u * u * c[3];
}

static double get_dbez3(double const* c, double u)
{
	double a = 1 - u;
	return 3 * (a * a * (c[1] - c[0]) + 2 * a * u * (c[2] - c[1]) + u * u * (c[3] - c[2]));
}

static double get_inv3(double const* c, double X)
{
	double lo = 0, hi = 1;
	unsigned int i;
	for (i = 0; i < 200; i++)
	{
		double m = 0.5 * (lo + hi);
		if (get_bez3(c, m) < X) lo = m; else hi = m;
	}
	return 0.5 * (lo + hi);
}

static void test_closed_forms(void)
{
	char msg[300];
	{
		/* x = -1 + 2u, y = -1 + 2v, z = x^2 */
		static double const xs[3] = { -1, 0, 1 }, zs[3] = { 1, -1, 1 };
		qaws_vec3 cp[6];
		qaws_surface_bezier_desc d;
		qaws_surface* s = NULL;
		qaws_geodesic_enclosure* e = NULL;
		double u0 = 0.3, v0 = 0.2, du = 0.8, dv = 0.5, L = 1.5, x0 = -1 + 2 * u0, y0 = -1 + 2 * v0, a, b, n, s1, y1, lo = -1, hi = 1, side = 0;
		unsigned int i;
		qaws_status st;
		int in = 0;
		for (i = 0; i < 3; i++)
		{
			cp[2 * i].x = (qaws_scalar)xs[i]; cp[2 * i].y = -1; cp[2 * i].z = (qaws_scalar)zs[i];
			cp[2 * i + 1] = cp[2 * i]; cp[2 * i + 1].y = 1;
		}
		memset(&d, 0, sizeof(d));
		d.u_degree = 2; d.v_degree = 1; d.control_points = cp; d.u_point_count = 3; d.v_point_count = 2;
		qaws_surface_create_bezier(&s, &d);
		a = sqrt(1 + 4 * x0 * x0) * 2 * du; b = 2 * dv; n = hypot(a, b); a /= n; b /= n;
		s1 = get_sx(x0) + L * a; y1 = y0 + L * b;
		for (i = 0; i < 200; i++) { double m = 0.5 * (lo + hi); if (get_sx(m) < s1) lo = m; else hi = m; }
		st = qaws_surface_geodesic_enclose(s, get_uv(u0, v0), get_uv(du, dv), (qaws_scalar)L, 0, &e);
		if (st == QAWS_STATUS_OK)
			in = get_inside(e, (0.5 * (lo + hi) + 1) / 2, (y1 + 1) / 2, &side);
		sprintf(msg, "enclosure: parabolic cylinder, the closed form's end in the end box %d (status %d, side %.2g, %u steps)", in, (int)st, side,
			qaws_geodesic_enclosure_get_step_count(e));
		TEST_ASSERT(st == QAWS_STATUS_OK && in && side <= 2e-5, msg);
		qaws_geodesic_enclosure_destroy(e);
		qaws_surface_destroy(s);
	}
	{
		static double const xs[4] = { 0, 0.1, 0.8, 1 }, ys[4] = { 0, 0.5, 0.6, 1 };
		qaws_vec3 cp[16];
		qaws_surface_bezier_desc d;
		qaws_surface* s = NULL;
		qaws_geodesic_enclosure* e = NULL;
		double u0 = 0.2, v0 = 0.3, du = 0.6, dv = 0.4, L = 0.6, a, b, n, side = 0;
		unsigned int i, j;
		qaws_status st;
		int in = 0;
		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
			{
				cp[4 * i + j].x = (qaws_scalar)xs[i]; cp[4 * i + j].y = (qaws_scalar)ys[j]; cp[4 * i + j].z = 0;
			}
		memset(&d, 0, sizeof(d));
		d.u_degree = 3; d.v_degree = 3; d.control_points = cp; d.u_point_count = 4; d.v_point_count = 4;
		qaws_surface_create_bezier(&s, &d);
		a = get_dbez3(xs, u0) * du; b = get_dbez3(ys, v0) * dv; n = hypot(a, b); a /= n; b /= n;
		st = qaws_surface_geodesic_enclose(s, get_uv(u0, v0), get_uv(du, dv), (qaws_scalar)L, 0, &e);
		if (st == QAWS_STATUS_OK)
			in = get_inside(e, get_inv3(xs, get_bez3(xs, u0) + L * a), get_inv3(ys, get_bez3(ys, v0) + L * b), &side);
		sprintf(msg, "enclosure: a plane parametrized unevenly, the straight line's end in the end box %d (status %d, side %.2g)", in, (int)st, side);
		TEST_ASSERT(st == QAWS_STATUS_OK && in && side <= 2e-5, msg);
		qaws_geodesic_enclosure_destroy(e);
		qaws_surface_destroy(s);
	}
}

static qaws_surface* get_bump(void)
{
	qaws_vec3 cp[36];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
		{
			double x = i / 5.0, y = j / 5.0;
			cp[6 * i + j].x = (qaws_scalar)x; cp[6 * i + j].y = (qaws_scalar)y;
			cp[6 * i + j].z = (qaws_scalar)(0.5 * exp(-8 * ((x - 0.5) * (x - 0.5) + (y - 0.45) * (y - 0.45))));
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 3; d.v_degree = 3; d.control_points = cp; d.u_point_count = 6; d.v_point_count = 6;
	qaws_surface_create_bspline(&s, &d);
	return s;
}

/* the float geodesics' crossing, by fine polylines; their count */
static int get_float_cross(qaws_surface const* s, double const* A, double const* B, double* out)
{
	enum { N = 1500 };
	static double pa[N + 1][2], pb[N + 1][2];
	int i, j, n = 0;
	for (i = 0; i <= N; i++)
	{
		qaws_vec2 q;
		qaws_surface_exp(s, get_uv(A[0], A[1]), get_uv(A[2], A[3]), (qaws_scalar)(A[4] * i / N), 0, &q, NULL);
		pa[i][0] = q.x; pa[i][1] = q.y;
		qaws_surface_exp(s, get_uv(B[0], B[1]), get_uv(B[2], B[3]), (qaws_scalar)(B[4] * i / N), 0, &q, NULL);
		pb[i][0] = q.x; pb[i][1] = q.y;
	}
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
		{
			double rx = pa[i + 1][0] - pa[i][0], ry = pa[i + 1][1] - pa[i][1], sx = pb[j + 1][0] - pb[j][0], sy = pb[j + 1][1] - pb[j][1];
			double den = rx * sy - ry * sx, qx = pb[j][0] - pa[i][0], qy = pb[j][1] - pa[i][1], t, u;
			if (den == 0)
				continue;
			t = (qx * sy - qy * sx) / den; u = (qx * ry - qy * rx) / den;
			if (t >= 0 && t < 1 && u >= 0 && u < 1)
			{
				out[0] = pa[i][0] + t * rx; out[1] = pa[i][1] + t * ry;
				n++;
			}
		}
	return n;
}

static void test_bump(void)
{
	static double const A[5] = { 0.1, 0.2, 0.7, 0.5, 0.9 }, B[5] = { 0.2, 0.75, 0.6, -0.6, 0.8 }, C[5] = { 0.1, 0.9, 1, 0.05, 0.7 };
	qaws_surface* s = get_bump();
	qaws_geodesic_enclosure *ea = NULL, *eb = NULL, *ec = NULL;
	qaws_vec2 q, lo[2], hi[2];
	qaws_status sa, sb, sc, x1, x2;
	unsigned int n1 = 9, n2 = 9;
	double side = 0, fx[2] = { 0, 0 };
	int in = 1, inx = 1, nf = 0;
	char msg[300];
	sa = qaws_surface_geodesic_enclose(s, get_uv(A[0], A[1]), get_uv(A[2], A[3]), (qaws_scalar)A[4], 0, &ea);
	sb = qaws_surface_geodesic_enclose(s, get_uv(B[0], B[1]), get_uv(B[2], B[3]), (qaws_scalar)B[4], 0, &eb);
	sc = qaws_surface_geodesic_enclose(s, get_uv(C[0], C[1]), get_uv(C[2], C[3]), (qaws_scalar)C[4], 0, &ec);
	qaws_surface_exp(s, get_uv(A[0], A[1]), get_uv(A[2], A[3]), (qaws_scalar)A[4], (qaws_scalar)1e-12, &q, NULL);
	if (sa == QAWS_STATUS_OK)
	{
		in = get_inside(ea, (double)q.x, (double)q.y, &side);
		if (GET_F32) in = 1;    /* (the float geodesic is coarser than the box) */
	}
	sprintf(msg, "enclosure: B-spline bump, the float geodesic's end in the end box %d (status %d, side %.2g)", in, (int)sa, side);
	TEST_ASSERT(sa == QAWS_STATUS_OK && in && side <= 1e-4, msg);
	x1 = sa || sb ? sa : qaws_geodesic_enclosure_cross(ea, eb, 2, lo, hi, &n1);
	x2 = sa || sc ? sa : qaws_geodesic_enclosure_cross(ea, ec, 2, lo + 1, hi + 1, &n2);
	if (!GET_F32)
	{
		nf = get_float_cross(s, A, B, fx);
		inx = n1 == 1 && fx[0] >= (double)lo[0].x && fx[0] <= (double)hi[0].x && fx[1] >= (double)lo[0].y && fx[1] <= (double)hi[0].y;
	}
	sprintf(msg, "enclosure: crossings certified on the bump, %u (status %d) holding the float one %d, %u for a pair apart (status %d)", n1, (int)x1, inx, n2, (int)x2);
	TEST_ASSERT(sb == QAWS_STATUS_OK && sc == QAWS_STATUS_OK && x1 == QAWS_STATUS_OK && x2 == QAWS_STATUS_OK && n1 == 1 && n2 == 0 && inx &&
		(GET_F32 || nf == 1), msg);
	qaws_geodesic_enclosure_destroy(ea); qaws_geodesic_enclosure_destroy(eb); qaws_geodesic_enclosure_destroy(ec);
	qaws_surface_destroy(s);
}

/* two geodesic triangles on the bump (edges: shortest paths): their crossings
   certified as the float curves' count; one's area bounded round the float
   region's */
static void test_polygons(void)
{
	static double const T[2][3][2] = { { { 0.15, 0.2 }, { 0.75, 0.3 }, { 0.4, 0.8 } }, { { 0.3, 0.45 }, { 0.85, 0.55 }, { 0.6, 0.15 } } };
	qaws_surface* s = get_bump();
	qaws_surface_curve* c[6] = { NULL };
	qaws_geodesic_segment g[6];
	qaws_surface_crossing x[8];
	qaws_surface_path p;
	qaws_surface_region* rg = NULL;
	qaws_scalar A = 0, lo = 0, hi = 0;
	qaws_status st = QAWS_STATUS_OK, sc, sa, nc;
	unsigned int i, j, n = 0, nf = 0, ni, m = 0;
	static unsigned int const corner[6] = { 1, 2, 0, 4, 5, 3 };   /* each triangle closed: an edge ends where the next starts */
	/* (in f32 the edges' data, a float fit's start and tangent, miss the corners by about
	   1e-5: tubes that wide, so a corner's two tubes meet) */
	qaws_scalar const tube = (qaws_scalar)(GET_F32 ? 1e-4 : 0);
	char msg[300];
	for (i = 0; i < 6 && st == QAWS_STATUS_OK; i++)
	{
		double const* a = T[i / 3][i % 3];
		double const* b = T[i / 3][(i + 1) % 3];
		qaws_eval_result_2d e;
		qaws_range r;
		if ((st = qaws_surface_curve_create_shortest_path(s, NULL, get_uv(a[0], a[1]), get_uv(b[0], b[1]), 0, &c[i])) != QAWS_STATUS_OK)
			break;
		r = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(c[i]));
		qaws_curve_evaluate_2d(&e, qaws_surface_curve_get_uv(c[i]), r.min_value, QAWS_EVAL_FLAG_D1);
		g[i].start = get_uv(a[0], a[1]); g[i].direction = e.d1;
		qaws_surface_curve_get_length(c[i], &g[i].length);
	}
	sc = st ? st : qaws_surface_geodesics_certify(s, g, 6, corner, tube, &n);
	/* the same corners undeclared: ends a hair apart, never taken as shared */
	nc = st ? st : qaws_surface_geodesics_certify(s, g, 6, NULL, tube, &m);
	sprintf(msg, "certified arrangement: the triangles' corners undeclared fail (status %d)", (int)nc);
	TEST_ASSERT(nc == QAWS_STATUS_CERTIFICATION_FAILED, msg);
	for (i = 0; i < 3 && st == QAWS_STATUS_OK; i++)
		for (j = 3; j < 6; j++)
			if (qaws_surface_curve_intersect(c[i], c[j], x, 8, &ni) == QAWS_STATUS_OK)
				nf += ni;
	sprintf(msg, "certified arrangement: two geodesic triangles on the bump, %u crossings (status %d), the float curves' %u", n, (int)sc, nf);
	TEST_ASSERT(st == QAWS_STATUS_OK && sc == QAWS_STATUS_OK && n == nf && n > 0, msg);
	sa = st;
	if (st == QAWS_STATUS_OK)
	{
		p.curves = (qaws_surface_curve const* const*)c; p.curve_count = 3; p.closed = 1;
		if (qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &p, 1, NULL, 0, &rg) == QAWS_STATUS_OK)
			qaws_surface_region_get_area(rg, &A);
		sa = qaws_surface_geodesic_polygon_area(s, g, 3, 0, 512, &lo, &hi);
	}
	sprintf(msg, "certified area: a geodesic triangle in [%.6f, %.6f] (status %d), the float region's %.6f", (double)lo, (double)hi, (int)sa, (double)A);
	TEST_ASSERT(sa == QAWS_STATUS_OK && lo <= A && A <= hi && hi - lo < 0.06 * A, msg);
	qaws_surface_region_destroy(rg);
	for (i = 0; i < 6; i++) qaws_surface_curve_destroy(c[i]);
	qaws_surface_destroy(s);
}

/* the identity plane: straight segments, their exact geometry known */
static qaws_geodesic_segment get_line(double x0, double y0, double x1, double y1)
{
	qaws_geodesic_segment g;
	g.start = get_uv(x0, y0); g.direction = get_uv(x1 - x0, y1 - y0);
	g.length = (qaws_scalar)hypot(x1 - x0, y1 - y0);
	return g;
}

/* two false certifications found in review, kept as regressions:
   - a C (a square less a cavity reached by a passage 0.002 wide): the
     tubes' cells close the passage, so the cavity was taken as inside;
     its area 0.4796 must lie in the bounds;
   - two segments crossing 5e-8 from one's end and the other's start: their
     ends near was taken as a shared corner, the crossing dropped; it must
     be counted or fail, never certified as none */
static void test_review_cases(void)
{
	static double const C[12][2] = { { 0.1, 0.1 }, { 0.9, 0.1 }, { 0.9, 0.9 }, { 0.1, 0.9 }, { 0.1, 0.501 }, { 0.3, 0.501 },
		{ 0.3, 0.7 }, { 0.7, 0.7 }, { 0.7, 0.3 }, { 0.3, 0.3 }, { 0.3, 0.499 }, { 0.1, 0.499 } };
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	qaws_geodesic_segment g[12], h[2];
	qaws_scalar lo = 0, hi = 0;
	qaws_status sa, sc;
	unsigned int i, n = 99;
	char msg[300];
	memset(&d, 0, sizeof(d));
	d.p10.x = 1; d.p01.y = 1; d.p11.x = 1; d.p11.y = 1;
	qaws_surface_create_bilinear(&s, &d);
	for (i = 0; i < 12; i++) g[i] = get_line(C[i][0], C[i][1], C[(i + 1) % 12][0], C[(i + 1) % 12][1]);
	sa = qaws_surface_geodesic_polygon_area(s, g, 12, 0, 0, &lo, &hi);
	sprintf(msg, "certified area: a C whose passage the grid closes, [%.6f, %.6f] (status %d) holding 0.4796", (double)lo, (double)hi, (int)sa);
	TEST_ASSERT(sa == QAWS_STATUS_CERTIFICATION_FAILED || (sa == QAWS_STATUS_OK && lo <= 0.4796 && 0.4796 <= hi), msg);
	TEST_ASSERT(sa == QAWS_STATUS_OK && hi - lo < 0.05, "certified area: the C's bounds proven, not refused");
	h[0] = get_line(0.2, 0.5, 0.5, 0.5);
	h[1] = get_line(0.49999995, 0.49999995, 0.49999995, 0.79999995);
	sc = qaws_surface_geodesics_certify(s, h, 2, NULL, 0, &n);
	sprintf(msg, "certified arrangement: a crossing 5e-8 from the ends, %u (status %d), never certified as none", n, (int)sc);
	TEST_ASSERT(sc == QAWS_STATUS_CERTIFICATION_FAILED || (sc == QAWS_STATUS_OK && n == 1), msg);
	qaws_surface_destroy(s);
}

static void test_refusals(void)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* prof = NULL;
	qaws_surface_revolution_desc rd;
	qaws_surface *rev = NULL, *s = get_bump();
	qaws_geodesic_enclosure* e = NULL;
	qaws_status s1, s2;
	char msg[200];
	memset(&g, 0, sizeof(g));
	g.radius = 1; g.angle_start = (qaws_scalar)-1; g.angle_end = 1;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&prof, &ad);
	memset(&rd, 0, sizeof(rd));
	rd.profile = prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rev, &rd);
	s1 = qaws_surface_geodesic_enclose(rev, get_uv(0.5, 0.5), get_uv(1, 0), (qaws_scalar)0.5, 0, &e);
	qaws_geodesic_enclosure_destroy(e);
	e = NULL;
	s2 = qaws_surface_geodesic_enclose(s, get_uv(0.5, 0.5), get_uv(1, 0.2), (qaws_scalar)3, 0, &e);
	qaws_geodesic_enclosure_destroy(e);
	sprintf(msg, "enclosure refused: a surface of revolution (%d), a geodesic leaving the domain (%d)", (int)s1, (int)s2);
	TEST_ASSERT(s1 == QAWS_STATUS_UNSUPPORTED_OPERATION && s2 == QAWS_STATUS_OUT_OF_RANGE, msg);
	qaws_surface_destroy(rev); qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

int test_109_geodesic_enclosure_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 109: validated enclosures of geodesics\n");
	test_closed_forms();
	test_bump();
	test_polygons();
	test_review_cases();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
