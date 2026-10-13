/*
 * Test 102: regions on surfaces, Booleans and metric areas
 * (qaws_surface_region.h)
 *
 *   - a sphere: geodesic circles' areas 2 pi R^2 (1 - cos(r / R)) and
 *     perimeters 2 pi R sin(r / R); a geodesic triangle's area by Girard,
 *     R^2 (angle sum - pi), the angles from the great circles' planes; two
 *     circles' union, intersection, difference and xor by inclusion and
 *     exclusion
 *   - a flat parallelogram patch: a geodesic circle's area pi r^2, and a
 *     circle with a geodesic square against the 2D engine on the plane's own
 *     coordinates
 *   - a cylinder (developable): a geodesic circle round the seam has the
 *     area pi r^2 of the unrolled plane
 *   - a wavy bicubic patch: Gauss-Bonnet on a geodesic triangle (the
 *     integral of the Gaussian curvature equals the angle sum less pi) and
 *     on a geodesic circle (with its geodesic curvature)
 *   - containment, and the region as a trimmed surface
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SRT_PI 3.14159265358979323846
#define SRT_F32 (sizeof(qaws_scalar) == 4)

static qaws_vec2 srt_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static qaws_surface* srt_sphere(double R, qaws_curve** prof)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)R; g.angle_start = (qaws_scalar)(-SRT_PI / 2); g.angle_end = (qaws_scalar)(SRT_PI / 2);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, prof);
	memset(&rd, 0, sizeof(rd));
	rd.profile = *prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rd, &s);
	return s;
}

static double srt_area(qaws_surface const* s, qaws_clip_type type, qaws_surface_path const* a, unsigned int na,
	qaws_surface_path const* b, unsigned int nb, qaws_status* st)
{
	qaws_surface_region* r = NULL;
	qaws_scalar area = 0;
	*st = qaws_surface_clip(s, type, QAWS_FILL_NON_ZERO, a, na, b, nb, &r);
	if (*st == QAWS_STATUS_OK)
		*st = qaws_surface_region_get_area(r, &area);
	qaws_surface_region_destroy(r);
	return (double)area;
}

static void test_sphere(void)
{
	double R = 1.3, tol = SRT_F32 ? 2e-4 : 1e-9;
	qaws_curve* prof = NULL;
	qaws_surface* s = srt_sphere(R, &prof);
	qaws_surface_curve *ca = NULL, *cb = NULL;
	qaws_surface_curve const *pa[1], *pb[1];
	qaws_surface_path A, B;
	double ra = 0.6, rb = 0.5, aA, aB, aU, aI, aD, aX, worst_cf;
	qaws_status st, s1 = QAWS_STATUS_OK;
	char msg[400];
	qaws_surface_curve_create_geodesic_arc(s, srt_uv(0.3, 0.55), (qaws_scalar)ra, 0, (qaws_scalar)(2 * SRT_PI), 0, &ca);
	qaws_surface_curve_create_geodesic_arc(s, srt_uv(0.38, 0.5), (qaws_scalar)rb, 0, (qaws_scalar)(2 * SRT_PI), 0, &cb);
	pa[0] = ca; pb[0] = cb;
	A.curves = pa; A.curve_count = 1; A.closed = 1;
	B.curves = pb; B.curve_count = 1; B.closed = 1;
	aA = srt_area(s, QAWS_CLIP_NONE, &A, 1, NULL, 0, &st); if (st) s1 = st;
	aB = srt_area(s, QAWS_CLIP_NONE, &B, 1, NULL, 0, &st); if (st) s1 = st;
	aU = srt_area(s, QAWS_CLIP_UNION, &A, 1, &B, 1, &st); if (st) s1 = st;
	aI = srt_area(s, QAWS_CLIP_INTERSECTION, &A, 1, &B, 1, &st); if (st) s1 = st;
	aD = srt_area(s, QAWS_CLIP_DIFFERENCE, &A, 1, &B, 1, &st); if (st) s1 = st;
	aX = srt_area(s, QAWS_CLIP_XOR, &A, 1, &B, 1, &st); if (st) s1 = st;
	worst_cf = fabs(aA - 2 * SRT_PI * R * R * (1 - cos(ra / R)));
	if (fabs(aB - 2 * SRT_PI * R * R * (1 - cos(rb / R))) > worst_cf) worst_cf = fabs(aB - 2 * SRT_PI * R * R * (1 - cos(rb / R)));
	sprintf(msg, "sphere: two geodesic circles' areas against 2 pi R^2 (1 - cos(r / R)) within %.2g; union + intersection - both %.2g, difference + intersection - A %.2g, xor - union + intersection %.2g",
		worst_cf, aU + aI - aA - aB, aD + aI - aA, aX - aU + aI);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && worst_cf <= tol && fabs(aU + aI - aA - aB) <= tol && fabs(aD + aI - aA) <= tol && fabs(aX - aU + aI) <= tol && aI > 0.1, msg);
	{
		qaws_surface_region* r = NULL;
		qaws_scalar per = 0;
		int in1 = 0, in2 = 0;
		qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &A, 1, NULL, 0, &r);
		qaws_surface_region_get_perimeter(r, &per);
		qaws_surface_region_contains(r, srt_uv(0.3, 0.55), &in1);
		qaws_surface_region_contains(r, srt_uv(0.7, 0.55), &in2);
		sprintf(msg, "sphere: a geodesic circle's perimeter against 2 pi R sin(r / R) within %.2g; its centre inside (%d), the far side outside (%d)",
			fabs((double)per - 2 * SRT_PI * R * sin(ra / R)), in1, in2);
		TEST_ASSERT(fabs((double)per - 2 * SRT_PI * R * sin(ra / R)) <= tol && in1 == 1 && in2 == 0, msg);
		qaws_surface_region_destroy(r);
	}
	/* a geodesic triangle: Girard */
	{
		double tri[3][2] = { { 0.1, 0.6 }, { 0.3, 0.68 }, { 0.22, 0.35 } }, P[3][3], angle = 0, area;
		qaws_surface_curve* e[3] = { 0 };
		qaws_surface_curve const* pe[3];
		qaws_surface_path T;
		unsigned int i;
		for (i = 0; i < 3; i++)
		{
			double lon = 2 * SRT_PI * tri[i][0], lat = -SRT_PI / 2 + SRT_PI * tri[i][1];
			P[i][0] = cos(lat) * cos(lon); P[i][1] = cos(lat) * sin(lon); P[i][2] = sin(lat);
			qaws_surface_curve_create_geodesic(s, srt_uv(tri[i][0], tri[i][1]), srt_uv(tri[(i + 1) % 3][0], tri[(i + 1) % 3][1]), 0, &e[i]);
			pe[i] = e[i];
		}
		/* the angle at each vertex between the planes through it and its neighbours */
		for (i = 0; i < 3; i++)
		{
			double const* a = P[i];
			double const* b = P[(i + 1) % 3];
			double const* c = P[(i + 2) % 3];
			double n1[3], n2[3], d;
			n1[0] = a[1] * b[2] - a[2] * b[1]; n1[1] = a[2] * b[0] - a[0] * b[2]; n1[2] = a[0] * b[1] - a[1] * b[0];
			n2[0] = a[1] * c[2] - a[2] * c[1]; n2[1] = a[2] * c[0] - a[0] * c[2]; n2[2] = a[0] * c[1] - a[1] * c[0];
			d = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / sqrt((n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]) * (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]));
			angle += acos(d);
		}
		T.curves = pe; T.curve_count = 3; T.closed = 1;
		area = srt_area(s, QAWS_CLIP_NONE, &T, 1, NULL, 0, &st);
		sprintf(msg, "sphere: a geodesic triangle's area against Girard's R^2 (angles - pi) within %.2g", fabs(area - R * R * (angle - SRT_PI)));
		TEST_ASSERT(st == QAWS_STATUS_OK && fabs(area - R * R * (angle - SRT_PI)) <= tol, msg);
		for (i = 0; i < 3; i++) qaws_surface_curve_destroy(e[i]);
	}
	qaws_surface_curve_destroy(ca); qaws_surface_curve_destroy(cb);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

static void test_plane(void)
{
	/* the plane through o spanned by e1 (unit) and e2 (unit, orthogonal), scaled 2 and 3 along u and v, sheared */
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	qaws_surface_curve *circle = NULL, *sq[4] = { 0 };
	qaws_surface_curve const *pc[1], *ps[4];
	qaws_surface_path C, Q;
	double r = 0.45, tol = SRT_F32 ? 1e-4 : 1e-9, aC, aU;
	qaws_vec2 corner[4];
	unsigned int i;
	qaws_status st, s1 = QAWS_STATUS_OK;
	char msg[300];
	memset(&d, 0, sizeof(d));
	d.p00.x = 0; d.p00.y = 0; d.p00.z = 0;
	d.p10.x = 2; d.p10.y = 0; d.p10.z = 0;
	d.p01.x = (qaws_scalar)0.6; d.p01.y = 3; d.p01.z = 0;
	d.p11.x = (qaws_scalar)2.6; d.p11.y = 3; d.p11.z = 0;
	qaws_surface_create_bilinear(&d, &s);
	qaws_surface_curve_create_geodesic_arc(s, srt_uv(0.5, 0.5), (qaws_scalar)r, 0, (qaws_scalar)(2 * SRT_PI), 0, &circle);
	corner[0] = srt_uv(0.3, 0.3); corner[1] = srt_uv(0.62, 0.3); corner[2] = srt_uv(0.62, 0.58); corner[3] = srt_uv(0.3, 0.58);
	for (i = 0; i < 4; i++) { qaws_surface_curve_create_geodesic(s, corner[i], corner[(i + 1) % 4], 0, &sq[i]); ps[i] = sq[i]; }
	pc[0] = circle;
	C.curves = pc; C.curve_count = 1; C.closed = 1;
	Q.curves = ps; Q.curve_count = 4; Q.closed = 1;
	aC = srt_area(s, QAWS_CLIP_NONE, &C, 1, NULL, 0, &st); if (st) s1 = st;
	aU = srt_area(s, QAWS_CLIP_UNION, &C, 1, &Q, 1, &st); if (st) s1 = st;
	{
		/* the 2D engine on the plane's own coordinates (x, y): a circle and the square's parallelogram */
		qaws_arc_segment g;
		qaws_arc_desc ad;
		qaws_curve *cc = NULL, *pq = NULL;
		qaws_curve const *cv1[1], *cv2[1];
		qaws_path_2d p1, p2;
		qaws_clip_result* res = NULL;
		qaws_scalar xy[8];
		double a2 = 0;
		unsigned int k;
		memset(&g, 0, sizeof(g));
		g.center[0] = (qaws_scalar)(0.5 * 2 + 0.5 * 0.6); g.center[1] = (qaws_scalar)(0.5 * 3); g.radius = (qaws_scalar)r;
		g.angle_start = 0; g.angle_end = (qaws_scalar)(2 * SRT_PI);
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
		qaws_curve_create_arc(&ad, &cc);
		for (k = 0; k < 4; k++)
		{
			xy[2 * k] = (qaws_scalar)(corner[k].x * 2 + corner[k].y * 0.6);
			xy[2 * k + 1] = (qaws_scalar)(corner[k].y * 3);
		}
		qaws_curve_create_polyline_2d(xy, 4, 1, &pq);
		cv1[0] = cc; cv2[0] = pq;
		p1.curves = cv1; p1.curve_count = 1; p1.closed = 1;
		p2.curves = cv2; p2.curve_count = 1; p2.closed = 1;
		qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &p1, 1, &p2, 1, &res);
		for (k = 0; res && k < qaws_clip_result_get_path_count(res); k++)
		{
			qaws_path_2d q;
			qaws_scalar a;
			qaws_clip_result_get_path(res, k, &q);
			qaws_path_compute_area_2d(&q, &a);
			a2 += (double)a;
		}
		sprintf(msg, "flat sheared patch: a geodesic circle's area against pi r^2 within %.2g; its union with a geodesic square against the 2D engine within %.2g",
			fabs(aC - SRT_PI * r * r), fabs(aU - a2));
		TEST_ASSERT(s1 == QAWS_STATUS_OK && fabs(aC - SRT_PI * r * r) <= tol && fabs(aU - a2) <= tol && aU > aC, msg);
		qaws_clip_result_destroy(res);
		qaws_curve_destroy(cc); qaws_curve_destroy(pq);
	}
	qaws_surface_curve_destroy(circle);
	for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
	qaws_surface_destroy(s);
}

static void test_cylinder(void)
{
	double R = 0.8, H = 2.5, r = 0.5, area = 0;
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* circle = NULL;
	qaws_surface_extrusion_desc ed;
	qaws_surface* s = NULL;
	qaws_surface_curve* c = NULL;
	qaws_surface_curve const* pc[1];
	qaws_surface_path P;
	qaws_range ru;
	qaws_status st;
	char msg[300];
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)R; g.angle_start = 0; g.angle_end = (qaws_scalar)(2 * SRT_PI);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, &circle);
	memset(&ed, 0, sizeof(ed));
	ed.profile = circle; ed.direction.z = 1; ed.length = (qaws_scalar)H;
	qaws_surface_create_extrusion(&ed, &s);
	ru = qaws_surface_get_u_range(s);
	/* centred near the middle of the domain, clear of the seam (step 4 cuts at seams) */
	st = qaws_surface_curve_create_geodesic_arc(s, srt_uv(0.5 * ru.max_value, 0.5), (qaws_scalar)r, 0, (qaws_scalar)(2 * SRT_PI), 0, &c);
	pc[0] = c;
	P.curves = pc; P.curve_count = 1; P.closed = 1;
	if (st == QAWS_STATUS_OK)
		area = srt_area(s, QAWS_CLIP_NONE, &P, 1, NULL, 0, &st);
	sprintf(msg, "cylinder: a geodesic circle's area against the unrolled pi r^2 within %.2g", fabs(area - SRT_PI * r * r));
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs(area - SRT_PI * r * r) <= (SRT_F32 ? 1e-4 : 1e-9), msg);
	qaws_surface_curve_destroy(c);
	qaws_surface_destroy(s);
	qaws_curve_destroy(circle);
}

/* the Gaussian curvature at (u, v) */
static void srt_curvature(void* user, qaws_scalar u, qaws_scalar v, qaws_vec3 p, qaws_scalar* out)
{
	qaws_surface const* s = (qaws_surface const*)user;
	qaws_surface_eval_result r;
	double n[3], nl, E, F, G, L, M, N;
	(void)p;
	qaws_surface_evaluate(s, u, v, QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV, &r);
	n[0] = (double)r.du.y * r.dv.z - (double)r.du.z * r.dv.y;
	n[1] = (double)r.du.z * r.dv.x - (double)r.du.x * r.dv.z;
	n[2] = (double)r.du.x * r.dv.y - (double)r.du.y * r.dv.x;
	nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	E = (double)r.du.x * r.du.x + (double)r.du.y * r.du.y + (double)r.du.z * r.du.z;
	F = (double)r.du.x * r.dv.x + (double)r.du.y * r.dv.y + (double)r.du.z * r.dv.z;
	G = (double)r.dv.x * r.dv.x + (double)r.dv.y * r.dv.y + (double)r.dv.z * r.dv.z;
	L = (r.duu.x * n[0] + r.duu.y * n[1] + r.duu.z * n[2]) / nl;
	M = (r.duv.x * n[0] + r.duv.y * n[1] + r.duv.z * n[2]) / nl;
	N = (r.dvv.x * n[0] + r.dvv.y * n[1] + r.dvv.z * n[2]) / nl;
	out[0] = (qaws_scalar)((L * N - M * M) / (E * G - F * F));
}

/* the angle between (u, v) directions a and b at p, by the metric */
static double srt_angle(qaws_surface const* s, qaws_vec2 p, double const* a, double const* b)
{
	qaws_scalar g[3];
	double ab, aa, bb;
	qaws_surface_get_metric(s, p.x, p.y, g, NULL);
	ab = g[0] * a[0] * b[0] + g[1] * (a[0] * b[1] + a[1] * b[0]) + g[2] * a[1] * b[1];
	aa = g[0] * a[0] * a[0] + 2 * g[1] * a[0] * a[1] + g[2] * a[1] * a[1];
	bb = g[0] * b[0] * b[0] + 2 * g[1] * b[0] * b[1] + g[2] * b[1] * b[1];
	return acos(ab / sqrt(aa * bb));
}

static void test_gauss_bonnet(void)
{
	static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
	qaws_vec3 cp[16];
	qaws_surface_bezier_desc bd;
	qaws_surface* s = NULL;
	qaws_vec2 tri[3];
	qaws_surface_curve* e[3] = { 0 };
	qaws_surface_curve const* pe[3];
	qaws_surface_path T;
	qaws_surface_region* r = NULL;
	qaws_scalar Kint = 0;
	double angles = 0, tol = SRT_F32 ? 2e-4 : 1e-8;
	unsigned int i;
	qaws_status st;
	char msg[300];
	for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
	memset(&bd, 0, sizeof(bd));
	bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
	qaws_surface_create_bezier(&bd, &s);
	tri[0] = srt_uv(0.15, 0.2); tri[1] = srt_uv(0.85, 0.3); tri[2] = srt_uv(0.4, 0.85);
	for (i = 0; i < 3; i++) { qaws_surface_curve_create_geodesic(s, tri[i], tri[(i + 1) % 3], 0, &e[i]); pe[i] = e[i]; }
	/* the interior angle at each vertex: between the incoming edge reversed and the outgoing edge */
	for (i = 0; i < 3; i++)
	{
		qaws_curve const* in = qaws_surface_curve_get_uv(e[(i + 2) % 3]);
		qaws_curve const* out = qaws_surface_curve_get_uv(e[i]);
		qaws_range ri = qaws_curve_get_parameter_range(in), ro = qaws_curve_get_parameter_range(out);
		qaws_eval_result_2d a, b;
		double da[2], db[2];
		qaws_curve_evaluate_2d(in, ri.max_value, QAWS_EVAL_FLAG_D1, &a);
		qaws_curve_evaluate_2d(out, ro.min_value, QAWS_EVAL_FLAG_D1, &b);
		da[0] = -a.d1.x; da[1] = -a.d1.y; db[0] = b.d1.x; db[1] = b.d1.y;
		angles += srt_angle(s, tri[i], da, db);
	}
	T.curves = pe; T.curve_count = 3; T.closed = 1;
	st = qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &T, 1, NULL, 0, &r);
	if (st == QAWS_STATUS_OK)
		st = qaws_surface_region_integrate(r, srt_curvature, s, 1, &Kint);
	sprintf(msg, "wavy patch, Gauss-Bonnet on a geodesic triangle: the integral of K (%.6f) against the angle sum less pi (%.6f): %.2g",
		(double)Kint, angles - SRT_PI, fabs((double)Kint - (angles - SRT_PI)));
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs((double)Kint - (angles - SRT_PI)) <= tol && fabs((double)Kint) > 1e-3, msg);
	{
		qaws_surface* tr = NULL;
		int in = 0;
		qaws_status s2 = qaws_surface_region_create_trimmed(r, &tr), s3 = QAWS_STATUS_INVALID_ARGUMENT;
		if (s2 == QAWS_STATUS_OK)
			s3 = qaws_surface_trim_contains(tr, (qaws_scalar)0.45, (qaws_scalar)0.45, &in);
		sprintf(msg, "the triangle as a trimmed surface (status %d), its inside point kept (%d)", (int)s2, in);
		TEST_ASSERT(s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_OK && in == 1, msg);
		qaws_surface_destroy(tr);
	}
	qaws_surface_region_destroy(r);
	for (i = 0; i < 3; i++) qaws_surface_curve_destroy(e[i]);
	qaws_surface_destroy(s);
}

int test_102_surface_region_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 102: regions on surfaces, Booleans and metric areas\n");
	test_sphere();
	test_plane();
	test_cylinder();
	test_gauss_bonnet();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
