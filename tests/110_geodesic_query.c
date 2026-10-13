/*
 * Test 110: geodesic queries on a surface (qaws_surface_distance,
 * qaws_surface_curve_closest_point, qaws_surface_curve_intersect,
 * qaws_surface_curve_closest_points, qaws_surface_region_distance)
 *
 *   - a sphere: the distance between two points against the great circle;
 *     a point's distance to a geodesic against the angle to its plane; two
 *     geodesics' crossing against the great circles' meeting (closest
 *     points then 0); a point's signed distance to a geodesic disk outside
 *     and inside against the distance to its centre less the radius
 *   - a cylinder: two geodesic segments apart, their closest distance
 *     against the unrolled segments'
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include "qaws_surface_revolution.h"
#include "qaws_bezier.h"
#include <math.h>
#include <string.h>

#define GQ_PI 3.14159265358979323846
#define GQ_F32 (sizeof(qaws_scalar) == 4)

static qaws_vec2 gq_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static void gq_at(qaws_surface const* s, double u, double v, double* x)
{
	qaws_surface_eval_result r;
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION, &r);
	x[0] = r.position.x; x[1] = r.position.y; x[2] = r.position.z;
}

static double gq_dot(double const* a, double const* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

static void gq_cross(double const* a, double const* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1]; c[1] = a[2] * b[0] - a[0] * b[2]; c[2] = a[0] * b[1] - a[1] * b[0];
}

/* the point's distance to the segment ab in the plane */
static double gq_seg(double const* a, double const* b, double const* p)
{
	double dx = b[0] - a[0], dy = b[1] - a[1], t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy);
	if (t < 0) t = 0;
	if (t > 1) t = 1;
	return hypot(a[0] + t * dx - p[0], a[1] + t * dy - p[1]);
}

static void test_sphere(void)
{
	double R = 1.3, tol = GQ_F32 ? 1e-5 : 1e-9;
	double A[2] = { 0.1, 0.4 }, B[2] = { 0.3, 0.65 }, X[2] = { 0.25, 0.35 }, C[2] = { 0.12, 0.62 }, D[2] = { 0.3, 0.38 };
	double a3[3], b3[3], x3[3], c3[3], d3[3], n[3], m[3], k[3], kn, p3[3], ref, off;
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_curve* prof = NULL;
	qaws_surface* s = NULL;
	qaws_surface_curve *ab = NULL, *cd = NULL;
	qaws_surface_crossing xs[4];
	qaws_scalar d = 0, d2 = 1;
	unsigned int nx = 0;
	qaws_status st;
	char msg[300];
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)R; g.angle_start = (qaws_scalar)(-GQ_PI / 2); g.angle_end = (qaws_scalar)(GQ_PI / 2);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, &prof);
	memset(&rd, 0, sizeof(rd));
	rd.profile = prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rd, &s);
	gq_at(s, A[0], A[1], a3); gq_at(s, B[0], B[1], b3); gq_at(s, X[0], X[1], x3); gq_at(s, C[0], C[1], c3); gq_at(s, D[0], D[1], d3);
	st = qaws_surface_distance(s, gq_uv(A[0], A[1]), gq_uv(B[0], B[1]), 0, &d);
	ref = R * acos(gq_dot(a3, b3) / (R * R));
	sprintf(msg, "geodesic distance on a sphere against the great circle: %.2g (status %d)", fabs((double)d - ref), (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs((double)d - ref) <= tol, msg);
	qaws_surface_curve_create_geodesic(s, gq_uv(A[0], A[1]), gq_uv(B[0], B[1]), 0, &ab);
	qaws_surface_curve_create_geodesic(s, gq_uv(C[0], C[1]), gq_uv(D[0], D[1]), 0, &cd);
	gq_cross(a3, b3, n);
	st = qaws_surface_curve_closest_point(ab, gq_uv(X[0], X[1]), 0, NULL, NULL, &d);
	ref = R * asin(fabs(gq_dot(n, x3)) / (sqrt(gq_dot(n, n)) * R));
	sprintf(msg, "a point's distance to a geodesic on a sphere against the angle to its plane: %.2g (status %d)", fabs((double)d - ref), (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs((double)d - ref) <= tol, msg);
	st = qaws_surface_curve_intersect(ab, cd, xs, 4, &nx);
	gq_cross(c3, d3, m); gq_cross(n, m, k); kn = sqrt(gq_dot(k, k));
	off = 1;
	if (nx)
	{
		gq_at(s, xs[0].point.x, xs[0].point.y, p3);
		off = fmin(sqrt(pow(p3[0] - R * k[0] / kn, 2) + pow(p3[1] - R * k[1] / kn, 2) + pow(p3[2] - R * k[2] / kn, 2)),
			sqrt(pow(p3[0] + R * k[0] / kn, 2) + pow(p3[1] + R * k[1] / kn, 2) + pow(p3[2] + R * k[2] / kn, 2)));
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_surface_curve_closest_points(ab, cd, 0, NULL, NULL, NULL, NULL, &d2);
	sprintf(msg, "two geodesics on a sphere: %u crossing, off the great circles' meeting %.2g, closest distance %.2g (status %d)", nx, off, (double)d2, (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && nx == 1 && off <= (GQ_F32 ? 1e-4 : 1e-8) && d2 == 0, msg);
	{
		/* a geodesic disk about A, radius 0.3 */
		qaws_surface_curve* cc = NULL;
		qaws_surface_curve const* pc[1];
		qaws_surface_path pp;
		qaws_surface_region* rg = NULL;
		double Y[2] = { 0.2, 0.5 }, Z[2] = { 0.11, 0.42 }, y3[3], z3[3], ry, rz;
		qaws_scalar dy = 0, dz = 0;
		qaws_status s1, s2;
		qaws_surface_curve_create_geodesic_arc(s, gq_uv(A[0], A[1]), (qaws_scalar)0.3, 0, (qaws_scalar)(2 * GQ_PI), 0, &cc);
		pc[0] = cc;
		pp.curves = pc; pp.curve_count = 1; pp.closed = 1;
		qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &pp, 1, NULL, 0, &rg);
		gq_at(s, Y[0], Y[1], y3); gq_at(s, Z[0], Z[1], z3);
		s1 = qaws_surface_region_distance(rg, gq_uv(Y[0], Y[1]), 0, &dy, NULL);
		s2 = qaws_surface_region_distance(rg, gq_uv(Z[0], Z[1]), 0, &dz, NULL);
		ry = R * acos(gq_dot(a3, y3) / (R * R)) - 0.3;
		rz = R * acos(gq_dot(a3, z3) / (R * R)) - 0.3;
		sprintf(msg, "signed distance to a geodesic disk, outside %.2g and inside %.2g off (status %d %d)", fabs((double)dy - ry), fabs((double)dz - rz), (int)s1, (int)s2);
		TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && dy > 0 && dz < 0 && fabs((double)dy - ry) <= (GQ_F32 ? 1e-5 : 1e-8) &&
			fabs((double)dz - rz) <= (GQ_F32 ? 2e-5 : 1e-8), msg);
		qaws_surface_region_destroy(rg);
		qaws_surface_curve_destroy(cc);
	}
	qaws_surface_curve_destroy(ab); qaws_surface_curve_destroy(cd);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

static void test_cylinder(void)
{
	/* radius 0.5, height pi: unrolled, (pi u, pi v) up to a reflection */
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_surface_revolution_desc rd;
	qaws_curve* line = NULL;
	qaws_surface* cy = NULL;
	qaws_surface_curve *e = NULL, *f = NULL;
	double a[2] = { 0.1, 0.2 }, b[2] = { 0.45, 0.35 }, c[2] = { 0.2, 0.6 }, dd[2] = { 0.55, 0.45 }, ua[2], ub[2], uc[2], ud[2], ref;
	qaws_scalar d = 0;
	qaws_status st;
	unsigned int i;
	char msg[200];
	cp[0] = (qaws_scalar)0.5; cp[1] = 0; cp[2] = (qaws_scalar)0.5; cp[3] = (qaws_scalar)-GQ_PI;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	qaws_curve_create_bezier(&bd, &line);
	memset(&rd, 0, sizeof(rd));
	rd.profile = line; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rd, &cy);
	qaws_surface_curve_create_geodesic(cy, gq_uv(a[0], a[1]), gq_uv(b[0], b[1]), 0, &e);
	qaws_surface_curve_create_geodesic(cy, gq_uv(c[0], c[1]), gq_uv(dd[0], dd[1]), 0, &f);
	for (i = 0; i < 2; i++) { ua[i] = GQ_PI * a[i]; ub[i] = GQ_PI * b[i]; uc[i] = GQ_PI * c[i]; ud[i] = GQ_PI * dd[i]; }
	ref = fmin(fmin(gq_seg(ua, ub, uc), gq_seg(ua, ub, ud)), fmin(gq_seg(uc, ud, ua), gq_seg(uc, ud, ub)));
	st = qaws_surface_curve_closest_points(e, f, 0, NULL, NULL, NULL, NULL, &d);
	sprintf(msg, "two geodesic segments on a cylinder: closest distance against the unrolled segments' %.2g (status %d)", fabs((double)d - ref), (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs((double)d - ref) <= (GQ_F32 ? 1e-5 : 1e-9), msg);
	qaws_surface_curve_destroy(e); qaws_surface_curve_destroy(f);
	qaws_surface_destroy(cy);
	qaws_curve_destroy(line);
}

int test_110_geodesic_query_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 110: geodesic queries\n");
	test_sphere();
	test_cylinder();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
