/*
 * Test 101: geodesics and curves on surfaces (qaws_surface_region.h)
 *
 *   - a flat bilinear patch with a non-affine parameterisation: geodesics
 *     are straight in space (their length the distance, their midpoint the
 *     chord's)
 *   - a sphere (a revolution): a geodesic's length is the great circle's,
 *     also across the seam; every point of a geodesic circle is at the
 *     geodesic distance r from its centre, and its length is
 *     2 pi R sin(r / R); a meridian drawn in (u, v) has the length R dphi
 *   - a cylinder (an extrusion of a circle): geodesics are the unrolled
 *     plane's straight lines, round the seam too
 *   - a bicubic Bezier patch against Mathematica
 *     (tests/reference/101_surface_geodesic.wls): the exponential map, a
 *     geodesic between two points, a geodesic circle's points
 *   - the metric and Christoffel symbols on the sphere against closed forms
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include "reference/101_surface_geodesic.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SGT_PI 3.14159265358979323846
#define SGT_F32 (sizeof(qaws_scalar) == 4)

static qaws_surface* sgt_sphere(double R, qaws_curve** prof)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)R; g.angle_start = (qaws_scalar)(-SGT_PI / 2); g.angle_end = (qaws_scalar)(SGT_PI / 2);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, prof);
	memset(&rd, 0, sizeof(rd));
	rd.profile = *prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rd, &s);
	return s;
}

static qaws_vec3 sgt_point(qaws_surface const* s, qaws_vec2 uv)
{
	qaws_surface_eval_result r;
	qaws_range ru = qaws_surface_get_u_range(s);
	double p = (double)ru.max_value - (double)ru.min_value, u = fmod((double)uv.x - ru.min_value, p);
	if (u < 0) u += p;
	qaws_surface_evaluate(s, (qaws_scalar)(ru.min_value + u), uv.y, QAWS_SURFACE_EVAL_POSITION, &r);
	return r.position;
}

static qaws_vec2 sgt_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static qaws_vec2 sgt_end(qaws_surface_curve const* c, int last)
{
	qaws_vec2 uv;
	qaws_range rg = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(c));
	qaws_surface_curve_evaluate(c, last ? rg.max_value : rg.min_value, &uv, NULL);
	return uv;
}

static void test_plane(void)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	double pts[6][2] = { { 0.1, 0.2 }, { 0.9, 0.7 }, { 0.2, 0.9 }, { 0.8, 0.1 }, { 0.5, 0.05 }, { 0.45, 0.95 } };
	double worst_len = 0, worst_mid = 0, worst_end = 0;
	unsigned int k;
	char msg[300];
	/* a flat quadrilateral, not a parallelogram: u and v bend in space */
	memset(&d, 0, sizeof(d));
	d.p00.x = 0; d.p00.y = 0;
	d.p10.x = 3; d.p10.y = (qaws_scalar)0.4;
	d.p01.x = (qaws_scalar)0.5; d.p01.y = 2;
	d.p11.x = (qaws_scalar)2.2; d.p11.y = (qaws_scalar)2.6;
	qaws_surface_create_bilinear(&d, &s);
	for (k = 0; k < 3; k++)
	{
		qaws_surface_curve* c = NULL;
		qaws_scalar L = 0;
		qaws_vec3 a, b, m;
		qaws_range rg;
		qaws_status st = qaws_surface_curve_create_geodesic(s, sgt_uv(pts[2 * k][0], pts[2 * k][1]), sgt_uv(pts[2 * k + 1][0], pts[2 * k + 1][1]), 0, &c);
		if (st != QAWS_STATUS_OK) { worst_len = 1e300; continue; }
		qaws_surface_curve_get_length(c, &L);
		rg = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(c));
		qaws_surface_curve_evaluate(c, rg.min_value, NULL, &a);
		qaws_surface_curve_evaluate(c, rg.max_value, NULL, &b);
		qaws_surface_curve_evaluate(c, (qaws_scalar)(0.5 * (rg.min_value + rg.max_value)), NULL, &m);
		{
			double dl = fabs((double)L - sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z)));
			double dm = sqrt((m.x - 0.5 * (a.x + b.x)) * (m.x - 0.5 * (a.x + b.x)) + (m.y - 0.5 * (a.y + b.y)) * (m.y - 0.5 * (a.y + b.y)));
			qaws_vec2 e = sgt_end(c, 1);
			double de = hypot(e.x - pts[2 * k + 1][0], e.y - pts[2 * k + 1][1]);
			if (dl > worst_len) worst_len = dl;
			if (dm > worst_mid) worst_mid = dm;
			if (de > worst_end) worst_end = de;
		}
		qaws_surface_curve_destroy(c);
	}
	sprintf(msg, "flat bilinear patch: geodesics straight in space (length against the distance %.2g, midpoint against the chord's %.2g, end %.2g)",
		worst_len, worst_mid, worst_end);
	TEST_ASSERT(worst_len <= (SGT_F32 ? 1e-4 : 1e-9) && worst_mid <= (SGT_F32 ? 1e-4 : 1e-9) && worst_end <= (SGT_F32 ? 1e-4 : 1e-9), msg);
	qaws_surface_destroy(s);
}

static void test_sphere(void)
{
	double R = 1.3, tol = SGT_F32 ? 2e-4 : 1e-8;
	qaws_curve* prof = NULL;
	qaws_surface* s = sgt_sphere(R, &prof);
	char msg[300];
	/* great circles, one across the seam */
	{
		double ab[3][4] = { { 0.1, 0.6, 0.35, 0.3 }, { 0.95, 0.5, 1.1, 0.45 }, { 0.6, 0.2, 0.62, 0.8 } };
		double worst = 0;
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			qaws_surface_curve* c = NULL;
			qaws_scalar L = 0;
			qaws_vec3 a, b;
			if (qaws_surface_curve_create_geodesic(s, sgt_uv(ab[k][0], ab[k][1]), sgt_uv(ab[k][2], ab[k][3]), 0, &c) != QAWS_STATUS_OK)
			{
				worst = 1e300;
				continue;
			}
			qaws_surface_curve_get_length(c, &L);
			a = sgt_point(s, sgt_uv(ab[k][0], ab[k][1]));
			b = sgt_point(s, sgt_uv(ab[k][2], ab[k][3]));
			{
				double dot = ((double)a.x * b.x + (double)a.y * b.y + (double)a.z * b.z) / (R * R);
				if (dot > 1) dot = 1;
				if (fabs((double)L - R * acos(dot)) > worst) worst = fabs((double)L - R * acos(dot));
			}
			qaws_surface_curve_destroy(c);
		}
		sprintf(msg, "sphere: three geodesics (one across the seam) against the great circles within %.2g", worst);
		TEST_ASSERT(worst <= tol, msg);
	}
	/* a geodesic circle: every point at the angle r / R from the centre, its length 2 pi R sin(r / R) */
	{
		qaws_surface_curve* c = NULL;
		qaws_vec2 cc = sgt_uv(0.3, 0.55);
		double r = 0.7, worst = 0, dl = 1e300;
		qaws_vec3 pc = sgt_point(s, cc);
		qaws_status st = qaws_surface_curve_create_geodesic_arc(s, cc, (qaws_scalar)r, 0, (qaws_scalar)(2 * SGT_PI), 0, &c);
		if (st == QAWS_STATUS_OK)
		{
			qaws_range rg = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(c));
			qaws_scalar L = 0;
			unsigned int k;
			for (k = 0; k <= 100; k++)
			{
				qaws_vec3 p;
				double dot;
				qaws_surface_curve_evaluate(c, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / 100), NULL, &p);
				dot = ((double)p.x * pc.x + (double)p.y * pc.y + (double)p.z * pc.z) / (R * R);
				if (dot > 1) dot = 1;
				if (fabs(R * acos(dot) - r) > worst) worst = fabs(R * acos(dot) - r);
			}
			qaws_surface_curve_get_length(c, &L);
			dl = fabs((double)L - 2 * SGT_PI * R * sin(r / R));
		}
		sprintf(msg, "sphere: a geodesic circle (r 0.7): its points at distance r within %.2g, its length against 2 pi R sin(r / R) within %.2g", worst, dl);
		TEST_ASSERT(st == QAWS_STATUS_OK && worst <= tol && dl <= tol, msg);
		qaws_surface_curve_destroy(c);
	}
	/* a meridian drawn in (u, v): length R dphi */
	{
		qaws_bezier_desc bd;
		qaws_scalar cp[4] = { (qaws_scalar)0.2, (qaws_scalar)0.25, (qaws_scalar)0.2, (qaws_scalar)0.75 };
		qaws_curve* line = NULL;
		qaws_surface_curve* c = NULL;
		qaws_scalar L = 0;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
		qaws_curve_create_bezier(&bd, &line);
		qaws_surface_curve_create_uv(s, line, &c);
		qaws_surface_curve_get_length(c, &L);
		sprintf(msg, "sphere: a meridian drawn in (u, v): length against R pi / 2 within %.2g", fabs((double)L - R * SGT_PI / 2));
		TEST_ASSERT(c && qaws_surface_curve_get_kind(c) == QAWS_SURFACE_CURVE_UV && fabs((double)L - R * SGT_PI / 2) <= tol, msg);
		qaws_surface_curve_destroy(c);
		qaws_curve_destroy(line);
	}
	/* the metric: E = (2 pi R cos phi)^2, F = 0, G = (pi R)^2; G1_12 = -pi tan(phi) (dphi/dv = pi) */
	{
		qaws_scalar g[3], ch[6];
		double v = 0.3, phi = -SGT_PI / 2 + SGT_PI * v, e1, e2, e3;
		qaws_surface_get_metric(s, (qaws_scalar)0.4, (qaws_scalar)v, g, ch);
		e1 = fabs((double)g[0] - pow(2 * SGT_PI * R * cos(phi), 2)) / pow(2 * SGT_PI * R, 2);
		e2 = fabs((double)g[2] - pow(SGT_PI * R, 2)) / pow(SGT_PI * R, 2);
		e3 = fabs((double)ch[1] + SGT_PI * tan(phi));
		sprintf(msg, "sphere: the metric (E %.2g, G %.2g relative) and G1_12 = -pi tan(phi) (%.2g) against closed forms", e1, e2, e3);
		TEST_ASSERT(e1 <= (SGT_F32 ? 1e-5 : 1e-12) && e2 <= (SGT_F32 ? 1e-5 : 1e-12) && e3 <= (SGT_F32 ? 1e-4 : 1e-10), msg);
	}
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

static void test_cylinder(void)
{
	double R = 0.8, H = 2.5, worst = 0;
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* circle = NULL;
	qaws_surface_extrusion_desc ed;
	qaws_surface* s = NULL;
	qaws_range ru;
	double ab[3][4] = { { 0.1, 0.1, 0.4, 0.9 }, { 0.9, 0.2, 1.15, 0.7 }, { 0.3, 0.5, -0.2, 0.6 } };
	unsigned int k;
	char msg[300];
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)R; g.angle_start = 0; g.angle_end = (qaws_scalar)(2 * SGT_PI);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, &circle);
	memset(&ed, 0, sizeof(ed));
	ed.profile = circle; ed.direction.z = 1; ed.length = (qaws_scalar)H;
	qaws_surface_create_extrusion(&ed, &s);
	ru = qaws_surface_get_u_range(s);
	for (k = 0; k < 3; k++)
	{
		/* u is the circle's arc length; v the height fraction? measured from the surface itself */
		qaws_surface_curve* c = NULL;
		qaws_scalar L = 0;
		qaws_surface_eval_result ra, rb;
		qaws_vec2 a = sgt_uv(ab[k][0] * ru.max_value, ab[k][1]), b = sgt_uv(ab[k][2] * ru.max_value, ab[k][3]);
		double du, dh;
		if (qaws_surface_curve_create_geodesic(s, a, b, 0, &c) != QAWS_STATUS_OK) { worst = 1e300; continue; }
		qaws_surface_curve_get_length(c, &L);
		/* unrolled: the arc length round (u is unwrapped in the covering plane) and the height */
		qaws_surface_evaluate(s, (qaws_scalar)0, a.y, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU, &ra);
		qaws_surface_evaluate(s, (qaws_scalar)0, b.y, QAWS_SURFACE_EVAL_POSITION, &rb);
		du = ((double)b.x - a.x) * sqrt((double)ra.du.x * ra.du.x + (double)ra.du.y * ra.du.y + (double)ra.du.z * ra.du.z);
		dh = (double)rb.position.z - ra.position.z;
		if (fabs((double)L - hypot(du, dh)) > worst) worst = fabs((double)L - hypot(du, dh));
		qaws_surface_curve_destroy(c);
	}
	sprintf(msg, "cylinder: three geodesics (two round the seam) against the unrolled distances within %.2g", worst);
	TEST_ASSERT(worst <= (SGT_F32 ? 2e-4 : 1e-8), msg);
	qaws_surface_destroy(s);
	qaws_curve_destroy(circle);
}

static void test_patch(void)
{
	qaws_vec3 cp[16];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	double tol = SGT_F32 ? 5e-4 : 1e-9, worst_exp = 0, dlen = 1e300, worst_circle = 0;
	unsigned int i, k;
	char msg[300];
	for (i = 0; i < 16; i++)
	{
		cp[i].x = (qaws_scalar)ref_patch[3 * i]; cp[i].y = (qaws_scalar)ref_patch[3 * i + 1]; cp[i].z = (qaws_scalar)ref_patch[3 * i + 2];
	}
	memset(&d, 0, sizeof(d));
	d.u_degree = 3; d.v_degree = 3; d.control_points = cp; d.u_point_count = 4; d.v_point_count = 4;
	qaws_surface_create_bezier(&d, &s);
	/* the exponential map: the angle's direction (from S_u, a quarter turn about S_u x S_v) */
	for (k = 0; k < 3 && s; k++)
	{
		qaws_vec2 p = sgt_uv(ref_exp_in[4 * k], ref_exp_in[4 * k + 1]), out, dir;
		qaws_scalar g[3];
		double a = ref_exp_in[4 * k + 2], L = ref_exp_in[4 * k + 3], E, F, G, det;
		qaws_surface_get_metric(s, p.x, p.y, g, NULL);
		E = g[0]; F = g[1]; G = g[2]; det = E * G - F * F;
		dir.x = (qaws_scalar)(cos(a) / sqrt(E) - sin(a) * F / sqrt(E * det));
		dir.y = (qaws_scalar)(sin(a) * E / sqrt(E * det));
		if (qaws_surface_exp(s, p, dir, (qaws_scalar)L, 0, &out, NULL) != QAWS_STATUS_OK) { worst_exp = 1e300; continue; }
		if (hypot(out.x - ref_exp_out[2 * k], out.y - ref_exp_out[2 * k + 1]) > worst_exp)
			worst_exp = hypot(out.x - ref_exp_out[2 * k], out.y - ref_exp_out[2 * k + 1]);
	}
	{
		qaws_surface_curve* c = NULL;
		qaws_scalar L = 0;
		if (qaws_surface_curve_create_geodesic(s, sgt_uv(ref_geodesic_ends[0], ref_geodesic_ends[1]), sgt_uv(ref_geodesic_ends[2], ref_geodesic_ends[3]), 0, &c) == QAWS_STATUS_OK)
		{
			qaws_surface_curve_get_length(c, &L);
			dlen = fabs((double)L - ref_geodesic_length[0]);
		}
		qaws_surface_curve_destroy(c);
	}
	{
		/* the circle: the arc's parameter is uniform in the angle, so angle k pi / 4 is at k / 8 of its range */
		qaws_surface_curve* c = NULL;
		if (qaws_surface_curve_create_geodesic_arc(s, sgt_uv(ref_circle_centre[0], ref_circle_centre[1]), (qaws_scalar)ref_circle_centre[2],
			0, (qaws_scalar)(2 * SGT_PI), 0, &c) == QAWS_STATUS_OK)
		{
			qaws_range rg = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(c));
			for (k = 0; k < 8; k++)
			{
				qaws_vec2 uv;
				qaws_surface_curve_evaluate(c, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / 8), &uv, NULL);
				if (hypot(uv.x - ref_circle_points[2 * k], uv.y - ref_circle_points[2 * k + 1]) > worst_circle)
					worst_circle = hypot(uv.x - ref_circle_points[2 * k], uv.y - ref_circle_points[2 * k + 1]);
			}
		}
		else
			worst_circle = 1e300;
		qaws_surface_curve_destroy(c);
	}
	sprintf(msg, "bicubic Bezier patch against Mathematica: the exponential map %.2g, a geodesic's length %.2g, a geodesic circle's points %.2g",
		worst_exp, dlen, worst_circle);
	TEST_ASSERT(s && worst_exp <= tol && dlen <= tol && worst_circle <= tol, msg);
	qaws_surface_destroy(s);
}

int test_101_surface_geodesic_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 101: geodesics and curves on surfaces\n");
	test_plane();
	test_sphere();
	test_cylinder();
	test_patch();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
