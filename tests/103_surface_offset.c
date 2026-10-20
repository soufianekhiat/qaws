/*
 * Test 103: geodesic offsets of regions and curves on surfaces
 * (qaws_surface_offset, qaws_surface_curve_create_offset)
 *
 *   - a sphere: a geodesic circle of radius rho grown and shrunk by d is
 *     the circle of radius rho +- d (area and perimeter in closed form); a
 *     geodesic triangle grown by d has the area A + L R sin(d/R) + (2 pi -
 *     A/R^2) R^2 (1 - cos(d/R)) (strips along great circles, sectors at the
 *     corners whose angles sum to the turning, by Gauss-Bonnet); a geodesic
 *     segment's tube has 2 l R sin(d/R) + 2 pi R^2 (1 - cos(d/R))
 *   - a flat sheared patch: a parallelogram grown, A + P d + pi d^2, and
 *     shrunk, A - P d + d^2 sum cot(angle / 2)
 *   - a cylinder: a circle across the seam grown and shrunk, pi (rho +- d)^2
 *     unrolled; a band round it grown by d on both sides
 *   - a wavy bicubic patch: by the Gauss lemma a geodesic circle's offset
 *     is the circle of radius rho + d about the same centre
 *   - a torus: the same across both seams
 *   - the offset curve of a circle: its points at the radius rho + d
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SOT_PI 3.14159265358979323846
#define SOT_F32 (sizeof(qaws_scalar) == 4)

static qaws_vec2 sot_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static qaws_surface* sot_revolve(double cx, double r, double a0, double a1, qaws_curve** prof)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)cx; g.radius = (qaws_scalar)r; g.angle_start = (qaws_scalar)a0; g.angle_end = (qaws_scalar)a1;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(prof, &ad);
	memset(&rd, 0, sizeof(rd));
	rd.profile = *prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&s, &rd);
	return s;
}

/* the offset of one path: its area, and its perimeter if asked */
static double sot_offset(qaws_surface const* s, qaws_surface_curve* const* cv, unsigned int n, int closed, double d,
	double* per, qaws_status* st)
{
	qaws_surface_curve const* pc[8];
	qaws_surface_path P;
	qaws_surface_region* r = NULL;
	qaws_scalar a = 0, p = 0;
	unsigned int i;
	for (i = 0; i < n; i++) pc[i] = cv[i];
	P.curves = pc; P.curve_count = n; P.closed = closed;
	*st = qaws_surface_offset(&r, s, &P, 1, QAWS_FILL_NON_ZERO, (qaws_scalar)d, 0);
	if (*st == QAWS_STATUS_OK) *st = qaws_surface_region_get_area(&a, r);
	if (*st == QAWS_STATUS_OK && per) *st = qaws_surface_region_get_perimeter(&p, r);
	if (per) *per = (double)p;
	qaws_surface_region_destroy(r);
	return (double)a;
}

static double sot_area(qaws_surface const* s, qaws_surface_curve* const* cv, unsigned int n, double* per, qaws_status* st)
{
	qaws_surface_curve const* pc[8];
	qaws_surface_path P;
	qaws_surface_region* r = NULL;
	qaws_scalar a = 0, p = 0;
	unsigned int i;
	for (i = 0; i < n; i++) pc[i] = cv[i];
	P.curves = pc; P.curve_count = n; P.closed = 1;
	*st = qaws_surface_clip(&r, s, &P, NULL, 1, 0, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO);
	if (*st == QAWS_STATUS_OK) *st = qaws_surface_region_get_area(&a, r);
	if (*st == QAWS_STATUS_OK && per) *st = qaws_surface_region_get_perimeter(&p, r);
	if (per) *per = (double)p;
	qaws_surface_region_destroy(r);
	return (double)a;
}

static qaws_surface_curve* sot_circle(qaws_surface const* s, double u, double v, double r)
{
	qaws_surface_curve* c = NULL;
	qaws_surface_curve_create_geodesic_arc(&c, s, sot_uv(u, v), (qaws_scalar)r, 0, (qaws_scalar)(2 * SOT_PI), 0);
	return c;
}

static void test_sphere(void)
{
	double R = 1.3, rho = 0.5, d = 0.2, tol = SOT_F32 ? 2e-4 : 1e-9, cap1, cap2, a1, a2, p1, p2, at, cf, as, cs;
	qaws_curve* prof = NULL;
	qaws_surface* s = sot_revolve(0, R, -SOT_PI / 2, SOT_PI / 2, &prof);
	qaws_surface_curve* c = sot_circle(s, 0.4, 0.55, rho);
	qaws_surface_curve* t[3] = { 0 };
	qaws_vec2 tri[3];
	qaws_status q1, q2, q3, q4, q5;
	unsigned int i;
	char msg[400];
	cap1 = 2 * SOT_PI * R * R * (1 - cos((rho + d) / R));
	cap2 = 2 * SOT_PI * R * R * (1 - cos((rho - d) / R));
	a1 = sot_offset(s, &c, 1, 1, d, &p1, &q1);
	a2 = sot_offset(s, &c, 1, 1, -d, &p2, &q2);
	sprintf(msg, "sphere: a geodesic circle grown and shrunk by d against the circles of radius rho +- d: areas %.2g, %.2g; perimeters %.2g, %.2g",
		fabs(a1 - cap1), fabs(a2 - cap2), fabs(p1 - 2 * SOT_PI * R * sin((rho + d) / R)), fabs(p2 - 2 * SOT_PI * R * sin((rho - d) / R)));
	TEST_ASSERT(c && !q1 && !q2 && fabs(a1 - cap1) <= tol && fabs(a2 - cap2) <= tol &&
		fabs(p1 - 2 * SOT_PI * R * sin((rho + d) / R)) <= tol && fabs(p2 - 2 * SOT_PI * R * sin((rho - d) / R)) <= tol, msg);
	tri[0] = sot_uv(0.3, 0.42); tri[1] = sot_uv(0.5, 0.40); tri[2] = sot_uv(0.4, 0.68);
	for (i = 0; i < 3; i++) qaws_surface_curve_create_geodesic(&t[i], s, tri[i], tri[(i + 1) % 3], 0);
	{
		double A, L;
		qaws_scalar l0 = 0;
		A = sot_area(s, t, 3, &L, &q3);
		at = sot_offset(s, t, 3, 1, d, NULL, &q4);
		cf = A + L * R * sin(d / R) + (2 * SOT_PI - A / (R * R)) * R * R * (1 - cos(d / R));
		as = sot_offset(s, t, 1, 0, d, NULL, &q5);
		qaws_surface_curve_get_length(&l0, t[0]);
		cs = 2 * (double)l0 * R * sin(d / R) + 2 * SOT_PI * R * R * (1 - cos(d / R));
	}
	sprintf(msg, "sphere: a geodesic triangle grown against strips and corner sectors %.2g; a geodesic segment's tube %.2g",
		fabs(at - cf), fabs(as - cs));
	TEST_ASSERT(!q3 && !q4 && !q5 && fabs(at - cf) <= tol && fabs(as - cs) <= tol, msg);
	/* the same across the seam, a corner exactly on it */
	{
		qaws_surface_curve* ts[3] = { 0 };
		qaws_vec2 tv[3];
		double A, L, ag, cfs;
		qaws_status w1, w2;
		tv[0] = sot_uv(0.92, 0.36); tv[1] = sot_uv(1.12, 0.40); tv[2] = sot_uv(1.0, 0.66);
		for (i = 0; i < 3; i++) qaws_surface_curve_create_geodesic(&ts[i], s, tv[i], tv[(i + 1) % 3], 0);
		A = sot_area(s, ts, 3, &L, &w1);
		ag = sot_offset(s, ts, 3, 1, d, NULL, &w2);
		cfs = A + L * R * sin(d / R) + (2 * SOT_PI - A / (R * R)) * R * R * (1 - cos(d / R));
		sprintf(msg, "sphere: a geodesic triangle across the seam, a corner on it, grown against the closed form %.2g", fabs(ag - cfs));
		TEST_ASSERT(!w1 && !w2 && fabs(ag - cfs) <= tol, msg);
		for (i = 0; i < 3; i++) qaws_surface_curve_destroy(ts[i]);
	}
	/* the offset curve itself: every point at the radius rho + d from the centre */
	{
		qaws_surface_curve* o = NULL;
		qaws_range rg;
		double worst = 0, cx[3], ccl = cos(SOT_PI * 0.55 - SOT_PI / 2);
		unsigned int k;
		qaws_status st = qaws_surface_curve_create_offset(&o, c, (qaws_scalar)d, 0);
		/* the circle runs counter-clockwise: its right is outward */
		if (st == QAWS_STATUS_OK) { qaws_surface_curve_destroy(o); o = NULL; st = qaws_surface_curve_create_offset(&o, c, (qaws_scalar)-d, 0); }
		cx[0] = R * ccl * cos(2 * SOT_PI * 0.4); cx[1] = R * ccl * sin(2 * SOT_PI * 0.4); cx[2] = R * sin(SOT_PI * 0.55 - SOT_PI / 2);
		if (st == QAWS_STATUS_OK)
		{
			rg = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(o));
			for (k = 0; k <= 200; k++)
			{
				qaws_vec3 p;
				double dot, ang;
				qaws_surface_curve_evaluate(NULL, &p, o, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / 200));
				dot = (p.x * cx[0] + p.y * cx[1] + p.z * cx[2]) / (R * R);
				ang = R * acos(dot > 1 ? 1 : dot);
				if (fabs(ang - (rho + d)) > worst) worst = fabs(ang - (rho + d));
			}
		}
		sprintf(msg, "sphere: the offset curve of a geodesic circle at the radius rho + d within %.2g", worst);
		TEST_ASSERT(st == QAWS_STATUS_OK && worst <= (SOT_F32 ? 1e-4 : 1e-9) && qaws_surface_curve_get_kind(o) == QAWS_SURFACE_CURVE_GEODESIC_OFFSET, msg);
		qaws_surface_curve_destroy(o);
	}
	/* overlapping operands: the union of two circles grown is the union of
	   the larger circles (the offset follows the region's boundary, not the
	   paths inside it) */
	{
		qaws_surface_curve *e1 = sot_circle(s, 0.36, 0.52, 0.4), *e2 = sot_circle(s, 0.44, 0.5, 0.35);
		qaws_surface_curve *f1 = sot_circle(s, 0.36, 0.52, 0.4 + d), *f2 = sot_circle(s, 0.44, 0.5, 0.35 + d);
		qaws_surface_curve const *pe[2], *pf[2];
		qaws_surface_path E[2], F[2];
		qaws_surface_region *rg = NULL, *rf = NULL;
		qaws_scalar ag = 0, af = 0;
		qaws_status s1, s2;
		pe[0] = e1; pe[1] = e2; pf[0] = f1; pf[1] = f2;
		E[0].curves = pe; E[0].curve_count = 1; E[0].closed = 1;
		E[1].curves = pe + 1; E[1].curve_count = 1; E[1].closed = 1;
		F[0].curves = pf; F[0].curve_count = 1; F[0].closed = 1;
		F[1].curves = pf + 1; F[1].curve_count = 1; F[1].closed = 1;
		s1 = qaws_surface_offset(&rg, s, E, 2, QAWS_FILL_NON_ZERO, (qaws_scalar)d, 0);
		s2 = qaws_surface_clip(&rf, s, F, NULL, 2, 0, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO);
		if (!s1) qaws_surface_region_get_area(&ag, rg);
		if (!s2) qaws_surface_region_get_area(&af, rf);
		sprintf(msg, "sphere: two overlapping circles' union grown against the union of the circles of radius rho + d %.2g", fabs((double)ag - (double)af));
		TEST_ASSERT(!s1 && !s2 && fabs((double)ag - (double)af) <= tol && af > 1, msg);
		qaws_surface_region_destroy(rg); qaws_surface_region_destroy(rf);
		qaws_surface_curve_destroy(e1); qaws_surface_curve_destroy(e2); qaws_surface_curve_destroy(f1); qaws_surface_curve_destroy(f2);
	}
	for (i = 0; i < 3; i++) qaws_surface_curve_destroy(t[i]);
	qaws_surface_curve_destroy(c);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

static void test_plane(void)
{
	/* the sheared plane of test 102: x = 2u + 0.6v, y = 3v */
	qaws_surface_bilinear_desc bd;
	qaws_surface* s = NULL;
	qaws_surface_curve* sq[4] = { 0 };
	qaws_vec2 corner[4];
	double d = 0.1, tol = SOT_F32 ? 1e-4 : 1e-9, A, P, th, ag, ash, bl;
	qaws_status q1, q2;
	unsigned int i;
	char msg[300];
	memset(&bd, 0, sizeof(bd));
	bd.p10.x = 2; bd.p01.x = (qaws_scalar)0.6; bd.p01.y = 3; bd.p11.x = (qaws_scalar)2.6; bd.p11.y = 3;
	qaws_surface_create_bilinear(&s, &bd);
	corner[0] = sot_uv(0.3, 0.3); corner[1] = sot_uv(0.62, 0.3); corner[2] = sot_uv(0.62, 0.58); corner[3] = sot_uv(0.3, 0.58);
	for (i = 0; i < 4; i++) qaws_surface_curve_create_geodesic(&sq[i], s, corner[i], corner[(i + 1) % 4], 0);
	/* the parallelogram: sides (0.64, 0) and (0.168, 0.84) */
	bl = sqrt(0.168 * 0.168 + 0.84 * 0.84);
	A = 0.64 * 0.84;
	P = 2 * (0.64 + bl);
	th = acos(0.168 / bl);
	ag = sot_offset(s, sq, 4, 1, d, NULL, &q1);
	ash = sot_offset(s, sq, 4, 1, -d, NULL, &q2);
	sprintf(msg, "flat sheared patch: a parallelogram grown, A + P d + pi d^2, %.2g; shrunk, A - P d + d^2 sum cot(angle/2), %.2g",
		fabs(ag - (A + P * d + SOT_PI * d * d)), fabs(ash - (A - P * d + d * d * 4 / sin(th))));
	TEST_ASSERT(!q1 && !q2 && fabs(ag - (A + P * d + SOT_PI * d * d)) <= tol && fabs(ash - (A - P * d + d * d * 4 / sin(th))) <= tol, msg);
	for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
	qaws_surface_destroy(s);
}

static void test_cylinder(void)
{
	double Rc = 0.8, H = 2.5, rho = 0.45, d = 0.15, tol = SOT_F32 ? 2e-4 : 1e-9, a1, a2, ab;
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve *circle = NULL, *k1 = NULL, *k2 = NULL;
	qaws_surface_extrusion_desc ed;
	qaws_surface* cy = NULL;
	qaws_surface_curve *c, *band[2] = { 0 };
	qaws_range ru;
	qaws_status q1, q2, q3 = QAWS_STATUS_OK;
	char msg[300];
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)Rc; g.angle_end = (qaws_scalar)(2 * SOT_PI);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&circle, &ad);
	memset(&ed, 0, sizeof(ed));
	ed.profile = circle; ed.direction.z = 1; ed.length = (qaws_scalar)H;
	qaws_surface_create_extrusion(&cy, &ed);
	ru = qaws_surface_get_u_range(cy);
	c = sot_circle(cy, ru.min_value + 0.02 * (ru.max_value - ru.min_value), 0.5, rho);
	a1 = sot_offset(cy, &c, 1, 1, d, NULL, &q1);
	a2 = sot_offset(cy, &c, 1, 1, -d, NULL, &q2);
	/* a band between two loops round it, grown by d on both sides */
	{
		qaws_bezier_desc bz;
		qaws_scalar cp[4];
		qaws_surface_curve const *p0[1], *p1[1];
		qaws_surface_path B[2];
		qaws_surface_region* r = NULL;
		qaws_scalar a = 0;
		memset(&bz, 0, sizeof(bz));
		bz.dimension = QAWS_DIMENSION_2D; bz.degree = 1; bz.control_points = cp; bz.control_point_count = 2;
		cp[0] = ru.min_value; cp[1] = (qaws_scalar)0.4; cp[2] = ru.max_value; cp[3] = (qaws_scalar)0.4;
		qaws_curve_create_bezier(&k1, &bz);
		cp[0] = ru.max_value; cp[1] = (qaws_scalar)0.6; cp[2] = ru.min_value; cp[3] = (qaws_scalar)0.6;
		qaws_curve_create_bezier(&k2, &bz);
		qaws_surface_curve_create_uv(&band[0], cy, k1);
		qaws_surface_curve_create_uv(&band[1], cy, k2);
		p0[0] = band[0]; p1[0] = band[1];
		B[0].curves = p0; B[0].curve_count = 1; B[0].closed = 1;
		B[1].curves = p1; B[1].curve_count = 1; B[1].closed = 1;
		q3 = qaws_surface_offset(&r, cy, B, 2, QAWS_FILL_NON_ZERO, (qaws_scalar)d, 0);
		if (q3 == QAWS_STATUS_OK) q3 = qaws_surface_region_get_area(&a, r);
		qaws_surface_region_destroy(r);
		ab = (double)a;
	}
	sprintf(msg, "cylinder: a circle across the seam grown and shrunk against pi (rho +- d)^2: %.2g, %.2g; a band round it grown on both sides %.2g",
		fabs(a1 - SOT_PI * (rho + d) * (rho + d)), fabs(a2 - SOT_PI * (rho - d) * (rho - d)), fabs(ab - 2 * SOT_PI * Rc * (0.2 * H + 2 * d)));
	TEST_ASSERT(c && !q1 && !q2 && !q3 && fabs(a1 - SOT_PI * (rho + d) * (rho + d)) <= tol && fabs(a2 - SOT_PI * (rho - d) * (rho - d)) <= tol &&
		fabs(ab - 2 * SOT_PI * Rc * (0.2 * H + 2 * d)) <= tol, msg);
	qaws_surface_curve_destroy(c);
	qaws_surface_curve_destroy(band[0]); qaws_surface_curve_destroy(band[1]);
	qaws_curve_destroy(k1); qaws_curve_destroy(k2);
	qaws_surface_destroy(cy);
	qaws_curve_destroy(circle);
}

static void test_gauss_lemma(void)
{
	static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
	qaws_vec3 cp[16];
	qaws_surface_bezier_desc bd;
	qaws_surface* s = NULL;
	qaws_curve* prof = NULL;
	qaws_surface* tor;
	double rho = 0.35, d = 0.12, tol = SOT_F32 ? 2e-4 : 1e-8, ag, as, ref1, ref2, at, reft;
	qaws_surface_curve *c, *c1, *c2, *ct, *ct1;
	qaws_status q[6];
	unsigned int i;
	char msg[300];
	for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
	memset(&bd, 0, sizeof(bd));
	bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
	qaws_surface_create_bezier(&s, &bd);
	c = sot_circle(s, 0.5, 0.5, rho);
	c1 = sot_circle(s, 0.5, 0.5, rho + d);
	c2 = sot_circle(s, 0.5, 0.5, rho - d);
	ag = sot_offset(s, &c, 1, 1, d, NULL, &q[0]);
	as = sot_offset(s, &c, 1, 1, -d, NULL, &q[1]);
	ref1 = sot_area(s, &c1, 1, NULL, &q[2]);
	ref2 = sot_area(s, &c2, 1, NULL, &q[3]);
	tor = sot_revolve(1.0, 0.42, 0, 2 * SOT_PI, &prof);
	ct = sot_circle(tor, 0.0, 0.0, 0.25);
	ct1 = sot_circle(tor, 0.0, 0.0, 0.25 + 0.1);
	at = sot_offset(tor, &ct, 1, 1, 0.1, NULL, &q[4]);
	reft = sot_area(tor, &ct1, 1, NULL, &q[5]);
	sprintf(msg, "Gauss lemma: a geodesic circle's offset against the circle of radius rho +- d, on a wavy patch %.2g, %.2g; on a torus across both seams %.2g",
		fabs(ag - ref1), fabs(as - ref2), fabs(at - reft));
	TEST_ASSERT(c && c1 && c2 && ct && ct1 && !q[0] && !q[1] && !q[2] && !q[3] && !q[4] && !q[5] &&
		fabs(ag - ref1) <= tol && fabs(as - ref2) <= tol && fabs(at - reft) <= tol, msg);
	qaws_surface_curve_destroy(c); qaws_surface_curve_destroy(c1); qaws_surface_curve_destroy(c2);
	qaws_surface_curve_destroy(ct); qaws_surface_curve_destroy(ct1);
	qaws_surface_destroy(s);
	qaws_surface_destroy(tor);
	qaws_curve_destroy(prof);
}

int test_103_surface_offset_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 103: geodesic offsets on surfaces\n");
	test_sphere();
	test_plane();
	test_cylinder();
	test_gauss_lemma();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
