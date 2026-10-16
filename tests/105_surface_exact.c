/*
 * Test 105: Booleans of regions on surfaces with exact topology
 * (qaws_surface_clip_exact: the (u, v) images through qaws_clip_exact, the
 * cuts at seams and the unions there too)
 *
 *   - against the float path: two geodesic circles on a sphere (union,
 *     intersection, difference, xor), circles across the sphere's seam,
 *     a circle across both seams of a torus: the same areas
 *   - curved boundaries cut at the seam: five circles and their unions
 *     certified, with the float areas
 *   - degenerate inputs drawn in (u, v), where only exact decisions hold:
 *     a region with itself (union and intersection the region, difference
 *     empty), two squares sharing an edge (union one loop, its area the sum),
 *     two squares touching at one corner (union two loops), a square sharing
 *     part of an edge with another across the seam
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include "qaws_bezier.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SXT_PI 3.14159265358979323846
#define SXT_F32 (sizeof(qaws_scalar) == 4)

static qaws_vec2 sxt_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static qaws_surface* sxt_revolve(double cx, double r, double a0, double a1, qaws_curve** prof)
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
	qaws_surface_create_revolution(&rd, &s);
	return s;
}

static double sxt_area(qaws_surface const* s, int exact, qaws_clip_type type, qaws_surface_path const* a, unsigned int na,
	qaws_surface_path const* b, unsigned int nb, unsigned int* loops, qaws_status* st)
{
	qaws_surface_region* r = NULL;
	qaws_scalar area = 0;
	*st = exact ? qaws_surface_clip_exact(s, type, QAWS_FILL_NON_ZERO, a, na, b, nb, NULL, &r)
		: qaws_surface_clip(s, type, QAWS_FILL_NON_ZERO, a, na, b, nb, &r);
	if (*st == QAWS_STATUS_OK)
		*st = qaws_surface_region_get_area(r, &area);
	if (loops)
		*loops = r ? qaws_surface_region_get_path_count(r) : 0;
	qaws_surface_region_destroy(r);
	return (double)area;
}

/* a square drawn in (u, v): four lines */
typedef struct sxt_square
{
	qaws_curve* k[4];
	qaws_surface_curve* c[4];
	qaws_surface_curve const* pc[4];
	qaws_surface_path path;
} sxt_square;

static void sxt_square_make(sxt_square* q, qaws_surface const* s, double u0, double v0, double u1, double v1)
{
	double x[5] = { u0, u1, u1, u0, u0 }, y[5] = { v0, v0, v1, v1, v0 };
	unsigned int i;
	for (i = 0; i < 4; i++)
	{
		qaws_bezier_desc bz;
		qaws_scalar cp[4];
		cp[0] = (qaws_scalar)x[i]; cp[1] = (qaws_scalar)y[i]; cp[2] = (qaws_scalar)x[i + 1]; cp[3] = (qaws_scalar)y[i + 1];
		memset(&bz, 0, sizeof(bz));
		bz.dimension = QAWS_DIMENSION_2D; bz.degree = 1; bz.control_points = cp; bz.control_point_count = 2;
		q->k[i] = NULL; q->c[i] = NULL;
		qaws_curve_create_bezier(&q->k[i], &bz);
		qaws_surface_curve_create_uv(s, q->k[i], &q->c[i]);
		q->pc[i] = q->c[i];
	}
	q->path.curves = q->pc; q->path.curve_count = 4; q->path.closed = 1;
}

static void sxt_square_free(sxt_square* q)
{
	unsigned int i;
	for (i = 0; i < 4; i++) { qaws_surface_curve_destroy(q->c[i]); qaws_curve_destroy(q->k[i]); }
}

static void test_against_float(void)
{
	double R = 1.3, tol = SXT_F32 ? 2e-4 : 1e-8, worst = 0;
	qaws_curve* prof = NULL;
	qaws_surface* s = sxt_revolve(0, R, -SXT_PI / 2, SXT_PI / 2, &prof);
	static qaws_clip_type const ops[4] = { QAWS_CLIP_UNION, QAWS_CLIP_INTERSECTION, QAWS_CLIP_DIFFERENCE, QAWS_CLIP_XOR };
	qaws_surface_curve *ca = NULL, *cb = NULL, *cs = NULL, *ct = NULL;
	qaws_surface_curve const *pa[1], *pb[1], *ps[1];
	qaws_surface_path A, B, Sp;
	qaws_status st, bad = QAWS_STATUS_OK;
	unsigned int k;
	char msg[300];
	qaws_surface_curve_create_geodesic_arc(s, sxt_uv(0.3, 0.55), (qaws_scalar)0.6, 0, (qaws_scalar)(2 * SXT_PI), 0, &ca);
	qaws_surface_curve_create_geodesic_arc(s, sxt_uv(0.38, 0.5), (qaws_scalar)0.5, 0, (qaws_scalar)(2 * SXT_PI), 0, &cb);
	qaws_surface_curve_create_geodesic_arc(s, sxt_uv(1.3, 0.5), (qaws_scalar)0.5, 0, (qaws_scalar)(2 * SXT_PI), 0, &cs);
	pa[0] = ca; pb[0] = cb; ps[0] = cs;
	A.curves = pa; A.curve_count = 1; A.closed = 1;
	B.curves = pb; B.curve_count = 1; B.closed = 1;
	Sp.curves = ps; Sp.curve_count = 1; Sp.closed = 1;
	for (k = 0; k < 4; k++)
	{
		double af = sxt_area(s, 0, ops[k], &A, 1, &B, 1, NULL, &st), ae;
		if (st) bad = st;
		ae = sxt_area(s, 1, ops[k], &A, 1, &B, 1, NULL, &st);
		if (st) bad = st;
		if (fabs(af - ae) > worst) worst = fabs(af - ae);
	}
	{
		double af = sxt_area(s, 0, QAWS_CLIP_UNION, &Sp, 1, &A, 1, NULL, &st), ae;
		if (st) bad = st;
		ae = sxt_area(s, 1, QAWS_CLIP_UNION, &Sp, 1, &A, 1, NULL, &st);
		if (st) bad = st;
		if (fabs(af - ae) > worst) worst = fabs(af - ae);
	}
	{
		qaws_curve* tp = NULL;
		qaws_surface* tor = sxt_revolve(1.0, 0.42, 0, 2 * SXT_PI, &tp);
		qaws_surface_curve const* pt[1];
		qaws_surface_path T;
		double af, ae;
		qaws_surface_curve_create_geodesic_arc(tor, sxt_uv(0.4, 0.3), (qaws_scalar)0.3, 0, (qaws_scalar)(2 * SXT_PI), 0, &ct);
		pt[0] = ct;
		T.curves = pt; T.curve_count = 1; T.closed = 1;
		af = sxt_area(tor, 0, QAWS_CLIP_NONE, &T, 1, NULL, 0, NULL, &st);
		if (st) bad = st;
		ae = sxt_area(tor, 1, QAWS_CLIP_NONE, &T, 1, NULL, 0, NULL, &st);
		if (st) bad = st;
		if (fabs(af - ae) > worst) worst = fabs(af - ae);
		qaws_surface_curve_destroy(ct);
		qaws_surface_destroy(tor);
		qaws_curve_destroy(tp);
	}
	sprintf(msg, "exact against float: two circles' Booleans, a circle given a period on (moved, not cut), a circle on a torus agree to %.2g (status %d)",
		worst, (int)bad);
	TEST_ASSERT(bad == QAWS_STATUS_OK && worst <= tol, msg);
	/* curved boundaries cut at the seam (a cut piece's end a few ulps past
	   its curve's start is that start on the parameter lattice): each circle
	   and its union with a shifted copy certified, the float areas */
	{
		static double const cs[5][3] = { { 0.02, 0.5, 0.5 }, { 0.0, 0.5, 0.4 }, { 0.98, 0.45, 0.6 }, { 0.05, 0.3, 0.3 }, { 0.99, 0.6, 0.25 } };
		double off = 0;
		qaws_status sbad = QAWS_STATUS_OK;
		unsigned int z;
		for (z = 0; z < 5; z++)
		{
			qaws_surface_curve *cx = NULL, *cy = NULL;
			qaws_surface_curve const *px[1], *py[1];
			qaws_surface_path X, Y;
			double af, ae, uf, ue;
			qaws_status sf, se, tf, te;
			qaws_surface_curve_create_geodesic_arc(s, sxt_uv(cs[z][0], cs[z][1]), (qaws_scalar)cs[z][2], 0, (qaws_scalar)(2 * SXT_PI), 0, &cx);
			qaws_surface_curve_create_geodesic_arc(s, sxt_uv(cs[z][0] + 0.06, cs[z][1] + 0.03), (qaws_scalar)cs[z][2], 0, (qaws_scalar)(2 * SXT_PI), 0, &cy);
			px[0] = cx; py[0] = cy;
			X.curves = px; X.curve_count = 1; X.closed = 1;
			Y.curves = py; Y.curve_count = 1; Y.closed = 1;
			af = sxt_area(s, 0, QAWS_CLIP_NONE, &X, 1, NULL, 0, NULL, &sf);
			ae = sxt_area(s, 1, QAWS_CLIP_NONE, &X, 1, NULL, 0, NULL, &se);
			uf = sxt_area(s, 0, QAWS_CLIP_UNION, &X, 1, &Y, 1, NULL, &tf);
			ue = sxt_area(s, 1, QAWS_CLIP_UNION, &X, 1, &Y, 1, NULL, &te);
			if (sf || se || tf || te) sbad = se ? se : te ? te : sf ? sf : tf;
			if (fabs(ae - af) > off) off = fabs(ae - af);
			if (fabs(ue - uf) > off) off = fabs(ue - uf);
			qaws_surface_curve_destroy(cx);
			qaws_surface_curve_destroy(cy);
		}
		sprintf(msg, "exact, five circles across the seam and their unions with shifted copies: status %d, areas %.2g off the float path's", (int)sbad, off);
		TEST_ASSERT(sbad == QAWS_STATUS_OK && off <= tol, msg);
	}
	qaws_surface_curve_destroy(ca); qaws_surface_curve_destroy(cb); qaws_surface_curve_destroy(cs);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

static void test_degenerate(void)
{
	double R = 1.3, tol = SXT_F32 ? 2e-5 : 1e-12;
	qaws_curve* prof = NULL;
	qaws_surface* s = sxt_revolve(0, R, -SXT_PI / 2, SXT_PI / 2, &prof);
	sxt_square a, a2, b, c, d1, d2;
	qaws_status s1, s2, s3, s4, s5, s6;
	unsigned int l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0;
	double aa, au, ai, ad, ab, ac, aab, aac, aw, ad1, ad2, aseam;
	char msg[400];
	/* (u, v) squares: a; a again; b sharing a's right edge; c touching a's
	   corner; across the seam, d1 at the right edge sharing part of an edge
	   with d2 at the left (through the seam) */
	sxt_square_make(&a, s, 0.25, 0.375, 0.375, 0.5);
	sxt_square_make(&a2, s, 0.25, 0.375, 0.375, 0.5);
	sxt_square_make(&b, s, 0.375, 0.375, 0.5, 0.5);
	sxt_square_make(&c, s, 0.375, 0.5, 0.5, 0.625);
	sxt_square_make(&d1, s, 0.875, 0.25, 1.0, 0.375);
	sxt_square_make(&d2, s, 0.0, 0.3125, 0.125, 0.4375);
	aa = sxt_area(s, 1, QAWS_CLIP_NONE, &a.path, 1, NULL, 0, NULL, &s1);
	au = sxt_area(s, 1, QAWS_CLIP_UNION, &a.path, 1, &a2.path, 1, &l1, &s2);
	ai = sxt_area(s, 1, QAWS_CLIP_INTERSECTION, &a.path, 1, &a2.path, 1, NULL, &s3);
	ad = sxt_area(s, 1, QAWS_CLIP_DIFFERENCE, &a.path, 1, &a2.path, 1, &l2, &s4);
	sprintf(msg, "exact, a region with itself: union %.2g and intersection %.2g off its area, %u loop; difference %.2g with %u loops",
		fabs(au - aa), fabs(ai - aa), l1, ad, l2);
	TEST_ASSERT(!s1 && !s2 && !s3 && !s4 && fabs(au - aa) <= tol && fabs(ai - aa) <= tol && l1 == 1 && l2 == 0 && ad == 0, msg);
	ab = sxt_area(s, 1, QAWS_CLIP_NONE, &b.path, 1, NULL, 0, NULL, &s1);
	ac = sxt_area(s, 1, QAWS_CLIP_NONE, &c.path, 1, NULL, 0, NULL, &s2);
	aab = sxt_area(s, 1, QAWS_CLIP_UNION, &a.path, 1, &b.path, 1, &l3, &s3);
	aac = sxt_area(s, 1, QAWS_CLIP_UNION, &a.path, 1, &c.path, 1, &l4, &s4);
	aw = sxt_area(s, 1, QAWS_CLIP_INTERSECTION, &a.path, 1, &b.path, 1, &l5, &s5);
	sprintf(msg, "exact: squares sharing an edge unite into %u loop (area %.2g off the sum), meet in nothing (%u loops); squares touching at a corner unite into %u loops (%.2g off)",
		l3, fabs(aab - aa - ab), l5, l4, fabs(aac - aa - ac));
	TEST_ASSERT(!s1 && !s2 && !s3 && !s4 && !s5 && l3 == 1 && fabs(aab - aa - ab) <= tol && l5 == 0 && aw == 0 && l4 == 2 && fabs(aac - aa - ac) <= tol, msg);
	/* across the seam: d1 = [0.875, 1] x [0.25, 0.375] and d2 = [0, 0.125] x
	   [0.3125, 0.4375] share the seam's stretch [0.3125, 0.375] */
	ad1 = sxt_area(s, 1, QAWS_CLIP_NONE, &d1.path, 1, NULL, 0, NULL, &s1);
	ad2 = sxt_area(s, 1, QAWS_CLIP_NONE, &d2.path, 1, NULL, 0, NULL, &s2);
	aseam = sxt_area(s, 1, QAWS_CLIP_UNION, &d1.path, 1, &d2.path, 1, &l1, &s6);
	sprintf(msg, "exact: squares meeting across the seam along part of an edge unite into %u loops in (u, v) (area %.2g off the sum)", l1, fabs(aseam - ad1 - ad2));
	TEST_ASSERT(!s1 && !s2 && !s6 && fabs(aseam - ad1 - ad2) <= tol && l1 == 2, msg);
	sxt_square_free(&a); sxt_square_free(&a2); sxt_square_free(&b); sxt_square_free(&c); sxt_square_free(&d1); sxt_square_free(&d2);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

int test_105_surface_exact_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 105: Booleans on surfaces with exact topology\n");
	test_against_float();
	test_degenerate();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
