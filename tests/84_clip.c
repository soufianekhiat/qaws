/*
 * Test 84: Boolean operations on regions of curves (qaws_clip)
 *
 *   - two overlapping squares: every clip type, areas and path counts
 *   - squares sharing an edge (one rectangle, collinear corners dropped or
 *     kept), touching at a corner (two paths), identical (union one, the
 *     rest empty)
 *   - a pentagram under every fill rule
 *   - nesting: a square with a hole and an island in it (depths 0, 1, 2)
 *   - two circles: lens areas in closed form, pieces that are NURBS
 *   - a circle minus a square: mixed curves
 *   - open subjects: a line through a square, kept inside or outside
 *   - the z callback at crossings
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define CLT_PI 3.14159265358979323846

#if QAWS_SCALAR_IS_FLOAT
#define CLT_TOL 2e-3
#else
#define CLT_TOL 1e-9
#endif

static qaws_curve* clt_poly(double const* xy, unsigned int n)
{
	qaws_scalar p[64];
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	qaws_curve_create_polyline_2d(p, n, 1, &c);
	return c;
}

static qaws_curve* clt_rect(double x0, double y0, double x1, double y1, int ccw)
{
	double a[8] = { x0, y0, x1, y0, x1, y1, x0, y1 }, b[8] = { x0, y0, x0, y1, x1, y1, x1, y0 };
	return clt_poly(ccw ? a : b, 4);
}

static double clt_area(qaws_clip_result const* r)
{
	unsigned int i, n = qaws_clip_result_get_path_count(r);
	double sum = 0.0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d p;
		qaws_scalar a = 0;
		qaws_clip_result_get_path(r, i, &p);
		qaws_path_compute_area_2d(&p, &a);
		sum += a;
	}
	return sum;
}

static unsigned int clt_vertices(qaws_clip_result const* r, unsigned int i)
{
	qaws_clip_vertex const* v = NULL;
	unsigned int n = 0;
	qaws_clip_result_get_vertices(r, 0, i, &v, &n);
	return n;
}

static qaws_clip_result* clt_run(qaws_clip_type ct, qaws_fill_rule fr, qaws_curve const** subj, unsigned int ns,
	qaws_curve const** clip, unsigned int nc, unsigned int flags)
{
	qaws_path_2d sp[8], cp[8];
	qaws_clip_desc d;
	qaws_clip_result* r = NULL;
	unsigned int i;
	qaws_status s;
	for (i = 0; i < ns; i++) { sp[i].curves = &subj[i]; sp[i].curve_count = 1; sp[i].closed = 1; }
	for (i = 0; i < nc; i++) { cp[i].curves = &clip[i]; cp[i].curve_count = 1; cp[i].closed = 1; }
	memset(&d, 0, sizeof(d));
	d.subjects = sp; d.subject_count = ns;
	d.clips = cp; d.clip_count = nc;
	d.clip_type = ct; d.fill_rule = fr; d.flags = flags;
	s = qaws_clip_execute(&d, &r);
	if (s != QAWS_STATUS_OK)
		printf("    clip status %d\n", (int)s);
	return r;
}

static void test_squares(void)
{
	static char const* names[5] = { "none", "intersection", "union", "difference", "xor" };
	static double const areas[5] = { 0, 1, 7, 3, 6 };
	static unsigned int const counts[5] = { 0, 1, 1, 1, 2 };
	qaws_curve* a = clt_rect(0, 0, 2, 2, 1);
	qaws_curve* b = clt_rect(1, 1, 3, 3, 1);
	qaws_curve const* sa = a;
	qaws_curve const* sb = b;
	char msg[256];
	int ct;
	for (ct = 1; ct <= 4; ct++)
	{
		qaws_clip_result* r = clt_run((qaws_clip_type)ct, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, 0);
		double area = r ? clt_area(r) : -1;
		unsigned int n = r ? qaws_clip_result_get_path_count(r) : 99;
		sprintf(msg, "two squares, %s: area %g (got %.12g), %u path(s) (got %u)", names[ct], areas[ct], area, counts[ct], n);
		TEST_ASSERT(r && fabs(area - areas[ct]) < CLT_TOL && n == counts[ct], msg);
		if (r && ct == QAWS_CLIP_UNION)
		{
			qaws_path_2d p;
			qaws_clip_result_get_path(r, 0, &p);
			sprintf(msg, "union: one closed polyline of 8 corners (curves %u, kind %d, vertices %u)", p.curve_count,
				p.curve_count ? (int)qaws_curve_get_kind(p.curves[0]) : -1, clt_vertices(r, 0));
			TEST_ASSERT(p.curve_count == 1 && qaws_curve_get_kind(p.curves[0]) == QAWS_CURVE_KIND_BSPLINE && clt_vertices(r, 0) == 8, msg);
			TEST_ASSERT(qaws_clip_result_get_parent(r, 0) == QAWS_CLIP_NONE_INDEX && !qaws_clip_result_is_hole(r, 0), "union: a top-level outer path");
		}
		qaws_clip_result_destroy(r);
	}
	/* reversed solution: the same areas, negative */
	{
		qaws_clip_result* r = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, QAWS_CLIP_REVERSE_SOLUTION);
		double area = r ? clt_area(r) : 0;
		sprintf(msg, "reverse solution: union area -7 (%.12g)", area);
		TEST_ASSERT(r && fabs(area + 7) < CLT_TOL, msg);
		qaws_clip_result_destroy(r);
	}
	qaws_curve_destroy(a);
	qaws_curve_destroy(b);
}

static void test_degenerate(void)
{
	char msg[256];
	/* sharing an edge */
	{
		qaws_curve* a = clt_rect(0, 0, 1, 1, 1);
		qaws_curve* b = clt_rect(1, 0, 2, 1, 1);
		qaws_curve const* sa = a;
		qaws_curve const* sb = b;
		qaws_clip_result* r = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, 0);
		qaws_clip_result* rk = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, QAWS_CLIP_PRESERVE_COLLINEAR);
		qaws_clip_result* ri = clt_run(QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, 0);
		sprintf(msg, "squares sharing an edge: union one rectangle of area 2 (%u paths, %.12g), 4 corners (%u), 6 with collinear kept (%u); intersection empty (%u)",
			r ? qaws_clip_result_get_path_count(r) : 99, r ? clt_area(r) : 0, r ? clt_vertices(r, 0) : 0,
			rk ? clt_vertices(rk, 0) : 0, ri ? qaws_clip_result_get_path_count(ri) : 99);
		TEST_ASSERT(r && qaws_clip_result_get_path_count(r) == 1 && fabs(clt_area(r) - 2) < CLT_TOL && clt_vertices(r, 0) == 4 &&
			rk && clt_vertices(rk, 0) == 6 && ri && qaws_clip_result_get_path_count(ri) == 0, msg);
		qaws_clip_result_destroy(r);
		qaws_clip_result_destroy(rk);
		qaws_clip_result_destroy(ri);
		qaws_curve_destroy(a);
		qaws_curve_destroy(b);
	}
	/* touching at a corner, and identical */
	{
		qaws_curve* a = clt_rect(0, 0, 1, 1, 1);
		qaws_curve* b = clt_rect(1, 1, 2, 2, 1);
		qaws_curve* c = clt_rect(0, 0, 1, 1, 1);
		qaws_curve const* sa = a;
		qaws_curve const* sb = b;
		qaws_curve const* sc = c;
		qaws_clip_result* r = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, 0);
		qaws_clip_result* u = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &sa, 1, &sc, 1, 0);
		qaws_clip_result* dd = clt_run(QAWS_CLIP_DIFFERENCE, QAWS_FILL_NON_ZERO, &sa, 1, &sc, 1, 0);
		qaws_clip_result* x = clt_run(QAWS_CLIP_XOR, QAWS_FILL_NON_ZERO, &sa, 1, &sc, 1, 0);
		sprintf(msg, "corner touch: two paths (%u); identical squares: union one of area 1 (%u, %.12g), difference and xor empty (%u, %u)",
			r ? qaws_clip_result_get_path_count(r) : 99, u ? qaws_clip_result_get_path_count(u) : 99, u ? clt_area(u) : 0,
			dd ? qaws_clip_result_get_path_count(dd) : 99, x ? qaws_clip_result_get_path_count(x) : 99);
		TEST_ASSERT(r && qaws_clip_result_get_path_count(r) == 2 && u && qaws_clip_result_get_path_count(u) == 1 &&
			fabs(clt_area(u) - 1) < CLT_TOL && dd && qaws_clip_result_get_path_count(dd) == 0 &&
			x && qaws_clip_result_get_path_count(x) == 0, msg);
		qaws_clip_result_destroy(r); qaws_clip_result_destroy(u); qaws_clip_result_destroy(dd); qaws_clip_result_destroy(x);
		qaws_curve_destroy(a); qaws_curve_destroy(b); qaws_curve_destroy(c);
	}
}

static void test_fill_rules(void)
{
	/* pentagram through the 5 points of a unit pentagon taken every second */
	double star[10];
	unsigned int i;
	qaws_curve* s;
	qaws_curve const* ss;
	char msg[256];
	double R = 1, inner, area_star, area_pent;
	for (i = 0; i < 5; i++)
	{
		double a = CLT_PI / 2 + 2 * CLT_PI * (2 * i % 5) / 5;
		star[2 * i] = R * cos(a);
		star[2 * i + 1] = R * sin(a);
	}
	s = clt_poly(star, 5);
	ss = s;
	/* inner pentagon radius r = R cos(72) / cos(36); star = pentagon + 5 tips */
	inner = R * cos(2 * CLT_PI / 5) / cos(CLT_PI / 5);
	area_pent = 2.5 * inner * inner * sin(2 * CLT_PI / 5);
	area_star = 5 * 0.5 * R * inner * sin(CLT_PI / 5) * 2;
	{
		qaws_clip_result* nz = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &ss, 1, NULL, 0, 0);
		qaws_clip_result* eo = clt_run(QAWS_CLIP_UNION, QAWS_FILL_EVEN_ODD, &ss, 1, NULL, 0, 0);
		qaws_clip_result* po = clt_run(QAWS_CLIP_UNION, QAWS_FILL_POSITIVE, &ss, 1, NULL, 0, 0);
		qaws_clip_result* ne = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NEGATIVE, &ss, 1, NULL, 0, 0);
		double anz = nz ? clt_area(nz) : 0, aeo = eo ? clt_area(eo) : 0;
		sprintf(msg, "pentagram: non-zero fills the star (%.10g vs %.10g), even-odd leaves the centre out (%.10g vs %.10g, %u paths); positive %u paths, negative %u",
			anz, area_star, aeo, area_star - area_pent, eo ? qaws_clip_result_get_path_count(eo) : 0,
			po ? qaws_clip_result_get_path_count(po) : 0, ne ? qaws_clip_result_get_path_count(ne) : 0);
		TEST_ASSERT(fabs(anz - area_star) < CLT_TOL && fabs(aeo - (area_star - area_pent)) < CLT_TOL &&
			eo && qaws_clip_result_get_path_count(eo) == 5 && po && qaws_clip_result_get_path_count(po) == 1 &&
			ne && qaws_clip_result_get_path_count(ne) == 0, msg);
		qaws_clip_result_destroy(nz); qaws_clip_result_destroy(eo); qaws_clip_result_destroy(po); qaws_clip_result_destroy(ne);
	}
	qaws_curve_destroy(s);
}

static void test_nesting(void)
{
	qaws_curve* o = clt_rect(0, 0, 10, 10, 1);
	qaws_curve* h = clt_rect(2, 2, 8, 8, 1);
	qaws_curve* i = clt_rect(4, 4, 6, 6, 1);
	qaws_curve const* s[3];
	qaws_clip_result* r;
	char msg[256];
	unsigned int n, k, d0 = 0, d1 = 0, d2 = 0, ok_parent = 1;
	s[0] = o; s[1] = h; s[2] = i;
	r = clt_run(QAWS_CLIP_UNION, QAWS_FILL_EVEN_ODD, s, 3, NULL, 0, 0);
	n = r ? qaws_clip_result_get_path_count(r) : 0;
	for (k = 0; k < n; k++)
	{
		unsigned int d = qaws_clip_result_get_depth(r, k), p = qaws_clip_result_get_parent(r, k);
		if (d == 0) d0++;
		if (d == 1) d1++;
		if (d == 2) d2++;
		if (d > 0 && (p == QAWS_CLIP_NONE_INDEX || qaws_clip_result_get_depth(r, p) != d - 1)) ok_parent = 0;
		if ((d & 1) != (unsigned int)qaws_clip_result_is_hole(r, k)) ok_parent = 0;
	}
	sprintf(msg, "square, hole, island (even-odd): 3 paths at depths 0, 1, 2 with their parents (%u paths: %u %u %u), area 100 - 36 + 4 (%.12g)",
		n, d0, d1, d2, r ? clt_area(r) : 0);
	TEST_ASSERT(n == 3 && d0 == 1 && d1 == 1 && d2 == 1 && ok_parent && fabs(clt_area(r) - 68) < CLT_TOL * 10, msg);
	qaws_clip_result_destroy(r);
	qaws_curve_destroy(o); qaws_curve_destroy(h); qaws_curve_destroy(i);
}

static void test_circles(void)
{
	qaws_curve *a = NULL, *b = NULL, *sq;
	qaws_curve const* sa;
	qaws_curve const* sb;
	char msg[256];
	double r = 1, dd = 1.2, lens, half;
	qaws_curve_create_ellipse_2d(0, 0, 1, 1, 0, &a);
	qaws_curve_create_ellipse_2d((qaws_scalar)dd, 0, 1, 1, 0, &b);
	sa = a; sb = b;
	/* lens of two unit circles at distance d */
	lens = 2 * r * r * acos(dd / (2 * r)) - 0.5 * dd * sqrt(4 * r * r - dd * dd);
	half = 0;
	{
		qaws_clip_result* ri = clt_run(QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, 0);
		qaws_clip_result* ru = clt_run(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, 0);
		qaws_path_2d p;
		int nurbs = 0;
		if (ri && qaws_clip_result_get_path_count(ri) == 1)
		{
			unsigned int k;
			qaws_clip_result_get_path(ri, 0, &p);
			nurbs = p.curve_count == 2 || p.curve_count == 3;   /* 3 when an arc crosses a circle's seam */
			for (k = 0; k < p.curve_count; k++)
				nurbs = nurbs && qaws_curve_get_kind(p.curves[k]) == QAWS_CURVE_KIND_NURBS;
		}
		sprintf(msg, "two circles: lens %.12g (got %.12g) of NURBS pieces (%d); union 2 pi - lens (got %.12g)",
			lens, ri ? clt_area(ri) : 0, nurbs, ru ? clt_area(ru) : 0);
		TEST_ASSERT(ri && fabs(clt_area(ri) - lens) < CLT_TOL && nurbs && ru && fabs(clt_area(ru) - (2 * CLT_PI - lens)) < CLT_TOL, msg);
		qaws_clip_result_destroy(ri);
		qaws_clip_result_destroy(ru);
	}
	/* circle minus a square through its centre: three quarters of the disc */
	sq = clt_rect(0, 0, 2, 2, 1);
	sb = sq;
	{
		qaws_clip_result* rd = clt_run(QAWS_CLIP_DIFFERENCE, QAWS_FILL_NON_ZERO, &sa, 1, &sb, 1, 0);
		sprintf(msg, "circle minus a square at its centre: 3 pi / 4 (got %.12g)", rd ? clt_area(rd) : 0);
		TEST_ASSERT(rd && fabs(clt_area(rd) - 0.75 * CLT_PI) < CLT_TOL, msg);
		qaws_clip_result_destroy(rd);
	}
	(void)half;
	qaws_curve_destroy(a); qaws_curve_destroy(b); qaws_curve_destroy(sq);
}

static qaws_scalar clt_z(void* user, qaws_clip_vertex const* v)
{
	(*(unsigned int*)user)++;
	return v->position.x + v->position.y;
}

static void test_open_and_z(void)
{
	qaws_curve* sq = clt_rect(0, 0, 2, 2, 1);
	qaws_curve* other = clt_rect(1, 1, 3, 3, 1);
	qaws_scalar lp[4] = { -1, 1, 3, 1 };
	qaws_curve* line = NULL;
	qaws_curve const* cs = sq;
	qaws_curve const* co = other;
	qaws_curve const* cl;
	qaws_path_2d clip, open, subj;
	qaws_clip_desc d;
	qaws_clip_result* r = NULL;
	unsigned int zc = 0;
	char msg[256];
	qaws_curve_create_polyline_2d(lp, 2, 0, &line);
	cl = line;
	clip.curves = &cs; clip.curve_count = 1; clip.closed = 1;
	open.curves = &cl; open.curve_count = 1; open.closed = 0;
	memset(&d, 0, sizeof(d));
	d.open_subjects = &open; d.open_subject_count = 1;
	d.clips = &clip; d.clip_count = 1;
	d.clip_type = QAWS_CLIP_INTERSECTION; d.fill_rule = QAWS_FILL_NON_ZERO;
	qaws_clip_execute(&d, &r);
	{
		unsigned int ni = r ? qaws_clip_result_get_open_path_count(r) : 0, no;
		double len = 0;
		qaws_clip_result* r2 = NULL;
		if (ni == 1)
		{
			qaws_path_2d p;
			qaws_scalar l = 0;
			qaws_clip_result_get_open_path(r, 0, &p);
			qaws_path_compute_length_2d(&p, &l);
			len = l;
		}
		d.clip_type = QAWS_CLIP_DIFFERENCE;
		qaws_clip_execute(&d, &r2);
		no = r2 ? qaws_clip_result_get_open_path_count(r2) : 0;
		sprintf(msg, "a line through a square: intersection keeps the inside piece (%u, length %.12g), difference the two outside (%u)", ni, len, no);
		TEST_ASSERT(ni == 1 && fabs(len - 2) < CLT_TOL && no == 2, msg);
		qaws_clip_result_destroy(r2);
	}
	qaws_clip_result_destroy(r);
	/* z callback: two squares crossing twice */
	subj.curves = &co; subj.curve_count = 1; subj.closed = 1;
	memset(&d, 0, sizeof(d));
	d.subjects = &subj; d.subject_count = 1;
	d.clips = &clip; d.clip_count = 1;
	d.clip_type = QAWS_CLIP_UNION; d.fill_rule = QAWS_FILL_NON_ZERO;
	d.z_fn = clt_z; d.z_user = &zc;
	r = NULL;
	qaws_clip_execute(&d, &r);
	{
		qaws_clip_vertex const* v = NULL;
		unsigned int n = 0, k, zok = 0;
		if (r) qaws_clip_result_get_vertices(r, 0, 0, &v, &n);
		for (k = 0; k < n; k++)
			if (v[k].operand_b != QAWS_CLIP_NONE_INDEX && fabs(v[k].z - (v[k].position.x + v[k].position.y)) < 1e-12)
				zok++;
		sprintf(msg, "z callback: called at the 2 crossings (%u calls), their output vertices carry z (%u)", zc, zok);
		TEST_ASSERT(zc == 2 && zok == 2, msg);
	}
	qaws_clip_result_destroy(r);
	qaws_curve_destroy(sq); qaws_curve_destroy(other); qaws_curve_destroy(line);
}

int test_84_clip_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 84: Boolean operations on regions of curves\n");
	test_squares();
	test_degenerate();
	test_fill_rules();
	test_nesting();
	test_circles();
	test_open_and_z();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
