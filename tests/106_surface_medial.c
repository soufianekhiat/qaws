/*
 * Test 106: medial axes of regions on surfaces
 * (qaws_surface_region_medial_axis)
 *
 *   - a sphere: a band between two latitudes has the middle latitude for
 *     its axis, every point at the half-width; a geodesic triangle's
 *     junction and largest inscribed circle have the inradius,
 *     tan r = sqrt(sin(s - a) sin(s - b) sin(s - c) / sin s)
 *   - a cylinder (developable): a rectangle's axis, unrolled, lies on the
 *     plane rectangle's (a middle segment and four bisectors), each point at
 *     its distance to the nearest side
 *   - a wavy bicubic patch, by the geodesic offsets: an axis point at radius
 *     r is inside the region shrunk by r - d and outside it shrunk by r + d;
 *     shrunk by the largest inscribed radius less d a part is left, by it
 *     plus d nothing
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include "qaws_bezier.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SMT_PI 3.14159265358979323846
#define SMT_F32 (sizeof(qaws_scalar) == 4)

static qaws_vec2 smt_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static qaws_surface_curve* smt_line(qaws_surface const* s, double a, double b, double c, double d, qaws_curve** k)
{
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_surface_curve* r = NULL;
	cp[0] = (qaws_scalar)a; cp[1] = (qaws_scalar)b; cp[2] = (qaws_scalar)c; cp[3] = (qaws_scalar)d;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	*k = NULL;
	qaws_curve_create_bezier(k, &bd);
	qaws_surface_curve_create_uv(s, *k, &r);
	return r;
}

static qaws_surface* smt_sphere(double R, qaws_curve** prof)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)R; g.angle_start = (qaws_scalar)(-SMT_PI / 2); g.angle_end = (qaws_scalar)(SMT_PI / 2);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(prof, &ad);
	memset(&rd, 0, sizeof(rd));
	rd.profile = *prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rd, &s);
	return s;
}

static void test_sphere(void)
{
	double R = 1.3, tol = SMT_F32 ? 2e-4 : 1e-9;
	qaws_curve* prof = NULL;
	qaws_surface* s = smt_sphere(R, &prof);
	qaws_surface_medial_desc md;
	char msg[300];
	memset(&md, 0, sizeof(md));
	md.resolution = (qaws_scalar)0.06;
	/* the band v in [0.4, 0.6]: the equator, at 0.1 pi R */
	{
		qaws_curve *k1, *k2;
		qaws_surface_curve *a = smt_line(s, 0, 0.4, 1, 0.4, &k1), *b = smt_line(s, 1, 0.6, 0, 0.6, &k2);
		qaws_surface_curve const *pa[1], *pb[1];
		qaws_surface_path B[2];
		qaws_surface_region* r = NULL;
		qaws_surface_medial_axis* m = NULL;
		qaws_status st;
		double worst = 0;
		unsigned int i, n = 0;
		pa[0] = a; pb[0] = b;
		B[0].curves = pa; B[0].curve_count = 1; B[0].closed = 1;
		B[1].curves = pb; B[1].curve_count = 1; B[1].closed = 1;
		qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, B, 2, NULL, 0, &r);
		st = qaws_surface_region_medial_axis(r, &md, &m);
		if (st == QAWS_STATUS_OK)
			for (i = 0; i < (n = qaws_surface_medial_axis_get_segment_count(m)); i++)
			{
				qaws_vec2 p[2];
				qaws_scalar rr[2];
				unsigned int k;
				qaws_surface_medial_axis_get_segment(m, i, p, rr);
				for (k = 0; k < 2; k++)
				{
					if (fabs(p[k].y - 0.5) * SMT_PI * R > worst) worst = fabs(p[k].y - 0.5) * SMT_PI * R;
					if (fabs(rr[k] - 0.1 * SMT_PI * R) > worst) worst = fabs(rr[k] - 0.1 * SMT_PI * R);
				}
			}
		sprintf(msg, "sphere: a band between latitudes, its axis (%u segments) on the equator at the half-width within %.2g", n, worst);
		TEST_ASSERT(st == QAWS_STATUS_OK && n > 50 && worst <= tol, msg);
		qaws_surface_medial_axis_destroy(m);
		qaws_surface_region_destroy(r);
		qaws_surface_curve_destroy(a); qaws_surface_curve_destroy(b);
		qaws_curve_destroy(k1); qaws_curve_destroy(k2);
	}
	/* a geodesic triangle: the inradius */
	{
		qaws_vec2 v[3];
		qaws_surface_curve* t[3] = { 0 };
		qaws_surface_curve const* pt[3];
		qaws_surface_path T;
		qaws_surface_region* r = NULL;
		qaws_surface_medial_axis* m = NULL;
		qaws_scalar la[3], jr = 0, br = 0;
		qaws_vec2 jp, bp;
		double ss, inr, a[3];
		unsigned int i, nj = 0;
		qaws_status st;
		v[0] = smt_uv(0.3, 0.4); v[1] = smt_uv(0.45, 0.42); v[2] = smt_uv(0.37, 0.65);
		for (i = 0; i < 3; i++)
		{
			qaws_surface_curve_create_geodesic(s, v[i], v[(i + 1) % 3], 0, &t[i]);
			pt[i] = t[i];
			qaws_surface_curve_get_length(t[i], &la[i]);
			a[i] = (double)la[i] / R;
		}
		T.curves = pt; T.curve_count = 3; T.closed = 1;
		qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &T, 1, NULL, 0, &r);
		ss = 0.5 * (a[0] + a[1] + a[2]);
		inr = R * atan(sqrt(sin(ss - a[0]) * sin(ss - a[1]) * sin(ss - a[2]) / sin(ss)));
		st = qaws_surface_region_medial_axis(r, &md, &m);
		if (st == QAWS_STATUS_OK)
		{
			nj = qaws_surface_medial_axis_get_junction_count(m);
			if (nj) qaws_surface_medial_axis_get_junction(m, 0, &jp, &jr);
			qaws_surface_medial_axis_get_max_inscribed(m, &bp, &br);
		}
		sprintf(msg, "sphere: a geodesic triangle, %u junction at the inradius %.2g off, the largest inscribed circle %.2g off",
			nj, fabs((double)jr - inr), fabs((double)br - inr));
		TEST_ASSERT(st == QAWS_STATUS_OK && nj == 1 && fabs((double)jr - inr) <= tol && fabs((double)br - inr) <= tol, msg);
		qaws_surface_medial_axis_destroy(m);
		qaws_surface_region_destroy(r);
		for (i = 0; i < 3; i++) qaws_surface_curve_destroy(t[i]);
	}
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

/* the distance from (x, y) to the plane rectangle [0, W] x [0, H]'s medial axis */
static double smt_rect_axis(double x, double y, double W, double H)
{
	double h = 0.5 * H, best = 1e300, seg[5][4];
	unsigned int i;
	seg[0][0] = h; seg[0][1] = h; seg[0][2] = W - h; seg[0][3] = h;
	seg[1][0] = 0; seg[1][1] = 0; seg[1][2] = h; seg[1][3] = h;
	seg[2][0] = 0; seg[2][1] = H; seg[2][2] = h; seg[2][3] = h;
	seg[3][0] = W; seg[3][1] = 0; seg[3][2] = W - h; seg[3][3] = h;
	seg[4][0] = W; seg[4][1] = H; seg[4][2] = W - h; seg[4][3] = h;
	for (i = 0; i < 5; i++)
	{
		double dx = seg[i][2] - seg[i][0], dy = seg[i][3] - seg[i][1], t = ((x - seg[i][0]) * dx + (y - seg[i][1]) * dy) / (dx * dx + dy * dy), d;
		if (t < 0) t = 0;
		if (t > 1) t = 1;
		d = hypot(x - seg[i][0] - t * dx, y - seg[i][1] - t * dy);
		if (d < best) best = d;
	}
	return best;
}

static void test_cylinder(void)
{
	double Rc = 0.8, Hc = 2.0, tol = SMT_F32 ? 2e-4 : 1e-9, worst = 0, P;
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* circle = NULL;
	qaws_surface_extrusion_desc ed;
	qaws_surface* cy = NULL;
	qaws_curve* k[4];
	qaws_surface_curve* e[4];
	qaws_surface_curve const* pe[4];
	qaws_surface_path Q;
	qaws_surface_region* r = NULL;
	qaws_surface_medial_axis* m = NULL;
	qaws_surface_medial_desc md;
	qaws_range ru;
	double u0, u1, v0 = 0.25, v1 = 0.55, W, H;
	unsigned int i, n = 0;
	qaws_status st;
	char msg[300];
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)Rc; g.angle_end = (qaws_scalar)(2 * SMT_PI);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&circle, &ad);
	memset(&ed, 0, sizeof(ed));
	ed.profile = circle; ed.direction.z = 1; ed.length = (qaws_scalar)Hc;
	qaws_surface_create_extrusion(&ed, &cy);
	ru = qaws_surface_get_u_range(cy);
	P = (double)ru.max_value - ru.min_value;
	u0 = ru.min_value + 0.2 * P; u1 = ru.min_value + 0.55 * P;
	W = (u1 - u0) / P * 2 * SMT_PI * Rc;
	H = (v1 - v0) * Hc;
	e[0] = smt_line(cy, u0, v0, u1, v0, &k[0]); e[1] = smt_line(cy, u1, v0, u1, v1, &k[1]);
	e[2] = smt_line(cy, u1, v1, u0, v1, &k[2]); e[3] = smt_line(cy, u0, v1, u0, v0, &k[3]);
	for (i = 0; i < 4; i++) pe[i] = e[i];
	Q.curves = pe; Q.curve_count = 4; Q.closed = 1;
	qaws_surface_clip(cy, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &Q, 1, NULL, 0, &r);
	memset(&md, 0, sizeof(md));
	md.resolution = (qaws_scalar)0.05;
	st = qaws_surface_region_medial_axis(r, &md, &m);
	if (st == QAWS_STATUS_OK)
		for (i = 0; i < (n = qaws_surface_medial_axis_get_segment_count(m)); i++)
		{
			qaws_vec2 p[2];
			qaws_scalar rr[2];
			unsigned int q;
			qaws_surface_medial_axis_get_segment(m, i, p, rr);
			for (q = 0; q < 2; q++)
			{
				double x = ((double)p[q].x - u0) / P * 2 * SMT_PI * Rc, y = ((double)p[q].y - v0) * Hc;
				double dside = fmin(fmin(x, W - x), fmin(y, H - y));
				if (smt_rect_axis(x, y, W, H) > worst) worst = smt_rect_axis(x, y, W, H);
				if (fabs((double)rr[q] - dside) > worst) worst = fabs((double)rr[q] - dside);
			}
		}
	sprintf(msg, "cylinder: a rectangle's axis, unrolled (%u segments), on the plane rectangle's within %.2g, radii the distances to the sides", n, worst);
	TEST_ASSERT(st == QAWS_STATUS_OK && n > 20 && worst <= tol, msg);
	qaws_surface_medial_axis_destroy(m);
	qaws_surface_region_destroy(r);
	for (i = 0; i < 4; i++) { qaws_surface_curve_destroy(e[i]); qaws_curve_destroy(k[i]); }
	qaws_surface_destroy(cy);
	qaws_curve_destroy(circle);
}

/* the plane L [0, 3] x [0, 1] u [0, 1] x [0, 3]: its straight skeleton's
   segments, with the times at their ends (linear along each) */
static double const smt_lsk[8][6] = {
	{ 0, 0, 0, 0.5, 0.5, 0.5 }, { 1, 1, 0, 0.5, 0.5, 0.5 }, { 0.5, 0.5, 0.5, 2.5, 0.5, 0.5 }, { 3, 0, 0, 2.5, 0.5, 0.5 },
	{ 3, 1, 0, 2.5, 0.5, 0.5 }, { 0.5, 0.5, 0.5, 0.5, 2.5, 0.5 }, { 0, 3, 0, 0.5, 2.5, 0.5 }, { 1, 3, 0, 0.5, 2.5, 0.5 } };

/* (x, y)'s distance to the L's skeleton, and the time there */
static double smt_l_skeleton(double x, double y, double* time)
{
	double best = 1e300;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		double const* g = smt_lsk[i];
		double dx = g[3] - g[0], dy = g[4] - g[1], l = ((x - g[0]) * dx + (y - g[1]) * dy) / (dx * dx + dy * dy), d;
		l = l < 0 ? 0 : l > 1 ? 1 : l;
		d = hypot(x - g[0] - l * dx, y - g[1] - l * dy);
		if (d < best) { best = d; *time = g[2] + l * (g[5] - g[2]); }
	}
	return best;
}

/* the mitered skeleton of an L on a cylinder (developable): unrolled, on
   the plane L's straight skeleton, times the distances swept; the mitre
   from the reflex corner (not on the medial axis) found */
static void test_mitered(void)
{
	double Rc = 0.8, Hc = 2.0, w = 0.4, tol = SMT_F32 ? 2e-4 : 1e-8, worst = 0, cworst = 0, P, u0, v0 = 0.2, best = 1e300;
	static double const L[6][2] = { { 0, 0 }, { 3, 0 }, { 3, 1 }, { 1, 1 }, { 1, 3 }, { 0, 3 } };
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* circle = NULL;
	qaws_surface_extrusion_desc ed;
	qaws_surface* cy = NULL;
	qaws_curve* k[6];
	qaws_surface_curve* e[6];
	qaws_surface_curve const* pe[6];
	qaws_surface_path Q;
	qaws_surface_region* r = NULL;
	qaws_surface_medial_axis* m = NULL;
	qaws_surface_medial_desc md;
	qaws_range ru;
	unsigned int i, n = 0;
	qaws_status st;
	char msg[300];
	memset(&g, 0, sizeof(g));
	g.radius = (qaws_scalar)Rc; g.angle_end = (qaws_scalar)(2 * SMT_PI);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&circle, &ad);
	memset(&ed, 0, sizeof(ed));
	ed.profile = circle; ed.direction.z = 1; ed.length = (qaws_scalar)Hc;
	qaws_surface_create_extrusion(&ed, &cy);
	ru = qaws_surface_get_u_range(cy);
	P = (double)ru.max_value - ru.min_value;
	u0 = ru.min_value + 0.2 * P;
	/* (x, y) unrolled in units of w: u = u0 + x w P / (2 pi Rc), v = v0 + y w / Hc */
	for (i = 0; i < 6; i++)
	{
		double const* a = L[i];
		double const* b = L[(i + 1) % 6];
		e[i] = smt_line(cy, u0 + a[0] * w * P / (2 * SMT_PI * Rc), v0 + a[1] * w / Hc, u0 + b[0] * w * P / (2 * SMT_PI * Rc), v0 + b[1] * w / Hc, &k[i]);
		pe[i] = e[i];
	}
	Q.curves = pe; Q.curve_count = 6; Q.closed = 1;
	qaws_surface_clip(cy, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &Q, 1, NULL, 0, &r);
	memset(&md, 0, sizeof(md));
	md.resolution = (qaws_scalar)0.04;
	st = qaws_surface_region_mitered_skeleton(r, &md, &m);
	if (st == QAWS_STATUS_OK)
		for (i = 0; i < (n = qaws_surface_medial_axis_get_segment_count(m)); i++)
		{
			qaws_vec2 p[2];
			qaws_scalar rr[2];
			unsigned int q;
			qaws_surface_medial_axis_get_segment(m, i, p, rr);
			for (q = 0; q < 2; q++)
			{
				double x = ((double)p[q].x - u0) / P * 2 * SMT_PI * Rc / w, y = ((double)p[q].y - v0) * Hc / w, t = 0, d = smt_l_skeleton(x, y, &t);
				/* (an end at a corner, radius about 0: the feet so near a line settle to about 1e-6) */
				if ((double)rr[q] / w < 1e-4) cworst = fmax(cworst, d); else worst = fmax(worst, fmax(d, fabs((double)rr[q] / w - t)));
				best = fmin(best, hypot(x - 0.75, y - 0.75));
			}
		}
	sprintf(msg, "mitered skeleton: an L on a cylinder, unrolled (%u segments, status %d), on the plane L's straight skeleton within %.2g (times too), its reflex mitre met %.2g from its middle, ends at corners within %.2g",
		n, (int)st, worst, best, cworst);
	TEST_ASSERT(st == QAWS_STATUS_OK && n > 20 && worst <= tol && best <= 0.2 && cworst <= (SMT_F32 ? 1e-3 : 1e-5), msg);
	qaws_surface_medial_axis_destroy(m);
	qaws_surface_region_destroy(r);
	for (i = 0; i < 6; i++) { qaws_surface_curve_destroy(e[i]); qaws_curve_destroy(k[i]); }
	qaws_surface_destroy(cy);
	qaws_curve_destroy(circle);
}

static void test_wavy(void)
{
	static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
	qaws_vec3 cp[16];
	qaws_surface_bezier_desc bd;
	qaws_surface* s = NULL;
	qaws_vec2 corner[4];
	qaws_surface_curve* q[4] = { 0 };
	qaws_surface_curve const* pq[4];
	qaws_surface_path Q;
	qaws_surface_region* r = NULL;
	qaws_surface_medial_axis* m = NULL;
	qaws_surface_medial_desc md;
	double dlt = SMT_F32 ? 2e-3 : 1e-6;
	unsigned int i, n = 0, checked = 0, good = 0;
	qaws_status st;
	char msg[300];
	for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
	memset(&bd, 0, sizeof(bd));
	bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
	qaws_surface_create_bezier(&bd, &s);
	corner[0] = smt_uv(0.15, 0.2); corner[1] = smt_uv(0.85, 0.25); corner[2] = smt_uv(0.75, 0.7); corner[3] = smt_uv(0.2, 0.8);
	for (i = 0; i < 4; i++) { qaws_surface_curve_create_geodesic(s, corner[i], corner[(i + 1) % 4], 0, &q[i]); pq[i] = q[i]; }
	Q.curves = pq; Q.curve_count = 4; Q.closed = 1;
	qaws_surface_clip(s, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &Q, 1, NULL, 0, &r);
	memset(&md, 0, sizeof(md));
	md.resolution = (qaws_scalar)0.06;
	st = qaws_surface_region_medial_axis(r, &md, &m);
	if (st == QAWS_STATUS_OK)
	{
		n = qaws_surface_medial_axis_get_segment_count(m);
		/* a few axis points, each against the region shrunk to either side of
		   its radius */
		for (i = 0; i < n && checked < 5; i += 7)
		{
			qaws_vec2 p[2];
			qaws_scalar rr[2];
			qaws_surface_region *in = NULL, *out = NULL;
			int c1 = 0, c2 = 1;
			qaws_status s1, s2;
			qaws_surface_medial_axis_get_segment(m, i, p, rr);
			if (rr[0] < 1e-3)
				continue;      /* (an end at a corner, on the boundary) */
			s1 = qaws_surface_offset(s, QAWS_FILL_NON_ZERO, &Q, 1, -(rr[0] - (qaws_scalar)dlt), 0, &in);
			s2 = qaws_surface_offset(s, QAWS_FILL_NON_ZERO, &Q, 1, -(rr[0] + (qaws_scalar)dlt), 0, &out);
			if (s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK)
			{
				qaws_surface_region_contains(in, p[0], &c1);
				if (out) qaws_surface_region_contains(out, p[0], &c2);
				checked++;
				if (c1 == 1 && c2 == 0) good++;
			}
			qaws_surface_region_destroy(in); qaws_surface_region_destroy(out);
		}
	}
	sprintf(msg, "wavy patch: a geodesic quadrilateral's axis (%u segments): %u of %u points inside the region shrunk by r - %.0e and outside it shrunk by r + %.0e",
		n, good, checked, dlt, dlt);
	TEST_ASSERT(st == QAWS_STATUS_OK && checked == 5 && good == checked, msg);
	/* the largest inscribed circle: the region shrunk by just less keeps a
	   part, shrunk by just more vanishes */
	{
		qaws_vec2 bc;
		qaws_scalar br = 0, a1 = 0, a2 = 0;
		qaws_surface_region *in = NULL, *out = NULL;
		qaws_status s1 = QAWS_STATUS_INVALID_ARGUMENT, s2 = QAWS_STATUS_INVALID_ARGUMENT;
		if (st == QAWS_STATUS_OK)
		{
			qaws_surface_medial_axis_get_max_inscribed(m, &bc, &br);
			s1 = qaws_surface_offset(s, QAWS_FILL_NON_ZERO, &Q, 1, -(br - (qaws_scalar)dlt), 0, &in);
			s2 = qaws_surface_offset(s, QAWS_FILL_NON_ZERO, &Q, 1, -(br + (qaws_scalar)dlt), 0, &out);
			if (in) qaws_surface_region_get_area(in, &a1);
			if (out) qaws_surface_region_get_area(out, &a2);
		}
		sprintf(msg, "wavy patch: the region shrunk by the largest inscribed radius %.6f less %.0e keeps %.2g, by it plus %.0e keeps %.2g",
			(double)br, dlt, (double)a1, dlt, (double)a2);
		TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && a1 > 0 && a2 == 0, msg);
		qaws_surface_region_destroy(in); qaws_surface_region_destroy(out);
	}
	qaws_surface_medial_axis_destroy(m);
	qaws_surface_region_destroy(r);
	for (i = 0; i < 4; i++) qaws_surface_curve_destroy(q[i]);
	qaws_surface_destroy(s);
}

int test_106_surface_medial_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 106: medial axes of regions on surfaces\n");
	test_sphere();
	test_cylinder();
	test_wavy();
	test_mitered();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
