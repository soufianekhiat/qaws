/*
 * Test 100: mean value coordinates over regions bounded by curves
 * (qaws_mean_value)
 *
 *   - against Mathematica (tests/reference/100_mean_value.wls): a blob of
 *     four cubic Beziers, its twelve handles' coordinates at four points
 *     (one near a corner, one near an edge); the transfinite interpolant
 *     of a nonlinear function over a NURBS ellipse and over a disc with a
 *     square hole
 *   - polygons: an L with a hole against Floater's (and Hormann's) closed
 *     form, computed here from the angles; the same polygon with its edges
 *     as straight quadratic Beziers (integrated, not in closed form) gives
 *     the same coordinates; the interpolant of piecewise linear boundary
 *     data equals the coordinates' combination
 *   - every kind (lines, arcs, Beziers, rational Beziers, B-splines,
 *     NURBS, a Hermite fitted): the coordinates sum to 1 and reproduce the
 *     point, at points sampled over each region
 *   - on the boundary the coordinates are its basis functions, and they
 *     tend to them from inside
 *   - a loop's orientation changes nothing
 *   - deformations: the handles unmoved, translated or under an affine map
 *     move the points the same way; an executor gives the serial points
 */

#include "test_common.h"
#include "qaws_mean_value.h"
#include "qaws_region_sample.h"
#include "reference/100_mean_value.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MVT_PI 3.14159265358979323846
#define MVT_F32 (sizeof(qaws_scalar) == 4)

static qaws_curve* mvt_bezier(double const* xy, unsigned int n)
{
	qaws_scalar p[32];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = n - 1; d.control_points = p; d.control_point_count = n;
	qaws_curve_create_bezier(&c, &d);
	return c;
}

static qaws_curve* mvt_circle(double cx, double cy, double r, int ccw)
{
	qaws_arc_segment g;
	qaws_arc_desc d;
	qaws_curve* c = NULL;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)cx; g.center[1] = (qaws_scalar)cy; g.radius = (qaws_scalar)r;
	g.angle_start = ccw ? 0 : (qaws_scalar)(2 * MVT_PI);
	g.angle_end = ccw ? (qaws_scalar)(2 * MVT_PI) : 0;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
	qaws_curve_create_arc(&c, &d);
	return c;
}

static qaws_curve* mvt_polyline(double const* xy, unsigned int n)
{
	qaws_scalar p[64];
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	qaws_curve_create_polyline_2d(&c, p, n, 1);
	return c;
}

static double const mvt_blob[4][8] = {
	{ 0, 0, 1, -0.9, 2.4, -0.5, 2.2, 0.6 }, { 2.2, 0.6, 2, 1.6, 1.4, 1.1, 1, 1 },
	{ 1, 1, 0.6, 0.9, 0.2, 1.8, -0.5, 1.2 }, { -0.5, 1.2, -1.2, 0.6, -0.6, 0.2, 0, 0 } };

static qaws_status mvt_make(qaws_path_2d* p, unsigned int np, qaws_mean_value** out)
{
	qaws_mean_value_desc d;
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = np;
	return qaws_mean_value_create(&d, out);
}

/* the coordinates' sum less 1 and the reproduction error at x */
static qaws_status mvt_check(qaws_mean_value const* m, double x, double y, double* sum_err, double* rep_err)
{
	unsigned int n = qaws_mean_value_get_handle_count(m), i;
	qaws_scalar* w = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (n + 1));
	qaws_vec2 pt;
	double s = 0, sx = 0, sy = 0;
	qaws_status st;
	pt.x = (qaws_scalar)x; pt.y = (qaws_scalar)y;
	st = qaws_mean_value_coordinates(m, pt, w);
	for (i = 0; i < n; i++)
	{
		qaws_vec2 h;
		qaws_mean_value_get_handle(m, i, &h);
		s += w[i];
		sx += w[i] * h.x;
		sy += w[i] * h.y;
	}
	*sum_err = fabs(s - 1);
	*rep_err = hypot(sx - x, sy - y);
	free(w);
	return st;
}

/* ------------------------------------------------------------------------ */

static void mvt_ell_fn(void* user, unsigned int path, unsigned int curve, qaws_scalar t, qaws_vec2 p, qaws_scalar* out)
{
	double x = p.x, y = p.y;
	(void)user; (void)path; (void)curve; (void)t;
	out[0] = (qaws_scalar)(x * x - y + x * y * y * y);
}

static void mvt_disc_fn(void* user, unsigned int path, unsigned int curve, qaws_scalar t, qaws_vec2 p, qaws_scalar* out)
{
	(void)user; (void)path; (void)curve; (void)t;
	out[0] = (qaws_scalar)(sin(2 * (double)p.x) + (double)p.y * p.y);
}

static void test_reference(void)
{
	double tol = MVT_F32 ? 2e-4 : 1e-11;
	char msg[300];
	/* the blob's handles' coordinates */
	{
		qaws_curve* cs[4];
		qaws_curve const* cv[4];
		qaws_path_2d p;
		qaws_mean_value* m = NULL;
		unsigned int i, k, n;
		double worst = 0;
		qaws_status st;
		for (i = 0; i < 4; i++) { cs[i] = mvt_bezier(mvt_blob[i], 4); cv[i] = cs[i]; }
		p.curves = cv; p.curve_count = 4; p.closed = 1;
		st = mvt_make(&p, 1, &m);
		n = m ? qaws_mean_value_get_handle_count(m) : 0;
		for (k = 0; k < 4 && n == 12; k++)
		{
			qaws_scalar w[12];
			qaws_vec2 pt;
			pt.x = (qaws_scalar)ref_blob_points[2 * k]; pt.y = (qaws_scalar)ref_blob_points[2 * k + 1];
			if (qaws_mean_value_coordinates(m, pt, w) != QAWS_STATUS_OK)
				worst = 1e300;
			for (i = 0; i < 12; i++)
				if (fabs((double)w[i] - ref_blob_coords[12 * k + i]) > worst)
					worst = fabs((double)w[i] - ref_blob_coords[12 * k + i]);
		}
		sprintf(msg, "blob of four cubics: %u handles (12), the coordinates at four points against Mathematica within %.2g", n, worst);
		TEST_ASSERT(st == QAWS_STATUS_OK && n == 12 && worst <= tol, msg);
		qaws_mean_value_destroy(m);
		for (i = 0; i < 4; i++) qaws_curve_destroy(cs[i]);
	}
	/* the ellipse's interpolant */
	{
		qaws_curve* c = NULL;
		qaws_curve const* cv;
		qaws_path_2d p;
		qaws_mean_value* m = NULL;
		unsigned int k;
		double worst = 0;
		qaws_curve_create_ellipse_2d(&c, (qaws_scalar)0.1, (qaws_scalar)-0.2, (qaws_scalar)1.5, (qaws_scalar)0.8, (qaws_scalar)0.3);
		cv = c;
		p.curves = &cv; p.curve_count = 1; p.closed = 1;
		mvt_make(&p, 1, &m);
		for (k = 0; k < 3 && m; k++)
		{
			qaws_vec2 pt;
			qaws_scalar v = 0;
			pt.x = (qaws_scalar)ref_ellipse_points[2 * k]; pt.y = (qaws_scalar)ref_ellipse_points[2 * k + 1];
			if (qaws_mean_value_interpolate(m, pt, mvt_ell_fn, NULL, 1, &v) != QAWS_STATUS_OK)
				worst = 1e300;
			if (fabs((double)v - ref_ellipse_values[k]) > worst)
				worst = fabs((double)v - ref_ellipse_values[k]);
		}
		sprintf(msg, "NURBS ellipse: the interpolant of x^2 - y + x y^3 at three points against Mathematica within %.2g", worst);
		TEST_ASSERT(m && worst <= tol, msg);
		qaws_mean_value_destroy(m);
		qaws_curve_destroy(c);
	}
	/* the disc with a square hole */
	{
		double sq[8] = { -0.2, -0.2, -0.2, 0.2, 0.2, 0.2, 0.2, -0.2 };
		qaws_curve* cs[2];
		qaws_curve const* cv[2];
		qaws_path_2d p[2];
		qaws_mean_value* m = NULL;
		unsigned int k;
		double worst = 0;
		cs[0] = mvt_circle(0, 0, 1, 1); cs[1] = mvt_polyline(sq, 4);
		cv[0] = cs[0]; cv[1] = cs[1];
		p[0].curves = cv; p[0].curve_count = 1; p[0].closed = 1;
		p[1].curves = cv + 1; p[1].curve_count = 1; p[1].closed = 1;
		mvt_make(p, 2, &m);
		for (k = 0; k < 3 && m; k++)
		{
			qaws_vec2 pt;
			qaws_scalar v = 0;
			pt.x = (qaws_scalar)ref_disc_points[2 * k]; pt.y = (qaws_scalar)ref_disc_points[2 * k + 1];
			if (qaws_mean_value_interpolate(m, pt, mvt_disc_fn, NULL, 1, &v) != QAWS_STATUS_OK)
				worst = 1e300;
			if (fabs((double)v - ref_disc_values[k]) > worst)
				worst = fabs((double)v - ref_disc_values[k]);
		}
		sprintf(msg, "disc with a square hole: the interpolant of sin 2x + y^2 at three points against Mathematica within %.2g", worst);
		TEST_ASSERT(m && worst <= tol, msg);
		qaws_mean_value_destroy(m);
		qaws_curve_destroy(cs[0]); qaws_curve_destroy(cs[1]);
	}
}

/* ------------------------------------------------------------------------ */

/* Floater's weights of a loop at x, the loop turned to the region's side
   (outside counter-clockwise, holes clockwise): added to w */
static void mvt_floater(double const* v, unsigned int n, double x, double y, double* w)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		unsigned int a = (i + n - 1) % n, b = (i + 1) % n;
		double dix = v[2 * i] - x, diy = v[2 * i + 1] - y, r = hypot(dix, diy);
		double dax = v[2 * a] - x, day = v[2 * a + 1] - y, dbx = v[2 * b] - x, dby = v[2 * b + 1] - y;
		double t1 = (dax * diy - day * dix) / (hypot(dax, day) * r + dax * dix + day * diy);
		double t2 = (dix * dby - diy * dbx) / (r * hypot(dbx, dby) + dix * dbx + diy * dby);
		w[i] += (t1 + t2) / r;
	}
}

static double mvt_hv[16];
static void mvt_pl_fn(void* user, unsigned int path, unsigned int curve, qaws_scalar t, qaws_vec2 p, qaws_scalar* out)
{
	/* outer loop: 6 vertices (handles 0-5), hole: 4 (handles 6-9); a polyline's parameter counts its edges */
	unsigned int n = path ? 4 : 6, base = path ? 6 : 0, k = (unsigned int)floor((double)t);
	double u;
	(void)user; (void)curve; (void)p;
	if (k >= n) k = n - 1;
	u = (double)t - k;
	out[0] = (qaws_scalar)((1 - u) * mvt_hv[base + k] + u * mvt_hv[base + (k + 1) % n]);
}

static void test_polygons(void)
{
	/* an L counter-clockwise, a square hole counter-clockwise too (turned by its nesting) */
	double L[12] = { 0, 0, 2, 0, 2, 1, 1, 1, 1, 2, 0, 2 };
	double H[8] = { 0.3, 0.3, 0.7, 0.3, 0.7, 0.7, 0.3, 0.7 };
	double Hcw[8] = { 0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.3 };
	double pts[10] = { 0.5, 0.15, 1.5, 0.5, 0.5, 1.5, 0.85, 0.85, 0.15, 0.9 };
	double tol = MVT_F32 ? 2e-5 : 1e-13, worst_cf = 0, worst_quad = 0, worst_int = 0;
	qaws_curve *cs[2], *qs[10];
	qaws_curve const *cv[2], *qv[10];
	qaws_path_2d p[2];
	qaws_mean_value *m = NULL, *mq = NULL;
	unsigned int k, i;
	char msg[400];
	cs[0] = mvt_polyline(L, 6); cs[1] = mvt_polyline(H, 4);
	cv[0] = cs[0]; cv[1] = cs[1];
	p[0].curves = cv; p[0].curve_count = 1; p[0].closed = 1;
	p[1].curves = cv + 1; p[1].curve_count = 1; p[1].closed = 1;
	mvt_make(p, 2, &m);
	/* the same edges as quadratic Beziers, their middle control points halfway */
	for (i = 0; i < 10; i++)
	{
		double const* v = i < 6 ? L : H;
		unsigned int n = i < 6 ? 6 : 4, j = i < 6 ? i : i - 6;
		double q[6];
		q[0] = v[2 * j]; q[1] = v[2 * j + 1];
		q[4] = v[2 * ((j + 1) % n)]; q[5] = v[2 * ((j + 1) % n) + 1];
		q[2] = 0.5 * (q[0] + q[4]); q[3] = 0.5 * (q[1] + q[5]);
		qs[i] = mvt_bezier(q, 3);
		qv[i] = qs[i];
	}
	p[0].curves = qv; p[0].curve_count = 6;
	p[1].curves = qv + 6; p[1].curve_count = 4;
	mvt_make(p, 2, &mq);
	for (i = 0; i < 10; i++) mvt_hv[i] = sin(1.7 * i) + 0.3 * i;
	for (k = 0; k < 5 && m && mq; k++)
	{
		double ref[10] = { 0 }, s = 0, comb = 0;
		qaws_scalar w[10], wq[20], v = 0;
		qaws_vec2 pt;
		pt.x = (qaws_scalar)pts[2 * k]; pt.y = (qaws_scalar)pts[2 * k + 1];
		mvt_floater(L, 6, pts[2 * k], pts[2 * k + 1], ref);
		mvt_floater(Hcw, 4, pts[2 * k], pts[2 * k + 1], ref + 6);
		/* Hcw runs the other way round: its weights back in H's order */
		{
			double t = ref[7];
			ref[7] = ref[9];
			ref[9] = t;
		}
		for (i = 0; i < 10; i++) s += ref[i];
		qaws_mean_value_coordinates(m, pt, w);
		qaws_mean_value_coordinates(mq, pt, wq);
		for (i = 0; i < 10; i++)
		{
			/* a corner of the quadratic edges: itself and half of each middle point beside it */
			unsigned int n = i < 6 ? 6 : 4, base = i < 6 ? 0 : 12, j = i < 6 ? i : i - 6;
			double quad = (double)wq[base + 2 * j] + 0.5 * ((double)wq[base + 2 * j + 1] + (double)wq[base + 2 * ((j + n - 1) % n) + 1]);
			if (fabs((double)w[i] - ref[i] / s) > worst_cf) worst_cf = fabs((double)w[i] - ref[i] / s);
			if (fabs(quad - ref[i] / s) > worst_quad) worst_quad = fabs(quad - ref[i] / s);
			comb += (double)w[i] * mvt_hv[i];
		}
		qaws_mean_value_interpolate(m, pt, mvt_pl_fn, NULL, 1, &v);
		if (fabs((double)v - comb) > worst_int) worst_int = fabs((double)v - comb);
	}
	sprintf(msg, "an L with a hole: %u handles (10); against Floater's closed form %.2g; its edges as quadratic Beziers %.2g; the interpolant of piecewise linear data %.2g",
		m ? qaws_mean_value_get_handle_count(m) : 0, worst_cf, worst_quad, worst_int);
	TEST_ASSERT(m && mq && qaws_mean_value_get_handle_count(m) == 10 && qaws_mean_value_get_handle_count(mq) == 20 &&
		worst_cf <= tol && worst_quad <= (MVT_F32 ? 2e-4 : 1e-11) && worst_int <= (MVT_F32 ? 2e-4 : 1e-11), msg);
	/* inside a convex polygon every coordinate is positive */
	{
		double hex[12];
		qaws_curve* hc;
		qaws_curve const* hv;
		qaws_mean_value* mh = NULL;
		unsigned int neg = 0;
		for (i = 0; i < 6; i++) { hex[2 * i] = cos(MVT_PI / 3 * i + 0.1); hex[2 * i + 1] = 0.7 * sin(MVT_PI / 3 * i + 0.1); }
		hc = mvt_polyline(hex, 6); hv = hc;
		p[0].curves = &hv; p[0].curve_count = 1;
		mvt_make(p, 1, &mh);
		for (k = 0; k < 50 && mh; k++)
		{
			qaws_scalar w[6];
			qaws_vec2 pt;
			pt.x = (qaws_scalar)(0.8 * cos(2.4 * k) * (k % 7) / 7.0);
			pt.y = (qaws_scalar)(0.5 * sin(2.4 * k) * (k % 7) / 7.0);
			qaws_mean_value_coordinates(mh, pt, w);
			for (i = 0; i < 6; i++) if (!(w[i] > 0)) neg++;
		}
		sprintf(msg, "a convex hexagon: coordinates positive at 50 points (%u not)", neg);
		TEST_ASSERT(mh && neg == 0, msg);
		qaws_mean_value_destroy(mh);
		qaws_curve_destroy(hc);
	}
	qaws_mean_value_destroy(m);
	qaws_mean_value_destroy(mq);
	qaws_curve_destroy(cs[0]); qaws_curve_destroy(cs[1]);
	for (i = 0; i < 10; i++) qaws_curve_destroy(qs[i]);
}

/* ------------------------------------------------------------------------ */

/* the shapes of every kind, as paths into cs / cv (returns the path count) */
static unsigned int mvt_shape(unsigned int which, qaws_curve** cs, qaws_curve const** cv, qaws_path_2d* p, unsigned int* nc, char const** name)
{
	unsigned int i;
	*nc = 0;
	switch (which)
	{
	case 0:
		*name = "blob of cubic Beziers";
		for (i = 0; i < 4; i++) cs[(*nc)++] = mvt_bezier(mvt_blob[i], 4);
		for (i = 0; i < *nc; i++) cv[i] = cs[i];
		p[0].curves = cv; p[0].curve_count = 4; p[0].closed = 1;
		return 1;
	case 1:
		*name = "NURBS ellipse with a round hole";
		qaws_curve_create_ellipse_2d(&cs[(*nc)++], (qaws_scalar)0.1, (qaws_scalar)-0.2, (qaws_scalar)1.5, (qaws_scalar)0.8, (qaws_scalar)0.3);
		cs[(*nc)++] = mvt_circle(0.5, -0.1, 0.3, 0);
		cv[0] = cs[0]; cv[1] = cs[1];
		p[0].curves = cv; p[0].curve_count = 1; p[0].closed = 1;
		p[1].curves = cv + 1; p[1].curve_count = 1; p[1].closed = 1;
		return 2;
	case 2:
	{
		/* a lens: two rational quadratic Beziers (circular arcs, weight cos 60) */
		qaws_rational_bezier_desc d;
		qaws_scalar a[6] = { -1, 0, 0, (qaws_scalar)-1.7320508075688772, 1, 0 }, b[6] = { 1, 0, 0, (qaws_scalar)1.7320508075688772, -1, 0 };
		qaws_scalar w[3] = { 1, (qaws_scalar)0.5, 1 };
		*name = "lens of rational Beziers";
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_point_count = 3; d.weights = w; d.weight_count = 3;
		d.control_points = a; qaws_curve_create_rational_bezier(&cs[(*nc)++], &d);
		d.control_points = b; qaws_curve_create_rational_bezier(&cs[(*nc)++], &d);
		cv[0] = cs[0]; cv[1] = cs[1];
		p[0].curves = cv; p[0].curve_count = 2; p[0].closed = 1;
		return 1;
	}
	case 3:
	{
		/* a closed cubic B-spline (clamped, its ends on one point) */
		qaws_bspline_desc d;
		qaws_scalar cp[18] = { 0, 0, 1, -0.5, 2, 0, 2.2, 1, 1.4, 1.6, 0.6, 1.2, -0.4, 1.5, -0.6, 0.4, 0, 0 };
		qaws_scalar kn[13] = { 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6 };
		*name = "closed cubic B-spline";
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = 9; d.knots = kn; d.knot_count = 13;
		qaws_curve_create_bspline(&cs[(*nc)++], &d);
		cv[0] = cs[0];
		p[0].curves = cv; p[0].curve_count = 1; p[0].closed = 1;
		return 1;
	}
	case 4:
	{
		/* a D: an arc and a Hermite cubic (fitted, exactly) */
		qaws_arc_segment g;
		qaws_arc_desc ad;
		qaws_hermite_desc hd;
		qaws_scalar hp[4] = { 0, -1, 0, 1 }, ht[4] = { -1, 1, 1, 1 };
		*name = "an arc and a Hermite cubic";
		memset(&g, 0, sizeof(g));
		g.radius = 1; g.angle_start = (qaws_scalar)(MVT_PI / 2); g.angle_end = (qaws_scalar)(3 * MVT_PI / 2);
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
		qaws_curve_create_arc(&cs[(*nc)++], &ad);
		memset(&hd, 0, sizeof(hd));
		hd.dimension = QAWS_DIMENSION_2D; hd.degree = 3; hd.points = hp; hd.derivatives = ht; hd.point_count = 2; hd.derivative_count = 2;
		qaws_curve_create_hermite(&cs[(*nc)++], &hd);
		cv[0] = cs[0]; cv[1] = cs[1];
		p[0].curves = cv; p[0].curve_count = 2; p[0].closed = 1;
		return 1;
	}
	}
	return 0;
}

static void test_kinds(void)
{
	unsigned int which;
	char msg[300];
	for (which = 0; which < 5; which++)
	{
		qaws_curve* cs[8] = { 0 };
		qaws_curve const* cv[8];
		qaws_path_2d p[2];
		qaws_mean_value* m = NULL;
		qaws_region_sampler* rs = NULL;
		qaws_region_sampler_desc rd;
		qaws_vec2 pts[40];
		unsigned int nc = 0, np, k, bad = 0;
		char const* name = "";
		double ws = 0, wr = 0, size;
		qaws_status st;
		np = mvt_shape(which, cs, cv, p, &nc, &name);
		st = mvt_make(p, np, &m);
		memset(&rd, 0, sizeof(rd));
		rd.paths = p; rd.path_count = np; rd.fill_rule = QAWS_FILL_EVEN_ODD;
		qaws_region_sampler_create(&rd, &rs);
		size = rs ? sqrt((double)qaws_region_sampler_get_area(rs)) : 1;
		if (rs) qaws_region_sampler_generate(rs, 40, QAWS_REGION_SAMPLE_SOBOL, 7, pts);
		for (k = 0; k < 40 && m && rs; k++)
		{
			double se, re;
			if (mvt_check(m, pts[k].x, pts[k].y, &se, &re) != QAWS_STATUS_OK) bad++;
			if (se > ws) ws = se;
			if (re > wr) wr = re;
		}
		sprintf(msg, "%s: %u handles; at 40 points the coordinates sum to 1 within %.2g and reproduce the point within %.2g of the size (%u failed)",
			name, m ? qaws_mean_value_get_handle_count(m) : 0, ws, wr / size, bad);
		TEST_ASSERT(st == QAWS_STATUS_OK && rs && bad == 0 && ws <= (MVT_F32 ? 1e-5 : 1e-13) && wr <= (MVT_F32 ? 1e-4 : 1e-11) * size, msg);
		qaws_region_sampler_destroy(rs);
		qaws_mean_value_destroy(m);
		for (k = 0; k < nc; k++) qaws_curve_destroy(cs[k]);
	}
}

/* ------------------------------------------------------------------------ */

static void mvt_xy_fn(void* user, unsigned int path, unsigned int curve, qaws_scalar t, qaws_vec2 p, qaws_scalar* out)
{
	(void)user; (void)path; (void)curve; (void)t;
	out[0] = (qaws_scalar)(cos(3 * (double)p.x) * exp((double)p.y));
}

static void test_boundary(void)
{
	qaws_curve* cs[4];
	qaws_curve const* cv[4];
	qaws_path_2d p;
	qaws_mean_value* m = NULL;
	unsigned int i, k, n;
	char msg[300];
	for (i = 0; i < 4; i++) { cs[i] = mvt_bezier(mvt_blob[i], 4); cv[i] = cs[i]; }
	p.curves = cv; p.curve_count = 4; p.closed = 1;
	mvt_make(&p, 1, &m);
	n = m ? qaws_mean_value_get_handle_count(m) : 0;
	if (m && n == 12)
	{
		/* on curve 1 at t = 0.3: the Bernstein values of its handles 3..6 */
		double t = 0.3, b[4], worst_on = 0, gaps[4];
		qaws_eval_result_2d e, e2;
		qaws_vec2 pt, nrm;
		qaws_scalar w[12], v0 = 0, f0 = 0;
		double ds[4] = { 1e-2, 1e-4, 1e-6, 1e-8 };
		b[0] = (1 - t) * (1 - t) * (1 - t); b[1] = 3 * t * (1 - t) * (1 - t); b[2] = 3 * t * t * (1 - t); b[3] = t * t * t;
		qaws_curve_evaluate_2d(&e, cs[1], (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		qaws_mean_value_coordinates(m, e.position, w);
		for (i = 0; i < 12; i++)
		{
			double want = i >= 3 && i <= 6 ? b[i - 3] : 0;
			if (fabs((double)w[i] - want) > worst_on) worst_on = fabs((double)w[i] - want);
		}
		sprintf(msg, "on the boundary: the coordinates are the curve's Bernstein values within %.2g", worst_on);
		TEST_ASSERT(worst_on <= (MVT_F32 ? 1e-5 : 1e-13), msg);
		/* from inside, along the inward normal: the interpolant tends to f there */
		mvt_xy_fn(NULL, 0, 1, (qaws_scalar)t, e.position, &f0);
		qaws_mean_value_interpolate(m, e.position, mvt_xy_fn, NULL, 1, &v0);
		{
			double l = hypot(e.d1.x, e.d1.y);
			nrm.x = (qaws_scalar)(-e.d1.y / l); nrm.y = (qaws_scalar)(e.d1.x / l);
		}
		for (k = 0; k < 4; k++)
		{
			qaws_scalar v = 0;
			pt.x = (qaws_scalar)(e.position.x + ds[k] * nrm.x);
			pt.y = (qaws_scalar)(e.position.y + ds[k] * nrm.y);
			qaws_mean_value_interpolate(m, pt, mvt_xy_fn, NULL, 1, &v);
			gaps[k] = fabs((double)v - (double)f0);
		}
		(void)e2;
		sprintf(msg, "towards the boundary (1e-2 to 1e-8 inside): the interpolant's gap to f there %.2g %.2g %.2g %.2g (on it %.2g)",
			gaps[0], gaps[1], gaps[2], gaps[3], fabs((double)v0 - (double)f0));
		TEST_ASSERT(fabs((double)v0 - (double)f0) <= (MVT_F32 ? 1e-6 : 1e-14) && gaps[1] < gaps[0] && gaps[2] < gaps[1] &&
			gaps[2] <= (MVT_F32 ? 1e-3 : 1e-4) && gaps[3] <= (MVT_F32 ? 1e-3 : 1e-6), msg);
	}
	else
		TEST_ASSERT(0, "blob for the boundary checks");
	qaws_mean_value_destroy(m);
	for (i = 0; i < 4; i++) qaws_curve_destroy(cs[i]);
}

static void test_orientation(void)
{
	double sq[8] = { -0.2, -0.2, 0.2, -0.2, 0.2, 0.2, -0.2, 0.2 };
	double sqr[8] = { -0.2, -0.2, -0.2, 0.2, 0.2, 0.2, 0.2, -0.2 };
	double pts[6] = { 0.5, 0.1, -0.7, -0.6, 0, 0.9 };
	unsigned int o, k, i;
	double v[4][3], worst = 0;
	char msg[200];
	for (o = 0; o < 4; o++)
	{
		qaws_curve* cs[2];
		qaws_curve const* cv[2];
		qaws_path_2d p[2];
		qaws_mean_value* m = NULL;
		cs[0] = mvt_circle(0, 0, 1, (o & 1) == 0);
		cs[1] = mvt_polyline((o & 2) ? sq : sqr, 4);
		cv[0] = cs[0]; cv[1] = cs[1];
		p[0].curves = cv; p[0].curve_count = 1; p[0].closed = 1;
		p[1].curves = cv + 1; p[1].curve_count = 1; p[1].closed = 1;
		mvt_make(p, 2, &m);
		for (k = 0; k < 3; k++)
		{
			qaws_vec2 pt;
			qaws_scalar val = 0;
			pt.x = (qaws_scalar)pts[2 * k]; pt.y = (qaws_scalar)pts[2 * k + 1];
			qaws_mean_value_interpolate(m, pt, mvt_disc_fn, NULL, 1, &val);
			v[o][k] = (double)val;
		}
		qaws_mean_value_destroy(m);
		qaws_curve_destroy(cs[0]); qaws_curve_destroy(cs[1]);
	}
	for (o = 1; o < 4; o++)
		for (i = 0; i < 3; i++)
			if (fabs(v[o][i] - v[0][i]) > worst) worst = fabs(v[o][i] - v[0][i]);
	sprintf(msg, "a disc with a hole, both loops either way round: the interpolant the same within %.2g", worst);
	TEST_ASSERT(worst <= (MVT_F32 ? 1e-5 : 1e-13), msg);
}

/* ------------------------------------------------------------------------ */

static void mvt_backwards(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int b = count;
	(void)user;
	while (b > 0)
	{
		unsigned int a = b > 3 ? b - 3 : 0;
		task(ctx, a, b);
		b = a;
	}
}

static void test_deform(void)
{
	qaws_curve* cs[8] = { 0 };
	qaws_curve const* cv[8];
	qaws_path_2d p[2];
	qaws_mean_value* m = NULL;
	qaws_region_sampler* rs = NULL;
	qaws_region_sampler_desc rd;
	qaws_vec2 pts[30], out[30], out2[30], *hv = NULL, *moved = NULL;
	unsigned int nc = 0, np, k, n = 0, i;
	char const* name;
	double we = 0, wt = 0, wa = 0, wx = 0;
	char msg[300];
	np = mvt_shape(1, cs, cv, p, &nc, &name);
	mvt_make(p, np, &m);
	memset(&rd, 0, sizeof(rd));
	rd.paths = p; rd.path_count = np; rd.fill_rule = QAWS_FILL_EVEN_ODD;
	qaws_region_sampler_create(&rd, &rs);
	if (m && rs)
	{
		n = qaws_mean_value_get_handle_count(m);
		hv = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
		moved = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
		for (i = 0; i < n; i++) qaws_mean_value_get_handle(m, i, &hv[i]);
		qaws_region_sampler_generate(rs, 30, QAWS_REGION_SAMPLE_R2, 3, pts);
		/* unmoved */
		qaws_mean_value_deform(m, hv, pts, 30, out);
		for (k = 0; k < 30; k++) if (hypot(out[k].x - pts[k].x, out[k].y - pts[k].y) > we) we = hypot(out[k].x - pts[k].x, out[k].y - pts[k].y);
		/* translated */
		for (i = 0; i < n; i++) { moved[i].x = hv[i].x + (qaws_scalar)0.7; moved[i].y = hv[i].y - (qaws_scalar)0.2; }
		qaws_mean_value_deform(m, moved, pts, 30, out);
		for (k = 0; k < 30; k++) { double d = hypot(out[k].x - pts[k].x - 0.7, out[k].y - pts[k].y + 0.2); if (d > wt) wt = d; }
		/* affine */
		for (i = 0; i < n; i++)
		{
			moved[i].x = (qaws_scalar)(1.2 * hv[i].x - 0.4 * hv[i].y + 0.1);
			moved[i].y = (qaws_scalar)(0.3 * hv[i].x + 0.9 * hv[i].y - 0.3);
		}
		qaws_mean_value_deform(m, moved, pts, 30, out);
		for (k = 0; k < 30; k++)
		{
			double ex = 1.2 * pts[k].x - 0.4 * pts[k].y + 0.1, ey = 0.3 * pts[k].x + 0.9 * pts[k].y - 0.3, d = hypot(out[k].x - ex, out[k].y - ey);
			if (d > wa) wa = d;
		}
	}
	sprintf(msg, "deformations of the ellipse with a hole (%u handles): unmoved %.2g, translated %.2g, affine %.2g", n, we, wt, wa);
	TEST_ASSERT(m && rs && we <= (MVT_F32 ? 1e-4 : 1e-11) && wt <= (MVT_F32 ? 1e-4 : 1e-11) && wa <= (MVT_F32 ? 1e-4 : 1e-11), msg);
	if (m && rs)
	{
		/* a bend: the handles moved nonlinearly; the same through an executor */
		qaws_mean_value_desc d;
		qaws_mean_value* me = NULL;
		qaws_batch_executor ex;
		for (i = 0; i < n; i++)
		{
			moved[i].x = (qaws_scalar)(hv[i].x + 0.2 * hv[i].y * hv[i].y);
			moved[i].y = (qaws_scalar)(hv[i].y + 0.1 * sin(2 * hv[i].x));
		}
		qaws_mean_value_deform(m, moved, pts, 30, out);
		memset(&d, 0, sizeof(d));
		d.paths = p; d.path_count = np;
		ex.parallel_for = mvt_backwards; ex.user = NULL;
		d.executor = &ex;
		qaws_mean_value_create(&d, &me);
		if (me) qaws_mean_value_deform(me, moved, pts, 30, out2);
		for (k = 0; k < 30 && me; k++) if (hypot(out[k].x - out2[k].x, out[k].y - out2[k].y) > wx) wx = hypot(out[k].x - out2[k].x, out[k].y - out2[k].y);
		sprintf(msg, "a bend through an executor (chunks backwards): the serial points within %.2g", wx);
		TEST_ASSERT(me && wx == 0, msg);
		qaws_mean_value_destroy(me);
	}
	free(hv); free(moved);
	qaws_region_sampler_destroy(rs);
	qaws_mean_value_destroy(m);
	for (k = 0; k < nc; k++) qaws_curve_destroy(cs[k]);
}

int test_100_mean_value_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 100: mean value coordinates over regions bounded by curves\n");
	test_reference();
	test_polygons();
	test_kinds();
	test_boundary();
	test_orientation();
	test_deform();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
