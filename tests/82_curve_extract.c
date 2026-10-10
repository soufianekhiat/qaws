/*
 * Test 82: Exact pieces of curves (qaws_curve_extract), reverse, closedness
 *
 *   For one curve of every kind and several parameter pairs (forward,
 *   reversed, whole domain, span boundaries):
 *     - the piece starts at C(t0) and ends at C(t1);
 *     - every point of the piece lies on the source between t0 and t1, and
 *       every point of the source there lies on the piece (two-way distance,
 *       found by dense sampling and a golden-section refinement);
 *     - the piece has the expected kind.
 *   Also: reverse for every kind, split inside a Hermite span, and is_closed
 *   for closed Bezier, B-spline and NURBS loops.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#if QAWS_SCALAR_IS_FLOAT
#define CXT_TOL 2e-4
#else
#define CXT_TOL 1e-9
#endif

typedef struct cxt_case
{
	char const* name;
	qaws_curve* curve;
	qaws_curve_kind piece_kind;   /* INVALID: either Bezier or B-spline */
	double tol;                   /* relative; 0 = CXT_TOL */
	qaws_curve const* keep;       /* source of a reparameterized wrapper */
} cxt_case;

static void cxt_pos(qaws_curve const* c, double t, double* p)
{
	qaws_eval_result_2d r;
	qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION, &r);
	p[0] = r.position.x;
	p[1] = r.position.y;
}

static double cxt_dist2(qaws_curve const* c, double t, double const* q)
{
	double p[2];
	cxt_pos(c, t, p);
	return (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]);
}

/* distance from q to c restricted to [lo, hi] */
static double cxt_distance(qaws_curve const* c, double lo, double hi, double const* q)
{
	unsigned int const n = 600;
	unsigned int i, bi = 0;
	double best = 1e300, a, b, g = 0.6180339887498949;
	for (i = 0; i <= n; i++)
	{
		double d = cxt_dist2(c, lo + (hi - lo) * i / n, q);
		if (d < best) { best = d; bi = i; }
	}
	a = lo + (hi - lo) * (bi > 0 ? bi - 1 : 0) / n;
	b = lo + (hi - lo) * (bi < n ? bi + 1 : n) / n;
	for (i = 0; i < 80; i++)
	{
		double x1 = b - g * (b - a), x2 = a + g * (b - a);
		if (cxt_dist2(c, x1, q) < cxt_dist2(c, x2, q)) b = x2; else a = x1;
	}
	a = cxt_dist2(c, 0.5 * (a + b), q);
	return sqrt(a < best ? a : best);
}

static double cxt_extent(qaws_curve const* c)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, p[2];
	unsigned int i, k;
	for (i = 0; i <= 200; i++)
	{
		cxt_pos(c, r.min_value + (r.max_value - r.min_value) * i / 200.0, p);
		for (k = 0; k < 2; k++)
		{
			if (p[k] < lo[k]) lo[k] = p[k];
			if (p[k] > hi[k]) hi[k] = p[k];
		}
	}
	return sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]));
}

/* worst two-way distance between piece e and source c on [t0, t1], and the
   end-point errors */
static void cxt_compare(qaws_curve const* c, double t0, double t1, qaws_curve const* e,
	double* out_ends, double* out_inside)
{
	qaws_range re = qaws_curve_get_parameter_range(e);
	double lo = t0 < t1 ? t0 : t1, hi = t0 < t1 ? t1 : t0, p[2], q[2], w = 0.0;
	unsigned int i;
	cxt_pos(c, t0, p); cxt_pos(e, re.min_value, q);
	*out_ends = sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]));
	cxt_pos(c, t1, p); cxt_pos(e, re.max_value, q);
	w = sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]));
	if (w > *out_ends) *out_ends = w;
	w = 0.0;
	for (i = 0; i <= 48; i++)
	{
		double d;
		cxt_pos(e, re.min_value + (re.max_value - re.min_value) * i / 48.0, q);
		d = cxt_distance(c, lo, hi, q);
		if (d > w) w = d;
		cxt_pos(c, lo + (hi - lo) * i / 48.0, q);
		d = cxt_distance(e, re.min_value, re.max_value, q);
		if (d > w) w = d;
	}
	*out_inside = w;
}

/* -------------------------------------------------------------------------- */
/*  One curve of every kind                                                   */
/* -------------------------------------------------------------------------- */

static qaws_curve* cxt_bezier(void)
{
	qaws_scalar cp[] = { 0, 0, 1, 3, 3, -1, 4, 2, 5, 0, 6, 1 };
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 5; d.control_points = cp; d.control_point_count = 6;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

static qaws_curve* cxt_rational_bezier(void)
{
	qaws_scalar cp[] = { 2, 0, 2, 2, 0, 2 }, w[] = { 1, (qaws_scalar)0.70710678118654752, 1 };
	qaws_rational_bezier_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = cp; d.control_point_count = 3;
	d.weights = w; d.weight_count = 3;
	qaws_curve_create_rational_bezier(&d, &c);
	return c;
}

static qaws_curve* cxt_bspline(int uniform)
{
	qaws_scalar cp[] = { 0, 0, 1, 2, 2, -1, 3, 1, 4, 3, 5, 0, 6, 2 };
	qaws_scalar kn[] = { 0, 0, 0, 0, (qaws_scalar)0.3, (qaws_scalar)0.45, (qaws_scalar)0.8, 1, 1, 1, 1 };
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = 7;
	if (uniform)
		d.is_uniform = 1;
	else
	{
		d.knots = kn; d.knot_count = 11;
	}
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static qaws_curve* cxt_nurbs_circle(void)
{
	static double const ux[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, uy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	qaws_scalar cp[18], w[9], kn[12] = { 0, 0, 0, (qaws_scalar)0.25, (qaws_scalar)0.25, (qaws_scalar)0.5,
		(qaws_scalar)0.5, (qaws_scalar)0.75, (qaws_scalar)0.75, 1, 1, 1 };
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 9; i++)
	{
		cp[2 * i] = (qaws_scalar)(3 * ux[i] + 1);
		cp[2 * i + 1] = (qaws_scalar)(2 * uy[i]);
		w[i] = (qaws_scalar)(i % 2 ? sqrt(0.5) : 1);
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = cp; d.control_point_count = 9;
	d.knots = kn; d.knot_count = 12; d.weights = w; d.weight_count = 9;
	qaws_curve_create_nurbs(&d, &c);
	return c;
}

static qaws_curve* cxt_hermite(void)
{
	qaws_scalar p[] = { 0, 0, 2, 1, 3, -1, 5, 0 }, v[] = { 1, 2, 2, 0, 1, -2, 1, 1 };
	qaws_hermite_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.points = p; d.derivatives = v;
	d.point_count = 4; d.derivative_count = 4;
	qaws_curve_create_hermite(&d, &c);
	return c;
}

static qaws_curve* cxt_catmull_rom(qaws_parameterization par, int closed)
{
	qaws_scalar p[] = { 0, 0, 1, 2, 3, 2, 4, 0, 2, -1, 1, -0.5f };
	qaws_catmull_rom_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.control_points = p; d.control_point_count = 6;
	d.parameterization = par; d.closed = closed;
	qaws_curve_create_catmull_rom(&d, &c);
	return c;
}

static qaws_curve* cxt_trajectory(void)
{
	qaws_scalar p[] = { 0, 0, 2, 3, 5, 1, 6, 4 }, t[] = { 0, 1, (qaws_scalar)2.5, 3 };
	qaws_trajectory_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.key_positions = p; d.key_count = 4;
	d.key_times = t; d.key_time_count = 4;
	qaws_curve_create_trajectory(&d, &c);
	return c;
}

static qaws_curve* cxt_yuksel(qaws_yuksel_mode mode)
{
	qaws_scalar p[] = { 0, 0, 1, 2, 3, 2, 4, 0, 6, 1 };
	qaws_yuksel_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.control_points = p; d.control_point_count = 5; d.mode = mode;
	qaws_curve_create_yuksel(&d, &c);
	return c;
}

static qaws_curve* cxt_polynomial(void)
{
	qaws_scalar co[] = { 1, 0, 2, 1, -1, 3, (qaws_scalar)0.5, -2, (qaws_scalar)0.25, (qaws_scalar)0.5 };
	qaws_polynomial_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 4; d.coefficients = co; d.coefficient_count = 5;
	d.t_min = -1; d.t_max = 2;
	qaws_curve_create_polynomial(&d, &c);
	return c;
}

static qaws_curve* cxt_arc(void)
{
	qaws_arc_segment s[2];
	qaws_arc_desc d;
	qaws_curve* c = NULL;
	memset(s, 0, sizeof(s));
	s[0].center[0] = 0; s[0].center[1] = 0; s[0].radius = 2;
	s[0].angle_start = 0; s[0].angle_end = (qaws_scalar)1.5;
	s[1].center[0] = (qaws_scalar)(2 * cos(1.5) - cos(1.5)); s[1].center[1] = (qaws_scalar)(2 * sin(1.5) - sin(1.5));
	s[1].radius = 1; s[1].angle_start = (qaws_scalar)1.5; s[1].angle_end = (qaws_scalar)3.5;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.segments = s; d.segment_count = 2;
	qaws_curve_create_arc(&d, &c);
	return c;
}

static qaws_curve* cxt_clothoid(void)
{
	qaws_clothoid_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.origin_x = 1; d.origin_y = -1; d.start_angle = (qaws_scalar)0.3;
	d.start_curvature = (qaws_scalar)-0.2; d.end_curvature = (qaws_scalar)1.1; d.length = 6;
	qaws_curve_create_clothoid(&d, &c);
	return c;
}

static qaws_curve* cxt_subdivision(void)
{
	qaws_scalar p[] = { 0, 0, 3, 0, 4, 2, 2, 4, -1, 2 };
	qaws_subdivision_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_3;
	d.control_points = p; d.control_point_count = 5; d.closed = 1; d.refinement_levels = 3;
	qaws_curve_create_subdivision(&d, &c);
	return c;
}

static qaws_curve* cxt_composite(void)
{
	qaws_scalar l[] = { 0, 0, 2, 0 };
	qaws_scalar b[] = { 2, 0, 3, 0, 4, 1, 4, 2 };
	qaws_arc_segment s;
	qaws_bezier_desc bd;
	qaws_arc_desc ad;
	qaws_composite_desc d;
	qaws_curve* seg[3] = { NULL, NULL, NULL };
	qaws_curve* c = NULL;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = l; bd.control_point_count = 2;
	qaws_curve_create_bezier(&bd, &seg[0]);
	bd.degree = 3; bd.control_points = b; bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &seg[1]);
	memset(&s, 0, sizeof(s));
	s.center[0] = 3; s.center[1] = 2; s.radius = 1; s.angle_start = 0; s.angle_end = (qaws_scalar)3.14159265358979;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &s; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, &seg[2]);
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.segments = seg; d.segment_count = 3;
	qaws_curve_create_composite(&d, &c);
	return c;
}

static unsigned int cxt_cases(cxt_case* cs)
{
	unsigned int n = 0;
	qaws_curve* src;
#define CXT_ADD(nm, cv, kd, tl) do { cs[n].name = nm; cs[n].curve = cv; cs[n].piece_kind = kd; cs[n].tol = tl; cs[n].keep = NULL; n++; } while (0)
	CXT_ADD("bezier", cxt_bezier(), QAWS_CURVE_KIND_BEZIER, 0);
	CXT_ADD("rational bezier", cxt_rational_bezier(), QAWS_CURVE_KIND_RATIONAL_BEZIER, 0);
	CXT_ADD("bspline", cxt_bspline(0), QAWS_CURVE_KIND_BSPLINE, 0);
	CXT_ADD("uniform bspline", cxt_bspline(1), QAWS_CURVE_KIND_BSPLINE, 0);
	CXT_ADD("nurbs circle", cxt_nurbs_circle(), QAWS_CURVE_KIND_NURBS, 0);
	CXT_ADD("hermite", cxt_hermite(), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("catmull-rom uniform", cxt_catmull_rom(QAWS_PARAMETERIZATION_UNIFORM, 0), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("catmull-rom centripetal", cxt_catmull_rom(QAWS_PARAMETERIZATION_CENTRIPETAL, 0), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("catmull-rom chordal closed", cxt_catmull_rom(QAWS_PARAMETERIZATION_CHORDAL, 1), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("trajectory", cxt_trajectory(), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("yuksel bezier", cxt_yuksel(QAWS_YUKSEL_MODE_BEZIER), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("yuksel circular", cxt_yuksel(QAWS_YUKSEL_MODE_CIRCULAR), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("yuksel hybrid", cxt_yuksel(QAWS_YUKSEL_MODE_HYBRID), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("polynomial", cxt_polynomial(), QAWS_CURVE_KIND_POLYNOMIAL, 0);
	CXT_ADD("arc", cxt_arc(), QAWS_CURVE_KIND_ARC, 0);
	CXT_ADD("clothoid", cxt_clothoid(), QAWS_CURVE_KIND_CLOTHOID, 0);
	CXT_ADD("subdivision closed", cxt_subdivision(), QAWS_CURVE_KIND_INVALID, 0);
	CXT_ADD("composite", cxt_composite(), QAWS_CURVE_KIND_COMPOSITE, 0);
	src = cxt_bezier();
	{
		qaws_curve* rp = NULL;
		qaws_curve_reparameterize_arc_length(src, 0, &rp);
		CXT_ADD("reparameterized", rp, QAWS_CURVE_KIND_BEZIER, 0);
		cs[n - 1].keep = src;
	}
#undef CXT_ADD
	return n;
}

/* -------------------------------------------------------------------------- */

static void test_extract(cxt_case const* cs, unsigned int n)
{
	static double const pairs[][2] = {
		{ 0.1, 0.7 }, { 0.0, 0.35 }, { 0.42, 1.0 }, { 0.8, 0.2 }, { 0.0, 1.0 },
		{ 1.0 / 3.0, 2.0 / 3.0 }, { 0.5, 0.5000001 }, { 0.93, 0.02 }
	};
	unsigned int i, j, bad_status = 0, bad_ends = 0, bad_inside = 0, bad_kind = 0, total = 0;
	double worst = 0.0;
	char msg[256];
	for (i = 0; i < n; i++)
	{
		qaws_range r;
		double ext, tol;
		if (!cs[i].curve)
		{
			printf("  FAIL: %s not created\n", cs[i].name);
			bad_status++;
			continue;
		}
		r = qaws_curve_get_parameter_range(cs[i].curve);
		ext = cxt_extent(cs[i].curve);
		tol = (cs[i].tol > 0 ? cs[i].tol : CXT_TOL) * ext;
		for (j = 0; j < sizeof(pairs) / sizeof(pairs[0]); j++)
		{
			double t0 = r.min_value + (r.max_value - r.min_value) * pairs[j][0];
			double t1 = r.min_value + (r.max_value - r.min_value) * pairs[j][1];
			double ends, inside;
			qaws_curve* e = NULL;
			qaws_status s;
#if QAWS_SCALAR_IS_FLOAT
			if (j == 6) continue;   /* a 1e-7 piece is below float rounding */
#endif
			total++;
			s = qaws_curve_extract(cs[i].curve, (qaws_scalar)t0, (qaws_scalar)t1, &e);
			if (s != QAWS_STATUS_OK)
			{
				printf("    %s [%.3f, %.3f]: status %d\n", cs[i].name, pairs[j][0], pairs[j][1], (int)s);
				bad_status++;
				continue;
			}
			cxt_compare(cs[i].curve, t0, t1, e, &ends, &inside);
			if (ends > tol)
			{
				printf("    %s [%.3f, %.3f]: end points off by %.2e\n", cs[i].name, pairs[j][0], pairs[j][1], ends);
				bad_ends++;
			}
			if (inside > tol)
			{
				printf("    %s [%.3f, %.3f]: piece off the source by %.2e\n", cs[i].name, pairs[j][0], pairs[j][1], inside);
				bad_inside++;
			}
			if (inside / ext > worst) worst = inside / ext;
			if (cs[i].piece_kind == QAWS_CURVE_KIND_INVALID
				? (qaws_curve_get_kind(e) != QAWS_CURVE_KIND_BEZIER && qaws_curve_get_kind(e) != QAWS_CURVE_KIND_BSPLINE)
				: (qaws_curve_get_kind(e) != cs[i].piece_kind &&
				   !(cs[i].piece_kind == QAWS_CURVE_KIND_COMPOSITE)))
			{
				printf("    %s [%.3f, %.3f]: piece of kind %d\n", cs[i].name, pairs[j][0], pairs[j][1], (int)qaws_curve_get_kind(e));
				bad_kind++;
			}
			qaws_curve_destroy(e);
		}
	}
	printf("    %u pieces of %u kinds, worst distance to the source %.1e of the extent\n", total, n, worst);
	sprintf(msg, "every extract succeeds (%u failed)", bad_status);
	TEST_ASSERT(bad_status == 0, msg);
	TEST_ASSERT(bad_ends == 0, "every piece starts at C(t0) and ends at C(t1)");
	TEST_ASSERT(bad_inside == 0, "every piece lies on its source and covers it");
	TEST_ASSERT(bad_kind == 0, "every piece has the expected kind");
}

static void test_reverse(cxt_case const* cs, unsigned int n)
{
	unsigned int i, bad = 0;
	for (i = 0; i < n; i++)
	{
		qaws_curve* rv = NULL;
		qaws_range r, rr;
		double ends, inside, tol;
		if (!cs[i].curve) continue;
		if (qaws_curve_reverse(cs[i].curve, &rv) != QAWS_STATUS_OK)
		{
			printf("    reverse %s failed\n", cs[i].name);
			bad++;
			continue;
		}
		r = qaws_curve_get_parameter_range(cs[i].curve);
		rr = qaws_curve_get_parameter_range(rv);
		(void)rr;
		tol = (cs[i].tol > 0 ? cs[i].tol : CXT_TOL) * cxt_extent(cs[i].curve);
		cxt_compare(cs[i].curve, r.max_value, r.min_value, rv, &ends, &inside);
		if (ends > tol || inside > tol)
		{
			printf("    reverse %s: ends %.2e, inside %.2e\n", cs[i].name, ends, inside);
			bad++;
		}
		qaws_curve_destroy(rv);
	}
	TEST_ASSERT(bad == 0, "reverse works on every kind and traces the curve backward");
}

static void test_hermite_split(void)
{
	qaws_curve* c = cxt_hermite();
	qaws_curve *l = NULL, *r = NULL;
	double ends, inside, tol = CXT_TOL * cxt_extent(c);
	int ok = qaws_curve_split(c, (qaws_scalar)1.4, &l, &r) == QAWS_STATUS_OK;
	if (ok)
	{
		cxt_compare(c, 0.0, 1.4, l, &ends, &inside);
		ok = ends <= tol && inside <= tol;
		cxt_compare(c, 1.4, 3.0, r, &ends, &inside);
		ok = ok && ends <= tol && inside <= tol;
	}
	TEST_ASSERT(ok, "splitting inside a Hermite span keeps both halves on the curve");
	qaws_curve_destroy(l);
	qaws_curve_destroy(r);
	qaws_curve_destroy(c);
}

static void test_closed(void)
{
	qaws_scalar sq[] = { 0, 0, 1, 0, 1, 1, 0, 1, 0, 0 };
	qaws_scalar lp[] = { 0, 0, 3, 1, -1, 2, 0, 0 };
	qaws_bspline_desc bd;
	qaws_bezier_desc zd;
	qaws_curve *poly = NULL, *loop = NULL, *circle = cxt_nurbs_circle(), *open = cxt_bspline(0);
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = sq; bd.control_point_count = 5;
	qaws_curve_create_bspline(&bd, &poly);
	memset(&zd, 0, sizeof(zd));
	zd.dimension = QAWS_DIMENSION_2D; zd.degree = 3; zd.control_points = lp; zd.control_point_count = 4;
	qaws_curve_create_bezier(&zd, &loop);
	TEST_ASSERT(qaws_curve_is_closed(poly) == 1, "a closed degree-1 B-spline polygon is closed");
	TEST_ASSERT(qaws_curve_is_closed(loop) == 1, "a Bezier loop is closed");
	TEST_ASSERT(qaws_curve_is_closed(circle) == 1, "a NURBS circle is closed");
	TEST_ASSERT(qaws_curve_is_closed(open) == 0, "an open B-spline is not closed");
	qaws_curve_destroy(poly);
	qaws_curve_destroy(loop);
	qaws_curve_destroy(circle);
	qaws_curve_destroy(open);
}

int test_82_curve_extract_main(void)
{
	cxt_case cs[32];
	unsigned int n, i;
	g_pass = 0;
	g_fail = 0;
	printf("Test 82: Exact pieces of curves\n");
	n = cxt_cases(cs);
	test_extract(cs, n);
	test_reverse(cs, n);
	test_hermite_split();
	test_closed();
	for (i = 0; i < n; i++)
	{
		qaws_curve_destroy(cs[i].curve);
		if (cs[i].keep) qaws_curve_destroy((qaws_curve*)cs[i].keep);
	}
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
