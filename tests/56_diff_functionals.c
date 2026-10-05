/*
 * Test 56: Integral functionals with gradients and Hessian-vector products
 *
 * Curves: length, bending, curvature squared. Surfaces: area, thin plate,
 * Willmore.
 *   - exact values on simple shapes (segment length, quarter circle,
 *     unit square area)
 *   - gradient matches finite differences of the value
 *   - tangent equals <gradient, direction>
 *   - second tangent equals <direction, H direction>
 *   - HVP matches finite differences of the gradient, and is symmetric
 *   - HVP refused for families that are not linear in their fields
 */

#include "test_diff.h"

#define FN_MAX 64

static qaws_diff_views one_view(qaws_field_view* fv, qaws_scalar* data, unsigned int count, unsigned int comps)
{
	qaws_diff_views v;
	*fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, data, count, comps);
	v.fields = fv;
	v.field_count = 1;
	v.children = NULL;
	v.child_count = 0;
	return v;
}

/* ------------------------------------------------------------------ */
/*  Curves                                                            */
/* ------------------------------------------------------------------ */

static qaws_scalar const g_fn_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_curve* fn_curve(qaws_scalar const* p)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 6;
	d.knots = g_fn_knots;
	d.knot_count = 10;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static void test_curve_values(void)
{
	/* Straight Bezier: length = chord. */
	static qaws_scalar const line[12] = { 0, 0, 0, 1, 1, 0, 2, 2, 0, 3, 3, 0 };
	qaws_bezier_desc bd;
	qaws_curve* c = NULL;
	qaws_scalar len = 0, bend = 0;
	bd.dimension = QAWS_DIMENSION_3D;
	bd.degree = 3;
	bd.control_points = line;
	bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &c);
	qaws_curve_functional_eval(NULL, c, QAWS_FUNCTIONAL_LENGTH, 0, NULL, &len, NULL, NULL);
	qaws_curve_functional_eval(NULL, c, QAWS_FUNCTIONAL_BENDING, 0, NULL, &bend, NULL, NULL);
	TEST_ASSERT(approx_eq(len, (qaws_scalar)(3 * 1.4142135623730951)), "segment length");
	TEST_ASSERT(approx_eq(bend, 0), "a uniformly parameterized segment has no bending");
	qaws_curve_destroy(c);

	/* Rational quarter circle: length pi/2 and curvature energy pi/2. */
	{
		static qaws_scalar const cps[6] = { 1, 0, 1, 1, 0, 1 };
		qaws_scalar w[3] = { 1, (qaws_scalar)0.70710678118654752, 1 };
		qaws_rational_bezier_desc rd;
		qaws_scalar k2 = 0;
		rd.dimension = QAWS_DIMENSION_2D;
		rd.degree = 2;
		rd.control_points = cps;
		rd.control_point_count = 3;
		rd.weights = w;
		rd.weight_count = 3;
		qaws_curve_create_rational_bezier(&rd, &c);
		qaws_curve_functional_eval(NULL, c, QAWS_FUNCTIONAL_LENGTH, 8, NULL, &len, NULL, NULL);
		qaws_curve_functional_eval(NULL, c, QAWS_FUNCTIONAL_CURVATURE_SQUARED, 8, NULL, &k2, NULL, NULL);
		TEST_ASSERT(approx_eq_loose(len, (qaws_scalar)1.5707963267948966), "quarter circle length is pi/2");
		TEST_ASSERT(approx_eq_loose(k2, (qaws_scalar)1.5707963267948966), "unit circle curvature energy is pi/2");
		{
			qaws_scalar dir[6] = { 0 }, hv[6] = { 0 };
			qaws_field_view a, b;
			qaws_diff_views va = one_view(&a, dir, 3, 2), vb = one_view(&b, hv, 3, 2);
			TEST_ASSERT(qaws_curve_functional_hvp(NULL, c, QAWS_FUNCTIONAL_LENGTH, 0, &va, &vb)
				== QAWS_STATUS_UNSUPPORTED_OPERATION, "direct HVP refused for rational curves");
		}
		qaws_curve_destroy(c);
	}
}

static double curve_value(qaws_scalar const* p, qaws_curve_functional f)
{
	qaws_curve* c = fn_curve(p);
	qaws_scalar v = 0;
	qaws_curve_functional_eval(NULL, c, f, 0, NULL, &v, NULL, NULL);
	qaws_curve_destroy(c);
	return v;
}

static void curve_gradient(qaws_scalar const* p, qaws_curve_functional f, qaws_scalar* g)
{
	qaws_curve* c = fn_curve(p);
	qaws_field_view fv;
	qaws_diff_views v;
	memset(g, 0, sizeof(qaws_scalar) * 18);
	v = one_view(&fv, g, 6, 3);
	qaws_curve_functional_gradient(NULL, c, f, 0, &v, NULL);
	qaws_curve_destroy(c);
}

static void check_curve_functional(qaws_curve_functional f, char const* name)
{
	qaws_scalar p[18], d[18], e[18], g[18], hv[18], he[18], pp[18], pm[18], gp[18], gm[18];
	qaws_curve* c;
	qaws_field_view fv, fv2;
	qaws_diff_views vd, vh;
	qaws_scalar val, tan, tan2;
	double h = DIFF_FD_STEP, fd, gd;
	unsigned int i;
	int ok_hvp = 1;

	for (i = 0; i < 6; i++)
	{
		p[3 * i] = (qaws_scalar)i;
		p[3 * i + 1] = (qaws_scalar)(0.6 * sin(1.3 * i));
		p[3 * i + 2] = (qaws_scalar)(0.3 * cos(0.8 * i));
	}
	diff_rand_fill(d, 18);
	diff_rand_fill(e, 18);
	c = fn_curve(p);

	vd = one_view(&fv, d, 6, 3);
	qaws_curve_functional_eval(NULL, c, f, 0, &vd, &val, &tan, &tan2);
	curve_gradient(p, f, g);

	for (i = 0; i < 18; i++)
	{
		pp[i] = (qaws_scalar)(p[i] + h * d[i]);
		pm[i] = (qaws_scalar)(p[i] - h * d[i]);
	}
	fd = (curve_value(pp, f) - curve_value(pm, f)) / (2 * h);
	gd = diff_dot(g, d, 18);
	printf("    %s: value %.6g, tangent %.9g, <grad, d> %.9g, fd %.9g\n", name, (double)val, (double)tan, gd, fd);
	TEST_ASSERT(diff_close(fd, gd, DIFF_TOL * 300), "functional gradient matches finite differences");
	TEST_ASSERT(diff_close(tan, gd, DIFF_TOL * 10), "functional tangent equals <gradient, direction>");

	/* HVP: finite differences of the gradient, symmetry, and tangent2. */
	memset(hv, 0, sizeof(hv));
	vh = one_view(&fv2, hv, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_functional_hvp(NULL, c, f, 0, &vd, &vh));
	curve_gradient(pp, f, gp);
	curve_gradient(pm, f, gm);
	for (i = 0; i < 18; i++)
		if (!diff_close((gp[i] - gm[i]) / (2 * h), hv[i], DIFF_TOL * 300))
			ok_hvp = 0;
	TEST_ASSERT(ok_hvp, "HVP matches finite differences of the gradient");
	TEST_ASSERT(diff_close(tan2, diff_dot(d, hv, 18), DIFF_TOL * 10), "second tangent equals <d, H d>");
	{
		qaws_field_view fe, fh;
		qaws_diff_views ve = one_view(&fe, e, 6, 3);
		qaws_diff_views vhe;
		memset(he, 0, sizeof(he));
		vhe = one_view(&fh, he, 6, 3);
		qaws_curve_functional_hvp(NULL, c, f, 0, &ve, &vhe);
		TEST_ASSERT(diff_close(diff_dot(e, hv, 18), diff_dot(d, he, 18), DIFF_TOL * 10), "HVP is symmetric");
	}
	qaws_curve_destroy(c);
}

/* ------------------------------------------------------------------ */
/*  Surfaces                                                          */
/* ------------------------------------------------------------------ */

static qaws_surface* fn_surface(qaws_scalar const* p)
{
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)p;
	d.u_point_count = 4;
	d.v_point_count = 4;
	qaws_surface_create_bezier(&d, &s);
	return s;
}

static void fill_patch(qaws_scalar* p, double bump)
{
	unsigned int a, b;
	for (a = 0; a < 4; a++)
		for (b = 0; b < 4; b++)
		{
			qaws_scalar* q = &p[(a * 4 + b) * 3];
			q[0] = (qaws_scalar)(a / 3.0);
			q[1] = (qaws_scalar)(b / 3.0);
			q[2] = (qaws_scalar)(bump * sin(1.7 * a + 0.3) * cos(1.1 * b));
		}
}

static void test_surface_values(void)
{
	qaws_scalar p[48];
	qaws_surface* s;
	qaws_scalar area = 0, plate = 0, will = 0;
	fill_patch(p, 0);
	s = fn_surface(p);
	qaws_surface_functional_eval(NULL, s, QAWS_FUNCTIONAL_AREA, 0, NULL, &area, NULL, NULL);
	qaws_surface_functional_eval(NULL, s, QAWS_FUNCTIONAL_THIN_PLATE, 0, NULL, &plate, NULL, NULL);
	qaws_surface_functional_eval(NULL, s, QAWS_FUNCTIONAL_WILLMORE, 0, NULL, &will, NULL, NULL);
	TEST_ASSERT(approx_eq(area, 1), "unit square area");
	TEST_ASSERT(approx_eq(plate, 0) && approx_eq(will, 0), "a plane has no thin-plate or Willmore energy");
	qaws_surface_destroy(s);
}

static double surface_value(qaws_scalar const* p, qaws_surface_functional f)
{
	qaws_surface* s = fn_surface(p);
	qaws_scalar v = 0;
	qaws_surface_functional_eval(NULL, s, f, 4, NULL, &v, NULL, NULL);
	qaws_surface_destroy(s);
	return v;
}

static void surface_gradient(qaws_scalar const* p, qaws_surface_functional f, qaws_scalar* g)
{
	qaws_surface* s = fn_surface(p);
	qaws_field_view fv;
	qaws_diff_views v;
	memset(g, 0, sizeof(qaws_scalar) * 48);
	v = one_view(&fv, g, 16, 3);
	qaws_surface_functional_gradient(NULL, s, f, 4, &v, NULL);
	qaws_surface_destroy(s);
}

static void check_surface_functional(qaws_surface_functional f, char const* name)
{
	qaws_scalar p[48], d[48], g[48], hv[48], pp[48], pm[48], gp[48], gm[48];
	qaws_surface* s;
	qaws_field_view fv, fv2;
	qaws_diff_views vd, vh;
	qaws_scalar val, tan, tan2;
	double h = DIFF_FD_STEP, fd, gd;
	unsigned int i;
	int ok_hvp = 1;

	fill_patch(p, 0.35);
	diff_rand_fill(d, 48);
	s = fn_surface(p);
	vd = one_view(&fv, d, 16, 3);
	qaws_surface_functional_eval(NULL, s, f, 4, &vd, &val, &tan, &tan2);
	surface_gradient(p, f, g);
	for (i = 0; i < 48; i++)
	{
		pp[i] = (qaws_scalar)(p[i] + h * d[i]);
		pm[i] = (qaws_scalar)(p[i] - h * d[i]);
	}
	fd = (surface_value(pp, f) - surface_value(pm, f)) / (2 * h);
	gd = diff_dot(g, d, 48);
	printf("    %s: value %.6g, tangent %.9g, <grad, d> %.9g, fd %.9g\n", name, (double)val, (double)tan, gd, fd);
	TEST_ASSERT(diff_close(fd, gd, DIFF_TOL * 300), "surface functional gradient matches finite differences");
	TEST_ASSERT(diff_close(tan, gd, DIFF_TOL * 10), "surface functional tangent equals <gradient, direction>");

	memset(hv, 0, sizeof(hv));
	vh = one_view(&fv2, hv, 16, 3);
	TEST_ASSERT_STATUS(qaws_surface_functional_hvp(NULL, s, f, 4, &vd, &vh));
	surface_gradient(pp, f, gp);
	surface_gradient(pm, f, gm);
	for (i = 0; i < 48; i++)
		if (!diff_close((gp[i] - gm[i]) / (2 * h), hv[i], DIFF_TOL * 300))
			ok_hvp = 0;
	TEST_ASSERT(ok_hvp, "surface HVP matches finite differences of the gradient");
	TEST_ASSERT(diff_close(tan2, diff_dot(d, hv, 48), DIFF_TOL * 10), "surface second tangent equals <d, H d>");
	qaws_surface_destroy(s);
}

int test_56_diff_functionals_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 56: Integral functionals\n");
	diff_seed(5656);
	test_curve_values();
	check_curve_functional(QAWS_FUNCTIONAL_LENGTH, "length");
	check_curve_functional(QAWS_FUNCTIONAL_BENDING, "bending");
	check_curve_functional(QAWS_FUNCTIONAL_CURVATURE_SQUARED, "curvature squared");
	test_surface_values();
	check_surface_functional(QAWS_FUNCTIONAL_AREA, "area");
	check_surface_functional(QAWS_FUNCTIONAL_THIN_PLATE, "thin plate");
	check_surface_functional(QAWS_FUNCTIONAL_WILLMORE, "willmore");
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
