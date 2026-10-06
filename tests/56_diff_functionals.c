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
 *   - Mathematica ground truth (tests/reference/56_functionals.wls): value,
 *     tangent, second tangent, gradient and HVP of a cubic B-spline curve
 *     and of a bicubic B-spline surface, from symbolic derivatives of the
 *     integrands at 30 digits; the surface error is reported per cell count
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
/*  Mathematica reference (tests/reference/56_functionals.wls)        */
/* ------------------------------------------------------------------ */

#include "reference/56_functionals.h"

/* |a - b| / max(1, |a|, |b|), the measure diff_close bounds */
static double ref_err(double a, double b)
{
	double scale = 1.0;
	if (fabs(a) > scale) scale = fabs(a);
	if (fabs(b) > scale) scale = fabs(b);
	return fabs(a - b) / scale;
}

static double ref_max_err(qaws_scalar const* x, double const* ref, unsigned int n)
{
	double e = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
		if (ref_err(x[i], ref[i]) > e)
			e = ref_err(x[i], ref[i]);
	return e;
}

static qaws_scalar const g_fn_ref_p[18] = { 0, 0, 0, 1, 2, 0.5f, 2.5f, 2, -0.5f, 3.5f, 0.5f, 1, 5, 1, 0, 6, 3, 0.5f };
static double const g_fn_ref_v[18] = { 0.2, -0.3, 0.1, -0.25, 0.2, 0.3, 0.1, 0.1, -0.2, 0.3, -0.2, 0.25,
	-0.2, 0.3, -0.1, 0.1, -0.25, 0.2 };

/* Errors of value/tangent/tangent2, gradient and HVP with `quadrature`
   points per piece. vals = { value, tangent, tangent2 } along g_fn_ref_v. */
static void ref_curve_errs(qaws_curve_functional f, unsigned int quadrature, double const* vals, double const* grad,
	double const* hvp, double* e)
{
	qaws_scalar dir[18], g[18], hv[18], got[3];
	qaws_field_view fd, fg, fh;
	qaws_diff_views vd, vg, vh;
	qaws_curve* c = fn_curve(g_fn_ref_p);
	unsigned int i;
	for (i = 0; i < 18; i++)
		dir[i] = (qaws_scalar)g_fn_ref_v[i];
	memset(g, 0, sizeof(g));
	memset(hv, 0, sizeof(hv));
	vd = one_view(&fd, dir, 6, 3);
	vg = one_view(&fg, g, 6, 3);
	vh = one_view(&fh, hv, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_functional_eval(NULL, c, f, quadrature, &vd, &got[0], &got[1], &got[2]));
	TEST_ASSERT_STATUS(qaws_curve_functional_gradient(NULL, c, f, quadrature, &vg, NULL));
	TEST_ASSERT_STATUS(qaws_curve_functional_hvp(NULL, c, f, quadrature, &vd, &vh));
	e[0] = ref_max_err(got, vals, 3);
	e[1] = ref_max_err(g, grad, 18);
	e[2] = ref_max_err(hv, hvp, 18);
	qaws_curve_destroy(c);
}

/* The default rule (6 points x 8 pieces per span) is reported; the check
   uses the 8-point rule, the curvature squared integrand |C'|^-5 being the
   least smooth. */
static void test_reference_curve(qaws_curve_functional f, char const* name, double const* vals, double const* grad,
	double const* hvp)
{
	double tol = QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-11, e6[3], e8[3];
	char msg[160];
	ref_curve_errs(f, 0, vals, grad, hvp, e6);
	ref_curve_errs(f, 8, vals, grad, hvp, e8);
	printf("    reference %-18s value/tangents %.1e, gradient %.1e, HVP %.1e (8 points: %.1e, %.1e, %.1e)\n", name,
		e6[0], e6[1], e6[2], e8[0], e8[1], e8[2]);
	sprintf(msg, "reference (%s): value, tangent and second tangent", name);
	TEST_ASSERT(e8[0] <= tol, msg);
	sprintf(msg, "reference (%s): gradient", name);
	TEST_ASSERT(e8[1] <= tol, msg);
	sprintf(msg, "reference (%s): Hessian-vector product", name);
	TEST_ASSERT(e8[2] <= tol, msg);
}

/* Bicubic B-spline, 5 x 5 control points, knots { 0, 0, 0, 0, 1, 2, 2, 2, 2 }
   in both directions: four polynomial patches on [0, 2]^2. An even number
   of cells keeps the quadrature cells inside the patches. */
static qaws_surface* ref_surface(qaws_scalar const* p)
{
	static qaws_scalar const knots[9] = { 0, 0, 0, 0, 1, 2, 2, 2, 2 };
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)p;
	d.u_point_count = 5;
	d.v_point_count = 5;
	d.u_knots = knots;
	d.u_knot_count = 9;
	d.v_knots = knots;
	d.v_knot_count = 9;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

/* The exact rationals of the script, row-major P[i][j]. */
static void ref_surface_data(qaws_scalar* p, qaws_scalar* v)
{
	unsigned int i, j, c;
	for (i = 0; i < 5; i++)
		for (j = 0; j < 5; j++)
		{
			qaws_scalar* q = &p[(i * 5 + j) * 3];
			q[0] = (qaws_scalar)((10.0 * i + (double)((2 * i + 3 * j) % 5) - 2.0) / 20.0);
			q[1] = (qaws_scalar)((10.0 * j + (double)((3 * i + j) % 5) - 2.0) / 20.0);
			q[2] = (qaws_scalar)(((double)((3 * i + 5 * j) % 7) - 3.0) / 10.0);
			for (c = 0; c < 3; c++)
				v[(i * 5 + j) * 3 + c] = (qaws_scalar)(((double)((7 * i + 3 * j + 5 * c) % 11) - 5.0) / 10.0);
		}
}

/* Largest relative error over value, tangents, gradient and HVP with
   `cells` quadrature cells per direction. */
static double ref_surface_err(qaws_surface_functional f, unsigned int cells, double const* vals, double const* grad,
	double const* hvp, double* e_parts)
{
	qaws_scalar p[75], dir[75], g[75], hv[75], got[3];
	qaws_field_view fd, fg, fh;
	qaws_diff_views vd, vg, vh;
	qaws_surface* s;
	ref_surface_data(p, dir);
	s = ref_surface(p);
	memset(g, 0, sizeof(g));
	memset(hv, 0, sizeof(hv));
	vd = one_view(&fd, dir, 25, 3);
	vg = one_view(&fg, g, 25, 3);
	vh = one_view(&fh, hv, 25, 3);
	qaws_surface_functional_eval(NULL, s, f, cells, &vd, &got[0], &got[1], &got[2]);
	qaws_surface_functional_gradient(NULL, s, f, cells, &vg, NULL);
	qaws_surface_functional_hvp(NULL, s, f, cells, &vd, &vh);
	qaws_surface_destroy(s);
	e_parts[0] = ref_max_err(got, vals, 3);
	e_parts[1] = ref_max_err(g, grad, 75);
	e_parts[2] = ref_max_err(hv, hvp, 75);
	return e_parts[0] > e_parts[1] ? (e_parts[0] > e_parts[2] ? e_parts[0] : e_parts[2])
		: (e_parts[1] > e_parts[2] ? e_parts[1] : e_parts[2]);
}

static void test_reference_surface(qaws_surface_functional f, char const* name, double const* vals, double const* grad,
	double const* hvp)
{
	/* 4 x 4 Gauss points per cell: the error falls as cells^-8 */
	static unsigned int const cells[5] = { 8, 16, 32, 64, 128 };
	double e[5], parts[3], tol = QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-11;
	char msg[160];
	int k;
	for (k = 0; k < 5; k++)
		e[k] = ref_surface_err(f, cells[k], vals, grad, hvp, parts);
	printf("    reference %-10s cells 8: %.1e, 16: %.1e, 32: %.1e, 64: %.1e, 128: %.1e\n", name, e[0], e[1], e[2], e[3], e[4]);
	printf("      with 128 cells: value/tangents %.1e, gradient %.1e, HVP %.1e\n", parts[0], parts[1], parts[2]);
	sprintf(msg, "reference (%s): value, tangents, gradient and HVP with 128 x 128 cells", name);
	TEST_ASSERT(e[4] <= tol, msg);
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


/* ------------------------------------------------------------------ */
/*  Knots as parameters of curve functionals                          */
/*                                                                    */
/*  Unclamped knots so every knot moves freely; the quadratic case     */
/*  has a bending integrand that jumps at knots, so the moving span    */
/*  boundaries contribute.                                             */
/* ------------------------------------------------------------------ */

static qaws_curve* knot_curve(unsigned int degree, qaws_scalar const* cps, qaws_scalar const* knots)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = degree;
	d.control_points = cps;
	d.control_point_count = 6;
	d.knots = knots;
	d.knot_count = 6 + degree + 1;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static double knot_value(unsigned int degree, qaws_scalar const* cps, qaws_scalar const* knots,
	qaws_scalar const* dk, double h, qaws_curve_functional f)
{
	qaws_scalar k[16];
	qaws_scalar v = 0;
	unsigned int i;
	qaws_curve* c;
	for (i = 0; i < 6 + degree + 1; i++)
		k[i] = (qaws_scalar)(knots[i] + h * dk[i]);
	c = knot_curve(degree, cps, k);
	qaws_curve_functional_eval(NULL, c, f, 0, NULL, &v, NULL, NULL);
	qaws_curve_destroy(c);
	return v;
}

static void check_knot_functional(unsigned int degree, qaws_curve_functional f, char const* name)
{
	static qaws_scalar const k3[10] = { 0, 0.4f, 0.9f, 1.5f, 2.1f, 2.6f, 3.2f, 3.6f, 4.1f, 4.5f };
	static qaws_scalar const k2[9] = { 0, 0.5f, 1.1f, 1.6f, 2.3f, 2.9f, 3.4f, 4.0f, 4.4f };
	qaws_scalar const* knots = degree == 3 ? k3 : k2;
	unsigned int kc = 6 + degree + 1, i;
	qaws_scalar cps[18], dk[16], gk[16], gc[18], v = 0, tan = 0, tan2 = 0, tp = 0, tm = 0;
	qaws_field_view fv[2];
	qaws_diff_views dir, grad;
	qaws_curve* c;
	double h = DIFF_FD_STEP * 0.5, fd, fd2, gdot;
	char label[96];

	for (i = 0; i < 6; i++)
	{
		cps[3 * i] = (qaws_scalar)i;
		cps[3 * i + 1] = (qaws_scalar)(0.8 * sin(1.3 * i));
		cps[3 * i + 2] = (qaws_scalar)(0.3 * cos(0.9 * i));
	}
	diff_rand_fill(dk, kc);
	for (i = 0; i < kc; i++)
		dk[i] *= (qaws_scalar)0.3;
	c = knot_curve(degree, cps, knots);

	fv[0] = qaws_field_view_make(QAWS_FIELD_KNOTS, dk, kc, 1);
	dir.fields = fv;
	dir.field_count = 1;
	dir.children = NULL;
	dir.child_count = 0;
	TEST_ASSERT_STATUS(qaws_curve_functional_eval(NULL, c, f, 0, &dir, &v, &tan, &tan2));

	fd = (knot_value(degree, cps, knots, dk, h, f) - knot_value(degree, cps, knots, dk, -h, f)) / (2 * h);
	{
		/* second derivative from tangents at shifted knots */
		qaws_scalar ks[16];
		qaws_curve* cs;
		for (i = 0; i < kc; i++) ks[i] = (qaws_scalar)(knots[i] + h * dk[i]);
		cs = knot_curve(degree, cps, ks);
		qaws_curve_functional_eval(NULL, cs, f, 0, &dir, NULL, &tp, NULL);
		qaws_curve_destroy(cs);
		for (i = 0; i < kc; i++) ks[i] = (qaws_scalar)(knots[i] - h * dk[i]);
		cs = knot_curve(degree, cps, ks);
		qaws_curve_functional_eval(NULL, cs, f, 0, &dir, NULL, &tm, NULL);
		qaws_curve_destroy(cs);
		fd2 = ((double)tp - tm) / (2 * h);
	}

	memset(gk, 0, sizeof(gk));
	memset(gc, 0, sizeof(gc));
	fv[0] = qaws_field_view_make(QAWS_FIELD_KNOTS, gk, kc, 1);
	fv[1] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, gc, 6, 3);
	grad.fields = fv;
	grad.field_count = 2;
	grad.children = NULL;
	grad.child_count = 0;
	TEST_ASSERT_STATUS(qaws_curve_functional_gradient(NULL, c, f, 0, &grad, NULL));
	gdot = diff_dot(gk, dk, kc);

	printf("    degree %u %-18s tangent %.8f fd %.8f gradient %.8f | tangent2 %.6f fd %.6f\n",
		degree, name, (double)tan, fd, gdot, (double)tan2, fd2);
	sprintf(label, "knot tangent of %s matches finite differences", name);
	TEST_ASSERT(diff_close(tan, fd, DIFF_TOL * 50), label);
	sprintf(label, "knot gradient of %s matches its tangent", name);
	TEST_ASSERT(diff_close(tan, gdot, DIFF_TOL * 10), label);
	sprintf(label, "knot tangent2 of %s matches finite differences", name);
	TEST_ASSERT(diff_close(tan2, fd2, DIFF_TOL * 100), label);
	qaws_curve_destroy(c);
}

static void test_knot_functionals(void)
{
	printf("  knots as parameters\n");
	check_knot_functional(3, QAWS_FUNCTIONAL_LENGTH, "length");
	check_knot_functional(3, QAWS_FUNCTIONAL_BENDING, "bending");
	check_knot_functional(3, QAWS_FUNCTIONAL_CURVATURE_SQUARED, "curvature squared");
	check_knot_functional(2, QAWS_FUNCTIONAL_LENGTH, "length");
	check_knot_functional(2, QAWS_FUNCTIONAL_BENDING, "bending");
}

/* Surface knots: interior knots change the integrand, the end knots of the
   domain move the whole cell grid and scale its weights. */
#define SK_N 5
#define SK_KN (SK_N + 4)

static qaws_surface* sk_surface(qaws_scalar const* cps, qaws_scalar const* uk, qaws_scalar const* vk)
{
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = SK_N;
	d.v_point_count = SK_N;
	d.u_knots = uk;
	d.u_knot_count = SK_KN;
	d.v_knots = vk;
	d.v_knot_count = SK_KN;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

static void sk_shift(qaws_scalar const* k, qaws_scalar const* dk, double h, qaws_scalar* out)
{
	unsigned int i;
	for (i = 0; i < SK_KN; i++)
		out[i] = (qaws_scalar)(k[i] + h * dk[i]);
}

static void check_surface_knot_functional(qaws_surface_functional f, char const* name)
{
	static qaws_scalar const uk[SK_KN] = { 0, 0, 0, 0, 1.3f, 3, 3, 3, 3 };
	static qaws_scalar const vk[SK_KN] = { 0, 0, 0, 0, 1.1f, 2.5f, 2.5f, 2.5f, 2.5f };
	qaws_scalar cps[SK_N * SK_N * 3], duk[SK_KN], dvk[SK_KN], ukp[SK_KN], vkp[SK_KN], ukm[SK_KN], vkm[SK_KN];
	qaws_scalar gu[SK_KN], gv[SK_KN], value, t1, t2, vp, vm, tp, tm;
	qaws_field_view fv[2], fg[2];
	qaws_diff_views vd, vg;
	qaws_surface *s, *sp, *sm;
	double h = QAWS_SCALAR_IS_FLOAT ? 1e-3 : 1e-5, tol = QAWS_SCALAR_IS_FLOAT ? 3e-2 : 1e-6;
	unsigned int i, j;
	char msg[128];
	for (i = 0; i < SK_N; i++)
		for (j = 0; j < SK_N; j++)
		{
			cps[3 * (i * SK_N + j) + 0] = (qaws_scalar)(0.75 * j);
			cps[3 * (i * SK_N + j) + 1] = (qaws_scalar)(0.6 * i);
			cps[3 * (i * SK_N + j) + 2] = (qaws_scalar)(0.3 * sin(1.3 * i + 0.7 * j));
		}
	diff_rand_fill(duk, SK_KN);
	diff_rand_fill(dvk, SK_KN);
	for (i = 0; i < SK_KN; i++)
	{
		duk[i] *= (qaws_scalar)0.2;
		dvk[i] *= (qaws_scalar)0.2;
	}
	fv[0] = qaws_field_view_make(QAWS_FIELD_U_KNOTS, duk, SK_KN, 1);
	fv[1] = qaws_field_view_make(QAWS_FIELD_V_KNOTS, dvk, SK_KN, 1);
	vd.fields = fv; vd.field_count = 2; vd.children = NULL; vd.child_count = 0;
	memset(gu, 0, sizeof(gu));
	memset(gv, 0, sizeof(gv));
	fg[0] = qaws_field_view_make(QAWS_FIELD_U_KNOTS, gu, SK_KN, 1);
	fg[1] = qaws_field_view_make(QAWS_FIELD_V_KNOTS, gv, SK_KN, 1);
	vg.fields = fg; vg.field_count = 2; vg.children = NULL; vg.child_count = 0;

	s = sk_surface(cps, uk, vk);
	sk_shift(uk, duk, h, ukp); sk_shift(vk, dvk, h, vkp);
	sk_shift(uk, duk, -h, ukm); sk_shift(vk, dvk, -h, vkm);
	sp = sk_surface(cps, ukp, vkp);
	sm = sk_surface(cps, ukm, vkm);
	TEST_ASSERT_STATUS(qaws_surface_functional_eval(NULL, s, f, 6, &vd, &value, &t1, &t2));
	qaws_surface_functional_eval(NULL, sp, f, 6, &vd, &vp, &tp, NULL);
	qaws_surface_functional_eval(NULL, sm, f, 6, &vd, &vm, &tm, NULL);
	TEST_ASSERT_STATUS(qaws_surface_functional_gradient(NULL, s, f, 6, &vg, NULL));
	printf("    %s: tangent %.10g fd %.10g, tangent2 %.10g fd %.10g\n", name, t1, (vp - vm) / (2 * h), t2, (tp - tm) / (2 * h));
	sprintf(msg, "surface %s: knot tangent matches finite differences", name);
	TEST_ASSERT(diff_close(t1, (vp - vm) / (2 * h), tol), msg);
	sprintf(msg, "surface %s: knot second tangent matches differences of tangents", name);
	TEST_ASSERT(diff_close(t2, (tp - tm) / (2 * h), 10 * tol), msg);
	sprintf(msg, "surface %s: <knot gradient, direction> = tangent", name);
	TEST_ASSERT(diff_close(diff_dot(gu, duk, SK_KN) + diff_dot(gv, dvk, SK_KN), t1, QAWS_SCALAR_IS_FLOAT ? 2e-3 : 1e-10), msg);
	qaws_surface_destroy(s);
	qaws_surface_destroy(sp);
	qaws_surface_destroy(sm);
}

static void test_surface_knot_functionals(void)
{
	printf("  surface knots as parameters\n");
	check_surface_knot_functional(QAWS_FUNCTIONAL_AREA, "area");
	check_surface_knot_functional(QAWS_FUNCTIONAL_THIN_PLATE, "thin plate");
	check_surface_knot_functional(QAWS_FUNCTIONAL_WILLMORE, "willmore");
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
	test_knot_functionals();
	test_surface_knot_functionals();
	test_reference_curve(QAWS_FUNCTIONAL_LENGTH, "length", ref_curve_length, ref_curve_grad_length, ref_curve_hvp_length);
	test_reference_curve(QAWS_FUNCTIONAL_BENDING, "bending", ref_curve_bending, ref_curve_grad_bending,
		ref_curve_hvp_bending);
	test_reference_curve(QAWS_FUNCTIONAL_CURVATURE_SQUARED, "curvature squared", ref_curve_curv2, ref_curve_grad_curv2,
		ref_curve_hvp_curv2);
	test_reference_surface(QAWS_FUNCTIONAL_AREA, "area", ref_surface_area, ref_surface_grad_area, ref_surface_hvp_area);
	test_reference_surface(QAWS_FUNCTIONAL_THIN_PLATE, "thin plate", ref_surface_plate, ref_surface_grad_plate,
		ref_surface_hvp_plate);
	test_reference_surface(QAWS_FUNCTIONAL_WILLMORE, "willmore", ref_surface_willmore, ref_surface_grad_willmore,
		ref_surface_hvp_willmore);
	test_surface_values();
	check_surface_functional(QAWS_FUNCTIONAL_AREA, "area");
	check_surface_functional(QAWS_FUNCTIONAL_THIN_PLATE, "thin plate");
	check_surface_functional(QAWS_FUNCTIONAL_WILLMORE, "willmore");
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
