/*
 * Test 140: the evaluation's degree limits
 *
 *   - each constructor refuses the degrees its evaluation cannot hold
 *     (INVALID_DEGREE, where they used to overrun a stack buffer): Bezier
 *     curves above 16, B-spline and NURBS curves and surfaces above 31,
 *     rational Bezier curves above 47, Bezier surfaces above 63;
 *   - the largest accepted degrees evaluate with every derivative asked: a
 *     B-spline of degree 20 with D3 (its buffer held degree 16 before) and
 *     of degree 31, a Bezier of degree 16, a rational Bezier of degree 47;
 *     each derivative against central differences of the one below.
 */

#include "test_common.h"
#include "qaws_bezier.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_rational_bezier.h"
#include "qaws_surface_bezier.h"
#include "qaws_surface_bspline.h"
#include "qaws_surface_nurbs.h"
#include "qaws_eval.h"
#include "qaws_curve.h"
#include "qaws_surface.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* n points of a wavy space curve, deterministic */
static qaws_scalar* tl_points(unsigned int n)
{
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double s = (double)i / (n > 1 ? n - 1 : 1);
		p[3 * i] = (qaws_scalar)s;
		p[3 * i + 1] = (qaws_scalar)(0.3 * sin(7.1 * s + 0.4 * i));
		p[3 * i + 2] = (qaws_scalar)(0.2 * cos(5.3 * s + 0.7 * i));
	}
	return p;
}

/* the largest relative error of D1, D2, D3 at t against central differences of the derivative below */
static double tl_derivatives(qaws_curve const* c, double t, double h)
{
	qaws_eval_result_3d a, m, p, m2, p2;
	unsigned int const all = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;
	double worst = 0;
	unsigned int k, d;
	/* the parameters as the scalar holds them, so the differences divide by the step actually taken */
	t = (qaws_scalar)t;
	h = (qaws_scalar)(t + h) - t;
	if (qaws_curve_evaluate_3d(c, (qaws_scalar)t, all, &a) != QAWS_STATUS_OK || qaws_curve_evaluate_3d(c, (qaws_scalar)(t - h), all, &m) != QAWS_STATUS_OK
		|| qaws_curve_evaluate_3d(c, (qaws_scalar)(t + h), all, &p) != QAWS_STATUS_OK || qaws_curve_evaluate_3d(c, (qaws_scalar)(t - 2 * h), all, &m2) != QAWS_STATUS_OK
		|| qaws_curve_evaluate_3d(c, (qaws_scalar)(t + 2 * h), all, &p2) != QAWS_STATUS_OK || !(a.valid_flags & QAWS_EVAL_FLAG_D3))
		return HUGE_VAL;
	for (d = 1; d <= 3; d++)
	{
		qaws_vec3 const *an = d == 1 ? &a.d1 : d == 2 ? &a.d2 : &a.d3, *lo = d == 1 ? &m.position : d == 2 ? &m.d1 : &m.d2,
			*hi = d == 1 ? &p.position : d == 2 ? &p.d1 : &p.d2, *lo2 = d == 1 ? &m2.position : d == 2 ? &m2.d1 : &m2.d2,
			*hi2 = d == 1 ? &p2.position : d == 2 ? &p2.d1 : &p2.d2;
		double num[3], ana[3], scale = 0, err = 0;
		/* fourth-order central differences: truncation h^4 f^(5) / 30, so a degree-n curve takes h near 1 / (2n) in float */
		num[0] = (8 * ((double)hi->x - lo->x) - ((double)hi2->x - lo2->x)) / (12 * h);
		num[1] = (8 * ((double)hi->y - lo->y) - ((double)hi2->y - lo2->y)) / (12 * h);
		num[2] = (8 * ((double)hi->z - lo->z) - ((double)hi2->z - lo2->z)) / (12 * h);
		ana[0] = an->x;
		ana[1] = an->y;
		ana[2] = an->z;
		for (k = 0; k < 3; k++)
		{
			scale = fmax(scale, fabs(ana[k]));
			err = fmax(err, fabs(num[k] - ana[k]));
		}
		worst = fmax(worst, err / fmax(scale, 1));
	}
	return worst;
}

static void test_refusals(void)
{
	qaws_scalar* p = tl_points(70);
	qaws_scalar w[70];
	qaws_vec3 grid[4];
	qaws_curve* c = NULL;
	qaws_surface* s = NULL;
	qaws_status st[8];
	unsigned int i;
	char msg[400];
	for (i = 0; i < 70; i++)
		w[i] = 1;
	memset(grid, 0, sizeof grid);
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 17;
		d.control_points = p;
		d.control_point_count = 18;
		st[0] = qaws_curve_create_bezier(&d, &c);
	}
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 32;
		d.control_points = p;
		d.control_point_count = 40;
		st[1] = qaws_curve_create_bspline(&d, &c);
	}
	{
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 32;
		d.control_points = p;
		d.control_point_count = 40;
		d.weights = w;
		d.weight_count = 40;
		st[2] = qaws_curve_create_nurbs(&d, &c);
	}
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 48;
		d.control_points = p;
		d.control_point_count = 49;
		d.weights = w;
		d.weight_count = 49;
		st[3] = qaws_curve_create_rational_bezier(&d, &c);
	}
	{
		qaws_surface_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 64;
		d.v_degree = 1;
		d.control_points = grid;
		d.u_point_count = 65;
		d.v_point_count = 2;
		st[4] = qaws_surface_create_bezier(&d, &s);
	}
	{
		qaws_surface_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 1;
		d.v_degree = 32;
		d.control_points = grid;
		d.u_point_count = 2;
		d.v_point_count = 33;
		st[5] = qaws_surface_create_bspline(&d, &s);
	}
	{
		qaws_surface_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 32;
		d.v_degree = 1;
		d.control_points = grid;
		d.u_point_count = 33;
		d.v_point_count = 2;
		d.weights = w;
		st[6] = qaws_surface_create_nurbs(&d, &s);
	}
	sprintf(msg, "refused past their limits: Bezier 17 %s, B-spline 32 %s, NURBS 32 %s, rational Bezier 48 %s, Bezier surface 64 %s, "
		"B-spline surface 32 %s, NURBS surface 32 %s", qaws_status_to_string(st[0]), qaws_status_to_string(st[1]), qaws_status_to_string(st[2]),
		qaws_status_to_string(st[3]), qaws_status_to_string(st[4]), qaws_status_to_string(st[5]), qaws_status_to_string(st[6]));
	for (i = 0; i < 7 && st[i] == QAWS_STATUS_INVALID_DEGREE; i++) {}
	TEST_ASSERT(i == 7 && !c && !s, msg);
	free(p);
}

static void test_largest(void)
{
	qaws_scalar* p = tl_points(60);
	qaws_scalar w[60];
	qaws_curve *b16 = NULL, *s20 = NULL, *s31 = NULL, *r47 = NULL;
	qaws_status st[4];
	/* relative to the derivative's size. double: steps far below truncation, rounding eps / h ~ 1e-11 per difference.
	   float: a step of 1 / (2 degree) puts truncation (1/2)^4 / 30 ~ 2e-3 and rounding (eps times the degree, over n h) ~ 1e-5 */
	double e[4], tol = sizeof(qaws_scalar) == sizeof(double) ? 1e-4 : 1e-2;
	int const f32 = sizeof(qaws_scalar) != sizeof(double);
	unsigned int i;
	char msg[400];
	for (i = 0; i < 60; i++)
		w[i] = (qaws_scalar)(1 + 0.3 * sin(1.7 * i));
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 16;
		d.control_points = p;
		d.control_point_count = 17;
		st[0] = qaws_curve_create_bezier(&d, &b16);
	}
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 20;
		d.control_points = p;
		d.control_point_count = 26;
		st[1] = qaws_curve_create_bspline(&d, &s20);
		d.degree = 31;
		d.control_point_count = 34;
		st[2] = qaws_curve_create_bspline(&d, &s31);
	}
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 47;
		d.control_points = p;
		d.control_point_count = 48;
		d.weights = w;
		d.weight_count = 48;
		st[3] = qaws_curve_create_rational_bezier(&d, &r47);
	}
	e[0] = b16 ? tl_derivatives(b16, 0.37, f32 ? 0.5 / 16 : 1e-5) : HUGE_VAL;
	e[1] = s20 ? tl_derivatives(s20, 0.41, f32 ? 0.5 / 20 : 1e-5) : HUGE_VAL;
	e[2] = s31 ? tl_derivatives(s31, 0.53, f32 ? 0.5 / 31 : 1e-5) : HUGE_VAL;
	/* (degree 47: the third derivative is 47 * 46 * 45 times the points; a wider step keeps the differences above rounding) */
	e[3] = r47 ? tl_derivatives(r47, 0.29, f32 ? 0.5 / 47 : 1e-4) : HUGE_VAL;
	sprintf(msg, "largest degrees evaluated with D3: Bezier 16 %s (derivatives off by %.1e), B-spline 20 %s (%.1e), B-spline 31 %s (%.1e), "
		"rational Bezier 47 %s (%.1e), tolerance %.0e", qaws_status_to_string(st[0]), e[0], qaws_status_to_string(st[1]), e[1],
		qaws_status_to_string(st[2]), e[2], qaws_status_to_string(st[3]), e[3], tol);
	TEST_ASSERT(st[0] == QAWS_STATUS_OK && st[1] == QAWS_STATUS_OK && st[2] == QAWS_STATUS_OK && st[3] == QAWS_STATUS_OK && e[0] < tol && e[1] < tol
		&& e[2] < tol && e[3] < tol, msg);
	qaws_curve_destroy(b16);
	qaws_curve_destroy(s20);
	qaws_curve_destroy(s31);
	qaws_curve_destroy(r47);
	free(p);
}

int test_140_eval_limits_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 140: the evaluation's degree limits\n");
	test_refusals();
	test_largest();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
