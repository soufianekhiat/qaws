/*
 * Test 153: parameterizations of B-spline fitting (qaws_bspline_fit_desc.parameterization)
 *
 *   - interpolating fits (as many control points as data points) through
 *     unevenly spaced points: each passes through the data at the
 *     parameters of its parameterization (uniform steps, chord lengths,
 *     their square roots), and not at another's;
 *   - DEFAULT is chord lengths, the same control points to the bit;
 *   - an unknown parameterization refused.
 *   (The derivatives of the centripetal and uniform fits: test 55.)
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define PZ_M 9
/* float: eps |P| / h ~ 6e-5 per difference, times the collocation condition */
#define PZ_F32_TANGENT 2e-3

static double const pz_x[PZ_M] = { 0, 0.2, 0.3, 1.5, 3.0, 3.1, 3.3, 5.0, 5.4 };

static void pz_points(qaws_vec2* p)
{
	unsigned int i;
	for (i = 0; i < PZ_M; i++)
	{
		p[i].x = (qaws_scalar)pz_x[i];
		p[i].y = (qaws_scalar)sin(1.3 * pz_x[i]);
	}
}

/* the parameters a parameterization gives the points */
static void pz_params(qaws_vec2 const* p, double alpha, double* t)
{
	double total = 0;
	unsigned int i;
	t[0] = 0;
	for (i = 1; i < PZ_M; i++)
	{
		double l = hypot((double)p[i].x - p[i - 1].x, (double)p[i].y - p[i - 1].y);
		total += alpha == 0 ? 1 : pow(l, alpha);
		t[i] = total;
	}
	for (i = 1; i < PZ_M; i++)
		t[i] /= total;
}

static qaws_curve* pz_fit(qaws_vec2 const* p, qaws_parameterization param)
{
	qaws_bspline_fit_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.data_points = p;
	d.data_point_count = PZ_M;
	d.degree = 3;
	d.control_point_count = PZ_M;
	d.parameterization = param;
	return qaws_curve_fit_bspline(&d, &c) == QAWS_STATUS_OK ? c : NULL;
}

/* the largest distance from the curve at parameters t to the points */
static double pz_miss(qaws_curve const* c, qaws_vec2 const* p, double const* t)
{
	double worst = 0;
	unsigned int i;
	for (i = 0; i < PZ_M; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)t[i], QAWS_EVAL_FLAG_POSITION, &e);
		worst = fmax(worst, hypot((double)e.position.x - p[i].x, (double)e.position.y - p[i].y));
	}
	return worst;
}

static void test_parameterizations(void)
{
	static char const* names[3] = { "uniform", "chordal", "centripetal" };
	qaws_parameterization params[3] = { QAWS_PARAMETERIZATION_UNIFORM, QAWS_PARAMETERIZATION_CHORDAL, QAWS_PARAMETERIZATION_CENTRIPETAL };
	double alpha[3] = { 0, 1, 0.5 }, t[3][PZ_M];
	qaws_vec2 p[PZ_M];
	unsigned int k;
	char msg[200];
	pz_points(p);
	for (k = 0; k < 3; k++)
		pz_params(p, alpha[k], t[k]);
	for (k = 0; k < 3; k++)
	{
		qaws_curve* c = pz_fit(p, params[k]);
		double own = c ? pz_miss(c, p, t[k]) : 1, other = c ? pz_miss(c, p, t[(k + 1) % 3]) : 0;
		snprintf(msg, sizeof msg, "%s: through the points at its parameters (%.1e), not at %s ones (%.2f)", names[k], own, names[(k + 1) % 3], other);
		printf("  %s\n", msg);
		/* (an interpolating fit solves the collocation system: its rounding is the system's condition, a few tens here, times
		   eps times the data's size, 5) */
		TEST_ASSERT(c && own <= fmax(1e-10, 64 * TEST_EPS * 5) && other > 1e-3, msg);
		qaws_curve_destroy(c);
	}
}

static void test_default_and_refusal(void)
{
	qaws_vec2 p[PZ_M];
	qaws_curve *a, *b;
	qaws_scalar ca[PZ_M * 2], cb[PZ_M * 2];
	unsigned int got;
	qaws_bspline_fit_desc d;
	qaws_curve* c = NULL;
	pz_points(p);
	a = pz_fit(p, QAWS_PARAMETERIZATION_DEFAULT);
	b = pz_fit(p, QAWS_PARAMETERIZATION_CHORDAL);
	qaws_curve_read_field(a, QAWS_FIELD_CONTROL_POINTS, ca, PZ_M * 2, &got);
	qaws_curve_read_field(b, QAWS_FIELD_CONTROL_POINTS, cb, PZ_M * 2, &got);
	TEST_ASSERT(memcmp(ca, cb, sizeof ca) == 0, "DEFAULT is chord lengths, to the bit");
	qaws_curve_destroy(a);
	qaws_curve_destroy(b);
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.data_points = p; d.data_point_count = PZ_M; d.degree = 3; d.control_point_count = PZ_M;
	d.parameterization = (qaws_parameterization)7;
	TEST_ASSERT(qaws_curve_fit_bspline(&d, &c) == QAWS_STATUS_INVALID_ARGUMENT, "an unknown parameterization refused");
}

/* control points of the interpolating fit of the points moved by h along dir */
static void pz_fit_at(qaws_vec2 const* p, double const* dir, double h, qaws_parameterization param, qaws_scalar* cps)
{
	qaws_vec2 q[PZ_M];
	qaws_curve* c;
	unsigned int i, got;
	for (i = 0; i < PZ_M; i++)
	{
		q[i].x = (qaws_scalar)(p[i].x + h * dir[2 * i]);
		q[i].y = (qaws_scalar)(p[i].y + h * dir[2 * i + 1]);
	}
	c = pz_fit(q, param);
	qaws_curve_read_field(c, QAWS_FIELD_CONTROL_POINTS, cps, PZ_M * 2, &got);
	qaws_curve_destroy(c);
}

/* the interpolating fit's map (its knots average p parameters) against differences of the fit */
static void test_interpolation_map(void)
{
	qaws_vec2 p[PZ_M];
	double dir[PZ_M * 2], worst = 0;
	qaws_scalar tq[PZ_M * 2], tk[PZ_M + 4], qp[PZ_M * 2], qm[PZ_M * 2], qp2[PZ_M * 2], qm2[PZ_M * 2];
	qaws_field_view in_f, out_f[2];
	qaws_diff_views in_v, out_v;
	qaws_diff_views const* ins[1];
	qaws_diff_views* outs[1];
	qaws_bspline_fit_desc d;
	qaws_curve* c = NULL;
	qaws_diff_map* map = NULL;
	qaws_scalar vdir[PZ_M * 2];
	unsigned int i, k;
	char msg[200];
	pz_points(p);
	for (i = 0; i < PZ_M * 2; i++)
		vdir[i] = (qaws_scalar)(dir[i] = sin(1.7 * i + 0.3));
	for (k = 0; k < 2; k++)
	{
		qaws_parameterization param = k ? QAWS_PARAMETERIZATION_CENTRIPETAL : QAWS_PARAMETERIZATION_CHORDAL;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_2D; d.data_points = p; d.data_point_count = PZ_M; d.degree = 3; d.control_point_count = PZ_M;
		d.parameterization = param;
		TEST_ASSERT(qaws_curve_fit_bspline_diff(&d, &c, &map) == QAWS_STATUS_OK, "interpolating fit with its map");
		in_f = qaws_field_view_make(QAWS_FIELD_POINTS, vdir, PZ_M, 2);
		in_v.fields = &in_f; in_v.field_count = 1; in_v.children = NULL; in_v.child_count = 0;
		out_f[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, tq, PZ_M, 2);
		out_f[1] = qaws_field_view_make(QAWS_FIELD_KNOTS, tk, PZ_M + 4, 1);
		out_v.fields = out_f; out_v.field_count = 2; out_v.children = NULL; out_v.child_count = 0;
		ins[0] = &in_v;
		outs[0] = &out_v;
		TEST_ASSERT(qaws_diff_map_tangent(map, NULL, ins, 1, outs, 1) == QAWS_STATUS_OK, "its tangent");
		/* double: central differences with h = 1e-6. float: rounding eps |P| / h rules, so fourth-order differences
		   (truncation h^4) with h = 0.01, a tenth of the shortest chord so the centripetal square roots stay smooth */
		if (sizeof(qaws_scalar) == sizeof(double))
		{
			pz_fit_at(p, dir, 1e-6, param, qp);
			pz_fit_at(p, dir, -1e-6, param, qm);
			for (i = 0; i < PZ_M * 2; i++)
				worst = fmax(worst, fabs(tq[i] - ((double)qp[i] - qm[i]) / 2e-6) / (1 + fabs((double)tq[i])));
		}
		else
		{
			double const h = 0.01;
			pz_fit_at(p, dir, h, param, qp);
			pz_fit_at(p, dir, -h, param, qm);
			pz_fit_at(p, dir, 2 * h, param, qp2);
			pz_fit_at(p, dir, -2 * h, param, qm2);
			for (i = 0; i < PZ_M * 2; i++)
				worst = fmax(worst, fabs(tq[i] - (8 * ((double)qp[i] - qm[i]) - ((double)qp2[i] - qm2[i])) / (12 * h)) / (1 + fabs((double)tq[i])));
		}
		qaws_diff_map_destroy(map);
		qaws_curve_destroy(c);
	}
	snprintf(msg, sizeof msg, "interpolating fits, chordal and centripetal: the map's tangent = differences of the fit (%.1e)", worst);
	printf("  %s\n", msg);
	TEST_ASSERT(worst <= (sizeof(qaws_scalar) == sizeof(double) ? 1e-6 : PZ_F32_TANGENT), msg);
}

int test_153_parameterization_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 153: parameterizations of B-spline fitting\n");
	test_parameterizations();
	test_default_and_refusal();
	test_interpolation_map();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
