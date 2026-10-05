/*
 * Test 54: Differentiable solves (implicit function theorem)
 *
 * Closest point on curves and surfaces:
 *   - tangents with respect to the query and the shape parameters match
 *     finite differences of re-solved problems
 *   - adjoint identity
 *   - reports: valid interior solutions, ambiguous competing solutions,
 *     active endpoints / boundaries held fixed
 */

#include "test_diff.h"

static qaws_scalar const g_imp_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_curve* imp_curve(qaws_scalar const* p, unsigned int dim)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)dim;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 6;
	d.knots = g_imp_knots;
	d.knot_count = 10;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static qaws_diff_views imp_views(qaws_field_view* fv, qaws_scalar* data, unsigned int count, unsigned int comps)
{
	qaws_diff_views v;
	*fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, data, count, comps);
	v.fields = fv;
	v.field_count = 1;
	v.children = NULL;
	v.child_count = 0;
	return v;
}

static void check_curve_closest(unsigned int dim)
{
	qaws_scalar cps[18], dir[18], bar[18], pp[18], pm[18];
	qaws_curve *c, *cp, *cm;
	unsigned int i, n, nc = 6 * dim;
	int ok_fd = 1, ok_adj = 1, ok_valid = 1;
	double h = DIFF_FD_STEP;

	diff_seed(600 + dim);
	for (n = 0; n < 6; n++)
	{
		cps[dim * n + 0] = (qaws_scalar)n;
		cps[dim * n + 1] = (qaws_scalar)(0.8 * sin(1.1 * n));
		if (dim == 3)
			cps[dim * n + 2] = (qaws_scalar)(0.4 * cos(0.7 * n));
	}
	c = imp_curve(cps, dim);

	for (i = 0; i < 6; i++)
	{
		qaws_diff_context ctx;
		qaws_diff_report report;
		qaws_field_view fv;
		qaws_diff_views views;
		qaws_vec3 q = qaws_v3((qaws_scalar)(0.6 + 0.7 * i), (qaws_scalar)(1.3 + 0.1 * i), dim == 3 ? (qaws_scalar)0.9 : 0);
		qaws_vec3 qd = diff_rand_vec3(), qbar = qaws_v3_zero();
		qaws_curve_closest_point val, tan, vp, vm, abar;
		double lhs, rhs;

		if (dim == 2)
			qd.z = 0;
		diff_rand_fill(dir, nc);
		views = imp_views(&fv, dir, 6, dim);
		qaws_diff_context_init(&ctx);
		qaws_diff_report_reset(&report);
		ctx.report = &report;
		TEST_ASSERT_STATUS(qaws_curve_closest_point_tangent(&ctx, c, q, &qd, &views, &val, &tan));
		if (report.validity != QAWS_DIFF_VALID)
			ok_valid = 0;

		for (n = 0; n < nc; n++)
		{
			pp[n] = (qaws_scalar)(cps[n] + h * dir[n]);
			pm[n] = (qaws_scalar)(cps[n] - h * dir[n]);
		}
		cp = imp_curve(pp, dim);
		cm = imp_curve(pm, dim);
		qaws_curve_closest_point_tangent(NULL, cp, qaws_v3_axpy(q, qd, (qaws_scalar)h), NULL, NULL, &vp, NULL);
		qaws_curve_closest_point_tangent(NULL, cm, qaws_v3_axpy(q, qd, (qaws_scalar)-h), NULL, NULL, &vm, NULL);
		if (!diff_close((vp.t - vm.t) / (2 * h), tan.t, DIFF_TOL * 300) ||
		    !diff_close((vp.position.y - vm.position.y) / (2 * h), tan.position.y, DIFF_TOL * 300) ||
		    !diff_close((vp.distance - vm.distance) / (2 * h), tan.distance, DIFF_TOL * 300))
		{
			printf("    closest %uD sample %u: fd t %.8g vs %.8g, d %.8g vs %.8g\n", dim, i,
				(vp.t - vm.t) / (2 * h), (double)tan.t, (vp.distance - vm.distance) / (2 * h), (double)tan.distance);
			ok_fd = 0;
		}
		qaws_curve_destroy(cp);
		qaws_curve_destroy(cm);

		abar.t = diff_rand();
		abar.position = diff_rand_vec3();
		if (dim == 2)
			abar.position.z = 0;
		abar.distance = diff_rand();
		lhs = (double)abar.t * tan.t + qaws_v3_dot(abar.position, tan.position) + (double)abar.distance * tan.distance;
		memset(bar, 0, sizeof(bar));
		views = imp_views(&fv, bar, 6, dim);
		TEST_ASSERT_STATUS(qaws_curve_closest_point_adjoint(NULL, c, q, &abar, &views, &qbar));
		rhs = diff_dot(bar, dir, nc) + qaws_v3_dot(qbar, qd);
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
			ok_adj = 0;
	}
	printf("    curve %uD closest point: tangent %s, adjoint %s\n", dim, ok_fd ? "ok" : "NO", ok_adj ? "ok" : "NO");
	TEST_ASSERT(ok_valid, "interior closest points are valid");
	TEST_ASSERT(ok_fd, "closest point tangent matches re-solved finite differences");
	TEST_ASSERT(ok_adj, "closest point adjoint identity (parameters and query)");
	qaws_curve_destroy(c);
}

static void test_curve_reports(void)
{
	/* Symmetric arch: the apex query below it has two equally close feet. */
	static qaws_scalar const arch[18] = { 0, 0, 0, 1, 2, 0, 2, 3, 0, 3, 3, 0, 4, 2, 0, 5, 0, 0 };
	qaws_curve* c = imp_curve(arch, 3);
	qaws_diff_context ctx;
	qaws_diff_report report;
	qaws_curve_closest_point val, tan;

	qaws_diff_context_init(&ctx);
	qaws_diff_report_reset(&report);
	ctx.report = &report;
	qaws_curve_closest_point_tangent(&ctx, c, qaws_v3((qaws_scalar)2.5, (qaws_scalar)0.4, 0), NULL, NULL, &val, NULL);
	TEST_ASSERT(report.validity == QAWS_DIFF_AMBIGUOUS, "competing closest points are reported as ambiguous");
	TEST_ASSERT(report.branch_gap < (qaws_scalar)0.05, "branch gap is small for a symmetric query");

	qaws_diff_report_reset(&report);
	{
		qaws_scalar dir[18];
		qaws_field_view fv;
		qaws_diff_views views;
		qaws_vec3 qd = qaws_v3(1, 1, 0);
		diff_seed(5);
		diff_rand_fill(dir, 18);
		views = imp_views(&fv, dir, 6, 3);
		qaws_curve_closest_point_tangent(&ctx, c, qaws_v3(-2, -1, 0), &qd, &views, &val, &tan);
	}
	TEST_ASSERT(approx_eq(val.t, 0), "query beyond the start projects onto the endpoint");
	TEST_ASSERT(report.validity == QAWS_DIFF_VALID_LOCALLY && (report.frozen_used & QAWS_FREEZE_ACTIVE_SET),
		"active endpoint is held fixed and reported");
	TEST_ASSERT(tan.t == 0, "parameter tangent is zero on the active endpoint");
	qaws_curve_destroy(c);
}

/* ------------------------------------------------------------------ */

static qaws_scalar const g_imp_sknots[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

static qaws_surface* imp_surface(qaws_scalar const* cps, qaws_scalar const* w)
{
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = 4;
	d.v_point_count = 4;
	d.weights = w;
	d.u_knots = g_imp_sknots;
	d.u_knot_count = 8;
	d.v_knots = g_imp_sknots;
	d.v_knot_count = 8;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

static void test_surface_closest(void)
{
	qaws_scalar cps[48], w[16], dir[64], bar[64], pp[48], pm[48], wp[16], wm[16];
	qaws_surface *s, *sp, *sm;
	unsigned int i, a, b, n;
	int ok_fd = 1, ok_adj = 1;
	double h = DIFF_FD_STEP;

	diff_seed(700);
	for (a = 0; a < 4; a++)
		for (b = 0; b < 4; b++)
		{
			qaws_scalar* p = &cps[(a * 4 + b) * 3];
			p[0] = (qaws_scalar)a;
			p[1] = (qaws_scalar)b;
			p[2] = (qaws_scalar)(0.3 * sin(1.3 * a + 0.4) * cos(0.9 * b));
			w[a * 4 + b] = (qaws_scalar)1 + (qaws_scalar)0.2 * diff_rand();
		}
	s = imp_surface(cps, w);

	for (i = 0; i < 5; i++)
	{
		qaws_field_view fv[2];
		qaws_diff_views views;
		qaws_vec3 q = qaws_v3((qaws_scalar)(0.7 + 0.4 * i), (qaws_scalar)(2.1 - 0.3 * i), (qaws_scalar)0.8);
		qaws_vec3 qd = diff_rand_vec3(), qbar = qaws_v3_zero();
		qaws_surface_closest_point val, tan, vp, vm, abar;
		double lhs, rhs;

		diff_rand_fill(dir, 64);
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, 16, 3);
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, dir + 48, 16, 1);
		views.fields = fv;
		views.field_count = 2;
		views.children = NULL;
		views.child_count = 0;
		TEST_ASSERT_STATUS(qaws_surface_closest_point_tangent(NULL, s, q, &qd, &views, &val, &tan));

		for (n = 0; n < 48; n++) { pp[n] = (qaws_scalar)(cps[n] + h * dir[n]); pm[n] = (qaws_scalar)(cps[n] - h * dir[n]); }
		for (n = 0; n < 16; n++) { wp[n] = (qaws_scalar)(w[n] + h * dir[48 + n]); wm[n] = (qaws_scalar)(w[n] - h * dir[48 + n]); }
		sp = imp_surface(pp, wp);
		sm = imp_surface(pm, wm);
		qaws_surface_closest_point_tangent(NULL, sp, qaws_v3_axpy(q, qd, (qaws_scalar)h), NULL, NULL, &vp, NULL);
		qaws_surface_closest_point_tangent(NULL, sm, qaws_v3_axpy(q, qd, (qaws_scalar)-h), NULL, NULL, &vm, NULL);
		if (!diff_close((vp.u - vm.u) / (2 * h), tan.u, DIFF_TOL * 300) ||
		    !diff_close((vp.v - vm.v) / (2 * h), tan.v, DIFF_TOL * 300) ||
		    !diff_close((vp.position.z - vm.position.z) / (2 * h), tan.position.z, DIFF_TOL * 300) ||
		    !diff_close((vp.distance - vm.distance) / (2 * h), tan.distance, DIFF_TOL * 300))
		{
			printf("    surface closest sample %u: fd u %.8g vs %.8g\n", i, (vp.u - vm.u) / (2 * h), (double)tan.u);
			ok_fd = 0;
		}
		qaws_surface_destroy(sp);
		qaws_surface_destroy(sm);

		abar.u = diff_rand();
		abar.v = diff_rand();
		abar.position = diff_rand_vec3();
		abar.distance = diff_rand();
		lhs = (double)abar.u * tan.u + (double)abar.v * tan.v + qaws_v3_dot(abar.position, tan.position) + (double)abar.distance * tan.distance;
		memset(bar, 0, sizeof(bar));
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, bar, 16, 3);
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, bar + 48, 16, 1);
		TEST_ASSERT_STATUS(qaws_surface_closest_point_adjoint(NULL, s, q, &abar, &views, &qbar));
		rhs = diff_dot(bar, dir, 64) + qaws_v3_dot(qbar, qd);
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
			ok_adj = 0;
	}
	printf("    surface closest point: tangent %s, adjoint %s\n", ok_fd ? "ok" : "NO", ok_adj ? "ok" : "NO");
	TEST_ASSERT(ok_fd, "surface closest point tangent matches re-solved finite differences");
	TEST_ASSERT(ok_adj, "surface closest point adjoint identity (control points, weights, query)");

	{
		qaws_diff_context ctx;
		qaws_diff_report report;
		qaws_surface_closest_point val, tan;
		qaws_vec3 qd = qaws_v3(1, 0, 0);
		qaws_diff_context_init(&ctx);
		qaws_diff_report_reset(&report);
		ctx.report = &report;
		qaws_surface_closest_point_tangent(&ctx, s, qaws_v3(-1, (qaws_scalar)1.5, 0), &qd, NULL, &val, &tan);
		TEST_ASSERT(approx_eq(val.u, 0), "query beyond the u = 0 edge lands on the boundary");
		TEST_ASSERT((report.frozen_used & QAWS_FREEZE_ACTIVE_SET) && tan.u == 0, "boundary coordinate held fixed");
	}
	qaws_surface_destroy(s);
}

int test_54_diff_implicit_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 54: Differentiable solves\n");
	check_curve_closest(3);
	check_curve_closest(2);
	test_curve_reports();
	test_surface_closest();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
