/*
 * Test 54: Differentiable solves (implicit function theorem)
 *
 * Closest point on curves and surfaces:
 *   - tangents with respect to the query and the shape parameters match
 *     finite differences of re-solved problems
 *   - adjoint identity
 *   - reports: valid interior solutions, ambiguous competing solutions,
 *     active endpoints / boundaries held fixed
 *   - Mathematica ground truth (tests/reference/54_implicit.wls): solution,
 *     tangent and adjoint of curve and surface closest points, a curve
 *     pair, a plane crossing, an extremum, an inflection and a
 *     surface-curve piercing, by the implicit function theorem at 30 digits
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


/* ------------------------------------------------------------------ */
/*  Roots seeded by the discrete finders                              */
/* ------------------------------------------------------------------ */

static void root_curve_a(qaws_scalar* p, unsigned int dim)
{
	unsigned int n;
	for (n = 0; n < 6; n++)
	{
		p[dim * n + 0] = (qaws_scalar)n;
		p[dim * n + 1] = (qaws_scalar)(0.8 * sin(1.1 * n));
		if (dim == 3)
			p[dim * n + 2] = (qaws_scalar)(0.4 * cos(0.7 * n));
	}
}

static void root_shift(qaws_scalar* out, qaws_scalar const* p, qaws_scalar const* dir, unsigned int n, double h)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		out[i] = (qaws_scalar)(p[i] + h * dir[i]);
}

static double v3w(qaws_vec3 a, qaws_vec3 w)
{
	return (double)a.x * w.x + (double)a.y * w.y + (double)a.z * w.z;
}

static double pair_dot(qaws_curve_pair_point const* a, qaws_curve_pair_point const* w)
{
	return (double)a->t_a * w->t_a + (double)a->t_b * w->t_b + v3w(a->position_a, w->position_a) +
		v3w(a->position_b, w->position_b) + (double)a->distance * w->distance;
}

static double pair_diff(qaws_curve_pair_point const* p, qaws_curve_pair_point const* m, qaws_curve_pair_point const* w, double h)
{
	return (pair_dot(p, w) - pair_dot(m, w)) / (2 * h);
}

static void check_pair(unsigned int dim)
{
	qaws_scalar ca[18], cb[18], da[18], db[18], ba[18], bb[18], ap[18], am[18], bp[18], bm[18];
	qaws_curve *A, *B;
	qaws_scalar seeds[4][2];
	unsigned int ns = 0, i, n, nc = 6 * dim;
	int ok_fd = 1, ok_adj = 1, ok_valid = 1;
	double h = DIFF_FD_STEP;

	diff_seed(900 + dim);
	root_curve_a(ca, dim);
	for (n = 0; n < 6; n++)
	{
		cb[dim * n + 0] = (qaws_scalar)(0.4 + 0.9 * n);
		cb[dim * n + 1] = (qaws_scalar)(1.3 - 0.55 * n);
		if (dim == 3)
			cb[dim * n + 2] = (qaws_scalar)(0.9 - 0.1 * n);
	}
	A = imp_curve(ca, dim);
	B = imp_curve(cb, dim);
	if (dim == 2)
	{
		qaws_intersection_2d hits[4];
		qaws_curve_find_intersections_2d(A, B, hits, 4, &ns);
		for (i = 0; i < ns; i++)
		{
			seeds[i][0] = hits[i].parameter_a;
			seeds[i][1] = hits[i].parameter_b;
		}
	}
	else
	{
		/* closest approach seeds from a coarse grid */
		double best = 1e30;
		unsigned int a, b;
		for (a = 0; a <= 40; a++)
			for (b = 0; b <= 40; b++)
			{
				qaws_eval_result_3d ra, rb;
				double d;
				qaws_curve_evaluate_3d(A, (qaws_scalar)(3.0 * a / 40), QAWS_EVAL_FLAG_POSITION, &ra);
				qaws_curve_evaluate_3d(B, (qaws_scalar)(3.0 * b / 40), QAWS_EVAL_FLAG_POSITION, &rb);
				d = v3w(qaws_v3_sub(ra.position, rb.position), qaws_v3_sub(ra.position, rb.position));
				if (d < best)
				{
					best = d;
					seeds[0][0] = (qaws_scalar)(3.0 * a / 40);
					seeds[0][1] = (qaws_scalar)(3.0 * b / 40);
				}
			}
		ns = 1;
	}
	printf("    pair %uD: %u solutions\n", dim, ns);
	TEST_ASSERT(ns >= 1, "curve pair has a solution to differentiate");

	for (i = 0; i < ns; i++)
	{
		qaws_diff_context ctx;
		qaws_diff_report report;
		qaws_field_view fa, fb;
		qaws_diff_views va, vb;
		qaws_curve_pair_point val, tan, vp, vm, w;
		qaws_curve *Ap, *Am, *Bp, *Bm;
		double lhs, rhs;

		diff_rand_fill(da, nc);
		diff_rand_fill(db, nc);
		w.t_a = diff_rand();
		w.t_b = diff_rand();
		w.position_a = diff_rand_vec3();
		w.position_b = diff_rand_vec3();
		w.distance = dim == 3 ? diff_rand() : 0;   /* |r| is not differentiable at an intersection */
		if (dim == 2)
			w.position_a.z = w.position_b.z = 0;
		va = imp_views(&fa, da, 6, dim);
		vb = imp_views(&fb, db, 6, dim);
		qaws_diff_context_init(&ctx);
		qaws_diff_report_reset(&report);
		ctx.report = &report;
		TEST_ASSERT_STATUS(qaws_curve_pair_point_tangent(&ctx, A, B, seeds[i][0], seeds[i][1], &va, &vb, &val, &tan));
		if (report.validity != QAWS_DIFF_VALID)
			ok_valid = 0;

		root_shift(ap, ca, da, nc, h);
		root_shift(am, ca, da, nc, -h);
		root_shift(bp, cb, db, nc, h);
		root_shift(bm, cb, db, nc, -h);
		Ap = imp_curve(ap, dim);
		Am = imp_curve(am, dim);
		Bp = imp_curve(bp, dim);
		Bm = imp_curve(bm, dim);
		qaws_curve_pair_point_tangent(NULL, Ap, Bp, val.t_a, val.t_b, NULL, NULL, &vp, NULL);
		qaws_curve_pair_point_tangent(NULL, Am, Bm, val.t_a, val.t_b, NULL, NULL, &vm, NULL);
		lhs = pair_dot(&tan, &w);
		rhs = pair_diff(&vp, &vm, &w, h);
		if (!diff_close(lhs, rhs, DIFF_TOL * 50))
			ok_fd = 0;
		printf("      t_a %.4f t_b %.4f distance %.4g: tangent %.8f fd %.8f\n", (double)val.t_a, (double)val.t_b,
			(double)val.distance, lhs, rhs);

		memset(ba, 0, sizeof(ba));
		memset(bb, 0, sizeof(bb));
		va = imp_views(&fa, ba, 6, dim);
		vb = imp_views(&fb, bb, 6, dim);
		TEST_ASSERT_STATUS(qaws_curve_pair_point_adjoint(NULL, A, B, seeds[i][0], seeds[i][1], &w, &va, &vb));
		rhs = diff_dot(ba, da, nc) + diff_dot(bb, db, nc);
		if (!diff_close(lhs, rhs, DIFF_TOL))
			ok_adj = 0;
		qaws_curve_destroy(Ap);
		qaws_curve_destroy(Am);
		qaws_curve_destroy(Bp);
		qaws_curve_destroy(Bm);
	}
	TEST_ASSERT(ok_valid, "pair solutions report valid");
	TEST_ASSERT(ok_fd, "pair tangent matches finite differences of re-solved problems");
	TEST_ASSERT(ok_adj, "pair adjoint identity");
	qaws_curve_destroy(A);
	qaws_curve_destroy(B);
}

static void check_plane(void)
{
	qaws_scalar c[18], d[18], bar[18], cp[18], cm[18], seeds[8];
	qaws_vec3 hits_pos[8];
	qaws_curve* C;
	qaws_plane plane;
	unsigned int ns = 0, i;
	int ok_fd = 1, ok_adj = 1;
	double h = DIFF_FD_STEP;

	diff_seed(950);
	root_curve_a(c, 3);
	C = imp_curve(c, 3);
	plane.point = qaws_v3((qaws_scalar)2.2, (qaws_scalar)0.3, (qaws_scalar)0.1);
	plane.normal = qaws_v3((qaws_scalar)1, (qaws_scalar)0.3, (qaws_scalar)-0.2);
	qaws_curve_find_plane_intersections(C, &plane, seeds, hits_pos, 8, &ns);
	printf("    plane: %u crossings\n", ns);
	TEST_ASSERT(ns >= 1, "plane crossing found");
	for (i = 0; i < ns; i++)
	{
		qaws_field_view fv;
		qaws_diff_views v;
		qaws_plane pd, pp, pm, pbar;
		qaws_curve_plane_point val, tan, vp, vm, w;
		qaws_curve *Cp, *Cm;
		double lhs, rhs;

		diff_rand_fill(d, 18);
		pd.point = diff_rand_vec3();
		pd.normal = diff_rand_vec3();
		w.t = diff_rand();
		w.position = diff_rand_vec3();
		v = imp_views(&fv, d, 6, 3);
		TEST_ASSERT_STATUS(qaws_curve_plane_point_tangent(NULL, C, &plane, seeds[i], &pd, &v, &val, &tan));
		root_shift(cp, c, d, 18, h);
		root_shift(cm, c, d, 18, -h);
		Cp = imp_curve(cp, 3);
		Cm = imp_curve(cm, 3);
		pp.point = qaws_v3_axpy(plane.point, pd.point, (qaws_scalar)h);
		pp.normal = qaws_v3_axpy(plane.normal, pd.normal, (qaws_scalar)h);
		pm.point = qaws_v3_axpy(plane.point, pd.point, (qaws_scalar)-h);
		pm.normal = qaws_v3_axpy(plane.normal, pd.normal, (qaws_scalar)-h);
		qaws_curve_plane_point_tangent(NULL, Cp, &pp, val.t, NULL, NULL, &vp, NULL);
		qaws_curve_plane_point_tangent(NULL, Cm, &pm, val.t, NULL, NULL, &vm, NULL);
		lhs = tan.t * w.t + v3w(tan.position, w.position);
		rhs = (((double)vp.t - vm.t) * w.t + v3w(qaws_v3_sub(vp.position, vm.position), w.position)) / (2 * h);
		if (!diff_close(lhs, rhs, DIFF_TOL * 50))
			ok_fd = 0;
		printf("      t %.4f: tangent %.8f fd %.8f\n", (double)val.t, lhs, rhs);

		memset(bar, 0, sizeof(bar));
		memset(&pbar, 0, sizeof(pbar));
		v = imp_views(&fv, bar, 6, 3);
		TEST_ASSERT_STATUS(qaws_curve_plane_point_adjoint(NULL, C, &plane, seeds[i], &w, &v, &pbar));
		rhs = diff_dot(bar, d, 18) + v3w(pbar.point, pd.point) + v3w(pbar.normal, pd.normal);
		if (!diff_close(lhs, rhs, DIFF_TOL))
			ok_adj = 0;
		qaws_curve_destroy(Cp);
		qaws_curve_destroy(Cm);
	}
	TEST_ASSERT(ok_fd, "plane crossing tangent matches finite differences");
	TEST_ASSERT(ok_adj, "plane crossing adjoint identity");
	qaws_curve_destroy(C);
}

static void check_extremum_inflection(void)
{
	qaws_scalar c[12], d[12], bar[12], cp[12], cm[12], seeds[8];
	qaws_curve* C;
	unsigned int ns = 0, i;
	int ok_fd = 1, ok_adj = 1;
	double h = DIFF_FD_STEP;

	diff_seed(970);
	root_curve_a(c, 2);
	C = imp_curve(c, 2);

	qaws_curve_find_extrema(C, 1, seeds, 8, &ns);
	printf("    extrema in y: %u\n", ns);
	TEST_ASSERT(ns >= 1, "extremum found");
	for (i = 0; i < ns; i++)
	{
		qaws_field_view fv;
		qaws_diff_views v;
		qaws_vec3 e = qaws_v3(0, 1, 0), ed = diff_rand_vec3(), ebar = qaws_v3_zero();
		qaws_curve_extremum val, tan, vp, vm, w;
		qaws_curve *Cp, *Cm;
		double lhs, rhs;

		ed.z = 0;
		diff_rand_fill(d, 12);
		w.t = diff_rand();
		w.position = diff_rand_vec3();
		w.position.z = 0;
		w.value = diff_rand();
		v = imp_views(&fv, d, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_extremum_tangent(NULL, C, e, seeds[i], &ed, &v, &val, &tan));
		root_shift(cp, c, d, 12, h);
		root_shift(cm, c, d, 12, -h);
		Cp = imp_curve(cp, 2);
		Cm = imp_curve(cm, 2);
		qaws_curve_extremum_tangent(NULL, Cp, qaws_v3_axpy(e, ed, (qaws_scalar)h), val.t, NULL, NULL, &vp, NULL);
		qaws_curve_extremum_tangent(NULL, Cm, qaws_v3_axpy(e, ed, (qaws_scalar)-h), val.t, NULL, NULL, &vm, NULL);
		lhs = tan.t * w.t + v3w(tan.position, w.position) + tan.value * w.value;
		rhs = (((double)vp.t - vm.t) * w.t + v3w(qaws_v3_sub(vp.position, vm.position), w.position) +
			((double)vp.value - vm.value) * w.value) / (2 * h);
		if (!diff_close(lhs, rhs, DIFF_TOL * 50))
			ok_fd = 0;
		printf("      extremum t %.4f: tangent %.8f fd %.8f\n", (double)val.t, lhs, rhs);
		memset(bar, 0, sizeof(bar));
		v = imp_views(&fv, bar, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_extremum_adjoint(NULL, C, e, seeds[i], &w, &v, &ebar));
		rhs = diff_dot(bar, d, 12) + v3w(ebar, ed);
		if (!diff_close(lhs, rhs, DIFF_TOL))
			ok_adj = 0;
		qaws_curve_destroy(Cp);
		qaws_curve_destroy(Cm);
	}
	TEST_ASSERT(ok_fd, "extremum tangent matches finite differences");
	TEST_ASSERT(ok_adj, "extremum adjoint identity");

	ok_fd = ok_adj = 1;
	qaws_curve_find_inflection_points(C, seeds, 8, &ns);
	printf("    inflections: %u\n", ns);
	TEST_ASSERT(ns >= 1, "inflection found");
	for (i = 0; i < ns; i++)
	{
		qaws_field_view fv;
		qaws_diff_views v;
		qaws_curve_inflection val, tan, vp, vm, w;
		qaws_curve *Cp, *Cm;
		double lhs, rhs;

		diff_rand_fill(d, 12);
		w.t = diff_rand();
		w.position = diff_rand_vec3();
		w.position.z = 0;
		v = imp_views(&fv, d, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_inflection_tangent(NULL, C, seeds[i], &v, &val, &tan));
		root_shift(cp, c, d, 12, h);
		root_shift(cm, c, d, 12, -h);
		Cp = imp_curve(cp, 2);
		Cm = imp_curve(cm, 2);
		qaws_curve_inflection_tangent(NULL, Cp, val.t, NULL, &vp, NULL);
		qaws_curve_inflection_tangent(NULL, Cm, val.t, NULL, &vm, NULL);
		lhs = tan.t * w.t + v3w(tan.position, w.position);
		rhs = (((double)vp.t - vm.t) * w.t + v3w(qaws_v3_sub(vp.position, vm.position), w.position)) / (2 * h);
		if (!diff_close(lhs, rhs, DIFF_TOL * 50))
			ok_fd = 0;
		printf("      inflection t %.4f: tangent %.8f fd %.8f\n", (double)val.t, lhs, rhs);
		memset(bar, 0, sizeof(bar));
		v = imp_views(&fv, bar, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_inflection_adjoint(NULL, C, seeds[i], &w, &v));
		rhs = diff_dot(bar, d, 12);
		if (!diff_close(lhs, rhs, DIFF_TOL))
			ok_adj = 0;
		qaws_curve_destroy(Cp);
		qaws_curve_destroy(Cm);
	}
	TEST_ASSERT(ok_fd, "inflection tangent matches finite differences");
	TEST_ASSERT(ok_adj, "inflection adjoint identity");
	qaws_curve_destroy(C);
}

static void check_pierce(void)
{
	qaws_scalar sc[48], sw[16], sd[48], sbar[48], sp[48], sm[48];
	qaws_scalar c[18], d[18], bar[18], cp[18], cm[18];
	qaws_surface *S;
	qaws_curve* C;
	qaws_surface_curve_intersection hits[4];
	unsigned int ns = 0, i, a, b, n;
	int ok_fd = 1, ok_adj = 1;
	double h = DIFF_FD_STEP;

	diff_seed(990);
	for (a = 0; a < 4; a++)
		for (b = 0; b < 4; b++)
		{
			qaws_scalar* p = &sc[(a * 4 + b) * 3];
			p[0] = (qaws_scalar)a;
			p[1] = (qaws_scalar)b;
			p[2] = (qaws_scalar)(0.3 * sin(1.3 * a + 0.4) * cos(0.9 * b));
			sw[a * 4 + b] = (qaws_scalar)1 + (qaws_scalar)0.2 * diff_rand();
		}
	for (n = 0; n < 6; n++)
	{
		c[3 * n + 0] = (qaws_scalar)(0.6 + 0.4 * n);
		c[3 * n + 1] = (qaws_scalar)(0.9 + 0.3 * n);
		c[3 * n + 2] = (qaws_scalar)(-0.8 + 0.35 * n);
	}
	S = imp_surface(sc, sw);
	C = imp_curve(c, 3);
	qaws_surface_find_curve_intersections(S, C, hits, 4, &ns);
	printf("    surface-curve: %u piercings\n", ns);
	TEST_ASSERT(ns >= 1, "piercing found");
	for (i = 0; i < ns; i++)
	{
		qaws_field_view fs, fc;
		qaws_diff_views vs, vc;
		qaws_surface_curve_point val, tan, vp, vm, w;
		qaws_surface *Sp, *Sm;
		qaws_curve *Cp, *Cm;
		double lhs, rhs;

		diff_rand_fill(sd, 48);
		diff_rand_fill(d, 18);
		w.u = diff_rand();
		w.v = diff_rand();
		w.t = diff_rand();
		w.position = diff_rand_vec3();
		vs = imp_views(&fs, sd, 16, 3);
		vc = imp_views(&fc, d, 6, 3);
		TEST_ASSERT_STATUS(qaws_surface_curve_point_tangent(NULL, S, C, hits[i].u, hits[i].v, hits[i].t, &vs, &vc, &val, &tan));
		root_shift(sp, sc, sd, 48, h);
		root_shift(sm, sc, sd, 48, -h);
		root_shift(cp, c, d, 18, h);
		root_shift(cm, c, d, 18, -h);
		Sp = imp_surface(sp, sw);
		Sm = imp_surface(sm, sw);
		Cp = imp_curve(cp, 3);
		Cm = imp_curve(cm, 3);
		qaws_surface_curve_point_tangent(NULL, Sp, Cp, val.u, val.v, val.t, NULL, NULL, &vp, NULL);
		qaws_surface_curve_point_tangent(NULL, Sm, Cm, val.u, val.v, val.t, NULL, NULL, &vm, NULL);
		lhs = tan.u * w.u + tan.v * w.v + tan.t * w.t + v3w(tan.position, w.position);
		rhs = (((double)vp.u - vm.u) * w.u + ((double)vp.v - vm.v) * w.v + ((double)vp.t - vm.t) * w.t +
			v3w(qaws_v3_sub(vp.position, vm.position), w.position)) / (2 * h);
		if (!diff_close(lhs, rhs, DIFF_TOL * 50))
			ok_fd = 0;
		printf("      (u %.4f, v %.4f, t %.4f): tangent %.8f fd %.8f\n", (double)val.u, (double)val.v, (double)val.t, lhs, rhs);
		memset(sbar, 0, sizeof(sbar));
		memset(bar, 0, sizeof(bar));
		vs = imp_views(&fs, sbar, 16, 3);
		vc = imp_views(&fc, bar, 6, 3);
		TEST_ASSERT_STATUS(qaws_surface_curve_point_adjoint(NULL, S, C, hits[i].u, hits[i].v, hits[i].t, &w, &vs, &vc));
		rhs = diff_dot(sbar, sd, 48) + diff_dot(bar, d, 18);
		if (!diff_close(lhs, rhs, DIFF_TOL))
			ok_adj = 0;
		qaws_surface_destroy(Sp);
		qaws_surface_destroy(Sm);
		qaws_curve_destroy(Cp);
		qaws_curve_destroy(Cm);
	}
	TEST_ASSERT(ok_fd, "piercing tangent matches finite differences");
	TEST_ASSERT(ok_adj, "piercing adjoint identity");
	qaws_surface_destroy(S);
	qaws_curve_destroy(C);
}

static void test_roots(void)
{
	printf("  roots seeded by the finders\n");
	check_pair(2);
	check_pair(3);
	check_plane();
	check_extremum_inflection();
	check_pierce();
}

/* ------------------------------------------------------------------ */
/*  Mathematica reference (tests/reference/54_implicit.wls)           */
/*                                                                    */
/*  theta lists the parameters of each problem in the script's order; */
/*  the direction and the output weight are the script's rationals.   */
/* ------------------------------------------------------------------ */

#include "reference/54_implicit.h"

static double const g_ref_pa[18] = { 0, 0, 0, 1, 0.7, 0.3, 2, 0.8, -0.2, 3, 0.1, 0.25, 4, -0.6, -0.1, 5, -0.5, 0.2 };

/* theta_dot[k] = (((7 k + 3) mod 11) - 5) / 10 */
static void ref_dir(qaws_scalar* out, unsigned int first, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		out[i] = (qaws_scalar)(((double)((7 * (first + i) + 3) % 11) - 5.0) / 10.0);
}

/* ybar[k] = (((5 k + 2) mod 9) - 4) / 10 */
static qaws_scalar ref_bar(unsigned int k)
{
	return (qaws_scalar)(((double)((5 * k + 2) % 9) - 4.0) / 10.0);
}

static double ref_err(double a, double b)
{
	double scale = 1.0;
	if (fabs(a) > scale) scale = fabs(a);
	if (fabs(b) > scale) scale = fabs(b);
	return fabs(a - b) / scale;
}

static double ref_max_err(double const* x, double const* ref, unsigned int n)
{
	double e = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
		if (ref_err(x[i], ref[i]) > e)
			e = ref_err(x[i], ref[i]);
	return e;
}

static void ref_check(char const* name, double const* val, double const* tan, double const* adj, unsigned int ny,
	unsigned int nt, double const* rval, double const* rtan, double const* radj)
{
	double tol = QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-11;
	double ev = ref_max_err(val, rval, ny), et = ref_max_err(tan, rtan, ny), ea = ref_max_err(adj, radj, nt);
	char msg[160];
	printf("    reference %-22s value %.1e, tangent %.1e, adjoint %.1e\n", name, ev, et, ea);
	sprintf(msg, "reference (%s): solution", name);
	TEST_ASSERT(ev <= tol, msg);
	sprintf(msg, "reference (%s): tangent", name);
	TEST_ASSERT(et <= tol, msg);
	sprintf(msg, "reference (%s): adjoint", name);
	TEST_ASSERT(ea <= tol, msg);
}

/* 2D: the inflection of the xy projection of g_ref_pa sits on a knot */
static double const g_ref_pa2[12] = { 0, 0, 1, 0.7, 2, 0.8, 3, 0.2, 4, -0.6, 5, -0.5 };

static void ref_curve_points(qaws_scalar* p, unsigned int dim)
{
	unsigned int n;
	for (n = 0; n < 6 * dim; n++)
		p[n] = (qaws_scalar)(dim == 3 ? g_ref_pa[n] : g_ref_pa2[n]);
}

/* The NURBS patch of the script (also tests/reference/52_geometry.wls). */
static void ref_patch(qaws_scalar* cps, qaws_scalar* w)
{
	unsigned int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			qaws_scalar* p = &cps[(i * 4 + j) * 3];
			p[0] = (qaws_scalar)(i + ((double)((i + 2 * j) % 3) - 1.0) / 10.0);
			p[1] = (qaws_scalar)(j + ((double)((2 * i + j) % 3) - 1.0) / 10.0);
			p[2] = (qaws_scalar)((16.0 * (i - 1.5) * (i - 1.5) - 10.0 * (j - 1.5) * (j - 1.5) + 4.0 * ((double)((i + j) % 3) - 1.0))
				/ 40.0);
			w[i * 4 + j] = (qaws_scalar)(1.0 + (2.0 * ((i + 3 * j) % 4) - 3.0) / 10.0);
		}
}

static qaws_vec3 ref_v3(qaws_scalar const* x)
{
	return qaws_v3(x[0], x[1], x[2]);
}

static void ref_put_v3(double* out, qaws_vec3 v)
{
	out[0] = v.x;
	out[1] = v.y;
	out[2] = v.z;
}

static void ref_curve_closest(unsigned int k)
{
	static double const qs[2][3] = { { 1.5, 1.2, 0.9 }, { 2.8, -0.2, -0.5 } };
	qaws_scalar cps[18], dir[21], bar[18];
	qaws_field_view fv;
	qaws_diff_views views;
	qaws_curve* c;
	qaws_curve_closest_point val, tan, abar;
	qaws_vec3 q = qaws_v3((qaws_scalar)qs[k][0], (qaws_scalar)qs[k][1], (qaws_scalar)qs[k][2]), qd, qbar = qaws_v3_zero();
	double y[5], yt[5], adj[21];
	unsigned int n;
	ref_curve_points(cps, 3);
	ref_dir(dir, 0, 21);
	qd = ref_v3(dir + 18);
	c = imp_curve(cps, 3);
	views = imp_views(&fv, dir, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_closest_point_tangent(NULL, c, q, &qd, &views, &val, &tan));
	y[0] = val.t; ref_put_v3(y + 1, val.position); y[4] = val.distance;
	yt[0] = tan.t; ref_put_v3(yt + 1, tan.position); yt[4] = tan.distance;
	abar.t = ref_bar(0);
	abar.position = qaws_v3(ref_bar(1), ref_bar(2), ref_bar(3));
	abar.distance = ref_bar(4);
	memset(bar, 0, sizeof(bar));
	views = imp_views(&fv, bar, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_closest_point_adjoint(NULL, c, q, &abar, &views, &qbar));
	for (n = 0; n < 18; n++)
		adj[n] = bar[n];
	ref_put_v3(adj + 18, qbar);
	ref_check(k ? "curve closest point 2" : "curve closest point 1", y, yt, adj, 5, 21,
		k ? ref_cc1_value : ref_cc0_value, k ? ref_cc1_tangent : ref_cc0_tangent, k ? ref_cc1_adjoint : ref_cc0_adjoint);
	qaws_curve_destroy(c);
}

static void ref_surface_closest(void)
{
	qaws_scalar cps[48], w[16], dir[67], bar[64];
	qaws_field_view fv[2];
	qaws_diff_views views;
	qaws_surface* s;
	qaws_surface_closest_point val, tan, abar;
	qaws_vec3 q = qaws_v3((qaws_scalar)1.5, (qaws_scalar)1.4, (qaws_scalar)1.2), qd, qbar = qaws_v3_zero();
	double y[6], yt[6], adj[67];
	unsigned int n;
	ref_patch(cps, w);
	ref_dir(dir, 0, 67);
	qd = ref_v3(dir + 64);
	s = imp_surface(cps, w);
	fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, 16, 3);
	fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, dir + 48, 16, 1);
	views.fields = fv;
	views.field_count = 2;
	views.children = NULL;
	views.child_count = 0;
	TEST_ASSERT_STATUS(qaws_surface_closest_point_tangent(NULL, s, q, &qd, &views, &val, &tan));
	y[0] = val.u; y[1] = val.v; ref_put_v3(y + 2, val.position); y[5] = val.distance;
	yt[0] = tan.u; yt[1] = tan.v; ref_put_v3(yt + 2, tan.position); yt[5] = tan.distance;
	abar.u = ref_bar(0);
	abar.v = ref_bar(1);
	abar.position = qaws_v3(ref_bar(2), ref_bar(3), ref_bar(4));
	abar.distance = ref_bar(5);
	memset(bar, 0, sizeof(bar));
	fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, bar, 16, 3);
	fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, bar + 48, 16, 1);
	TEST_ASSERT_STATUS(qaws_surface_closest_point_adjoint(NULL, s, q, &abar, &views, &qbar));
	for (n = 0; n < 64; n++)
		adj[n] = bar[n];
	ref_put_v3(adj + 64, qbar);
	ref_check("surface closest point", y, yt, adj, 6, 67, ref_sc_value, ref_sc_tangent, ref_sc_adjoint);
	qaws_surface_destroy(s);
}

static void ref_pair(void)
{
	qaws_scalar ca[18], cb[18], dir[36], ba[18], bb[18];
	qaws_field_view fa, fb;
	qaws_diff_views va, vb;
	qaws_curve *A, *B;
	qaws_curve_pair_point val, tan, w;
	qaws_scalar sa = (qaws_scalar)ref_pair_value[0], sb = (qaws_scalar)ref_pair_value[1];
	double y[9], yt[9], adj[36];
	unsigned int n;
	ref_curve_points(ca, 3);
	for (n = 0; n < 6; n++)
	{
		cb[3 * n] = (qaws_scalar)((4.0 + 9.0 * n) / 10.0);
		cb[3 * n + 1] = (qaws_scalar)((26.0 - 11.0 * n) / 20.0);
		cb[3 * n + 2] = (qaws_scalar)((9.0 - n) / 10.0);
	}
	ref_dir(dir, 0, 36);
	A = imp_curve(ca, 3);
	B = imp_curve(cb, 3);
	va = imp_views(&fa, dir, 6, 3);
	vb = imp_views(&fb, dir + 18, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_pair_point_tangent(NULL, A, B, sa, sb, &va, &vb, &val, &tan));
	y[0] = val.t_a; y[1] = val.t_b; ref_put_v3(y + 2, val.position_a); ref_put_v3(y + 5, val.position_b); y[8] = val.distance;
	yt[0] = tan.t_a; yt[1] = tan.t_b; ref_put_v3(yt + 2, tan.position_a); ref_put_v3(yt + 5, tan.position_b);
	yt[8] = tan.distance;
	w.t_a = ref_bar(0);
	w.t_b = ref_bar(1);
	w.position_a = qaws_v3(ref_bar(2), ref_bar(3), ref_bar(4));
	w.position_b = qaws_v3(ref_bar(5), ref_bar(6), ref_bar(7));
	w.distance = ref_bar(8);
	memset(ba, 0, sizeof(ba));
	memset(bb, 0, sizeof(bb));
	va = imp_views(&fa, ba, 6, 3);
	vb = imp_views(&fb, bb, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_pair_point_adjoint(NULL, A, B, sa, sb, &w, &va, &vb));
	for (n = 0; n < 18; n++)
	{
		adj[n] = ba[n];
		adj[18 + n] = bb[n];
	}
	ref_check("curve pair", y, yt, adj, 9, 36, ref_pair_value, ref_pair_tangent, ref_pair_adjoint);
	qaws_curve_destroy(A);
	qaws_curve_destroy(B);
}

static void ref_plane(void)
{
	qaws_scalar cps[18], dir[24], bar[18];
	qaws_field_view fv;
	qaws_diff_views views;
	qaws_curve* c;
	qaws_plane plane, pd, pbar;
	qaws_curve_plane_point val, tan, w;
	double y[4], yt[4], adj[24];
	unsigned int n;
	ref_curve_points(cps, 3);
	ref_dir(dir, 0, 24);
	c = imp_curve(cps, 3);
	plane.point = qaws_v3((qaws_scalar)2.2, (qaws_scalar)0.3, (qaws_scalar)0.1);
	plane.normal = qaws_v3(1, (qaws_scalar)0.3, (qaws_scalar)-0.2);
	pd.point = ref_v3(dir + 18);
	pd.normal = ref_v3(dir + 21);
	views = imp_views(&fv, dir, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_plane_point_tangent(NULL, c, &plane, (qaws_scalar)ref_plane_value[0], &pd, &views, &val, &tan));
	y[0] = val.t; ref_put_v3(y + 1, val.position);
	yt[0] = tan.t; ref_put_v3(yt + 1, tan.position);
	w.t = ref_bar(0);
	w.position = qaws_v3(ref_bar(1), ref_bar(2), ref_bar(3));
	memset(bar, 0, sizeof(bar));
	memset(&pbar, 0, sizeof(pbar));
	views = imp_views(&fv, bar, 6, 3);
	TEST_ASSERT_STATUS(qaws_curve_plane_point_adjoint(NULL, c, &plane, (qaws_scalar)ref_plane_value[0], &w, &views, &pbar));
	for (n = 0; n < 18; n++)
		adj[n] = bar[n];
	ref_put_v3(adj + 18, pbar.point);
	ref_put_v3(adj + 21, pbar.normal);
	ref_check("plane crossing", y, yt, adj, 4, 24, ref_plane_value, ref_plane_tangent, ref_plane_adjoint);
	qaws_curve_destroy(c);
}

static void ref_extremum_inflection(void)
{
	qaws_scalar cps[12], dir[15], bar[12];
	qaws_field_view fv;
	qaws_diff_views views;
	qaws_curve* c;
	double y[4], yt[4], adj[15];
	unsigned int n;
	ref_curve_points(cps, 2);
	ref_dir(dir, 0, 15);
	c = imp_curve(cps, 2);
	{
		/* extremum of the height along e = (0, 1, 0) */
		qaws_curve_extremum val, tan, w;
		qaws_vec3 e = qaws_v3(0, 1, 0), ed = ref_v3(dir + 12), ebar = qaws_v3_zero();
		views = imp_views(&fv, dir, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_extremum_tangent(NULL, c, e, (qaws_scalar)ref_ext_value[0], &ed, &views, &val, &tan));
		y[0] = val.t; y[1] = val.position.x; y[2] = val.position.y; y[3] = val.value;
		yt[0] = tan.t; yt[1] = tan.position.x; yt[2] = tan.position.y; yt[3] = tan.value;
		w.t = ref_bar(0);
		w.position = qaws_v3(ref_bar(1), ref_bar(2), 0);
		w.value = ref_bar(3);
		memset(bar, 0, sizeof(bar));
		views = imp_views(&fv, bar, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_extremum_adjoint(NULL, c, e, (qaws_scalar)ref_ext_value[0], &w, &views, &ebar));
		for (n = 0; n < 12; n++)
			adj[n] = bar[n];
		ref_put_v3(adj + 12, ebar);
		ref_check("extremum", y, yt, adj, 4, 15, ref_ext_value, ref_ext_tangent, ref_ext_adjoint);
	}
	{
		qaws_curve_inflection val, tan, w;
		views = imp_views(&fv, dir, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_inflection_tangent(NULL, c, (qaws_scalar)ref_infl_value[0], &views, &val, &tan));
		y[0] = val.t; y[1] = val.position.x; y[2] = val.position.y;
		yt[0] = tan.t; yt[1] = tan.position.x; yt[2] = tan.position.y;
		w.t = ref_bar(0);
		w.position = qaws_v3(ref_bar(1), ref_bar(2), 0);
		memset(bar, 0, sizeof(bar));
		views = imp_views(&fv, bar, 6, 2);
		TEST_ASSERT_STATUS(qaws_curve_inflection_adjoint(NULL, c, (qaws_scalar)ref_infl_value[0], &w, &views));
		for (n = 0; n < 12; n++)
			adj[n] = bar[n];
		ref_check("inflection", y, yt, adj, 3, 12, ref_infl_value, ref_infl_tangent, ref_infl_adjoint);
	}
	qaws_curve_destroy(c);
}

static void ref_pierce(void)
{
	qaws_scalar sc[48], sw[16], c[18], dir[66], sbar[48], bar[18];
	qaws_field_view fs, fc;
	qaws_diff_views vs, vc;
	qaws_surface* S;
	qaws_curve* C;
	qaws_surface_curve_point val, tan, w;
	qaws_scalar su = (qaws_scalar)ref_pierce_value[0], sv = (qaws_scalar)ref_pierce_value[1];
	qaws_scalar st = (qaws_scalar)ref_pierce_value[2];
	double y[6], yt[6], adj[66];
	unsigned int n;
	ref_patch(sc, sw);
	for (n = 0; n < 6; n++)
	{
		c[3 * n] = (qaws_scalar)((3.0 + 2.0 * n) / 5.0);
		c[3 * n + 1] = (qaws_scalar)((9.0 + 3.0 * n) / 10.0);
		c[3 * n + 2] = (qaws_scalar)((-16.0 + 7.0 * n) / 20.0);
	}
	ref_dir(dir, 0, 66);
	S = imp_surface(sc, sw);
	C = imp_curve(c, 3);
	vs = imp_views(&fs, dir, 16, 3);
	vc = imp_views(&fc, dir + 48, 6, 3);
	TEST_ASSERT_STATUS(qaws_surface_curve_point_tangent(NULL, S, C, su, sv, st, &vs, &vc, &val, &tan));
	y[0] = val.u; y[1] = val.v; y[2] = val.t; ref_put_v3(y + 3, val.position);
	yt[0] = tan.u; yt[1] = tan.v; yt[2] = tan.t; ref_put_v3(yt + 3, tan.position);
	w.u = ref_bar(0);
	w.v = ref_bar(1);
	w.t = ref_bar(2);
	w.position = qaws_v3(ref_bar(3), ref_bar(4), ref_bar(5));
	memset(sbar, 0, sizeof(sbar));
	memset(bar, 0, sizeof(bar));
	vs = imp_views(&fs, sbar, 16, 3);
	vc = imp_views(&fc, bar, 6, 3);
	TEST_ASSERT_STATUS(qaws_surface_curve_point_adjoint(NULL, S, C, su, sv, st, &w, &vs, &vc));
	for (n = 0; n < 48; n++)
		adj[n] = sbar[n];
	for (n = 0; n < 18; n++)
		adj[48 + n] = bar[n];
	ref_check("surface-curve piercing", y, yt, adj, 6, 66, ref_pierce_value, ref_pierce_tangent, ref_pierce_adjoint);
	qaws_surface_destroy(S);
	qaws_curve_destroy(C);
}

static void test_reference(void)
{
	printf("  Mathematica reference\n");
	ref_curve_closest(0);
	ref_curve_closest(1);
	ref_surface_closest();
	ref_pair();
	ref_plane();
	ref_extremum_inflection();
	ref_pierce();
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
	test_roots();
	test_reference();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
