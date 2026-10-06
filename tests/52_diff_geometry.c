/*
 * Test 52: Differentiable local geometry
 *
 * Frenet frame, speed, curvature, torsion (2D and 3D) and surface normal,
 * fundamental forms, Gaussian / mean / principal curvatures:
 *   - values match the existing inspection functions
 *   - tangents through the full chain (control points -> jet -> quantity)
 *     match finite differences
 *   - adjoint identity through the full chain
 *   - second tangent matches finite differences of the tangent
 *   - validity reports at straight points and umbilics
 *   - Mathematica ground truth (tests/reference/52_geometry.wls): values,
 *     tangents, second tangents and adjoints through the chain for a 3D
 *     B-spline, a 2D Bezier and a NURBS patch (weights and (u, v) moving),
 *     from exact series expansions at 30 digits
 */

#include "test_diff.h"

#define GEOM_SAMPLES 6

/* ------------------------------------------------------------------ */
/*  Fixtures                                                          */
/* ------------------------------------------------------------------ */

static qaws_scalar const g_curve_cps[6 * 3] = {
	0, 0, 0,  1, 2, 0.5f,  2.5f, 1.5f, -1,  3, -1, 1,  4.5f, 0.5f, 2,  6, 2, 0 };
static qaws_scalar const g_curve_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_curve* make_curve_3d(qaws_scalar const* cps)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = 6;
	d.knots = g_curve_knots;
	d.knot_count = 10;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static qaws_curve* make_curve_2d(qaws_scalar const* cps)
{
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 4;
	d.control_points = cps;
	d.control_point_count = 5;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

static qaws_scalar const g_curve2_cps[5 * 2] = { 0, 0, 1, 2, 2.5f, -1, 4, 1.5f, 5, 0 };

/* NURBS patch with curvature in both directions. */
static qaws_scalar g_surf_cps[4 * 4 * 3];
static qaws_scalar g_surf_w[16];
static qaws_scalar const g_surf_knots[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

static void init_surface_net(void)
{
	unsigned int i, j;
	diff_seed(4242);
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			qaws_scalar x = (qaws_scalar)i, y = (qaws_scalar)j;
			qaws_scalar* p = &g_surf_cps[(i * 4 + j) * 3];
			p[0] = x + diff_rand() * (qaws_scalar)0.1;
			p[1] = y + diff_rand() * (qaws_scalar)0.1;
			p[2] = (qaws_scalar)0.4 * (x - (qaws_scalar)1.5) * (x - (qaws_scalar)1.5)
			     - (qaws_scalar)0.25 * (y - (qaws_scalar)1.5) * (y - (qaws_scalar)1.5) + diff_rand() * (qaws_scalar)0.1;
			g_surf_w[i * 4 + j] = (qaws_scalar)1 + (qaws_scalar)0.3 * diff_rand();
		}
}

static qaws_surface* make_surface(qaws_scalar const* cps, qaws_scalar const* w)
{
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = 4;
	d.v_point_count = 4;
	d.weights = w;
	d.u_knots = g_surf_knots;
	d.u_knot_count = 8;
	d.v_knots = g_surf_knots;
	d.v_knot_count = 8;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

static void one_view(qaws_field_view* v, qaws_diff_views* views, qaws_diff_field field,
	qaws_scalar* data, unsigned int count, unsigned int comps)
{
	*v = qaws_field_view_make(field, data, count, comps);
	views->fields = v;
	views->field_count = 1;
	views->children = NULL;
	views->child_count = 0;
}

static double curve3_dot(qaws_curve_geometry_3d const* a, qaws_curve_geometry_3d const* b)
{
	return (double)qaws_v3_dot(a->tangent, b->tangent) + qaws_v3_dot(a->normal, b->normal) +
	       qaws_v3_dot(a->binormal, b->binormal) + (double)a->speed * b->speed +
	       (double)a->curvature * b->curvature + (double)a->torsion * b->torsion;
}

static void rand_curve3(qaws_curve_geometry_3d* g)
{
	g->tangent = diff_rand_vec3();
	g->normal = diff_rand_vec3();
	g->binormal = diff_rand_vec3();
	g->speed = diff_rand();
	g->curvature = diff_rand();
	g->torsion = diff_rand();
}

/* Flatten a 3D geometry into 12 scalars for comparisons. */
static void flat_curve3(qaws_curve_geometry_3d const* g, double* out)
{
	out[0] = g->tangent.x; out[1] = g->tangent.y; out[2] = g->tangent.z;
	out[3] = g->normal.x; out[4] = g->normal.y; out[5] = g->normal.z;
	out[6] = g->binormal.x; out[7] = g->binormal.y; out[8] = g->binormal.z;
	out[9] = g->speed; out[10] = g->curvature; out[11] = g->torsion;
}

static void flat_surface(qaws_surface_geometry const* g, double* out)
{
	out[0] = g->normal.x; out[1] = g->normal.y; out[2] = g->normal.z;
	out[3] = g->E; out[4] = g->F; out[5] = g->G;
	out[6] = g->L; out[7] = g->M; out[8] = g->N;
	out[9] = g->gaussian; out[10] = g->mean; out[11] = g->kappa1; out[12] = g->kappa2;
}

static int flat_close(double const* a, double const* b, unsigned int n, double tol)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		if (!diff_close(a[i], b[i], tol))
			return 0;
	return 1;
}

static qaws_status curve3_geometry_at(qaws_curve const* c, qaws_scalar t, qaws_curve_geometry_3d* g)
{
	qaws_curve_jet_3d p, tg;
	qaws_status st = qaws_curve_eval_tangent_3d(NULL, c, t, 0, 0xF, NULL, &p, &tg);
	if (st != QAWS_STATUS_OK)
		return st;
	return qaws_curve_geometry_eval_3d(&p, NULL, NULL, g, NULL, NULL, NULL);
}

/* ------------------------------------------------------------------ */
/*  Curve 3D                                                          */
/* ------------------------------------------------------------------ */

static void test_curve_3d(void)
{
	qaws_curve* c = make_curve_3d(g_curve_cps);
	unsigned int i;
	int ok_values = 1, ok_fd = 1, ok_adj = 1, ok_t2 = 1;
	double h = DIFF_FD_STEP;
	diff_seed(11);

	for (i = 0; i < GEOM_SAMPLES; i++)
	{
		qaws_scalar t = (qaws_scalar)0.17 + (qaws_scalar)0.47 * (qaws_scalar)i;
		qaws_scalar tdot = diff_rand();
		qaws_scalar dir[18], pbar[18], cps_p[18], cps_m[18], tbar = 0;
		qaws_field_view v;
		qaws_diff_views views;
		qaws_curve_jet_3d p, tg, tt, ybar;
		qaws_curve_geometry_3d g, gt, gtt, gbar, gp, gm;
		qaws_diff_validity validity;
		qaws_curve *cp, *cm;
		unsigned int n;
		double a[12], b[12];

		/* Values against the inspection functions. */
		{
			qaws_scalar kappa, tau, speed;
			qaws_vec3 T, N, B;
			qaws_curve_eval_tangent_3d(NULL, c, t, 0, 0xF, NULL, &p, &tg);
			qaws_curve_geometry_eval_3d(&p, NULL, NULL, &g, NULL, NULL, &validity);
			qaws_curve_compute_curvature_3d(c, t, &kappa);
			qaws_curve_compute_torsion_3d(c, t, &tau);
			qaws_curve_compute_speed(c, t, &speed);
			qaws_curve_compute_frenet_frame_3d(c, t, &T, &N, &B);
			if (validity != QAWS_DIFF_VALID || !diff_close(g.curvature, kappa, 1e-4) ||
			    !diff_close(g.torsion, tau, 1e-4) || !diff_close(g.speed, speed, 1e-4) ||
			    !diff_close(qaws_v3_dot(g.tangent, T), 1, 1e-4) || !diff_close(qaws_v3_dot(g.normal, N), 1, 1e-4) ||
			    !diff_close(qaws_v3_dot(g.binormal, B), 1, 1e-4))
				ok_values = 0;
		}

		/* Chain tangent and second tangent. */
		diff_rand_fill(dir, 18);
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, dir, 6, 3);
		qaws_curve_eval_batch_tangent2_3d(NULL, c, &t, &tdot, 1, 0xF, &views, &p, &tg, &tt);
		qaws_curve_geometry_eval_3d(&p, &tg, &tt, &g, &gt, &gtt, &validity);

		for (n = 0; n < 18; n++)
		{
			cps_p[n] = (qaws_scalar)(g_curve_cps[n] + h * dir[n]);
			cps_m[n] = (qaws_scalar)(g_curve_cps[n] - h * dir[n]);
		}
		cp = make_curve_3d(cps_p);
		cm = make_curve_3d(cps_m);
		curve3_geometry_at(cp, (qaws_scalar)(t + h * tdot), &gp);
		curve3_geometry_at(cm, (qaws_scalar)(t - h * tdot), &gm);
		flat_curve3(&gp, a);
		flat_curve3(&gm, b);
		for (n = 0; n < 12; n++)
			a[n] = (a[n] - b[n]) / (2 * h);
		flat_curve3(&gt, b);
		if (!flat_close(a, b, 12, DIFF_TOL * 100))
			ok_fd = 0;

		/* Second tangent: finite difference of the chained tangent. */
		{
			qaws_curve_jet_3d pp, tgp, pm, tgm;
			qaws_curve_geometry_3d g0, gtp, gtm;
			qaws_curve_eval_tangent_3d(NULL, cp, (qaws_scalar)(t + h * tdot), tdot, 0xF, &views, &pp, &tgp);
			qaws_curve_eval_tangent_3d(NULL, cm, (qaws_scalar)(t - h * tdot), tdot, 0xF, &views, &pm, &tgm);
			qaws_curve_geometry_eval_3d(&pp, &tgp, NULL, &g0, &gtp, NULL, NULL);
			qaws_curve_geometry_eval_3d(&pm, &tgm, NULL, &g0, &gtm, NULL, NULL);
			flat_curve3(&gtp, a);
			flat_curve3(&gtm, b);
			for (n = 0; n < 12; n++)
				a[n] = (a[n] - b[n]) / (2 * h);
			flat_curve3(&gtt, b);
			if (!flat_close(a, b, 12, DIFF_TOL * 300))
				ok_t2 = 0;
		}
		qaws_curve_destroy(cp);
		qaws_curve_destroy(cm);

		/* Adjoint through the chain. */
		rand_curve3(&gbar);
		memset(&ybar, 0, sizeof(ybar));
		qaws_curve_geometry_adjoint_3d(&p, &gbar, &ybar, NULL);
		memset(pbar, 0, sizeof(pbar));
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, pbar, 6, 3);
		qaws_curve_eval_adjoint_3d(NULL, c, t, ybar.channels, &ybar, &views, &tbar);
		if (!diff_close(curve3_dot(&gbar, &gt), diff_dot(pbar, dir, 18) + (double)tbar * tdot, DIFF_TOL * 10))
			ok_adj = 0;
	}

	TEST_ASSERT(ok_values, "3D frame, curvature, torsion and speed match inspection");
	TEST_ASSERT(ok_fd, "3D geometry tangent matches finite differences through the chain");
	TEST_ASSERT(ok_t2, "3D geometry second tangent matches finite differences");
	TEST_ASSERT(ok_adj, "3D geometry adjoint identity through the chain");
	qaws_curve_destroy(c);
}

/* ------------------------------------------------------------------ */
/*  Curve 2D                                                          */
/* ------------------------------------------------------------------ */

static void test_curve_2d(void)
{
	qaws_curve* c = make_curve_2d(g_curve2_cps);
	unsigned int i;
	int ok_values = 1, ok_adj = 1;
	diff_seed(12);

	for (i = 0; i < GEOM_SAMPLES; i++)
	{
		qaws_scalar t = (qaws_scalar)0.08 + (qaws_scalar)0.16 * (qaws_scalar)i;
		qaws_scalar tdot = diff_rand(), tbar = 0, kappa;
		qaws_scalar dir[10], pbar[10];
		qaws_field_view v;
		qaws_diff_views views;
		qaws_curve_jet_2d p, tg, ybar;
		qaws_curve_geometry_2d g, gt, gbar;
		double lhs, rhs;

		diff_rand_fill(dir, 10);
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, dir, 5, 2);
		qaws_curve_eval_tangent_2d(NULL, c, t, tdot, 0xF, &views, &p, &tg);
		qaws_curve_geometry_eval_2d(&p, &tg, NULL, &g, &gt, NULL, NULL);
		qaws_curve_compute_curvature_2d(c, t, &kappa);
		if (!diff_close(g.curvature, kappa, 1e-4))
			ok_values = 0;

		gbar.tangent.x = diff_rand(); gbar.tangent.y = diff_rand();
		gbar.normal.x = diff_rand(); gbar.normal.y = diff_rand();
		gbar.speed = diff_rand(); gbar.curvature = diff_rand();
		lhs = (double)gbar.tangent.x * gt.tangent.x + (double)gbar.tangent.y * gt.tangent.y +
		      (double)gbar.normal.x * gt.normal.x + (double)gbar.normal.y * gt.normal.y +
		      (double)gbar.speed * gt.speed + (double)gbar.curvature * gt.curvature;

		memset(&ybar, 0, sizeof(ybar));
		qaws_curve_geometry_adjoint_2d(&p, &gbar, &ybar, NULL);
		memset(pbar, 0, sizeof(pbar));
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, pbar, 5, 2);
		qaws_curve_eval_adjoint_2d(NULL, c, t, ybar.channels, &ybar, &views, &tbar);
		rhs = diff_dot(pbar, dir, 10) + (double)tbar * tdot;
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
			ok_adj = 0;
	}
	TEST_ASSERT(ok_values, "2D signed curvature matches inspection");
	TEST_ASSERT(ok_adj, "2D geometry adjoint identity through the chain");
	qaws_curve_destroy(c);
}

static void test_curve_straight(void)
{
	static qaws_scalar const line[4 * 3] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 };
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	qaws_curve_jet_3d p, tg;
	qaws_curve_geometry_3d g;
	qaws_diff_validity validity = QAWS_DIFF_VALID;

	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = line;
	d.control_point_count = 4;
	qaws_curve_create_bezier(&d, &c);
	qaws_curve_eval_tangent_3d(NULL, c, (qaws_scalar)0.4, 0, 0xF, NULL, &p, &tg);
	qaws_curve_geometry_eval_3d(&p, NULL, NULL, &g, NULL, NULL, &validity);
	TEST_ASSERT(validity == QAWS_DIFF_ILL_CONDITIONED, "straight point reports ill-conditioned frame");
	TEST_ASSERT(g.normal.x == 0 && g.binormal.x == 0 && g.torsion == 0, "undefined frame returned as zero");
	TEST_ASSERT(approx_eq(g.curvature, 0) && g.speed > 0, "curvature value still reported");
	qaws_curve_destroy(c);
}

/* ------------------------------------------------------------------ */
/*  Surfaces                                                          */
/* ------------------------------------------------------------------ */

static qaws_status surface_geometry_at(qaws_surface const* s, qaws_scalar u, qaws_scalar v, qaws_surface_geometry* g)
{
	qaws_surface_jet j;
	qaws_status st = qaws_surface_eval_jet(s, u, v, QAWS_SJET_ORDER2, &j);
	if (st != QAWS_STATUS_OK)
		return st;
	return qaws_surface_geometry_eval(&j, NULL, NULL, g, NULL, NULL, NULL);
}

static void test_surface(void)
{
	qaws_surface* s;
	unsigned int i, n;
	int ok_values = 1, ok_fd = 1, ok_adj = 1, ok_t2 = 1;
	double h = DIFF_FD_STEP;

	init_surface_net();
	s = make_surface(g_surf_cps, g_surf_w);
	diff_seed(13);

	for (i = 0; i < GEOM_SAMPLES; i++)
	{
		qaws_scalar u = (qaws_scalar)0.1 + (qaws_scalar)0.15 * (qaws_scalar)i;
		qaws_scalar v = (qaws_scalar)0.85 - (qaws_scalar)0.13 * (qaws_scalar)i;
		qaws_scalar udot = diff_rand(), vdot = diff_rand(), ubar = 0, vbar = 0;
		qaws_scalar dir[48 + 16], pbar[48 + 16];
		qaws_scalar cps_p[48], cps_m[48], w_p[16], w_m[16];
		qaws_field_view fv[2];
		qaws_diff_views views;
		qaws_surface_jet p, tg, tt, ybar;
		qaws_surface_geometry g, gt, gtt, gbar, gp, gm;
		qaws_surface_curvature_result ref;
		qaws_diff_validity validity;
		qaws_surface *sp, *sm;
		double a[13], b[13], lhs, rhs;

		qaws_surface_eval_jet(s, u, v, QAWS_SJET_ORDER3, &p);
		qaws_surface_geometry_eval(&p, NULL, NULL, &g, NULL, NULL, &validity);
		qaws_surface_compute_principal_curvatures(s, u, v, &ref);
		if (validity != QAWS_DIFF_VALID || !diff_close(g.gaussian, ref.gaussian, 1e-4) ||
		    !diff_close(g.mean, ref.mean, 1e-4) ||
		    !diff_close(g.kappa1 + g.kappa2, ref.kappa1 + ref.kappa2, 1e-4) ||
		    !diff_close(g.kappa1 * g.kappa2, ref.kappa1 * ref.kappa2, 1e-4))
			ok_values = 0;

		diff_rand_fill(dir, 64);
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, 16, 3);
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, dir + 48, 16, 1);
		views.fields = fv;
		views.field_count = 2;
		views.children = NULL;
		views.child_count = 0;
		qaws_surface_eval_batch_tangent2(NULL, s, &u, &v, &udot, &vdot, 1, QAWS_SJET_ORDER2, &views, &p, &tg, &tt);
		qaws_surface_geometry_eval(&p, &tg, &tt, &g, &gt, &gtt, NULL);

		for (n = 0; n < 48; n++)
		{
			cps_p[n] = (qaws_scalar)(g_surf_cps[n] + h * dir[n]);
			cps_m[n] = (qaws_scalar)(g_surf_cps[n] - h * dir[n]);
		}
		for (n = 0; n < 16; n++)
		{
			w_p[n] = (qaws_scalar)(g_surf_w[n] + h * dir[48 + n]);
			w_m[n] = (qaws_scalar)(g_surf_w[n] - h * dir[48 + n]);
		}
		sp = make_surface(cps_p, w_p);
		sm = make_surface(cps_m, w_m);
		surface_geometry_at(sp, (qaws_scalar)(u + h * udot), (qaws_scalar)(v + h * vdot), &gp);
		surface_geometry_at(sm, (qaws_scalar)(u - h * udot), (qaws_scalar)(v - h * vdot), &gm);
		flat_surface(&gp, a);
		flat_surface(&gm, b);
		for (n = 0; n < 13; n++)
			a[n] = (a[n] - b[n]) / (2 * h);
		flat_surface(&gt, b);
		if (!flat_close(a, b, 13, DIFF_TOL * 100))
			ok_fd = 0;

		{
			qaws_surface_jet pp, tgp, pm, tgm;
			qaws_surface_geometry g0, gtp, gtm;
			qaws_surface_eval_tangent(NULL, sp, (qaws_scalar)(u + h * udot), (qaws_scalar)(v + h * vdot),
				udot, vdot, QAWS_SJET_ORDER2, &views, &pp, &tgp);
			qaws_surface_eval_tangent(NULL, sm, (qaws_scalar)(u - h * udot), (qaws_scalar)(v - h * vdot),
				udot, vdot, QAWS_SJET_ORDER2, &views, &pm, &tgm);
			qaws_surface_geometry_eval(&pp, &tgp, NULL, &g0, &gtp, NULL, NULL);
			qaws_surface_geometry_eval(&pm, &tgm, NULL, &g0, &gtm, NULL, NULL);
			flat_surface(&gtp, a);
			flat_surface(&gtm, b);
			for (n = 0; n < 13; n++)
				a[n] = (a[n] - b[n]) / (2 * h);
			flat_surface(&gtt, b);
			if (!flat_close(a, b, 13, DIFF_TOL * 300))
				ok_t2 = 0;
		}
		qaws_surface_destroy(sp);
		qaws_surface_destroy(sm);

		gbar.normal = diff_rand_vec3();
		gbar.E = diff_rand(); gbar.F = diff_rand(); gbar.G = diff_rand();
		gbar.L = diff_rand(); gbar.M = diff_rand(); gbar.N = diff_rand();
		gbar.gaussian = diff_rand(); gbar.mean = diff_rand();
		gbar.kappa1 = diff_rand(); gbar.kappa2 = diff_rand();
		flat_surface(&gbar, a);
		flat_surface(&gt, b);
		lhs = 0;
		for (n = 0; n < 13; n++)
			lhs += a[n] * b[n];

		memset(&ybar, 0, sizeof(ybar));
		qaws_surface_geometry_adjoint(&p, &gbar, &ybar, NULL);
		memset(pbar, 0, sizeof(pbar));
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, pbar, 16, 3);
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, pbar + 48, 16, 1);
		qaws_surface_eval_adjoint(NULL, s, u, v, ybar.channels, &ybar, &views, &ubar, &vbar);
		rhs = diff_dot(pbar, dir, 64) + (double)ubar * udot + (double)vbar * vdot;
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
			ok_adj = 0;
	}

	TEST_ASSERT(ok_values, "surface curvatures match inspection");
	TEST_ASSERT(ok_fd, "surface geometry tangent matches finite differences through the chain");
	TEST_ASSERT(ok_t2, "surface geometry second tangent matches finite differences");
	TEST_ASSERT(ok_adj, "surface geometry adjoint identity through the chain");
	qaws_surface_destroy(s);
}

static void test_surface_umbilic(void)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	qaws_surface_jet j;
	qaws_surface_geometry g;
	qaws_diff_validity validity = QAWS_DIFF_VALID;

	d.p00 = qaws_v3(0, 0, 0);
	d.p10 = qaws_v3(1, 0, 0);
	d.p01 = qaws_v3(0, 1, 0);
	d.p11 = qaws_v3(1, 1, 0);
	qaws_surface_create_bilinear(&d, &s);
	qaws_surface_eval_jet(s, (qaws_scalar)0.3, (qaws_scalar)0.6, QAWS_SJET_ORDER2, &j);
	qaws_surface_geometry_eval(&j, NULL, NULL, &g, NULL, NULL, &validity);
	TEST_ASSERT(validity == QAWS_DIFF_ILL_CONDITIONED, "plane is umbilic: principal curvatures flagged");
	TEST_ASSERT(approx_eq(g.normal.z, 1) && approx_eq(g.gaussian, 0) && approx_eq(g.mean, 0), "plane normal and curvatures");
	qaws_surface_destroy(s);
}

/* ------------------------------------------------------------------ */
/*  Mathematica reference (tests/reference/52_geometry.wls)           */
/* ------------------------------------------------------------------ */

#include "reference/52_geometry.h"

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

static void ref_scalars(qaws_scalar const* x, double* out, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		out[i] = x[i];
}

static double const* const g_ref_c3[4][3] = {
	{ ref_c3_value0, ref_c3_value1, ref_c3_value2 },
	{ ref_c3_tangent0, ref_c3_tangent1, ref_c3_tangent2 },
	{ ref_c3_tangent2_0, ref_c3_tangent2_1, ref_c3_tangent2_2 },
	{ ref_c3_adjoint0, ref_c3_adjoint1, ref_c3_adjoint2 }
};
static double const* const g_ref_c2[4][3] = {
	{ ref_c2_value0, ref_c2_value1, ref_c2_value2 },
	{ ref_c2_tangent0, ref_c2_tangent1, ref_c2_tangent2 },
	{ ref_c2_tangent2_0, ref_c2_tangent2_1, ref_c2_tangent2_2 },
	{ ref_c2_adjoint0, ref_c2_adjoint1, ref_c2_adjoint2 }
};
static double const* const g_ref_s[4][3] = {
	{ ref_s_value0, ref_s_value1, ref_s_value2 },
	{ ref_s_tangent0, ref_s_tangent1, ref_s_tangent2 },
	{ ref_s_tangent2_0, ref_s_tangent2_1, ref_s_tangent2_2 },
	{ ref_s_adjoint0, ref_s_adjoint1, ref_s_adjoint2 }
};

/* Reports and checks the largest error of value, tangent, tangent2 and
   adjoint over the samples. */
static void ref_report(char const* name, double const* e)
{
	double tol = QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-11;
	char msg[160];
	printf("    reference %-8s value %.1e, tangent %.1e, tangent2 %.1e, adjoint %.1e\n", name, e[0], e[1], e[2], e[3]);
	sprintf(msg, "reference (%s): values", name);
	TEST_ASSERT(e[0] <= tol, msg);
	sprintf(msg, "reference (%s): tangents through the chain", name);
	TEST_ASSERT(e[1] <= tol, msg);
	sprintf(msg, "reference (%s): second tangents through the chain", name);
	TEST_ASSERT(e[2] <= 10 * tol, msg);
	sprintf(msg, "reference (%s): adjoints through the chain", name);
	TEST_ASSERT(e[3] <= tol, msg);
}

static void ref_track(double* e, int k, double err)
{
	if (err > e[k])
		e[k] = err;
}

static void test_reference_curve_3d(void)
{
	static qaws_scalar const tdots[3] = { (qaws_scalar)0.5, (qaws_scalar)-0.3, (qaws_scalar)0.7 };
	static double const dir_d[18] = { 0.2, -0.3, 0.1, -0.25, 0.2, 0.3, 0.1, 0.1, -0.2, 0.3, -0.2, 0.25,
		-0.2, 0.3, -0.1, 0.1, -0.25, 0.2 };
	static double const gbar_d[12] = { 0.3, -0.2, 0.5, -0.4, 0.1, 0.25, 0.2, -0.3, 0.4, -0.5, 0.6, -0.25 };
	qaws_curve* c = make_curve_3d(g_curve_cps);
	qaws_scalar dir[18], pbar[19];
	qaws_field_view v;
	qaws_diff_views views;
	double e[4] = { 0, 0, 0, 0 }, a[19];
	unsigned int i, n;
	for (n = 0; n < 18; n++)
		dir[n] = (qaws_scalar)dir_d[n];
	for (i = 0; i < 3; i++)
	{
		qaws_scalar t = (qaws_scalar)(i == 0 ? 0.4 : (i == 1 ? 1.3 : 2.5));
		qaws_curve_jet_3d p, tg, tt, ybar;
		qaws_curve_geometry_3d g, gt, gtt, gbar;
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, dir, 6, 3);
		qaws_curve_eval_batch_tangent2_3d(NULL, c, &t, &tdots[i], 1, 0xF, &views, &p, &tg, &tt);
		qaws_curve_geometry_eval_3d(&p, &tg, &tt, &g, &gt, &gtt, NULL);
		flat_curve3(&g, a);
		ref_track(e, 0, ref_max_err(a, g_ref_c3[0][i], 12));
		flat_curve3(&gt, a);
		ref_track(e, 1, ref_max_err(a, g_ref_c3[1][i], 12));
		flat_curve3(&gtt, a);
		ref_track(e, 2, ref_max_err(a, g_ref_c3[2][i], 12));

		gbar.tangent = qaws_v3((qaws_scalar)gbar_d[0], (qaws_scalar)gbar_d[1], (qaws_scalar)gbar_d[2]);
		gbar.normal = qaws_v3((qaws_scalar)gbar_d[3], (qaws_scalar)gbar_d[4], (qaws_scalar)gbar_d[5]);
		gbar.binormal = qaws_v3((qaws_scalar)gbar_d[6], (qaws_scalar)gbar_d[7], (qaws_scalar)gbar_d[8]);
		gbar.speed = (qaws_scalar)gbar_d[9];
		gbar.curvature = (qaws_scalar)gbar_d[10];
		gbar.torsion = (qaws_scalar)gbar_d[11];
		memset(&ybar, 0, sizeof(ybar));
		qaws_curve_geometry_adjoint_3d(&p, &gbar, &ybar, NULL);
		memset(pbar, 0, sizeof(pbar));
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, pbar, 6, 3);
		qaws_curve_eval_adjoint_3d(NULL, c, t, ybar.channels, &ybar, &views, &pbar[18]);
		ref_scalars(pbar, a, 19);
		ref_track(e, 3, ref_max_err(a, g_ref_c3[3][i], 19));
	}
	ref_report("3D curve", e);
	qaws_curve_destroy(c);
}

static void test_reference_curve_2d(void)
{
	static qaws_scalar const tdots[3] = { (qaws_scalar)-0.4, (qaws_scalar)0.3, (qaws_scalar)0.5 };
	static double const dir_d[10] = { 0.3, -0.2, 0.1, 0.4, -0.25, 0.2, 0.1, -0.3, 0.2, 0.1 };
	static double const gbar_d[6] = { 0.2, -0.3, 0.4, 0.1, -0.5, 0.3 };
	qaws_curve* c = make_curve_2d(g_curve2_cps);
	qaws_scalar dir[10], pbar[11];
	qaws_field_view v;
	qaws_diff_views views;
	double e[4] = { 0, 0, 0, 0 }, a[11];
	unsigned int i, n;
	for (n = 0; n < 10; n++)
		dir[n] = (qaws_scalar)dir_d[n];
	for (i = 0; i < 3; i++)
	{
		qaws_scalar t = (qaws_scalar)(i == 0 ? 0.2 : (i == 1 ? 0.5 : 0.8));
		qaws_curve_jet_2d p, tg, tt, ybar;
		qaws_curve_geometry_2d g[3], gbar;
		unsigned int k;
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, dir, 5, 2);
		qaws_curve_eval_batch_tangent2_2d(NULL, c, &t, &tdots[i], 1, 0xF, &views, &p, &tg, &tt);
		qaws_curve_geometry_eval_2d(&p, &tg, &tt, &g[0], &g[1], &g[2], NULL);
		for (k = 0; k < 3; k++)
		{
			a[0] = g[k].tangent.x; a[1] = g[k].tangent.y;
			a[2] = g[k].normal.x; a[3] = g[k].normal.y;
			a[4] = g[k].speed; a[5] = g[k].curvature;
			ref_track(e, (int)k, ref_max_err(a, g_ref_c2[k][i], 6));
		}
		gbar.tangent.x = (qaws_scalar)gbar_d[0]; gbar.tangent.y = (qaws_scalar)gbar_d[1];
		gbar.normal.x = (qaws_scalar)gbar_d[2]; gbar.normal.y = (qaws_scalar)gbar_d[3];
		gbar.speed = (qaws_scalar)gbar_d[4]; gbar.curvature = (qaws_scalar)gbar_d[5];
		memset(&ybar, 0, sizeof(ybar));
		qaws_curve_geometry_adjoint_2d(&p, &gbar, &ybar, NULL);
		memset(pbar, 0, sizeof(pbar));
		one_view(&v, &views, QAWS_FIELD_CONTROL_POINTS, pbar, 5, 2);
		qaws_curve_eval_adjoint_2d(NULL, c, t, ybar.channels, &ybar, &views, &pbar[10]);
		ref_scalars(pbar, a, 11);
		ref_track(e, 3, ref_max_err(a, g_ref_c2[3][i], 11));
	}
	ref_report("2D curve", e);
	qaws_curve_destroy(c);
}

static void test_reference_surface(void)
{
	static double const uv[3][4] = { { 0.2, 0.75, 0.5, -0.3 }, { 0.5, 0.4, -0.2, 0.4 }, { 0.8, 0.1, 0.3, 0.1 } };
	static double const gbar_d[13] = { 0.2, -0.3, 0.5, 0.3, -0.2, 0.1, -0.4, 0.25, 0.3, -0.1, 0.2, 0.4, -0.3 };
	qaws_scalar cps[48], w[16], dir[64], pbar[66];
	qaws_field_view fv[2];
	qaws_diff_views views;
	qaws_surface* s;
	double e[4] = { 0, 0, 0, 0 }, a[66];
	unsigned int i, j, c;
	/* the exact rationals of the script */
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			qaws_scalar* p = &cps[(i * 4 + j) * 3];
			p[0] = (qaws_scalar)(i + ((double)((i + 2 * j) % 3) - 1.0) / 10.0);
			p[1] = (qaws_scalar)(j + ((double)((2 * i + j) % 3) - 1.0) / 10.0);
			p[2] = (qaws_scalar)((16.0 * (i - 1.5) * (i - 1.5) - 10.0 * (j - 1.5) * (j - 1.5) + 4.0 * ((double)((i + j) % 3) - 1.0))
				/ 40.0);
			w[i * 4 + j] = (qaws_scalar)(1.0 + (2.0 * ((i + 3 * j) % 4) - 3.0) / 10.0);
			for (c = 0; c < 3; c++)
				dir[(i * 4 + j) * 3 + c] = (qaws_scalar)(((double)((5 * i + 3 * j + 7 * c) % 9) - 4.0) / 10.0);
			dir[48 + i * 4 + j] = (qaws_scalar)(((double)((3 * i + j) % 5) - 2.0) / 20.0);
		}
	s = make_surface(cps, w);
	for (i = 0; i < 3; i++)
	{
		qaws_scalar u = (qaws_scalar)uv[i][0], v = (qaws_scalar)uv[i][1];
		qaws_scalar udot = (qaws_scalar)uv[i][2], vdot = (qaws_scalar)uv[i][3];
		qaws_surface_jet p, tg, tt, ybar;
		qaws_surface_geometry g, gt, gtt, gbar;
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, 16, 3);
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, dir + 48, 16, 1);
		views.fields = fv;
		views.field_count = 2;
		views.children = NULL;
		views.child_count = 0;
		qaws_surface_eval_batch_tangent2(NULL, s, &u, &v, &udot, &vdot, 1, QAWS_SJET_ORDER2, &views, &p, &tg, &tt);
		qaws_surface_geometry_eval(&p, &tg, &tt, &g, &gt, &gtt, NULL);
		flat_surface(&g, a);
		ref_track(e, 0, ref_max_err(a, g_ref_s[0][i], 13));
		flat_surface(&gt, a);
		ref_track(e, 1, ref_max_err(a, g_ref_s[1][i], 13));
		flat_surface(&gtt, a);
		ref_track(e, 2, ref_max_err(a, g_ref_s[2][i], 13));

		gbar.normal = qaws_v3((qaws_scalar)gbar_d[0], (qaws_scalar)gbar_d[1], (qaws_scalar)gbar_d[2]);
		gbar.E = (qaws_scalar)gbar_d[3]; gbar.F = (qaws_scalar)gbar_d[4]; gbar.G = (qaws_scalar)gbar_d[5];
		gbar.L = (qaws_scalar)gbar_d[6]; gbar.M = (qaws_scalar)gbar_d[7]; gbar.N = (qaws_scalar)gbar_d[8];
		gbar.gaussian = (qaws_scalar)gbar_d[9]; gbar.mean = (qaws_scalar)gbar_d[10];
		gbar.kappa1 = (qaws_scalar)gbar_d[11]; gbar.kappa2 = (qaws_scalar)gbar_d[12];
		memset(&ybar, 0, sizeof(ybar));
		qaws_surface_geometry_adjoint(&p, &gbar, &ybar, NULL);
		memset(pbar, 0, sizeof(pbar));
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, pbar, 16, 3);
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, pbar + 48, 16, 1);
		qaws_surface_eval_adjoint(NULL, s, u, v, ybar.channels, &ybar, &views, &pbar[64], &pbar[65]);
		ref_scalars(pbar, a, 66);
		ref_track(e, 3, ref_max_err(a, g_ref_s[3][i], 66));
	}
	ref_report("surface", e);
	qaws_surface_destroy(s);
}

int test_52_diff_geometry_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 52: Differentiable local geometry\n");
	test_curve_3d();
	test_curve_2d();
	test_curve_straight();
	test_surface();
	test_surface_umbilic();
	test_reference_curve_3d();
	test_reference_curve_2d();
	test_reference_surface();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
