/*
 * Test 58: Constant-speed (arc-length) sampling, first and second order,
 * forward and backward
 *
 *   - evenly spaced samples on a straight segment with a non-uniform
 *     parameterization
 *   - Mathematica ground truth (tests/reference/58_arc_length.wls): t, the
 *     samples, their first and second tangents, the gradient and a
 *     Hessian-vector product of a B-spline, from symbolic derivatives of
 *     the speed and the implicit function theorem at 30 digits
 *   - tangent and tangent2 against finite differences (B-spline, NURBS
 *     with weights, 2D), distance tangents included
 *   - adjoint identity <adjoint, tangent> = <gradient, direction>
 *   - HVP against finite differences of the gradient, symmetric, refused
 *     for rational curves
 */

#include "test_diff.h"
#include "qaws_diff_sampling.h"

#define SM_SAMPLES 3
#define SM_PARAMS 24     /* 6 control points x 3 + 6 weights */

static qaws_scalar const g_sm_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_arc_length_target const g_sm_targets[SM_SAMPLES] = {
	{ (qaws_scalar)0, (qaws_scalar)0.35 }, { (qaws_scalar)1.7, (qaws_scalar)0 }, { (qaws_scalar)-0.5, (qaws_scalar)1 }
};

static qaws_scalar const g_sm_p0[18] = { 0, 0, 0, 1, 2, 0.5f, 2.5f, 2, -0.5f, 3.5f, 0.5f, 1, 5, 1, 0, 6, 3, 0.5f };
static double const g_sm_v[18] = { 0.2, -0.3, 0.1, -0.25, 0.2, 0.3, 0.1, 0.1, -0.2, 0.3, -0.2, 0.25,
	-0.2, 0.3, -0.1, 0.1, -0.25, 0.2 };

/* curve kinds: 0 B-spline 3D, 1 NURBS 3D (weights in x[18..23]), 2 B-spline 2D */
static qaws_curve* sm_curve(int kind, qaws_scalar const* x)
{
	qaws_curve* c = NULL;
	if (kind == 1)
	{
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = x;
		d.control_point_count = 6;
		d.knots = g_sm_knots;
		d.knot_count = 10;
		d.weights = x + 18;
		d.weight_count = 6;
		qaws_curve_create_nurbs(&d, &c);
	}
	else
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = kind == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = x;
		d.control_point_count = 6;
		d.knots = g_sm_knots;
		d.knot_count = 10;
		qaws_curve_create_bspline(&d, &c);
	}
	return c;
}

static unsigned int sm_dim(int kind)
{
	return kind == 2 ? 2u : 3u;
}

static unsigned int sm_param_count(int kind)
{
	return kind == 1 ? 24u : 6u * sm_dim(kind);
}

static qaws_diff_views sm_views(int kind, qaws_field_view* fv, qaws_scalar* data)
{
	qaws_diff_views v;
	unsigned int dim = sm_dim(kind);
	fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, data, 6, dim);
	v.field_count = 1;
	if (kind == 1)
	{
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, data + 18, 6, 1);
		v.field_count = 2;
	}
	v.fields = fv;
	v.children = NULL;
	v.child_count = 0;
	return v;
}

static void sm_base(int kind, qaws_scalar* x)
{
	unsigned int i, k;
	if (kind == 2)
	{
		for (i = 0; i < 6; i++)
			for (k = 0; k < 2; k++)
				x[2 * i + k] = g_sm_p0[3 * i + k];
		return;
	}
	memcpy(x, g_sm_p0, sizeof(g_sm_p0));
	if (kind == 1)
		for (i = 0; i < 6; i++)
			x[18 + i] = (qaws_scalar)(1.0 + 0.3 * sin(1.7 * i));
}

static void sm_values(int kind, qaws_scalar const* x, qaws_scalar const* dist_shift, qaws_arc_length_sample* out)
{
	qaws_arc_length_target tg[SM_SAMPLES];
	qaws_curve* c = sm_curve(kind, x);
	unsigned int i;
	for (i = 0; i < SM_SAMPLES; i++)
	{
		tg[i] = g_sm_targets[i];
		if (dist_shift)
			tg[i].distance += dist_shift[i];
	}
	qaws_curve_arc_length_sample_tangent(NULL, c, tg, NULL, SM_SAMPLES, 0, NULL, out, NULL, NULL, NULL);
	qaws_curve_destroy(c);
}

static double sm_sample_dot(qaws_arc_length_sample const* a, qaws_arc_length_sample const* b, unsigned int n)
{
	double s = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
		s += (double)a[i].t * b[i].t + (double)a[i].position.x * b[i].position.x +
			(double)a[i].position.y * b[i].position.y + (double)a[i].position.z * b[i].position.z;
	return s;
}

static void sm_rand_adjoint(qaws_arc_length_sample* a, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		a[i].t = diff_rand();
		a[i].position = diff_rand_vec3();
	}
}

/* ------------------------------------------------------------------ */

static void test_straight(void)
{
	/* A straight cubic with clustered control points: the parameter is not
	   proportional to arc length, the samples must still be evenly spaced. */
	static qaws_scalar const cps[12] = { 0, 0, 0, 0.2f, 0, 0, 2.6f, 0, 0, 3, 0, 0 };
	qaws_bezier_desc bd;
	qaws_curve* c = NULL;
	qaws_arc_length_target tg[5];
	qaws_arc_length_sample s[5];
	qaws_scalar total = 0;
	unsigned int i;
	int ok = 1;
	bd.dimension = QAWS_DIMENSION_3D;
	bd.degree = 3;
	bd.control_points = cps;
	bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &c);
	for (i = 0; i < 5; i++)
	{
		tg[i].distance = 0;
		tg[i].fraction = (qaws_scalar)(i / 4.0);
	}
	TEST_ASSERT(qaws_curve_arc_length_sample_tangent(NULL, c, tg, NULL, 5, 0, NULL, s, NULL, NULL, &total) == QAWS_STATUS_OK,
		"straight segment samples");
	TEST_ASSERT(approx_eq(total, 3), "segment length");
	for (i = 0; i < 5; i++)
		ok &= diff_close(s[i].position.x, 0.75 * i, 1e3 * DIFF_TOL) && diff_close(s[i].position.y, 0, DIFF_TOL);
	TEST_ASSERT(ok, "constant-speed samples are evenly spaced");
	TEST_ASSERT(!diff_close(s[1].t, 0.25, 1e-2), "the parameter is not arc length");
	qaws_curve_destroy(c);
}

/* Mathematica reference (tests/reference/58_arc_length.wls). */
#include "reference/58_arc_length.h"

static void test_reference(void)
{
	qaws_scalar x[SM_PARAMS], dir[SM_PARAMS], grad[SM_PARAMS], hv[SM_PARAMS];
	qaws_field_view fv[2], fd[2], fg[2], fh[2];
	qaws_diff_views vx, vd, vg, vh;
	qaws_arc_length_sample val[SM_SAMPLES], tan1[SM_SAMPLES], tan2[SM_SAMPLES], adj[SM_SAMPLES];
	qaws_curve* c;
	qaws_scalar total = 0;
	double tol = QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-11;
	unsigned int i;
	int ok_v = 1, ok_1 = 1, ok_2 = 1, ok_g = 1, ok_h = 1;
	static double const tbar[SM_SAMPLES] = { 0.3, -0.2, 0.25 };
	static double const pbar[3 * SM_SAMPLES] = { 0.5, -0.25, 0.2, -0.3, 0.4, 0.1, 0.2, 0.2, -0.4 };
	sm_base(0, x);
	for (i = 0; i < 18; i++)
		dir[i] = (qaws_scalar)g_sm_v[i];
	vx = sm_views(0, fv, x);
	vd = sm_views(0, fd, dir);
	vg = sm_views(0, fg, grad);
	vh = sm_views(0, fh, hv);
	(void)vx;
	c = sm_curve(0, x);
	TEST_ASSERT(qaws_curve_arc_length_sample_tangent(NULL, c, g_sm_targets, NULL, SM_SAMPLES, 0, &vd, val, tan1, tan2, &total)
		== QAWS_STATUS_OK, "reference: forward");
	TEST_ASSERT(diff_close(total, ref_total[0], tol), "reference: total length");
	for (i = 0; i < SM_SAMPLES; i++)
	{
		ok_v &= diff_close(val[i].t, ref_t[i], tol) && diff_close(val[i].position.x, ref_p[3 * i], tol) &&
			diff_close(val[i].position.y, ref_p[3 * i + 1], tol) && diff_close(val[i].position.z, ref_p[3 * i + 2], tol);
		ok_1 &= diff_close(tan1[i].t, ref_t1[i], tol) && diff_close(tan1[i].position.x, ref_p1[3 * i], tol) &&
			diff_close(tan1[i].position.y, ref_p1[3 * i + 1], tol) && diff_close(tan1[i].position.z, ref_p1[3 * i + 2], tol);
		ok_2 &= diff_close(tan2[i].t, ref_t2[i], 10 * tol) && diff_close(tan2[i].position.x, ref_p2[3 * i], 10 * tol) &&
			diff_close(tan2[i].position.y, ref_p2[3 * i + 1], 10 * tol) && diff_close(tan2[i].position.z, ref_p2[3 * i + 2], 10 * tol);
		adj[i].t = (qaws_scalar)tbar[i];
		adj[i].position.x = (qaws_scalar)pbar[3 * i];
		adj[i].position.y = (qaws_scalar)pbar[3 * i + 1];
		adj[i].position.z = (qaws_scalar)pbar[3 * i + 2];
	}
	TEST_ASSERT(ok_v, "reference: parameters and samples");
	TEST_ASSERT(ok_1, "reference: first order tangents");
	TEST_ASSERT(ok_2, "reference: second order tangents");
	memset(grad, 0, sizeof(grad));
	memset(hv, 0, sizeof(hv));
	TEST_ASSERT(qaws_curve_arc_length_sample_adjoint(NULL, c, g_sm_targets, SM_SAMPLES, 0, adj, &vg, NULL) == QAWS_STATUS_OK,
		"reference: adjoint");
	TEST_ASSERT(qaws_curve_arc_length_sample_hvp(NULL, c, g_sm_targets, SM_SAMPLES, 0, adj, &vd, &vh) == QAWS_STATUS_OK,
		"reference: hvp");
	for (i = 0; i < 18; i++)
	{
		ok_g &= diff_close(grad[i], ref_grad[i], tol);
		ok_h &= diff_close(hv[i], ref_hvp[i], 10 * tol);
	}
	TEST_ASSERT(ok_g, "reference: gradient");
	TEST_ASSERT(ok_h, "reference: Hessian-vector product");
	qaws_curve_destroy(c);
}

static void test_finite_differences(int kind, char const* name)
{
	unsigned int np = sm_param_count(kind), i, k;
	qaws_scalar x[SM_PARAMS], dir[SM_PARAMS], xp[SM_PARAMS], xm[SM_PARAMS], grad[SM_PARAMS];
	qaws_scalar ddir[SM_SAMPLES], dadj[SM_SAMPLES], sp[SM_SAMPLES], smn[SM_SAMPLES];
	qaws_field_view fd[2], fg[2];
	qaws_diff_views vd, vg;
	qaws_arc_length_sample val[SM_SAMPLES], t1[SM_SAMPLES], t2[SM_SAMPLES], vp[SM_SAMPLES], vm[SM_SAMPLES];
	qaws_arc_length_sample t1p[SM_SAMPLES], t1m[SM_SAMPLES], adj[SM_SAMPLES];
	qaws_curve* c;
	double h = QAWS_SCALAR_IS_FLOAT ? 1e-2 : 1e-5, tol = QAWS_SCALAR_IS_FLOAT ? 3e-2 : 1e-6;
	int ok1 = 1, ok2 = 1;
	char msg[160];
	diff_seed(77u + (unsigned int)kind);
	sm_base(kind, x);
	diff_rand_fill(dir, np);
	if (kind == 1)
		for (i = 18; i < 24; i++)
			dir[i] *= (qaws_scalar)0.3;
	diff_rand_fill(ddir, SM_SAMPLES);
	vd = sm_views(kind, fd, dir);
	vg = sm_views(kind, fg, grad);
	c = sm_curve(kind, x);
	qaws_curve_arc_length_sample_tangent(NULL, c, g_sm_targets, ddir, SM_SAMPLES, 0, &vd, val, t1, t2, NULL);
	for (i = 0; i < np; i++)
	{
		xp[i] = (qaws_scalar)(x[i] + h * dir[i]);
		xm[i] = (qaws_scalar)(x[i] - h * dir[i]);
	}
	for (i = 0; i < SM_SAMPLES; i++)
	{
		sp[i] = (qaws_scalar)(h * ddir[i]);
		smn[i] = (qaws_scalar)(-h * ddir[i]);
	}
	sm_values(kind, xp, sp, vp);
	sm_values(kind, xm, smn, vm);
	{
		/* tangents at the shifted points, for the second order check */
		qaws_curve* cp = sm_curve(kind, xp);
		qaws_curve* cm = sm_curve(kind, xm);
		qaws_arc_length_target tp[SM_SAMPLES], tm[SM_SAMPLES];
		for (i = 0; i < SM_SAMPLES; i++)
		{
			tp[i] = g_sm_targets[i];
			tm[i] = g_sm_targets[i];
			tp[i].distance += sp[i];
			tm[i].distance += smn[i];
		}
		qaws_curve_arc_length_sample_tangent(NULL, cp, tp, ddir, SM_SAMPLES, 0, &vd, NULL, t1p, NULL, NULL);
		qaws_curve_arc_length_sample_tangent(NULL, cm, tm, ddir, SM_SAMPLES, 0, &vd, NULL, t1m, NULL, NULL);
		qaws_curve_destroy(cp);
		qaws_curve_destroy(cm);
	}
	for (i = 0; i < SM_SAMPLES; i++)
	{
		ok1 &= diff_close(t1[i].t, (vp[i].t - vm[i].t) / (2 * h), tol);
		ok1 &= diff_close(t1[i].position.x, (vp[i].position.x - vm[i].position.x) / (2 * h), tol);
		ok1 &= diff_close(t1[i].position.y, (vp[i].position.y - vm[i].position.y) / (2 * h), tol);
		ok1 &= diff_close(t1[i].position.z, (vp[i].position.z - vm[i].position.z) / (2 * h), tol);
		ok2 &= diff_close(t2[i].t, (t1p[i].t - t1m[i].t) / (2 * h), tol);
		ok2 &= diff_close(t2[i].position.x, (t1p[i].position.x - t1m[i].position.x) / (2 * h), tol);
		ok2 &= diff_close(t2[i].position.y, (t1p[i].position.y - t1m[i].position.y) / (2 * h), tol);
		ok2 &= diff_close(t2[i].position.z, (t1p[i].position.z - t1m[i].position.z) / (2 * h), tol);
	}
	sprintf(msg, "%s: first order tangents match finite differences", name);
	TEST_ASSERT(ok1, msg);
	sprintf(msg, "%s: second order tangents match finite differences of the tangents", name);
	TEST_ASSERT(ok2, msg);

	/* adjoint identity, distances included */
	sm_rand_adjoint(adj, SM_SAMPLES);
	if (kind == 2)
		for (i = 0; i < SM_SAMPLES; i++)
			adj[i].position.z = 0;
	memset(grad, 0, sizeof(grad));
	memset(dadj, 0, sizeof(dadj));
	qaws_curve_arc_length_sample_adjoint(NULL, c, g_sm_targets, SM_SAMPLES, 0, adj, &vg, dadj);
	{
		double lhs = sm_sample_dot(adj, t1, SM_SAMPLES);
		double rhs = diff_dot(grad, dir, np) + diff_dot(dadj, ddir, SM_SAMPLES);
		sprintf(msg, "%s: adjoint identity (%.12g vs %.12g)", name, lhs, rhs);
		TEST_ASSERT(diff_close(lhs, rhs, QAWS_SCALAR_IS_FLOAT ? 2e-3 : 1e-10), msg);
	}

	/* HVP */
	{
		qaws_scalar hv[SM_PARAMS], gp[SM_PARAMS], gm[SM_PARAMS];
		qaws_field_view fh[2], fp[2], fm[2];
		qaws_diff_views vh = sm_views(kind, fh, hv), vgp = sm_views(kind, fp, gp), vgm = sm_views(kind, fm, gm);
		qaws_status st;
		memset(hv, 0, sizeof(hv));
		st = qaws_curve_arc_length_sample_hvp(NULL, c, g_sm_targets, SM_SAMPLES, 0, adj, &vd, &vh);
		if (kind == 1)
		{
			sprintf(msg, "%s: HVP refused for a rational curve", name);
			TEST_ASSERT(st == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
		}
		else
		{
			qaws_curve* cp;
			qaws_curve* cm;
			int okh = 1;
			for (i = 0; i < np; i++)
			{
				xp[i] = (qaws_scalar)(x[i] + h * dir[i]);
				xm[i] = (qaws_scalar)(x[i] - h * dir[i]);
			}
			cp = sm_curve(kind, xp);
			cm = sm_curve(kind, xm);
			memset(gp, 0, sizeof(gp));
			memset(gm, 0, sizeof(gm));
			qaws_curve_arc_length_sample_adjoint(NULL, cp, g_sm_targets, SM_SAMPLES, 0, adj, &vgp, NULL);
			qaws_curve_arc_length_sample_adjoint(NULL, cm, g_sm_targets, SM_SAMPLES, 0, adj, &vgm, NULL);
			for (k = 0; k < np; k++)
				okh &= diff_close(hv[k], (gp[k] - gm[k]) / (2 * h), tol);
			sprintf(msg, "%s: HVP matches finite differences of the gradient", name);
			TEST_ASSERT(st == QAWS_STATUS_OK && okh, msg);
			qaws_curve_destroy(cp);
			qaws_curve_destroy(cm);
			{
				/* symmetry: <u, H v> = <v, H u> */
				qaws_scalar u[SM_PARAMS], hu[SM_PARAMS];
				qaws_field_view fu[2], fhu[2];
				qaws_diff_views vu = sm_views(kind, fu, u), vhu = sm_views(kind, fhu, hu);
				diff_rand_fill(u, np);
				memset(hu, 0, sizeof(hu));
				qaws_curve_arc_length_sample_hvp(NULL, c, g_sm_targets, SM_SAMPLES, 0, adj, &vu, &vhu);
				sprintf(msg, "%s: HVP is symmetric", name);
				TEST_ASSERT(diff_close(diff_dot(u, hv, np), diff_dot(dir, hu, np), QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-9), msg);
			}
		}
	}
	qaws_curve_destroy(c);
}

static void test_refusals(void)
{
	qaws_scalar x[SM_PARAMS], k[10];
	qaws_field_view fk;
	qaws_diff_views vk;
	qaws_arc_length_sample s[SM_SAMPLES];
	qaws_curve* c;
	sm_base(0, x);
	c = sm_curve(0, x);
	memset(k, 0, sizeof(k));
	fk = qaws_field_view_make(QAWS_FIELD_KNOTS, k, 10, 1);
	vk.fields = &fk;
	vk.field_count = 1;
	vk.children = NULL;
	vk.child_count = 0;
	TEST_ASSERT(qaws_curve_arc_length_sample_tangent(NULL, c, g_sm_targets, NULL, SM_SAMPLES, 0, &vk, s, s, NULL, NULL) ==
		QAWS_STATUS_UNSUPPORTED_OPERATION, "knots move the quadrature spans: refused");
	TEST_ASSERT(qaws_curve_arc_length_sample_tangent(NULL, NULL, g_sm_targets, NULL, SM_SAMPLES, 0, NULL, s, NULL, NULL, NULL) ==
		QAWS_STATUS_INVALID_ARGUMENT, "NULL curve");
	qaws_curve_destroy(c);
}

int test_58_diff_sampling_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 58: Constant-speed sampling derivatives\n");
	test_straight();
	test_reference();
	test_finite_differences(0, "B-spline");
	test_finite_differences(1, "NURBS");
	test_finite_differences(2, "2D B-spline");
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
