/*
 * Test 59: Inverse-CDF warping of the unit square onto surfaces, first and
 * second order, forward and backward
 *
 *   - samples of a flat square are the warped points themselves (area)
 *   - tangent and tangent2 against finite differences (Bezier patch, NURBS
 *     patch with weights; area and density measures), xi tangents included
 *   - adjoint identity <adjoint, tangent> = <gradient, direction> + xi part
 *   - HVP against finite differences of the gradient, and symmetric
 */

#include "test_diff.h"
#include "qaws_diff_sampling.h"

#define SS_SAMPLES 3
#define SS_CP 16
#define SS_PARAMS (3 * SS_CP + SS_CP)
#define SS_CELLS 4
#define SS_QUAD 6

static qaws_scalar const g_ss_xi[2 * SS_SAMPLES] = { 0.3f, 0.6f, 0.75f, 0.2f, 0.5f, 0.9f };
static qaws_scalar const g_ss_knots[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

/* rho = 1 + (x^2 + y^2) / 4 + z / 3 */
static qaws_scalar ss_density(qaws_vec3 p, void* user, qaws_vec3* g, qaws_scalar* H)
{
	(void)user;
	g->x = (qaws_scalar)(p.x / 2);
	g->y = (qaws_scalar)(p.y / 2);
	g->z = (qaws_scalar)(1.0 / 3);
	H[0] = (qaws_scalar)0.5; H[1] = 0; H[2] = 0;
	H[3] = (qaws_scalar)0.5; H[4] = 0; H[5] = 0;
	return (qaws_scalar)(1 + (p.x * p.x + p.y * p.y) / 4 + p.z / 3);
}

static qaws_sample_measure_desc const g_ss_density = { QAWS_MEASURE_DENSITY, 0, ss_density, NULL };
static qaws_sample_measure_desc const g_ss_curvature = { QAWS_MEASURE_CURVATURE, (qaws_scalar)0.5, NULL, NULL };

/* kind 0: Bezier bicubic, 1: NURBS bicubic (weights in x[48..63]) */
static qaws_surface* ss_surface(int kind, qaws_scalar const* x)
{
	qaws_surface* s = NULL;
	if (kind == 0)
	{
		qaws_surface_bezier_desc d;
		d.u_degree = 3;
		d.v_degree = 3;
		d.control_points = (qaws_vec3 const*)x;
		d.u_point_count = 4;
		d.v_point_count = 4;
		qaws_surface_create_bezier(&d, &s);
	}
	else
	{
		qaws_surface_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.u_degree = 3;
		d.v_degree = 3;
		d.control_points = (qaws_vec3 const*)x;
		d.u_point_count = 4;
		d.v_point_count = 4;
		d.weights = x + 48;
		d.u_knots = g_ss_knots;
		d.u_knot_count = 8;
		d.v_knots = g_ss_knots;
		d.v_knot_count = 8;
		qaws_surface_create_nurbs(&d, &s);
	}
	return s;
}

static unsigned int ss_param_count(int kind)
{
	return kind == 1 ? 64u : 48u;
}

static qaws_diff_views ss_views(int kind, qaws_field_view* fv, qaws_scalar* data)
{
	qaws_diff_views v;
	fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, data, SS_CP, 3);
	v.field_count = 1;
	if (kind == 1)
	{
		fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, data + 48, SS_CP, 1);
		v.field_count = 2;
	}
	v.fields = fv;
	v.children = NULL;
	v.child_count = 0;
	return v;
}

/* A gently curved patch over [0, 3] x [0, 3]. */
static void ss_base(int kind, qaws_scalar* x)
{
	unsigned int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			qaws_scalar* p = x + 3 * (i * 4 + j);
			p[0] = (qaws_scalar)j;
			p[1] = (qaws_scalar)i;
			p[2] = (qaws_scalar)(0.4 * sin(1.3 * i + 0.7 * j) + 0.2 * i * (3 - j) / 3.0);
		}
	if (kind == 1)
		for (i = 0; i < SS_CP; i++)
			x[48 + i] = (qaws_scalar)(1.0 + 0.25 * sin(2.1 * i));
}

static double ss_dot(qaws_surface_cdf_sample const* a, qaws_surface_cdf_sample const* b, unsigned int n)
{
	double s = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
		s += (double)a[i].u * b[i].u + (double)a[i].v * b[i].v + (double)a[i].position.x * b[i].position.x +
			(double)a[i].position.y * b[i].position.y + (double)a[i].position.z * b[i].position.z;
	return s;
}

static void test_flat(void)
{
	/* A flat, uniformly parameterized square: the area measure is uniform,
	   so the samples are the points xi scaled to the square. */
	qaws_scalar x[48];
	qaws_surface* s;
	qaws_surface_cdf_sample out[SS_SAMPLES];
	qaws_scalar total = 0;
	unsigned int i, j;
	int ok = 1;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			x[3 * (i * 4 + j)] = (qaws_scalar)j;
			x[3 * (i * 4 + j) + 1] = (qaws_scalar)i;
			x[3 * (i * 4 + j) + 2] = 0;
		}
	s = ss_surface(0, x);
	TEST_ASSERT(qaws_surface_cdf_sample_tangent(NULL, s, NULL, g_ss_xi, NULL, SS_SAMPLES, SS_CELLS, SS_QUAD, NULL, out, NULL, NULL,
		&total) == QAWS_STATUS_OK, "flat square: samples");
	TEST_ASSERT(diff_close(total, 9, 1e3 * DIFF_TOL), "flat square: area 9");
	for (i = 0; i < SS_SAMPLES; i++)
		ok &= diff_close(out[i].u, g_ss_xi[2 * i], 1e3 * DIFF_TOL) && diff_close(out[i].v, g_ss_xi[2 * i + 1], 1e3 * DIFF_TOL);
	TEST_ASSERT(ok, "flat square: (u, v) are the points xi");
	/* no curvature: the curvature measure is the floor times the area */
	TEST_ASSERT(qaws_surface_cdf_sample_tangent(NULL, s, &g_ss_curvature, g_ss_xi, NULL, SS_SAMPLES, SS_CELLS, SS_QUAD, NULL, out,
		NULL, NULL, &total) == QAWS_STATUS_OK, "flat square: curvature samples");
	TEST_ASSERT(diff_close(total, 9 * 0.5, 1e3 * DIFF_TOL), "flat square: curvature measure 0.5 x 9");
	qaws_surface_destroy(s);
}

/* Quarter cylinder of radius 2 and height 3 (exact NURBS): the principal
   curvatures are 1/2 and 0, so the curvature measure is sqrt(f^2 + 1/4)
   times the area 3 pi. */
static void test_cylinder(void)
{
	double const r = 2, hgt = 3, c = sqrt(0.5);
	qaws_vec3 cps[9];
	qaws_scalar w[9], total = 0;
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	qaws_surface_cdf_sample out[SS_SAMPLES];
	unsigned int i, j;
	double expect = sqrt(0.25 + 0.25) * 3 * 3.14159265358979323846;
	for (i = 0; i < 3; i++)       /* u: around the circle */
		for (j = 0; j < 3; j++)   /* v: along the axis */
		{
			qaws_vec3* p = &cps[i * 3 + j];
			p->x = (qaws_scalar)(i == 0 ? r : (i == 1 ? r : 0));
			p->y = (qaws_scalar)(i == 0 ? 0 : r);
			p->z = (qaws_scalar)(hgt * j / 2.0);
			w[i * 3 + j] = (qaws_scalar)(i == 1 ? c : 1.0);
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cps;
	d.u_point_count = 3;
	d.v_point_count = 3;
	d.weights = w;
	TEST_ASSERT_STATUS(qaws_surface_create_nurbs(&d, &s));
	TEST_ASSERT(qaws_surface_cdf_sample_tangent(NULL, s, &g_ss_curvature, g_ss_xi, NULL, SS_SAMPLES, 8, 8, NULL, out, NULL, NULL,
		&total) == QAWS_STATUS_OK, "quarter cylinder: curvature samples");
	{
		char msg[128];
		sprintf(msg, "quarter cylinder: curvature measure sqrt(1/2) x 3 pi (%.12g vs %.12g)", (double)total, expect);
		TEST_ASSERT(diff_close(total, expect, QAWS_SCALAR_IS_FLOAT ? 1e-4 : 1e-10), msg);
	}
	qaws_surface_destroy(s);
}

/* The forward-over-reverse HVP of a Bezier patch (linear family) against
   the polarized HVP of the same patch as a NURBS with unit weights. */
static void test_hvp_paths(qaws_sample_measure_desc const* m, char const* name)
{
	qaws_scalar x[SS_PARAMS], dir[48], hb[48], hn[SS_PARAMS];
	qaws_surface_cdf_sample adj[SS_SAMPLES];
	qaws_field_view fd, fb, fn[2];
	qaws_diff_views vd, vb, vn;
	qaws_surface* sb;
	qaws_surface* sn;
	unsigned int i;
	int ok = 1;
	char msg[160];
	diff_seed(77u);
	ss_base(0, x);
	for (i = 0; i < SS_CP; i++)
		x[48 + i] = 1;
	diff_rand_fill(dir, 48);
	for (i = 0; i < SS_SAMPLES; i++)
	{
		adj[i].u = diff_rand();
		adj[i].v = diff_rand();
		adj[i].position = diff_rand_vec3();
	}
	sb = ss_surface(0, x);
	sn = ss_surface(1, x);
	fd = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, SS_CP, 3);
	fb = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, hb, SS_CP, 3);
	vd.fields = &fd;
	vb.fields = &fb;
	vd.field_count = vb.field_count = 1;
	vd.children = vb.children = NULL;
	vd.child_count = vb.child_count = 0;
	vn = ss_views(1, fn, hn);
	memset(hb, 0, sizeof(hb));
	memset(hn, 0, sizeof(hn));
	TEST_ASSERT_STATUS(qaws_surface_cdf_sample_hvp(NULL, sb, m, g_ss_xi, SS_SAMPLES, SS_CELLS, SS_QUAD, adj, &vd, &vb));
	TEST_ASSERT_STATUS(qaws_surface_cdf_sample_hvp(NULL, sn, m, g_ss_xi, SS_SAMPLES, SS_CELLS, SS_QUAD, adj, &vd, &vn));
	for (i = 0; i < 48; i++)
		ok &= diff_close(hb[i], hn[i], QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-9);
	sprintf(msg, "%s: forward-over-reverse HVP matches the polarized HVP", name);
	TEST_ASSERT(ok, msg);
	qaws_surface_destroy(sb);
	qaws_surface_destroy(sn);
}

static void ss_values(int kind, qaws_sample_measure_desc const* m, qaws_scalar const* x, qaws_scalar const* xi,
	qaws_diff_views const* dir, qaws_scalar const* xid, qaws_surface_cdf_sample* val, qaws_surface_cdf_sample* t1)
{
	qaws_surface* s = ss_surface(kind, x);
	qaws_surface_cdf_sample_tangent(NULL, s, m, xi, xid, SS_SAMPLES, SS_CELLS, SS_QUAD, dir, val, t1, NULL, NULL);
	qaws_surface_destroy(s);
}

static void test_surface(int kind, qaws_sample_measure_desc const* m, char const* name)
{
	unsigned int np = ss_param_count(kind), i, k;
	qaws_scalar x[SS_PARAMS], dir[SS_PARAMS], xp[SS_PARAMS], xm[SS_PARAMS], grad[SS_PARAMS], hv[SS_PARAMS];
	qaws_scalar xid[2 * SS_SAMPLES], xip[2 * SS_SAMPLES], xim[2 * SS_SAMPLES], xadj[2 * SS_SAMPLES];
	qaws_field_view fd[2], fg[2], fh[2];
	qaws_diff_views vd, vg, vh;
	qaws_surface_cdf_sample val[SS_SAMPLES], t1[SS_SAMPLES], t2[SS_SAMPLES], vp[SS_SAMPLES], vm[SS_SAMPLES];
	qaws_surface_cdf_sample t1p[SS_SAMPLES], t1m[SS_SAMPLES], adj[SS_SAMPLES];
	qaws_surface* s;
	double h = QAWS_SCALAR_IS_FLOAT ? 1e-2 : 1e-5, tol = QAWS_SCALAR_IS_FLOAT ? 3e-2 : 1e-6;
	int ok1 = 1, ok2 = 1;
	char msg[160];
	diff_seed(5959u + (unsigned int)kind);
	ss_base(kind, x);
	diff_rand_fill(dir, np);
	for (i = 0; i < np; i++)
		dir[i] *= (qaws_scalar)(i >= 48 ? 0.2 : 0.3);
	diff_rand_fill(xid, 2 * SS_SAMPLES);
	vd = ss_views(kind, fd, dir);
	vg = ss_views(kind, fg, grad);
	vh = ss_views(kind, fh, hv);
	s = ss_surface(kind, x);
	TEST_ASSERT(qaws_surface_cdf_sample_tangent(NULL, s, m, g_ss_xi, xid, SS_SAMPLES, SS_CELLS, SS_QUAD, &vd, val, t1, t2, NULL) ==
		QAWS_STATUS_OK, name);
	for (i = 0; i < np; i++)
	{
		xp[i] = (qaws_scalar)(x[i] + h * dir[i]);
		xm[i] = (qaws_scalar)(x[i] - h * dir[i]);
	}
	for (i = 0; i < 2 * SS_SAMPLES; i++)
	{
		xip[i] = (qaws_scalar)(g_ss_xi[i] + h * xid[i]);
		xim[i] = (qaws_scalar)(g_ss_xi[i] - h * xid[i]);
	}
	ss_values(kind, m, xp, xip, &vd, xid, vp, t1p);
	ss_values(kind, m, xm, xim, &vd, xid, vm, t1m);
	for (i = 0; i < SS_SAMPLES; i++)
	{
		ok1 &= diff_close(t1[i].u, (vp[i].u - vm[i].u) / (2 * h), tol) && diff_close(t1[i].v, (vp[i].v - vm[i].v) / (2 * h), tol);
		ok1 &= diff_close(t1[i].position.x, (vp[i].position.x - vm[i].position.x) / (2 * h), tol);
		ok1 &= diff_close(t1[i].position.y, (vp[i].position.y - vm[i].position.y) / (2 * h), tol);
		ok1 &= diff_close(t1[i].position.z, (vp[i].position.z - vm[i].position.z) / (2 * h), tol);
		ok2 &= diff_close(t2[i].u, (t1p[i].u - t1m[i].u) / (2 * h), tol) && diff_close(t2[i].v, (t1p[i].v - t1m[i].v) / (2 * h), tol);
		ok2 &= diff_close(t2[i].position.x, (t1p[i].position.x - t1m[i].position.x) / (2 * h), tol);
		ok2 &= diff_close(t2[i].position.y, (t1p[i].position.y - t1m[i].position.y) / (2 * h), tol);
		ok2 &= diff_close(t2[i].position.z, (t1p[i].position.z - t1m[i].position.z) / (2 * h), tol);
	}
	sprintf(msg, "%s: first order tangents match finite differences", name);
	TEST_ASSERT(ok1, msg);
	sprintf(msg, "%s: second order tangents match finite differences of the tangents", name);
	TEST_ASSERT(ok2, msg);

	/* adjoint identity, xi included */
	for (i = 0; i < SS_SAMPLES; i++)
	{
		adj[i].u = diff_rand();
		adj[i].v = diff_rand();
		adj[i].position = diff_rand_vec3();
	}
	memset(grad, 0, sizeof(grad));
	memset(xadj, 0, sizeof(xadj));
	qaws_surface_cdf_sample_adjoint(NULL, s, m, g_ss_xi, SS_SAMPLES, SS_CELLS, SS_QUAD, adj, &vg, xadj);
	{
		double lhs = ss_dot(adj, t1, SS_SAMPLES), rhs = diff_dot(grad, dir, np) + diff_dot(xadj, xid, 2 * SS_SAMPLES);
		sprintf(msg, "%s: adjoint identity (%.12g vs %.12g)", name, lhs, rhs);
		TEST_ASSERT(diff_close(lhs, rhs, QAWS_SCALAR_IS_FLOAT ? 2e-3 : 1e-10), msg);
	}

	/* HVP against finite differences of the gradient */
	{
		qaws_scalar gp[SS_PARAMS], gm[SS_PARAMS];
		qaws_field_view fp[2], fm[2];
		qaws_diff_views vgp = ss_views(kind, fp, gp), vgm = ss_views(kind, fm, gm);
		qaws_surface* sp = ss_surface(kind, xp);
		qaws_surface* sm = ss_surface(kind, xm);
		int okh = 1;
		memset(hv, 0, sizeof(hv));
		memset(gp, 0, sizeof(gp));
		memset(gm, 0, sizeof(gm));
		TEST_ASSERT(qaws_surface_cdf_sample_hvp(NULL, s, m, g_ss_xi, SS_SAMPLES, SS_CELLS, SS_QUAD, adj, &vd, &vh) == QAWS_STATUS_OK,
			name);
		qaws_surface_cdf_sample_adjoint(NULL, sp, m, g_ss_xi, SS_SAMPLES, SS_CELLS, SS_QUAD, adj, &vgp, NULL);
		qaws_surface_cdf_sample_adjoint(NULL, sm, m, g_ss_xi, SS_SAMPLES, SS_CELLS, SS_QUAD, adj, &vgm, NULL);
		for (k = 0; k < np; k++)
			okh &= diff_close(hv[k], (gp[k] - gm[k]) / (2 * h), QAWS_SCALAR_IS_FLOAT ? 5e-2 : 1e-5);
		sprintf(msg, "%s: HVP matches finite differences of the gradient", name);
		TEST_ASSERT(okh, msg);
		qaws_surface_destroy(sp);
		qaws_surface_destroy(sm);
		{
			qaws_scalar u[SS_PARAMS], hu[SS_PARAMS];
			qaws_field_view fu[2], fhu[2];
			qaws_diff_views vu = ss_views(kind, fu, u), vhu = ss_views(kind, fhu, hu);
			diff_rand_fill(u, np);
			memset(hu, 0, sizeof(hu));
			qaws_surface_cdf_sample_hvp(NULL, s, m, g_ss_xi, SS_SAMPLES, SS_CELLS, SS_QUAD, adj, &vu, &vhu);
			sprintf(msg, "%s: HVP is symmetric", name);
			TEST_ASSERT(diff_close(diff_dot(u, hv, np), diff_dot(dir, hu, np), QAWS_SCALAR_IS_FLOAT ? 5e-2 : 1e-6), msg);
		}
	}
	qaws_surface_destroy(s);
}

int test_59_diff_surface_sampling_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 59: Surface inverse-CDF sampling derivatives\n");
	test_flat();
	test_surface(0, NULL, "Bezier patch, area");
	test_surface(0, &g_ss_density, "Bezier patch, density");
	test_surface(1, NULL, "NURBS patch, area");
	test_surface(1, &g_ss_density, "NURBS patch, density");
	test_surface(0, &g_ss_curvature, "Bezier patch, curvature");
	test_surface(1, &g_ss_curvature, "NURBS patch, curvature");
	test_cylinder();
	test_hvp_paths(NULL, "area");
	test_hvp_paths(&g_ss_density, "density");
	test_hvp_paths(&g_ss_curvature, "curvature");
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
