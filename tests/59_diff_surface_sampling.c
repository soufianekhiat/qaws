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
	qaws_surface_destroy(s);
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

/* ------------------------------------------------------------------ */
/*  Mathematica reference (tests/reference/59_surface_sampling.wls)   */
/* ------------------------------------------------------------------ */

#include "reference/59_surface_sampling.h"

typedef struct ss_ref
{
	char const* name;
	qaws_sample_measure_desc const* measure;
	double const* total;
	double const* value[SS_SAMPLES];
	double const* tangent[SS_SAMPLES];
	double const* tangent2[SS_SAMPLES];
	double const* grad[SS_SAMPLES];
	double const* hvp[SS_SAMPLES];
} ss_ref;

static ss_ref const g_ss_refs[2] = {
	{ "area", NULL, ref_area_total, { ref_area_value0, ref_area_value1, ref_area_value2 },
		{ ref_area_tangent0, ref_area_tangent1, ref_area_tangent2 },
		{ ref_area_tangent2_0, ref_area_tangent2_1, ref_area_tangent2_2 }, { ref_area_grad0, ref_area_grad1, ref_area_grad2 },
		{ ref_area_hvp0, ref_area_hvp1, ref_area_hvp2 } },
	{ "density", &g_ss_density, ref_dens_total, { ref_dens_value0, ref_dens_value1, ref_dens_value2 },
		{ ref_dens_tangent0, ref_dens_tangent1, ref_dens_tangent2 },
		{ ref_dens_tangent2_0, ref_dens_tangent2_1, ref_dens_tangent2_2 }, { ref_dens_grad0, ref_dens_grad1, ref_dens_grad2 },
		{ ref_dens_hvp0, ref_dens_hvp1, ref_dens_hvp2 } }
};

static double ref_err(double a, double b)
{
	double scale = 1.0;
	if (fabs(a) > scale) scale = fabs(a);
	if (fabs(b) > scale) scale = fabs(b);
	return fabs(a - b) / scale;
}

static double ref_sample_err(qaws_surface_cdf_sample const* s, double const* ref)
{
	double e = ref_err(s->u, ref[0]);
	e = fmax(e, ref_err(s->v, ref[1]));
	e = fmax(e, ref_err(s->position.x, ref[2]));
	e = fmax(e, ref_err(s->position.y, ref[3]));
	return fmax(e, ref_err(s->position.z, ref[4]));
}

/* Largest errors of { total and samples, tangents, tangent2, gradient, HVP }
   of the bicubic Bezier patch of the script with `cells` x `cells` cells. */
static void ref_errors(ss_ref const* r, unsigned int cells, unsigned int quadrature, double* e)
{
	static double const xi_d[2 * SS_SAMPLES] = { 0.3, 0.6, 0.75, 0.2, 0.5, 0.9 };
	qaws_scalar x[48], dir[48], xi[2 * SS_SAMPLES], grad[48], hv[48], total = 0;
	qaws_field_view fd[2], fg[2], fh[2];
	qaws_diff_views vd, vg, vh;
	qaws_surface_cdf_sample val[SS_SAMPLES], t1[SS_SAMPLES], t2[SS_SAMPLES], adj[SS_SAMPLES];
	qaws_surface* s;
	unsigned int i, j, k;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			x[3 * (i * 4 + j)] = (qaws_scalar)j;
			x[3 * (i * 4 + j) + 1] = (qaws_scalar)i;
			x[3 * (i * 4 + j) + 2] = (qaws_scalar)(((double)((i + 2 * j) % 4) - 1.5) / 5.0);
		}
	for (k = 0; k < 48; k++)
		dir[k] = (qaws_scalar)(((double)((7 * k + 3) % 11) - 5.0) / 10.0);
	for (k = 0; k < 2 * SS_SAMPLES; k++)
		xi[k] = (qaws_scalar)xi_d[k];
	/* ybar[k] = (((5 k + 2) mod 9) - 4) / 10 over (u, v, x, y, z) per sample */
	for (i = 0; i < SS_SAMPLES; i++)
	{
		qaws_scalar b[5];
		for (k = 0; k < 5; k++)
			b[k] = (qaws_scalar)(((double)((5 * (5 * i + k) + 2) % 9) - 4.0) / 10.0);
		adj[i].u = b[0];
		adj[i].v = b[1];
		adj[i].position = qaws_v3(b[2], b[3], b[4]);
	}
	vd = ss_views(0, fd, dir);
	vg = ss_views(0, fg, grad);
	vh = ss_views(0, fh, hv);
	s = ss_surface(0, x);
	memset(e, 0, sizeof(double) * 5);
	memset(grad, 0, sizeof(grad));
	memset(hv, 0, sizeof(hv));
	TEST_ASSERT_STATUS(qaws_surface_cdf_sample_tangent(NULL, s, r->measure, xi, NULL, SS_SAMPLES, cells, quadrature, &vd, val, t1, t2,
		&total));
	TEST_ASSERT_STATUS(qaws_surface_cdf_sample_adjoint(NULL, s, r->measure, xi, SS_SAMPLES, cells, quadrature, adj, &vg, NULL));
	TEST_ASSERT_STATUS(qaws_surface_cdf_sample_hvp(NULL, s, r->measure, xi, SS_SAMPLES, cells, quadrature, adj, &vd, &vh));
	e[0] = ref_err(total, r->total[0]);
	for (i = 0; i < SS_SAMPLES; i++)
	{
		e[0] = fmax(e[0], ref_sample_err(&val[i], r->value[i]));
		e[1] = fmax(e[1], ref_sample_err(&t1[i], r->tangent[i]));
		e[2] = fmax(e[2], ref_sample_err(&t2[i], r->tangent2[i]));
	}
	/* the reference splits the gradient and the HVP per sample */
	for (k = 0; k < 48; k++)
	{
		double gk = 0, hk = 0;
		for (i = 0; i < SS_SAMPLES; i++)
		{
			gk += r->grad[i][k];
			hk += r->hvp[i][k];
		}
		e[3] = fmax(e[3], ref_err(grad[k], gk));
		e[4] = fmax(e[4], ref_err(hv[k], hk));
	}
	qaws_surface_destroy(s);
}

static void test_reference(ss_ref const* r)
{
	static unsigned int const cells[3] = { 4, 8, 16 };
	double e[3][5], tol = QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-11;
	char msg[160];
	int c, k;
	for (c = 0; c < 3; c++)
		ref_errors(r, cells[c], 8, e[c]);
	printf("    reference %s, 8 x 8 Gauss points per cell (cells 4 / 8 / 16):\n", r->name);
	printf("      samples %.1e / %.1e / %.1e, tangents %.1e / %.1e / %.1e, tangent2 %.1e / %.1e / %.1e\n",
		e[0][0], e[1][0], e[2][0], e[0][1], e[1][1], e[2][1], e[0][2], e[1][2], e[2][2]);
	printf("      gradient %.1e / %.1e / %.1e, HVP %.1e / %.1e / %.1e\n", e[0][3], e[1][3], e[2][3], e[0][4], e[1][4], e[2][4]);
	for (k = 0; k < 5; k++)
	{
		static char const* const what[5] = { "total and samples", "tangents", "second tangents", "gradient",
			"Hessian-vector product" };
		sprintf(msg, "reference (%s, 16 x 16 cells): %s", r->name, what[k]);
		TEST_ASSERT(e[2][k] <= (k == 2 || k == 4 ? 10 : 1) * tol, msg);
	}
}

int test_59_diff_surface_sampling_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 59: Surface inverse-CDF sampling derivatives\n");
	test_flat();
	test_reference(&g_ss_refs[0]);
	test_reference(&g_ss_refs[1]);
	test_surface(0, NULL, "Bezier patch, area");
	test_surface(0, &g_ss_density, "Bezier patch, density");
	test_surface(1, NULL, "NURBS patch, area");
	test_surface(1, &g_ss_density, "NURBS patch, density");
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
