/*
 * Test 60: Backend-neutral sampling and functional kernels
 *
 * core/qaws_bspline_sampling_core.h compiled on the C backend, driven the
 * way a GPU pass would drive it (one call per span, per sample, per control
 * point), must reproduce the runtime:
 *   - span measures, totals and the solved sample parameters
 *   - first and second sample tangents along control point tangents and
 *     target rates (arc length and curvature measures)
 *   - the sampling adjoint as a gather over control points
 *   - length, bending and curvature-squared functionals: value, tangents,
 *     gradient
 */

#include "test_diff.h"
#include "qaws_diff_sampling.h"
#include "qaws_diff_functionals.h"
#include "core/qaws_bspline_sampling_core.h"

#define SC_CP 9
#define SC_DEG 3
#define SC_KC (SC_CP + SC_DEG + 1)
#define SC_SPANS 6
#define SC_N 7

static qaws_scalar const g_sc_knots[SC_KC] = { 0, 0, 0, 0, 0.7f, 1.1f, 2.0f, 2.4f, 3.3f, 4, 4, 4, 4 };

/* Span k (0..SC_SPANS-1) starts at knot SC_DEG + k. */
typedef struct sc_window
{
	qaws_scalar lk[QAWS_CORE_MAX_POINTS * 2], lcp[QAWS_CORE_MAX_POINTS * 3], ldir[QAWS_CORE_MAX_POINTS * 3];
	qaws_scalar a, b;
	int first;   /* global index of local control point 0 */
} sc_window;

static void sc_load(qaws_scalar const* cps, qaws_scalar const* dir, int k, sc_window* w)
{
	int span = SC_DEG + k, j, c;
	memset(w, 0, sizeof(*w));
	for (j = 0; j < 2 * (SC_DEG + 1); j++)
		w->lk[j] = g_sc_knots[span - SC_DEG + j];
	for (j = 0; j <= SC_DEG; j++)
		for (c = 0; c < 3; c++)
		{
			w->lcp[j * 3 + c] = cps[(span - SC_DEG + j) * 3 + c];
			w->ldir[j * 3 + c] = dir ? dir[(span - SC_DEG + j) * 3 + c] : 0;
		}
	w->a = g_sc_knots[span];
	w->b = g_sc_knots[span + 1];
	w->first = span - SC_DEG;
}

static qaws_curve* sc_curve(qaws_scalar const* cps)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = SC_DEG;
	d.control_points = cps;
	d.control_point_count = SC_CP;
	d.knots = g_sc_knots;
	d.knot_count = SC_KC;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static int sc_vclose(qaws_vec3 a, qaws_vec3 b, double tol)
{
	return diff_close(a.x, b.x, tol) && diff_close(a.y, b.y, tol) && diff_close(a.z, b.z, tol);
}

static qaws_diff_views sc_views(qaws_field_view* fv, qaws_scalar* data)
{
	qaws_diff_views v;
	*fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, data, SC_CP, 3);
	v.fields = fv;
	v.field_count = 1;
	v.children = NULL;
	v.child_count = 0;
	return v;
}

static void test_sampling(int id, qaws_sample_measure_desc const* measure, char const* name)
{
	static double const dist[SC_N] = { 0, 0.4, 1.1, -0.3, 0, 2.0, 0.05 };
	static double const frac[SC_N] = { 0, 0.1, 0.25, 0.5, 0.77, 0.3, 1 };
	qaws_scalar cps[SC_CP * 3], dir[SC_CP * 3], grad[SC_CP * 3], gat[SC_CP * 3], dt[SC_N], dt2[SC_N], dadj[SC_N];
	qaws_scalar flo = measure ? measure->curvature_floor : 0;
	qaws_cdf_target tg[SC_N];
	qaws_cdf_sample val[SC_N], t1[SC_N], t2[SC_N], adj[SC_N];
	qaws_scalar total = 0;
	qaws_field_view fd, fg;
	qaws_diff_views vd, vg;
	qaws_curve* c;
	sc_window w[SC_SPANS];
	double cum[SC_SPANS + 1], cum1[SC_SPANS + 1], cum2[SC_SPANS + 1], lambda[SC_N], wspan[SC_SPANS], F = 0;
	int span_of[SC_N];
	double ts[SC_N];
	double tol = QAWS_SCALAR_IS_FLOAT ? 2e-3 : 1e-10;
	int ok_tot, ok_t = 1, ok_1 = 1, ok_2 = 1, ok_g = 1, ok_d = 1, k, i, e, q;
	char msg[160];

	diff_seed(6060u + (unsigned int)id);
	diff_rand_fill(cps, SC_CP * 3);
	for (i = 0; i < SC_CP; i++)
		cps[3 * i] += (qaws_scalar)(1.5 * i);   /* a regular curve */
	diff_rand_fill(dir, SC_CP * 3);
	for (i = 0; i < SC_N; i++)
	{
		tg[i].distance = (qaws_scalar)dist[i];
		tg[i].fraction = (qaws_scalar)frac[i];
		dt[i] = diff_rand();
		dt2[i] = diff_rand();
		adj[i].t = diff_rand();
		adj[i].position = diff_rand_vec3();
	}
	c = sc_curve(cps);
	vd = sc_views(&fd, dir);
	vg = sc_views(&fg, grad);
	TEST_ASSERT_STATUS(qaws_curve_cdf_sample_tangent(NULL, c, measure, tg, dt, dt2, SC_N, 0, &vd, val, t1, t2, &total));

	/* one thread per span: measures and their tangents, then a prefix sum */
	cum[0] = cum1[0] = cum2[0] = 0;
	for (k = 0; k < SC_SPANS; k++)
	{
		qaws_dual1 m;
		sc_load(cps, dir, k, &w[k]);
		m = qaws_bspline_integrate(id, flo, w[k].lcp, w[k].ldir, w[k].lk, SC_DEG, SC_DEG, w[k].a, w[k].b);
		cum[k + 1] = cum[k] + m.v;
		cum1[k + 1] = cum1[k] + m.t;
		cum2[k + 1] = cum2[k] + m.tt;
	}
	ok_tot = diff_close(cum[SC_SPANS], total, tol);
	sprintf(msg, "%s: span kernels give the runtime total (%.12g vs %.12g)", name, cum[SC_SPANS], (double)total);
	TEST_ASSERT(ok_tot, msg);

	/* one thread per sample: solve, then tangents */
	for (i = 0; i < SC_N; i++)
	{
		double sigma = dist[i] + frac[i] * cum[SC_SPANS];
		qaws_dual1 part;
		qaws_scalar ot[2];
		qaws_vec3 op[2];
		if (sigma < 0) sigma = 0;
		if (sigma > cum[SC_SPANS]) sigma = cum[SC_SPANS];
		k = 0;
		while (k + 1 < SC_SPANS && cum[k + 1] <= sigma)
			k++;
		span_of[i] = k;
		ts[i] = qaws_bspline_cdf_solve(id, flo, w[k].lcp, w[k].lk, SC_DEG, SC_DEG, w[k].a, w[k].b, (qaws_scalar)(sigma - cum[k]));
		ok_t &= diff_close(ts[i], val[i].t, QAWS_SCALAR_IS_FLOAT ? 1e-3 : 1e-11);
		part = qaws_bspline_integrate(id, flo, w[k].lcp, w[k].ldir, w[k].lk, SC_DEG, SC_DEG, w[k].a, (qaws_scalar)ts[i]);
		qaws_bspline_cdf_tangent(id, flo, w[k].lcp, w[k].ldir, w[k].lk, SC_DEG, SC_DEG, (qaws_scalar)ts[i],
			(qaws_scalar)(dt[i] + frac[i] * cum1[SC_SPANS]), (qaws_scalar)(dt2[i] + frac[i] * cum2[SC_SPANS]),
			(qaws_scalar)(cum1[k] + part.t), (qaws_scalar)(cum2[k] + part.tt), ot, op);
		ok_1 &= diff_close(ot[0], t1[i].t, tol) && sc_vclose(op[0], t1[i].position, tol);
		ok_2 &= diff_close(ot[1], t2[i].t, tol) && sc_vclose(op[1], t2[i].position, tol);
	}
	sprintf(msg, "%s: solve kernel gives the runtime samples", name);
	TEST_ASSERT(ok_t, msg);
	sprintf(msg, "%s: first order tangent kernel matches the runtime", name);
	TEST_ASSERT(ok_1, msg);
	sprintf(msg, "%s: second order tangent kernel matches the runtime", name);
	TEST_ASSERT(ok_2, msg);

	/* adjoint: lambda per sample, span weights, then a gather per control point */
	memset(grad, 0, sizeof(grad));
	memset(dadj, 0, sizeof(dadj));
	TEST_ASSERT_STATUS(qaws_curve_cdf_sample_adjoint(NULL, c, measure, tg, SC_N, 0, adj, &vg, dadj));
	for (i = 0; i < SC_N; i++)
	{
		k = span_of[i];
		lambda[i] = qaws_bspline_cdf_lambda(id, flo, w[k].lcp, w[k].lk, SC_DEG, SC_DEG, (qaws_scalar)ts[i], adj[i].t, adj[i].position);
		ok_d &= diff_close(lambda[i], dadj[i], tol);
		F += lambda[i] * frac[i];
	}
	for (k = 0; k < SC_SPANS; k++)
	{
		wspan[k] = F;
		for (i = 0; i < SC_N; i++)
			if (span_of[i] > k)
				wspan[k] -= lambda[i];
	}
	for (e = 0; e < SC_CP; e++)
	{
		qaws_vec3 acc = qaws_v3_zero();
		for (k = 0; k < SC_SPANS; k++)
		{
			int j = e - w[k].first;
			double h = ((double)w[k].b - w[k].a) / QAWS_SAMPLING_PIECES;
			if (j < 0 || j > SC_DEG)
				continue;
			/* full span nodes */
			for (q = 0; q < 8 * QAWS_SAMPLING_PIECES; q++)
			{
				qaws_scalar xw[2];
				qaws_gauss8(q % 8, xw);
				acc = qaws_v3_add(acc, qaws_bspline_node_adjoint_cp(id, flo, w[k].lcp, w[k].lk, SC_DEG, SC_DEG,
					(qaws_scalar)(w[k].a + h * (q / 8) + 0.5 * h * (1 + xw[0])), j, (qaws_scalar)(wspan[k] * 0.5 * h * xw[1])));
			}
			/* partial spans and positions of the samples in this span */
			for (i = 0; i < SC_N; i++)
			{
				qaws_eval_3d yb;
				double hp = (ts[i] - w[k].a) / QAWS_SAMPLING_PIECES;
				if (span_of[i] != k)
					continue;
				for (q = 0; q < 8 * QAWS_SAMPLING_PIECES && hp > 0; q++)
				{
					qaws_scalar xw[2];
					qaws_gauss8(q % 8, xw);
					acc = qaws_v3_add(acc, qaws_bspline_node_adjoint_cp(id, flo, w[k].lcp, w[k].lk, SC_DEG, SC_DEG,
						(qaws_scalar)(w[k].a + hp * (q / 8) + 0.5 * hp * (1 + xw[0])), j, (qaws_scalar)(-lambda[i] * 0.5 * hp * xw[1])));
				}
				yb.position = adj[i].position;
				yb.d1 = qaws_v3_zero();
				yb.d2 = qaws_v3_zero();
				acc = qaws_v3_add(acc, qaws_bspline_adjoint_cp_3d(w[k].lk, SC_DEG, SC_DEG, (qaws_scalar)ts[i], j, yb));
			}
		}
		gat[3 * e] = acc.x;
		gat[3 * e + 1] = acc.y;
		gat[3 * e + 2] = acc.z;
		ok_g &= sc_vclose(acc, qaws_v3(grad[3 * e], grad[3 * e + 1], grad[3 * e + 2]), QAWS_SCALAR_IS_FLOAT ? 5e-3 : 1e-10);
	}
	sprintf(msg, "%s: lambda kernel gives the runtime distance adjoints", name);
	TEST_ASSERT(ok_d, msg);
	sprintf(msg, "%s: gather over control points gives the runtime adjoint (|g| %.6f vs %.6f)", name,
		sqrt(diff_dot(gat, gat, SC_CP * 3)), sqrt(diff_dot(grad, grad, SC_CP * 3)));
	printf("    %s: total %.15g (runtime %.15g), |g| %.15g (runtime %.15g)\n", name, cum[SC_SPANS], (double)total,
		sqrt(diff_dot(gat, gat, SC_CP * 3)), sqrt(diff_dot(grad, grad, SC_CP * 3)));
	TEST_ASSERT(ok_g, msg);
	qaws_curve_destroy(c);
}

static void test_functional(int id, qaws_curve_functional functional, char const* name)
{
	qaws_scalar cps[SC_CP * 3], dir[SC_CP * 3], grad[SC_CP * 3], value = 0, tan1 = 0, tan2 = 0, gv = 0;
	qaws_field_view fd, fg;
	qaws_diff_views vd, vg;
	qaws_curve* c;
	sc_window w[SC_SPANS];
	double v = 0, d1 = 0, d2 = 0, tol = QAWS_SCALAR_IS_FLOAT ? 2e-3 : 1e-10;
	int k, e, q, ok_g = 1;
	char msg[160];
	diff_seed(6161u);
	diff_rand_fill(cps, SC_CP * 3);
	for (e = 0; e < SC_CP; e++)
		cps[3 * e] += (qaws_scalar)(1.5 * e);
	diff_rand_fill(dir, SC_CP * 3);
	c = sc_curve(cps);
	vd = sc_views(&fd, dir);
	vg = sc_views(&fg, grad);
	memset(grad, 0, sizeof(grad));
	TEST_ASSERT_STATUS(qaws_curve_functional_eval(NULL, c, functional, 8, &vd, &value, &tan1, &tan2));
	TEST_ASSERT_STATUS(qaws_curve_functional_gradient(NULL, c, functional, 8, &vg, &gv));
	for (k = 0; k < SC_SPANS; k++)
	{
		qaws_dual1 m;
		sc_load(cps, dir, k, &w[k]);
		m = qaws_bspline_integrate(id, 0, w[k].lcp, w[k].ldir, w[k].lk, SC_DEG, SC_DEG, w[k].a, w[k].b);
		v += m.v;
		d1 += m.t;
		d2 += m.tt;
	}
	sprintf(msg, "%s: span kernels give the runtime value and tangents (%.12g vs %.12g)", name, v, (double)value);
	TEST_ASSERT(diff_close(v, value, tol) && diff_close(d1, tan1, tol) && diff_close(d2, tan2, tol), msg);
	for (e = 0; e < SC_CP; e++)
	{
		qaws_vec3 acc = qaws_v3_zero();
		for (k = 0; k < SC_SPANS; k++)
		{
			int j = e - w[k].first;
			double h = ((double)w[k].b - w[k].a) / QAWS_SAMPLING_PIECES;
			if (j < 0 || j > SC_DEG)
				continue;
			for (q = 0; q < 8 * QAWS_SAMPLING_PIECES; q++)
			{
				qaws_scalar xw[2];
				qaws_gauss8(q % 8, xw);
				acc = qaws_v3_add(acc, qaws_bspline_node_adjoint_cp(id, 0, w[k].lcp, w[k].lk, SC_DEG, SC_DEG,
					(qaws_scalar)(w[k].a + h * (q / 8) + 0.5 * h * (1 + xw[0])), j, (qaws_scalar)(0.5 * h * xw[1])));
			}
		}
		ok_g &= sc_vclose(acc, qaws_v3(grad[3 * e], grad[3 * e + 1], grad[3 * e + 2]), tol);
	}
	sprintf(msg, "%s: gather over control points gives the runtime gradient", name);
	TEST_ASSERT(ok_g, msg);
	qaws_curve_destroy(c);
}

int test_60_diff_sampling_core_main(void)
{
	static qaws_sample_measure_desc const curv = { QAWS_MEASURE_CURVATURE, (qaws_scalar)0.4, NULL, NULL };
	g_pass = 0;
	g_fail = 0;
	printf("Test 60: Backend-neutral sampling and functional kernels\n");
	test_sampling(QAWS_CORE_INTEGRAND_SPEED, NULL, "arc length");
	test_sampling(QAWS_CORE_INTEGRAND_CURVATURE, &curv, "curvature");
	test_functional(QAWS_CORE_INTEGRAND_SPEED, QAWS_FUNCTIONAL_LENGTH, "length");
	test_functional(QAWS_CORE_INTEGRAND_BENDING, QAWS_FUNCTIONAL_BENDING, "bending");
	test_functional(QAWS_CORE_INTEGRAND_CURVATURE_SQ, QAWS_FUNCTIONAL_CURVATURE_SQUARED, "curvature squared");
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
