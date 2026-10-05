/*
 * Test 49: Differentiation model
 *
 * Context and report defaults, parameter keys, field names, views, and the
 * second-order forward / adjoint kernels of core/qaws_dual_core.h.
 */

#include "test_diff.h"

static void test_context_defaults(void)
{
	qaws_diff_context ctx;
	qaws_diff_report report;

	qaws_diff_context_init(&ctx);
	TEST_ASSERT(ctx.order == 1, "default order is 1");
	TEST_ASSERT(ctx.accumulation == QAWS_ACCUMULATE_SCATTER, "default accumulation is scatter");
	TEST_ASSERT(ctx.report == NULL, "no report by default");

	qaws_diff_report_reset(&report);
	TEST_ASSERT(report.validity == QAWS_DIFF_VALID, "report starts valid");
	TEST_ASSERT(report.diff_class == QAWS_DIFF_SMOOTH, "report starts smooth");
	TEST_ASSERT(report.evaluation_count == 0, "report starts empty");
}

static void test_param_keys(void)
{
	char buf[128];
	qaws_param_key k, parsed;
	unsigned int n;

	k = qaws_param_key_make(QAWS_FIELD_CONTROL_POINTS, 4, 2, 3);
	TEST_ASSERT_STATUS(qaws_param_key_prepend_child(&k, 0));
	n = qaws_param_key_to_string(&k, buf, sizeof(buf));
	TEST_ASSERT(strcmp(buf, "child[0]/control_points/4/z") == 0, "vector key text");
	TEST_ASSERT(n == strlen(buf), "key length");

	TEST_ASSERT_STATUS(qaws_param_key_parse(buf, &parsed));
	TEST_ASSERT(qaws_param_key_compare(&k, &parsed) == 0, "vector key round trip");

	k = qaws_param_key_make(QAWS_FIELD_WEIGHTS, 17, 0, 1);
	qaws_param_key_to_string(&k, buf, sizeof(buf));
	TEST_ASSERT(strcmp(buf, "weights/17") == 0, "scalar key has no component");
	TEST_ASSERT_STATUS(qaws_param_key_parse(buf, &parsed));
	TEST_ASSERT(qaws_param_key_compare(&k, &parsed) == 0, "scalar key round trip");

	k = qaws_param_key_make(QAWS_FIELD_RADIUS, 0, 0, 1);
	TEST_ASSERT_STATUS(qaws_param_key_prepend_child(&k, 3));
	TEST_ASSERT_STATUS(qaws_param_key_prepend_child(&k, 1));
	qaws_param_key_to_string(&k, buf, sizeof(buf));
	TEST_ASSERT(strcmp(buf, "child[1]/child[3]/radius/0") == 0, "nested child path order");

	TEST_ASSERT(qaws_param_key_parse("bogus/1", &parsed) != QAWS_STATUS_OK, "unknown field rejected");
	TEST_ASSERT(qaws_param_key_parse("control_points", &parsed) != QAWS_STATUS_OK, "missing element rejected");
	TEST_ASSERT(qaws_param_key_parse("control_points/1/w", &parsed) != QAWS_STATUS_OK, "bad component rejected");

	/* Truncation still reports the full length. */
	k = qaws_param_key_make(QAWS_FIELD_CONTROL_POINTS, 123, 1, 3);
	n = qaws_param_key_to_string(&k, buf, 8);
	TEST_ASSERT(n == strlen("control_points/123/y"), "truncated length");
	TEST_ASSERT(strlen(buf) == 7, "truncated buffer is terminated");

	{
		qaws_param_key a = qaws_param_key_make(QAWS_FIELD_CONTROL_POINTS, 1, 0, 3);
		qaws_param_key b = qaws_param_key_make(QAWS_FIELD_CONTROL_POINTS, 2, 0, 3);
		TEST_ASSERT(qaws_param_key_compare(&a, &b) < 0, "keys order by element");
	}
}

static void test_field_names(void)
{
	unsigned int i;
	int ok = 1;
	for (i = 1; i < (unsigned int)QAWS_FIELD_COUNT; i++)
		if (qaws_diff_field_from_name(qaws_diff_field_name((qaws_diff_field)i)) != (qaws_diff_field)i)
			ok = 0;
	TEST_ASSERT(ok, "every field name round trips");
	TEST_ASSERT(qaws_diff_field_from_name("nope") == QAWS_FIELD_NONE, "unknown name");
}

static void test_views(void)
{
	qaws_scalar data[6] = { 1, 2, 3, 4, 5, 6 };
	unsigned char active[2] = { 1, 0 };
	qaws_field_view v = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, data, 2, 3);
	qaws_diff_views views;

	views.fields = &v;
	views.field_count = 1;
	views.children = NULL;
	views.child_count = 0;

	TEST_ASSERT(qaws_diff_views_find(&views, QAWS_FIELD_CONTROL_POINTS) == &v, "find present field");
	TEST_ASSERT(qaws_diff_views_find(&views, QAWS_FIELD_WEIGHTS) == NULL, "absent field is inactive");

	v.active = active;
	v.component_mask = 1u | 4u;
	qaws_diff_views_clear(&views);
	TEST_ASSERT(data[0] == 0 && data[1] == 2 && data[2] == 0, "clear respects component mask");
	TEST_ASSERT(data[3] == 4 && data[4] == 5 && data[5] == 6, "clear respects element mask");
}

/* ------------------------------------------------------------------ */
/*  Dual kernels                                                      */
/* ------------------------------------------------------------------ */

static qaws_dual3 dual3_line(qaws_vec3 v, qaws_vec3 t, qaws_vec3 tt)
{
	return qaws_dual3_make(v, t, tt);
}

/* Evaluate the kernel at x + e*t + e^2/2*tt for a scalar e (primal only). */
typedef qaws_vec3 (*vec_fn)(qaws_vec3 a, qaws_vec3 b);

static qaws_vec3 fn_cross(qaws_vec3 a, qaws_vec3 b) { return qaws_v3_cross(a, b); }
static qaws_vec3 fn_normalize(qaws_vec3 a, qaws_vec3 b)
{
	qaws_scalar l = QAWS_SQRT(qaws_v3_dot(a, a));
	(void)b;
	return qaws_v3_scale(a, QAWS_ONE / l);
}
static qaws_vec3 fn_div(qaws_vec3 a, qaws_vec3 b)
{
	/* b.x acts as the scalar denominator */
	return qaws_v3_scale(a, QAWS_ONE / b.x);
}

static qaws_vec3 path(qaws_vec3 v, qaws_vec3 t, qaws_vec3 tt, double e)
{
	qaws_vec3 r;
	r.x = (qaws_scalar)(v.x + e * t.x + 0.5 * e * e * tt.x);
	r.y = (qaws_scalar)(v.y + e * t.y + 0.5 * e * e * tt.y);
	r.z = (qaws_scalar)(v.z + e * t.z + 0.5 * e * e * tt.z);
	return r;
}

/* Second-order expansion error at step e. */
static double taylor_error(vec_fn f, qaws_dual3 a, qaws_dual3 b, qaws_dual3 out, double e)
{
	qaws_vec3 y = f(path(a.v, a.t, a.tt, e), path(b.v, b.t, b.tt, e));
	double ex = y.x - (out.v.x + e * out.t.x + 0.5 * e * e * out.tt.x);
	double ey = y.y - (out.v.y + e * out.t.y + 0.5 * e * e * out.tt.y);
	double ez = y.z - (out.v.z + e * out.t.z + 0.5 * e * e * out.tt.z);
	return sqrt(ex * ex + ey * ey + ez * ez);
}

static void check_taylor(char const* name, vec_fn f, qaws_dual3 a, qaws_dual3 b, qaws_dual3 out)
{
#if QAWS_SCALAR_IS_FLOAT
	/* f32 cannot resolve O(e^3); check the first-order term only. */
	double e = 1e-2;
	qaws_vec3 yp = f(path(a.v, a.t, a.tt, e), path(b.v, b.t, b.tt, e));
	qaws_vec3 ym = f(path(a.v, a.t, a.tt, -e), path(b.v, b.t, b.tt, -e));
	int ok = diff_close((yp.x - ym.x) / (2 * e), out.t.x, 5e-2) &&
	         diff_close((yp.y - ym.y) / (2 * e), out.t.y, 5e-2) &&
	         diff_close((yp.z - ym.z) / (2 * e), out.t.z, 5e-2);
	TEST_ASSERT(ok, name);
#else
	double e1 = 1e-2, e2 = 5e-3;
	double r1 = taylor_error(f, a, b, out, e1);
	double r2 = taylor_error(f, a, b, out, e2);
	/* O(e^3): halving e divides the error by ~8. */
	int ok = (r1 < 1e-12) || (r1 / (r2 > 1e-300 ? r2 : 1e-300) > 6.0);
	if (!ok)
		printf("    %s: taylor errors %.3e %.3e\n", name, r1, r2);
	TEST_ASSERT(ok, name);
#endif
}

static void test_dual_kernels(void)
{
	qaws_dual3 a, b, out;
	qaws_dual1 w;
	diff_seed(7);

	a = dual3_line(qaws_v3(1, 2, 0.5f), diff_rand_vec3(), diff_rand_vec3());
	b = dual3_line(qaws_v3(-0.5f, 1, 2), diff_rand_vec3(), diff_rand_vec3());

	out = qaws_dual3_cross(a, b);
	check_taylor("cross second-order expansion", fn_cross, a, b, out);

	out = qaws_dual3_normalize(a);
	check_taylor("normalize second-order expansion", fn_normalize, a, b, out);

	w = qaws_dual1_make(b.v.x, b.t.x, b.tt.x);
	b.v.x = 1.7f;
	w.v = b.v.x;
	out = qaws_dual3_div(a, w);
	check_taylor("divide second-order expansion", fn_div, a, b, out);

	{
		/* sqrt and length against closed form */
		qaws_dual1 x = qaws_dual1_make(4, 1, 0);
		qaws_dual1 s = qaws_dual1_sqrt(x);
		TEST_ASSERT(approx_eq(s.v, 2) && approx_eq(s.t, 0.25f) && approx_eq(s.tt, -1.0f / 32.0f), "sqrt rule");
	}

	{
		/* degenerate normalize returns zero, never NaN */
		qaws_dual3 z = qaws_dual3_const(qaws_v3(0, 0, 0));
		qaws_dual3 n = qaws_dual3_normalize(z);
		TEST_ASSERT(n.v.x == 0 && n.t.x == 0 && n.tt.x == 0, "degenerate normalize is zero");
	}
}

/* <y_bar, J x_dot> == <J^T y_bar, x_dot> for each adjoint kernel */
static void test_adjoint_kernels(void)
{
	int i, ok_cross = 1, ok_norm = 1, ok_div = 1, ok_dot = 1, ok_len = 1, ok_scale = 1;
	diff_seed(99);
	for (i = 0; i < 20; i++)
	{
		qaws_vec3 a = diff_rand_vec3(), b = diff_rand_vec3();
		qaws_vec3 da = diff_rand_vec3(), db = diff_rand_vec3();
		qaws_vec3 ybar = diff_rand_vec3();
		qaws_scalar sbar = diff_rand();
		qaws_scalar w = (qaws_scalar)1.5 + diff_rand() * (qaws_scalar)0.5;
		qaws_scalar dw = diff_rand();
		double lhs, rhs;

		{
			qaws_dual3 c = qaws_dual3_cross(qaws_dual3_make(a, da, qaws_v3_zero()), qaws_dual3_make(b, db, qaws_v3_zero()));
			qaws_vec3_pair g = qaws_cross_adjoint(a, b, ybar);
			lhs = qaws_v3_dot(ybar, c.t);
			rhs = (double)qaws_v3_dot(g.a, da) + qaws_v3_dot(g.b, db);
			if (!diff_close(lhs, rhs, DIFF_TOL)) ok_cross = 0;
		}
		{
			qaws_dual3 n = qaws_dual3_normalize(qaws_dual3_make(a, da, qaws_v3_zero()));
			qaws_vec3 g = qaws_normalize_adjoint(a, ybar);
			lhs = qaws_v3_dot(ybar, n.t);
			rhs = qaws_v3_dot(g, da);
			if (!diff_close(lhs, rhs, DIFF_TOL)) ok_norm = 0;
		}
		{
			qaws_dual3 p = qaws_dual3_div(qaws_dual3_make(a, da, qaws_v3_zero()), qaws_dual1_make(w, dw, 0));
			qaws_vec3_scalar g = qaws_div_adjoint(a, w, ybar);
			lhs = qaws_v3_dot(ybar, p.t);
			rhs = (double)qaws_v3_dot(g.x, da) + (double)g.s * dw;
			if (!diff_close(lhs, rhs, DIFF_TOL)) ok_div = 0;
		}
		{
			qaws_dual1 s = qaws_dual3_dot(qaws_dual3_make(a, da, qaws_v3_zero()), qaws_dual3_make(b, db, qaws_v3_zero()));
			qaws_vec3_pair g = qaws_dot_adjoint(a, b, sbar);
			lhs = (double)sbar * s.t;
			rhs = (double)qaws_v3_dot(g.a, da) + qaws_v3_dot(g.b, db);
			if (!diff_close(lhs, rhs, DIFF_TOL)) ok_dot = 0;
		}
		{
			qaws_dual1 l = qaws_dual3_length(qaws_dual3_make(a, da, qaws_v3_zero()));
			qaws_vec3 g = qaws_length_adjoint(a, sbar);
			lhs = (double)sbar * l.t;
			rhs = qaws_v3_dot(g, da);
			if (!diff_close(lhs, rhs, DIFF_TOL)) ok_len = 0;
		}
		{
			qaws_dual3 y = qaws_dual3_scale(qaws_dual3_make(a, da, qaws_v3_zero()), qaws_dual1_make(w, dw, 0));
			qaws_vec3_scalar g = qaws_scale_adjoint(a, w, ybar);
			lhs = qaws_v3_dot(ybar, y.t);
			rhs = (double)qaws_v3_dot(g.x, da) + (double)g.s * dw;
			if (!diff_close(lhs, rhs, DIFF_TOL)) ok_scale = 0;
		}
	}
	TEST_ASSERT(ok_cross, "cross adjoint identity");
	TEST_ASSERT(ok_norm, "normalize adjoint identity");
	TEST_ASSERT(ok_div, "divide adjoint identity");
	TEST_ASSERT(ok_dot, "dot adjoint identity");
	TEST_ASSERT(ok_len, "length adjoint identity");
	TEST_ASSERT(ok_scale, "scale adjoint identity");
}

int test_49_diff_model_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 49: Differentiation model\n");
	test_context_defaults();
	test_param_keys();
	test_field_names();
	test_views();
	test_dual_kernels();
	test_adjoint_kernels();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
