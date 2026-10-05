/*
 * Test 57: Backend-neutral B-spline derivative kernels
 *
 * core/qaws_bspline_diff_core.h compiled on the C backend must reproduce
 * the runtime: tangent jets of qaws_curve_eval_batch_tangent_3d, and a
 * gather pass (one "thread" per control point over the samples listed by
 * qaws_curve_build_support_index) must reproduce the batch adjoint.
 */

#include "test_diff.h"
#include "core/qaws_bspline_diff_core.h"

#define CORE_CP 9
#define CORE_DEG 3
#define CORE_KC (CORE_CP + CORE_DEG + 1)
#define CORE_SAMPLES 40

static qaws_scalar const g_core_knots[CORE_KC] = { 0, 0, 0, 0, 0.7f, 1.1f, 2.0f, 2.4f, 3.3f, 4, 4, 4, 4 };

/* Local window of a sample, as a shader would load it. */
static int core_window(qaws_scalar const* cps, qaws_scalar const* dir, qaws_scalar t,
	qaws_scalar* lk, qaws_scalar* lcp, qaws_scalar* ldir)
{
	qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2];
	int span, j;
	memset(knots, 0, sizeof(knots));
	memcpy(knots, g_core_knots, sizeof(g_core_knots));
	span = qaws_find_span(knots, CORE_DEG, CORE_CP, t);
	for (j = 0; j < 2 * (CORE_DEG + 1); j++)
		lk[j] = g_core_knots[span - CORE_DEG + j];
	for (j = 0; j <= CORE_DEG; j++)
	{
		int c;
		for (c = 0; c < 3; c++)
		{
			lcp[j * 3 + c] = cps[(span - CORE_DEG + j) * 3 + c];
			if (ldir)
				ldir[j * 3 + c] = dir[(span - CORE_DEG + j) * 3 + c];
		}
	}
	return span;
}

static int vclose(qaws_vec3 a, qaws_vec3 b, double tol)
{
	double s = 1 + fabs(b.x) + fabs(b.y) + fabs(b.z);
	return fabs(a.x - b.x) <= tol * s && fabs(a.y - b.y) <= tol * s && fabs(a.z - b.z) <= tol * s;
}

static void test_core_kernels(void)
{
	qaws_scalar cps[CORE_CP * 3], dir[CORE_CP * 3], ts[CORE_SAMPLES], tdots[CORE_SAMPLES];
	qaws_scalar lk[QAWS_CORE_MAX_POINTS * 2], lcp[QAWS_CORE_MAX_POINTS * 3], ldir[QAWS_CORE_MAX_POINTS * 3];
	qaws_curve_jet_3d tan[CORE_SAMPLES], ybar[CORE_SAMPLES];
	qaws_scalar ref_bar[CORE_CP * 3], ref_tbar[CORE_SAMPLES], gat_bar[CORE_CP * 3], gat_tbar[CORE_SAMPLES];
	unsigned int offsets[CORE_CP + 1], samples[CORE_SAMPLES * (CORE_DEG + 1)], entries = 0, i, e;
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	qaws_field_view fv;
	qaws_diff_views views;
	int ok_tan = 1, ok_cp = 1, ok_t = 1;

	diff_seed(5757);
	diff_rand_fill(cps, CORE_CP * 3);
	diff_rand_fill(dir, CORE_CP * 3);
	for (i = 0; i < CORE_SAMPLES; i++)
	{
		ts[i] = (qaws_scalar)(4.0 * (i + 0.5) / CORE_SAMPLES);
		tdots[i] = diff_rand();
		diff_rand_jet3(&ybar[i]);
		ybar[i].d[3] = qaws_v3_zero();
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = CORE_DEG;
	d.control_points = cps;
	d.control_point_count = CORE_CP;
	d.knots = g_core_knots;
	d.knot_count = CORE_KC;
	TEST_ASSERT_STATUS(qaws_curve_create_bspline(&d, &c));

	/* tangents */
	fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, CORE_CP, 3);
	views.fields = &fv; views.field_count = 1; views.children = NULL; views.child_count = 0;
	TEST_ASSERT_STATUS(qaws_curve_eval_batch_tangent_3d(NULL, c, ts, tdots, CORE_SAMPLES, 0x7, &views, NULL, tan));
	for (i = 0; i < CORE_SAMPLES; i++)
	{
		int span = core_window(cps, dir, ts[i], lk, lcp, ldir);
		qaws_eval_3d k = qaws_bspline_tangent_3d(lcp, ldir, lk, CORE_DEG, CORE_DEG, ts[i], tdots[i]);
		(void)span;
		if (!vclose(k.position, tan[i].d[0], 1e-5) || !vclose(k.d1, tan[i].d[1], 1e-5) || !vclose(k.d2, tan[i].d[2], 1e-5))
			ok_tan = 0;
	}
	TEST_ASSERT(ok_tan, "core tangent kernel matches the runtime tangent");

	/* runtime adjoint */
	memset(ref_bar, 0, sizeof(ref_bar));
	memset(ref_tbar, 0, sizeof(ref_tbar));
	fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, ref_bar, CORE_CP, 3);
	TEST_ASSERT_STATUS(qaws_curve_eval_batch_adjoint_3d(NULL, c, ts, CORE_SAMPLES, 0x7, ybar, &views, ref_tbar));

	/* gather: one pass per control point over its samples */
	TEST_ASSERT_STATUS(qaws_curve_build_support_index(c, ts, CORE_SAMPLES, QAWS_FIELD_CONTROL_POINTS,
		offsets, CORE_CP + 1, samples, CORE_SAMPLES * (CORE_DEG + 1), &entries));
	memset(gat_bar, 0, sizeof(gat_bar));
	for (e = 0; e < CORE_CP; e++)
	{
		qaws_vec3 acc = qaws_v3_zero();
		unsigned int k;
		for (k = offsets[e]; k < offsets[e + 1]; k++)
		{
			unsigned int s = samples[k];
			qaws_eval_3d yb;
			int span = core_window(cps, NULL, ts[s], lk, lcp, NULL);
			yb.position = ybar[s].d[0];
			yb.d1 = ybar[s].d[1];
			yb.d2 = ybar[s].d[2];
			acc = qaws_v3_add(acc, qaws_bspline_adjoint_cp_3d(lk, CORE_DEG, CORE_DEG, ts[s], (int)e - (span - CORE_DEG), yb));
		}
		gat_bar[e * 3 + 0] = acc.x;
		gat_bar[e * 3 + 1] = acc.y;
		gat_bar[e * 3 + 2] = acc.z;
		if (!vclose(acc, qaws_v3(ref_bar[e * 3], ref_bar[e * 3 + 1], ref_bar[e * 3 + 2]), 1e-5))
			ok_cp = 0;
	}
	for (i = 0; i < CORE_SAMPLES; i++)
	{
		qaws_eval_3d yb;
		core_window(cps, NULL, ts[i], lk, lcp, NULL);
		yb.position = ybar[i].d[0];
		yb.d1 = ybar[i].d[1];
		yb.d2 = ybar[i].d[2];
		gat_tbar[i] = qaws_bspline_adjoint_t_3d(lcp, lk, CORE_DEG, CORE_DEG, ts[i], yb);
		if (!diff_close(gat_tbar[i], ref_tbar[i], 1e-5))
			ok_t = 0;
	}
	printf("    %u support entries for %u samples; gather |bar| %.6f, runtime %.6f\n", entries, CORE_SAMPLES,
		sqrt(diff_dot(gat_bar, gat_bar, CORE_CP * 3)), sqrt(diff_dot(ref_bar, ref_bar, CORE_CP * 3)));
	TEST_ASSERT(ok_cp, "gather kernel over the support index matches the runtime control point adjoint");
	TEST_ASSERT(ok_t, "core parameter adjoint matches the runtime");
	qaws_curve_destroy(c);
}

int test_57_diff_core_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 57: Backend-neutral derivative kernels\n");
	test_core_kernels();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
