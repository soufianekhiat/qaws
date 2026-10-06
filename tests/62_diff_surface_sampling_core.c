/*
 * Test 62: Backend-neutral surface sampling kernels
 *
 * core/qaws_bspline_surface_sampling_core.h compiled on the C backend and
 * driven like a GPU pass (one call per u cell, per sample, per control
 * point) must reproduce qaws_surface_cdf_sample_* with the area measure:
 * total, samples, first and second tangents along control point and xi
 * tangents, xi adjoints, and the control point adjoint as a gather.
 */

#include "test_diff.h"
#include "qaws_diff_sampling.h"
#include "core/qaws_bspline_surface_sampling_core.h"

#define SW_NU 5
#define SW_NV 6
#define SW_DEG 3
#define SW_CELLS 4
#define SW_N 4

static qaws_scalar const g_sw_uk[SW_NU + SW_DEG + 1] = { 0, 0, 0, 0, 0.45f, 1, 1, 1, 1 };
static qaws_scalar const g_sw_vk[SW_NV + SW_DEG + 1] = { 0, 0, 0, 0, 0.3f, 0.6f, 1, 1, 1, 1 };
static qaws_scalar const g_sw_xi[2 * SW_N] = { 0.2f, 0.7f, 0.55f, 0.35f, 0.9f, 0.85f, 0.4f, 0.1f };

typedef struct sw_patch
{
	qaws_scalar cp[QAWS_CORE_SURFACE_GRID], dot[QAWS_CORE_SURFACE_GRID], zero[QAWS_CORE_SURFACE_GRID];
	qaws_scalar uk[QAWS_CORE_MAX_POINTS * 2], vk[QAWS_CORE_MAX_POINTS * 2];
} sw_patch;

#define SW_ARGS(P, D) (P)->cp, (D), (P)->uk, (P)->vk, SW_DEG, SW_DEG, SW_NU, SW_NV

static void test_surface_sampling_core(void)
{
	static sw_patch P;
	qaws_scalar cps[SW_NU * SW_NV * 3], dir[SW_NU * SW_NV * 3], grad[SW_NU * SW_NV * 3], xd[2 * SW_N], xadj[2 * SW_N];
	qaws_scalar pv[QAWS_CORE_MAX_CELLS + 1], pt[QAWS_CORE_MAX_CELLS + 1], ptt[QAWS_CORE_MAX_CELLS + 1], total = 0;
	qaws_surface_cdf_sample val[SW_N], t1[SW_N], t2[SW_N], adj[SW_N];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	qaws_field_view fd, fg;
	qaws_diff_views vd, vg;
	double tol = QAWS_SCALAR_IS_FLOAT ? 3e-3 : 1e-9, lam[SW_N], mu[SW_N], su[SW_N], sv[SW_N], cellw[SW_CELLS], gn = 0, kn = 0;
	int ok_s = 1, ok_1 = 1, ok_2 = 1, ok_x = 1, ok_g = 1, i, j, c, q, n;
	char msg[200];

	diff_seed(6262u);
	memset(&P, 0, sizeof(P));
	for (i = 0; i < SW_NU; i++)
		for (j = 0; j < SW_NV; j++)
		{
			qaws_scalar* p = &cps[(i * SW_NV + j) * 3];
			p[0] = (qaws_scalar)(0.8 * i + 0.15 * j * j / 5.0);
			p[1] = (qaws_scalar)(0.5 * j * j / 5.0 + 0.2 * j);
			p[2] = (qaws_scalar)(0.3 * sin(1.2 * i + 0.7 * j));
		}
	diff_rand_fill(dir, SW_NU * SW_NV * 3);
	diff_rand_fill(xd, 2 * SW_N);
	for (i = 0; i < SW_NU * SW_NV * 3; i++)
	{
		P.cp[i] = cps[i];
		P.dot[i] = dir[i];
	}
	memcpy(P.uk, g_sw_uk, sizeof(g_sw_uk));
	memcpy(P.vk, g_sw_vk, sizeof(g_sw_vk));

	memset(&d, 0, sizeof(d));
	d.u_degree = SW_DEG;
	d.v_degree = SW_DEG;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = SW_NU;
	d.v_point_count = SW_NV;
	d.u_knots = g_sw_uk;
	d.u_knot_count = SW_NU + SW_DEG + 1;
	d.v_knots = g_sw_vk;
	d.v_knot_count = SW_NV + SW_DEG + 1;
	TEST_ASSERT_STATUS(qaws_surface_create_bspline(&d, &s));
	fd = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, SW_NU * SW_NV, 3);
	fg = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, grad, SW_NU * SW_NV, 3);
	vd.fields = &fd;
	vg.fields = &fg;
	vd.field_count = vg.field_count = 1;
	vd.children = vg.children = NULL;
	vd.child_count = vg.child_count = 0;
	TEST_ASSERT_STATUS(qaws_surface_cdf_sample_tangent(NULL, s, NULL, g_sw_xi, xd, SW_N, SW_CELLS, 8, &vd, val, t1, t2, &total));

	/* one thread per u cell, then an exclusive prefix sum */
	pv[0] = pt[0] = ptt[0] = 0;
	for (c = 0; c < SW_CELLS; c++)
	{
		qaws_dual1 m = qaws_patch_cell_mass(SW_ARGS(&P, P.dot), QAWS_ONE, SW_CELLS, c);
		pv[c + 1] = pv[c] + m.v;
		pt[c + 1] = pt[c] + m.t;
		ptt[c + 1] = ptt[c] + m.tt;
	}
	sprintf(msg, "cell kernels give the runtime area (%.12g vs %.12g)", (double)pv[SW_CELLS], (double)total);
	TEST_ASSERT(diff_close(pv[SW_CELLS], total, tol), msg);

	/* one thread per sample */
	for (n = 0; n < SW_N; n++)
	{
		qaws_scalar uv[2], tg[4];
		qaws_vec3 tp[2];
		qaws_patch_cdf_solve(SW_ARGS(&P, P.dot), SW_CELLS, pv, g_sw_xi[2 * n], g_sw_xi[2 * n + 1], uv);
		su[n] = uv[0];
		sv[n] = uv[1];
		ok_s &= diff_close(uv[0], val[n].u, QAWS_SCALAR_IS_FLOAT ? 1e-3 : 1e-11) && diff_close(uv[1], val[n].v, QAWS_SCALAR_IS_FLOAT ? 1e-3 : 1e-11);
		qaws_patch_cdf_tangent(SW_ARGS(&P, P.dot), SW_CELLS, pv, pt, ptt, g_sw_xi[2 * n], g_sw_xi[2 * n + 1], xd[2 * n], xd[2 * n + 1],
			uv[0], uv[1], tg, tp);
		ok_1 &= diff_close(tg[0], t1[n].u, tol) && diff_close(tg[1], t1[n].v, tol) && diff_close(tp[0].x, t1[n].position.x, tol) &&
			diff_close(tp[0].y, t1[n].position.y, tol) && diff_close(tp[0].z, t1[n].position.z, tol);
		ok_2 &= diff_close(tg[2], t2[n].u, tol) && diff_close(tg[3], t2[n].v, tol) && diff_close(tp[1].x, t2[n].position.x, tol) &&
			diff_close(tp[1].y, t2[n].position.y, tol) && diff_close(tp[1].z, t2[n].position.z, tol);
	}
	TEST_ASSERT(ok_s, "solve kernel gives the runtime samples");
	TEST_ASSERT(ok_1, "first order tangent kernel matches the runtime (control points and xi)");
	TEST_ASSERT(ok_2, "second order tangent kernel matches the runtime (control points and xi)");

	/* adjoint: multipliers per sample, then a gather per control point */
	for (n = 0; n < SW_N; n++)
	{
		adj[n].u = diff_rand();
		adj[n].v = diff_rand();
		adj[n].position = diff_rand_vec3();
	}
	memset(grad, 0, sizeof(grad));
	memset(xadj, 0, sizeof(xadj));
	TEST_ASSERT_STATUS(qaws_surface_cdf_sample_adjoint(NULL, s, NULL, g_sw_xi, SW_N, SW_CELLS, 8, adj, &vg, xadj));
	for (c = 0; c < SW_CELLS; c++)
		cellw[c] = 0;
	for (n = 0; n < SW_N; n++)
	{
		qaws_scalar out[4];
		int ku;
		qaws_patch_cdf_multipliers(SW_ARGS(&P, P.zero), SW_CELLS, pv, g_sw_xi[2 * n], g_sw_xi[2 * n + 1], (qaws_scalar)su[n],
			(qaws_scalar)sv[n], adj[n].u, adj[n].v, adj[n].position, out);
		lam[n] = out[0];
		mu[n] = out[1];
		ok_x &= diff_close(out[2], xadj[2 * n], tol) && diff_close(out[3], xadj[2 * n + 1], tol);
		ku = (int)(su[n] * SW_CELLS);
		if (ku >= SW_CELLS) ku = SW_CELLS - 1;
		for (c = 0; c < SW_CELLS; c++)
			cellw[c] += lam[n] * g_sw_xi[2 * n] - (c < ku ? lam[n] : 0);
	}
	TEST_ASSERT(ok_x, "multiplier kernel gives the runtime xi adjoints");
	for (i = 0; i < SW_NU; i++)
		for (j = 0; j < SW_NV; j++)
		{
			qaws_vec3 acc = qaws_v3_zero();
			int k = (i * SW_NV + j) * 3;
			for (n = 0; n < SW_N; n++)
			{
				int ku = (int)(su[n] * SW_CELLS);
				double uk;
				qaws_scalar ders_u[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS], ders_v[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
				int us = qaws_find_span(P.uk, SW_DEG, SW_NU, (qaws_scalar)su[n]), vs = qaws_find_span(P.vk, SW_DEG, SW_NV, (qaws_scalar)sv[n]);
				if (ku >= SW_CELLS) ku = SW_CELLS - 1;
				uk = (double)ku / SW_CELLS;
				/* the position term N_ij(u, v) p_bar */
				if (i >= us - SW_DEG && i <= us && j >= vs - SW_DEG && j <= vs)
				{
					double w;
					qaws_patch_rows(P.uk, SW_DEG, us, (qaws_scalar)su[n], ders_u);
					qaws_patch_rows(P.vk, SW_DEG, vs, (qaws_scalar)sv[n], ders_v);
					w = ders_u[i - (us - SW_DEG)] * ders_v[j - (vs - SW_DEG)];
					acc = qaws_v3_axpy(acc, adj[n].position, (qaws_scalar)w);
				}
				/* conditional lines */
				acc = qaws_v3_add(acc, qaws_patch_line_adjoint_cp(SW_ARGS(&P, P.zero), SW_CELLS, (qaws_scalar)su[n], (qaws_scalar)sv[n], i, j,
					(qaws_scalar)(-mu[n])));
				acc = qaws_v3_add(acc, qaws_patch_line_adjoint_cp(SW_ARGS(&P, P.zero), SW_CELLS, (qaws_scalar)su[n], 1, i, j,
					(qaws_scalar)(mu[n] * g_sw_xi[2 * n + 1])));
				/* lines of the partial u cell */
				for (q = 0; q < 8; q++)
				{
					qaws_scalar xw[2];
					double len = su[n] - uk;
					qaws_gauss8(q, xw);
					acc = qaws_v3_add(acc, qaws_patch_line_adjoint_cp(SW_ARGS(&P, P.zero), SW_CELLS, (qaws_scalar)(uk + 0.5 * (1 + xw[0]) * len), 1,
						i, j, (qaws_scalar)(-lam[n] * 0.5 * xw[1] * len)));
				}
			}
			/* lines of the full u cells */
			for (c = 0; c < SW_CELLS; c++)
				for (q = 0; q < 8; q++)
				{
					qaws_scalar xw[2];
					double hu = 1.0 / SW_CELLS;
					qaws_gauss8(q, xw);
					acc = qaws_v3_add(acc, qaws_patch_line_adjoint_cp(SW_ARGS(&P, P.zero), SW_CELLS, (qaws_scalar)(hu * (c + 0.5 + 0.5 * xw[0])), 1,
						i, j, (qaws_scalar)(cellw[c] * 0.5 * hu * xw[1])));
				}
			ok_g &= diff_close(acc.x, grad[k], tol) && diff_close(acc.y, grad[k + 1], tol) && diff_close(acc.z, grad[k + 2], tol);
			gn += (double)grad[k] * grad[k] + (double)grad[k + 1] * grad[k + 1] + (double)grad[k + 2] * grad[k + 2];
			kn += (double)acc.x * acc.x + (double)acc.y * acc.y + (double)acc.z * acc.z;
		}
	printf("    area %.15g (runtime %.15g), |g| %.15g (runtime %.15g)\n", (double)pv[SW_CELLS], (double)total, sqrt(kn), sqrt(gn));
	TEST_ASSERT(ok_g, "gather over control points gives the runtime adjoint");
	qaws_surface_destroy(s);
}

int test_62_diff_surface_sampling_core_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 62: Backend-neutral surface sampling kernels\n");
	test_surface_sampling_core();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
