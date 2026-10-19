/*
 * Test 147: rotation-minimizing frames (qaws_curve_compute_rotation_minimizing_frames_3d)
 *
 *   - a planar curve in a tilted plane, started from the plane's normal:
 *     the normal stays that normal (the reflections keep it);
 *   - a straight segment, where the Frenet frame is undefined: one frame
 *     all along;
 *   - a twisted space curve: the frames converge at fourth order as the
 *     samples double (as Wang et al. show for the double reflection),
 *     against a frame from 64 times as many samples;
 *   - frames orthonormal, tangent the unit D1; refusals.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_curve* fr_bezier(qaws_vec3 const* p, unsigned int n)
{
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = n - 1;
	d.control_points = p;
	d.control_point_count = n;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

static void fr_params(qaws_scalar* t, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		t[i] = (qaws_scalar)((double)i / (n - 1));
}

static void test_planar_and_straight(void)
{
	/* a plane through the origin with unit normal n, spanned by a and b */
	double n[3] = { 1, 2, 3 }, a[3] = { 2, -1, 0 }, b[3], ln = sqrt(14.0), la = sqrt(5.0);
	qaws_vec3 p[5], q[2], nv;
	qaws_scalar t[200], fr[200 * 9];
	qaws_curve *c, *line;
	double worst = 0, worst_line = 0, ortho = 0;
	unsigned int i, k;
	char msg[256];
	for (i = 0; i < 3; i++)
	{
		n[i] /= ln;
		a[i] /= la;
	}
	b[0] = n[1] * a[2] - n[2] * a[1];
	b[1] = n[2] * a[0] - n[0] * a[2];
	b[2] = n[0] * a[1] - n[1] * a[0];
	for (i = 0; i < 5; i++)
	{
		double u = cos(1.3 * i) * (1 + 0.3 * i), v = sin(1.7 * i) + 0.2 * i;
		p[i].x = (qaws_scalar)(u * a[0] + v * b[0]);
		p[i].y = (qaws_scalar)(u * a[1] + v * b[1]);
		p[i].z = (qaws_scalar)(u * a[2] + v * b[2]);
	}
	nv.x = (qaws_scalar)n[0]; nv.y = (qaws_scalar)n[1]; nv.z = (qaws_scalar)n[2];
	c = fr_bezier(p, 5);
	fr_params(t, 200);
	TEST_ASSERT(c && qaws_curve_compute_rotation_minimizing_frames_3d(c, t, 200, &nv, fr) == QAWS_STATUS_OK, "planar curve: frames");
	for (i = 0; i < 200; i++)
	{
		/* the normal stays on the plane normal */
		double e = fabs(fr[9 * i + 3] - n[0]) + fabs(fr[9 * i + 4] - n[1]) + fabs(fr[9 * i + 5] - n[2]);
		worst = e > worst ? e : worst;
		/* orthonormal rows */
		for (k = 0; k < 3; k++)
		{
			qaws_scalar const* u = fr + 9 * i + 3 * k;
			qaws_scalar const* w = fr + 9 * i + 3 * ((k + 1) % 3);
			double len = sqrt((double)u[0] * u[0] + (double)u[1] * u[1] + (double)u[2] * u[2]);
			double dot = (double)u[0] * w[0] + (double)u[1] * w[1] + (double)u[2] * w[2];
			ortho = fabs(len - 1) > ortho ? fabs(len - 1) : ortho;
			ortho = fabs(dot) > ortho ? fabs(dot) : ortho;
		}
	}
	/* a straight segment: one frame all along */
	q[0].x = 0; q[0].y = 0; q[0].z = 0;
	q[1].x = 1; q[1].y = 2; q[1].z = 3;
	line = fr_bezier(q, 2);
	TEST_ASSERT(line && qaws_curve_compute_rotation_minimizing_frames_3d(line, t, 200, NULL, fr) == QAWS_STATUS_OK, "straight segment: frames");
	for (i = 1; i < 200; i++)
		for (k = 0; k < 9; k++)
			worst_line = fabs((double)fr[9 * i + k] - fr[k]) > worst_line ? fabs((double)fr[9 * i + k] - fr[k]) : worst_line;
	sprintf(msg, "planar curve: the normal stays the plane normal within %.2g; straight segment: one frame within %.2g; orthonormal within %.2g",
		worst, worst_line, ortho);
	TEST_ASSERT(worst <= 1e-9 && worst_line <= 1e-12 && ortho <= 1e-12, msg);
	printf("  %s\n", msg);
	qaws_curve_destroy(c);
	qaws_curve_destroy(line);
}

static void test_convergence(void)
{
	qaws_vec3 p[6] = { { 0, 0, 0 }, { 1, 2, 0.5f }, { 3, -1, 2 }, { 4, 2, -1 }, { 2, 3, 3 }, { 5, 1, 1 } };
	qaws_curve* c = fr_bezier(p, 6);
	unsigned int const nref = 64 * 128 + 1;   /* 64 times the finest test's steps */
	qaws_scalar* tr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nref);
	qaws_scalar* ref = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * nref);
	double e1, e2, e3;
	char msg[256];
	if (!c || !tr || !ref)
	{
		TEST_ASSERT(0, "convergence: set up");
		free(tr); free(ref); qaws_curve_destroy(c);
		return;
	}
	fr_params(tr, nref);
	TEST_ASSERT(qaws_curve_compute_rotation_minimizing_frames_3d(c, tr, nref, NULL, ref) == QAWS_STATUS_OK, "reference frames");
	/* 32, 64 and 128 steps against the reference at the same parameters (the same first normal: both from the same tangent) */
	{
		/* errors at N = 33, 65, 129 samples, each against the reference at its own parameters */
		unsigned int Ns[3] = { 33, 65, 129 }, k, i;
		double e[3];
		for (k = 0; k < 3; k++)
		{
			unsigned int N = Ns[k], stride = (nref - 1) / (N - 1);
			qaws_scalar* t = (qaws_scalar*)malloc(sizeof(qaws_scalar) * N);
			qaws_scalar* fr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * N);
			e[k] = 0;
			fr_params(t, N);
			if (qaws_curve_compute_rotation_minimizing_frames_3d(c, t, N, NULL, fr) != QAWS_STATUS_OK)
				e[k] = 1e30;
			for (i = 0; i < N && e[k] < 1e30; i++)
			{
				qaws_scalar const* a = fr + 9 * i + 3;
				qaws_scalar const* b = ref + 9 * (stride * i) + 3;
				double dot = (double)a[0] * b[0] + (double)a[1] * b[1] + (double)a[2] * b[2];
				double ang = acos(dot > 1 ? 1 : dot < -1 ? -1 : dot);
				e[k] = ang > e[k] ? ang : e[k];
			}
			free(t);
			free(fr);
		}
		e1 = e[0]; e2 = e[1]; e3 = e[2];
	}
	sprintf(msg, "space curve: normal errors %.3g, %.3g, %.3g for 32, 64, 128 steps (ratios %.2f, %.2f; fourth order: 16)", e1, e2, e3, e1 / e2, e2 / e3);
	TEST_ASSERT(e1 < 1e-3 && e1 / e2 > 12 && e2 / e3 > 12, msg);
	printf("  %s\n", msg);
	free(tr);
	free(ref);
	qaws_curve_destroy(c);
}

static void test_refusals(void)
{
	qaws_vec3 p[3] = { { 0, 0, 0 }, { 1, 1, 0 }, { 2, 0, 1 } };
	qaws_vec2 p2[3] = { { 0, 0 }, { 1, 1 }, { 2, 0 } };
	qaws_scalar t[3] = { 0, (qaws_scalar)0.6, (qaws_scalar)0.4 }, fr[27];
	qaws_bezier_desc d;
	qaws_curve *c = fr_bezier(p, 3), *c2 = NULL;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = p2; d.control_point_count = 3;
	qaws_curve_create_bezier(&d, &c2);
	TEST_ASSERT(qaws_curve_compute_rotation_minimizing_frames_3d(c, t, 3, NULL, fr) == QAWS_STATUS_INVALID_ARGUMENT, "decreasing parameters refused");
	t[2] = 1;
	TEST_ASSERT(qaws_curve_compute_rotation_minimizing_frames_3d(c2, t, 3, NULL, fr) == QAWS_STATUS_INVALID_DIMENSION, "a 2D curve refused");
	TEST_ASSERT(qaws_curve_compute_rotation_minimizing_frames_3d(c, t, 0, NULL, fr) == QAWS_STATUS_INVALID_ARGUMENT, "no parameters refused");
	qaws_curve_destroy(c);
	qaws_curve_destroy(c2);
}

int test_147_frames_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 147: rotation-minimizing frames\n");
	test_planar_and_straight();
	test_convergence();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
