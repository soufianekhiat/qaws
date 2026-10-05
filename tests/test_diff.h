/*
 * test_diff.h - Shared helpers for differentiability tests
 *
 * Checks used across the diff suites:
 *   - adjoint identity   <y_bar, J x_dot> == <J^T y_bar, x_dot>   (exact, f32 and f64)
 *   - finite differences (f64 only for tight tolerances)
 *   - Taylor remainder   second-order expansion error shrinks as O(e^3)
 */

#ifndef QAWS_TEST_DIFF_H
#define QAWS_TEST_DIFF_H

#include "test_common.h"
#include "qaws_diff.h"
#include "core/qaws_dual_core.h"

#if QAWS_SCALAR_IS_FLOAT
#define DIFF_TOL 2e-3
#define DIFF_FD_STEP 1e-2
#else
#define DIFF_TOL 1e-9
#define DIFF_FD_STEP 1e-5
#endif

static unsigned int g_diff_seed = 12345u;

static void diff_seed(unsigned int s)
{
	g_diff_seed = s;
}

/* Uniform in [-1, 1]. */
static qaws_scalar diff_rand(void)
{
	g_diff_seed = g_diff_seed * 1664525u + 1013904223u;
	return (qaws_scalar)((double)(g_diff_seed >> 8) / (double)(1u << 24) * 2.0 - 1.0);
}

static void diff_rand_fill(qaws_scalar* x, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		x[i] = diff_rand();
}

static double diff_dot(qaws_scalar const* a, qaws_scalar const* b, unsigned int n)
{
	double s = 0.0;
	unsigned int i;
	for (i = 0; i < n; i++)
		s += (double)a[i] * (double)b[i];
	return s;
}

/* |a - b| <= tol * max(1, |a|, |b|) */
static int diff_close(double a, double b, double tol)
{
	double scale = 1.0, d = a - b;
	if (fabs(a) > scale) scale = fabs(a);
	if (fabs(b) > scale) scale = fabs(b);
	if (d < 0) d = -d;
	return d <= tol * scale;
}

static double diff_jet3_dot(qaws_curve_jet_3d const* a, qaws_curve_jet_3d const* b, unsigned int channels)
{
	double s = 0.0;
	unsigned int k;
	for (k = 0; k < 4; k++)
	{
		if (!(channels & (1u << k)))
			continue;
		s += (double)a->d[k].x * b->d[k].x + (double)a->d[k].y * b->d[k].y + (double)a->d[k].z * b->d[k].z;
	}
	return s;
}

static double diff_jet2_dot(qaws_curve_jet_2d const* a, qaws_curve_jet_2d const* b, unsigned int channels)
{
	double s = 0.0;
	unsigned int k;
	for (k = 0; k < 4; k++)
	{
		if (!(channels & (1u << k)))
			continue;
		s += (double)a->d[k].x * b->d[k].x + (double)a->d[k].y * b->d[k].y;
	}
	return s;
}

static void diff_rand_jet3(qaws_curve_jet_3d* j)
{
	unsigned int k;
	for (k = 0; k < 4; k++)
	{
		j->d[k].x = diff_rand();
		j->d[k].y = diff_rand();
		j->d[k].z = diff_rand();
	}
	j->channels = 0xF;
}

static void diff_rand_jet2(qaws_curve_jet_2d* j)
{
	unsigned int k;
	for (k = 0; k < 4; k++)
	{
		j->d[k].x = diff_rand();
		j->d[k].y = diff_rand();
	}
	j->channels = 0xF;
}

static qaws_vec3 diff_rand_vec3(void)
{
	qaws_vec3 v;
	v.x = diff_rand();
	v.y = diff_rand();
	v.z = diff_rand();
	return v;
}

#endif /* QAWS_TEST_DIFF_H */
