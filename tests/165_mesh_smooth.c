/*
 * Test 165: smoothing and fairing (qaws_mesh_3d_smooth, qaws_mesh_3d_fair)
 *
 *   - a unit icosphere with radial noise: each method lowers the noise;
 *     Laplacian steps shrink it, Taubin's keep the volume within 2%, and
 *     QAWS_SMOOTH_KEEP_VOLUME keeps it to rounding;
 *   - a noisy square grid: the border stays exactly where it was;
 *   - a roof (two planes meeting along a crease) with noise: with the crease
 *     as a feature line its vertices stay on the crease's plane and it stays
 *     sharp; without, it rounds off;
 *   - fairing a disturbed disk in a flat grid gives the plane back (orders
 *     1 to 3); fairing a disturbed cap of the sphere: the thin plate (order
 *     2) follows the sphere far closer than the membrane, order 3 closer still;
 *   - refusals: no mesh, an unknown method, order 4.
 */

#include "test_common.h"
#include "qaws_mesh_smooth.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static double st_rand(unsigned int* s)
{
	*s = *s * 1664525u + 1013904223u;
	return ((*s >> 8) & 0xFFFFFF) / (double)0x1000000 - 0.5;
}

static unsigned int st_mid(unsigned long long* key, unsigned int* val, unsigned int cap, double* x, unsigned int* nv, unsigned int a, unsigned int b)
{
	unsigned long long k = a < b ? ((unsigned long long)a << 32 | b) + 1 : ((unsigned long long)b << 32 | a) + 1;
	unsigned int i = (unsigned int)((k * 0x9E3779B97F4A7C15ull) >> 40) & (cap - 1), j;
	double l;
	while (key[i] && key[i] != k)
		i = (i + 1) & (cap - 1);
	if (key[i])
		return val[i];
	key[i] = k;
	val[i] = *nv;
	for (j = 0; j < 3; j++)
		x[3 * *nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
	l = sqrt(x[3 * *nv] * x[3 * *nv] + x[3 * *nv + 1] * x[3 * *nv + 1] + x[3 * *nv + 2] * x[3 * *nv + 2]);
	for (j = 0; j < 3; j++)
		x[3 * *nv + j] /= l;
	return (*nv)++;
}

/* a unit icosphere, its vertices moved radially by noise (amplitude) */
static qaws_mesh_3d* st_icosphere(unsigned int levels, double noise, unsigned int seed)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2;
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int nt = 20, nv = 12, cap, i, j, lv, *t, *t2, *val;
	unsigned long long* key;
	double* x;
	qaws_scalar* xs;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (lv = 0; lv < levels; lv++)
		nt *= 4;
	x = (double*)malloc(sizeof(double) * 3 * (nt / 2 + 2));
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	for (cap = 16; cap < 4 * nt; cap *= 2)
		;
	key = (unsigned long long*)malloc(sizeof(unsigned long long) * cap);
	val = (unsigned int*)malloc(sizeof(unsigned int) * cap);
	for (i = 0; i < 12; i++)
	{
		double l = sqrt(v0[3 * i] * v0[3 * i] + v0[3 * i + 1] * v0[3 * i + 1] + v0[3 * i + 2] * v0[3 * i + 2]);
		for (j = 0; j < 3; j++)
			x[3 * i + j] = v0[3 * i + j] / l;
	}
	for (i = 0; i < 60; i++)
		t[i] = (unsigned int)ic[i];
	for (nt = 20, lv = 0; lv < levels; lv++)
	{
		memset(key, 0, sizeof(unsigned long long) * cap);
		for (i = 0; i < nt; i++)
		{
			unsigned int a = t[3 * i], b = t[3 * i + 1], c = t[3 * i + 2];
			unsigned int ab = st_mid(key, val, cap, x, &nv, a, b), bc = st_mid(key, val, cap, x, &nv, b, c), ca = st_mid(key, val, cap, x, &nv, c, a);
			unsigned int const q[12] = { a, ab, ca, ab, b, bc, ca, bc, c, ab, bc, ca };
			memcpy(t2 + 12 * i, q, sizeof q);
		}
		nt *= 4;
		memcpy(t, t2, sizeof(unsigned int) * 3 * nt);
	}
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	for (i = 0; i < nv; i++)
	{
		double r = 1 + noise * st_rand(&seed);
		for (j = 0; j < 3; j++)
			xs[3 * i + j] = (qaws_scalar)(r * x[3 * i + j]);
	}
	memset(&d, 0, sizeof d);
	d.positions = xs;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	free(x); free(t); free(t2); free(key); free(val); free(xs);
	return m;
}

/* an n x n grid over [0, 1]^2, z from roof (|x - 1/2| times it) plus noise at the inner vertices */
static qaws_mesh_3d* st_grid(unsigned int n, double roof, double noise, unsigned int seed)
{
	unsigned int nv = (n + 1) * (n + 1), i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			unsigned int v = j * (n + 1) + i;
			double u = (double)i / n, w = (double)j / n, z = roof * fabs(u - 0.5);
			if (i && j && i < n && j < n)
				z += noise * st_rand(&seed);
			x[3 * v] = (qaws_scalar)u;
			x[3 * v + 1] = (qaws_scalar)w;
			x[3 * v + 2] = (qaws_scalar)z;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int v = j * (n + 1) + i;
			/* (diagonals turned to the crease at the middle, so it is a line of edges either way) */
			int flip = i < n / 2;
			t[k++] = v; t[k++] = v + 1; t[k++] = flip ? v + n + 1 : v + n + 2;
			t[k++] = flip ? v + 1 : v; t[k++] = v + n + 2; t[k++] = v + n + 1;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = 2 * n * n;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

/* the radii's spread about their mean (root mean square: the noise, whatever the shrinking), the volume */
static void st_sphere_stats(qaws_mesh_3d const* m, double* dev, double* vol)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), i;
	double s = 0, s2 = 0;
	*vol = 0;
	for (i = 0; i < nv; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		double r = sqrt((double)p.x * p.x + (double)p.y * p.y + (double)p.z * p.z);
		s += r;
		s2 += r * r;
	}
	*dev = sqrt(fmax(s2 / nv - (s / nv) * (s / nv), 0));
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		*vol += ((double)a.x * (b.y * e.z - b.z * e.y) + (double)a.y * (b.z * e.x - b.x * e.z) + (double)a.z * (b.x * e.y - b.y * e.x)) / 6;
	}
}

static void test_sphere(void)
{
	static char const* const name[4] = { "Laplacian, 10 steps", "Taubin, 10 steps", "implicit, one step", "bi-Laplacian, one step" };
	unsigned int k;
	double dev0, vol0;
	char msg[256];
	{
		qaws_mesh_3d* m = st_icosphere(4, 0.04, 7);
		st_sphere_stats(m, &dev0, &vol0);
		qaws_mesh_3d_destroy(m);
	}
	for (k = 0; k < 5; k++)
	{
		qaws_mesh_3d* m = st_icosphere(4, 0.04, 7);
		qaws_smooth_desc d;
		double dev, vol;
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.method = (qaws_smooth_method)(k < 4 ? k : 0);
		if (k == 4)
			d.flags = QAWS_SMOOTH_KEEP_VOLUME;
		st = qaws_mesh_3d_smooth(m, &d);
		st_sphere_stats(m, &dev, &vol);
		sprintf(msg, "noisy sphere, %s%s: radius spread %.5f (was %.5f), volume %.6f (was %.6f)", k < 4 ? name[k] : name[0],
			k == 4 ? ", keeping the volume" : "", dev, dev0, vol, vol0);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && dev < 0.5 * dev0 && (k != 1 || fabs(vol / vol0 - 1) < 0.02) && (k != 0 || vol < 0.99 * vol0)
			&& (k != 4 || fabs(vol / vol0 - 1) < 1e-12), msg);
		qaws_mesh_3d_destroy(m);
	}
}

static void test_border(void)
{
	unsigned int n = 20, nv = (n + 1) * (n + 1), i, k, same = 1;
	qaws_mesh_3d* m = st_grid(n, 0, 0.05, 3);
	qaws_scalar* x0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_scalar* x1 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double z0 = 0, z1 = 0;
	char msg[200];
	qaws_mesh_3d_get_positions(x0, m);
	for (k = 0; k < 2; k++)
	{
		qaws_smooth_desc d;
		memset(&d, 0, sizeof d);
		d.method = k ? QAWS_SMOOTH_IMPLICIT : QAWS_SMOOTH_TAUBIN;
		TEST_ASSERT(qaws_mesh_3d_smooth(m, &d) == QAWS_STATUS_OK, "noisy grid smoothed");
	}
	qaws_mesh_3d_get_positions(x1, m);
	for (i = 0; i < nv; i++)
	{
		if (qaws_mesh_3d_is_boundary_vertex(m, i))
			same &= memcmp(x0 + 3 * i, x1 + 3 * i, 3 * sizeof(qaws_scalar)) == 0;
		z0 += x0[3 * i + 2] * x0[3 * i + 2];
		z1 += x1[3 * i + 2] * x1[3 * i + 2];
	}
	sprintf(msg, "noisy grid, Taubin then implicit: the border unchanged %d, height rms %.5f (was %.5f)", same, sqrt(z1 / nv), sqrt(z0 / nv));
	printf("  %s\n", msg);
	TEST_ASSERT(same && z1 < 0.1 * z0, msg);
	qaws_mesh_3d_destroy(m);
	free(x0);
	free(x1);
}

static void test_crease(void)
{
	unsigned int n = 24, nv = (n + 1) * (n + 1), j, k;
	char msg[256];
	for (k = 0; k < 2; k++)
	{
		qaws_mesh_3d* m = st_grid(n, 0.5, 0.01, 11);
		unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * n);
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_smooth_desc d;
		double xdev = 0, zdrop = 0;
		for (j = 0; j < n; j++)
		{
			e[2 * j] = j * (n + 1) + n / 2;
			e[2 * j + 1] = (j + 1) * (n + 1) + n / 2;
		}
		if (k)
			qaws_mesh_3d_set_feature_edges(m, e, n);
		memset(&d, 0, sizeof d);
		d.method = QAWS_SMOOTH_TAUBIN;
		d.iterations = 30;
		qaws_mesh_3d_smooth(m, &d);
		qaws_mesh_3d_get_positions(x, m);
		/* the crease's vertices: off its plane x = 1/2, and below the ridge height 0 (the roof is 0 there) */
		for (j = 1; j < n; j++)
		{
			unsigned int v = j * (n + 1) + n / 2;
			xdev = fmax(xdev, fabs(x[3 * v] - 0.5));
			zdrop += x[3 * v + 2] / (n - 1);
		}
		sprintf(msg, "noisy roof, Taubin 30 steps, %s: crease off its plane by %.2e, its mean height %.4f (sharp: 0)",
			k ? "the crease a feature line" : "no features", xdev, zdrop);
		printf("  %s\n", msg);
		if (k)
			TEST_ASSERT(xdev < 1e-12 && fabs(zdrop) < 0.004, msg);
		else
			TEST_ASSERT(zdrop > 0.004, msg);
		qaws_mesh_3d_destroy(m);
		free(e);
		free(x);
	}
}

static void test_fair_plane(void)
{
	unsigned int n = 24, nv = (n + 1) * (n + 1), i, order;
	char msg[200];
	for (order = 1; order <= 3; order++)
	{
		qaws_mesh_3d* m = st_grid(n, 0, 0, 1);
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		unsigned char* fr = (unsigned char*)calloc(nv, 1);
		unsigned int seed = 5;
		qaws_fair_desc d;
		double zmax = 0;
		qaws_status st;
		qaws_mesh_3d_get_positions(x, m);
		for (i = 0; i < nv; i++)
		{
			double u = x[3 * i] - 0.5, w = x[3 * i + 1] - 0.5;
			if (u * u + w * w < 0.09)
			{
				fr[i] = 1;
				x[3 * i + 2] = (qaws_scalar)(0.2 * st_rand(&seed));
			}
		}
		qaws_mesh_3d_set_positions(m, x);
		memset(&d, 0, sizeof d);
		d.order = order;
		st = qaws_mesh_3d_fair(m, fr, &d);
		qaws_mesh_3d_get_positions(x, m);
		for (i = 0; i < nv; i++)
			zmax = fmax(zmax, fabs(x[3 * i + 2]));
		sprintf(msg, "a disturbed disk in a flat grid faired at order %u: height left %.2e", order, zmax);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && zmax < 1e-12, msg);
		qaws_mesh_3d_destroy(m);
		free(x);
		free(fr);
	}
}

static void test_fair_cap(void)
{
	unsigned int order, i;
	double dev[4] = { 0, 0, 0, 0 };
	char msg[200];
	for (order = 1; order <= 3; order++)
	{
		qaws_mesh_3d* m = st_icosphere(4, 0, 1);
		unsigned int nv = qaws_mesh_3d_get_vertex_count(m), seed = 9;
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		unsigned char* fr = (unsigned char*)calloc(nv, 1);
		qaws_fair_desc d;
		qaws_status st;
		qaws_mesh_3d_get_positions(x, m);
		for (i = 0; i < nv; i++)
			if (x[3 * i + 2] > 0.55)
			{
				fr[i] = 1;
				x[3 * i] += (qaws_scalar)(0.1 * st_rand(&seed));
				x[3 * i + 1] += (qaws_scalar)(0.1 * st_rand(&seed));
				x[3 * i + 2] += (qaws_scalar)(0.1 * st_rand(&seed));
			}
		qaws_mesh_3d_set_positions(m, x);
		memset(&d, 0, sizeof d);
		d.order = order;
		st = qaws_mesh_3d_fair(m, fr, &d);
		qaws_mesh_3d_get_positions(x, m);
		for (i = 0; i < nv; i++)
			if (fr[i])
				dev[order] = fmax(dev[order], fabs(sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]) - 1));
		TEST_ASSERT(st == QAWS_STATUS_OK, "cap faired");
		qaws_mesh_3d_destroy(m);
		free(x);
		free(fr);
	}
	sprintf(msg, "a disturbed cap of the sphere faired: off the sphere by %.4f (membrane), %.4f (thin plate), %.4f (order 3)", dev[1], dev[2],
		dev[3]);
	printf("  %s\n", msg);
	TEST_ASSERT(dev[2] < 0.25 * dev[1] && dev[3] < dev[2] && dev[2] < 0.1, msg);
}

static void test_refusals(void)
{
	qaws_mesh_3d* m = st_grid(2, 0, 0, 1);
	qaws_smooth_desc d;
	qaws_fair_desc f;
	unsigned char fr[9] = { 0 };
	memset(&d, 0, sizeof d);
	memset(&f, 0, sizeof f);
	TEST_ASSERT(qaws_mesh_3d_smooth(NULL, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "no mesh refused");
	d.method = (qaws_smooth_method)9;
	TEST_ASSERT(qaws_mesh_3d_smooth(m, &d) == QAWS_STATUS_INVALID_ARGUMENT, "an unknown method refused");
	f.order = 4;
	TEST_ASSERT(qaws_mesh_3d_fair(m, fr, &f) == QAWS_STATUS_INVALID_ARGUMENT, "order 4 refused");
	qaws_mesh_3d_destroy(m);
}

int test_165_mesh_smooth_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 165: smoothing and fairing\n");
	test_sphere();
	test_border();
	test_crease();
	test_fair_plane();
	test_fair_cap();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
