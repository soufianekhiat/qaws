/*
 * Test 115: Barnes-Hut energies and gradients (qaws_repulsive.h, R4)
 *
 *   - theta = 0 is the direct evaluator, energy and gradient
 *   - the approximation converges to the direct values as theta falls
 *   - the gradient (the paper's scheme: each element's own derivative of
 *     both orders, far nodes through their averaged projectors) converges
 *     to the direct gradient as theta falls
 *   - curves (a 3D trefoil polygon) and meshes (a refined lumpy sphere);
 *     the time against the direct evaluator at the larger sizes
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define RB_F32 (sizeof(qaws_scalar) == 4)
#define RB_PI 3.14159265358979323846

typedef struct rb_case
{
	qaws_repulsive_curves* c;
	qaws_repulsive_mesh* m;
	unsigned int n;                 /* coordinates */
} rb_case;

static qaws_status rb_direct(rb_case const* k, qaws_scalar const* x, double* E, qaws_scalar* g)
{
	qaws_scalar e = 0;
	qaws_status st;
	memset(g, 0, sizeof(qaws_scalar) * k->n);
	st = k->c ? qaws_repulsive_curves_gradient(g, &e, k->c, x, 1) : qaws_repulsive_mesh_gradient(g, &e, k->m, x, 1);
	*E = (double)e;
	return st;
}

static qaws_status rb_bh(rb_case const* k, qaws_scalar const* x, double theta, double* E, qaws_scalar* g, qaws_repulsive_bh_stats* stats)
{
	qaws_scalar e = 0;
	qaws_status st;
	if (g) memset(g, 0, sizeof(qaws_scalar) * k->n);
	st = k->c ? qaws_repulsive_curves_bh(g, &e, stats, k->c, x, (qaws_scalar)theta, 1) : qaws_repulsive_mesh_bh(g, &e, stats, k->m, x, (qaws_scalar)theta, 1);
	*E = (double)e;
	return st;
}

static double rb_gerr(qaws_scalar const* a, qaws_scalar const* b, unsigned int n)
{
	double e = 0, m = 0;
	unsigned int i;
	for (i = 0; i < n; i++) { e = fmax(e, fabs((double)a[i] - b[i])); m = fmax(m, fabs((double)b[i])); }
	return e / m;
}

static void rb_check(rb_case const* k, qaws_scalar* x, char const* name)
{
	qaws_scalar *g0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * k->n), *g1 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * k->n);
	qaws_scalar* y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * k->n);
	double E0 = 0, E1 = 0, eE[3], eG[3], th[3] = { 0.5, 0.25, 0.1 }, fd = 0, an = 0;
	qaws_repulsive_bh_stats st0, st1[3];
	qaws_status s0, s1;
	unsigned int i, t;
	char msg[400];
	s0 = rb_direct(k, x, &E0, g0);
	s1 = rb_bh(k, x, 0, &E1, g1, &st0);
	sprintf(msg, "%s, theta 0: the direct energy %.2g and gradient %.2g off (%llu exact pairs, %llu far)",
		name, fabs(E1 - E0) / E0, rb_gerr(g1, g0, k->n), st0.near_pairs, st0.far_pairs);
	TEST_ASSERT(s0 == QAWS_STATUS_OK && s1 == QAWS_STATUS_OK && fabs(E1 - E0) / E0 <= (RB_F32 ? 1e-6 : 1e-12) &&
		rb_gerr(g1, g0, k->n) <= (RB_F32 ? 1e-5 : 1e-12) && st0.far_pairs == 0, msg);
	for (t = 0; t < 3; t++)
	{
		rb_bh(k, x, th[t], &E1, g1, &st1[t]);
		eE[t] = fabs(E1 - E0) / E0;
		eG[t] = rb_gerr(g1, g0, k->n);
	}
	sprintf(msg, "%s, theta 0.5 / 0.25 / 0.1: energy %.2g / %.2g / %.2g off, gradient %.2g / %.2g / %.2g off (max norm); far pairs %llu / %llu / %llu of %llu",
		name, eE[0], eE[1], eE[2], eG[0], eG[1], eG[2], st1[0].far_pairs, st1[1].far_pairs, st1[2].far_pairs, st0.near_pairs);
	TEST_ASSERT(eE[2] <= eE[1] && eE[1] <= eE[0] && eG[2] <= eG[1] && eG[1] <= eG[0] && eE[1] <= 1e-2 && eG[2] <= 1e-2 && st1[0].far_pairs > 0, msg);
	(void)y; (void)an; (void)fd; (void)i;
	free(g0); free(g1); free(y);
}

static void rb_trefoil(qaws_scalar* x, unsigned int* e, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * RB_PI * i / n;
		/* (a lumpy trefoil: a symmetric one ties edge centres on split axes, and
		   the median split then builds a different tree at x + h v and x - h v) */
		x[3 * i] = (qaws_scalar)(sin(t) + 2 * sin(2 * t) + 0.13 * cos(5 * t)); x[3 * i + 1] = (qaws_scalar)(cos(t) - 2 * cos(2 * t) + 0.07 * sin(4 * t)); x[3 * i + 2] = (qaws_scalar)(-0.6 * sin(3 * t));
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % n;
	}
}

static void test_curves(void)
{
	enum { N = 300 };
	qaws_scalar x[3 * N];
	unsigned int e[2 * N];
	qaws_repulsive_curves_desc d;
	rb_case k;
	memset(&k, 0, sizeof(k));
	rb_trefoil(x, e, N);
	memset(&d, 0, sizeof(d));
	d.dimension = 3; d.vertex_count = N; d.edges = e; d.edge_count = N;
	qaws_repulsive_curves_create(&k.c, &d);
	k.n = 3 * N;
	rb_check(&k, x, "trefoil (300 edges)");
	qaws_repulsive_curves_destroy(k.c);
}

/* an octahedron subdivided `level` times, lumpy */
static void rb_sphere(unsigned int level, qaws_scalar* xs, unsigned int* tris, unsigned int* nv, unsigned int* nf)
{
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[8][3] = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 }, { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
	unsigned int l, i, n = 6, f = 8, q;
	double* x = (double*)malloc(sizeof(double) * 3 * (2 + 4 * (1u << (2 * level))) * 2);
	unsigned int* t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 8 * ((size_t)1 << (2 * level)));
	memcpy(x, ov, sizeof(ov));
	memcpy(tris, of, sizeof(of));
	for (l = 0; l < level; l++)
	{
		unsigned int nf2 = 0;
		for (i = 0; i < f; i++)
		{
			unsigned int c[3], mid[3], kk;
			for (kk = 0; kk < 3; kk++) c[kk] = tris[3 * i + kk];
			for (kk = 0; kk < 3; kk++)
			{
				unsigned int a = c[kk], b = c[(kk + 1) % 3], j;
				double m[3], len;
				for (q = 0; q < 3; q++) m[q] = 0.5 * (x[3 * a + q] + x[3 * b + q]);
				len = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
				for (q = 0; q < 3; q++) m[q] /= len;
				for (j = 0; j < n && !(fabs(x[3 * j] - m[0]) < 1e-12 && fabs(x[3 * j + 1] - m[1]) < 1e-12 && fabs(x[3 * j + 2] - m[2]) < 1e-12); j++)
					;
				if (j == n) { for (q = 0; q < 3; q++) x[3 * n + q] = m[q]; n++; }
				mid[kk] = j;
			}
			t2[3 * nf2] = c[0]; t2[3 * nf2 + 1] = mid[0]; t2[3 * nf2 + 2] = mid[2]; nf2++;
			t2[3 * nf2] = mid[0]; t2[3 * nf2 + 1] = c[1]; t2[3 * nf2 + 2] = mid[1]; nf2++;
			t2[3 * nf2] = mid[2]; t2[3 * nf2 + 1] = mid[1]; t2[3 * nf2 + 2] = c[2]; nf2++;
			t2[3 * nf2] = mid[0]; t2[3 * nf2 + 1] = mid[1]; t2[3 * nf2 + 2] = mid[2]; nf2++;
		}
		memcpy(tris, t2, sizeof(unsigned int) * 3 * nf2);
		f = nf2;
	}
	for (i = 0; i < n; i++)
	{
		double r = 1 + 0.15 * sin(3 * x[3 * i]) * cos(2 * x[3 * i + 1]);
		for (q = 0; q < 3; q++) xs[3 * i + q] = (qaws_scalar)(x[3 * i + q] * r);
	}
	free(x); free(t2);
	*nv = n; *nf = f;
}

static void test_meshes(void)
{
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 1100);
	unsigned int* tris = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 2100);
	unsigned int nv, nf;
	qaws_repulsive_mesh_desc d;
	rb_case k;
	memset(&k, 0, sizeof(k));
	rb_sphere(3, x, tris, &nv, &nf);
	memset(&d, 0, sizeof(d));
	d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf;
	qaws_repulsive_mesh_create(&k.m, &d);
	k.n = 3 * nv;
	rb_check(&k, x, "lumpy sphere (512 faces)");
	qaws_repulsive_mesh_destroy(k.m);
	free(x); free(tris);
}

/* time: direct against Barnes-Hut (theta 0.25), energy and gradient */
static void test_timing(void)
{
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 2100);
	qaws_scalar* g = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 2100);
	unsigned int* tris = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 8200);
	unsigned int nv, nf;
	qaws_repulsive_mesh_desc d;
	rb_case k;
	double E0, E1, td, tb;
	clock_t c0;
	char msg[300];
	memset(&k, 0, sizeof(k));
	rb_sphere(4, x, tris, &nv, &nf);
	memset(&d, 0, sizeof(d));
	d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf;
	qaws_repulsive_mesh_create(&k.m, &d);
	k.n = 3 * nv;
	c0 = clock(); rb_direct(&k, x, &E0, g); td = (double)(clock() - c0) / CLOCKS_PER_SEC;
	c0 = clock(); rb_bh(&k, x, 0.25, &E1, g, NULL); tb = (double)(clock() - c0) / CLOCKS_PER_SEC;
	sprintf(msg, "time, a sphere of %u faces, energy and gradient: direct %.2f s, Barnes-Hut (theta 0.25) %.2f s, energy %.2g off", nf, td, tb, fabs(E1 - E0) / E0);
	/* (the times reported, not asserted: they vary with the machine's load) */
	TEST_ASSERT(fabs(E1 - E0) / E0 <= 2e-2, msg);
	qaws_repulsive_mesh_destroy(k.m);
	free(x); free(g); free(tris);
}

/* an executor that runs the chunks backwards: the results must not move */
static void rb_backwards(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int i;
	(void)user;
	for (i = count; i-- > 0;) task(ctx, i, i + 1);
}

static void test_executor(void)
{
	enum { N = 700 };
	qaws_scalar x[3 * N], g0[3 * N], g1[3 * N], E0 = 0, E1 = 0;
	unsigned int e[2 * N], i, same = 1;
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves *c0 = NULL, *c1 = NULL;
	qaws_batch_executor ex;
	char msg[200];
	rb_trefoil(x, e, N);
	memset(&d, 0, sizeof(d));
	d.dimension = 3; d.vertex_count = N; d.edges = e; d.edge_count = N;
	qaws_repulsive_curves_create(&c0, &d);
	ex.parallel_for = rb_backwards; ex.user = NULL;
	d.executor = &ex;
	qaws_repulsive_curves_create(&c1, &d);
	memset(g0, 0, sizeof(g0)); memset(g1, 0, sizeof(g1));
	qaws_repulsive_curves_bh(g0, &E0, NULL, c0, x, (qaws_scalar)0.25, 1);
	qaws_repulsive_curves_bh(g1, &E1, NULL, c1, x, (qaws_scalar)0.25, 1);
	for (i = 0; i < 3 * N; i++) same &= g0[i] == g1[i];
	sprintf(msg, "an executor running the chunks backwards: energy and gradient identical %d (E %.17g, %.17g)", same && E0 == E1, (double)E0, (double)E1);
	TEST_ASSERT(same && E0 == E1, msg);
	qaws_repulsive_curves_destroy(c0); qaws_repulsive_curves_destroy(c1);
}

int test_115_repulsive_bh_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 115: Barnes-Hut tangent-point energies and gradients\n");
	test_curves();
	test_meshes();
	test_timing();
	test_executor();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
