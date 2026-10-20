/*
 * Test 116: hierarchical metrics and the iterative solve (qaws_repulsive.h,
 * R4)
 *
 *   - the metric applied from Barnes-Hut lists: theta 0 is the dense
 *     metric's product (curves in 2D and 3D, several components; a mesh);
 *     the product converges to it as theta falls
 *   - the iterative direction (pins eliminated, the other rows by the
 *     Schur complement, GMRES preconditioned by the Laplacian sandwich)
 *     with every pair near is the dense direction; at theta 0.25 close to
 *     it; curves (pins, lengths, a held direction, 2D) and meshes (area,
 *     volume, pins)
 *   - hierarchical solves lower the energy with their constraints held;
 *     the time per step against the dense mode reported
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "internal/qaws_internal_repulsive.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define RH_PI 3.14159265358979323846

static void rh_sphere(unsigned int level, double* x, unsigned int* tris, unsigned int* nv, unsigned int* nf)
{
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[8][3] = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 }, { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
	unsigned int l, i, n = 6, f = 8, q;
	unsigned int* t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 4 * 8 * (1u << (2 * level)));
	memcpy(x, ov, sizeof(ov));
	memcpy(tris, of, sizeof(of));
	for (l = 0; l < level; l++)
	{
		unsigned int nf2 = 0;
		for (i = 0; i < f; i++)
		{
			unsigned int c[3], mid[3], k;
			for (k = 0; k < 3; k++) c[k] = tris[3 * i + k];
			for (k = 0; k < 3; k++)
			{
				unsigned int a = c[k], b = c[(k + 1) % 3], j;
				double m[3], len;
				for (q = 0; q < 3; q++) m[q] = 0.5 * (x[3 * a + q] + x[3 * b + q]);
				len = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
				for (q = 0; q < 3; q++) m[q] /= len;
				for (j = 0; j < n && !(fabs(x[3 * j] - m[0]) < 1e-12 && fabs(x[3 * j + 1] - m[1]) < 1e-12 && fabs(x[3 * j + 2] - m[2]) < 1e-12); j++)
					;
				if (j == n) { for (q = 0; q < 3; q++) x[3 * n + q] = m[q]; n++; }
				mid[k] = j;
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
		for (q = 0; q < 3; q++) x[3 * i + q] *= r;
	}
	free(t2);
	*nv = n; *nf = f;
}

/* a lumpy trefoil of n edges (3D), or two lumpy circles of n / 2 edges (2D) */
static void rh_curves(unsigned int dim, unsigned int n, double* x, unsigned int* e)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * RH_PI * i / n;
		if (dim == 3)
		{
			x[3 * i] = sin(t) + 2 * sin(2 * t) + 0.13 * cos(5 * t); x[3 * i + 1] = cos(t) - 2 * cos(2 * t) + 0.07 * sin(4 * t); x[3 * i + 2] = -0.6 * sin(3 * t);
			e[2 * i] = i; e[2 * i + 1] = (i + 1) % n;
		}
		else
		{
			unsigned int h = n / 2, k = i / h, j = i % h;
			double s = 2 * RH_PI * j / h;
			x[2 * i] = (1 + 0.1 * cos(3 * s)) * cos(s) + 2.6 * k; x[2 * i + 1] = (1 + 0.1 * cos(3 * s)) * sin(s) + 0.3 * k;
			e[2 * i] = i; e[2 * i + 1] = k * h + (j + 1) % h;
		}
	}
}

static double rh_relerr(double const* a, double const* b, unsigned int n)
{
	double e = 0, m = 0;
	unsigned int i;
	for (i = 0; i < n; i++) { e = fmax(e, fabs(a[i] - b[i])); m = fmax(m, fabs(b[i])); }
	return e / m;
}

/* theta 0 against the dense product, then theta 0.5 / 0.25 / 0.1 */
static void rh_check(qaws_internal_rbh_op* (*make)(void const*, double const*, double), void const* obj, double const* x, double const* A,
	unsigned int nv, char const* name)
{
	double *u = (double*)malloc(sizeof(double) * nv), *y0 = (double*)malloc(sizeof(double) * nv), *y = (double*)malloc(sizeof(double) * nv), e[4];
	static double const th[4] = { 0, 0.5, 0.25, 0.1 };
	unsigned int i, j, t;
	size_t nfar[4];
	char msg[300];
	for (i = 0; i < nv; i++) u[i] = sin(0.37 * i) + 0.5 * cos(1.3 * i);
	for (i = 0; i < nv; i++)
	{
		y0[i] = 0;
		for (j = 0; j < nv; j++) y0[i] += A[(size_t)i * nv + j] * u[j];
	}
	for (t = 0; t < 4; t++)
	{
		qaws_internal_rbh_op* op = make(obj, x, th[t]);
		e[t] = 1;
		nfar[t] = 0;
		if (!op)
			continue;
		qaws_internal_rbh_apply(y, op, u);
		e[t] = rh_relerr(y, y0, nv);
		qaws_internal_rbh_links(&nfar[t], op);
		qaws_internal_rbh_destroy(op);
	}
	sprintf(msg, "%s: theta 0 against the dense metric %.2g; theta 0.5 / 0.25 / 0.1: %.2g / %.2g / %.2g (far links %u / %u / %u)", name, e[0], e[1], e[2],
		e[3], (unsigned int)nfar[1], (unsigned int)nfar[2], (unsigned int)nfar[3]);
	TEST_ASSERT(e[0] <= 1e-12 && e[3] <= e[2] && e[2] <= e[1] && e[2] <= 5e-2 && nfar[1] > 0, msg);
	free(u); free(y0); free(y);
}

static qaws_internal_rbh_op* rh_make_curves(void const* c, double const* x, double theta)
{
	qaws_internal_rbh_op* op = NULL;
	qaws_internal_rbh_curves_metric(&op, (qaws_repulsive_curves const*)c, x, theta);
	return op;
}

static qaws_internal_rbh_op* rh_make_mesh(void const* m, double const* x, double theta)
{
	qaws_internal_rbh_op* op = NULL;
	qaws_internal_rbh_mesh_metric(&op, (qaws_repulsive_mesh const*)m, x, theta);
	return op;
}

static void test_metric(void)
{
	unsigned int dim;
	for (dim = 2; dim <= 3; dim++)
	{
		enum { N = 200 };
		double x[3 * N], *A = (double*)malloc(sizeof(double) * N * N);
		unsigned int e[2 * N];
		qaws_repulsive_curves_desc d;
		qaws_repulsive_curves* c = NULL;
		rh_curves(dim, N, x, e);
		memset(&d, 0, sizeof(d));
		d.dimension = dim; d.vertex_count = N; d.edges = e; d.edge_count = N;
		qaws_repulsive_curves_create(&c, &d);
		qaws_internal_repulsive_curves_metric(A, c, x);
		rh_check(rh_make_curves, c, x, A, N, dim == 3 ? "trefoil (200 edges)" : "two circles in 2D (200 edges)");
		qaws_repulsive_curves_destroy(c);
		free(A);
	}
	{
		double* x = (double*)malloc(sizeof(double) * 3 * 300);
		unsigned int* tris = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 520);
		unsigned int nv, nf;
		double* A;
		qaws_repulsive_mesh_desc d;
		qaws_repulsive_mesh* m = NULL;
		rh_sphere(3, x, tris, &nv, &nf);
		A = (double*)malloc(sizeof(double) * nv * nv);
		memset(&d, 0, sizeof(d));
		d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf;
		qaws_repulsive_mesh_create(&m, &d);
		qaws_internal_repulsive_mesh_metric(A, m, x);
		rh_check(rh_make_mesh, m, x, A, nv, "lumpy sphere (512 faces)");
		qaws_repulsive_mesh_destroy(m);
		free(x); free(tris); free(A);
	}
}

static double rh_dot0(qaws_scalar const* a, unsigned int n)
{
	double s = 0;
	unsigned int i;
	for (i = 0; i < n; i++) s += (double)a[i] * a[i];
	return s;
}

/* the iterative direction with every pair near (theta 1e-9) against the
   dense one; then theta 0.25 near it */
static void rh_dirs(void* solver_dense, void* solver_h, void* solver_q, int mesh, qaws_scalar const* x, unsigned int N, char const* name)
{
	qaws_scalar *g0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * N), *g1 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * N),
		*g2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * N);
	qaws_status s0, s1, s2;
	double e1 = 1, e2 = 1, m = 0, a = 0, b = 0;
	unsigned int i;
	char msg[300];
	if (mesh)
	{
		s0 = qaws_repulsive_mesh_solver_direction(g0, (qaws_repulsive_mesh_solver*)solver_dense, x);
		s1 = qaws_repulsive_mesh_solver_direction(g1, (qaws_repulsive_mesh_solver*)solver_h, x);
		s2 = qaws_repulsive_mesh_solver_direction(g2, (qaws_repulsive_mesh_solver*)solver_q, x);
	}
	else
	{
		s0 = qaws_repulsive_curves_solver_direction(g0, (qaws_repulsive_curves_solver*)solver_dense, x);
		s1 = qaws_repulsive_curves_solver_direction(g1, (qaws_repulsive_curves_solver*)solver_h, x);
		s2 = qaws_repulsive_curves_solver_direction(g2, (qaws_repulsive_curves_solver*)solver_q, x);
	}
	if (s0 == QAWS_STATUS_OK && s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK)
	{
		e1 = e2 = 0;
		for (i = 0; i < N; i++)
		{
			m = fmax(m, fabs((double)g0[i]));
			e1 = fmax(e1, fabs((double)g1[i] - g0[i]));
			a += (double)g2[i] * g0[i]; b += (double)g2[i] * g2[i];
		}
		e1 /= m;
		e2 = 1 - a / sqrt(b * rh_dot0(g0, N));
	}
	sprintf(msg, "%s: the iterative direction (all pairs near) against the dense, %.2g; at theta 0.25, 1 - cos %.2g", name, e1, e2);
	TEST_ASSERT(e1 <= (sizeof(qaws_scalar) == 4 ? 1e-4 : 1e-6) && e2 <= 1e-2, msg);
	free(g0); free(g1); free(g2);
}

static qaws_repulsive_curves* rh_curves_make(unsigned int dim, unsigned int n, unsigned int* e)
{
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = dim; d.vertex_count = n; d.edges = e; d.edge_count = n;
	qaws_repulsive_curves_create(&c, &d);
	return c;
}

static void rh_curve_dirs(unsigned int dim, qaws_repulsive_curves_solve_desc* sd, char const* name)
{
	enum { N = 80 };
	double xd[3 * N];
	qaws_scalar x[3 * N];
	unsigned int e[2 * N], i, t;
	qaws_repulsive_curves* c;
	qaws_repulsive_curves_solver* s[3] = { NULL, NULL, NULL };
	static double const th[3] = { 0, 1e-9, 0.25 };
	rh_curves(dim, N, xd, e);
	for (i = 0; i < dim * N; i++) x[i] = (qaws_scalar)xd[i];
	c = rh_curves_make(dim, N, e);
	for (t = 0; t < 3; t++)
	{
		sd->theta = (qaws_scalar)th[t];
		qaws_repulsive_curves_solver_create(&s[t], c, x, sd);
	}
	if (s[0] && s[1] && s[2])
		rh_dirs(s[0], s[1], s[2], 0, x, dim * N, name);
	else
		TEST_ASSERT(0, name);
	for (t = 0; t < 3; t++) qaws_repulsive_curves_solver_destroy(s[t]);
	qaws_repulsive_curves_destroy(c);
}

static void rh_mesh_dirs(qaws_repulsive_mesh_solve_desc* sd, char const* name)
{
	double* xd = (double*)malloc(sizeof(double) * 3 * 70);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 70);
	unsigned int* tris = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 130);
	unsigned int nv, nf, i, t;
	qaws_repulsive_mesh_desc d;
	qaws_repulsive_mesh* m = NULL;
	qaws_repulsive_mesh_solver* s[3] = { NULL, NULL, NULL };
	static double const th[3] = { 0, 1e-9, 0.25 };
	rh_sphere(2, xd, tris, &nv, &nf);
	for (i = 0; i < 3 * nv; i++) x[i] = (qaws_scalar)xd[i];
	memset(&d, 0, sizeof(d));
	d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf;
	qaws_repulsive_mesh_create(&m, &d);
	for (t = 0; t < 3; t++)
	{
		sd->theta = (qaws_scalar)th[t];
		qaws_repulsive_mesh_solver_create(&s[t], m, x, sd);
	}
	if (s[0] && s[1] && s[2])
		rh_dirs(s[0], s[1], s[2], 1, x, 3 * nv, name);
	else
		TEST_ASSERT(0, name);
	for (t = 0; t < 3; t++) qaws_repulsive_mesh_solver_destroy(s[t]);
	qaws_repulsive_mesh_destroy(m);
	free(xd); free(x); free(tris);
}

static void test_directions(void)
{
	static unsigned int const pins[2] = { 0, 40 }, tedge[1] = { 3 }, mpins[2] = { 0, 1 };
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_mesh_solve_desc md;
	memset(&sd, 0, sizeof(sd));
	rh_curve_dirs(3, &sd, "trefoil (80 edges), barycenter held");
	memset(&sd, 0, sizeof(sd));
	sd.pinned = pins; sd.pinned_count = 2; sd.hold_total_length = 1;
	rh_curve_dirs(3, &sd, "trefoil, two pins and the total length");
	memset(&sd, 0, sizeof(sd));
	sd.tangent_edges = tedge; sd.tangent_edge_count = 1; sd.hold_total_length = 1;
	rh_curve_dirs(2, &sd, "two circles in 2D, barycenter, a held direction and the total length");
	memset(&md, 0, sizeof(md));
	md.hold_area = 1; md.hold_volume = 1;
	rh_mesh_dirs(&md, "lumpy sphere (128 faces), barycenter, area and volume");
	memset(&md, 0, sizeof(md));
	md.pinned = mpins; md.pinned_count = 2; md.hold_area = 1;
	rh_mesh_dirs(&md, "lumpy sphere, two pins and the area");
}

/* hierarchical solves: the energy falls, the constraints held; the time
   per step against the dense mode */
static void test_solves(void)
{
	/* the solves at a third of the size by default (the timings they report are the heavy tier's point) */
	{
		unsigned int const N = test_heavy() ? 600 : 200;
		double* xd = (double*)malloc(sizeof(double) * 3 * N);
		qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * N), *x2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * N);
		unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * N);
		unsigned int i, t;
		qaws_repulsive_curves* c;
		qaws_repulsive_curves_solve_desc sd;
		qaws_repulsive_report r[2];
		double tm[2];
		char msg[300];
		rh_curves(3, N, xd, e);
		for (i = 0; i < 3 * N; i++) x[i] = x2[i] = (qaws_scalar)xd[i];
		c = rh_curves_make(3, N, e);
		for (t = 0; t < 2; t++)
		{
			clock_t c0 = clock();
			memset(&sd, 0, sizeof(sd));
			sd.hold_total_length = 1; sd.max_iterations = t ? 20 : 3; sd.theta = (qaws_scalar)(t ? 0.25 : 0);
			memset(&r[t], 0, sizeof(r[t]));
			qaws_repulsive_curves_solve(&r[t], t ? x : x2, c, &sd);
			tm[t] = (double)(clock() - c0) / CLOCKS_PER_SEC / (r[t].iterations ? r[t].iterations : 1);
		}
		sprintf(msg, "trefoil (%u edges), total length held: hierarchical (theta 0.25) E %.4g -> %.4g in %u steps, length off %.2g, reason %d; "
			"%.3f s per step against %.3f s dense", N,
			(double)r[1].initial_energy, (double)r[1].energy, r[1].iterations, (double)r[1].constraint_residual, (int)r[1].reason, tm[1], tm[0]);
		TEST_ASSERT(r[1].energy < 0.9 * r[1].initial_energy && r[1].iterations >= 10 && (double)r[1].constraint_residual <= (sizeof(qaws_scalar) == 4 ? 1e-4 : 1e-7), msg);
		qaws_repulsive_curves_destroy(c);
		free(xd); free(x); free(x2); free(e);
	}
	{
		double* xd = (double*)malloc(sizeof(double) * 3 * 1100);
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 1100);
		unsigned int* tris = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 2100);
		unsigned int nv, nf, i;
		qaws_repulsive_mesh_desc d;
		qaws_repulsive_mesh* m = NULL;
		qaws_repulsive_mesh_solve_desc md;
		qaws_repulsive_report r;
		clock_t c0;
		double tm;
		char msg[300];
		rh_sphere(test_heavy() ? 4 : 3, xd, tris, &nv, &nf);
		for (i = 0; i < 3 * nv; i++) x[i] = (qaws_scalar)xd[i];
		memset(&d, 0, sizeof(d));
		d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf;
		qaws_repulsive_mesh_create(&m, &d);
		memset(&md, 0, sizeof(md));
		md.hold_area = 1; md.hold_volume = 1; md.max_iterations = 8; md.theta = (qaws_scalar)0.25;
		memset(&r, 0, sizeof(r));
		c0 = clock();
		qaws_repulsive_mesh_solve(&r, x, m, &md);
		tm = (double)(clock() - c0) / CLOCKS_PER_SEC / (r.iterations ? r.iterations : 1);
		sprintf(msg, "lumpy sphere (%u faces), area and volume held: hierarchical (theta 0.25) E %.4g -> %.4g in %u steps, residual %.2g, reason %d; %.2f s per step",
			nf, (double)r.initial_energy, (double)r.energy, r.iterations, (double)r.constraint_residual, (int)r.reason, tm);
		TEST_ASSERT(r.energy < r.initial_energy && r.iterations >= 4 && (double)r.constraint_residual <= (sizeof(qaws_scalar) == 4 ? 1e-4 : 1e-7), msg);
		qaws_repulsive_mesh_destroy(m);
		free(xd); free(x); free(tris);
	}
}

int test_116_repulsive_hsolve_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 116: hierarchical metrics and the iterative solve\n");
	test_metric();
	test_directions();
	test_solves();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
