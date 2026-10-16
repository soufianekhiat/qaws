/*
 * Test 119: the backend-neutral tangent-point kernels
 *
 * core/qaws_repulsive_core.h on the C backend, driven as the GPU kernels
 * drive it (a pass per element over every other element, then a gather per
 * vertex over the elements holding it, no atomics), must reproduce the
 * library: the energy and gradient of qaws_repulsive_*_gradient, and the
 * metric's product (the dense metric of the solvers), for a closed 3D
 * curve with a second component, an open 2D curve (z = 0), and a sphere.
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "internal/qaws_internal_repulsive.h"
#include "core/qaws_repulsive_core.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RK_F32 (sizeof(qaws_scalar) == 4)
#define RK_PI 3.14159265358979323846

static double rk_rel(double const* a, double const* b, unsigned int n)
{
	double e = 0, m = 0;
	unsigned int i;
	for (i = 0; i < n; i++) { e = fmax(e, fabs(a[i] - b[i])); m = fmax(m, fabs(b[i])); }
	return e / fmax(m, 1e-300);
}

/* the vertex -> (element, slot) adjacency, CSR */
static void rk_adjacency(unsigned int const* idx, unsigned int ne, unsigned int k, unsigned int nv, unsigned int* off, unsigned int* ent)
{
	unsigned int i, a, *fill = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	memset(off, 0, sizeof(unsigned int) * (nv + 1));
	for (i = 0; i < k * ne; i++) off[idx[i] + 1]++;
	for (i = 0; i < nv; i++) off[i + 1] += off[i];
	for (i = 0; i < ne; i++)
		for (a = 0; a < k; a++) ent[off[idx[k * i + a]] + fill[idx[k * i + a]]++] = k * i + a;
	free(fill);
}

static void rk_load(qaws_scalar const* x, unsigned int dim, unsigned int const* idx, unsigned int k, unsigned int I, qaws_scalar* p)
{
	unsigned int a, c;
	for (a = 0; a < k; a++)
		for (c = 0; c < 3; c++) p[3 * a + c] = c < dim ? x[dim * idx[k * I + a] + c] : 0;
}

/* the kernels' passes: energy, gradient (dim per vertex), y = A u */
static void rk_kernels(qaws_scalar const* x, unsigned int dim, unsigned int nv, unsigned int const* idx, unsigned int ne, unsigned int k, double alpha, double beta,
	qaws_scalar const* u, double* E, double* grad, double* y)
{
	unsigned int* off = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
	unsigned int* ent = (unsigned int*)malloc(sizeof(unsigned int) * k * ne);
	qaws_scalar* blocks = (qaws_scalar*)calloc(9 * ne, sizeof(qaws_scalar));
	qaws_scalar* acc = (qaws_scalar*)calloc(4 * ne, sizeof(qaws_scalar));
	qaws_scalar* coef = (qaws_scalar*)calloc(9 * ne, sizeof(qaws_scalar));
	double* eb = (double*)calloc(ne, sizeof(double));
	double ex = k == 2 ? 2 * ((beta - 1) / alpha - 1) + 1 : 2 * (1 - 2 / alpha) + 2;
	unsigned int I, J, a, c, v;
	rk_adjacency(idx, ne, k, nv, off, ent);
	/* element pass */
	for (I = 0; I < ne; I++)
	{
		qaws_scalar p[9], q[9], g[9], e[1], kw[2], cI[9], gI[3] = { 0, 0, 0 }, uI = 0;
		rk_load(x, dim, idx, k, I, p);
		if (k == 2) qaws_repulsive_curve_coefficients(p, cI); else qaws_repulsive_mesh_coefficients(p, cI);
		for (a = 0; a < k; a++)
		{
			uI += u[idx[k * I + a]] / (qaws_scalar)k;
			for (c = 0; c < 3; c++) gI[c] += u[idx[k * I + a]] * cI[3 * a + c];
		}
		for (c = 0; c < 9; c++) coef[9 * I + c] = cI[c];
		for (J = 0; J < ne; J++)
		{
			int shared = 0;
			qaws_scalar gJ[3] = { 0, 0, 0 }, uJ = 0;
			if (J == I)
				continue;
			for (a = 0; a < k; a++)
				for (c = 0; c < k; c++) shared |= idx[k * I + a] == idx[k * J + c];
			if (k == 2 && shared)
				continue;
			rk_load(x, dim, idx, k, J, q);
			if (k == 2) qaws_repulsive_curve_pair(p, q, (qaws_scalar)alpha, (qaws_scalar)beta, 0, e, g); else qaws_repulsive_mesh_pair(p, q, (qaws_scalar)alpha, 0, e, g);
			eb[I] += (double)e[0];
			for (c = 0; c < 3 * k; c++) blocks[9 * I + c] += g[c];
			if (k == 2) qaws_repulsive_curve_pair(q, p, (qaws_scalar)alpha, (qaws_scalar)beta, 1, e, g); else qaws_repulsive_mesh_pair(q, p, (qaws_scalar)alpha, 1, e, g);
			for (c = 0; c < 3 * k; c++) blocks[9 * I + c] += g[c];
			/* the metric's row I */
			if (k == 2) qaws_repulsive_curve_metric_pair(p, q, (qaws_scalar)ex, kw); else qaws_repulsive_mesh_metric_pair(p, q, (qaws_scalar)ex, kw);
			{
				qaws_scalar cJ2[9];
				if (k == 2) qaws_repulsive_curve_coefficients(q, cJ2); else qaws_repulsive_mesh_coefficients(q, cJ2);
				for (a = 0; a < k; a++)
				{
					uJ += u[idx[k * J + a]] / (qaws_scalar)k;
					for (c = 0; c < 3; c++) gJ[c] += u[idx[k * J + a]] * cJ2[3 * a + c];
				}
			}
			for (c = 0; c < 3; c++) acc[4 * I + c] += kw[0] * (gI[c] - gJ[c]);
			acc[4 * I + 3] += kw[1] * (uI - uJ);
		}
	}
	/* vertex gather */
	*E = 0;
	for (I = 0; I < ne; I++) *E += eb[I];
	for (v = 0; v < nv; v++)
	{
		double s = 0, sg[3] = { 0, 0, 0 };
		for (a = off[v]; a < off[v + 1]; a++)
		{
			unsigned int el = ent[a] / k, slot = ent[a] % k;
			for (c = 0; c < 3; c++)
			{
				sg[c] += (double)blocks[9 * el + 3 * slot + c];
				s += (double)coef[9 * el + 3 * slot + c] * (double)acc[4 * el + c];
			}
			s += (double)acc[4 * el + 3] / k;
		}
		for (c = 0; c < dim; c++) grad[dim * v + c] = sg[c];
		y[v] = s;
	}
	free(off); free(ent); free(blocks); free(acc); free(coef); free(eb);
}

static void rk_check(char const* name, qaws_scalar const* x, unsigned int dim, unsigned int nv, unsigned int const* idx, unsigned int ne, unsigned int k)
{
	qaws_scalar* gl = (qaws_scalar*)calloc(dim * nv, sizeof(qaws_scalar));
	qaws_scalar* u = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	double *g = (double*)malloc(sizeof(double) * dim * nv), *gr = (double*)malloc(sizeof(double) * dim * nv), *y = (double*)malloc(sizeof(double) * nv);
	double *xd = (double*)malloc(sizeof(double) * dim * nv), *A = (double*)malloc(sizeof(double) * nv * nv), *yr = (double*)calloc(nv, sizeof(double));
	double E = 0, alpha, beta, eE, eG, eA;
	qaws_scalar El = 0;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_mesh* m = NULL;
	unsigned int i, j;
	char msg[300];
	for (i = 0; i < nv; i++) u[i] = (qaws_scalar)(sin(0.37 * i) + 0.5 * cos(1.3 * i));
	for (i = 0; i < dim * nv; i++) xd[i] = (double)x[i];
	if (k == 2)
	{
		qaws_repulsive_curves_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = dim; d.vertex_count = nv; d.edges = idx; d.edge_count = ne;
		qaws_repulsive_curves_create(&c, &d);
		qaws_repulsive_curves_gradient(gl, &El, c, x, 1);
		qaws_internal_repulsive_curves_metric(c, xd, A);
		alpha = 3; beta = 6;
	}
	else
	{
		qaws_repulsive_mesh_desc d;
		memset(&d, 0, sizeof(d));
		d.vertex_count = nv; d.triangles = idx; d.triangle_count = ne;
		qaws_repulsive_mesh_create(&m, &d);
		qaws_repulsive_mesh_gradient(gl, &El, m, x, 1);
		qaws_internal_repulsive_mesh_metric(m, xd, A);
		alpha = 6; beta = 12;
	}
	for (i = 0; i < nv; i++)
		for (j = 0; j < nv; j++) yr[i] += A[(size_t)i * nv + j] * (double)u[j];
	for (i = 0; i < dim * nv; i++) gr[i] = (double)gl[i];
	rk_kernels(x, dim, nv, idx, ne, k, alpha, beta, u, &E, g, y);
	eE = fabs(E - (double)El) / (double)El;
	eG = rk_rel(g, gr, dim * nv);
	eA = rk_rel(y, yr, nv);
	sprintf(msg, "%s: the kernels against the library: energy %.2g, gradient %.2g, metric product %.2g", name, eE, eG, eA);
	TEST_ASSERT(eE <= (RK_F32 ? 1e-5 : 1e-13) && eG <= (RK_F32 ? 1e-4 : 1e-12) && eA <= (RK_F32 ? 1e-4 : 1e-12), msg);
	qaws_repulsive_curves_destroy(c); qaws_repulsive_mesh_destroy(m);
	free(gl); free(u); free(g); free(gr); free(y); free(xd); free(A); free(yr);
}

static void test_curves(void)
{
	enum { N = 60, M = 12 };
	qaws_scalar x[3 * (N + M)], x2[2 * 20];
	unsigned int e[2 * (N + M)], e2[2 * 19], i;
	for (i = 0; i < N; i++)
	{
		double t = 2 * RK_PI * i / N;
		x[3 * i] = (qaws_scalar)(sin(t) + 2 * sin(2 * t)); x[3 * i + 1] = (qaws_scalar)(cos(t) - 2 * cos(2 * t)); x[3 * i + 2] = (qaws_scalar)(-0.6 * sin(3 * t));
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % N;
	}
	for (i = 0; i < M; i++)
	{
		double t = 2 * RK_PI * i / M;
		x[3 * (N + i)] = (qaws_scalar)(0.4 * cos(t)); x[3 * (N + i) + 1] = (qaws_scalar)(0.4 * sin(t)); x[3 * (N + i) + 2] = (qaws_scalar)1.5;
		e[2 * (N + i)] = N + i; e[2 * (N + i) + 1] = N + (i + 1) % M;
	}
	rk_check("trefoil and a ring (72 edges, 3D)", x, 3, N + M, e, N + M, 2);
	for (i = 0; i < 20; i++)
	{
		double t = (double)i / 19;
		x2[2 * i] = (qaws_scalar)(2 * t); x2[2 * i + 1] = (qaws_scalar)(0.5 * sin(9 * t));
		if (i < 19) { e2[2 * i] = i; e2[2 * i + 1] = i + 1; }
	}
	rk_check("an open wave in 2D (19 edges)", x2, 2, 20, e2, 19, 2);
}

static void test_mesh(void)
{
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[8][3] = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 }, { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
	double p[70][3];
	qaws_scalar x[3 * 70];
	unsigned int tris[3 * 128], t2[3 * 128], n = 6, f = 8, l, i, k, q;
	memcpy(p, ov, sizeof(ov));
	memcpy(tris, of, sizeof(of));
	for (l = 0; l < 2; l++)
	{
		unsigned int nf = 0;
		for (i = 0; i < f; i++)
		{
			unsigned int c[3], mid[3];
			for (k = 0; k < 3; k++) c[k] = tris[3 * i + k];
			for (k = 0; k < 3; k++)
			{
				unsigned int a = c[k], b = c[(k + 1) % 3], j;
				double m[3], len;
				for (q = 0; q < 3; q++) m[q] = 0.5 * (p[a][q] + p[b][q]);
				len = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
				for (q = 0; q < 3; q++) m[q] /= len;
				for (j = 0; j < n && !(fabs(p[j][0] - m[0]) < 1e-12 && fabs(p[j][1] - m[1]) < 1e-12 && fabs(p[j][2] - m[2]) < 1e-12); j++)
					;
				if (j == n) { memcpy(p[n], m, sizeof(m)); n++; }
				mid[k] = j;
			}
			t2[3 * nf] = c[0]; t2[3 * nf + 1] = mid[0]; t2[3 * nf + 2] = mid[2]; nf++;
			t2[3 * nf] = mid[0]; t2[3 * nf + 1] = c[1]; t2[3 * nf + 2] = mid[1]; nf++;
			t2[3 * nf] = mid[2]; t2[3 * nf + 1] = mid[1]; t2[3 * nf + 2] = c[2]; nf++;
			t2[3 * nf] = mid[0]; t2[3 * nf + 1] = mid[1]; t2[3 * nf + 2] = mid[2]; nf++;
		}
		memcpy(tris, t2, sizeof(unsigned int) * 3 * nf);
		f = nf;
	}
	for (i = 0; i < n; i++)
	{
		double r = 1 + 0.15 * sin(3 * p[i][0]) * cos(2 * p[i][1]);
		for (q = 0; q < 3; q++) x[3 * i + q] = (qaws_scalar)(p[i][q] * r);
	}
	rk_check("a lumpy sphere (128 faces)", x, 3, n, tris, f, 3);
}

int test_119_repulsive_core_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 119: backend-neutral tangent-point kernels\n");
	test_curves();
	test_mesh();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
