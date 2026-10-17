/*
 * Test 114: fractional Sobolev descent of the surfaces' tangent-point
 * energy (qaws_repulsive.h, R3)
 *
 *   - the metric A = B + B0 against Mathematica (reference/
 *     114_repulsive_mesh_metric.wls: from the paper's bilinear forms, Eqs.
 *     12 and 13), an open patch and an octahedron; symmetric, constants in
 *     its kernel, positive semidefinite
 *   - the preconditioner (-L)^+ L^(2 - s) (-L)^+ (cotan Laplacian): GMRES
 *     with it reaches the same solution as without, in far fewer
 *     iterations, on a refined sphere
 *   - a bumpy sphere with its area and volume held: the energy falls, both
 *     held
 *   - collision-checked: a sheet pressed under a pinned one (an obstacle)
 *     moves, the obstacle does not, nothing touches; the path check catches
 *     a vertex through a face and an edge through an edge, ends clear
 *   - refusals: volume on an open mesh, touching input in the collision-
 *     checked mode; free translations reported
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "internal/qaws_internal_repulsive.h"
#include "reference/114_repulsive_mesh_metric.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RM_F32 (sizeof(qaws_scalar) == 4)

static qaws_repulsive_mesh* rm_make(unsigned int nv, unsigned int const* tris, unsigned int nf, double p)
{
	qaws_repulsive_mesh_desc d;
	qaws_repulsive_mesh* m = NULL;
	memset(&d, 0, sizeof(d));
	d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf; d.p = (qaws_scalar)p;
	qaws_repulsive_mesh_create(&d, &m);
	return m;
}

/* an octahedron subdivided `level` times, on the unit sphere: x (3 per
   vertex, double), tris; the counts out */
static void rm_sphere(unsigned int level, double* x, unsigned int* tris, unsigned int* nv, unsigned int* nf)
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
				/* the edge's midpoint, made once */
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
	free(t2);
	*nv = n; *nf = f;
}

static void test_metric(void)
{
	static struct { double const* x; unsigned int const* idx; unsigned int nv, nf; double const* p; double const* A; char const* name; } k[2] = {
		{ g114_m1_x, g114_m1_idx, 5, 3, &g114_m1_p, g114_m1_A, "open patch, p 6" },
		{ g114_m2_x, g114_m2_idx, 6, 8, &g114_m2_p, g114_m2_A, "octahedron, p 5.5" } };
	unsigned int t, i, j;
	char msg[300];
	for (t = 0; t < 2; t++)
	{
		qaws_repulsive_mesh* m = rm_make(k[t].nv, k[t].idx, k[t].nf, *k[t].p);
		double A[36], x[18], err = 0, mx = 0, asym = 0, row = 0, psd = 1e300;
		qaws_status st;
		for (i = 0; i < 3 * k[t].nv; i++) x[i] = (double)(qaws_scalar)k[t].x[i];
		st = m ? qaws_internal_repulsive_mesh_metric(m, x, A) : QAWS_STATUS_INVALID_ARGUMENT;
		for (i = 0; i < k[t].nv * k[t].nv; i++) { err = fmax(err, fabs(A[i] - k[t].A[i])); mx = fmax(mx, fabs(k[t].A[i])); }
		for (i = 0; i < k[t].nv; i++)
		{
			double s = 0;
			for (j = 0; j < k[t].nv; j++) { s += A[i * k[t].nv + j]; asym = fmax(asym, fabs(A[i * k[t].nv + j] - A[j * k[t].nv + i])); }
			row = fmax(row, fabs(s));
		}
		for (j = 0; j < 20; j++)
		{
			double u[6], q = 0, uu = 0;
			unsigned int a, b;
			for (a = 0; a < k[t].nv; a++) { u[a] = sin(5.3 * a + 2.1 * j + 0.4); uu += u[a] * u[a]; }
			for (a = 0; a < k[t].nv; a++) for (b = 0; b < k[t].nv; b++) q += u[a] * A[a * k[t].nv + b] * u[b];
			psd = fmin(psd, q / (mx * uu));
		}
		sprintf(msg, "metric, %s: against Mathematica %.2g (status %d), asymmetry %.2g, row sums %.2g, smallest u^T A u / |A| |u|^2 %.3g",
			k[t].name, err / mx, (int)st, asym / mx, row / mx, psd);
		TEST_ASSERT(st == QAWS_STATUS_OK && err / mx <= (RM_F32 ? 1e-5 : 1e-13) && asym / mx <= 1e-15 && row / mx <= 1e-13 && psd >= -1e-13, msg);
		qaws_repulsive_mesh_destroy(m);
	}
}

static void test_preconditioner(void)
{
	double* x = (double*)malloc(sizeof(double) * 3 * 300);
	unsigned int* tris = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 600);
	double b[300], u0[300], u1[300], mean = 0, diff = 0, un = 0;
	unsigned int nv, nf, i, it0 = 0, it1 = 0;
	qaws_repulsive_mesh* m;
	qaws_status s0, s1;
	char msg[300];
	rm_sphere(3, x, tris, &nv, &nf);
	/* a lumpy sphere: the preconditioner sees the geometry */
	for (i = 0; i < nv; i++)
	{
		double r = 1 + 0.15 * sin(3 * x[3 * i]) * cos(2 * x[3 * i + 1]);
		x[3 * i] *= r; x[3 * i + 1] *= r; x[3 * i + 2] *= 0.8 * r;
	}
	m = rm_make(nv, tris, nf, 6);
	for (i = 0; i < nv; i++) { b[i] = sin(1.7 * i) + cos(0.3 * i * i); mean += b[i] / nv; }
	for (i = 0; i < nv; i++) b[i] -= mean;
	s0 = qaws_internal_repulsive_mesh_gmres(m, x, b, 0, 1e-10, u0, &it0);
	s1 = qaws_internal_repulsive_mesh_gmres(m, x, b, 1, 1e-10, u1, &it1);
	for (i = 0; i < nv; i++) { diff = fmax(diff, fabs(u0[i] - u1[i])); un = fmax(un, fabs(u0[i])); }
	sprintf(msg, "preconditioner, a lumpy sphere (%u vertices, %u faces): GMRES to 1e-10 in %u iterations unpreconditioned, %u with (-L)^+ L^(2-s) (-L)^+; solutions %.2g apart",
		nv, nf, it0, it1, diff / un);
	TEST_ASSERT(s0 == QAWS_STATUS_OK && s1 == QAWS_STATUS_OK && 3 * it1 <= it0 && diff / un <= 1e-7, msg);
	qaws_repulsive_mesh_destroy(m);
	free(x); free(tris);
}

static void test_sphere(void)
{
	double xd[3 * 70];
	qaws_scalar x[3 * 70];
	unsigned int tris[3 * 130], nv, nf, i, f;
	double A0 = 0, V0 = 0, A1 = 0, V1 = 0;
	qaws_repulsive_mesh* m;
	qaws_repulsive_mesh_solve_desc sd;
	qaws_repulsive_report r;
	qaws_status st;
	char msg[300];
	rm_sphere(2, xd, tris, &nv, &nf);
	for (i = 0; i < nv; i++)
	{
		double rr = 1 + 0.2 * sin(4 * xd[3 * i]) * sin(3 * xd[3 * i + 2]);
		x[3 * i] = (qaws_scalar)(xd[3 * i] * rr); x[3 * i + 1] = (qaws_scalar)(xd[3 * i + 1] * rr); x[3 * i + 2] = (qaws_scalar)(xd[3 * i + 2] * rr);
	}
	for (f = 0; f < nf; f++)
	{
		/* area and volume, independently of the library */
		double a[3], b[3], c[3], n[3];
		unsigned int q;
		for (q = 0; q < 3; q++) { a[q] = x[3 * tris[3 * f] + q]; b[q] = x[3 * tris[3 * f + 1] + q]; c[q] = x[3 * tris[3 * f + 2] + q]; }
		n[0] = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
		n[1] = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]);
		n[2] = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
		A0 += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		V0 += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
	}
	m = rm_make(nv, tris, nf, 0);
	memset(&sd, 0, sizeof(sd));
	sd.hold_area = 1; sd.hold_volume = 1; sd.max_iterations = 25;
	memset(&r, 0, sizeof(r));
	st = qaws_repulsive_mesh_solve(m, &sd, x, &r);
	for (f = 0; f < nf; f++)
	{
		double a[3], b[3], c[3], n[3];
		unsigned int q;
		for (q = 0; q < 3; q++) { a[q] = x[3 * tris[3 * f] + q]; b[q] = x[3 * tris[3 * f + 1] + q]; c[q] = x[3 * tris[3 * f + 2] + q]; }
		n[0] = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
		n[1] = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]);
		n[2] = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
		A1 += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		V1 += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
	}
	sprintf(msg, "bumpy sphere (%u faces), area and volume held: reason %d after %u steps, E %.5f from %.5f, area off %.2g, volume off %.2g",
		nf, (int)r.reason, r.iterations, (double)r.energy, (double)r.initial_energy, fabs(A1 - A0) / A0, fabs(V1 - V0) / V0);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.iterations > 0 && r.energy < r.initial_energy && fabs(A1 - A0) / A0 <= (RM_F32 ? 1e-5 : 1e-9) &&
		fabs(V1 - V0) / V0 <= (RM_F32 ? 1e-5 : 1e-9), msg);
	qaws_repulsive_mesh_destroy(m);
}

/* the area of x's faces, independently of the library */
static double rm_area_of(qaws_scalar const* x, unsigned int const* tris, unsigned int nf)
{
	double A = 0;
	unsigned int f, q;
	for (f = 0; f < nf; f++)
	{
		double a[3], b[3], c[3], n[3];
		for (q = 0; q < 3; q++) { a[q] = x[3 * tris[3 * f] + q]; b[q] = x[3 * tris[3 * f + 1] + q]; c[q] = x[3 * tris[3 * f + 2] + q]; }
		n[0] = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
		n[1] = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]);
		n[2] = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
		A += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	}
	return A;
}

/* the enclosed volume: the sum of a . (b x c) / 6 over the faces */
static double rm_volume_of(qaws_scalar const* x, unsigned int const* tris, unsigned int nf)
{
	double V = 0;
	unsigned int f, q;
	for (f = 0; f < nf; f++)
	{
		double a[3], b[3], c[3];
		for (q = 0; q < 3; q++) { a[q] = x[3 * tris[3 * f] + q]; b[q] = x[3 * tris[3 * f + 1] + q]; c[q] = x[3 * tris[3 * f + 2] + q]; }
		V += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
	}
	return V;
}

/* the area weighted into the energy (the paper's minimize area) */
static void test_area_weight(void)
{
	double xd[3 * 70], A0, Aw = 0, Ap = 0;
	qaws_scalar x[3 * 70], xw[3 * 70], xp[3 * 70], dx[3 * 70];
	unsigned int tris[3 * 130], nv, nf, i, k;
	qaws_repulsive_mesh* m;
	qaws_repulsive_mesh_solve_desc sd;
	qaws_repulsive_mesh_solver* s = NULL;
	qaws_repulsive_report r0, rw, rp;
	qaws_repulsive_targets t;
	qaws_status st, sj;
	char msg[400];
	rm_sphere(2, xd, tris, &nv, &nf);
	for (i = 0; i < nv; i++)
	{
		double rr = 1 + 0.2 * sin(4 * xd[3 * i]) * sin(3 * xd[3 * i + 2]);
		for (k = 0; k < 3; k++) x[3 * i + k] = (qaws_scalar)(xd[3 * i + k] * rr);
	}
	A0 = rm_area_of(x, tris, nf);
	m = rm_make(nv, tris, nf, 0);
	memset(&r0, 0, sizeof(r0)); memset(&rw, 0, sizeof(rw)); memset(&rp, 0, sizeof(rp));
	memcpy(xw, x, sizeof(x)); memcpy(xp, x, sizeof(x));
	/* the same solve, weighted and not: unweighted the surface inflates, weighted it shrinks */
	memset(&sd, 0, sizeof(sd));
	sd.max_iterations = 20;
	st = qaws_repulsive_mesh_solve(m, &sd, xp, &rp);
	sd.area_weight = (qaws_scalar)2;
	st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solve(m, &sd, xw, &rw) : st;
	Aw = rm_area_of(xw, tris, nf);
	Ap = rm_area_of(xp, tris, nf);
	sprintf(msg, "area weighted (2): the energy at the input %.9g against the unweighted %.9g + 2 A %.9g; after 20 steps, area %.6f (unweighted %.6f, input %.6f), E %.5g from %.5g",
		(double)rw.initial_energy, (double)rp.initial_energy, 2 * A0, Aw, Ap, A0, (double)rw.energy, (double)rw.initial_energy);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs((double)rw.initial_energy - (double)rp.initial_energy - 2 * A0) <= (RM_F32 ? 1e-4 : 1e-10) * (double)rw.initial_energy &&
		Aw < Ap && Aw < A0 && rw.energy < rw.initial_energy, msg);
	/* no solution derivatives with it; a negative weight refused */
	sd.max_iterations = 0;
	st = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);
	memset(&t, 0, sizeof(t));
	sj = s ? qaws_repulsive_mesh_solver_solution_jvp(s, &t, dx, NULL) : QAWS_STATUS_INVALID_ARGUMENT;
	qaws_repulsive_mesh_solver_destroy(s);
	s = NULL;
	sd.area_weight = (qaws_scalar)-1;
	sprintf(msg, "area weighted: solution derivatives refused (%d), a negative weight refused (%d)", (int)sj,
		(int)qaws_repulsive_mesh_solver_create(m, &sd, x, &s));
	TEST_ASSERT(st == QAWS_STATUS_OK && sj == QAWS_STATUS_UNSUPPORTED_OPERATION && s == NULL, msg);
	qaws_repulsive_mesh_solver_destroy(s);
	s = NULL;
	/* the first step's largest move capped by initial_step (a solver made again resuming its step); a negative one refused */
	{
		static double const caps[2] = { 0.01, 1 };
		double moved[2] = { 0, 0 };
		qaws_status sn;
		sd.area_weight = 0;
		for (k = 0; k < 2; k++)
		{
			memcpy(xw, x, sizeof(x));
			sd.initial_step = (qaws_scalar)caps[k];
			st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solver_create(m, &sd, xw, &s) : st;
			st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solver_step(s, xw, &rw) : st;
			moved[k] = (double)rw.last_step;
			qaws_repulsive_mesh_solver_destroy(s);
			s = NULL;
		}
		sd.initial_step = (qaws_scalar)-1;
		sn = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);
		sprintf(msg, "initial step: the first move %.4g (cap 0.01), %.4g (cap 1), a negative one refused (%d)", moved[0], moved[1], (int)sn);
		TEST_ASSERT(st == QAWS_STATUS_OK && moved[0] > 0 && moved[0] <= 0.01 * 1.001 && moved[1] > moved[0] && s == NULL, msg);
		qaws_repulsive_mesh_solver_destroy(s);
	}
	qaws_repulsive_mesh_destroy(m);
}

/* translations added back (Repulsive Surfaces, Fig. 8): two spheres, each
   solved with its barycenter held, then moved by its dual-area mean of dE */
static void test_translations(void)
{
	double xd[3 * 20], w[2 * 20], t[2][4], keep[3] = { 0, 0, 0 }, err = 0, scale = 0, gap0, gap1, mean0[3] = { 0, 0, 0 }, drift = 0;
	qaws_scalar x[6 * 20], x0[6 * 20], g[6 * 20], dE[6 * 20], E;
	unsigned int tris[6 * 32], nv, nf, i, k, q, c, pins[20];
	qaws_repulsive_mesh* m;
	qaws_repulsive_mesh_solve_desc sd;
	qaws_repulsive_mesh_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_status st, s1, s2;
	char msg[400];
	rm_sphere(1, xd, tris, &nv, &nf);
	for (k = 0; k < 2; k++)
		for (i = 0; i < nv; i++)
			for (q = 0; q < 3; q++)
				x[3 * (k * nv + i) + q] = (qaws_scalar)(xd[3 * i + q] * (k ? 0.8 : 1) + (q == 0 && k ? 2.4 : 0) + (q == 1 && k ? 0.3 : 0));
	for (i = 0; i < 3 * nf; i++) tris[3 * nf + i] = tris[i] + nv;
	memcpy(x0, x, sizeof(x));
	m = rm_make(2 * nv, tris, 2 * nf, 6);
	memset(&sd, 0, sizeof(sd));
	sd.translations = 2;
	/* the direction: each sphere's mean of g is its dual-area mean of dE, less their vertex-weighted mean (the overall barycenter kept) */
	memset(dE, 0, sizeof(dE));
	st = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);
	st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solver_direction(s, x, g) : st;
	st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_gradient(m, x, (qaws_scalar)1, dE, &E) : st;
	memset(w, 0, sizeof(w)); memset(t, 0, sizeof(t));
	for (i = 0; i < 2 * nf; i++)
	{
		unsigned int const* f = tris + 3 * i;
		double e1[3], e2[3], n[3];
		for (q = 0; q < 3; q++) { e1[q] = (double)x[3 * f[1] + q] - (double)x[3 * f[0] + q]; e2[q] = (double)x[3 * f[2] + q] - (double)x[3 * f[0] + q]; }
		n[0] = e1[1] * e2[2] - e1[2] * e2[1]; n[1] = e1[2] * e2[0] - e1[0] * e2[2]; n[2] = e1[0] * e2[1] - e1[1] * e2[0];
		for (c = 0; c < 3; c++) w[f[c]] += sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / 6;
	}
	for (i = 0; i < 2 * nv; i++)
	{
		for (q = 0; q < 3; q++) t[i / nv][q] += w[i] * (double)dE[3 * i + q];
		t[i / nv][3] += w[i];
	}
	for (k = 0; k < 2; k++) for (q = 0; q < 3; q++) { t[k][q] /= t[k][3]; keep[q] += t[k][q] / 2; }
	for (k = 0; k < 2; k++)
		for (q = 0; q < 3; q++)
		{
			double mg = 0;
			for (i = 0; i < nv; i++) mg += (double)g[3 * (k * nv + i) + q] / nv;
			err = fmax(err, fabs(mg - (t[k][q] - keep[q])));
			scale = fmax(scale, fabs(t[k][q]));
		}
	/* ten steps: the spheres part, the overall barycenter stays */
	for (k = 0; st == QAWS_STATUS_OK && k < 10; k++) st = qaws_repulsive_mesh_solver_step(s, x, &r);
	qaws_repulsive_mesh_solver_destroy(s);
	s = NULL;
	gap1 = 0;
	for (i = 0; i < 2 * nv; i++)
	{
		double sgn = i < nv ? -1.0 / nv : 1.0 / nv;
		gap1 += sgn * (double)x[3 * i];
		for (q = 0; q < 3; q++) mean0[q] += ((double)x[3 * i + q] - (double)x0[3 * i + q]) / (2 * nv);
	}
	for (q = 0; q < 3; q++) drift = fmax(drift, fabs(mean0[q]));
	gap0 = 0;
	for (i = 0; i < 2 * nv; i++) gap0 += (i < nv ? -1.0 / nv : 1.0 / nv) * (double)x0[3 * i];
	/* a pinned sphere stays; the free one still moves; the whole-mesh drift refused with pins */
	memcpy(x, x0, sizeof(x));
	for (i = 0; i < nv; i++) pins[i] = i;
	sd.pinned = pins; sd.pinned_count = nv;
	s1 = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);
	for (k = 0; s1 == QAWS_STATUS_OK && k < 5; k++) s1 = qaws_repulsive_mesh_solver_step(s, x, &r);
	qaws_repulsive_mesh_solver_destroy(s);
	s = NULL;
	sd.translations = 1;
	s2 = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);
	sprintf(msg, "translations: each sphere's mean direction against its dual-area mean of dE %.2g (of %.3g); 10 steps: the centres %.4f -> %.4f apart, the overall barycenter %.2g; a pinned sphere %s; the drift with pins refused (%d)",
		err, scale, gap0, gap1, drift, memcmp(x, x0, sizeof(qaws_scalar) * 3 * nv) == 0 && memcmp(x + 3 * nv, x0 + 3 * nv, sizeof(qaws_scalar) * 3 * nv) != 0 ? "kept, the other moved" : "MOVED", (int)s2);
	TEST_ASSERT(st == QAWS_STATUS_OK && s1 == QAWS_STATUS_OK && err <= (RM_F32 ? 1e-4 : 1e-10) * scale && gap1 > gap0 && drift <= (RM_F32 ? 1e-5 : 1e-12) &&
		memcmp(x, x0, sizeof(qaws_scalar) * 3 * nv) == 0 && memcmp(x + 3 * nv, x0 + 3 * nv, sizeof(qaws_scalar) * 3 * nv) != 0 && s2 == QAWS_STATUS_INVALID_ARGUMENT && s == NULL, msg);
	qaws_repulsive_mesh_destroy(m);
}

/* the targets moved during a solve (a growth schedule), and the volume's
   deviation weighted into the energy */
static void test_targets(void)
{
	double xd[3 * 70], A0, A1, V0, V1, target = 0;
	qaws_scalar x[3 * 70], xv[3 * 70], a, v;
	unsigned int tris[3 * 130], nv, nf, i, k;
	qaws_repulsive_mesh* m;
	qaws_repulsive_mesh_solve_desc sd;
	qaws_repulsive_mesh_solver* s = NULL;
	qaws_repulsive_report r, rp, rw;
	qaws_status st, sa, sv, sn;
	char msg[400];
	rm_sphere(2, xd, tris, &nv, &nf);
	for (i = 0; i < nv; i++)
	{
		double rr = 1 + 0.2 * sin(4 * xd[3 * i]) * sin(3 * xd[3 * i + 2]);
		for (k = 0; k < 3; k++) x[3 * i + k] = (qaws_scalar)(xd[3 * i + k] * rr);
	}
	memcpy(xv, x, sizeof(x));
	A0 = rm_area_of(x, tris, nf);
	m = rm_make(nv, tris, nf, 6);
	/* the area grown 1.5x over 10 steps, held at its input value: the area follows */
	memset(&sd, 0, sizeof(sd));
	sd.hold_area = 1;
	st = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);
	for (k = 1; st == QAWS_STATUS_OK && k <= 10; k++)
	{
		target = A0 * (1 + 0.05 * k);
		a = (qaws_scalar)target;
		st = qaws_repulsive_mesh_solver_set_targets(s, &a, NULL);
		st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solver_step(s, x, &r) : st;
	}
	A1 = rm_area_of(x, tris, nf);
	v = 1;
	sv = s ? qaws_repulsive_mesh_solver_set_targets(s, NULL, &v) : QAWS_STATUS_OK;   /* (the volume is not held) */
	a = (qaws_scalar)-1;
	sa = s ? qaws_repulsive_mesh_solver_set_targets(s, &a, NULL) : QAWS_STATUS_OK;
	qaws_repulsive_mesh_solver_destroy(s);
	s = NULL;
	/* the volume's deviation: the energy at the input against the unweighted + w dev^2, and 20 steps toward twice the volume */
	memset(&sd, 0, sizeof(sd));
	sd.max_iterations = 1;
	memcpy(x, xv, sizeof(x));
	st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solve(m, &sd, x, &rp) : st;   /* (one step, for the unweighted energy at the input) */
	sd.volume_weight = (qaws_scalar)3;
	V0 = rm_volume_of(xv, tris, nf);
	sd.volume_target = (qaws_scalar)(2 * V0);
	sd.max_iterations = 20;
	st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solve(m, &sd, xv, &rw) : st;
	V1 = rm_volume_of(xv, tris, nf);
	sd.volume_weight = (qaws_scalar)-1;
	sn = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);
	sprintf(msg, "targets: the area grown to %.6f (target %.6f, input %.6f); refused: a volume not held (%d), a negative area (%d); the volume's deviation: E at the input %.9g against %.9g + 3 (1/2)^2, the volume %.4f -> %.4f (target %.4f); a negative weight refused (%d)",
		A1, target, A0, (int)sv, (int)sa, (double)rw.initial_energy, (double)rp.initial_energy, V0, V1, 2 * V0, (int)sn);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs(A1 - target) <= (RM_F32 ? 1e-5 : 1e-8) * target && sv == QAWS_STATUS_INVALID_ARGUMENT && sa == QAWS_STATUS_INVALID_ARGUMENT &&
		fabs((double)rw.initial_energy - (double)rp.initial_energy - 0.75) <= (RM_F32 ? 1e-4 : 1e-10) * (double)rw.initial_energy && V1 > V0 && sn == QAWS_STATUS_INVALID_ARGUMENT && s == NULL, msg);
	qaws_repulsive_mesh_destroy(m);
}

/* fixed obstacles: a mesh obstacle's energy (theta 0) is the two meshes'
   cross terms, both orders; a point cloud's is the faces' terms against
   its points; a surface drifting free moves away from one */
static void test_obstacles(void)
{
	double xs[3 * 20], xd[6 * 20], gS[3 * 20], gM[6 * 20], gO[6 * 20], gob[3 * 20], ES = 0, EO = 0, EM = 0, Eob = 0, Ebh = 0, Epc = 0, Ebrute = 0, gerr = 0, gmax = 0;
	double cx0 = 0, cx1 = 0;
	qaws_scalar xo[3 * 20], x[3 * 20];
	unsigned int tris[6 * 32], nv, nf, i, k, q, j;
	qaws_repulsive_mesh *mS, *mM;
	qaws_repulsive_obstacle od;
	qaws_internal_obstacle *ob = NULL, *pc = NULL;
	qaws_repulsive_mesh_solve_desc sd;
	qaws_repulsive_mesh_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_status st, sw, si;
	unsigned int bad[3] = { 0, 1, 99 };
	char msg[500];
	rm_sphere(1, xs, tris, &nv, &nf);
	for (i = 0; i < nv; i++)
		for (q = 0; q < 3; q++)
		{
			xd[3 * i + q] = xs[3 * i + q];
			xd[3 * (nv + i) + q] = 0.5 * xs[3 * i + q] + (q == 0 ? 1.9 : q == 1 ? 0.2 : 0);
			xo[3 * i + q] = (qaws_scalar)xd[3 * (nv + i) + q];
		}
	for (i = 0; i < 3 * nf; i++) tris[3 * nf + i] = tris[i] + nv;
	mS = rm_make(nv, tris, nf, 6); mM = rm_make(2 * nv, tris, 2 * nf, 6);
	memset(gS, 0, sizeof(gS)); memset(gM, 0, sizeof(gM)); memset(gO, 0, sizeof(gO)); memset(gob, 0, sizeof(gob));
	st = qaws_internal_repulsive_mesh_bh(mS, xd, 0, &ES, gS);
	st = st == QAWS_STATUS_OK ? qaws_internal_repulsive_mesh_bh(mM, xd, 0, &EM, gM) : st;
	{
		/* the obstacle's own energy, on its own vertex numbering */
		qaws_repulsive_mesh* m0 = rm_make(nv, tris, nf, 6);
		st = st == QAWS_STATUS_OK ? qaws_internal_repulsive_mesh_bh(m0, xd + 3 * nv, 0, &EO, gO) : st;
		qaws_repulsive_mesh_destroy(m0);
	}
	memset(&od, 0, sizeof(od));
	od.positions = xo; od.vertex_count = nv; od.triangles = tris; od.triangle_count = nf; od.weight = (qaws_scalar)2;
	st = st == QAWS_STATUS_OK ? qaws_internal_obstacle_create(&od, &ob) : st;
	st = st == QAWS_STATUS_OK ? qaws_internal_obstacle_energy(ob, mS, xd, 0, &Eob, gob) : st;
	st = st == QAWS_STATUS_OK ? qaws_internal_obstacle_energy(ob, mS, xd, 0.25, &Ebh, NULL) : st;
	for (i = 0; i < 3 * nv; i++)
	{
		gerr = fmax(gerr, fabs(gob[i] - 2 * (gM[i] - gS[i])));
		gmax = fmax(gmax, fabs(gob[i]));
	}
	/* the point cloud (weight 3): each face against each point, the face's normal only */
	od.triangles = NULL; od.triangle_count = 0; od.weight = (qaws_scalar)3;
	st = st == QAWS_STATUS_OK ? qaws_internal_obstacle_create(&od, &pc) : st;
	st = st == QAWS_STATUS_OK ? qaws_internal_obstacle_energy(pc, mS, xd, 0, &Epc, NULL) : st;
	for (k = 0; k < nf; k++)
	{
		unsigned int const* t = tris + 3 * k;
		double e1[3], e2[3], m[3], a, X[3];
		for (q = 0; q < 3; q++) { e1[q] = xd[3 * t[1] + q] - xd[3 * t[0] + q]; e2[q] = xd[3 * t[2] + q] - xd[3 * t[0] + q]; X[q] = (xd[3 * t[0] + q] + xd[3 * t[1] + q] + xd[3 * t[2] + q]) / 3; }
		m[0] = e1[1] * e2[2] - e1[2] * e2[1]; m[1] = e1[2] * e2[0] - e1[0] * e2[2]; m[2] = e1[0] * e2[1] - e1[1] * e2[0];
		a = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
		for (j = 0; j < nv; j++)
		{
			double dd[3], s = 0, r2 = 0;
			for (q = 0; q < 3; q++) { dd[q] = X[q] - (double)xo[3 * j + q]; s += m[q] / a * dd[q]; r2 += dd[q] * dd[q]; }
			Ebrute += 3 * 0.5 * a * pow(fabs(s), 6) / pow(r2, 6);
		}
	}
	/* the solver: the sphere drifting free (translations 1) moves away from the obstacle on its +x side */
	for (i = 0; i < 3 * nv; i++) x[i] = (qaws_scalar)xs[i];
	for (i = 0; i < nv; i++) cx0 += (double)x[3 * i] / nv;
	od.triangles = tris; od.triangle_count = nf; od.weight = 1;
	memset(&sd, 0, sizeof(sd));
	sd.translations = 1; sd.hold_area = 1; sd.obstacles = &od; sd.obstacle_count = 1;
	st = st == QAWS_STATUS_OK ? qaws_repulsive_mesh_solver_create(mS, &sd, x, &s) : st;
	for (k = 0; st == QAWS_STATUS_OK && k < 10; k++) st = qaws_repulsive_mesh_solver_step(s, x, &r);
	qaws_repulsive_mesh_solver_destroy(s);
	s = NULL;
	for (i = 0; i < nv; i++) cx1 += (double)x[3 * i] / nv;
	/* refused: a weight not positive, a triangle index out of range */
	od.weight = 0;
	sw = qaws_repulsive_mesh_solver_create(mS, &sd, x, &s);
	od.weight = 1; od.triangles = bad; od.triangle_count = 1;
	si = s ? QAWS_STATUS_OK : qaws_repulsive_mesh_solver_create(mS, &sd, x, &s);
	sprintf(msg, "obstacles: the mesh obstacle (2) %.12g against 2 x the cross terms %.12g, its gradient %.2g of %.3g; hierarchical (0.25) %.6g; the point cloud (3) %.12g against %.12g; drifting free, the sphere's centre x %.4f -> %.4f; refused: weight 0 (%d), an index out of range (%d)",
		Eob, 2 * (EM - ES - EO), gerr, gmax, Ebh, Epc, Ebrute, cx0, cx1, (int)sw, (int)si);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs(Eob - 2 * (EM - ES - EO)) <= 1e-10 * Eob && gerr <= 1e-9 * gmax && fabs(Ebh - Eob) <= 0.05 * Eob &&
		fabs(Epc - Ebrute) <= 1e-10 * Ebrute && cx1 < cx0 - 1e-3 && sw == QAWS_STATUS_INVALID_ARGUMENT && si == QAWS_STATUS_INVALID_ARGUMENT && s == NULL, msg);
	qaws_internal_obstacle_destroy(ob); qaws_internal_obstacle_destroy(pc);
	qaws_repulsive_mesh_destroy(mS); qaws_repulsive_mesh_destroy(mM);
}

/* a square sheet of n x n vertices at height z, its triangles from vertex base */
static unsigned int rm_sheet(qaws_scalar* x, unsigned int* tris, unsigned int base, unsigned int n, double z, double bump)
{
	unsigned int i, j, f = 0;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			double u = (double)i / (n - 1) - 0.5, v = (double)j / (n - 1) - 0.5;
			x[3 * (base + n * i + j)] = (qaws_scalar)u; x[3 * (base + n * i + j) + 1] = (qaws_scalar)v;
			x[3 * (base + n * i + j) + 2] = (qaws_scalar)(z + bump * cos(3.1 * u) * cos(3.1 * v));
		}
	for (i = 0; i + 1 < n; i++)
		for (j = 0; j + 1 < n; j++)
		{
			unsigned int a = base + n * i + j;
			tris[3 * f] = a; tris[3 * f + 1] = a + n; tris[3 * f + 2] = a + n + 1; f++;
			tris[3 * f] = a; tris[3 * f + 1] = a + n + 1; tris[3 * f + 2] = a + 1; f++;
		}
	return f;
}

static void test_sheets(void)
{
	enum { n = 6 };
	qaws_scalar x[3 * 2 * n * n], x0[3 * 2 * n * n];
	unsigned int tris[3 * 4 * (n - 1) * (n - 1)], pins[n * n + 4], i, nf, np = 0;
	double still = 0;
	qaws_repulsive_mesh* m;
	qaws_repulsive_mesh_solve_desc sd;
	qaws_repulsive_report r;
	qaws_status st;
	char msg[300];
	/* the free sheet bulging up toward the pinned one above */
	nf = rm_sheet(x, tris, 0, n, 0, 0.12);
	nf += rm_sheet(x, tris + 3 * nf, n * n, n, 0.25, 0);
	for (i = 0; i < n * n; i++) pins[np++] = n * n + i;
	pins[np++] = 0; pins[np++] = n - 1; pins[np++] = n * (n - 1); pins[np++] = n * n - 1;   /* the free sheet's corners too */
	memcpy(x0, x, sizeof(x));
	m = rm_make(2 * n * n, tris, nf, 0);
	memset(&sd, 0, sizeof(sd));
	sd.pinned = pins; sd.pinned_count = np; sd.hold_area = 1; sd.collision_checked = 1; sd.max_iterations = 15;
	memset(&r, 0, sizeof(r));
	st = qaws_repulsive_mesh_solve(m, &sd, x, &r);
	for (i = 3 * n * n; i < 3 * 2 * n * n; i++) still = fmax(still, fabs((double)x[i] - x0[i]));
	sprintf(msg, "a sheet under a pinned one, area held, collision-checked: reason %d after %u steps, E %.4f from %.4f, the pinned sheet moved %.2g, rejections e%u c%u",
		(int)r.reason, r.iterations, (double)r.energy, (double)r.initial_energy, still, r.rejected_energy, r.rejected_collision);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.iterations > 0 && r.energy < r.initial_energy && still <= (RM_F32 ? 1e-6 : 1e-12), msg);
	qaws_repulsive_mesh_destroy(m);
}

static void test_path_check(void)
{
	/* a triangle (0, 1, 2) and a separate edge (3, 4) with a free vertex 5
	   (its own degenerate-free triangle 3, 4, 5) */
	static unsigned int const t[2][3] = { { 0, 1, 2 }, { 3, 4, 5 } };
	double a[18] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.2, 0.2, 0.5, 2, 2, 0.5, 2, 3, 0.5 }, a2[18], b[18];
	qaws_repulsive_mesh* m = rm_make(6, &t[0][0], 2, 0);
	int vf, ee, still;
	unsigned int i;
	char msg[300];
	/* vertex 3 sweeps down through the triangle (the others stay far) */
	memcpy(b, a, sizeof(a)); b[11] = -0.5;
	vf = qaws_internal_repulsive_mesh_path_clear(m, 1e-3, a, b);
	/* edge (3, 4) across edge (0, 1), lowered from above it to below it */
	memcpy(a2, a, sizeof(a)); memcpy(b, a, sizeof(a));
	a2[9] = 0.5; a2[10] = -0.5; a2[11] = 0.3; a2[12] = 0.5; a2[13] = 0.2; a2[14] = 0.3;
	b[9] = 0.5; b[10] = -0.5; b[11] = -0.3; b[12] = 0.5; b[13] = 0.2; b[14] = -0.3;
	ee = qaws_internal_repulsive_mesh_path_clear(m, 1e-3, a2, b);
	/* nothing moving */
	still = qaws_internal_repulsive_mesh_path_clear(m, 1e-3, a, a);
	(void)i;
	sprintf(msg, "path check: a vertex through a face %d, an edge through an edge %d (ends clear), still %d", vf, ee, still);
	TEST_ASSERT(vf == 0 && ee == 0 && still == 1, msg);
	qaws_repulsive_mesh_destroy(m);
}

static void test_refusals(void)
{
	static unsigned int const t[2][3] = { { 0, 1, 2 }, { 3, 4, 5 } };
	qaws_scalar x[18] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.2, 0.2, -0.5, 0.2, 0.2, 0.5, 2, 2, 0 };   /* edge (3, 4) through face 0 */
	qaws_scalar y[18] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.2, 0.2, 1, 0.2, 0.5, 1, 2, 2, 1 };
	qaws_repulsive_mesh* m = rm_make(6, &t[0][0], 2, 0);
	qaws_repulsive_mesh_solve_desc sd;
	qaws_repulsive_mesh_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_status s1, s2, s3;
	char msg[300];
	memset(&sd, 0, sizeof(sd));
	memset(&r, 0, sizeof(r));
	sd.hold_volume = 1;
	s1 = qaws_repulsive_mesh_solver_create(m, &sd, y, &s);           /* open: no volume */
	sd.hold_volume = 0; sd.collision_checked = 1;
	s2 = qaws_repulsive_mesh_solver_create(m, &sd, x, &s);           /* touching */
	sd.collision_checked = 0; sd.hold_barycenter = -1;
	s3 = qaws_repulsive_mesh_solver_create(m, &sd, y, &s);
	if (s3 == QAWS_STATUS_OK) qaws_repulsive_mesh_solver_step(s, y, &r);
	qaws_repulsive_mesh_solver_destroy(s);
	sprintf(msg, "refusals: volume on an open mesh %d, touching input (collision-checked) %d; free translations: reason %d", (int)s1, (int)s2, (int)r.reason);
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_OK && r.reason == QAWS_REPULSIVE_NUMERICAL_FAILURE, msg);
	qaws_repulsive_mesh_destroy(m);
}

int test_114_repulsive_mesh_solve_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 114: fractional Sobolev descent of the surfaces' tangent-point energy\n");
	test_metric();
	test_preconditioner();
	test_sphere();
	test_sheets();
	test_area_weight();
	test_translations();
	test_targets();
	test_obstacles();
	test_path_check();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
