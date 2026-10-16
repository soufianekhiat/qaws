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
	qaws_repulsive_mesh_destroy(m);
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
	test_path_check();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
