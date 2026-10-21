/*
 * Test 170: repulsive curves on a mesh, distances along it
 * (qaws_repulsive_mesh_curves.h), small meshes
 *
 *   - euclidean: the energy and its gradient are the curves' own at the
 *     same positions;
 *   - geodesic on a flat grid: the energy close to the curves' (the
 *     chord is kept there);
 *   - a strip folded over, a loop on each sheet facing the other across a
 *     small gap: the chords see them close, the surface far apart: the
 *     geodesic energy a small part of the euclidean one; its gradient
 *     against central differences of moves inside the triangles;
 *   - an ellipse on a flat grid, its length held: rounder, the energy
 *     lower, its vertices on the grid, walked across triangles;
 *   - refusals: a location off the triangles, theta set.
 */

#include "test_common.h"
#include "qaws_repulsive_mesh_curves.h"
#include "qaws_mesh_reference.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RC_PI 3.14159265358979323846

static qaws_mesh_3d* rc_mesh(qaws_scalar const* x, unsigned int nv, unsigned int const* t, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

/* a grid of nu x nw cells mapped by f(u in [0, 1], w in [0, 1]) */
static qaws_mesh_3d* rc_grid(unsigned int nu, unsigned int nw, void (*f)(double u, double w, double* p))
{
	unsigned int i, j, k = 0, nv = (nu + 1) * (nw + 1);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nw);
	qaws_mesh_3d* m;
	for (j = 0; j <= nw; j++)
		for (i = 0; i <= nu; i++)
		{
			double p[3];
			f((double)i / nu, (double)j / nw, p);
			x[3 * (j * (nu + 1) + i)] = (qaws_scalar)p[0];
			x[3 * (j * (nu + 1) + i) + 1] = (qaws_scalar)p[1];
			x[3 * (j * (nu + 1) + i) + 2] = (qaws_scalar)p[2];
		}
	for (j = 0; j < nw; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int v = j * (nu + 1) + i;
			t[k++] = v; t[k++] = v + 1; t[k++] = v + nu + 2;
			t[k++] = v; t[k++] = v + nu + 2; t[k++] = v + nu + 1;
		}
	m = rc_mesh(x, nv, t, 2 * nu * nw);
	free(x);
	free(t);
	return m;
}

static void rc_flat(double u, double w, double* p)
{
	p[0] = 2 * u - 1; p[1] = 2 * w - 1; p[2] = 0;
}

/* a strip 1 wide folded over: z = 0 from x = 0 to 2, a half circle of radius 0.075, back at z = 0.15 */
#define RC_GAP 0.15
static void rc_fold(double u, double w, double* p)
{
	double rr = 0.5 * RC_GAP, L = 4 + RC_PI * rr, s = u * L;
	p[1] = w;
	if (s <= 2) { p[0] = s; p[2] = 0; }
	else if (s <= 2 + RC_PI * rr) { double a = (s - 2) / rr; p[0] = 2 + rr * sin(a); p[2] = rr - rr * cos(a); }
	else { p[0] = 2 - (s - 2 - RC_PI * rr); p[2] = RC_GAP; }
}

/* a closed loop of n vertices (edges i, i + 1) */
static qaws_repulsive_curves* rc_loops(unsigned int n, unsigned int loops)
{
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * n * loops), l, i;
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	for (l = 0; l < loops; l++)
		for (i = 0; i < n; i++)
		{
			e[2 * (l * n + i)] = l * n + i;
			e[2 * (l * n + i) + 1] = l * n + (i + 1) % n;
		}
	memset(&d, 0, sizeof d);
	d.dimension = 3;
	d.vertex_count = n * loops;
	d.edges = e;
	d.edge_count = n * loops;
	qaws_repulsive_curves_create(&c, &d);
	free(e);
	return c;
}

/* the locations of points, projected onto the mesh's triangles */
static void rc_locate(qaws_mesh_location* out, qaws_mesh_3d const* m, qaws_scalar const* p, unsigned int n)
{
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	memset(&rd, 0, sizeof rd);
	rd.interpolant = QAWS_INTERPOLANT_LINEAR;
	qaws_mesh_reference_create(&r, m, &rd);
	qaws_mesh_reference_project(out, NULL, r, p, NULL, n);
	qaws_mesh_reference_destroy(r);
}

/* an ellipse of n vertices about c, radii a (x), b (y), lifted by z */
static void rc_ellipse(qaws_scalar* p, unsigned int n, double cx, double cy, double z, double a, double b)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * RC_PI * (i + 0.25) / n;
		p[3 * i] = (qaws_scalar)(cx + a * cos(t)); p[3 * i + 1] = (qaws_scalar)(cy + b * sin(t)); p[3 * i + 2] = (qaws_scalar)z;
	}
}

/* the gradient at L against central differences of moves inside the triangles of every step-th vertex: the worst relative error */
static double rc_fd(qaws_repulsive_curves const* c, qaws_mesh_location const* L, unsigned int n, unsigned int step, qaws_repulsive_mesh_curves_desc const* md)
{
	qaws_repulsive_mesh_curves_solver* t = NULL;
	qaws_mesh_location* Lp = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * 2 * n);
	qaws_mesh_location* Lm = Lp + n;
	qaws_scalar* gg = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * n);
	qaws_scalar *xp = gg + 3 * n, *xm = xp + 3 * n, E = 0;
	double worst = 0;
	unsigned int i, k;
	qaws_repulsive_mesh_curves_solver_create(&t, c, L, md, NULL);
	qaws_repulsive_mesh_curves_solver_energy(&E, gg, t);
	for (i = 0; i < n; i += step)
		for (k = 0; k < 2; k++)
		{
			qaws_repulsive_mesh_curves_solver *sp = NULL, *sm = NULL;
			qaws_scalar Ep = 0, Em = 0;
			double h = 1e-4, fd, an = 0;
			unsigned int q;
			memcpy(Lp, L, sizeof(qaws_mesh_location) * n); memcpy(Lm, L, sizeof(qaws_mesh_location) * n);
			if (k == 0) { Lp[i].a = (qaws_scalar)((double)Lp[i].a + h); Lm[i].a = (qaws_scalar)((double)Lm[i].a - h); }
			else { Lp[i].b = (qaws_scalar)((double)Lp[i].b + h); Lm[i].b = (qaws_scalar)((double)Lm[i].b - h); }
			qaws_repulsive_mesh_curves_solver_create(&sp, c, Lp, md, NULL);
			qaws_repulsive_mesh_curves_solver_create(&sm, c, Lm, md, NULL);
			qaws_repulsive_mesh_curves_solver_energy(&Ep, NULL, sp);
			qaws_repulsive_mesh_curves_solver_energy(&Em, NULL, sm);
			qaws_repulsive_mesh_curves_solver_get_positions(xp, sp);
			qaws_repulsive_mesh_curves_solver_get_positions(xm, sm);
			for (q = 0; q < 3; q++) an += (double)gg[3 * i + q] * ((double)xp[3 * i + q] - (double)xm[3 * i + q]);
			fd = (double)Ep - (double)Em;
			worst = fmax(worst, fabs(fd - an) / fmax(fabs(fd), 1e-3 * (double)E * 2 * h));
			qaws_repulsive_mesh_curves_solver_destroy(sp);
			qaws_repulsive_mesh_curves_solver_destroy(sm);
		}
	qaws_repulsive_mesh_curves_solver_destroy(t);
	free(Lp);
	free(gg);
	return worst;
}

static void test_flat(void)
{
	qaws_mesh_3d* m = rc_grid(16, 16, rc_flat);
	qaws_repulsive_curves* c = rc_loops(12, 1);
	qaws_mesh_location L[12];
	qaws_scalar p[36], x[36], g[36], ge[36], E = 0, Ec = 0, Eg = 0, gg[36];
	qaws_repulsive_mesh_curves_desc md;
	qaws_repulsive_mesh_curves_solver *s = NULL, *t = NULL;
	char msg[300];
	double err = 0, gmax = 0;
	unsigned int i;
	rc_ellipse(p, 12, 0.03, -0.02, 0, 0.5, 0.35);
	rc_locate(L, m, p, 12);
	memset(&md, 0, sizeof md);
	md.mesh = m;
	md.euclidean = 1;
	qaws_repulsive_mesh_curves_solver_create(&s, c, L, &md, NULL);
	qaws_repulsive_mesh_curves_solver_get_positions(x, s);
	qaws_repulsive_mesh_curves_solver_energy(&E, g, s);
	memset(ge, 0, sizeof ge);
	qaws_repulsive_curves_gradient(ge, &Ec, c, x, 1);
	for (i = 0; i < 36; i++) { err = fmax(err, fabs((double)g[i] - (double)ge[i])); gmax = fmax(gmax, fabs((double)ge[i])); }
	sprintf(msg, "euclidean: energy %.6g, the curves' %.6g; gradient error %.2g of %.3g", (double)E, (double)Ec, err, gmax);
	printf("  %s\n", msg);
	TEST_ASSERT(s && fabs((double)E - (double)Ec) <= 1e-9 * (double)Ec && err <= 1e-8 * gmax, msg);
	/* geodesic: the same energy up to the heat method's errors */
	md.euclidean = 0;
	qaws_repulsive_mesh_curves_solver_create(&t, c, L, &md, NULL);
	qaws_repulsive_mesh_curves_solver_energy(&Eg, gg, t);
	sprintf(msg, "geodesic on a flat grid: energy %.6g, the chords' %.6g", (double)Eg, (double)Ec);
	printf("  %s\n", msg);
	TEST_ASSERT(t && fabs((double)Eg - (double)Ec) <= 0.1 * (double)Ec, msg);
	qaws_repulsive_mesh_curves_solver_destroy(s);
	qaws_repulsive_mesh_curves_solver_destroy(t);
	qaws_repulsive_curves_destroy(c);
	qaws_mesh_3d_destroy(m);
}

static void test_fold(void)
{
	qaws_mesh_3d* m = rc_grid(56, 12, rc_fold);
	qaws_repulsive_curves* c = rc_loops(10, 2);
	qaws_mesh_location L[20];
	qaws_scalar p[60], Ee = 0, Eg = 0;
	qaws_repulsive_mesh_curves_desc md;
	qaws_repulsive_mesh_curves_solver *s = NULL, *t = NULL;
	char msg[300];
	rc_ellipse(p, 10, 1, 0.5, 0, 0.3, 0.3);
	rc_ellipse(p + 30, 10, 1, 0.5, RC_GAP, 0.3, 0.3);
	rc_locate(L, m, p, 20);
	memset(&md, 0, sizeof md);
	md.mesh = m;
	md.euclidean = 1;
	qaws_repulsive_mesh_curves_solver_create(&s, c, L, &md, NULL);
	qaws_repulsive_mesh_curves_solver_energy(&Ee, NULL, s);
	md.euclidean = 0;
	qaws_repulsive_mesh_curves_solver_create(&t, c, L, &md, NULL);
	qaws_repulsive_mesh_curves_solver_energy(&Eg, NULL, t);
	sprintf(msg, "two loops facing across a fold %.2f apart: euclidean energy %.5g, geodesic %.5g", RC_GAP, (double)Ee, (double)Eg);
	printf("  %s\n", msg);
	TEST_ASSERT(s && t && (double)Eg < 0.2 * (double)Ee, msg);
	{
		double worst = rc_fd(c, L, 20, 3, &md);
		sprintf(msg, "its geodesic gradient against central differences: worst relative error %.3f", worst);
		printf("  %s\n", msg);
		TEST_ASSERT(worst < 0.25, msg);
	}
	qaws_repulsive_mesh_curves_solver_destroy(s);
	qaws_repulsive_mesh_curves_solver_destroy(t);
	qaws_repulsive_curves_destroy(c);
	qaws_mesh_3d_destroy(m);
}

static double rc_aspect(qaws_scalar const* x, unsigned int n)
{
	double lo[2] = { 1e30, 1e30 }, hi[2] = { -1e30, -1e30 };
	unsigned int i, k;
	for (i = 0; i < n; i++)
		for (k = 0; k < 2; k++) { lo[k] = fmin(lo[k], (double)x[3 * i + k]); hi[k] = fmax(hi[k], (double)x[3 * i + k]); }
	return (hi[0] - lo[0]) / (hi[1] - lo[1]);
}

static void test_solve(void)
{
	qaws_mesh_3d* m = rc_grid(16, 16, rc_flat);
	qaws_repulsive_curves* c = rc_loops(16, 1);
	qaws_mesh_location L[16], L0[16];
	qaws_scalar p[48], x[48];
	qaws_repulsive_mesh_curves_desc md;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_mesh_curves_solver* s = NULL;
	qaws_repulsive_report r;
	char msg[400];
	double zmax = 0, a0, a1;
	unsigned int i, crossed = 0;
	rc_ellipse(p, 16, 0.02, 0.01, 0, 0.55, 0.25);
	a0 = rc_aspect(p, 16);
	rc_locate(L, m, p, 16);
	memcpy(L0, L, sizeof L);
	memset(&md, 0, sizeof md);
	md.mesh = m;
	memset(&sd, 0, sizeof sd);
	sd.hold_total_length = 1;
	sd.max_iterations = 30;
	memset(&r, 0, sizeof r);
	qaws_repulsive_mesh_curves_solve(&r, L, c, &md, &sd);
	qaws_repulsive_mesh_curves_solver_create(&s, c, L, &md, NULL);
	qaws_repulsive_mesh_curves_solver_get_positions(x, s);
	for (i = 0; i < 16; i++)
	{
		zmax = fmax(zmax, fabs((double)x[3 * i + 2]));
		crossed += L[i].patch != L0[i].patch;
	}
	a1 = rc_aspect(x, 16);
	sprintf(msg, "an ellipse, its length held: %u steps, energy %.5g to %.5g, aspect %.2f to %.2f, %u of 16 vertices in other triangles, |z| %.1g",
		r.iterations, (double)r.initial_energy, (double)r.energy, a0, a1, crossed, zmax);
	printf("  %s\n", msg);
	TEST_ASSERT(s && r.iterations > 0 && r.energy < r.initial_energy && a1 < 0.75 * a0 && crossed >= 8 && zmax < 1e-12 &&
		(double)r.constraint_residual < 1e-6, msg);
	qaws_repulsive_mesh_curves_solver_destroy(s);
	qaws_repulsive_curves_destroy(c);
	qaws_mesh_3d_destroy(m);
}

static void test_refusals(void)
{
	qaws_mesh_3d* m = rc_grid(4, 4, rc_flat);
	qaws_repulsive_curves* c = rc_loops(4, 1);
	qaws_mesh_location L[4];
	qaws_scalar p[12];
	qaws_repulsive_mesh_curves_desc md;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_mesh_curves_solver* s = NULL;
	qaws_status a, b;
	char msg[200];
	rc_ellipse(p, 4, 0, 0, 0, 0.5, 0.5);
	rc_locate(L, m, p, 4);
	memset(&md, 0, sizeof md);
	md.mesh = m;
	memset(&sd, 0, sizeof sd);
	sd.theta = 0.5;
	a = qaws_repulsive_mesh_curves_solver_create(&s, c, L, &md, &sd);
	L[2].patch = 1000;
	b = qaws_repulsive_mesh_curves_solver_create(&s, c, L, &md, NULL);
	sprintf(msg, "refused: theta set (%d), a location off the triangles (%d)", (int)a, (int)b);
	printf("  %s\n", msg);
	TEST_ASSERT(a == QAWS_STATUS_INVALID_ARGUMENT && b == QAWS_STATUS_INVALID_ARGUMENT && !s, msg);
	qaws_repulsive_curves_destroy(c);
	qaws_mesh_3d_destroy(m);
}

int test_170_repulsive_mesh_curves_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 170: repulsive curves on a mesh\n");
	test_flat();
	test_fold();
	test_solve();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
