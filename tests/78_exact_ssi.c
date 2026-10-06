/*
 * Test 78: Certified surface / surface intersection curves
 *
 *   - a plane and a dome (one closed loop), a plane and a shifted saddle
 *     (two open branches), two elliptic tubes, rational NURBS conics
 *     extruded, the thin one piercing the wide one (two closed loops,
 *     chained across patch boundaries)
 *   - every certified point lies on both surfaces; min_depth gives more
 *     points with the same topology
 *   - a plane tangent to the dome's top is refused
 */

#include "test_common.h"
#include "qaws_exact.h"
#include <math.h>
#include <string.h>

static qaws_exact_surface* ssi_patch(double const* z, unsigned int deg, double x0, double x1)
{
	qaws_vec3 cps[16];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	qaws_exact_surface* e = NULL;
	unsigned int i, j;
	for (i = 0; i <= deg; i++)
		for (j = 0; j <= deg; j++)
		{
			cps[i * (deg + 1) + j].x = (qaws_scalar)(x0 + (x1 - x0) * i / deg);
			cps[i * (deg + 1) + j].y = (qaws_scalar)(x0 + (x1 - x0) * j / deg);
			cps[i * (deg + 1) + j].z = (qaws_scalar)z[i * (deg + 1) + j];
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = deg;
	d.v_degree = deg;
	d.control_points = cps;
	d.u_point_count = deg + 1;
	d.v_point_count = deg + 1;
	qaws_surface_create_bezier(&d, &s);
	qaws_exact_surface_prepare(NULL, s, &e, NULL);
	qaws_surface_destroy(s);
	return e;
}

/* An elliptic tube: the 3-arc NURBS conic (scaled by k) extruded; axis 2: along z, axis 0: along x. */
static qaws_exact_surface* ssi_tube(double k, int along_x, double h)
{
	static double const ex[7] = { 2, 2, -1, -4, -1, 2, 2 }, ey[7] = { 0, 4, 2, 0, -2, -4, 0 }, ew[7] = { 2, 1, 2, 1, 2, 1, 2 };
	qaws_vec3 cps[14];
	qaws_scalar ws[14], ku[10] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 }, kv[4] = { 0, 0, 1, 1 };
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	qaws_exact_surface* e = NULL;
	unsigned int i, j;
	for (i = 0; i < 7; i++)
		for (j = 0; j < 2; j++)
		{
			double a = k * ex[i], b = k * ey[i], t = j ? h : -h;
			qaws_vec3* p = &cps[i * 2 + j];
			if (along_x)
			{
				p->x = (qaws_scalar)t;
				p->y = (qaws_scalar)a;
				p->z = (qaws_scalar)b;
			}
			else
			{
				p->x = (qaws_scalar)a;
				p->y = (qaws_scalar)b;
				p->z = (qaws_scalar)t;
			}
			ws[i * 2 + j] = (qaws_scalar)ew[i];
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 2;
	d.v_degree = 1;
	d.control_points = cps;
	d.u_point_count = 7;
	d.v_point_count = 2;
	d.weights = ws;
	d.u_knots = ku;
	d.u_knot_count = 10;
	d.v_knots = kv;
	d.v_knot_count = 4;
	qaws_surface_create_nurbs(&d, &s);
	qaws_exact_surface_prepare(NULL, s, &e, NULL);
	qaws_surface_destroy(s);
	return e;
}

static qaws_exact_ssi_point g_pts[8192];

/* Runs the intersection; checks every point on both surfaces; returns branch counts. */
static qaws_status ssi_run(char const* name, qaws_exact_surface* a, qaws_exact_surface* b, unsigned int min_depth, unsigned int* open, unsigned int* closed,
	unsigned int* npoints, double* gap)
{
	qaws_exact_ssi_branch br[64];
	unsigned int np = 0, nb = 0, i;
	qaws_status st = qaws_exact_surface_surface_hits(a, b, min_depth, g_pts, 8192, &np, br, 64, &nb);
	*open = *closed = 0;
	*gap = 0;
	for (i = 0; i < nb; i++)
	{
		if (br[i].closed) (*closed)++;
		else (*open)++;
	}
	for (i = 0; i < np; i++)
	{
		double pa[3], pb[3], g;
		qaws_exact_surface_evaluate(a, 0.5 * (g_pts[i].u1_lo + g_pts[i].u1_hi), 0.5 * (g_pts[i].v1_lo + g_pts[i].v1_hi), 0, pa, NULL);
		qaws_exact_surface_evaluate(b, 0.5 * (g_pts[i].u2_lo + g_pts[i].u2_hi), 0.5 * (g_pts[i].v2_lo + g_pts[i].v2_hi), 0, pb, NULL);
		g = fabs(pa[0] - pb[0]) + fabs(pa[1] - pb[1]) + fabs(pa[2] - pb[2]);
		if (g > *gap) *gap = g;
	}
	*npoints = np;
	printf("    %-28s status %d: %u points, %u open, %u closed branches, worst gap %.1e\n", name, (int)st, np, *open, *closed, *gap);
	return st;
}

static void test_ssi(void)
{
	double plane[4] = { 0, 0, 0, 0 }, top[4] = { 0.5, 0.5, 0.5, 0.5 };
	double dome[9] = { -1, 0.5, -1, 0.5, 2, 0.5, -1, 0.5, -1 }, saddle[9];
	unsigned int i, j, op, cl, np, np2;
	double gap;
	qaws_status st;
	qaws_exact_surface* P = ssi_patch(plane, 1, -1, 2);
	qaws_exact_surface* T = ssi_patch(top, 1, -1, 2);
	qaws_exact_surface* D = ssi_patch(dome, 2, 0, 1);
	qaws_exact_surface* S;
	qaws_exact_surface* A = ssi_tube(1, 0, 4);
	qaws_exact_surface* B = ssi_tube(0.5, 1, 6);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			saddle[i * 3 + j] = ((double)i - 1) * ((double)j - 1) - 0.125;
	S = ssi_patch(saddle, 2, 0, 1);
	st = ssi_run("plane x dome", P, D, 0, &op, &cl, &np, &gap);
	TEST_ASSERT(st == QAWS_STATUS_OK && op == 0 && cl == 1 && gap < 1e-6, "plane x dome: one closed loop, every point on both surfaces");
	st = ssi_run("plane x dome, min_depth 4", P, D, 4, &op, &cl, &np2, &gap);
	TEST_ASSERT(st == QAWS_STATUS_OK && op == 0 && cl == 1 && np2 > 4 * np && gap < 1e-6, "min_depth 4: the same loop, many more certified points");
	st = ssi_run("plane x saddle", P, S, 0, &op, &cl, &np, &gap);
	TEST_ASSERT(st == QAWS_STATUS_OK && op == 2 && cl == 0 && gap < 1e-6, "plane x shifted saddle: two open branches");
	st = ssi_run("tube x tube", A, B, 0, &op, &cl, &np, &gap);
	TEST_ASSERT(st == QAWS_STATUS_OK && op == 0 && cl == 2 && gap < 1e-6,
		"elliptic tubes (rational NURBS, 3 x 3 patch pairs): two closed loops chained across patches");
	st = ssi_run("tube x tube, min_depth 3", A, B, 3, &op, &cl, &np2, &gap);
	TEST_ASSERT(st == QAWS_STATUS_OK && op == 0 && cl == 2 && np2 > 4 * np && gap < 1e-6,
		"min_depth 3: the same two loops (they cross the seams and the symmetric mid-planes), many more points");
	st = ssi_run("plane tangent to the dome", T, D, 0, &op, &cl, &np, &gap);
	TEST_ASSERT(st == QAWS_STATUS_CERTIFICATION_FAILED, "a tangential contact is refused, not guessed");
	qaws_exact_surface_destroy(P);
	qaws_exact_surface_destroy(T);
	qaws_exact_surface_destroy(D);
	qaws_exact_surface_destroy(S);
	qaws_exact_surface_destroy(A);
	qaws_exact_surface_destroy(B);
}

int test_78_exact_ssi_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 78: Certified surface / surface intersection curves\n");
	test_ssi();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
