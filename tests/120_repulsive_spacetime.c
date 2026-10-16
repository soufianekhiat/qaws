/*
 * Test 120: space-time paths and rooms (Repulsive Curves, Sec. 8.5)
 *
 * A room's walls repel each vertex by w times the sum over walls of the
 * integral of rho^-2 along the wall (half the integral of r^-3 over the
 * walls extruded along time: the paper's mesh obstacle):
 *   - its value and gradient against Mathematica's integrals from the
 *     definition (reference/120_repulsive_room.wls), on the paper's H room,
 *     at points that take each branch (a point 1e-6 off a wall's line,
 *     beside the wall, takes the series);
 *   - the 2D integral against half the 3D one;
 *   - refusals: a point outside, a point on a wall; nothing added.
 * The solver: four agents crossing a room diagonally, their paths curves in
 * (x, t, z), time held, the ends pinned, the length weighted, the room's
 * walls repelling: every time held, every vertex inside after every step,
 * the agents' smallest distance in a time slice grown; refusals.
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include "reference/120_repulsive_room.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RT_F32 (sizeof(qaws_scalar) == 4)

static void test_room_term(void)
{
	unsigned int i, q;
	double emax = 0, gmax = 0, gscale = 0;
	char msg[300];
	for (i = 0; i < 5; i++)
	{
		double E = 0, g[2] = { 0, 0 };
		int ok = qaws_internal_rs_room(g120_room, 12, g120_w, 2, 0, 1, g120_pts + 2 * i, 1, &E, g);
		if (!ok) { emax = HUGE_VAL; break; }
		emax = fmax(emax, fabs(E - g120_e[3 * i]) / fabs(g120_e[3 * i]));
		for (q = 0; q < 2; q++)
		{
			gmax = fmax(gmax, fabs(g[q] - g120_e[3 * i + 1 + q]));
			gscale = fmax(gscale, fabs(g120_e[3 * i + 1 + q]));
		}
	}
	sprintf(msg, "room term, the H room at 5 points (one 1e-6 off a wall's line): against Mathematica, energy %.2g, gradient %.2g (of %.3g)", emax, gmax, gscale);
	TEST_ASSERT(emax <= 1e-13 && gmax <= 1e-11 * gscale, msg);
	/* the 2D integral against half the 3D one, one wall */
	{
		double sq[8] = { 1, -1, 1, 1, -1, 1, -1, -1 }, p[2] = { 0, 0 }, E = 0;   /* (a square: its 4 walls alike at the origin) */
		qaws_internal_rs_room(sq, 4, 1, 2, 0, 1, p, 1, &E, NULL);
		sprintf(msg, "room term against the paper's obstacle: a square's wall at the origin %.17g, half the 3D integral of r^-3 %.17g (Mathematica)",
			E / 4, g120_half3d);
		TEST_ASSERT(fabs(E / 4 - g120_wall2d) <= 1e-15 && fabs(g120_wall2d - g120_half3d) <= 1e-15, msg);
	}
	/* refusals: outside, on a wall; nothing added */
	{
		double out[2] = { 0, 2 }, on[2] = { 2, -3 }, E = 5, g[2] = { 7, 7 };
		int a = qaws_internal_rs_room(g120_room, 12, 1, 2, 0, 1, out, 1, &E, g), b = qaws_internal_rs_room(g120_room, 12, 1, 2, 0, 1, on, 1, &E, g);
		sprintf(msg, "room term refusals: a point outside (in the H's notch) %d, on a wall %d; nothing added %d", a, b, E == 5 && g[0] == 7 && g[1] == 7);
		TEST_ASSERT(!a && !b && E == 5 && g[0] == 7 && g[1] == 7, msg);
	}
}

/* four agents, each a path of NV vertices from (x0, z0) at t 0 to (x1, z1)
   at t 8, through the room's middle, nearly meeting there */
enum { RT_A = 4, RT_NV = 65 };

static void rt_paths(qaws_scalar* x, unsigned int* e)
{
	static double const corner[4][2] = { { -2, -2 }, { 2, -2 }, { 2, 2 }, { -2, 2 } };
	unsigned int a, i;
	for (a = 0; a < RT_A; a++)
		for (i = 0; i < RT_NV; i++)
		{
			double s = (double)i / (RT_NV - 1), b = 4 * s * (1 - s), th = 1.5707963267948966 * a + 0.3;
			double const *p0 = corner[a], *p1 = corner[(a + 2) % 4];
			unsigned int v = a * RT_NV + i;
			x[3 * v] = (qaws_scalar)(p0[0] + s * (p1[0] - p0[0]) + 0.15 * b * cos(th));
			x[3 * v + 1] = (qaws_scalar)(8 * s);
			x[3 * v + 2] = (qaws_scalar)(p0[1] + s * (p1[1] - p0[1]) + 0.15 * b * sin(th));
			if (i + 1 < RT_NV) { e[2 * (a * (RT_NV - 1) + i)] = v; e[2 * (a * (RT_NV - 1) + i) + 1] = v + 1; }
		}
}

/* the agents' smallest distance in a time slice (vertex i of each path shares its time) */
static double rt_closest(qaws_scalar const* x)
{
	double m = HUGE_VAL;
	unsigned int a, b, i;
	for (i = 0; i < RT_NV; i++)
		for (a = 0; a < RT_A; a++)
			for (b = a + 1; b < RT_A; b++)
			{
				unsigned int u = a * RT_NV + i, v = b * RT_NV + i;
				m = fmin(m, hypot((double)x[3 * u] - x[3 * v], (double)x[3 * u + 2] - x[3 * v + 2]));
			}
	return m;
}

static void test_agents(void)
{
	static qaws_scalar const room[8] = { -3, -3, 3, -3, 3, 3, -3, 3 };
	qaws_scalar x[3 * RT_A * RT_NV], t0[RT_A * RT_NV];
	unsigned int e[2 * RT_A * (RT_NV - 1)], pins[2 * RT_A], i, a, it, inside = 1;
	qaws_repulsive_curves_desc cd;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_curves_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_status st;
	double close0, close1, terr = 0;
	char msg[400];
	rt_paths(x, e);
	for (i = 0; i < RT_A * RT_NV; i++) t0[i] = x[3 * i + 1];
	for (a = 0; a < RT_A; a++) { pins[2 * a] = a * RT_NV; pins[2 * a + 1] = a * RT_NV + RT_NV - 1; }
	close0 = rt_closest(x);
	memset(&cd, 0, sizeof(cd));
	cd.dimension = 3; cd.vertex_count = RT_A * RT_NV; cd.edges = e; cd.edge_count = RT_A * (RT_NV - 1); cd.alpha = 3; cd.beta = 6;
	qaws_repulsive_curves_create(&cd, &c);
	memset(&d, 0, sizeof(d));
	d.pinned = pins; d.pinned_count = 2 * RT_A;
	d.hold_coordinate = 2; d.length_weight = 1; d.room = room; d.room_count = 4; d.room_weight = (qaws_scalar)0.1;
	d.tolerance = (qaws_scalar)(RT_F32 ? 1e-3 : 1e-4);
	st = qaws_repulsive_curves_solver_create(c, &d, x, &s);
	memset(&r, 0, sizeof(r));
	for (it = 0; st == QAWS_STATUS_OK && it < 300; it++)
	{
		qaws_repulsive_curves_solver_step(s, x, &r);
		for (i = 0; i < RT_A * RT_NV; i++)
			inside &= fabs((double)x[3 * i]) < 3 && fabs((double)x[3 * i + 2]) < 3;
		if (r.reason != QAWS_REPULSIVE_RUNNING)
			break;
	}
	for (i = 0; i < RT_A * RT_NV; i++) terr = fmax(terr, fabs((double)x[3 * i + 1] - t0[i]));
	close1 = rt_closest(x);
	sprintf(msg, "four agents crossing a room in space-time (status %d): reason %d after %u steps, E %.5g from %.5g; times held %.2g, inside %d; the closest agents in a time slice %.3f -> %.3f",
		(int)st, (int)r.reason, r.iterations, (double)r.energy, (double)r.initial_energy, terr, inside, close0, close1);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.reason == QAWS_REPULSIVE_CONVERGED && r.energy < r.initial_energy && terr <= (RT_F32 ? 1e-4 : 1e-7) &&
		inside && close0 < 0.3 && close1 > 1.5, msg);
	/* no solution derivatives with these terms */
	{
		qaws_repulsive_targets t;
		qaws_scalar dx[3 * RT_A * RT_NV];
		memset(&t, 0, sizeof(t));
		st = s ? qaws_repulsive_curves_solver_solution_jvp(s, &t, dx, NULL) : QAWS_STATUS_INVALID_ARGUMENT;
		sprintf(msg, "space-time solver: solution derivatives refused (%d)", (int)st);
		TEST_ASSERT(st == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	}
	qaws_repulsive_curves_solver_destroy(s);
	/* refusals */
	{
		qaws_repulsive_curves_solve_desc b = d;
		qaws_scalar y[3 * RT_A * RT_NV];
		qaws_status s1, s2, s3, s4;
		qaws_repulsive_curves_solver* t = NULL;
		rt_paths(y, e);
		b.hold_coordinate = 4; s1 = qaws_repulsive_curves_solver_create(c, &b, y, &t); qaws_repulsive_curves_solver_destroy(t); t = NULL;
		b = d; b.room_count = 2; s2 = qaws_repulsive_curves_solver_create(c, &b, y, &t); qaws_repulsive_curves_solver_destroy(t); t = NULL;
		b = d; b.length_weight = -1; s3 = qaws_repulsive_curves_solver_create(c, &b, y, &t); qaws_repulsive_curves_solver_destroy(t); t = NULL;
		y[3 * 5] = 4;   /* (a vertex out of the room) */
		b = d; s4 = qaws_repulsive_curves_solver_create(c, &b, y, &t); qaws_repulsive_curves_solver_destroy(t);
		sprintf(msg, "space-time refusals: coordinate 4 of 3 (%d), a room of 2 points (%d), a negative length weight (%d), a vertex outside (%d)",
			(int)s1, (int)s2, (int)s3, (int)s4);
		TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT &&
			s4 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	}
	qaws_repulsive_curves_destroy(c);
}

int test_120_repulsive_spacetime_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 120: space-time paths and rooms\n");
	test_room_term();
	test_agents();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
