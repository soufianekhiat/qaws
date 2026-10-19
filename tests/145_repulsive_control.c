/*
 * Test 145: host control of repulsive runs
 *
 *   - a step cancelled at each of its poll points in turn (dense and
 *     hierarchical): QAWS_STATUS_CANCELLED, the reason
 *     QAWS_REPULSIVE_CANCELLED, the positions untouched bit for bit; the
 *     solver steps on afterwards;
 *   - *_solve reports every accepted step (on_step: as many calls as
 *     iterations, energies not rising, the positions passed being the
 *     caller's); cancelled after a few steps it returns CANCELLED with the
 *     positions of the last reported step;
 *   - the solver's objective: with nothing added, the energy and gradient of
 *     qaws_repulsive_curves_gradient; with a length term, that plus the
 *     weighted length, its gradient against central differences;
 *   - a mesh solve cancelled after two steps.
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RC_N 24

typedef struct rc_state
{
	int polls, cancel_at;                      /* cancel when polls reaches cancel_at (< 0: never) */
	unsigned int steps, cancel_after_steps;    /* on_step count; cancel once it reaches cancel_after_steps (0: never) */
	qaws_scalar energy_prev;
	unsigned int rising, wrong_array;
	qaws_scalar const* positions;              /* the caller's array, to check on_step receives it */
	qaws_scalar last[3 * RC_N + 3 * 64];
	unsigned int count;
} rc_state;

static int rc_cancel(void* u)
{
	rc_state* s = (rc_state*)u;
	if (s->cancel_after_steps && s->steps >= s->cancel_after_steps)
		return 1;
	return s->cancel_at >= 0 && s->polls++ >= s->cancel_at;
}

static void rc_on_step(qaws_repulsive_report const* r, qaws_scalar const* positions, void* u)
{
	rc_state* s = (rc_state*)u;
	s->rising += s->steps && r->energy > s->energy_prev;
	s->energy_prev = r->energy;
	s->wrong_array += positions != s->positions;
	memcpy(s->last, positions, sizeof(qaws_scalar) * s->count);
	s->steps++;
}

/* a wobbly closed loop in the plane z = 0 */
static void rc_loop(qaws_scalar* x, unsigned int* edges)
{
	unsigned int i;
	for (i = 0; i < RC_N; i++)
	{
		double a = 2 * M_PI * i / RC_N, r = 1 + 0.25 * sin(3 * a);
		x[3 * i] = (qaws_scalar)(r * cos(a));
		x[3 * i + 1] = (qaws_scalar)(r * sin(a));
		x[3 * i + 2] = (qaws_scalar)(0.1 * sin(2 * a));
		edges[2 * i] = i;
		edges[2 * i + 1] = (i + 1) % RC_N;
	}
}

static qaws_repulsive_curves* rc_curves(unsigned int* edges)
{
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	memset(&d, 0, sizeof d);
	d.dimension = 3;
	d.vertex_count = RC_N;
	d.edges = edges;
	d.edge_count = RC_N;
	qaws_repulsive_curves_create(&d, &c);
	return c;
}

static void test_cancel_step(void)
{
	qaws_scalar x[3 * RC_N], x0[3 * RC_N];
	unsigned int edges[2 * RC_N], mode, k, bad = 0, wrong = 0, cases = 0;
	qaws_repulsive_curves* c;
	char msg[200];
	rc_loop(x, edges);
	c = rc_curves(edges);
	if (!c)
	{
		TEST_ASSERT(0, "repulsive curves created");
		return;
	}
	for (mode = 0; mode < 2; mode++)
		for (k = 0; k < 6; k++)
		{
			qaws_repulsive_curves_solve_desc d;
			qaws_repulsive_curves_solver* s = NULL;
			qaws_repulsive_report r;
			rc_state st;
			qaws_status res;
			memset(&d, 0, sizeof d);
			memset(&st, 0, sizeof st);
			st.cancel_at = (int)k;
			d.theta = mode ? (qaws_scalar)0.5 : 0;
			d.progress.cancel = rc_cancel;
			d.progress.user_data = &st;
			rc_loop(x, edges);
			if (qaws_repulsive_curves_solver_create(c, &d, x, &s) != QAWS_STATUS_OK)
			{
				wrong++;
				continue;
			}
			memcpy(x0, x, sizeof x);
			res = qaws_repulsive_curves_solver_step(s, x, &r);
			/* k at or past the step's poll count: the step completes */
			if (res == QAWS_STATUS_CANCELLED)
			{
				cases++;
				bad += r.reason != QAWS_REPULSIVE_CANCELLED || memcmp(x, x0, sizeof x) != 0;
			}
			else
				wrong += res != QAWS_STATUS_OK;
			/* the solver steps on once the host stops cancelling */
			st.cancel_at = -1;
			wrong += qaws_repulsive_curves_solver_step(s, x, &r) != QAWS_STATUS_OK || r.reason == QAWS_REPULSIVE_CANCELLED;
			qaws_repulsive_curves_solver_destroy(s);
		}
	sprintf(msg, "steps cancelled at each poll point (dense and hierarchical): %u cancelled, %u moved or misreported, %u other failures", cases,
		bad, wrong);
	TEST_ASSERT(cases >= 8 && bad == 0 && wrong == 0, msg);
	qaws_repulsive_curves_destroy(c);
}

static void test_solve_progress(void)
{
	qaws_scalar x[3 * RC_N];
	unsigned int edges[2 * RC_N];
	qaws_repulsive_curves* c;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_report r;
	rc_state st;
	qaws_status res;
	char msg[256];
	rc_loop(x, edges);
	c = rc_curves(edges);
	if (!c)
		return;
	memset(&d, 0, sizeof d);
	memset(&st, 0, sizeof st);
	st.cancel_at = -1;
	st.positions = x;
	st.count = 3 * RC_N;
	d.max_iterations = 12;
	d.progress.cancel = rc_cancel;
	d.progress.on_step = rc_on_step;
	d.progress.user_data = &st;
	res = qaws_repulsive_curves_solve(c, &d, x, &r);
	sprintf(msg, "solve reports each accepted step: %u calls for %u iterations, %u rising energies, %u other arrays (%s)", st.steps, r.iterations,
		st.rising, st.wrong_array, qaws_status_to_string(res));
	TEST_ASSERT(res == QAWS_STATUS_OK && st.steps == r.iterations && st.steps > 0 && st.rising == 0 && st.wrong_array == 0, msg);
	/* cancelled once three steps are reported */
	rc_loop(x, edges);
	memset(&st, 0, sizeof st);
	st.cancel_at = -1;
	st.cancel_after_steps = 3;
	st.positions = x;
	st.count = 3 * RC_N;
	res = qaws_repulsive_curves_solve(c, &d, x, &r);
	sprintf(msg, "solve cancelled after 3 steps: %s, reason %d, %u steps reported, positions those of the last reported step: %d",
		qaws_status_to_string(res), (int)r.reason, st.steps, memcmp(st.last, x, sizeof x) == 0);
	TEST_ASSERT(res == QAWS_STATUS_CANCELLED && r.reason == QAWS_REPULSIVE_CANCELLED && st.steps == 3 && memcmp(st.last, x, sizeof x) == 0, msg);
	qaws_repulsive_curves_destroy(c);
}

static void test_objective(void)
{
	qaws_scalar x[3 * RC_N], g[3 * RC_N], gref[3 * RC_N], E = 0, Eref = 0, El = 0;
	unsigned int edges[2 * RC_N], i;
	qaws_repulsive_curves* c;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_curves_solver* s = NULL;
	double worst = 0, fd = 0, len = 0, h = 1e-6;
	char msg[256];
	rc_loop(x, edges);
	c = rc_curves(edges);
	if (!c)
		return;
	memset(&d, 0, sizeof d);
	memset(gref, 0, sizeof gref);
	qaws_repulsive_curves_gradient(c, x, 1, gref, &Eref);
	if (qaws_repulsive_curves_solver_create(c, &d, x, &s) == QAWS_STATUS_OK && qaws_repulsive_curves_solver_objective(s, x, &E, g) == QAWS_STATUS_OK)
	{
		for (i = 0; i < 3 * RC_N; i++)
			worst = fabs((double)g[i] - gref[i]) > worst ? fabs((double)g[i] - gref[i]) : worst;
		sprintf(msg, "objective without added terms: energy %.12g against %.12g, gradient within %.2g", (double)E, (double)Eref, worst);
		TEST_ASSERT(fabs((double)E - Eref) <= 1e-9 * fabs((double)Eref) && worst <= 1e-9 * (1 + fabs((double)Eref)), msg);
	}
	else
		TEST_ASSERT(0, "objective evaluated");
	qaws_repulsive_curves_solver_destroy(s);
	s = NULL;
	/* a length term: E + w L, its gradient against central differences */
	d.length_weight = (qaws_scalar)0.7;
	for (i = 0; i < RC_N; i++)
	{
		double dx = x[3 * ((i + 1) % RC_N)] - x[3 * i], dy = x[3 * ((i + 1) % RC_N) + 1] - x[3 * i + 1], dz = x[3 * ((i + 1) % RC_N) + 2] - x[3 * i + 2];
		len += sqrt(dx * dx + dy * dy + dz * dz);
	}
	if (qaws_repulsive_curves_solver_create(c, &d, x, &s) == QAWS_STATUS_OK && qaws_repulsive_curves_solver_objective(s, x, &El, g) == QAWS_STATUS_OK)
	{
		for (i = 0; i < 3 * RC_N; i += 5)
		{
			qaws_scalar xp[3 * RC_N], Ep = 0, Em = 0;
			memcpy(xp, x, sizeof x);
			xp[i] = (qaws_scalar)(x[i] + h);
			qaws_repulsive_curves_solver_objective(s, xp, &Ep, NULL);
			xp[i] = (qaws_scalar)(x[i] - h);
			qaws_repulsive_curves_solver_objective(s, xp, &Em, NULL);
			fd = fabs(((double)Ep - Em) / (2 * h) - g[i]) > fd ? fabs(((double)Ep - Em) / (2 * h) - g[i]) : fd;
		}
		sprintf(msg, "objective with a length term: %.12g against E + w L = %.12g; gradient against differences %.2g", (double)El, Eref + 0.7 * len, fd);
		TEST_ASSERT(fabs((double)El - (Eref + 0.7 * len)) <= 1e-9 * fabs((double)El) && fd <= 1e-4 * (1 + fabs((double)El)), msg);
		printf("  %s\n", msg);
	}
	else
		TEST_ASSERT(0, "objective with a length term evaluated");
	qaws_repulsive_curves_solver_destroy(s);
	qaws_repulsive_curves_destroy(c);
}

/* a subdivided octahedron on the unit sphere, squashed */
static void test_mesh_cancel(void)
{
	static qaws_scalar const oct[6 * 3] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1 };
	static unsigned int const tris[8 * 3] = { 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5 };
	qaws_scalar x[6 * 3];
	qaws_repulsive_mesh_desc md;
	qaws_repulsive_mesh* m = NULL;
	qaws_repulsive_mesh_solve_desc d;
	qaws_repulsive_report r;
	rc_state st;
	qaws_status res;
	unsigned int i;
	char msg[256];
	for (i = 0; i < 18; i++)
		x[i] = oct[i] * (qaws_scalar)(i % 3 == 2 ? 0.6 : 1.0);
	memset(&md, 0, sizeof md);
	md.vertex_count = 6;
	md.triangles = tris;
	md.triangle_count = 8;
	if (qaws_repulsive_mesh_create(&md, &m) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "repulsive mesh created");
		return;
	}
	memset(&d, 0, sizeof d);
	memset(&st, 0, sizeof st);
	st.cancel_at = -1;
	st.cancel_after_steps = 2;
	st.positions = x;
	st.count = 18;
	d.max_iterations = 20;
	d.hold_area = 1;
	d.progress.cancel = rc_cancel;
	d.progress.on_step = rc_on_step;
	d.progress.user_data = &st;
	res = qaws_repulsive_mesh_solve(m, &d, x, &r);
	sprintf(msg, "mesh solve cancelled after 2 steps: %s, %u steps reported, positions those of the last: %d", qaws_status_to_string(res), st.steps,
		memcmp(st.last, x, sizeof x) == 0);
	TEST_ASSERT(res == QAWS_STATUS_CANCELLED && r.reason == QAWS_REPULSIVE_CANCELLED && st.steps == 2 && memcmp(st.last, x, sizeof x) == 0, msg);
	qaws_repulsive_mesh_destroy(m);
}

int test_145_repulsive_control_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 145: host control of repulsive runs\n");
	test_cancel_step();
	test_solve_progress();
	test_objective();
	test_mesh_cancel();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
