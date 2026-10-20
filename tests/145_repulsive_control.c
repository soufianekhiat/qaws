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
 *   - a mesh solve cancelled after two steps;
 *   - the line search and GMRES settings: their defaults given explicitly
 *     reproduce zeros bit for bit (dense and hierarchical); from a far first
 *     trial a faster backtrack refuses fewer trials, max_trials 2 stops the
 *     step after two; out-of-range values refused at creation;
 *   - _ex allocators: every block of the objects and of the calls on them
 *     from the given allocator, none from the thread allocator, all returned
 *     after destroy; with a scratch arena no call leaves a block in it.
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "qaws_alloc.h"
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
	/* the difference step: 1e-6 in double. float: fourth-order differences with h = 2^-8 (exact in float); measured, their
	   error falls as h^4 down to 3e-4 there and rises below it with rounding, within the same 1e-4 (1 + E) bound */
	double worst = 0, fd = 0, len = 0, h = sizeof(qaws_scalar) == sizeof(double) ? 1e-6 : 1.0 / 256;
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
			if (sizeof(qaws_scalar) == sizeof(double))
			{
				/* (divided by the step the scalar coordinates actually took) */
				double step = (double)(qaws_scalar)(x[i] + h) - (double)(qaws_scalar)(x[i] - h);
				fd = fabs(((double)Ep - Em) / step - g[i]) > fd ? fabs(((double)Ep - Em) / step - g[i]) : fd;
			}
			else
			{
				/* float: fourth-order differences (the energy's third derivative is large, so second order needs a step
				   too small for float's rounding) */
				qaws_scalar Ep2 = 0, Em2 = 0;
				double e;
				xp[i] = (qaws_scalar)(x[i] + 2 * h);
				qaws_repulsive_curves_solver_objective(s, xp, &Ep2, NULL);
				xp[i] = (qaws_scalar)(x[i] - 2 * h);
				qaws_repulsive_curves_solver_objective(s, xp, &Em2, NULL);
				e = fabs((8 * ((double)Ep - Em) - ((double)Ep2 - Em2)) / (12 * h) - g[i]);
				fd = e > fd ? e : fd;
			}
		}
		sprintf(msg, "objective with a length term: %.12g against E + w L = %.12g; gradient against differences %.2g", (double)El, Eref + 0.7 * len, fd);
		/* (the energy is a sum over the RC_N^2 pairs and RC_N edge lengths: recursive summation errs by up to n eps relative) */
		TEST_ASSERT(fabs((double)El - (Eref + 0.7 * len)) <= fmax(1e-9, RC_N * TEST_EPS) * fabs((double)El) && fd <= 1e-4 * (1 + fabs((double)El)), msg);
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

static void rc_octahedron(qaws_scalar* x, qaws_repulsive_mesh_desc* md, qaws_allocator const* allocator, qaws_repulsive_mesh** m)
{
	static qaws_scalar const oct[6 * 3] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1 };
	static unsigned int const tris[8 * 3] = { 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5 };
	unsigned int i;
	for (i = 0; i < 18; i++)
		x[i] = oct[i] * (qaws_scalar)(i % 3 == 2 ? 0.6 : 1.0);
	memset(md, 0, sizeof *md);
	md->vertex_count = 6;
	md->triangles = tris;
	md->triangle_count = 8;
	*m = NULL;
	qaws_repulsive_mesh_create_ex(md, allocator, m);
}

/* the line search and GMRES settings: explicit defaults reproduce zeros bit for bit; fewer trials, a faster backtrack;
   values out of range refused */
static void test_settings(void)
{
	qaws_scalar xa[3 * RC_N], xb[3 * RC_N], m0[18], m1[18];
	unsigned int edges[2 * RC_N], mode, k, refused = 0, rej[3] = { 0, 0, 0 };
	qaws_repulsive_curves* c;
	qaws_repulsive_mesh_desc md;
	qaws_repulsive_mesh* m = NULL;
	qaws_repulsive_report ra, rb, rm[3];
	qaws_status sa, sb, sm[3] = { QAWS_STATUS_INVALID_ARGUMENT, QAWS_STATUS_INVALID_ARGUMENT, QAWS_STATUS_INVALID_ARGUMENT };
	char msg[300];
	rc_loop(xa, edges);
	c = rc_curves(edges);
	if (!c)
		return;
	for (mode = 0; mode < 2; mode++)
	{
		qaws_repulsive_curves_solve_desc d;
		memset(&d, 0, sizeof d);
		d.max_iterations = 8;
		d.theta = mode ? (qaws_scalar)0.5 : 0;
		rc_loop(xa, edges);
		rc_loop(xb, edges);
		sa = qaws_repulsive_curves_solve(c, &d, xa, &ra);
		d.armijo = (qaws_scalar)0.1; d.backtrack = (qaws_scalar)0.5; d.max_trials = 40;
		d.gmres_max_iterations = 1000; d.gmres_restart = 40; d.gmres_tolerance = (qaws_scalar)1e-8;
		sb = qaws_repulsive_curves_solve(c, &d, xb, &rb);
		sprintf(msg, "%s curves, the settings' defaults given explicitly: %u and %u steps, positions bit for bit %d, reports equal %d",
			mode ? "hierarchical" : "dense", ra.iterations, rb.iterations, memcmp(xa, xb, sizeof xa) == 0, memcmp(&ra, &rb, sizeof ra) == 0);
		TEST_ASSERT(sa == QAWS_STATUS_OK && sb == QAWS_STATUS_OK && ra.iterations > 0 && memcmp(xa, xb, sizeof xa) == 0 && memcmp(&ra, &rb, sizeof ra) == 0, msg);
	}
	/* the hierarchical direction's scalar solves within the limits: one iteration cannot solve them (NUMERICAL_FAILURE);
	   a looser tolerance gives another step */
	{
		qaws_repulsive_curves_solve_desc d;
		qaws_repulsive_curves_solver* s = NULL;
		qaws_repulsive_report r1, r2;
		qaws_status s1 = QAWS_STATUS_INVALID_ARGUMENT, s2 = QAWS_STATUS_INVALID_ARGUMENT;
		memset(&d, 0, sizeof d);
		d.theta = (qaws_scalar)0.5;
		d.gmres_max_iterations = 1;
		rc_loop(xa, edges);
		memset(&r1, 0, sizeof r1); memset(&r2, 0, sizeof r2);
		if (qaws_repulsive_curves_solver_create(c, &d, xa, &s) == QAWS_STATUS_OK)
			s1 = qaws_repulsive_curves_solver_step(s, xa, &r1);
		qaws_repulsive_curves_solver_destroy(s);
		s = NULL;
		d.max_iterations = 8;
		d.gmres_max_iterations = 0;
		d.gmres_tolerance = (qaws_scalar)1e-3;
		rc_loop(xa, edges);
		s2 = qaws_repulsive_curves_solve(c, &d, xa, &r2);
		sprintf(msg, "hierarchical, gmres_max_iterations 1: step status %d, reason %d (NUMERICAL_FAILURE); gmres_tolerance 1e-3: %u steps, "
			"energy %.10g (default %.10g), positions differ %d", (int)s1, (int)r1.reason, r2.iterations, (double)r2.energy, (double)ra.energy,
			memcmp(xa, xb, sizeof xa) != 0);
		TEST_ASSERT(s1 == QAWS_STATUS_OK && r1.reason == QAWS_REPULSIVE_NUMERICAL_FAILURE && r1.iterations == 0 && s2 == QAWS_STATUS_OK &&
			r2.iterations > 0 && r2.energy < r2.initial_energy && memcmp(xa, xb, sizeof xa) != 0, msg);
		printf("  %s\n", msg);
	}
	/* one mesh step from a first trial 40 mean edge lengths long: refused trials until one is accepted */
	rc_octahedron(m0, &md, NULL, &m);
	for (k = 0; k < 3 && m; k++)
	{
		qaws_repulsive_mesh_solve_desc d;
		qaws_repulsive_mesh_solver* s = NULL;
		memset(&d, 0, sizeof d);
		d.hold_area = 1;
		d.initial_step = 40;
		if (k == 1) d.backtrack = (qaws_scalar)0.1;
		if (k == 2) d.max_trials = 2;
		memcpy(m1, m0, sizeof m0);
		memset(&rm[k], 0, sizeof rm[k]);
		sm[k] = qaws_repulsive_mesh_solver_create(m, &d, m1, &s);
		if (sm[k] == QAWS_STATUS_OK)
			sm[k] = qaws_repulsive_mesh_solver_step(s, m1, &rm[k]);
		rej[k] = rm[k].rejected_energy + rm[k].rejected_collision + rm[k].rejected_projection;
		if (k == 2) refused = memcmp(m1, m0, sizeof m0) == 0;
		qaws_repulsive_mesh_solver_destroy(s);
	}
	sprintf(msg, "a mesh step from a 40-edge first trial: backtrack 0.5 %u refused trials then accepted (%u), 0.1 %u then accepted (%u); "
		"max_trials 2: %u refused, reason %d, positions untouched %d", rej[0], rm[0].iterations, rej[1], rm[1].iterations, rej[2], (int)rm[2].reason, (int)refused);
	TEST_ASSERT(m && sm[0] == QAWS_STATUS_OK && sm[1] == QAWS_STATUS_OK && sm[2] == QAWS_STATUS_OK && rm[0].iterations == 1 && rm[1].iterations == 1 &&
		rej[0] > 2 && rej[1] < rej[0] && rm[2].iterations == 0 && rej[2] == 2 && refused &&
		(rm[2].reason == QAWS_REPULSIVE_STALLED || rm[2].reason == QAWS_REPULSIVE_INFEASIBLE || rm[2].reason == QAWS_REPULSIVE_COLLISION_LIMITED), msg);
	printf("  %s\n", msg);
	/* refusals at creation */
	refused = 0;
	for (k = 0; k < 6; k++)
	{
		qaws_repulsive_curves_solve_desc d;
		qaws_repulsive_mesh_solve_desc dm;
		qaws_repulsive_curves_solver* s = NULL;
		qaws_repulsive_mesh_solver* sv = NULL;
		qaws_scalar const bad[6] = { (qaws_scalar)0.5, (qaws_scalar)-0.1, 1, (qaws_scalar)-0.5, -1, (qaws_scalar)0.7 };
		memset(&d, 0, sizeof d);
		memset(&dm, 0, sizeof dm);
		if (k == 0 || k == 1 || k == 5) d.armijo = dm.armijo = bad[k];
		if (k == 2 || k == 3) d.backtrack = dm.backtrack = bad[k];
		if (k == 4) d.gmres_tolerance = dm.gmres_tolerance = bad[k];
		rc_loop(xa, edges);
		refused += qaws_repulsive_curves_solver_create(c, &d, xa, &s) == QAWS_STATUS_INVALID_ARGUMENT && !s;
		refused += m && qaws_repulsive_mesh_solver_create(m, &dm, m0, &sv) == QAWS_STATUS_INVALID_ARGUMENT && !sv;
		qaws_repulsive_curves_solver_destroy(s);
		qaws_repulsive_mesh_solver_destroy(sv);
	}
	sprintf(msg, "armijo 0.5, -0.1, 0.7, backtrack 1, -0.5, gmres_tolerance -1 refused by the curve and mesh solvers: %u of 12", refused);
	TEST_ASSERT(refused == 12, msg);
	qaws_repulsive_mesh_destroy(m);
	qaws_repulsive_curves_destroy(c);
}

/* a counting allocator (legacy interface) */
typedef struct rc_count
{
	unsigned long allocs;
	long live;
} rc_count;

static void* rc_count_alloc(unsigned long size, void* user)
{
	rc_count* c = (rc_count*)user;
	c->allocs++;
	c->live++;
	return malloc(size);
}

static void rc_count_free(void* p, void* user)
{
	rc_count* c = (rc_count*)user;
	if (p)
		c->live--;
	free(p);
}

static unsigned char g_rc_arena[1 << 22];

/* the objects' allocator: every block of objects made with _ex and of calls on them comes from it, none from the thread
   allocator (a second counter, installed as such); everything returned after destroy. Then with a scratch arena: no
   call leaves a block in it. */
static void test_allocator(void)
{
	qaws_scalar x[3 * RC_N], g[3 * RC_N], dir[3 * RC_N], mx[18], dx[18], xbar[18], E = 0;
	unsigned int edges[2 * RC_N], i, it, pass, held = 0, calls = 0;
	rc_count oc, tc;
	qaws_allocator oa, ta, prev, prev_scratch;
	qaws_arena arena;
	qaws_status bad = QAWS_STATUS_OK, st;
	char msg[300];
	memset(&oc, 0, sizeof oc); memset(&tc, 0, sizeof tc);
	memset(&oa, 0, sizeof oa); memset(&ta, 0, sizeof ta);
	oa.alloc = rc_count_alloc; oa.dealloc = rc_count_free; oa.user_data = &oc;
	ta.alloc = rc_count_alloc; ta.dealloc = rc_count_free; ta.user_data = &tc;
	qaws_arena_init(&arena, g_rc_arena, sizeof g_rc_arena);
	qaws_exchange_thread_allocator(&ta, &prev);
	for (pass = 0; pass < 2; pass++)
	{
		qaws_repulsive_curves_desc cd;
		qaws_repulsive_curves* c = NULL;
		qaws_repulsive_curves_solve_desc d;
		qaws_repulsive_curves_solver *s = NULL, *sh = NULL;
		qaws_repulsive_mesh_desc md;
		qaws_repulsive_mesh* m = NULL;
		qaws_repulsive_mesh_solve_desc dm;
		qaws_repulsive_mesh_solver* ms = NULL;
		qaws_repulsive_report r;
		qaws_repulsive_targets t, tb;
		qaws_repulsive_derivative_report rep;
#define RC_CALL(expr) do { st = (expr); if (st != QAWS_STATUS_OK && bad == QAWS_STATUS_OK) bad = st; calls++; held += pass && arena.live != 0; } while (0)
		if (pass)
			qaws_exchange_thread_scratch_allocator(&arena.allocator, &prev_scratch);
		rc_loop(x, edges);
		memset(&cd, 0, sizeof cd);
		cd.dimension = 3; cd.vertex_count = RC_N; cd.edges = edges; cd.edge_count = RC_N;
		RC_CALL(qaws_repulsive_curves_create_ex(&cd, &oa, &c));
		if (!c)
			break;
		memset(g, 0, sizeof g);
		for (i = 0; i < 3 * RC_N; i++) dir[i] = (qaws_scalar)sin(0.3 * i);
		RC_CALL(qaws_repulsive_curves_gradient(c, x, 1, g, &E));
		RC_CALL(qaws_repulsive_curves_eval(c, x, dir, &E, NULL, NULL));
		RC_CALL(qaws_repulsive_curves_hvp(c, x, dir, g));
		RC_CALL(qaws_repulsive_curves_bh(c, x, (qaws_scalar)0.5, 1, g, &E, NULL));
		memset(&d, 0, sizeof d);
		d.collision_checked = 1;
		RC_CALL(qaws_repulsive_curves_solver_create_ex(c, &d, x, &oa, &s));
		for (it = 0; s && it < 3; it++) RC_CALL(qaws_repulsive_curves_solver_step(s, x, &r));
		RC_CALL(qaws_repulsive_curves_solver_direction(s, x, dir));
		RC_CALL(qaws_repulsive_curves_solver_objective(s, x, &E, g));
		d.collision_checked = 0;
		d.theta = (qaws_scalar)0.5;
		RC_CALL(qaws_repulsive_curves_solver_create_ex(c, &d, x, &oa, &sh));
		for (it = 0; sh && it < 2; it++) RC_CALL(qaws_repulsive_curves_solver_step(sh, x, &r));
		d.theta = 0;
		d.max_iterations = 3;
		RC_CALL(qaws_repulsive_curves_solve(c, &d, x, &r));   /* (a solver in the curves' allocator) */
		/* a mesh converged, its derivatives */
		rc_octahedron(mx, &md, &oa, &m);
		memset(&dm, 0, sizeof dm);
		dm.hold_area = 1; dm.hold_volume = 1;
		RC_CALL(m ? qaws_repulsive_mesh_solver_create_ex(m, &dm, mx, &oa, &ms) : QAWS_STATUS_INVALID_ARGUMENT);
		memset(&r, 0, sizeof r);
		for (it = 0; ms && it < 500 && r.reason == QAWS_REPULSIVE_RUNNING; it++) RC_CALL(qaws_repulsive_mesh_solver_step(ms, mx, &r));
		memset(&t, 0, sizeof t);
		t.area = 1;
		RC_CALL(qaws_repulsive_mesh_solver_solution_jvp_checked(ms, &t, dx, NULL, &rep));
		for (i = 0; i < 18; i++) xbar[i] = (qaws_scalar)cos(0.4 * i);
		memset(&tb, 0, sizeof tb);
		RC_CALL(qaws_repulsive_mesh_solver_solution_vjp(ms, xbar, &tb, NULL));
		RC_CALL(qaws_repulsive_mesh_eval(m, mx, NULL, &E, NULL, NULL));
		qaws_repulsive_mesh_solver_destroy(ms);
		qaws_repulsive_mesh_destroy(m);
		qaws_repulsive_curves_solver_destroy(s);
		qaws_repulsive_curves_solver_destroy(sh);
		qaws_repulsive_curves_destroy(c);
		held += pass && arena.live != 0;
		if (pass)
			qaws_exchange_thread_scratch_allocator(&prev_scratch, NULL);
#undef RC_CALL
	}
	qaws_set_thread_allocator(&prev);
	sprintf(msg, "objects made with _ex: %lu blocks from their allocator, %ld left after destroy; %lu from the thread allocator; %u calls, first failure %s",
		oc.allocs, oc.live, tc.allocs, calls, qaws_status_to_string(bad));
	TEST_ASSERT(bad == QAWS_STATUS_OK && oc.allocs > 0 && oc.live == 0 && tc.allocs == 0, msg);
	printf("  %s\n", msg);
	sprintf(msg, "with a scratch arena: %u calls left blocks in it, peak %lu bytes used", held, (unsigned long)arena.peak);
	TEST_ASSERT(held == 0 && arena.peak > 0, msg);
	printf("  %s\n", msg);
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
	test_settings();
	test_allocator();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
