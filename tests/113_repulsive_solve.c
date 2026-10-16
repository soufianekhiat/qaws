/*
 * Test 113: fractional Sobolev descent of the curves' tangent-point energy
 * (qaws_repulsive.h, R2)
 *
 *   - the metric A = B + B0 against Mathematica (reference/
 *     113_repulsive_metric.wls: from the paper's bilinear form, not its
 *     increment rules); symmetric, constants in its kernel, positive
 *     semidefinite
 *   - the direction: tangent to the constraints (each held edge length's
 *     derivative along it 0) and a descent direction
 *   - a trefoil with its edge lengths held converges, the energy falling at
 *     every step, the lengths held
 *   - an open 2D curve with its ends pinned, their tangents and its total
 *     length held
 *   - collision-checked: a Hopf link relaxed keeps linking number 1 (a
 *     Gauss integral over its segments, independent of the solver); the
 *     path check catches a crossing whose ends are both clear
 *   - refusals (touching input in the collision-checked mode, redundant
 *     constraints), the translation nullspace reported, the iteration limit
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "qaws_arc.h"
#include "qaws_bezier.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_repulsive.h"
#include "reference/113_repulsive_metric.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RS_F32 (sizeof(qaws_scalar) == 4)
#define RS_PI 3.14159265358979323846

static qaws_repulsive_curves* rs_make(unsigned int dim, unsigned int nv, unsigned int const* edges, unsigned int ne, double alpha, double beta)
{
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = dim; d.vertex_count = nv; d.edges = edges; d.edge_count = ne;
	d.alpha = (qaws_scalar)alpha; d.beta = (qaws_scalar)beta;
	qaws_repulsive_curves_create(&c, &d);
	return c;
}

static void test_metric(void)
{
	static struct { double const* x; unsigned int const* idx; unsigned int nv, ne, dim; double const* par; double const* A; char const* name; } k[2] = {
		{ g113_c1_x, g113_c1_idx, 7, 7, 3, g113_c1_param, g113_c1_A, "closed 3D heptagon" },
		{ g113_c2_x, g113_c2_idx, 6, 5, 2, g113_c2_param, g113_c2_A, "open 2D curve (2, 4.5)" } };
	unsigned int t, i, j;
	char msg[300];
	for (t = 0; t < 2; t++)
	{
		qaws_repulsive_curves* c = rs_make(k[t].dim, k[t].nv, k[t].idx, k[t].ne, k[t].par[0], k[t].par[1]);
		double A[64], x[24], err = 0, m = 0, asym = 0, row = 0, psd = 1e300;
		qaws_status st;
		for (i = 0; i < k[t].nv * k[t].dim; i++) x[i] = (double)(qaws_scalar)k[t].x[i];
		st = c ? qaws_internal_repulsive_curves_metric(c, x, A) : QAWS_STATUS_INVALID_ARGUMENT;
		for (i = 0; i < k[t].nv * k[t].nv; i++) { err = fmax(err, fabs(A[i] - k[t].A[i])); m = fmax(m, fabs(k[t].A[i])); }
		for (i = 0; i < k[t].nv; i++)
		{
			double s = 0;
			for (j = 0; j < k[t].nv; j++) { s += A[i * k[t].nv + j]; asym = fmax(asym, fabs(A[i * k[t].nv + j] - A[j * k[t].nv + i])); }
			row = fmax(row, fabs(s));
		}
		for (j = 0; j < 20; j++)
		{
			/* u^T A u over random u, relative to |A| |u|^2 */
			double u[8], q = 0, uu = 0;
			unsigned int a, b;
			for (a = 0; a < k[t].nv; a++) { u[a] = sin(7.1 * a + 3.3 * j + 0.2); uu += u[a] * u[a]; }
			for (a = 0; a < k[t].nv; a++) for (b = 0; b < k[t].nv; b++) q += u[a] * A[a * k[t].nv + b] * u[b];
			psd = fmin(psd, q / (m * uu));
		}
		sprintf(msg, "metric, %s: against Mathematica %.2g (status %d), asymmetry %.2g, row sums %.2g, smallest u^T A u / |A| |u|^2 %.3g",
			k[t].name, err / m, (int)st, asym / m, row / m, psd);
		TEST_ASSERT(st == QAWS_STATUS_OK && err / m <= (RS_F32 ? 1e-5 : 1e-13) && asym == 0 && row / m <= 1e-13 && psd >= -1e-13, msg);
		qaws_repulsive_curves_destroy(c);
	}
}

/* a trefoil, n vertices, closed; edges into e */
static void rs_trefoil(qaws_scalar* x, unsigned int* e, unsigned int n, double lump)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * RS_PI * i / n;
		x[3 * i] = (qaws_scalar)(sin(t) + 2 * sin(2 * t) + lump * cos(5 * t));
		x[3 * i + 1] = (qaws_scalar)(cos(t) - 2 * cos(2 * t));
		x[3 * i + 2] = (qaws_scalar)(-0.6 * sin(3 * t));
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % n;
	}
}

static double rs_len(qaws_scalar const* x, unsigned int dim, unsigned int a, unsigned int b)
{
	double s = 0;
	unsigned int q;
	for (q = 0; q < dim; q++) s += ((double)x[dim * b + q] - x[dim * a + q]) * ((double)x[dim * b + q] - x[dim * a + q]);
	return sqrt(s);
}

static void test_trefoil(void)
{
	enum { N = 48 };
	qaws_scalar x[3 * N], g[3 * N], dE[3 * N], E0 = 0;
	unsigned int e[2 * N], i, it;
	double l0[N], lerr = 0, tang = 0, gd = 0, gn = 0, prev = 1e300;
	int mono = 1;
	qaws_repulsive_curves* c;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_curves_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_status st;
	char msg[300];
	memset(&r, 0, sizeof(r));
	rs_trefoil(x, e, N, 0);
	c = rs_make(3, N, e, N, 0, 0);
	for (i = 0; i < N; i++) l0[i] = rs_len(x, 3, e[2 * i], e[2 * i + 1]);
	memset(&sd, 0, sizeof(sd));
	sd.hold_edge_lengths = 1;
	st = qaws_repulsive_curves_solver_create(&s, c, x, &sd);
	/* the direction: each edge length's derivative along it is 0; it descends */
	memset(dE, 0, sizeof(dE));
	if (st == QAWS_STATUS_OK) st = qaws_repulsive_curves_solver_direction(g, s, x);
	qaws_repulsive_curves_gradient(dE, &E0, c, x, 1);
	for (i = 0; i < N; i++)
	{
		unsigned int a = e[2 * i], b = e[2 * i + 1], q;
		double dl = 0;
		for (q = 0; q < 3; q++) dl += ((double)x[3 * b + q] - x[3 * a + q]) * ((double)g[3 * b + q] - g[3 * a + q]);
		tang = fmax(tang, fabs(dl) / l0[i]);
	}
	for (i = 0; i < 3 * N; i++) { gd += (double)dE[i] * g[i]; gn = fmax(gn, fabs((double)g[i])); }
	sprintf(msg, "trefoil: the direction along the held lengths %.2g (of |g| %.3g), descending (<dE, g> = %.4g)", tang, gn, gd);
	TEST_ASSERT(st == QAWS_STATUS_OK && tang <= (RS_F32 ? 1e-5 : 1e-12) * gn && gd > 0, msg);
	for (it = 0; it < 100 && st == QAWS_STATUS_OK; it++)
	{
		st = qaws_repulsive_curves_solver_step(&r, x, s);
		if ((double)r.energy > prev) mono = 0;
		prev = (double)r.energy;
		if (r.reason != QAWS_REPULSIVE_RUNNING)
			break;
	}
	for (i = 0; i < N; i++) lerr = fmax(lerr, fabs(rs_len(x, 3, e[2 * i], e[2 * i + 1]) - l0[i]) / l0[i]);
	sprintf(msg, "trefoil, edge lengths held: reason %d after %u steps, E %.6f from %.6f (falling at every step %d), stationarity %.2g, lengths off %.2g",
		(int)r.reason, r.iterations, (double)r.energy, (double)r.initial_energy, mono, (double)r.stationarity, lerr);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.reason == QAWS_REPULSIVE_CONVERGED && mono && r.energy < r.initial_energy &&
		lerr <= (RS_F32 ? 2e-5 : 1e-8), msg);
	qaws_repulsive_curves_solver_destroy(s);
	qaws_repulsive_curves_destroy(c);
}

static void test_pinned(void)
{
	enum { N = 21 };
	qaws_scalar x[2 * N], x0[2 * N];
	unsigned int e[2 * (N - 1)], i, pins[2] = { 0, N - 1 }, tans[2] = { 0, N - 2 };
	double L0 = 0, L1 = 0, perr, terr = 0;
	qaws_repulsive_curves* c;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_report r;
	qaws_status st;
	char msg[300];
	/* a squeezed S */
	for (i = 0; i < N; i++)
	{
		double t = (double)i / (N - 1);
		x[2 * i] = (qaws_scalar)(t * 2); x[2 * i + 1] = (qaws_scalar)(0.25 * sin(2 * RS_PI * t) + 0.05 * sin(9 * t));
		if (i + 1 < N) { e[2 * i] = i; e[2 * i + 1] = i + 1; }
	}
	memcpy(x0, x, sizeof(x));
	for (i = 0; i + 1 < N; i++) L0 += rs_len(x, 2, i, i + 1);
	c = rs_make(2, N, e, N - 1, 0, 0);
	memset(&sd, 0, sizeof(sd));
	sd.pinned = pins; sd.pinned_count = 2; sd.tangent_edges = tans; sd.tangent_edge_count = 2;
	sd.hold_total_length = 1; sd.collision_checked = 1; sd.max_iterations = 60;
	st = qaws_repulsive_curves_solve(&r, x, c, &sd);
	for (i = 0; i + 1 < N; i++) L1 += rs_len(x, 2, i, i + 1);
	perr = fmax(fmax(fabs((double)x[0] - x0[0]), fabs((double)x[1] - x0[1])), fmax(fabs((double)x[2 * N - 2] - x0[2 * N - 2]), fabs((double)x[2 * N - 1] - x0[2 * N - 1])));
	for (i = 0; i < 2; i++)
	{
		/* the end edges' directions against the initial ones (cross products) */
		unsigned int a = e[2 * tans[i]], b = e[2 * tans[i] + 1];
		double u0 = (double)x0[2 * b] - x0[2 * a], v0 = (double)x0[2 * b + 1] - x0[2 * a + 1], u1 = (double)x[2 * b] - x[2 * a], v1 = (double)x[2 * b + 1] - x[2 * a + 1];
		terr = fmax(terr, fabs(u0 * v1 - v0 * u1) / (hypot(u0, v0) * hypot(u1, v1)));
	}
	sprintf(msg, "open 2D curve, ends pinned and their tangents and the total length held (collision-checked): reason %d after %u steps, E %.5f from %.5f, pins off %.2g, tangents off %.2g, length off %.2g",
		(int)r.reason, r.iterations, (double)r.energy, (double)r.initial_energy, perr, terr, fabs(L1 - L0) / L0);
	TEST_ASSERT(st == QAWS_STATUS_OK && (r.reason == QAWS_REPULSIVE_CONVERGED || r.reason == QAWS_REPULSIVE_ITERATION_LIMIT) && r.energy < r.initial_energy &&
		perr <= (RS_F32 ? 1e-5 : 1e-10) && terr <= (RS_F32 ? 1e-5 : 1e-9) && fabs(L1 - L0) / L0 <= (RS_F32 ? 2e-5 : 1e-8), msg);
	qaws_repulsive_curves_destroy(c);
}

/* the Gauss linking integral of two closed polygons (a: na vertices, b: nb),
   Gauss-Legendre on each pair of segments */
static double rs_linking(qaws_scalar const* a, unsigned int na, qaws_scalar const* b, unsigned int nb)
{
	static double const gx[8] = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
	static double const gw[8] = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
	double s = 0;
	unsigned int i, j, p, q, c;
	for (i = 0; i < na; i++)
		for (j = 0; j < nb; j++)
		{
			double da[3], db[3];
			for (c = 0; c < 3; c++) { da[c] = (double)a[3 * ((i + 1) % na) + c] - a[3 * i + c]; db[c] = (double)b[3 * ((j + 1) % nb) + c] - b[3 * j + c]; }
			for (p = 0; p < 8; p++)
				for (q = 0; q < 8; q++)
				{
					double sa = 0.5 * (gx[p] + 1), sb = 0.5 * (gx[q] + 1), r[3], cr[3], n;
					for (c = 0; c < 3; c++) r[c] = (double)a[3 * i + c] + sa * da[c] - (double)b[3 * j + c] - sb * db[c];
					cr[0] = da[1] * db[2] - da[2] * db[1]; cr[1] = da[2] * db[0] - da[0] * db[2]; cr[2] = da[0] * db[1] - da[1] * db[0];
					n = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
					s += 0.25 * gw[p] * gw[q] * (cr[0] * r[0] + cr[1] * r[1] + cr[2] * r[2]) / (n * n * n);
				}
		}
	return s / (4 * RS_PI);
}

static void test_hopf(void)
{
	enum { M = 32 };
	qaws_scalar x[3 * 2 * M];
	unsigned int e[2 * 2 * M], i;
	double lk0, lk1;
	qaws_repulsive_curves* c;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_report r;
	qaws_status st;
	char msg[300];
	/* two circles through each other, the second squeezed toward the first */
	for (i = 0; i < M; i++)
	{
		double t = 2 * RS_PI * i / M;
		x[3 * i] = (qaws_scalar)cos(t); x[3 * i + 1] = (qaws_scalar)sin(t); x[3 * i + 2] = 0;
		x[3 * (M + i)] = (qaws_scalar)(0.9 + 0.8 * cos(t)); x[3 * (M + i) + 1] = 0; x[3 * (M + i) + 2] = (qaws_scalar)(0.45 * sin(t));
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % M;
		e[2 * (M + i)] = M + i; e[2 * (M + i) + 1] = M + (i + 1) % M;
	}
	lk0 = rs_linking(x, M, x + 3 * M, M);
	c = rs_make(3, 2 * M, e, 2 * M, 0, 0);
	memset(&sd, 0, sizeof(sd));
	sd.hold_edge_lengths = 1; sd.collision_checked = 1; sd.max_iterations = 60;
	st = qaws_repulsive_curves_solve(&r, x, c, &sd);
	lk1 = rs_linking(x, M, x + 3 * M, M);
	sprintf(msg, "Hopf link, collision-checked: reason %d after %u steps, E %.5f from %.5f, linking number %.4f -> %.4f",
		(int)r.reason, r.iterations, (double)r.energy, (double)r.initial_energy, lk0, lk1);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.iterations > 0 && r.energy < r.initial_energy && fabs(fabs(lk0) - 1) < 0.02 && fabs(lk1 - lk0) < 0.02, msg);
	qaws_repulsive_curves_destroy(c);
}

static void test_path_check(void)
{
	/* segment (0, 1) still; segment (2, 3) sweeping across it: clear at both
	   ends of the motion, crossing in between */
	static unsigned int const e[2][2] = { { 0, 1 }, { 2, 3 } };
	double a[12] = { -1, 0, 0, 1, 0, 0, 0, -1, 0.5, 0, 1, 0.5 }, b[12] = { -1, 0, 0, 1, 0, 0, 0, -1, -0.5, 0, 1, -0.5 }, p[12];
	qaws_repulsive_curves* c = rs_make(3, 4, &e[0][0], 2, 0, 0);
	int cross, slide, near;
	unsigned int i;
	char msg[300];
	cross = qaws_internal_repulsive_curves_path_clear(c, 1e-3, a, b);
	/* sliding along, never closer than 0.5 */
	for (i = 0; i < 12; i++) p[i] = a[i] + (i >= 6 && i % 3 == 0 ? 0.7 : 0);
	slide = qaws_internal_repulsive_curves_path_clear(c, 1e-3, a, p);
	/* approaching to within 5e-4: under the clearance (and half the gap) */
	for (i = 0; i < 12; i++) p[i] = a[i] - (i >= 6 && i % 3 == 2 ? 0.4995 : 0);
	near = qaws_internal_repulsive_curves_path_clear(c, 1e-3, a, p);
	sprintf(msg, "path check: a sweep through the other segment %d (ends clear), sliding past %d, closing under the clearance %d", cross, slide, near);
	TEST_ASSERT(cross == 0 && slide == 1 && near == 0, msg);
	qaws_repulsive_curves_destroy(c);
}

/* qaws curves -> working polyline (with provenance) -> solve -> exact refit */
static void test_adapter(void)
{
	qaws_arc_segment as;
	qaws_arc_desc ad;
	qaws_bezier_desc bd;
	qaws_scalar bcp[8] = { -0.5f, -0.3f, -0.1f, 0.5f, 0.1f, -0.5f, 0.5f, 0.3f };
	qaws_curve* cv[2] = { NULL, NULL };
	qaws_curve* fit[2] = { NULL, NULL };
	qaws_repulsive_polyline_desc pd;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_report r;
	qaws_scalar *x = NULL, *t = NULL;
	unsigned int *e = NULL, *src = NULL, nv = 0, ne = 0, i, loop_v = 0, loop_e = 0;
	double longest = 0, at = 0, back = 0;
	qaws_status s0, s1, s2 = QAWS_STATUS_INVALID_ARGUMENT, s3 = QAWS_STATUS_INVALID_ARGUMENT;
	char msg[400];
	memset(&r, 0, sizeof(r));
	memset(&as, 0, sizeof(as));
	as.radius = 1; as.angle_end = (qaws_scalar)(2 * RS_PI);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &as; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, &cv[0]);
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 3; bd.control_points = bcp; bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &cv[1]);
	memset(&pd, 0, sizeof(pd));
	pd.curves = (qaws_curve const* const*)cv; pd.curve_count = 2; pd.dimension = 2; pd.edge_length = (qaws_scalar)0.15;
	s0 = qaws_repulsive_polyline_from_curves(NULL, NULL, NULL, NULL, &nv, &ne, 0, 0, &pd);   /* sizing */
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (nv + 1)); t = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nv + 1));
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1)); src = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
	s1 = qaws_repulsive_polyline_from_curves(x, e, src, t, &nv, &ne, nv, ne, &pd);
	for (i = 0; s1 == QAWS_STATUS_OK && i < nv; i++)
	{
		qaws_eval_result_2d ev;
		qaws_curve_evaluate_2d(cv[src[i]], t[i], QAWS_EVAL_FLAG_POSITION, &ev);
		at = fmax(at, hypot((double)ev.position.x - x[2 * i], (double)ev.position.y - x[2 * i + 1]));
		loop_v += src[i] == 0;
	}
	for (i = 0; s1 == QAWS_STATUS_OK && i < ne; i++)
	{
		longest = fmax(longest, rs_len(x, 2, e[2 * i], e[2 * i + 1]));
		loop_e += src[e[2 * i]] == 0;
	}
	sprintf(msg, "adapter: a circle and a cubic sampled (sizing %d, %d): %u vertices, %u edges, the circle a loop (%u vertices, %u edges), longest edge %.4f (of 0.15), vertices on their curves at their parameters %.2g",
		(int)s0, (int)s1, nv, ne, loop_v, loop_e, longest, at);
	TEST_ASSERT(s0 == QAWS_STATUS_BUFFER_TOO_SMALL && s1 == QAWS_STATUS_OK && loop_v == loop_e && loop_v > 3 && longest <= 0.15 && at == 0, msg);
	if (s1 == QAWS_STATUS_OK && (c = rs_make(2, nv, e, ne, 0, 0)) != NULL)
	{
		memset(&sd, 0, sizeof(sd));
		sd.hold_edge_lengths = 1; sd.collision_checked = 1; sd.max_iterations = 25;
		s2 = qaws_repulsive_curves_solve(&r, x, c, &sd);
		s3 = qaws_repulsive_polyline_to_curve(&fit[0], x, src, t, nv, 2, 0, cv[0]);
		if (s3 == QAWS_STATUS_OK) s3 = qaws_repulsive_polyline_to_curve(&fit[1], x, src, t, nv, 2, 1, cv[1]);
		/* the refits through the moved vertices at their source parameters */
		for (i = 0; s3 == QAWS_STATUS_OK && i < nv; i++)
		{
			qaws_eval_result_2d ev;
			qaws_curve_evaluate_2d(fit[src[i]], t[i], QAWS_EVAL_FLAG_POSITION, &ev);
			back = fmax(back, hypot((double)ev.position.x - x[2 * i], (double)ev.position.y - x[2 * i + 1]));
		}
		if (s3 == QAWS_STATUS_OK)
		{
			/* the loop closes at the source's end parameter */
			qaws_eval_result_2d ev;
			qaws_curve_evaluate_2d(fit[0], qaws_curve_get_parameter_range(cv[0]).max_value, QAWS_EVAL_FLAG_POSITION, &ev);
			back = fmax(back, hypot((double)ev.position.x - x[0], (double)ev.position.y - x[1]));
		}
	}
	sprintf(msg, "adapter: solved (status %d, reason %d, E %.4f from %.4f), refitted (status %d) through every vertex at its source parameter %.2g",
		(int)s2, (int)r.reason, (double)r.energy, (double)r.initial_energy, (int)s3, back);
	TEST_ASSERT(s2 == QAWS_STATUS_OK && r.energy < r.initial_energy && s3 == QAWS_STATUS_OK && back <= (RS_F32 ? 1e-6 : 1e-14), msg);
	qaws_repulsive_curves_destroy(c);
	for (i = 0; i < 2; i++) { qaws_curve_destroy(cv[i]); qaws_curve_destroy(fit[i]); }
	free(x); free(t); free(e); free(src);
}

/* an obstacle: a loop entirely pinned, every edge length held (its rows are
   dropped, they would repeat the pins); a strand pinned at its ends beside
   it moves, the loop does not */
static void test_obstacle(void)
{
	enum { M = 13, O = 16 };
	qaws_scalar x[2 * (M + O)], x0[2 * (M + O)];
	unsigned int e[2 * (M - 1 + O)], pins[2 + O], i, ne = 0, np = 0;
	double still = 0;
	qaws_repulsive_curves* c;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_report r;
	qaws_status st;
	char msg[300];
	for (i = 0; i < M; i++)
	{
		double t = (double)i / (M - 1);
		x[2 * i] = (qaws_scalar)(-1 + 2 * t); x[2 * i + 1] = (qaws_scalar)(0.3 + 0.1 * sin(6 * t));
		if (i + 1 < M) { e[2 * ne] = i; e[2 * ne + 1] = i + 1; ne++; }
	}
	pins[np++] = 0; pins[np++] = M - 1;
	for (i = 0; i < O; i++)
	{
		double t = 2 * RS_PI * i / O;
		x[2 * (M + i)] = (qaws_scalar)(0.25 * cos(t)); x[2 * (M + i) + 1] = (qaws_scalar)(0.25 * sin(t));
		e[2 * ne] = M + i; e[2 * ne + 1] = M + (i + 1) % O; ne++;
		pins[np++] = M + i;
	}
	memcpy(x0, x, sizeof(x));
	c = rs_make(2, M + O, e, ne, 0, 0);
	memset(&sd, 0, sizeof(sd));
	sd.pinned = pins; sd.pinned_count = np; sd.hold_edge_lengths = 1; sd.collision_checked = 1; sd.max_iterations = 20;
	memset(&r, 0, sizeof(r));
	st = qaws_repulsive_curves_solve(&r, x, c, &sd);
	for (i = 2 * M; i < 2 * (M + O); i++) still = fmax(still, fabs((double)x[i] - x0[i]));
	sprintf(msg, "obstacle: a pinned loop with its edge lengths held (redundant rows dropped): reason %d after %u steps, E %.4f from %.4f, the loop moved %.2g",
		(int)r.reason, r.iterations, (double)r.energy, (double)r.initial_energy, still);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.reason != QAWS_REPULSIVE_NUMERICAL_FAILURE && r.iterations > 0 && r.energy < r.initial_energy && still <= (RS_F32 ? 1e-6 : 1e-12), msg);
	qaws_repulsive_curves_destroy(c);
}

static void test_refusals(void)
{
	static unsigned int const e[3][2] = { { 0, 1 }, { 1, 2 }, { 3, 4 } };
	qaws_scalar x[15] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0.5, -1, 0, 0.5, 1, 0 };   /* (3, 4) crosses (0, 1) */
	qaws_scalar y[15] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0.5, -1, 1, 0.5, 1, 1 };
	qaws_repulsive_curves* c = rs_make(3, 5, &e[0][0], 3, 0, 0);
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_curves_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_status s1, s2, s3, s4;
	char msg[300];
	memset(&sd, 0, sizeof(sd));
	memset(&r, 0, sizeof(r));
	sd.collision_checked = 1;
	s1 = qaws_repulsive_curves_solver_create(&s, c, x, &sd);           /* touching input */
	sd.collision_checked = 0; sd.hold_total_length = 1; sd.hold_edge_lengths = 1;
	s2 = qaws_repulsive_curves_solver_create(&s, c, y, &sd);           /* redundant lengths */
	/* no barycenter and nothing pinned: translations free, the system singular */
	sd.hold_total_length = 0; sd.hold_edge_lengths = 0; sd.hold_barycenter = -1;
	s3 = qaws_repulsive_curves_solver_create(&s, c, y, &sd);
	if (s3 == QAWS_STATUS_OK) qaws_repulsive_curves_solver_step(&r, y, s);
	qaws_repulsive_curves_solver_destroy(s);
	/* the iteration budget */
	sd.hold_barycenter = 0; sd.hold_total_length = 1; sd.max_iterations = 2;
	{
		qaws_repulsive_report r2;
		memset(&r2, 0, sizeof(r2));
		s4 = qaws_repulsive_curves_solve(&r2, y, c, &sd);
		sprintf(msg, "refusals: touching input (collision-checked) %d, total and edge lengths %d; free translations: reason %d; a budget of 2: reason %d after %u",
			(int)s1, (int)s2, (int)r.reason, (int)r2.reason, r2.iterations);
		TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_OK &&
			r.reason == QAWS_REPULSIVE_NUMERICAL_FAILURE && s4 == QAWS_STATUS_OK && r2.reason == QAWS_REPULSIVE_ITERATION_LIMIT && r2.iterations == 2, msg);
	}
	qaws_repulsive_curves_destroy(c);
}

int test_113_repulsive_solve_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 113: fractional Sobolev descent of tangent-point energies\n");
	test_metric();
	test_trefoil();
	test_pinned();
	test_hopf();
	test_path_check();
	test_adapter();
	test_obstacle();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
