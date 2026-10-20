/*
 * Test 118: native control points and the derivatives of converged
 * solutions (qaws_repulsive.h, R5)
 *
 *   - control points: a cubic B-spline, a cubic Bezier and a NURBS (weights
 *     held), sampled once; the energy's gradient in the control points
 *     against central differences, the HVP against differences of the
 *     gradient, the directional derivatives against both; the solve lowers
 *     the energy (dense and Barnes-Hut), and the curves rebuilt from the
 *     result reproduce the solved samples
 *   - a curve solution's derivative: an open strand pinned at both ends,
 *     its length held, converged; d x* / d (a pin's position) and
 *     d x* / d (the length) against re-solved differences; the adjoint
 *     identity xbar . JVP(d) = VJP(xbar) . d
 *   - a closed loop (barycenter and length held, nothing pinned): three
 *     rotations gauged out, the derivative orthogonal to them and matching
 *     the re-solved difference less its rotation
 *   - a mesh solution's derivative: an octahedral sphere, area and volume
 *     held, against the re-solved difference less its rotation
 *   - refusals: before convergence NOT_CONVERGED; the hierarchical mode
 *     UNSUPPORTED_OPERATION
 *   - the validity report: on the converged strand the KKT residual under
 *     1e-7, converged, a finite condition estimate >= 1, the unchecked JVP
 *     bit for bit; a looser gmres_tolerance takes fewer iterations; two
 *     iterations at most: NUMERICAL_FAILURE, converged 0
 *   - a mesh solution's VJP (a pin, the barycenter, area and volume held):
 *     the adjoint identity over random pairs, and the JVP for the area
 *     target against re-solved central differences
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "qaws_bspline.h"
#include "qaws_bezier.h"
#include "qaws_nurbs.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RC_F32 (sizeof(qaws_scalar) == 4)
#define RC_PI 3.14159265358979323846

static double rc_maxrel(double const* a, double const* b, unsigned int n)
{
	double e = 0, m = 0;
	unsigned int i;
	for (i = 0; i < n; i++) { e = fmax(e, fabs(a[i] - b[i])); m = fmax(m, fabs(b[i])); }
	return e / fmax(m, 1e-300);
}

/* three curves: a cubic B-spline, a cubic Bezier, a NURBS */
static void rc_curves(qaws_curve** c)
{
	static qaws_scalar const knots[12] = { 0, 0, 0, 0, (qaws_scalar)0.2, (qaws_scalar)0.4, (qaws_scalar)0.6, (qaws_scalar)0.8, 1, 1, 1, 1 };
	static qaws_scalar const w[8] = { 1, (qaws_scalar)1.5, (qaws_scalar)0.7, 1, (qaws_scalar)2, (qaws_scalar)0.8, (qaws_scalar)1.2, 1 };
	qaws_scalar p[3][24];
	qaws_bspline_desc bd;
	qaws_bezier_desc zd;
	qaws_nurbs_desc nd;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		double t = (double)i / 7;
		p[0][3 * i] = (qaws_scalar)(2 * t); p[0][3 * i + 1] = (qaws_scalar)(0.5 * sin(5 * t)); p[0][3 * i + 2] = (qaws_scalar)(0.3 * cos(4 * t));
		p[2][3 * i] = (qaws_scalar)(2 * t); p[2][3 * i + 1] = (qaws_scalar)(1.4 + 0.4 * cos(6 * t)); p[2][3 * i + 2] = (qaws_scalar)(0.5 * sin(3 * t));
	}
	for (i = 0; i < 4; i++)
	{
		double t = (double)i / 3;
		p[1][3 * i] = (qaws_scalar)(0.2 + 1.6 * t); p[1][3 * i + 1] = (qaws_scalar)(0.7 + 0.3 * sin(4 * t)); p[1][3 * i + 2] = (qaws_scalar)(1 - 0.5 * t);
	}
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_3D; bd.degree = 3; bd.control_points = p[0]; bd.control_point_count = 8; bd.knots = knots; bd.knot_count = 12;
	qaws_curve_create_bspline(&c[0], &bd);
	memset(&zd, 0, sizeof(zd));
	zd.dimension = QAWS_DIMENSION_3D; zd.degree = 3; zd.control_points = p[1]; zd.control_point_count = 4;
	qaws_curve_create_bezier(&c[1], &zd);
	memset(&nd, 0, sizeof(nd));
	nd.dimension = QAWS_DIMENSION_3D; nd.degree = 3; nd.control_points = p[2]; nd.control_point_count = 8; nd.knots = knots; nd.knot_count = 12;
	nd.weights = w; nd.weight_count = 8;
	qaws_curve_create_nurbs(&c[2], &nd);
}

static void test_controls(void)
{
	qaws_curve* c[3] = { NULL, NULL, NULL };
	qaws_repulsive_controls_desc d;
	qaws_repulsive_controls* ctl = NULL;
	qaws_status st;
	unsigned int nc = 0, nv = 0, ne = 0, i, k;
	char msg[300];
	rc_curves(c);
	memset(&d, 0, sizeof(d));
	d.polyline.curves = (qaws_curve const* const*)c; d.polyline.curve_count = 3; d.polyline.dimension = 3; d.polyline.edge_length = (qaws_scalar)0.12;
	st = qaws_repulsive_controls_create(&ctl, &d);
	qaws_repulsive_controls_get_counts(&nc, &nv, &ne, ctl);
	sprintf(msg, "controls: created %d (%u control points, %u samples, %u edges)", (int)st, nc, nv, ne);
	TEST_ASSERT(st == QAWS_STATUS_OK && nc == 20 && nv > 40, msg);
	if (st != QAWS_STATUS_OK)
		goto done;
	{
		unsigned int n = 3 * nc;
		qaws_scalar *P = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n), *Q = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n), *g = (qaws_scalar*)calloc(n, sizeof(qaws_scalar));
		qaws_scalar *gp = (qaws_scalar*)calloc(n, sizeof(qaws_scalar)), *gm = (qaws_scalar*)calloc(n, sizeof(qaws_scalar)), *hv = (qaws_scalar*)calloc(n, sizeof(qaws_scalar));
		qaws_scalar *dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n), E = 0, t = 0, t2 = 0;
		double *fd = (double*)malloc(sizeof(double) * n), *ga = (double*)malloc(sizeof(double) * n), *hfd = (double*)malloc(sizeof(double) * n), *ha = (double*)malloc(sizeof(double) * n);
		double h = RC_F32 ? 2e-3 : 1e-5, gdot = 0, eg, eh, et;
		qaws_repulsive_controls_read(P, ctl);
		for (i = 0; i < n; i++) dir[i] = (qaws_scalar)sin(0.7 * i + 0.3);
		qaws_repulsive_controls_gradient(g, &E, ctl, P, 1);
		qaws_repulsive_controls_hvp(hv, ctl, P, dir);
		qaws_repulsive_controls_eval(NULL, &t, &t2, ctl, P, dir);
		for (i = 0; i < n; i++)
		{
			qaws_scalar ep = 0, em = 0;
			memcpy(Q, P, sizeof(qaws_scalar) * n);
			Q[i] = (qaws_scalar)((double)P[i] + h);
			qaws_repulsive_controls_eval(&ep, NULL, NULL, ctl, Q, NULL);
			Q[i] = (qaws_scalar)((double)P[i] - h);
			qaws_repulsive_controls_eval(&em, NULL, NULL, ctl, Q, NULL);
			fd[i] = ((double)ep - em) / (2 * h);
			ga[i] = (double)g[i];
			gdot += (double)g[i] * dir[i];
		}
		for (i = 0; i < n; i++) { Q[i] = (qaws_scalar)((double)P[i] + h * dir[i]); gp[i] = gm[i] = 0; }
		qaws_repulsive_controls_gradient(gp, NULL, ctl, Q, 1);
		for (i = 0; i < n; i++) Q[i] = (qaws_scalar)((double)P[i] - h * dir[i]);
		qaws_repulsive_controls_gradient(gm, NULL, ctl, Q, 1);
		for (i = 0; i < n; i++) { hfd[i] = ((double)gp[i] - gm[i]) / (2 * h); ha[i] = (double)hv[i]; }
		eg = rc_maxrel(ga, fd, n);
		eh = rc_maxrel(ha, hfd, n);
		{
			double dhd = 0;
			for (i = 0; i < n; i++) dhd += (double)dir[i] * hv[i];
			et = fmax(fabs((double)t - gdot) / fabs(gdot), fabs((double)t2 - dhd) / fabs(dhd));
		}
		sprintf(msg, "controls: the gradient against differences %.2g, the HVP against differences of the gradient %.2g, t and t2 against them %.2g", eg, eh, et);
		TEST_ASSERT(eg <= (RC_F32 ? 2e-2 : 1e-7) && eh <= (RC_F32 ? 5e-2 : 1e-6) && et <= (RC_F32 ? 1e-3 : 1e-10), msg);
		/* solves, dense and Barnes-Hut; the curves rebuilt reproduce the samples */
		for (k = 0; k < 2; k++)
		{
			qaws_repulsive_curves_solve_desc sd;
			qaws_repulsive_report r;
			qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			unsigned int* src = NULL;
			double rebuilt = 0;
			memcpy(Q, P, sizeof(qaws_scalar) * n);
			memset(&sd, 0, sizeof(sd));
			sd.hold_total_length = 1; sd.max_iterations = 15; sd.theta = (qaws_scalar)(k ? 0.25 : 0);
			memset(&r, 0, sizeof(r));
			st = qaws_repulsive_controls_solve(&r, Q, ctl, &sd);
			qaws_repulsive_controls_positions(x, ctl, Q);
			{
				/* the samples' sources, from the polyline again */
				qaws_scalar* ts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
				unsigned int nv2, ne2, ci;
				src = (unsigned int*)malloc(sizeof(unsigned int) * nv);
				qaws_repulsive_polyline_from_curves(NULL, NULL, NULL, NULL, &nv2, &ne2, 0, 0, &d.polyline);
				{
					qaws_scalar* pos = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv2);
					unsigned int* ed = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne2);
					qaws_repulsive_polyline_from_curves(pos, ed, src, ts, &nv2, &ne2, nv2, ne2, &d.polyline);
					free(pos); free(ed);
				}
				for (ci = 0; ci < 3; ci++)
				{
					qaws_curve* cc = NULL;
					qaws_repulsive_controls_to_curve(&cc, ctl, Q, ci);
					for (i = 0; cc && i < nv; i++)
						if (src[i] == ci)
						{
							qaws_eval_result_3d er;
							qaws_curve_evaluate_3d(&er, cc, ts[i], QAWS_EVAL_FLAG_POSITION);
							rebuilt = fmax(rebuilt, fabs((double)er.position.x - x[3 * i]) + fabs((double)er.position.y - x[3 * i + 1]) + fabs((double)er.position.z - x[3 * i + 2]));
						}
					if (!cc) rebuilt = 1;
					qaws_curve_destroy(cc);
				}
				free(ts);
			}
			sprintf(msg, "controls, solve (%s): status %d, reason %d after %u steps, E %.4g -> %.4g, length residual %.2g; rebuilt curves off the samples %.2g",
				k ? "Barnes-Hut, theta 0.25" : "dense", (int)st, (int)r.reason, r.iterations, (double)r.initial_energy, (double)r.energy, (double)r.constraint_residual, rebuilt);
			TEST_ASSERT(st == QAWS_STATUS_OK && r.iterations >= 5 && r.energy < 0.95 * r.initial_energy && (double)r.constraint_residual <= (RC_F32 ? 1e-4 : 1e-7) &&
				rebuilt <= (RC_F32 ? 1e-4 : 1e-10), msg);
			free(x); free(src);
		}
		free(P); free(Q); free(g); free(gp); free(gm); free(hv); free(dir); free(fd); free(ga); free(hfd); free(ha);
	}
done:
	qaws_repulsive_controls_destroy(ctl);
	for (i = 0; i < 3; i++) qaws_curve_destroy(c[i]);
}

/* ---- solution derivatives ---- */

/* solve to convergence (tight); the solver kept for derivatives */
static qaws_repulsive_curves_solver* rc_converge(qaws_repulsive_curves const* c, qaws_repulsive_curves_solve_desc const* d0, qaws_scalar* x, qaws_repulsive_report* r)
{
	qaws_repulsive_curves_solve_desc d = *d0;
	qaws_repulsive_curves_solver* s = NULL;
	unsigned int it;
	d.tolerance = (qaws_scalar)(RC_F32 ? 1e-4 : 1e-7);
	if (qaws_repulsive_curves_solver_create(&s, c, x, &d) != QAWS_STATUS_OK)
		return NULL;
	memset(r, 0, sizeof(*r));
	for (it = 0; it < 3000; it++)
	{
		qaws_repulsive_curves_solver_step(r, x, s);
		if (r->reason != QAWS_REPULSIVE_RUNNING)
			break;
	}
	return s;
}

/* rotations about the barycenter removed from d (least squares) */
static void rc_unrotate(double* d, qaws_scalar const* x, unsigned int nv)
{
	double b[3] = { 0, 0, 0 }, R[3][64 * 3], G[3][3], rhs[3], w[3];
	unsigned int i, q, a, c;
	for (i = 0; i < nv; i++)
		for (q = 0; q < 3; q++) b[q] += (double)x[3 * i + q] / nv;
	for (a = 0; a < 3; a++)
		for (i = 0; i < nv; i++)
		{
			double r[3];
			for (q = 0; q < 3; q++) r[q] = (double)x[3 * i + q] - b[q];
			R[a][3 * i + a] = 0; R[a][3 * i + (a + 1) % 3] = -r[(a + 2) % 3]; R[a][3 * i + (a + 2) % 3] = r[(a + 1) % 3];
		}
	for (a = 0; a < 3; a++)
	{
		rhs[a] = 0;
		for (i = 0; i < 3 * nv; i++) rhs[a] += R[a][i] * d[i];
		for (c = 0; c < 3; c++)
		{
			G[a][c] = 0;
			for (i = 0; i < 3 * nv; i++) G[a][c] += R[a][i] * R[c][i];
		}
	}
	/* 3x3 solve by Cramer */
	{
		double det = G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1]) - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0]) + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]);
		for (a = 0; a < 3; a++)
		{
			double M[3][3];
			memcpy(M, G, sizeof(M));
			for (c = 0; c < 3; c++) M[c][a] = rhs[c];
			w[a] = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / det;
		}
	}
	for (i = 0; i < 3 * nv; i++)
		for (a = 0; a < 3; a++) d[i] -= w[a] * R[a][i];
}

static void test_curve_solution(void)
{
	enum { N = 14 };
	qaws_scalar x0[3 * N], x[3 * N], xp[3 * N], xm[3 * N], dx[3 * N], xbar[3 * N], dpin[6];
	unsigned int e[2 * N], pins[2] = { 0, N - 1 }, i, gauge = 9;
	qaws_repulsive_curves_desc cd;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_curves_solver *s = NULL, *sp, *sm;
	qaws_repulsive_report r, rp, rm;
	qaws_repulsive_targets t, tb;
	qaws_status st;
	double h = RC_F32 ? 1e-2 : 1e-3, fd[3 * N], an[3 * N], err, lhs = 0, rhs = 0;
	double orth_s = 1;
	char msg[300];
	for (i = 0; i < N; i++)
	{
		double u = (double)i / (N - 1);
		x0[3 * i] = (qaws_scalar)(1.6 * u); x0[3 * i + 1] = (qaws_scalar)(0.35 * sin(3 * RC_PI * u)); x0[3 * i + 2] = (qaws_scalar)(0.2 * sin(2 * RC_PI * u));
		e[2 * i] = i; e[2 * i + 1] = i + 1;
	}
	memset(&cd, 0, sizeof(cd));
	cd.dimension = 3; cd.vertex_count = N; cd.edges = e; cd.edge_count = N - 1;
	qaws_repulsive_curves_create(&c, &cd);
	memset(&d, 0, sizeof(d));
	d.pinned = pins; d.pinned_count = 2; d.hold_edge_lengths = 1;
	/* before convergence: refused */
	memcpy(x, x0, sizeof(x));
	{
		qaws_repulsive_curves_solver* s0 = NULL;
		qaws_status s1, s2;
		memset(&t, 0, sizeof(t));
		qaws_repulsive_curves_solver_create(&s0, c, x, &d);
		s1 = qaws_repulsive_curves_solver_solution_jvp(dx, NULL, s0, &t);
		qaws_repulsive_curves_solver_destroy(s0);
		s0 = NULL;
		d.theta = (qaws_scalar)0.25;
		qaws_repulsive_curves_solver_create(&s0, c, x, &d);
		s2 = qaws_repulsive_curves_solver_solution_jvp(dx, NULL, s0, &t);
		qaws_repulsive_curves_solver_destroy(s0);
		d.theta = 0;
		sprintf(msg, "solution derivatives refused: not converged %d, hierarchical %d", (int)s1, (int)s2);
		TEST_ASSERT(s1 == QAWS_STATUS_NOT_CONVERGED && s2 == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	}
	s = rc_converge(c, &d, x, &r);
	/* a pin moved (the last, along y): its edge's length target moves with it */
	memcpy(xp, x, sizeof(x)); memcpy(xm, x, sizeof(x));
	xp[3 * (N - 1) + 1] = (qaws_scalar)((double)xp[3 * (N - 1) + 1] + h);
	xm[3 * (N - 1) + 1] = (qaws_scalar)((double)xm[3 * (N - 1) + 1] - h);
	sp = rc_converge(c, &d, xp, &rp);
	sm = rc_converge(c, &d, xm, &rm);
	{
		qaws_scalar dl[N];
		for (i = 0; i + 1 < N; i++)
		{
			unsigned int q;
			double b = 0, cc = 0;
			for (q = 0; q < 3; q++)
			{
				double v = (double)xp[3 * (i + 1) + q] - xp[3 * i + q], w = (double)xm[3 * (i + 1) + q] - xm[3 * i + q];
				b += v * v; cc += w * w;
			}
			dl[i] = (qaws_scalar)((sqrt(b) - sqrt(cc)) / (2 * h));
		}
		memset(&t, 0, sizeof(t));
		memset(dpin, 0, sizeof(dpin));
		dpin[4] = 1;
		t.pinned = dpin; t.edge_lengths = dl;
		st = qaws_repulsive_curves_solver_solution_jvp(dx, &gauge, s, &t);
	}
	for (i = 0; i < 3 * N; i++) { fd[i] = ((double)xp[i] - xm[i]) / (2 * h); an[i] = (double)dx[i]; }
	{
		/* (less the rotation about the pins' line, free between re-solves) */
		double a[3], R[3 * N], rr = 0, fr = 0, an2 = 0;
		unsigned int q;
		for (q = 0; q < 3; q++) a[q] = (double)x[3 * (N - 1) + q] - x[q];
		for (i = 0; i < N; i++)
		{
			double p[3];
			for (q = 0; q < 3; q++) p[q] = (double)x[3 * i + q] - x[q];
			R[3 * i] = a[1] * p[2] - a[2] * p[1]; R[3 * i + 1] = a[2] * p[0] - a[0] * p[2]; R[3 * i + 2] = a[0] * p[1] - a[1] * p[0];
		}
		for (i = 0; i < 3 * N; i++) { rr += R[i] * R[i]; fr += fd[i] * R[i]; an2 += an[i] * R[i]; }
		for (i = 0; i < 3 * N; i++) fd[i] -= fr / rr * R[i];
		orth_s = fabs(an2) / sqrt(rr);
	}
	err = rc_maxrel(an, fd, 3 * N);
	sprintf(msg, "strand, edge lengths held, a pin moved: converged %d %d %d, JVP status %d, gauge %u (the rotation about the pins' line), the JVP's part along it %.2g, against re-solved differences less it %.2g",
		(int)r.reason, (int)rp.reason, (int)rm.reason, (int)st, gauge, orth_s, err);
	TEST_ASSERT(st == QAWS_STATUS_OK && gauge == 1 && orth_s <= 1e-8 && r.reason == QAWS_REPULSIVE_CONVERGED && err <= (RC_F32 ? 5e-2 : 1e-4), msg);
	/* the adjoint identity */
	{
		qaws_scalar dl[N], lb[N];
		memset(&t, 0, sizeof(t));
		for (i = 0; i < 6; i++) dpin[i] = (qaws_scalar)cos(1.3 * i + 0.2);
		for (i = 0; i + 1 < N; i++) dl[i] = (qaws_scalar)sin(2.1 * i + 0.5);
		t.pinned = dpin; t.edge_lengths = dl;
		qaws_repulsive_curves_solver_solution_jvp(dx, NULL, s, &t);
		for (i = 0; i < 3 * N; i++) { xbar[i] = (qaws_scalar)sin(0.9 * i + 0.4); lhs += (double)xbar[i] * dx[i]; }
		{
			qaws_scalar pb[6];
			memset(&tb, 0, sizeof(tb));
			tb.pinned = pb; tb.edge_lengths = lb;
			st = qaws_repulsive_curves_solver_solution_vjp(&tb, NULL, s, xbar);
			for (i = 0; i < 6; i++) rhs += (double)pb[i] * dpin[i];
			for (i = 0; i + 1 < N; i++) rhs += (double)lb[i] * dl[i];
		}
		sprintf(msg, "strand, the adjoint identity: xbar . JVP(d) %.12g, VJP(xbar) . d %.12g (status %d)", lhs, rhs, (int)st);
		TEST_ASSERT(st == QAWS_STATUS_OK && fabs(lhs - rhs) <= (RC_F32 ? 1e-3 : 1e-9) * fabs(lhs), msg);
	}
	qaws_repulsive_curves_solver_destroy(s); qaws_repulsive_curves_solver_destroy(sp); qaws_repulsive_curves_solver_destroy(sm);
	qaws_repulsive_curves_destroy(c);
}

static void test_loop_solution(void)
{
	enum { N = 16 };
	qaws_scalar x[3 * N], xp[3 * N], xm[3 * N], dx[3 * N];
	unsigned int e[2 * N], i, gauge = 0;
	qaws_repulsive_curves_desc cd;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_curves_solver *s, *sp, *sm;
	qaws_repulsive_report r, rp, rm;
	qaws_repulsive_targets t;
	qaws_status st;
	double h = RC_F32 ? 1e-2 : 1e-4, fd[3 * N], an[3 * N], err, orth = 0, Lp = 0, Lm = 0;
	char msg[300];
	for (i = 0; i < N; i++)
	{
		double u = 2 * RC_PI * i / N;
		x[3 * i] = (qaws_scalar)(cos(u) + 0.2 * cos(2 * u)); x[3 * i + 1] = (qaws_scalar)(sin(u) - 0.15 * sin(3 * u)); x[3 * i + 2] = (qaws_scalar)(0.3 * sin(2 * u));
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % N;
	}
	memset(&cd, 0, sizeof(cd));
	cd.dimension = 3; cd.vertex_count = N; cd.edges = e; cd.edge_count = N;
	qaws_repulsive_curves_create(&c, &cd);
	memset(&d, 0, sizeof(d));
	d.hold_edge_lengths = 1;
	s = rc_converge(c, &d, x, &r);
	/* the re-solved difference: the loop scaled by 1 +- h / L about its barycenter */
	{
		double b[3] = { 0, 0, 0 }, L = 0;
		unsigned int q;
		for (i = 0; i < N; i++)
		{
			double a = 0;
			for (q = 0; q < 3; q++) { double u = (double)x[3 * ((i + 1) % N) + q] - x[3 * i + q]; a += u * u; b[q] += (double)x[3 * i + q] / N; }
			L += sqrt(a);
		}
		for (i = 0; i < N; i++)
			for (q = 0; q < 3; q++)
			{
				xp[3 * i + q] = (qaws_scalar)(b[q] + ((double)x[3 * i + q] - b[q]) * (1 + h / L));
				xm[3 * i + q] = (qaws_scalar)(b[q] + ((double)x[3 * i + q] - b[q]) * (1 - h / L));
			}
		Lp = L + h; Lm = L - h;
		{
			/* every edge scaled alike: d len_e / dL = len_e / L */
			qaws_scalar dl[N];
			for (i = 0; i < N; i++)
			{
				double a = 0;
				for (q = 0; q < 3; q++) { double u = (double)x[3 * ((i + 1) % N) + q] - x[3 * i + q]; a += u * u; }
				dl[i] = (qaws_scalar)(sqrt(a) / L);
			}
			memset(&t, 0, sizeof(t));
			t.edge_lengths = dl;
			st = qaws_repulsive_curves_solver_solution_jvp(dx, &gauge, s, &t);
		}
	}
	sp = rc_converge(c, &d, xp, &rp);
	sm = rc_converge(c, &d, xm, &rm);
	for (i = 0; i < 3 * N; i++) { fd[i] = ((double)xp[i] - xm[i]) / (Lp - Lm); an[i] = (double)dx[i]; }
	rc_unrotate(fd, x, N);
	{
		double tmp[3 * N];
		memcpy(tmp, an, sizeof(tmp));
		rc_unrotate(tmp, x, N);
		for (i = 0; i < 3 * N; i++) orth = fmax(orth, fabs(tmp[i] - an[i]));
	}
	err = rc_maxrel(an, fd, 3 * N);
	sprintf(msg, "loop, nothing pinned, edge lengths held and scaled: converged %d %d %d, gauge %u, the JVP's rotation part %.2g, against the re-solved difference less its rotation %.2g",
		(int)r.reason, (int)rp.reason, (int)rm.reason, gauge, orth, err);
	TEST_ASSERT(st == QAWS_STATUS_OK && gauge == 3 && orth <= 1e-8 && err <= (RC_F32 ? 5e-2 : 1e-4), msg);
	qaws_repulsive_curves_solver_destroy(s); qaws_repulsive_curves_solver_destroy(sp); qaws_repulsive_curves_solver_destroy(sm);
	qaws_repulsive_curves_destroy(c);
}

static void rc_measures(qaws_scalar const* x, unsigned int const* tris, unsigned int nf, double* A, double* V)
{
	unsigned int f, q;
	*A = *V = 0;
	for (f = 0; f < nf; f++)
	{
		double a[3], b[3], c[3], u[3], v[3], n[3];
		for (q = 0; q < 3; q++) { a[q] = x[3 * tris[3 * f] + q]; b[q] = x[3 * tris[3 * f + 1] + q]; c[q] = x[3 * tris[3 * f + 2] + q]; u[q] = b[q] - a[q]; v[q] = c[q] - a[q]; }
		n[0] = u[1] * v[2] - u[2] * v[1]; n[1] = u[2] * v[0] - u[0] * v[2]; n[2] = u[0] * v[1] - u[1] * v[0];
		*A += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		*V += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
	}
}

static qaws_repulsive_mesh_solver* rc_mesh_converge(qaws_repulsive_mesh const* m, qaws_scalar* x, qaws_repulsive_report* r)
{
	qaws_repulsive_mesh_solve_desc d;
	qaws_repulsive_mesh_solver* s = NULL;
	unsigned int it;
	memset(&d, 0, sizeof(d));
	d.hold_area = 1; d.hold_volume = 1; d.tolerance = (qaws_scalar)(RC_F32 ? 1e-4 : 1e-6);
	if (qaws_repulsive_mesh_solver_create(&s, m, x, &d) != QAWS_STATUS_OK)
		return NULL;
	memset(r, 0, sizeof(*r));
	for (it = 0; it < 3000; it++)
	{
		qaws_repulsive_mesh_solver_step(r, x, s);
		if (r->reason != QAWS_REPULSIVE_RUNNING)
			break;
	}
	return s;
}

static void test_mesh_solution(void)
{
	qaws_scalar x[18], xp[18], xm[18], dx[18];
	unsigned int tris[24], nv = 6, nf = 8, i, gauge = 0, q;
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[24] = { 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5 };
	qaws_repulsive_mesh_desc md;
	qaws_repulsive_mesh* m = NULL;
	qaws_repulsive_mesh_solver *s, *sp, *sm;
	qaws_repulsive_report r, rp, rm;
	qaws_repulsive_targets t;
	qaws_status st;
	double h = RC_F32 ? 1e-2 : 5e-3, fd[18], an[18], err, A, V, Ap, Vp, Am, Vm, b[3] = { 0, 0, 0 };
	char msg[300];
	/* a lumpy octahedron (no vertex free to slide within the surface) */
	for (i = 0; i < 6; i++)
		for (q = 0; q < 3; q++) x[3 * i + q] = (qaws_scalar)(ov[i][q] * (1 + 0.15 * sin(2.0 * i + q)));
	memcpy(tris, of, sizeof(of));
	memset(&md, 0, sizeof(md));
	md.vertex_count = nv; md.triangles = tris; md.triangle_count = nf;
	qaws_repulsive_mesh_create(&m, &md);
	s = rc_mesh_converge(m, x, &r);
	rc_measures(x, tris, nf, &A, &V);
	/* the re-solved difference: the solution pushed out along x - b, by +-h */
	for (i = 0; i < nv; i++)
		for (q = 0; q < 3; q++) b[q] += (double)x[3 * i + q] / nv;
	for (i = 0; i < nv; i++)
		for (q = 0; q < 3; q++)
		{
			double w = ((double)x[3 * i + q] - b[q]) * (1 + 0.3 * sin(3.0 * i));
			xp[3 * i + q] = (qaws_scalar)((double)x[3 * i + q] + h * w);
			xm[3 * i + q] = (qaws_scalar)((double)x[3 * i + q] - h * w);
		}
	rc_measures(xp, tris, nf, &Ap, &Vp);
	rc_measures(xm, tris, nf, &Am, &Vm);
	sp = rc_mesh_converge(m, xp, &rp);
	sm = rc_mesh_converge(m, xm, &rm);
	memset(&t, 0, sizeof(t));
	t.area = (qaws_scalar)((Ap - Am) / (2 * h)); t.volume = (qaws_scalar)((Vp - Vm) / (2 * h));
	for (i = 0; i < nv; i++)
		for (q = 0; q < 3; q++) t.barycenter[q] += (qaws_scalar)(((double)xp[3 * i + q] - xm[3 * i + q]) / (2 * h * nv));   /* (each re-solve holds its input's barycenter) */
	st = qaws_repulsive_mesh_solver_solution_jvp(dx, &gauge, s, &t);
	for (i = 0; i < 3 * nv; i++) { fd[i] = ((double)xp[i] - xm[i]) / (2 * h); an[i] = (double)dx[i]; }
	rc_unrotate(fd, x, nv);
	err = rc_maxrel(an, fd, 3 * nv);
	sprintf(msg, "mesh (a lumpy octahedron, %u faces), area and volume changed: converged %d %d %d, gauge %u, against the re-solved difference less its rotation %.2g (A %.4g, V %.4g)",
		nf, (int)r.reason, (int)rp.reason, (int)rm.reason, gauge, err, A, V);
	/* (the re-solves converge to about 1e-7 in stationarity: the differences reach 5e-4 at best, h 3e-3 to 1e-2) */
	TEST_ASSERT(st == QAWS_STATUS_OK && gauge == 3 && r.reason == QAWS_REPULSIVE_CONVERGED && err <= (RC_F32 ? 5e-2 : 2e-3), msg);
	qaws_repulsive_mesh_solver_destroy(s); qaws_repulsive_mesh_solver_destroy(sp); qaws_repulsive_mesh_solver_destroy(sm);
	qaws_repulsive_mesh_destroy(m);
}

/* the validity report of a derivative: the strand of test_curve_solution
   converged three times alike, its solvers differing only in their GMRES
   limits (the dense descent does not use them) */
static void test_derivative_report(void)
{
	enum { N = 14 };
	qaws_scalar x0[3 * N], x[3 * N], dx[3 * N], dx2[3 * N], xbar[3 * N], dpin[6], pb[6];
	unsigned int e[2 * N], pins[2] = { 0, N - 1 }, i, k, gauge = 9, same = 1;
	qaws_repulsive_curves_desc cd;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_curves_solver* s[3] = { NULL, NULL, NULL };
	qaws_repulsive_curves_solver* s0 = NULL;
	qaws_repulsive_report r[3];
	qaws_repulsive_targets t, tb;
	qaws_repulsive_derivative_report rep[3], repv, repn;
	qaws_status st[3], stu, stv, stn;
	char msg[400];
	for (i = 0; i < N; i++)
	{
		double u = (double)i / (N - 1);
		x0[3 * i] = (qaws_scalar)(1.6 * u); x0[3 * i + 1] = (qaws_scalar)(0.35 * sin(3 * RC_PI * u)); x0[3 * i + 2] = (qaws_scalar)(0.2 * sin(2 * RC_PI * u));
		e[2 * i] = i; e[2 * i + 1] = i + 1;
	}
	memset(&cd, 0, sizeof(cd));
	cd.dimension = 3; cd.vertex_count = N; cd.edges = e; cd.edge_count = N - 1;
	qaws_repulsive_curves_create(&c, &cd);
	memset(&d, 0, sizeof(d));
	d.pinned = pins; d.pinned_count = 2; d.hold_edge_lengths = 1;
	/* refused before convergence: the report zeroed */
	memset(&t, 0, sizeof(t));
	memset(dpin, 0, sizeof(dpin));
	dpin[4] = 1;
	t.pinned = dpin;
	memcpy(x, x0, sizeof(x));
	qaws_repulsive_curves_solver_create(&s0, c, x, &d);
	memset(&repn, 0x7f, sizeof(repn));
	stn = qaws_repulsive_curves_solver_solution_jvp_checked(s0, &t, dx, NULL, &repn);
	qaws_repulsive_curves_solver_destroy(s0);
	for (k = 0; k < 3; k++)
	{
		qaws_repulsive_curves_solve_desc dk = d;
		if (k == 1) dk.gmres_tolerance = (qaws_scalar)1e-4;
		if (k == 2) dk.gmres_max_iterations = 2;
		memcpy(x, x0, sizeof(x));
		s[k] = rc_converge(c, &dk, x, &r[k]);
		st[k] = s[k] ? qaws_repulsive_curves_solver_solution_jvp_checked(s[k], &t, dx, k ? NULL : &gauge, &rep[k]) : QAWS_STATUS_INVALID_ARGUMENT;
		if (k == 0) memcpy(dx2, dx, sizeof(dx));
	}
	stu = s[0] ? qaws_repulsive_curves_solver_solution_jvp(dx, NULL, s[0], &t) : QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * N; i++) same &= dx[i] == dx2[i];
	for (i = 0; i < 3 * N; i++) xbar[i] = (qaws_scalar)sin(0.9 * i + 0.4);
	memset(&tb, 0, sizeof(tb));
	tb.pinned = pb;
	stv = s[0] ? qaws_repulsive_curves_solver_solution_vjp_checked(s[0], xbar, &tb, NULL, &repv) : QAWS_STATUS_INVALID_ARGUMENT;
	sprintf(msg, "refused before convergence: status %d, report zeroed %d", (int)stn, repn.converged == 0 && repn.iterations == 0 && repn.residual == 0);
	TEST_ASSERT(stn == QAWS_STATUS_NOT_CONVERGED && repn.converged == 0 && repn.iterations == 0 && repn.residual == 0, msg);
	sprintf(msg, "converged strand, default limits: JVP status %d, residual %.3g, %u iterations, converged %d, condition estimate %.4g, gauge %u (out %u); "
		"the unchecked JVP bit for bit %d; the VJP status %d, residual %.3g, converged %d",
		(int)st[0], (double)rep[0].residual, rep[0].iterations, rep[0].converged, (double)rep[0].condition_estimate, rep[0].gauge, gauge, (int)same,
		(int)stv, (double)repv.residual, repv.converged);
	TEST_ASSERT(st[0] == QAWS_STATUS_OK && stu == QAWS_STATUS_OK && rep[0].converged == 1 && (double)rep[0].residual <= 1e-7 && rep[0].iterations > 0 &&
		(double)rep[0].condition_estimate >= 1 && (double)rep[0].condition_estimate < 1e30 && rep[0].gauge == 1 && gauge == 1 && same &&
		stv == QAWS_STATUS_OK && repv.converged == 1 && (double)repv.residual <= 1e-7, msg);
	printf("  %s\n", msg);
	sprintf(msg, "gmres_tolerance 1e-4: status %d, residual %.3g, %u iterations (default %u), converged %d",
		(int)st[1], (double)rep[1].residual, rep[1].iterations, rep[0].iterations, rep[1].converged);
	TEST_ASSERT(st[1] == QAWS_STATUS_OK && rep[1].converged == 1 && (double)rep[1].residual <= 1e-4 && rep[1].iterations < rep[0].iterations, msg);
	printf("  %s\n", msg);
	sprintf(msg, "gmres_max_iterations 2: status %d (NUMERICAL_FAILURE), residual %.3g, %u iterations, converged %d, condition estimate %.4g",
		(int)st[2], (double)rep[2].residual, rep[2].iterations, rep[2].converged, (double)rep[2].condition_estimate);
	TEST_ASSERT(st[2] == QAWS_STATUS_NUMERICAL_FAILURE && rep[2].converged == 0 && rep[2].iterations == 2 && (double)rep[2].residual > 1e-7 &&
		(double)rep[2].residual < 1e30 && (double)rep[2].condition_estimate >= 1, msg);
	printf("  %s\n", msg);
	for (k = 0; k < 3; k++) qaws_repulsive_curves_solver_destroy(s[k]);
	qaws_repulsive_curves_destroy(c);
}

/* the least-norm move d with C d = e_area: C the Jacobian of the mesh's barycenter (3 rows), the pinned vertex (3),
   the area and the volume, at x; d = C^T (C C^T)^-1 e_area by Gaussian elimination */
static void rc_area_move(qaws_scalar const* x, unsigned int const* tris, unsigned int nv, unsigned int nf, unsigned int pin, qaws_scalar* out)
{
	double C[8][3 * 64], G[8][9], y[8];
	unsigned int f, q, i, j, k, a;
	memset(C, 0, sizeof(C));
	for (i = 0; i < nv; i++)
		for (q = 0; q < 3; q++) C[q][3 * i + q] = 1.0 / nv;
	for (q = 0; q < 3; q++) C[3 + q][3 * pin + q] = 1;
	for (f = 0; f < nf; f++)
		for (a = 0; a < 3; a++)
		{
			/* dA/dx_a = n x (x_c - x_b) / 2, n the unit normal; dV/dx_a = x_b x x_c / 6 */
			unsigned int va = tris[3 * f + a], vb = tris[3 * f + (a + 1) % 3], vc = tris[3 * f + (a + 2) % 3];
			double p0[3], p1[3], p2[3], u[3], v[3], n[3], e[3], len;
			for (q = 0; q < 3; q++) { p0[q] = x[3 * tris[3 * f] + q]; p1[q] = x[3 * tris[3 * f + 1] + q]; p2[q] = x[3 * tris[3 * f + 2] + q]; }
			for (q = 0; q < 3; q++) { u[q] = p1[q] - p0[q]; v[q] = p2[q] - p0[q]; e[q] = (double)x[3 * vc + q] - x[3 * vb + q]; }
			n[0] = u[1] * v[2] - u[2] * v[1]; n[1] = u[2] * v[0] - u[0] * v[2]; n[2] = u[0] * v[1] - u[1] * v[0];
			len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
			for (q = 0; q < 3; q++) n[q] /= len;
			C[6][3 * va] += 0.5 * (n[1] * e[2] - n[2] * e[1]); C[6][3 * va + 1] += 0.5 * (n[2] * e[0] - n[0] * e[2]); C[6][3 * va + 2] += 0.5 * (n[0] * e[1] - n[1] * e[0]);
			{
				double xb[3], xc[3];
				for (q = 0; q < 3; q++) { xb[q] = x[3 * vb + q]; xc[q] = x[3 * vc + q]; }
				C[7][3 * va] += (xb[1] * xc[2] - xb[2] * xc[1]) / 6; C[7][3 * va + 1] += (xb[2] * xc[0] - xb[0] * xc[2]) / 6; C[7][3 * va + 2] += (xb[0] * xc[1] - xb[1] * xc[0]) / 6;
			}
		}
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			G[i][j] = 0;
			for (k = 0; k < 3 * nv; k++) G[i][j] += C[i][k] * C[j][k];
		}
		G[i][8] = i == 6;
	}
	for (k = 0; k < 8; k++)
	{
		unsigned int p = k;
		for (i = k + 1; i < 8; i++) if (fabs(G[i][k]) > fabs(G[p][k])) p = i;
		for (j = 0; j < 9; j++) { double t = G[k][j]; G[k][j] = G[p][j]; G[p][j] = t; }
		for (i = k + 1; i < 8; i++)
		{
			double m = G[i][k] / G[k][k];
			for (j = k; j < 9; j++) G[i][j] -= m * G[k][j];
		}
	}
	for (k = 8; k-- > 0;)
	{
		y[k] = G[k][8];
		for (j = k + 1; j < 8; j++) y[k] -= G[k][j] * y[j];
		y[k] /= G[k][k];
	}
	for (k = 0; k < 3 * nv; k++)
	{
		double s = 0;
		for (i = 0; i < 8; i++) s += C[i][k] * y[i];
		out[k] = (qaws_scalar)s;
	}
}

/* a mesh solution's VJP: pins, barycenter, area and volume held */
static void test_mesh_vjp(void)
{
	enum { NV = 6, NF = 8 };
	qaws_scalar x[3 * NV], xp[3 * NV], xm[3 * NV], dx[3 * NV], xbar[3 * NV], dpin[3], pb[3];
	unsigned int tris[3 * NF], pins[1] = { 4 }, i, q, k, it, gauge = 0, resolved = 0;
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[24] = { 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5 };
	qaws_repulsive_mesh_desc md;
	qaws_repulsive_mesh* m = NULL;
	qaws_repulsive_mesh_solve_desc d;
	qaws_repulsive_mesh_solver *s = NULL, *sp = NULL, *sm = NULL;
	qaws_repulsive_report r, rp, rm;
	qaws_repulsive_targets t, tb;
	qaws_status st = QAWS_STATUS_OK, sj, sv;
	double h = RC_F32 ? 2e-2 : 4e-3, worst = 0, fd[3 * NV], fdk[2][3 * NV], an[3 * NV], err, A, V, orth;
	char msg[400];
	for (i = 0; i < NV; i++)
		for (q = 0; q < 3; q++) x[3 * i + q] = (qaws_scalar)(ov[i][q] * (1 + 0.15 * sin(2.0 * i + q)));
	memcpy(tris, of, sizeof(of));
	memset(&md, 0, sizeof(md));
	md.vertex_count = NV; md.triangles = tris; md.triangle_count = NF;
	qaws_repulsive_mesh_create(&m, &md);
	memset(&d, 0, sizeof(d));
	d.pinned = pins; d.pinned_count = 1; d.hold_barycenter = 1; d.hold_area = 1; d.hold_volume = 1;
	d.tolerance = (qaws_scalar)(RC_F32 ? 1e-4 : 1e-6);
	memset(&r, 0, sizeof(r));
	if (qaws_repulsive_mesh_solver_create(&s, m, x, &d) == QAWS_STATUS_OK)
		for (it = 0; it < 3000 && r.reason == QAWS_REPULSIVE_RUNNING; it++) qaws_repulsive_mesh_solver_step(&r, x, s);
	/* the adjoint identity <VJP(xbar), dt> = <xbar, JVP(dt)> for random xbar and dt over every target */
	for (k = 0; k < 3; k++)
	{
		double lhs = 0, rhs = 0;
		memset(&t, 0, sizeof(t));
		for (q = 0; q < 3; q++)
		{
			dpin[q] = (qaws_scalar)sin(1.7 * k + 0.9 * q + 0.3);
			t.barycenter[q] = (qaws_scalar)cos(2.3 * k + 1.1 * q);
		}
		t.pinned = dpin;
		t.area = (qaws_scalar)sin(3.1 * k + 0.7); t.volume = (qaws_scalar)cos(1.9 * k + 0.2);
		for (i = 0; i < 3 * NV; i++) xbar[i] = (qaws_scalar)sin(0.77 * i + 1.3 * k + 0.1);
		sj = s ? qaws_repulsive_mesh_solver_solution_jvp(dx, &gauge, s, &t) : QAWS_STATUS_INVALID_ARGUMENT;
		memset(&tb, 0, sizeof(tb));
		tb.pinned = pb;
		sv = s ? qaws_repulsive_mesh_solver_solution_vjp(&tb, NULL, s, xbar) : QAWS_STATUS_INVALID_ARGUMENT;
		if (sj != QAWS_STATUS_OK) st = sj;
		if (sv != QAWS_STATUS_OK) st = sv;
		for (i = 0; i < 3 * NV; i++) lhs += (double)xbar[i] * dx[i];
		for (q = 0; q < 3; q++) rhs += (double)pb[q] * dpin[q] + (double)tb.barycenter[q] * t.barycenter[q];
		rhs += (double)tb.area * t.area + (double)tb.volume * t.volume;
		worst = fmax(worst, fabs(lhs - rhs) / fmax(fabs(lhs), 1e-300));
	}
	sprintf(msg, "mesh VJP (a pin, barycenter, area and volume held): converged %d, gauge %u (the rotation about the pin-barycenter line), "
		"the adjoint identity over 3 random pairs, worst relative gap %.3g", (int)r.reason, gauge, worst);
	/* (both sides come from one KKT matrix solved to 1e-11: the gap is GMRES's residual times the system's condition) */
	TEST_ASSERT(st == QAWS_STATUS_OK && r.reason == QAWS_REPULSIVE_CONVERGED && gauge == 1 && worst <= (RC_F32 ? 1e-3 : 1e-8), msg);
	printf("  %s\n", msg);
	/* the JVP for the area target against central differences of re-solved meshes. Each re-solve starts at x +- h
	   delta, delta the least-norm move changing the area target alone to first order (C delta = e_area, C the rows'
	   Jacobian: barycenter, pin, area, volume); its other targets then differ from the solution's by O(h^2), alike
	   for + and -, so they cancel in the central difference. (Not set_targets at x itself: there the direction is
	   about 0, so the line search stalls; nor from the far input: the two re-solves then drift apart along the free
	   rotation by more than the difference resolves.) The solution bends sharply with the area here: measured, the
	   central difference's error falls as h^2, from 0.14 at h = 2e-2 to 3e-3 at 2e-3. So the differences at h and
	   h / 2 are combined by Richardson extrapolation, (4 D(h / 2) - D(h)) / 3, leaving O(h^4) and the re-solves' own
	   error (they stop at stationarity 1e-6) over h. Measured 6e-4; the bound 2e-3 as in test_mesh_solution. */
	rc_measures(x, tris, NF, &A, &V);
	rc_area_move(x, tris, NV, NF, pins[0], dx);
	for (k = 0; k < 2; k++)
	{
		double hk = k ? h / 2 : h;
		for (i = 0; i < 3 * NV; i++) { xp[i] = (qaws_scalar)((double)x[i] + hk * dx[i]); xm[i] = (qaws_scalar)((double)x[i] - hk * dx[i]); }
		memset(&rp, 0, sizeof(rp)); memset(&rm, 0, sizeof(rm));
		if (qaws_repulsive_mesh_solver_create(&sp, m, xp, &d) == QAWS_STATUS_OK && qaws_repulsive_mesh_solver_create(&sm, m, xm, &d) == QAWS_STATUS_OK)
		{
			for (it = 0; it < 3000 && rp.reason == QAWS_REPULSIVE_RUNNING; it++) qaws_repulsive_mesh_solver_step(&rp, xp, sp);
			for (it = 0; it < 3000 && rm.reason == QAWS_REPULSIVE_RUNNING; it++) qaws_repulsive_mesh_solver_step(&rm, xm, sm);
		}
		resolved += rp.reason == QAWS_REPULSIVE_CONVERGED && rm.reason == QAWS_REPULSIVE_CONVERGED;
		for (i = 0; i < 3 * NV; i++) fdk[k][i] = ((double)xp[i] - xm[i]) / (2 * hk);
		qaws_repulsive_mesh_solver_destroy(sp); qaws_repulsive_mesh_solver_destroy(sm);
		sp = sm = NULL;
	}
	memset(&t, 0, sizeof(t));
	t.area = 1;
	sj = s ? qaws_repulsive_mesh_solver_solution_jvp(dx, NULL, s, &t) : QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * NV; i++) { fd[i] = (4 * fdk[1][i] - fdk[0][i]) / 3; an[i] = (double)dx[i]; }
	{
		/* (less the rotation about the pin-barycenter line, free between the re-solves; the JVP is orthogonal to it) */
		double b[3] = { 0, 0, 0 }, a[3], R[3 * NV], rr = 0, fr = 0, ar = 0;
		for (i = 0; i < NV; i++)
			for (q = 0; q < 3; q++) b[q] += (double)x[3 * i + q] / NV;
		for (q = 0; q < 3; q++) a[q] = b[q] - x[3 * pins[0] + q];
		for (i = 0; i < NV; i++)
		{
			double p[3];
			for (q = 0; q < 3; q++) p[q] = (double)x[3 * i + q] - x[3 * pins[0] + q];
			R[3 * i] = a[1] * p[2] - a[2] * p[1]; R[3 * i + 1] = a[2] * p[0] - a[0] * p[2]; R[3 * i + 2] = a[0] * p[1] - a[1] * p[0];
		}
		for (i = 0; i < 3 * NV; i++) { rr += R[i] * R[i]; fr += fd[i] * R[i]; ar += an[i] * R[i]; }
		for (i = 0; i < 3 * NV; i++) fd[i] -= fr / rr * R[i];
		orth = fabs(ar) / sqrt(rr);
	}
	err = rc_maxrel(an, fd, 3 * NV);
	sprintf(msg, "mesh JVP, the area target (A %.4g) moved by +-%.0e and +-%.0e: %u of 2 re-solved pairs converged, the JVP's part along "
		"the free rotation %.2g, against the extrapolated central difference less it %.3g", A, h, h / 2, resolved, orth, err);
	/* float: the re-solves stop at float's tolerances, and their difference over h carries that stop error; measured 0.17, so
	   in float this is a check of sign and size only (0.3); the gauge removal leaves rounding along the rotation, a few hundred float
	   eps (measured 1.3e-7) against 1e-8 in double; the derivative's precision being checked by the adjoint identity
	   above and by the double build */
	TEST_ASSERT(sj == QAWS_STATUS_OK && resolved == 2 && orth <= (RC_F32 ? 1e-5 : 1e-8) && err <= (RC_F32 ? 3e-1 : 2e-3), msg);
	printf("  %s\n", msg);
	qaws_repulsive_mesh_solver_destroy(s); qaws_repulsive_mesh_solver_destroy(sp); qaws_repulsive_mesh_solver_destroy(sm);
	qaws_repulsive_mesh_destroy(m);
}

int test_118_repulsive_controls_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 118: native control points and solution derivatives\n");
	test_controls();
	test_curve_solution();
	test_loop_solution();
	test_mesh_solution();
	test_derivative_report();
	test_mesh_vjp();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
