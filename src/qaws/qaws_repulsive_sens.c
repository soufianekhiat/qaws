/*
 * Derivatives of a converged constrained solution (qaws_repulsive.h, R5).
 *
 * At a converged x*, dE = C^T lambda and c(x*) = theta (each constraint row
 * held at its target). Moving the targets by dtheta moves the solution by
 * dx, from the KKT system of the Lagrangian E - lambda . c:
 *   [K C^T; C 0] [dx; mu] = [0; dtheta],  K = Hess E - sum_r lambda_r Hess c_r.
 * lambda is C's least-squares multiplier of dE. K is applied matrix-free
 * (the energy's HVP and the constraints' Hessians). The system is solved by
 * flexible GMRES, preconditioned by the descent's own dense system
 * [A~ C^T; C 0] (the Sobolev metric standing in for K: the same principal
 * part). The system is symmetric, so the adjoint (VJP) is the same solve
 * with [xbar; 0]: theta_bar = its multiplier part.
 *
 * Rank and nullspaces:
 *   - a rotation about the barycenter that every constraint row allows (no
 *     pins, no held directions) is a symmetry of the problem. The solution
 *     is then not isolated, so K is singular along it. Such rotations are
 *     detected (C R = 0) and gauged out: dx is orthogonal to them, and
 *     their count is returned;
 *   - a free translation (nothing pinned, no barycenter) is refused, as by
 *     the descent;
 *   - a singular preconditioner, or a solve that does not converge
 *     (dependent rows, or not a strict minimum on the constraint set), is a
 *     NUMERICAL_FAILURE;
 *   - a state not converged (the solver's last reason) is NOT_CONVERGED.
 *     The derivative of a fixed sequence of steps is not provided.
 */

#include "qaws_repulsive.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RSN_RESTART 80       /* the defaults of the desc's gmres_* settings */
#define RSN_MAXIT 2400
#define RSN_TOL 1e-11
#define RSN_ACCEPT 1e-7   /* the residual accepted once restarts stop gaining (the floor of rounding) */

typedef qaws_internal_rs_core rs_core;

typedef struct rsn_sys
{
	rs_core const* s;
	unsigned int N, k, ng, n;       /* n = N + k + ng */
	double *C, *G, *lambda, *M, *t1, *t2;
	unsigned int* piv;
} rsn_sys;

static double rsn_dot(double const* a, double const* b, unsigned int n)
{
	double r = 0;
	unsigned int i;
	for (i = 0; i < n; i++) r += a[i] * b[i];
	return r;
}

/* the symmetric m x m A diagonalized in place by cyclic Jacobi rotations,
   V (the identity on entry) its eigenvectors as columns */
static void rsn_jacobi(double* A, double* V, unsigned int m)
{
	unsigned int sweep, p, q, i;
	for (sweep = 0; sweep < 50; sweep++)
	{
		double off = 0;
		for (p = 0; p < m; p++)
			for (q = p + 1; q < m; q++) off += A[p * m + q] * A[p * m + q];
		if (off <= 1e-60)
			break;
		for (p = 0; p < m; p++)
			for (q = p + 1; q < m; q++)
			{
				double apq = A[p * m + q], th, t, c, s;
				if (fabs(apq) < 1e-300)
					continue;
				th = (A[q * m + q] - A[p * m + p]) / (2 * apq);
				t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1));
				c = 1 / sqrt(t * t + 1); s = t * c;
				for (i = 0; i < m; i++)
				{
					double aip = A[i * m + p], aiq = A[i * m + q];
					A[i * m + p] = c * aip - s * aiq; A[i * m + q] = s * aip + c * aiq;
				}
				for (i = 0; i < m; i++)
				{
					double api = A[p * m + i], aqi = A[q * m + i];
					A[p * m + i] = c * api - s * aqi; A[q * m + i] = s * api + c * aqi;
				}
				for (i = 0; i < m; i++)
				{
					double vip = V[i * m + p], viq = V[i * m + q];
					V[i * m + p] = c * vip - s * viq; V[i * m + q] = s * vip + c * viq;
				}
			}
	}
}

/* y = [K dx + C^T mu + G^T g; C dx; G dx] */
static qaws_status rsn_apply(rsn_sys* S, double const* v, double* y)
{
	unsigned int N = S->N, k = S->k, i, r;
	qaws_status st;
	memset(y, 0, sizeof(double) * S->n);
	if ((st = S->s->ops->hvp(S->s, S->s->x, v, y)) != QAWS_STATUS_OK)
		return st;
	memset(S->t1, 0, sizeof(double) * N);
	S->s->ops->constraints_hvp(S->s, S->s->x, S->lambda, v, S->t1);
	for (i = 0; i < N; i++) y[i] -= S->t1[i];
	for (r = 0; r < k; r++)
	{
		double const* c = S->C + (size_t)r * N;
		for (i = 0; i < N; i++) y[i] += c[i] * v[N + r];
		y[N + r] = rsn_dot(c, v, N);
	}
	for (r = 0; r < S->ng; r++)
	{
		double const* g = S->G + (size_t)r * N;
		for (i = 0; i < N; i++) y[i] += g[i] * v[N + k + r];
		y[N + k + r] = rsn_dot(g, v, N);
	}
	return QAWS_STATUS_OK;
}

/* S->M = [A~ C^T G^T; C 0 0; G 0 0], factored */
static qaws_status rsn_precond(rsn_sys* S, double const* A)
{
	unsigned int N = S->N, n = S->n, dim = S->s->dim, nv = S->s->nv, i, j, q, r;
	double* rows;
	memset(S->M, 0, sizeof(double) * n * n);
	for (i = 0; i < nv; i++)
		for (j = 0; j < nv; j++)
			for (q = 0; q < dim; q++) S->M[(dim * i + q) * n + dim * j + q] = A[(size_t)i * nv + j];
	for (r = 0; r < S->k + S->ng; r++)
	{
		rows = r < S->k ? S->C + (size_t)r * N : S->G + (size_t)(r - S->k) * N;
		for (i = 0; i < N; i++) { S->M[(N + r) * n + i] = rows[i]; S->M[i * n + N + r] = rows[i]; }
	}
	return qaws_internal_rs_lu(S->M, S->piv, n) ? QAWS_STATUS_OK : QAWS_STATUS_NUMERICAL_FAILURE;
}

/* the ratio of the largest to the smallest singular value of the m x m
   upper triangle R (row stride ld): the square roots of R^T R's extreme
   eigenvalues (HUGE_VAL when singular; 1 when out of memory: no estimate) */
static double rsn_cond(double const* R, size_t ld, unsigned int m)
{
	double *A = (double*)qaws_internal_scratch_zalloc(qaws_internal_mul3(2, m, m), sizeof(double)), *V, lo = HUGE_VAL, hi = 0;
	unsigned int i, j, k;
	if (!A)
		return 1;
	V = A + (size_t)m * m;
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
		{
			double sum = 0;
			for (k = 0; k <= (i < j ? i : j); k++) sum += R[k * ld + i] * R[k * ld + j];
			A[(size_t)i * m + j] = sum;
			V[(size_t)i * m + j] = i == j;
		}
	rsn_jacobi(A, V, m);
	for (i = 0; i < m; i++) { lo = fmin(lo, A[(size_t)i * m + i]); hi = fmax(hi, A[(size_t)i * m + i]); }
	qaws_internal_free(A);
	return lo > 0 ? fmax(sqrt(hi / lo), 1) : HUGE_VAL;
}

/* the system's solve by flexible GMRES, right-preconditioned, within the
   desc's gmres_* limits; rep (may be NULL): what it achieved */
static qaws_status rsn_solve(rsn_sys* S, double const* b, double* u, qaws_repulsive_derivative_report* rep)
{
	rs_core const* s = S->s;
	size_t R = s->gmres_restart ? s->gmres_restart : RSN_RESTART;
	unsigned int n = S->n, i, j, kk, outer, it = 0;
	unsigned int maxit = s->gmres_max ? s->gmres_max : RSN_MAXIT, maxouter = (unsigned int)(maxit / R + (maxit % R != 0));
	double *V = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, (R + 1))), *Z = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), n, (R + 1))), *H = (double*)qaws_internal_zalloc(qaws_internal_mul2((R + 1), R), sizeof(double));
	double *w = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), n)), *g = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), qaws_internal_add2(qaws_internal_mul2(4, R), 1)));
	double *cs = NULL, *sn = NULL, *gq = NULL, *y = NULL, bn = sqrt(rsn_dot(b, b, n)), prev = HUGE_VAL, beta = 0, cond = 1, tol = s->gmres_tol > 0 ? s->gmres_tol : RSN_TOL;
	qaws_status st = QAWS_STATUS_NUMERICAL_FAILURE;
	memset(u, 0, sizeof(double) * n);
	if (!V || !Z || !H || !w || !g)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	cs = g; sn = g + R; gq = g + 2 * R; y = g + 3 * R + 1;
	if (!(bn > 0))
	{
		st = QAWS_STATUS_OK;
		goto done;
	}
	for (outer = 0;; outer++)   /* (the defaults: 30 restarts of 80, then a last residual) */
	{
		qaws_status a;
		if ((a = rsn_apply(S, u, w)) != QAWS_STATUS_OK) { st = a; goto done; }
		for (i = 0; i < n; i++) w[i] = b[i] - w[i];
		beta = sqrt(rsn_dot(w, w, n));
		if (beta <= tol * bn || (beta <= RSN_ACCEPT * bn && beta > 0.1 * prev))
		{
			st = QAWS_STATUS_OK;
			break;
		}
		if (outer == maxouter)
			break;
		prev = beta;
		for (i = 0; i < n; i++) V[i] = w[i] / beta;
		memset(gq, 0, sizeof(double) * (R + 1));
		gq[0] = beta;
		for (j = 0; j < R; j++)
		{
			double* z = Z + (size_t)j * n;
			double hn;
			memcpy(z, V + (size_t)j * n, sizeof(double) * n);
			qaws_internal_rs_lu_solve(z, S->M, S->piv, n);
			if ((a = rsn_apply(S, z, w)) != QAWS_STATUS_OK) { st = a; goto done; }
			for (kk = 0; kk <= j; kk++)
			{
				double hk = rsn_dot(w, V + (size_t)kk * n, n);
				H[kk * R + j] = hk;
				for (i = 0; i < n; i++) w[i] -= hk * V[(size_t)kk * n + i];
			}
			hn = sqrt(rsn_dot(w, w, n));
			H[(j + 1) * R + j] = hn;
			if (hn > 0)
				for (i = 0; i < n; i++) V[(size_t)(j + 1) * n + i] = w[i] / hn;
			for (kk = 0; kk < j; kk++)
			{
				double a0 = H[kk * R + j], a1 = H[(kk + 1) * R + j];
				H[kk * R + j] = cs[kk] * a0 + sn[kk] * a1;
				H[(kk + 1) * R + j] = -sn[kk] * a0 + cs[kk] * a1;
			}
			{
				double a0 = H[j * R + j], a1 = H[(j + 1) * R + j], rr = hypot(a0, a1);
				cs[j] = rr > 0 ? a0 / rr : 1; sn[j] = rr > 0 ? a1 / rr : 0;
				H[j * R + j] = rr; H[(j + 1) * R + j] = 0;
				gq[j + 1] = -sn[j] * gq[j]; gq[j] *= cs[j];
			}
			it++;
			if (fabs(gq[j + 1]) <= tol * bn || j + 1 == R || !(hn > 0) || it == maxit)
			{
				unsigned int m = j + 1;
				for (kk = m; kk-- > 0;)
				{
					double sm = gq[kk];
					for (i = kk + 1; i < m; i++) sm -= H[kk * R + i] * y[i];
					y[kk] = H[kk * R + kk] != 0 ? sm / H[kk * R + kk] : 0;
				}
				for (kk = 0; kk < m; kk++)
					for (i = 0; i < n; i++) u[i] += y[kk] * Z[(size_t)kk * n + i];
				if (rep)
					cond = fmax(cond, rsn_cond(H, R, m));   /* (H is R of its QR now: the same singular values) */
				break;
			}
		}
	}
done:
	if (rep && st != QAWS_STATUS_ALLOCATION_FAILURE)
	{
		rep->residual = (qaws_scalar)(bn > 0 ? beta / bn : 0);
		rep->iterations = it;
		rep->converged = st == QAWS_STATUS_OK;
		rep->condition_estimate = (qaws_scalar)cond;
	}
	qaws_internal_free(V); qaws_internal_free(Z); qaws_internal_free(H); qaws_internal_free(w); qaws_internal_free(g);
	return st;
}

qaws_status qaws_internal_rs_sensitivity(double* out, unsigned int* out_gauge, qaws_repulsive_derivative_report* out_report, rs_core const* s, double const* drows,
	double const* xbar)
{
	rsn_sys S;
	unsigned int N, k, i, r, q, dim, nv;
	double *A = NULL, *dE = NULL, *phi = NULL, *rhs = NULL, *u = NULL, *CC = NULL, b[3] = { 0, 0, 0 }, E, cn = 0;
	unsigned int* cpiv = NULL;
	qaws_status st;
	if (out_gauge) *out_gauge = 0;
	if (!s || (!drows && !xbar) || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (s->theta > 0 || s->full || !s->ops->hvp || !s->ops->constraints_hvp)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (dense x-space solvers only) */
	if (s->rep.reason != QAWS_REPULSIVE_CONVERGED)
		return QAWS_STATUS_NOT_CONVERGED;
	if (s->translates)
		return QAWS_STATUS_NUMERICAL_FAILURE;
	N = s->N; k = s->k; dim = s->dim; nv = s->nv;
	memset(&S, 0, sizeof(S));
	S.s = s; S.N = N; S.k = k;
	A = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)nv * nv + 1)));
	dE = (double*)qaws_internal_scratch_zalloc(N + 1, sizeof(double));
	phi = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (k + 1)));
	S.C = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)k * N + 1)));
	S.G = (double*)qaws_internal_zalloc(6 * (size_t)N + 1, sizeof(double));
	S.lambda = (double*)qaws_internal_zalloc(k + 1, sizeof(double));
	S.t1 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (N + 1)));
	CC = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (k * k + 1)));
	cpiv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (k + 1)));
	if (!A || !dE || !phi || !S.C || !S.G || !S.lambda || !S.t1 || !CC || !cpiv)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if ((st = s->ops->energy(s, s->x, &E, dE)) != QAWS_STATUS_OK || (st = s->ops->metric(s, s->x, A)) != QAWS_STATUS_OK)
		goto done;
	s->ops->constraints(s, s->x, phi, S.C);
	/* lambda: (C C^T) lambda = C dE */
	if (k)
	{
		for (r = 0; r < k; r++)
		{
			for (q = 0; q < k; q++) CC[r * k + q] = rsn_dot(S.C + (size_t)r * N, S.C + (size_t)q * N, N);
			S.lambda[r] = rsn_dot(S.C + (size_t)r * N, dE, N);
			cn = fmax(cn, sqrt(CC[r * k + r]));
		}
		if (!qaws_internal_rs_lu(CC, cpiv, k))
		{
			st = QAWS_STATUS_NUMERICAL_FAILURE;   /* (dependent constraint rows) */
			goto done;
		}
		qaws_internal_rs_lu_solve(S.lambda, CC, cpiv, k);
	}
	/* the rigid motions every row allows: of the translations and the
	   rotations about the barycenter (unit vectors M_a), the combinations
	   with C M c = 0 (the rows normalized): the eigenvectors of
	   (C M)^T (C M) whose eigenvalue is about 0. A strand pinned at both
	   ends keeps the rotation about the line through them, a free loop all
	   three about its barycenter. */
	for (i = 0; i < nv; i++)
		for (q = 0; q < dim; q++) b[q] += s->x[dim * i + q] / nv;
	{
		unsigned int m = dim == 3 ? 6 : 3, a, c2, j;
		double *Mo = (double*)qaws_internal_scratch_zalloc(qaws_internal_add2(qaws_internal_mul2(m, N), 1), sizeof(double)), Gm[36], V[36], ev[6];
		if (!Mo)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (a = 0; a < m; a++)
		{
			double* g = Mo + (size_t)a * N, gn;
			for (i = 0; i < nv; i++)
			{
				double d[3] = { 0, 0, 0 };
				unsigned int c;
				for (c = 0; c < dim; c++) d[c] = s->x[dim * i + c] - b[c];
				if (a < dim)
					g[dim * i + a] = 1;
				else if (dim == 2) { g[2 * i] = -d[1]; g[2 * i + 1] = d[0]; }
				else
				{
					unsigned int r3 = a - 3;
					g[3 * i + (r3 + 1) % 3] = -d[(r3 + 2) % 3]; g[3 * i + (r3 + 2) % 3] = d[(r3 + 1) % 3];
				}
			}
			gn = sqrt(rsn_dot(g, g, N));
			for (i = 0; gn > 0 && i < N; i++) g[i] /= gn;
		}
		for (a = 0; a < m; a++)
			for (c2 = 0; c2 < m; c2++)
			{
				double sum = 0;
				for (r = 0; r < k; r++)
				{
					double const* cr = S.C + (size_t)r * N;
					double nr = rsn_dot(cr, cr, N);
					if (nr > 0) sum += rsn_dot(cr, Mo + (size_t)a * N, N) * rsn_dot(cr, Mo + (size_t)c2 * N, N) / nr;
				}
				Gm[a * m + c2] = sum;
				V[a * m + c2] = a == c2;
			}
		rsn_jacobi(Gm, V, m);
		for (a = 0; a < m; a++) ev[a] = Gm[a * m + a];
		for (a = 0; a < m; a++)
		{
			double* g = S.G + (size_t)S.ng * N, gn;
			if (!(ev[a] <= 1e-16))
				continue;   /* (|C M c| above 1e-8: a constraint holds this motion) */
			memset(g, 0, sizeof(double) * N);
			for (j = 0; j < m; j++)
				for (i = 0; i < N; i++) g[i] += V[j * m + a] * Mo[(size_t)j * N + i];
			gn = sqrt(rsn_dot(g, g, N));
			for (i = 0; gn > 0 && i < N; i++) g[i] /= gn;
			S.ng++;
		}
		qaws_internal_free(Mo);
	}
	(void)cn;
	S.n = N + k + S.ng;
	S.M = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)S.n * S.n + 1)));
	S.piv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (S.n + 1)));
	rhs = (double*)qaws_internal_scratch_zalloc(S.n + 1, sizeof(double));
	u = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (S.n + 1)));
	if (!S.M || !S.piv || !rhs || !u)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if ((st = rsn_precond(&S, A)) != QAWS_STATUS_OK)
		goto done;
	if (drows)
		for (r = 0; r < k; r++) rhs[N + r] = drows[r];
	else
		for (i = 0; i < N; i++) rhs[i] = xbar[i];
	if (out_report) memset(out_report, 0, sizeof(*out_report));
	st = rsn_solve(&S, rhs, u, out_report);
	if (out_report && (st == QAWS_STATUS_OK || st == QAWS_STATUS_NUMERICAL_FAILURE)) out_report->gauge = S.ng;
	if (st != QAWS_STATUS_OK)
		goto done;
	if (drows)
		for (i = 0; i < N; i++) out[i] = u[i];
	else
		for (r = 0; r < k; r++) out[r] = u[N + r];
	if (out_gauge) *out_gauge = S.ng;
done:
	qaws_internal_free(A); qaws_internal_free(dE); qaws_internal_free(phi); qaws_internal_free(S.C); qaws_internal_free(S.G); qaws_internal_free(S.lambda); qaws_internal_free(S.t1); qaws_internal_free(S.M); qaws_internal_free(S.piv); qaws_internal_free(rhs); qaws_internal_free(u); qaws_internal_free(CC); qaws_internal_free(cpiv);
	return st;
}
