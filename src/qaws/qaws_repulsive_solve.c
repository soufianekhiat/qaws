/*
 * Fractional Sobolev descent of tangent-point energies (qaws_repulsive.h;
 * Repulsive Curves, Sec. 4.2 and 5.2 to 5.4; Repulsive Surfaces, Sec. 3.2
 * and 3.3; the contracts in docs/repulsive_plan.md).
 *
 * A core shared by curves and meshes (internal/qaws_internal_repulsive_solve.h):
 * the direction from the constrained system [A~ C^T; C 0], the backtracking
 * step, the Newton projection back onto the constraints, the termination
 * reasons, and conservative advancement for the collision checks. Then the
 * curves' metric, constraints and checks. The dense mode: the metric is
 * |V| x |V|, the systems (dim |V| + k)^2, solved by LU with partial
 * pivoting; the hierarchical mode (theta > 0) is qaws_repulsive_hsolve.c.
 * Everything in double; a solver keeps its own double copy of the
 * positions.
 */

#include "qaws_repulsive.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if QAWS_SCALAR_IS_FLOAT
#define RS_CTOL 1e-5
#else
#define RS_CTOL 1e-8
#endif
#define RS_ARMIJO 1e-4
#define RS_TRIALS 40

typedef qaws_internal_rs_core rs_core;

/* ------------------------------------------------------------------------ */
/*  Small helpers                                                           */
/* ------------------------------------------------------------------------ */

static double rs_norm(double const* a, unsigned int dim)
{
	double s = 0;
	unsigned int c;
	for (c = 0; c < dim; c++) s += a[c] * a[c];
	return sqrt(s);
}

static int rs_shared(unsigned int const* a, unsigned int const* b)
{
	return a[0] == b[0] || a[0] == b[1] || a[1] == b[0] || a[1] == b[1];
}

/* LU with partial pivoting of the n x n row-major m, in place; 0 when a
   pivot is negligible against the matrix's largest entry (singular) */
int qaws_internal_rs_lu(double* m, unsigned int n, unsigned int* piv)
{
	unsigned int i, j, k;
	double big = 0;
	for (i = 0; i < n * n; i++) big = fmax(big, fabs(m[i]));
	for (k = 0; k < n; k++)
	{
		unsigned int p = k;
		double best = fabs(m[k * n + k]);
		for (i = k + 1; i < n; i++)
			if (fabs(m[i * n + k]) > best) { best = fabs(m[i * n + k]); p = i; }
		if (!(best > 1e-13 * big))
			return 0;
		piv[k] = p;
		if (p != k)
			for (j = 0; j < n; j++) { double t = m[k * n + j]; m[k * n + j] = m[p * n + j]; m[p * n + j] = t; }
		for (i = k + 1; i < n; i++)
		{
			double f = m[i * n + k] / m[k * n + k];
			m[i * n + k] = f;
			if (f != 0)
				for (j = k + 1; j < n; j++) m[i * n + j] -= f * m[k * n + j];
		}
	}
	return 1;
}

void qaws_internal_rs_lu_solve(double const* m, unsigned int n, unsigned int const* piv, double* b)
{
	unsigned int i, j, k;
	for (k = 0; k < n; k++)
		if (piv[k] != k) { double t = b[k]; b[k] = b[piv[k]]; b[piv[k]] = t; }
	for (i = 1; i < n; i++)
		for (j = 0; j < i; j++) b[i] -= m[i * n + j] * b[j];
	for (i = n; i-- > 0;)
	{
		for (j = i + 1; j < n; j++) b[i] -= m[i * n + j] * b[j];
		b[i] /= m[i * n + i];
	}
}

/* ------------------------------------------------------------------------ */
/*  The core                                                                */
/* ------------------------------------------------------------------------ */

static double rs_phi_max(rs_core const* s, double const* phi)
{
	double m = 0;
	unsigned int r;
	for (r = 0; r < s->k; r++) m = fmax(m, fabs(phi[r]));
	return m;
}

/* the metric's storage: nv x nv (A~ = A on each coordinate), or N x N */
static size_t rs_asize(rs_core const* s)
{
	return s->full ? (size_t)s->N * s->N : (size_t)s->nv * s->nv;
}

/* how far vertex v moves (in the embedding) for the change d at x */
static double rs_vmove(rs_core const* s, double const* x, double const* d, unsigned int v)
{
	return s->ops->move ? s->ops->move(s, x, d, v) : rs_norm(d + s->dim * v, s->dim);
}

/* the saddle-point matrix [A~ C^T; C 0] (A~ = A on each coordinate),
   factored; 0 when singular */
static int rs_kkt(rs_core const* s, double const* A, double const* C, double* M, unsigned int* piv)
{
	unsigned int dim = s->dim, N = s->N, n = N + s->k, i, j, q, r;
	memset(M, 0, sizeof(double) * n * n);
	if (s->full)
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) M[i * n + j] = A[(size_t)i * N + j];
	else
		for (i = 0; i < s->nv; i++)
			for (j = 0; j < s->nv; j++)
				for (q = 0; q < dim; q++)
					M[(dim * i + q) * n + dim * j + q] = A[i * s->nv + j];
	for (r = 0; r < s->k; r++)
		for (i = 0; i < N; i++)
		{
			M[(N + r) * n + i] = C[r * N + i];
			M[i * n + N + r] = C[r * N + i];
		}
	return qaws_internal_rs_lu(M, n, piv);
}

/* the direction at x: g (N + k: the multipliers after g), with A, dE (N)
   and E; the KKT factors left in M */
static qaws_status rs_direction(rs_core const* s, double const* x, double* A, double* C, double* phi, double* M,
	unsigned int* piv, double* dE, double* E, double* g)
{
	unsigned int i, n = s->N + s->k;
	qaws_status st;
	memset(dE, 0, sizeof(double) * s->N);
	if ((st = s->ops->energy(s, x, E, dE)) != QAWS_STATUS_OK)
		return st;
	if ((st = s->ops->metric(s, x, A)) != QAWS_STATUS_OK)
		return st;
	s->ops->constraints(s, x, phi, C);
	if (!rs_kkt(s, A, C, M, piv))
		return QAWS_STATUS_NUMERICAL_FAILURE;
	for (i = 0; i < n; i++) g[i] = i < s->N ? dE[i] : 0;
	qaws_internal_rs_lu_solve(M, n, piv, g);
	return QAWS_STATUS_OK;
}

qaws_status qaws_internal_rs_direction(rs_core const* s, qaws_scalar const* positions, qaws_scalar* out_direction)
{
	unsigned int n, i;
	double *buf, *x, *A, *C, *phi, *M, *dE, *g, E;
	unsigned int* piv;
	qaws_status st;
	if (!s || !positions || !out_direction)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (s->theta > 0)
	{
		qaws_internal_rs_hsys* H = NULL;
		buf = (double*)calloc(3 * (size_t)s->N + 1, sizeof(double));
		if (!buf)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		x = buf; dE = x + s->N; g = dE + s->N;
		for (i = 0; i < s->N; i++) x[i] = (double)positions[i];
		st = s->ops->energy(s, x, &E, dE);
		if (st == QAWS_STATUS_OK)
			st = qaws_internal_rs_hdirection(s, x, dE, g, &H);
		for (i = 0; st == QAWS_STATUS_OK && i < s->N; i++) out_direction[i] = (qaws_scalar)g[i];
		qaws_internal_rs_hsys_free(H);
		free(buf);
		return st;
	}
	n = s->N + s->k;
	buf = (double*)malloc(sizeof(double) * ((size_t)s->N * 3 + rs_asize(s) + (size_t)s->k * s->N + s->k + (size_t)n * n + n + 1));
	piv = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	if (!buf || !piv)
	{
		free(buf); free(piv);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	x = buf; dE = x + s->N; A = dE + s->N; C = A + rs_asize(s); phi = C + (size_t)s->k * s->N; M = phi + s->k; g = M + (size_t)n * n;
	for (i = 0; i < s->N; i++) x[i] = (double)positions[i];
	st = rs_direction(s, x, A, C, phi, M, piv, dE, &E, g);
	for (i = 0; st == QAWS_STATUS_OK && i < s->N; i++) out_direction[i] = (qaws_scalar)g[i];
	free(buf); free(piv);
	return st;
}

/* y back onto the constraints: Newton steps [A~ C(y)^T; C(y) 0] [dx; m] =
   [0; -phi(y)], the metric A fixed; 1 when |phi| <= ctol */
static int rs_project(rs_core const* s, double const* A, double* y, double* C, double* phi, double* M, unsigned int* piv, double* rhs)
{
	unsigned int it, i, n = s->N + s->k;
	if (s->k == 0)
		return 1;
	for (it = 0; it < 12; it++)
	{
		s->ops->constraints(s, y, phi, C);
		if (rs_phi_max(s, phi) <= s->ctol)
			return 1;
		if (!rs_kkt(s, A, C, M, piv))
			return 0;
		for (i = 0; i < n; i++) rhs[i] = i < s->N ? 0 : -phi[i - s->N];
		qaws_internal_rs_lu_solve(M, n, piv, rhs);
		for (i = 0; i < s->N; i++) y[i] += rhs[i];
	}
	s->ops->constraints(s, y, phi, NULL);
	return rs_phi_max(s, phi) <= s->ctol;
}

qaws_status qaws_internal_rs_step(rs_core* s, qaws_scalar* positions, qaws_repulsive_report* out_report)
{
	unsigned int n, i, v, trial, same = 1;
	double *buf, *A, *C, *phi, *M, *dE, *g, *y, *rhs, E, Ey = 0, gdot = 0, gmax = 0, h;
	unsigned int* piv;
	qaws_status st;
	int last = 0;   /* the last rejection: 1 energy, 2 collision, 3 projection */
	qaws_internal_rs_hsys* H = NULL;
	size_t dense;
	if (!s || !positions)
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = s->N + s->k;
	dense = s->theta > 0 ? 0 : rs_asize(s) + (size_t)s->k * s->N + (size_t)n * n + n;
	/* the caller's positions, unless they are this solver's own rounded */
	for (i = 0; i < s->N; i++) same &= positions[i] == (qaws_scalar)s->x[i];
	if (!same)
		for (i = 0; i < s->N; i++) s->x[i] = (double)positions[i];
	buf = (double*)malloc(sizeof(double) * ((size_t)s->N * 2 + s->k + n + dense + 1));
	piv = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	if (!buf || !piv)
	{
		free(buf); free(piv);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	dE = buf; y = dE + s->N; phi = y + s->N; g = phi + s->k;
	A = C = M = rhs = NULL;
	if (dense)
	{
		A = g + n; C = A + rs_asize(s); M = C + (size_t)s->k * s->N; rhs = M + (size_t)n * n;
		st = rs_direction(s, s->x, A, C, phi, M, piv, dE, &E, g);
	}
	else
	{
		memset(dE, 0, sizeof(double) * s->N);
		st = s->ops->energy(s, s->x, &E, dE);
		if (st == QAWS_STATUS_OK)
			st = qaws_internal_rs_hdirection(s, s->x, dE, g, &H);
	}
	if (st != QAWS_STATUS_OK)
	{
		s->rep.reason = QAWS_REPULSIVE_NUMERICAL_FAILURE;
		goto done;
	}
	for (i = 0; i < s->N; i++) gdot += dE[i] * g[i];
	for (v = 0; v < s->nv; v++) gmax = fmax(gmax, rs_vmove(s, s->x, g, v));
	s->rep.energy = (qaws_scalar)E;
	s->rep.stationarity = (qaws_scalar)(gmax > 0 && E > 0 ? gdot * s->lmean / (E * gmax) : 0);
	s->ops->constraints(s, s->x, phi, NULL);
	s->rep.constraint_residual = (qaws_scalar)(rs_phi_max(s, phi) / s->lmean);
	if ((double)s->rep.stationarity <= s->tol)
	{
		s->rep.reason = QAWS_REPULSIVE_CONVERGED;
		goto done;
	}
	/* backtracking on the largest vertex move h */
	h = s->h;
	for (trial = 0; trial < RS_TRIALS; trial++, h *= 0.5)
	{
		double tau = h / gmax;
		for (i = 0; i < s->N; i++) y[i] = s->x[i] - tau * g[i];
		if (dense ? !rs_project(s, A, y, C, phi, M, piv, rhs) : !qaws_internal_rs_hproject(H, y))
		{
			s->rep.rejected_projection++;
			last = 3;
			continue;
		}
		if (s->ops->energy(s, y, &Ey, NULL) != QAWS_STATUS_OK || !(Ey <= E - RS_ARMIJO * tau * gdot))
		{
			s->rep.rejected_energy++;
			last = 1;
			continue;
		}
		if (s->checked && !s->ops->path_ok(s, s->x, y))
		{
			s->rep.rejected_collision++;
			last = 2;
			continue;
		}
		break;
	}
	if (trial == RS_TRIALS)
	{
		s->rep.reason = last == 2 ? QAWS_REPULSIVE_COLLISION_LIMITED : last == 3 ? QAWS_REPULSIVE_INFEASIBLE : QAWS_REPULSIVE_STALLED;
		goto done;
	}
	/* accepted: the move, and the next trial twice this one */
	{
		double move = 0;
		for (i = 0; i < s->N; i++) dE[i] = y[i] - s->x[i];   /* (dE no longer needed) */
		for (v = 0; v < s->nv; v++) move = fmax(move, rs_vmove(s, s->x, dE, v));
		s->rep.last_step = (qaws_scalar)(move / s->lmean);
	}
	memcpy(s->x, y, sizeof(double) * s->N);
	s->h = 2 * h;
	s->rep.energy = (qaws_scalar)Ey;
	s->rep.iterations++;
	/* accepted but vanishing: pairs closing in under the collision rule
	   (each step may halve a gap, never close it), or nothing left to gain */
	s->rep.reason = last == 2 && (double)s->rep.last_step < 1e-8 ? QAWS_REPULSIVE_COLLISION_LIMITED :
		(double)s->rep.last_step < 1e-12 ? QAWS_REPULSIVE_STALLED : QAWS_REPULSIVE_RUNNING;
	s->ops->constraints(s, s->x, phi, NULL);
	s->rep.constraint_residual = (qaws_scalar)(rs_phi_max(s, phi) / s->lmean);
	for (i = 0; i < s->N; i++) positions[i] = (qaws_scalar)s->x[i];
done:
	qaws_internal_rs_hsys_free(H);
	free(buf); free(piv);
	if (out_report) *out_report = s->rep;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Collision checks                                                        */
/* ------------------------------------------------------------------------ */

double qaws_internal_rs_segdist(double const* p0, double const* p1, double const* q0, double const* q1, unsigned int dim)
{
	double d1[3] = { 0, 0, 0 }, d2[3] = { 0, 0, 0 }, r[3] = { 0, 0, 0 }, a = 0, e = 0, f = 0, b = 0, c = 0, s, t, den, out = 0;
	unsigned int q;
	for (q = 0; q < dim; q++) { d1[q] = p1[q] - p0[q]; d2[q] = q1[q] - q0[q]; r[q] = p0[q] - q0[q]; }
	for (q = 0; q < dim; q++) { a += d1[q] * d1[q]; e += d2[q] * d2[q]; f += d2[q] * r[q]; b += d1[q] * d2[q]; c += d1[q] * r[q]; }
	if (a <= 0 && e <= 0)
		s = t = 0;
	else if (a <= 0)
	{
		s = 0;
		t = fmin(fmax(f / e, 0), 1);
	}
	else if (e <= 0)
	{
		t = 0;
		s = fmin(fmax(-c / a, 0), 1);
	}
	else
	{
		den = a * e - b * b;
		s = den > 0 ? fmin(fmax((b * f - c * e) / den, 0), 1) : 0;
		t = (b * s + f) / e;
		if (t < 0) { t = 0; s = fmin(fmax(-c / a, 0), 1); }
		else if (t > 1) { t = 1; s = fmin(fmax((b - c) / a, 0), 1); }
	}
	for (q = 0; q < dim; q++)
	{
		double d = p0[q] + s * d1[q] - (q0[q] + t * d2[q]);
		out += d * d;
	}
	return sqrt(out);
}

/* over [t0, t1] the distance is at least (d(t0) + d(t1) - L (t1 - t0)) / 2;
   intervals not settled so are halved, down to 2^-24 (then a collision) */
int qaws_internal_rs_advance(double (*dist)(void const* user, double t), void const* user, double d0, double L, double req)
{
	double st[64][4];
	unsigned int top = 0;
	if (!(req > 0) || d0 < req)
		return 0;
	if (d0 - L >= req)
		return 1;
	st[top][0] = 0; st[top][1] = 1; st[top][2] = d0; st[top][3] = dist(user, 1); top++;
	while (top)
	{
		double t0, t1, a, b, tm, dm;
		top--;
		t0 = st[top][0]; t1 = st[top][1]; a = st[top][2]; b = st[top][3];
		if (0.5 * (a + b - L * (t1 - t0)) >= req)
			continue;
		if (a < req || b < req || t1 - t0 < 1.0 / 16777216 || top + 2 > 64)
			return 0;
		tm = 0.5 * (t0 + t1);
		dm = dist(user, tm);
		st[top][0] = t0; st[top][1] = tm; st[top][2] = a; st[top][3] = dm; top++;
		st[top][0] = tm; st[top][1] = t1; st[top][2] = dm; st[top][3] = b; top++;
	}
	return 1;
}

/* ------------------------------------------------------------------------ */
/*  Curves                                                                  */
/* ------------------------------------------------------------------------ */

struct qaws_repulsive_curves_solver
{
	rs_core core;                   /* (first: the solver is its core) */
	qaws_repulsive_curves const* c;
	unsigned int *pin, npin, *tedge, nt, nlen;
	unsigned char* fixed;           /* pinned vertices */
	double *pin_target, *tbasis, *len0, L0, bary0[3];
	int total_len, edge_len, bary;
};

/* A = B + B0 at x (nv x nv, row major), the paper's increments over every
   ordered pair of edges sharing no vertex (Repulsive Curves, Sec. 5.2) */
static qaws_status rs_metric(qaws_repulsive_curves const* c, double const* x, double* A)
{
	unsigned int dim = c->dim, nv = c->nv, I, J, o, a, b, i, j, q;
	double sigma = (c->beta - 1) / c->alpha - 1, ex = 2 * sigma + 1;
	memset(A, 0, sizeof(double) * nv * nv);
	for (I = 0; I < c->ne; I++)
		for (J = I + 1; J < c->ne; J++)
		{
			if (rs_shared(c->edges + 2 * I, c->edges + 2 * J))
				continue;
			for (o = 0; o < 2; o++)
			{
				unsigned int const* P = c->edges + 2 * (o ? J : I);
				unsigned int const* Q = c->edges + 2 * (o ? I : J);
				double eP[3] = { 0, 0, 0 }, eQ[3] = { 0, 0, 0 }, lP, lQ, w = 0, w0 = 0, dot = 0;
				for (q = 0; q < dim; q++) { eP[q] = x[dim * P[1] + q] - x[dim * P[0] + q]; eQ[q] = x[dim * Q[1] + q] - x[dim * Q[0] + q]; }
				lP = rs_norm(eP, dim); lQ = rs_norm(eQ, dim);
				if (!(lP > 0 && lQ > 0))
					return QAWS_STATUS_DEGENERATE_CURVE;
				for (q = 0; q < dim; q++) dot += eP[q] * eQ[q] / (lP * lQ);
				for (i = 0; i < 2; i++)
					for (j = 0; j < 2; j++)
					{
						double d[3] = { 0, 0, 0 }, r2 = 0, td = 0, cr2, inv;
						for (q = 0; q < dim; q++) { d[q] = x[dim * P[i] + q] - x[dim * Q[j] + q]; r2 += d[q] * d[q]; td += eP[q] * d[q] / lP; }
						if (!(r2 > 0))
							return QAWS_STATUS_NUMERICAL_FAILURE;
						inv = pow(r2, -0.5 * ex);
						cr2 = fmax(r2 - td * td, 0);   /* |T_P x d|^2 */
						w += inv;
						w0 += cr2 / (r2 * r2) * inv;
					}
				w *= 0.25 * lP * lQ;
				w0 *= 0.25 * lP * lQ;
				for (a = 0; a < 2; a++)
					for (b = 0; b < 2; b++)
					{
						double s = a == b ? 1 : -1;
						A[P[a] * nv + P[b]] += s * w / (lP * lP) + 0.25 * w0;
						A[Q[a] * nv + Q[b]] += s * w / (lQ * lQ) + 0.25 * w0;
						A[P[a] * nv + Q[b]] -= s * w * dot / (lP * lQ) + 0.25 * w0;
						A[Q[a] * nv + P[b]] -= s * w * dot / (lP * lQ) + 0.25 * w0;
					}
			}
		}
	return QAWS_STATUS_OK;
}

static qaws_status rs_c_energy(rs_core const* core, double const* x, double* E, double* grad)
{
	qaws_repulsive_curves const* c = ((qaws_repulsive_curves_solver const*)core)->c;
	return core->theta > 0 ? qaws_internal_repulsive_curves_bh(c, x, core->theta, E, grad) : qaws_internal_repulsive_curves_energy(c, x, E, grad);
}

static qaws_status rs_c_metric(rs_core const* core, double const* x, double* A)
{
	return rs_metric(((qaws_repulsive_curves_solver const*)core)->c, x, A);
}

/* the hierarchical operators: the metric; masses l/2 per edge end, -L by
   1/l per edge; L^(2 - s) by vertex quadrature, r^-(1 + 2 (2 - s)) */
static qaws_status rs_c_hier(rs_core const* core, double const* x, qaws_internal_rs_hier* h)
{
	qaws_repulsive_curves const* c = ((qaws_repulsive_curves_solver const*)core)->c;
	unsigned int e, q, dim = c->dim;
	double s = (c->beta - 1) / c->alpha;
	qaws_status st;
	memset(h, 0, sizeof(*h));
	h->lij = (unsigned int*)malloc(sizeof(unsigned int) * 2 * c->ne);
	h->lw = (double*)malloc(sizeof(double) * c->ne);
	h->mass = (double*)calloc(c->nv + 1, sizeof(double));
	if (!h->lij || !h->lw || !h->mass)
	{
		qaws_internal_rs_hier_free(h);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (e = 0; e < c->ne; e++)
	{
		double d[3] = { 0, 0, 0 }, l;
		for (q = 0; q < dim; q++) d[q] = x[dim * c->edges[2 * e + 1] + q] - x[dim * c->edges[2 * e] + q];
		l = rs_norm(d, dim);
		h->lij[2 * e] = c->edges[2 * e]; h->lij[2 * e + 1] = c->edges[2 * e + 1];
		h->lw[e] = 1 / l;
		h->mass[c->edges[2 * e]] += 0.5 * l; h->mass[c->edges[2 * e + 1]] += 0.5 * l;
	}
	h->nl = c->ne;
	if ((st = qaws_internal_rbh_curves_metric(c, x, core->theta, &h->A)) == QAWS_STATUS_OK)
		st = qaws_internal_rbh_fractional(x, h->mass, c->nv, dim, 1 + 2 * (2 - s), core->theta, &h->F, core->executor);
	if (st != QAWS_STATUS_OK)
		qaws_internal_rs_hier_free(h);
	return st;
}

/* phi (k) and, if C, its Jacobian (k x N, row major) at x */
static void rs_c_constraints(rs_core const* core, double const* x, double* phi, double* C)
{
	qaws_repulsive_curves_solver const* s = (qaws_repulsive_curves_solver const*)core;
	unsigned int dim = core->dim, N = core->N, r = 0, i, q, e;
	if (C) memset(C, 0, sizeof(double) * core->k * N);
	if (s->bary)
		for (q = 0; q < dim; q++, r++)
		{
			double m = 0;
			for (i = 0; i < core->nv; i++)
			{
				m += x[dim * i + q] / core->nv;
				if (C) C[r * N + dim * i + q] = 1.0 / core->nv;
			}
			phi[r] = m - s->bary0[q];
		}
	for (i = 0; core->theta <= 0 && i < s->npin; i++)   /* (hierarchical: eliminated) */
		for (q = 0; q < dim; q++, r++)
		{
			phi[r] = x[dim * s->pin[i] + q] - s->pin_target[dim * i + q];
			if (C) C[r * N + dim * s->pin[i] + q] = 1;
		}
	for (e = 0; e < s->nt; e++)
	{
		unsigned int const* E = s->c->edges + 2 * s->tedge[e];
		unsigned int m;
		for (m = 0; m + 1 < dim; m++, r++)
		{
			double const* n = s->tbasis + dim * ((dim - 1) * e + m);
			phi[r] = 0;
			for (q = 0; q < dim; q++)
			{
				phi[r] += n[q] * (x[dim * E[1] + q] - x[dim * E[0] + q]);
				if (C) { C[r * N + dim * E[1] + q] += n[q]; C[r * N + dim * E[0] + q] -= n[q]; }
			}
		}
	}
	if (s->total_len || s->edge_len)
	{
		double total = 0;
		unsigned int r0 = r;
		for (e = 0; e < s->c->ne; e++)
		{
			unsigned int const* E = s->c->edges + 2 * e;
			double d[3] = { 0, 0, 0 }, l;
			unsigned int row = s->total_len ? r0 : r;
			for (q = 0; q < dim; q++) d[q] = x[dim * E[1] + q] - x[dim * E[0] + q];
			l = rs_norm(d, dim);
			total += l;
			if (s->edge_len && s->fixed[E[0]] && s->fixed[E[1]])
				continue;   /* (both ends pinned: its length is held already) */
			if (s->edge_len) { phi[row] = l - s->len0[e]; r++; }
			if (C && l > 0)
				for (q = 0; q < dim; q++) { C[row * N + dim * E[1] + q] += d[q] / l; C[row * N + dim * E[0] + q] -= d[q] / l; }
		}
		if (s->total_len) phi[r0] = total - s->L0;
	}
}

/* a pair of edges along the straight motion x0 -> x1 */
typedef struct rs_epair
{
	unsigned int const *I, *J;
	double const *x0, *x1;
	unsigned int dim;
} rs_epair;

static double rs_epair_at(void const* user, double t)
{
	rs_epair const* p = (rs_epair const*)user;
	double v[4][3];
	unsigned int k, q;
	unsigned int const* id[4];
	id[0] = &p->I[0]; id[1] = &p->I[1]; id[2] = &p->J[0]; id[3] = &p->J[1];
	for (k = 0; k < 4; k++)
		for (q = 0; q < p->dim; q++) v[k][q] = (1 - t) * p->x0[p->dim * *id[k] + q] + t * p->x1[p->dim * *id[k] + q];
	return qaws_internal_rs_segdist(v[0], v[1], v[2], v[3], p->dim);
}

/* the largest move of a vertex along x0 -> x1 */
static double rs_move(double const* x0, double const* x1, unsigned int v, unsigned int dim)
{
	double m = 0;
	unsigned int q;
	for (q = 0; q < dim; q++) m += (x1[dim * v + q] - x0[dim * v + q]) * (x1[dim * v + q] - x0[dim * v + q]);
	return sqrt(m);
}

/* every pair of edges sharing no vertex stays at least min(clearance, half
   its distance at t = 0) apart along x0 -> x1 (conservative advancement:
   the distance is Lipschitz in t with the pair's largest endpoint moves) */
static int rs_curves_path_ok(qaws_repulsive_curves const* c, double clearance, double const* x0, double const* x1)
{
	unsigned int I, J;
	for (I = 0; I < c->ne; I++)
		for (J = I + 1; J < c->ne; J++)
		{
			rs_epair p;
			double L, d0;
			p.I = c->edges + 2 * I; p.J = c->edges + 2 * J; p.x0 = x0; p.x1 = x1; p.dim = c->dim;
			if (rs_shared(p.I, p.J))
				continue;
			L = fmax(rs_move(x0, x1, p.I[0], c->dim), rs_move(x0, x1, p.I[1], c->dim)) + fmax(rs_move(x0, x1, p.J[0], c->dim), rs_move(x0, x1, p.J[1], c->dim));
			d0 = rs_epair_at(&p, 0);
			if (!qaws_internal_rs_advance(rs_epair_at, &p, d0, L, fmin(clearance, 0.5 * d0)))
				return 0;
		}
	return 1;
}

static int rs_c_path_ok(rs_core const* core, double const* x0, double const* x1)
{
	return rs_curves_path_ok(((qaws_repulsive_curves_solver const*)core)->c, core->clearance, x0, x1);
}

static qaws_status rs_c_hvp(rs_core const* core, double const* x, double const* v, double* hv)
{
	return qaws_internal_repulsive_curves_hvp(((qaws_repulsive_curves_solver const*)core)->c, x, v, hv);
}

/* hv += sum_r lambda_r Hess(c_r) v: only the length rows are not linear;
   Hess(l_e) v = (I - T T^T)(v_1 - v_0) / l at the edge's end, minus at its start */
static void rs_c_constraints_hvp(rs_core const* core, double const* x, double const* lambda, double const* v, double* hv)
{
	qaws_repulsive_curves_solver const* s = (qaws_repulsive_curves_solver const*)core;
	unsigned int dim = core->dim, r, e, q;
	if (!s->total_len && !s->edge_len)
		return;
	r = (s->bary ? dim : 0) + (core->theta > 0 ? 0 : dim * s->npin) + (dim - 1) * s->nt;
	for (e = 0; e < s->c->ne; e++)
	{
		unsigned int const* E = s->c->edges + 2 * e;
		double d[3] = { 0, 0, 0 }, w[3] = { 0, 0, 0 }, l, tw = 0, lam;
		if (s->edge_len && s->fixed[E[0]] && s->fixed[E[1]])
			continue;
		lam = lambda[s->total_len ? r : r++];
		for (q = 0; q < dim; q++) { d[q] = x[dim * E[1] + q] - x[dim * E[0] + q]; w[q] = v[dim * E[1] + q] - v[dim * E[0] + q]; }
		l = rs_norm(d, dim);
		if (!(l > 0))
			continue;
		for (q = 0; q < dim; q++) tw += d[q] * w[q] / l;
		for (q = 0; q < dim; q++)
		{
			double h = lam * (w[q] - d[q] / l * tw) / l;
			hv[dim * E[1] + q] += h; hv[dim * E[0] + q] -= h;
		}
	}
}

static qaws_internal_rs_ops const rs_curve_ops = { rs_c_energy, rs_c_metric, rs_c_constraints, rs_c_path_ok, rs_c_hier, NULL, rs_c_hvp,
	rs_c_constraints_hvp };

/* the rows' targets from t (to rows, dir 0) or back (dir 1) */
static void rs_c_targets(qaws_repulsive_curves_solver const* s, qaws_repulsive_targets* t, double* rows, int dir)
{
	unsigned int dim = s->core.dim, r = 0, i, q, e;
	if (dir)
	{
		memset(t->barycenter, 0, sizeof(t->barycenter));
		t->length = t->area = t->volume = 0;
	}
	if (s->bary)
		for (q = 0; q < dim; q++, r++)
		{
			if (dir) t->barycenter[q] = (qaws_scalar)rows[r]; else rows[r] = (double)t->barycenter[q];
		}
	for (i = 0; s->core.theta <= 0 && i < s->npin; i++)   /* (hierarchical: pins are not rows) */
		for (q = 0; q < dim; q++, r++)
		{
			if (dir) { if (t->pinned) t->pinned[dim * i + q] = (qaws_scalar)rows[r]; }
			else rows[r] = t->pinned ? (double)t->pinned[dim * i + q] : 0;
		}
	for (i = 0; i < (dim - 1) * s->nt; i++, r++)
		if (!dir) rows[r] = 0;   /* (held directions: no target) */
	if (s->total_len)
	{
		if (dir) t->length = (qaws_scalar)rows[r]; else rows[r] = (double)t->length;
		r++;
	}
	for (e = 0; s->edge_len && e < s->c->ne; e++)
	{
		double val = 0;
		if (!(s->fixed[s->c->edges[2 * e]] && s->fixed[s->c->edges[2 * e + 1]]))
		{
			if (dir) val = rows[r]; else rows[r] = t->edge_lengths ? (double)t->edge_lengths[e] : 0;
			r++;
		}
		if (dir && t->edge_lengths) t->edge_lengths[e] = (qaws_scalar)val;
	}
}

qaws_status qaws_repulsive_curves_solver_solution_jvp(qaws_repulsive_curves_solver const* s, qaws_repulsive_targets const* d, qaws_scalar* out_dx,
	unsigned int* out_gauge)
{
	double *rows, *dx;
	unsigned int i;
	qaws_status st;
	if (!s || !d || !out_dx)
		return QAWS_STATUS_INVALID_ARGUMENT;
	rows = (double*)calloc(s->core.k + 1, sizeof(double));
	dx = (double*)malloc(sizeof(double) * (s->core.N + 1));
	if (!rows || !dx)
	{
		free(rows); free(dx);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	rs_c_targets(s, (qaws_repulsive_targets*)d, rows, 0);
	st = qaws_internal_rs_sensitivity(&s->core, rows, NULL, dx, out_gauge);
	for (i = 0; st == QAWS_STATUS_OK && i < s->core.N; i++) out_dx[i] = (qaws_scalar)dx[i];
	free(rows); free(dx);
	return st;
}

qaws_status qaws_repulsive_curves_solver_solution_vjp(qaws_repulsive_curves_solver const* s, qaws_scalar const* xbar, qaws_repulsive_targets* out_bar,
	unsigned int* out_gauge)
{
	double *rows, *xb;
	unsigned int i;
	qaws_status st;
	if (!s || !xbar || !out_bar)
		return QAWS_STATUS_INVALID_ARGUMENT;
	rows = (double*)calloc(s->core.k + 1, sizeof(double));
	xb = (double*)malloc(sizeof(double) * (s->core.N + 1));
	if (!rows || !xb)
	{
		free(rows); free(xb);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < s->core.N; i++) xb[i] = (double)xbar[i];
	st = qaws_internal_rs_sensitivity(&s->core, NULL, xb, rows, out_gauge);
	if (st == QAWS_STATUS_OK)
		rs_c_targets(s, out_bar, rows, 1);
	free(rows); free(xb);
	return st;
}

void qaws_repulsive_curves_solver_destroy(qaws_repulsive_curves_solver* s)
{
	if (!s)
		return;
	free(s->fixed); free(s->pin); free(s->pin_target); free(s->tedge); free(s->tbasis); free(s->len0); free(s->core.x);
	free(s);
}

qaws_status qaws_repulsive_curves_solver_create(qaws_repulsive_curves const* c, qaws_repulsive_curves_solve_desc const* desc,
	qaws_scalar const* positions, qaws_repulsive_curves_solver** out_solver)
{
	qaws_repulsive_curves_solver* s;
	rs_core* k;
	qaws_repulsive_curves_solve_desc d;
	unsigned int i, q, dim;
	double total = 0;
	if (!c || !positions || !out_solver)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_solver = NULL;
	memset(&d, 0, sizeof(d));
	if (desc) d = *desc;
	dim = c->dim;
	if ((d.pinned_count && !d.pinned) || (d.tangent_edge_count && !d.tangent_edges) || (d.hold_total_length && d.hold_edge_lengths) || c->ne == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < d.pinned_count; i++) if (d.pinned[i] >= c->nv) return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < d.tangent_edge_count; i++) if (d.tangent_edges[i] >= c->ne) return QAWS_STATUS_INVALID_ARGUMENT;
	s = (qaws_repulsive_curves_solver*)calloc(1, sizeof(*s));
	if (!s)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	k = &s->core;
	k->ops = &rs_curve_ops;
	s->c = c; k->dim = dim; k->nv = c->nv; k->N = dim * c->nv;
	s->npin = d.pinned_count;
	s->total_len = d.hold_total_length != 0; s->edge_len = d.hold_edge_lengths != 0;
	s->bary = d.hold_barycenter > 0 || (d.hold_barycenter == 0 && d.pinned_count == 0);
	k->checked = d.collision_checked != 0;
	k->x = (double*)malloc(sizeof(double) * k->N);
	s->pin = (unsigned int*)malloc(sizeof(unsigned int) * (s->npin + 1));
	s->pin_target = (double*)malloc(sizeof(double) * (dim * s->npin + 1));
	s->tedge = (unsigned int*)malloc(sizeof(unsigned int) * (d.tangent_edge_count + 1));
	s->tbasis = (double*)malloc(sizeof(double) * (dim * (dim - 1) * d.tangent_edge_count + 1));
	s->fixed = (unsigned char*)calloc(c->nv + 1, 1);
	s->len0 = (double*)malloc(sizeof(double) * (c->ne + 1));
	if (!k->x || !s->pin || !s->pin_target || !s->tedge || !s->tbasis || !s->len0 || !s->fixed)
	{
		qaws_repulsive_curves_solver_destroy(s);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < k->N; i++) k->x[i] = (double)positions[i];
	/* the targets */
	for (q = 0; q < dim; q++)
	{
		s->bary0[q] = 0;
		for (i = 0; i < k->nv; i++) s->bary0[q] += k->x[dim * i + q] / k->nv;
	}
	for (i = 0; i < s->npin; i++)
	{
		s->pin[i] = d.pinned[i];
		for (q = 0; q < dim; q++) s->pin_target[dim * i + q] = k->x[dim * d.pinned[i] + q];
		s->fixed[d.pinned[i]] = 1;
	}
	for (i = 0; i < c->ne; i++)
	{
		double e[3] = { 0, 0, 0 };
		for (q = 0; q < dim; q++) e[q] = k->x[dim * c->edges[2 * i + 1] + q] - k->x[dim * c->edges[2 * i] + q];
		s->len0[i] = rs_norm(e, dim);
		total += s->len0[i];
		if (!(s->len0[i] > 0))
		{
			qaws_repulsive_curves_solver_destroy(s);
			return QAWS_STATUS_DEGENERATE_CURVE;
		}
	}
	s->L0 = total;
	k->lmean = total / c->ne;
	for (i = 0; i < d.tangent_edge_count; i++)
	{
		/* the normals of the held direction: their products with the edge stay 0
		   (an edge with both ends pinned holds its direction already) */
		unsigned int const* E = c->edges + 2 * d.tangent_edges[i];
		double t[3] = { 0, 0, 0 }, l;
		double* n = s->tbasis + dim * (dim - 1) * s->nt;
		if (s->fixed[E[0]] && s->fixed[E[1]])
			continue;
		s->tedge[s->nt++] = d.tangent_edges[i];
		for (q = 0; q < dim; q++) t[q] = k->x[dim * E[1] + q] - k->x[dim * E[0] + q];
		l = rs_norm(t, dim);
		for (q = 0; q < dim; q++) t[q] /= l;
		if (dim == 2)
		{
			n[0] = -t[1]; n[1] = t[0];
		}
		else
		{
			double a[3] = { 0, 0, 0 }, m;
			unsigned int lo = fabs(t[0]) <= fabs(t[1]) && fabs(t[0]) <= fabs(t[2]) ? 0 : fabs(t[1]) <= fabs(t[2]) ? 1 : 2;
			a[lo] = 1;
			n[0] = t[1] * a[2] - t[2] * a[1]; n[1] = t[2] * a[0] - t[0] * a[2]; n[2] = t[0] * a[1] - t[1] * a[0];
			m = rs_norm(n, 3);
			for (q = 0; q < 3; q++) n[q] /= m;
			n[3] = t[1] * n[2] - t[2] * n[1]; n[4] = t[2] * n[0] - t[0] * n[2]; n[5] = t[0] * n[1] - t[1] * n[0];
		}
	}
	for (i = 0; s->edge_len && i < c->ne; i++) s->nlen += !(s->fixed[c->edges[2 * i]] && s->fixed[c->edges[2 * i + 1]]);
	k->theta = (double)d.theta;
	k->executor = c->executor;
	k->fixed = s->fixed;
	k->translates = s->npin == 0 && !s->bary;
	k->k = (s->bary ? dim : 0) + (k->theta > 0 ? 0 : dim * s->npin) + (dim - 1) * s->nt + (s->total_len ? 1 : 0) + s->nlen;
	if (!(d.theta >= 0) || (k->theta > 0 && k->k > QAWS_INTERNAL_RS_HROWS))
	{
		/* (hierarchical: a solve per row and coordinate; per-edge lengths are too many) */
		qaws_repulsive_curves_solver_destroy(s);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	k->clearance = d.clearance > 0 ? (double)d.clearance : 1e-3 * k->lmean;
	k->tol = d.tolerance > 0 ? (double)d.tolerance : 1e-4;
	k->ctol = (d.constraint_tolerance > 0 ? (double)d.constraint_tolerance : RS_CTOL) * k->lmean;
	k->h = 0.25 * k->lmean;
	k->rep.reason = QAWS_REPULSIVE_RUNNING;
	k->rep.stationarity = 1;
	if (k->ops->energy(k, k->x, &total, NULL) != QAWS_STATUS_OK)
	{
		qaws_repulsive_curves_solver_destroy(s);
		return QAWS_STATUS_NUMERICAL_FAILURE;
	}
	k->rep.initial_energy = k->rep.energy = (qaws_scalar)total;
	/* collision-checked: no two edges sharing no vertex may touch already */
	if (k->checked)
	{
		unsigned int I, J;
		for (I = 0; I < c->ne; I++)
			for (J = I + 1; J < c->ne; J++)
			{
				rs_epair p;
				p.I = c->edges + 2 * I; p.J = c->edges + 2 * J; p.x0 = p.x1 = k->x; p.dim = dim;
				if (!rs_shared(p.I, p.J) && !(rs_epair_at(&p, 0) > 1e-12 * k->lmean))
				{
					qaws_repulsive_curves_solver_destroy(s);
					return QAWS_STATUS_INVALID_ARGUMENT;
				}
			}
	}
	*out_solver = s;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_curves_solver_direction(qaws_repulsive_curves_solver const* s, qaws_scalar const* positions, qaws_scalar* out_direction)
{
	return s ? qaws_internal_rs_direction(&s->core, positions, out_direction) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_curves_solver_step(qaws_repulsive_curves_solver* s, qaws_scalar* positions, qaws_repulsive_report* out_report)
{
	return s ? qaws_internal_rs_step(&s->core, positions, out_report) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_curves_solve(qaws_repulsive_curves const* curves, qaws_repulsive_curves_solve_desc const* desc,
	qaws_scalar* positions, qaws_repulsive_report* out_report)
{
	qaws_repulsive_curves_solver* s = NULL;
	qaws_repulsive_report rep;
	unsigned int it, max = desc && desc->max_iterations ? desc->max_iterations : 100;
	qaws_status st = qaws_repulsive_curves_solver_create(curves, desc, positions, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(&rep, 0, sizeof(rep));
	for (it = 0; it < max && st == QAWS_STATUS_OK; it++)
	{
		st = qaws_repulsive_curves_solver_step(s, positions, &rep);
		if (rep.reason != QAWS_REPULSIVE_RUNNING)
			break;
	}
	if (st == QAWS_STATUS_OK && rep.reason == QAWS_REPULSIVE_RUNNING)
		rep.reason = QAWS_REPULSIVE_ITERATION_LIMIT;
	if (out_report) *out_report = rep;
	qaws_repulsive_curves_solver_destroy(s);
	return st;
}

/* ------------------------------------------------------------------------ */
/*  For tests                                                               */
/* ------------------------------------------------------------------------ */

qaws_status qaws_internal_repulsive_curves_metric(qaws_repulsive_curves const* curves, double const* x, double* A)
{
	return rs_metric(curves, x, A);
}

int qaws_internal_repulsive_curves_path_clear(qaws_repulsive_curves const* curves, double clearance, double const* x0, double const* x1)
{
	return rs_curves_path_ok(curves, clearance, x0, x1);
}
