/*
 * The hierarchical mode of the fractional Sobolev descent (qaws_repulsive.h,
 * R4; Repulsive Surfaces, Sec. 5.2 and 5.3), on the core of
 * qaws_repulsive_solve.c: nothing of size |V|^2 is stored.
 *
 * The constrained direction [A~ C^T; C 0] [g; lambda] = [dE; 0] is solved
 * with A applied from Barnes-Hut lists (qaws_repulsive_bh.c):
 *   - pinned vertices are eliminated (their rows and columns dropped, g 0
 *     there), so C holds the other rows only (barycenter, lengths, area,
 *     volume, held directions; at most QAWS_INTERNAL_RS_HROWS);
 *   - with nothing pinned, A's kernel (the constants, a translation) is
 *     lifted by A + tau 1 1^T / |V|; the barycenter rows keep g's mean 0, so
 *     the solution is unchanged;
 *   - z = A~^-1 dE and Y = A~^-1 C^T by flexible GMRES, one scalar solve per
 *     coordinate a right-hand side has; then the Schur complement
 *     (C Y) lambda = C z, g = z - Y lambda.
 * The preconditioner is the paper's (-L)^-1 L^(2 - s) (-L)^-1: L the
 * Laplacian (1/l weights on curves, cotan on meshes), shifted by a small
 * multiple of the lumped mass so that every component is invertible, and
 * solved by a profile Cholesky factored once per direction; L^(2 - s) the
 * fractional operator by vertex quadrature, applied hierarchically too.
 * The constants (on the free vertices) are handled apart: the
 * preconditioner returns their mean over tau.
 *
 * The Newton projection after a step reuses Y: dx = -Y m, (C(y) Y) m =
 * phi(y).
 */

#include "qaws_repulsive.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include "internal/qaws_internal_parallel.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RH_RESTART 40        /* the defaults of the desc's gmres_* settings */
#define RH_MAXIT 1000
#define RH_TOL 1e-8          /* each scalar solve, relative */
#define RH_ACCEPT 1e-6       /* accepted once restarts stop gaining: the floor rounding sets in the approximate products */

typedef qaws_internal_rs_core rs_core;

struct qaws_internal_rs_hsys
{
	rs_core const* s;
	qaws_internal_rs_hier h;
	double tau, shift;
	int aug;                        /* nothing pinned: A + tau 1 1^T / |V| */
	unsigned int nfree, k;
	double *Y, *C, *phi, *S;        /* k x N, k x N, k, k x k */
	unsigned int* piv;
	double *t0, *t1;                /* nv each */
	unsigned int *perm, *inv, *first;  /* the Laplacian's ordering and profile */
	size_t* start;
	double* L;
};

void qaws_internal_rs_hier_free(qaws_internal_rs_hier* h)
{
	qaws_internal_rbh_destroy(h->A);
	qaws_internal_rbh_destroy(h->F);
	qaws_internal_free(h->lij); qaws_internal_free(h->lw); qaws_internal_free(h->mass);
	memset(h, 0, sizeof(*h));
}

void qaws_internal_rs_hsys_free(qaws_internal_rs_hsys* H)
{
	if (!H)
		return;
	qaws_internal_rs_hier_free(&H->h);
	qaws_internal_free(H->Y); qaws_internal_free(H->C); qaws_internal_free(H->phi); qaws_internal_free(H->S); qaws_internal_free(H->piv);
	qaws_internal_free(H->t0); qaws_internal_free(H->t1);
	qaws_internal_free(H->perm); qaws_internal_free(H->inv); qaws_internal_free(H->first); qaws_internal_free(H->start); qaws_internal_free(H->L);
	qaws_internal_free(H);
}

static void rh_mask(rs_core const* s, double* v)
{
	unsigned int i;
	if (s->fixed)
		for (i = 0; i < s->nv; i++)
			if (s->fixed[i]) v[i] = 0;
}

static double rh_dot(double const* a, double const* b, unsigned int n)
{
	double r = 0;
	unsigned int i;
	for (i = 0; i < n; i++) r += a[i] * b[i];
	return r;
}

/* one scalar solve's workspace: GMRES's vectors, the products' scratch */
typedef struct rh_work
{
	double *V, *Z, *H, *cs, *sn, *gq, *y;
	unsigned int R;                 /* the restart */
	double *t0, *t3, *b, *r, *tmp;  /* nv each */
	qaws_internal_rbh_scratch *sa, *sf;
} rh_work;

static void rh_work_destroy(rh_work* W)
{
	if (!W)
		return;
	qaws_internal_free(W->V); qaws_internal_free(W->Z); qaws_internal_free(W->H); qaws_internal_free(W->cs); qaws_internal_free(W->sn); qaws_internal_free(W->gq); qaws_internal_free(W->y);
	qaws_internal_free(W->t0); qaws_internal_free(W->t3); qaws_internal_free(W->b); qaws_internal_free(W->r); qaws_internal_free(W->tmp);
	qaws_internal_rbh_scratch_destroy(W->sa);
	qaws_internal_rbh_scratch_destroy(W->sf);
	qaws_internal_free(W);
}

static rh_work* rh_work_create(qaws_internal_rs_hsys const* H)
{
	size_t nv = H->s->nv, R = H->s->gmres_restart ? H->s->gmres_restart : RH_RESTART;
	rh_work* W = (rh_work*)qaws_internal_zalloc(1, sizeof(rh_work));
	if (!W)
		return NULL;
	W->R = (unsigned int)R;
	W->y = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), R));
	W->V = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), nv, (R + 1)));
	W->Z = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), nv, (R + 1)));
	W->H = (double*)qaws_internal_zalloc(qaws_internal_mul2((R + 1), R), sizeof(double));
	W->cs = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), R)); W->sn = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), R));
	W->gq = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (R + 1)));
	W->t0 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nv)); W->t3 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nv));
	W->b = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nv)); W->r = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nv)); W->tmp = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nv));
	W->sa = qaws_internal_rbh_scratch_create(H->h.A);
	W->sf = qaws_internal_rbh_scratch_create(H->h.F);
	if (!W->y || !W->V || !W->Z || !W->H || !W->cs || !W->sn || !W->gq || !W->t0 || !W->t3 || !W->b || !W->r || !W->tmp || !W->sa || !W->sf)
	{
		rh_work_destroy(W);
		return NULL;
	}
	return W;
}

/* y = A u on the free vertices (u's pinned entries ignored) */
static void rh_apply(qaws_internal_rs_hsys const* H, rh_work* W, double const* u, double* y)
{
	rs_core const* s = H->s;
	unsigned int i;
	memcpy(W->t3, u, sizeof(double) * s->nv);
	rh_mask(s, W->t3);
	qaws_internal_rbh_apply_with(H->h.A, W->sa, W->t3, y);
	rh_mask(s, y);
	if (H->aug)
	{
		double m = 0;
		for (i = 0; i < s->nv; i++) m += W->t3[i] / s->nv;
		for (i = 0; i < s->nv; i++) y[i] += H->tau * m;
	}
}

/* -L + shift M on the free vertices (the pinned ones' rows the identity),
   factored once at x by profile Cholesky in a reverse Cuthill-McKee order.
   shortcut: the profile is |V| times the ordering's bandwidth (about
   sqrt|V| on a surface), fine to some 1e4 vertices; a nested-dissection
   ordering, or conjugate gradients, beyond */
static int rh_chol_make(qaws_internal_rs_hsys* H)
{
	rs_core const* s = H->s;
	unsigned int nv = s->nv, i, l, head = 0, tail = 0, n = 0;
	unsigned int *deg = (unsigned int*)qaws_internal_zalloc(nv + 1, sizeof(unsigned int)), *adj0 = (unsigned int*)qaws_internal_zalloc(nv + 2, sizeof(unsigned int)), *adj = NULL;
	unsigned char* seen = (unsigned char*)qaws_internal_scratch_zalloc(nv + 1, 1);
	size_t total = 0;
	int ok = 0;
	H->perm = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nv + 1)));
	H->inv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nv + 1)));
	H->first = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nv + 1)));
	H->start = (size_t*)qaws_internal_malloc(qaws_internal_mul2(sizeof(size_t), (nv + 1)));
	if (!deg || !adj0 || !seen || !H->perm || !H->inv || !H->first || !H->start)
		goto done;
	/* the adjacency (free vertices only) */
	for (l = 0; l < H->h.nl; l++)
	{
		unsigned int a = H->h.lij[2 * l], b = H->h.lij[2 * l + 1];
		if (a != b && !(s->fixed && (s->fixed[a] || s->fixed[b]))) { deg[a]++; deg[b]++; }
	}
	for (i = 0; i < nv; i++) adj0[i + 1] = adj0[i] + deg[i];
	adj = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (adj0[nv] + 1)));
	if (!adj)
		goto done;
	memset(deg, 0, sizeof(unsigned int) * nv);
	for (l = 0; l < H->h.nl; l++)
	{
		unsigned int a = H->h.lij[2 * l], b = H->h.lij[2 * l + 1];
		if (a != b && !(s->fixed && (s->fixed[a] || s->fixed[b]))) { adj[adj0[a] + deg[a]++] = b; adj[adj0[b] + deg[b]++] = a; }
	}
	/* Cuthill-McKee: breadth first from a least-degree vertex per component,
	   neighbors by increasing degree; then reversed */
	while (n < nv)
	{
		unsigned int best = nv;
		for (i = 0; i < nv; i++)
			if (!seen[i] && (best == nv || deg[i] < deg[best])) best = i;
		seen[best] = 1;
		H->inv[tail++] = best;
		n++;
		while (head < tail)
		{
			unsigned int v = H->inv[head++], from = tail, a, b;
			for (a = adj0[v]; a < adj0[v + 1]; a++)
				if (!seen[adj[a]]) { seen[adj[a]] = 1; H->inv[tail++] = adj[a]; n++; }
			for (a = from + 1; a < tail; a++)
				for (b = a; b > from && deg[H->inv[b - 1]] > deg[H->inv[b]]; b--) { unsigned int t = H->inv[b]; H->inv[b] = H->inv[b - 1]; H->inv[b - 1] = t; }
		}
	}
	for (i = 0; i < nv / 2; i++) { unsigned int t = H->inv[i]; H->inv[i] = H->inv[nv - 1 - i]; H->inv[nv - 1 - i] = t; }
	for (i = 0; i < nv; i++) H->perm[H->inv[i]] = i;
	/* the profile: row r from first[r] to r */
	for (i = 0; i < nv; i++)
	{
		unsigned int r = H->perm[i], a;
		H->first[r] = r;
		for (a = adj0[i]; a < adj0[i + 1]; a++) if (H->perm[adj[a]] < H->first[r]) H->first[r] = H->perm[adj[a]];
	}
	for (i = 0; i < nv; i++) { H->start[i] = total; total += i - H->first[i] + 1; }
	H->L = (double*)qaws_internal_zalloc(total + 1, sizeof(double));
	if (!H->L)
		goto done;
	/* assembled, then factored in place: L L^T */
	for (i = 0; i < nv; i++)
	{
		unsigned int r = H->perm[i];
		H->L[H->start[r] + r - H->first[r]] = s->fixed && s->fixed[i] ? 1 : H->shift * H->h.mass[i];
	}
	for (l = 0; l < H->h.nl; l++)
	{
		unsigned int a = H->h.lij[2 * l], b = H->h.lij[2 * l + 1], ra, rb;
		if (a == b || (s->fixed && (s->fixed[a] || s->fixed[b])))
		{
			/* (a pinned neighbor: its term stays on the free diagonal) */
			if (a != b && !s->fixed[a]) H->L[H->start[H->perm[a]] + H->perm[a] - H->first[H->perm[a]]] += H->h.lw[l];
			if (a != b && !s->fixed[b]) H->L[H->start[H->perm[b]] + H->perm[b] - H->first[H->perm[b]]] += H->h.lw[l];
			continue;
		}
		ra = H->perm[a]; rb = H->perm[b];
		H->L[H->start[ra] + ra - H->first[ra]] += H->h.lw[l];
		H->L[H->start[rb] + rb - H->first[rb]] += H->h.lw[l];
		if (ra < rb) { unsigned int t = ra; ra = rb; rb = t; }
		H->L[H->start[ra] + rb - H->first[ra]] -= H->h.lw[l];
	}
	for (i = 0; i < nv; i++)
	{
		unsigned int j, k;
		double* Li = H->L + H->start[i] - H->first[i];
		for (j = H->first[i]; j <= i; j++)
		{
			double* Lj = H->L + H->start[j] - H->first[j];
			double sum = Li[j];
			for (k = H->first[i] > H->first[j] ? H->first[i] : H->first[j]; k < j; k++) sum -= Li[k] * Lj[k];
			if (j < i)
				Li[j] = sum / Lj[j];
			else if (sum > 0)
				Li[i] = sqrt(sum);
			else
				goto done;   /* (not positive definite) */
		}
	}
	ok = 1;
done:
	qaws_internal_free(deg); qaws_internal_free(adj0); qaws_internal_free(adj); qaws_internal_free(seen);
	return ok;
}

/* (-L + shift M) z = b on the free vertices (W's t0 as scratch) */
static void rh_chol_solve(qaws_internal_rs_hsys const* H, rh_work* W, double const* b, double* z)
{
	unsigned int nv = H->s->nv, i, j;
	double* y = W->t0;
	for (i = 0; i < nv; i++) y[H->perm[i]] = b[i];
	for (i = 0; i < nv; i++)
	{
		double const* Li = H->L + H->start[i] - H->first[i];
		for (j = H->first[i]; j < i; j++) y[i] -= Li[j] * y[j];
		y[i] /= Li[i];
	}
	for (i = nv; i-- > 0;)
	{
		double const* Li = H->L + H->start[i] - H->first[i];
		y[i] /= Li[i];
		for (j = H->first[i]; j < i; j++) y[j] -= Li[j] * y[i];
	}
	for (i = 0; i < nv; i++) z[i] = y[H->perm[i]];
	rh_mask(H->s, z);
}

/* z = P v: the free constants by tau, the rest by the Laplacian sandwich */
static void rh_precond(qaws_internal_rs_hsys const* H, rh_work* W, double const* v, double* z)
{
	rs_core const* s = H->s;
	unsigned int i;
	double c = 0, *tmp = W->tmp;
	for (i = 0; i < s->nv; i++) tmp[i] = v[i];
	rh_mask(s, tmp);
	for (i = 0; i < s->nv; i++) c += tmp[i] / H->nfree;
	for (i = 0; i < s->nv; i++) if (!s->fixed || !s->fixed[i]) tmp[i] -= c;
	rh_chol_solve(H, W, tmp, z);
	qaws_internal_rbh_apply_with(H->h.F, W->sf, z, tmp);
	rh_mask(s, tmp);
	rh_chol_solve(H, W, tmp, z);
	for (i = 0; i < s->nv; i++) if (!s->fixed || !s->fixed[i]) z[i] += c / H->tau;
}

/* A u = b on the free vertices by flexible GMRES (restarted); 0 when it
   does not converge */
static int rh_solve(qaws_internal_rs_hsys const* H, rh_work* W, double const* b0, double* u)
{
	size_t R = W->R;
	unsigned int nv = H->s->nv, i, j, k, outer, it = 0;
	unsigned int maxit = H->s->gmres_max ? H->s->gmres_max : RH_MAXIT, maxouter = (unsigned int)(maxit / R + (maxit % R != 0));
	double *V = W->V, *Z = W->Z, *Hm = W->H, *b = W->b, *r = W->r, *y = W->y, bn, prev = HUGE_VAL, tol = H->s->gmres_tol > 0 ? H->s->gmres_tol : RH_TOL;
	int ok = 0;
	memset(u, 0, sizeof(double) * nv);
	memcpy(b, b0, sizeof(double) * nv);
	rh_mask(H->s, b);
	bn = sqrt(rh_dot(b, b, nv));
	if (!(bn > 0))
		return 1;
	for (outer = 0; outer < maxouter; outer++)   /* (the defaults: 25 restarts of 40) */
	{
		double beta;
		rh_apply(H, W, u, r);
		for (i = 0; i < nv; i++) r[i] = b[i] - r[i];
		beta = sqrt(rh_dot(r, r, nv));
		if (beta <= tol * bn || (beta <= RH_ACCEPT * bn && beta > 0.5 * prev))
		{
			ok = 1;
			break;
		}
		prev = beta;
		for (i = 0; i < nv; i++) V[i] = r[i] / beta;
		memset(W->gq, 0, sizeof(double) * (R + 1));
		W->gq[0] = beta;
		for (j = 0; j < R; j++)
		{
			double* z = Z + (size_t)j * nv;
			double hn;
			rh_precond(H, W, V + (size_t)j * nv, z);
			rh_apply(H, W, z, r);
			for (k = 0; k <= j; k++)
			{
				double hk = rh_dot(r, V + (size_t)k * nv, nv);
				Hm[k * R + j] = hk;
				for (i = 0; i < nv; i++) r[i] -= hk * V[(size_t)k * nv + i];
			}
			hn = sqrt(rh_dot(r, r, nv));
			Hm[(j + 1) * R + j] = hn;
			if (hn > 0)
				for (i = 0; i < nv; i++) V[(size_t)(j + 1) * nv + i] = r[i] / hn;
			for (k = 0; k < j; k++)
			{
				double a = Hm[k * R + j], c2 = Hm[(k + 1) * R + j];
				Hm[k * R + j] = W->cs[k] * a + W->sn[k] * c2;
				Hm[(k + 1) * R + j] = -W->sn[k] * a + W->cs[k] * c2;
			}
			{
				double a = Hm[j * R + j], c2 = Hm[(j + 1) * R + j], rr = hypot(a, c2);
				W->cs[j] = rr > 0 ? a / rr : 1; W->sn[j] = rr > 0 ? c2 / rr : 0;
				Hm[j * R + j] = rr; Hm[(j + 1) * R + j] = 0;
				W->gq[j + 1] = -W->sn[j] * W->gq[j]; W->gq[j] *= W->cs[j];
			}
			it++;
			if (fabs(W->gq[j + 1]) <= tol * bn || j + 1 == R || !(hn > 0) || it == maxit)
			{
				unsigned int n = j + 1;
				for (k = n; k-- > 0;)
				{
					double sm = W->gq[k];
					for (i = k + 1; i < n; i++) sm -= Hm[k * R + i] * y[i];
					y[k] = Hm[k * R + k] != 0 ? sm / Hm[k * R + k] : 0;
				}
				for (k = 0; k < n; k++)
					for (i = 0; i < nv; i++) u[i] += y[k] * Z[(size_t)k * nv + i];
				break;
			}
		}
	}
	rh_mask(H->s, u);
	return ok;
}

/* The scalar solves of a direction: for each N-vector source and each
   coordinate it has free entries in, A u = that coordinate. Identical
   right-hand sides (the barycenter rows, one per coordinate) are solved
   once. The distinct ones are independent: they run through the executor,
   each with its own workspace, and give the same result in any order. */
typedef struct rh_batch
{
	qaws_internal_rs_hsys const* H;
	double const* rhs;              /* the distinct right-hand sides, nv each */
	double* sol;
	int* ok;
} rh_batch;

static qaws_status rh_batch_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	rh_batch* b = (rh_batch*)user;
	unsigned int nv = b->H->s->nv, j;
	rh_work* W = rh_work_create(b->H);
	(void)chunk;
	if (!W)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (j = begin; j < end; j++) b->ok[j] = rh_solve(b->H, W, b->rhs + (size_t)j * nv, b->sol + (size_t)j * nv);
	rh_work_destroy(W);
	return QAWS_STATUS_OK;
}

static int rh_solve_all(qaws_internal_rs_hsys const* H, double const* const* src, double* const* dst, unsigned int count)
{
	rs_core const* s = H->s;
	unsigned int dim = s->dim, nv = s->nv, njob = 0, nuniq = 0, j, q, i, t;
	unsigned int* which = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (count * dim + 1)));   /* job -> distinct system */
	double *rhs = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)count * dim * nv + 1))), *sol = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)count * dim * nv + 1)));
	int* ok = (int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(int), (count * dim + 1)));
	int good = which && rhs && sol && ok;
	rh_batch b;
	for (j = 0; good && j < count; j++)
		for (q = 0; q < dim; q++)
		{
			double* v = rhs + (size_t)nuniq * nv;
			int any = 0;
			for (i = 0; i < nv; i++)
			{
				v[i] = src[j][dim * i + q];
				if (s->fixed && s->fixed[i]) v[i] = 0;
				any |= v[i] != 0;
			}
			if (!any)
			{
				which[njob++] = 0xffffffffu;   /* (nothing free: u = 0) */
				continue;
			}
			for (t = 0; t < nuniq && memcmp(rhs + (size_t)t * nv, v, sizeof(double) * nv) != 0; t++)
				;
			which[njob++] = t;
			if (t == nuniq)
				nuniq++;
		}
	if (good)
	{
		b.H = H; b.rhs = rhs; b.sol = sol; b.ok = ok;
		good = qaws_internal_parallel(s->executor, nuniq, 1, rh_batch_job, &b) == QAWS_STATUS_OK;
	}
	for (t = 0; good && t < nuniq; t++) good = ok[t];
	for (j = 0, njob = 0; good && j < count; j++)
		for (q = 0; q < dim; q++, njob++)
			for (i = 0; i < nv; i++) dst[j][dim * i + q] = which[njob] == 0xffffffffu ? 0 : sol[(size_t)which[njob] * nv + i];
	qaws_internal_free(which); qaws_internal_free(rhs); qaws_internal_free(sol); qaws_internal_free(ok);
	return good;
}

qaws_status qaws_internal_rs_hdirection(rs_core const* s, double const* x, double const* dE, double* g, qaws_internal_rs_hsys** out)
{
	qaws_internal_rs_hsys* H;
	unsigned int nv = s->nv, N = s->N, k = s->k, i, r, c, q;
	qaws_status st;
	double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 }, diag2 = 0;
	*out = NULL;
	if (s->translates)
		return QAWS_STATUS_NUMERICAL_FAILURE;   /* (as the dense system: singular) */
	H = (qaws_internal_rs_hsys*)qaws_internal_zalloc(1, sizeof(*H));
	if (!H)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	H->s = s; H->k = k;
	H->Y = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)k * N + 1)));
	H->C = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)k * N + 1)));
	H->phi = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (k + 1)));
	H->S = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (k * k + 1)));
	H->piv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (k + 1)));
	H->t0 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nv)); H->t1 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), nv));
	if (!H->Y || !H->C || !H->phi || !H->S || !H->piv || !H->t0 || !H->t1)
	{
		qaws_internal_rs_hsys_free(H);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	st = s->ops->hier(s, x, &H->h);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_rs_hsys_free(H);
		return st;
	}
	for (i = 0; i < nv; i++)
	{
		H->nfree += !(s->fixed && s->fixed[i]);
		for (q = 0; q < s->dim; q++) { lo[q] = fmin(lo[q], x[s->dim * i + q]); hi[q] = fmax(hi[q], x[s->dim * i + q]); }
	}
	for (q = 0; q < s->dim; q++) diag2 += (hi[q] - lo[q]) * (hi[q] - lo[q]);
	/* shortcut: the shift (a mode of the bounding box's size) is a guess
	   that keeps each component invertible; tune it if the iterations grow */
	H->shift = 1 / fmax(diag2, 1e-300);
	H->aug = H->nfree == nv;
	if (!rh_chol_make(H))
		goto fail;
	/* tau: A's scale on the smoothest mode left, the free constants (with
	   pins) or a coordinate function (without) */
	{
		double m = 0, num, den;
		H->tau = 1;
		for (i = 0; i < nv; i++) H->t0[i] = H->aug ? x[s->dim * i] : 1;
		rh_mask(s, H->t0);
		if (H->aug)
		{
			for (i = 0; i < nv; i++) m += H->t0[i] / nv;
			for (i = 0; i < nv; i++) H->t0[i] -= m;
		}
		rh_work* W0 = rh_work_create(H);
		if (!W0)
			goto fail;
		H->aug = 0;   /* (A alone for tau) */
		rh_apply(H, W0, H->t0, H->t1);
		H->aug = H->nfree == nv;
		rh_work_destroy(W0);
		num = rh_dot(H->t0, H->t1, nv); den = rh_dot(H->t0, H->t0, nv);
		if (num > 0 && den > 0) H->tau = num / den;
	}
	/* z and Y together (their solves run at once), then the Schur complement */
	if (k)
		s->ops->constraints(s, x, H->phi, H->C);
	{
		double const* src[QAWS_INTERNAL_RS_HROWS + 1];
		double* dst[QAWS_INTERNAL_RS_HROWS + 1];
		src[0] = dE; dst[0] = g;
		for (r = 0; r < k; r++) { src[r + 1] = H->C + (size_t)r * N; dst[r + 1] = H->Y + (size_t)r * N; }
		if (!rh_solve_all(H, src, dst, k + 1))
			goto fail;
	}
	if (k)
	{
		for (r = 0; r < k; r++)
		{
			for (c = 0; c < k; c++) H->S[r * k + c] = rh_dot(H->C + (size_t)r * N, H->Y + (size_t)c * N, N);
			H->phi[r] = rh_dot(H->C + (size_t)r * N, g, N);
		}
		if (!qaws_internal_rs_lu(H->S, k, H->piv))
		{
			qaws_internal_rs_hsys_free(H);
			return QAWS_STATUS_NUMERICAL_FAILURE;
		}
		qaws_internal_rs_lu_solve(H->S, k, H->piv, H->phi);
		for (r = 0; r < k; r++)
			for (i = 0; i < N; i++) g[i] -= H->phi[r] * H->Y[(size_t)r * N + i];
	}
	*out = H;
	return QAWS_STATUS_OK;
fail:
	qaws_internal_rs_hsys_free(H);
	return QAWS_STATUS_NUMERICAL_FAILURE;
}

int qaws_internal_rs_hproject(qaws_internal_rs_hsys* H, double* y)
{
	rs_core const* s = H->s;
	unsigned int it, r, c, i, k = H->k, N = s->N;
	double m;
	if (k == 0)
		return 1;
	for (it = 0; it < 12; it++)
	{
		s->ops->constraints(s, y, H->phi, H->C);
		m = 0;
		for (r = 0; r < k; r++) m = fmax(m, fabs(H->phi[r]));
		if (m <= s->ctol)
			return 1;
		for (r = 0; r < k; r++)
			for (c = 0; c < k; c++) H->S[r * k + c] = rh_dot(H->C + (size_t)r * N, H->Y + (size_t)c * N, N);
		if (!qaws_internal_rs_lu(H->S, k, H->piv))
			return 0;
		qaws_internal_rs_lu_solve(H->S, k, H->piv, H->phi);
		for (r = 0; r < k; r++)
			for (i = 0; i < N; i++) y[i] -= H->phi[r] * H->Y[(size_t)r * N + i];
	}
	s->ops->constraints(s, y, H->phi, NULL);
	m = 0;
	for (r = 0; r < k; r++) m = fmax(m, fabs(H->phi[r]));
	return m <= s->ctol;
}
