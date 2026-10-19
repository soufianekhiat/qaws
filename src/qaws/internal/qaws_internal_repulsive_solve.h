#ifndef QAWS_INTERNAL_REPULSIVE_SOLVE_H
#define QAWS_INTERNAL_REPULSIVE_SOLVE_H

#include "../qaws_repulsive.h"

/* The fractional Sobolev descent shared by curves and meshes: a core
   driven by four callbacks (qaws_repulsive_solve.c). A solver object
   starts with the core, so a pointer to it is a pointer to the core. */

typedef struct qaws_internal_rs_core qaws_internal_rs_core;

/* the hierarchical mode's operators at x (R4): the metric, the fractional
   operator of the preconditioner (-L)^-1 L^(2 - s) (-L)^-1, and -L as
   weighted vertex pairs (sum w (u_i - u_j)(v_i - v_j)) with lumped masses */
typedef struct qaws_internal_rs_hier
{
	struct qaws_internal_rbh_op *A, *F;
	unsigned int* lij;                   /* nl pairs */
	double *lw, *mass;
	unsigned int nl;
} qaws_internal_rs_hier;

typedef struct qaws_internal_rs_ops
{
	/* E at x, and grad (N, may be NULL) += its gradient */
	qaws_status (*energy)(qaws_internal_rs_core const* core, double const* x, double* E, double* grad);
	/* the metric A at x: nv x nv (A~ = A on each coordinate), or N x N if core->full */
	qaws_status (*metric)(qaws_internal_rs_core const* core, double const* x, double* A);
	/* phi (k) and, if C, the Jacobian (k x N, row major) at x */
	void (*constraints)(qaws_internal_rs_core const* core, double const* x, double* phi, double* C);
	/* 1 when the straight motion x0 -> x1 keeps every checked pair clear */
	int (*path_ok)(qaws_internal_rs_core const* core, double const* x0, double const* x1);
	/* the hierarchical mode's operators at x (theta > 0) */
	qaws_status (*hier)(qaws_internal_rs_core const* core, double const* x, qaws_internal_rs_hier* out);
	/* how far vertex v moves in the embedding for the change d at x (NULL: |d_v|) */
	double (*move)(qaws_internal_rs_core const* core, double const* x, double const* d, unsigned int v);
	/* for the solution's derivatives (NULL: not provided): hv += H v (the
	   energy's Hessian), and hv += sum_r lambda_r Hess(c_r) v */
	qaws_status (*hvp)(qaws_internal_rs_core const* core, double const* x, double const* v, double* hv);
	void (*constraints_hvp)(qaws_internal_rs_core const* core, double const* x, double const* lambda, double const* v, double* hv);
	/* optional: g += the translations added back after the direction's solve
	   (Repulsive Surfaces, Fig. 8: each group solved with its barycenter
	   held, then moved by its mean of dE) */
	void (*translate)(qaws_internal_rs_core const* core, double const* x, double const* dE, double* g);
} qaws_internal_rs_ops;

struct qaws_internal_rs_core
{
	qaws_internal_rs_ops const* ops;
	unsigned int dim, nv, N, k;          /* N = dim nv; k constraint rows (hierarchical: pins not among them) */
	int checked;
	int full;                            /* the metric is N x N (parameters: not one per coordinate) */
	double clearance, tol, ctol, lmean, h;
	double hmax;                         /* > 0: the largest vertex move a step may take (0: none) */
	double theta;                        /* 0: dense; > 0: hierarchical, pins eliminated */
	qaws_batch_executor const* executor; /* optional: the hierarchical solves run through it */
	int translates;                      /* nothing pinned, no barycenter held: a free translation */
	unsigned char const* fixed;          /* pinned vertices */
	double* x;                           /* the current positions, in double */
	qaws_repulsive_report rep;
	qaws_repulsive_progress progress;    /* the desc's: cancellation, per-step reports */
};

/* the host's cancel callback says stop */
#define QAWS_INTERNAL_RS_CANCELLED(core) ((core)->progress.cancel && (core)->progress.cancel((core)->progress.user_data))

/* steps until a termination reason or max steps, reporting each accepted one (the *_solve loop) */
qaws_status qaws_internal_rs_solve(qaws_internal_rs_core* core, unsigned int max_steps, qaws_scalar* positions, qaws_repulsive_report* out_report);

/* the objective at positions and its gradient (may be NULL) */
qaws_status qaws_internal_rs_objective(qaws_internal_rs_core const* core, qaws_scalar const* positions, qaws_scalar* out_energy,
	qaws_scalar* out_gradient);

/* the most constraint rows (pins aside) the hierarchical mode solves for */
#define QAWS_INTERNAL_RS_HROWS 32

void qaws_internal_rs_hier_free(qaws_internal_rs_hier* h);

/* the hierarchical direction g (N) at x from dE (qaws_repulsive_hsolve.c),
   its state kept for the projection; the projection of y back onto the
   constraints (1 when within ctol) */
typedef struct qaws_internal_rs_hsys qaws_internal_rs_hsys;
qaws_status qaws_internal_rs_hdirection(qaws_internal_rs_core const* core, double const* x, double const* dE, double* g, qaws_internal_rs_hsys** out);
int qaws_internal_rs_hproject(qaws_internal_rs_hsys* H, double* y);
void qaws_internal_rs_hsys_free(qaws_internal_rs_hsys* H);

/* the converged solution's derivative (qaws_repulsive_sens.c), dense mode:
   with drows (k: the change of each constraint row's target, C dx = drows)
   dx (N); or with xbar (N) the rows' adjoints (k) in out. *out_gauge: the
   free rotations removed (dx orthogonal to them). */
qaws_status qaws_internal_rs_sensitivity(qaws_internal_rs_core const* core, double const* drows, double const* xbar, double* out, unsigned int* out_gauge);

/* LU with partial pivoting of the n x n row-major m in place (0: singular),
   and its solve */
int qaws_internal_rs_lu(double* m, unsigned int n, unsigned int* piv);
void qaws_internal_rs_lu_solve(double const* m, unsigned int n, unsigned int const* piv, double* b);

/* the descent direction g at positions (N) */
qaws_status qaws_internal_rs_direction(qaws_internal_rs_core const* core, qaws_scalar const* positions, qaws_scalar* out_direction);

/* one step (positions in place) */
qaws_status qaws_internal_rs_step(qaws_internal_rs_core* core, qaws_scalar* positions, qaws_repulsive_report* out_report);

/* a room's walls on nv vertices (dim each; u and v the plane's coordinates;
   the polygon's count points, 2 each): E += w sum over vertices and walls of
   the integral along the wall of rho^-2, grad (may be NULL) += its
   gradient; 0 (nothing added) when a vertex is not strictly inside */
int qaws_internal_rs_room(double const* room, unsigned int count, double w, unsigned int dim, unsigned int u, unsigned int v,
	double const* x, unsigned int nv, double* E, double* grad);

/* the distance between segments p0 p1 and q0 q1 (dim 2 or 3) */
double qaws_internal_rs_segdist(double const* p0, double const* p1, double const* q0, double const* q1, unsigned int dim);

/* conservative advancement for one pair: dist(user, t) its distance at time
   t in [0, 1], d0 = dist(user, 0), L the Lipschitz constant (the sum of the
   pair's largest point moves); 1 when it stays at least req apart */
int qaws_internal_rs_advance(double (*dist)(void const* user, double t), void const* user, double d0, double L, double req);

#endif /* QAWS_INTERNAL_REPULSIVE_SOLVE_H */
