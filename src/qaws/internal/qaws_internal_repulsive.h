#ifndef QAWS_INTERNAL_REPULSIVE_H
#define QAWS_INTERNAL_REPULSIVE_H

#include "../qaws_repulsive.h"
#include <stddef.h>

/* the energy layer's state and its double-precision entry points, shared
   with the solver (qaws_repulsive_solve.c) */

struct qaws_repulsive_curves
{
	unsigned int dim, nv, ne;
	unsigned int* edges;
	double alpha, beta;
	qaws_batch_executor const* executor;
};

struct qaws_repulsive_mesh
{
	unsigned int nv, nf;
	unsigned int* tris;
	double p;
	qaws_batch_executor const* executor;
};

/* the mesh's energy at x (3 per vertex); grad (may be NULL) += its gradient */
qaws_status qaws_internal_repulsive_mesh_energy(qaws_repulsive_mesh const* mesh, double const* x, double* out_energy, double* grad);

/* the energy at x (dim per vertex); grad (may be NULL) += its gradient */
qaws_status qaws_internal_repulsive_curves_energy(qaws_repulsive_curves const* curves, double const* x, double* out_energy, double* grad);

/* the surfaces' metric A = B + B0 at x (vertex_count^2, row major) */
qaws_status qaws_internal_repulsive_mesh_metric(qaws_repulsive_mesh const* mesh, double const* x, double* A);

/* A u = b (b mean-zero; one coordinate) by GMRES to the relative residual
   tol, right-preconditioned by (-L)^+ L^(2 - s) (-L)^+ (L the cotan
   Laplacian) when precondition; u mean-zero; the iterations taken */
qaws_status qaws_internal_repulsive_mesh_gmres(qaws_repulsive_mesh const* mesh, double const* x, double const* b, int precondition, double tol,
	double* u, unsigned int* out_iterations);

/* 1 when along x0 -> x1 every vertex and face it is not a corner of, and
   every two edges sharing no vertex, stay at least min(clearance, half
   their distance at x0) apart */
int qaws_internal_repulsive_mesh_path_clear(qaws_repulsive_mesh const* mesh, double clearance, double const* x0, double const* x1);

/* the fractional Sobolev metric A = B + B0 at x (vertex_count^2, row major) */
qaws_status qaws_internal_repulsive_curves_metric(qaws_repulsive_curves const* curves, double const* x, double* A);

/* 1 when along the straight motion x0 -> x1 every pair of edges sharing no
   vertex stays at least min(clearance, half its distance at x0) apart */
int qaws_internal_repulsive_curves_path_clear(qaws_repulsive_curves const* curves, double clearance, double const* x0, double const* x1);

/* hv += H v at x (in double) */
qaws_status qaws_internal_repulsive_curves_hvp(qaws_repulsive_curves const* curves, double const* x, double const* v, double* hv);
qaws_status qaws_internal_repulsive_mesh_hvp(qaws_repulsive_mesh const* mesh, double const* x, double const* v, double* hv);

/* Barnes-Hut energy at x and, if grad, grad += its gradient (in double) */
qaws_status qaws_internal_repulsive_curves_bh(qaws_repulsive_curves const* curves, double const* x, double theta, double* E, double* grad);
qaws_status qaws_internal_repulsive_mesh_bh(qaws_repulsive_mesh const* mesh, double const* x, double theta, double* E, double* grad);

/* Hierarchical operators on vertex functions (qaws_repulsive_bh.c): the
   metric A = B + B0 at x, and the fractional operator sum over ordered
   pairs of distinct points of (u_i - u_j)(v_i - v_j) m_i m_j r^-ex; theta 0
   is the dense form. Built once at x, applied as y = A u (nv each). */
typedef struct qaws_internal_rbh_op qaws_internal_rbh_op;
qaws_status qaws_internal_rbh_curves_metric(qaws_repulsive_curves const* curves, double const* x, double theta, qaws_internal_rbh_op** out);
qaws_status qaws_internal_rbh_mesh_metric(qaws_repulsive_mesh const* mesh, double const* x, double theta, qaws_internal_rbh_op** out);
qaws_status qaws_internal_rbh_fractional(double const* points, double const* mass, unsigned int count, unsigned int dim, double ex, double theta,
	qaws_internal_rbh_op** out, qaws_batch_executor const* executor);
void qaws_internal_rbh_apply(qaws_internal_rbh_op* op, double const* u, double* y);
/* a product's scratch: products with one each may run at once on one operator */
typedef struct qaws_internal_rbh_scratch
{
	double *gi, *ui, *G, *U, *accB, *acc0;
} qaws_internal_rbh_scratch;
qaws_internal_rbh_scratch* qaws_internal_rbh_scratch_create(qaws_internal_rbh_op const* op);
void qaws_internal_rbh_scratch_destroy(qaws_internal_rbh_scratch* scratch);
void qaws_internal_rbh_apply_with(qaws_internal_rbh_op const* op, qaws_internal_rbh_scratch* scratch, double const* u, double* y);
/* the near links (and far, if asked) */
size_t qaws_internal_rbh_links(qaws_internal_rbh_op const* op, size_t* out_far);
void qaws_internal_rbh_destroy(qaws_internal_rbh_op* op);

/* every pair i < j of overlapping boxes that accept allows (NULL: all), visited once, by a
   self traversal of a median-split tree; a non-OK visit stops and is returned */
struct qaws_bp_box;
qaws_status qaws_internal_rbh_box_pairs(struct qaws_bp_box const* boxes, unsigned int count, unsigned int dim,
	int (*accept)(void* user, unsigned int i, unsigned int j), qaws_status (*visit)(void* user, unsigned int i, unsigned int j), void* user);

#endif /* QAWS_INTERNAL_REPULSIVE_H */
