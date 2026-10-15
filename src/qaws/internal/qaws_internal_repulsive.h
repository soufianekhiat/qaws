#ifndef QAWS_INTERNAL_REPULSIVE_H
#define QAWS_INTERNAL_REPULSIVE_H

#include "../qaws_repulsive.h"

/* the energy layer's state and its double-precision entry points, shared
   with the solver (qaws_repulsive_solve.c) */

struct qaws_repulsive_curves
{
	unsigned int dim, nv, ne;
	unsigned int* edges;
	double alpha, beta;
};

/* the energy at x (dim per vertex); grad (may be NULL) += its gradient */
qaws_status qaws_internal_repulsive_curves_energy(qaws_repulsive_curves const* curves, double const* x, double* out_energy, double* grad);

/* the fractional Sobolev metric A = B + B0 at x (vertex_count^2, row major) */
qaws_status qaws_internal_repulsive_curves_metric(qaws_repulsive_curves const* curves, double const* x, double* A);

/* 1 when along the straight motion x0 -> x1 every pair of edges sharing no
   vertex stays at least min(clearance, half its distance at x0) apart */
int qaws_internal_repulsive_curves_path_clear(qaws_repulsive_curves const* curves, double clearance, double const* x0, double const* x1);

#endif /* QAWS_INTERNAL_REPULSIVE_H */
