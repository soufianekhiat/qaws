#ifndef QAWS_REMESH_CVT_H
#define QAWS_REMESH_CVT_H

/* Remeshing by a centroidal Voronoi tessellation of the reference (QAWS_REMESH_CVT) */

#include "qaws_mesh_internal.h"

/* the working mesh's free vertices moved to a centroidal Voronoi tessellation of the reference (density size^-4), its
   creases' vertices first spread along them and locked, its edges kept Delaunay by flips; 0 on allocation failure */
int rv_optimize(qaws_hm* m, qaws_mesh_reference const* r, qaws_size_field const* size, double L, unsigned int size_attr, unsigned int region_attr,
	mr_scratch* s, unsigned int steps, qaws_mesh_report* rep);

/* what the last optimisation reached, written when rv_trace is set (the tests): the energy after the Lloyd steps and at
   the end, the energy evaluations; shortcut: a global probe, not for threads; when remeshing runs in threads, report the
   energies through qaws_mesh_report */
typedef struct rv_probe
{
	double lloyd_energy, energy;
	unsigned int evaluations;
} rv_probe;

extern rv_probe* rv_trace;

#endif /* QAWS_REMESH_CVT_H */
