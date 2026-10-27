#ifndef QAWS_INTERNAL_ARRAY_H
#define QAWS_INTERNAL_ARRAY_H

#include "../qaws_mesh_view.h"

/* element i's component c of a host array, read bytewise (qaws_view.c) */
double qaws_internal_array_get(qaws_array_view const* v, unsigned int i, unsigned int c);

/* triangle t's corner k of a mesh view (qaws_mesh_view.c) */
unsigned int qaws_internal_mesh_view_index(qaws_mesh_view const* v, unsigned int t, unsigned int k);

/* the view's positions as doubles into x (3 per vertex): 0 when one is not finite (qaws_mesh_view.c) */
int qaws_internal_mesh_view_positions(double* x, qaws_array_view const* positions);

#endif /* QAWS_INTERNAL_ARRAY_H */
