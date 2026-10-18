#ifndef QAWS_QUAD_IMPROVE_H
#define QAWS_QUAD_IMPROVE_H

/*
 * Quad mesh finishing: valence pairs cancelled by local quad operations
 * against a feature-aware ideal valence, then relaxation by angles guarded
 * by the scaled Jacobian, on the reference.
 */

#include "qaws_mesh_internal.h"

/* valence pairs cancelled and quads relaxed by shape on the reference; kind: per vertex 0 free, 1 on a line, 2 corner (updated);
   compacts m at the end; returns the operations applied, HM_NONE on allocation failure */
unsigned int qi_improve(qaws_hm* m, qaws_mesh_reference const* r, mr_scratch* s, unsigned char* kind, double rho, unsigned int rounds);

#endif /* QAWS_QUAD_IMPROVE_H */
