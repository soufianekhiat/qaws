#ifndef QAWS_CROSS_FIELD_H
#define QAWS_CROSS_FIELD_H

/*
 * Globally smoothest n-direction fields (Knoppel, Crane, Pinkall, Schroder
 * 2013) on a store's vertices: per vertex a complex u whose argument is n
 * times the direction's angle, angles measured from the vertex's first edge
 * and rescaled to sum to 2 pi round it (pi on the border), neighbours
 * compared through their edge's transport (the s = 0 connection: eq. 18
 * without its holonomy correction, fine on the dense meshes it serves). The
 * energy sum w |u_j - r_ij u_i|^2 (cotangent weights) is minimized with the
 * vertices along feature lines and the border held to their lines, corners
 * free; one sparse LDL^T factor serves the smoothest field (per connected
 * component one solve when lines hold it, else inverse power iteration)
 * and the one drawn to a guide (lambda_t = 0).
 */

#include "qaws_mesh_internal.h"

typedef struct cf_system cf_system;

/* the n-direction connection Laplacian of m's vertices (unit normals nrm), the held values from L0->con & MF_ALONG (but
   not MF_PINNED: corners free) along m's feature and border edges, factored once */
qaws_status cf_create(qaws_hm const* m, double const* nrm, mf_level const* L0, unsigned int n, cf_system** out);
/* the smoothest field: 3 per vertex, unit tangents (for n = 4 in the mf_level q convention) */
qaws_status cf_smoothest(cf_system* sys, double* out_q);
/* the field drawn to guide (3 per vertex, a direction scaled by its strength; 0: none) */
qaws_status cf_aligned(cf_system* sys, double const* guide, double* out_q);
/* the last field, on the components with held vertices, relaxed toward Ginzburg-Landau's at the scale eps (at most steps
   MBO steps, a second factor, kept for further calls at the same eps): its coinciding zeros parted */
qaws_status cf_relax(cf_system* sys, double eps, unsigned int steps, double* out_q);
/* per face the index of the last field computed, in 1/n turns */
void cf_index(cf_system const* sys, int* out_face_index);
/* the factor's entries and the lower triangle's of its matrix */
void cf_fill(cf_system const* sys, double* out_factor, double* out_matrix);
void cf_free(cf_system* sys);
/* level 0's q carried up the levels (each coarse vertex its children's, area-weighted); vertices along lines keep cq */
void cf_restrict(mf_hier* h);
/* qaws_mesh_3d_direction_field's SMOOTHEST: h from mf_build on the mesh, n its symmetry; aligned draws it to the
   principal curvature directions of the mesh's PN reference */
qaws_status cf_direction_field(qaws_mesh_3d const* mesh, mf_hier const* h, unsigned int n, int aligned, qaws_scalar* out_directions,
	int* out_face_index);

#endif /* QAWS_CROSS_FIELD_H */
