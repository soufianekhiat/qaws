#ifndef QAWS_REMESH_H
#define QAWS_REMESH_H

#include "qaws_mesh_reference.h"
#include "qaws_mesh_field.h"

/*
 * Surface remeshing: a new triangle mesh of the same surface, its edges
 * near a target length (uniform, or adaptive by a size field), its
 * vertices on the reference.
 *
 *   LOCAL  local operations (Botsch and Kobbelt 2004), repeated: edges
 *          longer than 4/3 of their target split at their midpoint on the
 *          reference; edges shorter than 4/5 collapsed when the link
 *          condition holds, no edge longer than 4/3 appears and no face
 *          turns over; edges flipped toward valence 6 (4 on the border);
 *          vertices moved toward their neighbours' centroid along the
 *          surface (walked on the reference, no projection). Features kept:
 *          feature corners never move, vertices on feature lines (and the
 *          border) slide along them, feature edges are never flipped or
 *          crossed (a feature line turning by more than 60 degrees has a
 *          corner there). The default.
 *   CVT    a centroidal Voronoi tessellation of the reference (Lloyd,
 *          then L-BFGS; Yan, Levy et al. 2009, Liu et al. 2009): LOCAL
 *          first (10 iterations), the vertex count set from the size (a
 *          hexagon of the target each), the vertices on feature lines
 *          spread along them and fixed, the others moved to the centroids
 *          of their Voronoi cells restricted to the reference (density
 *          size^-4, region by region: no cell crosses a crease), the edges
 *          kept Delaunay by flips, no move folding a face; solved on a
 *          coarse level first (a quarter of the vertices, then each
 *          triangle split in four) unless the surface bends tighter than
 *          twice the size.
 *          `iterations` counts the L-BFGS steps after 5 Lloyd steps, per
 *          level (0: 30). The most uniform and regular triangles, the
 *          edges on target; about five times slower than LOCAL.
 *   FIELD_ALIGNED  reserved (UNSUPPORTED_OPERATION today).
 *
 * The reference is the surface: by default the input's own (PN, its
 * features as creases). Every output vertex keeps its location on the
 * reference (qaws_mesh_3d_get_footpoints), and every attribute of the
 * reference's mesh is evaluated there by its rule: vertex attributes at the
 * vertex, corner attributes at the vertex as seen from the face's side (a
 * seam that is a feature keeps both sides), face attributes at the face's
 * centre. The input is not changed.
 */

typedef enum qaws_remesh_method
{
	QAWS_REMESH_LOCAL = 0,
	QAWS_REMESH_CVT,
	QAWS_REMESH_FIELD_ALIGNED
} qaws_remesh_method;

typedef struct qaws_remesh_desc
{
	qaws_remesh_method method;
	qaws_mesh_reference const* reference;   /* NULL: the input's own, PN */
	qaws_scalar target_length;              /* 0: the input's mean edge length */
	qaws_size_field const* size;            /* adaptive, over target_length; on `reference` when sampled on one */
	unsigned int iterations;                /* 0: 10 */
} qaws_remesh_desc;

/* out_report may be NULL */
qaws_status qaws_mesh_remesh(qaws_mesh_3d const* mesh, qaws_remesh_desc const* desc, qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report);

#endif /* QAWS_REMESH_H */
