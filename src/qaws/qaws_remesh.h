#ifndef QAWS_REMESH_H
#define QAWS_REMESH_H

#include "qaws_mesh_reference.h"
#include "qaws_mesh_field.h"

/*
 * Surface remeshing: a new triangle mesh of the same surface, its edges
 * near a target length (uniform, or adaptive by a size field), its
 * vertices on the reference.
 *
 *   LOCAL  local operations (Botsch and Kobbelt 2004, as pmp,
 *          geometry-central and Mmgs refine them), repeated until nothing
 *          changes: edges longer than 4/3 of their target split at their
 *          midpoint on the reference; edges shorter than 4/5 collapsed onto
 *          their midpoint on the reference (on the crease along a feature
 *          line; else at the kept end) when the link condition holds, no
 *          edge longer than 4/3 appears, no face turns over and no corner
 *          opens past 120 degrees or past the largest there was; edges
 *          flipped toward each vertex's ideal valence (6 inside; on
 *          features and the border, in each sector between two feature
 *          edges as many faces as the nearest whole number of 60-degree
 *          corners: 3 per side of a smooth crease or a straight border)
 *          when the new pair's smallest angle stays above the old pair's
 *          or 30 degrees; vertices moved to their ring's area
 *          centroid weighted by the size field (a graded mesh keeps its
 *          grading), walked on the reference (three sweeps in the last
 *          third of the iterations with a size field, angle bounds or a
 *          tolerance, one otherwise: INSTANT's dense mesh keeps fewer
 *          singularities with one); needles (under 15 degrees)
 *          collapsed, caps (over 120) flipped; Delaunay flips to finish.
 *          Features kept: feature corners never move, vertices on feature
 *          lines (and the border) slide along them, feature edges are never
 *          flipped or crossed (a feature line turning by more than 60
 *          degrees has a corner there). The default.
 *
 *          min_angle, max_angle: then an angle-bounded pass (Wang, Yan et
 *          al. 2019; their defaults 35 and 86 degrees), in rounds: the
 *          worst large angles cut by a flip, moved away (a pattern search
 *          over the face's corners) or split away (the longest edge split,
 *          the flip round the new vertex that settles best, relaxation),
 *          the worst small ones collapsed, flipped or moved away, valence
 *          flips and relaxation between; each operation kept only when the
 *          worst angle of its star gets better and no angle inside the
 *          bounds leaves them; until every angle is inside or it stops
 *          getting better (the report's quality says how far it got).
 *          INVALID_ARGUMENT for a negative value, min_angle of 60 degrees
 *          or more, max_angle of 60 or less, min_angle >= max_angle.
 *          Exempt: at a border, crease or corner vertex, the corners in a
 *          sector (between two feature or border edges) of angle theta when
 *          no whole number k >= 1 has theta / k inside the bounds (a
 *          20-degree wedge; not a cube's 90-degree corner, 2 corners of 45).
 *          tolerance (as Mmgs's hausd): the sizes shrink where the surface
 *          bends (an edge's chord gap at most 3/4 of it); no operation
 *          leaves a face whose centre or edge midpoints are farther than it
 *          from the reference (unless a face it replaces was farther
 *          already, in the loop), edges whose midpoint is farther are split
 *          whatever their length; after the loop the guard samples the
 *          edges' quarter points too, faces farther have their longest
 *          edge split, and every later operation (Delaunay flips, the
 *          angle pass) is held to it. Sampled only: a face may bulge
 *          farther between its samples. It wins over the angle bounds.
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
 * Pinned vertices (`pinned`, indices of the input's vertices) are kept as
 * feature corners are: never moved or collapsed away, at their input
 * positions, and they come out flagged as corners
 * (qaws_mesh_3d_is_feature_corner). INVALID_ARGUMENT for an index past the
 * input's vertices.
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
	qaws_scalar min_angle;                  /* radians, below 60 degrees; > 0: the angle-bounded pass lifts angles to it */
	qaws_scalar max_angle;                  /* radians, above 60 degrees; > 0: the angle-bounded pass brings angles under it */
	qaws_scalar tolerance;                  /* > 0: the distance from the reference no face's sampled points pass */
	unsigned int const* pinned;             /* input vertices kept where they are, as feature corners (LOCAL and CVT) */
	unsigned int pinned_count;
} qaws_remesh_desc;

/* out_report may be NULL */
qaws_status qaws_mesh_remesh(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* mesh, qaws_remesh_desc const* desc);

#endif /* QAWS_REMESH_H */
