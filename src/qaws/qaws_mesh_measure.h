#ifndef QAWS_MESH_MEASURE_H
#define QAWS_MESH_MEASURE_H

#include "qaws_mesh_reference.h"

/*
 * Measures of surface meshes in the forms remeshing papers report them:
 * shape (Hu et al. 2017 and Wang et al. 2019 for triangles, QuadWild for
 * quads) and the distance to a reference, bounded from both sides.
 *
 * Shape, one pass over the faces:
 *   quality     a triangle's Q = 2 sqrt(3) r / h_max (inradius over the
 *               longest edge, scaled to 1 when equilateral);
 *   angles      every face corner's angle; the shares under 30 and over 90
 *               degrees and the mean smallest angle over triangles;
 *   Jacobian    a quad's smallest corner scaled Jacobian (the corner's two
 *               edges' cross product along the quad's Newell normal, over
 *               their lengths: 1 for a square, negative where it folds);
 *   deviations  AD, the quads' corners' mean |angle - 90|; ED, the edges'
 *               standard deviation over their mean length;
 *   valence     a vertex's ideal valence from its corners: a smooth interior
 *               vertex (off the border and the features) wants 6 when
 *               triangles are the most, else 4; a border, crease or corner
 *               vertex wants as many faces as its corners' angles hold ideal
 *               corners (60 degrees for a triangle, 90 for a quad), within a
 *               half (a tie at a half rounds to the vertex's own count, one
 *               face at the least), plus 1 on the border; a vertex is irregular
 *               when its valence is not that ideal. A box's corners (3 quads),
 *               its edges (4), a quad grid's straight border (3 edges) and
 *               its corners (2) are all regular.
 *
 * Deviation, certified: the sup over the mesh of the distance to the
 * reference and the sup over the reference of the distance to the mesh, each
 * as a lower and an upper bound (branch and bound over the faces and the
 * patches: upper bounds from the Bezier nets of the difference of the two
 * surfaces over matched triangles, never from footpoints). The upper bounds
 * and the reference-to-mesh lower bound hold unconditionally; the
 * mesh-to-reference lower bound is a distance found by projection, exact on
 * LINEAR references and, on PN, up to the patch projection (the
 * reference's, then polished along the patches' sides and across them).
 *
 * Polygons are measured as their corner fans: a non-planar quad counts as
 * its two triangles (from its first vertex).
 */

typedef struct qaws_shape_quality
{
	qaws_scalar min_quality, mean_quality;   /* triangles: Q = 2 sqrt(3) r / h_max, 1 when equilateral */
	qaws_scalar min_angle, max_angle;        /* degrees, over every corner of every face */
	qaws_scalar mean_min_angle;              /* degrees, the mean of each triangle's smallest angle */
	qaws_scalar below_30, above_90;          /* shares of triangle corners */
	qaws_scalar min_jacobian, mean_jacobian; /* quads: per quad its smallest corner scaled Jacobian, then min and mean over quads */
	qaws_scalar angle_deviation;             /* quads: mean |corner angle - 90|, degrees (AD) */
	qaws_scalar edge_deviation;              /* std / mean over the unique edges (ED) */
	qaws_scalar regular_share;               /* triangle-majority meshes: share of smooth interior vertices of valence 5..7 (V567);
	                                            quad-majority: share of all vertices at their ideal valence */
	unsigned int irregular;                  /* vertices off their ideal valence, border, crease and corner vertices included */
	unsigned int valence_excess;             /* sum of |valence - rounded ideal| */
} qaws_shape_quality;

/* out_ideal (may be NULL): per vertex its rounded ideal valence (0 for an isolated vertex); irregular counts the vertices off it */
qaws_status qaws_mesh_3d_measure_shape(qaws_shape_quality* out_quality, unsigned int* out_ideal, qaws_mesh_3d const* mesh);

typedef struct qaws_deviation_desc
{
	qaws_scalar precision;   /* target gap (upper - lower) as a share of the reference's box diagonal; 0: 1e-4 */
	unsigned int budget;     /* subdivisions per face or patch past those that bring it down to the other side's mean
	                            triangle size (2 x their areas' ratio, the memory growing with it); 0: 256 */
} qaws_deviation_desc;

typedef struct qaws_mesh_deviation
{
	qaws_scalar to_reference[2];      /* lower, upper bound of the sup over the mesh of the distance to the reference */
	qaws_scalar from_reference[2];    /* the same, from the reference to the mesh */
	qaws_scalar hausdorff[2];         /* the larger of the two directions' bounds */
	qaws_scalar rms;                  /* mesh to reference, area-weighted (projected, not certified) */
	qaws_scalar diagonal;             /* the reference's bounding-box diagonal */
	qaws_scalar to_reference_at[3];   /* the point of the mesh where to_reference[0] is reached */
	qaws_scalar from_reference_at[3]; /* the point of the reference where from_reference[0] is reached */
	int converged;                    /* 1: both gaps within the precision; 0: a budget ran out (the bounds hold, only wider) */
} qaws_mesh_deviation;

/* desc NULL: the defaults. INVALID_ARGUMENT for a NULL mesh, reference or
   out, a mesh without faces, a reference without triangles, a negative
   precision. */
qaws_status qaws_mesh_3d_measure_deviation(qaws_mesh_deviation* out_deviation, qaws_mesh_reference const* reference, qaws_mesh_3d const* mesh,
	qaws_deviation_desc const* desc);

#endif /* QAWS_MESH_MEASURE_H */
