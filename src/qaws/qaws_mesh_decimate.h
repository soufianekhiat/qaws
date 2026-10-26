#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_DECIMATE_H_F32
#define QAWS_MESH_DECIMATE_H_F32
#define QAWS_MESH_DECIMATE_H_BODY
#endif
#else
#ifndef QAWS_MESH_DECIMATE_H_F64
#define QAWS_MESH_DECIMATE_H_F64
#define QAWS_MESH_DECIMATE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_DECIMATE_H_BODY
#undef QAWS_MESH_DECIMATE_H_BODY
#define QAWS_MESH_DECIMATE_H

#include "qaws_mesh_measure.h"

/*
 * Decimation by quadric error metrics (Garland, Heckbert 1997), in place:
 * edges collapsed cheapest first, each merged vertex placed where the sum
 * of its two ends' quadrics is least (the planes of the faces they had,
 * weighted by area). The error is read as a distance: the root of the
 * area-weighted mean squared distance to those planes.
 *
 * What stays:
 *   features   the mesh's feature lines and border, held near their place
 *              by planes standing on their edges, and its feature corners;
 *              a vertex on a line moves only along it. With crease_angle,
 *              the edges whose faces turn by more are added to the mesh's
 *              features first (and stay so), their corners found as
 *              qaws_mesh_3d_set_feature_edges finds them.
 *   silhouettes  an edge whose two faces turn opposite ways to a view is on
 *              its silhouette. For the given views (directions, or eye points
 *              with QAWS_DECIMATE_VIEW_POINTS) the plane through the edge and
 *              the view's ray is added to its ends' quadrics: moving in that
 *              plane keeps the outline, moving out of it costs (averaged
 *              over the views: one or a few views hold their outlines
 *              closely).
 *   views      octahedral_views^2 directions toward eyes all round, on an
 *              octahedral map's grid (as octahedral impostors sample the
 *              sphere; from the upper +z hemisphere only with
 *              QAWS_DECIMATE_HEMISPHERE): each face's plane weighs the share
 *              of the views it faces (a fiftieth at least), and 1 +
 *              silhouette_weight times its corner's share of the views it
 *              is on a silhouette of. Parts no view sees go first; the
 *              outlines get the faces. (Over the whole sphere every face
 *              faces half the views and every vertex is on some outline:
 *              the result is near plain decimation's. An outline seen along
 *              d is the one seen along -d, so for outlines alone the upper
 *              half is the whole sphere: the faces' facing makes the
 *              difference.)
 *   importance  per vertex, a factor on its quadric (above 1: kept
 *              longer); 0 or less: never removed nor moved. The border's
 *              vertices are held still with QAWS_DECIMATE_LOCK_BORDER.
 *   shape      a collapse turning a face by more than max_normal_angle is
 *              refused.
 * attributes  with attribute_weight, the per vertex LINEAR and NORMALIZED
 *              attributes join the quadric (Hoppe 1999): over each triangle
 *              each component is the linear function through its corners'
 *              values, and a vertex is charged its squared difference from
 *              it, a unit of attribute counting as attribute_weight of
 *              distance. A collapse places the position and the attributes
 *              together where the sum is least; collapses that would
 *              smear an attribute's change cost, so its edges are kept.
 *              The error then reads in the same units.
 * memoryless  with QAWS_DECIMATE_MEMORYLESS nothing is stored per vertex
 *              (Lindstrom and Turk 1998): each collapse is costed from the
 *              faces round its edge as they are. The volume under them is
 *              kept (one linear constraint on the new point), and the point
 *              is where the volumes the moving faces sweep, squared and
 *              summed, are least (with the border's and features' planes):
 *              a closed mesh keeps its volume to round-off. Where a collapse
 *              touches the border, Lindstrom and Turk's border constraint
 *              keeps the chain's area vector (sum of v_i x v_i+1) too: a
 *              planar region keeps its area. Views and
 *              attribute_weight are refused with it (they live in the
 *              stored quadrics); importance only locks.
 * Attributes follow the collapses by their rules (qaws_mesh.h); corner
 * attributes (of up to 4 components) by wedges (Hoppe): the corners round the
 * edge grouped (equal values at a vertex, to round-off, joined across the
 * edge's faces), each group's value the one at the new position on the
 * nearest of its corners' triangles before the collapse (LINEAR, VECTOR and
 * the unit rules interpolated, the others its nearest corner's). Their seams
 * stay: an inside edge whose faces' corner attributes differ at either end
 * (UVs cut apart, split normals) or whose faces' attributes differ
 * (materials) is made a feature first, kept as the lines are.
 *
 * It stops at the target face count, before a collapse whose error passes
 * max_error, or when no collapse is left. With QAWS_DECIMATE_MEASURE the
 * report holds the result's certified deviation from the input
 * (qaws_mesh_3d_measure_deviation against the input as a LINEAR
 * reference).
 */

#define QAWS_DECIMATE_LOCK_BORDER 1u        /* the border's vertices are neither removed nor moved */
#define QAWS_DECIMATE_MEASURE 2u            /* the report's deviation measured */
#define QAWS_DECIMATE_VIEW_POINTS 4u        /* views are eye points (perspective), not directions */
#define QAWS_DECIMATE_HEMISPHERE 8u         /* the octahedral views over the +z hemisphere only */
#define QAWS_DECIMATE_MEMORYLESS 16u        /* costs from the current faces, nothing stored (Lindstrom and Turk 1998) */

typedef struct qaws_decimate_desc
{
	unsigned int target_faces;          /* stop at this many faces; 0: no count */
	qaws_scalar max_error;              /* stop before an error past this distance; 0: no bound */
	unsigned int flags;                 /* QAWS_DECIMATE_* */
	qaws_scalar crease_angle;           /* radians: sharper edges become features; 0: the mesh's features only */
	qaws_scalar max_normal_angle;       /* radians: a collapse turning a face by more is refused; 0: 1.37 (about 78 degrees) */
	qaws_scalar const* importance;      /* per vertex; NULL: 1 */
	qaws_scalar const* views;           /* 3 each: directions, or eye points with QAWS_DECIMATE_VIEW_POINTS (the octahedral ones directions) */
	unsigned int view_count;
	unsigned int octahedral_views;      /* N: N x N directions on an octahedral map; 0: none */
	qaws_scalar silhouette_weight;      /* the given views' planes against the faces' (per squared edge length), the octahedral
	                                       views' factor; 0: 100 */
	qaws_scalar feature_weight;         /* the border's and features' planes, the same way; 0: 100 */
	qaws_scalar attribute_weight;       /* a unit of the vertex attributes (LINEAR, NORMALIZED) as this distance in the quadric
	                                       (Hoppe 1999); 0: the attributes not in it */
} qaws_decimate_desc;

typedef struct qaws_decimate_report
{
	unsigned int faces_before, faces_after;
	unsigned int collapses;
	qaws_scalar max_error;              /* the largest error of a collapse made */
	unsigned int creases;               /* edges made features by crease_angle */
	unsigned int silhouette_edges;      /* silhouette edges, summed over the views */
	qaws_mesh_deviation deviation;      /* with QAWS_DECIMATE_MEASURE */
	unsigned int seams;                 /* edges made features as attribute seams */
} qaws_decimate_report;

/* report may be NULL; desc NULL: half the faces */
qaws_status qaws_mesh_3d_decimate(qaws_decimate_report* report, qaws_mesh_3d* mesh, qaws_decimate_desc const* desc);

#endif /* QAWS_MESH_DECIMATE_H */
