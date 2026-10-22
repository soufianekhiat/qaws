#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_OFFSET_H_F32
#define QAWS_MESH_OFFSET_H_F32
#define QAWS_MESH_OFFSET_H_BODY
#endif
#else
#ifndef QAWS_MESH_OFFSET_H_F64
#define QAWS_MESH_OFFSET_H_F64
#define QAWS_MESH_OFFSET_H_BODY
#endif
#endif
#ifdef QAWS_MESH_OFFSET_H_BODY
#undef QAWS_MESH_OFFSET_H_BODY
#define QAWS_MESH_OFFSET_H

#include "qaws_mesh.h"
#include "qaws_mesh_intersect.h"

/*
 * Offsets of meshes: the surface at a signed distance, shells and
 * thickened sheets.
 *
 * Kinds (d the desc's distance):
 *   SURFACE  the surface at signed distance d from a closed mesh (outward
 *            for d > 0, inward for d < 0; negative inside, as
 *            qaws_mesh_bvh_signed_distance);
 *   SHELL    a closed mesh hollowed: the mesh and its inward offset by |d|
 *            turned inward, a wall |d| thick (two closed surfaces, one
 *            inside the other);
 *   THICKEN  any surface, open sheets too, made a solid |d| thick: the
 *            points within |d| / 2 of it.
 *
 * Methods:
 *   CONTOUR  the offset as a level set of the exact distance to the input
 *            (the BVH's, signed by winding numbers), contoured on a grid by
 *            dual contouring: crossings found along the grid's edges by
 *            regula falsi on the distance, their normals the distance's
 *            gradient (from the closest point), one vertex per cell at
 *            the least squares point of its crossings' planes. The true
 *            offset: rounded where it turns about the input's convex edges
 *            and corners (and a sheet's border), sharp where it meets
 *            itself (concave edges, inward offsets of corners). Detail below
 *            the cell is lost; resolution sets the cells along the longest
 *            side of the result's box. The result goes through
 *            qaws_mesh_repair to split non-manifold places.
 *   MOVE     the input's vertices moved: each along the direction m that
 *            moves its faces' planes by one (the least squares solution of
 *            n_f . m = 1 over its faces, weighted by their angles there),
 *            times d: flat faces offset exactly, edges and corners mitred
 *            (kept sharp: a cube stays a cube). The faces and their
 *            connectivity are kept; THICKEN joins the two sides along the
 *            border by quads. Where the moved surface overlaps or folds
 *            (inward offsets thinner than |d|, outward ones closing a gap)
 *            it is resolved by the self-union (qaws_mesh_boolean): the
 *            region its winding number covers by at least 1/2, inverted
 *            folds gone. The mitre is not the true offset at convex
 *            corners (farther from the input, up to |m| |d|).
 *   ROUND    the true offset of the triangles as they are: the points
 *            within |d| of the surface, the union of closed convex pieces
 *            (each triangle's slab |d| to either side, a cylinder round each
 *            edge whose faces turn or which is on the border, a ball at each
 *            of their vertices), united with the input outward, taken from it
 *            inward, by the mesh Booleans. Cylinders where the faces turn
 *            outward, spheres at convex corners, sharp where the offset meets
 *            itself; the circles split every round_angle at most, their
 *            points at |d| exactly and the chords between them inside by at
 *            most |d| (1 - cos(round_angle / 2)). SURFACE and SHELL want a
 *            closed mesh (the inside is read by winding numbers); THICKEN
 *            takes sheets, their border rounded.
 *
 * distances (MOVE): a distance per vertex, each vertex moved along its
 * mitre by its own; distance must still be nonzero and finite (its sign
 * no longer matters). A varying offset is not the true offset of anything:
 * where it changes fast the moved surface folds, which the self-union
 * resolves. CONTOUR with distances: INVALID_ARGUMENT (shortcut: a contour of
 * a varying offset needs a level per point of the distance field).
 *
 * The report measures the offset surface (not the input's copy in a
 * shell): its vertices' distances to the input, smallest and largest
 * (SURFACE and SHELL: d's size; THICKEN: half of it, exactly where it is
 * the true offset).
 */

#ifndef QAWS_ENUM_QAWS_MESH_OFFSET_KIND
#define QAWS_ENUM_QAWS_MESH_OFFSET_KIND
typedef enum qaws_mesh_offset_kind
{
	QAWS_MESH_OFFSET_SURFACE = 0,
	QAWS_MESH_OFFSET_SHELL,
	QAWS_MESH_OFFSET_THICKEN
} qaws_mesh_offset_kind;
#endif

#ifndef QAWS_ENUM_QAWS_MESH_OFFSET_METHOD
#define QAWS_ENUM_QAWS_MESH_OFFSET_METHOD
typedef enum qaws_mesh_offset_method
{
	QAWS_MESH_OFFSET_CONTOUR = 0,
	QAWS_MESH_OFFSET_MOVE,
	QAWS_MESH_OFFSET_ROUND
} qaws_mesh_offset_method;
#endif

typedef struct qaws_mesh_offset_desc
{
	qaws_mesh_offset_kind kind;
	qaws_mesh_offset_method method;
	qaws_scalar distance;               /* d (SHELL and THICKEN: its size) */
	unsigned int resolution;            /* CONTOUR: cells along the result's box's longest side; 0: 64 */
	qaws_scalar const* distances;       /* MOVE: per input vertex its own distance (signed for SURFACE, its size for SHELL and THICKEN); NULL: distance everywhere */
	qaws_scalar round_angle;            /* ROUND: the arcs' largest step, radians, 0.01 to 1.6; 0: pi / 8 */
} qaws_mesh_offset_desc;

typedef struct qaws_mesh_offset_report
{
	unsigned int grid[3];               /* CONTOUR: grid points per axis */
	qaws_scalar cell;                   /* CONTOUR: the cell size */
	unsigned int crossings;             /* CONTOUR: grid edges the offset crosses */
	unsigned int sharp;                 /* CONTOUR: cell vertices on an edge or a corner (their planes span 2 or 3 directions) */
	qaws_intersect_report intersect;    /* MOVE: the moved surface's self-intersections, resolved */
	qaws_scalar min_distance;           /* the offset surface's vertices' distances to the input */
	qaws_scalar max_distance;
	int closed;                         /* the result's directed edges all cancel */
	unsigned int evaluated;             /* CONTOUR: grid points the distance was evaluated at (only near the level: a block of
	                                       4^3 cells whose corner is farther than its diagonal holds no crossing) */
} qaws_mesh_offset_report;

/* report may be NULL. INVALID_ARGUMENT for a distance zero or not finite, an unknown kind or method, a mesh without faces.
   An empty result (an inward offset past the input's thickness): OK and out_mesh NULL. */
qaws_status qaws_mesh_3d_offset(qaws_mesh_3d** out_mesh, qaws_mesh_offset_report* report, qaws_mesh_3d const* mesh, qaws_mesh_offset_desc const* desc);

#endif /* QAWS_MESH_OFFSET_H */
