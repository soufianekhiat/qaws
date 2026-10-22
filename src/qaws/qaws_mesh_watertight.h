#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_WATERTIGHT_H_F32
#define QAWS_MESH_WATERTIGHT_H_F32
#define QAWS_MESH_WATERTIGHT_H_BODY
#endif
#else
#ifndef QAWS_MESH_WATERTIGHT_H_F64
#define QAWS_MESH_WATERTIGHT_H_F64
#define QAWS_MESH_WATERTIGHT_H_BODY
#endif
#endif
#ifdef QAWS_MESH_WATERTIGHT_H_BODY
#undef QAWS_MESH_WATERTIGHT_H_BODY
#define QAWS_MESH_WATERTIGHT_H

#include "qaws_mesh.h"

/*
 * Watertight surfaces from any triangle input: soups, scans with holes,
 * overlapping or self-intersecting parts, open sheets. The result is a
 * new closed, manifold, outward triangle surface that follows the input,
 * not a repair of it (qaws_mesh_repair keeps the input's faces).
 *
 * Holes are first filled by qaws_mesh_repair's fills (Liepa's triangulation,
 * refined, faired by the bi-Laplacian with the border's ring held: tangent
 * to the surface round them, continuing its curvature), unless
 * QAWS_WATERTIGHT_WINDING_CAPS asks for the winding number's own caps (its
 * 1/2 level across a hole: nearer a minimal surface, a ridge at the rim) or
 * an offset is given (a sheet's border is no hole there).
 *
 * The inside is where the generalized winding number is at least 1/2
 * (qaws_mesh_bvh.h): holes left are closed, overlaps united, inner faces
 * dropped. With an offset, points within it of the input are
 * inside too: open sheets and thin parts become closed shells of that
 * thickness (an envelope of the input).
 *
 * The surface is the inside's boundary, contoured on a grid by dual
 * contouring (Ju et al. 2002): one vertex per cell the boundary crosses,
 * placed at the least squares point of the crossings' tangent planes, so
 * sharp edges and corners of the input are kept. Crossings on the input
 * are found by rays along the grid's edges (the input's own triangles
 * give their points and normals); on an offset or across a hole, by
 * bisection of the inside test. Only grid points near the input are
 * evaluated; the others take their side from one winding number per
 * region they flood. The surface then goes through qaws_mesh_repair to
 * split the few non-manifold places a cell's two sheets make
 * (QAWS_MESH_REPAIR_SPLIT_NON_MANIFOLD) and drop slivers.
 *
 * Detail below the cell size is lost: resolution sets the cells along
 * the longest side of the input's box.
 */

#define QAWS_WATERTIGHT_WINDING_CAPS 1u     /* holes capped by the winding number's 1/2 level, not filled first */

typedef struct qaws_watertight_desc
{
	unsigned int resolution;           /* cells along the box's longest side; 0: 64 */
	qaws_scalar offset;                /* points within it of the input count as inside; 0: none */
	unsigned int flags;                /* QAWS_WATERTIGHT_* */
} qaws_watertight_desc;

typedef struct qaws_watertight_report
{
	unsigned int grid[3];              /* grid points per axis */
	unsigned int evaluated;            /* winding numbers computed (near points, one per far region) */
	unsigned int crossings;            /* grid edges the boundary crosses */
	unsigned int sharp;                /* cell vertices whose crossings' planes span 2 or 3 directions (on an edge or a corner) */
	qaws_scalar cell;                  /* the cell size */
	unsigned int holes_filled;         /* by the fills before contouring */
	int closed;                        /* the result's directed edges all cancel */
} qaws_watertight_report;

/* face_sizes NULL: triangles (polygons are fanned); desc and report may be NULL. INVALID_ARGUMENT for an index out of range,
   a face under 3 vertices, a position not finite, a negative offset, no faces. */
qaws_status qaws_mesh_make_watertight(qaws_mesh_3d** out_mesh, qaws_watertight_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_watertight_desc const* desc);

/* A watertight surface from an oriented point cloud (positions and normals,
   3 each; normals outward, any length): the inside is where the points'
   winding number is at least 1/2, each point a dipole of its share of the
   surface (Barill et al. 2018: pi r^2 / K, r the distance to its K = 10th
   neighbour), summed by the BVH's fast winding numbers; contoured on the
   grid by dual contouring, crossings by regula falsi on the winding number,
   their normals the nearest point's (the winding number's gradient there is
   all near field; the nearest, not a mean, keeps sampled edges sharp). desc's offset
   must be 0; its flags are not
   read. INVALID_ARGUMENT for fewer than 4 points, a value not finite;
   INVALID_TOPOLOGY for normals turned inward (sum of area n . (p - centre),
   three times the volume, not positive) or when nothing comes out inside.
   The report's holes_filled is 0, evaluated the grid's points. */
qaws_status qaws_mesh_make_watertight_points(qaws_mesh_3d** out_mesh, qaws_watertight_report* report, qaws_scalar const* positions,
	qaws_scalar const* normals, unsigned int count, qaws_watertight_desc const* desc);

#endif /* QAWS_MESH_WATERTIGHT_H */
