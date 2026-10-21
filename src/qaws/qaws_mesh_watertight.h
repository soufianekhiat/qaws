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
 * The inside is where the input's generalized winding number is at least
 * 1/2 (qaws_mesh_bvh.h): holes are closed smoothly, overlaps united,
 * inner faces dropped. With an offset, points within it of the input are
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

typedef struct qaws_watertight_desc
{
	unsigned int resolution;           /* cells along the box's longest side; 0: 64 */
	qaws_scalar offset;                /* points within it of the input count as inside; 0: none */
} qaws_watertight_desc;

typedef struct qaws_watertight_report
{
	unsigned int grid[3];              /* grid points per axis */
	unsigned int evaluated;            /* winding numbers computed (near points, one per far region) */
	unsigned int crossings;            /* grid edges the boundary crosses */
	unsigned int sharp;                /* cell vertices whose crossings' planes span 2 or 3 directions (on an edge or a corner) */
	qaws_scalar cell;                  /* the cell size */
	int closed;                        /* the result's directed edges all cancel */
} qaws_watertight_report;

/* face_sizes NULL: triangles (polygons are fanned); desc and report may be NULL. INVALID_ARGUMENT for an index out of range,
   a face under 3 vertices, a position not finite, a negative offset, no faces. */
qaws_status qaws_mesh_make_watertight(qaws_mesh_3d** out_mesh, qaws_watertight_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_watertight_desc const* desc);

#endif /* QAWS_MESH_WATERTIGHT_H */
