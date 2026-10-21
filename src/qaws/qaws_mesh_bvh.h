#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_BVH_H_F32
#define QAWS_MESH_BVH_H_F32
#define QAWS_MESH_BVH_H_BODY
#endif
#else
#ifndef QAWS_MESH_BVH_H_F64
#define QAWS_MESH_BVH_H_F64
#define QAWS_MESH_BVH_H_BODY
#endif
#endif
#ifdef QAWS_MESH_BVH_H_BODY
#undef QAWS_MESH_BVH_H_BODY
#define QAWS_MESH_BVH_H

#include "qaws_mesh.h"
#include "qaws_triangulate.h"
#include "qaws_volume_mesh.h"
#include "qaws_batch_executor.h"

/*
 * Spatial queries on meshes: a bounding volume hierarchy over simplices.
 *
 * A BVH holds simplices of one kind in 2D or 3D: segments (curves, a 2D
 * mesh's border), triangles (surfaces, 2D meshes' faces) or tetrahedra
 * (volumes). Each simplex belongs to an element of the input: a face fanned
 * from its first vertex, a cell split into tetrahedra from its first vertex
 * (a hexahedron into 6, a prism 3, a pyramid 2), or itself for arrays.
 * Soups are welcome: nothing is assumed of the topology.
 *
 * Build: binned SAH (8 bins per axis, as tinyBVH's default), leaves of at
 * most `leaf_size` simplices, collapsed into wide nodes (the largest child
 * opened first) whose float boxes, rounded outward, are tested together:
 * 8 per node with AVX2, 4 with SSE2, the path chosen once at creation from
 * what the processor has. Positions are held in double whatever the precision of the caller,
 * and every answer comes from the double tests: rays first screen a leaf's
 * triangles in float with a margin that only lets more through.
 *
 * Queries, each over a batch of points or rays, in chunks through the
 * caller's executor (NULL: the calling thread; the results do not depend on
 * it):
 *   closest      the nearest point on the simplices within max_distance;
 *                inside a triangle of a 2D BVH or a tetrahedron, the point
 *                itself at distance 0: point location, with the weights of
 *                the simplex's vertices for interpolation;
 *   ray          the first hit (or any hit) along o + t d, t in [0, tmax]:
 *                triangles in 3D by the watertight test of Woop, Benthin and
 *                Wald (2013), so a ray through a shared edge or vertex never
 *                slips between the triangles; segments in 2D;
 *   winding      the generalized winding number (Jacobson et al. 2013) of
 *                triangles in 3D (solid angles over 4 pi) or segments in 2D
 *                (angles over 2 pi): 1 inside a closed outward surface or
 *                counter-clockwise loop, 0 outside, a smooth value between
 *                for open, broken or overlapping input. Far nodes by their
 *                dipole expansion (Barill et al. 2018), exact sums near:
 *                a node is far when the point is past `winding_accuracy`
 *                times its radius from its centre;
 *   signed       the distance to the closest point, negative where the
 *                winding number is at least 1/2: robust on soups and holes.
 *
 * refit moves the vertices (same count, same simplices): the boxes and the
 * winding moments are recomputed, the tree kept.
 *
 * Thread-safe for concurrent queries; refit needs exclusive access.
 */

typedef struct qaws_mesh_bvh qaws_mesh_bvh;

typedef struct qaws_mesh_bvh_desc
{
	qaws_scalar const* positions;      /* dimension per vertex */
	unsigned int vertex_count;
	unsigned int dimension;            /* 2 or 3 */
	unsigned int const* simplices;     /* simplex_size per simplex */
	unsigned int simplex_size;         /* 2 segments, 3 triangles, 4 tetrahedra */
	unsigned int simplex_count;
	unsigned int const* elements;      /* per simplex its element; NULL: the simplex's own index */
	unsigned int leaf_size;            /* 0: 4 */
	qaws_scalar winding_accuracy;      /* 0: 2.25 (Barill et al.'s beta) */
} qaws_mesh_bvh_desc;

/* A query's answer. element QAWS_MESH_NONE: nothing within reach. */
typedef struct qaws_mesh_hit
{
	unsigned int element;              /* face, cell, segment or simplex (desc.elements) */
	unsigned int simplex;
	unsigned int vertices[4];          /* the simplex's vertices (QAWS_MESH_NONE past its size) */
	qaws_scalar weights[4];            /* barycentric weights of the point on them */
	qaws_scalar point[3];              /* z 0 in 2D */
	qaws_scalar distance;              /* closest: the distance; ray: t */
} qaws_mesh_hit;

qaws_status qaws_mesh_bvh_create(qaws_mesh_bvh** out_bvh, qaws_mesh_bvh_desc const* desc);
/* The faces (elements: face indices). desc may be NULL; its leaf_size and winding_accuracy are read. */
qaws_status qaws_mesh_bvh_create_mesh_3d(qaws_mesh_bvh** out_bvh, qaws_mesh_3d const* mesh, qaws_mesh_bvh_desc const* desc);
/* The faces (triangles, then quads), or with border the boundary edges (segments, region on their left). */
qaws_status qaws_mesh_bvh_create_mesh_2d(qaws_mesh_bvh** out_bvh, qaws_mesh_2d const* mesh, int border, qaws_mesh_bvh_desc const* desc);
/* The cells, as tetrahedra. */
qaws_status qaws_mesh_bvh_create_volume(qaws_mesh_bvh** out_bvh, qaws_volume_mesh const* mesh, qaws_mesh_bvh_desc const* desc);
void qaws_mesh_bvh_destroy(qaws_mesh_bvh* bvh);

/* positions: dimension per vertex, as many vertices as created with */
qaws_status qaws_mesh_bvh_refit(qaws_mesh_bvh* bvh, qaws_scalar const* positions);

/* points: dimension per point. max_distance 0: no limit. */
qaws_status qaws_mesh_bvh_closest(qaws_mesh_hit* out_hits, qaws_mesh_bvh const* bvh, qaws_scalar const* points, unsigned int count,
	qaws_scalar max_distance, qaws_batch_executor const* executor);

/* origins, directions: dimension per ray (directions need not be unit; t
   counts in their lengths). tmax 0: no limit. any_hit: stop at the first
   hit found (shadow rays). UNSUPPORTED_OPERATION but for triangles in 3D and
   segments in 2D. */
qaws_status qaws_mesh_bvh_ray(qaws_mesh_hit* out_hits, qaws_mesh_bvh const* bvh, qaws_scalar const* origins, qaws_scalar const* directions,
	unsigned int count, qaws_scalar tmax, int any_hit, qaws_batch_executor const* executor);

/* UNSUPPORTED_OPERATION but for triangles in 3D and segments in 2D */
qaws_status qaws_mesh_bvh_winding_number(qaws_scalar* out_winding, qaws_mesh_bvh const* bvh, qaws_scalar const* points, unsigned int count,
	qaws_batch_executor const* executor);

/* out_hits may be NULL. UNSUPPORTED_OPERATION as the winding number. */
qaws_status qaws_mesh_bvh_signed_distance(qaws_scalar* out_distances, qaws_mesh_hit* out_hits, qaws_mesh_bvh const* bvh, qaws_scalar const* points,
	unsigned int count, qaws_batch_executor const* executor);

#endif /* QAWS_MESH_BVH_H */
