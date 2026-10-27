#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_VIEW_H_F32
#define QAWS_MESH_VIEW_H_F32
#define QAWS_MESH_VIEW_H_BODY
#endif
#else
#ifndef QAWS_MESH_VIEW_H_F64
#define QAWS_MESH_VIEW_H_F64
#define QAWS_MESH_VIEW_H_BODY
#endif
#endif
#ifdef QAWS_MESH_VIEW_H_BODY
#undef QAWS_MESH_VIEW_H_BODY
#define QAWS_MESH_VIEW_H

#include "qaws_view.h"
#include "qaws_mesh.h"
#include "qaws_mesh_bvh.h"
#include "qaws_motion.h"

/*
 * Meshes from a host's buffers: a vertex buffer (f32 or f64, interleaved
 * with other attributes at any stride, or planar) and an index buffer (16-
 * or 32-bit), without a staging copy on the host's side.
 *
 * Derived data is built once and owned (a mesh's connectivity, a BVH's
 * tree); positions are read from the view at each use, converted as they
 * are read (no allocation to update them). The lifetime contract is the
 * one of qaws_view.h: qaws never writes, frees or keeps host memory, except
 * a borrowed vertex motion, which reads its two views at every evaluation
 * for its whole life (they must stay valid, and must not change while an
 * evaluation runs).
 *
 *   qaws_mesh_3d_create_from_view          an owned mesh (one pass)
 *   qaws_mesh_3d_set_positions_from_view   a frame's positions into it
 *   qaws_mesh_bvh_create_view              a BVH over the triangles
 *   qaws_mesh_bvh_refit_view               a frame's positions into it, the
 *                                          tree kept (a skinned or simulated
 *                                          mesh)
 *   qaws_motion_create_vertices_borrowed   each vertex on its chord from one
 *                                          positions view to another (last
 *                                          frame's and this one's): contact
 *                                          detection gives it exact times
 */

#ifndef QAWS_ENUM_QAWS_INDEX_TYPE
#define QAWS_ENUM_QAWS_INDEX_TYPE
typedef enum qaws_index_type
{
	QAWS_INDEX_U16 = 1,
	QAWS_INDEX_U32 = 2
} qaws_index_type;
#endif

typedef struct qaws_mesh_view
{
	qaws_array_view positions;          /* 3 components; count: the vertex count */
	void const* indices;                /* 3 per triangle, read bytewise (no alignment needed) */
	qaws_index_type index_type;
	size_t index_stride;                /* bytes from a triangle to the next; 0: tight (3 indices) */
	unsigned int triangle_count;
} qaws_mesh_view;

/* INVALID_ARGUMENT unless the positions are a valid view of 3 components, every position finite, a known index type, a stride
   of at least 3 indices (or 0), every index below the vertex count */
qaws_status qaws_mesh_view_validate(qaws_mesh_view const* view);

/* An owned mesh of the view's triangles. INVALID_ARGUMENT as qaws_mesh_view_validate, or for a vertex twice in a triangle. */
qaws_status qaws_mesh_3d_create_from_view(qaws_mesh_3d** out_mesh, qaws_mesh_view const* view);

/* The mesh's positions from a view of as many vertices (3 components, finite); no allocation. INVALID_ARGUMENT otherwise (the
   mesh unchanged). */
qaws_status qaws_mesh_3d_set_positions_from_view(qaws_mesh_3d* mesh, qaws_array_view const* positions);

/* A BVH over the view's triangles (elements: triangle indices). desc may be NULL; its leaf_size and winding_accuracy are read. */
qaws_status qaws_mesh_bvh_create_view(qaws_mesh_bvh** out_bvh, qaws_mesh_view const* view, qaws_mesh_bvh_desc const* desc);

/* positions: as many vertices as the BVH has, of its dimension (2 or 3 components); no allocation. INVALID_ARGUMENT for another
   count or dimension, a value not finite (the BVH unchanged). */
qaws_status qaws_mesh_bvh_refit_view(qaws_mesh_bvh* bvh, qaws_array_view const* positions);

/* A vertex motion moving each vertex on its chord at a constant rate from positions0 (t = 0) to positions1 (t = 1), both views
   of the same count, 3 components, borrowed for the motion's life (read at each evaluation). INVALID_ARGUMENT for invalid views,
   counts that differ, a value not finite. */
qaws_status qaws_motion_create_vertices_borrowed(qaws_motion** out_motion, qaws_array_view const* positions0, qaws_array_view const* positions1);

#endif /* QAWS_MESH_VIEW_H */
