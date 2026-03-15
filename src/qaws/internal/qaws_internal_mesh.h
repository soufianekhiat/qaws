#ifndef QAWS_INTERNAL_MESH_H
#define QAWS_INTERNAL_MESH_H

#include "qaws_internal_types.h"
#include <limits.h>

/*
 * Edge-adjacency triangle mesh built from tessellation output.
 * Transient: allocated, used, freed within a single API call.
 */

typedef struct qaws_internal_edge {
	unsigned int v[2];        /* vertex indices, v[0] < v[1] canonical */
	unsigned int tri[2];      /* adjacent triangles; tri[1]=UINT_MAX if boundary */
	unsigned int opposite[2]; /* vertex in tri[k] not on this edge */
} qaws_internal_edge;

typedef struct qaws_internal_tri {
	unsigned int v[3];    /* CCW vertex indices */
	unsigned int edge[3]; /* edge[i] opposite vertex v[i] */
} qaws_internal_tri;

typedef struct qaws_internal_mesh {
	qaws_scalar const* positions;  /* 3 * vertex_count floats (xyz) */
	qaws_scalar const* uvs;        /* 2 * vertex_count floats (uv) */
	unsigned int vertex_count;
	qaws_internal_tri* tris;
	unsigned int tri_count;
	qaws_internal_edge* edges;
	unsigned int edge_count;
	unsigned int* vert_edge_offset; /* CSR-like, vertex_count+1 */
	unsigned int* vert_edges;       /* edge indices per vertex */
} qaws_internal_mesh;

/* Build mesh from tessellation vertex/index buffers.
   positions: [vertex_count * 3], uvs: [vertex_count * 2],
   indices: [index_count] (triangles, index_count must be multiple of 3).
   Returns QAWS_STATUS_OK on success. Caller must call qaws_internal_mesh_destroy. */
qaws_status qaws_internal_mesh_build(
	qaws_internal_mesh* mesh,
	qaws_scalar const* positions,
	qaws_scalar const* uvs,
	unsigned int vertex_count,
	unsigned int const* indices,
	unsigned int index_count);

void qaws_internal_mesh_destroy(qaws_internal_mesh* mesh);

/* Dijkstra shortest path on mesh edges. Returns vertex path in out_path (reverse order).
   out_path must have capacity >= vertex_count. Returns path length in out_path_count.
   Returns QAWS_STATUS_OK on success, QAWS_STATUS_NUMERICAL_FAILURE if disconnected. */
qaws_status qaws_internal_mesh_dijkstra(
	qaws_internal_mesh const* mesh,
	unsigned int src,
	unsigned int dst,
	unsigned int* out_path,
	unsigned int* out_path_count);

/* Compute 3D edge length between two vertices */
qaws_scalar qaws_internal_mesh_edge_length(
	qaws_internal_mesh const* mesh,
	unsigned int v0,
	unsigned int v1);

/* Compute cotangent weight for an edge. Clamped to [-1e6, 1e6]. */
qaws_scalar qaws_internal_mesh_cotan_weight(
	qaws_internal_mesh const* mesh,
	unsigned int edge_idx);

/* Try to flip an edge. Returns 1 if flipped, 0 if not (boundary or non-convex). */
int qaws_internal_mesh_try_flip_edge(
	qaws_internal_mesh* mesh,
	unsigned int edge_idx);

/* Find the closest vertex to a given (u,v) parameter. */
unsigned int qaws_internal_mesh_closest_vertex_uv(
	qaws_internal_mesh const* mesh,
	qaws_scalar u,
	qaws_scalar v);

#endif /* QAWS_INTERNAL_MESH_H */
