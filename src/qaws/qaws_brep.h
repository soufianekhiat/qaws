#ifndef QAWS_BREP_H
#define QAWS_BREP_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"
#include "qaws_curve.h"

/* B-rep (Boundary Representation) shell: topological structure
   representing a solid or sheet body as a collection of faces.

   Each face is a (possibly trimmed) surface with oriented boundary edges.
   Edges are shared between exactly two faces (manifold) or are boundary edges.
   Vertices are shared between edges.

   The topology uses a half-edge data structure:
   - Each edge has two half-edges, one for each adjacent face.
   - Each half-edge knows its face, next/prev half-edges in the face loop,
     and its twin half-edge on the adjacent face.

   Surfaces and curves are borrowed (not owned by the shell). */

typedef struct qaws_brep_shell qaws_brep_shell;

/* Opaque handles for topological entities */
typedef unsigned int qaws_brep_face_id;
typedef unsigned int qaws_brep_edge_id;
typedef unsigned int qaws_brep_vertex_id;
typedef unsigned int qaws_brep_halfedge_id;

#define QAWS_BREP_INVALID_ID ((unsigned int)0xFFFFFFFF)

/* Face information */
typedef struct qaws_brep_face_desc
{
	qaws_surface const* surface;    /* the face geometry (borrowed) */
	int is_reversed;                /* 1 = normal points inward */
} qaws_brep_face_desc;

/* Edge information */
typedef struct qaws_brep_edge_desc
{
	qaws_brep_vertex_id v_start;
	qaws_brep_vertex_id v_end;
	qaws_curve const* curve_3d;     /* 3D edge curve (borrowed, may be NULL) */
} qaws_brep_edge_desc;

/* Shell creation */
qaws_status qaws_brep_create(qaws_brep_shell** out_shell);
void qaws_brep_destroy(qaws_brep_shell* shell);

/* Add topological entities */
qaws_status qaws_brep_add_vertex(
	qaws_brep_shell* shell,
	qaws_vec3 position,
	qaws_brep_vertex_id* out_id);

qaws_status qaws_brep_add_edge(
	qaws_brep_shell* shell,
	qaws_brep_edge_desc const* desc,
	qaws_brep_edge_id* out_id);

qaws_status qaws_brep_add_face(
	qaws_brep_shell* shell,
	qaws_brep_face_desc const* desc,
	qaws_brep_edge_id const* edge_ids,
	int const* edge_orientations,   /* +1 or -1 per edge */
	unsigned int edge_count,
	qaws_brep_face_id* out_id);

/* Query topology */
unsigned int qaws_brep_get_vertex_count(qaws_brep_shell const* shell);
unsigned int qaws_brep_get_edge_count(qaws_brep_shell const* shell);
unsigned int qaws_brep_get_face_count(qaws_brep_shell const* shell);

qaws_status qaws_brep_get_vertex_position(
	qaws_brep_shell const* shell,
	qaws_brep_vertex_id id,
	qaws_vec3* out_position);

qaws_status qaws_brep_get_face_surface(
	qaws_brep_shell const* shell,
	qaws_brep_face_id id,
	qaws_surface const** out_surface);

qaws_status qaws_brep_get_edge_vertices(
	qaws_brep_shell const* shell,
	qaws_brep_edge_id id,
	qaws_brep_vertex_id* out_v_start,
	qaws_brep_vertex_id* out_v_end);

qaws_status qaws_brep_get_edge_faces(
	qaws_brep_shell const* shell,
	qaws_brep_edge_id id,
	qaws_brep_face_id* out_face_left,
	qaws_brep_face_id* out_face_right);

qaws_status qaws_brep_get_face_edges(
	qaws_brep_shell const* shell,
	qaws_brep_face_id id,
	qaws_brep_edge_id* out_edges,
	unsigned int capacity,
	unsigned int* out_count);

/* Adjacency queries */
qaws_status qaws_brep_get_vertex_edges(
	qaws_brep_shell const* shell,
	qaws_brep_vertex_id id,
	qaws_brep_edge_id* out_edges,
	unsigned int capacity,
	unsigned int* out_count);

qaws_status qaws_brep_get_vertex_faces(
	qaws_brep_shell const* shell,
	qaws_brep_vertex_id id,
	qaws_brep_face_id* out_faces,
	unsigned int capacity,
	unsigned int* out_count);

/* Validation */
typedef enum qaws_brep_validation_flags
{
	QAWS_BREP_VALID            = 0,
	QAWS_BREP_NOT_MANIFOLD     = 1 << 0,  /* edge shared by more than 2 faces */
	QAWS_BREP_NOT_CLOSED       = 1 << 1,  /* boundary edges exist */
	QAWS_BREP_NOT_ORIENTABLE   = 1 << 2,  /* inconsistent face orientations */
	QAWS_BREP_HAS_DEGENERATE   = 1 << 3   /* zero-area face or zero-length edge */
} qaws_brep_validation_flags;

qaws_status qaws_brep_validate(
	qaws_brep_shell const* shell,
	unsigned int* out_flags);

/* Euler characteristic: V - E + F = 2 for closed genus-0 shell */
qaws_status qaws_brep_euler_characteristic(
	qaws_brep_shell const* shell,
	int* out_chi);

/* Compute shell volume (only valid for closed shells).
   Uses the divergence theorem on the face meshes. */
qaws_status qaws_brep_compute_volume(
	qaws_brep_shell const* shell,
	qaws_scalar* out_volume);

/* Compute total surface area (sum of all face areas). */
qaws_status qaws_brep_compute_surface_area(
	qaws_brep_shell const* shell,
	qaws_scalar* out_area);

#endif /* QAWS_BREP_H */
