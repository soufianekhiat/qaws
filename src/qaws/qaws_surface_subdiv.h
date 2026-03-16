#ifndef QAWS_SURFACE_SUBDIV_H
#define QAWS_SURFACE_SUBDIV_H

#include "qaws_surface_types.h"

/* Subdivision surface: limit surface from iterative mesh refinement.

   Supports two subdivision schemes:
   - Catmull-Clark: for quad-dominant meshes, produces C2 surfaces away from
     extraordinary vertices. Based on cubic B-spline refinement.
   - Loop: for triangle meshes, produces C2 surfaces away from extraordinary
     vertices. Based on quartic box spline refinement.

   The control mesh is defined by vertices and faces.
   Faces can be quads (Catmull-Clark) or triangles (Loop).
   Mixed meshes use Catmull-Clark (triangles are treated as degenerate quads).

   Evaluation maps (u,v) in [0,1]x[0,1] to the limit surface.
   Internally, the mesh is subdivided to the specified level before
   sampling, and the result is stored on a regular grid for fast evaluation
   via bilinear interpolation.

   The control mesh data is copied (not borrowed). */

typedef enum qaws_subdiv_scheme
{
	QAWS_SUBDIV_CATMULL_CLARK = 0,
	QAWS_SUBDIV_LOOP
} qaws_subdiv_scheme;

typedef struct qaws_surface_subdiv_desc
{
	qaws_vec3 const* vertices;          /* vertex positions */
	unsigned int vertex_count;
	unsigned int const* face_indices;   /* flat array of vertex indices per face */
	unsigned int const* face_sizes;     /* number of vertices per face (3 or 4) */
	unsigned int face_count;
	qaws_subdiv_scheme scheme;          /* subdivision scheme */
	unsigned int subdivision_level;     /* number of subdivision iterations (0 = default 3) */
} qaws_surface_subdiv_desc;

qaws_status qaws_surface_create_subdiv(
	qaws_surface_subdiv_desc const* desc,
	qaws_surface** out_surface);

/* Retrieve the subdivided triangle mesh for direct rendering.
   Returns pointers into the surface's internal data (valid until destroy).
   Triangles are CCW-wound index triples into the vertex array.
   Normals are per-vertex smooth normals. */
qaws_status qaws_surface_subdiv_get_mesh(
	qaws_surface const* surface,
	qaws_vec3 const** out_vertices,
	qaws_vec3 const** out_normals,
	unsigned int const** out_triangles,
	unsigned int* out_vertex_count,
	unsigned int* out_triangle_count);

#endif /* QAWS_SURFACE_SUBDIV_H */
