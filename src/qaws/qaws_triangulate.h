#ifndef QAWS_TRIANGULATE_H
#define QAWS_TRIANGULATE_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_mesh.h"

/*
 * Triangulating regions of curves.
 *
 * The region (closed paths read with the fill rule) is first cleaned by one
 * union (qaws_clip): nested loops of curves, no overlaps. Their curves are
 * sampled to the tolerance (chord deviation), every boundary vertex on its
 * curve, with the curve and its parameter kept. A constrained Delaunay
 * triangulation of the samples follows: incremental (Bowyer-Watson), the
 * boundary edges recovered by flips, the triangles inside kept (parity
 * across the boundary).
 *
 * Coordinates are snapped onto a lattice of 2^-48 of the region's extent
 * (a move far below any tolerance) so that orientation and in-circle tests
 * are exact: a double filter, then 384-bit integers.
 *
 * Refinement (Ruppert): a triangle with an angle below min_angle, an area
 * over max_area or an edge longer than the size field gets its circumcentre;
 * a circumcentre that encroaches on a boundary edge (lies in its diametral
 * circle) splits that edge instead, at the parameter midpoint ON the curve,
 * so the boundary converges to the true curves. min_angle up to 33 degrees
 * terminates for inputs without sharp corners; corners sharper than 60
 * degrees keep their small angles (their triangles are not refined below
 * the tolerance) and max_vertices bounds the work in any case.
 *
 * Constraints: curves the mesh must follow inside the region (clipped to
 * it; crossing ones meet at a vertex on both curves), kept as edges like
 * the boundary but without changing the side.
 *
 * Remeshing (Botsch-Kobbelt): edges split above 4/3 of the target length,
 * collapsed below 4/5, flipped to Delaunay, then smoothed; boundary and
 * constraint vertices stay on their curves (splits at parameter midpoints,
 * smoothing along the parameter), corners stay.
 *
 * Quadrangulation: triangle pairs into convex quads, then optionally all
 * quads (see qaws_mesh_2d_quadrangulate).
 *
 * Thread-safe on immutable curves; the mesh owns copies of the boundary
 * curves.
 */

typedef struct qaws_mesh_2d qaws_mesh_2d;

/* target edge length at a point (refinement and remeshing) */
typedef qaws_scalar (*qaws_size_fn)(void* user, qaws_vec2 point);

typedef struct qaws_triangulate_desc
{
	qaws_path_2d const* paths;
	unsigned int path_count;
	qaws_fill_rule fill_rule;
	qaws_path_2d const* constraints;  /* optional: curves the mesh must follow inside the region (clipped to it; they may meet, not cross) */
	unsigned int constraint_count;
	qaws_scalar tolerance;          /* boundary chord deviation; 0 = 1e-3 of the extent */
	qaws_scalar min_angle;          /* degrees, 0 = no refinement by angle (at most 33) */
	qaws_scalar max_area;           /* 0 = none */
	qaws_size_fn size_fn;           /* optional: longest edge allowed at a point */
	void* size_user;
	unsigned int max_vertices;      /* 0 = 4 million */
} qaws_triangulate_desc;

qaws_status qaws_triangulate_region_2d(qaws_mesh_2d** out_mesh, qaws_triangulate_desc const* desc);

/* Remeshing in place toward edges of target_length (or the size field when
   size_fn is given), iterations rounds of split / collapse / flip /
   smooth. */
qaws_status qaws_mesh_2d_remesh(qaws_mesh_2d* mesh, qaws_scalar target_length, qaws_size_fn size_fn, void* size_user,
	unsigned int iterations);

void qaws_mesh_2d_destroy(qaws_mesh_2d* mesh);

unsigned int qaws_mesh_2d_get_vertex_count(qaws_mesh_2d const* mesh);
qaws_vec2 qaws_mesh_2d_get_vertex(qaws_mesh_2d const* mesh, unsigned int vertex);
unsigned int qaws_mesh_2d_get_triangle_count(qaws_mesh_2d const* mesh);
/* counter-clockwise vertex indices */
void qaws_mesh_2d_get_triangle(unsigned int out_vertices[3], qaws_mesh_2d const* mesh, unsigned int triangle);

/* A vertex on the boundary or a constraint: its curve (index among
   qaws_mesh_2d_get_boundary_curve, which lists the boundary curves then the
   constraint pieces) and parameter there; returns 0 for a free vertex. A
   vertex where two boundary curves meet reports the curve it starts. */
int qaws_mesh_2d_get_vertex_source(unsigned int* out_curve, qaws_scalar* out_parameter, qaws_mesh_2d const* mesh, unsigned int vertex);
unsigned int qaws_mesh_2d_get_boundary_curve_count(qaws_mesh_2d const* mesh);
qaws_curve const* qaws_mesh_2d_get_boundary_curve(qaws_mesh_2d const* mesh, unsigned int curve);

/* The edges on the boundary (each once, oriented with the region on the
   left) */
unsigned int qaws_mesh_2d_get_boundary_edge_count(qaws_mesh_2d const* mesh);
void qaws_mesh_2d_get_boundary_edge(unsigned int out_vertices[2], qaws_mesh_2d const* mesh, unsigned int edge);

/* The edges along the constraints (each once) */
unsigned int qaws_mesh_2d_get_constraint_edge_count(qaws_mesh_2d const* mesh);
void qaws_mesh_2d_get_constraint_edge(unsigned int out_vertices[2], qaws_mesh_2d const* mesh, unsigned int edge);

/* Quadrangulation in place: triangles paired into convex quadrilaterals
   (the best shaped pairs first, never across a boundary or constraint
   edge); with all_quads, every quadrilateral then split into 4 and every
   triangle left into 3 through edge midpoints (boundary and constraint
   midpoints on their curves) and centroids, so the mesh is quadrilaterals
   only; then smoothing_iterations rounds of smoothing that keep every
   quadrilateral convex. The triangles left (none with all_quads) stay
   triangles. A quadrangulated mesh is not remeshed. */
qaws_status qaws_mesh_2d_quadrangulate(qaws_mesh_2d* mesh, int all_quads, unsigned int smoothing_iterations);
unsigned int qaws_mesh_2d_get_quad_count(qaws_mesh_2d const* mesh);
/* counter-clockwise vertex indices */
void qaws_mesh_2d_get_quad(unsigned int out_vertices[4], qaws_mesh_2d const* mesh, unsigned int quad);

/*
 * The mesh as a corner store (the meshing module's, qaws_mesh.h).
 *
 * Every planar mesh is also held in the corner store of the surface meshes,
 * at z = 0: qaws_mesh_2d_get_mesh gives it as a qaws_mesh_3d to read with
 * qaws_mesh.h (faces of any size and their adjacency, the dual, topology,
 * quality, features: the boundary and constraint edges; attributes).
 * Its vertex and face indices are the ones above (triangles, then quads),
 * but where the region pinches: a vertex two boundary loops touch is split
 * there, its copies appended. The view is rebuilt by remeshing and
 * quadrangulation (attributes on it dropped).
 *
 * Meshes can be made from arrays and edited here, each edit checked as
 * qaws_mesh.h's are (INVALID_TOPOLOGY when refused) and the getters above
 * kept in step. Such a mesh has no curves under its vertices
 * (qaws_mesh_2d_get_vertex_source reports none): remeshing and
 * quadrangulation, which follow the curves, return UNSUPPORTED_OPERATION on
 * it.
 */
typedef struct qaws_mesh_2d_desc
{
	qaws_scalar const* positions;           /* 2 per vertex */
	unsigned int vertex_count;
	unsigned int const* face_vertices;      /* each face's vertices counter-clockwise, in turn */
	unsigned int const* face_sizes;         /* NULL: triangles */
	unsigned int face_count;
} qaws_mesh_2d_desc;

qaws_status qaws_mesh_2d_create(qaws_mesh_2d** out_mesh, qaws_mesh_2d_desc const* desc);
qaws_mesh_3d const* qaws_mesh_2d_get_mesh(qaws_mesh_2d const* mesh);

/* positions: 2 (NULL: the edge's midpoint, the far vertex, the face's centroid) */
qaws_status qaws_mesh_2d_split_edge(unsigned int* out_vertex, qaws_mesh_2d* mesh, qaws_scalar const* position, unsigned int a, unsigned int b);
qaws_status qaws_mesh_2d_collapse_edge(qaws_mesh_2d* mesh, unsigned int a, unsigned int b, qaws_scalar const* position);
qaws_status qaws_mesh_2d_flip_edge(qaws_mesh_2d* mesh, unsigned int a, unsigned int b);
qaws_status qaws_mesh_2d_split_face(unsigned int* out_vertex, qaws_mesh_2d* mesh, qaws_scalar const* position, unsigned int face);

qaws_status qaws_mesh_2d_add_attribute(unsigned int* out_attribute, qaws_mesh_2d* mesh, qaws_scalar const* values, qaws_attribute_desc const* desc);
qaws_status qaws_mesh_2d_set_attribute(qaws_mesh_2d* mesh, unsigned int attribute, qaws_scalar const* values);
qaws_status qaws_mesh_2d_remove_attribute(qaws_mesh_2d* mesh, unsigned int attribute);

#endif /* QAWS_TRIANGULATE_H */
