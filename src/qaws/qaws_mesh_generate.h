#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_GENERATE_H_F32
#define QAWS_MESH_GENERATE_H_F32
#define QAWS_MESH_GENERATE_H_BODY
#endif
#else
#ifndef QAWS_MESH_GENERATE_H_F64
#define QAWS_MESH_GENERATE_H_F64
#define QAWS_MESH_GENERATE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_GENERATE_H_BODY
#undef QAWS_MESH_GENERATE_H_BODY
#define QAWS_MESH_GENERATE_H

#include "qaws_mesh_reference.h"
#include "qaws_shell_region.h"

/*
 * Meshes of exact geometry: a surface (NURBS, any qaws surface) or a shell
 * (a B-rep: faces, trimmed or not, meeting along shared edges).
 *
 *   REMESH  the shell meshed densely and conformally (qaws_shell_region_
 *           triangulate on the whole shell, half the target edge), the
 *           faces' shared edges and the free sides made features; that mesh's
 *           PN reference remeshed (triangles) or quadrangulated (quads) at
 *           the target size; every vertex then snapped onto its face's exact
 *           surface: its footpoint the face and its (u, v), its position the
 *           surface's there. Elements never straddle a shared edge. The
 *           default.
 *   STRUCTURED  an iso-parametric grid per face: each side split at
 *           equal lengths in space (so a shared edge is split alike from
 *           both faces), the nodes inside blending opposite sides' (u, v);
 *           the counts made the same across shared edges and opposite sides,
 *           each its sides' longest length over the target. With chord or
 *           angle tolerances (no target): one surface's iso-lines placed by
 *           need (each u-line where the integral of the largest |S_u| over
 *           the size along it, the size from the larger principal curvature,
 *           reaches the next whole number; the v-lines alike); a shell's
 *           counts doubled, per direction, until every cell's sides and
 *           centre meet them (uniform: the worst cell sets the spacing).
 *           Quads (QUADS and QUAD_DOMINANT; a cell at a pole a triangle) or
 *           triangles by the shorter diagonal; nodes welded where faces,
 *           seams and poles meet; every vertex exactly on its face, edges
 *           between faces features. A shell's faces must be untrimmed
 *           (UNSUPPORTED_OPERATION otherwise).
 *   PARAMETRIC  triangles meshed in each face's (u, v) and lifted: the
 *           boundary split in space (each side by equal steps of its length
 *           over the size: the target, or with tolerances the side curve's
 *           own chord and angle sizes, within [min_length, max_length]; the
 *           same points from both faces of an edge), the counts made the
 *           larger across shared edges; inside, LOCAL remeshing in (u, v)
 *           with the metric of the surface's first fundamental form over
 *           the size squared (the size from the larger principal curvature
 *           with tolerances), the boundary pinned. Every vertex exactly on
 *           its face; isotropic in space whatever the parametrization.
 *           Triangles only (UNSUPPORTED_OPERATION for quads); a shell's
 *           faces untrimmed.
 *   RESTRICTED_DELAUNAY  PARAMETRIC's mesh made Delaunay (edge flips,
 *           never across a feature or border edge) and refined as Chew and
 *           Ruppert do: a triangle with an angle under 30 degrees (but
 *           between two feature edges), a circumradius over 0.8 of the size
 *           or its centroid off the face by more than the chord tolerance
 *           gets its circumcentre, snapped onto its face; a point
 *           encroaching on a feature or border edge splits that edge at its
 *           middle on its curve instead. On return every edge off the
 *           features is Delaunay (its facing angles sum to pi at most) and
 *           every angle is at least 30 degrees but the input's between
 *           feature edges (at most 8 times the start's vertices added).
 *           Triangles only; untrimmed faces.
 *
 * Sizes: a uniform target length, or adaptive by tolerances (the chord and
 * angle tolerances of qaws_curvature_size_desc, clamped and graded) read
 * from the dense mesh's curvature. The output's footpoints are on the exact
 * geometry: patch = the face (0 for a surface), a = u, b = v. A face
 * attribute "face" holds each element's face. The report's deviation is how
 * far the snapping moved the vertices (the dense reference's error).
 */

#ifndef QAWS_ENUM_QAWS_GENERATE_METHOD
#define QAWS_ENUM_QAWS_GENERATE_METHOD
typedef enum qaws_generate_method
{
	QAWS_GENERATE_REMESH = 0,
	QAWS_GENERATE_PARAMETRIC,
	QAWS_GENERATE_RESTRICTED_DELAUNAY,
	QAWS_GENERATE_STRUCTURED
} qaws_generate_method;
#endif

#ifndef QAWS_ENUM_QAWS_ELEMENTS
#define QAWS_ENUM_QAWS_ELEMENTS
typedef enum qaws_elements
{
	QAWS_ELEMENTS_TRIANGLES = 0,
	QAWS_ELEMENTS_QUADS,                /* all quads */
	QAWS_ELEMENTS_QUAD_DOMINANT
} qaws_elements;
#endif

typedef struct qaws_generate_desc
{
	qaws_generate_method method;
	qaws_elements elements;
	qaws_scalar target_length;          /* uniform; 0: from the tolerances, else a 30th of the square root of the area */
	qaws_scalar chord_tolerance;        /* adaptive (with angle_tolerance): see qaws_curvature_size_desc */
	qaws_scalar angle_tolerance;
	qaws_scalar min_length, max_length;
	qaws_scalar gradation;
} qaws_generate_desc;

qaws_status qaws_mesh_generate_surface(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_surface const* surface,
	qaws_generate_desc const* desc);
qaws_status qaws_mesh_generate_shell(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_shell_atlas const* shell,
	qaws_generate_desc const* desc);

/*
 * Curved elements: the nodes of a generated mesh's faces at order 2 or 3
 * (P2, P3 triangles; Q2, Q3 quads), on the exact geometry it came from.
 *
 * Nodes: the mesh's vertices first (positions and footpoints as they are),
 * then order - 1 per edge, shared by the faces either side, then the faces'
 * inner nodes: a cubic triangle's centroid, a quad's (order - 1)^2; none in
 * a polygon of five corners or more. Each new node is its element's
 * straight point (the edge's at j / order, the triangle's barycentre, the
 * quad's bilinear point at (i / order, j / order) along its sides 0-1 and
 * 0-3) snapped onto the element's face (its "face" attribute), started from
 * its corners' (u, v); a node of a feature or border edge held on the
 * domain's side as its vertices are.
 *
 * Elements: per face in turn, its corners, then its edges' nodes from
 * corner k toward corner k + 1, then its inner nodes, i fastest: k * order
 * + inner per face of k corners (qaws_mesh_generate_node_counts).
 * out_locations (may be NULL): each node's face and (u, v).
 * The mesh must come from qaws_mesh_generate_surface on that surface, or
 * qaws_mesh_generate_shell on that shell (INVALID_ARGUMENT without its
 * "face" attribute). Order 1: the vertices and faces as they are.
 */
qaws_status qaws_mesh_generate_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_mesh_3d const* mesh,
	unsigned int order);
qaws_status qaws_mesh_generate_nodes_surface(qaws_scalar* out_positions, qaws_mesh_location* out_locations, unsigned int* out_element_nodes,
	qaws_mesh_3d const* mesh, qaws_surface const* surface, unsigned int order);
qaws_status qaws_mesh_generate_nodes_shell(qaws_scalar* out_positions, qaws_mesh_location* out_locations, unsigned int* out_element_nodes,
	qaws_mesh_3d const* mesh, qaws_shell_atlas const* shell, unsigned int order);

#endif /* QAWS_MESH_GENERATE_H */
