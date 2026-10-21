#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_INTERSECT_H_F32
#define QAWS_MESH_INTERSECT_H_F32
#define QAWS_MESH_INTERSECT_H_BODY
#endif
#else
#ifndef QAWS_MESH_INTERSECT_H_F64
#define QAWS_MESH_INTERSECT_H_F64
#define QAWS_MESH_INTERSECT_H_BODY
#endif
#endif
#ifdef QAWS_MESH_INTERSECT_H_BODY
#undef QAWS_MESH_INTERSECT_H_BODY
#define QAWS_MESH_INTERSECT_H

#include "qaws_mesh.h"

/*
 * Self-intersections resolved numerically: every pair of triangles meeting
 * beyond their shared corners and edges is cut along their intersection,
 * so the result's triangles meet only at shared vertices and edges (a
 * conforming arrangement: the curves where the input crossed itself are
 * edges now, held by the faces of both sides). Nothing is removed: what is
 * inside is left to Booleans and winding numbers.
 *
 * Fast and robust without rational arithmetic (farq makes exact
 * arrangements): distances to planes, lines and points under the tolerance
 * count as zero, so coplanar and nearly coplanar pairs are handled
 * together (their overlaps cut each other along their edges) and a point
 * of an intersection is snapped to the input's features near it: a vertex
 * before an edge before a face. Every new point is keyed by what made it,
 * an edge through a triangle's plane, three triangles' planes, two edges of
 * coplanar triangles, and computed once, so every face holding it takes
 * the same vertex and the cuts match on both sides. Where three triangles
 * meet, the cuts cross at their common point. Each cut triangle is
 * triangulated again in its plane (constrained Delaunay) through its cuts
 * and the points on its edges.
 *
 * The cut faces' pieces are checked by exact predicates (qaws_mesh_repair's
 * test; the faces not cut are input faces the pairs already tested):
 * report.remaining counts the pairs of faces still crossing, which only
 * snapping (vertices moved by up to the tolerance) or a cut the
 * triangulation could not keep (report.unresolved) leave.
 *
 * Face attribute "source": the input face each triangle came from (faces
 * of more than 3 vertices are fanned first). Where an intersection curve
 * holds four faces the mesh gives its vertices copies (qaws_mesh.h; the
 * shared vertex by qaws_mesh_3d_get_vertex_source).
 */

typedef struct qaws_intersect_desc
{
	qaws_scalar tolerance;             /* distances under it count as zero; 0: 1e-9 D (D the diagonal of the input's box) */
} qaws_intersect_desc;

typedef struct qaws_intersect_report
{
	unsigned int pairs;                /* pairs of triangles found meeting beyond their shared corners and edges */
	unsigned int coplanar_pairs;       /* of them, overlapping in one plane */
	unsigned int points;               /* vertices added on the intersection curves */
	unsigned int split_faces;          /* input triangles cut */
	unsigned int welded_vertices;      /* input vertices within the tolerance of another */
	unsigned int unresolved;           /* cuts the triangulation could not keep */
	unsigned int remaining;            /* pairs of result faces still crossing (exact predicates) */
	qaws_scalar tolerance;             /* the tolerance used */
} qaws_intersect_report;

/* face_sizes NULL: triangles; desc and report may be NULL. INVALID_ARGUMENT for an index out of range, a face under
   3 vertices, a position not finite, a negative tolerance, no faces. */
qaws_status qaws_mesh_resolve_intersections(qaws_mesh_3d** out_mesh, qaws_intersect_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_intersect_desc const* desc);

#endif /* QAWS_MESH_INTERSECT_H */
