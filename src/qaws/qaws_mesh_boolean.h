#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_BOOLEAN_H_F32
#define QAWS_MESH_BOOLEAN_H_F32
#define QAWS_MESH_BOOLEAN_H_BODY
#endif
#else
#ifndef QAWS_MESH_BOOLEAN_H_F64
#define QAWS_MESH_BOOLEAN_H_F64
#define QAWS_MESH_BOOLEAN_H_BODY
#endif
#endif
#ifdef QAWS_MESH_BOOLEAN_H_BODY
#undef QAWS_MESH_BOOLEAN_H_BODY
#define QAWS_MESH_BOOLEAN_H

#include "qaws_mesh_intersect.h"
#include "qaws_boolean_2d.h"   /* qaws_boolean_op */

/*
 * Mesh Booleans, numerical and robust (farq makes exact ones): union,
 * intersection and difference of two closed triangle meshes, or the
 * self-union of one (b NULL: the region its winding number covers by at
 * least 1/2, overlaps and inner faces gone).
 *
 * The operands are resolved together (qaws_mesh_resolve_intersections:
 * every crossing cut, the pieces conforming), then each patch (pieces joined
 * by edges of two faces, bounded by the cuts) is sampled once,
 * just in front of and behind it: the generalized winding numbers of A and
 * B there (qaws_mesh_bvh.h, the operands' own triangles) say whether each
 * side is in the result. A piece between in and out is kept, turned so its
 * front faces out; the others are dropped. Where B lies on A (coplanar
 * overlaps), A's pieces stand for both. Because winding numbers decide
 * inside, operands that are not quite closed (small holes, soups,
 * overlapping parts) still give a sensible result.
 *
 * Each input face was triangulated on its own, so a flat region keeps its
 * input faces' edges (a cube's face diagonal through a cut); edges inside
 * flat regions (two faces of one operand in one plane within the
 * tolerance) are then turned toward Delaunay. The shape does not change.
 *
 * Face attribute "source": the face of the soup A then B (A's faces first,
 * B's from a's face_count) each triangle came from. report.closed: the
 * result's directed edges all cancel; it is not checked for being a
 * manifold (touching operands meet at edges or vertices, which the mesh
 * splits, qaws_mesh.h).
 */

typedef struct qaws_boolean_report
{
	qaws_intersect_report intersect;   /* the operands' resolution */
	unsigned int patches;              /* pieces joined by edges of two faces: each sampled once */
	unsigned int kept;                 /* pieces in the result */
	unsigned int flipped;              /* of them, turned (B's pieces facing into A - B) */
	unsigned int dropped;              /* pieces inside or outside the result */
	unsigned int flat_flips;           /* edges inside flat regions turned toward Delaunay */
	int closed;                        /* every directed edge of the result cancels */
} qaws_boolean_report;

/* b NULL: the self-union of a (op ignored). desc and report may be NULL. INVALID_ARGUMENT as
   qaws_mesh_resolve_intersections, or an unknown op. An empty result: OK and out_mesh NULL. */
qaws_status qaws_mesh_boolean(qaws_mesh_3d** out_mesh, qaws_boolean_report* report, qaws_mesh_3d_desc const* a,
	qaws_mesh_3d_desc const* b, qaws_boolean_op op, qaws_intersect_desc const* desc);

/*
 * The same Boolean of two meshes, their attributes carried (normals, UVs,
 * colours, material ids, ...):
 *   - matched by name (unnamed ones never match; an input's "source" is
 *     replaced by the result's); an attribute of both operands must have the
 *     same element, components and rule (INVALID_ARGUMENT otherwise); one
 *     operand's alone reads zeros on the other's faces, a NORMAL one that
 *     face's normal;
 *   - face attributes: the source face's value;
 *   - vertex and corner attributes: each result corner's value interpolated
 *     by the rule (colours in their rule's space) at its place in its source
 *     face (the face's fan triangle holding it, barycentric, clamped). Both
 *     sides of a cut keep their own values, so they come out as CORNER
 *     attributes of the same name: a cut stays sharp in the normals, two
 *     solids' colours or UV charts do not bleed into each other;
 *   - NORMAL attributes are negated on the pieces turned over (B's in A - B),
 *     TANGENT ones' handedness (w) too;
 *   - edge attributes: a result edge on an input edge (a piece of it) takes
 *     its value, the others zeros;
 *   - feature edges: those on the operands' feature edges stay features, and
 *     the cuts (edges between pieces of two source faces, on no input edge,
 *     not coplanar: where the operands meet) become features.
 * With attributes, edges inside flat regions are turned only between pieces
 * of one source face (values are never carried across faces). Footpoints
 * are not carried.
 */
qaws_status qaws_mesh_3d_boolean(qaws_mesh_3d** out_mesh, qaws_boolean_report* report, qaws_mesh_3d const* a, qaws_mesh_3d const* b,
	qaws_boolean_op op, qaws_intersect_desc const* desc);

#endif /* QAWS_MESH_BOOLEAN_H */
