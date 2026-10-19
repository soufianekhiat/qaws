#ifndef QAWS_MESH_REPAIR_H
#define QAWS_MESH_REPAIR_H

#include "qaws_mesh.h"

/*
 * Repair for farq: messy input (triangle and polygon soups, scans, meshes
 * exported from games and modelling tools) made a closed, consistently
 * oriented triangle surface, the input farq's exact Booleans take.
 *
 * The split (docs/meshing_plan.md, decision 9): farq is exact and takes any
 * input with a defined inside (every directed edge cancels), self-
 * intersecting, non-manifold and nested input included; every tolerance and
 * quality judgement is qaws's, made here. So the repair closes and orients
 * the surface and leaves to farq what needs no tolerance: self-intersections
 * (farq's self-union resolves them exactly) and non-manifold edges.
 *
 * The policies, in this order. Tolerances are world lengths; 0 means the
 * default, relative to D, the bounding diagonal of the vertices the faces
 * use; a negative value turns the policy off.
 *
 *   validation  an index out of range, a face under 3 vertices or a
 *               position not finite: INVALID_ARGUMENT, the report naming
 *               the first such face or vertex. Faces with a vertex twice
 *               are dropped. Polygons are cut into triangles (ear clipping
 *               in their plane): the output is triangles.
 *   welding     points within weld_tolerance (default 1e-6 D; off: only
 *               identical points) merged, transitively (union-find over a
 *               grid's pairs, so a chain of close points may span more),
 *               each group at its mean; faces collapsed by it dropped.
 *   degenerate  faces lower than degenerate_tolerance (twice the area over
 *               the longest edge; default the weld tolerance, or 1e-6 D
 *               when welding is off): a needle (an edge shorter than the
 *               tolerance) has that edge collapsed; a cap (a corner on its
 *               long edge) is removed and the faces across the long edge
 *               are split at the corner, so no hole opens.
 *   duplicates  faces on the same three vertices: kept once, turned the
 *               majority's way; a pair turned opposite ways (a membrane
 *               whose two sides cancel) removed both. Run after each step.
 *   gaps        boundary vertices within gap_tolerance (default 1e-3 D) of
 *               another boundary vertex merged at their mean, or of a
 *               boundary edge: that edge's face split at the vertex (the
 *               T-junction resolved), until the boundaries meet. Vertices
 *               already joined by an edge are never merged.
 *   components  connected parts (sharing a vertex) with less area than
 *               min_component_area or fewer faces than min_component_faces
 *               removed (0: none).
 *   orientation each part made consistent by turning faces across its
 *               manifold edges (edges of two faces); a part that cannot be
 *               (a Moebius strip) is reported and left, its holes open.
 *   holes       boundary loops of perimeter up to hole_perimeter (default
 *               all) filled: the triangulation least in its largest
 *               dihedral angle, then in area, by dynamic programming (Liepa
 *               2003) for loops up to 500 vertices, a fan from the
 *               centroid for longer ones; then refined toward the
 *               surrounding edge lengths (Liepa's centroid splits and
 *               Delaunay flips) and faired (the bi-Laplacian with Liepa's
 *               scale-dependent weights, solved with the hole's border and
 *               its ring held, so a fill continues the surface's curvature).
 *               A loop through a vertex twice is filled as two loops. Faces
 *               the fills add have source_face -1.
 *   outward     each part then turned so its signed volume is positive,
 *               except a part inside another closed part (judged by the
 *               generalized winding number at its vertices): at an odd
 *               depth it is a cavity and turned inward, or removed with
 *               QAWS_NESTED_SHELLS_REMOVE.
 *
 * Left for farq, and reported: non-manifold edges (three faces or more),
 * self-intersections (not looked for), boundary edges of holes not filled.
 * `closed` says every directed edge cancels: what farq needs.
 *
 * The result is a qaws_mesh_3d of triangles, the given vertices welded,
 * unused ones dropped, with a face attribute "source_face" (the given face
 * each triangle comes from, -1 for fills; exact up to 2^24 in the f32
 * build). Its creation splits non-manifold edges and vertices into copies
 * at the same positions (qaws_mesh.h), so its faces' indices no longer
 * cancel there; farq reads the positions and its exact point keys merge
 * the copies back (or map the indices through
 * qaws_mesh_3d_get_vertex_source): the directed edges cancel again.
 */

typedef enum qaws_nested_shells
{
	QAWS_NESTED_SHELLS_KEEP = 0,        /* a part at an odd depth inside others turned inward: a cavity */
	QAWS_NESTED_SHELLS_REMOVE           /* parts inside another closed part removed */
} qaws_nested_shells;

#define QAWS_MESH_REPAIR_FLAT_FILLS 1u      /* fills refined, not faired */
#define QAWS_MESH_REPAIR_COARSE_FILLS 2u    /* fills left as the minimal triangulation */

typedef struct qaws_mesh_repair_desc
{
	qaws_scalar weld_tolerance;         /* 0: 1e-6 D; negative: identical points only */
	qaws_scalar degenerate_tolerance;   /* 0: the weld tolerance (1e-6 D when it is off); negative: off */
	qaws_scalar gap_tolerance;          /* 0: 1e-3 D; negative: off */
	qaws_scalar hole_perimeter;         /* 0: every hole; negative: none */
	qaws_scalar min_component_area;     /* 0: none removed */
	unsigned int min_component_faces;   /* 0: none removed */
	qaws_nested_shells nested_shells;
	unsigned int flags;                 /* QAWS_MESH_REPAIR_*_FILLS */
} qaws_mesh_repair_desc;

typedef struct qaws_mesh_repair_report
{
	unsigned int invalid_face;          /* INVALID_ARGUMENT: the first face out of range or under 3 vertices; else QAWS_MESH_NONE */
	unsigned int invalid_vertex;        /* INVALID_ARGUMENT: the first position not finite; else QAWS_MESH_NONE */
	unsigned int dropped_faces;         /* a vertex twice, given or after welding */
	unsigned int welded_vertices;       /* merged into another by the weld */
	unsigned int needles, caps;         /* degenerate faces: edges collapsed, caps removed (their neighbours split) */
	unsigned int duplicate_faces;       /* removed */
	unsigned int stitched_vertices;     /* boundary vertices merged or snapped onto an edge across a gap */
	unsigned int removed_components, removed_faces;   /* small parts */
	unsigned int non_orientable;        /* parts with no consistent orientation, left */
	unsigned int holes_filled, holes_left;
	unsigned int fill_faces;            /* faces the fills added */
	unsigned int flipped_faces;         /* given faces turned */
	unsigned int nested_shells;         /* parts inside another closed part (turned inward at odd depths, or removed) */
	unsigned int non_manifold_edges;    /* left for farq */
	unsigned int boundary_edges;        /* edges of one face left */
	unsigned int open_edges;            /* edges whose directed uses do not cancel */
	int closed;                         /* every directed edge cancels (open_edges 0, at least one face) */
} qaws_mesh_repair_report;

/* face_sizes NULL: triangles; desc NULL: the defaults; report may be NULL */
qaws_status qaws_mesh_repair(qaws_scalar const* positions, unsigned int vertex_count, unsigned int const* face_vertices,
	unsigned int const* face_sizes, unsigned int face_count, qaws_mesh_repair_desc const* desc, qaws_mesh_3d** out_mesh,
	qaws_mesh_repair_report* report);

#endif /* QAWS_MESH_REPAIR_H */
