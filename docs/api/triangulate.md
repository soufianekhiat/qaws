# API Reference: qaws_triangulate.h

Meshing regions of curves: constrained Delaunay triangulation, Ruppert
refinement, isotropic remeshing and quadrangulation. Boundary and
constraint vertices stay on the true curves throughout.

```c
typedef struct qaws_triangulate_desc {
	qaws_path_2d const* paths; unsigned int path_count;   /* the region */
	qaws_fill_rule fill_rule;
	qaws_path_2d const* constraints; unsigned int constraint_count;   /* curves to follow inside */
	qaws_scalar tolerance;          /* boundary chord deviation; 0 = 1e-3 of the extent */
	qaws_scalar min_angle;          /* degrees (at most 33); 0 = none */
	qaws_scalar max_area;           /* 0 = none */
	qaws_size_fn size_fn; void* size_user;   /* longest edge allowed at a point */
	unsigned int max_vertices;      /* 0 = 4 million */
} qaws_triangulate_desc;

qaws_status qaws_triangulate_region_2d(qaws_mesh_2d** out_mesh, desc);
qaws_status qaws_mesh_2d_remesh(mesh, target_length, size_fn, size_user, iterations);
qaws_status qaws_mesh_2d_quadrangulate(mesh, all_quads, smoothing_iterations);
void qaws_mesh_2d_destroy(mesh);
```

The mesh's accessors:

| What | Count | Item |
|---|---|---|
| vertices | `get_vertex_count` | `get_vertex` |
| triangles (counter-clockwise) | `get_triangle_count` | `get_triangle` |
| quads (counter-clockwise) | `get_quad_count` | `get_quad` |
| boundary edges (region on the left) | `get_boundary_edge_count` | `get_boundary_edge` |
| constraint edges | `get_constraint_edge_count` | `get_constraint_edge` |
| a vertex's curve and parameter | | `get_vertex_source` (0 for a free vertex) |
| the boundary curves, then the constraint pieces | `get_boundary_curve_count` | `get_boundary_curve` |

## Triangulation

1. **The region.** It is cleaned by one union (`qaws_clip`) with its fill
   rule. That gives nested loops of curves, with self-crossings and
   overlaps resolved.
2. **Constraints.** These are clipped to the region (`qaws_clip`, open
   subjects against the clean loops). A piece that ends on the boundary
   forces a boundary sample at that same point, so the two share a vertex.
3. **Samples.** Each curve is sampled span by span, which keeps the
   corners of polylines and the knots of splines. Sampling is adaptive to
   the chord tolerance and to the size field. Every vertex lies on its
   curve, and its curve and parameter are kept.
4. **Lattice.** Coordinates are snapped onto a lattice of 2^-48 of the
   extent. Orientation and in-circle tests are then exact: a double filter
   first, then 384-bit integers.
5. **Constrained Delaunay.** Points are inserted by Bowyer-Watson (walking
   location, cavities that do not cross constrained edges). Segments are
   recovered with Anglada's method: the crossed triangles are removed and
   the two pseudo-polygons on either side are retriangulated, Delaunay.
   - A vertex lying exactly on a segment splits that segment.
   - Two constraints that cross meet at a new vertex: the chords' crossing
     point, polished onto both curves by Newton. Both are split there.
6. **Inside and outside.** Triangles are classified by parity across
   boundary segments, starting from the super triangle. Constraint
   segments do not change the side.

## Refinement (Ruppert)

A triangle is bad when:
- its smallest angle is below `min_angle`;
- its area is over `max_area`;
- or its longest edge is over the size field at its centroid.

A bad triangle gets its circumcentre. If the circumcentre encroaches on a
boundary or constraint segment (it lies inside the segment's diametral
circle), that segment is split instead. The split point is the segment's
parameter midpoint on its curve, so the boundary converges to the true
curve.

A skinny triangle whose two edges at its small angle are both segments
keeps its angle: that is the input's own angle. Up to 33 degrees,
refinement terminates on inputs without sharp corners. Near angles of the
input below 60 degrees, triangles cluster (the classic limit), and
`max_vertices` bounds the work.

## Remeshing (Botsch-Kobbelt)

Remeshing works toward edges of `target_length`, or of the size field.
Each iteration runs four passes:
1. **Split** edges longer than 4/3 of the target. A boundary or constraint
   edge is split at its parameter midpoint on its curve.
2. **Collapse** edges shorter than 4/5 of the target, under these guards:
   - the link condition holds;
   - no triangle inverts;
   - no new edge is longer than 4/3 of the target;
   - a feature vertex moves only along its own feature edge.
3. **Flip** edges to Delaunay, never a feature edge.
4. **Smooth.** Free vertices move halfway to their neighbours' centroid.
   Feature vertices move halfway to their neighbours' parameter midpoint,
   along the curve. No triangle may invert.

Corners never move. A corner is any of these:
- a feature vertex with other than two feature edges;
- a junction of two curves;
- the seam of a closed curve;
- a kink, where the tangent jumps at a span boundary.

## Quadrangulation

1. **Candidates.** Each pair of triangles across a free edge (not a
   boundary or constraint edge) is a candidate if together they form a
   convex quad. A candidate is scored by its worst corner off 90 degrees.
2. **Greedy pairing**, best shaped first.
3. **Augmenting paths** (Blossom-quad's idea, depth 8, without
   blossoms). They re-pair chains of matched triangles so that unmatched
   ones find partners.
4. **All quads** (`all_quads`). Every quad is split into 4 and every
   triangle left into 3, through edge midpoints and centroids. Boundary
   and constraint midpoints lie on their curves.
5. **Smoothing** that keeps every face convex and never worsens its worst
   corner.

A quadrangulated mesh is not remeshed.

## Numbers

Test 91 checks these in f64 and f32:

| Case | What is checked |
|---|---|
| Square | 2 triangles, area exact |
| NURBS disc | every boundary vertex on the circle (1e-12); every free edge locally Delaunay; the area the inscribed polygon's |
| Annulus | no triangle in the hole; Euler characteristic 0 |
| Degenerate inputs | collinear and repeated points, a bow-tie under both rules, squares touching at a corner |
| Refinement | 30 degrees reached; a maximum area kept; a size field kept |
| Remeshing | disc: mean edge within 20% of the target, smallest angle above 25, area kept; annulus under a size field: topology kept |
| Constrained square | a line and a circle: no triangle across the line, its edges exactly 4 long, plain, refined and remeshed; all quads, constrained: every quad convex, none across the line |
| Crossing constraints | both diagonals: none across either, edges 8√2 long |
| Quads on a disc | quad-dominant: at most a fifth triangles; all quads; every quad convex; boundary on the circle |

Figure 8 of `qaws_better2d_showcase` draws nine of these meshes.
