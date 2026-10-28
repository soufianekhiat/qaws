# Attribute promotion

`qaws_mesh_promote.h`, after Houdini's Attribute Promote: an attribute moved from one element of a mesh to another, the
values that meet at each target element combined by a method.

```c
qaws_attribute_promote_desc d = { 0 };
d.name = "color";                         /* a face attribute */
d.element = QAWS_MESH_VERTEX;             /* to the vertices */
d.method = QAWS_PROMOTE_AVERAGE;          /* by the attribute's rule: Oklab, sRGB, normals renormalized */
d.new_name = "color_vertex";              /* NULL: the same name, the original replaced */
qaws_mesh_3d_promote_attribute(mesh, &d);
```

**What meets what.**

| To | Its sources |
|---|---|
| a vertex | its corners, the faces about it, the edges at it |
| a corner | its vertex, its face, its two edges on that face |
| a face | its vertices, its corners, its edges |
| an edge | its two vertices, the corners at its ends on its faces, its faces |

**Methods.** `MAXIMUM`, `MINIMUM`, `AVERAGE`, `MODE` (the most common value, a tie the least), `MEDIAN` (an even count:
the higher middle), `SUM`, `SUM_OF_SQUARES`, `RMS`, `FIRST` and `LAST` (the source of least and greatest index). All work
per component except `AVERAGE`, which follows the attribute's rule (a `NORMAL` renormalized, colours averaged in their
space, quaternions on one side, `NEAREST` the heaviest). Houdini's Array of All and Unique Values need array attributes,
which qaws does not have.

**Options.**
- `piece`: a one-component attribute of the target element; its elements of equal value form a piece, the method runs
  over all the sources of the piece's elements and the result goes to each of them. The source and target element may
  then be the same (Houdini's "piece attribute").
- `index_name`: an attribute (NEAREST, one component) receiving the source element that gave each value (`MAXIMUM` and
  `MINIMUM`: the first component's; `MODE`, `FIRST`, `LAST`); -1 for the others.
- `new_name`, `delete_original`: the original is removed when asked, and always when the names are one.
- `weighted`: each source counts by its measure in `AVERAGE`, `SUM`, `SUM_OF_SQUARES` and `RMS`: a face its area, a
  corner its face's over the face's corner count, a vertex the sum of its corners', an edge its length. The others are
  unweighted.

**The whole mesh.** Houdini's detail class has no element here: `qaws_mesh_3d_reduce_attribute(values, &index, mesh, name,
method)` writes the method over every element into the caller's buffer; `qaws_mesh_3d_reduce_attribute_weighted` weighs
them as `weighted` does: `SUM` is the attribute's integral over the surface (along the edges for an edge attribute),
exact for values linear on each triangle, and `AVERAGE` its mean by area.

**Refused** (INVALID_ARGUMENT): no such attribute, an unknown element or method, the same element without a piece, a
piece that is not a one-component attribute of the target element, a new name another attribute holds.

**Test 186** (a cube of 12 triangles, each case against a brute force): a vertex height to faces by `AVERAGE` (exact),
face numbers to vertices by `MAXIMUM` with the index (exact), corner values to vertices by `AVERAGE` (1.8e-15), the
height to edges by `SUM` (exact), face numbers averaged per side through a piece (exact); reduce of 1 2 2 3 3 3 7 7:
`MODE` 3 from vertex 3, `MEDIAN` 3, `RMS` sqrt(134 / 8), `MAXIMUM` 7 from vertex 6; a vertex `NORMAL` averaged to faces
stays unit; the refusals; weighted: x integrated over the cube's surface 3 and its mean 0.5 (exact), the area 6, the
edges' length 12 + 6 sqrt 2, a face 1 summed onto the vertices (their faces' area) and x averaged onto the faces by area
shares (1.1e-16).

**Figure 75** (key `!`, `qaws_meshing_showcase`): a ball with a colour per face in six regions, promoted to vertices by
`AVERAGE` (the borders blend) and by `MODE` (each vertex one region's colour), and a field averaged per region through a
piece and written back to every face of it.
