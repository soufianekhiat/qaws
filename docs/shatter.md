# Seams and Shatter

This page covers `qaws_shatter.h`. It works in 3D on curved surfaces and in 2D on regions of curves. Both versions
answer two questions:

- **Seams:** where do the operands cross?
- **Shatter:** what are the pieces of operand A when it is cut along those crossings, and is each piece inside or
  outside operand B?

These match the "Seam" and "Shatter" modes of Houdini's Boolean node. Signatures are in
[api/shatter.md](api/shatter.md).

> Names are given without their precision suffix. Each one exists as `name_f32` and `name_f64`; see
> [precision.md](precision.md).

## Surfaces

### Inputs

`qaws_surface_seam_desc` holds:

- **A and B.** Each is a list of surfaces, trimmed or not, such as the faces of a body or a sheet. B may be empty when
  only A-A seams are wanted.
- **`kinds`.** Which seams to compute: `QAWS_SEAM_A_B`, `QAWS_SEAM_A_A`, `QAWS_SEAM_B_B`.
- **`resolution`.** The number of grid cells along each side of a surface's (u, v) domain (64 by default).
- **`tolerance`.** The largest gap allowed between the two surfaces at a seam point (1e-9 of the extent by default,
  1e-5 in f32).

### How seams are found

1. Each surface is laid on a `resolution x resolution` grid, with two triangles per cell.
2. A grid of boxes pairs the triangles of surfaces that may meet. Each pair's crossing segment is found.
3. The crossings are chained into curves through the triangle edges they cross.
4. Newton's method on S1(u1, v1) = S2(u2, v2) moves every point onto both true surfaces, taking the least change of
   the four parameters for as long as the gap shrinks.
5. A point outside a trimmed surface's trim is dropped, and the curve is broken there.

### Seam output

Each seam is a curve that either runs between two ends (a surface border, or a point where three surfaces meet) or is
closed. A closed seam does not repeat its first point. `qaws_seam_info` gives the seam's kind, its two surfaces and
their operands, and its point count. Each `qaws_seam_point` carries a position and its (u, v) on both surfaces.
Surface 1 is from the lower operand, then the lower index. A-A and B-B seams leave out crossings along a border that
the two surfaces share.

### Shatter

`qaws_surface_shatter_create` cuts every surface of A along its A-B seams. Each piece is a trimmed surface (see
[brep.md](brep.md)) over the A surface's base. Its trim loops are the original trim and the seams, flattened into (u, v)
polygons.

- **Inside or outside.** With `b_is_solid`, B's surfaces must bound a solid with outward normals. Each piece is then
  marked inside or outside B by the winding number of B's triangles at one of its points (inside from 1/2).
  `qaws_surface_shatter_get` reports 1 for inside, 0 for outside, or -1 when B is not a solid.
- **Ownership.** The shatter owns its pieces. A piece stays valid while the shatter and A's surfaces live.

### Robustness

- **Exact coincidences.** Cases such as a grid vertex lying on the other operand, a border lying in it, or coplanar
  faces are decided the same way for every triangle pair. The library does this by moving each operand's grid a little
  along a direction of its own while the crossings are found: 1e-8 of the extent, or 1e-5 in f32.
- **Touches.** A seam no larger than a hundred such moves, such as a touch at a point, is left out.
- **Surfaces on B's boundary.** A surface lying on B's boundary is cut along B's edges there. Its part on the boundary
  is marked inside or outside as the move decides.
- **Limits.** Tangential contacts and crossings finer than the grid are not resolved. Raise `resolution` for them.

## Regions

### Inputs

`qaws_region_seam_desc` works on paths (see [regions_2d.md](regions_2d.md)):

- **A** is a region: closed paths read with `a_rule`.
- **B** is a region read with `b_rule`. With `b_is_region = 0`, B is instead a set of curves with no inside (open or
  closed paths), the way a surface is in 3D.
- **`flatness`** bounds the chords' deviation from the curves (1e-3 of the extent by default).
- **`tolerance`** bounds a crossing's gap between its two curves (1e-12 of the extent by default, 1e-6 in f32).

### How it works

Every curve is flattened into chords within `flatness`. The chords of both operands are arranged into one planar
graph, and Newton's method moves every crossing onto both true curves.

### Seams

`qaws_region_seams_2d` returns crossing points, not curves. Each point carries:

- its kind (`A_B`, `A_A`, `B_B`);
- the operand, path, curve and parameter on both curves.

Points where a path's own curves join are not seams. The call counts the points when `out_points` is NULL, and returns
`BUFFER_TOO_SMALL` when there are more than `capacity`.

### Shatter

`qaws_region_shatter_create_2d` cuts A's region along B's curves into pieces:

- **What cuts A.** A's own curves do not cut it where its region lies on both sides of them, such as crossings of its
  paths inside it.
- **Piece shape.** Each piece is a region of closed paths. The outer path runs counter-clockwise and holes run
  clockwise.
- **Exact curves.** A piece's curves are exact `qaws_curve_extract` pieces of A's and B's curves, reversed where they
  are walked backwards.
- **Inside or outside.** When B is a region, each piece is marked inside or outside B. The winding numbers carried from
  face to face across the curves decide this, not sampling, with each operand's fill rule applied on its own.
- **Ownership.** The shatter owns the paths and their curves.

### Degenerate meetings

- Collinear chords split one another where they overlap. Boundaries that run together (shared sides, identical
  regions) therefore share their edges, and every operand's winding change across them counts.
- Seams leave out a meeting where the curves run together on either side, because an overlap is not a crossing.
- A meeting where the curves only touch (a tangency, or a corner on a side) gives one seam point.
- Crossings closer together than `flatness` are not resolved.

## Errors and threading

`INVALID_ARGUMENT` is returned for no surface, a NULL surface, a path without curves, or `kinds == 0`. Shatter also
returns it when B has no surface or path. Region shatter ignores `kinds`.

All of these functions only read their inputs, so they are safe to call from several threads on shared surfaces and
curves. They take no executor and no cancel callback. Seams and shatter objects do not change after they are created,
so several threads may read them at once.
