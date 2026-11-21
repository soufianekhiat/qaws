# 2D Paths, Regions and Clipping

This page explains the 2D region model: paths made of curves, regions read with a fill rule, the measures and point
location defined on them, and the clip types used by Boolean operations on regions. API pages:
[path](api/path.md), [clip](api/clip.md), [boolean_2d](api/boolean_2d.md).

> Names are given without their precision suffix. Each one exists as `name_f32` and `name_f64`; see
> [precision.md](precision.md).

## Paths

A `qaws_path_2d` is a list of 2D curves of any kind joined end to end: each curve ends where the next one starts. In a
closed path, the last curve also ends where the first one starts. The path borrows its curve pointers, and the curves
stay owned by the caller.

These constructors cover the common cases:

- `qaws_curve_create_polyline_2d` makes a degree-1 B-spline with parameter i at point i. A polygon is a closed path of
  one such curve.
- `qaws_curve_create_ellipse_2d` makes an exact rational quadratic NURBS, counter-clockwise.
- `qaws_curve_transform_2d` applies an affine map, keeping the curve's kind wherever that stays exact (see
  [operations.md](operations.md)).

## Regions and fill rules

A region is a list of closed paths read with a fill rule, as in Clipper2. A point is inside the region when the
winding number of all the paths around it passes the rule:

| Rule | Inside when the winding number is |
|---|---|
| `QAWS_FILL_EVEN_ODD` | odd |
| `QAWS_FILL_NON_ZERO` | non-zero |
| `QAWS_FILL_POSITIVE` | positive |
| `QAWS_FILL_NEGATIVE` | negative |

`qaws_fill_rule_inside(rule, winding)` applies a rule to a winding number. Orientation follows the usual convention
with x to the right and y up: a counter-clockwise path has positive area and winds +1 around the points it encloses.

## Measures and location

All of these work on the true curves, not on a polygon approximation.

- **Signed area.** `qaws_path_compute_area_2d` integrates x dy around a closed path, with Gauss-Legendre quadrature per
  span. The result is exact for polynomial spans and adaptive for the others. `qaws_region_compute_area_2d` sums the
  areas of a region's paths. `qaws_path_is_positive_2d` reports whether the area is non-negative.
- **Length.** `qaws_path_compute_length_2d` uses adaptive Gauss-Legendre quadrature of the speed.
- **Bounds.** `qaws_path_compute_bounds_2d` returns tight bounds: the end points, plus every point where the derivative
  of a coordinate vanishes.
- **Winding number.** `qaws_path_compute_winding_2d` returns the winding number of a closed path around a point, or
  `INVALID_ARGUMENT` when the point lies on the path.
- **Point location.** `qaws_region_locate_point_2d` returns `ON` when the point is within `tolerance` of a boundary.
  The default tolerance is 1e-10 of the region's extent (1e-5 in float builds). Otherwise it returns `INSIDE` or
  `OUTSIDE`, by the fill rule applied to the total winding number.

## Polyline utilities

Clipper2's path utilities, on flat arrays of 2D points. Each function writes at most `count` points, and the output
may be the input array.

- `qaws_polyline_strip_duplicates_2d` removes consecutive repeated points.
- `qaws_polyline_trim_collinear_2d` removes collinear points and spikes, using an exact orientation test.
- `qaws_polyline_simplify_2d` is Clipper2's SimplifyPath. It removes the point nearest the line through its neighbours
  while that distance is at most epsilon.
- `qaws_polyline_rdp_2d` is Ramer-Douglas-Peucker simplification.

## Clipping and Boolean operations

`qaws_clip.h` defines the model of Boolean operations on regions of curves, and `qaws_boolean_2d.h` defines the
`qaws_boolean_op` enum (union, intersection, difference). qaws itself provides the types and the result accessors,
but no function that runs a clip. Those types exist so that libraries built on qaws can share one description of
inputs and results.

### Inputs

A `qaws_clip_desc` holds:

- the **subjects** and **clips**: regions of closed paths, read with one fill rule;
- **open subjects**: open paths that are clipped by the closed regions;
- the **clip type**: what the result keeps.

| Clip type | A face is in the result when it is |
|---|---|
| `QAWS_CLIP_INTERSECTION` | in subject and in clip |
| `QAWS_CLIP_UNION` | in subject or in clip |
| `QAWS_CLIP_DIFFERENCE` | in subject and not in clip |
| `QAWS_CLIP_XOR` | in exactly one |

- **Tolerance.** Points closer than this become one vertex. The default is 6.4e-9 of the extent (1e-4 in float
  builds).
- **z callback.** `z_fn` is Clipper2's Z callback. It gives a z value to every output vertex where two inputs meet.
- **Executor.** `executor` lets the intersection pass run in parallel.

### The model

The header describes how the operation works:

1. All input curves go through one batched intersection, which finds crossings, touches and shared stretches.
2. Curves are cut into edges at every vertex. Edges that run along each other are merged into one edge, which carries
   the winding steps of all of them.
3. The edges form a planar graph. Each face gets its subject and clip winding numbers by walking across edges from the
   unbounded face.
4. The fill rule and the clip type decide which faces are inside the result.

The boundary between inside and outside faces, walked with the inside on the left, gives the result:

- Outer paths run counter-clockwise and holes clockwise. `QAWS_CLIP_REVERSE_SOLUTION` swaps these directions.
- Two parts that touch at a point come out as two paths.
- Nesting comes from the faces, not from point-in-polygon tests.

### The result

`qaws_clip_result` holds the result:

- **Closed paths.** Each has a parent (Clipper2's PolyTree, `QAWS_CLIP_NONE_INDEX` at the top level) and a depth.
  Holes are at odd depths.
- **Open paths.** These are the pieces of the open subjects.
- **Vertices.** Each `qaws_clip_vertex` records the input curves and parameters it lies on, and its z value.

Every output curve is a new curve that lies exactly on an input, as a `qaws_curve_extract` piece. Without
`QAWS_CLIP_PRESERVE_COLLINEAR`, runs of collinear line pieces are fused. The path views, curves and vertices belong to
the result and stay valid until `qaws_clip_result_destroy`.

## Seams and shatter of regions

To cut one region by another and classify the pieces, see [shatter.md](shatter.md).
`qaws_region_seams_2d` and `qaws_region_shatter_create_2d` work on the same paths and fill rules.

## Threading

The path and region functions only read their curves, so they are safe to call from several threads on shared curves.
