# Inspection and Analysis

`qaws_inspect.h` queries curve and surface properties and runs geometric
analysis. Names below are the unsuffixed forms ([precision.md](precision.md)).
Every function is thread-safe on immutable inputs. Intersections of many
curves at once are in `qaws_curve_batch.h`
([views_and_batches.md](views_and_batches.md)); exact predicates and
certified exact queries are in `qaws_exact.h`
([exact_kernel.md](exact_kernel.md)).

## Generic curve properties

```c
qaws_curve_kind kind     = qaws_curve_get_kind(curve);
qaws_dimension  dim      = qaws_curve_get_dimension(curve);
unsigned int    degree   = qaws_curve_get_degree(curve);
unsigned int    spans    = qaws_curve_get_span_count(curve);
qaws_range      range    = qaws_curve_get_parameter_range(curve);
int             closed   = qaws_curve_is_closed(curve);
int             periodic = qaws_curve_is_periodic(curve);
int             rational = qaws_curve_is_rational(curve);
qaws_continuity cont     = qaws_curve_get_continuity(curve);
```

`kind` is any `qaws_curve_kind` (the 13 families, `REPARAMETERIZED`,
`VIEW`). `qaws_curve_get_span_continuity(&c, curve, boundary_index)` gives the
continuity at one span boundary (0 to span_count - 2).

## Arc length

```c
qaws_scalar length, t;
qaws_curve_compute_arc_length(&length, curve, range.min_value, range.max_value);
qaws_curve_compute_parameter_at_arc_length(&t, curve, range.min_value, length / 2);
```

The length is adaptive Gauss-Kronrod 7-15 on each span, halving until the
Gauss and Kronrod estimates agree to `QAWS_EPSILON` of the length. The inverse
uses safeguarded Newton on the same measure; it returns
`QAWS_STATUS_OUT_OF_RANGE` for a negative length or one past the end.

## Bounds

`qaws_curve_compute_bounds_2d/3d` return the curve's true box for the
families with Bezier spans (Bezier, rational Bezier, B-spline, NURBS, views
over whole spans): each span's ends and the points where a coordinate's
derivative vanishes, found as the roots of that derivative's Bernstein
polynomial and evaluated on the curve (exact up to rounding). The other
families sample 64 points per span (at least 64) and return the samples'
box. For ranges guaranteed through rounding use the enclosures below.
`qaws_surface_compute_bounds` gives a surface's box.

## Closest points

```c
qaws_scalar t, dist, gap;
qaws_curve_find_closest_parameter_2d(&t, curve, point);
qaws_curve_find_closest_point_certified_2d(&t, &dist, &gap, curve, point);
```

- **Certified** (`qaws_curve_find_closest_point_certified_2d/3d`): for curves
  with Bezier spans (Bezier, rational Bezier, B-spline, NURBS, borrowed
  splines), branch and bound over the spans' Bezier forms with convex-hull and
  fat-line bounds. The global minimum distance lies in
  `[dist - gap, dist]`; the gap is within 1e-9 of the curve's extent
  (`QAWS_EPSILON` in f32) unless the work bound (4096 halvings) is reached.
  `QAWS_STATUS_UNSUPPORTED_OPERATION` for other families and non-positive
  weights.
- `qaws_curve_find_closest_parameter_2d/3d` use the certified search where it
  applies; the other families use sampling followed by Newton refinement.
- Surfaces: `qaws_surface_find_closest_point` (an 8 x 8 grid and Newton,
  then certified where the surface has Bezier patches) and
  `qaws_surface_find_closest_point_certified` (Bezier, B-spline, NURBS and
  borrowed surfaces of degree at most 15; work bound 8192 splits).

## Differential geometry of curves

| Function | Result |
|---|---|
| `qaws_curve_compute_tangent_2d/3d` | Unit tangent |
| `qaws_curve_compute_normal_2d` | Unit normal, the tangent turned 90 degrees counter-clockwise |
| `qaws_curve_compute_curvature_2d/3d` | Signed curvature in 2D, `|d1 x d2| / |d1|^3` in 3D |
| `qaws_curve_compute_torsion_3d` | Torsion |
| `qaws_curve_compute_speed` | `|d1|` |
| `qaws_curve_compute_frenet_frame_3d` | Unit tangent, principal normal, binormal |
| `qaws_curve_compute_rotation_minimizing_frames_3d` | Frames at many parameters |

`qaws_curve_compute_rotation_minimizing_frames_3d(out_frames, curve,
parameters, initial_normal, count)` writes 9 scalars per non-decreasing
parameter (unit tangent, normal, binormal). Each normal follows from the
previous by double reflection (`core/qaws_frame_core.h`), fourth-order
accurate in the step; the first is `initial_normal` made perpendicular to the
tangent (NULL: chosen from the tangent). It is defined where the Frenet frame
is not (straight pieces, inflections); a closed curve's frame need not close.
`QAWS_STATUS_DEGENERATE_CURVE` where the speed vanishes. No allocation.

## Features and topology

```c
qaws_scalar params[32];
unsigned int n;
qaws_curve_find_inflection_points(params, &n, curve, 32);
qaws_curve_find_extrema(params, &n, curve, 32, /* axis 0 x, 1 y, 2 z */ 0);

qaws_curvature_sample_2d comb[100];   /* position, signed curvature, normal */
qaws_curve_compute_curvature_comb_2d(comb, curve, 100, 100);

int winding;
qaws_curve_compute_winding_number_2d(&winding, curve, point);   /* 0: outside */

qaws_scalar h;
qaws_curve_compute_hausdorff_distance_2d(&h, curve_a, curve_b, 256);  /* by samples */
```

The comb and Hausdorff functions have `_3d` forms. For winding numbers and
point location in regions of several paths, see `qaws_path.h`
([regions_2d.md](regions_2d.md)).

## Enclosures (guaranteed ranges)

```c
double lo[9], hi[9];
qaws_curve_enclose(lo, hi, curve, t0, t1, /* order */ 2);
qaws_surface_enclose(lo, hi, surface, u0, u1, v0, v1, /* order */ 1);
```

For every parameter in the interval (box), the value lies in `[lo, hi]`, and
with `order` 1 or more the derivatives too (curves up to the second, surfaces
up to the third). The ranges come from the Bezier nets of the interval, every
operation rounded outward: wider than the true ranges, never narrower, and
shrinking with the interval.

- Curves: Bezier, rational Bezier, B-spline, NURBS (positive weights),
  Hermite, uniform Catmull-Rom and borrowed splines of degree at most 15.
- Surfaces: Bezier, B-spline, NURBS (positive weights) and borrowed spline
  surfaces of degree at most 15. Second and third derivatives are unbounded
  across an interior knot line where the surface is less than C2.
- `QAWS_STATUS_UNSUPPORTED_OPERATION` for other kinds,
  `QAWS_STATUS_NUMERICAL_FAILURE` where a weight's range reaches 0,
  `QAWS_STATUS_INVALID_ARGUMENT` for an empty, non-finite or out-of-range
  interval or an order too high.

## Surface analysis

| Function | Result |
|---|---|
| `qaws_surface_compute_area` | Area |
| `qaws_surface_compute_gaussian_curvature`, `qaws_surface_compute_mean_curvature` | Curvature at (u, v) |
| `qaws_surface_compute_principal_curvatures` | `qaws_surface_curvature_result`: gaussian, mean, kappa1, kappa2 |
| `qaws_surface_project_curve` | A 3D curve projected to a 2D B-spline in (u, v) by closest points (`sample_count` 0: 128) |

Surface tessellation (`qaws_surface_tessellate`,
`qaws_surface_tessellate_adaptive`) is described in
[sampling.md](sampling.md).

## Family-specific inspection

```c
qaws_vec2 cps[10];
unsigned int count;
qaws_curve_get_control_points(cps, &count, curve, 10);   /* Bezier, rational Bezier, B-spline, NURBS */
qaws_bezier_get_control_points(cps, &count, curve, 10);  /* Bezier only */

qaws_scalar knots[64], weights[32];
unsigned int knot_count, weight_count;
qaws_bspline_get_knots(knots, &knot_count, curve, 64);
qaws_nurbs_get_weights(weights, &weight_count, curve, 32);
```

A curve of another kind gives `QAWS_STATUS_UNSUPPORTED_OPERATION`; a short
buffer `QAWS_STATUS_BUFFER_TOO_SMALL` with the count. Borrowed curves refuse
these: the host holds the arrays.

## Types declared here

- `qaws_plane` (a point and a normal): the plane of the differentiable
  curve-plane intersection in `qaws_diff_ops.h`
  ([differentiation.md](differentiation.md)).
- `qaws_intersection_2d/3d` (`parameter_a`, `parameter_b`, `position`) and
  `qaws_surface_curve_intersection` (`u`, `v`, `t`, `position`): declared
  here; no public function in this header returns them. Curve intersections
  are found by the batch API in `qaws_curve_batch.h`.
