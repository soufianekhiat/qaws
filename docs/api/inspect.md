# API Reference: qaws_inspect.h

Curve and surface inspection: properties, arc length, bounds, closest points (sampled and certified), differential
geometry (tangents, curvature, torsion, frames), feature detection, control data, guaranteed enclosures,
Hausdorff distance, surface tessellation and curve projection. All functions are thread-safe on immutable curves
and surfaces.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Related: [eval.md](eval.md), [surface.md](surface.md), concepts in [../inspection.md](../inspection.md).

---

## Generic inspection

```c
qaws_curve_kind qaws_curve_get_kind(qaws_curve const* curve);
qaws_dimension  qaws_curve_get_dimension(qaws_curve const* curve);
unsigned int    qaws_curve_get_degree(qaws_curve const* curve);
unsigned int    qaws_curve_get_span_count(qaws_curve const* curve);
qaws_range      qaws_curve_get_parameter_range(qaws_curve const* curve);
int             qaws_curve_is_closed(qaws_curve const* curve);
int             qaws_curve_is_periodic(qaws_curve const* curve);
int             qaws_curve_is_rational(qaws_curve const* curve);
qaws_continuity qaws_curve_get_continuity(qaws_curve const* curve);
```

| Function | Returns | NULL curve |
|---|---|---|
| `qaws_curve_get_kind` | the curve family (`QAWS_CURVE_KIND_*`) | `QAWS_CURVE_KIND_INVALID` |
| `qaws_curve_get_dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` | 0 |
| `qaws_curve_get_degree` | the polynomial degree | 0 |
| `qaws_curve_get_span_count` | the number of spans | 0 |
| `qaws_curve_get_parameter_range` | the parameter domain | `{0, 0}` |
| `qaws_curve_is_closed` | non-zero if the curve forms a closed loop | 0 |
| `qaws_curve_is_periodic` | non-zero if the curve is periodic | 0 |
| `qaws_curve_is_rational` | non-zero for weighted (rational) curves | 0 |
| `qaws_curve_get_continuity` | the overall continuity class (`QAWS_CONTINUITY_C0` .. `C3`) | `QAWS_CONTINUITY_C0` |

```c
qaws_status qaws_curve_get_span_continuity(
    qaws_continuity* out_continuity,
    qaws_curve const* curve,
    unsigned int boundary_index);
```

Continuity at a span boundary, `boundary_index` from 0 to span_count - 2 (`QAWS_STATUS_INVALID_ARGUMENT` otherwise,
or for a null pointer). It reports the curve's overall continuity class; it does not measure the derivatives at
that boundary.

---

## Arc length

```c
qaws_status qaws_curve_compute_arc_length(
    qaws_scalar* out_length,
    qaws_curve const* curve,
    qaws_scalar parameter_min,
    qaws_scalar parameter_max);

qaws_status qaws_curve_compute_parameter_at_arc_length(
    qaws_scalar* out_parameter,
    qaws_curve const* curve,
    qaws_scalar parameter_min,
    qaws_scalar length);
```

- `compute_arc_length`: the length between two parameters, by adaptive Gauss-Kronrod (7-15) on each span. An
  interval is halved until its Gauss and Kronrod estimates agree to `QAWS_EPSILON` of the length. Does not allocate.
- `compute_parameter_at_arc_length`: the parameter at arc length `length` after `parameter_min`.
  - Method: safeguarded Newton on the length (the speed is its derivative) inside a bracket every step narrows; a
    step that leaves the bracket or does not halve the previous one is replaced by a bisection. Cusps and
    stationary spans are handled (inside a stationary span any of its parameters may be returned).
  - A `length` of 0 returns `parameter_min`; the whole remaining length returns the range's end.
  - `QAWS_STATUS_OUT_OF_RANGE` for a negative length, a length past the end of the range, or a `parameter_min`
    outside the range.
  - Each call integrates from `parameter_min`; for increasing lengths, start each call at the previous parameter
    and pass the difference.

---

## Bounds

```c
qaws_status qaws_curve_compute_bounds_2d(
    qaws_vec2* out_min,
    qaws_vec2* out_max,
    qaws_curve const* curve);

qaws_status qaws_curve_compute_bounds_3d(
    qaws_vec3* out_min,
    qaws_vec3* out_max,
    qaws_curve const* curve);

qaws_status qaws_surface_compute_bounds(
    qaws_vec3* out_min,
    qaws_vec3* out_max,
    qaws_surface const* surface);
```

Axis-aligned boxes. Curves with Bezier spans (Bezier, rational Bezier, B-spline, NURBS, views over whole spans): the
true box, from each span's ends and the points where a coordinate's derivative vanishes (the roots of the derivative's
Bernstein polynomial, `X'` or `X' W - X W'`, evaluated on the curve). Other curves: 64 samples per span (at least 64)
uniformly over the domain; surfaces on a 32 x 32 grid. A sampled box can be slightly smaller than the true one; for
ranges guaranteed through rounding use `qaws_curve_enclose` / `qaws_surface_enclose` below. Returns `QAWS_STATUS_INVALID_ARGUMENT`,
`QAWS_STATUS_INVALID_DIMENSION` (curves), or an evaluation status.

---

## Closest points on curves

```c
qaws_status qaws_curve_find_closest_parameter_2d(
    qaws_scalar* out_parameter,
    qaws_curve const* curve,
    qaws_vec2 point);

qaws_status qaws_curve_find_closest_parameter_3d(
    qaws_scalar* out_parameter,
    qaws_curve const* curve,
    qaws_vec3 point);

qaws_status qaws_curve_find_closest_point_certified_2d(
    qaws_scalar* out_parameter,
    qaws_scalar* out_distance,
    qaws_scalar* out_gap,
    qaws_curve const* curve,
    qaws_vec2 point);

qaws_status qaws_curve_find_closest_point_certified_3d(
    qaws_scalar* out_parameter,
    qaws_scalar* out_distance,
    qaws_scalar* out_gap,
    qaws_curve const* curve,
    qaws_vec3 point);
```

- `find_closest_parameter`: 32 uniform samples, then Newton refinement. For curves with Bezier spans (Bezier,
  rational Bezier, B-spline, NURBS and borrowed splines) the result is then certified as below; other families keep
  the sampled answer, which can be a local minimum.
- `find_closest_point_certified`: the global minimum distance lies in `[*out_distance - *out_gap, *out_distance]`.
  - Method: branch and bound over each span's Bezier form (blossomed from the span's B-spline window). A piece is
    dropped when a lower bound of its distance (to the box of its control points, or to its fat line: the chord's
    segment widened by the points' largest offset) is no less than the best distance found. Pieces are halved by
    de Casteljau until their box is under 1e-9 of the curve's extent (`QAWS_EPSILON` in float); a Newton step on
    `(C - p) . C'` finishes. Allocation-free.
  - Work bound: after 4096 halvings the remaining pieces end the search and widen the gap. A point at a circle's
    center ends this way: the distance is right and the gap says the parameter is not unique.
  - `QAWS_STATUS_UNSUPPORTED_OPERATION` for other families (Hermite, Catmull-Rom, polynomial, arc, clothoid,
    Yuksel, trajectory, composite, subdivision, wrappers) and for non-positive weights;
    `QAWS_STATUS_INVALID_ARGUMENT` for a null output.

---

## Differential geometry of curves

```c
qaws_status qaws_curve_compute_tangent_2d(qaws_vec2* out_tangent, qaws_curve const* curve, qaws_scalar parameter);
qaws_status qaws_curve_compute_tangent_3d(qaws_vec3* out_tangent, qaws_curve const* curve, qaws_scalar parameter);
qaws_status qaws_curve_compute_curvature_2d(qaws_scalar* out_curvature, qaws_curve const* curve, qaws_scalar parameter);
qaws_status qaws_curve_compute_curvature_3d(qaws_scalar* out_curvature, qaws_curve const* curve, qaws_scalar parameter);
qaws_status qaws_curve_compute_torsion_3d(qaws_scalar* out_torsion, qaws_curve const* curve, qaws_scalar parameter);
qaws_status qaws_curve_compute_speed(qaws_scalar* out_speed, qaws_curve const* curve, qaws_scalar parameter);
qaws_status qaws_curve_compute_normal_2d(qaws_vec2* out_normal, qaws_curve const* curve, qaws_scalar parameter);
qaws_status qaws_curve_compute_frenet_frame_3d(
    qaws_vec3* out_tangent,
    qaws_vec3* out_normal,
    qaws_vec3* out_binormal,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

| Function | Result | Where the speed is below 1e-12 |
|---|---|---|
| `compute_tangent_2d/3d` | unit `d1` | zero vector |
| `compute_curvature_2d` | signed `(d1.x d2.y - d1.y d2.x) / abs(d1)^3`, positive turning counter-clockwise | 0 |
| `compute_curvature_3d` | `abs(d1 x d2) / abs(d1)^3` | 0 |
| `compute_torsion_3d` | `((d1 x d2) . d3) / abs(d1 x d2)^2` | 0 when `abs(d1 x d2)^2 < 1e-24` |
| `compute_speed` | `abs(d1)`, 2D or 3D | - |
| `compute_normal_2d` | unit tangent turned 90 degrees counter-clockwise | zero vector |
| `compute_frenet_frame_3d` | unit T, N = B x T, B = unit(d1 x d2) | all zero |

On straight pieces and inflections (`d1 x d2` near zero) the Frenet frame picks a binormal perpendicular to the
tangent from the coordinate axis least aligned with it. All return `QAWS_STATUS_INVALID_ARGUMENT` for a null
pointer and `QAWS_STATUS_INVALID_DIMENSION` for the wrong dimension (`compute_speed` takes either).

### Rotation-minimizing frames

```c
qaws_status qaws_curve_compute_rotation_minimizing_frames_3d(
    qaws_scalar* out_frames,
    qaws_curve const* curve,
    qaws_scalar const* parameters,
    qaws_vec3 const* initial_normal,
    unsigned int count);
```

Rotation-minimizing (parallel-transport, Bishop) frames at `count` non-decreasing parameters. `out_frames` gets 9
scalars per parameter: the unit tangent, the normal, then the binormal (tangent x normal), row i for
`parameters[i]`.

- Defined where the Frenet frame is not (straight pieces, inflections) and free of its twist; suited to sweeps,
  tubes and camera rails.
- Each normal follows from the previous one by double reflection (Wang, Juttler, Zheng and Liu 2008; the step is
  `qaws_double_reflection` in `core/qaws_frame_core.h`, see [core.md](core.md)). Fourth order in the parameter
  step: halving the step divides the error by about 16.
- The first normal is `initial_normal` made perpendicular to the tangent; with NULL, or a normal along the tangent,
  a perpendicular is chosen from the tangent alone.
- A closed curve's frame need not close (its holonomy).
- `QAWS_STATUS_DEGENERATE_CURVE` where the speed vanishes. No allocation.

---

## Feature detection

```c
qaws_status qaws_curve_find_inflection_points(
    qaws_scalar* out_parameters,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int parameter_capacity);

qaws_status qaws_curve_find_extrema(
    qaws_scalar* out_parameters,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int parameter_capacity,
    unsigned int axis);
```

Both scan 32 samples per span and bisect each bracketed root.

- Inflections: in 2D, sign changes of `d1 x d2`; in 3D, places where `abs(d1 x d2)` crosses 1e-12.
- Extrema: sign changes of the derivative's `axis` component (0 = x, 1 = y, 2 = z); `QAWS_STATUS_INVALID_ARGUMENT`
  for an axis not below the dimension.
- At most `parameter_capacity` parameters are written; `*out_count` is the number found. When it exceeds the
  capacity the call returns `QAWS_STATUS_BUFFER_TOO_SMALL` (the first `parameter_capacity` are written).

---

## Curvature comb

```c
typedef struct qaws_curvature_sample_2d {
    qaws_vec2 position;
    qaws_scalar curvature;
    qaws_vec2 normal;
} qaws_curvature_sample_2d;

typedef struct qaws_curvature_sample_3d {
    qaws_vec3 position;
    qaws_scalar curvature;
    qaws_vec3 normal;
} qaws_curvature_sample_3d;

qaws_status qaws_curve_compute_curvature_comb_2d(
    qaws_curvature_sample_2d* out_samples,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int sample_capacity);

qaws_status qaws_curve_compute_curvature_comb_3d(
    qaws_curvature_sample_3d* out_samples,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int sample_capacity);
```

Position, curvature and normal at `sample_count` uniform parameters (ends included), for drawing curvature combs.
Curvature and normal are as in the table above (signed in 2D, Frenet normal in 3D; zero where the speed vanishes).
`sample_count` must be at least 2 (`QAWS_STATUS_INVALID_ARGUMENT`); `QAWS_STATUS_BUFFER_TOO_SMALL` when
`sample_capacity < sample_count`.

---

## Winding number

```c
qaws_status qaws_curve_compute_winding_number_2d(
    int* out_winding_number,
    qaws_curve const* curve,
    qaws_vec2 point);
```

How many times a closed 2D curve winds around `point` (positive counter-clockwise), from the angle summed over 256
uniform samples and rounded. `QAWS_STATUS_INVALID_ARGUMENT` for a curve that is not closed.

---

## Intersection records

```c
typedef struct qaws_intersection_2d {
    qaws_scalar parameter_a;
    qaws_scalar parameter_b;
    qaws_vec2 position;
} qaws_intersection_2d;

typedef struct qaws_intersection_3d {
    qaws_scalar parameter_a;
    qaws_scalar parameter_b;
    qaws_vec3 position;
} qaws_intersection_3d;

/* Plane defined by a point and a normal vector. */
typedef struct qaws_plane {
    qaws_vec3 point;
    qaws_vec3 normal;
} qaws_plane;

/* A point where a 3D curve pierces a surface: S(u, v) = C(t). */
typedef struct qaws_surface_curve_intersection {
    qaws_scalar u;
    qaws_scalar v;
    qaws_scalar t;
    qaws_vec3 position;
} qaws_surface_curve_intersection;
```

Shared record types: a crossing of two curves (one parameter on each), a plane (used by the curve-plane functions
in [diff.md](diff.md)), and a curve-surface piercing point. No function in this header fills them.

---

## Control data

```c
qaws_status qaws_bezier_get_control_points(
    void* out_control_points,
    unsigned int* out_point_count,
    qaws_curve const* curve,
    unsigned int point_capacity);

qaws_status qaws_curve_get_control_points(
    void* out_control_points,
    unsigned int* out_point_count,
    qaws_curve const* curve,
    unsigned int point_capacity);

qaws_status qaws_bspline_get_knots(
    qaws_scalar* out_knots,
    unsigned int* out_knot_count,
    qaws_curve const* curve,
    unsigned int knot_capacity);

qaws_status qaws_nurbs_get_weights(
    qaws_scalar* out_weights,
    unsigned int* out_weight_count,
    qaws_curve const* curve,
    unsigned int weight_capacity);
```

- Points are written as `dimension` scalars each (`qaws_vec2` / `qaws_vec3` layout); capacities count points,
  knots or weights.
- `qaws_bezier_get_control_points`: Bezier curves only.
- `qaws_curve_get_control_points`: Bezier, rational Bezier, B-spline and NURBS curves, Cartesian points (weights
  through `qaws_nurbs_get_weights` or the rational Bezier desc).
- `qaws_bspline_get_knots`: B-spline and NURBS curves. `qaws_nurbs_get_weights`: NURBS curves.
- The count is always written; `QAWS_STATUS_BUFFER_TOO_SMALL` when the capacity is below it (nothing copied).
  `QAWS_STATUS_UNSUPPORTED_OPERATION` for other kinds, `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer.

---

## Surface measures and curvature

```c
qaws_status qaws_surface_compute_area(
    qaws_scalar* out_area,
    qaws_surface const* surface);

qaws_status qaws_surface_compute_gaussian_curvature(
    qaws_scalar* out_curvature,
    qaws_surface const* surface,
    qaws_scalar u,
    qaws_scalar v);

qaws_status qaws_surface_compute_mean_curvature(
    qaws_scalar* out_curvature,
    qaws_surface const* surface,
    qaws_scalar u,
    qaws_scalar v);

typedef struct qaws_surface_curvature_result {
    qaws_scalar gaussian;
    qaws_scalar mean;
    qaws_scalar kappa1;
    qaws_scalar kappa2;
} qaws_surface_curvature_result;

qaws_status qaws_surface_compute_principal_curvatures(
    qaws_surface_curvature_result* out_result,
    qaws_surface const* surface,
    qaws_scalar u,
    qaws_scalar v);
```

- Area: composite Simpson's rule of `abs(S_u x S_v)` on a 65 x 65 grid over the domain.
- Curvatures from the first (E, F, G) and second (L, M, N) fundamental forms with the unit normal
  `unit(S_u x S_v)`: Gaussian `K = (LN - M^2) / (EG - F^2)`, mean `H = (EN + GL - 2FM) / (2(EG - F^2))`, principal
  `kappa1, kappa2 = H +- sqrt(max(H^2 - K, 0))` (`kappa1 >= kappa2`). At a degenerate point (singular first form or
  vanishing normal) every value is 0.
- `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer, otherwise the evaluation status.

---

## Hausdorff distance

```c
qaws_status qaws_curve_compute_hausdorff_distance_2d(
    qaws_scalar* out_distance,
    qaws_curve const* curve_a,
    qaws_curve const* curve_b,
    unsigned int sample_count);

qaws_status qaws_curve_compute_hausdorff_distance_3d(
    qaws_scalar* out_distance,
    qaws_curve const* curve_a,
    qaws_curve const* curve_b,
    unsigned int sample_count);
```

The symmetric Hausdorff distance, estimated: `sample_count` uniform samples on each curve, each projected onto the
other with `qaws_curve_find_closest_parameter_*`; the result is the larger of the two directed maxima.
`sample_count` must be at least 2 (`QAWS_STATUS_INVALID_ARGUMENT`); both curves must have the function's dimension.

---

## Closest points on surfaces

```c
qaws_status qaws_surface_find_closest_point(
    qaws_scalar* out_u,
    qaws_scalar* out_v,
    qaws_vec3* out_closest_point,
    qaws_surface const* surface,
    qaws_vec3 point);

qaws_status qaws_surface_find_closest_point_certified(
    qaws_scalar* out_u,
    qaws_scalar* out_v,
    qaws_vec3* out_closest_point,
    qaws_scalar* out_distance,
    qaws_scalar* out_gap,
    qaws_surface const* surface,
    qaws_vec3 point);
```

- `find_closest_point`: an 8 x 8 grid and Newton on the gradient of the squared distance, then certified as below
  where the surface allows it; other surfaces keep the sampled answer, which can be a local minimum.
- `find_closest_point_certified`: the global minimum distance lies in `[*out_distance - *out_gap, *out_distance]`.
  - Method: branch and bound over the Bezier patches (each blossomed from its B-spline window). A piece is dropped
    when a lower bound of its distance is no less than the best found: the distance to the axis-aligned box of its
    net, or to its prism in the patch's own frame (the parallelogram of the corner edges in the tangent plane,
    times the net's spread along the normal). Pieces are quartered until their box is under 1e-9 of the extent
    (`QAWS_EPSILON` in float); a Newton step finishes where the Hessian is positive definite. Allocation-free.
  - Work bound: after 8192 splits the remaining pieces widen the gap.
  - Bezier, B-spline, NURBS and borrowed surfaces of degree at most 15; `QAWS_STATUS_UNSUPPORTED_OPERATION` for
    other surfaces and non-positive weights; `QAWS_STATUS_INVALID_ARGUMENT` for a null `out_distance` / `out_gap`.

---

## Guaranteed enclosures

```c
qaws_status qaws_surface_enclose(
    double* out_lo,
    double* out_hi,
    qaws_surface const* surface,
    double u0,
    double u1,
    double v0,
    double v1,
    unsigned int order);

qaws_status qaws_curve_enclose(
    double* out_lo,
    double* out_hi,
    qaws_curve const* curve,
    double t0,
    double t1,
    unsigned int order);
```

Ranges that hold for every parameter in the box or interval: wider than the true ranges, never narrower, and
shrinking with the box. Computed from the Bezier nets of the box (each knot span's B-spline window, or its
derivatives' windows, blossomed over the part of the box in the span), every operation rounded outward:
Bernstein's hull bounds the polynomial parts, interval arithmetic and the quotient rule the rational ones.

- `qaws_surface_enclose`, slots of `out_lo` / `out_hi` (3 doubles each): `[0..2]` S; order >= 1: `[3..5]` S_u,
  `[6..8]` S_v; order >= 2: `[9..11]` S_uu, `[12..14]` S_uv, `[15..17]` S_vv; order 3: `[18..20]` S_uuu,
  `[21..23]` S_uuv, `[24..26]` S_uvv, `[27..29]` S_vvv. Size the arrays for the order (up to 30).
  Second and higher derivatives are unbounded (`-INFINITY`, `INFINITY`) when a knot line strictly inside the box
  has the surface less than C2 across it (a knot of multiplicity past the degree less 2).
  Bezier, B-spline, NURBS (positive weights) and borrowed spline surfaces of degree at most 15.
- `qaws_curve_enclose`, slots of `dim` doubles: `[0 .. dim-1]` C; order >= 1: `[dim .. 2dim-1]` C'; order 2:
  `[2dim .. 3dim-1]` C''. Hermite and uniform Catmull-Rom spans use their Bezier points, in intervals.
  Bezier, rational Bezier, B-spline, NURBS (positive weights), Hermite, uniform Catmull-Rom and borrowed spline
  curves of degree at most 15.
- Errors: `QAWS_STATUS_UNSUPPORTED_OPERATION` for other surfaces / curves; `QAWS_STATUS_NUMERICAL_FAILURE` where a
  weight's range reaches 0; `QAWS_STATUS_INVALID_ARGUMENT` for a box or interval that is empty, not finite or past
  the domain, or an order past 3 (surfaces) or 2 (curves).

---

## Surface tessellation

```c
typedef struct qaws_tessellation_vertex {
    qaws_vec3 position;
    qaws_vec3 normal;
    qaws_scalar u;
    qaws_scalar v;
} qaws_tessellation_vertex;
```

Both tessellators write vertices with position, unit normal and (u, v), and `unsigned int` triangle indices
(3 per triangle) into caller buffers.

### Adaptive by tolerances

```c
typedef struct qaws_surface_tessellation_desc
{
    qaws_scalar chord_tolerance;        /* distance; 0: unused */
    qaws_scalar angle_tolerance;        /* radians; 0: unused (one of the two is needed) */
    unsigned int min_depth;             /* quarterings everywhere first; 0: 2 */
    unsigned int max_depth;             /* quarterings at most; 0: 10, at most 20 */
} qaws_surface_tessellation_desc;

qaws_status qaws_surface_tessellate_adaptive(
    qaws_tessellation_vertex* out_vertices,
    unsigned int* out_vertex_count,
    unsigned int* out_indices,
    unsigned int* out_index_count,
    qaws_surface const* surface,
    unsigned int vertex_capacity,
    unsigned int index_capacity,
    qaws_surface_tessellation_desc const* desc);
```

A crack-free triangle mesh within the tolerances, without allocation.

- Subdivision: a quadtree over (u, v). A quad splits below `min_depth`, never at `max_depth`, and otherwise when a
  test point (side midpoints, center, and the centroids in (u, v) of its two triangles) lies farther than the chord
  tolerance from those triangles, or a normal on its 3 x 3 points turns from the center's by more than half the
  angle tolerance.
- The split decision depends only on the quad, so the tree is never stored and every call rebuilds the same one.
- No cracks: a quad's sides carry the vertices of finer neighbours; a quad without such vertices is two triangles,
  one with them a fan around its center. Triangles are counter-clockwise in (u, v).
- Counts first: call with NULL buffers for the counts, then with buffers of those sizes.
  `QAWS_STATUS_BUFFER_TOO_SMALL` (with the counts) when a capacity is short.
- Vertices are written per quad: a point shared by quads repeats, with the same bits (corners sit on dyadic
  parameters). Weld by (u, v) when one index per point is needed.
- Limits: the tests sample points, so they are estimates; there is no 2:1 grading, so a coarse quad next to much
  finer ones fans out into thin triangles.
- `QAWS_STATUS_INVALID_ARGUMENT` for a null `surface`, `desc` or count pointer, a negative tolerance, or both
  tolerances 0.

### Fixed thresholds

```c
typedef struct qaws_tessellation_desc {
    unsigned int max_depth;            /* 0: 5 */
    qaws_scalar curvature_threshold;   /* 0: 0.1 */
    qaws_scalar max_edge_length;       /* 0: disabled */
} qaws_tessellation_desc;

qaws_status qaws_surface_tessellate(
    qaws_tessellation_vertex* out_vertices,
    unsigned int* out_vertex_count,
    unsigned int* out_indices,
    unsigned int* out_index_count,
    qaws_surface const* surface,
    unsigned int vertex_capacity,
    unsigned int index_capacity,
    qaws_tessellation_desc const* desc);
```

A quadtree over (u, v): a quad splits (up to `max_depth`) while its center lies farther than `curvature_threshold`
(a distance) from the average of its corners, or, when `max_edge_length > 0`, a side is longer than it. Leaves are
fanned around their centers, picking up the vertices of finer neighbours, so the mesh has no T-junctions; vertices
are shared by (u, v). Uses scratch memory. All buffers are required. When either buffer fills, the call returns
`QAWS_STATUS_BUFFER_TOO_SMALL` with the counts written so far, not the counts needed (shared vertices are found among
those written, so the total is not known in advance); size the buffers generously or use the adaptive version,
which counts first.

---

## Curve projection onto a surface

```c
qaws_status qaws_surface_project_curve(
    qaws_curve** out_uv_curve,
    qaws_surface const* surface,
    qaws_curve const* curve,
    unsigned int sample_count);
```

Projects a 3D curve onto a surface as a 2D curve in parameter space (x = u, y = v): `sample_count` uniform samples
(0: 128) are each mapped to the closest (u, v) with `qaws_surface_find_closest_point`, and a cubic B-spline is
fitted through them (`sample_count / 4` control points, at least 6). The new curve uses the input curve's
allocator; the caller destroys it with `qaws_curve_destroy`. `QAWS_STATUS_INVALID_DIMENSION` for a 2D input
curve; `*out_uv_curve` is NULL on failure.
