# API Reference: qaws_path.h

2D paths and regions of curves: construction helpers, measures, point
location, and Clipper2's path utilities on point arrays. Paths are the
operands of [clip.md](clip.md) and the 2D [shatter.md](shatter.md).

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)). Thread-safe on immutable curves.

---

## Types

```c
typedef struct qaws_path_2d {
	qaws_curve const* const* curves;
	unsigned int curve_count;
	int closed;
} qaws_path_2d;
```

A **path** is a list of 2D curves joined end to end: each curve ends where
the next one starts, and a closed path's last curve ends where its first
one starts. A polygon is a path of one degree-1 B-spline. The path borrows
its curves.

A **region** is a list of closed paths read with a fill rule, as in
Clipper2: a point is inside when the total winding number of all the paths
around it passes the rule.

```c
typedef enum qaws_fill_rule {
	QAWS_FILL_EVEN_ODD = 0,   /* odd winding numbers are inside */
	QAWS_FILL_NON_ZERO,       /* non-zero winding numbers are inside */
	QAWS_FILL_POSITIVE,       /* positive winding numbers are inside */
	QAWS_FILL_NEGATIVE        /* negative winding numbers are inside */
} qaws_fill_rule;

typedef enum qaws_point_location {
	QAWS_POINT_OUTSIDE = 0,
	QAWS_POINT_INSIDE,
	QAWS_POINT_ON             /* within the tolerance of a boundary */
} qaws_point_location;
```

```c
int qaws_fill_rule_inside(qaws_fill_rule rule, int winding);
```

1 when winding number `winding` is inside under `rule`, else 0.

---

## Construction

```c
qaws_status qaws_curve_create_polyline_2d(qaws_curve** out_curve, qaws_scalar const* points,
	unsigned int point_count, int closed);
qaws_status qaws_curve_create_ellipse_2d(qaws_curve** out_curve, qaws_scalar center_x, qaws_scalar center_y,
	qaws_scalar radius_x, qaws_scalar radius_y, qaws_scalar rotation);
qaws_status qaws_curve_transform_2d(qaws_curve** out_curve, qaws_curve const* curve, qaws_scalar const m[6]);
```

Each creates a new curve the caller frees with `qaws_curve_destroy`.

- **Polyline**: a degree-1 B-spline through the points (2 scalars each),
  with parameter i at point i. `closed` = 1 adds the edge back to the first
  point, so do not repeat it.
- **Ellipse**: an exact closed rational quadratic NURBS of four quarter
  arcs, counter-clockwise from center + rx (cos rot, sin rot). rx == ry
  gives a circle.
- **Transform**: the curve under `x' = m[0] x + m[1] y + m[2]`,
  `y' = m[3] x + m[4] y + m[5]`.
  - Bezier, rational Bezier, B-spline, NURBS and polynomial curves keep
    their kind.
  - Arcs and clothoids keep theirs under a similarity (rotation, uniform
    scale, reflection, translation).
  - Under other maps an arc becomes an exact NURBS and a clothoid a cubic
    B-spline fitted to 1e-9 of its extent.
  - Other kinds become their exact `qaws_curve_extract` form first.

## Measures

```c
qaws_status qaws_path_compute_area_2d(qaws_scalar* out_area, qaws_path_2d const* path);
qaws_status qaws_region_compute_area_2d(qaws_scalar* out_area, qaws_path_2d const* paths, unsigned int path_count);
int         qaws_path_is_positive_2d(qaws_path_2d const* path);
qaws_status qaws_path_compute_length_2d(qaws_scalar* out_length, qaws_path_2d const* path);
qaws_status qaws_path_compute_bounds_2d(qaws_vec2* out_min, qaws_vec2* out_max, qaws_path_2d const* path);
```

- **Area**: the signed integral of x dy around the path, positive
  counter-clockwise (x right, y up). An open path is closed by the straight
  line back to its start. Gauss-Legendre per span: exact for polynomial
  spans, adaptive for the others.
- **Region area**: the sum of the signed areas of the paths.
- **Is positive**: 1 when the closed path's signed area is >= 0 (Clipper2's
  IsPositive).
- **Length**: arc length by adaptive Gauss-Legendre per span (exact for
  lines, to rounding for the other kinds).
- **Bounds**: tight, from the end points and every point where a
  coordinate's derivative vanishes (bracketed by sampling each span, then
  refined).
- Statuses: `QAWS_STATUS_INVALID_ARGUMENT` on a NULL pointer or a malformed
  path.

## Point location

```c
qaws_status qaws_path_compute_winding_2d(int* out_winding, qaws_path_2d const* path, qaws_vec2 point);
qaws_status qaws_region_locate_point_2d(qaws_point_location* out_location, qaws_path_2d const* paths,
	unsigned int path_count, qaws_fill_rule fill_rule, qaws_vec2 point, qaws_scalar tolerance);
```

- **Winding**: the winding number of a closed path around the point.
  `QAWS_STATUS_INVALID_ARGUMENT` when the point lies on the path (within
  the default tolerance).
- **Locate**: `QAWS_POINT_ON` within `tolerance` of a boundary (0 = 1e-10 of
  the region's extent, 1e-5 in float builds), else `QAWS_POINT_INSIDE` or
  `QAWS_POINT_OUTSIDE` by the fill rule on the total winding number.

The winding number comes from flattened chords: a chord close enough to the
point to pass on the wrong side of it is split until it cannot, so the
angle the chords sweep around the point is the angle the curves sweep. An
open path is closed by its chord back to the start.

## Polylines

Clipper2's path utilities on arrays of 2D points. Each takes `count` points
(2 scalars each) and writes at most `count` points to `out`, which may equal
`points`, with the kept count in `*out_count`.

```c
qaws_status qaws_polyline_strip_duplicates_2d(qaws_scalar* out, unsigned int* out_count,
	qaws_scalar const* points, unsigned int count, int closed);
qaws_status qaws_polyline_trim_collinear_2d(qaws_scalar* out, unsigned int* out_count,
	qaws_scalar const* points, unsigned int count, int closed);
qaws_status qaws_polyline_simplify_2d(qaws_scalar* out, unsigned int* out_count,
	qaws_scalar const* points, unsigned int count, int closed, qaws_scalar epsilon);
qaws_status qaws_polyline_rdp_2d(qaws_scalar* out, unsigned int* out_count,
	qaws_scalar const* points, unsigned int count, qaws_scalar epsilon);
```

- **Strip duplicates**: consecutive repeated points removed; `closed` also
  drops a last point equal to the first.
- **Trim collinear**: points on the line through their neighbours removed,
  and spikes too (a point whose neighbours lie on the same side along the
  line), by an exact orientation test. Fewer than 3 points left (2 for
  open) gives 0.
- **Simplify** (Clipper2's SimplifyPath): repeatedly removes the point
  nearest the line through its neighbours while that distance is at most
  `epsilon`.
- **RDP** (Ramer-Douglas-Peucker): keeps the ends, and recursively the point
  farthest from the chord while it is farther than `epsilon`.
