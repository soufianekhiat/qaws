# API Reference: qaws_path.h

2D paths and regions of curves, and Clipper2's path utilities.

A **path** is a list of 2D curves joined end to end:

```c
typedef struct qaws_path_2d {
	qaws_curve const* const* curves;
	unsigned int curve_count;
	int closed;
} qaws_path_2d;
```

Each curve ends where the next one starts. A closed path's last curve also
ends where its first one starts. A polygon is a path of one degree-1 B-spline.

A **region** is a list of closed paths read with a fill rule, as in Clipper2:

| Rule | A point is inside when the total winding number w is |
|---|---|
| `QAWS_FILL_EVEN_ODD` | odd |
| `QAWS_FILL_NON_ZERO` | not 0 |
| `QAWS_FILL_POSITIVE` | > 0 |
| `QAWS_FILL_NEGATIVE` | < 0 |

---

## Construction

```c
qaws_status qaws_curve_create_polyline_2d(qaws_scalar const* points, unsigned int point_count, int closed,
	qaws_curve** out_curve);
qaws_status qaws_curve_create_ellipse_2d(qaws_scalar center_x, qaws_scalar center_y,
	qaws_scalar radius_x, qaws_scalar radius_y, qaws_scalar rotation, qaws_curve** out_curve);
qaws_status qaws_curve_transform_2d(qaws_curve const* curve, qaws_scalar const m[6], qaws_curve** out_curve);
```

- **Polyline**: a degree-1 B-spline through the points, with parameter i at
  point i. With `closed`, the edge back to the first point is added, so
  do not repeat that point.
- **Ellipse**: an exact closed rational quadratic NURBS made of four
  quarter arcs, counter-clockwise.
- **Transform**: maps the curve through `x' = m0 x + m1 y + m2`,
  `y' = m3 x + m4 y + m5`.
  - Bezier, rational Bezier, B-spline, NURBS and polynomial curves keep
    their kind, since only their control points or coefficients move.
  - Arcs and clothoids keep their kind under a similarity: rotation,
    uniform scale, reflection or translation.
  - Under any other map, an arc becomes an exact NURBS and a clothoid
    becomes a cubic B-spline fitted to 1e-9 of its extent.
  - Other kinds are first turned into their exact `qaws_curve_extract`
    form.

## Measures

```c
qaws_status qaws_path_compute_area_2d(qaws_path_2d const* path, qaws_scalar* out_area);
qaws_status qaws_region_compute_area_2d(qaws_path_2d const* paths, unsigned int path_count, qaws_scalar* out_area);
int         qaws_path_is_positive_2d(qaws_path_2d const* path);
qaws_status qaws_path_compute_length_2d(qaws_path_2d const* path, qaws_scalar* out_length);
qaws_status qaws_path_compute_bounds_2d(qaws_path_2d const* path, qaws_vec2* out_min, qaws_vec2* out_max);
```

- **Area** is the signed integral of x dy around the path. It is positive
  counter-clockwise, with x pointing right and y up. An open path is closed
  by its chord. Each span is integrated by 8- and 16-point Gauss-Legendre,
  subdividing until the two agree. This is exact for polynomial spans, and
  rational spans and arcs converge in a few steps.
- **Bounds** are tight. They come from the end points and from every point
  where a coordinate's derivative changes sign: 32 samples per span bracket
  the sign changes, and bisection refines each one.

## Point location

```c
qaws_status qaws_region_locate_point_2d(qaws_path_2d const* paths, unsigned int path_count,
	qaws_fill_rule fill_rule, qaws_vec2 point, qaws_scalar tolerance, qaws_point_location* out_location);
qaws_status qaws_path_compute_winding_2d(qaws_path_2d const* path, qaws_vec2 point, int* out_winding);
```

The location is one of:

- `QAWS_POINT_ON`: within `tolerance` of a boundary. A tolerance of 0
  means 1e-10 of the region's extent, or 1e-5 in float builds.
- `QAWS_POINT_INSIDE` or `QAWS_POINT_OUTSIDE`: decided by the fill rule
  applied to the total winding number.

How the winding number is computed:

1. The curves are flattened into chords.
2. A chord within four times its deviation of the point is split, until no
   chord can pass on the wrong side of the point.
3. The angle the chords sweep around the point is then the angle the curves
   sweep.

Test 83 resolves a point 1e-7 inside and 1e-7 outside a circle of radius 2.

## Polylines

These are Clipper2's utilities, on arrays of 2D points:

```c
qaws_status qaws_polyline_strip_duplicates_2d(points, count, closed, out, out_count);
qaws_status qaws_polyline_trim_collinear_2d(points, count, closed, out, out_count);
qaws_status qaws_polyline_simplify_2d(points, count, closed, epsilon, out, out_count);
qaws_status qaws_polyline_rdp_2d(points, count, epsilon, out, out_count);
```

- **Strip duplicates** removes consecutive repeated points. A closed
  polyline also loses a last point equal to the first.
- **Trim collinear** removes points on the line through their neighbours,
  which includes spikes. It uses the exact orientation predicate in builds
  with the exact kernel.
- **Simplify** (Clipper2's SimplifyPath) repeatedly removes the point
  nearest the line through its neighbours, while that distance is at most
  epsilon.
- **RDP** (Ramer-Douglas-Peucker) keeps the farthest point from each chord,
  recursively, while it is farther than epsilon.

`out` may be the input array.
