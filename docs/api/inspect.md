# API Reference: qaws_inspect.h

Curve inspection, analysis, and family-specific queries.

---

## Generic inspection

### qaws_curve_get_kind

```c
qaws_curve_kind qaws_curve_get_kind(qaws_curve const* curve);
```

Returns the curve family. Returns `QAWS_CURVE_KIND_INVALID` if curve is NULL.

### qaws_curve_get_dimension

```c
qaws_dimension qaws_curve_get_dimension(qaws_curve const* curve);
```

Returns `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D`.

### qaws_curve_get_degree

```c
unsigned int qaws_curve_get_degree(qaws_curve const* curve);
```

Returns the polynomial degree.

### qaws_curve_get_span_count

```c
unsigned int qaws_curve_get_span_count(qaws_curve const* curve);
```

Returns the number of piecewise spans.

### qaws_curve_get_parameter_range

```c
qaws_range qaws_curve_get_parameter_range(qaws_curve const* curve);
```

Returns the valid parameter domain.

### qaws_curve_is_closed

```c
int qaws_curve_is_closed(qaws_curve const* curve);
```

Returns 1 if the curve forms a closed loop.

### qaws_curve_is_periodic

```c
int qaws_curve_is_periodic(qaws_curve const* curve);
```

Returns 1 if the curve is periodic (wraps around).

### qaws_curve_is_rational

```c
int qaws_curve_is_rational(qaws_curve const* curve);
```

Returns 1 for NURBS curves (weighted control points).

### qaws_curve_get_continuity

```c
qaws_continuity qaws_curve_get_continuity(qaws_curve const* curve);
```

Returns the overall continuity class (C0, C1, C2, C3).

---

## Analysis helpers

### qaws_curve_compute_arc_length

```c
qaws_status qaws_curve_compute_arc_length(
    qaws_scalar* out_length,
    qaws_curve const* curve,
    qaws_scalar parameter_min,
    qaws_scalar parameter_max);
```

Computes the arc length (distance) along the curve between two parameter values. Uses Gauss-Legendre quadrature internally.

### qaws_curve_compute_bounds_2d / 3d

```c
qaws_status qaws_curve_compute_bounds_2d(
    qaws_vec2* out_min,
    qaws_vec2* out_max,
    qaws_curve const* curve);

qaws_status qaws_curve_compute_bounds_3d(
    qaws_vec3* out_min,
    qaws_vec3* out_max,
    qaws_curve const* curve);
```

Computes an axis-aligned bounding box by dense sampling (64 samples per span).

### qaws_curve_find_closest_parameter_2d / 3d

```c
qaws_status qaws_curve_find_closest_parameter_2d(
    qaws_scalar* out_parameter,
    qaws_curve const* curve,
    qaws_vec2 point);

qaws_status qaws_curve_find_closest_parameter_3d(
    qaws_scalar* out_parameter,
    qaws_curve const* curve,
    qaws_vec3 point);
```

Finds the curve parameter closest to a given point. Uses coarse sampling (32 points) followed by Newton iteration refinement.

### qaws_curve_get_span_continuity

```c
qaws_status qaws_curve_get_span_continuity(
    qaws_continuity* out_continuity,
    qaws_curve const* curve,
    unsigned int boundary_index);
```

Queries continuity at a specific span boundary (index 0 to span_count-2).

---

## Derived geometric helpers

### qaws_curve_compute_tangent_2d / 3d

```c
qaws_status qaws_curve_compute_tangent_2d(
    qaws_vec2* out_tangent,
    qaws_curve const* curve,
    qaws_scalar parameter);

qaws_status qaws_curve_compute_tangent_3d(
    qaws_vec3* out_tangent,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

Computes the unit tangent vector at the given parameter. Returns a zero vector if the derivative magnitude is near zero.

### qaws_curve_compute_curvature_2d

```c
qaws_status qaws_curve_compute_curvature_2d(
    qaws_scalar* out_curvature,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

Computes signed curvature: `(d1.x * d2.y - d1.y * d2.x) / |d1|^3`. Positive for counter-clockwise bending, negative for clockwise.

### qaws_curve_compute_curvature_3d

```c
qaws_status qaws_curve_compute_curvature_3d(
    qaws_scalar* out_curvature,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

Computes unsigned curvature: `|d1 x d2| / |d1|^3`.

### qaws_curve_compute_torsion_3d

```c
qaws_status qaws_curve_compute_torsion_3d(
    qaws_scalar* out_torsion,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

Computes torsion (rate of change of the osculating plane): `(d1 x d2) . d3 / |d1 x d2|^2`. Only meaningful for 3D curves.

### qaws_curve_compute_speed

```c
qaws_status qaws_curve_compute_speed(
    qaws_scalar* out_speed,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

Computes the parametric speed `|d1|` (magnitude of the first derivative). Works for both 2D and 3D curves.

---

## Frenet frame

### qaws_curve_compute_normal_2d

```c
qaws_status qaws_curve_compute_normal_2d(
    qaws_vec2* out_normal,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

Computes the unit normal vector at the given parameter for a 2D curve. The normal is the tangent rotated 90 degrees counter-clockwise.

### qaws_curve_compute_frenet_frame_3d

```c
qaws_status qaws_curve_compute_frenet_frame_3d(
    qaws_vec3* out_tangent,
    qaws_vec3* out_normal,
    qaws_vec3* out_binormal,
    qaws_curve const* curve,
    qaws_scalar parameter);
```

Computes the Frenet frame (tangent, normal, binormal) at the given parameter for a 3D curve. Returns unit vectors forming an orthonormal basis.

---

## Inflection point detection

### qaws_curve_find_inflection_points

```c
qaws_status qaws_curve_find_inflection_points(
    qaws_scalar* out_parameters,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int parameter_capacity);
```

Finds parameter values where the curvature changes sign (inflection points).

**Parameters:**
- `curve` -- The curve to analyze.
- `out_parameters` -- Caller-allocated buffer for parameter values.
- `parameter_capacity` -- Size of the buffer.
- `out_count` -- Receives the number of inflection points found.

---

## Extrema detection

### qaws_curve_find_extrema

```c
qaws_status qaws_curve_find_extrema(
    qaws_scalar* out_parameters,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int parameter_capacity,
    unsigned int axis);
```

Finds parameter values where the curve reaches local extrema along a given axis.

**Parameters:**
- `curve` -- The curve to analyze.
- `axis` -- Coordinate axis (0 = x, 1 = y, 2 = z).
- `out_parameters` -- Caller-allocated buffer for parameter values.
- `parameter_capacity` -- Size of the buffer.
- `out_count` -- Receives the number of extrema found.

---

## Curvature comb

### qaws_curve_compute_curvature_comb_2d

```c
qaws_status qaws_curve_compute_curvature_comb_2d(
    qaws_curvature_sample_2d* out_samples,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int sample_capacity);
```

Computes curvature comb data (position, curvature, normal) at uniformly-spaced samples along a 2D curve. Useful for visualizing curvature distribution.

### qaws_curve_compute_curvature_comb_3d

```c
qaws_status qaws_curve_compute_curvature_comb_3d(
    qaws_curvature_sample_3d* out_samples,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int sample_capacity);
```

Same as `curvature_comb_2d` but for 3D curves.

---

## Winding number

### qaws_curve_compute_winding_number_2d

```c
qaws_status qaws_curve_compute_winding_number_2d(
    int* out_winding_number,
    qaws_curve const* curve,
    qaws_vec2 point);
```

Computes the winding number of a closed 2D curve around the given point. The winding number indicates how many times the curve wraps around the point (positive for counter-clockwise, negative for clockwise).

---

## Intersection detection

### qaws_curve_find_intersections_2d

```c
qaws_status qaws_curve_find_intersections_2d(
    qaws_intersection_2d* out_intersections,
    unsigned int* out_count,
    qaws_curve const* curve_a,
    qaws_curve const* curve_b,
    unsigned int intersection_capacity);
```

Finds all intersection points between two 2D curves. Results include parameter values on both curves and the intersection position.

### qaws_curve_find_intersections_3d

```c
qaws_status qaws_curve_find_intersections_3d(
    qaws_intersection_3d* out_intersections,
    unsigned int* out_count,
    qaws_curve const* curve_a,
    qaws_curve const* curve_b,
    unsigned int intersection_capacity);
```

Same as `find_intersections_2d` but for 3D curves.

---

## Self-intersection detection

### qaws_curve_find_self_intersections_2d

```c
qaws_status qaws_curve_find_self_intersections_2d(
    qaws_intersection_2d* out_intersections,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int intersection_capacity);
```

Finds all points where a 2D curve crosses itself. `parameter_a` and `parameter_b` in each result are distinct parameter values on the same curve that map to the same position.

### qaws_curve_find_self_intersections_3d

```c
qaws_status qaws_curve_find_self_intersections_3d(
    qaws_intersection_3d* out_intersections,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int intersection_capacity);
```

Same as `find_self_intersections_2d` but for 3D curves.

---

## Family-specific inspection

### qaws_bezier_get_control_points

```c
qaws_status qaws_bezier_get_control_points(
    void* out_control_points,
    unsigned int* out_point_count,
    qaws_curve const* curve,
    unsigned int point_capacity);
```

Extracts Bezier control points into a caller-provided buffer.

- `out_control_points` -- Buffer of `qaws_vec2` or `qaws_vec3` (matching the curve's dimension).
- `point_capacity` -- Size of the buffer in number of points.
- `out_point_count` -- Receives the actual number of control points.

Returns `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve is not a Bezier.

### qaws_bspline_get_knots

```c
qaws_status qaws_bspline_get_knots(
    qaws_scalar* out_knots,
    unsigned int* out_knot_count,
    qaws_curve const* curve,
    unsigned int knot_capacity);
```

Extracts the knot vector. Works for both B-Spline and NURBS curves.

Returns `QAWS_STATUS_UNSUPPORTED_OPERATION` for other curve types.

### qaws_nurbs_get_weights

```c
qaws_status qaws_nurbs_get_weights(
    qaws_scalar* out_weights,
    unsigned int* out_weight_count,
    qaws_curve const* curve,
    unsigned int weight_capacity);
```

Extracts NURBS weights.

Returns `QAWS_STATUS_UNSUPPORTED_OPERATION` if not a NURBS curve.
