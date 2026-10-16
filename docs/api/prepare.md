# API Reference: qaws_prepare.h

Precomputation / preparation utilities (C backend only).

Each function takes raw geometry input and outputs pre-computed flat arrays suitable for upload to any backend. The caller owns the output memory.

Available only when `QAWS_BACKEND == QAWS_BACKEND_C`.

---

## qaws_hermite_prepare_2d

```c
qaws_status qaws_hermite_prepare_2d(
    qaws_vec2* out_a,
    qaws_vec2* out_b,
    qaws_vec2* out_c,
    qaws_vec2* out_d,
    const qaws_vec2* points,
    const qaws_vec2* tangents,
    unsigned int point_count);
```

Converts Hermite points and tangents into polynomial coefficients (a, b, c, d) per span.

**Parameters:**
- `points` -- Array of interpolation points.
- `tangents` -- Array of tangent vectors at each point.
- `point_count` -- Number of points (>= 2).
- `out_a`, `out_b`, `out_c`, `out_d` -- Caller-allocated arrays of `point_count - 1` elements each. Receives the cubic coefficients per span.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_hermite_prepare_3d

```c
qaws_status qaws_hermite_prepare_3d(
    qaws_vec3* out_a,
    qaws_vec3* out_b,
    qaws_vec3* out_c,
    qaws_vec3* out_d,
    const qaws_vec3* points,
    const qaws_vec3* tangents,
    unsigned int point_count);
```

Same as `qaws_hermite_prepare_2d` but for 3D data.

---

## qaws_catmull_rom_prepare_2d

```c
qaws_status qaws_catmull_rom_prepare_2d(
    qaws_vec2* out_a,
    qaws_vec2* out_b,
    qaws_vec2* out_c,
    qaws_vec2* out_d,
    unsigned int* out_span_count,
    const qaws_vec2* points,
    unsigned int point_count,
    qaws_parameterization param,
    int closed);
```

Converts Catmull-Rom interpolation points into polynomial coefficients per span.

**Parameters:**
- `points` -- Array of interpolation points.
- `point_count` -- Number of points.
- `param` -- Parameterization (knot spacing).
- `closed` -- Non-zero for a closed curve.
- `out_a`, `out_b`, `out_c`, `out_d` -- Caller-allocated coefficient arrays.
- `out_span_count` -- Receives the number of spans produced.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_catmull_rom_prepare_3d

```c
qaws_status qaws_catmull_rom_prepare_3d(
    qaws_vec3* out_a,
    qaws_vec3* out_b,
    qaws_vec3* out_c,
    qaws_vec3* out_d,
    unsigned int* out_span_count,
    const qaws_vec3* points,
    unsigned int point_count,
    qaws_parameterization param,
    int closed);
```

Same as `qaws_catmull_rom_prepare_2d` but for 3D data.

---

## qaws_trajectory_prepare_2d

```c
qaws_status qaws_trajectory_prepare_2d(
    qaws_vec2* out_a,
    qaws_vec2* out_b,
    qaws_vec2* out_c,
    qaws_vec2* out_d,
    unsigned int* out_span_count,
    const qaws_vec2* positions,
    const qaws_scalar* times,
    unsigned int key_count,
    unsigned int degree,
    int closed);
```

Converts trajectory keyframes (positions + times) into polynomial coefficients per span.

**Parameters:**
- `positions` -- Array of keyframe positions.
- `times` -- Array of keyframe times.
- `key_count` -- Number of keyframes.
- `degree` -- Polynomial degree.
- `closed` -- Non-zero for a closed trajectory.
- `out_a`, `out_b`, `out_c`, `out_d` -- Caller-allocated coefficient arrays.
- `out_span_count` -- Receives the number of spans produced.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_trajectory_prepare_3d

```c
qaws_status qaws_trajectory_prepare_3d(
    qaws_vec3* out_a,
    qaws_vec3* out_b,
    qaws_vec3* out_c,
    qaws_vec3* out_d,
    unsigned int* out_span_count,
    const qaws_vec3* positions,
    const qaws_scalar* times,
    unsigned int key_count,
    unsigned int degree,
    int closed);
```

Same as `qaws_trajectory_prepare_2d` but for 3D data.

---

## qaws_bezier_prepare_deriv_2d

```c
qaws_status qaws_bezier_prepare_deriv_2d(
    qaws_vec2* out_d1_cp,
    qaws_vec2* out_d2_cp,
    qaws_vec2* out_d3_cp,
    const qaws_vec2* control_points,
    unsigned int degree);
```

Computes derivative control points (D1, D2, D3) from Bezier control points.

**Parameters:**
- `control_points` -- Array of `degree + 1` control points.
- `degree` -- Bezier degree.
- `out_d1_cp` -- Receives first-derivative control points (`degree` elements).
- `out_d2_cp` -- Receives second-derivative control points (`degree - 1` elements).
- `out_d3_cp` -- Receives third-derivative control points (`degree - 2` elements).

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_bezier_prepare_deriv_3d

```c
qaws_status qaws_bezier_prepare_deriv_3d(
    qaws_vec3* out_d1_cp,
    qaws_vec3* out_d2_cp,
    qaws_vec3* out_d3_cp,
    const qaws_vec3* control_points,
    unsigned int degree);
```

Same as `qaws_bezier_prepare_deriv_2d` but for 3D data.

---

## qaws_rational_bezier_prepare_2d

```c
qaws_status qaws_rational_bezier_prepare_2d(
    qaws_scalar* out_weighted_cp,
    const qaws_vec2* control_points,
    const qaws_scalar* weights,
    unsigned int degree);
```

Converts rational Bezier control points and weights into homogeneous weighted control points.

**Parameters:**
- `control_points` -- Array of `degree + 1` control points.
- `weights` -- Array of `degree + 1` positive weights.
- `degree` -- Bezier degree.
- `out_weighted_cp` -- Receives `(degree + 1) * 3` scalars (wx, wy, w per point).

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_rational_bezier_prepare_3d

```c
qaws_status qaws_rational_bezier_prepare_3d(
    qaws_scalar* out_weighted_cp,
    const qaws_vec3* control_points,
    const qaws_scalar* weights,
    unsigned int degree);
```

Same as `qaws_rational_bezier_prepare_2d` but for 3D. Output is `(degree + 1) * 4` scalars (wx, wy, wz, w per point).

---

## qaws_arc_prepare_2d

```c
qaws_status qaws_arc_prepare_2d(
    qaws_vec2* out_axis_u,
    qaws_vec2* out_axis_v,
    qaws_scalar* out_theta_start,
    qaws_scalar* out_theta_range,
    qaws_vec2 center,
    qaws_scalar radius_x,
    qaws_scalar radius_y,
    qaws_scalar start_angle,
    qaws_scalar end_angle);
```

Prepares arc parameters from center, radii, and angles into axis vectors and theta range.

**Parameters:**
- `center` -- Arc center.
- `radius_x` -- Radius along the x-axis (or semi-major axis).
- `radius_y` -- Radius along the y-axis (or semi-minor axis).
- `start_angle` -- Start angle in radians.
- `end_angle` -- End angle in radians.
- `out_axis_u` -- Receives the cosine-direction axis vector.
- `out_axis_v` -- Receives the sine-direction axis vector.
- `out_theta_start` -- Receives the starting theta.
- `out_theta_range` -- Receives the theta range.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_arc_prepare_3d

```c
qaws_status qaws_arc_prepare_3d(
    qaws_vec3* out_axis_u,
    qaws_vec3* out_axis_v,
    qaws_scalar* out_theta_start,
    qaws_scalar* out_theta_range,
    qaws_vec3 center,
    qaws_vec3 normal,
    qaws_scalar radius_x,
    qaws_scalar radius_y,
    qaws_scalar start_angle,
    qaws_scalar end_angle);
```

Same as `qaws_arc_prepare_2d` but for 3D arcs lying in a plane defined by `normal`.

**Parameters:**
- `center` -- Arc center in 3D.
- `normal` -- Plane normal vector.
- `radius_x`, `radius_y` -- Elliptical radii.
- `start_angle`, `end_angle` -- Angles in radians.
- `out_axis_u`, `out_axis_v` -- Receive the plane basis vectors.
- `out_theta_start`, `out_theta_range` -- Receive the theta parameters.

**Returns:** `QAWS_STATUS_OK` on success.
