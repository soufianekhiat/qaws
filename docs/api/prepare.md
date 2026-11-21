# API Reference: qaws_prepare.h

Preparation helpers: each takes raw geometry input and writes precomputed flat arrays (polynomial coefficients,
derivative control points, homogeneous points, arc axes) ready to upload to any backend, for example to the
shader-side evaluators in [core.md](core.md). The caller allocates and owns every output array; no curve object is
created.

The header declares its functions only when `QAWS_BACKEND == QAWS_BACKEND_C` (see [precision.md](precision.md)).
Names are shown unsuffixed; each exists as `_f32` and `_f64`.

Cubic coefficient layout used by the Hermite, Catmull-Rom and trajectory helpers: span `s` is
`p(u) = a[s] u^3 + b[s] u^2 + c[s] u + d[s]` for a local `u` in [0, 1] (so `d` is the span's start point and `c` its
start tangent scaled to the span).

---

## Hermite

```c
qaws_status qaws_hermite_prepare_2d(
    qaws_vec2* out_a,
    qaws_vec2* out_b,
    qaws_vec2* out_c,
    qaws_vec2* out_d,
    const qaws_vec2* points,
    const qaws_vec2* tangents,
    unsigned int point_count);

qaws_status qaws_hermite_prepare_3d(
    qaws_vec3* out_a,
    qaws_vec3* out_b,
    qaws_vec3* out_c,
    qaws_vec3* out_d,
    const qaws_vec3* points,
    const qaws_vec3* tangents,
    unsigned int point_count);
```

Cubic Hermite coefficients from points and per-point tangents: `a = 2p0 + m0 - 2p1 + m1`,
`b = -3p0 - 2m0 + 3p1 - m1`, `c = m0`, `d = p0`.

- `point_count >= 2`; the outputs hold `point_count - 1` spans each.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer or `point_count < 2`.

---

## Catmull-Rom

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

Catmull-Rom spans as cubic coefficients, one span per four consecutive points (the span runs between the middle
two).

- `param`: knot spacing exponent: `QAWS_PARAMETERIZATION_UNIFORM` 0, `QAWS_PARAMETERIZATION_CENTRIPETAL` 0.5,
  `QAWS_PARAMETERIZATION_CHORDAL` 1, `QAWS_PARAMETERIZATION_DEFAULT` (and anything else) 0.5. Knots are computed
  locally for each span from the distances between its four points.
- Open (`closed == 0`): at least 4 points, `point_count - 3` spans (the first and last points only shape the ends).
  Closed: at least 3 points, `point_count` spans, indices wrapping.
- Size the outputs for `point_count` spans; `*out_span_count` receives the count written.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` (null pointer) or `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT`.

---

## Trajectory

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

Cubic Hermite spans through timed keys. Velocities come from finite differences of the keys (one-sided at the ends,
central inside, time steps floored at 1e-15) and are scaled by each span's duration, so the local parameter of span
`s` is `(t - times[s]) / (times[s + 1] - times[s])`.

- `degree` must be 3 and `key_count >= 2` (`QAWS_STATUS_INVALID_ARGUMENT` otherwise). `closed` is accepted but not
  used.
- The outputs hold `key_count - 1` spans; `*out_span_count` receives that count.
- Past 64 keys a temporary velocity array is allocated (`QAWS_STATUS_ALLOCATION_FAILURE` if that fails).

---

## Bezier derivative control points

```c
qaws_status qaws_bezier_prepare_deriv_2d(
    qaws_vec2* out_d1_cp,
    qaws_vec2* out_d2_cp,
    qaws_vec2* out_d3_cp,
    const qaws_vec2* control_points,
    unsigned int degree);

qaws_status qaws_bezier_prepare_deriv_3d(
    qaws_vec3* out_d1_cp,
    qaws_vec3* out_d2_cp,
    qaws_vec3* out_d3_cp,
    const qaws_vec3* control_points,
    unsigned int degree);
```

Control points of the hodographs of a Bezier curve with `degree + 1` control points: `out_d1_cp` gets `degree`
points `n (P[i+1] - P[i])`, `out_d2_cp` `degree - 1` points (the second derivative), `out_d3_cp` `degree - 2`
points (the third).

- Each output may be NULL to skip it. In 2D the second derivative is built from the first, so `out_d2_cp` needs
  `out_d1_cp`, and `out_d3_cp` needs both; D2 needs degree >= 2 and D3 degree >= 3.
- The 3D version fills only `out_d1_cp`; `out_d2_cp` and `out_d3_cp` are not written.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` for null `control_points` or `degree < 1`.

---

## Rational Bezier homogeneous points

```c
qaws_status qaws_rational_bezier_prepare_2d(
    qaws_scalar* out_weighted_cp,
    const qaws_vec2* control_points,
    const qaws_scalar* weights,
    unsigned int degree);

qaws_status qaws_rational_bezier_prepare_3d(
    qaws_scalar* out_weighted_cp,
    const qaws_vec3* control_points,
    const qaws_scalar* weights,
    unsigned int degree);
```

Homogeneous control points from `degree + 1` Cartesian points and weights: `(wx, wy, w)` per point in 2D
(`(degree + 1) * 3` scalars), `(wx, wy, wz, w)` in 3D (`(degree + 1) * 4` scalars). Weights are not validated.
Returns `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer.

---

## Arc axes

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

Axis vectors for evaluating an elliptical arc as `center + cos(theta) axis_u + sin(theta) axis_v`, with
`theta_start = start_angle` and `theta_range = end_angle - start_angle` (radians). `center` is not used by the
computation (the caller keeps it).

- 2D: `axis_u = (radius_x, 0)`, `axis_v = (0, radius_y)`.
- 3D: an orthonormal basis of the plane perpendicular to `normal` (normalized), scaled by the radii;
  `QAWS_STATUS_INVALID_ARGUMENT` for a zero-length normal.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` for a null output.
