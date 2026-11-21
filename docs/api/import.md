# API Reference: qaws_import.h

Curves from raw point arrays (polylines): an interpolating Catmull-Rom spline or a cubic trajectory through the
points. Each function creates a new curve with the default allocator, owned by the caller (destroy it with
`qaws_curve_destroy`, see [curve.md](curve.md)).

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Related: [export.md](export.md) (the reverse direction), [catmull_rom.md](catmull_rom.md),
[trajectory.md](trajectory.md).

---

## Catmull-Rom from points

```c
qaws_status qaws_polyline_import_catmull_rom(
    qaws_curve** out_curve,
    void const* points,
    qaws_dimension dimension,
    unsigned int point_count,
    qaws_parameterization parameterization,
    int closed);
```

A Catmull-Rom spline interpolating `points`, the same as `qaws_curve_create_catmull_rom` with these fields.

- `points`: flat array of `point_count` `qaws_vec2` or `qaws_vec3`, as `dimension` says (`QAWS_DIMENSION_2D` or
  `QAWS_DIMENSION_3D`). The points are copied.
- `parameterization`: knot spacing; `QAWS_PARAMETERIZATION_DEFAULT` (0) is centripetal.
- `closed`: non-zero for a closed loop.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer or `point_count < 2`, otherwise the creation status.

---

## Trajectory from points

```c
qaws_status qaws_polyline_import_trajectory(
    qaws_curve** out_curve,
    void const* points,
    qaws_dimension dimension,
    unsigned int point_count,
    int closed);
```

A cubic trajectory through `points` with uniform timing: key `i` at time `i` (so the domain is
`[0, point_count - 1]`), created with `qaws_curve_create_trajectory`.

- `points`: flat array of `point_count` `qaws_vec2` or `qaws_vec3`. `closed`: non-zero for a closed loop.
- Past 64 points a temporary time array is allocated (`QAWS_STATUS_ALLOCATION_FAILURE` if that fails).
- Returns `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer or `point_count < 2`, otherwise the creation status.
