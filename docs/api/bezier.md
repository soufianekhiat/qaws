# qaws_bezier.h: Bezier curves

A single-span polynomial Bezier curve of any degree from 1 to 16. Each name exists as `name_f32` and `name_f64`
([precision.md](precision.md)). For a stack-held Bezier with no allocation, see [inline.md](inline.md).

## qaws_bezier_desc

```c
typedef struct qaws_bezier_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* control_points;
    unsigned int control_point_count;
} qaws_bezier_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | 1 (linear) to `QAWS_CORE_MAX_DEGREE` (16) |
| `control_points` | `control_point_count` points, `dimension` scalars each, packed (`qaws_vec2` / `qaws_vec3` array) |
| `control_point_count` | must be `degree + 1` |

## qaws_curve_create_bezier

```c
qaws_status qaws_curve_create_bezier(
    qaws_curve** out_curve,
    qaws_bezier_desc const* desc);
```

Creates the curve. The control points are copied; the caller's array may be freed afterwards. Parameter range
[0, 1], one span. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `control_points` NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree 0 or above 16 |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | count not `degree + 1` |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is written only on success. The allocator-taking form `qaws_curve_create_bezier_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
