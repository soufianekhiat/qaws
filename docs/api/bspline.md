# qaws_bspline.h: B-spline curves

A non-rational B-spline of degree 1 to 31, with a given knot vector or a generated clamped uniform one. Each name
exists as `name_f32` and `name_f64` ([precision.md](precision.md)).

## qaws_bspline_desc

```c
typedef struct qaws_bspline_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* control_points;
    unsigned int control_point_count;
    qaws_scalar const* knots;
    unsigned int knot_count;
    int is_uniform;
    int is_closed;
} qaws_bspline_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | 1 to 31 |
| `control_points` | `dimension` scalars per point, packed |
| `control_point_count` | at least `degree + 1` |
| `knots` | the knot vector, non-decreasing; NULL: a clamped uniform vector is generated |
| `knot_count` | with `knots`: must be `control_point_count + degree + 1` |
| `is_uniform` | non-zero with `knots` NULL asks for the generated vector (the same as `knots` NULL alone); ignored when `knots` is given |
| `is_closed` | not read by creation; whether the curve is closed is reported from its geometry (`qaws_curve_is_closed`, [inspect.md](inspect.md)) |

## qaws_curve_create_bspline, qaws_curve_create_bspline_ex

```c
qaws_status qaws_curve_create_bspline(
    qaws_curve** out_curve,
    qaws_bspline_desc const* desc);

qaws_status qaws_curve_create_bspline_ex(
    qaws_curve** out_curve,
    qaws_bspline_desc const* desc,
    qaws_allocator const* allocator);
```

Creates the curve. The `_ex` form takes its memory from `allocator` (NULL: the default, see [alloc.md](alloc.md)).

- Parameter range [knots[degree], knots[control_point_count]]. Each non-empty knot interval in that range is one span.
- The generated vector has `degree + 1` zeros, interior knots 1, 2, ..., and `degree + 1` copies of
  `control_point_count - degree`: the curve starts and ends at the end control points. For degree 3 and 6 points it is
  `{0, 0, 0, 0, 1, 2, 3, 3, 3, 3}`, range [0, 3], 3 spans.
- Continuity at a simple interior knot is C(degree - 1).
- Control points and knots are copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `control_points` NULL; no non-empty span |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree 0 or above 31 |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | fewer than `degree + 1` points, or a count too large to index |
| `QAWS_STATUS_INVALID_KNOT_COUNT` | `knot_count` not `control_point_count + degree + 1` (or 0) |
| `QAWS_STATUS_INVALID_KNOT_VECTOR` | a knot smaller than the one before it |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first.
