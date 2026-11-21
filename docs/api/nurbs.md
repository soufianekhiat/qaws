# qaws_nurbs.h: NURBS curves

A rational B-spline of degree 1 to 31 with an explicit knot vector and positive weights. Each name exists as
`name_f32` and `name_f64` ([precision.md](precision.md)).

## qaws_nurbs_desc

```c
typedef struct qaws_nurbs_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* control_points;
    unsigned int control_point_count;
    qaws_scalar const* knots;
    unsigned int knot_count;
    qaws_scalar const* weights;
    unsigned int weight_count;
    int is_closed;
} qaws_nurbs_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | 1 to 31 |
| `control_points` | `dimension` scalars per point, packed (not pre-multiplied by the weights) |
| `control_point_count` | at least `degree + 1` |
| `knots` | required, non-decreasing |
| `knot_count` | must be `control_point_count + degree + 1` |
| `weights` | required, one per control point, each > 0 |
| `weight_count` | must equal `control_point_count` |
| `is_closed` | not read by creation; whether the curve is closed is reported from its geometry (`qaws_curve_is_closed`, [inspect.md](inspect.md)) |

## qaws_curve_create_nurbs, qaws_curve_create_nurbs_ex

```c
qaws_status qaws_curve_create_nurbs(
    qaws_curve** out_curve,
    qaws_nurbs_desc const* desc);

qaws_status qaws_curve_create_nurbs_ex(
    qaws_curve** out_curve,
    qaws_nurbs_desc const* desc,
    qaws_allocator const* allocator);
```

Creates the curve. The `_ex` form takes its memory from `allocator` (NULL: the default, see [alloc.md](alloc.md)).

- Parameter range [knots[degree], knots[control_point_count]]; each non-empty knot interval in it is one span.
- `qaws_curve_is_rational` reports 1. With every weight 1 the curve is the B-spline of the same points and knots; a
  larger weight pulls the curve toward its point.
- Inputs are copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve`, `control_points`, `knots` or `weights` NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree 0 or above 31 |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | fewer than `degree + 1` points, or a count too large to index |
| `QAWS_STATUS_INVALID_KNOT_COUNT` | `knot_count` not `control_point_count + degree + 1` |
| `QAWS_STATUS_INVALID_KNOT_VECTOR` | a knot smaller than the one before it, or no non-empty span |
| `QAWS_STATUS_INVALID_WEIGHT_COUNT` | `weight_count != control_point_count`, or a weight <= 0 |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first.

## Example: quarter circle

```c
qaws_vec2 pts[] = { {1, 0}, {1, 1}, {0, 1} };
qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };
qaws_scalar weights[] = { 1, (qaws_scalar)0.70710678118654752, 1 };   /* 1/sqrt(2) */
qaws_nurbs_desc d = { QAWS_DIMENSION_2D, 2, pts, 3, knots, 6, weights, 3, 0 };
qaws_curve* c;
qaws_status st = qaws_curve_create_nurbs(&c, &d);   /* an exact 90-degree arc of the unit circle */
```

A full circle can be built from four such arcs; `tests/25_svg_curves.c` has one.
