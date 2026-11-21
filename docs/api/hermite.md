# qaws_hermite.h: Hermite curves

A piecewise cubic Hermite curve through given points with given derivatives. Each name exists as `name_f32` and
`name_f64` ([precision.md](precision.md)).

## qaws_hermite_desc

```c
typedef struct qaws_hermite_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* points;
    void const* derivatives;
    unsigned int point_count;
    unsigned int derivative_count;
} qaws_hermite_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | must be 3 |
| `points` | the points interpolated, `dimension` scalars each, packed |
| `derivatives` | the derivative (with respect to the curve parameter) at each point, same layout |
| `point_count` | at least 2 |
| `derivative_count` | must equal `point_count` |

## qaws_curve_create_hermite

```c
qaws_status qaws_curve_create_hermite(
    qaws_curve** out_curve,
    qaws_hermite_desc const* desc);
```

Creates the curve: `point_count - 1` cubic spans, span i on [i, i + 1], parameter range [0, point_count - 1]. The curve
passes through every point with the given derivative there, so it is C1 at the points. Inputs are copied. Destroy the
result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve`, `points` or `derivatives` NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree not 3 |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | fewer than 2 points, `derivative_count != point_count`, or a count too large to index |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_hermite_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
