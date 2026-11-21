# qaws_rational_bezier.h: rational Bezier curves

A single-span rational Bezier curve of degree 1 to 47; with suitable weights it represents conics (circles, ellipses,
parabolas, hyperbolas) exactly. Each name exists as `name_f32` and `name_f64` ([precision.md](precision.md)).

## qaws_rational_bezier_desc

```c
typedef struct qaws_rational_bezier_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* control_points;       /* dim_count * (degree+1) scalars */
    unsigned int control_point_count; /* must be degree+1 */
    qaws_scalar const* weights;       /* degree+1 positive weights */
    unsigned int weight_count;        /* must be degree+1 */
} qaws_rational_bezier_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | 1 to 47 |
| `control_points` | `degree + 1` points, `dimension` scalars each, packed |
| `control_point_count` | must be `degree + 1` |
| `weights` | `degree + 1` weights, each > 0 |
| `weight_count` | must be `degree + 1` |

## qaws_curve_create_rational_bezier

```c
qaws_status qaws_curve_create_rational_bezier(
    qaws_curve** out_curve,
    qaws_rational_bezier_desc const* desc);
```

Creates the curve. Parameter range [0, 1], one span. With every weight 1 it is the polynomial Bezier of the same
points. Points and weights are copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve`, `control_points` or `weights` NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree 0 or above 47 |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | `control_point_count != degree + 1` |
| `QAWS_STATUS_INVALID_WEIGHT_COUNT` | `weight_count != degree + 1`, or a weight <= 0 |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_rational_bezier_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
