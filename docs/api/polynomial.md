# qaws_polynomial.h: polynomial curves

A single-span polynomial curve in monomial form, `P(t) = c0 + c1*t + c2*t^2 + ... + cn*t^n`, over a chosen parameter
range. Each name exists as `name_f32` and `name_f64` ([precision.md](precision.md)). For a stack-held polynomial curve
with no allocation, see [inline.md](inline.md).

## qaws_polynomial_desc

```c
typedef struct qaws_polynomial_desc {
    qaws_dimension dimension;
    unsigned int degree;                /* polynomial degree (>= 1) */
    void const* coefficients;           /* (degree+1) * dim_count scalars */
    unsigned int coefficient_count;     /* must be degree+1 */
    qaws_scalar t_min;                  /* parameter range min (default 0) */
    qaws_scalar t_max;                  /* parameter range max (default 1) */
} qaws_polynomial_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | at least 1 |
| `coefficients` | each coefficient a point: `[c0_x, c0_y, (c0_z), c1_x, c1_y, (c1_z), ...]` |
| `coefficient_count` | must be `degree + 1` |
| `t_min`, `t_max` | the parameter range; both 0 means the default [0, 1]. `t` in `P(t)` is the curve parameter itself, not rescaled to [0, 1] |

## qaws_curve_create_polynomial

```c
qaws_status qaws_curve_create_polynomial(
    qaws_curve** out_curve,
    qaws_polynomial_desc const* desc);
```

Creates the curve with parameter range [t_min, t_max], one span. Coefficients are copied. Destroy the result with
`qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `coefficients` NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree 0, or too large to index |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | `coefficient_count != degree + 1` |
| `QAWS_STATUS_INVALID_PARAMETER_RANGE` | `t_max <= t_min` (after the default) |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_polynomial_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
