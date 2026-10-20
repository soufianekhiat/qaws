# API Reference: qaws_polynomial.h

Polynomial curve creation in monomial form.

---

## qaws_polynomial_desc

```c
typedef struct qaws_polynomial_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* coefficients;
    unsigned int coefficient_count;
    qaws_scalar t_min;
    qaws_scalar t_max;
} qaws_polynomial_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | Polynomial degree (>= 1) |
| `coefficients` | `(degree + 1) * dim_count` scalars. Layout: `[c0_x, c0_y, [c0_z], c1_x, c1_y, [c1_z], ...]` |
| `coefficient_count` | Must equal `degree + 1` |
| `t_min` | Parameter range minimum (default 0) |
| `t_max` | Parameter range maximum (default 1) |

The curve is defined as: `P(t) = c0 + c1*t + c2*t^2 + ... + cn*t^n`

Each coefficient is a point (2 or 3 scalars depending on dimension).

---

## qaws_curve_create_polynomial

```c
qaws_status qaws_curve_create_polynomial(
    qaws_curve** out_curve,
    qaws_polynomial_desc const* desc);
```

Creates a polynomial curve in monomial form.

**Parameters:**
- `desc` -- Curve descriptor. Must not be NULL.
- `out_curve` -- Receives the created curve. Set to NULL on error.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The descriptor data is copied. The caller may free `coefficients` after creation.
- Parameter range is [t_min, t_max].
