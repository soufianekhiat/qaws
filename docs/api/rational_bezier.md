# API Reference: qaws_rational_bezier.h

Rational Bezier curve creation.

---

## qaws_rational_bezier_desc

```c
typedef struct qaws_rational_bezier_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* control_points;
    unsigned int control_point_count;
    qaws_scalar const* weights;
    unsigned int weight_count;
} qaws_rational_bezier_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | Polynomial degree (>= 1) |
| `control_points` | Array of `qaws_vec2` or `qaws_vec3` (`degree + 1` points) |
| `control_point_count` | Must equal `degree + 1` |
| `weights` | Array of `degree + 1` positive weights |
| `weight_count` | Must equal `degree + 1` |

---

## qaws_curve_create_rational_bezier

```c
qaws_status qaws_curve_create_rational_bezier(
    qaws_curve** out_curve,
    qaws_rational_bezier_desc const* desc);
```

Creates a rational Bezier curve.

**Parameters:**
- `desc` -- Curve descriptor. Must not be NULL.
- `out_curve` -- Receives the created curve. Set to NULL on error.

**Returns:** `QAWS_STATUS_OK` on success.

**Errors:**
- `INVALID_ARGUMENT` -- NULL desc or out_curve, or NULL control_points/weights
- `INVALID_DIMENSION` -- dimension is not 2D or 3D
- `INVALID_DEGREE` -- degree is 0
- `INVALID_CONTROL_POINT_COUNT` -- control_point_count != degree + 1
- `INVALID_WEIGHT_COUNT` -- weight_count != degree + 1
- `ALLOCATION_FAILURE` -- malloc failed

**Notes:**
- The descriptor data is copied. The caller may free `control_points` and `weights` after creation.
- Parameter range is [0, 1].
- Rational Bezier curves can exactly represent conic sections (circles, ellipses, parabolas, hyperbolas) by choosing appropriate weights.
- With all weights equal to 1.0, the curve reduces to a standard Bezier.
