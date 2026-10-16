# API Reference: qaws_bezier.h

Bezier curve creation.

---

## qaws_bezier_desc

```c
typedef struct qaws_bezier_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* control_points;
    unsigned int control_point_count;
} qaws_bezier_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | Polynomial degree (1=linear, 2=quadratic, 3=cubic, ...) |
| `control_points` | Pointer to array of `qaws_vec2` or `qaws_vec3` |
| `control_point_count` | Number of control points. Must equal `degree + 1`. |

---

## qaws_curve_create_bezier

```c
qaws_status qaws_curve_create_bezier(
    qaws_curve** out_curve,
    qaws_bezier_desc const* desc);
```

Creates a single-span Bezier curve.

**Parameters:**
- `desc` -- Curve descriptor. Must not be NULL.
- `out_curve` -- Receives the created curve. Set to NULL on error.

**Returns:** `QAWS_STATUS_OK` on success.

**Errors:**
- `INVALID_ARGUMENT` -- NULL desc or out_curve, or NULL control_points
- `INVALID_DIMENSION` -- dimension is not 2D or 3D
- `INVALID_DEGREE` -- degree is 0
- `INVALID_CONTROL_POINT_COUNT` -- count != degree + 1
- `ALLOCATION_FAILURE` -- malloc failed

**Notes:**
- The descriptor data is copied. The caller may free `control_points` after creation.
- Parameter range is [0, 1].
- Span count is always 1.
