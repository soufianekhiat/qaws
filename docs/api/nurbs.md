# API Reference: qaws_nurbs.h

NURBS (Non-Uniform Rational B-Spline) curve creation.

---

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

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | Polynomial degree (commonly 2 for conics, 3 for general) |
| `control_points` | Array of `qaws_vec2` or `qaws_vec3` |
| `control_point_count` | Number of control points |
| `knots` | Knot vector. Must not be NULL. |
| `knot_count` | Must equal control_point_count + degree + 1 |
| `weights` | Per-control-point weights. Must not be NULL. |
| `weight_count` | Must equal control_point_count |
| `is_closed` | 0 for open, 1 for closed |

---

## qaws_curve_create_nurbs

```c
qaws_status qaws_curve_create_nurbs(
    qaws_curve** out_curve,
    qaws_nurbs_desc const* desc);
```

Creates a NURBS curve.

**Parameters:**
- `desc` -- Curve descriptor
- `out_curve` -- Receives the created curve

**Returns:** `QAWS_STATUS_OK` on success.

**Errors:**
- `INVALID_KNOT_COUNT` -- knot_count mismatch
- `INVALID_KNOT_VECTOR` -- not non-decreasing
- `INVALID_WEIGHT_COUNT` -- weight_count != control_point_count

**Notes:**
- `qaws_curve_is_rational()` returns 1 for NURBS curves.
- Weight of 1.0 means the control point has normal influence. Higher weights pull the curve closer.
- NURBS with all weights equal to 1.0 reduces to a standard B-Spline.

**Example: Quarter circle (degree 2)**

```c
qaws_vec2 pts[] = { {1,0}, {1,1}, {0,1} };
qaws_scalar knots[] = { 0,0,0, 1,1,1 };
qaws_scalar weights[] = { 1, 0.70710678118f, 1 };  // 1/sqrt(2)
```

This produces an exact 90-degree circular arc.

**Example: Full circle (4 quarter arcs)**

Use 4 NURBS quarter-arcs with appropriate control points and weights. See `tests/test_runner.c` `test_svg_circle` for a complete example.
