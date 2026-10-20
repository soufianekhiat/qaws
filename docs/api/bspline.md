# API Reference: qaws_bspline.h

B-Spline curve creation.

---

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

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | Polynomial degree (commonly 3) |
| `control_points` | Array of `qaws_vec2` or `qaws_vec3` |
| `control_point_count` | Number of control points |
| `knots` | Custom knot vector, or NULL for auto-generated |
| `knot_count` | Number of knots (must be control_point_count + degree + 1 if provided) |
| `is_uniform` | 1 for uniform clamped knots (auto-generates if knots is NULL) |
| `is_closed` | 0 for open, 1 for closed |

---

## qaws_curve_create_bspline

```c
qaws_status qaws_curve_create_bspline(
    qaws_curve** out_curve,
    qaws_bspline_desc const* desc);
```

Creates a B-Spline curve.

**Parameters:**
- `desc` -- Curve descriptor
- `out_curve` -- Receives the created curve

**Returns:** `QAWS_STATUS_OK` on success.

**Errors:**
- `INVALID_KNOT_COUNT` -- knot_count != control_point_count + degree + 1
- `INVALID_KNOT_VECTOR` -- knots are not non-decreasing

**Notes:**
- When `is_uniform = 1` and `knots = NULL`, a clamped uniform knot vector is generated automatically. Clamped means the first and last knots have multiplicity degree+1, causing the curve to interpolate the first and last control points.
- When `knots` is provided, it is used directly.
- Continuity at internal knots: up to C(degree-1) for simple knots.
- The parameter range spans the active portion of the knot vector.

**Knot vector reference:**

For a degree-3 B-Spline with 6 control points, a clamped uniform knot vector is:
```
{0, 0, 0, 0, 1, 2, 3, 3, 3, 3}
```
The active parameter range is [0, 3], with 3 spans.
