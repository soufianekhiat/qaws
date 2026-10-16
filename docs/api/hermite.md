# API Reference: qaws_hermite.h

Hermite curve creation.

---

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

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | Must be 3 (cubic Hermite) |
| `points` | Array of interpolation points (`qaws_vec2*` or `qaws_vec3*`) |
| `derivatives` | Array of tangent vectors at each point (same type as points) |
| `point_count` | Number of interpolation points. Must be >= 2. |
| `derivative_count` | Must equal `point_count` |

---

## qaws_curve_create_hermite

```c
qaws_status qaws_curve_create_hermite(
    qaws_curve** out_curve,
    qaws_hermite_desc const* desc);
```

Creates a piecewise cubic Hermite curve.

**Parameters:**
- `desc` -- Curve descriptor
- `out_curve` -- Receives the created curve

**Returns:** `QAWS_STATUS_OK` on success.

**Errors:**
- `INVALID_ARGUMENT` -- NULL pointers
- `INVALID_DIMENSION` -- not 2D or 3D
- `INVALID_DEGREE` -- degree != 3
- `INVALID_CONTROL_POINT_COUNT` -- point_count < 2, or derivative_count != point_count
- `ALLOCATION_FAILURE` -- malloc failed

**Notes:**
- The curve interpolates all points exactly.
- Tangents at each point are as specified in `derivatives`.
- Parameter range is [0, point_count-1]. Each span maps to a unit interval.
- Continuity: C1 at span boundaries.
