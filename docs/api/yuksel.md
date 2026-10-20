# API Reference: qaws_yuksel.h

Yuksel C2 interpolating spline creation.

Based on: Cem Yuksel, "A Class of C2 Interpolating Splines", ACM Transactions on Graphics, 2020.

---

## qaws_yuksel_mode

```c
typedef enum qaws_yuksel_mode {
    QAWS_YUKSEL_MODE_BEZIER = 0,
    QAWS_YUKSEL_MODE_CIRCULAR,
    QAWS_YUKSEL_MODE_ELLIPTICAL,
    QAWS_YUKSEL_MODE_HYBRID
} qaws_yuksel_mode;
```

| Mode | Description | Dimension |
|---|---|---|
| `BEZIER` | Quadratic Bezier sub-curves through point triplets | 2D, 3D |
| `CIRCULAR` | Circular arcs through point triplets | 2D only |
| `ELLIPTICAL` | Elliptical arcs with iterative fitting | 2D only |
| `HYBRID` | Circular when arc angle is small, elliptical fallback | 2D only |

---

## qaws_yuksel_desc

```c
typedef struct qaws_yuksel_desc {
    qaws_dimension dimension;
    void const* control_points;
    unsigned int control_point_count;
    qaws_yuksel_mode mode;
    int closed;
} qaws_yuksel_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `control_points` | Array of `qaws_vec2` or `qaws_vec3` |
| `control_point_count` | Number of control points. Must be >= 3. |
| `mode` | Interpolation mode (see above) |
| `closed` | 0 for open, 1 for closed |

---

## qaws_curve_create_yuksel

```c
qaws_status qaws_curve_create_yuksel(
    qaws_curve** out_curve,
    qaws_yuksel_desc const* desc);
```

Creates a Yuksel C2 interpolating spline.

**Parameters:**
- `desc` -- Curve descriptor
- `out_curve` -- Receives the created curve

**Returns:** `QAWS_STATUS_OK` on success.

**Errors:**
- `INVALID_ARGUMENT` -- NULL pointers, or non-Bezier mode with 3D dimension
- `INVALID_CONTROL_POINT_COUNT` -- fewer than 3 points
- `ALLOCATION_FAILURE` -- malloc failed

**Notes:**
- All control points are interpolated exactly.
- Continuity is C2 everywhere (achieved via cos^2/sin^2 blending).
- For open curves: span_count = control_point_count - 1.
- For closed curves: span_count = control_point_count.
- Parameter range is [0, span_count].
- The degree reported by `qaws_curve_get_degree()` is 2 (quadratic sub-curves).
- 3D is only supported with `QAWS_YUKSEL_MODE_BEZIER`. Circular, Elliptical, and Hybrid modes require 2D because they rely on plane-specific geometry (perpendicular bisectors, arc angles).

**Open curve boundary behavior:**
- The first span uses the sub-curve at P[1] unblended (only the left half).
- The last span uses the sub-curve at P[N-2] unblended (only the right half).
- Interior spans blend two adjacent sub-curves.
