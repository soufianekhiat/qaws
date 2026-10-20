# API Reference: qaws_subdivision.h

Subdivision curve creation.

---

## qaws_subdivision_scheme

```c
typedef enum qaws_subdivision_scheme {
    QAWS_SUBDIVISION_CHAIKIN = 0,
    QAWS_SUBDIVISION_LANE_RIESENFELD_3,
    QAWS_SUBDIVISION_LANE_RIESENFELD_4
} qaws_subdivision_scheme;
```

| Value | Description |
|---|---|
| `QAWS_SUBDIVISION_CHAIKIN` | Chaikin's algorithm. Quadratic B-spline limit (C1). |
| `QAWS_SUBDIVISION_LANE_RIESENFELD_3` | Lane-Riesenfeld degree 3. Cubic B-spline limit (C2). |
| `QAWS_SUBDIVISION_LANE_RIESENFELD_4` | Lane-Riesenfeld degree 4. Quartic B-spline limit (C3). |

---

## qaws_subdivision_desc

```c
typedef struct qaws_subdivision_desc {
    qaws_dimension dimension;
    qaws_subdivision_scheme scheme;
    void const* control_points;
    unsigned int control_point_count;
    int closed;
    unsigned int refinement_levels;
} qaws_subdivision_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `scheme` | Subdivision scheme (see `qaws_subdivision_scheme`) |
| `control_points` | Initial control polygon (`qaws_vec2*` or `qaws_vec3*`) |
| `control_point_count` | Number of control points (>= 3) |
| `closed` | 1 for a closed polygon, 0 for open |
| `refinement_levels` | Number of subdivision iterations (default 6) |

---

## qaws_curve_create_subdivision

```c
qaws_status qaws_curve_create_subdivision(
    qaws_curve** out_curve,
    qaws_subdivision_desc const* desc);
```

Creates a subdivision curve by iteratively refining a control polygon.

**Parameters:**
- `desc` -- Subdivision descriptor. Must not be NULL.
- `out_curve` -- Receives the created curve. Set to NULL on error.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The descriptor data is copied. The caller may free `control_points` after creation.
- Higher `refinement_levels` produce smoother curves at the cost of more internal points.
- The limit curve's continuity depends on the chosen scheme (C1, C2, or C3).
