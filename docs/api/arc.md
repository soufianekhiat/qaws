# API Reference: qaws_arc.h

Circular and elliptical arc creation.

---

## qaws_arc_segment

```c
typedef struct qaws_arc_segment {
    qaws_scalar center[3];
    qaws_scalar radius;
    qaws_scalar angle_start;
    qaws_scalar angle_end;
    qaws_scalar axis_u[3];
    qaws_scalar axis_v[3];
} qaws_arc_segment;
```

| Field | Description |
|---|---|
| `center` | Arc center (x, y, z). For 2D, only x and y are used. |
| `radius` | Arc radius. Must be > 0. |
| `angle_start` | Start angle in radians. For 2D, CCW from +x axis. |
| `angle_end` | End angle in radians. |
| `axis_u` | 3D only: unit vector for the cosine direction (ignored for 2D). |
| `axis_v` | 3D only: unit vector for the sine direction (ignored for 2D). |

---

## qaws_arc_desc

```c
typedef struct qaws_arc_desc {
    qaws_dimension dimension;
    qaws_arc_segment const* segments;
    unsigned int segment_count;
} qaws_arc_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `segments` | Array of arc segments |
| `segment_count` | Number of segments (>= 1) |

---

## qaws_curve_create_arc

```c
qaws_status qaws_curve_create_arc(
    qaws_curve** out_curve,
    qaws_arc_desc const* desc);
```

Creates an arc curve from one or more arc segments.

**Parameters:**
- `desc` -- Arc descriptor. Must not be NULL.
- `out_curve` -- Receives the created curve. Set to NULL on error.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- For 2D arcs, angles are in radians, measured counter-clockwise from the +x axis.
- For 3D arcs, the arc lies in a plane defined by the `axis_u` and `axis_v` vectors.
- Degenerate zero-radius segments are not supported.
- Multiple segments are chained end-to-end.
