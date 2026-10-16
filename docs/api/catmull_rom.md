# API Reference: qaws_catmull_rom.h

Catmull-Rom spline creation.

---

## qaws_catmull_rom_desc

```c
typedef struct qaws_catmull_rom_desc {
    qaws_dimension dimension;
    void const* control_points;
    unsigned int control_point_count;
    qaws_parameterization parameterization;
    int closed;
} qaws_catmull_rom_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `control_points` | Array of `qaws_vec2` or `qaws_vec3` |
| `control_point_count` | Number of control points. Must be >= 4 for open curves. |
| `parameterization` | `UNIFORM`, `CHORDAL`, or `CENTRIPETAL` |
| `closed` | 0 for open, 1 for closed |

---

## qaws_curve_create_catmull_rom

```c
qaws_status qaws_curve_create_catmull_rom(
    qaws_curve** out_curve,
    qaws_catmull_rom_desc const* desc);
```

Creates a Catmull-Rom interpolating spline.

**Parameters:**
- `desc` -- Curve descriptor
- `out_curve` -- Receives the created curve

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- For an open curve, the first and last points are ghost points that define tangent directions. The curve interpolates the interior points (index 1 to N-2).
- For a closed curve, all points are interpolated.
- `CENTRIPETAL` parameterization is recommended for general use.
- Degree is always 3 (cubic).
- Continuity: C1 at span boundaries.
