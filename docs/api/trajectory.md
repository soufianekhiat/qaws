# API Reference: qaws_trajectory.h

Trajectory curve creation.

---

## qaws_trajectory_desc

```c
typedef struct qaws_trajectory_desc {
    qaws_dimension dimension;
    unsigned int degree;
    void const* key_positions;
    unsigned int key_count;
    qaws_scalar const* key_times;
    unsigned int key_time_count;
    void const* key_velocities;
    unsigned int key_velocity_count;
    void const* key_accelerations;
    unsigned int key_acceleration_count;
    int closed;
} qaws_trajectory_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | Polynomial degree (typically 3) |
| `key_positions` | Array of keyframe positions (`qaws_vec2*` or `qaws_vec3*`) |
| `key_count` | Number of keyframes. Must be >= 2. |
| `key_times` | Time value for each keyframe |
| `key_time_count` | Must equal `key_count` |
| `key_velocities` | Optional velocity constraints at keys (NULL if unused) |
| `key_velocity_count` | Number of velocity constraints (0 if unused) |
| `key_accelerations` | Optional acceleration constraints at keys (NULL if unused) |
| `key_acceleration_count` | Number of acceleration constraints (0 if unused) |
| `closed` | 0 for open, 1 for closed |

---

## qaws_curve_create_trajectory

```c
qaws_status qaws_curve_create_trajectory(
    qaws_curve** out_curve,
    qaws_trajectory_desc const* desc);
```

Creates a time-keyed trajectory curve.

**Parameters:**
- `desc` -- Trajectory descriptor
- `out_curve` -- Receives the created curve

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The parameter domain matches the time domain (from first key_time to last key_time).
- Without velocity/acceleration constraints, the curve uses natural interpolation.
- Works naturally with the traversal system in `QAWS_TRAVERSAL_MODE_TIME`.
- Continuity: C1 at keyframe boundaries.
