# qaws_trajectory.h: time-keyed trajectories

A piecewise cubic curve through keyed positions whose parameter is time: the curve is at each key's position at that
key's time. Each name exists as `name_f32` and `name_f64` ([precision.md](precision.md)). Pairs with
`QAWS_TRAVERSAL_MODE_TIME` ([traversal.md](traversal.md)).

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

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `degree` | must be 3 (cubic Hermite spans) |
| `key_positions` | `key_count` positions, `dimension` scalars each, packed |
| `key_count` | at least 2 |
| `key_times` | one time per key, strictly increasing |
| `key_time_count` | must equal `key_count` |
| `key_velocities` | optional velocity (d position / d time) at each key, same layout as positions |
| `key_velocity_count` | used only when it equals `key_count` (with `key_velocities` non-NULL); otherwise velocities are estimated from the keys |
| `key_accelerations` | optional; copied and kept with the curve, not used to shape it |
| `key_acceleration_count` | number of accelerations given (0: none) |
| `closed` | non-zero marks the trajectory closed: reported by `qaws_curve_is_closed` and used when velocities are estimated; it adds no span |

## qaws_curve_create_trajectory

```c
qaws_status qaws_curve_create_trajectory(
    qaws_curve** out_curve,
    qaws_trajectory_desc const* desc);
```

Creates the curve: `key_count - 1` spans, span boundaries at the key times, parameter range
[key_times[0], key_times[key_count - 1]]. With velocities given, each span is the cubic Hermite with those end
velocities; without, the velocities are estimated by finite differences of the keys. C1 at the keys. Inputs are
copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve`, `key_positions` or `key_times` NULL; fewer than 2 keys; `key_time_count != key_count`; a count too large to index |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree not 3 |
| `QAWS_STATUS_INVALID_PARAMETER_RANGE` | times not strictly increasing |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_trajectory_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
