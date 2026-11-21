# API Reference: qaws_traversal.h

A traversal walks one curve, or several curves joined end to end, by parameter, arc length or time, with a motion
profile, easing and wrap mode. It holds an arc-length table built at creation.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Threading: the evaluate and map functions are thread-safe for concurrent reads. `qaws_traversal_advance_*` and
`qaws_traversal_reset` mutate state and need external synchronization across threads.

Related: [types.md](types.md) (the enums and `qaws_traversal_desc`), [motion.md](motion.md),
[sampling.md](sampling.md), concepts in [../traversal.md](../traversal.md).

---

## Types (declared in qaws_types.h)

```c
typedef struct qaws_traversal qaws_traversal;   /* opaque */

typedef enum qaws_traversal_mode {
    QAWS_TRAVERSAL_MODE_PARAMETER = 0,   /* input is the curve parameter */
    QAWS_TRAVERSAL_MODE_ARC_LENGTH,      /* input is a distance along the curve */
    QAWS_TRAVERSAL_MODE_TIME             /* input is a time, mapped through the motion profile */
} qaws_traversal_mode;

typedef enum qaws_motion_profile {
    QAWS_MOTION_PROFILE_NONE = 0,               /* distance = time */
    QAWS_MOTION_PROFILE_CONSTANT_SPEED,         /* speed * (t - start_time) */
    QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION,  /* speed * dt + acceleration * dt^2 / 2 */
    QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED,      /* ramp up at acceleration to max_speed, cruise, ramp down */
    QAWS_MOTION_PROFILE_S_CURVE,                /* seven-segment jerk-limited profile (jerk, acceleration, max_speed) */
    QAWS_MOTION_PROFILE_CUSTOM                  /* speed from custom_speed_fn(distance), integrated over time */
} qaws_motion_profile;

typedef qaws_scalar (*qaws_speed_fn)(qaws_scalar distance, void* user_data);

typedef enum qaws_easing {
    QAWS_EASING_LINEAR = 0,
    QAWS_EASING_QUAD_IN,  QAWS_EASING_QUAD_OUT,  QAWS_EASING_QUAD_IN_OUT,
    QAWS_EASING_CUBIC_IN, QAWS_EASING_CUBIC_OUT, QAWS_EASING_CUBIC_IN_OUT,
    QAWS_EASING_SINE_IN,  QAWS_EASING_SINE_OUT,  QAWS_EASING_SINE_IN_OUT
} qaws_easing;

typedef enum qaws_wrap_mode {
    QAWS_WRAP_CLAMP = 0,   /* distance clamped to [0, total length] */
    QAWS_WRAP_LOOP,        /* distance modulo total length */
    QAWS_WRAP_PING_PONG    /* back and forth over [0, total length] */
} qaws_wrap_mode;

typedef struct qaws_traversal_desc {
    qaws_traversal_mode traversal_mode;
    qaws_motion_profile motion_profile;
    qaws_scalar speed;
    qaws_scalar acceleration;
    qaws_scalar max_speed;
    qaws_scalar start_time;
    qaws_scalar end_time;
    int clamp_to_domain;
    qaws_easing easing;
    qaws_wrap_mode wrap_mode;
    qaws_scalar jerk;                    /* max jerk for the S-curve profile */
    qaws_speed_fn custom_speed_fn;       /* custom speed callback */
    void* custom_speed_fn_user_data;     /* user data for the custom speed callback */
} qaws_traversal_desc;
```

Notes on the profiles:
- Trapezoidal and S-curve use `start_time` / `end_time` as the move's duration; acceleration, max speed and jerk
  of 0 or less are taken as 1. When the duration is too short to reach `max_speed`, the ramps shrink to fit.
- Custom: `custom_speed_fn` is integrated (RK4) over `[start_time, end_time]` into a table at creation, capped at
  the total arc length. With no callback or a non-positive duration no table is built and the profile maps
  every time to distance 0.
- `clamp_to_domain` clamps the final parameter to the curve's domain.

---

## Lifetime

```c
qaws_status qaws_traversal_create(
    qaws_traversal** out_traversal,
    qaws_curve const* curve,
    qaws_traversal_desc const* desc);

qaws_status qaws_traversal_create_multi(
    qaws_traversal** out_traversal,
    qaws_curve const* const* curves,
    unsigned int curve_count,
    qaws_traversal_desc const* desc);

void qaws_traversal_destroy(qaws_traversal* traversal);
```

- `qaws_traversal_create` wraps one curve. It copies `desc`, builds the arc-length table (and the custom speed
  table when used) and starts at distance 0.
- `qaws_traversal_create_multi` joins `curve_count` curves end to end. All must have the same dimension
  (`QAWS_STATUS_INVALID_DIMENSION` otherwise). The combined parameter runs from 0 to the sum of the curves'
  parameter range lengths; the arc-length table covers the whole chain. The curve array itself is copied.
- The traversal does not own the curves: they must outlive it.
- `*out_traversal` is set to NULL on failure. Returns `QAWS_STATUS_INVALID_ARGUMENT` (null pointer,
  `curve_count == 0`, a null curve), `QAWS_STATUS_ALLOCATION_FAILURE`, or a status from building the table.
- `qaws_traversal_destroy` frees the traversal; NULL is a no-op.

---

## Evaluation

```c
qaws_status qaws_traversal_evaluate_2d(
    qaws_eval_result_2d* out_result,
    qaws_traversal const* traversal,
    qaws_scalar input_value,
    unsigned int eval_flags);

qaws_status qaws_traversal_evaluate_3d(
    qaws_eval_result_3d* out_result,
    qaws_traversal const* traversal,
    qaws_scalar input_value,
    unsigned int eval_flags);
```

Evaluates the curve at `input_value`, read according to `traversal_mode`:
- PARAMETER: used as the parameter (the combined parameter for a multi-curve traversal).
- ARC_LENGTH: a distance, wrapped by `wrap_mode`, then mapped to a parameter.
- TIME: eased (when `easing` is not linear and `end_time > start_time`, over the normalized `[start_time, end_time]`),
  mapped to a distance by the motion profile, wrapped, then mapped to a parameter.

For a multi-curve traversal the combined parameter is resolved to one curve and its local parameter.
`eval_flags` is a `QAWS_EVAL_FLAG_*` bitmask ([eval.md](eval.md)). Returns `QAWS_STATUS_INVALID_ARGUMENT` for a
null pointer, otherwise the status of the curve evaluation (for example `QAWS_STATUS_INVALID_DIMENSION`).

---

## Mapping

```c
qaws_status qaws_traversal_map_time_to_parameter(
    qaws_scalar* out_parameter,
    qaws_traversal const* traversal,
    qaws_scalar time_value);

qaws_status qaws_traversal_map_distance_to_parameter(
    qaws_scalar* out_parameter,
    qaws_traversal const* traversal,
    qaws_scalar distance_value);

qaws_status qaws_traversal_map_parameter_to_distance(
    qaws_scalar* out_distance,
    qaws_traversal const* traversal,
    qaws_scalar parameter);
```

- `map_time_to_parameter`: time -> distance through the motion profile -> parameter. Easing and wrap mode are not
  applied here (only in evaluation).
- `map_distance_to_parameter`: distance -> parameter through the arc-length table.
- `map_parameter_to_distance`: parameter -> distance through the arc-length table.
- The first two clamp to the domain of the traversal's curve when `clamp_to_domain` is set (for a multi-curve
  traversal, the first curve's domain). For a multi-curve traversal the parameters are combined parameters.
- Return `QAWS_STATUS_OK`, or `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer.

---

## Stepping

```c
qaws_status qaws_traversal_advance_2d(
    qaws_eval_result_2d* out_result,
    qaws_traversal* traversal,
    qaws_scalar delta_time,
    unsigned int eval_flags);

qaws_status qaws_traversal_advance_3d(
    qaws_eval_result_3d* out_result,
    qaws_traversal* traversal,
    qaws_scalar delta_time,
    unsigned int eval_flags);

qaws_status qaws_traversal_reset(
    qaws_traversal* traversal);
```

- `advance`: adds `delta_time * desc.speed` to the traversal's current distance, wraps it by `wrap_mode`, and
  evaluates there. The motion profile, easing and traversal mode are not used. Mutates the traversal.
- `reset`: sets the current distance back to 0.
- Return `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer; `advance` otherwise returns the curve evaluation status.
