# Traversal

A traversal wraps a curve (or a chain of curves) with a motion profile, so
you evaluate by distance or time instead of raw parameter. Names below are the
unsuffixed forms ([precision.md](precision.md)). For rigid-body and
per-vertex motions over t in [0, 1] (poses, frames along a curve), see
`qaws_motion.h` ([api/motion.md](api/motion.md)).

## Concept

A curve defines *shape*; a traversal defines *how to move along it*. The same
curve can be evaluated by parameter, traversed at constant speed, with an
acceleration profile, or by time with easing.

The traversal does **not** own its curves: they must outlive it. Creation
builds an arc-length table (and, for a custom speed profile, a time table);
evaluation then maps the input to a parameter and evaluates the curve there,
so the returned derivatives are the curve's, with respect to its parameter.

## Creating a traversal

```c
qaws_traversal_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    .motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED,
    .speed = 1.0f,
    .clamp_to_domain = 1
};
qaws_traversal* trav = NULL;
qaws_traversal_create(&trav, curve, &desc);
```

`qaws_traversal_desc` members:

| Member | Meaning |
|---|---|
| `traversal_mode` | How `input_value` is read (below) |
| `motion_profile` | Time to distance mapping in time mode (below) |
| `speed` | Constant speed; initial speed for constant acceleration |
| `acceleration` | Constant acceleration; ramp acceleration for trapezoidal and S-curve (0: 1) |
| `max_speed` | Cruise speed for trapezoidal and S-curve (S-curve 0: 1) |
| `start_time`, `end_time` | Time window of the motion (trapezoidal, S-curve, custom, easing) |
| `clamp_to_domain` | Clamp the resulting parameter to the curve's range |
| `easing` | `qaws_easing`, applied to normalized time in time mode |
| `wrap_mode` | `QAWS_WRAP_CLAMP`, `QAWS_WRAP_LOOP`, `QAWS_WRAP_PING_PONG` on the distance |
| `jerk` | Maximum jerk of the S-curve profile (0: 1) |
| `custom_speed_fn`, `custom_speed_fn_user_data` | Custom profile: speed as a function of distance |

## Traversal modes

| Mode | `input_value` means |
|---|---|
| `QAWS_TRAVERSAL_MODE_PARAMETER` | The curve parameter |
| `QAWS_TRAVERSAL_MODE_ARC_LENGTH` | Distance from the start (wrapped by `wrap_mode`) |
| `QAWS_TRAVERSAL_MODE_TIME` | Time: eased, mapped to distance by the profile, then wrapped |

## Motion profiles (time mode)

| Profile | Distance at time t (dt = t - start_time) |
|---|---|
| `QAWS_MOTION_PROFILE_NONE` | `t` itself |
| `QAWS_MOTION_PROFILE_CONSTANT_SPEED` | `speed * dt` |
| `QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION` | `speed * dt + acceleration * dt^2 / 2` |
| `QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED` | Accelerate to `max_speed`, cruise, decelerate, over `[start_time, end_time]` (ramps shortened when the window is short) |
| `QAWS_MOTION_PROFILE_S_CURVE` | Seven-segment jerk-limited profile from `max_speed`, `acceleration`, `jerk` |
| `QAWS_MOTION_PROFILE_CUSTOM` | Integrates `custom_speed_fn(distance, user)` over `[start_time, end_time]` into a table |

Easing (`QAWS_EASING_LINEAR`, `QUAD_IN/OUT/IN_OUT`, `CUBIC_IN/OUT/IN_OUT`,
`SINE_IN/OUT/IN_OUT`) remaps time normalized to `[start_time, end_time]`
before the profile.

## Evaluating and mapping

```c
qaws_eval_result_2d r;
qaws_traversal_evaluate_2d(&r, trav, 0.5f, QAWS_EVAL_FLAG_POSITION);

qaws_scalar param, distance;
qaws_traversal_map_time_to_parameter(&param, trav, 1.0f);
qaws_traversal_map_distance_to_parameter(&param, trav, 5.0f);
qaws_traversal_map_parameter_to_distance(&distance, trav, 0.5f);
```

`qaws_traversal_evaluate_3d` is the 3D form.

## Multi-curve traversal

```c
qaws_status qaws_traversal_create_multi(
    qaws_traversal** out_traversal,
    qaws_curve const* const* curves,
    unsigned int curve_count,
    qaws_traversal_desc const* desc);
```

The curves are treated as one path; the input (parameter, distance or time)
is mapped across their concatenated domain and resolved to the curve that
holds it. None of the curves is owned.

## Advance and reset

`qaws_traversal_advance_2d/3d(&r, trav, delta_time, flags)` steps the
traversal's own time forward by `delta_time` and evaluates there;
`qaws_traversal_reset(trav)` returns to the start. Both mutate the traversal:
synchronize them when several threads share one. Evaluate and map calls are
safe for concurrent reads ([threading.md](threading.md)).

## Lifecycle

```c
qaws_traversal_destroy(trav);
qaws_curve_destroy(curve);   /* the curve is destroyed separately */
```

## Example: constant-speed motion along a Bezier

```c
qaws_vec2 pts[] = { {0,0}, {0,2}, {2,2}, {2,0} };
qaws_bezier_desc bdesc = { .dimension = QAWS_DIMENSION_2D, .degree = 3,
                           .control_points = pts, .control_point_count = 4 };
qaws_curve* curve = NULL;
qaws_curve_create_bezier(&curve, &bdesc);

qaws_traversal_desc tdesc = { .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
                              .motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED,
                              .speed = 1.0f, .clamp_to_domain = 1 };
qaws_traversal* trav = NULL;
qaws_traversal_create(&trav, curve, &tdesc);

qaws_scalar total_len;
qaws_curve_compute_arc_length(&total_len, curve, 0, 1);
for (int i = 1; i <= 3; i++) {
    qaws_eval_result_2d r;
    qaws_traversal_evaluate_2d(&r, trav, total_len * (qaws_scalar)i / 4, QAWS_EVAL_FLAG_POSITION);
    /* r.position at 25%, 50%, 75% of the length */
}

qaws_traversal_destroy(trav);
qaws_curve_destroy(curve);
```
