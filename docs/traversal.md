# Traversal

Traversal wraps a curve with a motion profile, allowing you to evaluate by distance or time instead of raw parameter.

## Concept

A curve defines *shape*. A traversal defines *how to move along that shape*.

The same curve can be:
- evaluated directly by parameter
- traversed at constant speed
- traversed with acceleration profiles
- traversed by time

## Creating a traversal

```c
qaws_traversal_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    .motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED,
    .speed = 1.0f,
    .acceleration = 0,
    .max_speed = 0,
    .start_time = 0,
    .end_time = 0,
    .clamp_to_domain = 1
};

qaws_traversal* trav = NULL;
qaws_traversal_create(&trav, curve, &desc);
```

The traversal does **not** own the curve. The curve must remain alive for the traversal's lifetime.

## Traversal modes

| Mode | Input meaning |
|---|---|
| `QAWS_TRAVERSAL_MODE_PARAMETER` | Raw curve parameter |
| `QAWS_TRAVERSAL_MODE_ARC_LENGTH` | Distance along curve from start |
| `QAWS_TRAVERSAL_MODE_TIME` | Elapsed time (uses motion profile) |

## Motion profiles

| Profile | Description |
|---|---|
| `NONE` | No motion mapping, pass-through |
| `CONSTANT_SPEED` | Uniform speed along the curve |
| `CONSTANT_ACCELERATION` | Linearly increasing speed |
| `TRAPEZOIDAL_SPEED` | Accelerate, cruise, decelerate |

## Evaluating a traversal

```c
qaws_eval_result_2d result;
qaws_traversal_evaluate_2d(&result, trav,
    0.5f, QAWS_EVAL_FLAG_POSITION);
```

The `input_value` is interpreted according to the traversal mode:
- Arc-length mode: `0.5` means 0.5 units of distance along the curve
- Time mode: `0.5` means 0.5 seconds elapsed

## Parameter mapping

Convert between parameter spaces:

```c
qaws_scalar param;

// Time -> curve parameter
qaws_traversal_map_time_to_parameter(&param, trav, 1.0f);

// Distance -> curve parameter
qaws_traversal_map_distance_to_parameter(&param, trav, 5.0f);

// Curve parameter -> distance
qaws_scalar distance;
qaws_traversal_map_parameter_to_distance(&distance, trav, 0.5f);
```

## Lifecycle

```c
qaws_traversal_destroy(trav);
// curve is still valid -- destroy it separately
qaws_curve_destroy(curve);
```

## Multi-curve traversal

`qaws_traversal_create_multi` creates a traversal that spans multiple concatenated curves, treating them as a single continuous path. The input value (distance, time, or parameter) is mapped across the full concatenated domain.

```c
qaws_curve* curves[3] = { curve_a, curve_b, curve_c };

qaws_traversal_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    .motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED,
    .speed = 1.0f,
    .clamp_to_domain = 1
};

qaws_traversal* trav = NULL;
qaws_traversal_create_multi(&trav, curves, 3, &desc);

// Evaluate at a distance that may span across curve boundaries
qaws_eval_result_2d r;
qaws_traversal_evaluate_2d(&r, trav,
    5.0f, QAWS_EVAL_FLAG_POSITION);
```

Signature:

```c
qaws_status qaws_traversal_create_multi(
    qaws_traversal** out_traversal,
    qaws_curve const* const* curves,
    unsigned int curve_count,
    qaws_traversal_desc const* desc);
```

None of the input curves are owned by the traversal; they must all remain alive for the traversal's lifetime.

## Advance and reset

For incremental time-stepping (e.g., game loops or animations), use `qaws_traversal_advance_2d`/`3d` to step forward by a delta time from the current internal position. Use `qaws_traversal_reset` to return to the start.

These functions mutate the traversal's internal state and require external synchronization when used from multiple threads.

```c
qaws_traversal* trav = NULL;
qaws_traversal_create(&trav, curve, &desc);

// Game loop: advance by frame delta each tick
qaws_eval_result_2d r;
qaws_traversal_advance_2d(&r, trav,
    0.016f, QAWS_EVAL_FLAG_POSITION);
// r.position is the new position after advancing 16ms

// Reset to the beginning of the curve
qaws_traversal_reset(trav);
```

Signatures:

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

## Example: constant-speed motion along a Bezier

```c
qaws_vec2 pts[] = { {0,0}, {0,2}, {2,2}, {2,0} };
qaws_bezier_desc bdesc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = pts,
    .control_point_count = 4
};
qaws_curve* curve = NULL;
qaws_curve_create_bezier(&curve, &bdesc);

qaws_traversal_desc tdesc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    .motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED,
    .speed = 1.0f,
    .clamp_to_domain = 1
};
qaws_traversal* trav = NULL;
qaws_traversal_create(&trav, curve, &tdesc);

// Evaluate at 25%, 50%, 75% of total arc length
qaws_scalar total_len;
qaws_curve_compute_arc_length(&total_len, curve, 0, 1);

for (int i = 1; i <= 3; i++) {
    qaws_eval_result_2d r;
    qaws_traversal_evaluate_2d(&r,
        trav,
        total_len * (float)i / 4.0f, QAWS_EVAL_FLAG_POSITION);
    // r.position is evenly spaced along the curve
}

qaws_traversal_destroy(trav);
qaws_curve_destroy(curve);
```
