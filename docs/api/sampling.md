# API Reference: qaws_sampling.h

Curve sampling functions.

---

## qaws_curve_get_sample_count

```c
qaws_status qaws_curve_get_sample_count(
    unsigned int* out_sample_count,
    qaws_curve const* curve,
    qaws_sampling_desc const* desc);
```

Queries how many samples will be generated for a given sampling configuration. Use this to allocate a buffer before calling `sample_2d`/`sample_3d`.

---

## qaws_curve_sample_2d

```c
qaws_status qaws_curve_sample_2d(
    qaws_vec2* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int position_capacity,
    qaws_sampling_desc const* desc);
```

Samples a 2D curve into a caller-provided buffer of positions.

**Parameters:**
- `curve` -- The curve. Must be 2D.
- `desc` -- Sampling configuration.
- `out_positions` -- Buffer to receive sampled points.
- `position_capacity` -- Size of the buffer.
- `out_position_count` -- Receives the number of points written.

**Returns:** `QAWS_STATUS_OK` on success. `BUFFER_TOO_SMALL` if capacity is insufficient.

---

## qaws_curve_sample_3d

```c
qaws_status qaws_curve_sample_3d(
    qaws_vec3* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int position_capacity,
    qaws_sampling_desc const* desc);
```

Same as `sample_2d` but for 3D curves.

---

## qaws_curve_sample_eval_2d

```c
qaws_status qaws_curve_sample_eval_2d(
    qaws_eval_result_2d* out_samples,
    unsigned int* out_sample_count,
    qaws_curve const* curve,
    unsigned int sample_capacity,
    qaws_sampling_desc const* desc);
```

Samples a 2D curve into a buffer of full evaluation results (position + derivatives).

---

## qaws_curve_sample_eval_3d

```c
qaws_status qaws_curve_sample_eval_3d(
    qaws_eval_result_3d* out_samples,
    unsigned int* out_sample_count,
    qaws_curve const* curve,
    unsigned int sample_capacity,
    qaws_sampling_desc const* desc);
```

Same as `sample_eval_2d` but for 3D curves.

---

## qaws_curve_sample_curvature_weighted_2d

```c
qaws_status qaws_curve_sample_curvature_weighted_2d(
    qaws_vec2* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int position_capacity);
```

Distributes samples proportionally to curvature magnitude so high-curvature regions receive more points.

**Parameters:**
- `curve` -- The curve. Must be 2D.
- `sample_count` -- Desired number of samples.
- `out_positions` -- Buffer to receive sampled points.
- `position_capacity` -- Size of the buffer.
- `out_position_count` -- Receives the number of points written.

---

## qaws_curve_sample_curvature_weighted_3d

```c
qaws_status qaws_curve_sample_curvature_weighted_3d(
    qaws_vec3* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int position_capacity);
```

Same as `curvature_weighted_2d` but for 3D curves.

---

## qaws_curve_sample_feature_preserving_2d

```c
qaws_status qaws_curve_sample_feature_preserving_2d(
    qaws_vec2* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int base_sample_count,
    unsigned int position_capacity);
```

Inserts inflection points and extrema into a uniform parameter sample set so that important geometric features are never missed.

**Parameters:**
- `curve` -- The curve. Must be 2D.
- `base_sample_count` -- Number of uniform base samples before feature insertion.
- `out_positions` -- Buffer to receive sampled points.
- `position_capacity` -- Size of the buffer (should be larger than `base_sample_count` to accommodate inserted features).
- `out_position_count` -- Receives the actual number of points written.

---

## qaws_curve_sample_feature_preserving_3d

```c
qaws_status qaws_curve_sample_feature_preserving_3d(
    qaws_vec3* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int base_sample_count,
    unsigned int position_capacity);
```

Same as `feature_preserving_2d` but for 3D curves.

---

## Streaming / callback-based sampling

### qaws_sample_callback_2d / 3d

```c
typedef int (*qaws_sample_callback_2d)(
    qaws_scalar parameter,
    qaws_vec2 const* position,
    void* user_data);

typedef int (*qaws_sample_callback_3d)(
    qaws_scalar parameter,
    qaws_vec3 const* position,
    void* user_data);
```

Callback function types for streaming sampling. Return non-zero to continue, or 0 to stop early (not an error).

### qaws_curve_sample_streaming_2d

```c
qaws_status qaws_curve_sample_streaming_2d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_sample_callback_2d callback,
    void* user_data);
```

Samples a 2D curve and delivers each point through a callback instead of writing to a buffer. Useful for processing points one at a time without pre-allocating storage.

**Parameters:**
- `curve` -- The curve. Must be 2D.
- `desc` -- Sampling configuration.
- `callback` -- Called for each sample point. Return 0 to stop early.
- `user_data` -- Opaque pointer passed through to the callback.

---

### qaws_curve_sample_streaming_3d

```c
qaws_status qaws_curve_sample_streaming_3d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_sample_callback_3d callback,
    void* user_data);
```

Same as `streaming_2d` but for 3D curves.

---

## Adaptive tessellation

### qaws_curve_tessellate_2d / 3d

```c
typedef struct qaws_curve_tessellation_desc
{
    qaws_scalar chord_tolerance;        /* distance; 0: unused */
    qaws_scalar angle_tolerance;        /* radians; 0: unused (one of the two is needed) */
    unsigned int max_depth;             /* halvings per span at most; 0: 16, at most 30 */
} qaws_curve_tessellation_desc;

qaws_status qaws_curve_tessellate_2d(
    qaws_scalar* out_parameters,
    qaws_vec2* out_points,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int capacity,
    qaws_curve_tessellation_desc const* desc);

qaws_status qaws_curve_tessellate_3d(/* the same, qaws_vec3* out_points */);
```

A polyline that follows the curve within the given tolerances.

- **Subdivision.** Each span is halved until every segment passes both tests. The chord test measures the distance from the curve to the segment at 1/4, 1/2 and 3/4 of it. The angle test sums the tangent's turning through those points. The three test points catch S shapes, whose middle lies on the chord. Span boundaries are always vertices, and a kink there does not count as turning.
- **Determinism.** The vertices come out in parameter order and are the same on every call. Call once with `out_parameters` NULL for the count (`out_points` may be NULL), then again with buffers of that size.
- **Errors.** `QAWS_STATUS_BUFFER_TOO_SMALL` when `capacity` is below the count; `*out_count` then holds the count, and the buffers hold the first `capacity` vertices. Nothing is allocated.
- **Limits.** The tests sample three points, so they are estimates, not bounds (test 150: the true chord error stays under the tolerance on a circle and an S). `max_depth` limits the halvings near cusps, where the tangent reverses.
