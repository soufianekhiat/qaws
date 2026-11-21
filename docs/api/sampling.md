# API Reference: qaws_sampling.h

Turning a curve into points: uniform or arc-length sampling into caller buffers, sampling with derivatives,
curvature-weighted and feature-preserving sampling, streaming through a callback, and adaptive tessellation
to tolerances. All functions are thread-safe on immutable curves.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Related: [eval.md](eval.md), [traversal.md](traversal.md), concepts in [../sampling.md](../sampling.md).

---

## qaws_sampling_desc

Declared in `qaws_types.h` (see [types.md](types.md)); used by the sampling and streaming functions.

```c
typedef struct qaws_sampling_desc
{
    qaws_traversal_mode traversal_mode;
    unsigned int sample_count;
    qaws_scalar error_tolerance;
    int use_adaptive_sampling;
} qaws_sampling_desc;
```

- `traversal_mode`: `QAWS_TRAVERSAL_MODE_ARC_LENGTH` spaces samples evenly by arc length (through an internal
  arc-length table); any other mode spaces them evenly in parameter, from the domain's start to its end.
- `sample_count`: number of samples (the initial uniform count in adaptive mode, at least 2 there).
- `error_tolerance`: adaptive mode only; the largest distance from a chord midpoint to the curve
  (0 or less: 0.01).
- `use_adaptive_sampling`: non-zero halves each uniform interval recursively (at most 10 levels) until the curve
  point at the interval's middle is within `error_tolerance` of the chord midpoint. Adaptive mode samples in
  parameter, ignoring `traversal_mode`.

---

## Uniform and adaptive sampling into buffers

```c
qaws_status qaws_curve_get_sample_count(
    unsigned int* out_sample_count,
    qaws_curve const* curve,
    qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_2d(
    qaws_vec2* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int position_capacity,
    qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_3d(
    qaws_vec3* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int position_capacity,
    qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_eval_2d(
    qaws_eval_result_2d* out_samples,
    unsigned int* out_sample_count,
    qaws_curve const* curve,
    unsigned int sample_capacity,
    qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_eval_3d(
    qaws_eval_result_3d* out_samples,
    unsigned int* out_sample_count,
    qaws_curve const* curve,
    unsigned int sample_capacity,
    qaws_sampling_desc const* desc);
```

- `qaws_curve_get_sample_count` writes `desc->sample_count`: the count of the non-adaptive path. It does not
  predict the adaptive count; for that, call the sampler and read `*out_count` on `QAWS_STATUS_BUFFER_TOO_SMALL`.
- `sample_2d` / `sample_3d` write positions; `sample_eval_2d` / `sample_eval_3d` write full evaluation results
  (position and the derivatives the curve family provides, with `valid_flags`).
- Non-adaptive: the count is `sample_count` clamped to the capacity (the excess is not written, no error).
- Adaptive: when the result exceeds the capacity, returns `QAWS_STATUS_BUFFER_TOO_SMALL` with the required count
  in the count output and nothing written.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` (a null pointer), `QAWS_STATUS_INVALID_DIMENSION`,
  `QAWS_STATUS_ALLOCATION_FAILURE` (scratch tables, adaptive buffers), or an evaluation status.

---

## Curvature-weighted sampling

```c
qaws_status qaws_curve_sample_curvature_weighted_2d(
    qaws_vec2* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int position_capacity);

qaws_status qaws_curve_sample_curvature_weighted_3d(
    qaws_vec3* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    unsigned int position_capacity);
```

Distributes `sample_count` samples in parameter proportionally to curvature magnitude, so high-curvature regions
receive more points. Curvature is tabulated at uniform parameters, floored at a small epsilon (straight parts still
get samples), and integrated; samples are placed at equal steps of the integral. The first and last samples are the
curve's ends.

- `sample_count` must be at least 2 (`QAWS_STATUS_INVALID_ARGUMENT` otherwise).
- `QAWS_STATUS_BUFFER_TOO_SMALL` when `position_capacity < sample_count`; `*out_position_count` then holds
  `sample_count`.
- Also `QAWS_STATUS_INVALID_DIMENSION`, `QAWS_STATUS_ALLOCATION_FAILURE`, or a curvature/evaluation status.

---

## Feature-preserving sampling

```c
qaws_status qaws_curve_sample_feature_preserving_2d(
    qaws_vec2* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int base_sample_count,
    unsigned int position_capacity);

qaws_status qaws_curve_sample_feature_preserving_3d(
    qaws_vec3* out_positions,
    unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int base_sample_count,
    unsigned int position_capacity);
```

Takes `base_sample_count` uniform parameter samples and inserts the curve's inflection points and the per-axis
extrema (x, y, and z in 3D), found with `qaws_curve_find_inflection_points` and `qaws_curve_find_extrema`
([inspect.md](inspect.md)); at most 128 inflections and 128 extrema per axis. Parameters are kept sorted and
near-duplicates are merged. A feature query that fails is skipped.

- `base_sample_count` must be at least 2.
- `QAWS_STATUS_BUFFER_TOO_SMALL` when the merged count exceeds `position_capacity`; `*out_position_count` then
  holds the needed count and nothing is written. Size the buffer above `base_sample_count`.
- Also `QAWS_STATUS_INVALID_ARGUMENT`, `QAWS_STATUS_INVALID_DIMENSION`, `QAWS_STATUS_ALLOCATION_FAILURE`.

---

## Streaming sampling

```c
typedef int (*qaws_sample_callback_2d)(
    qaws_scalar parameter,
    qaws_vec2 const* position,
    void* user_data);

typedef int (*qaws_sample_callback_3d)(
    qaws_scalar parameter,
    qaws_vec3 const* position,
    void* user_data);

qaws_status qaws_curve_sample_streaming_2d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_sample_callback_2d callback,
    void* user_data);

qaws_status qaws_curve_sample_streaming_3d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_sample_callback_3d callback,
    void* user_data);
```

Delivers each sample to `callback` instead of a buffer. Uses `desc->sample_count` and `desc->traversal_mode`
(parameter or arc length); the adaptive fields are not used. The callback returns non-zero to continue, 0 to
stop early; stopping early returns `QAWS_STATUS_OK`. `position` is valid only during the call. `user_data` is
passed through.

- Returns `QAWS_STATUS_INVALID_ARGUMENT` (null `curve`, `desc` or `callback`), `QAWS_STATUS_INVALID_DIMENSION`,
  `QAWS_STATUS_ALLOCATION_FAILURE` (arc-length table), or an evaluation status.

---

## Adaptive tessellation

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

qaws_status qaws_curve_tessellate_3d(
    qaws_scalar* out_parameters,
    qaws_vec3* out_points,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int capacity,
    qaws_curve_tessellation_desc const* desc);
```

A polyline that follows the curve within the given tolerances.

- **Subdivision.** Each span is halved until every segment passes both enabled tests. The chord test measures the
  distance from the curve to the segment at 1/4, 1/2 and 3/4 of it (three points catch S shapes, whose middle lies
  on the chord). The angle test sums the tangent's turning through those points. Span boundaries are always
  vertices, and a kink there does not count as turning. `max_depth` limits the halvings, which matters near cusps.
- **Determinism.** Vertices come out in parameter order and are the same on every call. Call once with
  `out_parameters` NULL for the count (`out_points` may be NULL), then again with buffers of that size.
- **Errors.** `QAWS_STATUS_BUFFER_TOO_SMALL` when `capacity` is below the count; `*out_count` holds the count and
  the buffers hold the first `capacity` vertices. `QAWS_STATUS_INVALID_ARGUMENT` for a null `curve`, `desc` or
  `out_count`, a negative tolerance, or both tolerances 0. `QAWS_STATUS_INVALID_DIMENSION` for a dimension
  mismatch. Nothing is allocated.
- **Limits.** The tests sample three points per segment, so they are estimates, not bounds.
