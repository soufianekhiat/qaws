# Sampling

Sampling generates arrays of points (or full evaluation results) along a curve. This is the primary way to convert curves into renderable polylines.

## Sampling descriptor

```c
qaws_sampling_desc desc;
desc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;  // or ARC_LENGTH
desc.sample_count = 100;
desc.error_tolerance = 0;       // for adaptive sampling
desc.use_adaptive_sampling = 0; // 0 = uniform, 1 = adaptive
```

## Sampling modes

| Mode | Description |
|---|---|
| `QAWS_TRAVERSAL_MODE_PARAMETER` | Uniform spacing in parameter space |
| `QAWS_TRAVERSAL_MODE_ARC_LENGTH` | Uniform spacing in distance along curve |

Arc-length sampling produces visually even spacing regardless of parameterization speed.

## Sample positions only

```c
qaws_vec2 positions[256];
unsigned int count = 0;

qaws_sampling_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER,
    .sample_count = 256,
    .error_tolerance = 0,
    .use_adaptive_sampling = 0
};

qaws_curve_sample_2d(positions, &count, curve, 256, &desc);
// positions[0..count-1] now contain sampled points
```

## Sample with derivatives

When you need tangent or curvature information at each sample:

```c
qaws_eval_result_2d samples[256];
unsigned int count = 0;

qaws_curve_sample_eval_2d(samples, &count, curve, 256, &desc);
// samples[i].position, samples[i].d1, etc.
```

## Query sample count before allocation

```c
unsigned int count;
qaws_curve_get_sample_count(&count, curve, &desc);
// allocate buffer of size `count`, then call sample_2d/3d
```

## Adaptive sampling

When `use_adaptive_sampling` is set to 1, the sampling functions use recursive subdivision instead of uniform spacing. Segments where the chord-to-curve midpoint deviation exceeds `error_tolerance` are subdivided (up to 10 levels deep).

```c
qaws_sampling_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER,
    .sample_count = 16,           // initial uniform segments
    .error_tolerance = 0.01f,     // max chord-to-curve deviation
    .use_adaptive_sampling = 1
};

qaws_vec2 positions[4096];
unsigned int count = 0;

qaws_status s = qaws_curve_sample_2d(
    positions, &count, curve, 4096, &desc);
// count may be larger than sample_count due to subdivision
```

- `sample_count` sets the initial number of uniform segments to subdivide from.
- Tighter `error_tolerance` produces more points in high-curvature regions.
- Returns `QAWS_STATUS_BUFFER_TOO_SMALL` if the output buffer is insufficient.
- Works with both `sample_2d`/`sample_3d` (positions only). Adaptive mode is not used with `sample_eval_*` or arc-length traversal mode.

## Curvature-weighted sampling

Distributes samples proportionally to curvature magnitude so high-curvature regions receive more points. Unlike adaptive sampling (which subdivides recursively), curvature-weighted sampling takes a desired sample count and redistributes points according to curvature.

```c
qaws_vec2 positions[512];
unsigned int count = 0;

qaws_curve_sample_curvature_weighted_2d(
    // desired sample count
    positions,
    // buffer capacity
    &count,            curve,
    512,            512);
// count <= 512, concentrated in high-curvature areas
```

Signatures:

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

## Feature-preserving sampling

Inserts inflection points and extrema into a uniform parameter sample set. This ensures that sharp features and inflection points are never missed, regardless of the base sample count.

```c
qaws_vec2 positions[512];
unsigned int count = 0;

qaws_curve_sample_feature_preserving_2d(
    // base uniform sample count
    positions,
    // buffer capacity (may exceed base count)
    &count,             curve,
    64,            512);
// count >= 64, with inflection points and extrema injected
```

Signatures:

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

## Streaming sampling

Callback-based sampling that requires no output buffer allocation. A user-provided callback is invoked for each sample point. Return non-zero from the callback to continue, or 0 to stop early (not an error).

```c
int my_callback(qaws_scalar parameter,
                qaws_vec2 const* position,
                void* user_data) {
    // process sample point...
    return 1;  // continue (return 0 to stop early)
}

qaws_sampling_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER,
    .sample_count = 200,
    .error_tolerance = 0,
    .use_adaptive_sampling = 0
};

qaws_curve_sample_streaming_2d(curve, &desc, my_callback, NULL);
```

Callback types:

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

Signatures:

```c
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

## 3D variants

All sampling functions have 3D equivalents:

- `qaws_curve_sample_3d`
- `qaws_curve_sample_eval_3d`
- `qaws_curve_sample_curvature_weighted_3d`
- `qaws_curve_sample_feature_preserving_3d`
- `qaws_curve_sample_streaming_3d`
