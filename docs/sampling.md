# Sampling and Tessellation

Sampling turns a curve into arrays of points (or full evaluation results);
tessellation does the same by tolerances, for curves and surfaces, into caller
buffers. Names below are the unsuffixed forms ([precision.md](precision.md)).
Every function here is thread-safe on immutable curves and surfaces, and has a
2D and a 3D form where the geometry allows; a curve of the other dimension
gives `QAWS_STATUS_INVALID_DIMENSION`.

| Need | Function | Header |
|---|---|---|
| N points, uniform in parameter or arc length | `qaws_curve_sample_2d/3d` | `qaws_sampling.h` |
| N full results (position, D1, D2, D3) | `qaws_curve_sample_eval_2d/3d` | `qaws_sampling.h` |
| Points refined by chord deviation | `qaws_curve_sample_*` with `use_adaptive_sampling` | `qaws_sampling.h` |
| More points where curvature is high | `qaws_curve_sample_curvature_weighted_2d/3d` | `qaws_sampling.h` |
| Inflections and extrema never missed | `qaws_curve_sample_feature_preserving_2d/3d` | `qaws_sampling.h` |
| No output buffer (callback) | `qaws_curve_sample_streaming_2d/3d` | `qaws_sampling.h` |
| Curve vertices by chord and angle tolerance, deterministic counts | `qaws_curve_tessellate_2d/3d` | `qaws_sampling.h` |
| Surface position grid | `qaws_surface_sample` | `qaws_surface.h` |
| Surface triangle mesh (depth, curvature, edge length) | `qaws_surface_tessellate` | `qaws_inspect.h` |
| Surface triangle mesh by chord and angle tolerance, crack-free | `qaws_surface_tessellate_adaptive` | `qaws_inspect.h` |

Differentiable sampling (inverse-CDF samplers whose samples move with the
control points) is in `qaws_diff_sampling.h`
([differentiation.md](differentiation.md)).

## Sampling descriptor

```c
typedef struct qaws_sampling_desc {
    qaws_traversal_mode traversal_mode;   /* PARAMETER or ARC_LENGTH */
    unsigned int sample_count;            /* uniform: the count; adaptive: initial segments */
    qaws_scalar error_tolerance;          /* adaptive: chord-to-curve deviation */
    int use_adaptive_sampling;            /* 0 uniform, 1 adaptive */
} qaws_sampling_desc;
```

| Mode | Spacing |
|---|---|
| `QAWS_TRAVERSAL_MODE_PARAMETER` | Uniform in parameter |
| `QAWS_TRAVERSAL_MODE_ARC_LENGTH` | Uniform in distance along the curve (from an arc-length table) |

## Uniform sampling

```c
qaws_vec2 positions[256];
unsigned int count = 0;
qaws_sampling_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    .sample_count = 256,
    .error_tolerance = 0,
    .use_adaptive_sampling = 0
};
qaws_curve_sample_2d(positions, &count, curve, 256, &desc);
```

`qaws_curve_sample_eval_2d/3d` take the same descriptor and fill
`qaws_eval_result_2d/3d` records. In uniform mode the count written is
`sample_count`, cut to the capacity. `qaws_curve_get_sample_count` returns
`sample_count`, the uniform count; it does not predict adaptive counts.

## Adaptive sampling

With `use_adaptive_sampling = 1`, each of `sample_count` initial uniform
segments is halved recursively (at most 10 levels) while the curve's midpoint
deviates from the chord by more than `error_tolerance`. Adaptive mode works
in parameter space (`traversal_mode` is not used) and applies to
`qaws_curve_sample_*` and `qaws_curve_sample_eval_*`. When the result does not
fit, the call returns `QAWS_STATUS_BUFFER_TOO_SMALL` with the needed count in
the count output; the points are built in a temporary buffer.

```c
qaws_sampling_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER,
    .sample_count = 16,
    .error_tolerance = 0.01f,
    .use_adaptive_sampling = 1
};
qaws_vec2 positions[4096];
unsigned int count = 0;
qaws_status s = qaws_curve_sample_2d(positions, &count, curve, 4096, &desc);
```

## Curvature-weighted and feature-preserving sampling

```c
qaws_status qaws_curve_sample_curvature_weighted_2d(
    qaws_vec2* out_positions, unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int sample_count,        /* desired number of samples */
    unsigned int position_capacity);

qaws_status qaws_curve_sample_feature_preserving_2d(
    qaws_vec2* out_positions, unsigned int* out_position_count,
    qaws_curve const* curve,
    unsigned int base_sample_count,   /* uniform base samples */
    unsigned int position_capacity);  /* room for the inserted features too */
```

Curvature-weighted sampling places `sample_count` points with density
proportional to curvature magnitude. Feature-preserving sampling inserts the
inflection points and extrema into a uniform parameter sample set, so the
count can exceed `base_sample_count`. Both have `_3d` forms.

## Streaming sampling

```c
typedef int (*qaws_sample_callback_2d)(qaws_scalar parameter, qaws_vec2 const* position, void* user_data);
typedef int (*qaws_sample_callback_3d)(qaws_scalar parameter, qaws_vec3 const* position, void* user_data);

qaws_status qaws_curve_sample_streaming_2d(qaws_curve const* curve, qaws_sampling_desc const* desc,
    qaws_sample_callback_2d callback, void* user_data);
```

The callback runs once per sample; return non-zero to continue, 0 to stop
early (not an error).

## Curve tessellation by tolerances

```c
typedef struct qaws_curve_tessellation_desc {
    qaws_scalar chord_tolerance;   /* distance; 0: unused */
    qaws_scalar angle_tolerance;   /* radians; 0: unused (one of the two is needed) */
    unsigned int max_depth;        /* halvings per span at most; 0: 16, at most 30 */
} qaws_curve_tessellation_desc;

qaws_status qaws_curve_tessellate_2d(qaws_scalar* out_parameters, qaws_vec2* out_points,
    unsigned int* out_count, qaws_curve const* curve, unsigned int capacity,
    qaws_curve_tessellation_desc const* desc);
```

Each span is halved until every segment is within the chord tolerance (tested
at 1/4, 1/2 and 3/4 of it) and the tangent turns by at most the angle
tolerance along it. Span boundaries are always vertices, so a kink at a
repeated knot is kept. The result is deterministic: call with
`out_parameters` NULL for the count, then with buffers of that size.
`out_points` may be NULL. A short buffer gives `QAWS_STATUS_BUFFER_TOO_SMALL`
with the count, the buffers holding the first `capacity` vertices. No
allocation.

## Surfaces

- `qaws_surface_sample(out_positions, &count, surface, u_count, v_count, capacity)`:
  a `u_count x v_count` grid of positions, row major with v varying fastest.
- `qaws_surface_tessellate`: a triangle mesh of `qaws_tessellation_vertex`
  (position, normal, u, v) and indices into caller buffers, quads split by
  `qaws_tessellation_desc` (`max_depth` 0: 5, `curvature_threshold` 0: 0.1,
  `max_edge_length` 0: off).
- `qaws_surface_tessellate_adaptive`: a quadtree over (u, v) whose quads split
  until the surface at their 3 x 3 test points is within `chord_tolerance` of
  the bilinear patch of their corners and the normal turns by at most
  `angle_tolerance` (`min_depth` 0: 2; `max_depth` 0: 10, at most 20). No
  cracks, no allocation, deterministic: call with NULL buffers for the counts,
  then with buffers of those sizes. Triangles are counter-clockwise in
  (u, v). Vertices are written per quad, so a shared point repeats (with the
  same bits); there is no 2:1 grading between neighbouring quads.
