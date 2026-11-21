# qaws_composite.h: composite curves

One curve made of several existing curves of any kinds, chained in order. Each name exists as `name_f32` and
`name_f64` ([precision.md](precision.md)).

## qaws_composite_desc

```c
typedef struct qaws_composite_desc {
    qaws_dimension dimension;
    qaws_curve** segments;           /* array of segment curves; ownership transfers to composite */
    unsigned int segment_count;      /* >= 1 */
} qaws_composite_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D`; every segment must have it |
| `segments` | the segment curves, in order; none NULL |
| `segment_count` | at least 1 |

## qaws_curve_create_composite

```c
qaws_status qaws_curve_create_composite(
    qaws_curve** out_curve,
    qaws_composite_desc const* desc);
```

Creates the composite. Segment i occupies [i, i + 1] of the composite parameter (its own domain mapped onto that unit
interval); the range is [0, segment_count]. The degree reported is the largest segment degree. Segments are used as
given: they are not moved to meet.

Ownership: on success the composite owns every curve in `segments` and destroys them when it is destroyed; the caller
must not destroy them. The `segments` array itself is copied and stays the caller's. On failure the caller keeps the
curves. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `segments` NULL; `segment_count` 0; a NULL segment |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D, or a segment of another dimension |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is written only on success. The allocator-taking form `qaws_curve_create_composite_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
