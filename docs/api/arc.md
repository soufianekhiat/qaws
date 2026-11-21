# qaws_arc.h: circular arcs

A curve made of one or more circular arc segments, parameterized by arc length. Each name exists as `name_f32` and
`name_f64` ([precision.md](precision.md)).

## qaws_arc_segment

```c
typedef struct qaws_arc_segment {
    qaws_scalar center[3];
    qaws_scalar radius;
    qaws_scalar angle_start;    /* radians */
    qaws_scalar angle_end;      /* radians */
    qaws_scalar axis_u[3];      /* 3D only: unit vector, cos direction */
    qaws_scalar axis_v[3];      /* 3D only: unit vector, sin direction */
} qaws_arc_segment;
```

| Member | Meaning |
|---|---|
| `center` | the arc's center; in 2D only x and y are used |
| `radius` | must be > 0 |
| `angle_start`, `angle_end` | radians; in 2D counter-clockwise from +x. The sweep is `angle_end - angle_start` (negative: clockwise) and must not be zero |
| `axis_u`, `axis_v` | 3D only (ignored in 2D): the plane basis; the point at angle a is `center + radius * (cos(a) * axis_u + sin(a) * axis_v)` |

## qaws_arc_desc

```c
typedef struct qaws_arc_desc {
    qaws_dimension dimension;
    qaws_arc_segment const* segments;
    unsigned int segment_count;    /* >= 1 */
} qaws_arc_desc;
```

## qaws_curve_create_arc

```c
qaws_status qaws_curve_create_arc(
    qaws_curve** out_curve,
    qaws_arc_desc const* desc);
```

Creates the curve. Each segment is one span, in array order; the parameter is the cumulative arc length, range
[0, total length]. The segments are used as given: they are not moved to meet, so give end-to-end segments for a
connected curve. The degree reported is 2. Segments are copied. Destroy the result with `qaws_curve_destroy`
([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `segments` NULL; `segment_count` 0; a radius <= 0 |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_DEGENERATE_CURVE` | a segment with (near) zero sweep or length |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_arc_ex` is declared in `qaws_alloc.h`
([alloc.md](alloc.md)).
