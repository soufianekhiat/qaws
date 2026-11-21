# qaws_subdivision.h: subdivision curves

A curve made by refining a control polygon a number of times with a corner-cutting scheme, then evaluated as cubic
spans through the refined points. Each name exists as `name_f32` and `name_f64` ([precision.md](precision.md)).

## qaws_subdivision_scheme

| Value | Scheme | Limit curve (reported by `qaws_curve_get_continuity`) |
|---|---|---|
| `QAWS_SUBDIVISION_CHAIKIN` (0) | Chaikin corner cutting | quadratic B-spline, C1 |
| `QAWS_SUBDIVISION_LANE_RIESENFELD_3` | Lane-Riesenfeld, degree 3 | cubic B-spline, C2 |
| `QAWS_SUBDIVISION_LANE_RIESENFELD_4` | Lane-Riesenfeld, degree 4 | quartic B-spline, C3 |

## qaws_subdivision_desc

```c
typedef struct qaws_subdivision_desc {
    qaws_dimension dimension;
    qaws_subdivision_scheme scheme;
    void const* control_points;          /* initial control polygon */
    unsigned int control_point_count;    /* >= 3 */
    int closed;                          /* 1 for closed polygon */
    unsigned int refinement_levels;      /* number of subdivision iterations (default 6) */
} qaws_subdivision_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `scheme` | see above |
| `control_points` | the initial polygon, `dimension` scalars per point, packed |
| `control_point_count` | at least 3 |
| `closed` | non-zero: a closed polygon |
| `refinement_levels` | refinement steps; 0 means 6. Each step about doubles the point count |

## qaws_curve_create_subdivision

```c
qaws_status qaws_curve_create_subdivision(
    qaws_curve** out_curve,
    qaws_subdivision_desc const* desc);
```

Refines the polygon and creates the curve through the refined points: R - 1 spans open, R spans closed (R refined
points), parameter range [0, span_count]. More levels give a closer approximation of the limit curve and more spans.
The input points are copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `control_points` NULL; unknown scheme |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | fewer than 3 points, a count too large to index, or too few refined points (under 4 open, 3 closed) |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed, or the refined count would overflow |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_subdivision_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
