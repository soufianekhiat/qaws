# qaws_catmull_rom.h: Catmull-Rom splines

A cubic interpolating spline through control points, with uniform, chordal or centripetal knot spacing. Each name
exists as `name_f32` and `name_f64` ([precision.md](precision.md)).

## qaws_catmull_rom_desc

```c
typedef struct qaws_catmull_rom_desc {
    qaws_dimension dimension;
    void const* control_points;
    unsigned int control_point_count;
    qaws_parameterization parameterization;
    int closed;
} qaws_catmull_rom_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `control_points` | `dimension` scalars per point, packed |
| `control_point_count` | at least 4 open, at least 3 closed |
| `parameterization` | `UNIFORM` (alpha 0), `CHORDAL` (alpha 1), `CENTRIPETAL` (alpha 0.5); `DEFAULT` is centripetal |
| `closed` | non-zero: the last point joins the first |

## qaws_curve_create_catmull_rom

```c
qaws_status qaws_curve_create_catmull_rom(
    qaws_curve** out_curve,
    qaws_catmull_rom_desc const* desc);
```

Creates the curve. Degree 3, C1 at span joins.

- Open: the first and last points only shape the end tangents; the curve runs from point 1 to point N-2 in N-3 spans.
- Closed: every point is interpolated, N spans.
- Span i is on [i, i + 1]; parameter range [0, span_count] whatever the parameterization (the knot spacing shapes
  each span, it does not change the parameter domain).

Inputs are copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `control_points` NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | too few points (see above) or a count too large to index |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_catmull_rom_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
