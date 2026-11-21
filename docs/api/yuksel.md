# qaws_yuksel.h: Yuksel C2 interpolating splines

A spline through every control point, C2 everywhere, built by blending overlapping sub-curves (one per point triplet)
with cos^2 / sin^2 weights (Cem Yuksel, "A Class of C2 Interpolating Splines", ACM Transactions on Graphics, 2020).
Each name exists as `name_f32` and `name_f64` ([precision.md](precision.md)).

## qaws_yuksel_mode

The sub-curve fitted through each triplet.

| Value | Sub-curve |
|---|---|
| `QAWS_YUKSEL_MODE_BEZIER` (0) | quadratic Bezier |
| `QAWS_YUKSEL_MODE_CIRCULAR` | circular arc |
| `QAWS_YUKSEL_MODE_ELLIPTICAL` | elliptical arc |
| `QAWS_YUKSEL_MODE_HYBRID` | circular arc when its span is within pi/2, else elliptical |

All modes work in 2D and 3D; in 3D the arc modes fit each triplet in its own plane.

## qaws_yuksel_desc

```c
typedef struct qaws_yuksel_desc {
    qaws_dimension dimension;
    void const* control_points;
    unsigned int control_point_count;
    qaws_yuksel_mode mode;
    int closed;
} qaws_yuksel_desc;
```

| Member | Meaning |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `control_points` | `dimension` scalars per point, packed |
| `control_point_count` | at least 3 |
| `mode` | see above |
| `closed` | non-zero: the last point joins the first |

## qaws_curve_create_yuksel

```c
qaws_status qaws_curve_create_yuksel(
    qaws_curve** out_curve,
    qaws_yuksel_desc const* desc);
```

Creates the curve. It passes through every point; point i is at parameter i.

- Open: `control_point_count - 1` spans. The first and last spans use the end sub-curve alone (unblended); interior
  spans blend two neighbouring sub-curves.
- Closed: `control_point_count` spans.
- Parameter range [0, span_count]. `qaws_curve_get_degree` reports 2.
- Inputs are copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)). Reversing an arc-mode curve
  gives an extracted copy, not a Yuksel curve ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc`, `out_curve` or `control_points` NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | fewer than 3 points, or a count too large to index |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_yuksel_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
