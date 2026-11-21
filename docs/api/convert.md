# API Reference: qaws_convert.h

Conversions between curve families and Bezier degree changes. Each function creates a new curve, owned by the
caller (destroy it with `qaws_curve_destroy`, see [curve.md](curve.md)), allocated with the source curve's
allocator; the output is set to NULL on failure. All are thread-safe: they only read the input curve.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Common statuses: `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer, `QAWS_STATUS_UNSUPPORTED_OPERATION` when the
source is not of the required kind, `QAWS_STATUS_ALLOCATION_FAILURE`.

Related: [bezier.md](bezier.md), [bspline.md](bspline.md), [nurbs.md](nurbs.md), [hermite.md](hermite.md),
[catmull_rom.md](catmull_rom.md), [operations.md](operations.md) (`qaws_curve_extract` converts pieces of any
family).

---

## Family conversions

```c
qaws_status qaws_curve_convert_hermite_to_bezier(
    qaws_curve** out_bezier,
    qaws_curve const* curve,
    unsigned int span_index);

qaws_status qaws_curve_convert_catmull_rom_to_bezier(
    qaws_curve** out_bezier,
    qaws_curve const* curve,
    unsigned int span_index);

qaws_status qaws_curve_convert_bezier_to_bspline(
    qaws_curve** out_bspline,
    qaws_curve const* curve);

qaws_status qaws_curve_convert_bspline_to_nurbs(
    qaws_curve** out_nurbs,
    qaws_curve const* curve);
```

| Function | Source | Result |
|---|---|---|
| `convert_hermite_to_bezier` | Hermite | one span as a cubic Bezier: `P0`, `P0 + M0/3`, `P1 - M1/3`, `P1` |
| `convert_catmull_rom_to_bezier` | Catmull-Rom | one span (its cubic over [0, 1]) as a cubic Bezier |
| `convert_bezier_to_bspline` | Bezier | B-spline of the same degree and control points, clamped knots `[0 x (degree+1), 1 x (degree+1)]` |
| `convert_bspline_to_nurbs` | B-spline | NURBS with the same control points and knots, all weights 1 |

- The per-span conversions take a 0-based `span_index`; `QAWS_STATUS_OUT_OF_RANGE` when it is not below the span
  count. The Bezier's parameter runs over [0, 1] for that span.
- All four conversions are exact.

---

## Degree elevation and reduction

```c
qaws_status qaws_curve_elevate_degree(
    qaws_curve** out_elevated,
    qaws_curve const* curve);

qaws_status qaws_curve_reduce_degree(
    qaws_curve** out_reduced,
    qaws_curve const* curve);
```

Both accept Bezier curves only (`QAWS_STATUS_UNSUPPORTED_OPERATION` for other kinds).

- `elevate_degree`: degree n to n + 1, exact (`Q[i] = i/(n+1) P[i-1] + (1 - i/(n+1)) P[i]`).
- `reduce_degree`: degree n to n - 1, an approximation: the least-squares fit (L2 over [0, 1]) that keeps both end
  points. `QAWS_STATUS_INVALID_DEGREE` below degree 2; `QAWS_STATUS_NUMERICAL_FAILURE` if the normal equations are
  singular.
