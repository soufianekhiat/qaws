# qaws_inline.h: inline curves

Initialize a curve inside a caller-held `qaws_curve_inline` (see [types.md](types.md)) with no heap allocation. The
result works with every evaluation and inspection function through `qaws_curve_inline_get_curve(&ic)`. Each name
exists as `name_f32` and `name_f64` ([precision.md](precision.md)).

Constraints:

- Bezier and polynomial families only, one span (a Bezier has `degree + 1` control points).
- Degree 1 to 7 (8 control points x 3 scalars), 2D or 3D, within `QAWS_INLINE_STORAGE_BYTES` (512).
- Never call `qaws_curve_destroy` on it; it holds no resources, so it can simply go out of scope.
- It is self-contained: it may be copied with memcpy or assignment; call `qaws_curve_inline_get_curve(&copy)` once on
  the copy (it re-points the copy at its own storage) before other threads use it.

```c
qaws_curve_inline ic;
qaws_vec3 pts[] = { {0, 0, 0}, {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
qaws_bezier_desc desc = { QAWS_DIMENSION_3D, 3, pts, 4 };
if (qaws_curve_init_bezier_inline(&ic, &desc) == QAWS_STATUS_OK) {
    qaws_curve* c = qaws_curve_inline_get_curve(&ic);
    /* qaws_curve_evaluate_3d(c, ...) */
}
```

## qaws_curve_init_bezier_inline

```c
qaws_status qaws_curve_init_bezier_inline(
    qaws_curve_inline* out,
    qaws_bezier_desc const* desc);
```

Initializes `out` as the Bezier curve `desc` describes ([bezier.md](bezier.md)), parameter range [0, 1]. The control
points are copied into `out`.

## qaws_curve_init_polynomial_inline

```c
qaws_status qaws_curve_init_polynomial_inline(
    qaws_curve_inline* out,
    qaws_polynomial_desc const* desc);
```

Initializes `out` as the polynomial curve `desc` describes ([polynomial.md](polynomial.md)), parameter range
[t_min, t_max] (both 0: [0, 1]). The coefficients are copied into `out`. Unlike `qaws_curve_create_polynomial`, it
does not check that `t_max > t_min`.

## Statuses (both functions)

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | initialized |
| `QAWS_STATUS_INVALID_ARGUMENT` | `out`, `desc` or the point/coefficient array NULL |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D |
| `QAWS_STATUS_INVALID_DEGREE` | degree 0 or above 7 |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | count not `degree + 1` |
| `QAWS_STATUS_BUFFER_TOO_SMALL` | the layout does not fit the inline storage |
