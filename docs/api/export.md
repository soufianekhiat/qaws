# API Reference: qaws_export.h

Least-squares B-spline fitting, and export of curves and surfaces to SVG path data, polylines and Wavefront OBJ
text. Export writes into caller buffers and allocates nothing; fitting creates a new curve.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Related: [import.md](import.md) (the reverse direction), [bspline.md](bspline.md), [sampling.md](sampling.md),
[inspect.md](inspect.md) (`qaws_surface_tessellate_adaptive` for meshes by tolerance), [diff.md](diff.md)
(`qaws_curve_fit_bspline_diff`, the differentiable fit).

---

## B-spline fitting

```c
typedef struct qaws_bspline_fit_desc
{
    qaws_dimension dimension;
    void const* data_points;           /* m * dim_count scalars (flat array of qaws_vec2 or qaws_vec3) */
    unsigned int data_point_count;     /* m >= degree + 1 */
    unsigned int degree;               /* desired B-spline degree (typically 3) */
    unsigned int control_point_count;  /* n: desired number of CPs (degree < n <= m) */
    qaws_scalar const* parameters;     /* optional: m parameter values in [0,1], NULL = from the points as below */
    qaws_parameterization parameterization;  /* without parameters: DEFAULT and CHORDAL chord lengths, UNIFORM equal steps,
                                                CENTRIPETAL square roots of the chord lengths */
} qaws_bspline_fit_desc;

qaws_status qaws_curve_fit_bspline(
    qaws_curve** out_curve,
    qaws_bspline_fit_desc const* desc);
```

Fits a clamped B-spline of `degree` with `control_point_count` control points to `data_point_count` points by
least squares, on the domain [0, 1].

- Parameters: `desc->parameters` when given (m values in [0, 1]), otherwise derived from the points by
  `parameterization`: `QAWS_PARAMETERIZATION_DEFAULT` and `_CHORDAL` step by chord lengths, `_UNIFORM` by equal
  steps, `_CENTRIPETAL` by square roots of the chord lengths (steadier through sharp turns and uneven spacing).
  Other values are refused.
- Knots are placed by averaging the parameters (Piegl and Tiller). The first and last data points are
  interpolated exactly; with `control_point_count == data_point_count` the fit interpolates every point at its
  parameter.
- The new curve is owned by the caller (`qaws_curve_destroy`); `*out_curve` is NULL on failure.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` (null pointer, unknown parameterization, `m < degree + 1`),
  `QAWS_STATUS_INVALID_DIMENSION`, `QAWS_STATUS_INVALID_DEGREE` (degree 0),
  `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` (`n <= degree` or `n > m`), `QAWS_STATUS_NUMERICAL_FAILURE` (singular
  normal equations), `QAWS_STATUS_ALLOCATION_FAILURE`.

---

## SVG path data

```c
qaws_status qaws_curve_export_svg_path(
    char* out_path_data,
    unsigned int* out_length,
    qaws_curve const* curve,
    unsigned int capacity,
    unsigned int sample_count);
```

Writes a 2D curve as SVG path data (the `d` attribute), numbers with 6 significant digits.

| Family | Commands |
|---|---|
| Bezier of degree 1, 2, 3 | `M` then `L`, `Q` or `C`, exact |
| Bezier of higher degree | `M` then `4 x degree` line segments |
| Hermite, Catmull-Rom | one exact `C` per span |
| Rational Bezier | `M` then `sample_count` line segments |
| Every other family | `M` then `sample_count` (at least 4) cubic `C` segments built from positions and tangents at uniform parameters (approximate) |

- `sample_count`: 0 means 64.
- `capacity` is the buffer size in bytes including the terminator; the string is always NUL-terminated.
  `*out_length` is the length written; `*out_length >= capacity` means the string was truncated (the call still
  returns `QAWS_STATUS_OK`).
- Returns `QAWS_STATUS_INVALID_ARGUMENT` (null pointer), `QAWS_STATUS_INVALID_DIMENSION` (3D curve),
  `QAWS_STATUS_BUFFER_TOO_SMALL` (`capacity == 0`), or an evaluation status.

---

## Polylines

```c
typedef enum qaws_polyline_sampling
{
    QAWS_POLYLINE_UNIFORM   = 0,   /* uniform parameter spacing (default) */
    QAWS_POLYLINE_CURVATURE = 1    /* curvature-weighted: more points where the curve bends */
} qaws_polyline_sampling;

qaws_status qaws_polyline_export_2d(
    qaws_vec2* out_points,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    qaws_polyline_sampling sampling);

qaws_status qaws_polyline_export_3d(
    qaws_vec3* out_points,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    qaws_polyline_sampling sampling);
```

Samples any curve of the matching dimension into `sample_count` points (at least 2); `out_points` must hold
`sample_count` elements.

- Uniform: parameters evenly spaced, the last exactly at the domain's end. On an evaluation error `*out_count` is
  the number of points written before it.
- Curvature: `qaws_curve_sample_curvature_weighted_2d/3d` ([sampling.md](sampling.md)); `*out_count` is
  `sample_count` on success.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` (null pointer, `sample_count < 2`), `QAWS_STATUS_INVALID_DIMENSION`.

---

## OBJ surface mesh

```c
qaws_status qaws_surface_export_obj(
    char* out_obj_data,
    unsigned int* out_length,
    qaws_surface const* surface,
    unsigned int capacity,
    unsigned int u_samples,
    unsigned int v_samples,
    int include_normals);
```

Writes a surface as Wavefront OBJ text: the surface sampled on a `u_samples x v_samples` grid (each at least 2)
over its domain, `v` lines, then (with `include_normals` non-zero) `vn` lines, then two triangles per grid cell as
1-based `f` lines (`f a//a b//b c//c` with normals).

- With normals, the orientation is chosen so normals point outward: flipped (normals negated and winding reversed)
  when the surface closes in v, or when the center normal points toward the centroid of the four corners.
- `capacity` is the buffer size in bytes including the terminator; the text is always NUL-terminated and
  `*out_length` is the length written.
- Returns `QAWS_STATUS_BUFFER_TOO_SMALL` when the text did not fit (or `capacity == 0`),
  `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer or a grid below 2 x 2.
