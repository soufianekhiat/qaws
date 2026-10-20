# API Reference: qaws_export.h

Serialization, export, and curve fitting.

---

## qaws_bspline_fit_desc

```c
typedef struct qaws_bspline_fit_desc {
    qaws_dimension dimension;
    void const* data_points;
    unsigned int data_point_count;
    unsigned int degree;
    unsigned int control_point_count;
    qaws_scalar const* parameters;
    qaws_parameterization parameterization;
} qaws_bspline_fit_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `data_points` | Flat array of `qaws_vec2` or `qaws_vec3` (m points) |
| `data_point_count` | Number of data points (m). Must be >= `degree + 1`. |
| `degree` | Desired B-spline degree (typically 3) |
| `control_point_count` | Desired number of control points (n). Must satisfy `degree < n <= m`. |
| `parameters` | Optional: m parameter values in [0, 1]. Pass NULL to derive them from the points by `parameterization`. |
| `parameterization` | Without `parameters`: `QAWS_PARAMETERIZATION_DEFAULT` and `_CHORDAL` step by the chord lengths, `_UNIFORM` by equal steps, `_CENTRIPETAL` by the square roots of the chord lengths (steadier through sharp turns and uneven spacing). Others are refused. `qaws_curve_fit_bspline_diff` differentiates through all three (test 55). |

With `control_point_count` equal to `data_point_count` the fit interpolates the points at those parameters (test 153).

---

## qaws_curve_fit_bspline

```c
qaws_status qaws_curve_fit_bspline(
    qaws_curve** out_curve,
    qaws_bspline_fit_desc const* desc);
```

Fits a B-spline curve to a set of data points using least-squares approximation.

**Parameters:**
- `desc` -- Fitting descriptor.
- `out_curve` -- Receives the fitted B-spline curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_export_svg_path

```c
qaws_status qaws_curve_export_svg_path(
    char* out_path_data,
    unsigned int* out_length,
    qaws_curve const* curve,
    unsigned int capacity,
    unsigned int sample_count);
```

Converts a 2D curve to an SVG path data string (the `d` attribute).

**Parameters:**
- `curve` -- The 2D curve to export.
- `out_path_data` -- Caller-allocated buffer for the SVG path string.
- `capacity` -- Buffer size in bytes (including null terminator).
- `out_length` -- Receives the actual string length written (excluding null terminator). If `out_length >= capacity`, the string was truncated.
- `sample_count` -- Number of sample points for non-Bezier families. Pass 0 for the default (64).

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- For Bezier curves: maps directly to M/C/Q/L commands (exact).
- For Hermite/Catmull-Rom curves: converts each span to cubic Bezier (exact).
- For all other families: samples the curve and fits piecewise cubic Bezier segments (approximate, controlled by `sample_count`).
- The curve must be 2D.

---

## qaws_polyline_sampling

```c
typedef enum qaws_polyline_sampling {
    QAWS_POLYLINE_UNIFORM   = 0,
    QAWS_POLYLINE_CURVATURE = 1
} qaws_polyline_sampling;
```

| Value | Description |
|---|---|
| `QAWS_POLYLINE_UNIFORM` | Uniform parameter spacing (default) |
| `QAWS_POLYLINE_CURVATURE` | Curvature-weighted -- more points in high-curvature regions, fewer in straight segments |

---

## qaws_polyline_export_2d

```c
qaws_status qaws_polyline_export_2d(
    qaws_vec2* out_points,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    qaws_polyline_sampling sampling);
```

Samples a 2D curve and outputs positions as a flat point array (polyline).

**Parameters:**
- `curve` -- Any 2D curve (any family).
- `sample_count` -- Number of sample points (>= 2).
- `sampling` -- `QAWS_POLYLINE_UNIFORM` or `QAWS_POLYLINE_CURVATURE`.
- `out_points` -- Caller-allocated buffer. Must hold at least `sample_count` elements.
- `out_count` -- Receives the number of points actually written. For uniform sampling this equals `sample_count`. For curvature sampling this may be <= `sample_count`.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_polyline_export_3d

```c
qaws_status qaws_polyline_export_3d(
    qaws_vec3* out_points,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int sample_count,
    qaws_polyline_sampling sampling);
```

Same as `qaws_polyline_export_2d` but for 3D curves. The `out_points` buffer holds `qaws_vec3` elements.
