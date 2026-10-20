# API Reference: qaws_convert.h

Conversion between curve types and degree manipulation.

Thread-safe: operates on immutable input curves.

---

## qaws_curve_convert_hermite_to_bezier

```c
qaws_status qaws_curve_convert_hermite_to_bezier(
    qaws_curve** out_bezier,
    qaws_curve const* curve,
    unsigned int span_index);
```

Converts a single span of a Hermite curve to a Bezier curve.

**Parameters:**
- `curve` -- A Hermite curve.
- `span_index` -- 0-based index of the span to convert.
- `out_bezier` -- Receives the equivalent Bezier curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_convert_catmull_rom_to_bezier

```c
qaws_status qaws_curve_convert_catmull_rom_to_bezier(
    qaws_curve** out_bezier,
    qaws_curve const* curve,
    unsigned int span_index);
```

Converts a single span of a Catmull-Rom curve to a Bezier curve.

**Parameters:**
- `curve` -- A Catmull-Rom curve.
- `span_index` -- 0-based index of the span to convert.
- `out_bezier` -- Receives the equivalent Bezier curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_convert_bezier_to_bspline

```c
qaws_status qaws_curve_convert_bezier_to_bspline(
    qaws_curve** out_bspline,
    qaws_curve const* curve);
```

Converts a Bezier curve to a B-Spline curve.

**Parameters:**
- `curve` -- A Bezier curve.
- `out_bspline` -- Receives the equivalent B-Spline curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_convert_bspline_to_nurbs

```c
qaws_status qaws_curve_convert_bspline_to_nurbs(
    qaws_curve** out_nurbs,
    qaws_curve const* curve);
```

Converts a B-Spline curve to a NURBS curve (with all weights equal to 1.0).

**Parameters:**
- `curve` -- A B-Spline curve.
- `out_nurbs` -- Receives the equivalent NURBS curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_elevate_degree

```c
qaws_status qaws_curve_elevate_degree(
    qaws_curve** out_elevated,
    qaws_curve const* curve);
```

Elevates the degree of a curve by one (degree n to degree n+1).

**Parameters:**
- `curve` -- The input curve.
- `out_elevated` -- Receives the degree-elevated curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_reduce_degree

```c
qaws_status qaws_curve_reduce_degree(
    qaws_curve** out_reduced,
    qaws_curve const* curve);
```

Reduces the degree of a Bezier curve by one (degree n to degree n-1).

**Parameters:**
- `curve` -- A Bezier curve.
- `out_reduced` -- Receives the degree-reduced curve.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- Bezier only. Degree reduction is an approximation; the result may not exactly match the original curve.
