# API Reference: qaws_inline.h

Stack-allocated (inline) curve initialization without heap allocation.

---

## Overview

These functions initialize a `qaws_curve_inline` struct without any heap allocation. The resulting curve is usable with all evaluation and inspection functions via `&inline_curve.curve`.

**Constraints:**
- Single span only (degree+1 control points for Bezier).
- Maximum degree 7 (8 control points x 3 dimensions = 24 scalars).
- 2D or 3D.
- Bezier and Polynomial families only.

The `qaws_curve_inline` struct is self-contained. It can be copied with `memcpy` or assigned. Do NOT call `qaws_curve_destroy` on it.

**Example:**

```c
qaws_curve_inline ic;
qaws_vec3 pts[] = {{0,0,0}, {1,2,3}, {4,5,6}, {7,8,9}};
qaws_bezier_desc desc = { QAWS_DIMENSION_3D, 3, pts, 4 };
qaws_curve_init_bezier_inline(&ic, &desc);
// Use &ic.curve with qaws_curve_evaluate_3d, etc.
```

---

## qaws_curve_init_bezier_inline

```c
qaws_status qaws_curve_init_bezier_inline(
    qaws_curve_inline* out,
    qaws_bezier_desc const* desc);
```

Initializes an inline Bezier curve on the stack.

**Parameters:**
- `desc` -- Bezier curve descriptor. Must describe a single span (control_point_count == degree + 1) with degree <= 7.
- `out` -- Pointer to a `qaws_curve_inline` struct to initialize.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_init_polynomial_inline

```c
qaws_status qaws_curve_init_polynomial_inline(
    qaws_curve_inline* out,
    qaws_polynomial_desc const* desc);
```

Initializes an inline Polynomial curve on the stack.

**Parameters:**
- `desc` -- Polynomial curve descriptor. Must have degree <= 7.
- `out` -- Pointer to a `qaws_curve_inline` struct to initialize.

**Returns:** `QAWS_STATUS_OK` on success.
