# API Reference: qaws_alloc.h

Custom allocator support for curve and surface creation.

---

## Overview

The `_ex` variants of creation functions behave identically to their non-`_ex` counterparts but route all internal allocations through the provided `qaws_allocator`. Pass NULL for `allocator` to use the default `malloc`/`free`.

Curves and surfaces created with a custom allocator MUST be destroyed with `qaws_curve_destroy` or `qaws_surface_destroy` (which read the stored allocator). The allocator must remain valid until the object is destroyed.

---

## qaws_curve_create_bezier_ex

```c
qaws_status qaws_curve_create_bezier_ex(
    qaws_bezier_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_bezier`.

**Parameters:**
- `desc` -- Bezier curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_hermite_ex

```c
qaws_status qaws_curve_create_hermite_ex(
    qaws_hermite_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_hermite`.

**Parameters:**
- `desc` -- Hermite curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_catmull_rom_ex

```c
qaws_status qaws_curve_create_catmull_rom_ex(
    qaws_catmull_rom_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_catmull_rom`.

**Parameters:**
- `desc` -- Catmull-Rom curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_bspline_ex

```c
qaws_status qaws_curve_create_bspline_ex(
    qaws_bspline_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_bspline`.

**Parameters:**
- `desc` -- B-Spline curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_nurbs_ex

```c
qaws_status qaws_curve_create_nurbs_ex(
    qaws_nurbs_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_nurbs`.

**Parameters:**
- `desc` -- NURBS curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_trajectory_ex

```c
qaws_status qaws_curve_create_trajectory_ex(
    qaws_trajectory_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_trajectory`.

**Parameters:**
- `desc` -- Trajectory curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## The other curve families

```c
qaws_status qaws_curve_create_rational_bezier_ex(qaws_rational_bezier_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_yuksel_ex(qaws_yuksel_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_arc_ex(qaws_arc_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_polynomial_ex(qaws_polynomial_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_clothoid_ex(qaws_clothoid_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_composite_ex(qaws_composite_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_subdivision_ex(qaws_subdivision_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
```

Each is the allocator-aware variant of the constructor without `_ex`,
which calls it with NULL. Every allocation the constructor makes goes
through the allocator: the curve, its stored data, and its scratch (the
subdivision's refinement levels). A composite owns its segments only on
success; the segments keep the allocators they were made with.

**Checked by test 146.** Every curve and surface constructor here runs with
an allocator that fails at its k-th allocation, for every k. Each failure
returns `QAWS_STATUS_ALLOCATION_FAILURE` with no object and nothing left
allocated. A success, once destroyed, leaves nothing either.

---

## qaws_surface_create_bezier_ex

```c
qaws_status qaws_surface_create_bezier_ex(
    qaws_surface_bezier_desc const* desc,
    qaws_allocator const* allocator,
    qaws_surface** out_surface);
```

Allocator-aware variant of `qaws_surface_create_bezier`.

**Parameters:**
- `desc` -- Bezier surface descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_surface` -- Receives the created surface.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_surface_create_bspline_ex

```c
qaws_status qaws_surface_create_bspline_ex(
    qaws_surface_bspline_desc const* desc,
    qaws_allocator const* allocator,
    qaws_surface** out_surface);
```

Allocator-aware variant of `qaws_surface_create_bspline`.

**Parameters:**
- `desc` -- B-Spline surface descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_surface` -- Receives the created surface.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_surface_create_nurbs_ex

```c
qaws_status qaws_surface_create_nurbs_ex(
    qaws_surface_nurbs_desc const* desc,
    qaws_allocator const* allocator,
    qaws_surface** out_surface);
```

Allocator-aware variant of `qaws_surface_create_nurbs`.

**Parameters:**
- `desc` -- NURBS surface descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_surface` -- Receives the created surface.

**Returns:** `QAWS_STATUS_OK` on success.
