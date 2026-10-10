# API Reference: qaws_curve.h

Base curve lifecycle.

---

## qaws_curve_destroy

```c
void qaws_curve_destroy(qaws_curve* curve);
```

Destroys a curve object and frees all associated memory.

**Parameters:**
- `curve` -- Curve to destroy. May be NULL (no-op).

**Semantics:**
- Must be called for every curve created with `qaws_curve_create_*`.
- After destruction, the pointer must not be used.
- Any traversal objects referencing this curve become invalid.

---

## qaws_curve_reverse

```c
qaws_status qaws_curve_reverse(
    qaws_curve const* curve,
    qaws_curve** out_reversed);
```

Creates a new curve that traces the same path as the input but in the opposite direction.

**Parameters:**
- `curve` -- The source curve.
- `out_reversed` -- Receives the newly allocated reversed curve.

**Returns:** `QAWS_STATUS_OK` on success.

**Semantics:**
- The returned curve is a new allocation and must be destroyed with `qaws_curve_destroy`.
- Evaluating the reversed curve at parameter `t` yields the same position as evaluating the original at `t_max - t + t_min`.
- Every kind can be reversed. Bezier, Hermite, Catmull-Rom (closed ones keep
  their start point), B-spline, NURBS, trajectory, quadratic Yuksel, rational
  Bezier, polynomial, arc, clothoid and composite curves come back as their own
  kind. Subdivision, reparameterized and arc-mode Yuksel curves come back as
  the reversed `qaws_curve_extract` of their whole domain.
