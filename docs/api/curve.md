# qaws_curve.h: curve lifetime

Destroying a curve and making its reversed copy. Curves are created by the family headers (see
[types.md](types.md) for the list of kinds). Each name exists as `name_f32` and `name_f64`
([precision.md](precision.md)).

## qaws_curve_destroy

```c
void qaws_curve_destroy(qaws_curve* curve);
```

Destroys a curve and frees its memory with the allocator it was created with. NULL is a no-op.

- Call it once for every curve a `qaws_curve_create_*` (or other curve-returning) call gave you. Never call it on an
  inline curve ([inline.md](inline.md)).
- A composite curve destroys the segment curves it owns ([composite.md](composite.md)).
- The curve must not be in use by another thread, and traversals created on it must not be used afterwards.

## qaws_curve_reverse

```c
qaws_status qaws_curve_reverse(
    qaws_curve** out_reversed,
    qaws_curve const* curve);
```

Creates a new curve tracing the same path in the opposite direction: evaluating it at `t` gives the position of the
original at `t_min + t_max - t`.

- Bezier, Hermite, Catmull-Rom (a closed one keeps its start point), B-spline, NURBS, trajectory, quadratic-mode
  Yuksel, rational Bezier, polynomial, arc, clothoid and composite curves come back as their own kind. Arc-mode Yuksel
  (circular, elliptical, hybrid), subdivision, reparameterized and view curves come back as the reversed
  `qaws_curve_extract` of their whole domain ([operations.md](operations.md)).
- The result is a new allocation; destroy it with `qaws_curve_destroy`. The source is not changed.

Returns `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` for a NULL argument; `QAWS_STATUS_UNSUPPORTED_OPERATION` for an
unknown kind; `QAWS_STATUS_ALLOCATION_FAILURE`; or the status of the creation it delegates to. `*out_reversed` is set
to NULL first.
