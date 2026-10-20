# API Reference: qaws_clothoid.h

Clothoid (Euler spiral / Cornu spiral) curve creation.

---

## qaws_clothoid_desc

```c
typedef struct qaws_clothoid_desc {
    qaws_scalar origin_x;
    qaws_scalar origin_y;
    qaws_scalar start_angle;
    qaws_scalar start_curvature;
    qaws_scalar end_curvature;
    qaws_scalar length;
} qaws_clothoid_desc;
```

| Field | Description |
|---|---|
| `origin_x` | Starting point x-coordinate |
| `origin_y` | Starting point y-coordinate |
| `start_angle` | Initial heading in radians |
| `start_curvature` | Curvature at s=0 (kappa_0) |
| `end_curvature` | Curvature at s=L (kappa_1) |
| `length` | Total arc length L. Must be > 0. |

The curvature varies linearly with arc length:

```
kappa(s) = kappa_0 + (kappa_1 - kappa_0) * s / L
```

Position is computed by integrating:

```
x(s) = origin_x + integral_0^s cos(theta(u)) du
y(s) = origin_y + integral_0^s sin(theta(u)) du
```

where `theta(s) = start_angle + kappa_0*s + (kappa_1 - kappa_0)*s^2 / (2*L)`.

---

## qaws_curve_create_clothoid

```c
qaws_status qaws_curve_create_clothoid(
    qaws_curve** out_curve,
    qaws_clothoid_desc const* desc);
```

Creates a clothoid (Euler spiral) curve.

**Parameters:**
- `desc` -- Clothoid descriptor. Must not be NULL.
- `out_curve` -- Receives the created curve. Set to NULL on error.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- 2D only. 3D clothoids are not supported.
- The descriptor data is copied. The caller may free it after creation.
