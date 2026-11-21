# qaws_clothoid.h: clothoids

A 2D clothoid (Euler or Cornu spiral): curvature varies linearly with arc length. Each name exists as `name_f32` and
`name_f64` ([precision.md](precision.md)).

## qaws_clothoid_desc

```c
typedef struct qaws_clothoid_desc {
    qaws_scalar origin_x;
    qaws_scalar origin_y;
    qaws_scalar start_angle;        /* radians */
    qaws_scalar start_curvature;    /* kappa_0 */
    qaws_scalar end_curvature;      /* kappa_1 */
    qaws_scalar length;             /* arc length L > 0 */
} qaws_clothoid_desc;
```

| Member | Meaning |
|---|---|
| `origin_x`, `origin_y` | the start point |
| `start_angle` | the initial heading, radians |
| `start_curvature` | curvature at s = 0 (kappa_0) |
| `end_curvature` | curvature at s = L (kappa_1) |
| `length` | total arc length L, > 0 |

The curve:

```
kappa(s) = kappa_0 + (kappa_1 - kappa_0) * s / L
theta(s) = start_angle + kappa_0 * s + (kappa_1 - kappa_0) * s^2 / (2 * L)
x(s)     = origin_x + integral_0^s cos(theta(u)) du
y(s)     = origin_y + integral_0^s sin(theta(u)) du
```

## qaws_curve_create_clothoid

```c
qaws_status qaws_curve_create_clothoid(
    qaws_curve** out_curve,
    qaws_clothoid_desc const* desc);
```

Creates the curve: always 2D, one span, parameter = arc length, range [0, L]. The degree reported is 3. The descriptor
is copied. Destroy the result with `qaws_curve_destroy` ([curve.md](curve.md)).

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `desc` or `out_curve` NULL; `length <= 0` |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL first. The allocator-taking form `qaws_curve_create_clothoid_ex` is declared in
`qaws_alloc.h` ([alloc.md](alloc.md)).
