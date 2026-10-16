# Surface Types

Qaws supports five surface types. Each is created through a descriptor struct and a `qaws_surface_create_*` function. After creation, all types share the same evaluation, sampling, normal-query, and destruction APIs.

---

## Common API

**Headers:** `qaws_surface.h`, `qaws_surface_types.h`

All surfaces are immutable after creation. Multiple threads may concurrently evaluate or sample the same surface.

### Evaluation flags

Request specific outputs by combining flags:

| Flag | Value | Output field |
|---|---|---|
| `QAWS_SURFACE_EVAL_POSITION` | `1 << 0` | `position` |
| `QAWS_SURFACE_EVAL_DU` | `1 << 1` | `du` (partial derivative w.r.t. u) |
| `QAWS_SURFACE_EVAL_DV` | `1 << 2` | `dv` (partial derivative w.r.t. v) |
| `QAWS_SURFACE_EVAL_NORMAL` | `1 << 3` | `normal` (unit normal, du x dv) |
| `QAWS_SURFACE_EVAL_DUU` | `1 << 4` | `duu` (second partial d^2/du^2) |
| `QAWS_SURFACE_EVAL_DUV` | `1 << 5` | `duv` (mixed partial d^2/dudv) |
| `QAWS_SURFACE_EVAL_DVV` | `1 << 6` | `dvv` (second partial d^2/dv^2) |

### Evaluate

```c
qaws_surface_eval_result result;
qaws_surface_evaluate(&result, surface, u,
    v,
    QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
/* result.position and result.normal are now populated.
   result.valid_flags indicates which fields were computed. */
```

### Get normal

Convenience function that returns only the unit surface normal at (u, v):

```c
qaws_vec3 n;
qaws_surface_get_normal(&n, surface, u, v);
```

### Sample a grid

Produces a u_count x v_count grid of positions in row-major order (v varies fastest). The caller allocates the output array.

```c
unsigned int u_count = 32, v_count = 32;
qaws_vec3 positions[32 * 32];
unsigned int count;
qaws_surface_sample(positions, &count, surface,
    u_count, v_count, 32 * 32);
```

### Inspection

```c
qaws_surface_kind kind = qaws_surface_get_kind(surface);
unsigned int u_deg = qaws_surface_get_u_degree(surface);
unsigned int v_deg = qaws_surface_get_v_degree(surface);
qaws_range u_range = qaws_surface_get_u_range(surface);
qaws_range v_range = qaws_surface_get_v_range(surface);
int rational = qaws_surface_is_rational(surface);
```

### Destroy

```c
qaws_surface_destroy(surface);
```

---

## Bezier

**Header:** `qaws_surface_bezier.h`

Tensor-product Bezier patch defined by a (u_degree+1) x (v_degree+1) grid of 3D control points. The patch interpolates its four corner control points and is influenced by the interior ones.

### Descriptor fields

| Field | Type | Description |
|---|---|---|
| `u_degree` | `unsigned int` | Degree in the u direction |
| `v_degree` | `unsigned int` | Degree in the v direction |
| `control_points` | `qaws_vec3 const*` | Row-major grid of control points |
| `u_point_count` | `unsigned int` | Must be `u_degree + 1` |
| `v_point_count` | `unsigned int` | Must be `v_degree + 1` |

### Key properties

- **Degree:** arbitrary per direction (commonly 3x3 for bicubic)
- **Patches:** 1
- **Rational:** no
- **Domain:** u in [0,1], v in [0,1]
- **Control points:** `(u_degree+1) * (v_degree+1)`

**When to use:** Single-patch surfaces where you need direct control over shape via a control-point grid. Common in rendering, font outlines, and simple freeform modeling.

```c
/* Bilinear patch (degree 1x1) */
qaws_vec3 pts[] = {
    {0,0,0}, {1,0,0},
    {0,1,0}, {1,1,1}
};
qaws_surface_bezier_desc desc = {
    .u_degree = 1,
    .v_degree = 1,
    .control_points = pts,
    .u_point_count = 2,
    .v_point_count = 2
};
qaws_surface* s = NULL;
qaws_surface_create_bezier(&s, &desc);
```

---

## B-Spline

**Header:** `qaws_surface_bspline.h`

Tensor-product B-spline surface with independent knot vectors in u and v. Provides local control: moving one control point affects only a bounded region of the surface.

### Descriptor fields

| Field | Type | Description |
|---|---|---|
| `u_degree` | `unsigned int` | Degree in the u direction |
| `v_degree` | `unsigned int` | Degree in the v direction |
| `control_points` | `qaws_vec3 const*` | Row-major grid (u_point_count x v_point_count) |
| `u_point_count` | `unsigned int` | Number of control points in u |
| `v_point_count` | `unsigned int` | Number of control points in v |
| `u_knots` | `qaws_scalar const*` | Knot vector in u (NULL for auto uniform clamped) |
| `u_knot_count` | `unsigned int` | Length of u knot vector (0 for auto) |
| `v_knots` | `qaws_scalar const*` | Knot vector in v (NULL for auto uniform clamped) |
| `v_knot_count` | `unsigned int` | Length of v knot vector (0 for auto) |

### Key properties

- **Degree:** arbitrary per direction (commonly 3x3)
- **Patches:** determined by knot vectors and degrees
- **Rational:** no
- **Continuity:** up to C(degree-1) at internal knots in each direction
- **Local control:** yes

**When to use:** CAD surfaces, complex multi-patch shapes, situations where local control and custom knot spacing matter.

```c
qaws_vec3 pts[] = {
    {0,0,0}, {1,0,0}, {2,0,0}, {3,0,0},
    {0,1,0}, {1,1,1}, {2,1,1}, {3,1,0},
    {0,2,0}, {1,2,0}, {2,2,0}, {3,2,0}
};
qaws_surface_bspline_desc desc = {
    .u_degree = 2,
    .v_degree = 3,
    .control_points = pts,
    .u_point_count = 3,
    .v_point_count = 4,
    .u_knots = NULL,   /* auto uniform clamped */
    .u_knot_count = 0,
    .v_knots = NULL,
    .v_knot_count = 0
};
qaws_surface* s = NULL;
qaws_surface_create_bspline(&s, &desc);
```

---

## NURBS

**Header:** `qaws_surface_nurbs.h`

Non-Uniform Rational B-Spline surface. Extends B-spline surfaces with a per-control-point weight grid, enabling exact representation of spheres, cylinders, cones, tori, and other quadric surfaces.

### Descriptor fields

| Field | Type | Description |
|---|---|---|
| `u_degree` | `unsigned int` | Degree in the u direction |
| `v_degree` | `unsigned int` | Degree in the v direction |
| `control_points` | `qaws_vec3 const*` | Row-major grid (u_point_count x v_point_count) |
| `u_point_count` | `unsigned int` | Number of control points in u |
| `v_point_count` | `unsigned int` | Number of control points in v |
| `weights` | `qaws_scalar const*` | Weight grid, same layout as control points |
| `u_knots` | `qaws_scalar const*` | Knot vector in u (NULL for auto uniform clamped) |
| `u_knot_count` | `unsigned int` | Length of u knot vector (0 for auto) |
| `v_knots` | `qaws_scalar const*` | Knot vector in v (NULL for auto uniform clamped) |
| `v_knot_count` | `unsigned int` | Length of v knot vector (0 for auto) |

### Key properties

- **Degree:** arbitrary per direction (commonly 2 or 3)
- **Patches:** determined by knot vectors and degrees
- **Rational:** yes
- **Continuity:** up to C(degree-1) at internal knots in each direction
- **Local control:** yes

**When to use:** Exact quadric surfaces (spheres, cylinders, cones), CAD interchange formats (STEP, IGES), any situation where B-spline surfaces are insufficient due to the need for rational weighting.

```c
/* Rational biquadratic patch */
qaws_vec3 pts[] = {
    {1,0,0}, {1,1,0}, {0,1,0},
    {1,0,1}, {1,1,1}, {0,1,1},
    {1,0,2}, {1,1,2}, {0,1,2}
};
qaws_scalar weights[] = {
    1.0f, 0.707f, 1.0f,
    0.707f, 0.5f,  0.707f,
    1.0f, 0.707f, 1.0f
};
qaws_scalar u_knots[] = { 0,0,0, 1,1,1 };
qaws_scalar v_knots[] = { 0,0,0, 1,1,1 };
qaws_surface_nurbs_desc desc = {
    .u_degree = 2,
    .v_degree = 2,
    .control_points = pts,
    .u_point_count = 3,
    .v_point_count = 3,
    .weights = weights,
    .u_knots = u_knots,
    .u_knot_count = 6,
    .v_knots = v_knots,
    .v_knot_count = 6
};
qaws_surface* s = NULL;
qaws_surface_create_nurbs(&s, &desc);
```

---

## Swept

**Header:** `qaws_surface_swept.h`

A surface formed by sweeping a 2D cross-section profile along a 3D path curve. At each point along the path, the profile is oriented using the Frenet frame (normal and binormal):

```
S(u,v) = path(u) + profile_x(v) * N(u) + profile_y(v) * B(u)
```

where N and B are the Frenet normal and binormal of the path at parameter u.

### Descriptor fields

| Field | Type | Description |
|---|---|---|
| `path` | `qaws_curve const*` | 3D path curve (borrowed, not owned) |
| `profile` | `qaws_curve const*` | 2D cross-section curve (borrowed, not owned) |
| `scale` | `qaws_scalar` | Uniform scale of the cross-section (0 defaults to 1.0) |

### Key properties

- **u parameter:** follows the path curve domain
- **v parameter:** follows the profile curve domain
- **Rational:** depends on underlying curves
- **Ownership:** the surface borrows the path and profile curves; the caller must keep both curves alive for the lifetime of the surface

**When to use:** Tubes, pipes, extruded shapes, any geometry naturally described as a profile moving along a spine curve.

```c
/* Assume path_3d is a 3D Bezier spine curve,
   and circle_2d is a 2D arc describing a circular cross-section. */
qaws_surface_swept_desc desc = {
    .path = path_3d,
    .profile = circle_2d,
    .scale = 1.0f
};
qaws_surface* s = NULL;
qaws_surface_create_swept(&s, &desc);
/* path_3d and circle_2d must remain valid while s is in use. */
```

---

## Ruled

**Header:** `qaws_surface_ruled.h`

A surface formed by linearly interpolating between two 3D boundary curves:

```
S(u,v) = (1 - v) * curve_a(u) + v * curve_b(u)
```

At each u parameter value, the surface traces a straight line from curve_a to curve_b as v goes from 0 to 1.

### Descriptor fields

| Field | Type | Description |
|---|---|---|
| `curve_a` | `qaws_curve const*` | First boundary curve, 3D (borrowed, not owned) |
| `curve_b` | `qaws_curve const*` | Second boundary curve, 3D (borrowed, not owned) |

### Key properties

- **u parameter:** mapped to each curve's own domain [0,1]
- **v parameter:** blends linearly from curve_a (v=0) to curve_b (v=1)
- **Rational:** depends on underlying curves
- **Ownership:** the surface borrows both curves; the caller must keep them alive for the lifetime of the surface

**When to use:** Transition surfaces between two boundary curves, wing skins, lofted shapes with straight-line cross-sections, any geometry where rulings (straight lines) connect corresponding points on two curves.

```c
/* Assume top_curve and bottom_curve are 3D curves. */
qaws_surface_ruled_desc desc = {
    .curve_a = bottom_curve,
    .curve_b = top_curve
};
qaws_surface* s = NULL;
qaws_surface_create_ruled(&s, &desc);
/* bottom_curve and top_curve must remain valid while s is in use. */
```
