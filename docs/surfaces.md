# Surface Types

Each surface type is created from a descriptor by a `qaws_surface_create_*`
function. After creation every type shares the same evaluation, sampling,
normal, inspection and destruction calls. Names below are the unsuffixed
forms ([precision.md](precision.md)).

| Kind (`QAWS_SURFACE_KIND_*`) | Header | Constructor | Inputs |
|---|---|---|---|
| `BEZIER` | `qaws_surface_bezier.h` | `qaws_surface_create_bezier` | control grid (copied) |
| `BSPLINE` | `qaws_surface_bspline.h` | `qaws_surface_create_bspline` | control grid, knots (copied) |
| `NURBS` | `qaws_surface_nurbs.h` | `qaws_surface_create_nurbs` | control grid, weights, knots (copied) |
| `BILINEAR` | `qaws_surface_bilinear.h` | `qaws_surface_create_bilinear` | four corners |
| `BIQUADRATIC` | `qaws_surface_biquadratic.h` | `qaws_surface_create_biquadratic` | 3 x 3 control points |
| `SWEPT` | `qaws_surface_swept.h` | `qaws_surface_create_swept` | path and profile curves (borrowed) |
| `RULED` | `qaws_surface_ruled.h` | `qaws_surface_create_ruled` | two curves (borrowed) |
| `COONS` | `qaws_surface_coons.h` | `qaws_surface_create_coons` | four boundary curves (borrowed) |
| `EXTRUSION` | `qaws_surface_extrusion.h` | `qaws_surface_create_extrusion` | profile curve (borrowed), direction |
| `REVOLUTION` | `qaws_surface_revolution.h` | `qaws_surface_create_revolution` | 2D profile (borrowed), axis, angle |
| `PIPE` | `qaws_surface_pipe.h` | `qaws_surface_create_pipe` | path curve (borrowed), radii |
| `LOFT` | `qaws_surface_loft.h` | `qaws_surface_create_loft` | section curves (borrowed) |
| `GORDON` | `qaws_surface_gordon.h` | `qaws_surface_create_gordon` | curve network (borrowed) |
| `OFFSET` | `qaws_surface_offset.h` | `qaws_surface_create_offset` | base surface (borrowed), distance |
| `TRIMMED` | `qaws_surface_trim.h` | `qaws_surface_create_trimmed` | base surface and trim curves (borrowed) |
| `TSPLINE` | `qaws_surface_tspline.h` | `qaws_surface_create_tspline` | control points with local knot vectors |
| `SUBDIV` | `qaws_surface_subdiv.h` | `qaws_surface_create_subdiv` | control mesh (copied) |
| `VIEW` | `qaws_view.h` | `qaws_surface_create_borrowed` | host arrays (borrowed), see [views_and_batches.md](views_and_batches.md) |

`QAWS_SURFACE_KIND_FILLET` is declared but no public constructor makes it.
Every constructor except the borrowed one has an `_ex` form taking a
`qaws_allocator` (`qaws_alloc.h`, [memory.md](memory.md)). A surface that
borrows curves or another surface needs them alive, unchanged, until it is
destroyed.

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

### Batch evaluation

`qaws_surface_evaluate_batch` (`qaws_eval.h`) takes parallel arrays of `u`
and `v`; see [evaluation.md](evaluation.md). Bounds, area, curvatures,
closest points, enclosures, tessellation and curve projection are in
`qaws_inspect.h` ([inspection.md](inspection.md),
[sampling.md](sampling.md)).

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

---

## Bilinear and biquadratic

**Headers:** `qaws_surface_bilinear.h`, `qaws_surface_biquadratic.h`

- Bilinear: `S(u,v) = (1-u)(1-v) P00 + u(1-v) P10 + (1-u)v P01 + uv P11`
  from the corners `p00`, `p10`, `p01`, `p11` (u left to right, v bottom to
  top).
- Biquadratic: degree (2, 2) Bernstein patch from `control_points[9]`, row
  major: `cp[0..2]` the v = 0 row, `cp[3..5]` v = 1/2, `cp[6..8]` v = 1; u
  increases within a row.

---

## Coons

**Header:** `qaws_surface_coons.h`

Bilinearly blended patch from four 3D boundary curves:
`S = Lc + Ld - B`, `Lc` ruled between `c0(u)` (v = 0) and `c1(u)` (v = 1),
`Ld` ruled between `d0(v)` (u = 0) and `d1(v)` (u = 1), `B` the bilinear
patch of the corners. The corners must agree: `c0(0) = d0(0)`,
`c0(1) = d1(0)`, `c1(0) = d0(1)`, `c1(1) = d1(1)`. Curves are borrowed.

---

## Extrusion and revolution

**Headers:** `qaws_surface_extrusion.h`, `qaws_surface_revolution.h`

- Extrusion: `S(u,v) = profile(u) + v * direction`, u the profile's domain, v
  from 0 to `length` (0: `|direction|` as unit). A 2D profile is lifted to
  z = 0.
- Revolution: a 2D profile (x the radius, y the height) turned about
  `axis_origin` and `axis_direction` (normalized) by `angle` radians (0: a
  full turn); u in [0, 1] maps to the angle, v follows the profile.

Both borrow the profile curve.

---

## Pipe

**Header:** `qaws_surface_pipe.h`

`S(u,v) = path(u) + radius_x cos(2 pi v) N(u) + radius_y sin(2 pi v) B(u)`,
u in [0, 1] along the 3D path, v in [0, 1] around the section; `N`, `B` a
rotation-minimizing frame by double reflection. `radius_y` 0 means
`radius_x` (a circular pipe). The path is borrowed.

---

## Loft and Gordon

**Headers:** `qaws_surface_loft.h`, `qaws_surface_gordon.h`

- Loft: `section_count` (at least 2) 3D section curves; at
  `v = v_parameters[i]` (NULL: uniform) the surface passes through section
  `i`; between sections cubic B-spline blending, u along each section.
- Gordon: a network of `u_curve_count` u-curves sitting at `v_params` and
  `v_curve_count` v-curves sitting at `u_params` (at least 2 each);
  `S = L_u + L_v - T` (two lofts less their tensor product through the
  intersection points). Compatibility: `u_curve[i](v_params[j])` must equal
  `v_curve[j](u_params[i])`.

Curves are borrowed. Evaluation keeps per-section scratch on the stack up to
32 sections, beyond that from the surface's allocator.

---

## Offset

**Header:** `qaws_surface_offset.h`

`S_off(u,v) = S_base(u,v) + distance * N_base(u,v)`; positive distance along
the base normal. The base surface is borrowed.

---

## Trimmed

**Header:** `qaws_surface_trim.h`

A base surface restricted by trim loops in (u, v): each `qaws_trim_loop` is a
closed chain of 2D curves with `is_outer` 1 (kept region) or 0 (hole).
Without an outer loop the whole [0,1] x [0,1] domain is the outer boundary.
`qaws_surface_trim_contains(&inside, trimmed, u, v)` tells whether (u, v) is
kept. The base surface and the trim curves are borrowed. Trimmed pieces are
also what [shatter](shatter.md) produces.

---

## T-spline

**Header:** `qaws_surface_tspline.h`

Locally refinable surface: each `qaws_tspline_control_point` has a position,
a weight (1 for non-rational) and its own local knot vectors of length
`2 * u_degree + 2` and `2 * v_degree + 2`. The basis of point i is the product
of the two univariate B-spline bases of its local knots;
`S = sum w_i P_i B_i / sum w_i B_i` over the given `u_range` x `v_range`.
The local knot vectors are given explicitly, not inferred from a T-mesh.

---

## Subdivision

**Header:** `qaws_surface_subdiv.h`

Limit surface of a control mesh (`vertices`, `face_indices`, `face_sizes` of
3 or 4), by `QAWS_SUBDIV_CATMULL_CLARK` (quads; mixed meshes use it, triangles
as degenerate quads) or `QAWS_SUBDIV_LOOP` (triangles). The mesh is copied,
subdivided `subdivision_level` times (0: 3) and stored on a regular grid;
evaluation over [0,1] x [0,1] interpolates it bilinearly.
`qaws_surface_subdiv_get_mesh` returns the subdivided triangle mesh (CCW
triangles, smooth per-vertex normals) as pointers into the surface, valid
until it is destroyed.
