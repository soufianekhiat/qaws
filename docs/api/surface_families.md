# API Reference: surface families

Constructors for the tensor-product and constructed surface families. Every
constructor has the form

```c
qaws_status qaws_surface_create_<family>(qaws_surface** out_surface, qaws_surface_<family>_desc const* desc);
```

and returns a new `qaws_surface` that the caller frees with
`qaws_surface_destroy`. Evaluation, sampling and properties are the generic
functions of [surface.md](surface.md). Errors follow [status.md](status.md):
`QAWS_STATUS_INVALID_ARGUMENT` for a NULL or inconsistent descriptor,
`QAWS_STATUS_INVALID_DIMENSION` for a curve of the wrong dimension (for
example a 3D profile given to a revolution), `QAWS_STATUS_ALLOCATION_FAILURE`
when memory runs out.

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)).

**Ownership.** Control-point families (Bezier, B-spline, NURBS, bilinear,
biquadratic) copy their data. Families built from curves or surfaces
(swept, ruled, Coons, extrusion, revolution, pipe, loft, Gordon, offset)
borrow them: the curves and the base surface must outlive the new surface.

**Grid layout.** Control-point grids are row-major: row i holds the points
P[i][0..v_point_count-1], with i the u index.

| Header | Kind |
|---|---|
| qaws_surface_bezier.h | `QAWS_SURFACE_KIND_BEZIER` |
| qaws_surface_bspline.h | `QAWS_SURFACE_KIND_BSPLINE` |
| qaws_surface_nurbs.h | `QAWS_SURFACE_KIND_NURBS` |
| qaws_surface_bilinear.h | `QAWS_SURFACE_KIND_BILINEAR` |
| qaws_surface_biquadratic.h | `QAWS_SURFACE_KIND_BIQUADRATIC` |
| qaws_surface_swept.h | `QAWS_SURFACE_KIND_SWEPT` |
| qaws_surface_ruled.h | `QAWS_SURFACE_KIND_RULED` |
| qaws_surface_coons.h | `QAWS_SURFACE_KIND_COONS` |
| qaws_surface_extrusion.h | `QAWS_SURFACE_KIND_EXTRUSION` |
| qaws_surface_revolution.h | `QAWS_SURFACE_KIND_REVOLUTION` |
| qaws_surface_pipe.h | `QAWS_SURFACE_KIND_PIPE` |
| qaws_surface_loft.h | `QAWS_SURFACE_KIND_LOFT` |
| qaws_surface_gordon.h | `QAWS_SURFACE_KIND_GORDON` |
| qaws_surface_offset.h | `QAWS_SURFACE_KIND_OFFSET` |

---

## Control-point surfaces

### Bezier (qaws_surface_bezier.h)

```c
typedef struct qaws_surface_bezier_desc {
	unsigned int u_degree;
	unsigned int v_degree;
	qaws_vec3 const* control_points;   /* (u_degree+1) x (v_degree+1), row-major */
	unsigned int u_point_count;        /* must be u_degree + 1 */
	unsigned int v_point_count;        /* must be v_degree + 1 */
} qaws_surface_bezier_desc;

qaws_status qaws_surface_create_bezier(qaws_surface** out_surface, qaws_surface_bezier_desc const* desc);
```

A tensor-product Bezier patch.

### B-spline (qaws_surface_bspline.h)

```c
typedef struct qaws_surface_bspline_desc {
	unsigned int u_degree;
	unsigned int v_degree;
	qaws_vec3 const* control_points;   /* u_point_count x v_point_count, row-major */
	unsigned int u_point_count;
	unsigned int v_point_count;
	qaws_scalar const* u_knots;
	unsigned int u_knot_count;         /* 0 for auto uniform clamped */
	qaws_scalar const* v_knots;
	unsigned int v_knot_count;         /* 0 for auto uniform clamped */
} qaws_surface_bspline_desc;

qaws_status qaws_surface_create_bspline(qaws_surface** out_surface, qaws_surface_bspline_desc const* desc);
```

A tensor-product B-spline. A knot count of 0 in a direction makes a uniform
clamped knot vector there.

### NURBS (qaws_surface_nurbs.h)

```c
typedef struct qaws_surface_nurbs_desc {
	unsigned int u_degree;
	unsigned int v_degree;
	qaws_vec3 const* control_points;   /* u_point_count x v_point_count */
	unsigned int u_point_count;
	unsigned int v_point_count;
	qaws_scalar const* weights;        /* u_point_count x v_point_count, same layout as the points */
	qaws_scalar const* u_knots;
	unsigned int u_knot_count;         /* 0 for auto uniform clamped */
	qaws_scalar const* v_knots;
	unsigned int v_knot_count;         /* 0 for auto uniform clamped */
} qaws_surface_nurbs_desc;

qaws_status qaws_surface_create_nurbs(qaws_surface** out_surface, qaws_surface_nurbs_desc const* desc);
```

A rational tensor-product B-spline; `qaws_surface_is_rational` returns 1.

### Bilinear (qaws_surface_bilinear.h)

```c
typedef struct qaws_surface_bilinear_desc {
	qaws_vec3 p00;   /* corner at (u=0, v=0) */
	qaws_vec3 p10;   /* corner at (u=1, v=0) */
	qaws_vec3 p01;   /* corner at (u=0, v=1) */
	qaws_vec3 p11;   /* corner at (u=1, v=1) */
} qaws_surface_bilinear_desc;

qaws_status qaws_surface_create_bilinear(qaws_surface** out_surface, qaws_surface_bilinear_desc const* desc);
```

A degree (1, 1) patch through four corners:
`S(u,v) = (1-u)(1-v) P00 + u(1-v) P10 + (1-u)v P01 + uv P11`.

```
P01 --- P11
 |       |
P00 --- P10      u left to right, v bottom to top
```

### Biquadratic (qaws_surface_biquadratic.h)

```c
typedef struct qaws_surface_biquadratic_desc {
	qaws_vec3 control_points[9];   /* 3x3 grid, row-major */
} qaws_surface_biquadratic_desc;

qaws_status qaws_surface_create_biquadratic(qaws_surface** out_surface,
	qaws_surface_biquadratic_desc const* desc);
```

A degree (2, 2) Bernstein patch, B_i^2(t) = {(1-t)^2, 2t(1-t), t^2}. Note
the layout here is by v rows: `cp[0..2]` is the v = 0 row, `cp[3..5]` the
v = 1/2 row, `cp[6..8]` the v = 1 row; within a row u increases
(u = 0, 1/2, 1).

---

## Surfaces from curves

### Swept (qaws_surface_swept.h)

```c
typedef struct qaws_surface_swept_desc {
	qaws_curve const* path;      /* 3D path curve (borrowed) */
	qaws_curve const* profile;   /* 2D cross-section curve (borrowed) */
	qaws_scalar scale;           /* uniform scale of the cross-section; 0 = 1 */
} qaws_surface_swept_desc;

qaws_status qaws_surface_create_swept(qaws_surface** out_surface, qaws_surface_swept_desc const* desc);
```

A 2D profile swept along a 3D path:
`S(u,v) = path(u) + profile_x(v) N(u) + profile_y(v) B(u)`, with N and B the
Frenet normal and binormal of the path. u follows the path's domain, v the
profile's.

### Ruled (qaws_surface_ruled.h)

```c
typedef struct qaws_surface_ruled_desc {
	qaws_curve const* curve_a;   /* first boundary curve, 3D (borrowed) */
	qaws_curve const* curve_b;   /* second boundary curve, 3D (borrowed) */
} qaws_surface_ruled_desc;

qaws_status qaws_surface_create_ruled(qaws_surface** out_surface, qaws_surface_ruled_desc const* desc);
```

`S(u,v) = (1-v) curve_a(u) + v curve_b(u)`. u in [0, 1] is mapped to each
curve's own domain; v blends from `curve_a` (v = 0) to `curve_b` (v = 1).
Both curves must be 3D.

### Coons (qaws_surface_coons.h)

```c
typedef struct qaws_surface_coons_desc {
	qaws_curve const* c0;   /* bottom edge (v=0), u from 0 to 1, 3D (borrowed) */
	qaws_curve const* c1;   /* top edge (v=1), u from 0 to 1, 3D (borrowed) */
	qaws_curve const* d0;   /* left edge (u=0), v from 0 to 1, 3D (borrowed) */
	qaws_curve const* d1;   /* right edge (u=1), v from 0 to 1, 3D (borrowed) */
} qaws_surface_coons_desc;

qaws_status qaws_surface_create_coons(qaws_surface** out_surface, qaws_surface_coons_desc const* desc);
```

A bilinearly blended Coons patch, `S = Lc + Ld - B`: the ruled surface
between c0 and c1, plus the ruled surface between d0 and d1, minus the
bilinear patch of the four corners. The corners must agree:
c0(0) = d0(0), c0(1) = d1(0), c1(0) = d0(1), c1(1) = d1(1). All curves 3D.

### Extrusion (qaws_surface_extrusion.h)

```c
typedef struct qaws_surface_extrusion_desc {
	qaws_curve const* profile;   /* 2D or 3D (borrowed) */
	qaws_vec3 direction;         /* extrusion direction */
	qaws_scalar length;          /* extrusion length; 0 = use |direction| as the unit */
} qaws_surface_extrusion_desc;

qaws_status qaws_surface_create_extrusion(qaws_surface** out_surface,
	qaws_surface_extrusion_desc const* desc);
```

`S(u,v) = profile(u) + v direction`. u follows the profile, v the distance
along the direction (0 to length). A 2D profile is lifted to the z = 0
plane.

### Revolution (qaws_surface_revolution.h)

```c
typedef struct qaws_surface_revolution_desc {
	qaws_curve const* profile;    /* 2D profile curve (borrowed) */
	qaws_vec3 axis_origin;        /* point on the rotation axis */
	qaws_vec3 axis_direction;     /* rotation axis direction (normalized; zero length is invalid) */
	qaws_scalar angle;            /* rotation in radians; 0 = 2*pi */
} qaws_surface_revolution_desc;

qaws_status qaws_surface_create_revolution(qaws_surface** out_surface,
	qaws_surface_revolution_desc const* desc);
```

A 2D profile rotated about an axis. The profile's x is the radius r(v) and
its second coordinate the height h(v):
`S(u,v) = axis_origin + r(v) cos(theta) X + r(v) sin(theta) Y + h(v) Z`,
theta = u x angle, with Z along the axis. u in [0, 1] maps to the rotation
[0, angle]; v follows the profile's domain. The profile must be 2D.

### Pipe (qaws_surface_pipe.h)

```c
typedef struct qaws_surface_pipe_desc {
	qaws_curve const* path;   /* 3D path curve (borrowed) */
	qaws_scalar radius_x;     /* cross-section radius along the normal */
	qaws_scalar radius_y;     /* radius along the binormal; 0 = radius_x */
} qaws_surface_pipe_desc;

qaws_status qaws_surface_create_pipe(qaws_surface** out_surface, qaws_surface_pipe_desc const* desc);
```

A circular or elliptical section swept along a 3D path:
`S(u,v) = path(u) + radius_x cos(2 pi v) N(u) + radius_y sin(2 pi v) B(u)`.
N and B are a rotation-minimizing (parallel transport) frame computed by
double reflection. u in [0, 1] runs along the path, v in [0, 1] around the
section.

### Loft (qaws_surface_loft.h)

```c
typedef struct qaws_surface_loft_desc {
	qaws_curve const* const* sections;   /* N 3D cross-section curves (borrowed) */
	unsigned int section_count;          /* N >= 2 */
	qaws_scalar const* v_parameters;     /* N values in [0,1], one per section; NULL = uniform */
} qaws_surface_loft_desc;

qaws_status qaws_surface_create_loft(qaws_surface** out_surface, qaws_surface_loft_desc const* desc);
```

A skinning surface through the sections: at v = `v_parameters[i]` it
passes through `sections[i](u)`, and cubic B-spline blending interpolates
between sections. u in [0, 1] runs along each section, v in [0, 1] across
them.

### Gordon (qaws_surface_gordon.h)

```c
typedef struct qaws_surface_gordon_desc {
	qaws_curve const* const* u_curves;   /* M u-isoparameter curves, 3D (borrowed) */
	unsigned int u_curve_count;          /* M >= 2 */
	qaws_scalar const* v_params;         /* M values: where each u-curve sits in v */
	qaws_curve const* const* v_curves;   /* N v-isoparameter curves, 3D (borrowed) */
	unsigned int v_curve_count;          /* N >= 2 */
	qaws_scalar const* u_params;         /* N values: where each v-curve sits in u */
} qaws_surface_gordon_desc;

qaws_status qaws_surface_create_gordon(qaws_surface** out_surface, qaws_surface_gordon_desc const* desc);
```

Interpolates a curve network, generalizing the Coons patch:
`S = L_u + L_v - T`, the loft through the u-curves, plus the loft through the
v-curves, minus the tensor product through the intersection points
Q[i][j] = u_curve[i](v_param[j]). The network must be compatible:
u_curve[i] at v_param[j] equals v_curve[j] at u_param[i].

---

## Offset (qaws_surface_offset.h)

```c
typedef struct qaws_surface_offset_desc {
	qaws_surface const* base;   /* base surface (borrowed) */
	qaws_scalar distance;       /* positive = along the normal */
} qaws_surface_offset_desc;

qaws_status qaws_surface_create_offset(qaws_surface** out_surface, qaws_surface_offset_desc const* desc);
```

`S_off(u,v) = S_base(u,v) + distance N_base(u,v)`: the base moved a constant
distance along its unit normal.
