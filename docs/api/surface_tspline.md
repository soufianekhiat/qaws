# API Reference: qaws_surface_tspline.h

T-spline surfaces: locally refinable surfaces in which every control point
carries its own local knot vectors, so detail can be added in one area
without inserting whole rows or columns of control points. The result is a
`qaws_surface` of kind `QAWS_SURFACE_KIND_TSPLINE` (see
[surface.md](surface.md)).

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)).

---

## Model

Control point i has a position P_i, a weight w_i and local knot vectors of
length 2*u_degree+2 and 2*v_degree+2. Its basis function is the product of
the two univariate B-spline basis functions those vectors define:

```
B_i(u,v) = N[u_knots_i](u) * N[v_knots_i](v)
S(u,v)   = sum_i w_i P_i B_i(u,v) / sum_i w_i B_i(u,v)
```

With all weights 1 this is `sum_i P_i B_i / sum_i B_i`. In a T-spline the
local knot vectors follow from the T-mesh connectivity; this API takes them
explicitly from the caller.

## Types

```c
typedef struct qaws_tspline_control_point {
	qaws_vec3 position;
	qaws_scalar weight;           /* 1.0 for non-rational */
	qaws_scalar const* u_knots;   /* local u-knot vector, length 2*u_degree+2 */
	qaws_scalar const* v_knots;   /* local v-knot vector, length 2*v_degree+2 */
} qaws_tspline_control_point;

typedef struct qaws_surface_tspline_desc {
	qaws_tspline_control_point const* control_points;
	unsigned int control_point_count;
	unsigned int u_degree;        /* typically 3 */
	unsigned int v_degree;        /* typically 3 */
	qaws_range u_range;           /* parameter domain in u */
	qaws_range v_range;           /* parameter domain in v */
} qaws_surface_tspline_desc;
```

## Functions

```c
qaws_status qaws_surface_create_tspline(qaws_surface** out_surface, qaws_surface_tspline_desc const* desc);
```

Creates the T-spline over `u_range` x `v_range`. Free it with
`qaws_surface_destroy`.
