# API Reference: portable core headers

Evaluation and differentiation kernels that compile unchanged as C, HLSL,
GLSL and Halide: no pointers to the library, no allocation, fixed-size
arrays, branch-free selects where a value decides. A host includes them in
a shader, a Halide pipeline or plain C to evaluate curve and surface data
it has laid out itself; they reproduce the C runtime's evaluation of the
same families.

Installed headers:

- bootstrap: `qaws_hlsl.h`, `qaws_glsl.h`, `qaws_halide.h`;
- shared types: `qaws_core_types.h`;
- kernels, installed under `core/`: `qaws_decasteljau_core.h`,
  `qaws_cubic_poly_core.h`, `qaws_horner_core.h`,
  `qaws_rational_bezier_core.h`, `qaws_arc_core.h`,
  `qaws_bspline_basis_core.h`, `qaws_bspline_eval_core.h`,
  `qaws_nurbs_eval_core.h`, `qaws_surface_core.h`, `qaws_dual_core.h`,
  `qaws_bspline_diff_core.h`, `qaws_bspline_sampling_core.h`,
  `qaws_bspline_surface_sampling_core.h`.

Every function is `QAWS_INLINE` and takes `qaws_scalar`; the macros it uses
(`QAWS_OUT`, `QAWS_SELECT`, `QAWS_LITERAL`, ...) are described in
[precision.md](precision.md). On the C backend the function names are
precision-bearing (`name_f32` / `name_f64` through `QAWS_SFX`); each core
header has a single include guard, so a translation unit gets the
precision current at its first inclusion. The struct types the core
headers define (`qaws_eval_2d`, `qaws_eval_3d`, `qaws_surface_eval_3d`,
`qaws_dual1`, `qaws_dual3`, `qaws_vec3_pair`, `qaws_vec3_scalar`) have a
single, unsuffixed name. See also
[backends.md](../backends.md).

---

## Bootstrap headers

```c
#include "qaws_hlsl.h"     /* QAWS_BACKEND 1 (HLSL), unless already defined */
#include "qaws_glsl.h"     /* QAWS_BACKEND 2 (GLSL) */
#include "qaws_halide.h"   /* QAWS_BACKEND 3 (Halide) */
```

Each defines `QAWS_BACKEND` if it is not already defined, then includes
`qaws_platform.h` and `qaws_core_types.h`. Include one before any core
header.

## qaws_core_types.h

```c
QAWS_STRUCT(qaws_vec2) { qaws_scalar x, y; } QAWS_STRUCT_END(qaws_vec2);
QAWS_STRUCT(qaws_vec3) { qaws_scalar x, y, z; } QAWS_STRUCT_END(qaws_vec3);

QAWS_STRUCT(qaws_eval_2d) { qaws_vec2 position; qaws_vec2 d1; qaws_vec2 d2; } QAWS_STRUCT_END(qaws_eval_2d);
QAWS_STRUCT(qaws_eval_3d) { qaws_vec3 position; qaws_vec3 d1; qaws_vec3 d2; } QAWS_STRUCT_END(qaws_eval_3d);

#define QAWS_CORE_MAX_DEGREE 16                          /* may be defined before inclusion */
#define QAWS_CORE_MAX_POINTS (QAWS_CORE_MAX_DEGREE + 1)
```

- `qaws_vec2` and `qaws_vec3` are defined here for shader and Halide
  consumers; on the C backend they are skipped when `qaws_types.h` already
  defined them (the layouts are the same).
- `qaws_eval_2d` / `qaws_eval_3d` are the kernels' results: position, first
  and second derivative (guard `QAWS_EVAL_2D_DEFINED`).
- `QAWS_CORE_MAX_DEGREE` bounds every fixed-size array; array parameters
  below are sized from it.

## Conventions

- **Flat arrays.** Control points are flat scalars: `[x0, y0, x1, y1, ...]`
  in 2D, `[x0, y0, z0, ...]` in 3D.
- **Curve evaluation flags** (`int eval_flags`) use the values of
  `QAWS_EVAL_FLAG_*`: `0x2` first derivative, `0x4` second. The position is
  always computed; derivatives not asked for are zero.
- **Spline windows.** B-spline and NURBS kernels take `local_cp`, the
  degree + 1 control points of the span, `local_knots`, the knot window
  around it (2 (degree + 1) values), and `span`, the span index inside that
  window.

## qaws_decasteljau_core.h

Bezier curves by De Casteljau.

```c
qaws_vec2 qaws_bezier3_eval_2d(qaws_vec2 p0, qaws_vec2 p1, qaws_vec2 p2, qaws_vec2 p3, qaws_scalar t);
qaws_vec3 qaws_bezier3_eval_3d(qaws_vec3 p0, qaws_vec3 p1, qaws_vec3 p2, qaws_vec3 p3, qaws_scalar t);
qaws_vec2 qaws_decasteljau_2d(qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2], int degree, qaws_scalar t);
qaws_vec3 qaws_decasteljau_3d(qaws_scalar cp[QAWS_CORE_MAX_POINTS * 3], int degree, qaws_scalar t);
void qaws_bezier_deriv_points_2d(QAWS_OUT qaws_scalar out_dp[QAWS_CORE_MAX_POINTS * 2], qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2], int degree);
void qaws_bezier_deriv_points_3d(QAWS_OUT qaws_scalar out_dp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar cp[QAWS_CORE_MAX_POINTS * 3], int degree);
qaws_eval_2d qaws_decasteljau_eval_2d(qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2], int degree, qaws_scalar t, int eval_flags);
qaws_eval_3d qaws_decasteljau_eval_3d(qaws_scalar cp[QAWS_CORE_MAX_POINTS * 3], int degree, qaws_scalar t, int eval_flags);
```

| Function | Does |
|---|---|
| `qaws_bezier3_eval_*` | cubic Bezier position, four points by value |
| `qaws_decasteljau_*` | position of a Bezier curve of any degree up to `QAWS_CORE_MAX_DEGREE` |
| `qaws_bezier_deriv_points_*` | the degree - 1 derivative control points `degree (P[i+1] - P[i])` |
| `qaws_decasteljau_eval_*` | position, d1 (degree >= 1) and d2 (degree >= 2) |

## qaws_cubic_poly_core.h

```c
qaws_eval_2d qaws_cubic_eval_2d(qaws_vec2 a, qaws_vec2 b, qaws_vec2 c, qaws_vec2 d, qaws_scalar t, int eval_flags);
qaws_eval_3d qaws_cubic_eval_3d(qaws_vec3 a, qaws_vec3 b, qaws_vec3 c, qaws_vec3 d, qaws_scalar t, int eval_flags);
```

P(t) = a t^3 + b t^2 + c t + d with its derivatives: the evaluation shared
by Hermite, Catmull-Rom and trajectory spans (coefficients per span).

## qaws_horner_core.h

```c
qaws_vec2 qaws_horner_2d(qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 2], int degree, qaws_scalar t);
qaws_vec3 qaws_horner_3d(qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 3], int degree, qaws_scalar t);
qaws_eval_2d qaws_horner_eval_2d(qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 2], int degree, qaws_scalar t, int eval_flags);
qaws_eval_3d qaws_horner_eval_3d(qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 3], int degree, qaws_scalar t, int eval_flags);
```

P(t) = c0 + c1 t + ... + cn t^n by Horner; the `_eval_` forms add d1
(degree >= 1) and d2 (degree >= 2).

## qaws_rational_bezier_core.h

```c
qaws_eval_2d qaws_rational_bezier_eval_2d(qaws_scalar weighted_cp[QAWS_CORE_MAX_POINTS * 3], int degree, qaws_scalar t, int eval_flags);
qaws_eval_3d qaws_rational_bezier_eval_3d(qaws_scalar weighted_cp[QAWS_CORE_MAX_POINTS * 4], int degree, qaws_scalar t, int eval_flags);
```

Homogeneous control points `(w x, w y, w)` in 2D, `(w x, w y, w z, w)` in
3D; De Casteljau on the homogeneous points, then the quotient rule. The
second derivative is computed for degree 1 too (the quotient's does not
vanish).

## qaws_arc_core.h

```c
qaws_vec2 qaws_arc_eval_2d(qaws_vec2 center, qaws_vec2 axis_u, qaws_vec2 axis_v,
    qaws_scalar theta_start, qaws_scalar theta_range, qaws_scalar t);
qaws_vec3 qaws_arc_eval_3d(qaws_vec3 center, qaws_vec3 axis_u, qaws_vec3 axis_v,
    qaws_scalar theta_start, qaws_scalar theta_range, qaws_scalar t);
qaws_eval_2d qaws_arc_eval_full_2d(qaws_vec2 center, qaws_vec2 axis_u, qaws_vec2 axis_v,
    qaws_scalar theta_start, qaws_scalar theta_range, qaws_scalar t, int eval_flags);
qaws_eval_3d qaws_arc_eval_full_3d(qaws_vec3 center, qaws_vec3 axis_u, qaws_vec3 axis_v,
    qaws_scalar theta_start, qaws_scalar theta_range, qaws_scalar t, int eval_flags);
```

`center + cos(theta) axis_u + sin(theta) axis_v` with
`theta = theta_start + t theta_range`, t in [0, 1]. The axes carry the
radii. `_full_` adds d1 = (-sin u + cos v) theta_range and
d2 = (-cos u - sin v) theta_range^2.

## qaws_bspline_basis_core.h

```c
#define QAWS_CORE_MAX_DERIV 3

int qaws_find_span(qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2], int degree, int num_cp, qaws_scalar t);
          /* Halide: returns Halide::Expr */
void qaws_bspline_basis(QAWS_OUT qaws_scalar out_N[QAWS_CORE_MAX_POINTS], qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2],
    int degree, int span, qaws_scalar t);
void qaws_bspline_basis_derivs(qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2], int degree, int span, qaws_scalar t, int max_k,
    QAWS_OUT qaws_scalar out_ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS]);
```

- **`qaws_find_span`** finds i with `knots[i] <= t < knots[i + 1]`. In C,
  HLSL and GLSL: a bounded binary search (20 iterations); at the domain
  ends, the nearest non-empty span. In Halide: a branch-free linear scan
  returning an `Expr` evaluated per sample (at most `QAWS_CORE_MAX_POINTS`
  selects).
- **`qaws_bspline_basis`** writes the degree + 1 non-zero basis values at
  t (a zero denominator gives 0).
- **`qaws_bspline_basis_derivs`** writes basis values and derivatives up
  to order `max_k` (at most `QAWS_CORE_MAX_DERIV`), row-major
  `(max_k + 1) x (degree + 1)`, on stack buffers.

## qaws_bspline_eval_core.h

```c
qaws_eval_2d qaws_bspline_eval_2d(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 2], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree, int span, qaws_scalar t, int eval_flags);
qaws_eval_3d qaws_bspline_eval_3d(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree, int span, qaws_scalar t, int eval_flags);
```

B-spline position and derivatives on a span window. Derivatives above the
degree are zero.

## qaws_nurbs_eval_core.h

```c
qaws_eval_2d qaws_nurbs_eval_2d(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 2], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_weights[QAWS_CORE_MAX_POINTS], int degree, int span, qaws_scalar t, int eval_flags);
qaws_eval_3d qaws_nurbs_eval_3d(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_weights[QAWS_CORE_MAX_POINTS], int degree, int span, qaws_scalar t, int eval_flags);
```

Rational evaluation: weighted sums A^(k) and W^(k), then the quotient rule.
Cartesian control points with degree + 1 weights.

## qaws_surface_core.h

```c
#define QAWS_SURFACE_CORE_POSITION  0x1
#define QAWS_SURFACE_CORE_DU        0x2
#define QAWS_SURFACE_CORE_DV        0x4
#define QAWS_SURFACE_CORE_NORMAL    0x8
#define QAWS_SURFACE_CORE_DUU       0x10
#define QAWS_SURFACE_CORE_DUV       0x20
#define QAWS_SURFACE_CORE_DVV       0x40

QAWS_STRUCT(qaws_surface_eval_3d) {
    qaws_vec3 position;
    qaws_vec3 du, dv;
    qaws_vec3 duu, dvv, duv;
    qaws_vec3 normal;
} QAWS_STRUCT_END(qaws_surface_eval_3d);

qaws_surface_eval_3d qaws_surface_bezier_eval_3d(qaws_scalar cp[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3],
    int u_degree, int v_degree, int v_count, qaws_scalar u, qaws_scalar v, int eval_flags);
qaws_surface_eval_3d qaws_surface_bspline_eval_3d(qaws_scalar cp[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar u_knots[QAWS_CORE_MAX_POINTS * 2], qaws_scalar v_knots[QAWS_CORE_MAX_POINTS * 2],
    int u_degree, int v_degree, int u_count, int v_count, int u_span, int v_span,
    qaws_scalar u, qaws_scalar v, int eval_flags);
```

- The flag values are those of `QAWS_SURFACE_EVAL_*`
  ([surface.md](surface.md)): one flag word for the runtime and the cores.
- Control points row-major: `cp[(i * v_count + j) * 3 + k]`, row i along u.
- **Bezier**: `v_count = v_degree + 1`, u and v in [0, 1].
- **B-spline**: whole knot vectors and the precomputed spans `u_span`,
  `v_span`.
- The normal is `du x dv` normalized; below a length of 1e-12 it is
  (0, 0, 1). Asking for the normal also fills du and dv.

## qaws_dual_core.h

Second-order forward numbers and first-order adjoint rules, the building
blocks of the differential kernels. A dual value carries a primal v, its
derivative t along one tangent and the second derivative tt along the same
tangent (not divided by 2): x(e) = v + e t + e^2/2 tt + O(e^3).

```c
QAWS_STRUCT(qaws_dual1) { qaws_scalar v, t, tt; } QAWS_STRUCT_END(qaws_dual1);
QAWS_STRUCT(qaws_dual3) { qaws_vec3 v, t, tt; } QAWS_STRUCT_END(qaws_dual3);
QAWS_STRUCT(qaws_vec3_pair) { qaws_vec3 a, b; } QAWS_STRUCT_END(qaws_vec3_pair);
QAWS_STRUCT(qaws_vec3_scalar) { qaws_vec3 x; qaws_scalar s; } QAWS_STRUCT_END(qaws_vec3_scalar);
```

| Function | Signature and rule |
|---|---|
| `qaws_v3` | `qaws_vec3 qaws_v3(qaws_scalar x, qaws_scalar y, qaws_scalar z)` |
| `qaws_v3_add` | `qaws_vec3 qaws_v3_add(qaws_vec3 a, qaws_vec3 b)` |
| `qaws_v3_sub` | `qaws_vec3 qaws_v3_sub(qaws_vec3 a, qaws_vec3 b)` |
| `qaws_v3_scale` | `qaws_vec3 qaws_v3_scale(qaws_vec3 a, qaws_scalar s)` |
| `qaws_v3_axpy` | `qaws_vec3 qaws_v3_axpy(qaws_vec3 a, qaws_vec3 b, qaws_scalar s)`: a + b s |
| `qaws_v3_dot` | `qaws_scalar qaws_v3_dot(qaws_vec3 a, qaws_vec3 b)` |
| `qaws_v3_cross` | `qaws_vec3 qaws_v3_cross(qaws_vec3 a, qaws_vec3 b)` |
| `qaws_v3_zero` | `qaws_vec3 qaws_v3_zero(void)` |
| `qaws_dual1_make` | `qaws_dual1 qaws_dual1_make(qaws_scalar v, qaws_scalar t, qaws_scalar tt)` |
| `qaws_dual1_const` | `qaws_dual1 qaws_dual1_const(qaws_scalar v)`: t = tt = 0 |
| `qaws_dual3_make` | `qaws_dual3 qaws_dual3_make(qaws_vec3 v, qaws_vec3 t, qaws_vec3 tt)` |
| `qaws_dual3_const` | `qaws_dual3 qaws_dual3_const(qaws_vec3 v)` |
| `qaws_dual1_add` | `qaws_dual1 qaws_dual1_add(qaws_dual1 a, qaws_dual1 b)` |
| `qaws_dual1_sub` | `qaws_dual1 qaws_dual1_sub(qaws_dual1 a, qaws_dual1 b)` |
| `qaws_dual1_mul` | `qaws_dual1 qaws_dual1_mul(qaws_dual1 a, qaws_dual1 b)` |
| `qaws_dual1_div` | `qaws_dual1 qaws_dual1_div(qaws_dual1 a, qaws_dual1 b)` |
| `qaws_dual1_sqrt` | `qaws_dual1 qaws_dual1_sqrt(qaws_dual1 a)`: derivatives 0 at 0 |
| `qaws_dual3_add` | `qaws_dual3 qaws_dual3_add(qaws_dual3 a, qaws_dual3 b)` |
| `qaws_dual3_sub` | `qaws_dual3 qaws_dual3_sub(qaws_dual3 a, qaws_dual3 b)` |
| `qaws_dual3_scale` | `qaws_dual3 qaws_dual3_scale(qaws_dual3 x, qaws_dual1 s)`: x s |
| `qaws_dual3_dot` | `qaws_dual1 qaws_dual3_dot(qaws_dual3 a, qaws_dual3 b)` |
| `qaws_dual3_cross` | `qaws_dual3 qaws_dual3_cross(qaws_dual3 a, qaws_dual3 b)` |
| `qaws_dual3_div` | `qaws_dual3 qaws_dual3_div(qaws_dual3 x, qaws_dual1 w)`: x / w |
| `qaws_dual3_length` | `qaws_dual1 qaws_dual3_length(qaws_dual3 a)` |
| `qaws_dual3_mul_const` | `qaws_dual3 qaws_dual3_mul_const(qaws_dual3 a, qaws_scalar s)` |
| `qaws_dual3_normalize` | `qaws_dual3 qaws_dual3_normalize(qaws_dual3 a)`: zero dual below a length of 1e-30 |
| `qaws_dot_adjoint` | `qaws_vec3_pair qaws_dot_adjoint(qaws_vec3 a, qaws_vec3 b, qaws_scalar s_bar)`: s = a . b |
| `qaws_cross_adjoint` | `qaws_vec3_pair qaws_cross_adjoint(qaws_vec3 a, qaws_vec3 b, qaws_vec3 c_bar)`: a_bar = b x c_bar, b_bar = c_bar x a |
| `qaws_length_adjoint` | `qaws_vec3 qaws_length_adjoint(qaws_vec3 a, qaws_scalar l_bar)`: l = abs(a) |
| `qaws_normalize_adjoint` | `qaws_vec3 qaws_normalize_adjoint(qaws_vec3 a, qaws_vec3 n_bar)`: (n_bar - n (n . n_bar)) / abs(a) |
| `qaws_div_adjoint` | `qaws_vec3_scalar qaws_div_adjoint(qaws_vec3 x, qaws_scalar w, qaws_vec3 p_bar)`: x_bar = p_bar / w, w_bar = -(p_bar . p) / w |
| `qaws_scale_adjoint` | `qaws_vec3_scalar qaws_scale_adjoint(qaws_vec3 x, qaws_scalar s, qaws_vec3 y_bar)`: x_bar = s y_bar, s_bar = x . y_bar |

Adjoint rules take primal values and an output adjoint and return the input
adjoint contributions, which the caller adds.

## qaws_bspline_diff_core.h

B-spline derivative kernels on span windows; jets are `qaws_eval_3d`
(position, d1, d2). They reproduce `qaws_curve_eval_batch_tangent_3d` and
`_adjoint_3d` of the C runtime for B-spline curves ([diff.md](diff.md)).

```c
void qaws_bspline_diff_rows(qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2], int degree, int span, qaws_scalar t,
    QAWS_OUT qaws_scalar out_ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS]);
qaws_vec3 qaws_bspline_diff_row_sum(qaws_scalar values[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS], int degree, int k);
qaws_eval_3d qaws_bspline_tangent_3d(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2], int degree, int span, qaws_scalar t, qaws_scalar t_dot);
qaws_vec3 qaws_bspline_adjoint_cp_3d(qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2], int degree, int span, qaws_scalar t,
    int j, qaws_eval_3d ybar);
qaws_scalar qaws_bspline_adjoint_t_3d(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree, int span, qaws_scalar t, qaws_eval_3d ybar);
```

| Function | Does |
|---|---|
| `qaws_bspline_diff_rows` | basis rows 0..3, rows above the degree set to zero |
| `qaws_bspline_diff_row_sum` | sum_j N_j^(k) values[j] |
| `qaws_bspline_tangent_3d` | tangent jet: dC^(k) = sum N^(k) P_dot + t_dot sum N^(k+1) P, k = 0, 1, 2 |
| `qaws_bspline_adjoint_cp_3d` | sum_k N_j^(k)(t) ybar_k for local control point j: the per-element term of a gather kernel (one thread per control point over the samples it supports) |
| `qaws_bspline_adjoint_t_3d` | parameter adjoint sum_k ybar_k . C^(k+1)(t) of one sample |

## qaws_bspline_sampling_core.h

Inverse-CDF sampling and integral functionals of B-spline curves, matching
`qaws_curve_cdf_sample_*` and `qaws_curve_functional_*` with the default
quadrature: each span split into `QAWS_SAMPLING_PIECES` pieces of 8 Gauss
points ([sampling.md](sampling.md), [diff.md](diff.md)). Window
conventions as above; `local_cp_dot` holds the control point tangents
(zeros for none). Kernels run one thread per span, per sample or per
control point.

```c
#define QAWS_SAMPLING_PIECES 8

#define QAWS_CORE_INTEGRAND_SPEED        0   /* abs(C'): arc length */
#define QAWS_CORE_INTEGRAND_CURVATURE    1   /* sqrt(floor^2 abs(C')^2 + abs(C' x C'')^2 / abs(C')^4) */
#define QAWS_CORE_INTEGRAND_BENDING      2   /* abs(C'')^2 */
#define QAWS_CORE_INTEGRAND_CURVATURE_SQ 3   /* abs(C' x C'')^2 / abs(C')^5 */
```

```c
void qaws_gauss8(QAWS_OUT qaws_scalar xw[2], int k);
qaws_dual1 qaws_core_integrand(int id, qaws_scalar curvature_floor, qaws_dual3 d1, qaws_dual3 d2);
void qaws_bspline_jet_rows(QAWS_OUT qaws_vec3 c[4], QAWS_OUT qaws_vec3 c_dot[3],
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2], int degree, int span, qaws_scalar t);
qaws_dual1 qaws_bspline_integrand_at(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2], int id, qaws_scalar curvature_floor, int degree, int span, qaws_scalar t);
qaws_dual1 qaws_bspline_integrate(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2], int id, qaws_scalar curvature_floor, int degree, int span,
    qaws_scalar a, qaws_scalar b);
qaws_scalar qaws_bspline_cdf_solve(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int id, qaws_scalar curvature_floor, int degree, int span, qaws_scalar a, qaws_scalar b, qaws_scalar remaining);
void qaws_bspline_cdf_tangent(QAWS_OUT qaws_scalar out_t[2], QAWS_OUT qaws_vec3 out_p[2],
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2], int id, qaws_scalar curvature_floor, int degree, int span,
    qaws_scalar t, qaws_scalar sigma_dot, qaws_scalar sigma_dot2, qaws_scalar measure_dot, qaws_scalar measure_dot2);
qaws_scalar qaws_bspline_cdf_lambda(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int id, qaws_scalar curvature_floor, int degree, int span, qaws_scalar t, qaws_scalar t_bar, qaws_vec3 p_bar);
qaws_vec3 qaws_bspline_node_adjoint_cp(qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3], qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int id, qaws_scalar curvature_floor, int degree, int span, qaws_scalar t, int j, qaws_scalar c);
```

| Function | Does |
|---|---|
| `qaws_gauss8` | node `xw[0]` and weight `xw[1]` of point k (0..7) of the 8-point Gauss-Legendre rule on [-1, 1] |
| `qaws_core_integrand` | the integrand `id` on dual jets (C', C'') |
| `qaws_bspline_jet_rows` | C', C'', C''' at t in `c[1..3]` (`c[0]` the position) and the first rows along the control point tangents in `c_dot` |
| `qaws_bspline_integrand_at` | the integrand at t along the control point tangents (v, t, tt) |
| `qaws_bspline_integrate` | the integral over [a, b] inside one span with its first and second tangents (zero when b <= a); the span totals give the prefix measures |
| `qaws_bspline_cdf_solve` | the parameter in [a, b] at which the measure from a reaches `remaining`: safeguarded Newton, fixed iteration count, branch free |
| `qaws_bspline_cdf_tangent` | `out_t = (t', t'')` and `out_p` (position tangent, second tangent) of a sample, from the target rates `sigma_dot`, `sigma_dot2` and the measure rates at fixed t (implicit function theorem) |
| `qaws_bspline_cdf_lambda` | adjoint multiplier (t_bar + p_bar . C') / m of one sample, 0 where m vanishes |
| `qaws_bspline_node_adjoint_cp` | c (N_j' df/dC' + N_j'' df/dC'') at one quadrature node for local control point j: the per-element term of a gather (functional gradients, and the measure terms of the sampling adjoint) |

## qaws_bspline_surface_sampling_core.h

Inverse-CDF warping of the unit square onto a B-spline surface by its area
measure w = abs(S_u x S_v), matching `qaws_surface_cdf_sample_*` of the
runtime with the area measure, `cells` cells per direction and 8 Gauss
points per cell:

- A(u) = xi_u A_total, B(v; u) = xi_v B(v1; u),
- A(u) = integral from u0 to u of B(v1; .), B(b; u) = integral from v0 to
  b of w(u, .).

The whole control grid is passed (a v-line crosses several knot spans):
`cp[(i * v_count + j) * 3 + c]`, at most `QAWS_CORE_MAX_POINTS` per
direction; `cp_dot` holds the control point tangents and `dir` scales them
(0: value and parameter rates only). Duals carry first and second rates
along one path: control points moving at `cp_dot`, u and v at their own
rates. `prefix_v`, `prefix_t`, `prefix_tt` are the exclusive prefix sums of
the cell masses (value, first and second rate).

```c
#define QAWS_CORE_SURFACE_GRID (QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3)
#define QAWS_CORE_MAX_CELLS 32

#define QAWS_PATCH_PARAMS \
    qaws_scalar cp[QAWS_CORE_SURFACE_GRID], qaws_scalar cp_dot[QAWS_CORE_SURFACE_GRID], \
    qaws_scalar u_knots[QAWS_CORE_MAX_POINTS * 2], qaws_scalar v_knots[QAWS_CORE_MAX_POINTS * 2], \
    int u_degree, int v_degree, int u_count, int v_count
#define QAWS_PATCH_ARGS cp, cp_dot, u_knots, v_knots, u_degree, v_degree, u_count, v_count
```

`QAWS_PATCH_PARAMS` is the parameter list of the patch kernels below and
`QAWS_PATCH_ARGS` the matching argument list.

```c
qaws_vec3 qaws_patch_sum(qaws_scalar grid[QAWS_CORE_SURFACE_GRID],
    qaws_scalar nu[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS], qaws_scalar nv[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS],
    int u_degree, int v_degree, int v_count, int us, int vs, int a, int b);
void qaws_patch_rows(qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2], int degree, int span, qaws_scalar t,
    QAWS_OUT qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS]);
void qaws_patch_jets(QAWS_OUT qaws_dual3 y[3], QAWS_PATCH_PARAMS, qaws_scalar dir, qaws_dual1 u, qaws_dual1 v);
qaws_dual1 qaws_patch_area_at(QAWS_PATCH_PARAMS, qaws_scalar dir, qaws_dual1 u, qaws_dual1 v);
qaws_dual1 qaws_dual_affine(qaws_scalar a, qaws_dual1 x, qaws_scalar b);
qaws_dual1 qaws_patch_line(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, qaws_dual1 u, qaws_dual1 b);
qaws_dual1 qaws_patch_cell_mass(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, int c);
qaws_dual1 qaws_patch_marginal(qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar prefix_t[QAWS_CORE_MAX_CELLS + 1],
    qaws_scalar prefix_tt[QAWS_CORE_MAX_CELLS + 1], QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, qaws_dual1 u);
qaws_dual1 qaws_patch_F(qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar prefix_t[QAWS_CORE_MAX_CELLS + 1],
    qaws_scalar prefix_tt[QAWS_CORE_MAX_CELLS + 1], QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, qaws_dual1 u, qaws_dual1 xi);
qaws_dual1 qaws_patch_G(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, qaws_dual1 u, qaws_dual1 v, qaws_dual1 xi);
void qaws_patch_cdf_solve(QAWS_OUT qaws_scalar out_uv[2], qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1],
    QAWS_PATCH_PARAMS, int cells, qaws_scalar xu, qaws_scalar xv);
void qaws_patch_cdf_tangent(QAWS_OUT qaws_scalar result[4], QAWS_OUT qaws_vec3 out_p[2],
    qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar prefix_t[QAWS_CORE_MAX_CELLS + 1],
    qaws_scalar prefix_tt[QAWS_CORE_MAX_CELLS + 1], QAWS_PATCH_PARAMS, int cells, qaws_scalar xu, qaws_scalar xv,
    qaws_scalar xu_dot, qaws_scalar xv_dot, qaws_scalar u, qaws_scalar v);
void qaws_patch_cdf_multipliers(QAWS_OUT qaws_scalar result[4], qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1],
    QAWS_PATCH_PARAMS, int cells, qaws_scalar xu, qaws_scalar xv, qaws_scalar u,
    qaws_scalar v, qaws_scalar u_bar, qaws_scalar v_bar, qaws_vec3 p_bar);
qaws_vec3 qaws_patch_area_adjoint_cp(QAWS_PATCH_PARAMS, qaws_scalar u, qaws_scalar v, int i, int j, qaws_scalar c);
qaws_vec3 qaws_patch_line_adjoint_cp(QAWS_PATCH_PARAMS, int cells, qaws_scalar u, qaws_scalar b, int i, int j, qaws_scalar c);
```

| Function | Does |
|---|---|
| `qaws_patch_sum` | sum_ij N_i^(a)(u) N_j^(b)(v) P_ij over the spans (us, vs) |
| `qaws_patch_rows` | basis rows 0..3 at t, rows above the degree zero |
| `qaws_patch_jets` | dual (S, S_u, S_v) at the dual point (u, v), along the path |
| `qaws_patch_area_at` | w = abs(S_u x S_v) at the dual point |
| `qaws_dual_affine` | a x + b on a dual |
| `qaws_patch_line` | B(b; u): full v cells, then the partial cell whose nodes and weights move with b |
| `qaws_patch_cell_mass` | measure of u cell c with its rates (one thread per u cell); an exclusive prefix sum gives the prefixes |
| `qaws_patch_marginal` | A(u): the prefix plus the partial u cell moving with u |
| `qaws_patch_F` | F(u) = A(u) - xi_u A_total |
| `qaws_patch_G` | G(v; u) = B(v; u) - xi_v B(v1; u) |
| `qaws_patch_cdf_solve` | `out_uv = (u, v)` of the point (xu, xv): safeguarded Newton in the cell holding the target, fixed iteration count |
| `qaws_patch_cdf_tangent` | `result = (u', v', u'', v'')` and `out_p` (position tangent, second) along `cp_dot` and (xu_dot, xv_dot): one dual Newton step per order |
| `qaws_patch_cdf_multipliers` | `result = (lambda, mu, lambda A_total, mu B(v1; u))`: the adjoint multipliers and the xi_u, xi_v adjoints of one sample with adjoint (u_bar, v_bar, p_bar) |
| `qaws_patch_area_adjoint_cp` | c (N_i'(u) N_j(v) S_v x n + N_i(u) N_j'(v) n x S_u) at a node for global control point (i, j); zero outside the node's spans |
| `qaws_patch_line_adjoint_cp` | c times the gradient of B(b; u) for global control point (i, j): the per-line term of the gather |

## Other files in src/qaws/core

`qaws_curve_record_core.h`, `qaws_surface_record_core.h` and
`qaws_frame_core.h` are in the same directory but are not installed: they
are internal to the C runtime.
