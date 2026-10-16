# Evaluation

Qaws provides a unified evaluation model for all curve families. Every curve can be evaluated for position and up to three orders of derivatives.

## Evaluation flags

Request what you need via a bitmask:

```c
QAWS_EVAL_FLAG_POSITION  // position on the curve
QAWS_EVAL_FLAG_D1        // first derivative (tangent direction)
QAWS_EVAL_FLAG_D2        // second derivative (curvature-related)
QAWS_EVAL_FLAG_D3        // third derivative
```

Combine flags with bitwise OR:

```c
unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
```

Only requested quantities are computed. The `valid_flags` field in the result tells you what was actually filled.

## Evaluating by global parameter

```c
qaws_eval_result_2d result;
qaws_status s = qaws_curve_evaluate_2d(&result, curve,
    0.5f, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);

// result.position.x, result.position.y
// result.d1.x, result.d1.y
```

The parameter must be within the curve's domain. Query the domain with:

```c
qaws_range range = qaws_curve_get_parameter_range(curve);
// range.min_value ... range.max_value
```

## Evaluating by span and local parameter

For piecewise curves, you can evaluate a specific span directly:

```c
qaws_eval_result_2d result;
qaws_curve_evaluate_span_2d(&result,
    curve,       0,    // span index (0 to span_count-1)
    0.5f, // local parameter within span [0, 1]
    QAWS_EVAL_FLAG_POSITION);
```

This is useful when you already know which span you're in.

## 3D evaluation

All evaluation functions have 3D variants:

```c
qaws_eval_result_3d result;
qaws_curve_evaluate_3d(&result, curve,
    0.5f, QAWS_EVAL_FLAG_POSITION);
// result.position.x, .y, .z
```

## Derivative semantics

| Flag | Meaning | Typical use |
|---|---|---|
| `D1` | First derivative dP/dt | Tangent direction, velocity |
| `D2` | Second derivative d2P/dt2 | Curvature computation, acceleration |
| `D3` | Third derivative d3P/dt3 | Jerk, torsion in 3D |

Note: derivatives are with respect to the curve parameter, not arc length. For speed-independent derivatives, use traversal with arc-length mode.

## Computing derived quantities from derivatives

**Tangent vector (unit):**
```c
// normalize(d1)
float len = sqrtf(d1.x*d1.x + d1.y*d1.y);
float tx = d1.x / len, ty = d1.y / len;
```

**Normal vector (2D):**
```c
float nx = -ty, ny = tx;
```

**Curvature (2D):**
```c
// kappa = (d1.x*d2.y - d1.y*d2.x) / (d1.x^2 + d1.y^2)^(3/2)
```

**Speed:**
```c
float speed = sqrtf(d1.x*d1.x + d1.y*d1.y);
```

## Batch evaluation

When you need to evaluate a curve at many parameters (e.g. for rendering, sampling, or LUT construction), batch evaluation amortises per-call overhead and enables SIMD acceleration where available.

```c
// Evaluate 256 points along a 2D curve in one call
qaws_scalar params[256];
qaws_eval_result_2d results[256];
for (int i = 0; i < 256; i++)
    params[i] = (qaws_scalar)i / 255.0f;

qaws_status s = qaws_curve_evaluate_batch_2d(
    results, curve, params,
    256,
    QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
```

**Signatures** (from `qaws_eval.h`):

```c
qaws_status qaws_curve_evaluate_batch_2d(
    qaws_eval_result_2d* out_results,
    qaws_curve const* curve,
    qaws_scalar const* parameters,
    unsigned int count,
    unsigned int eval_flags);

qaws_status qaws_curve_evaluate_batch_3d(
    qaws_eval_result_3d* out_results,
    qaws_curve const* curve,
    qaws_scalar const* parameters,
    unsigned int count,
    unsigned int eval_flags);
```

`out_results` is a caller-allocated array of `count` elements. On partial failure the call evaluates as many parameters as possible and returns the first error status; entries where `out_results[i].valid_flags == 0` indicate that evaluation failed for that parameter.

## Surface evaluation

Surfaces are evaluated over a (u, v) parameter domain. The single-point and batch entry points mirror the curve API.

**Single point:**

```c
qaws_surface_eval_result result;
qaws_status s = qaws_surface_evaluate(
    surface, 0.5f, 0.25f,
    QAWS_EVAL_FLAG_POSITION, &result);
// result.position.x, .y, .z
```

**Batch:**

```c
qaws_scalar u_params[64], v_params[64];
qaws_surface_eval_result results[64];
// ... fill u_params, v_params ...
qaws_status s = qaws_surface_evaluate_batch(
    results, surface, u_params, v_params,
    64, QAWS_EVAL_FLAG_POSITION);
```

`qaws_surface_evaluate_batch` takes parallel arrays of `u` and `v` parameters plus a count, just like the curve batch functions.

Surface normals can also be obtained directly:

```c
qaws_vec3 normal;
qaws_surface_get_normal(surface, u, v, &normal);
```

## Core eval headers (multi-backend architecture)

The `src/qaws/core/` directory contains **header-only inline functions** that perform the mathematical evaluation kernels. These headers use only scalar math and fixed-size arrays, with no dependency on the `qaws_curve` runtime, `malloc`, or any C standard library calls. This makes them directly usable from:

- **C / C++** — included as normal headers
- **HLSL / GLSL** — the `QAWS_INLINE`, `QAWS_TYPE_DEF`, and scalar typedefs can be redefined to match shader semantics
- **Halide** — same mechanism, compiles into Halide pipelines

Available core headers:

| Header | Purpose |
|---|---|
| `qaws_horner_core.h` | Horner-form polynomial evaluation with derivatives |
| `qaws_cubic_poly_core.h` | Cubic polynomial P(t) = at^3 + bt^2 + ct + d |
| `qaws_decasteljau_core.h` | De Casteljau Bezier evaluation |
| `qaws_rational_bezier_core.h` | Rational Bezier (weighted/homogeneous) |
| `qaws_arc_core.h` | Circular/elliptical arc evaluation |
| `qaws_bspline_basis_core.h` | B-spline basis function computation |
| `qaws_bspline_eval_core.h` | B-spline curve evaluation |
| `qaws_nurbs_eval_core.h` | NURBS evaluation |
| `qaws_surface_core.h` | Tensor-product surface evaluation |

Example — evaluate a cubic polynomial directly without a `qaws_curve`:

```c
#include "core/qaws_cubic_poly_core.h"

qaws_vec2 a = {1, 0}, b = {0, 1}, c = {2, 3}, d = {0, 0};
qaws_eval_2d r = qaws_cubic_eval_2d(a, b, c, d, 0.5f,
    0x1 /* position only */);
```

## Prepare + eval architecture

The **prepare** functions (`qaws_prepare.h`) bridge high-level curve descriptions and the low-level core eval headers. They pre-compute flat coefficient arrays from user-facing descriptors (control points, tangents, knots, weights, etc.) so that core headers can evaluate from those arrays with no further interpretation.

This is the mechanism that enables non-C backends: on the CPU, `qaws_prepare_*` runs once to produce the coefficient arrays, which are then uploaded to a GPU buffer or passed to a Halide pipeline where the core eval kernels consume them.

**Workflow:**

1. Call a `qaws_*_prepare_*` function with the high-level geometry.
2. Receive flat coefficient arrays (e.g. `a, b, c, d` for cubics, weighted CPs for rational Beziers).
3. Feed those arrays to the matching `core/` eval function on any backend.

**Available prepare functions:**

| Prepare function | Outputs |
|---|---|
| `qaws_hermite_prepare_2d/3d` | a, b, c, d coefficients per span |
| `qaws_catmull_rom_prepare_2d/3d` | a, b, c, d + span count |
| `qaws_trajectory_prepare_2d/3d` | a, b, c, d + span count |
| `qaws_bezier_prepare_deriv_2d/3d` | derivative control points (D1, D2, D3) |
| `qaws_rational_bezier_prepare_2d/3d` | homogeneous weighted CPs |
| `qaws_arc_prepare_2d/3d` | axis vectors + theta range |

**Example — Hermite curve on GPU:**

```c
// CPU side: prepare coefficients
qaws_vec3 a[N], b[N], c[N], d[N];
qaws_hermite_prepare_3d(a, b, c,
    d, points, tangents, point_count);
// Upload a, b, c, d arrays to GPU buffer

// GPU side (HLSL/GLSL): evaluate using core header
// qaws_cubic_eval_3d(a[span], b[span], c[span], d[span], t, flags)
```

## Inline curves (stack-allocated evaluation)

For cases where you want full `qaws_curve` API compatibility but cannot use heap allocation (e.g. embedded systems, hot loops, shader-adjacent code), **inline curves** provide stack-allocated initialization.

```c
#include "qaws_inline.h"

qaws_curve_inline ic;
qaws_vec3 pts[] = {{0,0,0}, {1,2,3}, {4,5,6}, {7,8,9}};
qaws_bezier_desc desc = { QAWS_DIMENSION_3D, 3, pts, 4 };
qaws_curve_init_bezier_inline(&ic, &desc);

// Use via &ic.curve — works with all evaluation functions
qaws_eval_result_3d result;
qaws_curve_evaluate_3d(&result, &ic.curve,
    0.5f, QAWS_EVAL_FLAG_POSITION);
```

**Available inline constructors:**

| Function | Curve family |
|---|---|
| `qaws_curve_init_bezier_inline` | Bezier (single span, degree <= 7) |
| `qaws_curve_init_polynomial_inline` | Polynomial (single span, degree <= 7) |

Constraints:
- Single span only.
- Maximum degree 7 (8 control points).
- 2D or 3D.
- The `qaws_curve_inline` struct is self-contained and can be copied with `memcpy`. Do **not** call `qaws_curve_destroy` on it.
