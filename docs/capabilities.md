# Capabilities

What runs where, as the code stands, and what checks it. This page replaces
`feature_per_backends.md`, which listed every core as available on every
backend.

Status words:

- **tested**: a test compares results. The test is named.
- **compiled**: the compiler accepts the source. No results are compared.
  Checked by `tests/shaders/check_shaders.ps1` (dxc `-T cs_6_0` and
  glslangValidator `-V`, Vulkan SDK 1.4.363, 2026-10-19). This is a manual
  step, not part of the build.
- **written**: the source exists but has never been built.
- **absent**: the feature does not exist.

## Precision

| | CPU | HLSL | GLSL | Halide |
|---|---|---|---|---|
| computing precision | `qaws_scalar`: float or double per build (`QAWS_SCALAR_IS_FLOAT`) | float | float | float or double `Expr` |
| host data | packed `qaws_scalar` (owned objects); f32 or f64, strided or planar, borrowed (views, tests 142 and 143) | buffers of float | buffers of float | `Func` / buffers |

## Curves on the CPU

Every curve family evaluates the position and D1 to D3, in either
precision, and is tested by its own suite (01 to 10 and onwards). The other
columns differ by family.

| Family | Evaluation kernel | Borrowed / prepared | Derivatives w.r.t. data | Allocator | Degree limit |
|---|---|---|---|---|---|
| Bezier | core de Casteljau (to D2) + runtime D3 | view with n = degree + 1 (142) | yes | yes | 16 (140) |
| Rational Bezier | runtime | NURBS view with n = degree + 1 (142) | yes | malloc | 47 (140) |
| Hermite | runtime | absent | yes | yes | 3 |
| Catmull-Rom | runtime | absent | yes | yes | 3 |
| B-spline | internal basis | `qaws_curve_view`, prepared, borrowed (142) | yes | yes | 31 (140) |
| NURBS | internal basis | weights or homogeneous points (142) | yes | yes | 31 (140) |
| Trajectory | runtime | absent | absent | yes | |
| Yuksel | runtime | absent | yes | malloc | |
| Arc | core arc | absent | absent | malloc | |
| Polynomial | core Horner | absent | yes | malloc | |
| Clothoid | runtime | absent | yes | malloc | |
| Composite | its children | absent | yes | malloc | |
| Subdivision | runtime | absent | absent | malloc | |
| Arc-length / length-matched wrappers | chain rule over the source (144) | absent | absent | malloc | |

Borrowed curves: `qaws_view.h`. Evaluation does not allocate.
`qaws_curve_prepared_evaluate_range` fills caller buffers. Disjoint ranges
may run on host workers (142: four threads, bit-identical). A borrowed
`qaws_curve` works with every operation that evaluates curves. Extract,
split and reverse go through an owned copy (142).

## Surfaces on the CPU

| Family | Evaluation | Borrowed / prepared | Jets (closest point, curvature) | Derivatives w.r.t. data | Allocator |
|---|---|---|---|---|---|
| Bezier | runtime | absent | yes | yes | yes |
| B-spline | runtime | `qaws_surface_view` (143) | yes | yes | yes |
| NURBS | runtime | weights or homogeneous points (143) | yes | yes | yes |
| Ruled, Coons, extrusion, revolution, bilinear, biquadratic, loft, Gordon, offset | runtime over their curves | absent | yes | yes | malloc |
| Swept, pipe | finite differences | absent | absent | absent | malloc |
| Fillet, subdivision | runtime | absent | absent | absent | malloc |
| Borrowed (`QAWS_SURFACE_KIND_VIEW`) | prepared view | - | to third order (143) | absent | yes |

## Cores on each backend

The cores in `src/qaws/core/*.h` are header-only and do not allocate. They
use fixed arrays: degree at most 16, at most 34 knots, D2 at most for
evaluation.

| Core | Runtime caller | C | HLSL | GLSL | Halide |
|---|---|---|---|---|---|
| de Casteljau (Bezier) | Bezier curves | tested (01) | compiled | compiled | written (`examples/eval_bezier_halide.cpp`) |
| Rational Bezier | none | tested against the runtime (144) | compiled | compiled | written |
| Cubic (Hermite, Catmull-Rom, trajectory) | none | compiled | compiled | compiled | written |
| Horner (polynomial) | polynomial curves | tested (10) | compiled | compiled | written |
| Arc | arc curves | tested (09) | compiled | compiled | written |
| B-spline span and evaluation | none | compiled | compiled | compiled | written; needs the span at pipeline-build time |
| NURBS evaluation | none | compiled | compiled | compiled | written; same span limit |
| Bezier surface | none | tested against the runtime (144) | compiled | compiled | written |
| B-spline surface | none | compiled | compiled | compiled | written; the normal block uses C `if` on expressions |
| B-spline derivatives w.r.t. data | none | tested (57) | compiled (`diff_bspline_adjoint.hlsl`) | not checked | written |
| B-spline sampling (CDF) | none | tested (60) | compiled (`diff_cdf_sampling.hlsl`) | not checked | written; C `if` on expressions |
| Surface functionals | none | tested (61) | compiled (`diff_surface_functionals.hlsl`) | not checked | written |
| Surface sampling | none | tested (62) | compiled (`diff_surface_sampling.hlsl`) | not checked | written; C `if` on expressions |
| Repulsive pair kernels | none | tested (119) | run on a GPU once by hand, parity 1.1e-4 (`examples/repulsive_gpu_parity.c`) | not checked | written; C `if` on expressions |
| Dual numbers | differentiation | tested (49) | compiled (through the sampling and functional cores) | not checked | written |

"Compiled" says the shader source is accepted. Only the repulsive kernels
have been run on a GPU and compared with the CPU.

## Missing

- **Shader ABI for B-splines, NURBS and surfaces.** Today a shader copies a
  knot window into a 34-element array. Planned: a header per curve in a
  structured buffer, a bounded span search on the knot buffer, and the
  de Boor window gathered from buffers. Halide would take the span as a
  runtime `Expr`.
- **CPU and GPU result comparison.** Every eval core except the repulsive
  kernels is checked by compiling only.
- **Parallel-transport frames.** The kernel is private to pipe surfaces.
  Planned: a public kernel (double reflection) and a prepared frame table.
- **Arc length and its inverse.** Two implementations disagree: a 5-point
  Gauss rule with a 256-entry table, and the CDF sampler. Planned: one
  error-controlled measure with a safeguarded Newton inverse.
- **Certified closest point.** Today curves use 32 samples and
  Gauss-Newton; surfaces use an 8x8 grid and Newton. Both can miss the
  global minimum. Planned: subdivision with bounds.
- **Adaptive tessellation** with chord and angle tolerances into caller
  buffers. Today the quadtree uses an absolute threshold.
- **Sweeps and lofts with analytic derivatives** and without allocation
  per evaluation. Swept and pipe surfaces use finite differences. The loft
  and Gordon fallbacks allocate.
- **Views of the other families** (Hermite, Catmull-Rom, polynomial, arc
  and others), and owned objects rebuilt on views.
- **The allocator through every entry point** marked malloc above, and
  through operations, sampling, traversal, batches and repulsive.
