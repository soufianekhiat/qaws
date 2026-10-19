# Capabilities

What runs where, as the code stands, and what checks it. This page replaces
`feature_per_backends.md`, which listed every core as available on every
backend.

Status words:

- **tested**: a test compares results. The test is named.
- **run on a GPU**: a shader ran on a GPU, and its results were compared with
  the CPU by a tool run by hand. The largest difference is given (below).
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
| Rational Bezier | runtime | NURBS view with n = degree + 1 (142) | yes | yes (146) | 47 (140) |
| Hermite | runtime | `QAWS_VIEW_HERMITE` (142) | yes | yes | 3 |
| Catmull-Rom | runtime | `QAWS_VIEW_CATMULL_ROM` (142) | yes | yes | 3 |
| B-spline | internal basis | `qaws_curve_view`, prepared, borrowed (142) | yes | yes | 31 (140) |
| NURBS | internal basis | weights or homogeneous points (142) | yes | yes | 31 (140) |
| Trajectory | runtime | absent | absent | yes | |
| Yuksel | runtime | absent | yes | yes (146) | |
| Arc | core arc | absent | absent | yes (146) | |
| Polynomial | core Horner | `QAWS_VIEW_POLYNOMIAL` (142) | yes | yes (146) | |
| Clothoid | runtime | absent | yes | yes (146) | |
| Composite | its children | absent | yes | yes (146) | |
| Subdivision | runtime | absent | absent | yes (146) | |
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
| Ruled, Coons, extrusion, revolution, bilinear, biquadratic, loft, Gordon, offset | runtime over their curves | absent | yes | yes | yes (146) |
| Swept, pipe | finite differences | absent | absent | absent | yes (146) |
| Fillet, subdivision, trimmed, T-spline | runtime | absent | absent | absent | yes (146) |
| Borrowed (`QAWS_SURFACE_KIND_VIEW`) | prepared view | - | to third order (143) | absent | yes |

## Cores on each backend

The cores in `src/qaws/core/*.h` are header-only and do not allocate. They
use fixed arrays: degree at most 16, at most 34 knots, D2 at most for
evaluation.

| Core | Runtime caller | C | HLSL | GLSL | Halide |
|---|---|---|---|---|---|
| de Casteljau (Bezier) | Bezier curves | tested (01) | run on a GPU, 1.7e-6 (eval_gpu_parity) | compiled | written (`examples/eval_bezier_halide.cpp`) |
| Rational Bezier | none | tested against the runtime (144) | run on a GPU, 1.9e-6 (eval_gpu_parity) | compiled | written |
| Cubic (Hermite, Catmull-Rom, trajectory) | none | tested against the runtime (144) | run on a GPU, 1.3e-7 (eval_gpu_parity) | compiled | written |
| Horner (polynomial) | polynomial curves | tested (10) | run on a GPU, 2.1e-7 (eval_gpu_parity) | compiled | written |
| Arc | arc curves | tested (09) | run on a GPU, 2e-5 (eval_gpu_parity) | compiled | written |
| B-spline span and evaluation | none | tested against the runtime (144) | run on a GPU, 7.4e-7 (eval_gpu_parity) | compiled | written; needs the span at pipeline-build time |
| NURBS evaluation | none | tested against the runtime (144) | run on a GPU, 2.1e-6 (eval_gpu_parity) | compiled | written; same span limit |
| Bezier surface | none | tested against the runtime (144) | run on a GPU, 2.2e-6 (eval_gpu_parity) | compiled | written |
| B-spline surface | none | tested against the runtime (144) | run on a GPU, 1.9e-6 (eval_gpu_parity) | compiled | written; the normal block uses C `if` on expressions |
| Curve records: B-spline and NURBS read from buffers, any length (`qaws_curve_record_core.h`) | none | tested against the runtime (144) | run on a GPU, 4.9e-7 (eval_gpu_parity) | compiled | absent |
| Surface records: B-spline and NURBS surfaces read from buffers, any number of points (`qaws_surface_record_core.h`) | none | tested against the runtime (144) | run on a GPU, 7.6e-6 (eval_gpu_parity) | compiled | absent |
| B-spline derivatives w.r.t. data | none | tested (57) | compiled (`diff_bspline_adjoint.hlsl`) | not checked | written |
| B-spline sampling (CDF) | none | tested (60) | compiled (`diff_cdf_sampling.hlsl`) | not checked | written; C `if` on expressions |
| Surface functionals | none | tested (61) | compiled (`diff_surface_functionals.hlsl`) | not checked | written |
| Surface sampling | none | tested (62) | compiled (`diff_surface_sampling.hlsl`) | not checked | written; C `if` on expressions |
| Repulsive pair kernels | none | tested (119) | run on a GPU once by hand, parity 1.1e-4 (`examples/repulsive_gpu_parity.c`) | not checked | written; C `if` on expressions |
| Dual numbers | differentiation | tested (49) | compiled (through the sampling and functional cores) | not checked | written |

"Compiled" says the shader source is accepted. "Run on a GPU" gives the
largest difference, relative to 1 + |value|, between the core on a GPU
(float, Direct3D 11, fxc cs_5_0) and the same core compiled as C (double)
on the same inputs, over 256 parameters per family: position, D1 and D2,
and for surfaces du, dv, the normal and the second partials.
`examples/eval_gpu_parity.c` with `examples/eval_parity.hlsl` measures it
by hand on Windows (2026-10-19). With test 144 (C core = runtime) this
closes the chain GPU = C core = runtime. The repulsive kernels have their
own tool, `examples/repulsive_gpu_parity.c`.

## Missing

- **Halide spans at run time.** Curves and surfaces read from buffers have a
  shader interface (`qaws_curve_record`, `qaws_surface_record`): a bounded
  span search on the knot buffers, and only the de Boor window loaded. Halide
  still needs the span when the pipeline is built.
- **CPU and GPU comparison in a build step.** Both GPU tools run by hand on
  Windows only. GLSL is compiled, never run. Halide is never built.
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
- **Sweeps and lofts with analytic derivatives.** Swept and pipe surfaces use
  finite differences. Loft and Gordon evaluation scratch is on the stack up
  to 32 sections, beyond that from the surface's allocator.
- **Views of the remaining families** (arc, clothoid, Yuksel, trajectory,
  subdivision, composite). Hermite, Catmull-Rom and polynomial curves have
  views; owned objects are not rebuilt on views (their parity is tested).
- **The allocator through the other entry points.** Every curve and surface
  constructor takes one, and curves derived from a source curve use its
  allocator (test 146). Operation scratch, sampling, traversal, batches and
  repulsive still use malloc.
