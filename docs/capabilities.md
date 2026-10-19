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
| Clothoid | runtime | nothing to borrow: six scalars | yes | yes (146) | |
| Composite | its children | through borrowed children (142) | yes | yes (146) | |
| Subdivision | runtime | absent | absent | yes (146) | |
| Arc-length / length-matched wrappers | chain rule over the source (144) | absent | absent | malloc | |

Borrowed curves: `qaws_view.h`. Evaluation does not allocate.
`qaws_curve_prepared_evaluate_range` fills caller buffers. Disjoint ranges
may run on host workers (142: four threads, bit-identical). A borrowed
`qaws_curve` works with every operation that evaluates curves. Extract,
split and reverse go through an owned copy (142).

Arc length (`qaws_curve_compute_arc_length`): adaptive Gauss-Kronrod 7-15
on each span, to QAWS_EPSILON of the length. Its inverse
(`qaws_curve_compute_parameter_at_arc_length`) uses safeguarded Newton.
Both match Mathematica to double precision (148). Sampling, traversal and
the arc-length wrappers build their 256-entry tables from the same measure
and interpolate linearly between entries. The differentiable CDF sampler
(58, 60) keeps its fixed composite rule, because it differentiates through
that rule.

Closest points (`qaws_curve_find_closest_point_certified_2d` and `_3d`,
also behind `qaws_curve_find_closest_parameter_*`) are certified for
Bezier, rational Bezier, B-spline, NURBS and borrowed splines. They run
branch and bound over the spans' Bezier forms. The global minimum lies
within a returned gap below the distance (149: a spike between the
samples, which sampling misses). `qaws_surface_find_closest_point_certified`
(also behind `qaws_surface_find_closest_point`) does the same over the
Bezier patches of Bezier, B-spline, NURBS and borrowed surfaces.

## Surfaces on the CPU

| Family | Evaluation | Borrowed / prepared | Jets (closest point, curvature) | Derivatives w.r.t. data | Allocator |
|---|---|---|---|---|---|
| Bezier | runtime | absent | yes | yes | yes |
| B-spline | runtime | `qaws_surface_view` (143) | yes | yes | yes |
| NURBS | runtime | weights or homogeneous points (143) | yes | yes | yes |
| Ruled, Coons, extrusion, revolution, bilinear, biquadratic, loft, Gordon, offset | runtime over their curves | absent | yes | yes | yes (146) |
| Swept, pipe | analytic du, dv and the mixed derivative (Frenet-Serret; the rotation-minimizing frame law on a double-reflection frame), duu from the analytic du (152) | absent | absent | absent | yes (146) |
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
| Rotation-minimizing frame step (double reflection, `qaws_frame_core.h`) | `qaws_curve_compute_rotation_minimizing_frames_3d` | tested (147: fourth order) | compiled | compiled | written; C `if` on expressions |
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
- **Certified closest point on the other families.** Curves other than
  Bezier, rational Bezier, B-spline, NURBS and borrowed splines keep
  sampling with Newton. Surfaces other than Bezier, B-spline, NURBS and
  borrowed ones (degree at most 15) keep an 8 × 8 grid with Newton.
- **Graded surface tessellation.** `qaws_surface_tessellate_adaptive`
  (151) meets chord and angle tolerances without cracks into caller
  buffers, but has no 2:1 grading (fans of thin triangles beside much
  finer quads), and repeats shared vertices instead of welding them.
- **Third-order sweep derivatives.** Swept and pipe surfaces take duu from
  a difference of their analytic du (the torsion's rate would need the
  path's fourth derivative). Lofts and Gordon surfaces are analytic through
  their sections' jets; sections without jets fall back to differences.
  Their evaluation scratch is on the stack up to 32 sections, beyond that
  from the surface's allocator.
- **Views of the remaining families** (arc, Yuksel, trajectory,
  subdivision). Composites take borrowed children; clothoids hold no arrays. Hermite, Catmull-Rom and polynomial curves have
  views; owned objects are not rebuilt on views (their parity is tested).
