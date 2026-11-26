# Capabilities

What runs where, as the code stands, and what checks it. Numbers in
parentheses are test suites in `tests/` (`142` is `tests/142_curve_view.c`).

Status words:

- **tested**: a test suite compares results. The suite is named.
- **run on a GPU**: a shader ran on a GPU and its results were compared with
  the C cores by a tool run by hand (`examples/eval_gpu_parity.c`, Direct3D 11,
  Windows; it fails a family past 1e-4 relative to 1 + |value|).
- **compiled**: the compiler accepts the source; no results are compared.
  Checked by hand with `tests/shaders/check_shaders.ps1` (dxc `-T cs_6_0`,
  glslangValidator `-V`); not part of the build.
- **written**: the source exists; no build or tool in the repository
  compiles it.

## Precision

| | CPU | HLSL | GLSL | Halide |
|---|---|---|---|---|
| computing precision | both f32 and f64 in one library; the caller's unsuffixed names follow `QAWS_SCALAR_IS_FLOAT` (155) | float | float | float or double `Expr` |
| host data | owned objects copy into `qaws_scalar`; views read f32 or f64, packed, strided or planar, in place (142, 143) | buffers of float | buffers of float | `Func` / buffers |

## Curves on the CPU

Every family evaluates the position and D1 to D3 in either precision and has
its own suite (01 to 13). Degree limits are tested in 140.

| Family | Evaluation kernel | Borrowed view (142) | Derivatives w.r.t. data | Degree limit (140) |
|---|---|---|---|---|
| Bezier | core de Casteljau | `QAWS_VIEW_SPLINE` with n = degree + 1 | yes | 16 |
| Rational Bezier | runtime | `QAWS_VIEW_SPLINE` with weights, n = degree + 1 | yes | 47 |
| Hermite | core cubic | `QAWS_VIEW_HERMITE` | yes | 3 |
| Catmull-Rom | core cubic | `QAWS_VIEW_CATMULL_ROM` | yes | 3 |
| B-spline | runtime basis | `QAWS_VIEW_SPLINE` | yes | 31 |
| NURBS | runtime basis | `QAWS_VIEW_SPLINE`, weights or homogeneous points | yes | 31 |
| Trajectory | core cubic | `QAWS_VIEW_TRAJECTORY` | no | |
| Yuksel | runtime | `QAWS_VIEW_YUKSEL` | yes | |
| Arc | core arc | `QAWS_VIEW_ARC` | no | |
| Polynomial | core Horner | `QAWS_VIEW_POLYNOMIAL` | yes | |
| Clothoid | runtime | none (six scalars) | yes | |
| Composite | its segments | through borrowed segments | yes | |
| Subdivision | runtime | `QAWS_VIEW_SUBDIVISION` | no | |
| Arc-length and length-matched wrappers | chain rule over the source (144) | none | no | |

"Derivatives w.r.t. data" is the per-family support of `qaws_diff*.h`
([differentiation.md](differentiation.md)). Every family has an `_ex`
creator taking an allocator.

Borrowed curves (`qaws_view.h`): evaluation does not allocate;
`qaws_curve_prepared_evaluate_range` fills caller buffers, and disjoint ranges
may run on host threads (142: four threads, bit for bit the serial outputs).
A borrowed `qaws_curve` works with every operation that evaluates curves;
extract, split and reverse go through an owned copy.

Arc length (`qaws_curve_compute_arc_length`): adaptive Gauss-Kronrod 7-15 on
each span to `QAWS_EPSILON` of the length; the inverse uses safeguarded
Newton. Both match Mathematica references (148). Sampling, traversal and the
arc-length wrappers build their tables from the same measure. The
differentiable CDF sampler keeps its fixed composite rule, because it
differentiates through that rule (60).

Closest points (`qaws_curve_find_closest_point_certified_2d/3d`, also behind
`qaws_curve_find_closest_parameter_*`) are certified for Bezier, rational
Bezier, B-spline, NURBS and borrowed splines: branch and bound over the spans'
Bezier forms, the global minimum within a returned gap (149: a spike between
32 samples that sampling misses). `qaws_surface_find_closest_point_certified`
(also behind `qaws_surface_find_closest_point`) does the same over the Bezier
patches of Bezier, B-spline, NURBS and borrowed surfaces (149). Surface closest
points with their derivatives give bit-identical answers from four threads
(141).

Rotation-minimizing frames (`qaws_curve_compute_rotation_minimizing_frames_3d`)
converge at fourth order (147). Adaptive curve tessellation (150), guaranteed
curve ranges (157) and B-spline fit parameterizations (153) have their own
suites.

## Surfaces on the CPU

| Family | Evaluation | Borrowed view | Derivatives w.r.t. data |
|---|---|---|---|
| Bezier | runtime | none | yes |
| B-spline | runtime | `qaws_surface_view` (143) | yes |
| NURBS | runtime | weights or homogeneous points (143) | yes |
| Ruled, Coons, extrusion, revolution, bilinear, biquadratic, loft, Gordon, offset | runtime over their curves | none | yes |
| Swept, pipe | analytic du, dv and the mixed derivative; duu from the analytic du (152) | none | no |
| Subdivision, trimmed, T-spline | runtime | none | no |
| Borrowed (`QAWS_SURFACE_KIND_VIEW`) | prepared view, partials to third order for closest points and curvatures (143) | - | no |

Adaptive tessellation (151) and guaranteed ranges over parameter boxes (156)
have their own suites.

## Cores on each backend

The cores in `src/qaws/core/*.h` are header-only and do not allocate. They
use fixed arrays (degree at most 16) except the record cores, which read any
number of points and knots from buffers.

| Core | Runtime caller | C | HLSL | GLSL | Halide |
|---|---|---|---|---|---|
| de Casteljau (Bezier) | Bezier curves | tested (01) | run on a GPU | compiled | written (`examples/eval_bezier_halide.cpp`) |
| Rational Bezier | none | tested against the runtime (144) | run on a GPU | compiled | written |
| Cubic (Hermite, Catmull-Rom, trajectory) | Hermite, Catmull-Rom, trajectory spans; borrowed views | tested against the runtime (144) | run on a GPU | compiled | written |
| Horner (polynomial) | polynomial curves; borrowed views | tested (10) | run on a GPU | compiled | written |
| Arc | arc curves | tested (09) | run on a GPU | compiled | written |
| B-spline span and evaluation | none | tested against the runtime (144) | run on a GPU | compiled | written |
| NURBS evaluation | none | tested against the runtime (144) | run on a GPU | compiled | written |
| Bezier and B-spline surfaces | none | tested against the runtime (144) | run on a GPU | compiled | written |
| Curve records (`qaws_curve_record_core.h`) | none | tested against the runtime (144) | run on a GPU | compiled | run by Halide's JIT against the C core by hand (`examples/eval_record_halide.cpp`, fails past 1e-12 relative) |
| Surface records (`qaws_surface_record_core.h`) | none | tested against the runtime (144) | run on a GPU | compiled | run by Halide's JIT against the C core by hand (same tool) |
| Rotation-minimizing frame step (`qaws_frame_core.h`) | `qaws_curve_compute_rotation_minimizing_frames_3d`, pipe surfaces | tested (147) | compiled (`cores_check.hlsl`) | compiled (`cores_check.comp`) | written |
| B-spline derivatives w.r.t. data (`qaws_bspline_diff_core.h`) | none | tested against the runtime (57) | compiled (`diff_bspline_adjoint.hlsl`) | not checked | written |
| B-spline CDF sampling (`qaws_bspline_sampling_core.h`) | none | tested against the runtime (60) | compiled (`diff_cdf_sampling.hlsl`) | not checked | written |
| B-spline surface sampling (`qaws_bspline_surface_sampling_core.h`) | none | tested against the runtime (62) | compiled (`diff_surface_sampling.hlsl`) | not checked | written |
| Dual numbers (`qaws_dual_core.h`) | differentiation, several curve and surface families | tested (49) | compiled (through the sampling cores) | not checked | written |

"Run on a GPU" compares the core on a GPU (float, Direct3D 11) with the same
core compiled as C (double) on the same inputs, over position, D1 and D2,
and for surfaces du, dv, the normal and the second partials
(`examples/eval_gpu_parity.c` with `examples/eval_parity.hlsl`). With test 144
(C core = runtime) this closes the chain GPU = C core = runtime.

## Memory behaviour of the main operations

Every allocation goes through the `_ex` allocator argument, else the
thread's allocator, else the default ([memory.md](memory.md)). Columns:

- **Allocates:** memory taken during the call; "object" means only the
  object it returns.
- **Borrows input:** the call reads host arrays in place.
- **Caller output:** results written into caller buffers.
- **Prepared:** a reusable object for repeated calls.
- **Parallel:** calls on shared immutable inputs may run concurrently; the
  library has no threads of its own. "Executor" means the caller can pass a
  `qaws_batch_executor` ([threading.md](threading.md)).
- **Results outlive scratch:** returned objects own their memory.

| Operation | Allocates | Borrows input | Caller output | Prepared | Parallel | Results outlive scratch |
|---|---|---|---|---|---|---|
| Curve and surface evaluation (owned, borrowed, prepared ranges) | no (loft and Gordon past 32 sections: from the surface's allocator) | views (142, 143) | yes | `qaws_curve_prepare`, `qaws_curve_prepare_into` | yes (142) | n/a |
| Curve and surface construction | object | through views | object | - | yes | yes |
| Arc length, inverse, frames, certified closest points | no | any curve, borrowed included | yes | - | yes | n/a |
| Curve and surface tessellation by tolerances | no | yes | yes, counts first | - | yes | n/a |
| Sampling (`qaws_curve_sample_*`) | arc-length tables, adaptive buffers | yes | yes | - | yes | n/a |
| B-spline fitting | normal equations, parameters, knots | data read in place | object | - | yes | yes |
| Offsets (`qaws_curve_offset_2d`) | yes | yes | objects | - | yes | yes |
| Traversal | tables at creation | yes | object | the traversal | evaluate and map: yes | yes |
| Curve and surface batches | yes | yes | caller capacity and counts | prepared sets | executor | yes |

The batches, the exact closest-point batches, and seams and shatter take a cancel and progress callback
([api/progress.md](api/progress.md)). Prepared views publish
their size (`qaws_curve_prepared_size`, `qaws_surface_prepared_size`); other
operations publish no scratch-size bound.

## Limits

- **CPU and GPU comparison is manual.** Both GPU tools run by hand on
  Windows only. GLSL is compiled, never run. Halide runs the records by hand;
  its other cores are not built by any tool in the repository.
- **Certified closest point on the other families.** Curves other than
  Bezier, rational Bezier, B-spline, NURBS and borrowed splines use sampling
  with Newton. Surfaces other than Bezier, B-spline, NURBS and borrowed ones
  (degree at most 15) use an 8 x 8 grid with Newton.
- **Surface tessellation grading.** `qaws_surface_tessellate_adaptive` meets
  chord and angle tolerances without cracks into caller buffers, but has no
  2:1 grading (fans of thin triangles beside much finer quads) and repeats
  shared vertices instead of welding them.
- **Third-order sweep derivatives.** Swept and pipe surfaces take duu from a
  difference of their analytic du. Lofts and Gordon surfaces are analytic
  through their sections' jets; sections without jets fall back to
  differences.
- **SIMD kernels.** `src/qaws/batch/` has no public header; no public entry
  point and no test calls it. The curve batches use `src/qaws/simd` directly
  (B-spline basis functions a register of parameters at once).
