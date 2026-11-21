# Precision and Robustness

This page covers four topics:

- how qaws offers single and double precision side by side;
- the tolerances its functions use;
- what the floating-point API promises;
- how that compares to the exact kernel.

Headers: `qaws_platform.h`, `qaws_precision.h`, `qaws_names.h`, `qaws_both.h`. API page:
[api/precision.md](api/precision.md).

## Both precisions in one library

The C library always contains both precisions. Every function and type that carries a scalar exists twice, as
`name_f32` and `name_f64`. For example, `qaws_curve_evaluate_2d_f32` takes `qaws_f32` (float) and
`qaws_curve_evaluate_2d_f64` takes `qaws_f64` (double). The two sets do not collide and can be linked into one
program.

The unsuffixed names in the documentation are macros:

- `qaws_names.h` (generated) maps each precision-bearing name to `QAWS_SFX(name)`.
- `qaws_precision.h` defines `QAWS_SFX` as `_f32` or `_f64` according to `QAWS_SCALAR_IS_FLOAT`.

| `QAWS_SCALAR_IS_FLOAT` | unsuffixed names mean | `qaws_scalar` | `QAWS_EPSILON` | `QAWS_LITERAL(1.0)` |
|---|---|---|---|---|
| 1 (the default) | `_f32` | `float` | `1e-6f` | `1.0f` |
| 0 | `_f64` | `double` | `1e-12` | `1.0` |

`QAWS_SCALAR_IS_FLOAT` only chooses what the caller's unsuffixed names mean; it does not decide what is compiled into
the library. The CMake option of the same name sets the default for consumers of the package. A translation unit can
also define it before including `qaws.h`.

- **Both in one file.** `qaws_both.h` declares both precisions in one translation unit. The unsuffixed names keep the
  caller's precision, and the suffixed names of the other precision are available too.
- **Switching inside one file.** `qaws_precision.h` has no include guard, so including it again after changing
  `QAWS_SCALAR_IS_FLOAT` switches the meaning of `qaws_scalar`, the scalar macros (`QAWS_SQRT`, `QAWS_FABS`, ...) and
  the unsuffixed names.
- **Precision-free parts.** Some parts are declared and defined only once: the allocator and arena
  (`qaws_allocator.h`), the batch executor, `qaws_status` and the backend macros. A default allocator or an arena
  therefore serves the f32 and f64 functions alike.

### Shader and array backends

`qaws_platform.h` detects the backend: C, HLSL, GLSL or Halide. HLSL and GLSL are float only. Halide follows
`QAWS_SCALAR_IS_FLOAT`. See [backends.md](backends.md).

## Precision of data that is not qaws_scalar

- **Views** (`qaws_array_view`) hold f32 or f64 data whatever `qaws_scalar` is. Values are converted as they are read
  (see [views_and_batches.md](views_and_batches.md)).
- **The exact kernel** always takes `double` inputs, in both precisions. Its results are converted to `double` once,
  correctly rounded. See [exact_kernel.md](exact_kernel.md).

## Tolerances

There is no global tolerance setting. Each function documents its own, and most defaults scale with the extent of the
input so that results do not depend on units. The f32 forms use looser defaults. Passing 0 selects the default.

| Where | Default | f32 |
|---|---|---|
| `QAWS_EPSILON` (generic small value, e.g. B-rep degenerate edges) | 1e-12 | 1e-6 |
| Extract of a Yuksel curve; transform of a clothoid under a non-similarity | fitted to 1e-9 of the extent | extract: 1e-5 |
| Point location in a region (`ON` band) | 1e-10 of the extent | 1e-5 |
| Clip vertex merging | 6.4e-9 of the extent | 1e-4 |
| Surface seams: gap at a seam point | 1e-9 of the extent | 1e-5 |
| Surface seams: grid offset that resolves coincidences | 1e-8 of the extent | 1e-5 |
| Region seams: gap at a crossing / chord flatness | 1e-12 / 1e-3 of the extent | 1e-6 / 1e-3 |
| Curve batches: flatness | 2^-10 of the scene extent | same |
| Surface batches: flatness | 2^-9 of the scene extent | same |

Functions that take an explicit tolerance (`qaws_curve_offset_cubic_2d`, `qaws_curve_merge_chain`, adaptive
tessellation, polyline simplification) use it in the curve's own units.

## The robustness model of the floating-point API

The float API (f32 and f64) computes in floating point. It does not certify its answers.

- **Exactness by construction.** Results that are exact in the family's own arithmetic stay exact up to rounding.
  Examples are extract, split, conversions and transforms, which produce curves that lie on their source (see
  [operations.md](operations.md)).
- **Bounded by sampling.** Searches (closest points, intersections, seams, batches) start from a flattening or a grid
  and refine with Newton's method. Their accuracy at a found answer is set by the solver's tolerance. Whether an answer
  is found at all is set by the flattening. Contacts that are tangential, branches closer than the flatness, and
  features finer than the grid can be missed. Each header states which of these it does not resolve.
- **Decided degeneracies.** Where a decision must be consistent rather than exact, the float code fixes one answer for
  every query. Surface shatter does this by moving each operand by a small offset of its own. Region shatter splits
  collinear chords against each other. A few predicates use exact arithmetic, such as the orientation test in
  `qaws_polyline_trim_collinear_2d`.
- **Failures are reported.** Failures that a function detects come back as a status (`NUMERICAL_FAILURE`,
  `NOT_CONVERGED`, `DEGENERATE_CURVE`, ...). Check the status before using the output. See
  [error_handling.md](error_handling.md).

f32 gives about 7 significant digits, which is enough for display, tessellation and animation. f64 gives about 16,
which is what modelling queries, Boolean-style cuts and fitting need.

## The exact kernel

`qaws_exact.h` (built with `QAWS_ENABLE_EXACT`, on by default) is a separate model with a different promise:

- **Predicates are certified.** A predicate's sign is never wrong. The kernel first evaluates in f64 with a proven
  error bound. Only when that bound does not prove the sign does it fall back to exact multi-limb integer arithmetic.
- **Quantized once, then exact.** Curves and surfaces are rounded once onto integer lattices that the caller chooses.
  After that, evaluation does not round. Queries return enclosures, certified results, or a status that says
  certification failed.

Use the float API for speed and the exact kernel where a wrong sign or a missed root is not acceptable. The exact
kernel's design, lattices, statuses and certified batches are described in [exact_kernel.md](exact_kernel.md).
