# Qaws Overview

**Qaws** (Arabic for arc) is a dependency-free C11 library for parametric
geometry in 2D and 3D: curves and surfaces, their evaluation, sampling,
traversal and inspection, operations on them, 2D regions and Booleans,
surface seams and shatter, a B-rep shell, an exact (certified) kernel,
differentiation with respect to the geometry's data, and motions.

Every precision-bearing name exists as `name_f32` and `name_f64`; these pages
use the unsuffixed names (see [precision.md](precision.md)).

## What the library models

| Area | Headers | Concept page |
|---|---|---|
| Curves: 13 families created from descriptors (Bezier, rational Bezier, Hermite, Catmull-Rom, B-spline, NURBS, trajectory, Yuksel, arc, polynomial, clothoid, subdivision, composite), interpolating splines, polylines | `qaws_bezier.h` ... `qaws_composite.h`, `qaws_interpolating.h`, `qaws_path.h` | [curve_families.md](curve_families.md) |
| Surfaces: Bezier, B-spline, NURBS, bilinear, biquadratic, swept, ruled, Coons, extrusion, revolution, pipe, loft, Gordon, offset, trimmed, T-spline, subdivision | `qaws_surface_*.h` | [surfaces.md](surfaces.md) |
| Evaluation: position and derivatives at a parameter, per span, in batches; portable cores | `qaws_eval.h`, `qaws_prepare.h`, `qaws_inline.h`, `core/*.h` | [evaluation.md](evaluation.md) |
| Sampling and tessellation | `qaws_sampling.h`, `qaws_inspect.h` | [sampling.md](sampling.md) |
| Traversal: motion along a curve by distance or time | `qaws_traversal.h` | [traversal.md](traversal.md) |
| Inspection: arc length, bounds, closest points (certified), frames, curvature, enclosures | `qaws_inspect.h` | [inspection.md](inspection.md) |
| Operations: split, extract, join, offsets, fillets, conversion, export and import | `qaws_operations.h`, `qaws_convert.h`, `qaws_export.h`, `qaws_import.h` | [operations.md](operations.md) |
| Borrowed views and batches: host-owned arrays evaluated in place; many curves or surfaces at once | `qaws_view.h`, `qaws_curve_batch.h`, `qaws_surface_batch.h`, `qaws_batch_executor.h` | [views_and_batches.md](views_and_batches.md) |
| 2D paths, regions, clipping and Booleans on curves | `qaws_path.h`, `qaws_clip.h`, `qaws_boolean_2d.h` | [regions_2d.md](regions_2d.md) |
| Seams and shatter of curved surfaces | `qaws_shatter.h` | [shatter.md](shatter.md) |
| B-rep shell: faces, edges, vertices over borrowed surfaces | `qaws_brep.h` | [brep.md](brep.md) |
| Exact kernel: certified predicates, exact curves and surfaces | `qaws_exact.h` | [exact_kernel.md](exact_kernel.md) |
| Differentiation with respect to control points, knots, weights and parameters | `qaws_diff*.h` | [differentiation.md](differentiation.md) |
| Motions over t in [0, 1]: rigid (poses, along a curve) and per vertex | `qaws_motion.h` | [api/motion.md](api/motion.md) |
| Memory: allocators, thread and scratch allocators, arena | `qaws_alloc.h`, `qaws_allocator.h` | [memory.md](memory.md) |

`qaws.h` includes every public header. The function-by-function reference is
[api/index.md](api/index.md).

## Design principles

1. **C first.** Plain structs, enums and typedefs, opaque handles
   (`qaws_curve*`, `qaws_surface*`, ...), explicit ownership: a descriptor is
   read during the call, an object is destroyed by its own `_destroy`.
2. **Status codes.** Every fallible function returns `qaws_status`
   ([error_handling.md](error_handling.md)).
3. **One evaluation model.** Position and D1 to D3 for every curve family;
   position, partials to second order and the normal for every surface.
   Piecewise curves are the normal case.
4. **Geometry and motion are separate.** A curve is shape; a traversal or a
   motion says how something moves along it.
5. **Immutable objects.** Curves and surfaces do not change after creation:
   any number of threads may read one at once. The library has no threads of
   its own; batched work can be fanned out through a caller's
   `qaws_batch_executor` ([threading.md](threading.md)).
6. **Both precisions in one library.** f32 and f64 forms of every
   precision-bearing function and type; `QAWS_SCALAR_IS_FLOAT` (default 1)
   only picks what the unsuffixed names mean ([precision.md](precision.md)).
7. **Caller-controlled memory.** `_ex` constructors take a `qaws_allocator`;
   otherwise the thread's allocator, else the default. Many queries write into
   caller buffers and report counts ([memory.md](memory.md)).
8. **Portable cores.** The evaluation kernels in `src/qaws/core/` are
   header-only and compile as C, HLSL, GLSL and Halide
   ([backends.md](backends.md)).
9. **Features, not implementations, in names.** Internals (caches, solvers,
   basis code, the mesh machinery) stay behind the opaque handles.

## How the parts fit

- A **descriptor** goes into `qaws_curve_create_*` or `qaws_surface_create_*`
  and an owned object comes out. A **view** (`qaws_view.h`) instead wraps
  host arrays without copying them; a borrowed curve or surface works with the
  same evaluation, sampling and inspection calls.
- **Evaluation** is the base every other module uses: sampling, traversal,
  inspection, operations, Booleans and differentiation all evaluate through
  `qaws_curve_evaluate_*` / `qaws_surface_evaluate` or the same span code.
- **Composite objects borrow or own their inputs** as each header states:
  a composite curve owns its segments; swept, ruled, Coons, loft, Gordon,
  offset, trimmed surfaces and B-rep shells borrow theirs; traversals borrow
  their curves.
- **Batches** (`qaws_curve_batch.h`, `qaws_surface_batch.h`, the exact kernel)
  take many inputs at once and an optional executor; results do not depend on
  the executor.
- The **exact kernel** is C only and built with `QAWS_ENABLE_EXACT` (CMake,
  default ON); its predicates decide signs exactly on double inputs.

## Quick start

```c
#include "qaws.h"

qaws_vec2 points[] = { {0,0}, {0,1}, {1,1}, {1,0} };
qaws_bezier_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = points,
    .control_point_count = 4
};

qaws_curve* curve = NULL;
if (qaws_curve_create_bezier(&curve, &desc) == QAWS_STATUS_OK) {
    qaws_eval_result_2d r;
    qaws_curve_evaluate_2d(&r, curve, (qaws_scalar)0.5,
        QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
    /* r.position = (0.5, 0.75), r.d1 = first derivative at t = 0.5 */
    qaws_curve_destroy(curve);
}
```

## Building

Sharpmake is the primary build (Visual Studio solution); CMake mirrors the
public API for other platforms, installs and packaging; the library can also
be dropped into any project or used as a two-file amalgamation. See
[building.md](building.md) and [integration.md](integration.md).

## Repository layout

```
src/qaws/            public headers (qaws_*.h, qaws.h) and the C implementation
src/qaws/core/       portable header-only evaluation kernels (C, HLSL, GLSL, Halide)
src/qaws/exact/      exact kernel implementation (internal headers)
src/qaws/internal/   internal numerical machinery (not public)
src/qaws/simd/       SIMD abstraction (AVX2, SSE2, NEON, scalar)
src/qaws/batch/      SIMD evaluation kernels (built with SIMD on)
tests/               test suites (one per numbered file) and shader checks
examples/            C showcases, HLSL/GLSL shaders, Halide generators
bench/               benchmarks
tools/               API naming tool (qaws_names.h generator and checks)
buildsystem/         Sharpmake sources, bootstrap scripts, amalgamate.cmake, CMake templates
docs/                concept pages (this directory)
docs/api/            API reference, one page per header area
working_dir/         test and example output
```

## Other concept pages

[building.md](building.md), [integration.md](integration.md),
[curve_families.md](curve_families.md), [surfaces.md](surfaces.md),
[evaluation.md](evaluation.md), [sampling.md](sampling.md),
[traversal.md](traversal.md), [inspection.md](inspection.md),
[operations.md](operations.md), [views_and_batches.md](views_and_batches.md),
[regions_2d.md](regions_2d.md), [shatter.md](shatter.md), [brep.md](brep.md),
[exact_kernel.md](exact_kernel.md), [differentiation.md](differentiation.md),
[precision.md](precision.md), [memory.md](memory.md),
[threading.md](threading.md), [error_handling.md](error_handling.md),
[backends.md](backends.md), [capabilities.md](capabilities.md).
