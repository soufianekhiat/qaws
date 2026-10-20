# Qaws Overview

**Qaws** (Arabic: قوس, "arc") is a dependency-free C11 library for curve and surface evaluation, traversal, and inspection in 2D and 3D.

## What Qaws does

Qaws lets you:

- **Define** curves from control points, knots, weights, or keyframes
- **Evaluate** position and derivatives at any parameter -- on CPU, GPU, or DSP
- **Sample** curves into point buffers for rendering or export
- **Traverse** curves with constant-speed, acceleration, or time-based motion
- **Inspect** curve properties: tangent, curvature, torsion, arc length, bounds, closest point, intersections, and more
- **Operate** on curves: split, join, offset, reparameterize, reverse
- **Convert** between curve families and export to SVG, polyline, or B-spline fits
- **Evaluate surfaces**: Bezier, B-spline, NURBS, swept, and ruled patches

All curve families share one evaluation interface. A Bezier curve and a NURBS curve are evaluated the same way.

## Design principles

1. **C first** -- plain structs, enums, typedefs, explicit memory ownership, opaque handles
2. **Geometry and motion are separate** -- curves define shape; traversal objects define motion
3. **One evaluation model** -- position + D1/D2/D3 derivatives, uniform across all families
4. **Piecewise is normal** -- multi-span curves are first-class, not a special case
5. **Multi-backend** -- core evaluation math is isolated in portable headers reusable from C, HLSL, GLSL, and Halide
6. **Public API is stable** -- internal caches, solvers, and basis math are hidden

## Supported curve families

| Family | Header | Continuity | Interpolating | Closed | 3D |
|---|---|---|---|---|---|
| Bezier | `qaws_bezier.h` | C0 | endpoints only | no | yes |
| Rational Bezier | `qaws_rational_bezier.h` | C0 | endpoints only | no | yes |
| Hermite | `qaws_hermite.h` | C1 | yes | no | yes |
| Catmull-Rom | `qaws_catmull_rom.h` | C1 | yes | yes | yes |
| B-Spline | `qaws_bspline.h` | up to C(p-1) | no | yes | yes |
| NURBS | `qaws_nurbs.h` | up to C(p-1) | no | yes | yes |
| Trajectory | `qaws_trajectory.h` | C1 | yes | yes | yes |
| Yuksel | `qaws_yuksel.h` | C2 | yes | yes | yes |
| Arc | `qaws_arc.h` | C2 | endpoints | no | yes |
| Polynomial | `qaws_polynomial.h` | C-inf | no | no | yes |
| Clothoid | `qaws_clothoid.h` | C2 | no | no | yes |
| Subdivision | `qaws_subdivision.h` | varies | varies | yes | yes |
| Composite | `qaws_composite.h` | user-defined | varies | yes | yes |
| Inline variants | `qaws_inline.h` | (matches family) | (matches family) | (matches family) | yes |

Inline curves (`qaws_inline.h`) are stack-allocated versions of the above families -- no heap, no handle, no destroy call.

## Surfaces

| Surface type | Header |
|---|---|
| Bezier patch | `qaws_surface_bezier.h` |
| B-spline patch | `qaws_surface_bspline.h` |
| NURBS patch | `qaws_surface_nurbs.h` |
| Swept surface | `qaws_surface_swept.h` |
| Ruled surface | `qaws_surface_ruled.h` |

Surfaces share a common evaluation model (`qaws_surface.h`) with position and partial derivatives at any (u, v) parameter.

## Multi-backend support

Qaws separates portable evaluation math from the C runtime API:

- **C runtime** (`src/qaws/*.c`) -- full API with memory management, traversal, inspection, and operations
- **Core eval headers** (`src/qaws/core/`) -- standalone, dependency-free headers containing the raw math (de Casteljau, Horner, B-spline basis, NURBS, arc, surface evaluation). These can be included directly from HLSL, GLSL, Halide, or any C-like language.
- **Bootstrap headers** -- `qaws_hlsl.h`, `qaws_glsl.h`, `qaws_halide.h` provide type and macro scaffolding so core headers compile in each target language
- **SIMD batch evaluation** (`src/qaws/simd/`, `src/qaws/batch/`) -- AVX2, SSE2, and NEON backends for evaluating many parameters in one call

Example GPU shaders and Halide generators live in `examples/`.

## Operations, conversion, and inspection

**Operations** (`qaws_operations.h`): split, join, offset, reparameterize, reverse.

**Conversion** (`qaws_convert.h`): Hermite-to-Bezier, Catmull-Rom-to-Bezier, Bezier-to-B-spline, B-spline-to-NURBS, degree elevation, degree reduction.

**Export** (`qaws_export.h`): SVG path string, polyline sampling, B-spline fitting.

**Import** (`qaws_import.h`): polyline to Catmull-Rom or Trajectory.

**Inspection** (`qaws_inspect.h`): tangent, curvature, torsion, speed, normal, Frenet frame, inflection points, extrema, curvature comb, winding number, intersection, self-intersection, arc length, bounds, closest point.

**Prepare** (`qaws_prepare.h`): pre-compute coefficients and lookup tables for repeated evaluation -- avoids redundant basis/knot work in hot loops.

## Custom allocators

All `_create` functions have `_ex` variants that accept a `qaws_allocator` for custom memory management. The default path uses `malloc`/`free`.

## Configurable scalar type

Qaws supports `float` (f32) or `double` (f64) via the `QAWS_SCALAR_IS_FLOAT` compile option. Default is `float`.

## Quick start

```c
#include "qaws.h"

// Create a cubic Bezier
qaws_vec2 points[] = { {0,0}, {0,1}, {1,1}, {1,0} };
qaws_bezier_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = points,
    .control_point_count = 4
};

qaws_curve* curve = NULL;
qaws_curve_create_bezier(&curve, &desc);

// Evaluate at midpoint
qaws_eval_result_2d result;
qaws_curve_evaluate_2d(&result, curve,
    0.5f, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);

// result.position = {0.5, 0.75}
// result.d1       = first derivative at t=0.5

qaws_curve_destroy(curve);
```

## Three workflows

**A. Geometry-first** -- create a curve, evaluate it directly.

**B. Motion-on-top** -- wrap a curve with a traversal for constant-speed or accelerated motion.

**C. Tooling and inspection** -- query properties, convert between families, export for rendering.

## Building

### Sharpmake (primary)

```bash
cd buildsystem
bootstrap.bat
```

Generates Visual Studio solutions. This is the primary build system.

### CMake (backup)

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

## Repository structure

```
src/qaws/            Public headers and implementations
src/qaws/core/       Portable eval headers (usable from C, HLSL, GLSL, Halide)
src/qaws/simd/       SIMD abstraction (AVX2, SSE2, NEON, scalar fallback)
src/qaws/batch/      Batch evaluation using SIMD backends
src/qaws/internal/   Internal numerical machinery
tests/               Test suite with SVG visual output
docs/                This documentation
docs/api/            Per-header API reference
examples/            GPU compute shaders (HLSL, GLSL) and Halide generators
buildsystem/         Sharpmake project generation
```
