# Qaws

| Platform | f32 | f64 |
|---|---|---|
| **Linux** | [![Linux f32](https://github.com/soufianekhiat/qaws/actions/workflows/linux-f32.yml/badge.svg?branch=main)](https://github.com/soufianekhiat/qaws/actions/workflows/linux-f32.yml) | [![Linux f64](https://github.com/soufianekhiat/qaws/actions/workflows/linux-f64.yml/badge.svg?branch=main)](https://github.com/soufianekhiat/qaws/actions/workflows/linux-f64.yml) |
| **macOS** | [![macOS f32](https://github.com/soufianekhiat/qaws/actions/workflows/macos-f32.yml/badge.svg?branch=main)](https://github.com/soufianekhiat/qaws/actions/workflows/macos-f32.yml) | [![macOS f64](https://github.com/soufianekhiat/qaws/actions/workflows/macos-f64.yml/badge.svg?branch=main)](https://github.com/soufianekhiat/qaws/actions/workflows/macos-f64.yml) |
| **Windows** | [![Windows f32](https://github.com/soufianekhiat/qaws/actions/workflows/windows-f32.yml/badge.svg?branch=main)](https://github.com/soufianekhiat/qaws/actions/workflows/windows-f32.yml) | [![Windows f64](https://github.com/soufianekhiat/qaws/actions/workflows/windows-f64.yml/badge.svg?branch=main)](https://github.com/soufianekhiat/qaws/actions/workflows/windows-f64.yml) |

**Qaws** (قوس) — Arabic for *arc*, *bow*, *curve*. A dependency-free C11
library for creating, evaluating, sampling, traversing, and inspecting
parametric curves and surfaces in 2D and 3D.

One source, four backends: **C** &middot; **HLSL** &middot; **GLSL** &middot; **Halide**

---

## Quick start

```c
#include "qaws.h"

qaws_vec2 pts[] = { {0,0}, {0.3f,1}, {0.7f,1}, {1,0} };
qaws_bezier_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = pts,
    .control_point_count = 4
};

qaws_curve* curve = NULL;
qaws_curve_create_bezier(&curve, &desc);

qaws_eval_result_2d r;
qaws_curve_evaluate_2d(&r, curve, 0.5f, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
// r.position = point on curve at t=0.5
// r.d1       = first derivative (tangent direction)

qaws_curve_destroy(curve);
```

## Curve families

Every family is created through a descriptor and shares the same
evaluation, sampling, traversal, and inspection API after creation.

| Family | Degree | Spans | Continuity | Closed | Rational | Dim | Header |
|---|:---:|:---:|:---:|:---:|:---:|:---:|---|
| **Bezier** | 1..N | 1 | -- | -- | -- | 2D 3D | `qaws_bezier.h` |
| **Hermite** | 3 | n-1 | C1 | -- | -- | 2D 3D | `qaws_hermite.h` |
| **Catmull-Rom** | 3 | n-3 / n | C1 | yes | -- | 2D 3D | `qaws_catmull_rom.h` |
| **B-Spline** | 1..N | knot-dep. | C(p-1) | yes | -- | 2D 3D | `qaws_bspline.h` |
| **NURBS** | 1..N | knot-dep. | C(p-1) | yes | yes | 2D 3D | `qaws_nurbs.h` |
| **Trajectory** | 3 | k-1 | C1 | yes | -- | 2D 3D | `qaws_trajectory.h` |
| **Yuksel C2** | 2 | n-1 / n | C2 | yes | -- | 2D 3D | `qaws_yuksel.h` |
| **Rational Bezier** | 1..N | 1 | -- | -- | yes | 2D 3D | `qaws_rational_bezier.h` |
| **Arc** | analytic | seg. | C-inf | opt. | -- | 2D 3D | `qaws_arc.h` |
| **Polynomial** | 1..N | 1 | C-inf | -- | -- | 2D 3D | `qaws_polynomial.h` |
| **Clothoid** | transcend. | 1 | C-inf | -- | -- | 2D | `qaws_clothoid.h` |
| **Subdivision** | scheme-dep. | refined | C1..C3 | yes | -- | 2D 3D | `qaws_subdivision.h` |
| **Composite** | per-seg. | seg. | user-dep. | opt. | -- | 2D 3D | `qaws_composite.h` |

## Surfaces

| Surface | Header |
|---|---|
| Bezier patch | `qaws_surface_bezier.h` |
| B-Spline | `qaws_surface_bspline.h` |
| NURBS | `qaws_surface_nurbs.h` |
| Swept (profile along path) | `qaws_surface_swept.h` |
| Ruled (linear blend of two curves) | `qaws_surface_ruled.h` |

## What you can do

### Create

13 curve families + 5 surface types. See the [curve families table](#curve-families) above.

### Evaluate

| Capability | Function |
|---|---|
| Position, D1, D2, D3 | `qaws_curve_evaluate_2d/3d` |
| Per-span evaluation | `qaws_curve_evaluate_span_2d/3d` |
| Batch (SIMD-accelerated) | `qaws_curve_evaluate_batch_2d/3d` |
| Surface point + partials + normal | `qaws_surface_evaluate` |

### Sample

| Mode | Function |
|---|---|
| Uniform / adaptive | `qaws_curve_sample_2d/3d` |
| Curvature-weighted | `qaws_curve_sample_curvature_weighted_2d/3d` |
| Feature-preserving | `qaws_curve_sample_feature_preserving_2d/3d` |
| Streaming (callback) | `qaws_curve_sample_streaming_2d/3d` |

### Traverse

| Feature | |
|---|---|
| Modes | parameter, arc-length, timed |
| Motion profiles | constant speed, constant accel., trapezoidal, S-curve, custom |
| Easing | 10 functions (quad, cubic, sine -- in/out/in-out) |
| Wrap modes | clamp, loop, ping-pong |
| Multi-curve | `qaws_traversal_create_multi` |

### Inspect

| Category | Functions |
|---|---|
| Geometry | arc length, bounding box, closest point |
| Differential | tangent, normal, curvature, torsion, speed, Frenet frame |
| Topology | inflection points, extrema, winding number |
| Intersection | curve-curve, self-intersection |
| Visualization | curvature comb |

### Convert and operate

| Conversion | Operations | Export / Import |
|---|---|---|
| Hermite -> Bezier | Split at parameter | SVG path data |
| Catmull-Rom -> Bezier | Join two curves | Polyline export |
| Bezier -> B-Spline | Offset (2D) | B-spline fitting |
| B-Spline -> NURBS | Arc-length reparam. | Polyline import |
| Degree elevation | Reverse | |
| Degree reduction | | |

## Multi-backend

The core evaluation math lives in header-only files under `src/qaws/core/`.
These compile on all four backends from the same source:

| Backend | What you get | Scalar type |
|---|---|---|
| **C** | Full runtime + prepare + SIMD batch | `float` or `double` |
| **HLSL** | 22 core eval functions | `float` |
| **GLSL** | 22 core eval functions | `float` |
| **Halide** | 22 core eval functions (symbolic IR) | `Halide::Expr` |

```hlsl
// HLSL compute shader
#include "qaws_hlsl.h"
#include "core/qaws_decasteljau_core.h"

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID) {
    float t = (float)id.x / (float)(sample_count - 1);
    qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);
    results[id.x * 2]     = pos.x;
    results[id.x * 2 + 1] = pos.y;
}
```

Backend detection is automatic (`__HLSL_VERSION`, `GL_core_profile`,
`HALIDE_HALIDERUNTIME_H`) or forced via bootstrap headers (`qaws_hlsl.h`,
`qaws_glsl.h`, `qaws_halide.h`).

## Building

Sharpmake (primary) or CMake (backup).

```bash
# CMake
cmake -B build -DQAWS_SCALAR_IS_FLOAT=ON
cmake --build build --config Release
ctest --build-config Release --output-on-failure

# Sharpmake
sharpmake /sources(@'buildsystem/sharpmake/src/main.sharpmake.cs')
```

### Options

| Option | Default | Effect |
|---|:---:|---|
| `QAWS_SCALAR_IS_FLOAT` | `ON` | `float` (32-bit) or `double` (64-bit) |
| `QAWS_ENABLE_SIMD` | `OFF` | AVX2 / SSE2 / NEON batch eval |
| `QAWS_ENABLE_WARNINGS` | `ON` | Compiler warnings |
| `QAWS_WARNINGS_AS_ERRORS` | `ON`(f64) / `OFF`(f32) | `-Werror` / `/WX` |

## Design

- **Zero dependencies** -- C11 standard library only (`<math.h>`, `<stdlib.h>`, `<string.h>`)
- **Immutable after creation** -- `qaws_curve` and `qaws_surface` are thread-safe for concurrent reads
- **Opaque handles** -- all internal state is hidden behind `qaws_curve*` / `qaws_surface*`
- **Status codes everywhere** -- every fallible function returns `qaws_status`
- **Custom allocators** -- `_ex` variants accept `qaws_allocator` for arena/pool/tracking
- **Stack-allocated curves** -- `qaws_curve_inline` for hot paths (no malloc, degree <= 7)

## Repository layout

```
src/qaws/            Public headers and C implementation
src/qaws/core/       Header-only eval kernels (all backends)
src/qaws/internal/   Private implementation details
src/qaws/simd/       SIMD abstraction (AVX2, SSE2, NEON, scalar)
src/qaws/batch/      SIMD batch evaluation
tests/               Test runner
examples/            HLSL, GLSL, Halide shader examples
docs/                Documentation
buildsystem/         Sharpmake project files
```

## Documentation

Full docs are in [`docs/`](docs/index.md):

- [Overview](docs/overview.md) -- design principles, quick start
- [Backends](docs/backends.md) -- C / HLSL / GLSL / Halide architecture
- [Curve Families](docs/curve_families.md) -- all 13 families in detail
- [Evaluation](docs/evaluation.md) -- eval flags, batch eval, core headers, prepare+eval
- [Sampling](docs/sampling.md) -- uniform, adaptive, curvature-weighted, streaming
- [Traversal](docs/traversal.md) -- motion profiles, easing, wrap modes
- [Inspection](docs/inspection.md) -- geometry queries, intersection, analysis
- [Building](docs/building.md) -- CMake and Sharpmake configuration
- [Capabilities](docs/capabilities.md) -- what runs on which backend, and which test or check backs it

## License

See [LICENSE](LICENSE).
