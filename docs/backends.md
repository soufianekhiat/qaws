# Multi-backend architecture

Qaws core evaluation code compiles on four backends from a single set of
headers. The backend is selected at compile time — there is no runtime
dispatch and no build-time option to choose a backend.

## Backends at a glance

| Backend | Language | `QAWS_BACKEND` | Typical use |
|---|---|:---:|---|
| **C** | C11 | 0 | CPU runtime, tools, tests, SIMD batch eval |
| **HLSL** | HLSL (SM 5.0+) | 1 | DirectX compute / pixel shaders |
| **GLSL** | GLSL 450 / ES 3.1 | 2 | Vulkan / OpenGL compute shaders |
| **Halide** | C++ (Halide DSL) | 3 | AOT image-processing pipelines |

## How detection works

`qaws_platform.h` auto-detects the backend from compiler-defined macros:

```c
#if defined(__HLSL_VERSION)
#  define QAWS_BACKEND QAWS_BACKEND_HLSL
#elif defined(GL_core_profile) || defined(GL_es_profile)
#  define QAWS_BACKEND QAWS_BACKEND_GLSL
#elif defined(HALIDE_HALIDERUNTIME_H)
#  define QAWS_BACKEND QAWS_BACKEND_HALIDE
#else
#  define QAWS_BACKEND QAWS_BACKEND_C
#endif
```

You can override detection by defining `QAWS_BACKEND` before including any
Qaws header, or by including one of the bootstrap headers which set it for
you.

## Bootstrap headers

Each non-C backend has a one-line bootstrap header that forces the backend
and pulls in the platform layer:

| Header | Sets backend to | Include path |
|---|---|---|
| `qaws_hlsl.h` | `QAWS_BACKEND_HLSL` (1) | `#include "qaws_hlsl.h"` |
| `qaws_glsl.h` | `QAWS_BACKEND_GLSL` (2) | `#include "qaws_glsl.h"` |
| `qaws_halide.h` | `QAWS_BACKEND_HALIDE` (3) | `#include "qaws_halide.h"` |

After the bootstrap header, include the core eval headers you need:

```hlsl
#include "qaws_hlsl.h"
#include "core/qaws_decasteljau_core.h"
```

## What each backend gets

### C backend

The C backend is the full-featured runtime. It provides:

- **All core eval headers** — inline functions in `src/qaws/core/`
- **Prepare functions** — `qaws_prepare.h` pre-computes flat coefficient
  arrays from high-level descriptors
- **Full runtime API** — `qaws.h` with opaque `qaws_curve*` / `qaws_surface*`
  objects, sampling, traversal, inspection, operations, conversion, export,
  import
- **SIMD batch evaluation** — AVX2, SSE2, NEON accelerated paths in
  `src/qaws/batch/`
- **Scalar precision** — `float` or `double` via `QAWS_SCALAR_IS_FLOAT`

### HLSL backend

Shader-only. No runtime, no pointers, no malloc.

- **All 22 core eval functions** from `src/qaws/core/`
- **`QAWS_OUT` / `QAWS_INOUT`** expand to `out` / `inout` for output array
  parameters
- **`QAWS_INLINE`** expands to `inline`
- **`QAWS_TYPE_DEF`** is empty (HLSL struct tags are type names)
- **Math maps to HLSL builtins**: `sqrt`, `abs`, `floor`, `sin`, `cos`,
  `atan2`, `pow`, `fmod`
- **`qaws_lerp`** maps to `lerp`
- Scalar is always `float`

### GLSL backend

Same scope as HLSL with minor mapping differences:

- **All 22 core eval functions**
- **`QAWS_OUT` / `QAWS_INOUT`** expand to `out` / `inout`
- **`QAWS_INLINE`** is empty (GLSL has no `inline` keyword)
- **`QAWS_TYPE_DEF`** is empty
- **`qaws_lerp`** maps to `mix`
- **`QAWS_ATAN2`** maps to `atan` (GLSL uses `atan(y,x)`)
- **`QAWS_FMOD`** maps to `mod`
- Scalar is always `float`

### Halide backend

C++ meta-programming backend. `qaws_scalar` is `Halide::Expr`.

- **All 22 core eval functions** — loop bodies unroll at pipeline definition
  time and build symbolic expression graphs
- **`qaws_find_span`** has a dedicated implementation that builds a chain of
  `Halide::select()` calls returning a symbolic integer `Halide::Expr`
- **`QAWS_SELECT`** maps to `Halide::select` (evaluates both branches — safe
  for Qaws because discarded values are never used)
- **`QAWS_INLINE`** expands to `inline`
- **Math maps to Halide namespace**: `Halide::sqrt`, `Halide::sin`, etc.
- **`QAWS_LITERAL`** builds `Halide::Internal::make_const(Halide::Float(bits), x)`
- Scalar precision — `float` or `double` via `QAWS_SCALAR_IS_FLOAT`
  (controls `QAWS_HALIDE_FLOAT_BITS`)

## Platform abstraction macros

All core headers are written using these macros from `qaws_platform.h` so
that a single source compiles on every backend:

| Macro | C (float) | HLSL | GLSL | Halide |
|---|---|---|---|---|
| `qaws_scalar` | `float` | `float` | `float` | `Halide::Expr` |
| `QAWS_LITERAL(x)` | `x##f` | `(x)` | `(x)` | `make_const(...)` |
| `QAWS_ZERO` / `QAWS_ONE` | `0.0f` / `1.0f` | `(0.0)` / `(1.0)` | `(0.0)` / `(1.0)` | Halide const expr |
| `QAWS_INLINE` | `static inline` | `inline` | *(empty)* | `inline` |
| `QAWS_TYPE_DEF` | `typedef` | *(empty)* | *(empty)* | `typedef` |
| `QAWS_OUT` | *(empty)* | `out` | `out` | *(empty)* |
| `QAWS_INOUT` | *(empty)* | `inout` | `inout` | *(empty)* |
| `QAWS_SELECT(c,t,f)` | `(c)?(t):(f)` | `(c)?(t):(f)` | `(c)?(t):(f)` | `Halide::select` |
| `QAWS_SQRT(x)` | `sqrtf(x)` | `sqrt(x)` | `sqrt(x)` | `Halide::sqrt(x)` |
| `qaws_min/max/clamp/lerp` | C functions | `min/max/clamp/lerp` | `min/max/clamp/mix` | `Halide::min/max/clamp/lerp` |

## Typical workflows

### GPU shader (HLSL / GLSL)

1. On the CPU, call prepare functions to build flat arrays:
   ```c
   qaws_hermite_prepare_2d(points, tangents, count, a, b, c, d);
   ```
2. Upload the flat arrays to a GPU buffer.
3. In the shader, include the bootstrap header and the core eval header:
   ```hlsl
   #include "qaws_hlsl.h"
   #include "core/qaws_cubic_poly_core.h"

   // Load a, b, c, d from buffer for the target span
   qaws_eval_2d result = qaws_cubic_eval_2d(a, b, c, d, local_t, 0x7);
   ```

### Halide AOT pipeline

1. Include the bootstrap header and core headers in your Halide generator:
   ```cpp
   #include "qaws_halide.h"
   #include "core/qaws_decasteljau_core.h"
   ```
2. Build control points from `Halide::ImageParam` or `Halide::Param`:
   ```cpp
   qaws_vec2 p0 = { cp_buf(0), cp_buf(1) };
   // ...
   qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);
   ```
3. The function calls build Halide IR. Schedule and compile as usual.
   The expression graph for high-degree splines can be large but this is
   a one-time AOT compilation cost.

### Direct C usage (no runtime)

For hot paths where the `qaws_curve*` abstraction overhead matters, include
the core headers directly in C code:

```c
#include "qaws_platform.h"
#include "core/qaws_decasteljau_core.h"

qaws_eval_2d r = qaws_decasteljau_eval_2d(flat_cp, degree, t, 0x7);
```

## Example shaders

Complete working examples are in the `examples/` directory:

| File | Backend | What it does |
|---|---|---|
| `eval_bezier.hlsl` | HLSL | Cubic Bezier compute shader |
| `eval_bezier.comp` | GLSL | Cubic Bezier compute shader (Vulkan) |
| `eval_bezier_halide.cpp` | Halide | Cubic Bezier AOT generator |
| `eval_cubic_poly.hlsl` | HLSL | Hermite/CatRom cubic polynomial |
| `eval_cubic_poly.comp` | GLSL | Hermite/CatRom cubic polynomial |
| `eval_bspline.hlsl` | HLSL | B-spline with knot span search |
| `eval_bspline.comp` | GLSL | B-spline with knot span search |
| `diff_bspline_adjoint.hlsl` | HLSL | B-spline tangents, gather adjoint (one thread per control point) and parameter adjoints |
| `diff_cdf_sampling.hlsl` | HLSL | Inverse-CDF sampling of a B-spline (arc length or curvature measure): span measures, solved samples with first and second tangents, adjoint multipliers and a gather adjoint per control point |

## See also

- [Feature availability per backend](feature_per_backends.md) — complete
  function-by-function matrix
- [Evaluation](evaluation.md) — core eval headers and prepare+eval
  architecture in detail
- [Building](building.md) — SIMD options and install rules for core headers
