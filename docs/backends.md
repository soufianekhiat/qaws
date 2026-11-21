# Multi-backend Architecture

The evaluation cores in `src/qaws/core/` compile on four backends from one
set of headers. The backend is selected at compile time, when the headers are
included: there is no runtime dispatch and no build option for it. What runs
and is checked on each backend: [capabilities.md](capabilities.md).

## Backends at a glance

| Backend | Language | `QAWS_BACKEND` | Typical use |
|---|---|:---:|---|
| **C** | C11 | `QAWS_BACKEND_C` (0) | The full runtime library, tools, tests |
| **HLSL** | HLSL | `QAWS_BACKEND_HLSL` (1) | Direct3D compute and pixel shaders |
| **GLSL** | GLSL | `QAWS_BACKEND_GLSL` (2) | Vulkan and OpenGL compute shaders |
| **Halide** | C++ (Halide DSL) | `QAWS_BACKEND_HALIDE` (3) | Halide pipelines (AOT or JIT) |

## Detection

`qaws_platform.h` picks the backend from compiler-defined macros unless
`QAWS_BACKEND` is already defined:

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

## Bootstrap headers

Each non-C backend has a bootstrap header that sets `QAWS_BACKEND` (when not
already set) and includes `qaws_platform.h` and `qaws_core_types.h`:

| Header | Sets |
|---|---|
| `qaws_hlsl.h` | `QAWS_BACKEND` 1 |
| `qaws_glsl.h` | `QAWS_BACKEND` 2 |
| `qaws_halide.h` | `QAWS_BACKEND` 3 |

Then include the cores you need:

```hlsl
#include "qaws_hlsl.h"
#include "core/qaws_decasteljau_core.h"
```

`qaws_core_types.h` defines `qaws_vec2`, `qaws_vec3` (on C only when
`qaws_types.h` has not), the core result structs `qaws_eval_2d` /
`qaws_eval_3d` (position, d1, d2) and `QAWS_CORE_MAX_DEGREE` (16).

## What each backend gets

### C

- The whole runtime (`qaws.h`): opaque curves and surfaces, sampling,
  traversal, inspection, operations, regions, shatter, B-rep, exact kernel,
  differentiation, views and batches.
- Every core header, compiled as `static inline` functions.
- `qaws_prepare.h`: flat coefficient arrays from high-level inputs, for
  upload to the other backends ([evaluation.md](evaluation.md)).
- Both precisions in one library; `qaws_scalar` and the math macros follow
  `QAWS_SCALAR_IS_FLOAT` through `qaws_precision.h`
  ([precision.md](precision.md)).
- The SIMD kernels in `src/qaws/batch/` when built with SIMD
  ([building.md](building.md)).

### HLSL

Shader only: no runtime, no pointers, no allocation. `qaws_scalar` is
`float`. `QAWS_OUT` / `QAWS_INOUT` expand to `out` / `inout`, `QAWS_INLINE`
to `inline`, math to the HLSL intrinsics, `qaws_lerp` to `lerp`.

### GLSL

Same scope as HLSL. `qaws_scalar` is `float`. `QAWS_INLINE` and
`QAWS_TYPE_DEF` are empty (no `inline`, no `typedef`), structs are named
structs, `QAWS_ATAN2` maps to `atan(y, x)`, `QAWS_FMOD` to `mod`, `qaws_lerp`
to `mix`, int and scalar conversions use constructor casts.

### Halide

C++: `qaws_scalar` is `qaws_halide_scalar`, a `Halide::Expr` that also
converts from doubles and buffer elements. `QAWS_SCALAR_IS_FLOAT` (default 1)
picks 32- or 64-bit floats (`QAWS_HALIDE_FLOAT_BITS`). Core calls build
expression graphs at pipeline-definition time: loops unroll, `QAWS_SELECT`
is `Halide::select` (both branches evaluated), and `qaws_index` (a knot span,
a window start) is an `Expr` evaluated per element, so span search
(`qaws_find_span`) is a branchless scan returning a symbolic index.

## Platform macros

| Macro | C | HLSL | GLSL | Halide |
|---|---|---|---|---|
| `qaws_scalar` | `qaws_f32` or `qaws_f64` | `float` | `float` | `qaws_halide_scalar` |
| `QAWS_LITERAL(x)` | `x##f` (f32) or `x` (f64) | `(x)` | `(x)` | `make_const(Float(bits), x)` |
| `QAWS_ZERO` / `QAWS_ONE` | `QAWS_LITERAL(0.0)` / `(1.0)` | same | same | same |
| `QAWS_EPSILON` | `1e-6f` (f32), `1e-12` (f64) | `1e-6f` | `1e-6` | by precision |
| `QAWS_INLINE` | `static inline` | `inline` | empty | `inline` |
| `QAWS_CONSTEXPR` | `static const` | `static const` | `const` | `static const` |
| `QAWS_TYPE_DEF` | `typedef` | `typedef` | empty | `typedef` |
| `QAWS_STRUCT(n)` / `QAWS_STRUCT_END(n)` | `typedef struct` / `n` | same | `struct n` / empty | same as C |
| `QAWS_OUT`, `QAWS_INOUT` | empty | `out`, `inout` | `out`, `inout` | empty |
| `QAWS_SELECT(c, t, f)` | `(c) ? (t) : (f)` | same | same | `Halide::select` |
| `QAWS_ITOF(n)`, `QAWS_FTOI(x)` | C casts | C casts | `float(n)`, `int(x)` | `qaws_scalar(double(n))`, `Halide::cast<int>` |
| `qaws_index` | `int` | `int` | `int` | `Halide::Expr` |
| `QAWS_UNUSED(x)` | `(void)(x)` | same | empty | same |
| `QAWS_SQRT`, `QAWS_FABS`, `QAWS_FLOOR`, `QAWS_CEIL`, `QAWS_SIN`, `QAWS_COS`, `QAWS_ATAN2`, `QAWS_POW`, `QAWS_FMOD` | `<math.h>` (float or double forms) | intrinsics | intrinsics | `Halide::` functions |
| `qaws_min`, `qaws_max`, `qaws_clamp`, `qaws_lerp` | inline functions per precision | `min`, `max`, `clamp`, `lerp` | `min`, `max`, `clamp`, `mix` | `Halide::min`, `max`, `clamp`, `lerp` |

## Typical workflows

### GPU shader (HLSL / GLSL)

1. On the CPU, build flat arrays with a prepare function, or read records
   straight from buffers with `qaws_curve_record_core.h` /
   `qaws_surface_record_core.h`:
   ```c
   qaws_hermite_prepare_2d(a, b, c, d, points, tangents, count);
   ```
2. Upload the arrays to a GPU buffer.
3. In the shader, include the bootstrap header and the core:
   ```hlsl
   #include "qaws_hlsl.h"
   #include "core/qaws_cubic_poly_core.h"

   qaws_eval_2d r = qaws_cubic_eval_2d(a, b, c, d, local_t, 0x7);
   ```

### Halide pipeline

```cpp
#include "qaws_halide.h"
#include "core/qaws_decasteljau_core.h"

qaws_vec2 p0 = { cp_buf(0), cp_buf(1) };
/* ... */
qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);
```

The calls build Halide IR; schedule and compile as usual. High-degree
splines give large expression graphs, a one-time compilation cost.

### Direct C use (no runtime object)

```c
#include "qaws.h"
#include "qaws_core_types.h"
#include "core/qaws_decasteljau_core.h"

qaws_eval_2d r = qaws_decasteljau_eval_2d(flat_cp, degree, t, 0x7);
```

## Examples

| File | Backend | What it does |
|---|---|---|
| `eval_bezier.hlsl` / `eval_bezier.comp` | HLSL / GLSL | Cubic Bezier compute shader |
| `eval_cubic_poly.hlsl` / `eval_cubic_poly.comp` | HLSL / GLSL | Hermite / Catmull-Rom cubic polynomial |
| `eval_bspline.hlsl` / `eval_bspline.comp` | HLSL / GLSL | B-spline with knot span search |
| `eval_bezier_halide.cpp` | Halide | Cubic Bezier generator |
| `eval_record_halide.cpp` (+ `eval_record_halide_ref.c`) | Halide | B-spline and NURBS curve and surface records run by Halide's JIT against the C cores |
| `eval_parity.hlsl` + `eval_gpu_parity.c` | HLSL + C | GPU against C-core parity tool (Direct3D 11, Windows) |
| `diff_bspline_adjoint.hlsl` | HLSL | B-spline tangents, gather adjoint (one thread per control point), parameter adjoints |
| `diff_cdf_sampling.hlsl` | HLSL | Inverse-CDF sampling of a B-spline (arc length or curvature measure) with tangents and adjoints |
| `diff_surface_sampling.hlsl` | HLSL | Inverse-CDF warp of the unit square onto a B-spline surface by area |

`tests/shaders/check_shaders.ps1` compiles the HLSL (dxc `-T cs_6_0`) and
GLSL (glslangValidator `-V`) sources, including `tests/shaders/cores_check.*`,
by hand; it is not part of the build.

## See also

- [Capabilities](capabilities.md): what runs on which backend and what checks it
- [Evaluation](evaluation.md): core headers and prepare + eval
- [Building](building.md): SIMD options and install rules for the core headers
