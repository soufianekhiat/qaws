# API Reference: qaws_platform.h, qaws_precision.h, qaws_names.h, qaws_both.h

The backend, the scalar type and the precision of every name.

- `qaws_platform.h`: backend detection, the scalar type on shader and
  Halide backends, qualifier, struct, math, select and helper macros.
- `qaws_precision.h`: on the C backend, the names of the current precision
  (f32 or f64).
- `qaws_names.h`: the generated table mapping every precision-bearing name
  to its suffixed form.
- `qaws_both.h`: both precisions declared in one translation unit.

See [precision.md](../precision.md) for the precision model and
[backends.md](../backends.md) for the backends.

---

## Backends

```c
#define QAWS_BACKEND_C      0
#define QAWS_BACKEND_HLSL   1
#define QAWS_BACKEND_GLSL   2
#define QAWS_BACKEND_HALIDE 3
```

`QAWS_BACKEND` may be defined before inclusion. Otherwise it is detected:
HLSL when `__HLSL_VERSION` is defined, GLSL when `GL_core_profile` or
`GL_es_profile` is, Halide when `HALIDE_HALIDERUNTIME_H` is, else C. The
bootstrap headers `qaws_hlsl.h`, `qaws_glsl.h` and `qaws_halide.h` set it
explicitly ([core.md](core.md)).

## Scalar type

| Backend | `qaws_scalar` | `QAWS_SCALAR_IS_FLOAT` | `QAWS_LITERAL(x)` | `QAWS_EPSILON` |
|---|---|---|---|---|
| C | `qaws_f32` or `qaws_f64` (from `qaws_precision.h`) | default 1, settable | `x##f` (f32) or `x` (f64) | `1e-6f` or `1e-12` |
| HLSL | `float` (typedef) | 1 | `(x)` | `1e-6f` |
| GLSL | `float` (macro) | 1 | `(x)` | `1e-6` |
| Halide | `qaws_halide_scalar` | default 1, settable | a Halide constant of 32 or 64 bits | `1e-6f` or `1e-12` |

On the C backend `qaws_platform.h` declares both precisions:

```c
typedef float  qaws_f32;
typedef double qaws_f64;
```

On Halide, `qaws_halide_scalar` derives from `Halide::Expr` and is built
from a `double` (a constant of `QAWS_HALIDE_FLOAT_BITS`, 32 or 64), an
`Expr` or a `FuncRef` (a buffer's element).

Other scalar macros:

| Macro | Meaning |
|---|---|
| `QAWS_ZERO`, `QAWS_ONE` | `QAWS_LITERAL(0.0)`, `QAWS_LITERAL(1.0)` |
| `QAWS_ITOF(n)`, `QAWS_FTOI(x)` | int to scalar and scalar to int (constructor casts in GLSL, `Halide::cast` in Halide) |
| `qaws_index` | an index the data decide (a knot span): `int`, or `Halide::Expr` evaluated per element |
| `QAWS_SELECT(cond, t, f)` | branchless select: `?:`, or `Halide::select` |
| `QAWS_UNUSED(x)` | marks a parameter unused (empty in GLSL) |

## Qualifiers and structs

| Macro | C | HLSL | GLSL | Halide |
|---|---|---|---|---|
| `QAWS_INLINE` | `static inline` | `inline` | (empty) | `inline` |
| `QAWS_CONSTEXPR` | `static const` | `static const` | `const` | `static const` |
| `QAWS_TYPE_DEF` | `typedef` | `typedef` | (empty) | `typedef` |
| `QAWS_OUT` | (empty) | `out` | `out` | (empty) |
| `QAWS_INOUT` | (empty) | `inout` | `inout` | (empty) |

```c
QAWS_STRUCT(name) { members } QAWS_STRUCT_END(name);
```

A struct type every backend names the same way: a typedef of an anonymous
struct in C, HLSL and Halide; a named struct in GLSL, which has no typedef.

## Math and helpers

`QAWS_SQRT`, `QAWS_FABS`, `QAWS_FLOOR`, `QAWS_CEIL`, `QAWS_SIN`, `QAWS_COS`,
`QAWS_ATAN2(y, x)`, `QAWS_POW(x, y)`, `QAWS_FMOD(x, y)` map to the
backend's functions: `sqrtf`/`sqrt` and so on in C by precision, the
intrinsics in HLSL and GLSL (`atan(y, x)` and `mod` in GLSL), `Halide::`
functions in Halide (`QAWS_FMOD` as `x - floor(x / y) * y`).

```c
qaws_scalar qaws_min(qaws_scalar a, qaws_scalar b);
qaws_scalar qaws_max(qaws_scalar a, qaws_scalar b);
qaws_scalar qaws_clamp(qaws_scalar x, qaws_scalar lo, qaws_scalar hi);
qaws_scalar qaws_lerp(qaws_scalar a, qaws_scalar b, qaws_scalar t);   /* (1 - t) a + t b */
```

In C these are `QAWS_INLINE` functions, defined once per precision by
`qaws_precision.h` (precision-bearing: `qaws_min_f32`, `qaws_min_f64`, ...).
In HLSL, GLSL and Halide they are macros over `min`, `max`, `clamp` and
`lerp` (`mix` in GLSL).

## The precision of names (C backend)

The library holds both precisions: every scalar-bearing function and type
exists as `name_f32` and `name_f64`. Programs write the unsuffixed name.

```c
#define QAWS_SFX(name) name##_f32   /* QAWS_SCALAR_IS_FLOAT 1 */
#define QAWS_SFX(name) name##_f64   /* QAWS_SCALAR_IS_FLOAT 0 */
```

- **`QAWS_SCALAR_IS_FLOAT`** (default 1) picks which of the two the
  unsuffixed names, `qaws_scalar` and the scalar macros mean. Define it
  before including any qaws header.
- **`qaws_precision.h`** has no include guard. Each inclusion redefines
  `QAWS_SFX`, `qaws_scalar`, `QAWS_LITERAL`, `QAWS_EPSILON` and the math
  macros for the current `QAWS_SCALAR_IS_FLOAT`, includes `qaws_names.h`,
  and declares the scalar helpers for that precision once. Included again
  after `QAWS_SCALAR_IS_FLOAT` changes, it switches every precision-bearing
  name to the other form. `qaws_platform.h` includes it at its end,
  outside its own guard.
- **`qaws_names.h`** is generated (do not edit). It is one line per
  precision-bearing name, about 900 of them, each of the form

  ```c
  #define qaws_curve_create_bspline QAWS_SFX(qaws_curve_create_bspline)
  ```

  Naming rule: every function, struct and opaque type declared by a
  precision-guarded header is precision-bearing and gets the suffix, even
  when it holds no scalar itself (`qaws_curve`, `qaws_array_view`), since
  it belongs to one precision's API. The functions of the portable core
  headers are suffixed too; their small struct types are not (see
  [core.md](core.md)). Names declared once for both
  precisions have a single form: the precision-free headers
  (`qaws_allocator.h`, `qaws_batch_executor.h`, `qaws_status.h`) and the
  plain enums guarded by `QAWS_ENUM_*` (`qaws_element_type`,
  `qaws_curve_view_family`, ...). The mapping goes through
  `QAWS_SFX`, so it follows whatever precision is current where the name
  is used.
- **Public headers** are guarded once per precision
  (`QAWS_<HEADER>_H_F32` / `_F64`), so including a header again after a
  precision switch declares its other form.

Linking: one library holds both forms, so f32 and f64 code can live in one
program, each translation unit using its own `QAWS_SCALAR_IS_FLOAT`.

## qaws_both.h

```c
#include "qaws_both.h"
```

Both precisions in one translation unit: includes `qaws.h` once with
`QAWS_SCALAR_IS_FLOAT` 1 and once with 0, so every public function and type
is declared as `name_f32` and `name_f64`, then restores the caller's
precision. The unsuffixed names, `qaws_scalar` and the scalar macros keep
meaning the caller's precision (`QAWS_SCALAR_IS_FLOAT`, default 1).

```c
#define QAWS_SCALAR_IS_FLOAT 0
#include "qaws_both.h"

qaws_curve* c = NULL;                           /* f64: qaws_curve_f64 */
qaws_curve_f32* preview = NULL;                 /* the f32 form, named explicitly */
qaws_curve_create_bspline_f32(&preview, &desc32);
```
