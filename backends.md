# Qaws Multi-Backend Architecture Proposal

Inspired by the [Alwan](../alwan/) color library's compile-time backend selection across C, HLSL, GLSL, and Halide.

---

## Goal

Make qaws curve evaluation usable as a **header-include library** in HLSL, GLSL, and Halide — the same way Alwan's `*_core.h` headers are included directly in shaders. No vtables, no pointers, no heap allocation. Just `#include` and call.

The architecture is **prepare + eval**: prepare functions (C only) output flat arrays from raw geometry input; eval functions (all backends, header-only) do the math. The existing C runtime API (`qaws_curve*`, vtable dispatch, allocators) stays as an optional convenience wrapper on top.

---

## How Alwan Does It

Alwan's architecture has three layers:

| Layer | Files | Backends |
|-------|-------|----------|
| **Platform** | `alwan_platform.h` | All — defines `alwan_scalar`, math macros, qualifiers |
| **Core** | `core/*_core.h` (header-only) | All — pure math, value-returning, no pointers |
| **API** | `api/*.c` | C only — pointer wrappers around core functions |

Key mechanisms in `alwan_platform.h`:

- **Backend detection** via compiler predefined macros (`__HLSL_VERSION`, `GL_core_profile`, `HALIDE_HALIDERUNTIME_H`)
- **`ALWAN_TYPE_DEF`** — expands to `typedef` in C/Halide, nothing in HLSL/GLSL (where the struct tag itself is the type name)
- **`ALWAN_INLINE`** — `static inline` in C, `inline` in HLSL/Halide, empty in GLSL
- **`ALWAN_CONSTEXPR`** — `static const` in C/HLSL/Halide, `const` in GLSL
- **`ALWAN_LITERAL(x)`** — `x` or `x##f` in C, `(x)` in HLSL/GLSL, `Halide::Internal::make_const(...)` in Halide
- **`ALWAN_SELECT(cond, t, f)`** — ternary in C/HLSL/GLSL, `Halide::select()` in Halide
- **`ALWAN_REF(p)` / `ALWAN_ADDR(x)`** — dereference/address in C/Halide, pass-by-value in HLSL/GLSL
- **Math macros** (`ALWAN_SQRT`, `ALWAN_CBRT`, etc.) — dispatch to `sqrtf`/`sqrt`/`sqrt`/`Halide::sqrt`
- **Utility functions** (`alwan_min`, `alwan_lerp`, etc.) — inline functions in C, macros to intrinsics in HLSL/GLSL/Halide

**Halide scalar type**: A custom struct `alwan_halide_scalar` inherits from `Halide::Expr` and adds an implicit constructor from `double`, so bare numeric literals work naturally in struct initializers.

Bootstrap headers (`alwan_hlsl.h`, `alwan_glsl.h`, `alwan_halide.h`) force the backend define and include platform + types. Users then include whichever `*_core.h` modules they need.

---

## Proposed Architecture for Qaws

### Directory Layout

```
src/qaws/
├── (existing files unchanged — C runtime API)
│
├── qaws_platform.h              # Backend detection, scalar type, math macros
├── qaws_core_types.h            # Cross-backend type definitions (vec2, vec3, eval results)
├── qaws_prepare.h               # Prepare functions — C only, outputs flat arrays
├── qaws_prepare.c               # Prepare function implementations
├── qaws_hlsl.h                  # HLSL bootstrap: force backend + include platform + types
├── qaws_glsl.h                  # GLSL bootstrap
├── qaws_halide.h                # Halide bootstrap
│
├── core/                        # Header-only, backend-agnostic curve math
│   ├── qaws_decasteljau_core.h  # De Casteljau evaluation (any degree, 2D/3D)
│   ├── qaws_cox_deboor_core.h   # Cox-de Boor B-spline basis
│   ├── qaws_nurbs_core.h        # Weighted rational evaluation
│   ├── qaws_hermite_core.h      # Cubic Hermite polynomial
│   ├── qaws_horner_core.h       # Horner polynomial evaluation
│   ├── qaws_arc_core.h          # Circular/elliptical arc math
│   ├── qaws_fresnel_core.h      # Fresnel integrals (clothoid)
│   └── qaws_surface_core.h     # Tensor-product surface evaluation
│
├── simd/                        # SIMD acceleration (C backend only)
│   ├── qaws_simd.h              # ISA detection + dispatch
│   ├── qaws_simd_sse2.h         # SSE2 intrinsics (x86)
│   ├── qaws_simd_avx2.h         # AVX2 intrinsics (x86)
│   ├── qaws_simd_neon.h         # NEON intrinsics (ARM)
│   ├── qaws_simd_scalar.h       # Scalar fallback
│   └── qaws_simd_types.h        # Lane-width-agnostic type aliases
│
└── batch/                       # SIMD batch evaluation (C backend only)
    ├── qaws_batch_bezier.c      # SIMD De Casteljau for N parameters
    ├── qaws_batch_bspline.c     # SIMD Cox-de Boor
    ├── qaws_batch_nurbs.c       # SIMD weighted sums
    └── qaws_batch_surface.c     # SIMD tensor-product grid eval
```

---

## `qaws_platform.h` — The Backend Abstraction

Single header, following Alwan's pattern exactly.

### 1.1 Backend Detection

```c
#define QAWS_BACKEND_C      0
#define QAWS_BACKEND_HLSL   1
#define QAWS_BACKEND_GLSL   2
#define QAWS_BACKEND_HALIDE 3

#ifndef QAWS_BACKEND
# if defined(__HLSL_VERSION)
#   define QAWS_BACKEND QAWS_BACKEND_HLSL
# elif defined(GL_core_profile) || defined(GL_es_profile)
#   define QAWS_BACKEND QAWS_BACKEND_GLSL
# elif defined(HALIDE_HALIDERUNTIME_H)
#   define QAWS_BACKEND QAWS_BACKEND_HALIDE
# else
#   define QAWS_BACKEND QAWS_BACKEND_C
# endif
#endif
```

### 1.2 Scalar Type

```c
#if QAWS_BACKEND == QAWS_BACKEND_C
# include <math.h>
# if QAWS_SCALAR_IS_FLOAT
    typedef float qaws_scalar;
#   define QAWS_LITERAL(x) x##f
#   define QAWS_EPSILON 1e-6f
# else
    typedef double qaws_scalar;
#   define QAWS_LITERAL(x) x
#   define QAWS_EPSILON 1e-12
# endif

#elif QAWS_BACKEND == QAWS_BACKEND_HLSL
# define QAWS_SCALAR_IS_FLOAT 1
# define QAWS_LITERAL(x) (x)
# define QAWS_EPSILON 1e-6f
  /* qaws_scalar is just float in HLSL — no typedef needed */

#elif QAWS_BACKEND == QAWS_BACKEND_GLSL
# define QAWS_SCALAR_IS_FLOAT 1
# define QAWS_LITERAL(x) (x)
# define QAWS_EPSILON 1e-6

#elif QAWS_BACKEND == QAWS_BACKEND_HALIDE
# include <Halide.h>
# ifndef QAWS_SCALAR_IS_FLOAT
#   define QAWS_SCALAR_IS_FLOAT 1
# endif
# if QAWS_SCALAR_IS_FLOAT
#   define QAWS_HALIDE_FLOAT_BITS 32
#   define QAWS_EPSILON 1e-6f
#   define QAWS_LITERAL(x) Halide::Internal::make_const(Halide::Float(32), (x))
# else
#   define QAWS_HALIDE_FLOAT_BITS 64
#   define QAWS_EPSILON 1e-12
#   define QAWS_LITERAL(x) Halide::Internal::make_const(Halide::Float(64), (x))
# endif
  /* Wrapper struct: inherits Halide::Expr, adds implicit double constructor
   * so bare numeric literals work in array initializers without QAWS_LITERAL() */
  struct qaws_halide_scalar : Halide::Expr {
      using Halide::Expr::Expr;
      qaws_halide_scalar() = default;
      qaws_halide_scalar(double v)
          : Halide::Expr(Halide::Internal::make_const(
                Halide::Float(QAWS_HALIDE_FLOAT_BITS), v)) {}
      qaws_halide_scalar(const Halide::Expr& e) : Halide::Expr(e) {}
  };
  typedef qaws_halide_scalar qaws_scalar;
#endif

#define QAWS_ZERO QAWS_LITERAL(0.0)
#define QAWS_ONE  QAWS_LITERAL(1.0)
```

### 1.3 Qualifier Macros

```c
#if QAWS_BACKEND == QAWS_BACKEND_C
# define QAWS_INLINE     static inline
# define QAWS_CONSTEXPR  static const
# define QAWS_TYPE_DEF   typedef

#elif QAWS_BACKEND == QAWS_BACKEND_HLSL
# define QAWS_INLINE     inline
# define QAWS_CONSTEXPR  static const
# define QAWS_TYPE_DEF

#elif QAWS_BACKEND == QAWS_BACKEND_GLSL
# define QAWS_INLINE
# define QAWS_CONSTEXPR  const
# define QAWS_TYPE_DEF

#elif QAWS_BACKEND == QAWS_BACKEND_HALIDE
# define QAWS_INLINE     inline
# define QAWS_CONSTEXPR  static const
# define QAWS_TYPE_DEF   typedef
#endif
```

`QAWS_TYPE_DEF` is the key trick from Alwan: in C and Halide, `QAWS_TYPE_DEF struct { ... } name;` creates a typedef. In HLSL/GLSL, it expands to nothing so the struct tag IS the type name.

### 1.4 Math Macros

```c
#if QAWS_BACKEND == QAWS_BACKEND_C
  /* C: dispatch float vs double */
# if QAWS_SCALAR_IS_FLOAT
#   define QAWS_SQRT(x)     sqrtf(x)
#   define QAWS_FABS(x)     fabsf(x)
#   define QAWS_FLOOR(x)    floorf(x)
#   define QAWS_CEIL(x)     ceilf(x)
#   define QAWS_SIN(x)      sinf(x)
#   define QAWS_COS(x)      cosf(x)
#   define QAWS_ATAN2(y,x)  atan2f(y,x)
#   define QAWS_POW(x,y)    powf(x,y)
#   define QAWS_FMOD(x,y)   fmodf(x,y)
# else
#   define QAWS_SQRT(x)     sqrt(x)
#   define QAWS_FABS(x)     fabs(x)
#   define QAWS_FLOOR(x)    floor(x)
#   define QAWS_CEIL(x)     ceil(x)
#   define QAWS_SIN(x)      sin(x)
#   define QAWS_COS(x)      cos(x)
#   define QAWS_ATAN2(y,x)  atan2(y,x)
#   define QAWS_POW(x,y)    pow(x,y)
#   define QAWS_FMOD(x,y)   fmod(x,y)
# endif

#elif QAWS_BACKEND == QAWS_BACKEND_HLSL
# define QAWS_SQRT(x)     sqrt(x)
# define QAWS_FABS(x)     abs(x)
# define QAWS_FLOOR(x)    floor(x)
# define QAWS_CEIL(x)     ceil(x)
# define QAWS_SIN(x)      sin(x)
# define QAWS_COS(x)      cos(x)
# define QAWS_ATAN2(y,x)  atan2(y,x)
# define QAWS_POW(x,y)    pow(x,y)
# define QAWS_FMOD(x,y)   fmod(x,y)

#elif QAWS_BACKEND == QAWS_BACKEND_GLSL
# define QAWS_SQRT(x)     sqrt(x)
# define QAWS_FABS(x)     abs(x)
# define QAWS_FLOOR(x)    floor(x)
# define QAWS_CEIL(x)     ceil(x)
# define QAWS_SIN(x)      sin(x)
# define QAWS_COS(x)      cos(x)
# define QAWS_ATAN2(y,x)  atan(y,x)    /* GLSL: atan(y,x) not atan2 */
# define QAWS_POW(x,y)    pow(x,y)
# define QAWS_FMOD(x,y)   mod(x,y)     /* GLSL: mod() not fmod() */

#elif QAWS_BACKEND == QAWS_BACKEND_HALIDE
# define QAWS_SQRT(x)     Halide::sqrt(x)
# define QAWS_FABS(x)     Halide::abs(x)
# define QAWS_FLOOR(x)    Halide::floor(x)
# define QAWS_CEIL(x)     Halide::ceil(x)
# define QAWS_SIN(x)      Halide::sin(x)
# define QAWS_COS(x)      Halide::cos(x)
# define QAWS_ATAN2(y,x)  Halide::atan2(y,x)
# define QAWS_POW(x,y)    Halide::pow(x,y)
# define QAWS_FMOD(x,y)   ((x) - Halide::floor((x) / (y)) * (y))
#endif
```

### 1.5 Branchless Select

```c
#if QAWS_BACKEND == QAWS_BACKEND_HALIDE
# define QAWS_SELECT(cond, t, f) Halide::select((cond), (t), (f))
#else
# define QAWS_SELECT(cond, t, f) ((cond) ? (t) : (f))
#endif
```

### 1.6 Utility Functions

```c
#if QAWS_BACKEND == QAWS_BACKEND_C
  QAWS_INLINE qaws_scalar qaws_min(qaws_scalar a, qaws_scalar b) { return (a < b) ? a : b; }
  QAWS_INLINE qaws_scalar qaws_max(qaws_scalar a, qaws_scalar b) { return (a > b) ? a : b; }
  QAWS_INLINE qaws_scalar qaws_clamp(qaws_scalar x, qaws_scalar lo, qaws_scalar hi) {
      return (x < lo) ? lo : (x > hi) ? hi : x;
  }
  QAWS_INLINE qaws_scalar qaws_lerp(qaws_scalar a, qaws_scalar b, qaws_scalar t) {
      return (QAWS_ONE - t) * a + t * b;
  }

#elif QAWS_BACKEND == QAWS_BACKEND_HLSL
# define qaws_min(a, b)        min(a, b)
# define qaws_max(a, b)        max(a, b)
# define qaws_clamp(x, lo, hi) clamp(x, lo, hi)
# define qaws_lerp(a, b, t)    lerp(a, b, t)

#elif QAWS_BACKEND == QAWS_BACKEND_GLSL
# define qaws_min(a, b)        min(a, b)
# define qaws_max(a, b)        max(a, b)
# define qaws_clamp(x, lo, hi) clamp(x, lo, hi)
# define qaws_lerp(a, b, t)    mix(a, b, t)

#elif QAWS_BACKEND == QAWS_BACKEND_HALIDE
# define qaws_min(a, b)        Halide::min(a, b)
# define qaws_max(a, b)        Halide::max(a, b)
# define qaws_clamp(x, lo, hi) Halide::clamp(x, lo, hi)
# define qaws_lerp(a, b, t)    Halide::lerp(a, b, t)
#endif
```

---

## `qaws_core_types.h` — Cross-Backend Types

Using `QAWS_TYPE_DEF` to define the same `qaws_vec2`, `qaws_vec3` names on all backends. No `qaws_core_` prefix — the types are the same everywhere.

```c
#ifndef QAWS_CORE_TYPES_H
#define QAWS_CORE_TYPES_H

#include "qaws_platform.h"

/* 2D point/vector */
QAWS_TYPE_DEF struct {
    qaws_scalar x, y;
} qaws_vec2;

/* 3D point/vector */
QAWS_TYPE_DEF struct {
    qaws_scalar x, y, z;
} qaws_vec3;

/* Evaluation result: position + derivatives (2D) */
QAWS_TYPE_DEF struct {
    qaws_vec2 position;
    qaws_vec2 d1;
    qaws_vec2 d2;
} qaws_eval_2d;

/* Evaluation result: position + derivatives (3D) */
QAWS_TYPE_DEF struct {
    qaws_vec3 position;
    qaws_vec3 d1;
    qaws_vec3 d2;
} qaws_eval_3d;

/* Maximum supported degree in core headers (bounds stack arrays) */
#ifndef QAWS_CORE_MAX_DEGREE
# define QAWS_CORE_MAX_DEGREE 16
#endif

#define QAWS_CORE_MAX_POINTS (QAWS_CORE_MAX_DEGREE + 1)

#endif /* QAWS_CORE_TYPES_H */
```

In C and Halide this creates `typedef struct { ... } qaws_vec2;`.
In HLSL/GLSL this creates `struct qaws_vec2 { float x, y; };` — the struct tag is the type name.

The existing C runtime `qaws_types.h` already defines `qaws_vec2`/`qaws_vec3` with `typedef`. On the C backend, `qaws_core_types.h` is not included — the runtime uses `qaws_types.h` as before. The core types header is for shader/Halide consumers only, and its definitions are binary-compatible with the C runtime types.

---

## Bootstrap Headers

Following Alwan's pattern, each backend gets a one-file bootstrap header.

### `qaws_hlsl.h`

```hlsl
#ifndef QAWS_HLSL_H
#define QAWS_HLSL_H

#ifndef QAWS_BACKEND
#define QAWS_BACKEND 1
#endif

#include "qaws_platform.h"
#include "qaws_core_types.h"

/* Then include whichever core modules you need:
 *   #include "core/qaws_decasteljau_core.h"
 *   #include "core/qaws_hermite_core.h"
 */

#endif
```

### `qaws_glsl.h`

```glsl
#ifndef QAWS_GLSL_H
#define QAWS_GLSL_H

#ifndef QAWS_BACKEND
#define QAWS_BACKEND 2
#endif

#include "qaws_platform.h"
#include "qaws_core_types.h"

#endif
```

### `qaws_halide.h`

```cpp
#ifndef QAWS_HALIDE_H
#define QAWS_HALIDE_H

/* Force Halide backend */
#ifndef QAWS_BACKEND
#define QAWS_BACKEND 3
#endif

/* Platform: includes <Halide.h>, typedefs qaws_scalar = Halide::Expr wrapper,
 * defines QAWS_LITERAL, math macros, QAWS_SELECT,
 * utility functions (qaws_min, qaws_max, qaws_clamp, ...) */
#include "qaws_platform.h"

/* Types: defines qaws_vec2/vec3 and eval result types */
#include "qaws_core_types.h"

#endif /* QAWS_HALIDE_H */
```

---

## Core Headers — Examples

Header-only, backend-agnostic. Compile as C, HLSL, GLSL, and Halide.

### `qaws_decasteljau_core.h` — Bezier evaluation

```c
#ifndef QAWS_DECASTELJAU_CORE_H
#define QAWS_DECASTELJAU_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ----------------------------------------------------------------
 * Cubic Bezier: position only (2D)
 *
 * Matches the inner loop of qaws_bezier.c:bezier_eval_span_2d()
 * for degree==3. Four CPs passed by value — no arrays needed.
 * ---------------------------------------------------------------- */

QAWS_INLINE qaws_vec2 qaws_bezier3_eval_2d(
    qaws_vec2 p0, qaws_vec2 p1,
    qaws_vec2 p2, qaws_vec2 p3,
    qaws_scalar t)
{
    qaws_scalar u = QAWS_ONE - t;
    qaws_vec2 q0, q1, q2, r0, r1, result;

    q0.x = u * p0.x + t * p1.x;  q0.y = u * p0.y + t * p1.y;
    q1.x = u * p1.x + t * p2.x;  q1.y = u * p1.y + t * p2.y;
    q2.x = u * p2.x + t * p3.x;  q2.y = u * p2.y + t * p3.y;

    r0.x = u * q0.x + t * q1.x;  r0.y = u * q0.y + t * q1.y;
    r1.x = u * q1.x + t * q2.x;  r1.y = u * q1.y + t * q2.y;

    result.x = u * r0.x + t * r1.x;
    result.y = u * r0.y + t * r1.y;
    return result;
}

/* ----------------------------------------------------------------
 * Arbitrary-degree De Casteljau (2D, position only)
 *
 * Control points in flat array: [x0,y0, x1,y1, ..., xn,yn]
 * degree must be <= QAWS_CORE_MAX_DEGREE.
 *
 * Matches qaws_internal_decasteljau(cp, degree, 2, t, buf)
 * ---------------------------------------------------------------- */

QAWS_INLINE qaws_vec2 qaws_decasteljau_2d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    qaws_scalar t)
{
    qaws_scalar work[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar u = QAWS_ONE - t;
    qaws_vec2 result;
    int r, i;

    for (i = 0; i <= degree; i++) {
        work[i * 2]     = cp[i * 2];
        work[i * 2 + 1] = cp[i * 2 + 1];
    }

    for (r = 1; r <= degree; r++) {
        for (i = 0; i <= degree - r; i++) {
            work[i * 2]     = u * work[i * 2]     + t * work[(i + 1) * 2];
            work[i * 2 + 1] = u * work[i * 2 + 1] + t * work[(i + 1) * 2 + 1];
        }
    }

    result.x = work[0];
    result.y = work[1];
    return result;
}

/* ----------------------------------------------------------------
 * Bezier derivative control points (2D)
 *
 * Given degree-N CPs, writes degree-(N-1) differenced CPs.
 * out_dp must be QAWS_CORE_MAX_POINTS * 2 in size.
 *
 * Matches qaws_internal_bezier_derivative_points(cp, degree, 2, out)
 * ---------------------------------------------------------------- */

QAWS_INLINE void qaws_bezier_deriv_points_2d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    qaws_scalar out_dp[QAWS_CORE_MAX_POINTS * 2])
{
    qaws_scalar deg = (qaws_scalar)degree;  /* QAWS_LITERAL needed for Halide */
    int i;
    for (i = 0; i < degree; i++) {
        out_dp[i * 2]     = deg * (cp[(i + 1) * 2]     - cp[i * 2]);
        out_dp[i * 2 + 1] = deg * (cp[(i + 1) * 2 + 1] - cp[i * 2 + 1]);
    }
}

/* 3D variants follow the same pattern with .z component added */
/* ... qaws_bezier3_eval_3d, qaws_decasteljau_3d, etc. ... */

#endif /* QAWS_DECASTELJAU_CORE_H */
```

### `qaws_cubic_poly_core.h` — Shared by Hermite, Catmull-Rom, Trajectory

```c
#ifndef QAWS_CUBIC_POLY_CORE_H
#define QAWS_CUBIC_POLY_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ----------------------------------------------------------------
 * Cubic polynomial: P(t) = a*t³ + b*t² + c*t + d
 *
 * This is the exact eval code currently duplicated in:
 *   - qaws_hermite.c:hermite_eval_span_2d()
 *   - qaws_catmull_rom.c:catmull_rom_eval_span_2d()
 *   - qaws_trajectory.c:trajectory_eval_span_2d()
 *
 * All three pre-compute a,b,c,d per span per dimension at creation,
 * then this function is the shared eval hot path.
 * ---------------------------------------------------------------- */

QAWS_INLINE qaws_eval_2d qaws_cubic_eval_2d(
    qaws_vec2 a, qaws_vec2 b,
    qaws_vec2 c, qaws_vec2 d,
    qaws_scalar t, int eval_flags)
{
    qaws_scalar t2 = t * t;
    qaws_scalar t3 = t2 * t;
    qaws_eval_2d result;

    result.position.x = a.x * t3 + b.x * t2 + c.x * t + d.x;
    result.position.y = a.y * t3 + b.y * t2 + c.y * t + d.y;

    result.d1.x = QAWS_LITERAL(3.0) * a.x * t2
                + QAWS_LITERAL(2.0) * b.x * t + c.x;
    result.d1.y = QAWS_LITERAL(3.0) * a.y * t2
                + QAWS_LITERAL(2.0) * b.y * t + c.y;

    result.d2.x = QAWS_LITERAL(6.0) * a.x * t + QAWS_LITERAL(2.0) * b.x;
    result.d2.y = QAWS_LITERAL(6.0) * a.y * t + QAWS_LITERAL(2.0) * b.y;

    return result;
}

/* 3D variant adds .z component */
/* ... qaws_cubic_eval_3d ... */

#endif /* QAWS_CUBIC_POLY_CORE_H */
```

### `qaws_bspline_basis_core.h` — B-spline basis (replaces malloc)

This header compiles on C, HLSL, and GLSL. **Not Halide** — the Cox-de Boor recurrence has data-dependent indexing patterns that don't map to Halide's functional model. For Halide B-spline evaluation, either pre-compute basis values on CPU or use a specialized Halide generator with `RDom`.

```c
#ifndef QAWS_BSPLINE_BASIS_CORE_H
#define QAWS_BSPLINE_BASIS_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

#if QAWS_BACKEND != QAWS_BACKEND_HALIDE

/* Maximum derivative order we support */
#define QAWS_CORE_MAX_DERIV 3

/* ----------------------------------------------------------------
 * B-spline basis functions + derivatives
 *
 * Matches qaws_internal_bspline_basis_derivs() but uses a
 * fixed-size stack buffer instead of malloc.
 *
 * out_ders: (max_k+1) x (degree+1) row-major, max size
 *           (QAWS_CORE_MAX_DERIV+1) * QAWS_CORE_MAX_POINTS
 * ---------------------------------------------------------------- */

QAWS_INLINE void qaws_bspline_basis_derivs(
    qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    int max_k,
    qaws_scalar out_ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS])
{
    /* Cox-de Boor recurrence with derivative computation.
     * Stack-allocated — no malloc needed.
     * See The NURBS Book, Algorithm A2.3 */
    /* ... implementation ... */
}

#endif /* !QAWS_BACKEND_HALIDE */

#endif /* QAWS_BSPLINE_BASIS_CORE_H */
```

### Usage from Each Backend

**C (existing runtime, unchanged):**
```c
#include "qaws.h"                         /* full C API as before */
#include "core/qaws_decasteljau_core.h"   /* optional: direct core access */
#include "core/qaws_cubic_poly_core.h"    /* optional: Hermite/CatmullRom/Trajectory eval */
```

**HLSL (compute shader — Bezier evaluation):**
```hlsl
#include "qaws_hlsl.h"
#include "core/qaws_decasteljau_core.h"

StructuredBuffer<float> cp : register(t0);     /* 8 floats: 4 CPs × 2D */
RWStructuredBuffer<float> results : register(u0);

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID) {
    float t = id.x / 1024.0;
    qaws_vec2 p0 = { cp[0], cp[1] };
    qaws_vec2 p1 = { cp[2], cp[3] };
    qaws_vec2 p2 = { cp[4], cp[5] };
    qaws_vec2 p3 = { cp[6], cp[7] };
    qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);
    results[id.x * 2]     = pos.x;
    results[id.x * 2 + 1] = pos.y;
}
```

**HLSL (compute shader — Hermite/CatmullRom/Trajectory span evaluation):**
```hlsl
#include "qaws_hlsl.h"
#include "core/qaws_cubic_poly_core.h"

/* Pre-computed span coefficients uploaded from C side */
StructuredBuffer<float> span_coeffs : register(t0);  /* span_count * 8 floats (2D) */
RWTexture1D<float2> results : register(u0);

cbuffer Params : register(b0) {
    uint span_count;
};

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID) {
    /* Map thread to span + local_t */
    float global_t = id.x / 1024.0 * (float)span_count;
    uint span = min((uint)global_t, span_count - 1);
    float local_t = global_t - (float)span;

    /* Load pre-computed a,b,c,d coefficients for this span */
    uint base = span * 8;  /* 2 dims × 4 coeffs */
    qaws_vec2 a = { span_coeffs[base + 0], span_coeffs[base + 4] };
    qaws_vec2 b = { span_coeffs[base + 1], span_coeffs[base + 5] };
    qaws_vec2 c = { span_coeffs[base + 2], span_coeffs[base + 6] };
    qaws_vec2 d = { span_coeffs[base + 3], span_coeffs[base + 7] };

    qaws_eval_2d result = qaws_cubic_eval_2d(a, b, c, d, local_t, 1);
    results[id.x] = float2(result.position.x, result.position.y);
}
```

**GLSL (compute shader):**
```glsl
#include "qaws_glsl.h"
#include "core/qaws_decasteljau_core.h"

layout(std430, binding = 0) readonly buffer CPBuffer { float cp[]; };
layout(std430, binding = 1) buffer OutBuffer { float results[]; };

layout(local_size_x = 256) in;
void main() {
    float t = float(gl_GlobalInvocationID.x) / 1024.0;
    qaws_vec2 p0 = qaws_vec2(cp[0], cp[1]);
    qaws_vec2 p1 = qaws_vec2(cp[2], cp[3]);
    qaws_vec2 p2 = qaws_vec2(cp[4], cp[5]);
    qaws_vec2 p3 = qaws_vec2(cp[6], cp[7]);
    qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);
    results[gl_GlobalInvocationID.x * 2]     = pos.x;
    results[gl_GlobalInvocationID.x * 2 + 1] = pos.y;
}
```

**Halide (C++ DSL):**
```cpp
#include "qaws_halide.h"
#include "core/qaws_decasteljau_core.h"

// Buffers = ImageParam, scalars = Param<>
Halide::ImageParam cp_buf(Halide::Float(32), 1, "cp");
Halide::Param<float> t_start("t_start");
Halide::Param<float> t_step("t_step");

Halide::Func bezier("bezier");
Halide::Var i("i");

// qaws_scalar = Halide::Expr, so this builds an expression graph
qaws_vec2 p0 = { cp_buf(0), cp_buf(1) };
qaws_vec2 p1 = { cp_buf(2), cp_buf(3) };
qaws_vec2 p2 = { cp_buf(4), cp_buf(5) };
qaws_vec2 p3 = { cp_buf(6), cp_buf(7) };
qaws_scalar t = t_start + Halide::cast<float>(i) * t_step;

qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);
bezier(i) = Halide::Tuple(pos.x, pos.y);
bezier.vectorize(i, 8);
```

---

## Core Headers to Create

Based on auditing the actual eval functions in the codebase, the inner math falls into these distinct patterns:

| Header | Source | Algorithm | Shader priority |
|--------|--------|-----------|-----------------|
| `qaws_decasteljau_core.h` | `qaws_internal_basis.c`, `qaws_bezier.c` | De Casteljau triangular reduction, derivative control points. Cubic specialization + arbitrary degree, 2D/3D. | High |
| `qaws_cubic_poly_core.h` | `qaws_hermite.c`, `qaws_catmull_rom.c`, `qaws_trajectory.c` | Evaluate `a*t³ + b*t² + c*t + d` per dimension with up to D3. All three families pre-compute these coefficients at creation time, then eval is identical. | High |
| `qaws_bspline_basis_core.h` | `qaws_internal_basis.c` | B-spline basis functions + derivatives (Cox-de Boor). Knot span search. **Note**: current `bspline_eval_span` uses `malloc` for `ders_buf` — must be replaced with a stack buffer bounded by `(QAWS_CORE_MAX_DEGREE+1)*(max_deriv+1)`. | High |
| `qaws_bspline_eval_core.h` | `qaws_bspline.c` | Weighted sum: `sum(basis[j] * cp[span-degree+j])`. Builds on basis core. 2D/3D, position + D1/D2/D3. | High |
| `qaws_nurbs_eval_core.h` | `qaws_nurbs.c` | NURBS quotient rule: `A^(k)`, `W^(k)` accumulation then `C' = (A' - C*W') / W`. Builds on basis core. | High |
| `qaws_rational_bezier_core.h` | `qaws_rational_bezier.c` | Homogeneous De Casteljau (dim+1 components) + quotient rule for derivatives. Builds on decasteljau core. | High |
| `qaws_horner_core.h` | `qaws_polynomial.c` | Horner scheme with stride for monomial-basis polynomials. `poly_horner`, `poly_horner_d1`, `poly_horner_d2`, `poly_horner_d3`. | Medium |
| `qaws_arc_core.h` | `qaws_arc.c` | `cos(θ), sin(θ)` with plane basis vectors (axis_u, axis_v). Fixed per-segment data, no variable buffers. | Medium |
| `qaws_surface_core.h` | `qaws_surface_bezier.c`, `qaws_surface_bspline.c` | Tensor-product: evaluate basis along u then v. | High |

### Family Classification by Buffer Needs

Not all families need the buffer abstraction. This is critical for deciding what to extract:

| Category | Families | Buffer access at eval time? |
|----------|----------|---------------------------|
| **No buffer** | Bezier, Hermite, Catmull-Rom, Trajectory, Polynomial, Arc, Clothoid | No. All data is either passed by value (Bezier CPs) or pre-computed into per-span coefficients at creation. Core functions take small local arrays only. |
| **Local window of global buffer** | B-spline, NURBS | Yes. At eval time, must find knot span, then index into control points, knots, and weights arrays. But only `degree+1` CPs and `2*(degree+1)` knots are accessed per eval. |
| **Homogeneous buffer** | Rational Bezier | Sort of. Uses homogeneous weighted points (pre-computed at creation). Single span, so all points are local. Fits in a small array for `degree ≤ QAWS_CORE_MAX_DEGREE`. |
| **Nested dispatch** | Composite, Subdivision, Yuksel | These delegate to other curve types or have complex per-subcurve structures. Not candidates for core headers — remain C-only. |
| **Analytic** | Clothoid | Fresnel integral approximation via numerical integration. Could be a core header but low priority. |

**Takeaway**: Only **B-spline** and **NURBS** require the two-layer buffer pattern. All other families can be pure value-in/value-out core headers with no buffer concerns.

### Constraints for Shader Compatibility

Core headers must follow these rules:

- **No pointers.** Pass and return by value. Use structs for multi-value returns.
- **No dynamic arrays.** Fixed-size stack arrays bounded by `QAWS_CORE_MAX_DEGREE`.
- **No heap allocation.** No `malloc`/`free`. **NOTE**: Current `bspline_eval_span_2d/3d` mallocs `ders_buf` — this must be changed to a stack buffer.
- **No recursion.** Iterative loops with fixed upper bounds.
- **No `memcpy`, `memset`, `printf`, `assert`.** Use explicit assignment loops or `QAWS_ZERO` initialization.
- **No `unsigned int` in GLSL.** Use `int` for loop counters (GLSL has limited `uint` support in older versions).
- **Use `QAWS_SELECT` for Halide.** Halide cannot use C ternary operators on `Halide::Expr` for conditional selection — must use `Halide::select()`.
- **Use `QAWS_LITERAL` for all float constants.** Bare `1.0` works in C/HLSL/GLSL but not Halide (needs `make_const`).
- **Fixed-size array arguments.** In HLSL/GLSL, array parameters must have compile-time-known sizes. Use `QAWS_CORE_MAX_POINTS` as the array dimension, with `degree` passed separately to bound the loop.

---

## Architecture: Prepare + Eval

The architecture collapses to **two layers**. No vtable, no opaque objects, no lifecycle management for GPU/shader users.

| Layer | Where | What |
|-------|-------|------|
| **Prepare functions** | C only (`qaws_prepare.h` / `qaws_prepare.c`) | Validate input, compute derived data (span coefficients, homogeneous points, knot params). Output flat arrays the caller stores however they want. |
| **Eval functions** | All backends (`core/*_core.h`, header-only) | Pure math. Takes local data by value or small fixed-size arrays, returns result. No pointers, no heap, no lifecycle. |

The existing C runtime API (`qaws_curve*`, vtable dispatch, allocators) stays as an **optional convenience wrapper** — it calls prepare internally at creation time and eval internally at evaluation time. Users who only want GPU evaluation bypass it entirely.

### Prepare Functions (C Only)

Each family has a prepare function that takes raw geometry input and outputs pre-computed flat data suitable for upload to any backend.

```c
#include "qaws_prepare.h"

/* --- Hermite --- */
/* Input: points + tangents. Output: a,b,c,d coefficient vectors per span. */
qaws_status qaws_hermite_prepare_2d(
    const qaws_vec2* points,
    const qaws_vec2* tangents,
    unsigned int point_count,
    qaws_vec2* out_a,          /* [span_count] */
    qaws_vec2* out_b,          /* [span_count] */
    qaws_vec2* out_c,          /* [span_count] */
    qaws_vec2* out_d);         /* [span_count] */

/* --- Catmull-Rom --- */
/* Input: interpolation points. Output: same a,b,c,d coefficients. */
qaws_status qaws_catmull_rom_prepare_2d(
    const qaws_vec2* points,
    unsigned int point_count,
    qaws_parameterization param,
    int closed,
    qaws_vec2* out_a,
    qaws_vec2* out_b,
    qaws_vec2* out_c,
    qaws_vec2* out_d,
    unsigned int* out_span_count);

/* --- Trajectory --- */
/* Input: key positions + times. Output: same a,b,c,d coefficients. */
qaws_status qaws_trajectory_prepare_2d(
    const qaws_vec2* positions,
    const qaws_scalar* times,
    unsigned int key_count,
    unsigned int degree,
    int closed,
    qaws_vec2* out_a,
    qaws_vec2* out_b,
    qaws_vec2* out_c,
    qaws_vec2* out_d,
    unsigned int* out_span_count);

/* --- Bezier --- */
/* No prepare needed for position-only eval (CPs are the data).
   For derivative eval, pre-compute derivative control points. */
qaws_status qaws_bezier_prepare_deriv_2d(
    const qaws_vec2* control_points,
    unsigned int degree,
    qaws_vec2* out_d1_cp,      /* [degree] */
    qaws_vec2* out_d2_cp,      /* [degree-1] */
    qaws_vec2* out_d3_cp);     /* [degree-2] */

/* --- Rational Bezier --- */
/* Output: homogeneous weighted points (w*x, w*y, w) */
qaws_status qaws_rational_bezier_prepare_2d(
    const qaws_vec2* control_points,
    const qaws_scalar* weights,
    unsigned int degree,
    qaws_scalar* out_weighted_cp);  /* [(degree+1) * 3] */

/* --- B-spline / NURBS --- */
/* No prepare needed — knots, CPs, and weights are already flat arrays.
   The eval function takes them directly. */

/* --- Polynomial --- */
/* No prepare needed — coefficients are already flat. */

/* --- Arc --- */
/* Output: per-segment axis vectors and angle range */
qaws_status qaws_arc_prepare_2d(
    qaws_vec2 center,
    qaws_scalar radius_x,
    qaws_scalar radius_y,
    qaws_scalar start_angle,
    qaws_scalar end_angle,
    qaws_vec2* out_axis_u,
    qaws_vec2* out_axis_v,
    qaws_scalar* out_theta_start,
    qaws_scalar* out_theta_range);
```

All prepare functions output **flat arrays of scalars or vec2/vec3**. No opaque structs, no pointers-to-pointers. The caller owns the memory and decides where to store it (CPU array, GPU buffer, shared memory, etc.).

### Eval Functions (All Backends)

These are the `core/*_core.h` headers described above. They take pre-computed data by value or in small fixed-size arrays and return results. They compile on C, HLSL, GLSL, and Halide.

```c
/* core/qaws_cubic_poly_core.h — shared by Hermite, Catmull-Rom, Trajectory */
QAWS_INLINE qaws_eval_2d qaws_cubic_eval_2d(
    qaws_vec2 a, qaws_vec2 b, qaws_vec2 c, qaws_vec2 d,
    qaws_scalar t, int eval_flags);

/* core/qaws_decasteljau_core.h — Bezier */
QAWS_INLINE qaws_vec2 qaws_bezier3_eval_2d(
    qaws_vec2 p0, qaws_vec2 p1, qaws_vec2 p2, qaws_vec2 p3,
    qaws_scalar t);

/* core/qaws_bspline_eval_core.h — B-spline */
QAWS_INLINE qaws_vec3 qaws_bspline_eval_3d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree, int span, qaws_scalar t);

/* core/qaws_arc_core.h — Arc */
QAWS_INLINE qaws_vec2 qaws_arc_eval_2d(
    qaws_vec2 axis_u, qaws_vec2 axis_v,
    qaws_scalar theta_start, qaws_scalar theta_range,
    qaws_scalar t);
```

### Workflow

```
  C side (prepare)              Any backend (eval)
  ──────────────────            ────────────────────────────
  qaws_catmull_rom_prepare_2d()
      → a[], b[], c[], d[]  ──upload──→  StructuredBuffer / SSBO / ImageParam
                                              │
                                              ▼
                               thread: load a,b,c,d for span
                                       qaws_cubic_eval_2d(a,b,c,d,t)
                                       write result
```

1. **CPU (C only)**: Call prepare function → get flat arrays
2. **Upload**: Copy flat arrays to GPU buffer — or keep on CPU
3. **Eval (any backend)**: Call eval function with local data

### How the Existing C Runtime Wraps This

The existing `qaws_curve*` API stays unchanged for CPU-only users. Internally, it:

1. **At creation** (`qaws_curve_create_*`): Calls the prepare function, stores the result in the opaque impl struct, sets up the vtable.
2. **At eval** (`qaws_curve_evaluate_*`): Finds the span, extracts local data from the impl, calls the core eval function.
3. **At destruction** (`qaws_curve_destroy`): Frees the impl struct.

This is a convenience wrapper for users who want managed lifecycle and polymorphic dispatch. It adds zero overhead beyond the indirection. GPU/shader users bypass it entirely.

### Before / After Examples

#### Hermite / Catmull-Rom / Trajectory

All three families pre-compute `a,b,c,d` coefficients per span. Today this code is duplicated 3 times in the C runtime. With the new architecture, the prepare step outputs flat arrays and the eval is a single shared function.

**CPU side** (prepare once, upload):
```c
qaws_vec2 a[MAX_SPANS], b[MAX_SPANS], c[MAX_SPANS], d[MAX_SPANS];
unsigned int span_count;
qaws_catmull_rom_prepare_2d(points, 10, QAWS_PARAMETERIZATION_CENTRIPETAL, 0,
                            a, b, c, d, &span_count);
/* Upload a,b,c,d to GPU as a StructuredBuffer of span_count * 8 floats */
```

**GPU side** (eval per thread):
```hlsl
#include "qaws_hlsl.h"
#include "core/qaws_cubic_poly_core.h"

StructuredBuffer<float> span_data : register(t0);  /* span_count * 8 floats (2D) */
[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID) {
    uint span = /* ... map thread to span ... */;
    float local_t = /* ... */;
    uint base = span * 8;
    qaws_vec2 a = { span_data[base+0], span_data[base+4] };
    qaws_vec2 b = { span_data[base+1], span_data[base+5] };
    qaws_vec2 c = { span_data[base+2], span_data[base+6] };
    qaws_vec2 d = { span_data[base+3], span_data[base+7] };
    qaws_eval_2d result = qaws_cubic_eval_2d(a, b, c, d, local_t, 1);
    /* write result.position */
}
```

#### B-spline (malloc bug fix included)

**Current problem** (`qaws_bspline.c` today — allocates with `malloc` in the hot eval path):
```c
ders_buf = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (max_deriv_order + 1) * stride);
/* ... eval ... */
free(ders_buf);
```

This is unnecessary — `(max_deriv+1) * (degree+1)` is bounded by `4 * (QAWS_CORE_MAX_DEGREE+1)` = 68 scalars = 272 bytes on the stack.

**After** — the core eval function uses a stack buffer and works on all backends:
```c
/* core/qaws_bspline_eval_core.h */
QAWS_INLINE qaws_eval_2d qaws_bspline_eval_2d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree, int span, qaws_scalar t, int eval_flags)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS]; /* stack */
    /* ... Cox-de Boor + weighted sum ... */
}
```

B-spline and NURBS need no prepare function — their raw data (knots, CPs, weights) is already flat. The caller uploads those arrays directly and does buffer extraction + span finding per-backend (see the Buffer Abstraction section below).

---

## SIMD Layer (C Backend Only)

SIMD acceleration is a C-only concern, invisible to HLSL/GLSL/Halide users. It accelerates `qaws_curve_evaluate_batch_2d/3d` by processing multiple parameter values in parallel.

### ISA Detection (`simd/qaws_simd.h`)

```c
#if defined(__AVX2__)
#  include "qaws_simd_avx2.h"
#  define QAWS_SIMD_WIDTH_F32 8
#  define QAWS_SIMD_WIDTH_F64 4
#elif defined(__SSE2__) || defined(_M_X64)
#  include "qaws_simd_sse2.h"
#  define QAWS_SIMD_WIDTH_F32 4
#  define QAWS_SIMD_WIDTH_F64 2
#elif defined(__aarch64__) || defined(__ARM_NEON)
#  include "qaws_simd_neon.h"
#  define QAWS_SIMD_WIDTH_F32 4
#  define QAWS_SIMD_WIDTH_F64 2
#else
#  include "qaws_simd_scalar.h"
#  define QAWS_SIMD_WIDTH_F32 1
#  define QAWS_SIMD_WIDTH_F64 1
#endif
```

### Type-Generic Aliases (following Alwan's `alwan_map_internal.h` pattern)

```c
#if QAWS_SCALAR_IS_FLOAT
# define QAWS_SIMD_WIDTH   QAWS_SIMD_WIDTH_F32
# define qaws_simd         qaws_simd_f32
# define qaws_simd_set1    qaws_simd_f32_set1
# define qaws_simd_load    qaws_simd_f32_load
# define qaws_simd_store   qaws_simd_f32_store
# define qaws_simd_add     qaws_simd_f32_add
# define qaws_simd_sub     qaws_simd_f32_sub
# define qaws_simd_mul     qaws_simd_f32_mul
# define qaws_simd_fma     qaws_simd_f32_fma
#else
# define QAWS_SIMD_WIDTH   QAWS_SIMD_WIDTH_F64
# define qaws_simd         qaws_simd_f64
  /* ... same pattern ... */
#endif
```

The batch `.c` files use these generic aliases so the same code compiles for f32 or f64 SIMD.

---

## Implementation Plan

### Phase 1: Platform + Types `[S]`

Create `qaws_platform.h` and `qaws_core_types.h`. The existing `qaws_types.h` already defines `qaws_scalar`, `qaws_vec2`, `qaws_vec3` — the `qaws_core_types.h` mirrors these using `QAWS_TYPE_DEF` for cross-backend compatibility. The C runtime continues using `qaws_types.h`; shaders use `qaws_core_types.h`.

Also, replace `(qaws_scalar)3.0`, `(qaws_scalar)sqrt(...)`, etc. patterns scattered across existing `.c` files with `QAWS_LITERAL(3.0)`, `QAWS_SQRT(...)`. Pure refactor, no functional change.

### Phase 2: Core Eval Headers `[M]`

Extract inner math from existing family eval functions into `core/*_core.h`. Order by dependency:

1. **`qaws_decasteljau_core.h`** — from `qaws_internal_decasteljau()` + `qaws_internal_bezier_derivative_points()`. Serves: Bezier, Rational Bezier.

2. **`qaws_cubic_poly_core.h`** — the `a*t³ + b*t² + c*t + d` eval duplicated in Hermite, Catmull-Rom, Trajectory. Immediate code deduplication.

3. **`qaws_bspline_basis_core.h`** — Cox-de Boor basis + derivatives. Stack buffer instead of `malloc`. Serves: B-spline, NURBS, surfaces.

4. **`qaws_bspline_eval_core.h`** — weighted sum loop. Depends on (3).

5. **`qaws_nurbs_eval_core.h`** — quotient rule. Depends on (3).

6. **`qaws_horner_core.h`** — Horner scheme for monomial-basis polynomials.

7. **`qaws_arc_core.h`** — parametric `cos(θ)/sin(θ)` with plane basis vectors.

8. **`qaws_surface_core.h`** — tensor-product evaluation.

After extraction, the existing `.c` files `#include` the core headers and call them internally.

### Phase 3: Prepare Functions `[M]`

Create `qaws_prepare.h` / `qaws_prepare.c` with prepare functions for each family:

- `qaws_hermite_prepare_2d/3d` — compute `a,b,c,d` span coefficients from points + tangents
- `qaws_catmull_rom_prepare_2d/3d` — compute span coefficients from interpolation points
- `qaws_trajectory_prepare_2d/3d` — compute span coefficients from key positions + times
- `qaws_bezier_prepare_deriv_2d/3d` — compute derivative control points
- `qaws_rational_bezier_prepare_2d/3d` — compute homogeneous weighted points
- `qaws_arc_prepare_2d/3d` — compute axis vectors and angle range

Each prepare function extracts the computation currently done inside `qaws_curve_create_*()`. The existing create functions are refactored to call prepare internally, so the C runtime API is unchanged.

### Phase 4: Bootstrap Headers + Tests `[S]`

Create `qaws_hlsl.h`, `qaws_glsl.h`, `qaws_halide.h`. Write example shaders:
- `examples/eval_bezier.hlsl` — Bezier evaluation compute shader (DXC compilation test)
- `examples/eval_bezier.comp` — GLSL compute shader (glslangValidator test)
- `examples/eval_bspline.hlsl` — B-spline with buffer extraction pattern
- `examples/eval_cubic_poly.hlsl` — Hermite/CatmullRom/Trajectory with pre-computed coefficients
- CMake targets: `add_test(NAME hlsl_compile COMMAND dxc ...)`, `add_test(NAME glsl_compile COMMAND glslangValidator ...)`

### Phase 5: SIMD Batch Evaluation `[M]`

Implement `simd/qaws_simd_*.h` + `batch/qaws_batch_*.c`. Accelerate `qaws_curve_evaluate_batch_2d/3d` for Bezier and B-Spline families. Other families fall back to existing scalar loop.

### Priority

**Phase 1 → 2 → 3 → 4** delivers the prepare+eval shader library.
**Phase 5** delivers CPU SIMD acceleration.

All phases are additive — the existing C API is never broken. Phase 2 fixes the `malloc` in B-spline eval (performance improvement even without multi-backend). Phase 3 enables GPU workflows without requiring the C runtime.

---

## Build System Integration

```cmake
option(QAWS_SCALAR_IS_FLOAT "Use float instead of double" OFF)

# Core headers are header-only, always installed
install(FILES
    src/qaws/qaws_platform.h
    src/qaws/qaws_core_types.h
    src/qaws/qaws_hlsl.h
    src/qaws/qaws_glsl.h
    src/qaws/qaws_halide.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qaws)

install(DIRECTORY src/qaws/core/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qaws/core)

# SIMD is C-only, compile-time ISA detection
if(QAWS_ENABLE_SIMD)
    include(CheckCCompilerFlag)
    check_c_compiler_flag("-mavx2" QAWS_HAS_AVX2)
    if(QAWS_HAS_AVX2)
        target_compile_options(qaws PRIVATE -mavx2 -mfma)
    endif()
endif()
```

Users who only want the shader library can just copy `qaws_platform.h`, `qaws_core_types.h`, `qaws_hlsl.h`/`qaws_glsl.h`/`qaws_halide.h`, and `core/*.h` into their shader include path. No linking required.

---

## Buffer / Variable-Size Data Abstraction

### The Problem

Alwan avoids this entirely — all core functions are value-in/value-out on fixed-size 3-component color structs. No core function ever reads from or writes to a variable-length buffer.

Qaws cannot do this. Curve algorithms need variable-length data: control points, knot vectors, weights. A B-spline with 100 control points and a B-spline with 5 control points use the same algorithm, but the data size varies. Each backend represents "a buffer of N elements" differently:

| Backend | Read-only buffer | Read-write buffer |
|---------|-----------------|-------------------|
| **C**       | `const qaws_scalar*` (pointer) | `qaws_scalar*` |
| **HLSL**    | `StructuredBuffer<T>`, `Buffer<T>`, `Texture1D<T>` | `RWStructuredBuffer<T>`, `RWBuffer<T>`, `RWTexture1D<T>`, `RWTexture2D<T>` — where `T` is `float`, `float2`, `float3`, or `float4` |
| **GLSL**    | `readonly buffer { float data[]; }` (SSBO), `samplerBuffer` (TBO) | `buffer { float data[]; }` |
| **Halide**  | `Halide::Func` (pure), `Halide::ImageParam` | `Halide::Func` (update definition) |

There is no single abstraction that maps cleanly to all of these.

### Key Insight: Most Families Don't Need Buffers At All

From auditing the actual eval functions:

| Family | Eval-time data access | Buffer needed? |
|--------|----------------------|----------------|
| **Bezier** | `impl->control_points` — single span, all CPs are local | No (degree ≤ 16 → ≤ 51 scalars) |
| **Hermite** | `impl->span_coeffs[span_index * dim * 4]` — 4 pre-computed coefficients per dim per span | No (8 or 12 scalars) |
| **Catmull-Rom** | `impl->segment_coeffs[span_index * dim * 4]` — identical to Hermite | No (8 or 12 scalars) |
| **Trajectory** | `impl->span_coeffs[span_index * dim * 4]` — identical to Hermite | No (8 or 12 scalars) |
| **Polynomial** | `impl->coefficients` — single span, Horner over all coefficients | No (degree ≤ 16 → ≤ 51 scalars) |
| **Arc** | `impl->segments[span_index]` — fixed-size `qaws_arc_segment` struct | No (all scalars) |
| **Rational Bezier** | `impl->weighted_points` — single span, all homogeneous CPs | No (degree ≤ 16 → ≤ 68 scalars) |
| **B-spline** | `impl->control_points[knot_span - degree + j]` + `impl->knots` — indexes into potentially large arrays | **Yes** — window of `degree+1` CPs + `2*(degree+1)` knots |
| **NURBS** | Same as B-spline + `impl->weights[idx]` | **Yes** — window of CPs + knots + weights |

**Only B-spline and NURBS require the buffer extraction pattern.** All other families work with small local data that fits entirely on the stack. This means the two-layer architecture is only needed for 2 out of 9 eval-capable families.

For B-spline/NURBS, the local window is bounded by `QAWS_CORE_MAX_DEGREE`:
- `degree+1` CPs × 3 components = 51 scalars max
- `2*(degree+1)` knots = 34 scalars max
- `degree+1` weights = 17 scalars max
- **Total**: ~102 scalars = 408 bytes — easily fits in shader registers/local memory

### Solution: Two-Layer Architecture

**Layer 1 — Core functions (all backends):** Operate on small, fixed-size local data only. No buffer access at all. This is what compiles everywhere.

```c
/* Core: always works with local data passed by value or in small arrays */
QAWS_INLINE qaws_vec3 qaws_bspline_eval_3d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],  /* pre-extracted */
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t)
{
    /* Pure math — no buffer access, no pointers beyond the stack arrays */
}
```

**Layer 2 — Buffer access (backend-specific):** Extracts the local window from whatever storage the backend uses, then calls the core function.

This layer is NOT shared across backends. Each backend writes its own extraction code in its own idiom. The core headers provide a **span-finding function** (pure math, no buffer access) so the caller knows which indices to extract.

### Span Finding

The knot span search (`knots[span] <= t < knots[span+1]`) is needed only for B-spline and NURBS. Its implementation differs by backend:

**C / HLSL / GLSL** — binary search with a `while` loop (or bounded `for` loop for GLSL):

```c
QAWS_INLINE int qaws_find_span(
    qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2],
    int knot_count,
    int degree,
    int num_cp,
    qaws_scalar t)
{
    int low = degree;
    int high = num_cp;
    int mid, iter;

    if (t >= knots[high]) return high - 1;
    if (t <= knots[low])  return low;

    /* Bounded binary search — log2(QAWS_CORE_MAX_POINTS*2) ≈ 6 iterations max.
     * Uses a for loop with fixed upper bound (required for GLSL). */
    for (iter = 0; iter < 20 && high - low > 1; iter++) {
        mid = (low + high) / 2;
        if (t < knots[mid])
            high = mid;
        else
            low = mid;
    }
    return low;
}
```

Note: uses a bounded `for` loop instead of `while` — GLSL requires it. The bound of 20 iterations is sufficient for `log2(QAWS_CORE_MAX_POINTS * 2)` = ~6 iterations.

**Halide** — cannot use imperative loops at all. Two options:

1. **Pre-compute on CPU** (recommended): The C runtime already finds the knot span in `qaws_internal_find_knot_span()`. Pass the span index as a `Halide::Param<int>` to the shader. This is what you'd do for single-curve evaluation.

2. **RDom reduction** (for batch evaluation across many `t` values): Use Halide's `RDom` over the knot array with `argmin`/`select`, similar to the Mandelbrot `argmin(magnitude(...) < 4)` pattern:

```cpp
// Find span: last knot index where knots[k] <= t
Halide::ImageParam knot_buf(Halide::Float(32), 1, "knots");
Halide::Param<int> num_knots("num_knots");
Halide::Param<int> degree("degree");

Halide::Func span("span");
Halide::Var i("i");  // index over t values to evaluate

Halide::RDom k(0, num_knots);
// For each t(i), find the largest k where knots[k] <= t(i) and k >= degree
Halide::Expr valid = (knot_buf(k) <= t_values(i)) & (k >= degree);
span(i) = argmax(select(valid, k, -1))[0];
// Clamp: if span >= num_cp, use num_cp - 1
span(i) = clamp(span(i), degree, num_cp - 1);
```

The RDom approach is O(n) per eval (linear scan) rather than O(log n) (binary search), but Halide's vectorization and parallelism compensate. For most practical curves (< 100 knots), the linear scan is fast enough.

### Per-Backend Buffer Access Examples

#### C Backend

```c
#include "core/qaws_cox_deboor_core.h"

/* C: just pointer indexing */
void eval_bspline_3d(
    const qaws_scalar* all_cp,      /* pointer to N*3 scalars */
    const qaws_scalar* all_knots,   /* pointer to M scalars */
    int degree, int num_cp,
    qaws_scalar t,
    qaws_vec3* out)
{
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3];
    int num_knots = num_cp + degree + 1;
    int i, span;

    /* Load all knots (for small curves) or a window */
    for (i = 0; i < num_knots && i < QAWS_CORE_MAX_POINTS * 2; i++)
        local_knots[i] = all_knots[i];

    span = qaws_find_span(local_knots, num_knots, degree, t);

    /* Extract local control points: span-degree .. span */
    for (i = 0; i <= degree; i++) {
        int idx = (span - degree + i) * 3;
        local_cp[i * 3]     = all_cp[idx];
        local_cp[i * 3 + 1] = all_cp[idx + 1];
        local_cp[i * 3 + 2] = all_cp[idx + 2];
    }

    *out = qaws_bspline_eval_3d(local_cp, local_knots, degree, span, t);
}
```

#### HLSL Backend

```hlsl
#include "qaws_hlsl.h"
#include "core/qaws_cox_deboor_core.h"

StructuredBuffer<float>  g_cp     : register(t0);  /* N*3 floats */
StructuredBuffer<float>  g_knots  : register(t1);  /* M floats */
RWStructuredBuffer<float> g_out   : register(u0);  /* output positions */

cbuffer CurveParams : register(b0) {
    int degree;
    int num_cp;
    int num_knots;
};

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID)
{
    float t = id.x / float(1024);

    /* Load knots into local array */
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2];
    int k;
    for (k = 0; k < num_knots && k < QAWS_CORE_MAX_POINTS * 2; k++)
        local_knots[k] = g_knots[k];

    int span = qaws_find_span(local_knots, num_knots, degree, t);

    /* Extract local control points from StructuredBuffer */
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3];
    int i;
    for (i = 0; i <= degree; i++) {
        int idx = (span - degree + i) * 3;
        local_cp[i * 3]     = g_cp[idx];
        local_cp[i * 3 + 1] = g_cp[idx + 1];
        local_cp[i * 3 + 2] = g_cp[idx + 2];
    }

    qaws_vec3 pos = qaws_bspline_eval_3d(
        local_cp, local_knots, degree, span, t);

    g_out[id.x * 3]     = pos.x;
    g_out[id.x * 3 + 1] = pos.y;
    g_out[id.x * 3 + 2] = pos.z;
}
```

#### GLSL Backend

```glsl
#include "qaws_glsl.h"
#include "core/qaws_cox_deboor_core.h"

layout(std430, binding = 0) readonly buffer CPBuffer   { float cp[];   };
layout(std430, binding = 1) readonly buffer KnotBuffer { float knots[]; };
layout(std430, binding = 2) buffer OutBuffer           { float out_pos[]; };

uniform int degree;
uniform int num_cp;
uniform int num_knots;

layout(local_size_x = 256) in;
void main()
{
    float t = float(gl_GlobalInvocationID.x) / 1024.0;

    /* Load knots into local array */
    float local_knots[QAWS_CORE_MAX_POINTS * 2];
    for (int k = 0; k < num_knots && k < QAWS_CORE_MAX_POINTS * 2; k++)
        local_knots[k] = knots[k];

    int span = qaws_find_span(local_knots, num_knots, degree, t);

    /* Extract local control points from SSBO */
    float local_cp[QAWS_CORE_MAX_POINTS * 3];
    for (int i = 0; i <= degree; i++) {
        int idx = (span - degree + i) * 3;
        local_cp[i * 3]     = cp[idx];
        local_cp[i * 3 + 1] = cp[idx + 1];
        local_cp[i * 3 + 2] = cp[idx + 2];
    }

    qaws_vec3 pos = qaws_bspline_eval_3d(
        local_cp, local_knots, degree, span, t);

    out_pos[gl_GlobalInvocationID.x * 3]     = pos.x;
    out_pos[gl_GlobalInvocationID.x * 3 + 1] = pos.y;
    out_pos[gl_GlobalInvocationID.x * 3 + 2] = pos.z;
}
```

#### Halide Backend

Halide is fundamentally different from C/HLSL/GLSL: code builds expression DAGs, not imperative loops. The core functions still work because `qaws_scalar` becomes `Halide::Expr`, so `qaws_bezier3_eval_2d()` builds a Halide expression graph rather than computing a value.

For families without buffer access (Bezier, Hermite, etc.), this works directly:

```cpp
#include "qaws_halide.h"
#include "core/qaws_decasteljau_core.h"

// Buffers are ImageParam (1D arrays), scalars are Param<>
Halide::ImageParam cp_buf(Halide::Float(32), 1, "cp");  // 4*2 floats for cubic 2D
Halide::Param<float> t_start("t_start");  // scalar → Param<>
Halide::Param<float> t_step("t_step");    // scalar → Param<>

Halide::Func eval("eval");
Halide::Var i("i");

// Build control points from buffer — these are Halide::Expr, not float values
qaws_vec2 p0 = { cp_buf(0), cp_buf(1) };
qaws_vec2 p1 = { cp_buf(2), cp_buf(3) };
qaws_vec2 p2 = { cp_buf(4), cp_buf(5) };
qaws_vec2 p3 = { cp_buf(6), cp_buf(7) };

// t is a Halide::Expr dependent on the Var i
qaws_scalar t = t_start + Halide::cast<float>(i) * t_step;

// This call builds a Halide expression graph — no actual computation yet
qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);

eval(i) = Halide::Tuple(pos.x, pos.y);

// Schedule and compile
eval.vectorize(i, 8);
eval.compile_jit();
```

For B-spline/NURBS (buffer access needed), the local array pattern does **not** translate to Halide — you cannot fill a Halide array in a loop. Instead, unroll the buffer extraction for bounded degree, or use Halide `RDom` (reduction domain) for the basis function accumulation. This is Halide-specific glue code, not shared.

### Why Not a `QAWS_READ(buf, index)` Macro?

It's tempting to define a single macro:
```c
#if QAWS_BACKEND == QAWS_BACKEND_C
#  define QAWS_READ(buf, i) ((buf)[i])
#elif QAWS_BACKEND == QAWS_BACKEND_HLSL
#  define QAWS_READ(buf, i) ((buf)[i])  /* StructuredBuffer supports [] */
/* ... */
#endif
```

This fails because:

1. **Buffer declaration differs radically.** C uses `const float*`, HLSL uses `StructuredBuffer<float>`, GLSL uses `layout(std430) readonly buffer { float data[]; }`. A macro can hide the read, but not the declaration.
2. **HLSL has multiple buffer types** (`StructuredBuffer`, `Buffer`, `Texture1D`, `ByteAddressBuffer`) with different performance characteristics. The user must choose.
3. **GLSL SSBOs require layout qualifiers** that cannot be macro'd away.
4. **Halide buffer reads build expression DAGs**, not values. Putting `Halide::Func` reads inside a loop creates IR nodes, not iterations. The entire loop structure must be rethought for Halide.
5. **Core functions would become buffer-coupled.** Currently they are pure math on local data — testable, composable, and trivially correct. Adding buffer reads would break this.

The two-layer approach keeps core functions pure and lets each backend handle buffer access in its native idiom. The only shared code is the math.

### Summary

| Layer | Where | What it does |
|-------|-------|-------------|
| **Prepare functions** | C only | Validate input, compute derived data. Output flat arrays. |
| **Eval functions** | All backends (header-only) | Pure math on small local data. No buffer access, no heap, no lifecycle. |
| **Span finding** | All backends (except Halide: CPU pre-compute or RDom) | Binary search on local knot array. B-spline/NURBS only. |
| **Buffer extraction** | Per-backend (user code) | Load local window from native storage, call core eval. B-spline/NURBS only. |
| **C runtime wrapper** | C only (optional) | `qaws_curve*`, vtable dispatch, allocators, `create`/`destroy`. Convenience layer on top of prepare+eval. |

**7 of 9 families** (Bezier, Hermite, Catmull-Rom, Trajectory, Rational Bezier, Polynomial, Arc) need no buffer extraction at eval time — just call prepare on CPU, upload flat data, call eval in the shader.

**2 of 9 families** (B-spline, NURBS) need per-backend buffer extraction (~10 lines). We provide **example shaders** (`examples/eval_bspline.hlsl`, `examples/eval_bspline.comp`) to copy and adapt.

GPU/shader users never touch `qaws_curve*`. CPU-only users can continue using the full C runtime API unchanged — it wraps prepare+eval internally.
