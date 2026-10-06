#ifndef QAWS_PLATFORM_H
#define QAWS_PLATFORM_H

/* ===================================================================
 * Backend detection
 * =================================================================== */

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

/* ===================================================================
 * Scalar type and literal macro
 * =================================================================== */

#if QAWS_BACKEND == QAWS_BACKEND_C
# include <math.h>
# ifndef QAWS_SCALAR_IS_FLOAT
#   define QAWS_SCALAR_IS_FLOAT 1
# endif
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
  typedef float qaws_scalar;
# define QAWS_SCALAR_IS_FLOAT 1
# define QAWS_LITERAL(x) (x)
# define QAWS_EPSILON 1e-6f

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

/* int <-> scalar conversions (GLSL has constructor casts only) */
#if QAWS_BACKEND == QAWS_BACKEND_GLSL
# define QAWS_ITOF(n) float(n)
# define QAWS_FTOI(x) int(x)
#elif QAWS_BACKEND == QAWS_BACKEND_HALIDE
# define QAWS_ITOF(n) qaws_scalar((double)(n))
# define QAWS_FTOI(x) Halide::cast<int>(x)
#else
# define QAWS_ITOF(n) ((qaws_scalar)(n))
# define QAWS_FTOI(x) ((int)(x))
#endif

/* ===================================================================
 * Qualifier macros
 * =================================================================== */

#if QAWS_BACKEND == QAWS_BACKEND_C
# define QAWS_INLINE     static inline
# define QAWS_CONSTEXPR  static const
# define QAWS_TYPE_DEF   typedef
# define QAWS_OUT
# define QAWS_INOUT

#elif QAWS_BACKEND == QAWS_BACKEND_HLSL
# define QAWS_INLINE     inline
# define QAWS_CONSTEXPR  static const
# define QAWS_TYPE_DEF   typedef
# define QAWS_OUT        out
# define QAWS_INOUT      inout

#elif QAWS_BACKEND == QAWS_BACKEND_GLSL
# define QAWS_INLINE
# define QAWS_CONSTEXPR  const
# define QAWS_TYPE_DEF
# define QAWS_OUT        out
# define QAWS_INOUT      inout

#elif QAWS_BACKEND == QAWS_BACKEND_HALIDE
# define QAWS_INLINE     inline
# define QAWS_CONSTEXPR  static const
# define QAWS_TYPE_DEF   typedef
# define QAWS_OUT
# define QAWS_INOUT
#endif

/* ===================================================================
 * Math macros
 * =================================================================== */

#if QAWS_BACKEND == QAWS_BACKEND_C
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
# define QAWS_ATAN2(y,x)  atan(y,x)
# define QAWS_POW(x,y)    pow(x,y)
# define QAWS_FMOD(x,y)   mod(x,y)

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

/* ===================================================================
 * Branchless select
 * =================================================================== */

#if QAWS_BACKEND == QAWS_BACKEND_HALIDE
# define QAWS_SELECT(cond, t, f) Halide::select((cond), (t), (f))
#else
# define QAWS_SELECT(cond, t, f) ((cond) ? (t) : (f))
#endif

/* ===================================================================
 * Utility functions
 * =================================================================== */

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

#endif /* QAWS_PLATFORM_H */
