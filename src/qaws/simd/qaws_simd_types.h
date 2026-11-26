#ifndef QAWS_SIMD_TYPES_H
#define QAWS_SIMD_TYPES_H

/*
 * Width-generic type aliases.
 *
 * Each ISA header defines:
 *   qaws_simdf  — SIMD float register type
 *   qaws_simdd  — SIMD double register type (if width > 1)
 *
 *   qaws_simd_set1_f(v)    — broadcast float
 *   qaws_simd_load_f(ptr)  — aligned load
 *   qaws_simd_store_f(ptr, v) — aligned store
 *   qaws_simd_add_f(a, b)  — add
 *   qaws_simd_sub_f(a, b)  — sub
 *   qaws_simd_mul_f(a, b)  — mul
 *   qaws_simd_fmadd_f(a, b, c) — fused multiply-add: a*b + c
 *   qaws_simd_div_f(a, b)  - division, correctly rounded (each lane as the scalar a / b)
 *   qaws_simd_neg_f(a)     - negation (the sign flipped, zeros included)
 *
 * Scalar type selection based on QAWS_SCALAR_IS_FLOAT:
 */

#include "../qaws_platform.h"

#if QAWS_SCALAR_IS_FLOAT
  typedef qaws_simdf qaws_simd;
# define QAWS_SIMD_WIDTH QAWS_SIMD_WIDTH_F32
# define qaws_simd_set1    qaws_simd_set1_f
# define qaws_simd_load    qaws_simd_load_f
# define qaws_simd_store   qaws_simd_store_f
# define qaws_simd_add     qaws_simd_add_f
# define qaws_simd_sub     qaws_simd_sub_f
# define qaws_simd_mul     qaws_simd_mul_f
# define qaws_simd_fmadd   qaws_simd_fmadd_f
#else
  typedef qaws_simdd qaws_simd;
# define QAWS_SIMD_WIDTH QAWS_SIMD_WIDTH_F64
# define qaws_simd_set1    qaws_simd_set1_d
# define qaws_simd_load    qaws_simd_load_d
# define qaws_simd_store   qaws_simd_store_d
# define qaws_simd_add     qaws_simd_add_d
# define qaws_simd_sub     qaws_simd_sub_d
# define qaws_simd_mul     qaws_simd_mul_d
# define qaws_simd_fmadd   qaws_simd_fmadd_d
#endif

#endif /* QAWS_SIMD_TYPES_H */
