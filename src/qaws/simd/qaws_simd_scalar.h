#ifndef QAWS_SIMD_SCALAR_H
#define QAWS_SIMD_SCALAR_H

/*
 * Scalar fallback — width 1.
 * All operations are plain C arithmetic.
 */

typedef float  qaws_simdf;
typedef double qaws_simdd;

#define qaws_simd_set1_f(v)       (v)
#define qaws_simd_load_f(ptr)     (*(ptr))
#define qaws_simd_store_f(ptr, v) (*(ptr) = (v))
#define qaws_simd_add_f(a, b)     ((a) + (b))
#define qaws_simd_sub_f(a, b)     ((a) - (b))
#define qaws_simd_mul_f(a, b)     ((a) * (b))
#define qaws_simd_fmadd_f(a, b, c) ((a) * (b) + (c))
#define qaws_simd_div_f(a, b)     ((a) / (b))
#define qaws_simd_neg_f(a)        (-(a))

#define qaws_simd_set1_d(v)       (v)
#define qaws_simd_load_d(ptr)     (*(ptr))
#define qaws_simd_store_d(ptr, v) (*(ptr) = (v))
#define qaws_simd_add_d(a, b)     ((a) + (b))
#define qaws_simd_sub_d(a, b)     ((a) - (b))
#define qaws_simd_mul_d(a, b)     ((a) * (b))
#define qaws_simd_fmadd_d(a, b, c) ((a) * (b) + (c))
#define qaws_simd_div_d(a, b)     ((a) / (b))
#define qaws_simd_neg_d(a)        (-(a))

#endif /* QAWS_SIMD_SCALAR_H */
