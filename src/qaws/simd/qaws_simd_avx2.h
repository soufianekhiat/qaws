#ifndef QAWS_SIMD_AVX2_H
#define QAWS_SIMD_AVX2_H

#include <immintrin.h>

typedef __m256  qaws_simdf;
typedef __m256d qaws_simdd;

#define qaws_simd_set1_f(v)       _mm256_set1_ps(v)
#define qaws_simd_load_f(ptr)     _mm256_load_ps(ptr)
#define qaws_simd_store_f(ptr, v) _mm256_store_ps(ptr, v)
#define qaws_simd_add_f(a, b)     _mm256_add_ps(a, b)
#define qaws_simd_sub_f(a, b)     _mm256_sub_ps(a, b)
#define qaws_simd_mul_f(a, b)     _mm256_mul_ps(a, b)
#ifdef __FMA__
# define qaws_simd_fmadd_f(a, b, c) _mm256_fmadd_ps(a, b, c)
#else
# define qaws_simd_fmadd_f(a, b, c) _mm256_add_ps(_mm256_mul_ps(a, b), c)
#endif

#define qaws_simd_div_f(a, b)     _mm256_div_ps(a, b)
#define qaws_simd_neg_f(a)        _mm256_xor_ps(a, _mm256_set1_ps(-0.0f))

#define qaws_simd_set1_d(v)       _mm256_set1_pd(v)
#define qaws_simd_load_d(ptr)     _mm256_load_pd(ptr)
#define qaws_simd_store_d(ptr, v) _mm256_store_pd(ptr, v)
#define qaws_simd_add_d(a, b)     _mm256_add_pd(a, b)
#define qaws_simd_sub_d(a, b)     _mm256_sub_pd(a, b)
#define qaws_simd_mul_d(a, b)     _mm256_mul_pd(a, b)
#ifdef __FMA__
# define qaws_simd_fmadd_d(a, b, c) _mm256_fmadd_pd(a, b, c)
#else
# define qaws_simd_fmadd_d(a, b, c) _mm256_add_pd(_mm256_mul_pd(a, b), c)
#endif
#define qaws_simd_div_d(a, b)     _mm256_div_pd(a, b)
#define qaws_simd_neg_d(a)        _mm256_xor_pd(a, _mm256_set1_pd(-0.0))

#endif /* QAWS_SIMD_AVX2_H */
