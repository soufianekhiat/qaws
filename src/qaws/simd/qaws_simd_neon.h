#ifndef QAWS_SIMD_NEON_H
#define QAWS_SIMD_NEON_H

#include <arm_neon.h>

typedef float32x4_t qaws_simdf;
/* NEON f64 only available on AArch64; fallback to scalar for double */
typedef double qaws_simdd;

#define qaws_simd_set1_f(v)       vdupq_n_f32(v)
#define qaws_simd_load_f(ptr)     vld1q_f32(ptr)
#define qaws_simd_store_f(ptr, v) vst1q_f32(ptr, v)
#define qaws_simd_add_f(a, b)     vaddq_f32(a, b)
#define qaws_simd_sub_f(a, b)     vsubq_f32(a, b)
#define qaws_simd_mul_f(a, b)     vmulq_f32(a, b)
#define qaws_simd_fmadd_f(a, b, c) vmlaq_f32(c, a, b)
#define qaws_simd_neg_f(a)        vnegq_f32(a)
#if defined(__aarch64__) || defined(_M_ARM64)
# define qaws_simd_div_f(a, b)    vdivq_f32(a, b)
#else
/* ARMv7 has no vector division: lane by lane, correctly rounded (a reciprocal estimate refined is not) */
static inline float32x4_t qaws_simd_neon_div_f(float32x4_t a, float32x4_t b)
{
	float x[4], y[4];
	int i;
	vst1q_f32(x, a);
	vst1q_f32(y, b);
	for (i = 0; i < 4; i++)
		x[i] /= y[i];
	return vld1q_f32(x);
}
# define qaws_simd_div_f(a, b)    qaws_simd_neon_div_f(a, b)
#endif

#define qaws_simd_set1_d(v)       (v)
#define qaws_simd_load_d(ptr)     (*(ptr))
#define qaws_simd_store_d(ptr, v) (*(ptr) = (v))
#define qaws_simd_add_d(a, b)     ((a) + (b))
#define qaws_simd_sub_d(a, b)     ((a) - (b))
#define qaws_simd_mul_d(a, b)     ((a) * (b))
#define qaws_simd_fmadd_d(a, b, c) ((a) * (b) + (c))
#define qaws_simd_div_d(a, b)     ((a) / (b))
#define qaws_simd_neg_d(a)        (-(a))

#endif /* QAWS_SIMD_NEON_H */
