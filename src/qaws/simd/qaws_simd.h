#ifndef QAWS_SIMD_H
#define QAWS_SIMD_H

/*
 * ISA detection and dispatch.
 * Selects the widest available SIMD instruction set.
 */

#if defined(__AVX2__)
# include "qaws_simd_avx2.h"
# define QAWS_SIMD_WIDTH_F32 8
# define QAWS_SIMD_WIDTH_F64 4
# define QAWS_SIMD_ISA "AVX2"
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64)
# include "qaws_simd_sse2.h"
# define QAWS_SIMD_WIDTH_F32 4
# define QAWS_SIMD_WIDTH_F64 2
# define QAWS_SIMD_ISA "SSE2"
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
# include "qaws_simd_neon.h"
# define QAWS_SIMD_WIDTH_F32 4
# define QAWS_SIMD_WIDTH_F64 1
# define QAWS_SIMD_ISA "NEON"
#else
# include "qaws_simd_scalar.h"
# define QAWS_SIMD_WIDTH_F32 1
# define QAWS_SIMD_WIDTH_F64 1
# define QAWS_SIMD_ISA "scalar"
#endif

#include "qaws_simd_types.h"

/* a buffer aligned for the widest register's aligned loads and stores */
#if defined(_MSC_VER)
# define QAWS_SIMD_ALIGN __declspec(align(32))
#elif defined(__GNUC__) || defined(__clang__)
# define QAWS_SIMD_ALIGN __attribute__((aligned(32)))
#else
# define QAWS_SIMD_ALIGN
#endif

#endif /* QAWS_SIMD_H */
