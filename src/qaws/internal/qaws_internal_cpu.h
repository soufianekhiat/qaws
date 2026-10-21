#ifndef QAWS_INTERNAL_CPU_H
#define QAWS_INTERNAL_CPU_H

/*
 * The processor's instruction sets, probed once (CPUID, and XGETBV for the
 * operating system's saving of the YMM and ZMM registers), for the kernels
 * that pick a path at run time. Each wide kernel is compiled for its target
 * per function, never the whole library.
 */

#define QAWS_CPU_SSE41 1u
#define QAWS_CPU_AVX 2u
#define QAWS_CPU_AVX2 4u
#define QAWS_CPU_FMA 8u
#define QAWS_CPU_AVX512F 16u

/* the QAWS_CPU_* bits the processor and the operating system both support, less those masked off by the override */
unsigned int qaws_internal_cpu_features(void);

/* keeps only the features in mask from now on (~0u: all again); for tests comparing paths, not thread-safe */
void qaws_internal_cpu_override(unsigned int mask);

#endif /* QAWS_INTERNAL_CPU_H */
