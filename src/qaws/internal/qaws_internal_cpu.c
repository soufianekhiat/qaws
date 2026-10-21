/*
 * The processor's instruction sets (qaws_internal_cpu.h).
 */

#include "qaws_internal_cpu.h"

/* (precision-free, one state: defined by the f64 pass only, like the memory module) */
#if !QAWS_SCALAR_IS_FLOAT

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <intrin.h>
#include <immintrin.h>
#define CPU_X86 1
static void cpu_id(unsigned int leaf, unsigned int sub, unsigned int* r)
{
	int v[4];
	__cpuidex(v, (int)leaf, (int)sub);
	r[0] = (unsigned int)v[0];
	r[1] = (unsigned int)v[1];
	r[2] = (unsigned int)v[2];
	r[3] = (unsigned int)v[3];
}
static unsigned long long cpu_xcr0(void)
{
	return _xgetbv(0);
}
#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <cpuid.h>
#define CPU_X86 1
static void cpu_id(unsigned int leaf, unsigned int sub, unsigned int* r)
{
	if (!__get_cpuid_count(leaf, sub, r, r + 1, r + 2, r + 3))
		r[0] = r[1] = r[2] = r[3] = 0;
}
static unsigned long long cpu_xcr0(void)
{
	unsigned int lo, hi;
	__asm__ volatile("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
	return ((unsigned long long)hi << 32) | lo;
}
#else
#define CPU_X86 0
#endif

static unsigned int cpu_mask = ~0u;
static int cpu_done;
static unsigned int cpu_bits;

static unsigned int cpu_probe(void)
{
	unsigned int f = 0;
#if CPU_X86
	unsigned int r[4], max, ecx1;
	cpu_id(0, 0, r);
	max = r[0];
	if (max < 1)
		return 0;
	cpu_id(1, 0, r);
	ecx1 = r[2];
	if (ecx1 & (1u << 19))
		f |= QAWS_CPU_SSE41;
	/* the wide registers only when the operating system saves them (OSXSAVE, then XCR0: SSE and YMM state, and the 3 AVX-512 states) */
	if ((ecx1 & (1u << 27)) && (ecx1 & (1u << 28)))
	{
		unsigned long long x = cpu_xcr0();
		if ((x & 6) == 6)
		{
			f |= QAWS_CPU_AVX;
			if (ecx1 & (1u << 12))
				f |= QAWS_CPU_FMA;
			if (max >= 7)
			{
				cpu_id(7, 0, r);
				if (r[1] & (1u << 5))
					f |= QAWS_CPU_AVX2;
				if ((r[1] & (1u << 16)) && (x & 0xe6) == 0xe6)
					f |= QAWS_CPU_AVX512F;
			}
		}
	}
#endif
	return f;
}

unsigned int qaws_internal_cpu_features(void)
{
	/* (a race only writes the same bits twice) */
	if (!cpu_done)
	{
		cpu_bits = cpu_probe();
		cpu_done = 1;
	}
	return cpu_bits & cpu_mask;
}

void qaws_internal_cpu_override(unsigned int mask)
{
	cpu_mask = mask;
}

#else
typedef int qaws_internal_cpu_f32_unit;   /* (no empty unit in the f32 pass) */
#endif /* !QAWS_SCALAR_IS_FLOAT */
