#ifndef QAWS_EXACT_INT_H
#define QAWS_EXACT_INT_H

/*
 * Fixed-capacity signed integers for the exact kernel (C only).
 *
 * Sign-magnitude, little-endian 32-bit limbs (limb products fit in
 * uint64_t: portable C11, no compiler intrinsics), at most
 * QAWS_EXACT_MAX_BITS bits of magnitude. `size` counts the used limbs, so
 * small values stay cheap. Zero has size 0 and sign 0.
 *
 * Every operation that can grow a value is checked: a result that does not
 * fit returns QAWS_STATUS_EXACT_RANGE_EXCEEDED and leaves the destination
 * unspecified; nothing wraps. Destinations may alias operands.
 */

#include "../qaws_status.h"
#include <stdint.h>
#include <stddef.h>

#define QAWS_EXACT_MAX_BITS 2048
#define QAWS_EXACT_LIMBS (QAWS_EXACT_MAX_BITS / 32)

typedef struct qaws_exact_int
{
	uint32_t limb[QAWS_EXACT_LIMBS];
	int size;   /* used limbs; limb[size - 1] != 0 */
	int sign;   /* -1, 0, +1 */
} qaws_exact_int;

void qaws_exact_int_zero(qaws_exact_int* x);
void qaws_exact_int_from_i64(qaws_exact_int* x, int64_t v);

/* Decimal text with an optional sign. */
qaws_status qaws_exact_int_from_text(qaws_exact_int* x, char const* text);
/* Decimal text; BUFFER_TOO_SMALL when it does not fit (cap includes the 0). */
qaws_status qaws_exact_int_to_text(qaws_exact_int const* x, char* buf, size_t cap);

int qaws_exact_int_sign(qaws_exact_int const* x);
int qaws_exact_int_is_zero(qaws_exact_int const* x);
/* -1, 0, +1 */
int qaws_exact_int_cmp(qaws_exact_int const* a, qaws_exact_int const* b);
int qaws_exact_int_cmp_abs(qaws_exact_int const* a, qaws_exact_int const* b);
/* Bits of the magnitude (0 for zero). */
unsigned int qaws_exact_int_bits(qaws_exact_int const* x);

void qaws_exact_int_neg(qaws_exact_int* r, qaws_exact_int const* a);
qaws_status qaws_exact_int_add(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b);
qaws_status qaws_exact_int_sub(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b);
qaws_status qaws_exact_int_mul(qaws_exact_int* r, qaws_exact_int const* a, qaws_exact_int const* b);
qaws_status qaws_exact_int_mul_i64(qaws_exact_int* r, qaws_exact_int const* a, int64_t s);
/* a * 2^k */
qaws_status qaws_exact_int_shl(qaws_exact_int* r, qaws_exact_int const* a, unsigned int k);
/* sign(a) * floor(|a| / 2^k) */
void qaws_exact_int_shr(qaws_exact_int* r, qaws_exact_int const* a, unsigned int k);
/* a / d for a divisor d != 0 that divides a exactly; INTERNAL_ERROR otherwise. */
qaws_status qaws_exact_int_divexact_u32(qaws_exact_int* r, qaws_exact_int const* a, uint32_t d);

/* Nearest double (ties to even); +-HUGE_VAL past the double range. */
double qaws_exact_int_to_double(qaws_exact_int const* x);

/* Nearest double of num / den (den != 0), correctly rounded through a
   64-bit long division with a sticky remainder bit (subnormal results may
   round twice). */
double qaws_exact_ratio_to_double(qaws_exact_int const* num, qaws_exact_int const* den);

/*
 * A finite double is exactly m * 2^e with integer m: m odd (or zero), |m|
 * below 2^53. Returns 0 for infinities and NaN.
 */
int qaws_exact_split_double(double d, int64_t* m, int* e);

#endif /* QAWS_EXACT_INT_H */
