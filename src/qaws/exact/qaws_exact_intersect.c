#include "qaws_exact_curve.h"
#include "qaws_exact_roots.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

/* ------------------------------------------------------------------ */
/*  Hits                                                               */
/* ------------------------------------------------------------------ */

/* Trailing zero bits of a non-zero magnitude. */
static unsigned int ctz_int(qaws_exact_int const* x)
{
	unsigned int i, b;
	for (i = 0; i < (unsigned int)x->size; i++)
		if (x->limb[i])
		{
			for (b = 0; !((x->limb[i] >> b) & 1u); b++)
				;
			return 32 * i + b;
		}
	return 0;
}

/* t = (a 2^d + (b - a) index) 2^-(d + shift): nearest double, and whether it is exact. */
qaws_status qaws_exact_span_param_to_double(qaws_exact_span const* sp, int shift, uint64_t index, int depth, double* out, int* exact)
{
	qaws_exact_int num, t;
	qaws_status st;
	unsigned int bits;
	qaws_exact_int_from_i64(&num, sp->a);
	TRY(qaws_exact_int_shl(&num, &num, (unsigned int)depth));
	/* index < 2^62: split to stay in int64 */
	qaws_exact_int_from_i64(&t, (int64_t)(index >> 31));
	TRY(qaws_exact_int_mul_i64(&t, &t, sp->b - sp->a));
	TRY(qaws_exact_int_shl(&t, &t, 31));
	TRY(qaws_exact_int_add(&num, &num, &t));
	qaws_exact_int_from_i64(&t, (int64_t)(index & 0x7FFFFFFFu));
	TRY(qaws_exact_int_mul_i64(&t, &t, sp->b - sp->a));
	TRY(qaws_exact_int_add(&num, &num, &t));
	*out = ldexp(qaws_exact_int_to_double(&num), -(depth + shift));
	bits = qaws_exact_int_bits(&num);
	*exact = qaws_exact_int_is_zero(&num) || (bits - ctz_int(&num) <= 53 && fabs(*out) >= 2.2250738585072014e-308);
	return QAWS_STATUS_OK;
}
