#include "qaws_exact_curve.h"
#include "qaws_exact_roots.h"
#include "qaws_exact_poly.h"
#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define CC_MAX_IMPLICIT 6
#define CC_COEF QAWS_EXACT_POLY_COEF
#define CC_MAX_ROOTS QAWS_EXACT_SOLVE_MAX_ROOTS

void qaws_exact_span_box(qaws_exact_span const* sp, unsigned int dim, double lo[3], double hi[3])
{
	unsigned int D = dim + 1, j, k;
	for (k = 0; k < 3; k++)
	{
		lo[k] = k < dim ? HUGE_VAL : 0;
		hi[k] = k < dim ? -HUGE_VAL : 0;
	}
	for (j = 0; j <= sp->degree; j++)
		if (qaws_exact_int_sign(&sp->h[j * D + dim]) <= 0)
		{
			for (k = 0; k < dim; k++)
			{
				lo[k] = HUGE_VAL;
				hi[k] = -HUGE_VAL;
			}
			return;
		}
	for (j = 0; j <= sp->degree; j++)
		for (k = 0; k < dim; k++)
		{
			/* correctly rounded (twice at worst, for subnormals): two ulps out are sound */
			double v = qaws_exact_ratio_to_double(&sp->h[j * D + k], &sp->h[j * D + dim]);
			double l = nextafter(nextafter(v, -HUGE_VAL), -HUGE_VAL), h = nextafter(nextafter(v, HUGE_VAL), HUGE_VAL);
			if (l < lo[k]) lo[k] = l;
			if (h > hi[k]) hi[k] = h;
		}
}
