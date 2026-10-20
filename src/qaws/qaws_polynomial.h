#ifndef QAWS_POLYNOMIAL_H
#define QAWS_POLYNOMIAL_H

#include "qaws_types.h"
#include "qaws_status.h"

/* Polynomial curve in monomial form: C(t) = sum_{i=0}^{degree} coefficients[i] * t^i
 * Each coefficient is a point (2 or 3 scalars depending on dimension).
 * coefficients layout: [c0_x, c0_y, [c0_z], c1_x, c1_y, [c1_z], ...]
 * So the curve is: P(t) = c0 + c1*t + c2*t^2 + ... + cn*t^n */
typedef struct qaws_polynomial_desc
{
	qaws_dimension dimension;
	unsigned int degree;                /* polynomial degree (>= 1) */
	void const* coefficients;           /* (degree+1) * dim_count scalars */
	unsigned int coefficient_count;     /* must be degree+1 */
	qaws_scalar t_min;                  /* parameter range min (default 0) */
	qaws_scalar t_max;                  /* parameter range max (default 1) */
} qaws_polynomial_desc;

qaws_status qaws_curve_create_polynomial(
	qaws_curve** out_curve,
	qaws_polynomial_desc const* desc);

#endif /* QAWS_POLYNOMIAL_H */
