#ifndef QAWS_INTERNAL_KINDS_H
#define QAWS_INTERNAL_KINDS_H

/* Impl structs of the curve kinds whose public headers carry their own types. */

#include "qaws_internal_types.h"
#include "../qaws_arc.h"
#include "../qaws_subdivision.h"

typedef struct qaws_arc_impl
{
	qaws_arc_segment* segments;
	unsigned int segment_count;
} qaws_arc_impl;

typedef struct qaws_polynomial_impl
{
	qaws_scalar* coefficients;         /* (degree+1) * dim_count scalars */
	unsigned int coefficient_count;
} qaws_polynomial_impl;

typedef struct qaws_clothoid_impl
{
	qaws_scalar origin_x, origin_y;
	qaws_scalar start_angle;
	qaws_scalar kappa_0;             /* start curvature */
	qaws_scalar kappa_1;             /* end curvature */
	qaws_scalar length;              /* total arc length */
	qaws_scalar origin[2];           /* origin, contiguous (CENTER field) */
	qaws_scalar rate;                /* (kappa_1 - kappa_0) / length (CURVATURE_RATE field) */
} qaws_clothoid_impl;

typedef struct qaws_subdivision_impl
{
	qaws_scalar* refined_points;      /* refined_count * dim_count */
	unsigned int refined_count;
	int closed;
	qaws_scalar* segment_coeffs;      /* uniform Catmull-Rom cubic coefficients */
	qaws_subdivision_scheme scheme;
} qaws_subdivision_impl;

/* Arc-length wrappers (QAWS_CURVE_KIND_REPARAMETERIZED): the source curve and
   its parameter at wrapper parameter t. Returns 0 when curve is not one. */
int qaws_internal_reparam_source(
	qaws_curve const* curve,
	qaws_scalar t,
	qaws_curve const** out_source,
	qaws_scalar* out_source_t);

/* 1 when the curve's two ends coincide, relative to the size of the curve. */
int qaws_internal_curve_ends_meet(qaws_curve const* curve);

#endif /* QAWS_INTERNAL_KINDS_H */
