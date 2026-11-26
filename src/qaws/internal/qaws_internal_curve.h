#ifndef QAWS_INTERNAL_CURVE_H
#define QAWS_INTERNAL_CURVE_H

#include "qaws_internal_types.h"
#include <stddef.h>

qaws_curve* qaws_internal_curve_alloc(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable);

qaws_curve* qaws_internal_curve_alloc_ex(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable,
	qaws_allocator const* allocator);

void qaws_internal_curve_free(qaws_curve* curve);

/* the allocator-aware helpers (qaws_internal_alloc, _dealloc, _calloc, _scratch): qaws_internal_memory.h, precision-free */
#include "qaws_internal_memory.h"

/* qaws_curve_evaluate_batch_2d/3d on a B-spline of degree 1 to 5: the basis functions of a SIMD register's worth of
   parameters at once (qaws_internal_bspline_basis_derivs_lanes), each result as the one-point evaluation makes it, bit for
   bit; QAWS_STATUS_UNSUPPORTED_OPERATION for the other curves and degrees and without SIMD (the caller loops). */
qaws_status qaws_internal_bspline_eval_batch(
	void* out_results,
	qaws_curve const* curve,
	qaws_scalar const* parameters,
	unsigned int count,
	unsigned int eval_flags);

#endif /* QAWS_INTERNAL_CURVE_H */
