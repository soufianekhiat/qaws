#ifndef QAWS_INTERNAL_CURVE_H
#define QAWS_INTERNAL_CURVE_H

#include "qaws_internal_types.h"

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

/* Allocator-aware helpers for impl allocations. */
void* qaws_internal_alloc(unsigned long size, qaws_allocator const* allocator);
void  qaws_internal_dealloc(void* ptr, qaws_allocator const* allocator);

#endif /* QAWS_INTERNAL_CURVE_H */
