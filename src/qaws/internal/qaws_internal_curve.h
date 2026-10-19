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

/* Allocator-aware helpers for impl allocations. */
void* qaws_internal_alloc(qaws_allocator const* allocator, size_t size);
/* qaws_internal_mul2..4 for allocation sizes */
#include "qaws_internal_memory.h"
void  qaws_internal_dealloc(qaws_allocator const* allocator, void* ptr);
/* count * size zeroed bytes (NULL on overflow or failure) */
void* qaws_internal_calloc(qaws_allocator const* allocator, size_t count, size_t size);

/* Evaluation scratch: `stack` (stack_bytes, the caller's array) when bytes fit, else bytes from the allocator (an object's own, so
   evaluation never calls malloc); released by qaws_internal_scratch_release with the same stack pointer (NULL is ignored). */
#define QAWS_INTERNAL_SCRATCH 32
void* qaws_internal_scratch(qaws_allocator const* allocator, void* stack, size_t stack_bytes, size_t bytes);
void qaws_internal_scratch_release(qaws_allocator const* allocator, void* stack, void* scratch);

#endif /* QAWS_INTERNAL_CURVE_H */
