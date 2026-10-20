#ifndef QAWS_INTERNAL_MEMORY_H
#define QAWS_INTERNAL_MEMORY_H

#include "../qaws_types.h"
#include <stddef.h>

/*
 * Every allocation of the library goes through these. A block carries a header with its size and, from a custom
 * allocator, a copy of that allocator's descriptor, so it is freed (and grown) by it whatever the settings or the
 * caller's descriptor when it is released.
 *
 * Without an explicit allocator, memory comes from the calling thread's allocator (qaws_set_thread_allocator), else the
 * process default (qaws_set_default_allocator), else malloc.
 */

/* the allocator a block without an explicit one comes from now (NULL: malloc) */
qaws_allocator const* qaws_internal_current_allocator(void);   /* a copy the library holds */

void* qaws_internal_alloc_from(size_t size, qaws_allocator const* allocator);   /* NULL: the current allocator */
void* qaws_internal_malloc(size_t size);
void* qaws_internal_zalloc(size_t count, size_t size);
void* qaws_internal_realloc(void* ptr, size_t size);                           /* by the block's own allocator */
void qaws_internal_free(void* ptr);

/* a temporary of one call (freed before it returns, never part of a result): from the thread's scratch allocator when
   one is set, else as qaws_internal_malloc */
void* qaws_internal_scratch_malloc(size_t size);
void* qaws_internal_scratch_zalloc(size_t count, size_t size);

/* checked size products: SIZE_MAX when the product does not fit, which every allocation refuses (a size that wraps
   would allocate a short block the caller then overruns) */
size_t qaws_internal_mul2(size_t a, size_t b);
size_t qaws_internal_mul3(size_t a, size_t b, size_t c);
size_t qaws_internal_mul4(size_t a, size_t b, size_t c, size_t d);
size_t qaws_internal_add2(size_t a, size_t b);   /* SIZE_MAX when the sum does not fit */
/* the next capacity of a growing list: first when cap is 0, else cap doubled until it reaches need (UINT_MAX when doubling would pass it); 0 when need cannot be held (need > UINT_MAX) */
unsigned int qaws_internal_grow(unsigned int cap, size_t need, unsigned int first);

/* the user size of a block (0 for NULL) */
size_t qaws_internal_block_size(void const* ptr);

#endif /* QAWS_INTERNAL_MEMORY_H */
