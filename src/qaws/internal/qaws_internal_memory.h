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

void* qaws_internal_alloc_from(qaws_allocator const* allocator, size_t size);   /* NULL: the current allocator */
void* qaws_internal_malloc(size_t size);
void* qaws_internal_zalloc(size_t count, size_t size);
void* qaws_internal_realloc(void* ptr, size_t size);                           /* by the block's own allocator */
void qaws_internal_free(void* ptr);

/* the user size of a block (0 for NULL) */
size_t qaws_internal_block_size(void const* ptr);

#endif /* QAWS_INTERNAL_MEMORY_H */
