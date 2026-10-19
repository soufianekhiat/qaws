#ifndef QAWS_INTERNAL_MEMORY_H
#define QAWS_INTERNAL_MEMORY_H

#include "../qaws_types.h"
#include <stddef.h>

/*
 * Every allocation of the library goes through these. A block carries a 16-byte header: the allocator that gave it and
 * its size, so it is freed (and grown) by that allocator whatever the settings when it is released.
 *
 * Without an explicit allocator, memory comes from the calling thread's allocator (qaws_set_thread_allocator), else the
 * process default (qaws_set_default_allocator), else malloc.
 */

/* the allocator a block without an explicit one comes from now (NULL: malloc) */
qaws_allocator const* qaws_internal_current_allocator(void);

void* qaws_internal_alloc_from(qaws_allocator const* allocator, size_t size);   /* NULL: the current allocator */
void* qaws_internal_malloc(size_t size);
void* qaws_internal_zalloc(size_t count, size_t size);
void* qaws_internal_realloc(void* ptr, size_t size);                           /* by the block's own allocator */
void qaws_internal_free(void* ptr);

#endif /* QAWS_INTERNAL_MEMORY_H */
