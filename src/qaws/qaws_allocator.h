#ifndef QAWS_ALLOCATOR_H
#define QAWS_ALLOCATOR_H

/*
 * The allocator, its settings and the arena. Precision-free: declared once
 * and defined once in the library (internal/qaws_internal_memory.c), so the
 * process default, the thread allocators and the arenas are the same for the
 * f32 and the f64 functions.
 */

#include "qaws_status.h"
#include <stddef.h>

/*
 * Custom allocator. Pass to _ex creation functions or install it (qaws_alloc.h). If allocation returns NULL, the call
 * fails with QAWS_STATUS_ALLOCATION_FAILURE.
 *
 * The library copies this descriptor (callbacks and user_data) into every block it allocates and into the objects that
 * allocate later, so the descriptor itself may be a short-lived local; the state user_data points to must outlive every
 * block the allocator gave.
 *
 * Two interfaces. A legacy descriptor sets the first three members only; the others are never read. An extended one sets
 * alloc to qaws_allocator_extended (the marker) and the sized members:
 *   legacy: alloc(size, user_data), dealloc(ptr, user_data); sizes in unsigned long (32 bits on Windows: larger
 *           requests fail);
 *   sized:  alloc_sized(size, alignment, user_data), dealloc_sized(ptr, size, alignment, user_data) with the size and
 *           alignment of the allocation; resize (optional) grows or shrinks a block in place and returns non-zero, or
 *           returns 0 to let the library allocate, copy and free.
 * Blocks must be aligned to QAWS_MEMORY_ALIGNMENT (16); misaligned blocks are given back and the call fails.
 */
typedef struct qaws_allocator
{
	void* (*alloc)(unsigned long size, void* user_data);
	void  (*dealloc)(void* ptr, void* user_data);
	void* user_data;
	void* (*alloc_sized)(size_t size, size_t alignment, void* user_data);
	void  (*dealloc_sized)(void* ptr, size_t size, size_t alignment, void* user_data);
	int   (*resize)(void* ptr, size_t old_size, size_t new_size, size_t alignment, void* user_data);
} qaws_allocator;

#define QAWS_MEMORY_ALIGNMENT 16

/* the marker of an extended descriptor (alloc = qaws_allocator_extended); never allocates when called */
void* qaws_allocator_extended(unsigned long size, void* user_data);

/*
 * Where the library's memory comes from. Without an explicit allocator (a NULL argument to the _ex functions, and all
 * internal scratch and results), memory comes from the calling thread's allocator if set, else the process default if
 * set, else the C allocator (aligned to QAWS_MEMORY_ALIGNMENT).
 *
 * Lifetimes: every setting and every _ex argument is copied (callbacks and user_data). A descriptor may be a local that
 * goes out of scope; each block keeps the copy it was made with and is freed and grown by it, whatever the settings at
 * that time and on whichever thread. What user_data points to must outlive every block of that allocator.
 *
 * The default is write-once in effect: it can be set (or cleared with NULL) until the library first allocates from it;
 * after that the call returns QAWS_STATUS_UNSUPPORTED_OPERATION and nothing changes. Reading it then is a plain load,
 * so allocations take no lock.
 *
 * A thread's allocator applies to that thread only. qaws_exchange_thread_allocator installs one and returns the one it
 * replaces (zeroed: none) so a host can nest settings and restore them, on every return path, with
 * qaws_set_thread_allocator(&previous). The library never changes it.
 */
qaws_status qaws_set_default_allocator(qaws_allocator const* allocator);   /* NULL: the C allocator */
void qaws_set_thread_allocator(qaws_allocator const* allocator);            /* this thread; NULL: the default */
void qaws_get_thread_allocator(qaws_allocator* out_allocator);              /* zeroed: none */
void qaws_exchange_thread_allocator(qaws_allocator* out_previous, qaws_allocator const* allocator);

/*
 * Scratch apart from results: with a scratch allocator installed on a thread, the temporaries of adaptive sampling,
 * tessellation, meshing, clipping and the repulsive solvers (blocks allocated and released within one call) come from
 * it; results and objects still come from the _ex argument, the thread allocator or the default as above, so a frame
 * arena can serve as scratch and be reset after each call while results stay. Other operations take their scratch as
 * before. Installs it (NULL: none) and returns the one it replaces (zeroed: none); restore it on every return path.
 */
void qaws_exchange_thread_scratch_allocator(qaws_allocator* out_previous, qaws_allocator const* allocator);

/*
 * An allocator over caller memory (sized interface): a mapped file, pinned memory, a host frame arena. The memory must
 * be CPU-accessible and coherent for the CPU; placing the library's data there does not move any computation to a GPU.
 *
 * Blocks are cut in order, aligned to QAWS_MEMORY_ALIGNMENT, each after a small header. The latest block grows and shrinks
 * in place. A freed block is marked, and the freed blocks at the top give their space back, so blocks released in any
 * last-in-first-out order (nested scratch) are reclaimed at once; others wait for the blocks above them, or the reset. It is a hard bound: when full, allocation fails and
 * calls return QAWS_STATUS_ALLOCATION_FAILURE; nothing falls back to the heap. One thread at a time (one arena per
 * worker with qaws_set_thread_allocator). Pass &arena.allocator wherever an allocator is taken.
 *
 * Lifetime: objects and results allocated in an arena live in it. qaws_arena_reset refuses (QAWS_STATUS_INVALID_ARGUMENT)
 * while blocks are still held (live > 0): destroy every arena-backed object and finish every call first. Freeing a block
 * the arena does not hold (another arena, an earlier generation, twice) is counted in invalid_frees and asserts in
 * debug builds.
 */
typedef struct qaws_arena
{
	qaws_allocator allocator;
	unsigned char* base;
	size_t size;
	size_t used;                        /* bytes up to the end of the latest block */
	size_t peak;                        /* most bytes used since init */
	size_t last;                        /* the latest block's header offset ((size_t)-1: none) */
	size_t live;                        /* blocks held */
	size_t failures;                    /* allocations refused for room */
	size_t invalid_frees;
	unsigned int generation;            /* resets so far */
} qaws_arena;

void qaws_arena_init(qaws_arena* arena, void* memory, size_t size);
qaws_status qaws_arena_reset(qaws_arena* arena);

#endif /* QAWS_ALLOCATOR_H */
