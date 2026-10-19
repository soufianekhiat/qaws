/*
 * The library's one allocation point (qaws_internal_memory.h) and the public allocator settings and arena (qaws_alloc.h).
 * The only calls to the C allocator in the library are here.
 *
 * A block is [header][user bytes]. The header ends with a tail (user size, kind, magic) right before the user bytes;
 * blocks from a custom allocator also hold a copy of its descriptor at the start of the header, so the block is freed by
 * those callbacks even after the caller's descriptor is gone. Header lengths are multiples of QAWS_MEMORY_ALIGNMENT
 * derived from sizeof, so the user bytes keep the block's alignment on 32- and 64-bit targets.
 */

#include "qaws_internal_memory.h"
#include "../qaws_alloc.h"
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
#include <malloc.h>
#include <intrin.h>
#define QAWS_THREAD_LOCAL __declspec(thread)
#else
#include <sched.h>
#define QAWS_THREAD_LOCAL _Thread_local
#endif

#define MEM_ALIGN ((size_t)QAWS_MEMORY_ALIGNMENT)
#define MEM_ROUND(n) ((((n) + MEM_ALIGN - 1) / MEM_ALIGN) * MEM_ALIGN)
#define MEM_MAGIC 0x51A75A11u

enum { MEM_SYSTEM = 1, MEM_LEGACY = 2, MEM_SIZED = 3 };

typedef struct mem_tail
{
	size_t size;                        /* user bytes */
	unsigned int kind;
	unsigned int magic;
} mem_tail;

/* header lengths: system blocks the tail only, custom blocks the descriptor copy and the tail */
#define MEM_HEAD_SYSTEM MEM_ROUND(sizeof(mem_tail))
#define MEM_HEAD_CUSTOM MEM_ROUND(sizeof(qaws_allocator) + sizeof(mem_tail))

/* ---- the C allocator, aligned to MEM_ALIGN on every target ---------------- */

static void* sys_alloc(size_t n)
{
#if defined(_MSC_VER)
	return _aligned_malloc(n, MEM_ALIGN);
#else
	if (_Alignof(max_align_t) >= MEM_ALIGN)
		return malloc(n);
	return aligned_alloc(MEM_ALIGN, MEM_ROUND(n));
#endif
}

static void sys_free(void* p)
{
#if defined(_MSC_VER)
	_aligned_free(p);
#else
	free(p);
#endif
}

/* NULL on failure, the old block kept */
static void* sys_realloc(void* p, size_t old_n, size_t n)
{
#if defined(_MSC_VER)
	(void)old_n;
	return _aligned_realloc(p, n, MEM_ALIGN);
#else
	void* q;
	if (_Alignof(max_align_t) >= MEM_ALIGN)
		return realloc(p, n);
	q = sys_alloc(n);
	if (q)
	{
		memcpy(q, p, old_n < n ? old_n : n);
		sys_free(p);
	}
	return q;
#endif
}

/* ---- atomics for the default's state ------------------------------------ */

#if defined(_MSC_VER)
static long mem_cas(long volatile* p, long expect, long want) { return _InterlockedCompareExchange(p, want, expect) == expect; }
static long mem_load(long volatile* p) { return _InterlockedOr(p, 0); }
static void mem_store(long volatile* p, long v) { _InterlockedExchange(p, v); }
#if defined(_M_ARM64) || defined(_M_ARM)
static void mem_yield(void) { __yield(); }
#else
static void mem_yield(void) { _mm_pause(); }
#endif
#else
static long mem_cas(long volatile* p, long expect, long want) { return __atomic_compare_exchange_n(p, &expect, want, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST); }
static long mem_load(long volatile* p) { return __atomic_load_n(p, __ATOMIC_SEQ_CST); }
static void mem_store(long volatile* p, long v) { __atomic_store_n(p, v, __ATOMIC_SEQ_CST); }
static void mem_yield(void) { sched_yield(); }
#endif

/* the default: settable until the library first reads it, then sealed */
enum { DEF_OPEN = 0, DEF_WRITING = 1, DEF_OPEN_SET = 2, DEF_SEALED = 3, DEF_SEALED_SET = 4 };
static long volatile g_default_state;
static qaws_allocator g_default;

static QAWS_THREAD_LOCAL qaws_allocator g_thread;
static QAWS_THREAD_LOCAL int g_thread_set;

void* qaws_allocator_extended(unsigned long size, void* user_data)
{
	(void)size;
	(void)user_data;
	return NULL;
}

static int mem_sized(qaws_allocator const* a)
{
	return a->alloc == qaws_allocator_extended;
}

static int mem_custom(qaws_allocator const* a)
{
	return a && a->alloc && (!mem_sized(a) || a->alloc_sized);
}

/* the members a descriptor defines: a legacy one is never read past its first three */
static void mem_copy(qaws_allocator* to, qaws_allocator const* from)
{
	memset(to, 0, sizeof *to);
	to->alloc = from->alloc;
	to->dealloc = from->dealloc;
	to->user_data = from->user_data;
	if (mem_sized(from))
	{
		to->alloc_sized = from->alloc_sized;
		to->dealloc_sized = from->dealloc_sized;
		to->resize = from->resize;
	}
}

qaws_status qaws_set_default_allocator(qaws_allocator const* allocator)
{
	for (;;)
	{
		long s = mem_load(&g_default_state);
		if (s == DEF_SEALED || s == DEF_SEALED_SET)
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
		if (s == DEF_WRITING)
		{
			mem_yield();
			continue;
		}
		if (!mem_cas(&g_default_state, s, DEF_WRITING))
			continue;
		if (mem_custom(allocator))
			mem_copy(&g_default, allocator);
		else
			memset(&g_default, 0, sizeof g_default);
		mem_store(&g_default_state, mem_custom(allocator) ? DEF_OPEN_SET : DEF_OPEN);
		return QAWS_STATUS_OK;
	}
}

/* the default as the library reads it, sealing it */
static qaws_allocator const* mem_default(void)
{
	for (;;)
	{
		long s = mem_load(&g_default_state);
		if (s == DEF_SEALED_SET)
			return &g_default;
		if (s == DEF_SEALED)
			return NULL;
		if (s == DEF_WRITING)
		{
			mem_yield();
			continue;
		}
		if (mem_cas(&g_default_state, s, s == DEF_OPEN_SET ? DEF_SEALED_SET : DEF_SEALED))
			return s == DEF_OPEN_SET ? &g_default : NULL;
	}
}

void qaws_set_thread_allocator(qaws_allocator const* allocator)
{
	g_thread_set = mem_custom(allocator);
	if (g_thread_set)
		mem_copy(&g_thread, allocator);
	else
		memset(&g_thread, 0, sizeof g_thread);
}

void qaws_get_thread_allocator(qaws_allocator* out_allocator)
{
	if (out_allocator)
		*out_allocator = g_thread;
}

void qaws_exchange_thread_allocator(qaws_allocator const* allocator, qaws_allocator* out_previous)
{
	qaws_allocator previous = g_thread;
	qaws_set_thread_allocator(allocator);
	if (out_previous)
		*out_previous = previous;
}

qaws_allocator const* qaws_internal_current_allocator(void)
{
	return g_thread_set ? &g_thread : mem_default();
}

/* ---- blocks ------------------------------------------------------------- */

static void mem_read_tail(void const* ptr, mem_tail* t)
{
	memcpy(t, (unsigned char const*)ptr - sizeof(mem_tail), sizeof *t);
}

/* the raw block of a user pointer; 0 when it is not one of ours */
static int mem_block(void const* ptr, mem_tail* t, unsigned char** raw, qaws_allocator* a)
{
	mem_read_tail(ptr, t);
	assert(t->magic == MEM_MAGIC && "qaws: a pointer the library did not allocate");
	if (t->magic != MEM_MAGIC || t->kind < MEM_SYSTEM || t->kind > MEM_SIZED)
		return 0;
	*raw = (unsigned char*)ptr - (t->kind == MEM_SYSTEM ? MEM_HEAD_SYSTEM : MEM_HEAD_CUSTOM);
	if (t->kind != MEM_SYSTEM)
		memcpy(a, *raw, sizeof *a);
	return 1;
}

static void mem_release(unsigned char* raw, unsigned int kind, qaws_allocator const* a, size_t total)
{
	if (kind == MEM_SYSTEM)
		sys_free(raw);
	else if (kind == MEM_SIZED)
	{
		if (a->dealloc_sized)
			a->dealloc_sized(raw, total, MEM_ALIGN, a->user_data);
	}
	else if (a->dealloc)
		a->dealloc(raw, a->user_data);
}

static void* mem_finish(unsigned char* raw, unsigned int kind, qaws_allocator const* a, size_t size)
{
	size_t head = kind == MEM_SYSTEM ? MEM_HEAD_SYSTEM : MEM_HEAD_CUSTOM;
	mem_tail t;
	if (kind != MEM_SYSTEM)
	{
		qaws_allocator c;
		mem_copy(&c, a);
		memcpy(raw, &c, sizeof c);
	}
	t.size = size;
	t.kind = kind;
	t.magic = MEM_MAGIC;
	memcpy(raw + head - sizeof t, &t, sizeof t);
	return raw + head;
}

void* qaws_internal_alloc_from(qaws_allocator const* allocator, size_t size)
{
	unsigned char* raw;
	unsigned int kind;
	size_t head, total;
	if (!mem_custom(allocator))
		allocator = qaws_internal_current_allocator();
	kind = !mem_custom(allocator) ? MEM_SYSTEM : mem_sized(allocator) ? MEM_SIZED : MEM_LEGACY;
	head = kind == MEM_SYSTEM ? MEM_HEAD_SYSTEM : MEM_HEAD_CUSTOM;
	if (size > SIZE_MAX - head - MEM_ALIGN)
		return NULL;
	total = size + head;
	if (kind == MEM_SYSTEM)
		raw = (unsigned char*)sys_alloc(total);
	else if (kind == MEM_SIZED)
		raw = (unsigned char*)allocator->alloc_sized(total, MEM_ALIGN, allocator->user_data);
	else
	{
		/* the legacy interface counts in unsigned long: larger requests fail, never truncated */
		if (total > (size_t)ULONG_MAX)
			return NULL;
		raw = (unsigned char*)allocator->alloc((unsigned long)total, allocator->user_data);
	}
	if (!raw)
		return NULL;
	if ((uintptr_t)raw % MEM_ALIGN)
	{
		/* a misaligned block is given back: the call fails rather than misalign the library's data */
		assert(!"qaws: the allocator returned a block not aligned to QAWS_MEMORY_ALIGNMENT");
		mem_release(raw, kind, allocator, total);
		return NULL;
	}
	return mem_finish(raw, kind, allocator, size);
}

void* qaws_internal_malloc(size_t size)
{
	return qaws_internal_alloc_from(NULL, size);
}

void* qaws_internal_zalloc(size_t count, size_t size)
{
	void* p;
	if (size && count > SIZE_MAX / size)
		return NULL;
	p = qaws_internal_alloc_from(NULL, count * size);
	if (p)
		memset(p, 0, count * size);
	return p;
}

void qaws_internal_free(void* ptr)
{
	mem_tail t;
	unsigned char* raw;
	qaws_allocator a;
	if (!ptr || !mem_block(ptr, &t, &raw, &a))
		return;
	mem_release(raw, t.kind, &a, t.size + (t.kind == MEM_SYSTEM ? MEM_HEAD_SYSTEM : MEM_HEAD_CUSTOM));
}

void* qaws_internal_realloc(void* ptr, size_t size)
{
	mem_tail t;
	unsigned char* raw;
	qaws_allocator a;
	size_t head;
	void* p;
	if (!ptr)
		return qaws_internal_malloc(size);
	if (size == 0)
	{
		qaws_internal_free(ptr);
		return NULL;
	}
	if (!mem_block(ptr, &t, &raw, &a))
		return NULL;
	head = t.kind == MEM_SYSTEM ? MEM_HEAD_SYSTEM : MEM_HEAD_CUSTOM;
	if (size > SIZE_MAX - head - MEM_ALIGN)
		return NULL;
	if (t.kind == MEM_SYSTEM)
	{
		unsigned char* q = (unsigned char*)sys_realloc(raw, t.size + head, size + head);
		return q ? mem_finish(q, MEM_SYSTEM, NULL, size) : NULL;
	}
	if (t.kind == MEM_SIZED && a.resize && a.resize(raw, t.size + head, size + head, MEM_ALIGN, a.user_data))
		return mem_finish(raw, MEM_SIZED, &a, size);
	/* a new block from the same allocator; the old one kept on failure */
	p = qaws_internal_alloc_from(&a, size);
	if (!p)
		return NULL;
	memcpy(p, ptr, t.size < size ? t.size : size);
	qaws_internal_free(ptr);
	return p;
}

size_t qaws_internal_block_size(void const* ptr)
{
	mem_tail t;
	if (!ptr)
		return 0;
	mem_read_tail(ptr, &t);
	return t.magic == MEM_MAGIC ? t.size : 0;
}

/* ---- arena over caller memory ------------------------------------------- */

/* Each block has a small header before it: the offset of the block allocated before it and its size with a freed flag.
   A freed block is marked; the freed blocks at the top are then reclaimed, so blocks released in any last-in-first-out
   order (nested scratch) give their space back at once, and the rest at qaws_arena_reset. */
typedef struct arena_head
{
	size_t prev;                        /* the previous block's header offset, ARENA_NONE for the first */
	size_t size;                        /* the block's bytes */
	size_t freed;
} arena_head;

#define ARENA_NONE ((size_t)-1)
#define ARENA_HEAD MEM_ROUND(sizeof(arena_head))

static arena_head arena_read(qaws_arena const* a, size_t at)
{
	arena_head h;
	memcpy(&h, a->base + at, sizeof h);
	return h;
}

static void arena_write(qaws_arena* a, size_t at, arena_head const* h)
{
	memcpy(a->base + at, h, sizeof *h);
}

/* a block after an aligned header; NULL past the end (no arithmetic overflows) */
static void* arena_alloc(size_t size, size_t alignment, void* user)
{
	qaws_arena* a = (qaws_arena*)user;
	size_t align = alignment > MEM_ALIGN ? alignment : MEM_ALIGN;
	size_t pad = (align - a->used % align) % align;
	arena_head h;
	if (pad > a->size - a->used || ARENA_HEAD > a->size - a->used - pad || size > a->size - a->used - pad - ARENA_HEAD)
	{
		a->failures++;
		return NULL;
	}
	h.prev = a->last;
	h.size = size;
	h.freed = 0;
	a->last = a->used + pad;
	arena_write(a, a->last, &h);
	a->used = a->last + ARENA_HEAD + size;
	a->live++;
	if (a->used > a->peak)
		a->peak = a->used;
	return a->base + a->last + ARENA_HEAD;
}

/* the header offset of a block it holds, ARENA_NONE otherwise */
static size_t arena_find(qaws_arena const* a, void const* ptr)
{
	unsigned char const* p = (unsigned char const*)ptr;
	size_t at;
	if (!p || p < a->base + ARENA_HEAD || p >= a->base + a->used)
		return ARENA_NONE;
	at = (size_t)(p - a->base) - ARENA_HEAD;
	return at % MEM_ALIGN == 0 ? at : ARENA_NONE;
}

static void arena_dealloc(void* ptr, size_t size, size_t alignment, void* user)
{
	qaws_arena* a = (qaws_arena*)user;
	size_t at = arena_find(a, ptr);
	arena_head h;
	(void)size;
	(void)alignment;
	if (at == ARENA_NONE || (h = arena_read(a, at)).freed || a->live == 0)
	{
		/* a block from another generation, another arena, or freed twice */
		assert(!"qaws_arena: a block it does not hold");
		a->invalid_frees++;
		return;
	}
	h.freed = 1;
	arena_write(a, at, &h);
	a->live--;
	/* the freed blocks at the top give their space back */
	while (a->last != ARENA_NONE && (h = arena_read(a, a->last)).freed)
	{
		a->used = a->last;
		a->last = h.prev;
	}
}

/* the latest block grows or shrinks in place */
static int arena_resize(void* ptr, size_t old_size, size_t new_size, size_t alignment, void* user)
{
	qaws_arena* a = (qaws_arena*)user;
	size_t at = arena_find(a, ptr);
	arena_head h;
	(void)old_size;
	(void)alignment;
	if (at == ARENA_NONE || at != a->last || new_size > a->size - at - ARENA_HEAD)
		return 0;
	h = arena_read(a, at);
	h.size = new_size;
	arena_write(a, at, &h);
	a->used = at + ARENA_HEAD + new_size;
	if (a->used > a->peak)
		a->peak = a->used;
	return 1;
}

void qaws_arena_init(qaws_arena* arena, void* memory, size_t size)
{
	size_t skip = (size_t)((MEM_ALIGN - ((uintptr_t)memory % MEM_ALIGN)) % MEM_ALIGN);
	memset(arena, 0, sizeof *arena);
	arena->base = (unsigned char*)memory + (skip < size ? skip : size);
	arena->size = skip < size ? size - skip : 0;
	arena->last = ARENA_NONE;
	arena->allocator.alloc = qaws_allocator_extended;
	arena->allocator.alloc_sized = arena_alloc;
	arena->allocator.dealloc_sized = arena_dealloc;
	arena->allocator.resize = arena_resize;
	arena->allocator.user_data = arena;
}

qaws_status qaws_arena_reset(qaws_arena* arena)
{
	/* blocks still held (live objects, results) would be overwritten */
	if (arena->live)
		return QAWS_STATUS_INVALID_ARGUMENT;
	arena->used = 0;
	arena->last = ARENA_NONE;
	arena->generation++;
	return QAWS_STATUS_OK;
}
