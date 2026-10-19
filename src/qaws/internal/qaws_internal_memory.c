/*
 * The library's one allocation point (qaws_internal_memory.h) and the public allocator settings and arena (qaws_alloc.h).
 * The only calls to malloc, realloc and free in the library are here.
 */

#include "qaws_internal_memory.h"
#include "../qaws_alloc.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
#define QAWS_THREAD_LOCAL __declspec(thread)
#else
#define QAWS_THREAD_LOCAL _Thread_local
#endif

/* the header in front of every block: its allocator (NULL: malloc) and its size; 16 bytes keep 16-byte alignment */
#define QAWS_BLOCK_HEADER 16

static qaws_allocator const* g_default_allocator;
static QAWS_THREAD_LOCAL qaws_allocator const* g_thread_allocator;

void qaws_set_default_allocator(qaws_allocator const* allocator)
{
	g_default_allocator = allocator;
}

void qaws_set_thread_allocator(qaws_allocator const* allocator)
{
	g_thread_allocator = allocator;
}

qaws_allocator const* qaws_internal_current_allocator(void)
{
	return g_thread_allocator ? g_thread_allocator : g_default_allocator;
}

static void mem_header(unsigned char* raw, qaws_allocator const* allocator, size_t size)
{
	memcpy(raw, &allocator, sizeof allocator);
	memcpy(raw + 8, &size, sizeof size);
}

static qaws_allocator const* mem_owner(void const* ptr, size_t* size)
{
	qaws_allocator const* allocator;
	unsigned char const* raw = (unsigned char const*)ptr - QAWS_BLOCK_HEADER;
	memcpy(&allocator, raw, sizeof allocator);
	if (size)
		memcpy(size, raw + 8, sizeof *size);
	return allocator;
}

void* qaws_internal_alloc_from(qaws_allocator const* allocator, size_t size)
{
	unsigned char* raw;
	if (!allocator)
		allocator = qaws_internal_current_allocator();
	if (allocator && !allocator->alloc)
		allocator = NULL;
	if (size > (size_t)-1 - QAWS_BLOCK_HEADER)
		return NULL;
	if (allocator)
	{
		/* the allocator interface counts in unsigned long */
		if (size + QAWS_BLOCK_HEADER > (size_t)ULONG_MAX)
			return NULL;
		raw = (unsigned char*)allocator->alloc((unsigned long)(size + QAWS_BLOCK_HEADER), allocator->user_data);
	}
	else
		raw = (unsigned char*)malloc(size + QAWS_BLOCK_HEADER);
	if (!raw)
		return NULL;
	mem_header(raw, allocator, size);
	return raw + QAWS_BLOCK_HEADER;
}

void* qaws_internal_malloc(size_t size)
{
	return qaws_internal_alloc_from(NULL, size);
}

void* qaws_internal_zalloc(size_t count, size_t size)
{
	void* p;
	if (size && count > (size_t)-1 / size)
		return NULL;
	p = qaws_internal_alloc_from(NULL, count * size);
	if (p)
		memset(p, 0, count * size);
	return p;
}

void qaws_internal_free(void* ptr)
{
	qaws_allocator const* allocator;
	unsigned char* raw;
	if (!ptr)
		return;
	allocator = mem_owner(ptr, NULL);
	raw = (unsigned char*)ptr - QAWS_BLOCK_HEADER;
	if (allocator)
	{
		if (allocator->dealloc)
			allocator->dealloc(raw, allocator->user_data);
	}
	else
		free(raw);
}

void* qaws_internal_realloc(void* ptr, size_t size)
{
	qaws_allocator const* allocator;
	size_t old;
	void* p;
	if (!ptr)
		return qaws_internal_malloc(size);
	if (size == 0)
	{
		qaws_internal_free(ptr);
		return NULL;
	}
	allocator = mem_owner(ptr, &old);
	if (!allocator)
	{
		unsigned char* raw;
		if (size > (size_t)-1 - QAWS_BLOCK_HEADER)
			return NULL;
		raw = (unsigned char*)realloc((unsigned char*)ptr - QAWS_BLOCK_HEADER, size + QAWS_BLOCK_HEADER);
		if (!raw)
			return NULL;
		mem_header(raw, NULL, size);
		return raw + QAWS_BLOCK_HEADER;
	}
	/* grown in the same allocator: a new block, the contents, the old block released */
	p = qaws_internal_alloc_from(allocator, size);
	if (!p)
		return NULL;
	memcpy(p, ptr, old < size ? old : size);
	qaws_internal_free(ptr);
	return p;
}

/* ---- arena over caller memory ------------------------------------------- */

static void* arena_alloc(unsigned long size, void* user)
{
	qaws_arena* a = (qaws_arena*)user;
	size_t at = (a->used + 15) & ~(size_t)15;
	if (at > a->size || (size_t)size > a->size - at)
		return NULL;
	a->last = at;
	a->used = at + (size_t)size;
	if (a->used > a->peak)
		a->peak = a->used;
	return a->base + at;
}

/* releases only the latest block (growing arrays reuse their space); the rest at qaws_arena_reset */
static void arena_dealloc(void* ptr, void* user)
{
	qaws_arena* a = (qaws_arena*)user;
	if (ptr && (unsigned char*)ptr == a->base + a->last)
		a->used = a->last;
}

void qaws_arena_init(qaws_arena* arena, void* memory, size_t size)
{
	size_t skip = (size_t)((16 - ((size_t)memory & 15)) & 15);
	memset(arena, 0, sizeof *arena);
	arena->base = (unsigned char*)memory + (skip < size ? skip : size);
	arena->size = skip < size ? size - skip : 0;
	arena->allocator.alloc = arena_alloc;
	arena->allocator.dealloc = arena_dealloc;
	arena->allocator.user_data = arena;
}

void qaws_arena_reset(qaws_arena* arena)
{
	arena->used = 0;
	arena->last = 0;
}
