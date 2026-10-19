/*
 * Test 154: where the library's memory comes from (qaws_alloc.h, internal/qaws_internal_memory.h)
 *
 *   - a descriptor that was a stack local: the object outlives it and is
 *     destroyed through the copy, after the stack is reused;
 *   - the default: settable until first used, then refused; a mixed workload
 *     sees all its blocks come back;
 *   - thread allocators: before the default, exchanged and restored when
 *     nested, blocks made on one thread freed on another under other settings;
 *   - misaligned allocators refused (release builds: debug builds assert);
 *   - sizes near SIZE_MAX refused; a failed reallocation keeps the old block;
 *     5 GB requested through the sized interface without truncation, refused
 *     by the legacy one where unsigned long is 32 bits;
 *   - the arena: the workload inside a caller buffer, in-place growth of the
 *     latest block, reset refused while blocks are held, a full arena failing
 *     cleanly, frees it does not hold counted.
 */

#include "test_common.h"
#include "test_threads.h"
#include "internal/qaws_internal_memory.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

typedef struct mm_count
{
	long allocs, frees;
	size_t last_size;
} mm_count;

static void* mm_alloc(unsigned long size, void* user)
{
	((mm_count*)user)->allocs++;
	return malloc(size);
}

static void mm_free(void* p, void* user)
{
	((mm_count*)user)->frees++;
	free(p);
}

static qaws_allocator mm_legacy(mm_count* c)
{
	qaws_allocator a;
	memset(&a, 0, sizeof a);
	a.alloc = mm_alloc;
	a.dealloc = mm_free;
	a.user_data = c;
	return a;
}

static qaws_curve* mm_square(double x, double y, double s)
{
	qaws_vec2 p[5];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	p[0].x = (qaws_scalar)x; p[0].y = (qaws_scalar)y;
	p[1].x = (qaws_scalar)(x + s); p[1].y = (qaws_scalar)y;
	p[2].x = (qaws_scalar)(x + s); p[2].y = (qaws_scalar)(y + s);
	p[3].x = (qaws_scalar)x; p[3].y = (qaws_scalar)(y + s);
	p[4] = p[0];
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = p; d.control_point_count = 5;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

/* a mixed workload; its failures counted */
static unsigned int mm_workload(void)
{
	static qaws_vec2 samples[4096];
	static qaws_tessellation_vertex tv[4096];
	static unsigned int ti[24000];
	qaws_vec2 wave[6] = { { 0, 0 }, { 1, 2 }, { 2, -2 }, { 3, 2 }, { 4, -2 }, { 5, 0 } };
	qaws_vec3 net[16];
	qaws_bspline_desc d;
	qaws_curve *a, *b, *w = NULL, *offs[16];
	qaws_curve const* ca[1];
	qaws_curve const* cb[1];
	qaws_path_2d pa, pb;
	qaws_clip_result* r = NULL;
	qaws_sampling_desc sd;
	qaws_surface_bezier_desc bd;
	qaws_tessellation_desc td;
	qaws_surface* s = NULL;
	unsigned int failures = 0, n = 0, nv = 0, ni = 0, i, j;

	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = wave; d.control_point_count = 6;
	failures += qaws_curve_create_bspline(&d, &w) != QAWS_STATUS_OK;
	if (w)
	{
		memset(&sd, 0, sizeof sd);
		sd.sample_count = 4;
		sd.error_tolerance = (qaws_scalar)1e-4;
		sd.use_adaptive_sampling = 1;
		failures += qaws_curve_sample_2d(w, &sd, samples, 4096, &n) != QAWS_STATUS_OK;
		n = 0;
		failures += qaws_curve_offset_2d(w, (qaws_scalar)0.3, 1, offs, 16, &n) != QAWS_STATUS_OK;
		for (i = 0; i < n; i++)
			qaws_curve_destroy(offs[i]);
		qaws_curve_destroy(w);
	}
	a = mm_square(0, 0, 2);
	b = mm_square(1, 1, 2);
	if (a && b)
	{
		ca[0] = a; cb[0] = b;
		pa.curves = ca; pa.curve_count = 1; pa.closed = 1;
		pb.curves = cb; pb.curve_count = 1; pb.closed = 1;
		failures += qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &pa, 1, &pb, 1, &r) != QAWS_STATUS_OK;
		qaws_clip_result_destroy(r);
	}
	else
		failures++;
	qaws_curve_destroy(a);
	qaws_curve_destroy(b);
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			net[i * 4 + j].x = (qaws_scalar)i;
			net[i * 4 + j].y = (qaws_scalar)j;
			net[i * 4 + j].z = (qaws_scalar)sin(1.3 * i + 0.7 * j);
		}
	memset(&bd, 0, sizeof bd);
	bd.u_degree = 3; bd.v_degree = 3; bd.control_points = net; bd.u_point_count = 4; bd.v_point_count = 4;
	failures += qaws_surface_create_bezier(&bd, &s) != QAWS_STATUS_OK;
	if (s)
	{
		memset(&td, 0, sizeof td);
		td.max_depth = 4;
		td.curvature_threshold = (qaws_scalar)0.01;
		failures += qaws_surface_tessellate(s, &td, tv, 4096, &nv, ti, 24000, &ni) != QAWS_STATUS_OK;
		qaws_surface_destroy(s);
	}
	return failures;
}

/* ---- descriptor lifetime --------------------------------------------- */

static mm_count g_scoped;

/* the descriptor is this frame's local; the curve outlives it */
static qaws_curve* mm_make_scoped(void)
{
	qaws_allocator local = mm_legacy(&g_scoped);
	qaws_vec2 p[4] = { { 0, 0 }, { 1, 1 }, { 2, 0 }, { 3, 1 } };
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p; d.control_point_count = 4;
	qaws_curve_create_bezier_ex(&d, &local, &c);
	return c;
}

/* fills the stack where the descriptor was */
static void mm_clobber(void)
{
	volatile unsigned char junk[512];
	unsigned int i;
	for (i = 0; i < sizeof junk; i++)
		junk[i] = 0xCD;
}

static void test_descriptor_lifetime(void)
{
	qaws_curve* c = mm_make_scoped();
	qaws_curve* r = NULL;
	qaws_eval_result_2d e;
	mm_clobber();
	TEST_ASSERT(c && g_scoped.allocs > 0, "a curve made with a stack-local descriptor");
	/* allocating later through the object's allocator (a reversed copy) uses the object's copy too */
	TEST_ASSERT(qaws_curve_reverse(c, &r) == QAWS_STATUS_OK, "derived from it after the descriptor's scope");
	qaws_curve_evaluate_2d(r, (qaws_scalar)0.25, QAWS_EVAL_FLAG_POSITION, &e);
	qaws_curve_destroy(r);
	qaws_curve_destroy(c);
	TEST_ASSERT(g_scoped.frees == g_scoped.allocs, "destroyed through the copy: every block back to its allocator");
}

/* ---- default and thread allocators ------------------------------------- */

static mm_count g_default_count;
static qaws_allocator g_default_desc;

static void test_default(void)
{
	unsigned int failures;
	char msg[200];
	qaws_status st;
	g_default_desc = mm_legacy(&g_default_count);
	st = qaws_set_default_allocator(&g_default_desc);
	if (st == QAWS_STATUS_OK)
	{
		failures = mm_workload();
		snprintf(msg, sizeof msg, "default allocator: %ld blocks, %ld returned (%u failures)", g_default_count.allocs, g_default_count.frees,
			failures);
		printf("  %s\n", msg);
		TEST_ASSERT(failures == 0 && g_default_count.allocs > 50 && g_default_count.frees == g_default_count.allocs, msg);
		TEST_ASSERT(qaws_set_default_allocator(NULL) == QAWS_STATUS_UNSUPPORTED_OPERATION, "the default is sealed once used");
	}
	else
	{
		/* the whole runner in one process: an earlier suite allocated first */
		TEST_ASSERT(st == QAWS_STATUS_UNSUPPORTED_OPERATION, "the default is sealed once used (by an earlier suite)");
	}
}

static void test_thread_nesting(void)
{
	mm_count outer = { 0, 0, 0 }, inner = { 0, 0, 0 };
	qaws_allocator ao = mm_legacy(&outer), ai = mm_legacy(&inner), prev1, prev2, now;
	long base_default = g_default_count.allocs;
	qaws_curve* c;
	qaws_exchange_thread_allocator(&ao, &prev1);
	TEST_ASSERT(prev1.alloc == NULL, "no thread allocator to begin with");
	qaws_exchange_thread_allocator(&ai, &prev2);
	TEST_ASSERT(prev2.alloc == mm_alloc && prev2.user_data == &outer, "the nested exchange returns the outer one");
	c = mm_square(0, 0, 1);
	qaws_curve_destroy(c);
	qaws_set_thread_allocator(&prev2);
	qaws_get_thread_allocator(&now);
	TEST_ASSERT(now.user_data == &outer, "restored to the outer one");
	c = mm_square(0, 0, 1);
	qaws_curve_destroy(c);
	qaws_set_thread_allocator(&prev1);
	qaws_get_thread_allocator(&now);
	TEST_ASSERT(now.alloc == NULL, "restored to none");
	TEST_ASSERT(inner.allocs > 0 && inner.frees == inner.allocs && outer.allocs > 0 && outer.frees == outer.allocs
		&& g_default_count.allocs == base_default, "each level took its blocks; none from the default");
}

/* blocks made on worker threads with their allocators, freed on another thread under other settings */
typedef struct mm_job
{
	mm_count count;
	qaws_allocator a;
	qaws_curve* curves[8];
} mm_job;

static void mm_make(void* p)
{
	mm_job* j = (mm_job*)p;
	unsigned int i;
	qaws_set_thread_allocator(&j->a);
	for (i = 0; i < 8; i++)
		j->curves[i] = mm_square(i, 0, 1);
	qaws_set_thread_allocator(NULL);
}


/* thread 0 (the caller) idles, thread 1 destroys every job's curves under its own allocator */
static void mm_destroy(void* p)
{
	mm_job* jobs = (mm_job*)p;
	mm_count other = { 0, 0, 0 };
	qaws_allocator o;
	unsigned int t, i;
	if (!jobs)
		return;
	o = mm_legacy(&other);
	qaws_set_thread_allocator(&o);
	for (t = 0; t < 4; t++)
		for (i = 0; i < 8; i++)
			qaws_curve_destroy(jobs[t].curves[i]);
	qaws_set_thread_allocator(NULL);
}

static void test_cross_thread(void)
{
	static mm_job jobs[4];
	void* args[4];
	void* d[2];
	unsigned int t, ok = 1;
	for (t = 0; t < 4; t++)
	{
		memset(&jobs[t], 0, sizeof jobs[t]);
		jobs[t].a = mm_legacy(&jobs[t].count);
		args[t] = &jobs[t];
	}
	test_run_threads(4, mm_make, args);
	d[0] = NULL;
	d[1] = jobs;
	test_run_threads(2, mm_destroy, d);
	for (t = 0; t < 4; t++)
		ok &= jobs[t].count.allocs > 0 && jobs[t].count.frees == jobs[t].count.allocs;
	TEST_ASSERT(ok, "four workers' blocks freed on another thread under another allocator, each back to its own");
}

/* ---- alignment, sizes, reallocation ------------------------------------ */

static mm_count g_shift;

/* blocks 8 bytes off alignment */
static void* mm_shifted_alloc(unsigned long size, void* user)
{
	unsigned char* p = (unsigned char*)malloc((size_t)size + 32);
	((mm_count*)user)->allocs++;
	return p ? p + 8 + (16 - ((uintptr_t)p & 15)) % 16 : NULL;
}

static void mm_shifted_free(void* p, void* user)
{
	/* blocks are only ever given back here unused: the test keeps them (a few bytes) */
	(void)p;
	((mm_count*)user)->frees++;
}

typedef struct mm_sized
{
	size_t requested;
	long calls;
} mm_sized;

/* records the request, never allocates */
static void* mm_record(size_t size, size_t alignment, void* user)
{
	mm_sized* s = (mm_sized*)user;
	(void)alignment;
	s->requested = size;
	s->calls++;
	return NULL;
}

static void mm_nofree(void* p, size_t size, size_t alignment, void* user)
{
	(void)p; (void)size; (void)alignment; (void)user;
}

static void test_sizes_and_alignment(void)
{
	qaws_allocator shifted, legacy, sized;
	mm_count lc = { 0, 0, 0 };
	mm_sized rec = { 0, 0 };
	size_t five_gb = (size_t)5 << 30;
	unsigned char* p;
	unsigned int i;
	char msg[200];

#ifdef NDEBUG
	memset(&shifted, 0, sizeof shifted);
	shifted.alloc = mm_shifted_alloc;
	shifted.dealloc = mm_shifted_free;
	shifted.user_data = &g_shift;
	TEST_ASSERT(qaws_internal_alloc_from(&shifted, 64) == NULL && g_shift.allocs == 1 && g_shift.frees == 1,
		"a misaligned block is given back and the allocation fails");
#else
	(void)shifted;
	printf("  (misaligned allocators assert in debug builds; checked in release)\n");
#endif

	TEST_ASSERT(qaws_internal_malloc(SIZE_MAX) == NULL && qaws_internal_malloc(SIZE_MAX - 8) == NULL, "sizes near SIZE_MAX refused");
	TEST_ASSERT(qaws_internal_zalloc(SIZE_MAX / 2, 3) == NULL && qaws_internal_zalloc(SIZE_MAX, SIZE_MAX) == NULL, "count times size overflow refused");
	p = (unsigned char*)qaws_internal_malloc(100);
	for (i = 0; i < 100; i++)
		p[i] = (unsigned char)i;
	TEST_ASSERT(qaws_internal_realloc(p, SIZE_MAX - 4) == NULL, "reallocation near SIZE_MAX refused");
	for (i = 0; i < 100 && p[i] == (unsigned char)i; i++)
		;
	TEST_ASSERT(i == 100 && qaws_internal_block_size(p) == 100, "the old block kept, contents and size");
	TEST_ASSERT(((uintptr_t)p % QAWS_MEMORY_ALIGNMENT) == 0, "blocks aligned to QAWS_MEMORY_ALIGNMENT");
	qaws_internal_free(p);

	/* 5 GB: the sized interface sees it whole; the legacy one where unsigned long cannot hold it is not called */
	memset(&sized, 0, sizeof sized);
	sized.alloc = qaws_allocator_extended;
	sized.alloc_sized = mm_record;
	sized.dealloc_sized = mm_nofree;
	sized.user_data = &rec;
	if (sizeof(size_t) > 4)
	{
		TEST_ASSERT(qaws_internal_alloc_from(&sized, five_gb) == NULL && rec.calls == 1 && rec.requested > five_gb
			&& rec.requested - five_gb <= 256, "5 GB through the sized interface, not truncated");
		legacy = mm_legacy(&lc);
		qaws_internal_alloc_from(&legacy, five_gb);
		snprintf(msg, sizeof msg, "5 GB through the legacy interface (unsigned long of %u bytes): %s", (unsigned)sizeof(unsigned long),
			sizeof(unsigned long) < 8 ? "refused without calling it" : "passed whole");
		printf("  %s\n", msg);
		TEST_ASSERT(sizeof(unsigned long) < 8 ? lc.allocs == 0 : lc.allocs == 1, msg);
	}
}

/* ---- the arena ----------------------------------------------------------- */

static void test_arena(void)
{
	static double big[1 << 18];       /* 2 MB, 16-byte aligned through double */
	static double tiny_buffer[256];
	qaws_arena arena;
	unsigned int failures;
	char msg[200];
	qaws_allocator prev;

	qaws_arena_init(&arena, big, sizeof big);
	qaws_exchange_thread_allocator(&arena.allocator, &prev);
	failures = mm_workload();
	qaws_set_thread_allocator(&prev);
	snprintf(msg, sizeof msg, "arena over a 2 MB caller buffer: the workload ran (%u failures), peak %lu bytes, %lu blocks held after",
		failures, (unsigned long)arena.peak, (unsigned long)arena.live);
	printf("  %s\n", msg);
	TEST_ASSERT(failures == 0 && arena.peak > 0 && arena.live == 0, msg);
	TEST_ASSERT(qaws_arena_reset(&arena) == QAWS_STATUS_OK && arena.used == 0 && arena.generation == 1, "reset with nothing held");

	/* in-place growth of the latest block */
	{
		unsigned char *a, *b;
		a = (unsigned char*)qaws_internal_alloc_from(&arena.allocator, 64);
		b = (unsigned char*)qaws_internal_realloc(a, 4096);
		TEST_ASSERT(a && b == a, "the latest block grows in place");
		qaws_internal_free(b);
	}

	/* a failed reallocation keeps the block (not the latest: no room to grow in place) */
	{
		qaws_arena s2;
		unsigned char *a, *b;
		unsigned int i;
		qaws_arena_init(&s2, tiny_buffer, sizeof tiny_buffer);
		a = (unsigned char*)qaws_internal_alloc_from(&s2.allocator, 256);
		b = (unsigned char*)qaws_internal_alloc_from(&s2.allocator, 256);
		for (i = 0; i < 256; i++)
			a[i] = (unsigned char)(i * 7);
		TEST_ASSERT(a && b && qaws_internal_realloc(a, 4000) == NULL, "a reallocation the arena cannot hold fails");
		for (i = 0; i < 256 && a[i] == (unsigned char)(i * 7); i++)
			;
		TEST_ASSERT(i == 256, "and keeps the old block whole");
		qaws_internal_free(b);
		qaws_internal_free(a);
		TEST_ASSERT(s2.live == 0 && s2.failures >= 1, "blocks released, the refusal counted");
	}

	/* reset refused while an object lives in the arena */
	{
		qaws_vec2 p[4] = { { 0, 0 }, { 1, 1 }, { 2, 0 }, { 3, 1 } };
		qaws_bezier_desc bd;
		qaws_curve* c = NULL;
		qaws_eval_result_2d e;
		memset(&bd, 0, sizeof bd);
		bd.dimension = QAWS_DIMENSION_2D; bd.degree = 3; bd.control_points = p; bd.control_point_count = 4;
		TEST_ASSERT(qaws_curve_create_bezier_ex(&bd, &arena.allocator, &c) == QAWS_STATUS_OK
			&& (unsigned char*)c >= (unsigned char*)big && (unsigned char*)c < (unsigned char*)big + sizeof big, "a curve inside the caller's buffer");
		TEST_ASSERT(qaws_arena_reset(&arena) == QAWS_STATUS_INVALID_ARGUMENT, "reset refused while the curve lives");
		qaws_curve_evaluate_2d(c, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &e);
		TEST_ASSERT(fabs(e.position.x - 1.5) <= 1e-12, "and the curve still evaluates");
		qaws_curve_destroy(c);
		TEST_ASSERT(qaws_arena_reset(&arena) == QAWS_STATUS_OK, "reset once it is destroyed");
	}

	/* too tiny_buffer: failures reported, nothing taken from the heap */
	qaws_arena_init(&arena, tiny_buffer, sizeof tiny_buffer);
	qaws_exchange_thread_allocator(&arena.allocator, &prev);
	failures = mm_workload();
	qaws_set_thread_allocator(&prev);
	TEST_ASSERT(failures > 0 && arena.peak <= sizeof tiny_buffer && arena.live == 0, "a 2 KB arena: allocation failures reported cleanly, nothing held");

	/* overflow-safe arithmetic in the arena itself */
	qaws_arena_init(&arena, tiny_buffer, sizeof tiny_buffer);
	TEST_ASSERT(arena.allocator.alloc_sized(SIZE_MAX, 16, &arena) == NULL && arena.allocator.alloc_sized(SIZE_MAX - 15, 16, &arena) == NULL
		&& arena.used == 0, "arena requests near SIZE_MAX refused");
}

/* ---- concurrent work, one arena per worker ------------------------------ */

typedef struct mm_worker
{
	qaws_arena arena;
	unsigned int failures;
} mm_worker;

static void mm_work(void* p)
{
	mm_worker* w = (mm_worker*)p;
	qaws_allocator prev;
	unsigned int k;
	qaws_exchange_thread_allocator(&w->arena.allocator, &prev);
	/* a frame arena: reset between workloads, refused if anything were still held */
	for (k = 0; k < 3; k++)
	{
		w->failures += mm_workload();
		w->failures += qaws_arena_reset(&w->arena) != QAWS_STATUS_OK;
	}
	qaws_set_thread_allocator(&prev);
}

static void test_concurrent(void)
{
	static double buffers[4][1 << 16];      /* 512 KB each */
	static double tiny[4][64];
	static mm_worker w[4];
	void* args[4];
	unsigned int t, ok = 1, failed = 1;
	for (t = 0; t < 4; t++)
	{
		qaws_arena_init(&w[t].arena, buffers[t], sizeof buffers[t]);
		w[t].failures = 0;
		args[t] = &w[t];
	}
	test_run_threads(4, mm_work, args);
	for (t = 0; t < 4; t++)
		ok &= w[t].failures == 0 && w[t].arena.live == 0 && w[t].arena.invalid_frees == 0;
	TEST_ASSERT(ok, "four workers, each in its own frame arena reset between three workloads, at once: all ran, nothing held");
	for (t = 0; t < 4; t++)
	{
		qaws_arena_init(&w[t].arena, tiny[t], sizeof tiny[t]);
		w[t].failures = 0;
	}
	test_run_threads(4, mm_work, args);
	for (t = 0; t < 4; t++)
		failed &= w[t].failures > 0 && w[t].arena.live == 0 && w[t].arena.invalid_frees == 0;
	TEST_ASSERT(failed, "four workers failing for memory at once: clean failures, nothing held");
}

int test_154_memory_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 154: allocators and caller memory\n");
	test_default();
	test_descriptor_lifetime();
	test_thread_nesting();
	test_cross_thread();
	test_sizes_and_alignment();
	test_arena();
	test_concurrent();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
