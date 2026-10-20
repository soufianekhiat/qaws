/*
 * Test 154: where the library's memory comes from (qaws_alloc.h, internal/qaws_internal_memory.h)
 *
 *   - a descriptor that was a stack local: the object outlives it and is
 *     destroyed through the copy, after the stack is reused;
 *   - the default: settable until first used, then refused; a mixed workload
 *     sees all its blocks come back; set against the first allocation of
 *     seven workers in 24 fresh processes (installed or sealed first, never
 *     a partly written descriptor);
 *   - a scratch arena apart from results: sampling, offsets, clipping,
 *     tessellation, meshing and a repulsive solve leave nothing in it, their
 *     results live on the heap past its resets;
 *   - checked size products and sums, list capacities, an oversized net
 *     refused before its points are read;
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
	qaws_curve_create_bspline(&c, &d);
	return c;
}

/* the parts mm_workload runs: sampling, offsets, clipping, tessellation */
enum { MM_SAMPLE = 1, MM_OFFSET = 2, MM_CLIP = 4, MM_TESS = 8, MM_MESH = 16, MM_REPULSIVE = 32, MM_CLIP64 = 64, MM_ALL = 15 };
static unsigned int g_parts = MM_ALL;
static unsigned int g_published;   /* results a failed call left in its outputs */

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
	if (g_parts & (MM_SAMPLE | MM_OFFSET))
		failures += qaws_curve_create_bspline(&w, &d) != QAWS_STATUS_OK;
	if (w)
	{
		memset(&sd, 0, sizeof sd);
		sd.sample_count = 4;
		sd.error_tolerance = (qaws_scalar)1e-4;
		sd.use_adaptive_sampling = 1;
		if (g_parts & MM_SAMPLE)
			failures += qaws_curve_sample_2d(samples, &n, w, 4096, &sd) != QAWS_STATUS_OK;
		n = 0;
		if (g_parts & MM_OFFSET)
			failures += qaws_curve_offset_2d(offs, &n, w, 16, (qaws_scalar)0.3, 1) != QAWS_STATUS_OK;
		for (i = 0; i < n; i++)
			qaws_curve_destroy(offs[i]);
		qaws_curve_destroy(w);
	}
	a = g_parts & MM_CLIP ? mm_square(0, 0, 2) : NULL;
	b = g_parts & MM_CLIP ? mm_square(1, 1, 2) : NULL;
	if (!(g_parts & MM_CLIP))
		;
	else if (a && b)
	{
		ca[0] = a; cb[0] = b;
		pa.curves = ca; pa.curve_count = 1; pa.closed = 1;
		pb.curves = cb; pb.curve_count = 1; pb.closed = 1;
		failures += qaws_clip_boolean(&r, &pa, &pb, 1, 1, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO) != QAWS_STATUS_OK;
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
	if (g_parts & (MM_TESS | MM_MESH))
		failures += qaws_surface_create_bezier(&s, &bd) != QAWS_STATUS_OK;
	if (s && (g_parts & MM_TESS))
	{
		memset(&td, 0, sizeof td);
		td.max_depth = 4;
		td.curvature_threshold = (qaws_scalar)0.01;
		failures += qaws_surface_tessellate(tv, &nv, ti, &ni, s, 4096, 24000, &td) != QAWS_STATUS_OK;
	}
	if (s && (g_parts & MM_MESH))
	{
		qaws_generate_desc gd;
		qaws_mesh_3d* m = NULL;
		qaws_status st;
		memset(&gd, 0, sizeof gd);
		gd.method = QAWS_GENERATE_REMESH;
		gd.target_length = (qaws_scalar)0.5;
		st = qaws_mesh_generate_surface(&m, NULL, s, &gd);
		failures += st != QAWS_STATUS_OK;
		g_published += st != QAWS_STATUS_OK && m != NULL;   /* a failed call publishes no mesh */
		qaws_mesh_3d_destroy(m);
	}
	qaws_surface_destroy(s);
	if (g_parts & MM_CLIP64)
	{
		/* integer clipping: two overlapping squares and an open path through them */
		static int64_t const sa[8] = { 0, 0, 2000, 0, 2000, 2000, 0, 2000 }, sb[8] = { 1000, 1000, 3000, 1000, 3000, 3000, 1000, 3000 };
		static int64_t const op[4] = { -500, 1500, 3500, 1700 };
		qaws_path64 pa64 = { sa, 4 }, pb64 = { sb, 4 }, po64 = { op, 2 };
		qaws_clip64_desc cd;
		qaws_clip64_result* r64 = NULL;
		qaws_status st;
		memset(&cd, 0, sizeof cd);
		cd.subjects = &pa64; cd.subject_count = 1; cd.open_subjects = &po64; cd.open_subject_count = 1;
		cd.clips = &pb64; cd.clip_count = 1; cd.clip_type = QAWS_CLIP_UNION; cd.fill_rule = QAWS_FILL_NON_ZERO;
		st = qaws_clip64_execute(&r64, &cd);
		failures += st != QAWS_STATUS_OK;
		g_published += st != QAWS_STATUS_OK && r64 != NULL;
		qaws_clip64_result_destroy(r64);
	}
	if (g_parts & MM_REPULSIVE)
	{
		qaws_repulsive_curves_desc rd;
		qaws_repulsive_curves_solve_desc sd;
		qaws_repulsive_curves* rc = NULL;
		qaws_scalar pos[48];
		unsigned int edges[32];
		qaws_status st;
		for (i = 0; i < 16; i++)
		{
			double t = 6.283185307179586 * i / 16;
			pos[3 * i] = (qaws_scalar)cos(t);
			pos[3 * i + 1] = (qaws_scalar)(sin(t) * (1 + 0.2 * cos(3 * t)));
			pos[3 * i + 2] = (qaws_scalar)(0.2 * sin(2 * t));
			edges[2 * i] = i;
			edges[2 * i + 1] = (i + 1) % 16;
		}
		memset(&rd, 0, sizeof rd);
		rd.dimension = 3; rd.vertex_count = 16; rd.edges = edges; rd.edge_count = 16;
		st = qaws_repulsive_curves_create(&rc, &rd);
		failures += st != QAWS_STATUS_OK;
		g_published += st != QAWS_STATUS_OK && rc != NULL;
		if (rc)
		{
			memset(&sd, 0, sizeof sd);
			sd.max_iterations = 2;
			failures += qaws_repulsive_curves_solve(NULL, pos, rc, &sd) != QAWS_STATUS_OK;
			qaws_repulsive_curves_destroy(rc);
		}
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
	qaws_curve_create_bezier_ex(&c, &d, &local);
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
	TEST_ASSERT(qaws_curve_reverse(&r, c) == QAWS_STATUS_OK, "derived from it after the descriptor's scope");
	qaws_curve_evaluate_2d(&e, r, (qaws_scalar)0.25, QAWS_EVAL_FLAG_POSITION);
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

/* the first allocation against qaws_set_default_allocator: in a fresh process (a child of this runner), seven workers
   allocate and free from the moment thread 0 installs the default. The default is installed (OK) or sealed by a worker
   first (UNSUPPORTED_OPERATION); either way every block returns to the allocator that gave it, and the default's
   callbacks never see a partly written descriptor (their user_data checked on every call). */
#define MM_RACE_MAGIC 0x5EA1u
typedef struct mm_race
{
	unsigned int magic;
	long volatile allocs, frees, bad;
} mm_race;
static mm_race g_race = { MM_RACE_MAGIC, 0, 0, 0 };
static long volatile g_race_go;
/* the child's index mod 3: thread 0 installs at once (0), once a worker is about to allocate (1), once one has allocated (2) */
static unsigned int g_race_mode;
static long volatile g_race_entered, g_race_done;

static long mm_inc(long volatile* p)
{
#ifdef _WIN32
	return InterlockedIncrement(p);
#else
	return __atomic_add_fetch(p, 1, __ATOMIC_SEQ_CST);
#endif
}

static void* mm_race_alloc(unsigned long size, void* user)
{
	mm_race* r = (mm_race*)user;
	if (r != &g_race || r->magic != MM_RACE_MAGIC)
		mm_inc(&g_race.bad);
	mm_inc(&g_race.allocs);
	return malloc(size);
}

static void mm_race_free(void* p, void* user)
{
	mm_race* r = (mm_race*)user;
	if (r != &g_race || r->magic != MM_RACE_MAGIC)
		mm_inc(&g_race.bad);
	mm_inc(&g_race.frees);
	free(p);
}

static qaws_status g_race_status = QAWS_STATUS_INTERNAL_ERROR;
static long volatile g_race_failed;

static void mm_race_worker(void* p)
{
	unsigned int i;
	if (p)
	{
		/* thread 0: release the workers, then install the default */
		qaws_allocator a;
		memset(&a, 0, sizeof a);
		a.alloc = mm_race_alloc;
		a.dealloc = mm_race_free;
		a.user_data = &g_race;
		mm_inc(&g_race_go);
		while ((g_race_mode == 1 && !g_race_entered) || (g_race_mode == 2 && !g_race_done)) {}
		g_race_status = qaws_set_default_allocator(&a);
		memset(&a, 0xCD, sizeof a);   /* the library holds its own copy */
		return;
	}
	while (!g_race_go) {}
	for (i = 0; i < 2000; i++)
	{
		void* b;
		if (i == 0)
			mm_inc(&g_race_entered);
		b = qaws_internal_malloc(16 + i % 200);
		if (i == 0)
			mm_inc(&g_race_done);
		if (!b)
			mm_inc(&g_race_failed);
		else
			memset(b, 0x5A, 16);
		qaws_internal_free(b);
	}
}

static int test_default_race_child(void)
{
	void* args[8] = { &g_race, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
	int sealed_first;
	g_race_mode = ((unsigned int)atoi(getenv("QAWS_TEST_RACE")) - 1) % 3;
	test_run_threads(8, mm_race_worker, args);
	sealed_first = g_race_status == QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (g_race_failed || g_race.bad || g_race.allocs != g_race.frees || (g_race_status != QAWS_STATUS_OK && !sealed_first)
		|| (sealed_first && g_race.allocs != 0))
	{
		printf("  race: status %s, %ld blocks, %ld returned, %ld bad descriptors, %ld failures\n", qaws_status_to_string(g_race_status),
			g_race.allocs, g_race.frees, g_race.bad, g_race_failed);
		return 1;
	}
	return sealed_first ? 100 : 0;
}

static void test_default_race(void)
{
#ifdef _WIN32
	char path[MAX_PATH];
	char saved[300];
	unsigned int k, installed = 0, sealed = 0, bad = 0;
	char msg[200];
	if (!GetModuleFileNameA(NULL, path, sizeof path))
	{
		TEST_ASSERT(0, "default race: the runner's path");
		return;
	}
	{
		char const* was = getenv("QAWS_TEST");
		snprintf(saved, sizeof saved, "QAWS_TEST=%s", was ? was : "");
	}
	_putenv("QAWS_TEST=154");
	for (k = 0; k < 24; k++)
	{
		char var[40];
		intptr_t code;
		snprintf(var, sizeof var, "QAWS_TEST_RACE=%u", k + 1);
		_putenv(var);
		code = _spawnl(_P_WAIT, path, path, NULL);
		if (code == 0)
			installed++;
		else if (code == 100)
			sealed++;
		else
			bad++;
	}
	_putenv("QAWS_TEST_RACE=");
	_putenv(saved);   /* the runner's own setting back */
	snprintf(msg, sizeof msg, "default against first allocation, 24 fresh processes: installed in %u, sealed by a worker first in %u, wrong in %u",
		installed, sealed, bad);
	printf("  %s\n", msg);
	TEST_ASSERT(bad == 0 && installed + sealed == 24, msg);
#else
	/* shortcut: the child-process race runs on Windows only; on other targets add a fork/exec path when they are tested */
#endif
}

/* ---- scratch apart from results ------------------------------------------ */

/* after a reset: what the arena held is overwritten (a result still pointing into it reads 0xCD), then a fresh arena */
static void mm_poison(qaws_arena* arena, unsigned char* memory, size_t size)
{
	memset(memory, 0xCD, arena->peak < size ? arena->peak + 64 : size);
	qaws_arena_init(arena, memory, size);
}

/* each operation with a frame arena as the thread's scratch and a counting heap as its allocator: when the call returns
   the arena holds nothing (the temporaries came from it and went back), so it resets while the results live on the heap;
   the results then go back to the heap. */
static void test_scratch_split(void)
{
	static unsigned char memory[1 << 22];
	static qaws_vec2 samples[4096];
	static qaws_tessellation_vertex tv[8192];
	static unsigned int ti[49152];
	qaws_vec2 wave[6] = { { 0, 0 }, { 1, 2 }, { 2, -2 }, { 3, 2 }, { 4, -2 }, { 5, 0 } };
	qaws_vec3 net[16];
	mm_count heap = { 0, 0, 0 };
	qaws_allocator ha = mm_legacy(&heap), prev, prev_scratch;
	qaws_arena arena;
	qaws_bspline_desc d;
	qaws_curve *w = NULL, *a, *b, *offs[16];
	qaws_curve const* ca[1];
	qaws_curve const* cb[1];
	qaws_path_2d pa, pb;
	qaws_clip_result* r = NULL;
	qaws_surface* s = NULL;
	qaws_mesh_3d* m = NULL;
	qaws_repulsive_curves* rc = NULL;
	qaws_mesh_reference* ref = NULL;
	qaws_size_field* field = NULL;
	unsigned int n = 0, noff = 0, nv = 0, ni = 0, i, j, ok = 1, held = 0, edges[48];
	size_t scratch_peak[6];
	qaws_scalar pos[72];
	char msg[300];
	qaws_arena_init(&arena, memory, sizeof memory);
	qaws_exchange_thread_allocator(&prev, &ha);
	qaws_exchange_thread_scratch_allocator(&prev_scratch, &arena.allocator);
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = wave; d.control_point_count = 6;
	ok &= qaws_curve_create_bspline(&w, &d) == QAWS_STATUS_OK;
	/* adaptive sampling */
	{
		qaws_sampling_desc sd;
		memset(&sd, 0, sizeof sd);
		sd.sample_count = 4; sd.error_tolerance = (qaws_scalar)1e-5; sd.use_adaptive_sampling = 1;
		ok &= qaws_curve_sample_2d(samples, &n, w, 4096, &sd) == QAWS_STATUS_OK;
		held += arena.live != 0; scratch_peak[0] = arena.peak; ok &= qaws_arena_reset(&arena) == QAWS_STATUS_OK; mm_poison(&arena, memory, sizeof memory);
	}
	/* offsets: curves kept past the reset */
	ok &= qaws_curve_offset_2d(offs, &noff, w, 16, (qaws_scalar)0.3, 1) == QAWS_STATUS_OK;
	held += arena.live != 0; scratch_peak[1] = arena.peak; ok &= qaws_arena_reset(&arena) == QAWS_STATUS_OK; mm_poison(&arena, memory, sizeof memory);
	/* clipping: the result kept */
	a = mm_square(0, 0, 2);
	b = mm_square(1, 1, 2);
	ca[0] = a; cb[0] = b;
	pa.curves = ca; pa.curve_count = 1; pa.closed = 1;
	pb.curves = cb; pb.curve_count = 1; pb.closed = 1;
	ok &= qaws_clip_boolean(&r, &pa, &pb, 1, 1, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO) == QAWS_STATUS_OK;
	held += arena.live != 0; scratch_peak[2] = arena.peak; ok &= qaws_arena_reset(&arena) == QAWS_STATUS_OK; mm_poison(&arena, memory, sizeof memory);
	/* tessellation */
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			net[i * 4 + j].x = (qaws_scalar)i;
			net[i * 4 + j].y = (qaws_scalar)j;
			net[i * 4 + j].z = (qaws_scalar)sin(1.3 * i + 0.7 * j);
		}
	{
		qaws_surface_bezier_desc bd;
		qaws_tessellation_desc td;
		memset(&bd, 0, sizeof bd);
		bd.u_degree = 3; bd.v_degree = 3; bd.control_points = net; bd.u_point_count = 4; bd.v_point_count = 4;
		ok &= qaws_surface_create_bezier(&s, &bd) == QAWS_STATUS_OK;
		memset(&td, 0, sizeof td);
		td.max_depth = 5; td.curvature_threshold = (qaws_scalar)0.003;
		ok &= qaws_surface_tessellate(tv, &nv, ti, &ni, s, 8192, 49152, &td) == QAWS_STATUS_OK;
		held += arena.live != 0; scratch_peak[3] = arena.peak; ok &= qaws_arena_reset(&arena) == QAWS_STATUS_OK; mm_poison(&arena, memory, sizeof memory);
	}
	/* meshing: the mesh kept */
	{
		qaws_generate_desc gd;
		memset(&gd, 0, sizeof gd);
		gd.method = QAWS_GENERATE_REMESH;
		gd.target_length = (qaws_scalar)0.3;
		ok &= qaws_mesh_generate_surface(&m, NULL, s, &gd) == QAWS_STATUS_OK;
		held += arena.live != 0; scratch_peak[4] = arena.peak; ok &= qaws_arena_reset(&arena) == QAWS_STATUS_OK; mm_poison(&arena, memory, sizeof memory);
	}
	/* a curvature size field on that mesh: kept past the reset, read below against one made without scratch */
	if (m)
	{
		qaws_curvature_size_desc cd;
		memset(&cd, 0, sizeof cd);
		cd.chord_tolerance = (qaws_scalar)0.01;
		cd.gradation = (qaws_scalar)0.5;
		ok &= qaws_mesh_reference_create(&ref, m, NULL) == QAWS_STATUS_OK;
		ok &= ref && qaws_size_field_create_curvature(&field, ref, &cd) == QAWS_STATUS_OK;
		held += arena.live != 0; ok &= qaws_arena_reset(&arena) == QAWS_STATUS_OK; mm_poison(&arena, memory, sizeof memory);
	}
	/* a repulsive solve: the solver object kept */
	{
		qaws_repulsive_curves_desc rd;
		qaws_repulsive_curves_solve_desc sd;
		for (i = 0; i < 24; i++)
		{
			double t = 6.283185307179586 * i / 24;
			pos[3 * i] = (qaws_scalar)(cos(t) * (1 + 0.3 * cos(5 * t)));
			pos[3 * i + 1] = (qaws_scalar)(sin(t) * (1 + 0.3 * cos(5 * t)));
			pos[3 * i + 2] = (qaws_scalar)(0.2 * sin(3 * t));
			edges[2 * i] = i;
			edges[2 * i + 1] = (i + 1) % 24;
		}
		memset(&rd, 0, sizeof rd);
		rd.dimension = 3; rd.vertex_count = 24; rd.edges = edges; rd.edge_count = 24;
		ok &= qaws_repulsive_curves_create(&rc, &rd) == QAWS_STATUS_OK;
		memset(&sd, 0, sizeof sd);
		sd.max_iterations = 3;
		sd.hold_total_length = 1;
		ok &= qaws_repulsive_curves_solve(NULL, pos, rc, &sd) == QAWS_STATUS_OK;
		held += arena.live != 0; scratch_peak[5] = arena.peak; ok &= qaws_arena_reset(&arena) == QAWS_STATUS_OK; mm_poison(&arena, memory, sizeof memory);
	}
	qaws_exchange_thread_scratch_allocator(NULL, &prev_scratch);
	snprintf(msg, sizeof msg, "scratch arena apart from results: sampling, offsets, clipping, tessellation, meshing, repulsive; held after a call %u; "
		"scratch peaks %lu, %lu, %lu, %lu, %lu, %lu bytes", held, (unsigned long)scratch_peak[0],
		(unsigned long)scratch_peak[1], (unsigned long)scratch_peak[2], (unsigned long)scratch_peak[3], (unsigned long)scratch_peak[4],
		(unsigned long)scratch_peak[5]);
	printf("  %s\n", msg);
	TEST_ASSERT(ok && held == 0 && arena.invalid_frees == 0 && arena.failures == 0 && scratch_peak[5] > 0 && noff > 0 && r && m && rc, msg);
	/* the results are whole on the heap: used after the resets, then returned */
	ok = qaws_clip_result_get_path_count(r) > 0 && qaws_mesh_3d_get_face_count(m) > 0;
	{
		static qaws_scalar xyz[3 * 65536];
		qaws_path_2d path;
		ok &= qaws_mesh_3d_get_vertex_count(m) <= 65536 && qaws_mesh_3d_get_positions(xyz, m) == QAWS_STATUS_OK
			&& qaws_clip_result_get_path(&path, r, 0) == QAWS_STATUS_OK && path.curve_count > 0;
		/* the size field against one made without scratch, at the corners of the first triangles */
		if (ok && ref && field)
		{
			qaws_curvature_size_desc cd;
			qaws_size_field* plain = NULL;
			qaws_mesh_location loc[16];
			qaws_scalar at[48], s1[16], s2[16];
			unsigned int t, tv3[3], cnt = qaws_mesh_reference_get_triangle_count(ref) < 16 ? qaws_mesh_reference_get_triangle_count(ref) : 16;
			memset(&cd, 0, sizeof cd);
			cd.chord_tolerance = (qaws_scalar)0.01;
			cd.gradation = (qaws_scalar)0.5;
			ok &= qaws_size_field_create_curvature(&plain, ref, &cd) == QAWS_STATUS_OK;
			for (t = 0; t < cnt; t++)
			{
				qaws_mesh_reference_get_triangle(tv3, ref, t);
				loc[t].patch = t; loc[t].a = 1; loc[t].b = 0;
				at[3 * t] = xyz[3 * tv3[0]]; at[3 * t + 1] = xyz[3 * tv3[0] + 1]; at[3 * t + 2] = xyz[3 * tv3[0] + 2];
			}
			ok &= plain && cnt > 0 && qaws_size_field_evaluate(s1, field, loc, at, cnt) == QAWS_STATUS_OK
				&& qaws_size_field_evaluate(s2, plain, loc, at, cnt) == QAWS_STATUS_OK && memcmp(s1, s2, cnt * sizeof s1[0]) == 0;
			qaws_size_field_destroy(plain);
		}
	}
	for (i = 0; i < noff; i++)
		qaws_curve_destroy(offs[i]);
	qaws_clip_result_destroy(r);
	qaws_size_field_destroy(field);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(m);
	qaws_repulsive_curves_destroy(rc);
	qaws_surface_destroy(s);
	qaws_curve_destroy(a);
	qaws_curve_destroy(b);
	qaws_curve_destroy(w);
	qaws_set_thread_allocator(&prev);
	TEST_ASSERT(ok && heap.allocs > 0 && heap.frees == heap.allocs && arena.live == 0, "the results used after the resets and returned to the heap");
}

static void test_thread_nesting(void)
{
	mm_count outer = { 0, 0, 0 }, inner = { 0, 0, 0 };
	qaws_allocator ao = mm_legacy(&outer), ai = mm_legacy(&inner), prev1, prev2, now;
	long base_default = g_default_count.allocs;
	qaws_curve* c;
	qaws_exchange_thread_allocator(&prev1, &ao);
	TEST_ASSERT(prev1.alloc == NULL, "no thread allocator to begin with");
	qaws_exchange_thread_allocator(&prev2, &ai);
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
	qaws_allocator local = j->a;   /* a stack descriptor, gone before the blocks are freed on another thread */
	unsigned int i;
	qaws_set_thread_allocator(&local);
	for (i = 0; i < 8; i++)
		j->curves[i] = mm_square(i, 0, 1);
	qaws_set_thread_allocator(NULL);
	memset(&local, 0xCD, sizeof local);
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

/* records the pointer it is given back (the test keeps the few bytes) */
static void* g_shift_last;
static void mm_shifted_free(void* p, void* user)
{
	g_shift_last = p;
	((mm_count*)user)->frees++;
}

/* the sized interface, 8 bytes off the alignment it is asked */
static void* mm_shifted_sized(size_t size, size_t alignment, void* user)
{
	unsigned char* p = (unsigned char*)malloc(size + 64);
	((mm_count*)user)->allocs++;
	return p ? p + 8 + (alignment - ((uintptr_t)p % alignment)) % alignment : NULL;
}

static void mm_shifted_sized_free(void* p, size_t size, size_t alignment, void* user)
{
	(void)size; (void)alignment;
	g_shift_last = p;
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

	/* the legacy interface cannot be told the alignment: blocks 8 bytes off are aligned inside, freed at the address given */
	memset(&shifted, 0, sizeof shifted);
	shifted.alloc = mm_shifted_alloc;
	shifted.dealloc = mm_shifted_free;
	shifted.user_data = &g_shift;
	{
		unsigned char* q = (unsigned char*)qaws_internal_alloc_from(64, &shifted);
		TEST_ASSERT(q && (uintptr_t)q % QAWS_MEMORY_ALIGNMENT == 0, "a legacy allocator off alignment: the block aligned inside");
		memset(q, 0x5A, 64);
		qaws_internal_free(q);
		TEST_ASSERT(g_shift.frees == 1 && g_shift_last && (unsigned char*)g_shift_last < q && q - (unsigned char*)g_shift_last < 256,
			"and given back at the address the allocator returned");
	}
#ifdef NDEBUG
	/* the sized interface is told the alignment: a block off it is given back and the allocation fails */
	{
		qaws_allocator ss;
		mm_count sc = { 0, 0, 0 };
		memset(&ss, 0, sizeof ss);
		ss.alloc = qaws_allocator_extended;
		ss.alloc_sized = mm_shifted_sized;
		ss.dealloc_sized = mm_shifted_sized_free;
		ss.user_data = &sc;
		TEST_ASSERT(qaws_internal_alloc_from(64, &ss) == NULL && sc.allocs == 1 && sc.frees == 1,
			"a sized allocator off the alignment asked: refused, the block given back");
	}
#else
	printf("  (a sized allocator off its alignment asserts in debug builds; checked in release)\n");
#endif

	TEST_ASSERT(qaws_internal_malloc(SIZE_MAX) == NULL && qaws_internal_malloc(SIZE_MAX - 8) == NULL, "sizes near SIZE_MAX refused");
	TEST_ASSERT(qaws_internal_zalloc(SIZE_MAX / 2, 3) == NULL && qaws_internal_zalloc(SIZE_MAX, SIZE_MAX) == NULL, "count times size overflow refused");
	/* the checked size arithmetic: exact below the limit, SIZE_MAX (refused by every allocation) past it */
	TEST_ASSERT(qaws_internal_mul2(SIZE_MAX / 2, 2) == SIZE_MAX - 1 && qaws_internal_mul2(SIZE_MAX / 2 + 1, 2) == SIZE_MAX
		&& qaws_internal_mul2(0, SIZE_MAX) == 0 && qaws_internal_mul3(SIZE_MAX / 4 + 1, 2, 2) == SIZE_MAX
		&& qaws_internal_mul4(65536, 65536, 65536, 65536) == SIZE_MAX
		&& qaws_internal_mul3(0xFFFFFFFFu, 3, 8) == (sizeof(size_t) > 4 ? (size_t)0xFFFFFFFFu * 24 : SIZE_MAX)
		&& qaws_internal_add2(SIZE_MAX - 1, 1) == SIZE_MAX && qaws_internal_add2(SIZE_MAX, 1) == SIZE_MAX
		&& qaws_internal_malloc(qaws_internal_mul3(SIZE_MAX / 3, 2, 2)) == NULL, "checked size products and sums");
	TEST_ASSERT(qaws_internal_grow(0, 1, 16) == 16 && qaws_internal_grow(16, 17, 16) == 32 && qaws_internal_grow(16, 100, 16) == 128
		&& qaws_internal_grow(0x80000000u, 0x80000001u, 16) == 0xFFFFFFFFu && (sizeof(size_t) == 4 || qaws_internal_grow(0xFFFFFFFFu, (size_t)0xFFFFFFFFu + 1, 16) == 0),   /* (32-bit: no list reaches UINT_MAX elements) */
		"list capacities: doubled to the need, held at UINT_MAX, 0 past it");
	/* oversized topology refused before the control points are read (they are 16 here) */
	{
		qaws_vec3 net16[16];
		qaws_surface_bspline_desc sd;
		qaws_surface* s = NULL;
		qaws_status st;
		memset(net16, 0, sizeof net16);
		memset(&sd, 0, sizeof sd);
		sd.u_degree = 3; sd.v_degree = 3; sd.control_points = net16; sd.u_point_count = 0x10000; sd.v_point_count = 0x10000;
		st = qaws_surface_create_bspline(&s, &sd);
		snprintf(msg, sizeof msg, "a 65536 x 65536 net (2^32 points, past unsigned int): %s", qaws_status_to_string(st));
		TEST_ASSERT(st == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT && !s, msg);
		if (sizeof(size_t) == 4)
		{
			/* 32-bit: 2^29 points of three scalars overflow size_t; refused before they are read */
			qaws_bspline_desc bd;
			qaws_curve* c = NULL;
			memset(&bd, 0, sizeof bd);
			bd.dimension = QAWS_DIMENSION_3D; bd.degree = 3; bd.control_points = net16; bd.control_point_count = 0x20000000;
			st = qaws_curve_create_bspline(&c, &bd);
			snprintf(msg, sizeof msg, "32-bit: 2^29 control points of three scalars: %s", qaws_status_to_string(st));
			TEST_ASSERT(st != QAWS_STATUS_OK && !c, msg);
		}
	}
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
		TEST_ASSERT(qaws_internal_alloc_from(five_gb, &sized) == NULL && rec.calls == 1 && rec.requested > five_gb
			&& rec.requested - five_gb <= 256, "5 GB through the sized interface, not truncated");
		legacy = mm_legacy(&lc);
		qaws_internal_alloc_from(five_gb, &legacy);
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
	qaws_exchange_thread_allocator(&prev, &arena.allocator);
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
		a = (unsigned char*)qaws_internal_alloc_from(64, &arena.allocator);
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
		a = (unsigned char*)qaws_internal_alloc_from(256, &s2.allocator);
		b = (unsigned char*)qaws_internal_alloc_from(256, &s2.allocator);
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
		TEST_ASSERT(qaws_curve_create_bezier_ex(&c, &bd, &arena.allocator) == QAWS_STATUS_OK
			&& (unsigned char*)c >= (unsigned char*)big && (unsigned char*)c < (unsigned char*)big + sizeof big, "a curve inside the caller's buffer");
		TEST_ASSERT(qaws_arena_reset(&arena) == QAWS_STATUS_INVALID_ARGUMENT, "reset refused while the curve lives");
		qaws_curve_evaluate_2d(&e, c, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(fabs(e.position.x - 1.5) <= 1e-12, "and the curve still evaluates");
		qaws_curve_destroy(c);
		TEST_ASSERT(qaws_arena_reset(&arena) == QAWS_STATUS_OK, "reset once it is destroyed");
	}

	/* too tiny_buffer: failures reported, nothing taken from the heap */
	qaws_arena_init(&arena, tiny_buffer, sizeof tiny_buffer);
	qaws_exchange_thread_allocator(&prev, &arena.allocator);
	failures = mm_workload();
	qaws_set_thread_allocator(&prev);
	TEST_ASSERT(failures > 0 && arena.peak <= sizeof tiny_buffer && arena.live == 0, "a 2 KB arena: allocation failures reported cleanly, nothing held");

	/* every part in arenas from 256 bytes up (allocations fail at many different points): nothing may stay held */
	{
		static char const* names[7] = { "sampling", "offsets", "clipping", "tessellation", "meshing", "repulsive", "integer clipping" };
		static double sweep[1 << 16];
		unsigned int part, leaks = 0, runs = 0;
		size_t sz;
		g_published = 0;
		for (part = 0; part < 7; part++)
		{
			unsigned int part_leaks = 0;
			g_parts = 1u << part;
			for (sz = 256; sz <= sizeof sweep; sz = sz * 9 / 8 + 16)
			{
				qaws_arena_init(&arena, sweep, sz);
				qaws_exchange_thread_allocator(&prev, &arena.allocator);
				mm_workload();
				qaws_set_thread_allocator(&prev);
				runs++;
				if (arena.live || arena.invalid_frees)
				{
					if (!part_leaks)
						printf("    %s: %lu blocks held after failing in a %lu-byte arena\n", names[part], (unsigned long)arena.live, (unsigned long)sz);
					part_leaks++;
				}
			}
			leaks += part_leaks;
		}
		g_parts = MM_ALL;
		snprintf(msg, sizeof msg, "failure paths: %u runs in arenas from 256 bytes, %u left blocks held, %u published a result on failure", runs,
			leaks, g_published);
		printf("  %s\n", msg);
		TEST_ASSERT(leaks == 0 && g_published == 0, msg);
	}

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
	qaws_exchange_thread_allocator(&prev, &w->arena.allocator);
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
	{
		char const* race = getenv("QAWS_TEST_RACE");
		if (race && *race)
			exit(test_default_race_child());   /* the child's verdict as its exit code (the runner would fold it) */
	}
	printf("Test 154: allocators and caller memory\n");
	test_default();
	test_default_race();
	test_descriptor_lifetime();
	test_thread_nesting();
	test_cross_thread();
	test_scratch_split();
	test_sizes_and_alignment();
	test_arena();
	test_concurrent();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
