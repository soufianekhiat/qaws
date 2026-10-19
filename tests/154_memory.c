/*
 * Test 154: where the library's memory comes from (qaws_alloc.h)
 *
 *   - a counting default allocator sees every block of a mixed workload
 *     (curves, adaptive sampling, offsets, polygon clipping, tessellation),
 *     and every block comes back to it;
 *   - a thread allocator takes over from the default, on its thread only;
 *   - a block freed after the settings change returns to the allocator that
 *     gave it;
 *   - the workload inside an arena over a caller buffer: all of it there, the
 *     buffer's peak; too small a buffer: a clean allocation failure.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef struct mm_count
{
	long allocs, frees, live;
} mm_count;

static void* mm_alloc(unsigned long size, void* user)
{
	mm_count* c = (mm_count*)user;
	c->allocs++;
	c->live++;
	return malloc(size);
}

static void mm_free(void* p, void* user)
{
	mm_count* c = (mm_count*)user;
	c->frees++;
	c->live--;
	free(p);
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

/* a mixed workload; its failures counted (allocation failures under a small arena) */
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

static void test_default_and_thread(void)
{
	mm_count cd = { 0, 0, 0 }, ct = { 0, 0, 0 };
	qaws_allocator ad, at;
	unsigned int failures;
	char msg[200];
	ad.alloc = mm_alloc; ad.dealloc = mm_free; ad.user_data = &cd;
	at.alloc = mm_alloc; at.dealloc = mm_free; at.user_data = &ct;

	qaws_set_default_allocator(&ad);
	failures = mm_workload();
	qaws_set_default_allocator(NULL);
	snprintf(msg, sizeof msg, "default allocator: %ld blocks, %ld returned, the workload ran (%u failures)", cd.allocs, cd.frees, failures);
	printf("  %s\n", msg);
	TEST_ASSERT(failures == 0 && cd.allocs > 50 && cd.live == 0, msg);

	/* the thread allocator first */
	cd.allocs = cd.frees = 0;
	qaws_set_default_allocator(&ad);
	qaws_set_thread_allocator(&at);
	mm_workload();
	qaws_set_thread_allocator(NULL);
	qaws_set_default_allocator(NULL);
	TEST_ASSERT(ct.allocs > 50 && ct.live == 0 && cd.allocs == 0, "thread allocator before the default, all returned");

	/* made under one setting, destroyed under another: back to its own allocator */
	{
		qaws_curve* c;
		long before;
		qaws_set_default_allocator(&ad);
		c = mm_square(0, 0, 1);
		qaws_set_default_allocator(&at);
		before = cd.frees;
		qaws_curve_destroy(c);
		qaws_set_default_allocator(NULL);
		TEST_ASSERT(cd.frees > before && cd.live == 0 && ct.live == 0, "a block returns to the allocator that gave it");
	}
}

static void test_arena(void)
{
	static double big[1 << 18];       /* 2 MB, 16-byte aligned through double */
	static double small[256];
	qaws_arena arena;
	unsigned int failures;
	char msg[200];
	qaws_arena_init(&arena, big, sizeof big);
	qaws_set_thread_allocator(&arena.allocator);
	failures = mm_workload();
	qaws_set_thread_allocator(NULL);
	snprintf(msg, sizeof msg, "arena over a 2 MB caller buffer: the workload ran (%u failures), peak %lu bytes", failures, (unsigned long)arena.peak);
	printf("  %s\n", msg);
	TEST_ASSERT(failures == 0 && arena.peak > 0 && arena.peak <= sizeof big, msg);
	qaws_arena_reset(&arena);
	TEST_ASSERT(arena.used == 0, "reset");

	/* too small: failures reported, nothing crashes, nothing taken from the heap */
	qaws_arena_init(&arena, small, sizeof small);
	qaws_set_thread_allocator(&arena.allocator);
	failures = mm_workload();
	qaws_set_thread_allocator(NULL);
	TEST_ASSERT(failures > 0 && arena.peak <= sizeof small, "a 2 KB arena: allocation failures reported cleanly");

	/* an owned curve placed in caller memory explicitly */
	{
		qaws_vec2 p[4] = { { 0, 0 }, { 1, 1 }, { 2, 0 }, { 3, 1 } };
		qaws_bezier_desc bd;
		qaws_curve* c = NULL;
		qaws_eval_result_2d e;
		qaws_arena_init(&arena, big, sizeof big);
		memset(&bd, 0, sizeof bd);
		bd.dimension = QAWS_DIMENSION_2D; bd.degree = 3; bd.control_points = p; bd.control_point_count = 4;
		TEST_ASSERT(qaws_curve_create_bezier_ex(&bd, &arena.allocator, &c) == QAWS_STATUS_OK
			&& (unsigned char*)c >= (unsigned char*)big && (unsigned char*)c < (unsigned char*)big + sizeof big, "a curve inside the caller's buffer");
		qaws_curve_evaluate_2d(c, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &e);
		TEST_ASSERT(fabs(e.position.x - 1.5) <= 1e-12, "evaluated there");
		qaws_curve_destroy(c);
	}
}

int test_154_memory_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 154: allocators and caller memory\n");
	test_default_and_thread();
	test_arena();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
