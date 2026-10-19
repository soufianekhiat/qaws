/*
 * memory_bench.c - heap against arena allocation
 *
 * For each operation: the time per run with the C allocator and inside a qaws_arena (the thread allocator), the blocks
 * and the most live bytes on the heap (a counting sized allocator), and the arena's peak: what growing arrays cost an
 * arena that releases only its latest block (peak over live).
 */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct mb_count
{
	size_t live, peak;
	long blocks;
} mb_count;

static void* mb_alloc(size_t size, size_t alignment, void* user)
{
	mb_count* c = (mb_count*)user;
	void* p;
	(void)alignment;
#if defined(_MSC_VER)
	p = _aligned_malloc(size, 16);
#else
	p = aligned_alloc(16, (size + 15) / 16 * 16);
#endif
	if (p)
	{
		c->blocks++;
		c->live += size;
		if (c->live > c->peak)
			c->peak = c->live;
	}
	return p;
}

static void mb_free(void* p, size_t size, size_t alignment, void* user)
{
	mb_count* c = (mb_count*)user;
	(void)alignment;
	c->live -= size;
#if defined(_MSC_VER)
	_aligned_free(p);
#else
	free(p);
#endif
}

static qaws_curve* mb_wave(void)
{
	qaws_vec2 p[40];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 40; i++)
	{
		p[i].x = (qaws_scalar)i;
		p[i].y = (qaws_scalar)(2 * sin(0.9 * i) + 0.5 * cos(2.3 * i));
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p; d.control_point_count = 40;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static qaws_surface* mb_patch(void)
{
	qaws_vec3 n[12 * 12];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 12; i++)
		for (j = 0; j < 12; j++)
		{
			n[i * 12 + j].x = (qaws_scalar)i;
			n[i * 12 + j].y = (qaws_scalar)j;
			n[i * 12 + j].z = (qaws_scalar)(1.2 * sin(0.8 * i) * cos(0.7 * j));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3; d.v_degree = 3; d.control_points = n; d.u_point_count = 12; d.v_point_count = 12;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

static qaws_curve* mb_square(double x, double y, double s)
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

/* ---- the operations: inputs made once, the run allocates its own --------- */

static qaws_curve* g_wave;
static qaws_surface* g_patch;
static qaws_curve* g_squares[64];

static int op_sampling(void)
{
	static qaws_vec2 out[200000];
	qaws_sampling_desc sd;
	unsigned int n = 0;
	memset(&sd, 0, sizeof sd);
	sd.sample_count = 8;
	sd.error_tolerance = (qaws_scalar)1e-6;
	sd.use_adaptive_sampling = 1;
	return qaws_curve_sample_2d(g_wave, &sd, out, 200000, &n) == QAWS_STATUS_OK;
}

static int op_tessellation(void)
{
	static qaws_tessellation_vertex v[200000];
	static unsigned int idx[600000];
	qaws_tessellation_desc td;
	unsigned int nv = 0, ni = 0;
	memset(&td, 0, sizeof td);
	td.max_depth = 7;
	td.curvature_threshold = (qaws_scalar)1e-3;
	return qaws_surface_tessellate(g_patch, &td, v, 200000, &nv, idx, 600000, &ni) == QAWS_STATUS_OK;
}

static int op_meshing(void)
{
	qaws_generate_desc gd;
	qaws_mesh_3d* m = NULL;
	int ok;
	memset(&gd, 0, sizeof gd);
	gd.method = QAWS_GENERATE_REMESH;
	gd.target_length = (qaws_scalar)0.25;
	ok = qaws_mesh_generate_surface(g_patch, &gd, &m, NULL) == QAWS_STATUS_OK;
	qaws_mesh_3d_destroy(m);
	return ok;
}

static int op_clipping(void)
{
	qaws_path_2d paths[64];
	qaws_curve const* curves[64];
	qaws_clip_result* r = NULL;
	unsigned int i;
	int ok;
	for (i = 0; i < 64; i++)
	{
		curves[i] = g_squares[i];
		paths[i].curves = &curves[i];
		paths[i].curve_count = 1;
		paths[i].closed = 1;
	}
	ok = qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, paths, 64, NULL, 0, &r) == QAWS_STATUS_OK;
	qaws_clip_result_destroy(r);
	return ok;
}

static int op_offsets(void)
{
	qaws_curve* out[64];
	unsigned int n = 0, i;
	int ok = qaws_curve_offset_2d(g_wave, (qaws_scalar)0.4, 1, out, 64, &n) == QAWS_STATUS_OK;
	for (i = 0; i < n; i++)
		qaws_curve_destroy(out[i]);
	return ok;
}

static int op_repulsive(void)
{
	enum { N = 48 };
	unsigned int edges[2 * N], i;
	qaws_scalar pos[3 * N];
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_curves* rc = NULL;
	int ok;
	for (i = 0; i < N; i++)
	{
		double a = 6.283185307179586 * i / N;
		pos[3 * i] = (qaws_scalar)(cos(a) * (1 + 0.3 * cos(5 * a)));
		pos[3 * i + 1] = (qaws_scalar)(sin(a) * (1 + 0.3 * cos(5 * a)));
		pos[3 * i + 2] = (qaws_scalar)(0.2 * sin(3 * a));
		edges[2 * i] = i;
		edges[2 * i + 1] = (i + 1) % N;
	}
	memset(&d, 0, sizeof d);
	d.dimension = 3; d.vertex_count = N; d.edges = edges; d.edge_count = N;
	if (qaws_repulsive_curves_create(&d, &rc) != QAWS_STATUS_OK)
		return 0;
	memset(&sd, 0, sizeof sd);
	sd.max_iterations = 5;
	sd.hold_total_length = 1;
	ok = qaws_repulsive_curves_solve(rc, &sd, pos, NULL) == QAWS_STATUS_OK;
	qaws_repulsive_curves_destroy(rc);
	return ok;
}

typedef struct mb_op
{
	char const* name;
	int (*run)(void);
	int reps;
} mb_op;

static double mb_seconds(int (*run)(void), int reps, int* ok)
{
	clock_t t0 = clock();
	int r;
	*ok = 1;
	for (r = 0; r < reps; r++)
		*ok &= run();
	return (double)(clock() - t0) / CLOCKS_PER_SEC / reps;
}

int main(void)
{
	static mb_op ops[] = {
		{ "adaptive sampling", op_sampling, 50 },
		{ "tessellation", op_tessellation, 5 },
		{ "meshing (remesh)", op_meshing, 2 },
		{ "clipping (64 squares)", op_clipping, 20 },
		{ "offsets", op_offsets, 20 },
		{ "repulsive (48 vertices, 5 steps)", op_repulsive, 3 },
	};
	size_t arena_bytes = (size_t)512 << 20;
	unsigned char* memory = (unsigned char*)malloc(arena_bytes);
	unsigned int i, k;
	if (!memory)
		return 1;
	g_wave = mb_wave();
	g_patch = mb_patch();
	for (k = 0; k < 64; k++)
		g_squares[k] = mb_square((k % 8) * 0.7, (k / 8) * 0.7, 1.0);
	printf("%-34s %12s %12s %10s %14s %14s %8s\n", "operation", "heap ms", "arena ms", "blocks", "heap live peak", "arena peak", "ratio");
	for (i = 0; i < sizeof ops / sizeof ops[0]; i++)
	{
		mb_count c;
		qaws_allocator counting, prev;
		qaws_arena arena;
		double th, ta;
		int ok1, ok2, ok3;
		/* heap: the C allocator, then a counting sized allocator for the live peak */
		th = mb_seconds(ops[i].run, ops[i].reps, &ok1);
		memset(&c, 0, sizeof c);
		memset(&counting, 0, sizeof counting);
		counting.alloc = qaws_allocator_extended;
		counting.alloc_sized = mb_alloc;
		counting.dealloc_sized = mb_free;
		counting.user_data = &c;
		qaws_exchange_thread_allocator(&counting, &prev);
		ok2 = ops[i].run();
		qaws_set_thread_allocator(&prev);
		/* arena */
		qaws_arena_init(&arena, memory, arena_bytes);
		qaws_exchange_thread_allocator(&arena.allocator, &prev);
		ta = mb_seconds(ops[i].run, 1, &ok3);
		qaws_set_thread_allocator(&prev);
		printf("%-34s %12.3f %12.3f %10ld %14lu %14lu %8.2f%s\n", ops[i].name, th * 1e3, ta * 1e3, c.blocks, (unsigned long)c.peak,
			(unsigned long)arena.peak, c.peak ? (double)arena.peak / (double)c.peak : 0.0, ok1 && ok2 && ok3 ? "" : "  (a run failed)");
	}
	for (k = 0; k < 64; k++)
		qaws_curve_destroy(g_squares[k]);
	qaws_curve_destroy(g_wave);
	qaws_surface_destroy(g_patch);
	free(memory);
	return 0;
}
