/*
 * Test 182: synthetic stress of the hot calls, timed
 *
 *   Each case runs one call over deterministic synthetic inputs (a seeded generator) and prints its
 *   operation count, wall and thread CPU time, the time per operation and a hash of every output
 *   bit, so a faster build is checked against the slower one's hash at the same precision.
 *
 *   QAWS_STRESS=k (k >= 1) scales the counts by k; unset, each case runs a small count (a check
 *   that it runs, the times not meaningful). QAWS_STRESS_ONLY=substring runs the matching cases.
 *
 *   Cases: curve evaluation (one point, batches), arc length and its inverse, closest parameter,
 *   bounds, tessellation, winding numbers, rigorous enclosures, offsets, curve sets' closest
 *   points; surface evaluation (one point, batches), closest point, tessellation, area, rigorous
 *   enclosures, surface sets' closest points; region point location, polyline simplification;
 *   the exact orientation predicates, winding and batch closest points.
 */
#include "test_common.h"
#include "qaws_curve_batch.h"
#include "qaws_surface_batch.h"
#include "qaws_exact.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
#else
# include <time.h>
#endif

static unsigned int g_scale;          /* 0: smoke counts */
static char const* g_only;
static unsigned long long g_rng;

static double st_rand(void)           /* [0, 1) */
{
	g_rng = g_rng * 6364136223846793005ull + 1442695040888963407ull;
	return (double)(g_rng >> 11) * (1.0 / 9007199254740992.0);
}

static unsigned int st_count(unsigned int full)
{
	unsigned int n = g_scale ? full * g_scale : full / 1000;
	return n ? n : 1;
}

/* every bit of the outputs, 8 bytes a step (a byte a step would cost as much as a cheap call) */
static unsigned long long st_hash(unsigned long long h, void const* data, size_t bytes)
{
	unsigned char const* p = (unsigned char const*)data;
	size_t i;
	for (i = 0; i + 8 <= bytes; i += 8)
	{
		unsigned long long w;
		memcpy(&w, p + i, 8);
		h = (h ^ w) * 0x100000001B3ull;
		h ^= h >> 29;
	}
	for (; i < bytes; i++)
		h = (h ^ p[i]) * 1099511628211ull;
	return h;
}

typedef struct st_clock
{
	double wall, cpu;
} st_clock;

static st_clock st_now(void)
{
	st_clock c;
#ifdef _WIN32
	LARGE_INTEGER f, t;
	FILETIME a, b, k, u;
	QueryPerformanceFrequency(&f);
	QueryPerformanceCounter(&t);
	c.wall = (double)t.QuadPart * 1000.0 / (double)f.QuadPart;
	GetThreadTimes(GetCurrentThread(), &a, &b, &k, &u);
	c.cpu = (double)(((unsigned long long)u.dwHighDateTime << 32 | u.dwLowDateTime)
	                 + ((unsigned long long)k.dwHighDateTime << 32 | k.dwLowDateTime)) * 1e-4;
#else
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	c.wall = (double)t.tv_sec * 1000.0 + (double)t.tv_nsec * 1e-6;
	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
	c.cpu = (double)t.tv_sec * 1000.0 + (double)t.tv_nsec * 1e-6;
#endif
	return c;
}

static int st_want(char const* name)
{
	return !g_only || strstr(name, g_only) != NULL;
}

static void st_report(char const* name, unsigned int ops, st_clock t0, unsigned long long hash, int ok)
{
	st_clock t1 = st_now();
	double wall = t1.wall - t0.wall, cpu = t1.cpu - t0.cpu;
	char msg[160];
	printf("  STRESS %-34s %9u ops %10.1f ms wall %10.1f ms cpu %10.1f ns/op  hash %016llx\n", name, ops, wall, cpu,
	       cpu * 1e6 / (double)ops, hash);
	sprintf(msg, "%s: every call succeeded", name);
	TEST_ASSERT(ok, msg);
}

/* ---- the inputs ---- */

#define ST_CP 64

typedef struct st_inputs
{
	qaws_curve* bspline3;     /* cubic 3D B-spline, 64 points on a noisy helix */
	qaws_curve* nurbs3;       /* cubic 3D NURBS on the same points, weights 0.5..2 */
	qaws_curve* bezier2;      /* quintic 2D Bezier */
	qaws_curve* loop2;        /* closed cubic 2D B-spline, a 32-point star */
	qaws_curve* chain2[16];   /* cubic 2D Beziers round a wavy circle */
	qaws_surface* surface;    /* bicubic B-spline, 16 x 16 wavy points */
	qaws_exact_curve* exact2[16];
} st_inputs;

static void st_make(st_inputs* in)
{
	qaws_vec3 p3[ST_CP];
	qaws_scalar w[ST_CP], knots[ST_CP + 4], p2[2 * 32], b2[12];
	qaws_vec3 grid[16 * 16];
	qaws_bspline_desc bd;
	qaws_nurbs_desc nd;
	qaws_bezier_desc zd;
	qaws_surface_bspline_desc sd;
	qaws_exact_desc ed;
	unsigned int i, j;
	memset(in, 0, sizeof *in);
	g_rng = 12345;
	for (i = 0; i < ST_CP; i++)
	{
		double a = 0.35 * i;
		p3[i].x = (qaws_scalar)(cos(a) + 0.1 * st_rand());
		p3[i].y = (qaws_scalar)(sin(a) + 0.1 * st_rand());
		p3[i].z = (qaws_scalar)(0.05 * i);
		w[i] = (qaws_scalar)(0.5 + 1.5 * st_rand());
	}
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_3D, bd.degree = 3, bd.control_points = p3, bd.control_point_count = ST_CP, bd.is_uniform = 1;
	qaws_curve_create_bspline(&in->bspline3, &bd);
	for (i = 0; i < ST_CP + 4; i++)
		knots[i] = (qaws_scalar)(i < 4 ? 0 : i >= ST_CP ? ST_CP - 3 : i - 3);
	memset(&nd, 0, sizeof nd);
	nd.dimension = QAWS_DIMENSION_3D, nd.degree = 3, nd.control_points = p3, nd.control_point_count = ST_CP;
	nd.knots = knots, nd.knot_count = ST_CP + 4, nd.weights = w, nd.weight_count = ST_CP;
	qaws_curve_create_nurbs(&in->nurbs3, &nd);
	for (i = 0; i < 6; i++)
		b2[2 * i] = (qaws_scalar)i, b2[2 * i + 1] = (qaws_scalar)((i & 1) ? 2 : -1);
	memset(&zd, 0, sizeof zd);
	zd.dimension = QAWS_DIMENSION_2D, zd.degree = 5, zd.control_points = b2, zd.control_point_count = 6;
	qaws_curve_create_bezier(&in->bezier2, &zd);
	for (i = 0; i < 32; i++)
	{
		double a = 6.283185307179586 * i / 32, r = (i & 1) ? 0.5 : 1.0;
		p2[2 * i] = (qaws_scalar)(r * cos(a)), p2[2 * i + 1] = (qaws_scalar)(r * sin(a));
	}
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_2D, bd.degree = 3, bd.control_points = p2, bd.control_point_count = 32, bd.is_uniform = 1, bd.is_closed = 1;
	qaws_curve_create_bspline(&in->loop2, &bd);
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = 0;
	for (i = 0; i < 16; i++)
	{
		qaws_scalar c[8];
		double a0 = 6.283185307179586 * i / 16, a1 = 6.283185307179586 * (i + 1) / 16, h = 0.5523 * (a1 - a0) / 1.5708;
		double r0 = 0.8 + 0.1 * cos(3 * a0), r1 = 0.8 + 0.1 * cos(3 * a1);
		c[0] = (qaws_scalar)(r0 * cos(a0)), c[1] = (qaws_scalar)(r0 * sin(a0));
		c[2] = (qaws_scalar)(r0 * (cos(a0) - h * sin(a0))), c[3] = (qaws_scalar)(r0 * (sin(a0) + h * cos(a0)));
		c[4] = (qaws_scalar)(r1 * (cos(a1) + h * sin(a1))), c[5] = (qaws_scalar)(r1 * (sin(a1) - h * cos(a1)));
		c[6] = (qaws_scalar)(r1 * cos(a1)), c[7] = (qaws_scalar)(r1 * sin(a1));
		if (i == 15)
			c[6] = (qaws_scalar)(0.9 * cos(0.0)), c[7] = 0;   /* the loop closes on piece 0's start */
		memset(&zd, 0, sizeof zd);
		zd.dimension = QAWS_DIMENSION_2D, zd.degree = 3, zd.control_points = c, zd.control_point_count = 4;
		qaws_curve_create_bezier(&in->chain2[i], &zd);
		qaws_exact_curve_prepare(&in->exact2[i], NULL, in->chain2[i], &ed);
	}
	for (i = 0; i < 16; i++)
		for (j = 0; j < 16; j++)
		{
			grid[16 * i + j].x = (qaws_scalar)i, grid[16 * i + j].y = (qaws_scalar)j;
			grid[16 * i + j].z = (qaws_scalar)(0.8 * sin(0.7 * i) * cos(0.5 * j) + 0.2 * st_rand());
		}
	memset(&sd, 0, sizeof sd);
	sd.u_degree = 3, sd.v_degree = 3, sd.control_points = grid, sd.u_point_count = 16, sd.v_point_count = 16;
	qaws_surface_create_bspline(&in->surface, &sd);
}

static void st_free(st_inputs* in)
{
	unsigned int i;
	qaws_curve_destroy(in->bspline3);
	qaws_curve_destroy(in->nurbs3);
	qaws_curve_destroy(in->bezier2);
	qaws_curve_destroy(in->loop2);
	for (i = 0; i < 16; i++)
	{
		qaws_curve_destroy(in->chain2[i]);
		qaws_exact_curve_destroy(in->exact2[i]);
	}
	qaws_surface_destroy(in->surface);
}

/* ---- curves ---- */

static void st_curve_eval(char const* name, qaws_curve const* c, unsigned int flags, unsigned int full)
{
	unsigned int n = st_count(full), i, ok = 1;
	qaws_range r = qaws_curve_get_parameter_range(c);
	unsigned long long h = 14695981039346656037ull;
	st_clock t0;
	if (!st_want(name))
		return;
	t0 = st_now();
	for (i = 0; i < n; i++)
	{
		qaws_eval_result_3d e;
		qaws_scalar t = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)((i * 0.6180339887498949) - floor(i * 0.6180339887498949));
		memset(&e, 0, sizeof e);
		ok &= qaws_curve_evaluate_3d(&e, c, t, flags) == QAWS_STATUS_OK;
		h = st_hash(h, &e.position, sizeof e.position * 3);
	}
	st_report(name, n, t0, h, ok);
}

static void st_curve_batch(char const* name, qaws_curve const* c, int two_d, unsigned int flags, unsigned int full)
{
	unsigned int n = st_count(full), i, k, ok = 1, block = 4096;
	qaws_range r = qaws_curve_get_parameter_range(c);
	qaws_scalar* t = (qaws_scalar*)malloc(sizeof(qaws_scalar) * block);
	qaws_eval_result_3d* e3 = (qaws_eval_result_3d*)calloc(block, sizeof *e3);
	qaws_eval_result_2d* e2 = (qaws_eval_result_2d*)calloc(block, sizeof *e2);
	unsigned long long h = 14695981039346656037ull;
	st_clock t0;
	if (!st_want(name) || !t || !e3 || !e2)
		goto done;
	for (i = 0; i < block; i++)
		t[i] = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / (qaws_scalar)(block - 1);
	t0 = st_now();
	for (k = 0; k < n; k += block)
	{
		unsigned int m = n - k < block ? n - k : block;
		if (two_d)
		{
			ok &= qaws_curve_evaluate_batch_2d(e2, c, t, m, flags) == QAWS_STATUS_OK;
			h = st_hash(h, e2, sizeof *e2 * m);
		}
		else
		{
			ok &= qaws_curve_evaluate_batch_3d(e3, c, t, m, flags) == QAWS_STATUS_OK;
			h = st_hash(h, e3, sizeof *e3 * m);
		}
	}
	st_report(name, n, t0, h, ok);
done:
	free(t);
	free(e3);
	free(e2);
}

static void st_curve_measures(st_inputs const* in)
{
	qaws_curve const* c = in->bspline3;
	qaws_range r = qaws_curve_get_parameter_range(c);
	unsigned long long h;
	unsigned int n, i, ok;
	st_clock t0;
	if (st_want("curve arc_length"))
	{
		n = st_count(5000), ok = 1, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_scalar len = 0, a = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)((i % 97) / 200.0);
			ok &= qaws_curve_compute_arc_length(&len, c, a, a + (r.max_value - r.min_value) * (qaws_scalar)0.5) == QAWS_STATUS_OK;
			h = st_hash(h, &len, sizeof len);
		}
		st_report("curve arc_length", n, t0, h, ok);
	}
	if (st_want("curve parameter_at_arc_length"))
	{
		qaws_scalar total = 0;
		qaws_curve_compute_arc_length(&total, c, r.min_value, r.max_value);
		n = st_count(2000), ok = 1, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_scalar t = 0;
			ok &= qaws_curve_compute_parameter_at_arc_length(&t, c, r.min_value, total * (qaws_scalar)((i % 1000) / 1000.0)) == QAWS_STATUS_OK;
			h = st_hash(h, &t, sizeof t);
		}
		st_report("curve parameter_at_arc_length", n, t0, h, ok);
	}
	if (st_want("curve closest_parameter"))
	{
		n = st_count(10000), ok = 1, h = 14695981039346656037ull, g_rng = 7, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_vec3 p;
			qaws_scalar t = 0;
			p.x = (qaws_scalar)(2.4 * st_rand() - 1.2), p.y = (qaws_scalar)(2.4 * st_rand() - 1.2), p.z = (qaws_scalar)(3.2 * st_rand());
			ok &= qaws_curve_find_closest_parameter_3d(&t, c, p) == QAWS_STATUS_OK;
			h = st_hash(h, &t, sizeof t);
		}
		st_report("curve closest_parameter", n, t0, h, ok);
	}
	if (st_want("curve bounds"))
	{
		n = st_count(2000), ok = 1, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_vec3 lo, hi;
			ok &= qaws_curve_compute_bounds_3d(&lo, &hi, (i & 1) ? in->nurbs3 : c) == QAWS_STATUS_OK;
			h = st_hash(h, &lo, sizeof lo), h = st_hash(h, &hi, sizeof hi);
		}
		st_report("curve bounds", n, t0, h, ok);
	}
	if (st_want("curve tessellate"))
	{
		qaws_scalar* tp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 65536);
		qaws_vec3* pp = (qaws_vec3*)malloc(sizeof(qaws_vec3) * 65536);
		qaws_curve_tessellation_desc td;
		memset(&td, 0, sizeof td);
		td.chord_tolerance = (qaws_scalar)1e-4, td.angle_tolerance = (qaws_scalar)0.05, td.max_depth = 16;
		n = st_count(500), ok = tp && pp, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; ok && i < n; i++)
		{
			unsigned int cnt = 0;
			ok &= qaws_curve_tessellate_3d(tp, pp, &cnt, (i & 1) ? in->nurbs3 : c, 65536, &td) == QAWS_STATUS_OK;
			h = st_hash(h, pp, sizeof *pp * cnt);
		}
		st_report("curve tessellate", n, t0, h, ok);
		free(tp);
		free(pp);
	}
	if (st_want("curve winding"))
	{
		n = st_count(200000), ok = 1, h = 14695981039346656037ull, g_rng = 11, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_vec2 p;
			int wn = 0;
			p.x = (qaws_scalar)(2.4 * st_rand() - 1.2), p.y = (qaws_scalar)(2.4 * st_rand() - 1.2);
			qaws_curve_compute_winding_number_2d(&wn, in->loop2, p);   /* (on the curve: refused, counted as 0) */
			h = st_hash(h, &wn, sizeof wn);
		}
		st_report("curve winding", n, t0, h, ok);
	}
	if (st_want("curve enclose"))
	{
		n = st_count(100000), ok = 1, h = 14695981039346656037ull, g_rng = 13, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			/* order 2: C, C' and C'' (3 x 3 values each way) */
			double lo[9], hi[9], a = r.min_value + (r.max_value - r.min_value) * st_rand() * 0.99;
			double b = fmin(a + 0.01 * (r.max_value - r.min_value), (double)r.max_value);
			ok &= qaws_curve_enclose(lo, hi, c, a, b, 2) == QAWS_STATUS_OK;
			h = st_hash(h, lo, sizeof lo), h = st_hash(h, hi, sizeof hi);
		}
		st_report("curve enclose", n, t0, h, ok);
	}
	if (st_want("curve offset"))
	{
		qaws_curve* out[256];
		n = st_count(50), ok = 1, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			unsigned int cnt = 0, k;
			ok &= qaws_curve_offset_2d(out, &cnt, in->loop2, 256, (qaws_scalar)0.05, 1) == QAWS_STATUS_OK;
			for (k = 0; k < cnt; k++)
			{
				qaws_eval_result_2d e;
				memset(&e, 0, sizeof e);
				qaws_curve_evaluate_2d(&e, out[k], qaws_curve_get_parameter_range(out[k]).min_value, QAWS_EVAL_FLAG_POSITION);
				h = st_hash(h, &e.position, sizeof e.position);
				qaws_curve_destroy(out[k]);
			}
		}
		st_report("curve offset", n, t0, h, ok);
	}
}

static void st_curve_set(st_inputs const* in)
{
	qaws_curve const* cs[18];
	qaws_curve_batch_desc d;
	qaws_curve_set* set = NULL;
	unsigned int n = st_count(100000), i, ok;
	qaws_scalar* pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	qaws_closest_point* out = (qaws_closest_point*)calloc(n, sizeof *out);
	unsigned long long h = 14695981039346656037ull;
	st_clock t0;
	if (!st_want("curve_set closest") || !pts || !out)
		goto done;
	cs[0] = in->bspline3, cs[1] = in->nurbs3;
	for (i = 0; i < 16; i++)
		cs[2 + i] = in->chain2[i];
	memset(&d, 0, sizeof d);
	d.curves = cs, d.curve_count = 2;   /* the 3D ones */
	g_rng = 17;
	for (i = 0; i < 3 * n; i++)
		pts[i] = (qaws_scalar)(i % 3 == 2 ? 3.2 * st_rand() : 2.4 * st_rand() - 1.2);
	t0 = st_now();
	ok = qaws_curve_set_create(&set, &d) == QAWS_STATUS_OK;
	ok = ok && qaws_curve_set_find_closest(out, NULL, set, pts, n, 0) == QAWS_STATUS_OK;
	for (i = 0; i < n; i++)
		h = st_hash(h, &out[i].parameter, sizeof out[i].parameter), h = st_hash(h, &out[i].curve, sizeof out[i].curve);
	st_report("curve_set closest", n, t0, h, ok);
	qaws_curve_set_destroy(set);
done:
	free(pts);
	free(out);
}

/* ---- surfaces ---- */

static void st_surfaces(st_inputs const* in)
{
	qaws_surface const* s = in->surface;
	unsigned long long h;
	unsigned int n, i, ok;
	st_clock t0;
	qaws_range ur, vr;
	ur = qaws_surface_get_u_range(s);
	vr = qaws_surface_get_v_range(s);
	if (st_want("surface eval"))
	{
		n = st_count(1000000), ok = 1, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_surface_eval_result e;
			double fu = i * 0.6180339887498949, fv = i * 0.7548776662466927;
			memset(&e, 0, sizeof e);
			ok &= qaws_surface_evaluate(&e, s, ur.min_value + (ur.max_value - ur.min_value) * (qaws_scalar)(fu - floor(fu)),
			                            vr.min_value + (vr.max_value - vr.min_value) * (qaws_scalar)(fv - floor(fv)),
			                            QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL)
			      == QAWS_STATUS_OK;
			h = st_hash(h, &e, sizeof(qaws_vec3) * 4);
		}
		st_report("surface eval", n, t0, h, ok);
	}
	if (st_want("surface eval_batch"))
	{
		unsigned int block = 4096, k;
		qaws_scalar* u = (qaws_scalar*)malloc(sizeof(qaws_scalar) * block);
		qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * block);
		qaws_surface_eval_result* e = (qaws_surface_eval_result*)calloc(block, sizeof *e);
		n = st_count(1000000), ok = u && v && e, h = 14695981039346656037ull;
		for (i = 0; ok && i < block; i++)
			u[i] = ur.min_value + (ur.max_value - ur.min_value) * (qaws_scalar)(i % 64) / 63,
			v[i] = vr.min_value + (vr.max_value - vr.min_value) * (qaws_scalar)(i / 64) / 63;
		t0 = st_now();
		for (k = 0; ok && k < n; k += block)
		{
			unsigned int m = n - k < block ? n - k : block;
			ok &= qaws_surface_evaluate_batch(e, s, u, v, m, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL)
			      == QAWS_STATUS_OK;
			h = st_hash(h, e, sizeof *e * m);
		}
		st_report("surface eval_batch", n, t0, h, ok);
		free(u);
		free(v);
		free(e);
	}
	if (st_want("surface closest"))
	{
		n = st_count(2000), ok = 1, h = 14695981039346656037ull, g_rng = 19, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_vec3 p, q;
			qaws_scalar u = 0, v = 0;
			p.x = (qaws_scalar)(15 * st_rand()), p.y = (qaws_scalar)(15 * st_rand()), p.z = (qaws_scalar)(2 * st_rand() - 1);
			ok &= qaws_surface_find_closest_point(&u, &v, &q, s, p) == QAWS_STATUS_OK;
			h = st_hash(h, &u, sizeof u), h = st_hash(h, &v, sizeof v);
		}
		st_report("surface closest", n, t0, h, ok);
	}
	if (st_want("surface tessellate"))
	{
		unsigned int vcap = 1u << 18, icap = 3u << 19;
		qaws_tessellation_vertex* vx = (qaws_tessellation_vertex*)malloc(sizeof *vx * vcap);
		unsigned int* ix = (unsigned int*)malloc(sizeof *ix * icap);
		qaws_tessellation_desc td;
		memset(&td, 0, sizeof td);
		td.max_depth = 8, td.curvature_threshold = (qaws_scalar)0.02, td.max_edge_length = (qaws_scalar)0.25;
		n = st_count(4), ok = vx && ix, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; ok && i < n; i++)
		{
			unsigned int nv = 0, ni = 0;
			ok &= qaws_surface_tessellate(vx, &nv, ix, &ni, s, vcap, icap, &td) == QAWS_STATUS_OK;
			h = st_hash(h, vx, sizeof *vx * nv), h = st_hash(h, ix, sizeof *ix * ni);
		}
		st_report("surface tessellate", n, t0, h, ok);
		free(vx);
		free(ix);
	}
	if (st_want("surface area"))
	{
		n = st_count(500), ok = 1, h = 14695981039346656037ull, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_scalar a = 0;
			ok &= qaws_surface_compute_area(&a, s) == QAWS_STATUS_OK;
			h = st_hash(h, &a, sizeof a);
		}
		st_report("surface area", n, t0, h, ok);
	}
	if (st_want("surface enclose"))
	{
		n = st_count(20000), ok = 1, h = 14695981039346656037ull, g_rng = 23, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			double lo[18], hi[18];   /* order 2: S, S_u, S_v, S_uu, S_uv, S_vv */
			double u0 = ur.min_value + (ur.max_value - ur.min_value) * st_rand() * 0.98, v0 = vr.min_value + (vr.max_value - vr.min_value) * st_rand() * 0.98;
			double du = 0.02 * (ur.max_value - ur.min_value), dv = 0.02 * (vr.max_value - vr.min_value);
			ok &= qaws_surface_enclose(lo, hi, s, u0, fmin(u0 + du, (double)ur.max_value), v0, fmin(v0 + dv, (double)vr.max_value), 2) == QAWS_STATUS_OK;
			h = st_hash(h, lo, sizeof lo), h = st_hash(h, hi, sizeof hi);
		}
		st_report("surface enclose", n, t0, h, ok);
	}
	if (st_want("surface_set closest"))
	{
		qaws_surface_batch_desc d;
		qaws_surface_set* set = NULL;
		qaws_scalar* pts;
		qaws_surface_batch_closest* out;
		n = st_count(50000);
		pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		out = (qaws_surface_batch_closest*)calloc(n, sizeof *out);
		ok = pts && out, h = 14695981039346656037ull, g_rng = 29;
		for (i = 0; ok && i < 3 * n; i++)
			pts[i] = (qaws_scalar)(i % 3 == 2 ? 2 * st_rand() - 1 : 15 * st_rand());
		memset(&d, 0, sizeof d);
		d.surfaces = &s, d.surface_count = 1;
		t0 = st_now();
		ok = ok && qaws_surface_set_create(&set, &d) == QAWS_STATUS_OK;
		ok = ok && qaws_surface_set_find_closest(out, NULL, set, pts, n, 0) == QAWS_STATUS_OK;
		for (i = 0; ok && i < n; i++)
			h = st_hash(h, &out[i].u, sizeof out[i].u), h = st_hash(h, &out[i].v, sizeof out[i].v);
		st_report("surface_set closest", n, t0, h, ok);
		qaws_surface_set_destroy(set);
		free(pts);
		free(out);
	}
}

/* ---- regions, polylines, exact ---- */

static void st_paths(st_inputs const* in)
{
	unsigned long long h;
	unsigned int n, i, ok;
	st_clock t0;
	if (st_want("region locate"))
	{
		qaws_path_2d path;
		path.curves = (qaws_curve const* const*)in->chain2, path.curve_count = 16, path.closed = 1;
		n = st_count(20000), ok = 1, h = 14695981039346656037ull, g_rng = 31, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			qaws_point_location loc = QAWS_POINT_OUTSIDE;
			qaws_vec2 p;
			p.x = (qaws_scalar)(2.4 * st_rand() - 1.2), p.y = (qaws_scalar)(2.4 * st_rand() - 1.2);
			ok &= qaws_region_locate_point_2d(&loc, &path, 1, QAWS_FILL_NON_ZERO, p, 0) == QAWS_STATUS_OK;
			h = st_hash(h, &loc, sizeof loc);
		}
		st_report("region locate", n, t0, h, ok);
	}
	if (st_want("region prepared"))
	{
		/* the same points through a prepared region: the same hash as "region locate" */
		qaws_path_2d path;
		qaws_region_2d* r = NULL;
		qaws_vec2* pts;
		qaws_point_location* loc;
		path.curves = (qaws_curve const* const*)in->chain2, path.curve_count = 16, path.closed = 1;
		n = st_count(20000);
		pts = (qaws_vec2*)malloc(sizeof *pts * n);
		loc = (qaws_point_location*)malloc(sizeof *loc * n);
		ok = pts && loc, h = 14695981039346656037ull, g_rng = 31;
		for (i = 0; ok && i < n; i++)
			pts[i].x = (qaws_scalar)(2.4 * st_rand() - 1.2), pts[i].y = (qaws_scalar)(2.4 * st_rand() - 1.2);
		t0 = st_now();
		ok = ok && qaws_region_create_2d(&r, &path, 1, QAWS_FILL_NON_ZERO, 0) == QAWS_STATUS_OK;
		ok = ok && qaws_region_locate_points_2d(loc, r, NULL, pts, n, NULL) == QAWS_STATUS_OK;
		for (i = 0; ok && i < n; i++)
			h = st_hash(h, &loc[i], sizeof loc[i]);
		st_report("region prepared", n, t0, h, ok);
		qaws_region_destroy_2d(r);
		free(pts);
		free(loc);
	}
	if (st_want("polyline simplify"))
	{
		unsigned int m = 20000;
		qaws_scalar* pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * m);
		qaws_scalar* out = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * m);
		n = st_count(2), ok = pts && out, h = 14695981039346656037ull, g_rng = 37;
		for (i = 0; ok && i < m; i++)
		{
			double a = 6.283185307179586 * i / m;
			pts[2 * i] = (qaws_scalar)((1 + 0.01 * st_rand()) * cos(a)), pts[2 * i + 1] = (qaws_scalar)((1 + 0.01 * st_rand()) * sin(a));
		}
		t0 = st_now();
		for (i = 0; ok && i < n; i++)
		{
			unsigned int cnt = 0;
			ok &= qaws_polyline_simplify_2d(out, &cnt, pts, m, 1, (qaws_scalar)0.005) == QAWS_STATUS_OK;
			h = st_hash(h, out, sizeof(qaws_scalar) * 2 * cnt);
		}
		st_report("polyline simplify", n * m, t0, h, ok);
		free(pts);
		free(out);
	}
	if (st_want("exact orient2d"))
	{
		n = st_count(2000000), ok = 1, h = 14695981039346656037ull, g_rng = 41, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			double a[2], b[2], c[2];
			qaws_exact_sign sg;
			a[0] = st_rand(), a[1] = st_rand(), b[0] = st_rand(), b[1] = st_rand();
			c[0] = a[0] + (b[0] - a[0]) * 0.5, c[1] = a[1] + (b[1] - a[1]) * 0.5;   /* near-collinear: the exact path */
			if (i & 1)
				c[0] = st_rand(), c[1] = st_rand();
			ok &= qaws_exact_orient2d(&sg, NULL, a, b, c) == QAWS_STATUS_OK;
			h = st_hash(h, &sg, sizeof sg);
		}
		st_report("exact orient2d", n, t0, h, ok);
	}
	if (st_want("exact orient3d"))
	{
		n = st_count(1000000), ok = 1, h = 14695981039346656037ull, g_rng = 43, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			double a[3], b[3], c[3], d[3];
			qaws_exact_sign sg;
			unsigned int k;
			for (k = 0; k < 3; k++)
				a[k] = st_rand(), b[k] = st_rand(), c[k] = st_rand();
			for (k = 0; k < 3; k++)
				d[k] = (i & 1) ? st_rand() : (a[k] + b[k] + c[k]) / 3;   /* coplanar half the time */
			ok &= qaws_exact_orient3d(&sg, NULL, a, b, c, d) == QAWS_STATUS_OK;
			h = st_hash(h, &sg, sizeof sg);
		}
		st_report("exact orient3d", n, t0, h, ok);
	}
	if (st_want("exact winding"))
	{
		n = st_count(50000), ok = 1, h = 14695981039346656037ull, g_rng = 47, t0 = st_now();
		for (i = 0; i < n; i++)
		{
			double p[2];
			int wn = 0;
			p[0] = 2.4 * st_rand() - 1.2, p[1] = 2.4 * st_rand() - 1.2;
			qaws_exact_winding_2d(&wn, (qaws_exact_curve const* const*)in->exact2, p, 16);
			h = st_hash(h, &wn, sizeof wn);
		}
		st_report("exact winding", n, t0, h, ok);
	}
	if (st_want("exact batch_closest"))
	{
		qaws_exact_batch_desc d;
		double* pts;
		qaws_exact_closest_point* out;
		n = st_count(20000);
		pts = (double*)malloc(sizeof(double) * 2 * n);
		out = (qaws_exact_closest_point*)calloc(n, sizeof *out);
		ok = pts && out, h = 14695981039346656037ull, g_rng = 53;
		for (i = 0; ok && i < 2 * n; i++)
			pts[i] = 2.4 * st_rand() - 1.2;
		memset(&d, 0, sizeof d);
		d.curves = (qaws_exact_curve const* const*)in->exact2, d.curve_count = 16;
		t0 = st_now();
		ok = ok && qaws_exact_curve_batch_closest(out, NULL, pts, n, &d) == QAWS_STATUS_OK;
		for (i = 0; ok && i < n; i++)
			h = st_hash(h, &out[i], sizeof out[i]);
		st_report("exact batch_closest", n, t0, h, ok);
		free(pts);
		free(out);
	}
}

int test_182_stress_main(void)
{
	st_inputs in;
	char const* s = getenv("QAWS_STRESS");
	unsigned int all = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;
	g_pass = 0;
	g_fail = 0;
	g_scale = s ? (unsigned int)atoi(s) : 0;
	g_only = getenv("QAWS_STRESS_ONLY");
	printf("Test 182: synthetic stress of the hot calls (%s, scale %u%s)\n", QAWS_SCALAR_IS_FLOAT ? "f32" : "f64", g_scale,
	       g_scale ? "" : ": smoke counts");
	st_make(&in);
	TEST_ASSERT(in.bspline3 && in.nurbs3 && in.bezier2 && in.loop2 && in.surface && in.exact2[0] && in.exact2[15], "the inputs made");
	if (in.bspline3 && in.nurbs3 && in.bezier2 && in.loop2 && in.surface && in.exact2[0] && in.exact2[15])
	{
		st_curve_eval("curve eval bspline", in.bspline3, all, 1000000);
		st_curve_eval("curve eval nurbs", in.nurbs3, all, 1000000);
		st_curve_batch("curve eval_batch bspline", in.bspline3, 0, all, 2000000);
		st_curve_batch("curve eval_batch nurbs", in.nurbs3, 0, all, 2000000);
		st_curve_batch("curve eval_batch bezier2d", in.bezier2, 1, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, 2000000);
		st_curve_measures(&in);
		st_curve_set(&in);
		st_surfaces(&in);
		st_paths(&in);
	}
	st_free(&in);
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
