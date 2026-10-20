/*
 * Test 141: evaluation from many host threads
 *
 *   - closest points on a bumpy surface with their derivatives
 *     (qaws_surface_closest_point_tangent, whose global check kept its
 *     sample grid in a static array): four threads each solving every query
 *     several times give the serial answers bit for bit.
 */

#include "test_common.h"
#include "test_threads.h"
#include "qaws_surface_bezier.h"
#include "qaws_surface.h"
#include "qaws_diff_ops.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define TH_QUERIES 48
#define TH_THREADS 4
#define TH_ROUNDS 6

typedef struct th_job
{
	qaws_surface const* s;
	qaws_vec3 const* q;
	qaws_surface_closest_point const* serial;
	unsigned int mismatches, failures;
} th_job;

static void th_solve(void* p)
{
	th_job* j = (th_job*)p;
	unsigned int r, i;
	for (r = 0; r < TH_ROUNDS; r++)
		for (i = 0; i < TH_QUERIES; i++)
		{
			qaws_surface_closest_point v;
			if (qaws_surface_closest_point_tangent(&v, NULL, j->s, NULL, NULL, j->q[i], NULL) != QAWS_STATUS_OK)
				j->failures++;
			else if (memcmp(&v, &j->serial[i], sizeof v) != 0)
				j->mismatches++;
		}
}

static void test_surface_closest_points(void)
{
	qaws_vec3 grid[36], q[TH_QUERIES];
	qaws_surface_closest_point serial[TH_QUERIES];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	th_job job[TH_THREADS];
	void* args[TH_THREADS];
	unsigned int i, j, bad = 0, mism = 0, fail = 0;
	char msg[300];
	/* a quintic patch with bumps, so the global check finds other basins */
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
		{
			grid[i * 6 + j].x = (qaws_scalar)i / 5;
			grid[i * 6 + j].y = (qaws_scalar)j / 5;
			grid[i * 6 + j].z = (qaws_scalar)(0.6 * sin(2.3 * i + 1.1) * cos(1.9 * j + 0.5));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 5;
	d.v_degree = 5;
	d.control_points = grid;
	d.u_point_count = 6;
	d.v_point_count = 6;
	if (qaws_surface_create_bezier(&s, &d) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "surface: created");
		return;
	}
	for (i = 0; i < TH_QUERIES; i++)
	{
		q[i].x = (qaws_scalar)(0.5 + 0.6 * sin(0.7 * i));
		q[i].y = (qaws_scalar)(0.5 + 0.6 * cos(1.3 * i));
		q[i].z = (qaws_scalar)(0.8 * sin(2.1 * i));
		bad += qaws_surface_closest_point_tangent(&serial[i], NULL, s, NULL, NULL, q[i], NULL) != QAWS_STATUS_OK;
	}
	for (i = 0; i < TH_THREADS; i++)
	{
		job[i].s = s;
		job[i].q = q;
		job[i].serial = serial;
		job[i].mismatches = job[i].failures = 0;
		args[i] = &job[i];
	}
	test_run_threads(TH_THREADS, th_solve, args);
	for (i = 0; i < TH_THREADS; i++)
	{
		mism += job[i].mismatches;
		fail += job[i].failures;
	}
	sprintf(msg, "surface closest points on %u threads, %u queries %u times each: %u differ from the serial answers, %u failed (serial %u)",
		TH_THREADS, TH_QUERIES, TH_ROUNDS, mism, fail, bad);
	TEST_ASSERT(bad == 0 && mism == 0 && fail == 0, msg);
	qaws_surface_destroy(s);
}

int test_141_eval_threads_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 141: evaluation from many host threads\n");
	test_surface_closest_points();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
