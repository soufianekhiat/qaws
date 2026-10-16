/*
 * Test 95: many independent jobs in one call (qaws_clip_execute_batch,
 * qaws_clip64_execute_batch, qaws_offset_execute_batch,
 * qaws_offset64_execute_batch): the results equal one call each, bit for
 * bit, serially and on executors running the jobs in other orders.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define J95_JOBS 48

static unsigned int j95_seed = 12345;
static double j95_rand(void)
{
	j95_seed = j95_seed * 1103515245u + 12345u;
	return (double)((j95_seed >> 8) & 0xFFFFFF) / 16777216.0;
}

/* executors: backwards, and odd chunks before even ones */
static void j95_backwards(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int i;
	(void)user;
	for (i = count; i-- > 0;) task(ctx, i, i + 1);
}

static void j95_interleaved(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int i;
	(void)user;
	for (i = 1; i < count; i += 2) task(ctx, i, i + 1);
	for (i = 0; i < count; i += 2) task(ctx, i, i + 1);
}

/* two float results the same, vertex for vertex */
static int j95_same(qaws_clip_result const* a, qaws_clip_result const* b)
{
	unsigned int i, k;
	if (!a || !b || qaws_clip_result_get_path_count(a) != qaws_clip_result_get_path_count(b)) return 0;
	for (i = 0; i < qaws_clip_result_get_path_count(a); i++)
	{
		qaws_clip_vertex const *va = NULL, *vb = NULL;
		unsigned int na = 0, nb = 0;
		qaws_clip_result_get_vertices(&va, &na, a, 0, i);
		qaws_clip_result_get_vertices(&vb, &nb, b, 0, i);
		if (na != nb) return 0;
		for (k = 0; k < na; k++)
			if (va[k].position.x != vb[k].position.x || va[k].position.y != vb[k].position.y) return 0;
	}
	return 1;
}

static int j95_same64(qaws_clip64_result const* a, qaws_clip64_result const* b)
{
	unsigned int i;
	if (!a || !b || qaws_clip64_result_get_path_count(a) != qaws_clip64_result_get_path_count(b)) return 0;
	for (i = 0; i < qaws_clip64_result_get_path_count(a); i++)
	{
		int64_t const *pa, *pb;
		unsigned int na = 0, nb = 0;
		qaws_clip64_result_get_path(&pa, &na, a, i);
		qaws_clip64_result_get_path(&pb, &nb, b, i);
		if (na != nb || memcmp(pa, pb, sizeof(int64_t) * 2 * na)) return 0;
	}
	return 1;
}

int test_95_jobs_main(void)
{
	qaws_curve* curves[2 * J95_JOBS];
	qaws_curve const* views[2 * J95_JOBS];
	qaws_path_2d paths[2 * J95_JOBS];
	int64_t pts[2 * J95_JOBS][2 * 12];
	qaws_path64 p64[2 * J95_JOBS];
	qaws_clip_desc cd[J95_JOBS];
	qaws_clip64_desc cd64[J95_JOBS];
	qaws_offset_group og[J95_JOBS];
	qaws_offset_desc od[J95_JOBS];
	qaws_offset64_group og64[J95_JOBS];
	qaws_offset64_desc od64[J95_JOBS];
	qaws_batch_executor ex[3];
	unsigned int i, e;
	char msg[300];
	g_pass = 0;
	g_fail = 0;
	printf("Test 95: batch jobs\n");
	memset(ex, 0, sizeof(ex));
	ex[1].parallel_for = j95_backwards;
	ex[2].parallel_for = j95_interleaved;
	/* jobs: pairs of random star polygons (float as polylines and circles, int64) */
	for (i = 0; i < 2 * J95_JOBS; i++)
	{
		qaws_scalar xy[24];
		unsigned int k, n = 5 + (unsigned int)(j95_rand() * 7);
		double cx = 10 * j95_rand(), cy = 10 * j95_rand();
		for (k = 0; k < n; k++)
		{
			double a = 2 * 3.14159265358979 * k / n, r = 2 + 3 * j95_rand();
			xy[2 * k] = (qaws_scalar)(cx + r * cos(a));
			xy[2 * k + 1] = (qaws_scalar)(cy + r * sin(a));
			pts[i][2 * k] = (int64_t)floor(1000 * (cx + r * cos(a)));
			pts[i][2 * k + 1] = (int64_t)floor(1000 * (cy + r * sin(a)));
		}
		curves[i] = NULL;
		if (i % 5 == 4)
			qaws_curve_create_ellipse_2d(&curves[i], (qaws_scalar)cx, (qaws_scalar)cy, (qaws_scalar)(2 + 2 * j95_rand()), (qaws_scalar)(1 + j95_rand()), 0);
		else
			qaws_curve_create_polyline_2d(&curves[i], xy, n, 1);
		views[i] = curves[i];
		paths[i].curves = &views[i]; paths[i].curve_count = 1; paths[i].closed = 1;
		p64[i].points = pts[i]; p64[i].point_count = n;
	}
	for (i = 0; i < J95_JOBS; i++)
	{
		static qaws_clip_type const ct[4] = { QAWS_CLIP_INTERSECTION, QAWS_CLIP_UNION, QAWS_CLIP_DIFFERENCE, QAWS_CLIP_XOR };
		memset(&cd[i], 0, sizeof(cd[i]));
		cd[i].subjects = &paths[2 * i]; cd[i].subject_count = 1; cd[i].clips = &paths[2 * i + 1]; cd[i].clip_count = 1;
		cd[i].clip_type = ct[i % 4]; cd[i].fill_rule = QAWS_FILL_NON_ZERO;
		memset(&cd64[i], 0, sizeof(cd64[i]));
		cd64[i].subjects = &p64[2 * i]; cd64[i].subject_count = 1; cd64[i].clips = &p64[2 * i + 1]; cd64[i].clip_count = 1;
		cd64[i].clip_type = ct[i % 4]; cd64[i].fill_rule = QAWS_FILL_NON_ZERO;
		og[i].paths = &paths[2 * i]; og[i].path_count = 2; og[i].join_type = (qaws_join_type)(i % 4); og[i].end_type = QAWS_END_POLYGON;
		memset(&od[i], 0, sizeof(od[i]));
		od[i].groups = &og[i]; od[i].group_count = 1; od[i].delta = (qaws_scalar)(i % 2 ? 0.4 : -0.3);
		og64[i].paths = &p64[2 * i]; og64[i].path_count = 2; og64[i].join_type = (qaws_join_type)(i % 4); og64[i].end_type = QAWS_END_POLYGON;
		memset(&od64[i], 0, sizeof(od64[i]));
		od64[i].groups = &og64[i]; od64[i].group_count = 1; od64[i].delta = i % 2 ? 400 : -300;
	}
	/* one call each, then batches on three executors */
	{
		qaws_clip_result *r1[J95_JOBS], *r2[J95_JOBS], *o1[J95_JOBS], *o2[J95_JOBS];
		qaws_clip64_result *s1[J95_JOBS], *s2[J95_JOBS], *q1[J95_JOBS], *q2[J95_JOBS];
		for (i = 0; i < J95_JOBS; i++)
		{
			r1[i] = NULL; qaws_clip_execute(&r1[i], &cd[i]);
			s1[i] = NULL; qaws_clip64_execute(&s1[i], &cd64[i]);
			o1[i] = NULL; qaws_offset_execute(&o1[i], &od[i]);
			q1[i] = NULL; qaws_offset64_execute(&q1[i], &od64[i]);
		}
		for (e = 0; e < 3; e++)
		{
			unsigned int same_c = 1, same_c64 = 1, same_o = 1, same_o64 = 1;
			qaws_batch_executor const* x = e ? &ex[e] : NULL;
			qaws_status a = qaws_clip_execute_batch(r2, J95_JOBS, cd, x);
			qaws_status b = qaws_clip64_execute_batch(s2, J95_JOBS, cd64, x);
			qaws_status c = qaws_offset_execute_batch(o2, J95_JOBS, od, x);
			qaws_status d = qaws_offset64_execute_batch(q2, J95_JOBS, od64, x);
			for (i = 0; i < J95_JOBS; i++)
			{
				same_c = same_c && j95_same(r1[i], r2[i]);
				same_c64 = same_c64 && j95_same64(s1[i], s2[i]);
				same_o = same_o && j95_same(o1[i], o2[i]);
				same_o64 = same_o64 && j95_same64(q1[i], q2[i]);
				qaws_clip_result_destroy(r2[i]); qaws_clip64_result_destroy(s2[i]);
				qaws_clip_result_destroy(o2[i]); qaws_clip64_result_destroy(q2[i]);
			}
			sprintf(msg, "%u jobs of each kind, %s: clip %d, clip64 %d, offset %d, offset64 %d the same as one call each",
				J95_JOBS, e == 0 ? "serial" : (e == 1 ? "executor backwards" : "executor interleaved"), same_c, same_c64, same_o, same_o64);
			TEST_ASSERT(a == QAWS_STATUS_OK && b == QAWS_STATUS_OK && c == QAWS_STATUS_OK && d == QAWS_STATUS_OK && same_c && same_c64 && same_o && same_o64, msg);
		}
		for (i = 0; i < J95_JOBS; i++)
		{
			qaws_clip_result_destroy(r1[i]); qaws_clip64_result_destroy(s1[i]);
			qaws_clip_result_destroy(o1[i]); qaws_clip64_result_destroy(q1[i]);
		}
	}
	/* a failing job: nothing kept */
	{
		qaws_clip_result* r[3];
		qaws_clip_desc bad[3];
		qaws_status s;
		bad[0] = cd[0]; bad[1] = cd[1]; bad[2] = cd[2];
		bad[1].subjects = NULL; bad[1].subject_count = 3;
		s = qaws_clip_execute_batch(r, 3, bad, NULL);
		sprintf(msg, "a failing job: its status, no result kept (%d, %p %p %p)", (int)s, (void*)r[0], (void*)r[1], (void*)r[2]);
		TEST_ASSERT(s != QAWS_STATUS_OK && !r[0] && !r[1] && !r[2], msg);
	}
	for (i = 0; i < 2 * J95_JOBS; i++) qaws_curve_destroy(curves[i]);
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
