/*
 * A long call's side of qaws_progress (qaws_progress.h): its stages by name,
 * the work done in the current one, cancel polled and reports throttled. All
 * on the calling thread. A call without progress pays one untaken branch per
 * QAWS_PG_STEP.
 *
 *   qaws_pg_init(g, progress, names, count)   before the work
 *   qaws_pg_iteration(g, it, max)             an iterative call's next iteration (its stages start over)
 *   qaws_pg_stage(g, index, total)            the current stage ends, stage `index` starts
 *   QAWS_PG_STEP(g, done)                     done items in the stage: cancel every 64th call,
 *                                             a report once ~50 ms passed
 *   qaws_pg_chunk(g, done)                    as QAWS_PG_STEP for a step worth a poll of its own
 *                                             (a chunk of work): cancel polled on every call
 *   qaws_pg_end(g)                            the current stage's end, when no other follows
 *   qaws_pg_stop(g)                           as QAWS_PG_STEP without done, as a stop callback (void* g, may be NULL)
 *   qaws_pg_stop_each(g)                      the same polling on every call (a costly step each)
 *
 * The stage, iteration and step calls return non-zero once cancelled (and go
 * on returning it): the caller unwinds with QAWS_STATUS_CANCELLED.
 * qaws_internal_parallel_pg (qaws_internal_parallel.h) runs chunks with the
 * stage's done counting the items of the chunks finished.
 */

#ifndef QAWS_INTERNAL_PROGRESS_H
#define QAWS_INTERNAL_PROGRESS_H

#include "../qaws_progress.h"
#include <stddef.h>
#include <time.h>

typedef struct qaws_pg
{
	qaws_progress const* p;          /* NULL: none (or neither callback set) */
	char const* const* names;
	unsigned int count, index, iteration, max_iterations, calls;
	unsigned long long done, total, reported;
	long long last_ns;
	int started, cancelled;
} qaws_pg;

static inline long long qaws_pg_now_ns(void)
{
	struct timespec ts;
	if (timespec_get(&ts, TIME_UTC) != TIME_UTC)
		return 0;
	return (long long)ts.tv_sec * 1000000000LL + ts.tv_nsec;
}

static inline void qaws_pg_report(qaws_pg* g)
{
	qaws_progress_report r;
	if (!g->p->report || !g->started)
		return;
	if (g->done < g->reported)
		g->done = g->reported;
	g->reported = g->done;
	r.stage = g->index < g->count ? g->names[g->index] : "";
	r.stage_index = g->index;
	r.stage_count = g->count;
	r.iteration = g->iteration;
	r.max_iterations = g->max_iterations;
	r.done = g->done;
	r.total = g->total > g->done ? g->total : g->done;
	r.bytes_held = 0;
	g->last_ns = qaws_pg_now_ns();
	g->p->report(&r, g->p->user_data);
}

static inline int qaws_pg_poll(qaws_pg* g)
{
	if (!g->cancelled && g->p->cancel && g->p->cancel(g->p->user_data))
		g->cancelled = 1;
	return g->cancelled;
}

static inline void qaws_pg_init(qaws_pg* g, qaws_progress const* p, char const* const* names, unsigned int count)
{
	g->p = p && (p->cancel || p->report) ? p : NULL;
	g->names = names;
	g->count = count;
	g->index = g->iteration = g->max_iterations = g->calls = 0;
	g->done = g->total = g->reported = 0;
	g->last_ns = 0;
	g->started = g->cancelled = 0;
}

/* (the current stage's end: done raised to its total) */
static inline void qaws_pg_finish_(qaws_pg* g)
{
	if (!g->started)
		return;
	if (g->done < g->total)
		g->done = g->total;
	g->total = g->done;
	qaws_pg_report(g);
}

static inline int qaws_pg_stage(qaws_pg* g, unsigned int index, unsigned long long total)
{
	if (!g->p)
		return 0;
	if (g->cancelled || qaws_pg_poll(g))
		return 1;
	qaws_pg_finish_(g);
	g->index = index;
	g->done = g->reported = 0;
	g->total = total;
	g->started = 1;
	qaws_pg_report(g);
	return 0;
}

static inline int qaws_pg_iteration(qaws_pg* g, unsigned int iteration, unsigned int max_iterations)
{
	if (!g->p)
		return 0;
	if (g->cancelled || qaws_pg_poll(g))
		return 1;
	/* (the iteration before ends first: its last stage's end reported under it) */
	qaws_pg_finish_(g);
	g->started = 0;
	g->iteration = iteration;
	g->max_iterations = max_iterations;
	return 0;
}

static inline int qaws_pg_step_(qaws_pg* g, unsigned long long done)
{
	if (g->cancelled)
		return 1;
	if (done > g->done)
		g->done = done;
	if ((++g->calls & 63u) != 0)
		return 0;
	if (qaws_pg_poll(g))
		return 1;
	if (g->p->report && qaws_pg_now_ns() - g->last_ns >= 50000000LL)
		qaws_pg_report(g);
	return 0;
}

#define QAWS_PG_STEP(g, d) ((g)->p && qaws_pg_step_((g), (d)))

static inline int qaws_pg_chunk(qaws_pg* g, unsigned long long done)
{
	if (!g->p)
		return 0;
	g->calls |= 63u;   /* (the step's count then polls) */
	return qaws_pg_step_(g, done);
}

static inline int qaws_pg_stop(void* g)
{
	qaws_pg* p = (qaws_pg*)g;
	return p && p->p && qaws_pg_step_(p, p->done);
}

static inline int qaws_pg_stop_each(void* g)
{
	qaws_pg* p = (qaws_pg*)g;
	return p && qaws_pg_chunk(p, p->done);
}

static inline void qaws_pg_end(qaws_pg* g)
{
	if (g->p && !g->cancelled)
	{
		qaws_pg_finish_(g);
		g->started = 0;
	}
}

#endif
