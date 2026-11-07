#include "qaws_internal_parallel.h"
#include "qaws_internal_memory.h"
#include <stdlib.h>

#define PAR_MAX_CHUNKS 256

typedef struct par_wrap
{
	unsigned int count, chunks;
	qaws_status (*job)(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end);
	void* ctx;
	qaws_status* status;
	unsigned int base;              /* the wave's first chunk (waves: par_wave_task) */
} par_wrap;

unsigned int qaws_internal_chunk_count(unsigned int count, unsigned int grain)
{
	unsigned int k;
	if (count == 0)
		return 0;
	if (grain == 0)
		grain = 1;
	k = (count + grain - 1) / grain;
	return k > PAR_MAX_CHUNKS ? PAR_MAX_CHUNKS : k;
}

static void par_range(par_wrap const* w, unsigned int k, unsigned int* b, unsigned int* e)
{
	*b = (unsigned int)((unsigned long long)w->count * k / w->chunks);
	*e = (unsigned int)((unsigned long long)w->count * (k + 1) / w->chunks);
}

/* the executor's task: chunks [begin, end) */
static void par_task(void* ctx, unsigned int begin, unsigned int end)
{
	par_wrap* w = (par_wrap*)ctx;
	unsigned int k;
	for (k = begin; k < end && k < w->chunks; k++)
	{
		unsigned int b, e;
		par_range(w, k, &b, &e);
		w->status[k] = w->job(w->ctx, k, b, e);
	}
}

/* a wave's task: its chunks [begin, end) from the wave's base */
static void par_wave_task(void* ctx, unsigned int begin, unsigned int end)
{
	par_wrap* w = (par_wrap*)ctx;
	par_task(ctx, w->base + begin, w->base + end);
}

qaws_status qaws_internal_parallel(unsigned int count, unsigned int grain, qaws_status (*job)(void* user, unsigned int chunk, unsigned int begin, unsigned int end),
	void* user, qaws_batch_executor const* executor)
{
	return qaws_internal_parallel_stop(count, grain, job, user, executor, NULL, NULL);
}

qaws_status qaws_internal_parallel_stop(unsigned int count, unsigned int grain,
	qaws_status (*job)(void* user, unsigned int chunk, unsigned int begin, unsigned int end), void* user, qaws_batch_executor const* executor,
	int (*stop)(void* stop_user), void* stop_user)
{
	par_wrap w;
	unsigned int k;
	qaws_status s = QAWS_STATUS_OK;
	w.count = count;
	w.chunks = qaws_internal_chunk_count(count, grain);
	w.job = job;
	w.ctx = user;
	w.base = 0;
	if (w.chunks == 0)
		return QAWS_STATUS_OK;
	w.status = (qaws_status*)qaws_internal_malloc(qaws_internal_mul2(w.chunks, sizeof(qaws_status)));
	if (!w.status)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (k = 0; k < w.chunks; k++)
		w.status[k] = QAWS_STATUS_INTERNAL_ERROR;   /* a chunk the executor never ran */
	if (stop && stop(stop_user))
	{
		qaws_internal_free(w.status);
		return QAWS_STATUS_CANCELLED;
	}
	if (executor && executor->parallel_for && stop)
	{
		/* the chunks in waves, stop called between them on this thread (never on the executor's): at least 32 chunks a
		   wave (its threads kept busy), at most 16 waves */
		unsigned int wave = w.chunks < 32 ? w.chunks : 32;
		if (wave < (w.chunks + 15) / 16)
			wave = (w.chunks + 15) / 16;
		for (k = 0; k < w.chunks; k += wave)
		{
			unsigned int nk = w.chunks - k < wave ? w.chunks - k : wave;
			if (k && stop(stop_user))
			{
				qaws_internal_free(w.status);
				return QAWS_STATUS_CANCELLED;
			}
			executor->parallel_for(executor->user, nk, par_wave_task, &w);
			w.base += nk;
		}
	}
	else if (executor && executor->parallel_for)
		executor->parallel_for(executor->user, w.chunks, par_task, &w);
	else if (stop)
		for (k = 0; k < w.chunks; k++)
		{
			if (k && stop(stop_user))
			{
				qaws_internal_free(w.status);
				return QAWS_STATUS_CANCELLED;
			}
			par_task(&w, k, k + 1);
		}
	else
		par_task(&w, 0, w.chunks);
	for (k = 0; k < w.chunks && s == QAWS_STATUS_OK; k++)
		s = w.status[k];
	qaws_internal_free(w.status);
	return s;
}
