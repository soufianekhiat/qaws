#include "qaws_internal_parallel.h"
#include <stdlib.h>

#define PAR_MAX_CHUNKS 256

typedef struct par_wrap
{
	unsigned int count, chunks;
	qaws_status (*job)(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end);
	void* ctx;
	qaws_status* status;
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

qaws_status qaws_internal_parallel(qaws_batch_executor const* executor, unsigned int count, unsigned int grain,
	qaws_status (*job)(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end), void* ctx)
{
	par_wrap w;
	unsigned int k;
	qaws_status s = QAWS_STATUS_OK;
	w.count = count;
	w.chunks = qaws_internal_chunk_count(count, grain);
	w.job = job;
	w.ctx = ctx;
	if (w.chunks == 0)
		return QAWS_STATUS_OK;
	w.status = (qaws_status*)malloc(w.chunks * sizeof(qaws_status));
	if (!w.status)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (k = 0; k < w.chunks; k++)
		w.status[k] = QAWS_STATUS_INTERNAL_ERROR;   /* a chunk the executor never ran */
	if (executor && executor->parallel_for)
		executor->parallel_for(executor->user, w.chunks, par_task, &w);
	else
		par_task(&w, 0, w.chunks);
	for (k = 0; k < w.chunks && s == QAWS_STATUS_OK; k++)
		s = w.status[k];
	free(w.status);
	return s;
}
