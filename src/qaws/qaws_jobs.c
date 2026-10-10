/*
 * Many independent clip and offset jobs in one call, in parallel through an
 * executor: each job is solved on its own, results in order, the same as
 * one call each.
 */

#include "qaws_clip.h"
#include "qaws_clip64.h"
#include "qaws_offset.h"
#include "internal/qaws_internal_parallel.h"
#include <stdlib.h>
#include <string.h>

typedef enum jb_kind { JB_CLIP, JB_CLIP64, JB_OFFSET, JB_OFFSET64 } jb_kind;

typedef struct jb_ctx
{
	jb_kind kind;
	void const* descs;
	void** out;
} jb_ctx;

static qaws_status jb_job(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	jb_ctx const* j = (jb_ctx const*)ctx;
	unsigned int i;
	(void)chunk;
	for (i = begin; i < end; i++)
	{
		qaws_status s;
		switch (j->kind)
		{
		case JB_CLIP:
		{
			/* the job runs serially: the batch is the parallel part */
			qaws_clip_desc d = ((qaws_clip_desc const*)j->descs)[i];
			d.executor = NULL;
			s = qaws_clip_execute(&d, (qaws_clip_result**)&j->out[i]);
			break;
		}
		case JB_CLIP64:
			s = qaws_clip64_execute(&((qaws_clip64_desc const*)j->descs)[i], (qaws_clip64_result**)&j->out[i]);
			break;
		case JB_OFFSET:
		{
			qaws_offset_desc d = ((qaws_offset_desc const*)j->descs)[i];
			d.executor = NULL;
			s = qaws_offset_execute(&d, (qaws_clip_result**)&j->out[i]);
			break;
		}
		default:
			s = qaws_offset64_execute(&((qaws_offset64_desc const*)j->descs)[i], (qaws_clip64_result**)&j->out[i]);
			break;
		}
		if (s != QAWS_STATUS_OK) return s;
	}
	return QAWS_STATUS_OK;
}

static qaws_status jb_run(jb_kind kind, void const* descs, unsigned int count, qaws_batch_executor const* executor, void** out)
{
	jb_ctx c;
	qaws_status s;
	unsigned int i;
	if ((count && (!descs || !out)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++) out[i] = NULL;
	c.kind = kind;
	c.descs = descs;
	c.out = out;
	s = qaws_internal_parallel(executor, count, 1, jb_job, &c);
	if (s != QAWS_STATUS_OK)
	{
		/* all or nothing */
		for (i = 0; i < count; i++)
		{
			if (!out[i]) continue;
			if (kind == JB_CLIP || kind == JB_OFFSET) qaws_clip_result_destroy((qaws_clip_result*)out[i]);
			else qaws_clip64_result_destroy((qaws_clip64_result*)out[i]);
			out[i] = NULL;
		}
	}
	return s;
}

qaws_status qaws_clip_execute_batch(qaws_clip_desc const* descs, unsigned int count, qaws_batch_executor const* executor,
	qaws_clip_result** out_results)
{
	return jb_run(JB_CLIP, descs, count, executor, (void**)out_results);
}

qaws_status qaws_clip64_execute_batch(qaws_clip64_desc const* descs, unsigned int count, qaws_batch_executor const* executor,
	qaws_clip64_result** out_results)
{
	return jb_run(JB_CLIP64, descs, count, executor, (void**)out_results);
}

qaws_status qaws_offset_execute_batch(qaws_offset_desc const* descs, unsigned int count, qaws_batch_executor const* executor,
	qaws_clip_result** out_results)
{
	return jb_run(JB_OFFSET, descs, count, executor, (void**)out_results);
}

qaws_status qaws_offset64_execute_batch(qaws_offset64_desc const* descs, unsigned int count, qaws_batch_executor const* executor,
	qaws_clip64_result** out_results)
{
	return jb_run(JB_OFFSET64, descs, count, executor, (void**)out_results);
}
