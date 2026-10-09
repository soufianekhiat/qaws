#ifndef QAWS_INTERNAL_PARALLEL_H
#define QAWS_INTERNAL_PARALLEL_H

#include "../qaws_curve_batch.h"

/*
 * Work split into chunks the batch owns: [0, count) in chunk_count ranges of
 * about `grain` items. job(ctx, chunk, begin, end) runs once per chunk,
 * through the executor when one is given (any order, any thread), serially
 * otherwise; each chunk has its own index for its scratch and statistics.
 * Returns the first failing chunk's status.
 */
unsigned int qaws_internal_chunk_count(unsigned int count, unsigned int grain);

qaws_status qaws_internal_parallel(qaws_batch_executor const* executor, unsigned int count, unsigned int grain,
	qaws_status (*job)(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end), void* ctx);

#endif /* QAWS_INTERNAL_PARALLEL_H */
