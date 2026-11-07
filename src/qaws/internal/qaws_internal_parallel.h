#ifndef QAWS_INTERNAL_PARALLEL_H
#define QAWS_INTERNAL_PARALLEL_H

#include "../qaws_curve_batch.h"

/*
 * Work split into chunks the batch owns: [0, count) in chunk_count ranges of
 * about `grain` items. job(user, chunk, begin, end) runs once per chunk,
 * through the executor when one is given (any order, any thread), serially
 * otherwise; each chunk has its own index for its scratch and statistics.
 * Returns the first failing chunk's status.
 */
unsigned int qaws_internal_chunk_count(unsigned int count, unsigned int grain);

qaws_status qaws_internal_parallel(unsigned int count, unsigned int grain, qaws_status (*job)(void* user, unsigned int chunk, unsigned int begin, unsigned int end),
	void* user, qaws_batch_executor const* executor);
/* The same, stop (may be NULL) called with stop_user on this thread: before each chunk when the chunks run serially (no
   executor); through one, between waves of them (at least 32 chunks a wave, at most 16 waves; never on the executor's
   threads). Non-zero ends it with QAWS_STATUS_CANCELLED, the chunks after left undone */
qaws_status qaws_internal_parallel_stop(unsigned int count, unsigned int grain,
	qaws_status (*job)(void* user, unsigned int chunk, unsigned int begin, unsigned int end), void* user, qaws_batch_executor const* executor,
	int (*stop)(void* stop_user), void* stop_user);

#endif /* QAWS_INTERNAL_PARALLEL_H */
