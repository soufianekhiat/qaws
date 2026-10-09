#ifndef QAWS_BATCH_EXECUTOR_H
#define QAWS_BATCH_EXECUTOR_H

/*
 * Parallel execution of the batched operations (float and exact): qaws has
 * no threads of its own. A batch splits its independent work (query
 * points, rays, grid cells) into chunks it owns and hands them to the
 * caller's executor: parallel_for(user, count, task, ctx) must call
 * task(ctx, begin, end) over ranges covering [0, count) exactly once, on
 * any threads in any order, and return when all are done. The results do
 * not depend on the executor. Without one the work runs on the calling
 * thread.
 */
typedef void (*qaws_batch_task_fn)(void* ctx, unsigned int begin, unsigned int end);

typedef struct qaws_batch_executor
{
	void (*parallel_for)(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx);
	void* user;
} qaws_batch_executor;

#endif /* QAWS_BATCH_EXECUTOR_H */
