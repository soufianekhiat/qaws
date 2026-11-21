# API Reference: qaws_batch_executor.h

Parallel execution of the batched operations, float and exact. qaws has no
threads of its own: a batch splits its independent work (query points,
rays, grid cells) into chunks it owns and hands them to the caller's
executor. Without an executor the work runs on the calling thread.

This header is precision-free: the names have no `_f32` / `_f64` forms.
See [threading.md](../threading.md) for the threading model.

---

## Types

```c
typedef void (*qaws_batch_task_fn)(void* ctx, unsigned int begin, unsigned int end);

typedef struct qaws_batch_executor {
    void (*parallel_for)(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx);
    void* user;
} qaws_batch_executor;
```

- `parallel_for(user, count, task, ctx)` must call `task(ctx, begin, end)`
  over ranges covering `[0, count)` exactly once, on any threads in any
  order, and return when all are done. `user` is the executor's own
  pointer; `ctx` belongs to the batch and is passed back unchanged.
- `count` is a number of chunks, at most 256. The executor may hand a task
  one chunk or a range of them.
- A chunk the executor never runs makes the batch return
  `QAWS_STATUS_INTERNAL_ERROR`.
- Each chunk has its own scratch, hit lists and statistics, merged in
  chunk order afterwards, so the results do not depend on the executor.
- A NULL executor, or one with a NULL `parallel_for`, runs the chunks on
  the calling thread.
- With cancellation, the chunks go to the executor in waves (at least 32
  chunks a wave, at most 16 waves), and the stop check runs between waves
  on the calling thread, never on the executor's threads.

## Where it is taken

An optional `qaws_batch_executor const* executor` member of:

- the curve batch descriptors, [curve_batch.md](curve_batch.md)
  (`qaws_curve_batch_desc`, `qaws_level_crossing_desc`,
  `qaws_closest_desc`); prepared curve sets keep the one of their desc;
- the surface batch descriptors, [surface_batch.md](surface_batch.md);
- the clipping descriptor, for its batched intersection,
  [clip.md](clip.md);
- the exact kernel's batched descriptors, [exact.md](exact.md).

## Example

```c
static void run_all(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
    my_pool* pool = (my_pool*)user;
    my_pool_for(pool, count, task, ctx);   /* calls task(ctx, b, e) over [0, count), waits */
}

qaws_batch_executor ex = { run_all, &pool };
desc.executor = &ex;
```
