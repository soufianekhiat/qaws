# Threading

qaws creates no threads and keeps no worker pool. Its only lock is a short one in the memory layer. It is designed to
run inside a host's own threading: objects are immutable once created, so any number of threads may read them, and
batched operations can hand their independent chunks to the host's scheduler. API page:
[batch_executor](api/batch_executor.md).

## What may be shared

| Object | Concurrent use | Needs external synchronization |
|---|---|---|
| `qaws_curve`, `qaws_surface` (owned, borrowed, inline after re-pointing) | evaluate, sample, inspect, measure, query, derive new curves | create and destroy; never destroy one in use |
| prepared curve / surface (`qaws_curve_prepared`, ...) | evaluate, including disjoint ranges of one range call | the host must not change the arrays the object depends on during an evaluation |
| `qaws_curve_set`, `qaws_surface_set` | queries | destroy |
| seams, shatter, clip result | reads (counts, pieces, paths, vertices) | destroy |
| `qaws_motion` | reads | destroy |
| `qaws_traversal` | evaluate and map queries | `qaws_traversal_advance_*` and `qaws_traversal_reset` change its state |
| `qaws_brep_shell` | queries, validation, measures | `qaws_brep_add_*` |
| differentiation adjoint views | | adjoint functions write into caller storage, so threads sharing one adjoint view must synchronize ([differentiation.md](differentiation.md)) |
| `qaws_arena` | | one thread at a time |

Functions that create a new object from shared inputs (operations, conversions, shatter, sets) only read those inputs.
Several threads may run them at once on the same curves and surfaces.

## No hidden mutable state

Curves and surfaces hold no lazily filled caches. Everything an object needs is computed when it is created, and
evaluation writes only into the caller's output and stack scratch. The only process-wide state lives in the memory
layer (see [memory.md](memory.md)):

- **Default allocator.** It is write-once in effect. Set it before the library first allocates. After that, reading it
  is a plain atomic load.
- **Thread allocators.** The thread allocator and the thread scratch allocator are thread-local. A setting on one
  thread never affects another.
- **Kept memory.** The pool of kept blocks is shared by all threads and guarded by a short spin lock taken only when a
  block is taken or returned. `qaws_set_kept_memory_limit` and `qaws_release_kept_memory` are thread-safe.
- **Allocator calls.** Each block records the allocator that made it, so a block freed on another thread still goes
  back to its own allocator. Your allocator callbacks must be thread-safe if objects made with them are created or
  destroyed on several threads.

## The batch executor

Batched operations take an optional `qaws_batch_executor`. This covers curve and surface closest points, prepared sets,
and the certified batches of the exact kernel.

```c
typedef void (*qaws_batch_task_fn)(void* ctx, unsigned int begin, unsigned int end);

typedef struct qaws_batch_executor
{
	void (*parallel_for)(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx);
	void* user;
} qaws_batch_executor;
```

Your executor and the batch each have obligations:

- **Your executor.** `parallel_for` must call `task(ctx, begin, end)` over ranges that cover `[0, count)` exactly once.
  It may use any threads, in any order, and must return only when every range is done. If a chunk is never run, the
  batch returns `QAWS_STATUS_INTERNAL_ERROR`.
- **Chunking.** The batch splits its independent work (query points, rays, cells) into at most 256 chunks. Each chunk
  has its own scratch, hits and statistics, and these are merged in chunk order on the calling thread. Results are
  therefore identical whatever the executor does, including running everything serially.
- **Calling thread.** Grid building and the final merge stay on the calling thread. Without an executor, all work runs
  there.
- **Allocation in tasks.** Allocations made inside tasks follow the allocator rules of the thread that runs the task.
  If you install per-worker arenas, install them on the executor's threads.

An executor can be adapted from any task system, such as a thread pool, OpenMP, TBB or a job system. qaws needs
nothing else from it.

## Cancellation and progress

The public float API takes no cancel or progress callbacks. Long calls run to completion.

- **Bounding the work.** To limit how long a call takes, size the work: use coarser flatness or resolution, use fewer
  query points per call, or split one batch into several calls.
- **Streaming sampling.** `qaws_curve_sample_streaming_*` calls back for each sample. A return value of 0 stops early,
  and stopping early is not an error.
- **`QAWS_STATUS_CANCELLED`.** This status exists in `qaws_status` for operations that are stopped by a host's cancel
  callback. Its outputs are as that operation documents. The library's internal parallel helper supports such a stop:
  it checks between chunks or waves, on the calling thread, never on the executor's threads. This lets libraries built
  on qaws report a cancellation with the shared status.

## Practical rules

- Create objects on any thread, then share them read-only. Destroy an object only after every reader is done.
- Give each thread its own traversal, or synchronize its `advance` and `reset` calls.
- Copy an inline curve per thread, and call `qaws_curve_inline_get_curve` on the copy before sharing it.
- Keep arenas per thread.
- Set the default allocator once at start-up, before any other qaws call.
