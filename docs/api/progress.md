# API Reference: qaws_progress.h

Cancellation and progress reports for qaws's long calls. The header is precision-free: one type for both precisions,
no `_f32` / `_f64` forms. It is included by `qaws_batch_executor.h`, so every desc with an executor sees it.

See [threading.md](../threading.md) for the threading model.

---

## Types

```c
typedef struct qaws_progress_report {
    char const* stage;               /* the stage's name; static */
    unsigned int stage_index;        /* from 0 */
    unsigned int stage_count;
    unsigned int iteration;          /* iterative calls; else 0 */
    unsigned int max_iterations;     /* iterative calls; else 0 */
    unsigned long long done;         /* the stage's work items done */
    unsigned long long total;        /* the stage's work items known */
    unsigned long long bytes_held;   /* when known; else 0 */
} qaws_progress_report;

typedef struct qaws_progress {
    int (*cancel)(void* user_data);
    void (*report)(qaws_progress_report const* report, void* user_data);
    void* user_data;
    int keep_partial;                /* iterative calls that offer it */
} qaws_progress;
```

## Contract

- **Optional.** A NULL progress, or one with both callbacks NULL, costs one untaken branch per check.
- **Calling thread.** `cancel` and `report` run on the thread that made the call, never on an executor's threads. With
  an executor the work goes out in waves (at least 32 chunks a wave, at most 16 waves) and the checks run between them.
  The host side can be an atomic read and a few stores, with no locking.
- **cancel.** Polled between chunks of work. Non-zero stops the call: it returns `QAWS_STATUS_CANCELLED`, frees what it
  allocated, and leaves its inputs untouched. A create that is cancelled leaves its out pointer NULL. Filled buffers
  (the closest points) hold unspecified values. A time budget is a cancel that reads a clock.
- **report.** At the start and the end of each stage, and at most every ~50 ms in between. A stage's last report has
  `done == total`. `done` never decreases within a stage; `total` may grow once known, never below `done`. No report
  follows a cancel or a failure. The report is valid only during the callback.
- **Prepared sets** do not keep the progress they were created with: their queries take their own (the `_ex` forms).
- **Latency.** Cancel is checked between chunks (serial) or between waves (executor). A huge job split into 16 waves is
  checked every sixteenth of its run.

## Where it is taken

| Call | Where | Stages |
|---|---|---|
| `qaws_curve_batch_find_closest` | `qaws_closest_desc.progress` | `flatten` (curves), `query` (points) |
| `qaws_curve_set_create` | `qaws_curve_batch_desc.progress` | `flatten` |
| `qaws_curve_set_find_closest_ex` | argument | `query` |
| `qaws_surface_batch_find_closest` | `qaws_surface_batch_closest_desc.progress` | `flatten` (surfaces), `query` (points) |
| `qaws_surface_set_create` | `qaws_surface_batch_desc.progress` | `flatten` |
| `qaws_surface_set_find_closest_ex` | argument | `query` |
| `qaws_exact_curve_batch_closest` | `qaws_exact_batch_desc.progress` | `query` |
| `qaws_exact_surface_batch_closest` | `qaws_exact_ssi_batch_desc.progress` | `query` |
| `qaws_surface_seams_create` | `qaws_surface_seam_desc.progress` | `tessellate` (surfaces), `pairs` (triangle pairs tested) |
| `qaws_surface_shatter_create` | `qaws_surface_seam_desc.progress` | `tessellate`, `pairs`, `pieces` (A's surfaces) |
| `qaws_region_seams_2d` | `qaws_region_seam_desc.progress` | `flatten`, `pairs` (chord pairs tested) |
| `qaws_region_shatter_create_2d` | `qaws_region_seam_desc.progress` | `flatten`, `arrange` |

A `pairs` stage has no total known in advance: its `total` follows `done`.

## Example

```c
static int my_cancel(void* user) { return atomic_load(&((my_job*)user)->stop); }
static void my_report(qaws_progress_report const* r, void* user)
{
    my_job* j = (my_job*)user;
    j->stage = r->stage_index;
    j->fraction = r->total ? (double)r->done / (double)r->total : 0;
}

qaws_progress pg = { my_cancel, my_report, &job, 0 };
desc.progress = &pg;
st = qaws_curve_batch_find_closest(out, NULL, &desc);   /* QAWS_STATUS_CANCELLED once job.stop is set */
```
