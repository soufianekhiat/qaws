#ifndef QAWS_PROGRESS_H
#define QAWS_PROGRESS_H

/*
 * Cancellation and progress reports for qaws's long calls (the batches, the
 * prepared sets, the exact batches, seams and shatters). Precision-free: one
 * type for both precisions.
 *
 * A long call takes a qaws_progress (a desc's `progress` member, or an _ex
 * form's argument); NULL or a zeroed one costs one untaken branch per check.
 * The call invokes cancel and report on the thread that called it, never on
 * an executor's threads (an executor's work is split into waves, the checks
 * made between them), so the host side may be as simple as an atomic read
 * and a few stored numbers.
 *
 *   cancel   polled between chunks of work (with an executor: between its
 *            waves). Non-zero stops the call: QAWS_STATUS_CANCELLED,
 *            everything it allocated freed, its inputs untouched, its
 *            outputs as each header says (by default undefined but freeable:
 *            nothing the caller must release). A time budget is a cancel
 *            reading a clock.
 *   report   at most every ~50 ms, plus once at the start and the end of
 *            each stage (the end's done equal to its total; none after a
 *            cancel or a failure). done never decreases within a stage;
 *            total may grow once known, never below done. The report is
 *            valid only during the call.
 *
 * Prepared objects (sets) do not keep the progress they were created with:
 * their queries take their own.
 */

typedef struct qaws_progress_report
{
	char const* stage;               /* the stage's name ("flatten", "query", ...); static */
	unsigned int stage_index;        /* its place among the call's stages, from 0 */
	unsigned int stage_count;
	unsigned int iteration;          /* iterative calls: the iteration under way, from 0; else 0 */
	unsigned int max_iterations;     /* iterative calls: their budget; else 0 */
	unsigned long long done;         /* the stage's work items done */
	unsigned long long total;        /* the stage's work items known */
	unsigned long long bytes_held;   /* the memory the call holds, when known; else 0 */
} qaws_progress_report;

typedef struct qaws_progress
{
	int (*cancel)(void* user_data);
	void (*report)(qaws_progress_report const* report, void* user_data);
	void* user_data;
	int keep_partial;                /* iterative calls that offer it: a cancel keeps the last complete iteration */
} qaws_progress;

#endif /* QAWS_PROGRESS_H */
