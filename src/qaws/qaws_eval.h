#ifndef QAWS_EVAL_H
#define QAWS_EVAL_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"

/* All evaluation functions are thread-safe on immutable curves. */

qaws_status qaws_curve_evaluate_2d(
	qaws_eval_result_2d* out_result,
	qaws_curve const* curve,
	qaws_scalar parameter,
	unsigned int eval_flags);

qaws_status qaws_curve_evaluate_3d(
	qaws_eval_result_3d* out_result,
	qaws_curve const* curve,
	qaws_scalar parameter,
	unsigned int eval_flags);

qaws_status qaws_curve_evaluate_span_2d(
	qaws_eval_result_2d* out_result,
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_parameter,
	unsigned int eval_flags);

qaws_status qaws_curve_evaluate_span_3d(
	qaws_eval_result_3d* out_result,
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_parameter,
	unsigned int eval_flags);

/*
 * Batch evaluation: evaluate a curve at multiple parameters in one call.
 * parameters: array of count parameter values.
 * out_results: array of count result structs (caller-allocated).
 * Returns QAWS_STATUS_OK if all evaluations succeed.
 * On partial failure, evaluates as many as possible and returns the first
 * error status; out_results[i].valid_flags == 0 indicates that entry failed.
 */
qaws_status qaws_curve_evaluate_batch_2d(
	qaws_eval_result_2d* out_results,
	qaws_curve const* curve,
	qaws_scalar const* parameters,
	unsigned int count,
	unsigned int eval_flags);

qaws_status qaws_curve_evaluate_batch_3d(
	qaws_eval_result_3d* out_results,
	qaws_curve const* curve,
	qaws_scalar const* parameters,
	unsigned int count,
	unsigned int eval_flags);

/*
 * Batch surface evaluation: evaluate a surface at multiple (u,v) pairs.
 * u_params, v_params: arrays of count parameter values.
 * out_results: array of count result structs (caller-allocated).
 */
qaws_status qaws_surface_evaluate_batch(
	qaws_surface_eval_result* out_results,
	qaws_surface const* surface,
	qaws_scalar const* u_params,
	qaws_scalar const* v_params,
	unsigned int count,
	unsigned int eval_flags);

#endif /* QAWS_EVAL_H */
