#include "qaws_eval.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_span.h"
#include <string.h>

qaws_status qaws_curve_evaluate_2d(
	qaws_eval_result_2d* out_result,
	qaws_curve const* curve,
	qaws_scalar parameter,
	unsigned int eval_flags)
{
	qaws_scalar local_t;
	unsigned int span_index;

	if (!curve || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	span_index = qaws_internal_find_span(curve, parameter, &local_t);
	return curve->vtable->eval_span_2d(curve, span_index, local_t, eval_flags, out_result);
}

qaws_status qaws_curve_evaluate_3d(
	qaws_eval_result_3d* out_result,
	qaws_curve const* curve,
	qaws_scalar parameter,
	unsigned int eval_flags)
{
	qaws_scalar local_t;
	unsigned int span_index;

	if (!curve || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	span_index = qaws_internal_find_span(curve, parameter, &local_t);
	return curve->vtable->eval_span_3d(curve, span_index, local_t, eval_flags, out_result);
}

qaws_status qaws_curve_evaluate_span_2d(
	qaws_eval_result_2d* out_result,
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_parameter,
	unsigned int eval_flags)
{
	if (!curve || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (span_index >= curve->span_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	if (local_parameter < (qaws_scalar)0) local_parameter = (qaws_scalar)0;
	if (local_parameter > (qaws_scalar)1) local_parameter = (qaws_scalar)1;

	return curve->vtable->eval_span_2d(curve, span_index, local_parameter, eval_flags, out_result);
}

qaws_status qaws_curve_evaluate_span_3d(
	qaws_eval_result_3d* out_result,
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_parameter,
	unsigned int eval_flags)
{
	if (!curve || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (span_index >= curve->span_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	if (local_parameter < (qaws_scalar)0) local_parameter = (qaws_scalar)0;
	if (local_parameter > (qaws_scalar)1) local_parameter = (qaws_scalar)1;

	return curve->vtable->eval_span_3d(curve, span_index, local_parameter, eval_flags, out_result);
}

/* --- Batch evaluation --- */

qaws_status qaws_curve_evaluate_batch_2d(
	qaws_eval_result_2d* out_results,
	qaws_curve const* curve,
	qaws_scalar const* parameters,
	unsigned int count,
	unsigned int eval_flags)
{
	unsigned int i;
	qaws_status first_err = QAWS_STATUS_OK;

	if (!curve || !parameters || !out_results)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	for (i = 0; i < count; i++)
	{
		qaws_scalar local_t;
		unsigned int span_index;
		qaws_status s;

		memset(&out_results[i], 0, sizeof(out_results[i]));
		span_index = qaws_internal_find_span(curve, parameters[i], &local_t);
		s = curve->vtable->eval_span_2d(curve, span_index, local_t, eval_flags, &out_results[i]);
		if (s != QAWS_STATUS_OK && first_err == QAWS_STATUS_OK)
			first_err = s;
	}
	return first_err;
}

qaws_status qaws_curve_evaluate_batch_3d(
	qaws_eval_result_3d* out_results,
	qaws_curve const* curve,
	qaws_scalar const* parameters,
	unsigned int count,
	unsigned int eval_flags)
{
	unsigned int i;
	qaws_status first_err = QAWS_STATUS_OK;

	if (!curve || !parameters || !out_results)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	for (i = 0; i < count; i++)
	{
		qaws_scalar local_t;
		unsigned int span_index;
		qaws_status s;

		memset(&out_results[i], 0, sizeof(out_results[i]));
		span_index = qaws_internal_find_span(curve, parameters[i], &local_t);
		s = curve->vtable->eval_span_3d(curve, span_index, local_t, eval_flags, &out_results[i]);
		if (s != QAWS_STATUS_OK && first_err == QAWS_STATUS_OK)
			first_err = s;
	}
	return first_err;
}

qaws_status qaws_surface_evaluate_batch(
	qaws_surface_eval_result* out_results,
	qaws_surface const* surface,
	qaws_scalar const* u_params,
	qaws_scalar const* v_params,
	unsigned int count,
	unsigned int eval_flags)
{
	unsigned int i;
	qaws_status first_err = QAWS_STATUS_OK;

	if (!surface || !u_params || !v_params || !out_results)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (i = 0; i < count; i++)
	{
		qaws_status s;
		memset(&out_results[i], 0, sizeof(out_results[i]));
		s = qaws_surface_evaluate(&out_results[i], surface, u_params[i],
			v_params[i], eval_flags);
		if (s != QAWS_STATUS_OK && first_err == QAWS_STATUS_OK)
			first_err = s;
	}
	return first_err;
}
