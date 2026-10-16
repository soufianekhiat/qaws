#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_TRAVERSAL_H_F32
#define QAWS_TRAVERSAL_H_F32
#define QAWS_TRAVERSAL_H_BODY
#endif
#else
#ifndef QAWS_TRAVERSAL_H_F64
#define QAWS_TRAVERSAL_H_F64
#define QAWS_TRAVERSAL_H_BODY
#endif
#endif
#ifdef QAWS_TRAVERSAL_H_BODY
#undef QAWS_TRAVERSAL_H_BODY
#define QAWS_TRAVERSAL_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Traversal evaluate and map functions are thread-safe for concurrent reads.
 * qaws_traversal_advance_* and qaws_traversal_reset mutate state and require
 * external synchronization when used from multiple threads.
 */

qaws_status qaws_traversal_create(
	qaws_curve const* curve,
	qaws_traversal_desc const* desc,
	qaws_traversal** out_traversal);

qaws_status qaws_traversal_create_multi(
	qaws_curve const* const* curves,
	unsigned int curve_count,
	qaws_traversal_desc const* desc,
	qaws_traversal** out_traversal);

void qaws_traversal_destroy(qaws_traversal* traversal);

qaws_status qaws_traversal_evaluate_2d(
	qaws_traversal const* traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_2d* out_result);

qaws_status qaws_traversal_evaluate_3d(
	qaws_traversal const* traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_3d* out_result);

qaws_status qaws_traversal_map_time_to_parameter(
	qaws_traversal const* traversal,
	qaws_scalar time_value,
	qaws_scalar* out_parameter);

qaws_status qaws_traversal_map_distance_to_parameter(
	qaws_traversal const* traversal,
	qaws_scalar distance_value,
	qaws_scalar* out_parameter);

qaws_status qaws_traversal_map_parameter_to_distance(
	qaws_traversal const* traversal,
	qaws_scalar parameter,
	qaws_scalar* out_distance);

qaws_status qaws_traversal_advance_2d(
	qaws_traversal* traversal,
	qaws_scalar delta_time,
	unsigned int eval_flags,
	qaws_eval_result_2d* out_result);

qaws_status qaws_traversal_advance_3d(
	qaws_traversal* traversal,
	qaws_scalar delta_time,
	unsigned int eval_flags,
	qaws_eval_result_3d* out_result);

qaws_status qaws_traversal_reset(
	qaws_traversal* traversal);

#endif /* QAWS_TRAVERSAL_H */
