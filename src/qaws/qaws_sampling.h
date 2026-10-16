#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SAMPLING_H_F32
#define QAWS_SAMPLING_H_F32
#define QAWS_SAMPLING_H_BODY
#endif
#else
#ifndef QAWS_SAMPLING_H_F64
#define QAWS_SAMPLING_H_F64
#define QAWS_SAMPLING_H_BODY
#endif
#endif
#ifdef QAWS_SAMPLING_H_BODY
#undef QAWS_SAMPLING_H_BODY
#define QAWS_SAMPLING_H

#include "qaws_types.h"
#include "qaws_status.h"

/* All sampling functions are thread-safe on immutable curves. */

qaws_status qaws_curve_get_sample_count(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	unsigned int* out_sample_count);

qaws_status qaws_curve_sample_2d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_vec2* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

qaws_status qaws_curve_sample_3d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_vec3* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

qaws_status qaws_curve_sample_eval_2d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_eval_result_2d* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_sample_count);

qaws_status qaws_curve_sample_eval_3d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_eval_result_3d* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_sample_count);

/* Curvature-weighted sampling: distribute samples proportionally to
   curvature magnitude so high-curvature regions receive more points. */

qaws_status qaws_curve_sample_curvature_weighted_2d(
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_vec2* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

qaws_status qaws_curve_sample_curvature_weighted_3d(
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_vec3* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

/* Feature-preserving sampling: inserts inflection points and extrema
   into a uniform parameter sample set. */

qaws_status qaws_curve_sample_feature_preserving_2d(
	qaws_curve const* curve,
	unsigned int base_sample_count,
	qaws_vec2* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

qaws_status qaws_curve_sample_feature_preserving_3d(
	qaws_curve const* curve,
	unsigned int base_sample_count,
	qaws_vec3* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

/* Streaming / callback-based sampling.
   The callback is invoked for each sample point.  Return non-zero from
   the callback to continue, or 0 to stop early (not an error). */

typedef int (*qaws_sample_callback_2d)(
	qaws_scalar parameter,
	qaws_vec2 const* position,
	void* user_data);

typedef int (*qaws_sample_callback_3d)(
	qaws_scalar parameter,
	qaws_vec3 const* position,
	void* user_data);

qaws_status qaws_curve_sample_streaming_2d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_sample_callback_2d callback,
	void* user_data);

qaws_status qaws_curve_sample_streaming_3d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_sample_callback_3d callback,
	void* user_data);

#endif /* QAWS_SAMPLING_H */
