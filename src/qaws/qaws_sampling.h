#ifndef QAWS_SAMPLING_H
#define QAWS_SAMPLING_H

#include "qaws_types.h"
#include "qaws_status.h"

/* All sampling functions are thread-safe on immutable curves. */

qaws_status qaws_curve_get_sample_count(
	unsigned int* out_sample_count,
	qaws_curve const* curve,
	qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_2d(
	qaws_vec2* out_positions,
	unsigned int* out_position_count,
	qaws_curve const* curve,
	unsigned int position_capacity,
	qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_3d(
	qaws_vec3* out_positions,
	unsigned int* out_position_count,
	qaws_curve const* curve,
	unsigned int position_capacity,
	qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_eval_2d(
	qaws_eval_result_2d* out_samples,
	unsigned int* out_sample_count,
	qaws_curve const* curve,
	unsigned int sample_capacity,
	qaws_sampling_desc const* desc);

qaws_status qaws_curve_sample_eval_3d(
	qaws_eval_result_3d* out_samples,
	unsigned int* out_sample_count,
	qaws_curve const* curve,
	unsigned int sample_capacity,
	qaws_sampling_desc const* desc);

/* Curvature-weighted sampling: distribute samples proportionally to
   curvature magnitude so high-curvature regions receive more points. */

qaws_status qaws_curve_sample_curvature_weighted_2d(
	qaws_vec2* out_positions,
	unsigned int* out_position_count,
	qaws_curve const* curve,
	unsigned int sample_count,
	unsigned int position_capacity);

qaws_status qaws_curve_sample_curvature_weighted_3d(
	qaws_vec3* out_positions,
	unsigned int* out_position_count,
	qaws_curve const* curve,
	unsigned int sample_count,
	unsigned int position_capacity);

/* Feature-preserving sampling: inserts inflection points and extrema
   into a uniform parameter sample set. */

qaws_status qaws_curve_sample_feature_preserving_2d(
	qaws_vec2* out_positions,
	unsigned int* out_position_count,
	qaws_curve const* curve,
	unsigned int base_sample_count,
	unsigned int position_capacity);

qaws_status qaws_curve_sample_feature_preserving_3d(
	qaws_vec3* out_positions,
	unsigned int* out_position_count,
	qaws_curve const* curve,
	unsigned int base_sample_count,
	unsigned int position_capacity);

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

/* Adaptive tessellation into caller buffers: each span is halved until every
   segment is within the chord tolerance of the curve (measured at 1/4, 1/2
   and 3/4 of it) and the tangent turns by at most the angle tolerance along
   it (summed over those points). Span boundaries are always vertices. The
   subdivision is deterministic: call with out_parameters NULL for the count,
   then with buffers of that size. */
typedef struct qaws_curve_tessellation_desc
{
	qaws_scalar chord_tolerance;        /* distance; 0: unused */
	qaws_scalar angle_tolerance;        /* radians; 0: unused (one of the two is needed) */
	unsigned int max_depth;             /* halvings per span at most; 0: 16, at most 30 */
} qaws_curve_tessellation_desc;

/* out_points may be NULL. BUFFER_TOO_SMALL when capacity < the count (in
   *out_count); the buffers then hold the first capacity vertices. */
qaws_status qaws_curve_tessellate_2d(
	qaws_curve const* curve,
	qaws_curve_tessellation_desc const* desc,
	qaws_scalar* out_parameters,
	qaws_vec2* out_points,
	unsigned int capacity,
	unsigned int* out_count);

qaws_status qaws_curve_tessellate_3d(
	qaws_curve const* curve,
	qaws_curve_tessellation_desc const* desc,
	qaws_scalar* out_parameters,
	qaws_vec3* out_points,
	unsigned int capacity,
	unsigned int* out_count);

#endif /* QAWS_SAMPLING_H */
