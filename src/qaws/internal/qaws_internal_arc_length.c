#include "qaws_internal_arc_length.h"
#include "qaws_internal_types.h"
#include "qaws_internal_span.h"
#include <math.h>

/* 5-point Gauss-Legendre quadrature nodes and weights on [-1, 1] */
static const qaws_scalar gl_nodes[5] = {
	(qaws_scalar)-0.906179845938664,
	(qaws_scalar)-0.538469310105683,
	(qaws_scalar)0.0,
	(qaws_scalar)0.538469310105683,
	(qaws_scalar)0.906179845938664
};

static const qaws_scalar gl_weights[5] = {
	(qaws_scalar)0.236926885056189,
	(qaws_scalar)0.478628670499366,
	(qaws_scalar)0.568888888888889,
	(qaws_scalar)0.478628670499366,
	(qaws_scalar)0.236926885056189
};

static qaws_scalar arc_length_segment(
	qaws_curve const* curve,
	qaws_scalar t_min,
	qaws_scalar t_max)
{
	qaws_scalar half_len = (t_max - t_min) * (qaws_scalar)0.5;
	qaws_scalar mid = (t_min + t_max) * (qaws_scalar)0.5;
	qaws_scalar sum = (qaws_scalar)0;
	unsigned int dim = (unsigned int)curve->dimension;

	for (unsigned int i = 0; i < 5; i++) {
		qaws_scalar t = mid + half_len * gl_nodes[i];

		qaws_scalar local_t;
		unsigned int span_index = qaws_internal_find_span(&local_t, curve, t);

		qaws_scalar speed = (qaws_scalar)0;

		if (dim == 2) {
			qaws_eval_result_2d res;
			curve->vtable->eval_span_2d(
				curve, span_index, local_t,
				QAWS_EVAL_FLAG_D1, &res);
			speed = (qaws_scalar)sqrt(
				(double)(res.d1.x * res.d1.x + res.d1.y * res.d1.y));
		} else if (dim == 3) {
			qaws_eval_result_3d res;
			curve->vtable->eval_span_3d(
				curve, span_index, local_t,
				QAWS_EVAL_FLAG_D1, &res);
			speed = (qaws_scalar)sqrt(
				(double)(res.d1.x * res.d1.x + res.d1.y * res.d1.y + res.d1.z * res.d1.z));
		}

		sum += gl_weights[i] * speed;
	}

	return half_len * sum;
}

qaws_status qaws_internal_arc_length(
	qaws_scalar* out_length,
	qaws_curve const* curve,
	qaws_scalar t_min,
	qaws_scalar t_max)
{
	if (!curve || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (t_min >= t_max) {
		*out_length = (qaws_scalar)0;
		return QAWS_STATUS_OK;
	}

	/* Use span boundaries for subdivision to avoid integrating across discontinuities */
	qaws_scalar total = (qaws_scalar)0;
	qaws_scalar seg_start = t_min;

	for (unsigned int s = 0; s < curve->span_count; s++) {
		qaws_scalar boundary = curve->span_boundaries[s + 1];
		if (boundary <= t_min) continue;
		if (boundary >= t_max) break;

		total += arc_length_segment(curve, seg_start, boundary);
		seg_start = boundary;
	}

	total += arc_length_segment(curve, seg_start, t_max);

	*out_length = total;
	return QAWS_STATUS_OK;
}

qaws_status qaws_internal_build_arc_length_table(
	qaws_scalar* table_params,
	qaws_scalar* table_distances,
	qaws_curve const* curve,
	unsigned int table_size)
{
	if (!curve || !table_params || !table_distances || table_size < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	qaws_scalar t_min = curve->parameter_range.min_value;
	qaws_scalar t_max = curve->parameter_range.max_value;

	for (unsigned int i = 0; i < table_size; i++) {
		table_params[i] = t_min + (t_max - t_min) *
			(qaws_scalar)i / (qaws_scalar)(table_size - 1);
	}

	table_distances[0] = (qaws_scalar)0;

	for (unsigned int i = 1; i < table_size; i++) {
		qaws_scalar seg_len;
		qaws_internal_arc_length(&seg_len, curve, table_params[i - 1], table_params[i]);
		table_distances[i] = table_distances[i - 1] + seg_len;
	}

	return QAWS_STATUS_OK;
}

qaws_scalar qaws_internal_distance_to_parameter(
	qaws_scalar const* table_params,
	qaws_scalar const* table_distances,
	unsigned int table_size,
	qaws_scalar distance)
{
	if (distance <= table_distances[0])
		return table_params[0];
	if (distance >= table_distances[table_size - 1])
		return table_params[table_size - 1];

	unsigned int low = 0;
	unsigned int high = table_size - 1;

	while (high - low > 1) {
		unsigned int mid = (low + high) / 2;
		if (table_distances[mid] <= distance) {
			low = mid;
		} else {
			high = mid;
		}
	}

	qaws_scalar d_range = table_distances[high] - table_distances[low];
	if (d_range <= (qaws_scalar)0)
		return table_params[low];

	qaws_scalar alpha = (distance - table_distances[low]) / d_range;
	return table_params[low] + alpha * (table_params[high] - table_params[low]);
}

qaws_scalar qaws_internal_parameter_to_distance(
	qaws_scalar const* table_params,
	qaws_scalar const* table_distances,
	unsigned int table_size,
	qaws_scalar parameter)
{
	if (parameter <= table_params[0])
		return table_distances[0];
	if (parameter >= table_params[table_size - 1])
		return table_distances[table_size - 1];

	unsigned int low = 0;
	unsigned int high = table_size - 1;

	while (high - low > 1) {
		unsigned int mid = (low + high) / 2;
		if (table_params[mid] <= parameter) {
			low = mid;
		} else {
			high = mid;
		}
	}

	qaws_scalar t_range = table_params[high] - table_params[low];
	if (t_range <= (qaws_scalar)0)
		return table_distances[low];

	qaws_scalar alpha = (parameter - table_params[low]) / t_range;
	return table_distances[low] + alpha * (table_distances[high] - table_distances[low]);
}
