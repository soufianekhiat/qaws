#include "qaws_sampling.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_cache.h"
#include <stdlib.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Adaptive sampling helpers
 * ------------------------------------------------------------------------- */

#define ADAPTIVE_MAX_DEPTH 10

struct adaptive_buf_2d
{
	qaws_vec2 *points;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_buf_2d_push(struct adaptive_buf_2d *buf, qaws_vec2 pt)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_vec2 *new_ptr = (qaws_vec2 *)qaws_internal_realloc(
			buf->points, qaws_internal_mul2(new_cap, sizeof(qaws_vec2)));
		if (!new_ptr)
			return 0;
		buf->points = new_ptr;
		buf->capacity = new_cap;
	}
	buf->points[buf->count++] = pt;
	return 1;
}

static qaws_status adaptive_subdivide_2d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec2 p0,
	qaws_scalar t1, qaws_vec2 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_buf_2d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar mx, my, dx, dy, dist;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		if (!adaptive_buf_2d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_2d(
		&result, curve, tmid, QAWS_EVAL_FLAG_POSITION);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Distance from midpoint of chord to actual curve point */
	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));

	if (dist > tolerance)
	{
		status = adaptive_subdivide_2d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_subdivide_2d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		if (!adaptive_buf_2d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

struct adaptive_buf_3d
{
	qaws_vec3 *points;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_buf_3d_push(struct adaptive_buf_3d *buf, qaws_vec3 pt)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_vec3 *new_ptr = (qaws_vec3 *)qaws_internal_realloc(
			buf->points, qaws_internal_mul2(new_cap, sizeof(qaws_vec3)));
		if (!new_ptr)
			return 0;
		buf->points = new_ptr;
		buf->capacity = new_cap;
	}
	buf->points[buf->count++] = pt;
	return 1;
}

static qaws_status adaptive_subdivide_3d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec3 p0,
	qaws_scalar t1, qaws_vec3 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_buf_3d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar mx, my, mz, dx, dy, dz, dist;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		if (!adaptive_buf_3d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_3d(
		&result, curve, tmid, QAWS_EVAL_FLAG_POSITION);
	if (status != QAWS_STATUS_OK)
		return status;

	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	mz = (p0.z + p1.z) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dz = result.position.z - mz;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy + dz * dz));

	if (dist > tolerance)
	{
		status = adaptive_subdivide_3d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_subdivide_3d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		if (!adaptive_buf_3d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_2d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_vec2 *out_positions,
	unsigned int position_capacity,
	unsigned int *out_position_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_2d result;
	qaws_status status;
	struct adaptive_buf_2d buf;
	qaws_scalar t0, t1;
	qaws_vec2 p0, p1;
	unsigned int copy_count;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.points = (qaws_vec2 *)qaws_internal_malloc(qaws_internal_mul2(buf.capacity, sizeof(qaws_vec2)));
	if (!buf.points)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Evaluate first point */
	status = qaws_curve_evaluate_2d(
		&result, curve, range_min, QAWS_EVAL_FLAG_POSITION);
	if (status != QAWS_STATUS_OK)
	{
		qaws_internal_free(buf.points);
		return status;
	}
	p0 = result.position;
	if (!adaptive_buf_2d_push(&buf, p0))
	{
		qaws_internal_free(buf.points);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_2d(
			&result, curve, t1, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.points);
			return status;
		}
		p1 = result.position;

		status = adaptive_subdivide_2d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.points);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > position_capacity)
	{
		qaws_internal_free(buf.points);
		*out_position_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_positions[i] = buf.points[i];

	*out_position_count = copy_count;
	qaws_internal_free(buf.points);
	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_3d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_vec3 *out_positions,
	unsigned int position_capacity,
	unsigned int *out_position_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_3d result;
	qaws_status status;
	struct adaptive_buf_3d buf;
	qaws_scalar t0, t1;
	qaws_vec3 p0, p1;
	unsigned int copy_count;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.points = (qaws_vec3 *)qaws_internal_malloc(qaws_internal_mul2(buf.capacity, sizeof(qaws_vec3)));
	if (!buf.points)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	status = qaws_curve_evaluate_3d(
		&result, curve, range_min, QAWS_EVAL_FLAG_POSITION);
	if (status != QAWS_STATUS_OK)
	{
		qaws_internal_free(buf.points);
		return status;
	}
	p0 = result.position;
	if (!adaptive_buf_3d_push(&buf, p0))
	{
		qaws_internal_free(buf.points);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_3d(
			&result, curve, t1, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.points);
			return status;
		}
		p1 = result.position;

		status = adaptive_subdivide_3d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.points);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > position_capacity)
	{
		qaws_internal_free(buf.points);
		*out_position_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_positions[i] = buf.points[i];

	*out_position_count = copy_count;
	qaws_internal_free(buf.points);
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Adaptive eval sampling helpers
 * ------------------------------------------------------------------------- */

struct adaptive_eval_buf_2d
{
	qaws_eval_result_2d *samples;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_eval_buf_2d_push(struct adaptive_eval_buf_2d *buf,
	qaws_eval_result_2d sample)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_eval_result_2d *new_ptr = (qaws_eval_result_2d *)qaws_internal_realloc(
			buf->samples, qaws_internal_mul2(new_cap, sizeof(qaws_eval_result_2d)));
		if (!new_ptr)
			return 0;
		buf->samples = new_ptr;
		buf->capacity = new_cap;
	}
	buf->samples[buf->count++] = sample;
	return 1;
}

static qaws_status adaptive_eval_subdivide_2d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec2 p0,
	qaws_scalar t1, qaws_vec2 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_eval_buf_2d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_2d result;
	qaws_eval_result_2d end_result;
	qaws_status status;
	qaws_scalar mx, my, dx, dy, dist;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		status = qaws_curve_evaluate_2d(
			&end_result, curve, t1, eval_flags);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_2d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_2d(
		&result, curve, tmid, eval_flags);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Distance from midpoint of chord to actual curve point */
	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));

	if (dist > tolerance)
	{
		status = adaptive_eval_subdivide_2d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_eval_subdivide_2d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		status = qaws_curve_evaluate_2d(
			&end_result, curve, t1, eval_flags);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_2d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

struct adaptive_eval_buf_3d
{
	qaws_eval_result_3d *samples;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_eval_buf_3d_push(struct adaptive_eval_buf_3d *buf,
	qaws_eval_result_3d sample)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_eval_result_3d *new_ptr = (qaws_eval_result_3d *)qaws_internal_realloc(
			buf->samples, qaws_internal_mul2(new_cap, sizeof(qaws_eval_result_3d)));
		if (!new_ptr)
			return 0;
		buf->samples = new_ptr;
		buf->capacity = new_cap;
	}
	buf->samples[buf->count++] = sample;
	return 1;
}

static qaws_status adaptive_eval_subdivide_3d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec3 p0,
	qaws_scalar t1, qaws_vec3 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_eval_buf_3d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_3d result;
	qaws_eval_result_3d end_result;
	qaws_status status;
	qaws_scalar mx, my, mz, dx, dy, dz, dist;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		status = qaws_curve_evaluate_3d(
			&end_result, curve, t1, eval_flags);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_3d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_3d(
		&result, curve, tmid, eval_flags);
	if (status != QAWS_STATUS_OK)
		return status;

	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	mz = (p0.z + p1.z) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dz = result.position.z - mz;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy + dz * dz));

	if (dist > tolerance)
	{
		status = adaptive_eval_subdivide_3d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_eval_subdivide_3d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		status = qaws_curve_evaluate_3d(
			&end_result, curve, t1, eval_flags);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_3d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_eval_2d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_eval_result_2d *out_samples,
	unsigned int sample_capacity,
	unsigned int *out_sample_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_2d result;
	qaws_status status;
	struct adaptive_eval_buf_2d buf;
	qaws_scalar t0, t1;
	qaws_vec2 p0, p1;
	unsigned int copy_count;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.samples = (qaws_eval_result_2d *)qaws_internal_malloc(
		qaws_internal_mul2(buf.capacity, sizeof(qaws_eval_result_2d)));
	if (!buf.samples)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Evaluate first point */
	status = qaws_curve_evaluate_2d(
		&result, curve, range_min, eval_flags);
	if (status != QAWS_STATUS_OK)
	{
		qaws_internal_free(buf.samples);
		return status;
	}
	p0 = result.position;
	if (!adaptive_eval_buf_2d_push(&buf, result))
	{
		qaws_internal_free(buf.samples);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_2d(
			&result, curve, t1, eval_flags);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.samples);
			return status;
		}
		p1 = result.position;

		status = adaptive_eval_subdivide_2d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.samples);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > sample_capacity)
	{
		qaws_internal_free(buf.samples);
		*out_sample_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_samples[i] = buf.samples[i];

	*out_sample_count = copy_count;
	qaws_internal_free(buf.samples);
	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_eval_3d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_eval_result_3d *out_samples,
	unsigned int sample_capacity,
	unsigned int *out_sample_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_3d result;
	qaws_status status;
	struct adaptive_eval_buf_3d buf;
	qaws_scalar t0, t1;
	qaws_vec3 p0, p1;
	unsigned int copy_count;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.samples = (qaws_eval_result_3d *)qaws_internal_malloc(
		qaws_internal_mul2(buf.capacity, sizeof(qaws_eval_result_3d)));
	if (!buf.samples)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Evaluate first point */
	status = qaws_curve_evaluate_3d(
		&result, curve, range_min, eval_flags);
	if (status != QAWS_STATUS_OK)
	{
		qaws_internal_free(buf.samples);
		return status;
	}
	p0 = result.position;
	if (!adaptive_eval_buf_3d_push(&buf, result))
	{
		qaws_internal_free(buf.samples);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_3d(
			&result, curve, t1, eval_flags);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.samples);
			return status;
		}
		p1 = result.position;

		status = adaptive_eval_subdivide_3d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(buf.samples);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > sample_capacity)
	{
		qaws_internal_free(buf.samples);
		*out_sample_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_samples[i] = buf.samples[i];

	*out_sample_count = copy_count;
	qaws_internal_free(buf.samples);
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * get_sample_count
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_get_sample_count(
	unsigned int *out_sample_count,
	qaws_curve const *curve,
	qaws_sampling_desc const *desc)
{
	if (!curve || !desc || !out_sample_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_sample_count = desc->sample_count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_2d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_2d(
	qaws_vec2 *out_positions,
	unsigned int *out_position_count,
	qaws_curve const *curve,
	unsigned int position_capacity,
	qaws_sampling_desc const *desc)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_2d result;
	qaws_status status;

	if (!curve || !desc || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_2d(
			curve, desc, out_positions, position_capacity, out_position_count);

	count = desc->sample_count;
	if (count > position_capacity)
		count = position_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_distances)
		{
			qaws_internal_free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			table_params, table_distances, curve, table_size);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(table_params);
			qaws_internal_free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_2d(
				&result, curve, param, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return status;
			}

			out_positions[i] = result.position;
		}

		qaws_internal_free(table_params);
		qaws_internal_free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_2d(
				&result, curve, t, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
				return status;

			out_positions[i] = result.position;
		}
	}

	*out_position_count = count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_3d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_3d(
	qaws_vec3 *out_positions,
	unsigned int *out_position_count,
	qaws_curve const *curve,
	unsigned int position_capacity,
	qaws_sampling_desc const *desc)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_3d result;
	qaws_status status;

	if (!curve || !desc || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_3d(
			curve, desc, out_positions, position_capacity, out_position_count);

	count = desc->sample_count;
	if (count > position_capacity)
		count = position_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_distances)
		{
			qaws_internal_free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			table_params, table_distances, curve, table_size);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(table_params);
			qaws_internal_free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_3d(
				&result, curve, param, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return status;
			}

			out_positions[i] = result.position;
		}

		qaws_internal_free(table_params);
		qaws_internal_free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_3d(
				&result, curve, t, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
				return status;

			out_positions[i] = result.position;
		}
	}

	*out_position_count = count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_eval_2d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_eval_2d(
	qaws_eval_result_2d *out_samples,
	unsigned int *out_sample_count,
	qaws_curve const *curve,
	unsigned int sample_capacity,
	qaws_sampling_desc const *desc)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	unsigned int eval_flags;
	qaws_status status;

	if (!curve || !desc || !out_samples || !out_sample_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_eval_2d(
			curve, desc, out_samples, sample_capacity, out_sample_count);

	count = desc->sample_count;
	if (count > sample_capacity)
		count = sample_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_distances)
		{
			qaws_internal_free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			table_params, table_distances, curve, table_size);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(table_params);
			qaws_internal_free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_2d(
				&out_samples[i], curve, param, eval_flags);
			if (status != QAWS_STATUS_OK)
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return status;
			}
		}

		qaws_internal_free(table_params);
		qaws_internal_free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_2d(
				&out_samples[i], curve, t, eval_flags);
			if (status != QAWS_STATUS_OK)
				return status;
		}
	}

	*out_sample_count = count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_eval_3d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_eval_3d(
	qaws_eval_result_3d *out_samples,
	unsigned int *out_sample_count,
	qaws_curve const *curve,
	unsigned int sample_capacity,
	qaws_sampling_desc const *desc)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	unsigned int eval_flags;
	qaws_status status;

	if (!curve || !desc || !out_samples || !out_sample_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_eval_3d(
			curve, desc, out_samples, sample_capacity, out_sample_count);

	count = desc->sample_count;
	if (count > sample_capacity)
		count = sample_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_distances)
		{
			qaws_internal_free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			table_params, table_distances, curve, table_size);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(table_params);
			qaws_internal_free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_3d(
				&out_samples[i], curve, param, eval_flags);
			if (status != QAWS_STATUS_OK)
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return status;
			}
		}

		qaws_internal_free(table_params);
		qaws_internal_free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_3d(
				&out_samples[i], curve, t, eval_flags);
			if (status != QAWS_STATUS_OK)
				return status;
		}
	}

	*out_sample_count = count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Curvature-weighted sampling (Task 26)
 * ------------------------------------------------------------------------- */

#define CURV_WEIGHT_TABLE_SIZE 256

qaws_status qaws_curve_sample_curvature_weighted_2d(
	qaws_vec2 *out_positions,
	unsigned int *out_position_count,
	qaws_curve const *curve,
	unsigned int sample_count,
	unsigned int position_capacity)
{
	unsigned int n;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar *t_table;
	qaws_scalar *curv_table;
	qaws_scalar *cum_weight;
	qaws_scalar epsilon;
	qaws_scalar total;
	qaws_scalar target;
	qaws_scalar w;
	qaws_status status;
	qaws_eval_result_2d result;

	if (!curve || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (position_capacity < sample_count)
	{
		*out_position_count = sample_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	n = CURV_WEIGHT_TABLE_SIZE;
	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;
	epsilon = (qaws_scalar)1e-6;

	t_table = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(n, sizeof(qaws_scalar)));
	if (!t_table)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	curv_table = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(n, sizeof(qaws_scalar)));
	if (!curv_table)
	{
		qaws_internal_free(t_table);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	cum_weight = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(n, sizeof(qaws_scalar)));
	if (!cum_weight)
	{
		qaws_internal_free(t_table);
		qaws_internal_free(curv_table);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Step 1: sample curvature at uniform parameter points */
	for (i = 0; i < n; ++i)
	{
		qaws_scalar kappa;

		t_table[i] = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(n - 1);

		status = qaws_curve_compute_curvature_2d(
			&kappa, curve, t_table[i]);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(t_table);
			qaws_internal_free(curv_table);
			qaws_internal_free(cum_weight);
			return status;
		}

		w = kappa;
		if (w < (qaws_scalar)0.0)
			w = -w;
		if (w < epsilon)
			w = epsilon;
		curv_table[i] = w;
	}

	/* Step 2: build cumulative curvature-weighted table */
	cum_weight[0] = (qaws_scalar)0.0;
	for (i = 1; i < n; ++i)
	{
		cum_weight[i] = cum_weight[i - 1]
			+ curv_table[i] * (t_table[i] - t_table[i - 1]);
	}

	total = cum_weight[n - 1];

	/* Step 3: for each output sample, binary search in cum_weight */
	for (i = 0; i < sample_count; ++i)
	{
		unsigned int lo, hi, mid;
		qaws_scalar frac, t_val;

		if (sample_count == 1)
			target = (qaws_scalar)0.0;
		else
			target = (qaws_scalar)i * total / (qaws_scalar)(sample_count - 1);

		/* Binary search */
		lo = 0;
		hi = n - 1;
		while (lo < hi - 1)
		{
			mid = (lo + hi) / 2;
			if (cum_weight[mid] <= target)
				lo = mid;
			else
				hi = mid;
		}

		/* Linear interpolation between lo and hi */
		if (cum_weight[hi] - cum_weight[lo] > (qaws_scalar)0.0)
			frac = (target - cum_weight[lo])
				/ (cum_weight[hi] - cum_weight[lo]);
		else
			frac = (qaws_scalar)0.0;

		t_val = t_table[lo] + frac * (t_table[hi] - t_table[lo]);

		status = qaws_curve_evaluate_2d(
			&result, curve, t_val, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(t_table);
			qaws_internal_free(curv_table);
			qaws_internal_free(cum_weight);
			return status;
		}

		out_positions[i] = result.position;
	}

	qaws_internal_free(t_table);
	qaws_internal_free(curv_table);
	qaws_internal_free(cum_weight);

	*out_position_count = sample_count;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_sample_curvature_weighted_3d(
	qaws_vec3 *out_positions,
	unsigned int *out_position_count,
	qaws_curve const *curve,
	unsigned int sample_count,
	unsigned int position_capacity)
{
	unsigned int n;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar *t_table;
	qaws_scalar *curv_table;
	qaws_scalar *cum_weight;
	qaws_scalar epsilon;
	qaws_scalar total;
	qaws_scalar target;
	qaws_scalar w;
	qaws_status status;
	qaws_eval_result_3d result;

	if (!curve || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (position_capacity < sample_count)
	{
		*out_position_count = sample_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	n = CURV_WEIGHT_TABLE_SIZE;
	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;
	epsilon = (qaws_scalar)1e-6;

	t_table = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(n, sizeof(qaws_scalar)));
	if (!t_table)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	curv_table = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(n, sizeof(qaws_scalar)));
	if (!curv_table)
	{
		qaws_internal_free(t_table);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	cum_weight = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(n, sizeof(qaws_scalar)));
	if (!cum_weight)
	{
		qaws_internal_free(t_table);
		qaws_internal_free(curv_table);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Step 1: sample curvature at uniform parameter points */
	for (i = 0; i < n; ++i)
	{
		qaws_scalar kappa;

		t_table[i] = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(n - 1);

		status = qaws_curve_compute_curvature_3d(
			&kappa, curve, t_table[i]);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(t_table);
			qaws_internal_free(curv_table);
			qaws_internal_free(cum_weight);
			return status;
		}

		w = kappa;
		if (w < (qaws_scalar)0.0)
			w = -w;
		if (w < epsilon)
			w = epsilon;
		curv_table[i] = w;
	}

	/* Step 2: build cumulative curvature-weighted table */
	cum_weight[0] = (qaws_scalar)0.0;
	for (i = 1; i < n; ++i)
	{
		cum_weight[i] = cum_weight[i - 1]
			+ curv_table[i] * (t_table[i] - t_table[i - 1]);
	}

	total = cum_weight[n - 1];

	/* Step 3: for each output sample, binary search in cum_weight */
	for (i = 0; i < sample_count; ++i)
	{
		unsigned int lo, hi, mid;
		qaws_scalar frac, t_val;

		if (sample_count == 1)
			target = (qaws_scalar)0.0;
		else
			target = (qaws_scalar)i * total / (qaws_scalar)(sample_count - 1);

		/* Binary search */
		lo = 0;
		hi = n - 1;
		while (lo < hi - 1)
		{
			mid = (lo + hi) / 2;
			if (cum_weight[mid] <= target)
				lo = mid;
			else
				hi = mid;
		}

		/* Linear interpolation between lo and hi */
		if (cum_weight[hi] - cum_weight[lo] > (qaws_scalar)0.0)
			frac = (target - cum_weight[lo])
				/ (cum_weight[hi] - cum_weight[lo]);
		else
			frac = (qaws_scalar)0.0;

		t_val = t_table[lo] + frac * (t_table[hi] - t_table[lo]);

		status = qaws_curve_evaluate_3d(
			&result, curve, t_val, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(t_table);
			qaws_internal_free(curv_table);
			qaws_internal_free(cum_weight);
			return status;
		}

		out_positions[i] = result.position;
	}

	qaws_internal_free(t_table);
	qaws_internal_free(curv_table);
	qaws_internal_free(cum_weight);

	*out_position_count = sample_count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Feature-preserving sampling (Task 27)
 * ------------------------------------------------------------------------- */

#define FEATURE_MAX_INFLECTIONS 128
#define FEATURE_MAX_EXTREMA_PER_AXIS 128
#define FEATURE_DEDUP_TOL ((qaws_scalar)1e-10)

/* Insert a parameter into a sorted array, maintaining sort order and
   skipping duplicates within FEATURE_DEDUP_TOL. Returns updated count. */
static unsigned int feature_insert_unique(
	qaws_scalar *params,
	unsigned int count,
	qaws_scalar value)
{
	unsigned int i, j;
	qaws_scalar diff;

	/* Check for duplicates */
	for (i = 0; i < count; ++i)
	{
		diff = params[i] - value;
		if (diff < (qaws_scalar)0.0)
			diff = -diff;
		if (diff < FEATURE_DEDUP_TOL)
			return count; /* duplicate */
	}

	/* Find insertion position */
	i = count;
	while (i > 0 && params[i - 1] > value)
	{
		--i;
	}

	/* Shift elements right */
	for (j = count; j > i; --j)
		params[j] = params[j - 1];

	params[i] = value;
	return count + 1;
}

qaws_status qaws_curve_sample_feature_preserving_2d(
	qaws_vec2 *out_positions,
	unsigned int *out_position_count,
	qaws_curve const *curve,
	unsigned int base_sample_count,
	unsigned int position_capacity)
{
	unsigned int max_params;
	qaws_scalar *params;
	unsigned int param_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_status status;
	qaws_eval_result_2d result;

	/* Temporary buffers for feature detection */
	qaws_scalar inflections[FEATURE_MAX_INFLECTIONS];
	unsigned int inflection_count;
	qaws_scalar extrema_buf[FEATURE_MAX_EXTREMA_PER_AXIS];
	unsigned int extrema_count;

	if (!curve || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (base_sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Allocate parameter buffer: base + inflections + extrema for 2 axes */
	max_params = base_sample_count
		+ FEATURE_MAX_INFLECTIONS
		+ FEATURE_MAX_EXTREMA_PER_AXIS * 2;

	params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(max_params, sizeof(qaws_scalar)));
	if (!params)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Step 1: generate uniform base parameters */
	param_count = 0;
	for (i = 0; i < base_sample_count; ++i)
	{
		qaws_scalar t;

		if (base_sample_count == 1)
			t = range_min;
		else
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)(base_sample_count - 1);

		params[param_count++] = t;
	}

	/* Step 2: find inflection points */
	inflection_count = 0;
	status = qaws_curve_find_inflection_points(
		inflections, &inflection_count, curve, FEATURE_MAX_INFLECTIONS);
	if (status == QAWS_STATUS_OK)
	{
		for (i = 0; i < inflection_count; ++i)
			param_count = feature_insert_unique(
				params, param_count, inflections[i]);
	}

	/* Step 3: find extrema for axis 0 (x) and axis 1 (y) */
	{
		unsigned int axis;
		for (axis = 0; axis < 2; ++axis)
		{
			extrema_count = 0;
			status = qaws_curve_find_extrema(
				extrema_buf, &extrema_count, curve,
				FEATURE_MAX_EXTREMA_PER_AXIS, axis);
			if (status == QAWS_STATUS_OK)
			{
				for (i = 0; i < extrema_count; ++i)
					param_count = feature_insert_unique(
						params, param_count, extrema_buf[i]);
			}
		}
	}

	/* Check capacity */
	if (param_count > position_capacity)
	{
		*out_position_count = param_count;
		qaws_internal_free(params);
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	/* Step 4: evaluate positions at all merged parameters */
	for (i = 0; i < param_count; ++i)
	{
		status = qaws_curve_evaluate_2d(
			&result, curve, params[i], QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(params);
			return status;
		}
		out_positions[i] = result.position;
	}

	*out_position_count = param_count;
	qaws_internal_free(params);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_sample_feature_preserving_3d(
	qaws_vec3 *out_positions,
	unsigned int *out_position_count,
	qaws_curve const *curve,
	unsigned int base_sample_count,
	unsigned int position_capacity)
{
	unsigned int max_params;
	qaws_scalar *params;
	unsigned int param_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_status status;
	qaws_eval_result_3d result;

	/* Temporary buffers for feature detection */
	qaws_scalar inflections[FEATURE_MAX_INFLECTIONS];
	unsigned int inflection_count;
	qaws_scalar extrema_buf[FEATURE_MAX_EXTREMA_PER_AXIS];
	unsigned int extrema_count;

	if (!curve || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (base_sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Allocate parameter buffer: base + inflections + extrema for 3 axes */
	max_params = base_sample_count
		+ FEATURE_MAX_INFLECTIONS
		+ FEATURE_MAX_EXTREMA_PER_AXIS * 3;

	params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(max_params, sizeof(qaws_scalar)));
	if (!params)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Step 1: generate uniform base parameters */
	param_count = 0;
	for (i = 0; i < base_sample_count; ++i)
	{
		qaws_scalar t;

		if (base_sample_count == 1)
			t = range_min;
		else
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)(base_sample_count - 1);

		params[param_count++] = t;
	}

	/* Step 2: find inflection points */
	inflection_count = 0;
	status = qaws_curve_find_inflection_points(
		inflections, &inflection_count, curve, FEATURE_MAX_INFLECTIONS);
	if (status == QAWS_STATUS_OK)
	{
		for (i = 0; i < inflection_count; ++i)
			param_count = feature_insert_unique(
				params, param_count, inflections[i]);
	}

	/* Step 3: find extrema for axis 0 (x), 1 (y), 2 (z) */
	{
		unsigned int axis;
		for (axis = 0; axis < 3; ++axis)
		{
			extrema_count = 0;
			status = qaws_curve_find_extrema(
				extrema_buf, &extrema_count, curve,
				FEATURE_MAX_EXTREMA_PER_AXIS, axis);
			if (status == QAWS_STATUS_OK)
			{
				for (i = 0; i < extrema_count; ++i)
					param_count = feature_insert_unique(
						params, param_count, extrema_buf[i]);
			}
		}
	}

	/* Check capacity */
	if (param_count > position_capacity)
	{
		*out_position_count = param_count;
		qaws_internal_free(params);
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	/* Step 4: evaluate positions at all merged parameters */
	for (i = 0; i < param_count; ++i)
	{
		status = qaws_curve_evaluate_3d(
			&result, curve, params[i], QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(params);
			return status;
		}
		out_positions[i] = result.position;
	}

	*out_position_count = param_count;
	qaws_internal_free(params);
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Streaming / callback-based sampling (Task 28)
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_streaming_2d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_sample_callback_2d callback,
	void *user_data)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar t;
	qaws_eval_result_2d result;
	qaws_status status;

	if (!curve || !desc || !callback)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	count = desc->sample_count;
	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_distances)
		{
			qaws_internal_free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			table_params, table_distances, curve, table_size);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(table_params);
			qaws_internal_free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length
					/ (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_2d(
				&result, curve, param, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return status;
			}

			if (!callback(param, &result.position, user_data))
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return QAWS_STATUS_OK;
			}
		}

		qaws_internal_free(table_params);
		qaws_internal_free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min)
					/ (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_2d(
				&result, curve, t, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
				return status;

			if (!callback(t, &result.position, user_data))
				return QAWS_STATUS_OK;
		}
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_sample_streaming_3d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_sample_callback_3d callback,
	void *user_data)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar t;
	qaws_eval_result_3d result;
	qaws_status status;

	if (!curve || !desc || !callback)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	count = desc->sample_count;
	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)qaws_internal_scratch_malloc(qaws_internal_mul2(table_size, sizeof(qaws_scalar)));
		if (!table_distances)
		{
			qaws_internal_free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			table_params, table_distances, curve, table_size);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(table_params);
			qaws_internal_free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length
					/ (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_3d(
				&result, curve, param, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return status;
			}

			if (!callback(param, &result.position, user_data))
			{
				qaws_internal_free(table_params);
				qaws_internal_free(table_distances);
				return QAWS_STATUS_OK;
			}
		}

		qaws_internal_free(table_params);
		qaws_internal_free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min)
					/ (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_3d(
				&result, curve, t, QAWS_EVAL_FLAG_POSITION);
			if (status != QAWS_STATUS_OK)
				return status;

			if (!callback(t, &result.position, user_data))
				return QAWS_STATUS_OK;
		}
	}

	return QAWS_STATUS_OK;
}
