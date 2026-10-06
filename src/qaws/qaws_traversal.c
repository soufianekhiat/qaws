#include "qaws_traversal.h"
#include "qaws_diff_sampling.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_cache.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Easing
 * ------------------------------------------------------------------------- */

#define QAWS_PI (qaws_scalar)3.14159265358979323846
#define QAWS_CUSTOM_SPEED_TABLE_SIZE 256

static qaws_scalar apply_easing(qaws_easing easing, qaws_scalar t)
{
	qaws_scalar inv;

	switch (easing)
	{
	case QAWS_EASING_LINEAR:
		return t;

	case QAWS_EASING_QUAD_IN:
		return t * t;

	case QAWS_EASING_QUAD_OUT:
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - inv * inv;

	case QAWS_EASING_QUAD_IN_OUT:
		if (t < (qaws_scalar)0.5)
			return (qaws_scalar)2.0 * t * t;
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - (qaws_scalar)2.0 * inv * inv;

	case QAWS_EASING_CUBIC_IN:
		return t * t * t;

	case QAWS_EASING_CUBIC_OUT:
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - inv * inv * inv;

	case QAWS_EASING_CUBIC_IN_OUT:
		if (t < (qaws_scalar)0.5)
			return (qaws_scalar)4.0 * t * t * t;
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - (qaws_scalar)4.0 * inv * inv * inv;

	case QAWS_EASING_SINE_IN:
		return (qaws_scalar)1.0 - (qaws_scalar)cos((double)(t * QAWS_PI / (qaws_scalar)2.0));

	case QAWS_EASING_SINE_OUT:
		return (qaws_scalar)sin((double)(t * QAWS_PI / (qaws_scalar)2.0));

	case QAWS_EASING_SINE_IN_OUT:
		return ((qaws_scalar)1.0 - (qaws_scalar)cos((double)(t * QAWS_PI))) / (qaws_scalar)2.0;

	default:
		return t;
	}
}

/* ---------------------------------------------------------------------------
 * Wrap mode
 * ------------------------------------------------------------------------- */

static qaws_scalar apply_wrap(qaws_wrap_mode mode, qaws_scalar distance, qaws_scalar total_length)
{
	qaws_scalar wrapped;
	qaws_scalar cycle;

	if (total_length <= (qaws_scalar)0.0)
		return (qaws_scalar)0.0;

	switch (mode)
	{
	case QAWS_WRAP_CLAMP:
		if (distance < (qaws_scalar)0.0)
			return (qaws_scalar)0.0;
		if (distance > total_length)
			return total_length;
		return distance;

	case QAWS_WRAP_LOOP:
		wrapped = (qaws_scalar)fmod((double)distance, (double)total_length);
		if (wrapped < (qaws_scalar)0.0)
			wrapped += total_length;
		return wrapped;

	case QAWS_WRAP_PING_PONG:
		/* Fold into [0, 2*total_length) then mirror */
		cycle = (qaws_scalar)2.0 * total_length;
		wrapped = (qaws_scalar)fmod((double)distance, (double)cycle);
		if (wrapped < (qaws_scalar)0.0)
			wrapped += cycle;
		if (wrapped > total_length)
			wrapped = cycle - wrapped;
		return wrapped;

	default:
		return distance;
	}
}

/* ---------------------------------------------------------------------------
 * Motion profile: map time to distance
 * ------------------------------------------------------------------------- */

/* S-curve profile: distance, speed and acceleration at a time. */
static qaws_scalar scurve_eval(
	qaws_traversal_desc const *desc,
	qaws_scalar time,
	qaws_scalar *out_speed,
	qaws_scalar *out_accel)
{
	qaws_scalar dt;
	qaws_scalar total_time;
	qaws_scalar V;
	qaws_scalar A;
	qaws_scalar J;
	qaws_scalar t_j;
	qaws_scalar t_a;
	qaws_scalar t_v;
	qaws_scalar t_total_accel;
	qaws_scalar seg_times[7];
	qaws_scalar seg_v[8]; /* velocities at segment boundaries */
	qaws_scalar seg_d[8]; /* cumulative distances at segment boundaries */
	qaws_scalar local_t = (qaws_scalar)0.0;
	qaws_scalar v_at_start;
	qaws_scalar d_at_start;
	qaws_scalar lt;
	int seg;

	*out_speed = (qaws_scalar)0.0;
	*out_accel = (qaws_scalar)0.0;
	dt = time - desc->start_time;
	total_time = desc->end_time - desc->start_time;

	if (total_time <= (qaws_scalar)0.0 || dt <= (qaws_scalar)0.0)
		return (qaws_scalar)0.0;

	if (dt > total_time)
		dt = total_time;

	V = desc->max_speed > (qaws_scalar)0.0 ? desc->max_speed : (qaws_scalar)1.0;
	A = desc->acceleration > (qaws_scalar)0.0 ? desc->acceleration : (qaws_scalar)1.0;
	J = desc->jerk > (qaws_scalar)0.0 ? desc->jerk : (qaws_scalar)1.0;

	/* Compute jerk time and constant-acceleration time */
	t_j = A / J;
	t_a = V / A - t_j;

	if (t_a < (qaws_scalar)0.0)
	{
		/* Triangular profile: can't reach max acceleration */
		t_j = (qaws_scalar)sqrt((double)(V / J));
		t_a = (qaws_scalar)0.0;
	}

	/* Total time for acceleration phase (segments 1-3) */
	t_total_accel = (qaws_scalar)2.0 * t_j + t_a;

	/* Cruise time */
	t_v = total_time - (qaws_scalar)2.0 * t_total_accel;
	if (t_v < (qaws_scalar)0.0)
	{
		/* Not enough time: scale down proportionally */
		qaws_scalar scale = total_time / ((qaws_scalar)2.0 * t_total_accel);
		t_j *= scale;
		t_a *= scale;
		t_total_accel *= scale;
		t_v = (qaws_scalar)0.0;
		/* Recompute effective A and V after scaling */
		A = J * t_j;
		V = A * (t_j + t_a);
	}

	/* Segment durations:
	 * 1: [0, t_j]         increasing acceleration (jerk = +J)
	 * 2: [t_j, t_j+t_a]  constant acceleration (a = A)
	 * 3: [t_j+t_a, 2*t_j+t_a] decreasing acceleration (jerk = -J)
	 * 4: cruise at V for t_v
	 * 5-7: symmetric deceleration */
	seg_times[0] = t_j;
	seg_times[1] = t_a;
	seg_times[2] = t_j;
	seg_times[3] = t_v;
	seg_times[4] = t_j;
	seg_times[5] = t_a;
	seg_times[6] = t_j;

	/* Compute velocity and distance at each segment boundary */
	seg_v[0] = (qaws_scalar)0.0;
	seg_d[0] = (qaws_scalar)0.0;

	/* Segment 1: jerk = +J, a(t) = J*t, v(t) = J*t^2/2, d(t) = J*t^3/6 */
	{
		qaws_scalar t1 = seg_times[0];
		seg_v[1] = J * t1 * t1 / (qaws_scalar)2.0;
		seg_d[1] = seg_d[0] + J * t1 * t1 * t1 / (qaws_scalar)6.0;
	}

	/* Segment 2: constant acceleration A, v(t) = v1 + A*t, d(t) = d1 + v1*t + A*t^2/2 */
	{
		qaws_scalar t2 = seg_times[1];
		seg_v[2] = seg_v[1] + A * t2;
		seg_d[2] = seg_d[1] + seg_v[1] * t2 + A * t2 * t2 / (qaws_scalar)2.0;
	}

	/* Segment 3: jerk = -J, a(t) = A - J*t, v(t) = v2 + A*t - J*t^2/2,
	 * d(t) = d2 + v2*t + A*t^2/2 - J*t^3/6 */
	{
		qaws_scalar t3 = seg_times[2];
		seg_v[3] = seg_v[2] + A * t3 - J * t3 * t3 / (qaws_scalar)2.0;
		seg_d[3] = seg_d[2] + seg_v[2] * t3 + A * t3 * t3 / (qaws_scalar)2.0
			- J * t3 * t3 * t3 / (qaws_scalar)6.0;
	}

	/* Segment 4: cruise at V */
	{
		qaws_scalar t4 = seg_times[3];
		seg_v[4] = seg_v[3]; /* should be ~V */
		seg_d[4] = seg_d[3] + seg_v[3] * t4;
	}

	/* Segment 5: jerk = -J (start deceleration), a(t) = -J*t,
	 * v(t) = V - J*t^2/2, d(t) = d4 + V*t - J*t^3/6 */
	{
		qaws_scalar t5 = seg_times[4];
		seg_v[5] = seg_v[4] - J * t5 * t5 / (qaws_scalar)2.0;
		seg_d[5] = seg_d[4] + seg_v[4] * t5 - J * t5 * t5 * t5 / (qaws_scalar)6.0;
	}

	/* Segment 6: constant deceleration -A, v(t) = v5 - A*t,
	 * d(t) = d5 + v5*t - A*t^2/2 */
	{
		qaws_scalar t6 = seg_times[5];
		seg_v[6] = seg_v[5] - A * t6;
		seg_d[6] = seg_d[5] + seg_v[5] * t6 - A * t6 * t6 / (qaws_scalar)2.0;
	}

	/* Segment 7: jerk = +J (end deceleration), a(t) = -A + J*t,
	 * v(t) = v6 - A*t + J*t^2/2, d(t) = d6 + v6*t - A*t^2/2 + J*t^3/6 */
	{
		qaws_scalar t7 = seg_times[6];
		seg_v[7] = seg_v[6] - A * t7 + J * t7 * t7 / (qaws_scalar)2.0;
		seg_d[7] = seg_d[6] + seg_v[6] * t7 - A * t7 * t7 / (qaws_scalar)2.0
			+ J * t7 * t7 * t7 / (qaws_scalar)6.0;
	}

	/* Find which segment dt falls into */
	{
		qaws_scalar cumulative = (qaws_scalar)0.0;
		seg = -1;
		for (int i = 0; i < 7; i++)
		{
			if (dt <= cumulative + seg_times[i])
			{
				seg = i;
				local_t = dt - cumulative;
				break;
			}
			cumulative += seg_times[i];
		}
		if (seg < 0)
		{
			/* Past end, return total distance */
			return seg_d[7];
		}
	}

	v_at_start = seg_v[seg];
	d_at_start = seg_d[seg];
	lt = local_t;

	switch (seg)
	{
	case 0: /* Segment 1: jerk = +J */
		*out_speed = J * lt * lt / (qaws_scalar)2.0;
		*out_accel = J * lt;
		return d_at_start + J * lt * lt * lt / (qaws_scalar)6.0;

	case 1: /* Segment 2: constant acceleration A */
		*out_speed = v_at_start + A * lt;
		*out_accel = A;
		return d_at_start + v_at_start * lt + A * lt * lt / (qaws_scalar)2.0;

	case 2: /* Segment 3: jerk = -J */
		*out_speed = v_at_start + A * lt - J * lt * lt / (qaws_scalar)2.0;
		*out_accel = A - J * lt;
		return d_at_start + v_at_start * lt + A * lt * lt / (qaws_scalar)2.0
			- J * lt * lt * lt / (qaws_scalar)6.0;

	case 3: /* Segment 4: cruise */
		*out_speed = v_at_start;
		return d_at_start + v_at_start * lt;

	case 4: /* Segment 5: jerk = -J (begin decel) */
		*out_speed = v_at_start - J * lt * lt / (qaws_scalar)2.0;
		*out_accel = -J * lt;
		return d_at_start + v_at_start * lt - J * lt * lt * lt / (qaws_scalar)6.0;

	case 5: /* Segment 6: constant deceleration -A */
		*out_speed = v_at_start - A * lt;
		*out_accel = -A;
		return d_at_start + v_at_start * lt - A * lt * lt / (qaws_scalar)2.0;

	case 6: /* Segment 7: jerk = +J (end decel) */
		*out_speed = v_at_start - A * lt + J * lt * lt / (qaws_scalar)2.0;
		*out_accel = -A + J * lt;
		return d_at_start + v_at_start * lt - A * lt * lt / (qaws_scalar)2.0
			+ J * lt * lt * lt / (qaws_scalar)6.0;

	default:
		return seg_d[7];
	}
}

static qaws_scalar scurve_map_time_to_distance(
	qaws_traversal_desc const *desc,
	qaws_scalar time)
{
	qaws_scalar v, a;
	return scurve_eval(desc, time, &v, &a);
}

static qaws_scalar custom_speed_map_time_to_distance(
	qaws_traversal const *traversal,
	qaws_scalar time)
{
	unsigned int low;
	unsigned int high;
	unsigned int mid;
	qaws_scalar t_range;
	qaws_scalar alpha;

	if (!traversal->custom_time_table || traversal->custom_table_size == 0)
		return (qaws_scalar)0.0;

	if (time <= traversal->custom_time_table[0])
		return traversal->custom_dist_table[0];
	if (time >= traversal->custom_time_table[traversal->custom_table_size - 1])
		return traversal->custom_dist_table[traversal->custom_table_size - 1];

	/* Binary search */
	low = 0;
	high = traversal->custom_table_size - 1;

	while (high - low > 1)
	{
		mid = (low + high) / 2;
		if (traversal->custom_time_table[mid] <= time)
			low = mid;
		else
			high = mid;
	}

	t_range = traversal->custom_time_table[high] - traversal->custom_time_table[low];
	if (t_range <= (qaws_scalar)0.0)
		return traversal->custom_dist_table[low];

	alpha = (time - traversal->custom_time_table[low]) / t_range;
	return traversal->custom_dist_table[low]
		+ alpha * (traversal->custom_dist_table[high] - traversal->custom_dist_table[low]);
}

static qaws_scalar traversal_map_time_to_distance(
	qaws_traversal_desc const *desc,
	qaws_scalar time)
{
	qaws_scalar dt;
	qaws_scalar total_time;
	qaws_scalar accel;
	qaws_scalar t_ramp;
	qaws_scalar d_accel;
	qaws_scalar d_cruise;
	qaws_scalar d_total;
	qaws_scalar dt_decel;

	switch (desc->motion_profile)
	{
	case QAWS_MOTION_PROFILE_CONSTANT_SPEED:
		dt = time - desc->start_time;
		return desc->speed * dt;

	case QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION:
		dt = time - desc->start_time;
		return desc->speed * dt + (qaws_scalar)0.5 * desc->acceleration * dt * dt;

	case QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED:
		dt = time - desc->start_time;
		total_time = desc->end_time - desc->start_time;

		if (total_time <= (qaws_scalar)0.0)
			return (qaws_scalar)0.0;

		accel = desc->acceleration > (qaws_scalar)0.0 ? desc->acceleration : (qaws_scalar)1.0;
		t_ramp = desc->max_speed / accel;

		if ((qaws_scalar)2.0 * t_ramp > total_time)
			t_ramp = total_time / (qaws_scalar)2.0;

		if (dt <= t_ramp)
		{
			/* Acceleration phase */
			return (qaws_scalar)0.5 * accel * dt * dt;
		}
		else if (dt <= total_time - t_ramp)
		{
			/* Cruise phase */
			d_accel = (qaws_scalar)0.5 * accel * t_ramp * t_ramp;
			return d_accel + desc->max_speed * (dt - t_ramp);
		}
		else
		{
			/* Deceleration phase */
			dt_decel = total_time - dt;
			d_accel = (qaws_scalar)0.5 * accel * t_ramp * t_ramp;
			d_cruise = desc->max_speed * (total_time - (qaws_scalar)2.0 * t_ramp);
			d_total = d_accel + d_cruise + d_accel;
			return d_total - (qaws_scalar)0.5 * accel * dt_decel * dt_decel;
		}

	case QAWS_MOTION_PROFILE_S_CURVE:
		return scurve_map_time_to_distance(desc, time);

	case QAWS_MOTION_PROFILE_NONE:
	default:
		return time;
	}
}

/* ---------------------------------------------------------------------------
 * Multi-curve chain helpers
 * ------------------------------------------------------------------------- */

static void resolve_chain_parameter(
	qaws_traversal const *traversal,
	qaws_scalar combined_param,
	qaws_curve const **out_curve,
	qaws_scalar *out_local_param)
{
	unsigned int i;
	qaws_scalar total_range;

	/* Find which curve this parameter falls into */
	for (i = 0; i < traversal->chain_count - 1; i++)
	{
		if (combined_param < traversal->chain_param_offsets[i + 1])
			break;
	}

	/* Compute local parameter for this curve */
	{
		qaws_scalar offset = traversal->chain_param_offsets[i];
		qaws_scalar range_i = traversal->chain_param_offsets[i + 1] - offset;
		qaws_scalar local_frac;
		qaws_range curve_range;

		if (range_i <= (qaws_scalar)0.0)
		{
			*out_curve = traversal->chain_curves[i];
			*out_local_param = traversal->chain_curves[i]->parameter_range.min_value;
			return;
		}

		local_frac = (combined_param - offset) / range_i;
		if (local_frac < (qaws_scalar)0.0)
			local_frac = (qaws_scalar)0.0;
		if (local_frac > (qaws_scalar)1.0)
			local_frac = (qaws_scalar)1.0;

		curve_range = traversal->chain_curves[i]->parameter_range;
		*out_curve = traversal->chain_curves[i];
		*out_local_param = curve_range.min_value
			+ local_frac * (curve_range.max_value - curve_range.min_value);
	}

	(void)total_range;
}

/* ---------------------------------------------------------------------------
 * Custom speed table builder
 * ------------------------------------------------------------------------- */

static qaws_status build_custom_speed_table(
	qaws_traversal *traversal)
{
	unsigned int n;
	unsigned int i;
	qaws_scalar total_time;
	qaws_scalar dt_step;
	qaws_scalar s;
	qaws_scalar t;

	if (traversal->desc.motion_profile != QAWS_MOTION_PROFILE_CUSTOM
		|| !traversal->desc.custom_speed_fn)
	{
		traversal->custom_time_table = NULL;
		traversal->custom_dist_table = NULL;
		traversal->custom_table_size = 0;
		return QAWS_STATUS_OK;
	}

	n = QAWS_CUSTOM_SPEED_TABLE_SIZE;
	total_time = traversal->desc.end_time - traversal->desc.start_time;

	if (total_time <= (qaws_scalar)0.0)
	{
		traversal->custom_time_table = NULL;
		traversal->custom_dist_table = NULL;
		traversal->custom_table_size = 0;
		return QAWS_STATUS_OK;
	}

	traversal->custom_time_table = (qaws_scalar *)malloc(n * sizeof(qaws_scalar));
	if (!traversal->custom_time_table)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	traversal->custom_dist_table = (qaws_scalar *)malloc(n * sizeof(qaws_scalar));
	if (!traversal->custom_dist_table)
	{
		free(traversal->custom_time_table);
		traversal->custom_time_table = NULL;
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	traversal->custom_table_size = n;
	dt_step = total_time / (qaws_scalar)(n - 1);

	/* Numerically integrate ds/dt = v(s) using RK4 */
	s = (qaws_scalar)0.0;

	for (i = 0; i < n; i++)
	{
		t = traversal->desc.start_time + dt_step * (qaws_scalar)i;
		traversal->custom_time_table[i] = t;
		traversal->custom_dist_table[i] = s;

		if (i < n - 1)
		{
			/* RK4 step: ds/dt = v(s) */
			qaws_scalar k1 = traversal->desc.custom_speed_fn(
				s, traversal->desc.custom_speed_fn_user_data);
			qaws_scalar k2 = traversal->desc.custom_speed_fn(
				s + (qaws_scalar)0.5 * dt_step * k1,
				traversal->desc.custom_speed_fn_user_data);
			qaws_scalar k3 = traversal->desc.custom_speed_fn(
				s + (qaws_scalar)0.5 * dt_step * k2,
				traversal->desc.custom_speed_fn_user_data);
			qaws_scalar k4 = traversal->desc.custom_speed_fn(
				s + dt_step * k3,
				traversal->desc.custom_speed_fn_user_data);
			s += dt_step * (k1 + (qaws_scalar)2.0 * k2
				+ (qaws_scalar)2.0 * k3 + k4) / (qaws_scalar)6.0;

			/* Cap at total arc length */
			if (s > traversal->total_arc_length)
				s = traversal->total_arc_length;
		}
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Create / Destroy
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_create(
	qaws_curve const *curve,
	qaws_traversal_desc const *desc,
	qaws_traversal **out_traversal)
{
	qaws_traversal *traversal;
	unsigned int table_size;
	qaws_status status;

	if (!curve || !desc || !out_traversal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_traversal = NULL;

	traversal = (qaws_traversal *)malloc(sizeof(qaws_traversal));
	if (!traversal)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	traversal->curve = curve;
	traversal->desc = *desc;

	/* Initialize multi-curve chain fields */
	traversal->chain_curves = NULL;
	traversal->chain_count = 0;
	traversal->chain_param_offsets = NULL;

	/* Initialize custom speed table fields */
	traversal->custom_time_table = NULL;
	traversal->custom_dist_table = NULL;
	traversal->custom_table_size = 0;

	table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
	traversal->table_size = table_size;

	traversal->table_params = (qaws_scalar *)malloc(
		table_size * sizeof(qaws_scalar));
	if (!traversal->table_params)
	{
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	traversal->table_distances = (qaws_scalar *)malloc(
		table_size * sizeof(qaws_scalar));
	if (!traversal->table_distances)
	{
		free(traversal->table_params);
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	status = qaws_internal_build_arc_length_table(
		curve, table_size, traversal->table_params, traversal->table_distances);
	if (status != QAWS_STATUS_OK)
	{
		free(traversal->table_distances);
		free(traversal->table_params);
		free(traversal);
		return status;
	}

	traversal->total_arc_length = traversal->table_distances[table_size - 1];
	traversal->current_distance = (qaws_scalar)0.0;

	/* Build custom speed table if needed */
	status = build_custom_speed_table(traversal);
	if (status != QAWS_STATUS_OK)
	{
		free(traversal->table_distances);
		free(traversal->table_params);
		free(traversal);
		return status;
	}

	*out_traversal = traversal;
	return QAWS_STATUS_OK;
}

qaws_status qaws_traversal_create_multi(
	qaws_curve const *const *curves,
	unsigned int curve_count,
	qaws_traversal_desc const *desc,
	qaws_traversal **out_traversal)
{
	qaws_traversal *traversal;
	unsigned int table_size;
	qaws_scalar *param_offsets;
	qaws_curve const **chain;
	unsigned int i;
	qaws_scalar combined_range;
	qaws_status status;

	if (!curves || curve_count == 0 || !desc || !out_traversal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_traversal = NULL;

	/* Validate: all curves must have same dimension */
	{
		qaws_dimension dim0 = curves[0]->dimension;
		for (i = 1; i < curve_count; i++)
		{
			if (!curves[i])
				return QAWS_STATUS_INVALID_ARGUMENT;
			if (curves[i]->dimension != dim0)
				return QAWS_STATUS_INVALID_DIMENSION;
		}
	}

	traversal = (qaws_traversal *)malloc(sizeof(qaws_traversal));
	if (!traversal)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	traversal->curve = curves[0]; /* default curve for single-curve fallback */
	traversal->desc = *desc;

	/* Initialize custom speed table fields */
	traversal->custom_time_table = NULL;
	traversal->custom_dist_table = NULL;
	traversal->custom_table_size = 0;

	/* Build chain arrays */
	chain = (qaws_curve const **)malloc(curve_count * sizeof(qaws_curve const *));
	if (!chain)
	{
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	param_offsets = (qaws_scalar *)malloc((curve_count + 1) * sizeof(qaws_scalar));
	if (!param_offsets)
	{
		free(chain);
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Build cumulative parameter offsets */
	param_offsets[0] = (qaws_scalar)0.0;
	for (i = 0; i < curve_count; i++)
	{
		qaws_range r = curves[i]->parameter_range;
		chain[i] = curves[i];
		param_offsets[i + 1] = param_offsets[i] + (r.max_value - r.min_value);
	}

	traversal->chain_curves = chain;
	traversal->chain_count = curve_count;
	traversal->chain_param_offsets = param_offsets;

	combined_range = param_offsets[curve_count];

	/* Build combined arc-length table */
	table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
	traversal->table_size = table_size;

	traversal->table_params = (qaws_scalar *)malloc(
		table_size * sizeof(qaws_scalar));
	if (!traversal->table_params)
	{
		free(param_offsets);
		free(chain);
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	traversal->table_distances = (qaws_scalar *)malloc(
		table_size * sizeof(qaws_scalar));
	if (!traversal->table_distances)
	{
		free(traversal->table_params);
		free(param_offsets);
		free(chain);
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Fill combined parameter table with uniform spacing over combined range */
	for (i = 0; i < table_size; i++)
	{
		traversal->table_params[i] = combined_range
			* (qaws_scalar)i / (qaws_scalar)(table_size - 1);
	}

	/* Build cumulative arc-length table over the combined parameter space */
	traversal->table_distances[0] = (qaws_scalar)0.0;

	for (i = 1; i < table_size; i++)
	{
		qaws_scalar prev_combined = traversal->table_params[i - 1];
		qaws_scalar curr_combined = traversal->table_params[i];
		qaws_scalar seg_length = (qaws_scalar)0.0;

		/* This interval might span multiple curves, so accumulate piecewise */
		{
			qaws_scalar remaining_start = prev_combined;
			qaws_scalar remaining_end = curr_combined;
			unsigned int ci;

			for (ci = 0; ci < curve_count; ci++)
			{
				qaws_scalar off_min = param_offsets[ci];
				qaws_scalar off_max = param_offsets[ci + 1];
				qaws_scalar curve_range_len = off_max - off_min;
				qaws_scalar seg_start;
				qaws_scalar seg_end;
				qaws_scalar local_start;
				qaws_scalar local_end;
				qaws_range cr;
				qaws_scalar seg_len;

				if (remaining_start >= off_max)
					continue;
				if (remaining_end <= off_min)
					break;

				seg_start = remaining_start < off_min ? off_min : remaining_start;
				seg_end = remaining_end > off_max ? off_max : remaining_end;

				/* Map to local curve parameter */
				cr = curves[ci]->parameter_range;
				if (curve_range_len > (qaws_scalar)0.0)
				{
					local_start = cr.min_value
						+ (seg_start - off_min) / curve_range_len
						* (cr.max_value - cr.min_value);
					local_end = cr.min_value
						+ (seg_end - off_min) / curve_range_len
						* (cr.max_value - cr.min_value);
				}
				else
				{
					local_start = cr.min_value;
					local_end = cr.min_value;
				}

				seg_len = (qaws_scalar)0.0;
				qaws_internal_arc_length(curves[ci], local_start, local_end, &seg_len);
				seg_length += seg_len;
			}
		}

		traversal->table_distances[i] = traversal->table_distances[i - 1] + seg_length;
	}

	traversal->total_arc_length = traversal->table_distances[table_size - 1];
	traversal->current_distance = (qaws_scalar)0.0;

	/* Build custom speed table if needed */
	status = build_custom_speed_table(traversal);
	if (status != QAWS_STATUS_OK)
	{
		free(traversal->table_distances);
		free(traversal->table_params);
		free(param_offsets);
		free(chain);
		free(traversal);
		return status;
	}

	*out_traversal = traversal;
	return QAWS_STATUS_OK;
}

void qaws_traversal_destroy(qaws_traversal *traversal)
{
	if (!traversal)
		return;

	free(traversal->table_params);
	free(traversal->table_distances);

	/* Free multi-curve chain arrays */
	if (traversal->chain_curves)
		free(traversal->chain_curves);
	if (traversal->chain_param_offsets)
		free(traversal->chain_param_offsets);

	/* Free custom speed tables */
	if (traversal->custom_time_table)
		free(traversal->custom_time_table);
	if (traversal->custom_dist_table)
		free(traversal->custom_dist_table);

	free(traversal);
}

/* ---------------------------------------------------------------------------
 * Evaluate
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_evaluate_2d(
	qaws_traversal const *traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_scalar parameter;
	qaws_scalar distance;
	qaws_range range;
	qaws_scalar normalized_time;
	qaws_scalar eased_time;
	qaws_scalar time_span;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	switch (traversal->desc.traversal_mode)
	{
	case QAWS_TRAVERSAL_MODE_PARAMETER:
		parameter = input_value;
		break;

	case QAWS_TRAVERSAL_MODE_ARC_LENGTH:
		distance = input_value;
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	case QAWS_TRAVERSAL_MODE_TIME:
		eased_time = input_value;
		if (traversal->desc.easing != QAWS_EASING_LINEAR)
		{
			time_span = traversal->desc.end_time - traversal->desc.start_time;
			if (time_span > (qaws_scalar)0.0)
			{
				normalized_time = (input_value - traversal->desc.start_time) / time_span;
				if (normalized_time < (qaws_scalar)0.0)
					normalized_time = (qaws_scalar)0.0;
				if (normalized_time > (qaws_scalar)1.0)
					normalized_time = (qaws_scalar)1.0;
				normalized_time = apply_easing(traversal->desc.easing, normalized_time);
				eased_time = traversal->desc.start_time + normalized_time * time_span;
			}
		}
		if (traversal->desc.motion_profile == QAWS_MOTION_PROFILE_CUSTOM
			&& traversal->custom_table_size > 0)
		{
			distance = custom_speed_map_time_to_distance(traversal, eased_time);
		}
		else
		{
			distance = traversal_map_time_to_distance(
				&traversal->desc, eased_time);
		}
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	default:
		parameter = input_value;
		break;
	}

	/* Multi-curve chain: resolve combined parameter to specific curve */
	if (traversal->chain_count > 0)
	{
		qaws_curve const *resolved_curve;
		qaws_scalar local_param;

		resolve_chain_parameter(traversal, parameter, &resolved_curve, &local_param);

		if (traversal->desc.clamp_to_domain)
		{
			range = resolved_curve->parameter_range;
			if (local_param < range.min_value)
				local_param = range.min_value;
			if (local_param > range.max_value)
				local_param = range.max_value;
		}

		return qaws_curve_evaluate_2d(
			resolved_curve, local_param, eval_flags, out_result);
	}

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_2d(
		traversal->curve, parameter, eval_flags, out_result);
}

qaws_status qaws_traversal_evaluate_3d(
	qaws_traversal const *traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_scalar parameter;
	qaws_scalar distance;
	qaws_range range;
	qaws_scalar normalized_time;
	qaws_scalar eased_time;
	qaws_scalar time_span;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	switch (traversal->desc.traversal_mode)
	{
	case QAWS_TRAVERSAL_MODE_PARAMETER:
		parameter = input_value;
		break;

	case QAWS_TRAVERSAL_MODE_ARC_LENGTH:
		distance = input_value;
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	case QAWS_TRAVERSAL_MODE_TIME:
		eased_time = input_value;
		if (traversal->desc.easing != QAWS_EASING_LINEAR)
		{
			time_span = traversal->desc.end_time - traversal->desc.start_time;
			if (time_span > (qaws_scalar)0.0)
			{
				normalized_time = (input_value - traversal->desc.start_time) / time_span;
				if (normalized_time < (qaws_scalar)0.0)
					normalized_time = (qaws_scalar)0.0;
				if (normalized_time > (qaws_scalar)1.0)
					normalized_time = (qaws_scalar)1.0;
				normalized_time = apply_easing(traversal->desc.easing, normalized_time);
				eased_time = traversal->desc.start_time + normalized_time * time_span;
			}
		}
		if (traversal->desc.motion_profile == QAWS_MOTION_PROFILE_CUSTOM
			&& traversal->custom_table_size > 0)
		{
			distance = custom_speed_map_time_to_distance(traversal, eased_time);
		}
		else
		{
			distance = traversal_map_time_to_distance(
				&traversal->desc, eased_time);
		}
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	default:
		parameter = input_value;
		break;
	}

	/* Multi-curve chain: resolve combined parameter to specific curve */
	if (traversal->chain_count > 0)
	{
		qaws_curve const *resolved_curve;
		qaws_scalar local_param;

		resolve_chain_parameter(traversal, parameter, &resolved_curve, &local_param);

		if (traversal->desc.clamp_to_domain)
		{
			range = resolved_curve->parameter_range;
			if (local_param < range.min_value)
				local_param = range.min_value;
			if (local_param > range.max_value)
				local_param = range.max_value;
		}

		return qaws_curve_evaluate_3d(
			resolved_curve, local_param, eval_flags, out_result);
	}

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_3d(
		traversal->curve, parameter, eval_flags, out_result);
}

/* ---------------------------------------------------------------------------
 * Mapping utilities
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_map_time_to_parameter(
	qaws_traversal const *traversal,
	qaws_scalar time_value,
	qaws_scalar *out_parameter)
{
	qaws_scalar distance;
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (traversal->desc.motion_profile == QAWS_MOTION_PROFILE_CUSTOM
		&& traversal->custom_table_size > 0)
	{
		distance = custom_speed_map_time_to_distance(traversal, time_value);
	}
	else
	{
		distance = traversal_map_time_to_distance(&traversal->desc, time_value);
	}

	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance);

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	*out_parameter = parameter;
	return QAWS_STATUS_OK;
}

qaws_status qaws_traversal_map_distance_to_parameter(
	qaws_traversal const *traversal,
	qaws_scalar distance_value,
	qaws_scalar *out_parameter)
{
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance_value);

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	*out_parameter = parameter;
	return QAWS_STATUS_OK;
}

qaws_status qaws_traversal_map_parameter_to_distance(
	qaws_traversal const *traversal,
	qaws_scalar parameter,
	qaws_scalar *out_distance)
{
	if (!traversal || !out_distance)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_distance = qaws_internal_parameter_to_distance(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, parameter);

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Advance (iterator) functions
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_advance_2d(
	qaws_traversal *traversal,
	qaws_scalar delta_time,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_scalar distance;
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	traversal->current_distance += delta_time * traversal->desc.speed;
	distance = apply_wrap(traversal->desc.wrap_mode,
		traversal->current_distance, traversal->total_arc_length);

	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance);

	/* Multi-curve chain: resolve combined parameter to specific curve */
	if (traversal->chain_count > 0)
	{
		qaws_curve const *resolved_curve;
		qaws_scalar local_param;

		resolve_chain_parameter(traversal, parameter, &resolved_curve, &local_param);

		if (traversal->desc.clamp_to_domain)
		{
			range = resolved_curve->parameter_range;
			if (local_param < range.min_value)
				local_param = range.min_value;
			if (local_param > range.max_value)
				local_param = range.max_value;
		}

		return qaws_curve_evaluate_2d(
			resolved_curve, local_param, eval_flags, out_result);
	}

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_2d(
		traversal->curve, parameter, eval_flags, out_result);
}

qaws_status qaws_traversal_advance_3d(
	qaws_traversal *traversal,
	qaws_scalar delta_time,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_scalar distance;
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	traversal->current_distance += delta_time * traversal->desc.speed;
	distance = apply_wrap(traversal->desc.wrap_mode,
		traversal->current_distance, traversal->total_arc_length);

	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance);

	/* Multi-curve chain: resolve combined parameter to specific curve */
	if (traversal->chain_count > 0)
	{
		qaws_curve const *resolved_curve;
		qaws_scalar local_param;

		resolve_chain_parameter(traversal, parameter, &resolved_curve, &local_param);

		if (traversal->desc.clamp_to_domain)
		{
			range = resolved_curve->parameter_range;
			if (local_param < range.min_value)
				local_param = range.min_value;
			if (local_param > range.max_value)
				local_param = range.max_value;
		}

		return qaws_curve_evaluate_3d(
			resolved_curve, local_param, eval_flags, out_result);
	}

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_3d(
		traversal->curve, parameter, eval_flags, out_result);
}

qaws_status qaws_traversal_reset(qaws_traversal *traversal)
{
	if (!traversal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	traversal->current_distance = (qaws_scalar)0.0;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Differentiable sampling: traversal inputs as arc-length targets
 * ------------------------------------------------------------------------- */

/* Easing g(n) with g'(n) and g''(n) on [0, 1]. */
static double easing_eval(qaws_easing easing, double n, double* g1, double* g2)
{
	double h = 0.5 * 3.14159265358979323846, inv = 1 - n;
	switch (easing)
	{
	case QAWS_EASING_QUAD_IN:
		*g1 = 2 * n; *g2 = 2;
		return n * n;
	case QAWS_EASING_QUAD_OUT:
		*g1 = 2 * inv; *g2 = -2;
		return 1 - inv * inv;
	case QAWS_EASING_QUAD_IN_OUT:
		if (n < 0.5) { *g1 = 4 * n; *g2 = 4; return 2 * n * n; }
		*g1 = 4 * inv; *g2 = -4;
		return 1 - 2 * inv * inv;
	case QAWS_EASING_CUBIC_IN:
		*g1 = 3 * n * n; *g2 = 6 * n;
		return n * n * n;
	case QAWS_EASING_CUBIC_OUT:
		*g1 = 3 * inv * inv; *g2 = -6 * inv;
		return 1 - inv * inv * inv;
	case QAWS_EASING_CUBIC_IN_OUT:
		if (n < 0.5) { *g1 = 12 * n * n; *g2 = 24 * n; return 4 * n * n * n; }
		*g1 = 12 * inv * inv; *g2 = -24 * inv;
		return 1 - 4 * inv * inv * inv;
	case QAWS_EASING_SINE_IN:
		*g1 = h * sin(n * h); *g2 = h * h * cos(n * h);
		return 1 - cos(n * h);
	case QAWS_EASING_SINE_OUT:
		*g1 = h * cos(n * h); *g2 = -h * h * sin(n * h);
		return sin(n * h);
	case QAWS_EASING_SINE_IN_OUT:
		*g1 = h * sin(2 * h * n); *g2 = 2 * h * h * cos(2 * h * n);
		return 0.5 * (1 - cos(2 * h * n));
	default:
		*g1 = 1; *g2 = 0;
		return n;
	}
}

/* Distance of the motion profile at a time, with its speed and acceleration
   (the same branches as traversal_map_time_to_distance). */
static double profile_eval(qaws_traversal const* traversal, qaws_scalar time, double* v, double* a)
{
	qaws_traversal_desc const* desc = &traversal->desc;
	double dt = (double)time - desc->start_time, total_time = (double)desc->end_time - desc->start_time;
	*v = 0;
	*a = 0;
	if (desc->motion_profile == QAWS_MOTION_PROFILE_CUSTOM && traversal->custom_table_size > 0)
	{
		qaws_scalar const* tt = traversal->custom_time_table;
		qaws_scalar const* dd = traversal->custom_dist_table;
		unsigned int n = traversal->custom_table_size, lo = 0, hi = n - 1;
		if (time <= tt[0] || time >= tt[n - 1])
			return custom_speed_map_time_to_distance(traversal, time);
		while (hi - lo > 1)
		{
			unsigned int mid = (lo + hi) / 2;
			if (tt[mid] <= time) lo = mid; else hi = mid;
		}
		if (tt[hi] > tt[lo])
			*v = ((double)dd[hi] - dd[lo]) / ((double)tt[hi] - tt[lo]);
		return custom_speed_map_time_to_distance(traversal, time);
	}
	switch (desc->motion_profile)
	{
	case QAWS_MOTION_PROFILE_CONSTANT_SPEED:
		*v = desc->speed;
		return desc->speed * dt;
	case QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION:
		*v = desc->speed + desc->acceleration * dt;
		*a = desc->acceleration;
		return desc->speed * dt + 0.5 * desc->acceleration * dt * dt;
	case QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED:
	{
		double accel = desc->acceleration > 0 ? desc->acceleration : 1.0, t_ramp;
		if (total_time <= 0)
			return 0;
		t_ramp = desc->max_speed / accel;
		if (2 * t_ramp > total_time)
			t_ramp = total_time / 2;
		if (dt <= t_ramp)
		{
			*v = accel * dt;
			*a = accel;
			return 0.5 * accel * dt * dt;
		}
		if (dt <= total_time - t_ramp)
		{
			*v = desc->max_speed;
			return 0.5 * accel * t_ramp * t_ramp + desc->max_speed * (dt - t_ramp);
		}
		{
			double dd2 = total_time - dt, d_accel = 0.5 * accel * t_ramp * t_ramp;
			*v = accel * dd2;
			*a = -accel;
			return 2 * d_accel + desc->max_speed * (total_time - 2 * t_ramp) - 0.5 * accel * dd2 * dd2;
		}
	}
	case QAWS_MOTION_PROFILE_S_CURVE:
	{
		qaws_scalar sv, sa, d = scurve_eval(desc, time, &sv, &sa);
		*v = sv;
		*a = sa;
		return d;
	}
	default:
		*v = 1;
		return time;
	}
}

qaws_status qaws_traversal_cdf_targets(
	qaws_traversal const* traversal,
	qaws_scalar const* inputs,
	unsigned int count,
	qaws_cdf_target* out_targets,
	qaws_scalar* out_rate,
	qaws_scalar* out_rate2)
{
	unsigned int i;
	double L;
	if (!traversal || ((!inputs || !out_targets) && count))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (traversal->desc.traversal_mode == QAWS_TRAVERSAL_MODE_PARAMETER || traversal->chain_count > 0)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	L = traversal->total_arc_length;
	for (i = 0; i < count; i++)
	{
		double s, r1 = 1, r2 = 0, dist, frac = 0, sign = 1;
		if (traversal->desc.traversal_mode == QAWS_TRAVERSAL_MODE_TIME)
		{
			qaws_traversal_desc const* desc = &traversal->desc;
			double te = inputs[i], e1 = 1, e2 = 0, span = (double)desc->end_time - desc->start_time, v, a;
			if (desc->easing != QAWS_EASING_LINEAR && span > 0)
			{
				double n = ((double)inputs[i] - desc->start_time) / span, g1, g2, g;
				if (n < 0 || n > 1)
				{
					n = n < 0 ? 0 : 1;
					g = easing_eval(desc->easing, n, &g1, &g2);
					e1 = e2 = 0;
				}
				else
				{
					g = easing_eval(desc->easing, n, &g1, &g2);
					e1 = g1;
					e2 = g2 / span;
				}
				te = desc->start_time + g * span;
			}
			s = profile_eval(traversal, (qaws_scalar)te, &v, &a);
			r1 = v * e1;
			r2 = a * e1 * e1 + v * e2;
		}
		else
			s = inputs[i];
		dist = s;
		if (L > 0)
			switch (traversal->desc.wrap_mode)
			{
			case QAWS_WRAP_CLAMP:
				if (s < 0) { dist = 0; r1 = r2 = 0; }
				else if (s > L) { dist = 0; frac = 1; r1 = r2 = 0; }
				break;
			case QAWS_WRAP_LOOP:
				frac = -floor(s / L);
				break;
			case QAWS_WRAP_PING_PONG:
			{
				double k = floor(s / (2 * L)), w = s - 2 * k * L;
				if (w <= L)
					frac = -2 * k;
				else
				{
					dist = -s;
					frac = 2 * k + 2;
					sign = -1;
				}
				break;
			}
			default:
				break;
			}
		out_targets[i].distance = (qaws_scalar)dist;
		out_targets[i].fraction = (qaws_scalar)frac;
		if (out_rate) out_rate[i] = (qaws_scalar)(sign * r1);
		if (out_rate2) out_rate2[i] = (qaws_scalar)(sign * r2);
	}
	return QAWS_STATUS_OK;
}
