#include "qaws_internal_arc_length.h"
#include "qaws_internal_types.h"
#include "qaws_internal_span.h"
#include "../qaws_platform.h"
#include <math.h>

/* Gauss-Kronrod 7-15 on [-1, 1] (QUADPACK qk15): the Kronrod nodes, the Gauss nodes the odd ones */
static const double gk_nodes[8] = {
	0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
	0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
	0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
	0.207784955007898467600689403773245, 0.0
};

static const double gk_weights[8] = {
	0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
	0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
	0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
	0.204432940075298892414161999234649, 0.209482141084727828012999174891714
};

static const double g7_weights[4] = {
	0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
	0.381830050505118944950369775488975, 0.417959183673469387755102040816327
};

/* |C'(t)| */
static double curve_speed(qaws_curve const* curve, qaws_scalar t)
{
	qaws_scalar local_t;
	unsigned int span_index = qaws_internal_find_span(&local_t, curve, t);
	if (curve->dimension == QAWS_DIMENSION_2D) {
		qaws_eval_result_2d res;
		curve->vtable->eval_span_2d(curve, span_index, local_t, QAWS_EVAL_FLAG_D1, &res);
		return sqrt((double)res.d1.x * res.d1.x + (double)res.d1.y * res.d1.y);
	} else {
		qaws_eval_result_3d res;
		curve->vtable->eval_span_3d(curve, span_index, local_t, QAWS_EVAL_FLAG_D1, &res);
		return sqrt((double)res.d1.x * res.d1.x + (double)res.d1.y * res.d1.y + (double)res.d1.z * res.d1.z);
	}
}

/* The Kronrod estimate over [a, b] and its distance to the Gauss one */
static double gk15(qaws_curve const* curve, double a, double b, double* err)
{
	double h = 0.5 * (b - a), mid = 0.5 * (a + b);
	double f0 = curve_speed(curve, (qaws_scalar)mid);
	double k = gk_weights[7] * f0, g = g7_weights[3] * f0;
	for (unsigned int i = 0; i < 7; i++) {
		double f = curve_speed(curve, (qaws_scalar)(mid - h * gk_nodes[i])) + curve_speed(curve, (qaws_scalar)(mid + h * gk_nodes[i]));
		k += gk_weights[i] * f;
		if (i & 1)
			g += g7_weights[i / 2] * f;
	}
	*err = fabs(h * (k - g));
	return h * k;
}

#define ARC_STACK 40

/* Adaptive on [a, b] (one span): halves until the Gauss and Kronrod estimates agree to QAWS_EPSILON of the length (or
   a 1e-30 floor); deeper than ARC_STACK - 1 halvings keeps the Kronrod estimate. */
static qaws_scalar arc_length_segment(
	qaws_curve const* curve,
	qaws_scalar t_min,
	qaws_scalar t_max)
{
	double stack_a[ARC_STACK], stack_b[ARC_STACK];
	unsigned int top;
	double err, whole = gk15(curve, t_min, t_max, &err), total = 0.0;
	double tol = (double)QAWS_EPSILON * fabs(whole) + 1e-30, m = 0.5 * ((double)t_min + (double)t_max);
	if (err <= tol)
		return (qaws_scalar)whole;
	stack_a[0] = m; stack_b[0] = t_max;
	stack_a[1] = t_min; stack_b[1] = m;
	top = 2;
	while (top > 0) {
		double a = stack_a[top - 1], b = stack_b[top - 1], piece;
		top--;
		piece = gk15(curve, a, b, &err);
		m = 0.5 * (a + b);
		/* each piece gets its share of the tolerance */
		if (err <= tol * (b - a) / ((double)t_max - (double)t_min) || top + 2 > ARC_STACK || !(m > a && m < b)) {
			total += piece;
			continue;
		}
		stack_a[top] = m; stack_b[top] = b; top++;
		stack_a[top] = a; stack_b[top] = m; top++;
	}
	return (qaws_scalar)total;
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

/* rtsafe-style: a Newton step on s(t) - length with s' = |C'|, a bisection when the step leaves the bracket or does
   not halve the previous one; s(t) carried incrementally from the last iterate. */
qaws_status qaws_internal_parameter_at_arc_length(
	qaws_scalar* out_parameter,
	qaws_curve const* curve,
	qaws_scalar t_start,
	qaws_scalar length)
{
	qaws_scalar t_end, total, piece;
	double lo, hi, t, s, step_old, step;
	unsigned int it;

	if (!curve || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;
	t_end = curve->parameter_range.max_value;
	if (!(t_start >= curve->parameter_range.min_value && t_start <= t_end) || !(length >= (qaws_scalar)0))
		return QAWS_STATUS_OUT_OF_RANGE;
	qaws_internal_arc_length(&total, curve, t_start, t_end);
	if (length > total * (1 + 4 * QAWS_EPSILON))
		return QAWS_STATUS_OUT_OF_RANGE;
	if (length >= total || length == (qaws_scalar)0) {
		*out_parameter = length == (qaws_scalar)0 ? t_start : t_end;
		return QAWS_STATUS_OK;
	}

	lo = t_start;
	hi = t_end;
	t = lo + (hi - lo) * ((double)length / (double)total);
	qaws_internal_arc_length(&piece, curve, t_start, (qaws_scalar)t);
	s = piece;
	step_old = hi - lo;
	step = step_old;
	for (it = 0; it < 100; it++) {
		double f = s - (double)length, v, t_next;
		if (fabs(f) <= (double)QAWS_EPSILON * (double)total)
			break;
		if (f < 0)
			lo = t;
		else
			hi = t;
		if (hi - lo <= (double)QAWS_EPSILON * ((double)t_end - (double)t_start))
			break;
		v = curve_speed(curve, (qaws_scalar)t);
		t_next = v > 0 ? t - f / v : lo - 1;
		if (!(t_next > lo && t_next < hi) || fabs(2 * f) > fabs(step_old * v)) {
			step_old = step;
			t_next = 0.5 * (lo + hi);
		} else {
			step_old = step;
		}
		step = t_next - t;
		if (t_next > t) {
			qaws_internal_arc_length(&piece, curve, (qaws_scalar)t, (qaws_scalar)t_next);
			s += piece;
		} else {
			qaws_internal_arc_length(&piece, curve, (qaws_scalar)t_next, (qaws_scalar)t);
			s -= piece;
		}
		t = t_next;
	}
	if (it == 100)
		return QAWS_STATUS_NOT_CONVERGED;
	*out_parameter = (qaws_scalar)t;
	return QAWS_STATUS_OK;
}
