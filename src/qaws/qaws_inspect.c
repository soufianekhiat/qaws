#include "qaws_alloc.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_inspect.h"
#include "qaws_eval.h"
#include "qaws_surface.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_closest.h"
#include "internal/qaws_internal_solver.h"
#include "internal/qaws_internal_fit.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  Generic inspection functions                                       */
/* ------------------------------------------------------------------ */

qaws_curve_kind qaws_curve_get_kind(qaws_curve const *curve)
{
	if (!curve)
		return QAWS_CURVE_KIND_INVALID;
	return curve->kind;
}

qaws_dimension qaws_curve_get_dimension(qaws_curve const *curve)
{
	if (!curve)
		return (qaws_dimension)0;
	return curve->dimension;
}

unsigned int qaws_curve_get_degree(qaws_curve const *curve)
{
	if (!curve)
		return 0;
	return curve->degree;
}

unsigned int qaws_curve_get_span_count(qaws_curve const *curve)
{
	if (!curve)
		return 0;
	return curve->span_count;
}

qaws_range qaws_curve_get_parameter_range(qaws_curve const *curve)
{
	qaws_range range;
	if (!curve)
	{
		range.min_value = (qaws_scalar)0.0;
		range.max_value = (qaws_scalar)0.0;
		return range;
	}
	return curve->parameter_range;
}

int qaws_curve_is_closed(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->is_closed)
		return 0;
	return curve->vtable->is_closed(curve);
}

int qaws_curve_is_periodic(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->is_periodic)
		return 0;
	return curve->vtable->is_periodic(curve);
}

int qaws_curve_is_rational(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->is_rational)
		return 0;
	return curve->vtable->is_rational(curve);
}

qaws_continuity qaws_curve_get_continuity(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->get_continuity)
		return QAWS_CONTINUITY_C0;
	return curve->vtable->get_continuity(curve);
}

/* ------------------------------------------------------------------ */
/*  Analysis helpers                                                    */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_arc_length(
	qaws_scalar *out_length,
	qaws_curve const *curve,
	qaws_scalar parameter_min,
	qaws_scalar parameter_max)
{
	if (!curve || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;

	return qaws_internal_arc_length(
		out_length, curve, parameter_min, parameter_max);
}

qaws_status qaws_curve_compute_parameter_at_arc_length(
	qaws_scalar *out_parameter,
	qaws_curve const *curve,
	qaws_scalar parameter_min,
	qaws_scalar length)
{
	return qaws_internal_parameter_at_arc_length(
		out_parameter, curve, parameter_min, length);
}

qaws_status qaws_curve_compute_bounds_2d(
	qaws_vec2 *out_min,
	qaws_vec2 *out_max,
	qaws_curve const *curve)
{
	unsigned int sample_count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_2d result;
	qaws_status status;

	if (!curve || !out_min || !out_max)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	sample_count = 64 * curve->span_count;
	if (sample_count < 64)
		sample_count = 64;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Evaluate first point to initialize bounds */
	status = qaws_curve_evaluate_2d(
		&result, curve, range_min, QAWS_EVAL_FLAG_POSITION);
	if (status != QAWS_STATUS_OK)
		return status;

	out_min->x = result.position.x;
	out_min->y = result.position.y;
	out_max->x = result.position.x;
	out_max->y = result.position.y;

	for (i = 1; i < sample_count; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_2d(
			&result, curve, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		if (result.position.x < out_min->x)
			out_min->x = result.position.x;
		if (result.position.y < out_min->y)
			out_min->y = result.position.y;
		if (result.position.x > out_max->x)
			out_max->x = result.position.x;
		if (result.position.y > out_max->y)
			out_max->y = result.position.y;
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_bounds_3d(
	qaws_vec3 *out_min,
	qaws_vec3 *out_max,
	qaws_curve const *curve)
{
	unsigned int sample_count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_3d result;
	qaws_status status;

	if (!curve || !out_min || !out_max)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	sample_count = 64 * curve->span_count;
	if (sample_count < 64)
		sample_count = 64;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Evaluate first point to initialize bounds */
	status = qaws_curve_evaluate_3d(
		&result, curve, range_min, QAWS_EVAL_FLAG_POSITION);
	if (status != QAWS_STATUS_OK)
		return status;

	out_min->x = result.position.x;
	out_min->y = result.position.y;
	out_min->z = result.position.z;
	out_max->x = result.position.x;
	out_max->y = result.position.y;
	out_max->z = result.position.z;

	for (i = 1; i < sample_count; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_3d(
			&result, curve, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		if (result.position.x < out_min->x)
			out_min->x = result.position.x;
		if (result.position.y < out_min->y)
			out_min->y = result.position.y;
		if (result.position.z < out_min->z)
			out_min->z = result.position.z;
		if (result.position.x > out_max->x)
			out_max->x = result.position.x;
		if (result.position.y > out_max->y)
			out_max->y = result.position.y;
		if (result.position.z > out_max->z)
			out_max->z = result.position.z;
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Closest-point helpers                                              */
/* ------------------------------------------------------------------ */

struct closest_ctx_2d
{
	qaws_curve const *curve;
};

static void closest_eval_2d(
	void const *ctx,
	qaws_scalar t,
	qaws_scalar *out_pos,
	qaws_scalar *out_d1,
	unsigned int dim_count)
{
	struct closest_ctx_2d const *c = (struct closest_ctx_2d const *)ctx;
	qaws_eval_result_2d result;

	(void)dim_count;

	qaws_curve_evaluate_2d(&result, c->curve,
		t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);

	out_pos[0] = result.position.x;
	out_pos[1] = result.position.y;
	out_d1[0] = result.d1.x;
	out_d1[1] = result.d1.y;
}

static qaws_status closest_sampled_2d(
	qaws_curve const *curve,
	qaws_vec2 point,
	qaws_scalar *out_parameter)
{
	unsigned int i;
	qaws_scalar t;
	qaws_scalar best_t;
	qaws_scalar best_dist_sq;
	qaws_scalar dx;
	qaws_scalar dy;
	qaws_scalar dist_sq;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_2d result;
	qaws_status status;
	struct closest_ctx_2d ctx;
	qaws_scalar target_pos[2];

	if (!curve || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Initial coarse search: sample at 32 points */
	best_t = range_min;
	best_dist_sq = (qaws_scalar)1.0e30;

	for (i = 0; i < 32; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)31.0;

		status = qaws_curve_evaluate_2d(
			&result, curve, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result.position.x - point.x;
		dy = result.position.y - point.y;
		dist_sq = dx * dx + dy * dy;

		if (dist_sq < best_dist_sq)
		{
			best_dist_sq = dist_sq;
			best_t = t;
		}
	}

	/* Refine with Newton iteration */
	ctx.curve = curve;
	target_pos[0] = point.x;
	target_pos[1] = point.y;

	best_t = qaws_internal_closest_point_newton(
		target_pos, 2, closest_eval_2d, range_min,
		range_max, best_t, 32, &ctx);

	*out_parameter = best_t;
	return QAWS_STATUS_OK;
}

/* the sampled answer, then certified where the family has Bezier spans */
static qaws_status closest_2d(
	qaws_curve const *curve,
	qaws_vec2 point,
	qaws_scalar *out_parameter,
	double *out_distance,
	double *out_gap)
{
	double x[3] = { point.x, point.y, 0 };
	qaws_status status = closest_sampled_2d(curve, point, out_parameter);
	if (status != QAWS_STATUS_OK)
		return status;
	return qaws_internal_closest_certify(out_distance, out_gap, out_parameter, curve, x);
}

qaws_status qaws_curve_find_closest_parameter_2d(
	qaws_scalar *out_parameter,
	qaws_curve const *curve,
	qaws_vec2 point)
{
	double distance, gap;
	qaws_status status = closest_2d(curve, point, out_parameter, &distance, &gap);
	return status == QAWS_STATUS_UNSUPPORTED_OPERATION ? QAWS_STATUS_OK : status;
}

qaws_status qaws_curve_find_closest_point_certified_2d(
	qaws_scalar *out_parameter,
	qaws_scalar *out_distance,
	qaws_scalar *out_gap,
	qaws_curve const *curve,
	qaws_vec2 point)
{
	double distance, gap;
	qaws_scalar t;
	qaws_status status;
	if (!out_parameter || !out_distance || !out_gap)
		return QAWS_STATUS_INVALID_ARGUMENT;
	status = closest_2d(curve, point, &t, &distance, &gap);
	if (status != QAWS_STATUS_OK)
		return status;
	*out_parameter = t;
	*out_distance = (qaws_scalar)distance;
	*out_gap = (qaws_scalar)gap;
	return QAWS_STATUS_OK;
}

struct closest_ctx_3d
{
	qaws_curve const *curve;
};

static void closest_eval_3d(
	void const *ctx,
	qaws_scalar t,
	qaws_scalar *out_pos,
	qaws_scalar *out_d1,
	unsigned int dim_count)
{
	struct closest_ctx_3d const *c = (struct closest_ctx_3d const *)ctx;
	qaws_eval_result_3d result;

	(void)dim_count;

	qaws_curve_evaluate_3d(&result, c->curve,
		t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);

	out_pos[0] = result.position.x;
	out_pos[1] = result.position.y;
	out_pos[2] = result.position.z;
	out_d1[0] = result.d1.x;
	out_d1[1] = result.d1.y;
	out_d1[2] = result.d1.z;
}

static qaws_status closest_sampled_3d(
	qaws_curve const *curve,
	qaws_vec3 point,
	qaws_scalar *out_parameter)
{
	unsigned int i;
	qaws_scalar t;
	qaws_scalar best_t;
	qaws_scalar best_dist_sq;
	qaws_scalar dx;
	qaws_scalar dy;
	qaws_scalar dz;
	qaws_scalar dist_sq;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_3d result;
	qaws_status status;
	struct closest_ctx_3d ctx;
	qaws_scalar target_pos[3];

	if (!curve || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Initial coarse search: sample at 32 points */
	best_t = range_min;
	best_dist_sq = (qaws_scalar)1.0e30;

	for (i = 0; i < 32; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)31.0;

		status = qaws_curve_evaluate_3d(
			&result, curve, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result.position.x - point.x;
		dy = result.position.y - point.y;
		dz = result.position.z - point.z;
		dist_sq = dx * dx + dy * dy + dz * dz;

		if (dist_sq < best_dist_sq)
		{
			best_dist_sq = dist_sq;
			best_t = t;
		}
	}

	/* Refine with Newton iteration */
	ctx.curve = curve;
	target_pos[0] = point.x;
	target_pos[1] = point.y;
	target_pos[2] = point.z;

	best_t = qaws_internal_closest_point_newton(
		target_pos, 3, closest_eval_3d, range_min,
		range_max, best_t, 32, &ctx);

	*out_parameter = best_t;
	return QAWS_STATUS_OK;
}

/* the sampled answer, then certified where the family has Bezier spans */
static qaws_status closest_3d(
	qaws_curve const *curve,
	qaws_vec3 point,
	qaws_scalar *out_parameter,
	double *out_distance,
	double *out_gap)
{
	double x[3] = { point.x, point.y, point.z };
	qaws_status status = closest_sampled_3d(curve, point, out_parameter);
	if (status != QAWS_STATUS_OK)
		return status;
	return qaws_internal_closest_certify(out_distance, out_gap, out_parameter, curve, x);
}

qaws_status qaws_curve_find_closest_parameter_3d(
	qaws_scalar *out_parameter,
	qaws_curve const *curve,
	qaws_vec3 point)
{
	double distance, gap;
	qaws_status status = closest_3d(curve, point, out_parameter, &distance, &gap);
	return status == QAWS_STATUS_UNSUPPORTED_OPERATION ? QAWS_STATUS_OK : status;
}

qaws_status qaws_curve_find_closest_point_certified_3d(
	qaws_scalar *out_parameter,
	qaws_scalar *out_distance,
	qaws_scalar *out_gap,
	qaws_curve const *curve,
	qaws_vec3 point)
{
	double distance, gap;
	qaws_scalar t;
	qaws_status status;
	if (!out_parameter || !out_distance || !out_gap)
		return QAWS_STATUS_INVALID_ARGUMENT;
	status = closest_3d(curve, point, &t, &distance, &gap);
	if (status != QAWS_STATUS_OK)
		return status;
	*out_parameter = t;
	*out_distance = (qaws_scalar)distance;
	*out_gap = (qaws_scalar)gap;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Span continuity                                                    */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_get_span_continuity(
	qaws_continuity *out_continuity,
	qaws_curve const *curve,
	unsigned int boundary_index)
{
	if (!curve || !out_continuity)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->span_count == 0 || boundary_index >= curve->span_count - 1)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/*
	 * For v1, return the overall curve continuity.
	 * A more precise implementation would check actual derivatives
	 * at the span boundary.
	 */
	if (!curve->vtable || !curve->vtable->get_continuity)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_continuity = curve->vtable->get_continuity(curve);
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Family-specific inspection: Bezier                                 */
/* ------------------------------------------------------------------ */

qaws_status qaws_bezier_get_control_points(
	void *out_control_points,
	unsigned int *out_point_count,
	qaws_curve const *curve,
	unsigned int point_capacity)
{
	qaws_bezier_impl const *impl;
	unsigned int count;
	unsigned int dim_count;

	if (!curve || !out_control_points || !out_point_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->kind != QAWS_CURVE_KIND_BEZIER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bezier_impl const *)curve->impl;
	count = impl->control_point_count;
	dim_count = (unsigned int)curve->dimension;

	*out_point_count = count;

	if (point_capacity < count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	memcpy(out_control_points, impl->control_points,
		(size_t)count * (size_t)dim_count * sizeof(qaws_scalar));

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_get_control_points(
	void *out_control_points,
	unsigned int *out_point_count,
	qaws_curve const *curve,
	unsigned int point_capacity)
{
	qaws_scalar const *cps;
	unsigned int count;

	if (!curve || !out_control_points || !out_point_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	switch (curve->kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
		cps = ((qaws_bezier_impl const *)curve->impl)->control_points;
		count = ((qaws_bezier_impl const *)curve->impl)->control_point_count;
		break;
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
		cps = ((qaws_rational_bezier_impl const *)curve->impl)->control_points;
		count = ((qaws_rational_bezier_impl const *)curve->impl)->control_point_count;
		break;
	case QAWS_CURVE_KIND_BSPLINE:
		cps = ((qaws_bspline_impl const *)curve->impl)->control_points;
		count = ((qaws_bspline_impl const *)curve->impl)->control_point_count;
		break;
	case QAWS_CURVE_KIND_NURBS:
		cps = ((qaws_nurbs_impl const *)curve->impl)->control_points;
		count = ((qaws_nurbs_impl const *)curve->impl)->control_point_count;
		break;
	default:
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	}

	*out_point_count = count;
	if (point_capacity < count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	memcpy(out_control_points, cps,
		(size_t)count * (size_t)curve->dimension * sizeof(qaws_scalar));
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Family-specific inspection: B-Spline knots                         */
/* ------------------------------------------------------------------ */

qaws_status qaws_bspline_get_knots(
	qaws_scalar *out_knots,
	unsigned int *out_knot_count,
	qaws_curve const *curve,
	unsigned int knot_capacity)
{
	if (!curve || !out_knots || !out_knot_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->kind == QAWS_CURVE_KIND_BSPLINE)
	{
		qaws_bspline_impl const *impl =
			(qaws_bspline_impl const *)curve->impl;
		unsigned int knot_count = impl->knot_count;

		*out_knot_count = knot_count;

		if (knot_capacity < knot_count)
			return QAWS_STATUS_BUFFER_TOO_SMALL;

		memcpy(out_knots, impl->knots,
			(size_t)knot_count * sizeof(qaws_scalar));

		return QAWS_STATUS_OK;
	}
	else if (curve->kind == QAWS_CURVE_KIND_NURBS)
	{
		qaws_nurbs_impl const *impl =
			(qaws_nurbs_impl const *)curve->impl;
		unsigned int knot_count = impl->knot_count;

		*out_knot_count = knot_count;

		if (knot_capacity < knot_count)
			return QAWS_STATUS_BUFFER_TOO_SMALL;

		memcpy(out_knots, impl->knots,
			(size_t)knot_count * sizeof(qaws_scalar));

		return QAWS_STATUS_OK;
	}

	return QAWS_STATUS_UNSUPPORTED_OPERATION;
}

/* ------------------------------------------------------------------ */
/*  Family-specific inspection: NURBS weights                          */
/* ------------------------------------------------------------------ */

qaws_status qaws_nurbs_get_weights(
	qaws_scalar *out_weights,
	unsigned int *out_weight_count,
	qaws_curve const *curve,
	unsigned int weight_capacity)
{
	qaws_nurbs_impl const *impl;
	unsigned int weight_count;

	if (!curve || !out_weights || !out_weight_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->kind != QAWS_CURVE_KIND_NURBS)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_nurbs_impl const *)curve->impl;
	weight_count = impl->weight_count;

	*out_weight_count = weight_count;

	if (weight_capacity < weight_count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	memcpy(out_weights, impl->weights,
		(size_t)weight_count * sizeof(qaws_scalar));

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Derived geometric helpers                                          */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_tangent_2d(
	qaws_vec2 *out_tangent,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		&result, curve, parameter, QAWS_EVAL_FLAG_D1);
	if (status != QAWS_STATUS_OK)
		return status;

	len = (qaws_scalar)sqrt(
		(double)(result.d1.x * result.d1.x + result.d1.y * result.d1.y));

	if (len < (qaws_scalar)1.0e-12)
	{
		out_tangent->x = (qaws_scalar)0.0;
		out_tangent->y = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	out_tangent->x = result.d1.x / len;
	out_tangent->y = result.d1.y / len;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_tangent_3d(
	qaws_vec3 *out_tangent,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		&result, curve, parameter, QAWS_EVAL_FLAG_D1);
	if (status != QAWS_STATUS_OK)
		return status;

	len = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x +
		result.d1.y * result.d1.y +
		result.d1.z * result.d1.z));

	if (len < (qaws_scalar)1.0e-12)
	{
		out_tangent->x = (qaws_scalar)0.0;
		out_tangent->y = (qaws_scalar)0.0;
		out_tangent->z = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	out_tangent->x = result.d1.x / len;
	out_tangent->y = result.d1.y / len;
	out_tangent->z = result.d1.z / len;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_curvature_2d(
	qaws_scalar *out_curvature,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar speed;
	qaws_scalar cross;

	if (!curve || !out_curvature)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		&result, curve, parameter, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
	if (status != QAWS_STATUS_OK)
		return status;

	speed = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x + result.d1.y * result.d1.y));

	if (speed < (qaws_scalar)1.0e-12)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* Signed curvature: (d1.x * d2.y - d1.y * d2.x) / |d1|^3 */
	cross = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	*out_curvature = cross / (speed * speed * speed);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_curvature_3d(
	qaws_scalar *out_curvature,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar speed;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len;

	if (!curve || !out_curvature)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		&result, curve, parameter, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
	if (status != QAWS_STATUS_OK)
		return status;

	speed = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x +
		result.d1.y * result.d1.y +
		result.d1.z * result.d1.z));

	if (speed < (qaws_scalar)1.0e-12)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* |d1 x d2| / |d1|^3 */
	cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
	cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
	cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	cross_len = (qaws_scalar)sqrt((double)(cx * cx + cy * cy + cz * cz));

	*out_curvature = cross_len / (speed * speed * speed);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_torsion_3d(
	qaws_scalar *out_torsion,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len_sq;
	qaws_scalar dot;

	if (!curve || !out_torsion)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		&result, curve,
		parameter,
		QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3);
	if (status != QAWS_STATUS_OK)
		return status;

	/* (d1 x d2) . d3 / |d1 x d2|^2 */
	cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
	cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
	cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	cross_len_sq = cx * cx + cy * cy + cz * cz;

	if (cross_len_sq < (qaws_scalar)1.0e-24)
	{
		*out_torsion = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	dot = cx * result.d3.x + cy * result.d3.y + cz * result.d3.z;
	*out_torsion = dot / cross_len_sq;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_speed(
	qaws_scalar *out_speed,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_status status;

	if (!curve || !out_speed)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d result;
		status = qaws_curve_evaluate_2d(
			&result, curve, parameter, QAWS_EVAL_FLAG_D1);
		if (status != QAWS_STATUS_OK)
			return status;

		*out_speed = (qaws_scalar)sqrt((double)(
			result.d1.x * result.d1.x +
			result.d1.y * result.d1.y));
	}
	else
	{
		qaws_eval_result_3d result;
		status = qaws_curve_evaluate_3d(
			&result, curve, parameter, QAWS_EVAL_FLAG_D1);
		if (status != QAWS_STATUS_OK)
			return status;

		*out_speed = (qaws_scalar)sqrt((double)(
			result.d1.x * result.d1.x +
			result.d1.y * result.d1.y +
			result.d1.z * result.d1.z));
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Frenet frame                                                       */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_normal_2d(
	qaws_vec2 *out_normal,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_normal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		&result, curve, parameter, QAWS_EVAL_FLAG_D1);
	if (status != QAWS_STATUS_OK)
		return status;

	len = (qaws_scalar)sqrt(
		(double)(result.d1.x * result.d1.x + result.d1.y * result.d1.y));

	if (len < (qaws_scalar)1.0e-12)
	{
		out_normal->x = (qaws_scalar)0.0;
		out_normal->y = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* Rotate normalized tangent 90 degrees CCW: (-ty, tx) */
	out_normal->x = -result.d1.y / len;
	out_normal->y =  result.d1.x / len;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_frenet_frame_3d(
	qaws_vec3 *out_tangent,
	qaws_vec3 *out_normal,
	qaws_vec3 *out_binormal,
	qaws_curve const *curve,
	qaws_scalar parameter)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar d1_len;
	qaws_scalar tx, ty, tz;
	qaws_scalar bx, by, bz;
	qaws_scalar b_len;
	qaws_scalar ax, ay, az;

	if (!curve || !out_tangent || !out_normal || !out_binormal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		&result, curve,
		parameter, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
	if (status != QAWS_STATUS_OK)
		return status;

	/* T = normalize(D1) */
	d1_len = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x +
		result.d1.y * result.d1.y +
		result.d1.z * result.d1.z));

	if (d1_len < (qaws_scalar)1.0e-12)
	{
		out_tangent->x  = (qaws_scalar)0.0;
		out_tangent->y  = (qaws_scalar)0.0;
		out_tangent->z  = (qaws_scalar)0.0;
		out_normal->x   = (qaws_scalar)0.0;
		out_normal->y   = (qaws_scalar)0.0;
		out_normal->z   = (qaws_scalar)0.0;
		out_binormal->x = (qaws_scalar)0.0;
		out_binormal->y = (qaws_scalar)0.0;
		out_binormal->z = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	tx = result.d1.x / d1_len;
	ty = result.d1.y / d1_len;
	tz = result.d1.z / d1_len;

	/* B = normalize(D1 x D2) */
	bx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
	by = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
	bz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	b_len = (qaws_scalar)sqrt((double)(bx * bx + by * by + bz * bz));

	if (b_len < (qaws_scalar)1.0e-12)
	{
		/*
		 * D1 x D2 is near zero (straight line or inflection).
		 * Pick an arbitrary vector not parallel to T, then
		 * use cross products to build a perpendicular normal.
		 */
		ax = (qaws_scalar)0.0;
		ay = (qaws_scalar)0.0;
		az = (qaws_scalar)0.0;

		/* Choose the axis least aligned with T */
		if (tx * tx <= ty * ty && tx * tx <= tz * tz)
			ax = (qaws_scalar)1.0;
		else if (ty * ty <= tz * tz)
			ay = (qaws_scalar)1.0;
		else
			az = (qaws_scalar)1.0;

		/* B = normalize(T x arbitrary) */
		bx = ty * az - tz * ay;
		by = tz * ax - tx * az;
		bz = tx * ay - ty * ax;
		b_len = (qaws_scalar)sqrt((double)(bx * bx + by * by + bz * bz));

		bx /= b_len;
		by /= b_len;
		bz /= b_len;
	}
	else
	{
		bx /= b_len;
		by /= b_len;
		bz /= b_len;
	}

	/* N = B x T */
	out_tangent->x = tx;
	out_tangent->y = ty;
	out_tangent->z = tz;

	out_normal->x = by * tz - bz * ty;
	out_normal->y = bz * tx - bx * tz;
	out_normal->z = bx * ty - by * tx;

	out_binormal->x = bx;
	out_binormal->y = by;
	out_binormal->z = bz;

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Inflection point detection                                         */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_find_inflection_points(
	qaws_scalar *out_parameters,
	unsigned int *out_count,
	qaws_curve const *curve,
	unsigned int parameter_capacity)
{
	unsigned int span_count;
	unsigned int samples_per_span;
	unsigned int total_samples;
	unsigned int found;
	unsigned int i;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_scalar prev_val;
	qaws_scalar prev_t;

	if (!curve || !out_parameters || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_count = 0;
	found = 0;

	span_count = curve->span_count;
	if (span_count == 0)
		span_count = 1;
	samples_per_span = 32;
	total_samples = span_count * samples_per_span;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d result;
		qaws_status status;
		qaws_scalar t;
		qaws_scalar val;

		/* Evaluate first sample */
		status = qaws_curve_evaluate_2d(
			&result, curve,
			range_min, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		if (status != QAWS_STATUS_OK)
			return status;

		prev_val = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
		prev_t = range_min;

		for (i = 1; i <= total_samples; ++i)
		{
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)total_samples;

			status = qaws_curve_evaluate_2d(
				&result, curve,
				t, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
			if (status != QAWS_STATUS_OK)
				return status;

			val = result.d1.x * result.d2.y - result.d1.y * result.d2.x;

			/* Sign change detected - bisect */
			if ((prev_val > (qaws_scalar)0.0 && val < (qaws_scalar)0.0) ||
				(prev_val < (qaws_scalar)0.0 && val > (qaws_scalar)0.0))
			{
				qaws_scalar lo = prev_t;
				qaws_scalar hi = t;
				qaws_scalar lo_val = prev_val;
				qaws_scalar mid;
				qaws_scalar mid_val;
				unsigned int iter;

				for (iter = 0; iter < 32; ++iter)
				{
					mid = (lo + hi) * (qaws_scalar)0.5;
					if ((hi - lo) < (qaws_scalar)1.0e-12)
						break;

					status = qaws_curve_evaluate_2d(
						&result, curve,
						mid, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
					if (status != QAWS_STATUS_OK)
						return status;

					mid_val = result.d1.x * result.d2.y
						- result.d1.y * result.d2.x;

					if ((lo_val > (qaws_scalar)0.0 && mid_val > (qaws_scalar)0.0) ||
						(lo_val < (qaws_scalar)0.0 && mid_val < (qaws_scalar)0.0))
					{
						lo = mid;
						lo_val = mid_val;
					}
					else
					{
						hi = mid;
					}
				}

				mid = (lo + hi) * (qaws_scalar)0.5;
				if (found < parameter_capacity)
					out_parameters[found] = mid;
				++found;
			}

			prev_val = val;
			prev_t = t;
		}
	}
	else if (curve->dimension == QAWS_DIMENSION_3D)
	{
		qaws_eval_result_3d result;
		qaws_status status;
		qaws_scalar t;
		qaws_scalar cx, cy, cz;
		qaws_scalar val;
		qaws_scalar prev_mag;

		/* Evaluate first sample */
		status = qaws_curve_evaluate_3d(
			&result, curve,
			range_min, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		if (status != QAWS_STATUS_OK)
			return status;

		cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
		cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
		cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
		prev_mag = (qaws_scalar)sqrt((double)(cx * cx + cy * cy + cz * cz));
		prev_t = range_min;

		for (i = 1; i <= total_samples; ++i)
		{
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)total_samples;

			status = qaws_curve_evaluate_3d(
				&result, curve,
				t, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
			if (status != QAWS_STATUS_OK)
				return status;

			cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
			cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
			cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
			val = (qaws_scalar)sqrt((double)(cx * cx + cy * cy + cz * cz));

			/*
			 * For 3D, detect when magnitude approaches zero.
			 * If one side is above tolerance and the other is below,
			 * or if both bracket a minimum near zero, bisect.
			 */
			if ((prev_mag > (qaws_scalar)1.0e-12 && val < (qaws_scalar)1.0e-12) ||
				(prev_mag < (qaws_scalar)1.0e-12 && val > (qaws_scalar)1.0e-12))
			{
				qaws_scalar lo = prev_t;
				qaws_scalar hi = t;
				qaws_scalar mid;
				qaws_scalar mid_mag;
				unsigned int iter;

				for (iter = 0; iter < 32; ++iter)
				{
					mid = (lo + hi) * (qaws_scalar)0.5;
					if ((hi - lo) < (qaws_scalar)1.0e-12)
						break;

					status = qaws_curve_evaluate_3d(
						&result, curve,
						mid, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
					if (status != QAWS_STATUS_OK)
						return status;

					cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
					cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
					cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
					mid_mag = (qaws_scalar)sqrt(
						(double)(cx * cx + cy * cy + cz * cz));

					if (mid_mag < (qaws_scalar)1.0e-12)
					{
						/* Found the zero region */
						break;
					}

					/* Bisect toward the side with smaller magnitude */
					if (prev_mag < val)
						hi = mid;
					else
						lo = mid;
				}

				mid = (lo + hi) * (qaws_scalar)0.5;
				if (found < parameter_capacity)
					out_parameters[found] = mid;
				++found;
			}

			prev_mag = val;
			prev_t = t;
		}
	}
	else
	{
		return QAWS_STATUS_INVALID_DIMENSION;
	}

	*out_count = found;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Extrema detection                                                  */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_find_extrema(
	qaws_scalar *out_parameters,
	unsigned int *out_count,
	qaws_curve const *curve,
	unsigned int parameter_capacity,
	unsigned int axis)
{
	unsigned int span_count;
	unsigned int samples_per_span;
	unsigned int total_samples;
	unsigned int found;
	unsigned int i;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_scalar prev_val;
	qaws_scalar prev_t;

	if (!curve || !out_parameters || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (axis >= (unsigned int)curve->dimension)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_count = 0;
	found = 0;

	span_count = curve->span_count;
	if (span_count == 0)
		span_count = 1;
	samples_per_span = 32;
	total_samples = span_count * samples_per_span;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d result;
		qaws_status status;
		qaws_scalar t;
		qaws_scalar val;

		/* Evaluate first sample */
		status = qaws_curve_evaluate_2d(
			&result, curve, range_min, QAWS_EVAL_FLAG_D1);
		if (status != QAWS_STATUS_OK)
			return status;

		prev_val = (axis == 0) ? result.d1.x : result.d1.y;
		prev_t = range_min;

		for (i = 1; i <= total_samples; ++i)
		{
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)total_samples;

			status = qaws_curve_evaluate_2d(
				&result, curve, t, QAWS_EVAL_FLAG_D1);
			if (status != QAWS_STATUS_OK)
				return status;

			val = (axis == 0) ? result.d1.x : result.d1.y;

			/* Sign change detected - bisect */
			if ((prev_val > (qaws_scalar)0.0 && val < (qaws_scalar)0.0) ||
				(prev_val < (qaws_scalar)0.0 && val > (qaws_scalar)0.0))
			{
				qaws_scalar lo = prev_t;
				qaws_scalar hi = t;
				qaws_scalar lo_val = prev_val;
				qaws_scalar mid;
				qaws_scalar mid_val;
				unsigned int iter;

				for (iter = 0; iter < 32; ++iter)
				{
					mid = (lo + hi) * (qaws_scalar)0.5;
					if ((hi - lo) < (qaws_scalar)1.0e-12)
						break;

					status = qaws_curve_evaluate_2d(
						&result, curve, mid, QAWS_EVAL_FLAG_D1);
					if (status != QAWS_STATUS_OK)
						return status;

					mid_val = (axis == 0) ? result.d1.x : result.d1.y;

					if ((lo_val > (qaws_scalar)0.0 && mid_val > (qaws_scalar)0.0) ||
						(lo_val < (qaws_scalar)0.0 && mid_val < (qaws_scalar)0.0))
					{
						lo = mid;
						lo_val = mid_val;
					}
					else
					{
						hi = mid;
					}
				}

				mid = (lo + hi) * (qaws_scalar)0.5;
				if (found < parameter_capacity)
					out_parameters[found] = mid;
				++found;
			}

			prev_val = val;
			prev_t = t;
		}
	}
	else if (curve->dimension == QAWS_DIMENSION_3D)
	{
		qaws_eval_result_3d result;
		qaws_status status;
		qaws_scalar t;
		qaws_scalar val;

		/* Evaluate first sample */
		status = qaws_curve_evaluate_3d(
			&result, curve, range_min, QAWS_EVAL_FLAG_D1);
		if (status != QAWS_STATUS_OK)
			return status;

		if (axis == 0) prev_val = result.d1.x;
		else if (axis == 1) prev_val = result.d1.y;
		else prev_val = result.d1.z;
		prev_t = range_min;

		for (i = 1; i <= total_samples; ++i)
		{
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)total_samples;

			status = qaws_curve_evaluate_3d(
				&result, curve, t, QAWS_EVAL_FLAG_D1);
			if (status != QAWS_STATUS_OK)
				return status;

			if (axis == 0) val = result.d1.x;
			else if (axis == 1) val = result.d1.y;
			else val = result.d1.z;

			/* Sign change detected - bisect */
			if ((prev_val > (qaws_scalar)0.0 && val < (qaws_scalar)0.0) ||
				(prev_val < (qaws_scalar)0.0 && val > (qaws_scalar)0.0))
			{
				qaws_scalar lo = prev_t;
				qaws_scalar hi = t;
				qaws_scalar lo_val = prev_val;
				qaws_scalar mid;
				qaws_scalar mid_val;
				unsigned int iter;

				for (iter = 0; iter < 32; ++iter)
				{
					mid = (lo + hi) * (qaws_scalar)0.5;
					if ((hi - lo) < (qaws_scalar)1.0e-12)
						break;

					status = qaws_curve_evaluate_3d(
						&result, curve, mid, QAWS_EVAL_FLAG_D1);
					if (status != QAWS_STATUS_OK)
						return status;

					if (axis == 0) mid_val = result.d1.x;
					else if (axis == 1) mid_val = result.d1.y;
					else mid_val = result.d1.z;

					if ((lo_val > (qaws_scalar)0.0 && mid_val > (qaws_scalar)0.0) ||
						(lo_val < (qaws_scalar)0.0 && mid_val < (qaws_scalar)0.0))
					{
						lo = mid;
						lo_val = mid_val;
					}
					else
					{
						hi = mid;
					}
				}

				mid = (lo + hi) * (qaws_scalar)0.5;
				if (found < parameter_capacity)
					out_parameters[found] = mid;
				++found;
			}

			prev_val = val;
			prev_t = t;
		}
	}
	else
	{
		return QAWS_STATUS_INVALID_DIMENSION;
	}

	*out_count = found;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Curvature comb data                                                */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_curvature_comb_2d(
	qaws_curvature_sample_2d *out_samples,
	qaws_curve const *curve,
	unsigned int sample_count,
	unsigned int sample_capacity)
{
	unsigned int i;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_scalar t;
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar speed;
	qaws_scalar cross;
	qaws_scalar len;

	if (!curve || !out_samples)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (sample_capacity < sample_count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	for (i = 0; i < sample_count; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_2d(
			&result, curve,
			t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		if (status != QAWS_STATUS_OK)
			return status;

		out_samples[i].position = result.position;

		/* Curvature: (d1.x * d2.y - d1.y * d2.x) / |d1|^3 */
		speed = (qaws_scalar)sqrt((double)(
			result.d1.x * result.d1.x + result.d1.y * result.d1.y));

		if (speed < (qaws_scalar)1.0e-12)
		{
			out_samples[i].curvature = (qaws_scalar)0.0;
			out_samples[i].normal.x = (qaws_scalar)0.0;
			out_samples[i].normal.y = (qaws_scalar)0.0;
			continue;
		}

		cross = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
		out_samples[i].curvature = cross / (speed * speed * speed);

		/* Normal: rotate unit tangent 90 degrees CCW */
		len = speed;
		out_samples[i].normal.x = -result.d1.y / len;
		out_samples[i].normal.y =  result.d1.x / len;
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_curvature_comb_3d(
	qaws_curvature_sample_3d *out_samples,
	qaws_curve const *curve,
	unsigned int sample_count,
	unsigned int sample_capacity)
{
	unsigned int i;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_scalar t;
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar speed;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len;
	qaws_scalar d1_len;
	qaws_scalar tx, ty, tz;
	qaws_scalar bx, by, bz;
	qaws_scalar b_len;
	qaws_scalar ax, ay, az;

	if (!curve || !out_samples)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (sample_capacity < sample_count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	for (i = 0; i < sample_count; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_3d(
			&result, curve,
			t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		if (status != QAWS_STATUS_OK)
			return status;

		out_samples[i].position = result.position;

		/* Curvature: |d1 x d2| / |d1|^3 */
		speed = (qaws_scalar)sqrt((double)(
			result.d1.x * result.d1.x +
			result.d1.y * result.d1.y +
			result.d1.z * result.d1.z));

		if (speed < (qaws_scalar)1.0e-12)
		{
			out_samples[i].curvature = (qaws_scalar)0.0;
			out_samples[i].normal.x = (qaws_scalar)0.0;
			out_samples[i].normal.y = (qaws_scalar)0.0;
			out_samples[i].normal.z = (qaws_scalar)0.0;
			continue;
		}

		cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
		cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
		cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
		cross_len = (qaws_scalar)sqrt(
			(double)(cx * cx + cy * cy + cz * cz));

		out_samples[i].curvature = cross_len / (speed * speed * speed);

		/* Normal from Frenet frame: N = B x T where B = normalize(D1 x D2) */
		d1_len = speed;
		tx = result.d1.x / d1_len;
		ty = result.d1.y / d1_len;
		tz = result.d1.z / d1_len;

		b_len = cross_len;

		if (b_len < (qaws_scalar)1.0e-12)
		{
			/* Degenerate: pick arbitrary perpendicular normal */
			ax = (qaws_scalar)0.0;
			ay = (qaws_scalar)0.0;
			az = (qaws_scalar)0.0;

			if (tx * tx <= ty * ty && tx * tx <= tz * tz)
				ax = (qaws_scalar)1.0;
			else if (ty * ty <= tz * tz)
				ay = (qaws_scalar)1.0;
			else
				az = (qaws_scalar)1.0;

			bx = ty * az - tz * ay;
			by = tz * ax - tx * az;
			bz = tx * ay - ty * ax;
			b_len = (qaws_scalar)sqrt(
				(double)(bx * bx + by * by + bz * bz));

			bx /= b_len;
			by /= b_len;
			bz /= b_len;
		}
		else
		{
			bx = cx / b_len;
			by = cy / b_len;
			bz = cz / b_len;
		}

		/* N = B x T */
		out_samples[i].normal.x = by * tz - bz * ty;
		out_samples[i].normal.y = bz * tx - bx * tz;
		out_samples[i].normal.z = bx * ty - by * tx;
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Winding number                                                     */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_winding_number_2d(
	int *out_winding_number,
	qaws_curve const *curve,
	qaws_vec2 point)
{
	unsigned int sample_count;
	unsigned int i;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_scalar t;
	qaws_scalar angle_sum;
	qaws_scalar prev_dx;
	qaws_scalar prev_dy;
	qaws_scalar curr_dx;
	qaws_scalar curr_dy;
	qaws_scalar cross_val;
	qaws_scalar dot_val;
	qaws_scalar delta_angle;
	qaws_eval_result_2d result;
	qaws_status status;
	double pi2;

	if (!curve || !out_winding_number)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (!qaws_curve_is_closed(curve))
		return QAWS_STATUS_INVALID_ARGUMENT;

	sample_count = 256;
	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Evaluate first sample */
	status = qaws_curve_evaluate_2d(
		&result, curve, range_min, QAWS_EVAL_FLAG_POSITION);
	if (status != QAWS_STATUS_OK)
		return status;

	prev_dx = result.position.x - point.x;
	prev_dy = result.position.y - point.y;
	angle_sum = (qaws_scalar)0.0;

	for (i = 1; i <= sample_count; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)sample_count;

		status = qaws_curve_evaluate_2d(
			&result, curve, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		curr_dx = result.position.x - point.x;
		curr_dy = result.position.y - point.y;

		cross_val = prev_dx * curr_dy - prev_dy * curr_dx;
		dot_val = prev_dx * curr_dx + prev_dy * curr_dy;
		delta_angle = (qaws_scalar)atan2((double)cross_val, (double)dot_val);

		angle_sum += delta_angle;

		prev_dx = curr_dx;
		prev_dy = curr_dy;
	}

	/* Divide by 2*pi and round to nearest integer */
	pi2 = 6.283185307179586476925286766559;
	*out_winding_number = (int)floor((double)angle_sum / pi2 + 0.5);

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Intersection detection via sampling + Newton refinement            */
/* ------------------------------------------------------------------ */

#define ISECT_SAMPLE_COUNT   256
#define ISECT_NEWTON_ITERS   20
#define ISECT_DEDUP_TOLERANCE ((qaws_scalar)1.0e-3)

#if QAWS_SCALAR_IS_FLOAT
#define ISECT_POS_TOLERANCE  ((qaws_scalar)1.0e-3f)
#define ISECT_CONVERGE_TOL   ((qaws_scalar)1.0e-5f)
#define ISECT_CANDIDATE_TOL  ((qaws_scalar)1.0e-3f)
#else
#define ISECT_POS_TOLERANCE  ((qaws_scalar)1.0e-6)
#define ISECT_CONVERGE_TOL   ((qaws_scalar)1.0e-10)
#define ISECT_CANDIDATE_TOL  ((qaws_scalar)1.0e-3)
#endif

/* ------------------------------------------------------------------ */
/*  Surface inspection                                                 */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_compute_bounds(
	qaws_vec3 *out_min,
	qaws_vec3 *out_max,
	qaws_surface const *surface)
{
	unsigned int nu;
	unsigned int nv;
	unsigned int i;
	unsigned int j;
	qaws_scalar u_min;
	qaws_scalar u_max;
	qaws_scalar v_min;
	qaws_scalar v_max;
	qaws_scalar u;
	qaws_scalar v;
	qaws_range u_range;
	qaws_range v_range;
	qaws_surface_eval_result result;
	qaws_status status;

	if (!surface || !out_min || !out_max)
		return QAWS_STATUS_INVALID_ARGUMENT;

	nu = 32;
	nv = 32;

	u_range = qaws_surface_get_u_range(surface);
	v_range = qaws_surface_get_v_range(surface);
	u_min = u_range.min_value;
	u_max = u_range.max_value;
	v_min = v_range.min_value;
	v_max = v_range.max_value;

	/* Evaluate first point to initialize bounds */
	status = qaws_surface_evaluate(
		&result, surface, u_min, v_min, QAWS_SURFACE_EVAL_POSITION);
	if (status != QAWS_STATUS_OK)
		return status;

	out_min->x = result.position.x;
	out_min->y = result.position.y;
	out_min->z = result.position.z;
	out_max->x = result.position.x;
	out_max->y = result.position.y;
	out_max->z = result.position.z;

	for (i = 0; i < nu; ++i)
	{
		u = u_min + (qaws_scalar)i * (u_max - u_min)
			/ (qaws_scalar)(nu - 1);

		for (j = 0; j < nv; ++j)
		{
			/* Skip the (0,0) point already evaluated */
			if (i == 0 && j == 0)
				continue;

			v = v_min + (qaws_scalar)j * (v_max - v_min)
				/ (qaws_scalar)(nv - 1);

			status = qaws_surface_evaluate(
				&result, surface, u, v, QAWS_SURFACE_EVAL_POSITION);
			if (status != QAWS_STATUS_OK)
				return status;

			if (result.position.x < out_min->x)
				out_min->x = result.position.x;
			if (result.position.y < out_min->y)
				out_min->y = result.position.y;
			if (result.position.z < out_min->z)
				out_min->z = result.position.z;
			if (result.position.x > out_max->x)
				out_max->x = result.position.x;
			if (result.position.y > out_max->y)
				out_max->y = result.position.y;
			if (result.position.z > out_max->z)
				out_max->z = result.position.z;
		}
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_compute_area(
	qaws_scalar *out_area,
	qaws_surface const *surface)
{
	unsigned int nu;
	unsigned int nv;
	unsigned int i;
	unsigned int j;
	qaws_scalar u_min;
	qaws_scalar u_max;
	qaws_scalar v_min;
	qaws_scalar v_max;
	qaws_scalar hu;
	qaws_scalar hv;
	qaws_scalar u;
	qaws_scalar v;
	qaws_scalar area;
	qaws_scalar wu;
	qaws_scalar wv;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len;
	qaws_range u_range;
	qaws_range v_range;
	qaws_surface_eval_result result;
	qaws_status status;

	if (!surface || !out_area)
		return QAWS_STATUS_INVALID_ARGUMENT;

	nu = 64;
	nv = 64;

	u_range = qaws_surface_get_u_range(surface);
	v_range = qaws_surface_get_v_range(surface);
	u_min = u_range.min_value;
	u_max = u_range.max_value;
	v_min = v_range.min_value;
	v_max = v_range.max_value;

	hu = (u_max - u_min) / (qaws_scalar)(nu - 1);
	hv = (v_max - v_min) / (qaws_scalar)(nv - 1);

	area = (qaws_scalar)0.0;

	/* Composite Simpson's rule: nu and nv must be odd (64 is even, use 65) */
	/* Actually, for simplicity and to match the stated 64x64 grid,
	   we use nu=nv=65 intervals (65 points, 64 intervals). Simpson's
	   rule requires even number of intervals, so 64 is fine. */
	nu = 65;
	nv = 65;
	hu = (u_max - u_min) / (qaws_scalar)(nu - 1);
	hv = (v_max - v_min) / (qaws_scalar)(nv - 1);

	for (i = 0; i < nu; ++i)
	{
		u = u_min + (qaws_scalar)i * hu;

		/* Simpson weight for u */
		if (i == 0 || i == nu - 1)
			wu = (qaws_scalar)1.0;
		else if (i % 2 == 1)
			wu = (qaws_scalar)4.0;
		else
			wu = (qaws_scalar)2.0;

		for (j = 0; j < nv; ++j)
		{
			v = v_min + (qaws_scalar)j * hv;

			/* Simpson weight for v */
			if (j == 0 || j == nv - 1)
				wv = (qaws_scalar)1.0;
			else if (j % 2 == 1)
				wv = (qaws_scalar)4.0;
			else
				wv = (qaws_scalar)2.0;

			status = qaws_surface_evaluate(
				&result, surface, u,
				v,
				QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV);
			if (status != QAWS_STATUS_OK)
				return status;

			/* Cross product: du x dv */
			cx = result.du.y * result.dv.z - result.du.z * result.dv.y;
			cy = result.du.z * result.dv.x - result.du.x * result.dv.z;
			cz = result.du.x * result.dv.y - result.du.y * result.dv.x;
			cross_len = (qaws_scalar)sqrt(
				(double)(cx * cx + cy * cy + cz * cz));

			area += wu * wv * cross_len;
		}
	}

	area *= hu * hv / (qaws_scalar)9.0;
	*out_area = area;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_compute_gaussian_curvature(
	qaws_scalar *out_curvature,
	qaws_surface const *surface,
	qaws_scalar u,
	qaws_scalar v)
{
	qaws_surface_eval_result result;
	qaws_status status;
	qaws_scalar E, F, G;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len;
	qaws_scalar nx, ny, nz;
	qaws_scalar L, M, N_coeff;
	qaws_scalar denom;

	if (!surface || !out_curvature)
		return QAWS_STATUS_INVALID_ARGUMENT;

	status = qaws_surface_evaluate(
		&result, surface, u,
		v,
		QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
		| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV
		| QAWS_SURFACE_EVAL_DVV);
	if (status != QAWS_STATUS_OK)
		return status;

	/* First fundamental form */
	E = result.du.x * result.du.x + result.du.y * result.du.y
		+ result.du.z * result.du.z;
	F = result.du.x * result.dv.x + result.du.y * result.dv.y
		+ result.du.z * result.dv.z;
	G = result.dv.x * result.dv.x + result.dv.y * result.dv.y
		+ result.dv.z * result.dv.z;

	denom = E * G - F * F;
	if (denom < (qaws_scalar)1.0e-24 && denom > (qaws_scalar)-1.0e-24)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* Unit normal: N = (du x dv) / |du x dv| */
	cx = result.du.y * result.dv.z - result.du.z * result.dv.y;
	cy = result.du.z * result.dv.x - result.du.x * result.dv.z;
	cz = result.du.x * result.dv.y - result.du.y * result.dv.x;
	cross_len = (qaws_scalar)sqrt((double)(cx * cx + cy * cy + cz * cz));

	if (cross_len < (qaws_scalar)1.0e-12)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	nx = cx / cross_len;
	ny = cy / cross_len;
	nz = cz / cross_len;

	/* Second fundamental form */
	L = result.duu.x * nx + result.duu.y * ny + result.duu.z * nz;
	M = result.duv.x * nx + result.duv.y * ny + result.duv.z * nz;
	N_coeff = result.dvv.x * nx + result.dvv.y * ny + result.dvv.z * nz;

	/* Gaussian curvature: K = (L*N - M^2) / (E*G - F^2) */
	*out_curvature = (L * N_coeff - M * M) / denom;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_compute_mean_curvature(
	qaws_scalar *out_curvature,
	qaws_surface const *surface,
	qaws_scalar u,
	qaws_scalar v)
{
	qaws_surface_eval_result result;
	qaws_status status;
	qaws_scalar E, F, G;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len;
	qaws_scalar nx, ny, nz;
	qaws_scalar L, M, N_coeff;
	qaws_scalar denom;

	if (!surface || !out_curvature)
		return QAWS_STATUS_INVALID_ARGUMENT;

	status = qaws_surface_evaluate(
		&result, surface, u,
		v,
		QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
		| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV
		| QAWS_SURFACE_EVAL_DVV);
	if (status != QAWS_STATUS_OK)
		return status;

	/* First fundamental form */
	E = result.du.x * result.du.x + result.du.y * result.du.y
		+ result.du.z * result.du.z;
	F = result.du.x * result.dv.x + result.du.y * result.dv.y
		+ result.du.z * result.dv.z;
	G = result.dv.x * result.dv.x + result.dv.y * result.dv.y
		+ result.dv.z * result.dv.z;

	denom = E * G - F * F;
	if (denom < (qaws_scalar)1.0e-24 && denom > (qaws_scalar)-1.0e-24)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* Unit normal: N = (du x dv) / |du x dv| */
	cx = result.du.y * result.dv.z - result.du.z * result.dv.y;
	cy = result.du.z * result.dv.x - result.du.x * result.dv.z;
	cz = result.du.x * result.dv.y - result.du.y * result.dv.x;
	cross_len = (qaws_scalar)sqrt((double)(cx * cx + cy * cy + cz * cz));

	if (cross_len < (qaws_scalar)1.0e-12)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	nx = cx / cross_len;
	ny = cy / cross_len;
	nz = cz / cross_len;

	/* Second fundamental form */
	L = result.duu.x * nx + result.duu.y * ny + result.duu.z * nz;
	M = result.duv.x * nx + result.duv.y * ny + result.duv.z * nz;
	N_coeff = result.dvv.x * nx + result.dvv.y * ny + result.dvv.z * nz;

	/* Mean curvature: H = (E*N + G*L - 2*F*M) / (2*(E*G - F^2)) */
	*out_curvature = (E * N_coeff + G * L - (qaws_scalar)2.0 * F * M)
		/ ((qaws_scalar)2.0 * denom);
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_compute_principal_curvatures(
	qaws_surface_curvature_result *out_result,
	qaws_surface const *surface,
	qaws_scalar u,
	qaws_scalar v)
{
	qaws_surface_eval_result result;
	qaws_status status;
	qaws_scalar E, F, G;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len;
	qaws_scalar nx, ny, nz;
	qaws_scalar L, M, N_coeff;
	qaws_scalar denom;
	qaws_scalar K, H;
	qaws_scalar disc;

	if (!surface || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	status = qaws_surface_evaluate(
		&result, surface, u,
		v,
		QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
		| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV
		| QAWS_SURFACE_EVAL_DVV);
	if (status != QAWS_STATUS_OK)
		return status;

	/* First fundamental form */
	E = result.du.x * result.du.x + result.du.y * result.du.y
		+ result.du.z * result.du.z;
	F = result.du.x * result.dv.x + result.du.y * result.dv.y
		+ result.du.z * result.dv.z;
	G = result.dv.x * result.dv.x + result.dv.y * result.dv.y
		+ result.dv.z * result.dv.z;

	denom = E * G - F * F;
	if (denom < (qaws_scalar)1.0e-24 && denom > (qaws_scalar)-1.0e-24)
	{
		out_result->gaussian = (qaws_scalar)0.0;
		out_result->mean = (qaws_scalar)0.0;
		out_result->kappa1 = (qaws_scalar)0.0;
		out_result->kappa2 = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* Unit normal */
	cx = result.du.y * result.dv.z - result.du.z * result.dv.y;
	cy = result.du.z * result.dv.x - result.du.x * result.dv.z;
	cz = result.du.x * result.dv.y - result.du.y * result.dv.x;
	cross_len = (qaws_scalar)sqrt((double)(cx * cx + cy * cy + cz * cz));

	if (cross_len < (qaws_scalar)1.0e-12)
	{
		out_result->gaussian = (qaws_scalar)0.0;
		out_result->mean = (qaws_scalar)0.0;
		out_result->kappa1 = (qaws_scalar)0.0;
		out_result->kappa2 = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	nx = cx / cross_len;
	ny = cy / cross_len;
	nz = cz / cross_len;

	/* Second fundamental form */
	L = result.duu.x * nx + result.duu.y * ny + result.duu.z * nz;
	M = result.duv.x * nx + result.duv.y * ny + result.duv.z * nz;
	N_coeff = result.dvv.x * nx + result.dvv.y * ny + result.dvv.z * nz;

	/* Gaussian and mean curvature */
	K = (L * N_coeff - M * M) / denom;
	H = (E * N_coeff + G * L - (qaws_scalar)2.0 * F * M)
		/ ((qaws_scalar)2.0 * denom);

	out_result->gaussian = K;
	out_result->mean = H;

	/* Principal curvatures: kappa = H +/- sqrt(H^2 - K) */
	disc = H * H - K;
	if (disc < (qaws_scalar)0.0)
		disc = (qaws_scalar)0.0;
	disc = (qaws_scalar)sqrt((double)disc);

	out_result->kappa1 = H + disc;
	out_result->kappa2 = H - disc;

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Hausdorff distance                                                 */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_hausdorff_distance_2d(
	qaws_scalar *out_distance,
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	unsigned int sample_count)
{
	unsigned int i;
	qaws_scalar range_a_min;
	qaws_scalar range_a_max;
	qaws_scalar range_b_min;
	qaws_scalar range_b_max;
	qaws_scalar t;
	qaws_scalar param;
	qaws_scalar dx, dy;
	qaws_scalar dist;
	qaws_scalar max_dist_ab;
	qaws_scalar max_dist_ba;
	qaws_eval_result_2d result_a;
	qaws_eval_result_2d result_b;
	qaws_status status;

	if (!curve_a || !curve_b || !out_distance)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve_a->dimension != QAWS_DIMENSION_2D
		|| curve_b->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_a_min = curve_a->parameter_range.min_value;
	range_a_max = curve_a->parameter_range.max_value;
	range_b_min = curve_b->parameter_range.min_value;
	range_b_max = curve_b->parameter_range.max_value;

	/* Directed distance: h(a, b) = max over samples on a of min distance to b */
	max_dist_ab = (qaws_scalar)0.0;

	for (i = 0; i < sample_count; ++i)
	{
		qaws_vec2 pt;

		t = range_a_min + (qaws_scalar)i * (range_a_max - range_a_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_2d(
			&result_a, curve_a, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_a.position;

		/* Find closest point on curve_b */
		status = qaws_curve_find_closest_parameter_2d(
			&param, curve_b, pt);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_2d(
			&result_b, curve_b, param, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result_b.position.x - pt.x;
		dy = result_b.position.y - pt.y;
		dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));

		if (dist > max_dist_ab)
			max_dist_ab = dist;
	}

	/* Directed distance: h(b, a) = max over samples on b of min distance to a */
	max_dist_ba = (qaws_scalar)0.0;

	for (i = 0; i < sample_count; ++i)
	{
		qaws_vec2 pt;

		t = range_b_min + (qaws_scalar)i * (range_b_max - range_b_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_2d(
			&result_b, curve_b, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_b.position;

		/* Find closest point on curve_a */
		status = qaws_curve_find_closest_parameter_2d(
			&param, curve_a, pt);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_2d(
			&result_a, curve_a, param, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result_a.position.x - pt.x;
		dy = result_a.position.y - pt.y;
		dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));

		if (dist > max_dist_ba)
			max_dist_ba = dist;
	}

	/* Hausdorff distance = max(h(a,b), h(b,a)) */
	*out_distance = max_dist_ab > max_dist_ba ? max_dist_ab : max_dist_ba;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_hausdorff_distance_3d(
	qaws_scalar *out_distance,
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	unsigned int sample_count)
{
	unsigned int i;
	qaws_scalar range_a_min;
	qaws_scalar range_a_max;
	qaws_scalar range_b_min;
	qaws_scalar range_b_max;
	qaws_scalar t;
	qaws_scalar param;
	qaws_scalar dx, dy, dz;
	qaws_scalar dist;
	qaws_scalar max_dist_ab;
	qaws_scalar max_dist_ba;
	qaws_eval_result_3d result_a;
	qaws_eval_result_3d result_b;
	qaws_status status;

	if (!curve_a || !curve_b || !out_distance)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve_a->dimension != QAWS_DIMENSION_3D
		|| curve_b->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_a_min = curve_a->parameter_range.min_value;
	range_a_max = curve_a->parameter_range.max_value;
	range_b_min = curve_b->parameter_range.min_value;
	range_b_max = curve_b->parameter_range.max_value;

	/* Directed distance: h(a, b) */
	max_dist_ab = (qaws_scalar)0.0;

	for (i = 0; i < sample_count; ++i)
	{
		qaws_vec3 pt;

		t = range_a_min + (qaws_scalar)i * (range_a_max - range_a_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_3d(
			&result_a, curve_a, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_a.position;

		/* Find closest point on curve_b */
		status = qaws_curve_find_closest_parameter_3d(
			&param, curve_b, pt);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_3d(
			&result_b, curve_b, param, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result_b.position.x - pt.x;
		dy = result_b.position.y - pt.y;
		dz = result_b.position.z - pt.z;
		dist = (qaws_scalar)sqrt(
			(double)(dx * dx + dy * dy + dz * dz));

		if (dist > max_dist_ab)
			max_dist_ab = dist;
	}

	/* Directed distance: h(b, a) */
	max_dist_ba = (qaws_scalar)0.0;

	for (i = 0; i < sample_count; ++i)
	{
		qaws_vec3 pt;

		t = range_b_min + (qaws_scalar)i * (range_b_max - range_b_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_3d(
			&result_b, curve_b, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_b.position;

		/* Find closest point on curve_a */
		status = qaws_curve_find_closest_parameter_3d(
			&param, curve_a, pt);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_3d(
			&result_a, curve_a, param, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result_a.position.x - pt.x;
		dy = result_a.position.y - pt.y;
		dz = result_a.position.z - pt.z;
		dist = (qaws_scalar)sqrt(
			(double)(dx * dx + dy * dy + dz * dz));

		if (dist > max_dist_ba)
			max_dist_ba = dist;
	}

	/* Hausdorff distance = max(h(a,b), h(b,a)) */
	*out_distance = max_dist_ab > max_dist_ba ? max_dist_ab : max_dist_ba;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Closest point on surface                                           */
/* ------------------------------------------------------------------ */

/* Newton from (u, v) on the squared distance to point, clamped to the domain; the
   result's squared distance */
static qaws_status surface_closest_refine(
	qaws_surface const *surface,
	qaws_vec3 point,
	qaws_scalar *io_u,
	qaws_scalar *io_v,
	qaws_vec3 *out_point,
	qaws_scalar *out_dist_sq)
{
	unsigned int iter;
	qaws_scalar u = *io_u, v = *io_v;
	qaws_scalar g_u, g_v, h_uu, h_uv, h_vv, det, delta_u, delta_v, delta_len;
	qaws_range u_range = qaws_surface_get_u_range(surface);
	qaws_range v_range = qaws_surface_get_v_range(surface);
	qaws_surface_eval_result result;
	qaws_status status;
	qaws_vec3 r;

	for (iter = 0; iter < 20; ++iter)
	{
		status = qaws_surface_evaluate(
			&result, surface, u,
			v,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU
			| QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DUU
			| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV);
		if (status != QAWS_STATUS_OK)
			return status;

		r.x = result.position.x - point.x;
		r.y = result.position.y - point.y;
		r.z = result.position.z - point.z;

		/* Gradient */
		g_u = r.x * result.du.x + r.y * result.du.y + r.z * result.du.z;
		g_v = r.x * result.dv.x + r.y * result.dv.y + r.z * result.dv.z;

		/* Hessian approximation */
		h_uu = result.du.x * result.du.x + result.du.y * result.du.y
			+ result.du.z * result.du.z
			+ r.x * result.duu.x + r.y * result.duu.y
			+ r.z * result.duu.z;
		h_uv = result.du.x * result.dv.x + result.du.y * result.dv.y
			+ result.du.z * result.dv.z
			+ r.x * result.duv.x + r.y * result.duv.y
			+ r.z * result.duv.z;
		h_vv = result.dv.x * result.dv.x + result.dv.y * result.dv.y
			+ result.dv.z * result.dv.z
			+ r.x * result.dvv.x + r.y * result.dvv.y
			+ r.z * result.dvv.z;

		/* Solve 2x2 system: [h_uu h_uv; h_uv h_vv] * [du; dv] = [g_u; g_v] */
		det = h_uu * h_vv - h_uv * h_uv;
		if (QAWS_FABS(det) < QAWS_LITERAL(1.0e-30))
			break;

		delta_u = (h_vv * g_u - h_uv * g_v) / det;
		delta_v = (h_uu * g_v - h_uv * g_u) / det;

		u -= delta_u;
		v -= delta_v;

		/* Clamp to parameter domain */
		if (u < u_range.min_value) u = u_range.min_value;
		if (u > u_range.max_value) u = u_range.max_value;
		if (v < v_range.min_value) v = v_range.min_value;
		if (v > v_range.max_value) v = v_range.max_value;

		delta_len = QAWS_SQRT(delta_u * delta_u + delta_v * delta_v);
		if (delta_len < QAWS_LITERAL(1.0e-8))
			break;
	}

	/* Evaluate final position */
	status = qaws_surface_evaluate(
		&result, surface, u, v, QAWS_SURFACE_EVAL_POSITION);
	if (status != QAWS_STATUS_OK)
		return status;
	r.x = result.position.x - point.x;
	r.y = result.position.y - point.y;
	r.z = result.position.z - point.z;
	*io_u = u;
	*io_v = v;
	*out_point = result.position;
	*out_dist_sq = r.x * r.x + r.y * r.y + r.z * r.z;
	return QAWS_STATUS_OK;
}

/* an 8 x 8 grid's nearest samples, each refined by Newton, the nearest result kept:
   several seeds, as a closed surface's seam samples (u at both ends, one point) tie
   and the one taken by rounding can stall against the domain's end */
static qaws_status surface_closest_sampled(
	qaws_surface const *surface,
	qaws_vec3 point,
	qaws_scalar *out_u,
	qaws_scalar *out_v,
	qaws_vec3 *out_closest_point)
{
	unsigned int i, j, k, ns = 0;
	qaws_scalar u, v, dx, dy, dz, dist_sq, best = QAWS_LITERAL(1.0e30);
	qaws_scalar su[4], sv[4], sd[4];
	qaws_scalar u_min, u_max, v_min, v_max;
	qaws_range u_range;
	qaws_range v_range;
	qaws_surface_eval_result result;
	qaws_status status;

	if (!surface || !out_u || !out_v || !out_closest_point)
		return QAWS_STATUS_INVALID_ARGUMENT;

	u_range = qaws_surface_get_u_range(surface);
	v_range = qaws_surface_get_v_range(surface);
	u_min = u_range.min_value;
	u_max = u_range.max_value;
	v_min = v_range.min_value;
	v_max = v_range.max_value;

	/* Coarse search: 8x8 grid, its 4 nearest samples in order */
	for (i = 0; i < 8; ++i)
	{
		u = u_min + (qaws_scalar)i * (u_max - u_min)
			/ QAWS_LITERAL(7.0);

		for (j = 0; j < 8; ++j)
		{
			v = v_min + (qaws_scalar)j * (v_max - v_min)
				/ QAWS_LITERAL(7.0);

			status = qaws_surface_evaluate(
				&result, surface, u, v, QAWS_SURFACE_EVAL_POSITION);
			if (status != QAWS_STATUS_OK)
				return status;

			dx = result.position.x - point.x;
			dy = result.position.y - point.y;
			dz = result.position.z - point.z;
			dist_sq = dx * dx + dy * dy + dz * dz;

			if (ns < 4 || dist_sq < sd[ns - 1])
			{
				k = ns < 4 ? ns++ : ns - 1;
				for (; k > 0 && sd[k - 1] > dist_sq; --k)
				{
					su[k] = su[k - 1];
					sv[k] = sv[k - 1];
					sd[k] = sd[k - 1];
				}
				su[k] = u;
				sv[k] = v;
				sd[k] = dist_sq;
			}
		}
	}

	/* Newton refinement from each, the nearest result kept (ties: the earlier seed) */
	for (k = 0; k < ns; ++k)
	{
		qaws_vec3 p;
		u = su[k];
		v = sv[k];
		status = surface_closest_refine(surface, point, &u, &v, &p, &dist_sq);
		if (status != QAWS_STATUS_OK)
			return status;
		if (dist_sq < best)
		{
			best = dist_sq;
			*out_u = u;
			*out_v = v;
			*out_closest_point = p;
		}
	}
	return QAWS_STATUS_OK;
}

/* the sampled answer, then certified where the surface has Bezier patches */
static qaws_status surface_closest(
	qaws_surface const *surface,
	qaws_vec3 point,
	qaws_scalar *out_u,
	qaws_scalar *out_v,
	qaws_vec3 *out_closest_point,
	double *out_distance,
	double *out_gap)
{
	double x[3] = { point.x, point.y, point.z };
	qaws_surface_eval_result result;
	qaws_status status = surface_closest_sampled(surface, point, out_u, out_v, out_closest_point);
	if (status != QAWS_STATUS_OK)
		return status;
	status = qaws_internal_surface_closest_certify(out_distance, out_gap, out_u, out_v, surface, x);
	if (status != QAWS_STATUS_OK)
		return status;
	status = qaws_surface_evaluate(&result, surface, *out_u, *out_v, QAWS_SURFACE_EVAL_POSITION);
	if (status == QAWS_STATUS_OK)
		*out_closest_point = result.position;
	return status;
}

qaws_status qaws_surface_find_closest_point(
	qaws_scalar *out_u,
	qaws_scalar *out_v,
	qaws_vec3 *out_closest_point,
	qaws_surface const *surface,
	qaws_vec3 point)
{
	double distance, gap;
	qaws_status status = surface_closest(surface, point, out_u, out_v, out_closest_point, &distance, &gap);
	return status == QAWS_STATUS_UNSUPPORTED_OPERATION ? QAWS_STATUS_OK : status;
}

qaws_status qaws_surface_find_closest_point_certified(
	qaws_scalar *out_u,
	qaws_scalar *out_v,
	qaws_vec3 *out_closest_point,
	qaws_scalar *out_distance,
	qaws_scalar *out_gap,
	qaws_surface const *surface,
	qaws_vec3 point)
{
	double distance, gap;
	qaws_status status;
	if (!out_distance || !out_gap)
		return QAWS_STATUS_INVALID_ARGUMENT;
	status = surface_closest(surface, point, out_u, out_v, out_closest_point, &distance, &gap);
	if (status != QAWS_STATUS_OK)
		return status;
	*out_distance = (qaws_scalar)distance;
	*out_gap = (qaws_scalar)gap;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Adaptive tessellation                                              */
/* ------------------------------------------------------------------ */

/* Internal quad descriptor for the tessellation stack */
typedef struct qaws_tess_quad {
	qaws_scalar u0, v0;
	qaws_scalar u1, v1;
	unsigned int depth;
} qaws_tess_quad;

/* Leaf quad collected during subdivision */
typedef struct qaws_tess_leaf {
	qaws_scalar u0, v0;
	qaws_scalar u1, v1;
} qaws_tess_leaf;

/* Find existing vertex at (u,v).  Returns index or (unsigned)-1. */
static unsigned int qaws_tess_find_vertex(
	qaws_tessellation_vertex const *vertices,
	unsigned int vertex_count,
	qaws_scalar u,
	qaws_scalar v)
{
	unsigned int i;
	for (i = 0; i < vertex_count; ++i)
	{
		if (QAWS_FABS(vertices[i].u - u) < QAWS_LITERAL(1.0e-10)
			&& QAWS_FABS(vertices[i].v - v) < QAWS_LITERAL(1.0e-10))
			return i;
	}
	return (unsigned int)-1;
}

/* Evaluate surface and add vertex, with deduplication by (u,v). */
static unsigned int qaws_tess_eval_add(
	qaws_surface const *surface,
	qaws_tessellation_vertex *vertices,
	unsigned int *vertex_count,
	unsigned int vertex_capacity,
	qaws_scalar u,
	qaws_scalar v,
	qaws_status *out_status)
{
	unsigned int idx;
	unsigned int n;
	qaws_surface_eval_result r;

	idx = qaws_tess_find_vertex(vertices, *vertex_count, u, v);
	if (idx != (unsigned int)-1)
		return idx;

	memset(&r, 0, sizeof(r));
	*out_status = qaws_surface_evaluate(&r, surface, u,
		v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
	if (*out_status != QAWS_STATUS_OK)
		return 0;

	n = *vertex_count;
	if (n < vertex_capacity)
	{
		vertices[n].position = r.position;
		vertices[n].normal = r.normal;
		vertices[n].u = u;
		vertices[n].v = v;
		*vertex_count = n + 1;
		return n;
	}
	return n > 0 ? n - 1 : 0;
}

/* Collect all vertex indices on a horizontal edge (v = v_const,
   u in [u_lo, u_hi]).  Output is sorted by ascending u. */
static unsigned int qaws_tess_collect_h(
	qaws_tessellation_vertex const *verts,
	unsigned int vert_count,
	qaws_scalar v_const,
	qaws_scalar u_lo,
	qaws_scalar u_hi,
	unsigned int *out,
	unsigned int max_out)
{
	unsigned int n = 0, i, j;
	unsigned int tmp;

	for (i = 0; i < vert_count && n < max_out; ++i)
	{
		if (QAWS_FABS(verts[i].v - v_const) < QAWS_LITERAL(1.0e-10)
			&& verts[i].u >= u_lo - QAWS_LITERAL(1.0e-10)
			&& verts[i].u <= u_hi + QAWS_LITERAL(1.0e-10))
		{
			out[n++] = i;
		}
	}
	/* Insertion sort by u */
	for (i = 1; i < n; ++i)
		for (j = i; j > 0 && verts[out[j]].u < verts[out[j - 1]].u; --j)
		{
			tmp = out[j]; out[j] = out[j - 1]; out[j - 1] = tmp;
		}
	return n;
}

/* Collect all vertex indices on a vertical edge (u = u_const,
   v in [v_lo, v_hi]).  Output is sorted by ascending v. */
static unsigned int qaws_tess_collect_v(
	qaws_tessellation_vertex const *verts,
	unsigned int vert_count,
	qaws_scalar u_const,
	qaws_scalar v_lo,
	qaws_scalar v_hi,
	unsigned int *out,
	unsigned int max_out)
{
	unsigned int n = 0, i, j;
	unsigned int tmp;

	for (i = 0; i < vert_count && n < max_out; ++i)
	{
		if (QAWS_FABS(verts[i].u - u_const) < QAWS_LITERAL(1.0e-10)
			&& verts[i].v >= v_lo - QAWS_LITERAL(1.0e-10)
			&& verts[i].v <= v_hi + QAWS_LITERAL(1.0e-10))
		{
			out[n++] = i;
		}
	}
	/* Insertion sort by v */
	for (i = 1; i < n; ++i)
		for (j = i; j > 0 && verts[out[j]].v < verts[out[j - 1]].v; --j)
		{
			tmp = out[j]; out[j] = out[j - 1]; out[j - 1] = tmp;
		}
	return n;
}

qaws_status qaws_surface_tessellate(
	qaws_tessellation_vertex *out_vertices,
	unsigned int *out_vertex_count,
	unsigned int *out_indices,
	unsigned int *out_index_count,
	qaws_surface const *surface,
	unsigned int vertex_capacity,
	unsigned int index_capacity,
	qaws_tessellation_desc const *desc)
{
	unsigned int max_depth;
	qaws_scalar curv_thresh;
	qaws_scalar max_edge;
	unsigned int vert_count;
	unsigned int idx_count;
	unsigned int stack_top;
	unsigned int li;
	qaws_scalar u_min, u_max, v_min, v_max;
	qaws_range u_range;
	qaws_range v_range;
	qaws_status status;

	/* Work structures */
	unsigned int stack_capacity;
	qaws_tess_quad *stack;
	unsigned int leaf_capacity;
	qaws_tess_leaf *leaves;
	unsigned int leaf_count;

	if (!surface || !desc || !out_vertices || !out_vertex_count
		|| !out_indices || !out_index_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Apply defaults */
	max_depth = desc->max_depth > 0 ? desc->max_depth : 5;
	curv_thresh = desc->curvature_threshold > QAWS_LITERAL(0.0)
		? desc->curvature_threshold : QAWS_LITERAL(0.1);
	max_edge = desc->max_edge_length;

	u_range = qaws_surface_get_u_range(surface);
	v_range = qaws_surface_get_v_range(surface);
	u_min = u_range.min_value;
	u_max = u_range.max_value;
	v_min = v_range.min_value;
	v_max = v_range.max_value;

	vert_count = 0;
	idx_count = 0;

	/* Allocate work structures */
	stack_capacity = 4096;
	leaf_capacity = 4096;
	stack = (qaws_tess_quad *)qaws_internal_scratch_malloc(
		qaws_internal_mul2(stack_capacity, sizeof(qaws_tess_quad)));
	leaves = (qaws_tess_leaf *)qaws_internal_scratch_malloc(
		qaws_internal_mul2(leaf_capacity, sizeof(qaws_tess_leaf)));
	if (!stack || !leaves)
	{
		qaws_internal_free(stack);
		qaws_internal_free(leaves);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}

	/* ---- Phase 1: Quadtree subdivision, collect leaf quads ---- */
	leaf_count = 0;
	stack[0].u0 = u_min;
	stack[0].v0 = v_min;
	stack[0].u1 = u_max;
	stack[0].v1 = v_max;
	stack[0].depth = 0;
	stack_top = 1;

	while (stack_top > 0)
	{
		qaws_tess_quad quad;
		qaws_scalar u_mid, v_mid;
		qaws_surface_eval_result r00, r10, r01, r11, r_center;
		qaws_vec3 bilinear_mid;
		qaws_scalar flat_dx, flat_dy, flat_dz;
		qaws_scalar flatness;
		int subdivide;

		--stack_top;
		quad = stack[stack_top];

		u_mid = QAWS_LITERAL(0.5) * (quad.u0 + quad.u1);
		v_mid = QAWS_LITERAL(0.5) * (quad.v0 + quad.v1);

		/* Evaluate 4 corners + center */
		status = qaws_surface_evaluate(&r00, surface, quad.u0,
			quad.v0,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(stack); qaws_internal_free(leaves); return status;
		}

		status = qaws_surface_evaluate(&r10, surface, quad.u1,
			quad.v0,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(stack); qaws_internal_free(leaves); return status;
		}

		status = qaws_surface_evaluate(&r01, surface, quad.u0,
			quad.v1,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(stack); qaws_internal_free(leaves); return status;
		}

		status = qaws_surface_evaluate(&r11, surface, quad.u1,
			quad.v1,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(stack); qaws_internal_free(leaves); return status;
		}

		status = qaws_surface_evaluate(&r_center, surface, u_mid,
			v_mid,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
		if (status != QAWS_STATUS_OK)
		{
			qaws_internal_free(stack); qaws_internal_free(leaves); return status;
		}

		/* Check subdivision criteria */
		subdivide = 0;

		if (quad.depth < max_depth)
		{
			/* Flatness: compare center to bilinear interpolation */
			bilinear_mid.x = QAWS_LITERAL(0.25)
				* (r00.position.x + r10.position.x
					+ r01.position.x + r11.position.x);
			bilinear_mid.y = QAWS_LITERAL(0.25)
				* (r00.position.y + r10.position.y
					+ r01.position.y + r11.position.y);
			bilinear_mid.z = QAWS_LITERAL(0.25)
				* (r00.position.z + r10.position.z
					+ r01.position.z + r11.position.z);

			flat_dx = r_center.position.x - bilinear_mid.x;
			flat_dy = r_center.position.y - bilinear_mid.y;
			flat_dz = r_center.position.z - bilinear_mid.z;
			flatness = QAWS_SQRT(
				flat_dx * flat_dx + flat_dy * flat_dy
				+ flat_dz * flat_dz);

			if (flatness > curv_thresh)
				subdivide = 1;

			/* Edge length check */
			if (max_edge > QAWS_LITERAL(0.0) && !subdivide)
			{
				qaws_scalar ex, ey, ez, elen;

				ex = r10.position.x - r00.position.x;
				ey = r10.position.y - r00.position.y;
				ez = r10.position.z - r00.position.z;
				elen = QAWS_SQRT(ex * ex + ey * ey + ez * ez);
				if (elen > max_edge) subdivide = 1;

				ex = r11.position.x - r10.position.x;
				ey = r11.position.y - r10.position.y;
				ez = r11.position.z - r10.position.z;
				elen = QAWS_SQRT(ex * ex + ey * ey + ez * ez);
				if (elen > max_edge) subdivide = 1;

				ex = r11.position.x - r01.position.x;
				ey = r11.position.y - r01.position.y;
				ez = r11.position.z - r01.position.z;
				elen = QAWS_SQRT(ex * ex + ey * ey + ez * ez);
				if (elen > max_edge) subdivide = 1;

				ex = r01.position.x - r00.position.x;
				ey = r01.position.y - r00.position.y;
				ez = r01.position.z - r00.position.z;
				elen = QAWS_SQRT(ex * ex + ey * ey + ez * ez);
				if (elen > max_edge) subdivide = 1;
			}
		}

		if (subdivide && stack_top + 4 > stack_capacity)
		{
			unsigned int new_cap = stack_capacity * 2;
			qaws_tess_quad* new_stack = (qaws_tess_quad*)qaws_internal_realloc(
				stack, qaws_internal_mul2(new_cap, sizeof(qaws_tess_quad)));
			if (!new_stack)
			{
				qaws_internal_free(stack); qaws_internal_free(leaves);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			stack = new_stack;
			stack_capacity = new_cap;
		}
		if (subdivide)
		{
			/* Split into 4 sub-quads */
			stack[stack_top].u0 = quad.u0;
			stack[stack_top].v0 = quad.v0;
			stack[stack_top].u1 = u_mid;
			stack[stack_top].v1 = v_mid;
			stack[stack_top].depth = quad.depth + 1;
			++stack_top;

			stack[stack_top].u0 = u_mid;
			stack[stack_top].v0 = quad.v0;
			stack[stack_top].u1 = quad.u1;
			stack[stack_top].v1 = v_mid;
			stack[stack_top].depth = quad.depth + 1;
			++stack_top;

			stack[stack_top].u0 = quad.u0;
			stack[stack_top].v0 = v_mid;
			stack[stack_top].u1 = u_mid;
			stack[stack_top].v1 = quad.v1;
			stack[stack_top].depth = quad.depth + 1;
			++stack_top;

			stack[stack_top].u0 = u_mid;
			stack[stack_top].v0 = v_mid;
			stack[stack_top].u1 = quad.u1;
			stack[stack_top].v1 = quad.v1;
			stack[stack_top].depth = quad.depth + 1;
			++stack_top;
		}
		else
		{
			/* Store as leaf quad, grow buffer if needed */
			if (leaf_count >= leaf_capacity)
			{
				unsigned int new_cap = leaf_capacity * 2;
				qaws_tess_leaf* new_leaves = (qaws_tess_leaf*)qaws_internal_realloc(
					leaves, qaws_internal_mul2(new_cap, sizeof(qaws_tess_leaf)));
				if (!new_leaves)
				{
					qaws_internal_free(stack); qaws_internal_free(leaves);
					return QAWS_STATUS_ALLOCATION_FAILURE;
				}
				leaves = new_leaves;
				leaf_capacity = new_cap;
			}
			leaves[leaf_count].u0 = quad.u0;
			leaves[leaf_count].v0 = quad.v0;
			leaves[leaf_count].u1 = quad.u1;
			leaves[leaf_count].v1 = quad.v1;
			++leaf_count;
		}
	}

	qaws_internal_free(stack);

	/* ---- Phase 2: Add all leaf corner vertices ---- */
	status = QAWS_STATUS_OK;
	for (li = 0; li < leaf_count; ++li)
	{
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u0, leaves[li].v0, &status);
		if (status != QAWS_STATUS_OK) { qaws_internal_free(leaves); return status; }
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u1, leaves[li].v0, &status);
		if (status != QAWS_STATUS_OK) { qaws_internal_free(leaves); return status; }
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u0, leaves[li].v1, &status);
		if (status != QAWS_STATUS_OK) { qaws_internal_free(leaves); return status; }
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u1, leaves[li].v1, &status);
		if (status != QAWS_STATUS_OK) { qaws_internal_free(leaves); return status; }
	}

	/* ---- Phase 3: Triangulate each leaf with T-junction stitching ----
	   For each leaf, walk its 4 edges and collect ALL vertices that
	   lie on each edge (including midpoints from finer neighbours).
	   Fan-triangulate from the quad center to the perimeter ring.
	   This eliminates T-junctions: every edge vertex referenced by
	   a finer neighbour is also in the coarser quad's triangulation. */
	for (li = 0; li < leaf_count; ++li)
	{
		qaws_scalar u0, v0, u1, v1, u_mid, v_mid;
		unsigned int perim[128];
		unsigned int perim_count;
		unsigned int edge_buf[64];
		unsigned int edge_n;
		unsigned int i_center;
		unsigned int k;
		unsigned int tmp;

		u0 = leaves[li].u0;
		v0 = leaves[li].v0;
		u1 = leaves[li].u1;
		v1 = leaves[li].v1;
		u_mid = QAWS_LITERAL(0.5) * (u0 + u1);
		v_mid = QAWS_LITERAL(0.5) * (v0 + v1);
		perim_count = 0;

		/* Bottom edge: u ascending at v=v0 */
		edge_n = qaws_tess_collect_h(out_vertices, vert_count,
			v0, u0, u1, edge_buf, 64);
		for (k = 0; k < edge_n; ++k)
			perim[perim_count++] = edge_buf[k];

		/* Right edge: v ascending at u=u1, skip first (=last of bottom) */
		edge_n = qaws_tess_collect_v(out_vertices, vert_count,
			u1, v0, v1, edge_buf, 64);
		for (k = 1; k < edge_n; ++k)
			perim[perim_count++] = edge_buf[k];

		/* Top edge: u DESCENDING at v=v1, skip first (=last of right) */
		edge_n = qaws_tess_collect_h(out_vertices, vert_count,
			v1, u0, u1, edge_buf, 64);
		/* Reverse */
		for (k = 0; k < edge_n / 2; ++k)
		{
			tmp = edge_buf[k];
			edge_buf[k] = edge_buf[edge_n - 1 - k];
			edge_buf[edge_n - 1 - k] = tmp;
		}
		for (k = 1; k < edge_n; ++k)
			perim[perim_count++] = edge_buf[k];

		/* Left edge: v DESCENDING at u=u0, skip first and last */
		edge_n = qaws_tess_collect_v(out_vertices, vert_count,
			u0, v0, v1, edge_buf, 64);
		/* Reverse */
		for (k = 0; k < edge_n / 2; ++k)
		{
			tmp = edge_buf[k];
			edge_buf[k] = edge_buf[edge_n - 1 - k];
			edge_buf[edge_n - 1 - k] = tmp;
		}
		for (k = 1; k + 1 < edge_n; ++k)
			perim[perim_count++] = edge_buf[k];

		if (perim_count < 3)
			continue;

		/* Add center vertex for fan */
		i_center = qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, u_mid, v_mid, &status);
		if (status != QAWS_STATUS_OK) { qaws_internal_free(leaves); return status; }

		/* Emit fan triangles: (center, perim[k], perim[k+1]) */
		for (k = 0; k < perim_count; ++k)
		{
			unsigned int next = (k + 1) % perim_count;

			if (idx_count + 3 > index_capacity)
				break;
			out_indices[idx_count++] = i_center;
			out_indices[idx_count++] = perim[k];
			out_indices[idx_count++] = perim[next];
		}
	}

	qaws_internal_free(leaves);

	*out_vertex_count = vert_count;
	*out_index_count = idx_count;
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Curve projection onto surface (#20)                                       */
/* ========================================================================== */

qaws_status qaws_surface_project_curve(
	qaws_curve** out_uv_curve,
	qaws_surface const* surface,
	qaws_curve const* curve,
	unsigned int sample_count)
{
	unsigned int n_samples;
	qaws_scalar* params = NULL;
	qaws_scalar* uv_coords = NULL;
	qaws_range range;
	unsigned int i;
	unsigned int n_cp;
	qaws_status status;

	if (!surface || !curve || !out_uv_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(curve) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_uv_curve = NULL;
	n_samples = sample_count > 0 ? sample_count : 128;
	range = qaws_curve_get_parameter_range(curve);

	params = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)n_samples));
	uv_coords = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), (size_t)n_samples, 2));
	if (!params || !uv_coords)
	{
		qaws_internal_free(params); qaws_internal_free(uv_coords);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 0; i < n_samples; i++)
	{
		qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
			(qaws_scalar)i / (qaws_scalar)(n_samples - 1);
		qaws_eval_result_3d er;
		qaws_scalar u, v;
		qaws_vec3 closest;

		memset(&er, 0, sizeof(er));
		status = qaws_curve_evaluate_3d(&er, curve, t, QAWS_EVAL_FLAG_POSITION);
		if (status != QAWS_STATUS_OK) { qaws_internal_free(params); qaws_internal_free(uv_coords); return status; }

		status = qaws_surface_find_closest_point(&u, &v, &closest, surface, er.position);
		if (status != QAWS_STATUS_OK) { qaws_internal_free(params); qaws_internal_free(uv_coords); return status; }

		params[i] = t;
		uv_coords[i * 2 + 0] = u;
		uv_coords[i * 2 + 1] = v;
	}

	/* Fit a 2D B-spline through the (u,v) samples */
	n_cp = n_samples / 4;
	if (n_cp < 6) n_cp = 6;
	if (n_cp > n_samples) n_cp = n_samples;

	status = qaws_internal_fit_bspline_ex(out_uv_curve, 
		params, uv_coords,
		QAWS_DIMENSION_2D, n_samples,
		n_cp, 3,
		curve->allocator);

	qaws_internal_free(params);
	qaws_internal_free(uv_coords);
	return status;
}
