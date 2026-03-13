#include "qaws_inspect.h"
#include "qaws_eval.h"
#include "qaws_surface.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
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
	qaws_curve const *curve,
	qaws_scalar parameter_min,
	qaws_scalar parameter_max,
	qaws_scalar *out_length)
{
	if (!curve || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;

	return qaws_internal_arc_length(
		curve, parameter_min, parameter_max, out_length);
}

qaws_status qaws_curve_compute_bounds_2d(
	qaws_curve const *curve,
	qaws_vec2 *out_min,
	qaws_vec2 *out_max)
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
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
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
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
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
	qaws_curve const *curve,
	qaws_vec3 *out_min,
	qaws_vec3 *out_max)
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
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
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
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
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

	qaws_curve_evaluate_2d(c->curve, t,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);

	out_pos[0] = result.position.x;
	out_pos[1] = result.position.y;
	out_d1[0] = result.d1.x;
	out_d1[1] = result.d1.y;
}

qaws_status qaws_curve_find_closest_parameter_2d(
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
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
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
		closest_eval_2d, &ctx, target_pos, 2,
		range_min, range_max, best_t, 32);

	*out_parameter = best_t;
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

	qaws_curve_evaluate_3d(c->curve, t,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);

	out_pos[0] = result.position.x;
	out_pos[1] = result.position.y;
	out_pos[2] = result.position.z;
	out_d1[0] = result.d1.x;
	out_d1[1] = result.d1.y;
	out_d1[2] = result.d1.z;
}

qaws_status qaws_curve_find_closest_parameter_3d(
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
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
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
		closest_eval_3d, &ctx, target_pos, 3,
		range_min, range_max, best_t, 32);

	*out_parameter = best_t;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Span continuity                                                    */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_get_span_continuity(
	qaws_curve const *curve,
	unsigned int boundary_index,
	qaws_continuity *out_continuity)
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
	qaws_curve const *curve,
	void *out_control_points,
	unsigned int point_capacity,
	unsigned int *out_point_count)
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

/* ------------------------------------------------------------------ */
/*  Family-specific inspection: B-Spline knots                         */
/* ------------------------------------------------------------------ */

qaws_status qaws_bspline_get_knots(
	qaws_curve const *curve,
	qaws_scalar *out_knots,
	unsigned int knot_capacity,
	unsigned int *out_knot_count)
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
	qaws_curve const *curve,
	qaws_scalar *out_weights,
	unsigned int weight_capacity,
	unsigned int *out_weight_count)
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec2 *out_tangent)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		curve, parameter, QAWS_EVAL_FLAG_D1, &result);
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec3 *out_tangent)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		curve, parameter, QAWS_EVAL_FLAG_D1, &result);
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_curvature)
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
		curve, parameter, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_curvature)
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
		curve, parameter, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_torsion)
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
		curve, parameter,
		QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3,
		&result);
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_speed)
{
	qaws_status status;

	if (!curve || !out_speed)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d result;
		status = qaws_curve_evaluate_2d(
			curve, parameter, QAWS_EVAL_FLAG_D1, &result);
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
			curve, parameter, QAWS_EVAL_FLAG_D1, &result);
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec2 *out_normal)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_normal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		curve, parameter, QAWS_EVAL_FLAG_D1, &result);
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
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec3 *out_tangent,
	qaws_vec3 *out_normal,
	qaws_vec3 *out_binormal)
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
		curve, parameter,
		QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
	qaws_curve const *curve,
	qaws_scalar *out_parameters,
	unsigned int parameter_capacity,
	unsigned int *out_count)
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
			curve, range_min,
			QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		prev_val = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
		prev_t = range_min;

		for (i = 1; i <= total_samples; ++i)
		{
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)total_samples;

			status = qaws_curve_evaluate_2d(
				curve, t,
				QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
						curve, mid,
						QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
			curve, range_min,
			QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
				curve, t,
				QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
						curve, mid,
						QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
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
	qaws_curve const *curve,
	unsigned int axis,
	qaws_scalar *out_parameters,
	unsigned int parameter_capacity,
	unsigned int *out_count)
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
			curve, range_min, QAWS_EVAL_FLAG_D1, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		prev_val = (axis == 0) ? result.d1.x : result.d1.y;
		prev_t = range_min;

		for (i = 1; i <= total_samples; ++i)
		{
			t = range_min + (qaws_scalar)i * (range_max - range_min)
				/ (qaws_scalar)total_samples;

			status = qaws_curve_evaluate_2d(
				curve, t, QAWS_EVAL_FLAG_D1, &result);
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
						curve, mid, QAWS_EVAL_FLAG_D1, &result);
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
			curve, range_min, QAWS_EVAL_FLAG_D1, &result);
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
				curve, t, QAWS_EVAL_FLAG_D1, &result);
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
						curve, mid, QAWS_EVAL_FLAG_D1, &result);
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
	qaws_curve const *curve,
	unsigned int sample_count,
	qaws_curvature_sample_2d *out_samples,
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
			curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
			&result);
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
	qaws_curve const *curve,
	unsigned int sample_count,
	qaws_curvature_sample_3d *out_samples,
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
			curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
			&result);
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
	qaws_curve const *curve,
	qaws_vec2 point,
	int *out_winding_number)
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
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
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
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
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

/* Deduplication helper: checks both parameter proximity and position proximity */
static int isect_2d_is_duplicate(
	qaws_intersection_2d const *buf, unsigned int count, unsigned int capacity,
	qaws_scalar ta, qaws_scalar tb)
{
	unsigned int i;
	unsigned int check = count < capacity ? count : capacity;
	for (i = 0; i < check; ++i)
	{
		qaws_scalar da = buf[i].parameter_a - ta;
		qaws_scalar db = buf[i].parameter_b - tb;
		if (da < 0) da = -da;
		if (db < 0) db = -db;
		if (da < ISECT_DEDUP_TOLERANCE && db < ISECT_DEDUP_TOLERANCE)
			return 1;
	}
	return 0;
}

/* Position-based dedup for self-intersections: different parameter pairs
   can converge to the same geometric crossing on multi-span curves. */
static int isect_2d_is_pos_duplicate(
	qaws_intersection_2d const *buf, unsigned int count, unsigned int capacity,
	qaws_vec2 pos)
{
	unsigned int i;
	unsigned int check = count < capacity ? count : capacity;
	for (i = 0; i < check; ++i)
	{
		qaws_scalar dx = buf[i].position.x - pos.x;
		qaws_scalar dy = buf[i].position.y - pos.y;
		if (dx < 0) dx = -dx;
		if (dy < 0) dy = -dy;
		if (dx < ISECT_POS_TOLERANCE * (qaws_scalar)10.0 &&
			dy < ISECT_POS_TOLERANCE * (qaws_scalar)10.0)
			return 1;
	}
	return 0;
}

static int isect_3d_is_duplicate(
	qaws_intersection_3d const *buf, unsigned int count, unsigned int capacity,
	qaws_scalar ta, qaws_scalar tb)
{
	unsigned int i;
	unsigned int check = count < capacity ? count : capacity;
	for (i = 0; i < check; ++i)
	{
		qaws_scalar da = buf[i].parameter_a - ta;
		qaws_scalar db = buf[i].parameter_b - tb;
		if (da < 0) da = -da;
		if (db < 0) db = -db;
		if (da < ISECT_DEDUP_TOLERANCE && db < ISECT_DEDUP_TOLERANCE)
			return 1;
	}
	return 0;
}

/*
 * Newton-Raphson refinement for 2D curve-curve intersection.
 * Solves C_a(ta) - C_b(tb) = 0 using the 2x2 Jacobian [C_a', -C_b'].
 * Returns 1 if converged, 0 otherwise. ta/tb are updated in place.
 */
static int newton_refine_2d(
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	qaws_scalar *ta, qaws_scalar *tb,
	qaws_scalar a_min, qaws_scalar a_max,
	qaws_scalar b_min, qaws_scalar b_max)
{
	unsigned int iter;
	for (iter = 0; iter < ISECT_NEWTON_ITERS; ++iter)
	{
		qaws_eval_result_2d ra, rb;
		qaws_scalar fx, fy;
		qaws_scalar j00, j01, j10, j11;
		qaws_scalar det, inv_det;
		qaws_scalar dta, dtb;

		if (qaws_curve_evaluate_2d(curve_a, *ta,
				QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &ra) != QAWS_STATUS_OK)
			return 0;
		if (qaws_curve_evaluate_2d(curve_b, *tb,
				QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &rb) != QAWS_STATUS_OK)
			return 0;

		fx = ra.position.x - rb.position.x;
		fy = ra.position.y - rb.position.y;

		if (fx < 0) { if (-fx < ISECT_CONVERGE_TOL && (fy < 0 ? -fy : fy) < ISECT_CONVERGE_TOL) return 1; }
		else        { if ( fx < ISECT_CONVERGE_TOL && (fy < 0 ? -fy : fy) < ISECT_CONVERGE_TOL) return 1; }

		/* Jacobian: [da.x, -db.x; da.y, -db.y] */
		j00 = ra.d1.x;  j01 = -rb.d1.x;
		j10 = ra.d1.y;  j11 = -rb.d1.y;

		det = j00 * j11 - j01 * j10;
		if (det < 0) det = -det;
		if (det < (qaws_scalar)1.0e-30)
			return 0;

		det = j00 * j11 - j01 * j10;
		inv_det = (qaws_scalar)1.0 / det;

		dta = ( j11 * fx - j01 * fy) * inv_det;
		dtb = (-j10 * fx + j00 * fy) * inv_det;

		*ta -= dta;
		*tb -= dtb;

		/* Clamp to domain */
		if (*ta < a_min) *ta = a_min;
		if (*ta > a_max) *ta = a_max;
		if (*tb < b_min) *tb = b_min;
		if (*tb > b_max) *tb = b_max;
	}

	/* Check final residual */
	{
		qaws_eval_result_2d ra, rb;
		qaws_scalar dx, dy;
		if (qaws_curve_evaluate_2d(curve_a, *ta, QAWS_EVAL_FLAG_POSITION, &ra) != QAWS_STATUS_OK)
			return 0;
		if (qaws_curve_evaluate_2d(curve_b, *tb, QAWS_EVAL_FLAG_POSITION, &rb) != QAWS_STATUS_OK)
			return 0;
		dx = ra.position.x - rb.position.x;
		dy = ra.position.y - rb.position.y;
		if (dx < 0) dx = -dx;
		if (dy < 0) dy = -dy;
		return (dx < ISECT_POS_TOLERANCE && dy < ISECT_POS_TOLERANCE) ? 1 : 0;
	}
}

static int newton_refine_3d(
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	qaws_scalar *ta, qaws_scalar *tb,
	qaws_scalar a_min, qaws_scalar a_max,
	qaws_scalar b_min, qaws_scalar b_max)
{
	unsigned int iter;
	for (iter = 0; iter < ISECT_NEWTON_ITERS; ++iter)
	{
		qaws_eval_result_3d ra, rb;
		qaws_scalar fx, fy, fz;
		qaws_scalar ax, ay, az, bx, by, bz;
		qaws_scalar ata, atb, dta, dtb;

		if (qaws_curve_evaluate_3d(curve_a, *ta,
				QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &ra) != QAWS_STATUS_OK)
			return 0;
		if (qaws_curve_evaluate_3d(curve_b, *tb,
				QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &rb) != QAWS_STATUS_OK)
			return 0;

		fx = ra.position.x - rb.position.x;
		fy = ra.position.y - rb.position.y;
		fz = ra.position.z - rb.position.z;

		{
			qaws_scalar afx = fx < 0 ? -fx : fx;
			qaws_scalar afy = fy < 0 ? -fy : fy;
			qaws_scalar afz = fz < 0 ? -fz : fz;
			if (afx < ISECT_CONVERGE_TOL && afy < ISECT_CONVERGE_TOL
				&& afz < ISECT_CONVERGE_TOL)
				return 1;
		}

		/* 3D overdetermined: least-squares via normal equations
		   J = [da, -db] (3x2), solve J^T J [dta;dtb] = J^T f */
		ax = ra.d1.x; ay = ra.d1.y; az = ra.d1.z;
		bx = -rb.d1.x; by = -rb.d1.y; bz = -rb.d1.z;

		ata = ax * ax + ay * ay + az * az;
		atb = ax * bx + ay * by + az * bz;
		{
			qaws_scalar btb = bx * bx + by * by + bz * bz;
			qaws_scalar rhs_a = ax * fx + ay * fy + az * fz;
			qaws_scalar rhs_b = bx * fx + by * fy + bz * fz;
			qaws_scalar det = ata * btb - atb * atb;
			qaws_scalar abs_det = det < 0 ? -det : det;
			qaws_scalar inv_det;

			if (abs_det < (qaws_scalar)1.0e-30)
				return 0;

			inv_det = (qaws_scalar)1.0 / det;
			dta = ( btb * rhs_a - atb * rhs_b) * inv_det;
			dtb = (-atb * rhs_a + ata * rhs_b) * inv_det;
		}

		*ta -= dta;
		*tb -= dtb;

		if (*ta < a_min) *ta = a_min;
		if (*ta > a_max) *ta = a_max;
		if (*tb < b_min) *tb = b_min;
		if (*tb > b_max) *tb = b_max;
	}

	{
		qaws_eval_result_3d ra, rb;
		qaws_scalar dx, dy, dz;
		if (qaws_curve_evaluate_3d(curve_a, *ta, QAWS_EVAL_FLAG_POSITION, &ra) != QAWS_STATUS_OK)
			return 0;
		if (qaws_curve_evaluate_3d(curve_b, *tb, QAWS_EVAL_FLAG_POSITION, &rb) != QAWS_STATUS_OK)
			return 0;
		dx = ra.position.x - rb.position.x;
		dy = ra.position.y - rb.position.y;
		dz = ra.position.z - rb.position.z;
		if (dx < 0) dx = -dx;
		if (dy < 0) dy = -dy;
		if (dz < 0) dz = -dz;
		return (dx < ISECT_POS_TOLERANCE && dy < ISECT_POS_TOLERANCE
			&& dz < ISECT_POS_TOLERANCE) ? 1 : 0;
	}
}

/* ------------------------------------------------------------------ */
/*  Self-intersection detection                                        */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_find_self_intersections_2d(
	qaws_curve const *curve,
	qaws_intersection_2d *out_intersections,
	unsigned int intersection_capacity,
	unsigned int *out_count)
{
	unsigned int i, j, n, found;
	qaws_scalar range_min, range_max;
	qaws_scalar *params = NULL;
	qaws_vec2 *pts = NULL;
	qaws_status s;

	if (!curve || !out_intersections || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_count = 0;
	found = 0;

	n = ISECT_SAMPLE_COUNT;
	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	params = (qaws_scalar *)malloc(n * sizeof(qaws_scalar));
	pts = (qaws_vec2 *)malloc(n * sizeof(qaws_vec2));
	if (!params || !pts) { free(params); free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }

	/* Dense sampling */
	for (i = 0; i < n; ++i)
	{
		qaws_eval_result_2d r;
		params[i] = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(n - 1);
		s = qaws_curve_evaluate_2d(curve, params[i], QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) { free(params); free(pts); return s; }
		pts[i] = r.position;
	}

	/* Find candidate pairs: samples far apart in parameter but close in position.
	   Minimum index gap to avoid trivial near-neighbors. */
	{
		unsigned int min_gap = n / 8;
		if (min_gap < 4) min_gap = 4;

		for (i = 0; i < n; ++i)
		{
			for (j = i + min_gap; j < n; ++j)
			{
				qaws_scalar dx = pts[i].x - pts[j].x;
				qaws_scalar dy = pts[i].y - pts[j].y;
				qaws_scalar dist2;
				if (dx < 0) dx = -dx;
				if (dy < 0) dy = -dy;

				/* Coarse distance threshold: ~2x the sampling chord length */
				dist2 = dx * dx + dy * dy;
				if (dist2 > ISECT_CANDIDATE_TOL * ISECT_CANDIDATE_TOL * (qaws_scalar)10000.0)
					continue;

				/* Newton refinement */
				{
					qaws_scalar ta = params[i];
					qaws_scalar tb = params[j];

					if (newton_refine_2d(curve, curve, &ta, &tb,
						range_min, range_max, range_min, range_max))
					{
						/* Ensure ta < tb for self-intersection */
						if (ta > tb)
						{
							qaws_scalar tmp = ta; ta = tb; tb = tmp;
						}
						/* Parameter gap must be significant (not a trivial point) */
						if ((tb - ta) > ISECT_DEDUP_TOLERANCE * (qaws_scalar)10.0)
						{
							if (!isect_2d_is_duplicate(out_intersections, found, intersection_capacity, ta, tb))
							{
								qaws_eval_result_2d rp;
								qaws_curve_evaluate_2d(curve, ta,
									QAWS_EVAL_FLAG_POSITION, &rp);
								/* Also check position-based dedup for multi-span curves */
								if (!isect_2d_is_pos_duplicate(out_intersections, found, intersection_capacity, rp.position))
								{
									if (found < intersection_capacity)
									{
										out_intersections[found].parameter_a = ta;
										out_intersections[found].parameter_b = tb;
										out_intersections[found].position = rp.position;
									}
									++found;
								}
							}
						}
					}
				}
			}
		}
	}

	free(params);
	free(pts);
	*out_count = found;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_find_self_intersections_3d(
	qaws_curve const *curve,
	qaws_intersection_3d *out_intersections,
	unsigned int intersection_capacity,
	unsigned int *out_count)
{
	unsigned int i, j, n, found;
	qaws_scalar range_min, range_max;
	qaws_scalar *params = NULL;
	qaws_vec3 *pts = NULL;
	qaws_status s;

	if (!curve || !out_intersections || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_count = 0;
	found = 0;

	n = ISECT_SAMPLE_COUNT;
	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	params = (qaws_scalar *)malloc(n * sizeof(qaws_scalar));
	pts = (qaws_vec3 *)malloc(n * sizeof(qaws_vec3));
	if (!params || !pts) { free(params); free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }

	for (i = 0; i < n; ++i)
	{
		qaws_eval_result_3d r;
		params[i] = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(n - 1);
		s = qaws_curve_evaluate_3d(curve, params[i], QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) { free(params); free(pts); return s; }
		pts[i] = r.position;
	}

	{
		unsigned int min_gap = n / 8;
		if (min_gap < 4) min_gap = 4;

		for (i = 0; i < n; ++i)
		{
			for (j = i + min_gap; j < n; ++j)
			{
				qaws_scalar dx = pts[i].x - pts[j].x;
				qaws_scalar dy = pts[i].y - pts[j].y;
				qaws_scalar dz = pts[i].z - pts[j].z;
				qaws_scalar dist2;

				dist2 = dx * dx + dy * dy + dz * dz;
				if (dist2 > ISECT_CANDIDATE_TOL * ISECT_CANDIDATE_TOL * (qaws_scalar)10000.0)
					continue;

				{
					qaws_scalar ta = params[i];
					qaws_scalar tb = params[j];

					if (newton_refine_3d(curve, curve, &ta, &tb,
						range_min, range_max, range_min, range_max))
					{
						if (ta > tb)
						{
							qaws_scalar tmp = ta; ta = tb; tb = tmp;
						}
						if ((tb - ta) > ISECT_DEDUP_TOLERANCE * (qaws_scalar)10.0)
						{
							if (!isect_3d_is_duplicate(out_intersections, found, intersection_capacity, ta, tb))
							{
								if (found < intersection_capacity)
								{
									qaws_eval_result_3d rp;
									qaws_curve_evaluate_3d(curve, ta,
										QAWS_EVAL_FLAG_POSITION, &rp);
									out_intersections[found].parameter_a = ta;
									out_intersections[found].parameter_b = tb;
									out_intersections[found].position = rp.position;
								}
								++found;
							}
						}
					}
				}
			}
		}
	}

	free(params);
	free(pts);
	*out_count = found;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Curve-curve intersection                                           */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_find_intersections_2d(
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	qaws_intersection_2d *out_intersections,
	unsigned int intersection_capacity,
	unsigned int *out_count)
{
	unsigned int na, nb, i, j, found;
	qaws_scalar a_min, a_max, b_min, b_max;
	qaws_scalar *params_a = NULL, *params_b = NULL;
	qaws_vec2 *pts_a = NULL, *pts_b = NULL;
	qaws_status s;

	if (!curve_a || !curve_b || !out_intersections || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve_a->dimension != QAWS_DIMENSION_2D
		|| curve_b->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_count = 0;
	found = 0;

	na = ISECT_SAMPLE_COUNT;
	nb = ISECT_SAMPLE_COUNT;
	a_min = curve_a->parameter_range.min_value;
	a_max = curve_a->parameter_range.max_value;
	b_min = curve_b->parameter_range.min_value;
	b_max = curve_b->parameter_range.max_value;

	params_a = (qaws_scalar *)malloc(na * sizeof(qaws_scalar));
	pts_a = (qaws_vec2 *)malloc(na * sizeof(qaws_vec2));
	params_b = (qaws_scalar *)malloc(nb * sizeof(qaws_scalar));
	pts_b = (qaws_vec2 *)malloc(nb * sizeof(qaws_vec2));
	if (!params_a || !pts_a || !params_b || !pts_b)
	{
		free(params_a); free(pts_a); free(params_b); free(pts_b);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Sample both curves */
	for (i = 0; i < na; ++i)
	{
		qaws_eval_result_2d r;
		params_a[i] = a_min + (qaws_scalar)i * (a_max - a_min) / (qaws_scalar)(na - 1);
		s = qaws_curve_evaluate_2d(curve_a, params_a[i], QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) goto cleanup_2d;
		pts_a[i] = r.position;
	}
	for (j = 0; j < nb; ++j)
	{
		qaws_eval_result_2d r;
		params_b[j] = b_min + (qaws_scalar)j * (b_max - b_min) / (qaws_scalar)(nb - 1);
		s = qaws_curve_evaluate_2d(curve_b, params_b[j], QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) goto cleanup_2d;
		pts_b[j] = r.position;
	}

	/* Find candidates: nearby sample pairs from different curves */
	for (i = 0; i < na; ++i)
	{
		for (j = 0; j < nb; ++j)
		{
			qaws_scalar dx = pts_a[i].x - pts_b[j].x;
			qaws_scalar dy = pts_a[i].y - pts_b[j].y;
			qaws_scalar dist2 = dx * dx + dy * dy;

			if (dist2 > ISECT_CANDIDATE_TOL * ISECT_CANDIDATE_TOL * (qaws_scalar)10000.0)
				continue;

			{
				qaws_scalar ta = params_a[i];
				qaws_scalar tb = params_b[j];

				if (newton_refine_2d(curve_a, curve_b, &ta, &tb,
					a_min, a_max, b_min, b_max))
				{
					if (!isect_2d_is_duplicate(out_intersections, found, intersection_capacity, ta, tb))
					{
						qaws_eval_result_2d rp;
						qaws_curve_evaluate_2d(curve_a, ta,
							QAWS_EVAL_FLAG_POSITION, &rp);
						if (!isect_2d_is_pos_duplicate(out_intersections, found, intersection_capacity, rp.position))
						{
							if (found < intersection_capacity)
							{
								out_intersections[found].parameter_a = ta;
								out_intersections[found].parameter_b = tb;
								out_intersections[found].position = rp.position;
							}
							++found;
						}
					}
				}
			}
		}
	}

	s = QAWS_STATUS_OK;

cleanup_2d:
	free(params_a); free(pts_a); free(params_b); free(pts_b);
	*out_count = found;
	return s;
}

qaws_status qaws_curve_find_intersections_3d(
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	qaws_intersection_3d *out_intersections,
	unsigned int intersection_capacity,
	unsigned int *out_count)
{
	unsigned int na, nb, i, j, found;
	qaws_scalar a_min, a_max, b_min, b_max;
	qaws_scalar *params_a = NULL, *params_b = NULL;
	qaws_vec3 *pts_a = NULL, *pts_b = NULL;
	qaws_status s;

	if (!curve_a || !curve_b || !out_intersections || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve_a->dimension != QAWS_DIMENSION_3D
		|| curve_b->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_count = 0;
	found = 0;

	na = ISECT_SAMPLE_COUNT;
	nb = ISECT_SAMPLE_COUNT;
	a_min = curve_a->parameter_range.min_value;
	a_max = curve_a->parameter_range.max_value;
	b_min = curve_b->parameter_range.min_value;
	b_max = curve_b->parameter_range.max_value;

	params_a = (qaws_scalar *)malloc(na * sizeof(qaws_scalar));
	pts_a = (qaws_vec3 *)malloc(na * sizeof(qaws_vec3));
	params_b = (qaws_scalar *)malloc(nb * sizeof(qaws_scalar));
	pts_b = (qaws_vec3 *)malloc(nb * sizeof(qaws_vec3));
	if (!params_a || !pts_a || !params_b || !pts_b)
	{
		free(params_a); free(pts_a); free(params_b); free(pts_b);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 0; i < na; ++i)
	{
		qaws_eval_result_3d r;
		params_a[i] = a_min + (qaws_scalar)i * (a_max - a_min) / (qaws_scalar)(na - 1);
		s = qaws_curve_evaluate_3d(curve_a, params_a[i], QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) goto cleanup_3d;
		pts_a[i] = r.position;
	}
	for (j = 0; j < nb; ++j)
	{
		qaws_eval_result_3d r;
		params_b[j] = b_min + (qaws_scalar)j * (b_max - b_min) / (qaws_scalar)(nb - 1);
		s = qaws_curve_evaluate_3d(curve_b, params_b[j], QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) goto cleanup_3d;
		pts_b[j] = r.position;
	}

	for (i = 0; i < na; ++i)
	{
		for (j = 0; j < nb; ++j)
		{
			qaws_scalar dx = pts_a[i].x - pts_b[j].x;
			qaws_scalar dy = pts_a[i].y - pts_b[j].y;
			qaws_scalar dz = pts_a[i].z - pts_b[j].z;
			qaws_scalar dist2 = dx * dx + dy * dy + dz * dz;

			if (dist2 > ISECT_CANDIDATE_TOL * ISECT_CANDIDATE_TOL * (qaws_scalar)10000.0)
				continue;

			{
				qaws_scalar ta = params_a[i];
				qaws_scalar tb = params_b[j];

				if (newton_refine_3d(curve_a, curve_b, &ta, &tb,
					a_min, a_max, b_min, b_max))
				{
					if (!isect_3d_is_duplicate(out_intersections, found, intersection_capacity, ta, tb))
					{
						if (found < intersection_capacity)
						{
							qaws_eval_result_3d rp;
							qaws_curve_evaluate_3d(curve_a, ta,
								QAWS_EVAL_FLAG_POSITION, &rp);
							out_intersections[found].parameter_a = ta;
							out_intersections[found].parameter_b = tb;
							out_intersections[found].position = rp.position;
						}
						++found;
					}
				}
			}
		}
	}

	s = QAWS_STATUS_OK;

cleanup_3d:
	free(params_a); free(pts_a); free(params_b); free(pts_b);
	*out_count = found;
	return s;
}

/* ------------------------------------------------------------------ */
/*  Surface inspection                                                 */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_compute_bounds(
	qaws_surface const *surface,
	qaws_vec3 *out_min,
	qaws_vec3 *out_max)
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
		surface, u_min, v_min, QAWS_SURFACE_EVAL_POSITION, &result);
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
				surface, u, v, QAWS_SURFACE_EVAL_POSITION, &result);
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
	qaws_surface const *surface,
	qaws_scalar *out_area)
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
				surface, u, v,
				QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV,
				&result);
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
	qaws_surface const *surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_scalar *out_curvature)
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
		surface, u, v,
		QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
		| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV
		| QAWS_SURFACE_EVAL_DVV,
		&result);
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
	qaws_surface const *surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_scalar *out_curvature)
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
		surface, u, v,
		QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
		| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV
		| QAWS_SURFACE_EVAL_DVV,
		&result);
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
	qaws_surface const *surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_surface_curvature_result *out_result)
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
		surface, u, v,
		QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
		| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV
		| QAWS_SURFACE_EVAL_DVV,
		&result);
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
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	unsigned int sample_count,
	qaws_scalar *out_distance)
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
			curve_a, t, QAWS_EVAL_FLAG_POSITION, &result_a);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_a.position;

		/* Find closest point on curve_b */
		status = qaws_curve_find_closest_parameter_2d(
			curve_b, pt, &param);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_2d(
			curve_b, param, QAWS_EVAL_FLAG_POSITION, &result_b);
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
			curve_b, t, QAWS_EVAL_FLAG_POSITION, &result_b);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_b.position;

		/* Find closest point on curve_a */
		status = qaws_curve_find_closest_parameter_2d(
			curve_a, pt, &param);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_2d(
			curve_a, param, QAWS_EVAL_FLAG_POSITION, &result_a);
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
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	unsigned int sample_count,
	qaws_scalar *out_distance)
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
			curve_a, t, QAWS_EVAL_FLAG_POSITION, &result_a);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_a.position;

		/* Find closest point on curve_b */
		status = qaws_curve_find_closest_parameter_3d(
			curve_b, pt, &param);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_3d(
			curve_b, param, QAWS_EVAL_FLAG_POSITION, &result_b);
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
			curve_b, t, QAWS_EVAL_FLAG_POSITION, &result_b);
		if (status != QAWS_STATUS_OK)
			return status;

		pt = result_b.position;

		/* Find closest point on curve_a */
		status = qaws_curve_find_closest_parameter_3d(
			curve_a, pt, &param);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_3d(
			curve_a, param, QAWS_EVAL_FLAG_POSITION, &result_a);
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

qaws_status qaws_surface_find_closest_point(
	qaws_surface const *surface,
	qaws_vec3 point,
	qaws_scalar *out_u,
	qaws_scalar *out_v,
	qaws_vec3 *out_closest_point)
{
	unsigned int i;
	unsigned int j;
	unsigned int iter;
	qaws_scalar u;
	qaws_scalar v;
	qaws_scalar best_u;
	qaws_scalar best_v;
	qaws_scalar best_dist_sq;
	qaws_scalar dx, dy, dz;
	qaws_scalar dist_sq;
	qaws_scalar u_min, u_max, v_min, v_max;
	qaws_scalar g_u, g_v;
	qaws_scalar h_uu, h_uv, h_vv;
	qaws_scalar det;
	qaws_scalar delta_u, delta_v;
	qaws_scalar delta_len;
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

	/* Coarse search: 8x8 grid */
	best_u = u_min;
	best_v = v_min;
	best_dist_sq = QAWS_LITERAL(1.0e30);

	for (i = 0; i < 8; ++i)
	{
		u = u_min + (qaws_scalar)i * (u_max - u_min)
			/ QAWS_LITERAL(7.0);

		for (j = 0; j < 8; ++j)
		{
			v = v_min + (qaws_scalar)j * (v_max - v_min)
				/ QAWS_LITERAL(7.0);

			status = qaws_surface_evaluate(
				surface, u, v, QAWS_SURFACE_EVAL_POSITION, &result);
			if (status != QAWS_STATUS_OK)
				return status;

			dx = result.position.x - point.x;
			dy = result.position.y - point.y;
			dz = result.position.z - point.z;
			dist_sq = dx * dx + dy * dy + dz * dz;

			if (dist_sq < best_dist_sq)
			{
				best_dist_sq = dist_sq;
				best_u = u;
				best_v = v;
			}
		}
	}

	/* Newton refinement */
	u = best_u;
	v = best_v;

	for (iter = 0; iter < 20; ++iter)
	{
		qaws_vec3 r;

		status = qaws_surface_evaluate(
			surface, u, v,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU
			| QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DUU
			| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV,
			&result);
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
		if (u < u_min) u = u_min;
		if (u > u_max) u = u_max;
		if (v < v_min) v = v_min;
		if (v > v_max) v = v_max;

		delta_len = QAWS_SQRT(delta_u * delta_u + delta_v * delta_v);
		if (delta_len < QAWS_LITERAL(1.0e-8))
			break;
	}

	/* Evaluate final position */
	status = qaws_surface_evaluate(
		surface, u, v, QAWS_SURFACE_EVAL_POSITION, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	*out_u = u;
	*out_v = v;
	*out_closest_point = result.position;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Curve-plane intersection                                           */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_find_plane_intersections(
	qaws_curve const *curve,
	qaws_plane const *plane,
	qaws_scalar *out_parameters,
	qaws_vec3 *out_positions,
	unsigned int capacity,
	unsigned int *out_count)
{
	unsigned int i;
	unsigned int count;
	unsigned int bisect_iter;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_scalar t;
	qaws_scalar t_prev;
	qaws_scalar d_prev;
	qaws_scalar d_curr;
	qaws_scalar t_lo, t_hi;
	qaws_scalar d_lo;
	qaws_scalar t_mid, d_mid;
	qaws_scalar nx, ny, nz;
	qaws_eval_result_3d result;
	qaws_status status;
	unsigned int n_samples;

	if (!curve || !plane || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;
	nx = plane->normal.x;
	ny = plane->normal.y;
	nz = plane->normal.z;
	n_samples = 128;
	count = 0;

	/* Evaluate first sample */
	memset(&result, 0, sizeof(result));
	status = qaws_curve_evaluate_3d(
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	d_prev = (result.position.x - plane->point.x) * nx
		+ (result.position.y - plane->point.y) * ny
		+ (result.position.z - plane->point.z) * nz;
	t_prev = range_min;

	for (i = 1; i <= n_samples; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)n_samples;

		memset(&result, 0, sizeof(result));
		status = qaws_curve_evaluate_3d(
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		d_curr = (result.position.x - plane->point.x) * nx
			+ (result.position.y - plane->point.y) * ny
			+ (result.position.z - plane->point.z) * nz;

		/* Exact hit: d_curr is on the plane */
		if (QAWS_FABS(d_curr) < QAWS_LITERAL(1.0e-8))
		{
			if (count < capacity
				&& QAWS_FABS(d_prev) >= QAWS_LITERAL(1.0e-8))
			{
				if (out_parameters)
					out_parameters[count] = t;
				if (out_positions)
					out_positions[count] = result.position;
				++count;
			}
		}
		/* Sign change: root lies between t_prev and t */
		else if ((d_prev > QAWS_LITERAL(0.0) && d_curr < QAWS_LITERAL(0.0))
			|| (d_prev < QAWS_LITERAL(0.0) && d_curr > QAWS_LITERAL(0.0)))
		{
			if (count >= capacity)
				break;

			/* Bisection refinement */
			t_lo = t_prev;
			t_hi = t;
			d_lo = d_prev;

			for (bisect_iter = 0; bisect_iter < 50; ++bisect_iter)
			{
				t_mid = QAWS_LITERAL(0.5) * (t_lo + t_hi);

				status = qaws_curve_evaluate_3d(
					curve, t_mid, QAWS_EVAL_FLAG_POSITION, &result);
				if (status != QAWS_STATUS_OK)
					return status;

				d_mid = (result.position.x - plane->point.x) * nx
					+ (result.position.y - plane->point.y) * ny
					+ (result.position.z - plane->point.z) * nz;

				if (QAWS_FABS(d_mid) < QAWS_LITERAL(1.0e-8))
					break;

				if ((d_lo > QAWS_LITERAL(0.0) && d_mid > QAWS_LITERAL(0.0))
					|| (d_lo < QAWS_LITERAL(0.0)
						&& d_mid < QAWS_LITERAL(0.0)))
				{
					t_lo = t_mid;
					d_lo = d_mid;
				}
				else
				{
					t_hi = t_mid;
				}

				if ((t_hi - t_lo) < QAWS_LITERAL(1.0e-8))
					break;
			}

			t_mid = QAWS_LITERAL(0.5) * (t_lo + t_hi);

			/* Evaluate final position */
			status = qaws_curve_evaluate_3d(
				curve, t_mid, QAWS_EVAL_FLAG_POSITION, &result);
			if (status != QAWS_STATUS_OK)
				return status;

			if (out_parameters)
				out_parameters[count] = t_mid;
			if (out_positions)
				out_positions[count] = result.position;
			++count;
		}

		d_prev = d_curr;
		t_prev = t;
	}

	*out_count = count;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Surface-curve intersection                                         */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_find_curve_intersections(
	qaws_surface const *surface,
	qaws_curve const *curve,
	qaws_surface_curve_intersection *out_intersections,
	unsigned int capacity,
	unsigned int *out_count)
{
	unsigned int i;
	unsigned int gi, gj;
	unsigned int iter;
	unsigned int count;
	unsigned int k;
	int duplicate;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_scalar u_min, u_max, v_min, v_max;
	qaws_scalar u, v;
	qaws_scalar best_u, best_v;
	qaws_scalar best_dist_sq;
	qaws_scalar dx, dy, dz;
	qaws_scalar dist_sq;
	qaws_scalar fx, fy, fz;
	qaws_scalar det;
	qaws_scalar du, dv, dt;
	qaws_scalar delta_len;
	qaws_range u_range;
	qaws_range v_range;
	qaws_eval_result_3d curve_result;
	qaws_surface_eval_result surf_result;
	qaws_status status;
	unsigned int n_curve_samples;

	/* Jacobian columns and Cramer's rule temporaries */
	qaws_scalar j00, j01, j02;
	qaws_scalar j10, j11, j12;
	qaws_scalar j20, j21, j22;
	qaws_scalar c0x, c0y, c0z;
	qaws_scalar c1x, c1y, c1z;
	qaws_scalar c2x, c2y, c2z;

	if (!surface || !curve || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;
	u_range = qaws_surface_get_u_range(surface);
	v_range = qaws_surface_get_v_range(surface);
	u_min = u_range.min_value;
	u_max = u_range.max_value;
	v_min = v_range.min_value;
	v_max = v_range.max_value;

	n_curve_samples = 64;
	count = 0;

	for (i = 0; i < n_curve_samples; ++i)
	{
		if (count >= capacity)
			break;

		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(n_curve_samples - 1);

		status = qaws_curve_evaluate_3d(
			curve, t, QAWS_EVAL_FLAG_POSITION, &curve_result);
		if (status != QAWS_STATUS_OK)
			return status;

		/* Quick 8x8 grid search on surface for closest point */
		best_u = u_min;
		best_v = v_min;
		best_dist_sq = QAWS_LITERAL(1.0e30);

		for (gi = 0; gi < 8; ++gi)
		{
			u = u_min + (qaws_scalar)gi * (u_max - u_min)
				/ QAWS_LITERAL(7.0);

			for (gj = 0; gj < 8; ++gj)
			{
				v = v_min + (qaws_scalar)gj * (v_max - v_min)
					/ QAWS_LITERAL(7.0);

				status = qaws_surface_evaluate(
					surface, u, v, QAWS_SURFACE_EVAL_POSITION,
					&surf_result);
				if (status != QAWS_STATUS_OK)
					return status;

				dx = surf_result.position.x - curve_result.position.x;
				dy = surf_result.position.y - curve_result.position.y;
				dz = surf_result.position.z - curve_result.position.z;
				dist_sq = dx * dx + dy * dy + dz * dz;

				if (dist_sq < best_dist_sq)
				{
					best_dist_sq = dist_sq;
					best_u = u;
					best_v = v;
				}
			}
		}

		/* Only proceed to Newton if coarse distance is small enough.
		   Use generous threshold; Newton convergence is checked after. */
		if (best_dist_sq > QAWS_LITERAL(4.0))
			continue;

		/* Newton iteration: solve S(u,v) - C(t) = 0 */
		u = best_u;
		v = best_v;

		for (iter = 0; iter < 20; ++iter)
		{
			status = qaws_surface_evaluate(
				surface, u, v,
				QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU
				| QAWS_SURFACE_EVAL_DV,
				&surf_result);
			if (status != QAWS_STATUS_OK)
				return status;

			status = qaws_curve_evaluate_3d(
				curve, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1,
				&curve_result);
			if (status != QAWS_STATUS_OK)
				return status;

			/* F = S(u,v) - C(t) */
			fx = surf_result.position.x - curve_result.position.x;
			fy = surf_result.position.y - curve_result.position.y;
			fz = surf_result.position.z - curve_result.position.z;

			delta_len = QAWS_SQRT(fx * fx + fy * fy + fz * fz);
			if (delta_len < QAWS_LITERAL(1.0e-8))
				break;

			/* Jacobian: J = [S_u, S_v, -C'(t)] */
			j00 = surf_result.du.x;
			j10 = surf_result.du.y;
			j20 = surf_result.du.z;
			j01 = surf_result.dv.x;
			j11 = surf_result.dv.y;
			j21 = surf_result.dv.z;
			j02 = -curve_result.d1.x;
			j12 = -curve_result.d1.y;
			j22 = -curve_result.d1.z;

			/* Determinant via Cramer's rule */
			c0x = j11 * j22 - j21 * j12;
			c0y = j21 * j02 - j01 * j22;
			c0z = j01 * j12 - j11 * j02;

			det = j00 * c0x + j10 * c0y + j20 * c0z;
			if (QAWS_FABS(det) < QAWS_LITERAL(1.0e-30))
				break;

			/* Columns for right-hand side substitution */
			c1x = j10 * j22 - j20 * j12;
			c1y = j20 * j02 - j00 * j22;
			c1z = j00 * j12 - j10 * j02;

			c2x = j10 * j21 - j20 * j11;
			c2y = j20 * j01 - j00 * j21;
			c2z = j00 * j11 - j10 * j01;

			du = (fx * c0x + fy * c0y + fz * c0z) / det;
			dv = -(fx * c1x + fy * c1y + fz * c1z) / det;
			dt = (fx * c2x + fy * c2y + fz * c2z) / det;

			u -= du;
			v -= dv;
			t -= dt;

			/* Clamp to parameter domains */
			if (u < u_min) u = u_min;
			if (u > u_max) u = u_max;
			if (v < v_min) v = v_min;
			if (v > v_max) v = v_max;
			if (t < range_min) t = range_min;
			if (t > range_max) t = range_max;

			delta_len = QAWS_SQRT(du * du + dv * dv + dt * dt);
			if (delta_len < QAWS_LITERAL(1.0e-8))
				break;
		}

		/* Check convergence */
		status = qaws_surface_evaluate(
			surface, u, v, QAWS_SURFACE_EVAL_POSITION, &surf_result);
		if (status != QAWS_STATUS_OK)
			return status;

		status = qaws_curve_evaluate_3d(
			curve, t, QAWS_EVAL_FLAG_POSITION, &curve_result);
		if (status != QAWS_STATUS_OK)
			return status;

		fx = surf_result.position.x - curve_result.position.x;
		fy = surf_result.position.y - curve_result.position.y;
		fz = surf_result.position.z - curve_result.position.z;
		dist_sq = fx * fx + fy * fy + fz * fz;

		if (dist_sq > QAWS_LITERAL(1.0e-6))
			continue;

		/* Deduplicate: merge results closer than 1e-4 in parameter space */
		duplicate = 0;
		for (k = 0; k < count; ++k)
		{
			qaws_scalar du2, dv2, dt2;
			du2 = u - out_intersections[k].u;
			dv2 = v - out_intersections[k].v;
			dt2 = t - out_intersections[k].t;

			if (QAWS_SQRT(du2 * du2 + dv2 * dv2 + dt2 * dt2)
				< QAWS_LITERAL(1.0e-4))
			{
				duplicate = 1;
				break;
			}
		}

		if (!duplicate && count < capacity)
		{
			out_intersections[count].u = u;
			out_intersections[count].v = v;
			out_intersections[count].t = t;
			out_intersections[count].position = surf_result.position;
			++count;
		}
	}

	*out_count = count;
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
	*out_status = qaws_surface_evaluate(surface, u, v,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
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
	qaws_surface const *surface,
	qaws_tessellation_desc const *desc,
	qaws_tessellation_vertex *out_vertices,
	unsigned int vertex_capacity,
	unsigned int *out_vertex_count,
	unsigned int *out_indices,
	unsigned int index_capacity,
	unsigned int *out_index_count)
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
	stack = (qaws_tess_quad *)malloc(
		stack_capacity * sizeof(qaws_tess_quad));
	leaves = (qaws_tess_leaf *)malloc(
		leaf_capacity * sizeof(qaws_tess_leaf));
	if (!stack || !leaves)
	{
		free(stack);
		free(leaves);
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
		status = qaws_surface_evaluate(surface, quad.u0, quad.v0,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL,
			&r00);
		if (status != QAWS_STATUS_OK)
		{
			free(stack); free(leaves); return status;
		}

		status = qaws_surface_evaluate(surface, quad.u1, quad.v0,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL,
			&r10);
		if (status != QAWS_STATUS_OK)
		{
			free(stack); free(leaves); return status;
		}

		status = qaws_surface_evaluate(surface, quad.u0, quad.v1,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL,
			&r01);
		if (status != QAWS_STATUS_OK)
		{
			free(stack); free(leaves); return status;
		}

		status = qaws_surface_evaluate(surface, quad.u1, quad.v1,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL,
			&r11);
		if (status != QAWS_STATUS_OK)
		{
			free(stack); free(leaves); return status;
		}

		status = qaws_surface_evaluate(surface, u_mid, v_mid,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL,
			&r_center);
		if (status != QAWS_STATUS_OK)
		{
			free(stack); free(leaves); return status;
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

		if (subdivide && stack_top + 4 <= stack_capacity)
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
			/* Store as leaf quad */
			if (leaf_count < leaf_capacity)
			{
				leaves[leaf_count].u0 = quad.u0;
				leaves[leaf_count].v0 = quad.v0;
				leaves[leaf_count].u1 = quad.u1;
				leaves[leaf_count].v1 = quad.v1;
				++leaf_count;
			}
		}
	}

	free(stack);

	/* ---- Phase 2: Add all leaf corner vertices ---- */
	status = QAWS_STATUS_OK;
	for (li = 0; li < leaf_count; ++li)
	{
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u0, leaves[li].v0, &status);
		if (status != QAWS_STATUS_OK) { free(leaves); return status; }
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u1, leaves[li].v0, &status);
		if (status != QAWS_STATUS_OK) { free(leaves); return status; }
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u0, leaves[li].v1, &status);
		if (status != QAWS_STATUS_OK) { free(leaves); return status; }
		qaws_tess_eval_add(surface, out_vertices, &vert_count,
			vertex_capacity, leaves[li].u1, leaves[li].v1, &status);
		if (status != QAWS_STATUS_OK) { free(leaves); return status; }
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
		if (status != QAWS_STATUS_OK) { free(leaves); return status; }

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

	free(leaves);

	*out_vertex_count = vert_count;
	*out_index_count = idx_count;
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Curve projection onto surface (#20)                                       */
/* ========================================================================== */

qaws_status qaws_surface_project_curve(
	qaws_surface const* surface,
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_curve** out_uv_curve)
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

	params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samples);
	uv_coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samples * 2);
	if (!params || !uv_coords)
	{
		free(params); free(uv_coords);
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
		status = qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION, &er);
		if (status != QAWS_STATUS_OK) { free(params); free(uv_coords); return status; }

		status = qaws_surface_find_closest_point(surface, er.position, &u, &v, &closest);
		if (status != QAWS_STATUS_OK) { free(params); free(uv_coords); return status; }

		params[i] = t;
		uv_coords[i * 2 + 0] = u;
		uv_coords[i * 2 + 1] = v;
	}

	/* Fit a 2D B-spline through the (u,v) samples */
	n_cp = n_samples / 4;
	if (n_cp < 6) n_cp = 6;
	if (n_cp > n_samples) n_cp = n_samples;

	status = qaws_internal_fit_bspline(
		QAWS_DIMENSION_2D, 3,
		params, uv_coords,
		n_samples, n_cp,
		out_uv_curve);

	free(params);
	free(uv_coords);
	return status;
}

/* ========================================================================== */
/*  Geodesic curves on surfaces (#24)                                         */
/* ========================================================================== */

/* Compute Christoffel symbols from the first and second fundamental forms.
   Given surface eval with position, du, dv, duu, duv, dvv:
     E = du.du, F = du.dv, G = dv.dv
     e = N.duu, f = N.duv, g = N.dvv  (second fundamental form)
   Christoffel symbols (first kind):
     Gamma^1_{ij} from [E F; F G]^{-1} * coefficients of first fundamental form derivatives.
   We use the simpler formulation with geodesic equations in terms of metric. */

static void compute_christoffel(
	qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v,
	qaws_scalar gamma[2][2][2])
{
	qaws_surface_eval_result r;
	qaws_scalar E, F, G, det, inv_det;
	qaws_scalar dEdu, dEdv, dFdu, dFdv, dGdu, dGdv;
	qaws_surface_eval_result r_du, r_dv;
	qaws_scalar h = QAWS_LITERAL(1e-5);
	unsigned int eval_flags = QAWS_SURFACE_EVAL_POSITION
		| QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV;
	int i, j, k;

	memset(gamma, 0, sizeof(qaws_scalar) * 8);
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, u, v, eval_flags, &r);

	E = r.du.x * r.du.x + r.du.y * r.du.y + r.du.z * r.du.z;
	F = r.du.x * r.dv.x + r.du.y * r.dv.y + r.du.z * r.dv.z;
	G = r.dv.x * r.dv.x + r.dv.y * r.dv.y + r.dv.z * r.dv.z;

	det = E * G - F * F;
	if (QAWS_FABS(det) < QAWS_LITERAL(1e-20)) return;
	inv_det = QAWS_ONE / det;

	/* Numerical derivatives of metric coefficients */
	memset(&r_du, 0, sizeof(r_du));
	memset(&r_dv, 0, sizeof(r_dv));

	{
		qaws_surface_eval_result rp, rm;
		qaws_scalar Ep, Fp, Gp, Em, Fm, Gm;
		qaws_scalar u_clamp;

		/* dE/du, dF/du, dG/du via central difference in u */
		u_clamp = u + h;
		if (u_clamp > QAWS_ONE) u_clamp = QAWS_ONE;
		memset(&rp, 0, sizeof(rp));
		qaws_surface_evaluate(surface, u_clamp, v, eval_flags, &rp);

		u_clamp = u - h;
		if (u_clamp < QAWS_ZERO) u_clamp = QAWS_ZERO;
		memset(&rm, 0, sizeof(rm));
		qaws_surface_evaluate(surface, u_clamp, v, eval_flags, &rm);

		Ep = rp.du.x * rp.du.x + rp.du.y * rp.du.y + rp.du.z * rp.du.z;
		Fp = rp.du.x * rp.dv.x + rp.du.y * rp.dv.y + rp.du.z * rp.dv.z;
		Gp = rp.dv.x * rp.dv.x + rp.dv.y * rp.dv.y + rp.dv.z * rp.dv.z;
		Em = rm.du.x * rm.du.x + rm.du.y * rm.du.y + rm.du.z * rm.du.z;
		Fm = rm.du.x * rm.dv.x + rm.du.y * rm.dv.y + rm.du.z * rm.dv.z;
		Gm = rm.dv.x * rm.dv.x + rm.dv.y * rm.dv.y + rm.dv.z * rm.dv.z;

		dEdu = (Ep - Em) / (QAWS_LITERAL(2.0) * h);
		dFdu = (Fp - Fm) / (QAWS_LITERAL(2.0) * h);
		dGdu = (Gp - Gm) / (QAWS_LITERAL(2.0) * h);

		/* dE/dv, dF/dv, dG/dv via central difference in v */
		u_clamp = v + h;
		if (u_clamp > QAWS_ONE) u_clamp = QAWS_ONE;
		memset(&rp, 0, sizeof(rp));
		qaws_surface_evaluate(surface, u, u_clamp, eval_flags, &rp);

		u_clamp = v - h;
		if (u_clamp < QAWS_ZERO) u_clamp = QAWS_ZERO;
		memset(&rm, 0, sizeof(rm));
		qaws_surface_evaluate(surface, u, u_clamp, eval_flags, &rm);

		Ep = rp.du.x * rp.du.x + rp.du.y * rp.du.y + rp.du.z * rp.du.z;
		Fp = rp.du.x * rp.dv.x + rp.du.y * rp.dv.y + rp.du.z * rp.dv.z;
		Gp = rp.dv.x * rp.dv.x + rp.dv.y * rp.dv.y + rp.dv.z * rp.dv.z;
		Em = rm.du.x * rm.du.x + rm.du.y * rm.du.y + rm.du.z * rm.du.z;
		Fm = rm.du.x * rm.dv.x + rm.du.y * rm.dv.y + rm.du.z * rm.dv.z;
		Gm = rm.dv.x * rm.dv.x + rm.dv.y * rm.dv.y + rm.dv.z * rm.dv.z;

		dEdv = (Ep - Em) / (QAWS_LITERAL(2.0) * h);
		dFdv = (Fp - Fm) / (QAWS_LITERAL(2.0) * h);
		dGdv = (Gp - Gm) / (QAWS_LITERAL(2.0) * h);
	}

	/* Christoffel symbols of the second kind:
	   Gamma^k_{ij} = 0.5 * g^{km} * (dg_{mi}/dq^j + dg_{mj}/dq^i - dg_{ij}/dq^m)
	   where g = [E F; F G], g^{-1} = inv_det * [G -F; -F E]
	   and q^0 = u, q^1 = v, g_{00}=E, g_{01}=g_{10}=F, g_{11}=G */
	{
		/* dg_{ij}/dq^k stored as dg[i][j][k] */
		qaws_scalar dg[2][2][2];
		qaws_scalar ginv[2][2];

		dg[0][0][0] = dEdu; dg[0][0][1] = dEdv;
		dg[0][1][0] = dFdu; dg[0][1][1] = dFdv;
		dg[1][0][0] = dFdu; dg[1][0][1] = dFdv;
		dg[1][1][0] = dGdu; dg[1][1][1] = dGdv;

		ginv[0][0] = G * inv_det;
		ginv[0][1] = -F * inv_det;
		ginv[1][0] = -F * inv_det;
		ginv[1][1] = E * inv_det;

		for (k = 0; k < 2; k++)
			for (i = 0; i < 2; i++)
				for (j = 0; j < 2; j++)
				{
					int m;
					qaws_scalar val = QAWS_ZERO;
					for (m = 0; m < 2; m++)
						val += ginv[k][m] *
							(dg[m][i][j] + dg[m][j][i] - dg[i][j][m]);
					gamma[k][i][j] = val * QAWS_LITERAL(0.5);
				}
	}
}

qaws_status qaws_surface_compute_geodesic(
	qaws_surface const* surface,
	qaws_scalar start_u,
	qaws_scalar start_v,
	qaws_scalar end_u,
	qaws_scalar end_v,
	unsigned int max_iterations,
	unsigned int step_count,
	qaws_curve** out_curve)
{
	unsigned int max_iter;
	unsigned int n_steps;
	qaws_scalar* path_coords = NULL;
	qaws_scalar* path_params = NULL;
	unsigned int iter;
	qaws_scalar du_target, dv_target, dist_target;
	qaws_scalar best_angle, best_err;
	qaws_scalar angle_lo, angle_hi;
	qaws_status status;
	unsigned int n_cp;

	if (!surface || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;
	max_iter = max_iterations > 0 ? max_iterations : 20;
	n_steps = step_count > 0 ? step_count : 200;

	du_target = end_u - start_u;
	dv_target = end_v - start_v;
	dist_target = QAWS_SQRT(du_target * du_target + dv_target * dv_target);
	if (dist_target < QAWS_LITERAL(1e-10))
		return QAWS_STATUS_DEGENERATE_CURVE;

	path_coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(n_steps + 1) * 3);
	path_params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(n_steps + 1));
	if (!path_coords || !path_params)
	{
		free(path_coords); free(path_params);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Initial shooting angle: direction from start to end in parameter space */
	best_angle = QAWS_ATAN2(dv_target, du_target);
	best_err = QAWS_LITERAL(1e30);
	angle_lo = best_angle - QAWS_LITERAL(1.0);
	angle_hi = best_angle + QAWS_LITERAL(1.0);

	/* Shooting method with bisection on the initial angle */
	for (iter = 0; iter < max_iter; iter++)
	{
		qaws_scalar angle = (angle_lo + angle_hi) * QAWS_LITERAL(0.5);
		qaws_scalar speed = dist_target; /* initial speed in param space */
		qaws_scalar dt = QAWS_ONE / (qaws_scalar)n_steps;
		qaws_scalar u = start_u, v = start_v;
		qaws_scalar udot = speed * QAWS_COS(angle);
		qaws_scalar vdot = speed * QAWS_SIN(angle);
		unsigned int step;
		qaws_scalar final_du, final_dv, err;
		qaws_scalar cross;

		/* Integrate geodesic ODE with RK4:
		   u'' = -Gamma^0_{ij} u'^i u'^j
		   v'' = -Gamma^1_{ij} u'^i u'^j */
		for (step = 0; step <= n_steps; step++)
		{
			/* Store current position */
			if (step == n_steps || iter == max_iter - 1)
			{
				qaws_surface_eval_result sr;
				memset(&sr, 0, sizeof(sr));
				qaws_surface_evaluate(surface, qaws_clamp(u, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(v, QAWS_ZERO, QAWS_ONE),
					QAWS_SURFACE_EVAL_POSITION, &sr);
				path_coords[step * 3 + 0] = sr.position.x;
				path_coords[step * 3 + 1] = sr.position.y;
				path_coords[step * 3 + 2] = sr.position.z;
				path_params[step] = (qaws_scalar)step / (qaws_scalar)n_steps;
			}

			if (step == n_steps) break;

			/* RK4 step */
			{
				qaws_scalar gamma[2][2][2];
				qaws_scalar k1u, k1v, k1ud, k1vd;
				qaws_scalar k2u, k2v, k2ud, k2vd;
				qaws_scalar k3u, k3v, k3ud, k3vd;
				qaws_scalar k4u, k4v, k4ud, k4vd;
				qaws_scalar tu, tv, tud, tvd;
				qaws_scalar acc_u, acc_v;

				/* k1 */
				compute_christoffel(surface,
					qaws_clamp(u, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(v, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * udot * udot +
					QAWS_LITERAL(2.0) * gamma[0][0][1] * udot * vdot +
					gamma[0][1][1] * vdot * vdot);
				acc_v = -(gamma[1][0][0] * udot * udot +
					QAWS_LITERAL(2.0) * gamma[1][0][1] * udot * vdot +
					gamma[1][1][1] * vdot * vdot);
				k1u = udot * dt;
				k1v = vdot * dt;
				k1ud = acc_u * dt;
				k1vd = acc_v * dt;

				/* k2 */
				tu = u + k1u * QAWS_LITERAL(0.5);
				tv = v + k1v * QAWS_LITERAL(0.5);
				tud = udot + k1ud * QAWS_LITERAL(0.5);
				tvd = vdot + k1vd * QAWS_LITERAL(0.5);
				compute_christoffel(surface,
					qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * tud * tud +
					QAWS_LITERAL(2.0) * gamma[0][0][1] * tud * tvd +
					gamma[0][1][1] * tvd * tvd);
				acc_v = -(gamma[1][0][0] * tud * tud +
					QAWS_LITERAL(2.0) * gamma[1][0][1] * tud * tvd +
					gamma[1][1][1] * tvd * tvd);
				k2u = tud * dt;
				k2v = tvd * dt;
				k2ud = acc_u * dt;
				k2vd = acc_v * dt;

				/* k3 */
				tu = u + k2u * QAWS_LITERAL(0.5);
				tv = v + k2v * QAWS_LITERAL(0.5);
				tud = udot + k2ud * QAWS_LITERAL(0.5);
				tvd = vdot + k2vd * QAWS_LITERAL(0.5);
				compute_christoffel(surface,
					qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * tud * tud +
					QAWS_LITERAL(2.0) * gamma[0][0][1] * tud * tvd +
					gamma[0][1][1] * tvd * tvd);
				acc_v = -(gamma[1][0][0] * tud * tud +
					QAWS_LITERAL(2.0) * gamma[1][0][1] * tud * tvd +
					gamma[1][1][1] * tvd * tvd);
				k3u = tud * dt;
				k3v = tvd * dt;
				k3ud = acc_u * dt;
				k3vd = acc_v * dt;

				/* k4 */
				tu = u + k3u;
				tv = v + k3v;
				tud = udot + k3ud;
				tvd = vdot + k3vd;
				compute_christoffel(surface,
					qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * tud * tud +
					QAWS_LITERAL(2.0) * gamma[0][0][1] * tud * tvd +
					gamma[0][1][1] * tvd * tvd);
				acc_v = -(gamma[1][0][0] * tud * tud +
					QAWS_LITERAL(2.0) * gamma[1][0][1] * tud * tvd +
					gamma[1][1][1] * tvd * tvd);
				k4u = tud * dt;
				k4v = tvd * dt;
				k4ud = acc_u * dt;
				k4vd = acc_v * dt;

				/* Update state */
				u += (k1u + QAWS_LITERAL(2.0) * k2u + QAWS_LITERAL(2.0) * k3u + k4u) / QAWS_LITERAL(6.0);
				v += (k1v + QAWS_LITERAL(2.0) * k2v + QAWS_LITERAL(2.0) * k3v + k4v) / QAWS_LITERAL(6.0);
				udot += (k1ud + QAWS_LITERAL(2.0) * k2ud + QAWS_LITERAL(2.0) * k3ud + k4ud) / QAWS_LITERAL(6.0);
				vdot += (k1vd + QAWS_LITERAL(2.0) * k2vd + QAWS_LITERAL(2.0) * k3vd + k4vd) / QAWS_LITERAL(6.0);
			}
		}

		/* Check how close we got to the endpoint */
		final_du = u - end_u;
		final_dv = v - end_v;
		err = QAWS_SQRT(final_du * final_du + final_dv * final_dv);

		if (err < best_err)
		{
			best_err = err;
			best_angle = angle;
		}

		/* Bisect: use cross product to determine which side of target we landed */
		cross = du_target * final_dv - dv_target * final_du;
		if (cross > QAWS_ZERO)
			angle_hi = angle;
		else
			angle_lo = angle;

		if (err < dist_target * QAWS_LITERAL(0.01))
			break;
	}

	/* Do one final integration at the best angle to get the path */
	{
		qaws_scalar dt = QAWS_ONE / (qaws_scalar)n_steps;
		qaws_scalar u = start_u, v = start_v;
		qaws_scalar udot = dist_target * QAWS_COS(best_angle);
		qaws_scalar vdot = dist_target * QAWS_SIN(best_angle);
		unsigned int step;

		for (step = 0; step <= n_steps; step++)
		{
			qaws_surface_eval_result sr;
			memset(&sr, 0, sizeof(sr));
			qaws_surface_evaluate(surface, qaws_clamp(u, QAWS_ZERO, QAWS_ONE),
				qaws_clamp(v, QAWS_ZERO, QAWS_ONE),
				QAWS_SURFACE_EVAL_POSITION, &sr);
			path_coords[step * 3 + 0] = sr.position.x;
			path_coords[step * 3 + 1] = sr.position.y;
			path_coords[step * 3 + 2] = sr.position.z;
			path_params[step] = (qaws_scalar)step / (qaws_scalar)n_steps;

			if (step == n_steps) break;

			{
				qaws_scalar gamma[2][2][2];
				qaws_scalar acc_u, acc_v;
				qaws_scalar k1u, k1v, k1ud, k1vd;
				qaws_scalar k2u, k2v, k2ud, k2vd;
				qaws_scalar k3u, k3v, k3ud, k3vd;
				qaws_scalar k4u, k4v, k4ud, k4vd;
				qaws_scalar tu, tv, tud, tvd;

				compute_christoffel(surface,
					qaws_clamp(u, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(v, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * udot * udot +
					QAWS_LITERAL(2.0) * gamma[0][0][1] * udot * vdot +
					gamma[0][1][1] * vdot * vdot);
				acc_v = -(gamma[1][0][0] * udot * udot +
					QAWS_LITERAL(2.0) * gamma[1][0][1] * udot * vdot +
					gamma[1][1][1] * vdot * vdot);
				k1u = udot * dt; k1v = vdot * dt;
				k1ud = acc_u * dt; k1vd = acc_v * dt;

				tu = u + k1u * QAWS_LITERAL(0.5); tv = v + k1v * QAWS_LITERAL(0.5);
				tud = udot + k1ud * QAWS_LITERAL(0.5); tvd = vdot + k1vd * QAWS_LITERAL(0.5);
				compute_christoffel(surface, qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * tud * tud + QAWS_LITERAL(2.0) * gamma[0][0][1] * tud * tvd + gamma[0][1][1] * tvd * tvd);
				acc_v = -(gamma[1][0][0] * tud * tud + QAWS_LITERAL(2.0) * gamma[1][0][1] * tud * tvd + gamma[1][1][1] * tvd * tvd);
				k2u = tud * dt; k2v = tvd * dt; k2ud = acc_u * dt; k2vd = acc_v * dt;

				tu = u + k2u * QAWS_LITERAL(0.5); tv = v + k2v * QAWS_LITERAL(0.5);
				tud = udot + k2ud * QAWS_LITERAL(0.5); tvd = vdot + k2vd * QAWS_LITERAL(0.5);
				compute_christoffel(surface, qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * tud * tud + QAWS_LITERAL(2.0) * gamma[0][0][1] * tud * tvd + gamma[0][1][1] * tvd * tvd);
				acc_v = -(gamma[1][0][0] * tud * tud + QAWS_LITERAL(2.0) * gamma[1][0][1] * tud * tvd + gamma[1][1][1] * tvd * tvd);
				k3u = tud * dt; k3v = tvd * dt; k3ud = acc_u * dt; k3vd = acc_v * dt;

				tu = u + k3u; tv = v + k3v;
				tud = udot + k3ud; tvd = vdot + k3vd;
				compute_christoffel(surface, qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
					qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
				acc_u = -(gamma[0][0][0] * tud * tud + QAWS_LITERAL(2.0) * gamma[0][0][1] * tud * tvd + gamma[0][1][1] * tvd * tvd);
				acc_v = -(gamma[1][0][0] * tud * tud + QAWS_LITERAL(2.0) * gamma[1][0][1] * tud * tvd + gamma[1][1][1] * tvd * tvd);
				k4u = tud * dt; k4v = tvd * dt; k4ud = acc_u * dt; k4vd = acc_v * dt;

				u += (k1u + QAWS_LITERAL(2.0) * k2u + QAWS_LITERAL(2.0) * k3u + k4u) / QAWS_LITERAL(6.0);
				v += (k1v + QAWS_LITERAL(2.0) * k2v + QAWS_LITERAL(2.0) * k3v + k4v) / QAWS_LITERAL(6.0);
				udot += (k1ud + QAWS_LITERAL(2.0) * k2ud + QAWS_LITERAL(2.0) * k3ud + k4ud) / QAWS_LITERAL(6.0);
				vdot += (k1vd + QAWS_LITERAL(2.0) * k2vd + QAWS_LITERAL(2.0) * k3vd + k4vd) / QAWS_LITERAL(6.0);
			}
		}
	}

	/* Fit a 3D B-spline through the path */
	n_cp = (n_steps + 1) / 4;
	if (n_cp < 6) n_cp = 6;
	if (n_cp > n_steps + 1) n_cp = n_steps + 1;

	status = qaws_internal_fit_bspline(
		QAWS_DIMENSION_3D, 3,
		path_params, path_coords,
		n_steps + 1, n_cp,
		out_curve);

	free(path_coords);
	free(path_params);
	return status;
}
