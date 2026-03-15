#include "qaws_inspect.h"
#include "qaws_eval.h"
#include "qaws_surface.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_solver.h"
#include "internal/qaws_internal_fit.h"
#include "internal/qaws_internal_mesh.h"
#include "internal/qaws_internal_sparse.h"
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

		if (subdivide && stack_top + 4 > stack_capacity)
		{
			unsigned int new_cap = stack_capacity * 2;
			qaws_tess_quad* new_stack = (qaws_tess_quad*)realloc(
				stack, new_cap * sizeof(qaws_tess_quad));
			if (!new_stack)
			{
				free(stack); free(leaves);
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
				qaws_tess_leaf* new_leaves = (qaws_tess_leaf*)realloc(
					leaves, new_cap * sizeof(qaws_tess_leaf));
				if (!new_leaves)
				{
					free(stack); free(leaves);
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

/* Fit a cubic B-spline through path points with exact endpoint
   interpolation.  After least-squares fitting, overwrites the first
   and last control points to guarantee the curve passes through them
   (clamped knot vector ensures C(0) = CP[0], C(1) = CP[n-1]). */
static qaws_status create_path_bspline(
	qaws_scalar const* path_coords,
	unsigned int point_count,
	unsigned int max_cp,
	qaws_curve** out_curve)
{
	unsigned int n_cp, i;
	qaws_status status;
	qaws_bspline_impl* impl;
	unsigned int last_cp;
	qaws_scalar* params = NULL;

	if (point_count < 2)
		return QAWS_STATUS_DEGENERATE_CURVE;

	n_cp = point_count / 4;
	if (n_cp < 6) n_cp = 6;
	if (n_cp > max_cp) n_cp = max_cp;
	if (n_cp > point_count) n_cp = point_count;
	/* B-spline of degree 3 needs at least 4 control points */
	if (n_cp < 4) n_cp = 4;
	if (n_cp > point_count) n_cp = point_count;
	if (point_count < 4)
	{
		/* Too few points for degree-3 B-spline; duplicate interior points */
		/* (This path is a fallback and shouldn't happen in normal operation) */
		return QAWS_STATUS_DEGENERATE_CURVE;
	}

	/* Generate uniform parameters [0, 1] */
	params = (qaws_scalar*)malloc(point_count * sizeof(qaws_scalar));
	if (!params) return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < point_count; i++)
		params[i] = (qaws_scalar)i / (qaws_scalar)(point_count - 1);

	status = qaws_internal_fit_bspline(
		QAWS_DIMENSION_3D, 3,
		params, path_coords,
		point_count, n_cp,
		out_curve);
	free(params);

	if (status != QAWS_STATUS_OK)
		return status;

	/* Force exact endpoint interpolation by overwriting first/last CPs.
	   With clamped knot vector, C(t_min) = CP[0] and C(t_max) = CP[n-1]. */
	impl = (qaws_bspline_impl*)(*out_curve)->impl;
	last_cp = impl->control_point_count - 1;

	impl->control_points[0] = path_coords[0];
	impl->control_points[1] = path_coords[1];
	impl->control_points[2] = path_coords[2];

	impl->control_points[last_cp * 3 + 0] = path_coords[(point_count - 1) * 3 + 0];
	impl->control_points[last_cp * 3 + 1] = path_coords[(point_count - 1) * 3 + 1];
	impl->control_points[last_cp * 3 + 2] = path_coords[(point_count - 1) * 3 + 2];

	return QAWS_STATUS_OK;
}

/* Single RK4 step for geodesic ODE integration.
   Updates (u, v, udot, vdot) in place. */
static void geodesic_rk4_step(
	qaws_surface const* surface,
	qaws_scalar* u, qaws_scalar* v,
	qaws_scalar* udot, qaws_scalar* vdot,
	qaws_scalar dt)
{
	qaws_scalar gamma[2][2][2];
	qaws_scalar k1u, k1v, k1ud, k1vd;
	qaws_scalar k2u, k2v, k2ud, k2vd;
	qaws_scalar k3u, k3v, k3ud, k3vd;
	qaws_scalar k4u, k4v, k4ud, k4vd;
	qaws_scalar tu, tv, tud, tvd, au, av;

	compute_christoffel(surface,
		qaws_clamp(*u, QAWS_ZERO, QAWS_ONE),
		qaws_clamp(*v, QAWS_ZERO, QAWS_ONE), gamma);
	au = -(gamma[0][0][0] * *udot * *udot +
		QAWS_LITERAL(2.0) * gamma[0][0][1] * *udot * *vdot +
		gamma[0][1][1] * *vdot * *vdot);
	av = -(gamma[1][0][0] * *udot * *udot +
		QAWS_LITERAL(2.0) * gamma[1][0][1] * *udot * *vdot +
		gamma[1][1][1] * *vdot * *vdot);
	k1u = *udot * dt; k1v = *vdot * dt;
	k1ud = au * dt; k1vd = av * dt;

	tu = *u + k1u * QAWS_LITERAL(0.5);
	tv = *v + k1v * QAWS_LITERAL(0.5);
	tud = *udot + k1ud * QAWS_LITERAL(0.5);
	tvd = *vdot + k1vd * QAWS_LITERAL(0.5);
	compute_christoffel(surface,
		qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
		qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
	au = -(gamma[0][0][0]*tud*tud + QAWS_LITERAL(2.0)*gamma[0][0][1]*tud*tvd + gamma[0][1][1]*tvd*tvd);
	av = -(gamma[1][0][0]*tud*tud + QAWS_LITERAL(2.0)*gamma[1][0][1]*tud*tvd + gamma[1][1][1]*tvd*tvd);
	k2u = tud * dt; k2v = tvd * dt;
	k2ud = au * dt; k2vd = av * dt;

	tu = *u + k2u * QAWS_LITERAL(0.5);
	tv = *v + k2v * QAWS_LITERAL(0.5);
	tud = *udot + k2ud * QAWS_LITERAL(0.5);
	tvd = *vdot + k2vd * QAWS_LITERAL(0.5);
	compute_christoffel(surface,
		qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
		qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
	au = -(gamma[0][0][0]*tud*tud + QAWS_LITERAL(2.0)*gamma[0][0][1]*tud*tvd + gamma[0][1][1]*tvd*tvd);
	av = -(gamma[1][0][0]*tud*tud + QAWS_LITERAL(2.0)*gamma[1][0][1]*tud*tvd + gamma[1][1][1]*tvd*tvd);
	k3u = tud * dt; k3v = tvd * dt;
	k3ud = au * dt; k3vd = av * dt;

	tu = *u + k3u; tv = *v + k3v;
	tud = *udot + k3ud; tvd = *vdot + k3vd;
	compute_christoffel(surface,
		qaws_clamp(tu, QAWS_ZERO, QAWS_ONE),
		qaws_clamp(tv, QAWS_ZERO, QAWS_ONE), gamma);
	au = -(gamma[0][0][0]*tud*tud + QAWS_LITERAL(2.0)*gamma[0][0][1]*tud*tvd + gamma[0][1][1]*tvd*tvd);
	av = -(gamma[1][0][0]*tud*tud + QAWS_LITERAL(2.0)*gamma[1][0][1]*tud*tvd + gamma[1][1][1]*tvd*tvd);
	k4u = tud * dt; k4v = tvd * dt;
	k4ud = au * dt; k4vd = av * dt;

	*u += (k1u + QAWS_LITERAL(2.0)*k2u + QAWS_LITERAL(2.0)*k3u + k4u) / QAWS_LITERAL(6.0);
	*v += (k1v + QAWS_LITERAL(2.0)*k2v + QAWS_LITERAL(2.0)*k3v + k4v) / QAWS_LITERAL(6.0);
	*udot += (k1ud + QAWS_LITERAL(2.0)*k2ud + QAWS_LITERAL(2.0)*k3ud + k4ud) / QAWS_LITERAL(6.0);
	*vdot += (k1vd + QAWS_LITERAL(2.0)*k2vd + QAWS_LITERAL(2.0)*k3vd + k4vd) / QAWS_LITERAL(6.0);
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
	qaws_scalar du_target, dv_target, dist_target;
	qaws_scalar best_angle, best_err;
	qaws_status status;

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

	path_coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(n_steps + 2) * 3);
	if (!path_coords)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Shooting: grid search ±60° from direct angle, then golden-section
	   refinement.  Score = UV distance at closest approach to target.
	   After finding best angle, integrate final path and trim. */
	best_angle = QAWS_ATAN2(dv_target, du_target);
	best_err = QAWS_LITERAL(1e30);

	/* Phase 1: grid search over ±90° from direct angle */
	{
		unsigned int gi;
		unsigned int n_grid = 64;
		qaws_scalar base = QAWS_ATAN2(dv_target, du_target);
		qaws_scalar dt_g = QAWS_ONE / (qaws_scalar)n_steps;
		qaws_scalar half_range = QAWS_LITERAL(1.5708); /* π/2 */

		for (gi = 0; gi <= n_grid; gi++)
		{
			qaws_scalar angle = base - half_range +
				QAWS_LITERAL(2.0) * half_range * (qaws_scalar)gi / (qaws_scalar)n_grid;
			qaws_scalar u = start_u, v = start_v;
			qaws_scalar udot = dist_target * QAWS_COS(angle);
			qaws_scalar vdot = dist_target * QAWS_SIN(angle);
			unsigned int step;
			qaws_scalar min_d2 = QAWS_LITERAL(1e30);

			for (step = 0; step < n_steps; step++)
			{
				qaws_scalar du2, dv2, d2;
				geodesic_rk4_step(surface, &u, &v, &udot, &vdot, dt_g);
				du2 = u - end_u;
				dv2 = v - end_v;
				d2 = du2 * du2 + dv2 * dv2;
				if (d2 < min_d2) min_d2 = d2;
				if (u < -QAWS_LITERAL(0.05) || u > QAWS_LITERAL(1.05) ||
					v < -QAWS_LITERAL(0.05) || v > QAWS_LITERAL(1.05))
					break;
			}

			if (min_d2 < best_err * best_err)
			{
				best_err = QAWS_SQRT(min_d2);
				best_angle = angle;
			}
		}
	}

	/* Phase 2: golden-section refinement around best angle */
	{
		qaws_scalar a = best_angle - QAWS_LITERAL(0.15);
		qaws_scalar b = best_angle + QAWS_LITERAL(0.15);
		qaws_scalar gr = QAWS_LITERAL(0.6180339887);
		unsigned int ri;

		for (ri = 0; ri < max_iter; ri++)
		{
			qaws_scalar c = b - gr * (b - a);
			qaws_scalar d = a + gr * (b - a);
			qaws_scalar fc, fd;

			/* Evaluate angle c */
			{
				qaws_scalar u = start_u, v = start_v;
				qaws_scalar udot = dist_target * QAWS_COS(c);
				qaws_scalar vdot = dist_target * QAWS_SIN(c);
				qaws_scalar dt_g = QAWS_ONE / (qaws_scalar)n_steps;
				unsigned int step;
				fc = QAWS_LITERAL(1e30);
				for (step = 0; step < n_steps; step++)
				{
					qaws_scalar du2, dv2, d2;
					geodesic_rk4_step(surface, &u, &v, &udot, &vdot, dt_g);
					du2 = u - end_u; dv2 = v - end_v;
					d2 = du2 * du2 + dv2 * dv2;
					if (d2 < fc) fc = d2;
					if (u < -QAWS_LITERAL(0.05) || u > QAWS_LITERAL(1.05) ||
						v < -QAWS_LITERAL(0.05) || v > QAWS_LITERAL(1.05))
						break;
				}
			}

			/* Evaluate angle d */
			{
				qaws_scalar u = start_u, v = start_v;
				qaws_scalar udot = dist_target * QAWS_COS(d);
				qaws_scalar vdot = dist_target * QAWS_SIN(d);
				qaws_scalar dt_g = QAWS_ONE / (qaws_scalar)n_steps;
				unsigned int step;
				fd = QAWS_LITERAL(1e30);
				for (step = 0; step < n_steps; step++)
				{
					qaws_scalar du2, dv2, d2;
					geodesic_rk4_step(surface, &u, &v, &udot, &vdot, dt_g);
					du2 = u - end_u; dv2 = v - end_v;
					d2 = du2 * du2 + dv2 * dv2;
					if (d2 < fd) fd = d2;
					if (u < -QAWS_LITERAL(0.05) || u > QAWS_LITERAL(1.05) ||
						v < -QAWS_LITERAL(0.05) || v > QAWS_LITERAL(1.05))
						break;
				}
			}

			if (fc < fd)
				b = d;
			else
				a = c;
		}

		best_angle = (a + b) * QAWS_LITERAL(0.5);
	}

	/* Final integration at best angle.  Store both UV and 3D positions,
	   then trim at closest approach in UV to the target. */
	{
		qaws_scalar dt = QAWS_ONE / (qaws_scalar)n_steps;
		qaws_scalar u = start_u, v = start_v;
		qaws_scalar udot = dist_target * QAWS_COS(best_angle);
		qaws_scalar vdot = dist_target * QAWS_SIN(best_angle);
		unsigned int step;
		unsigned int actual_count = 0;
		unsigned int best_pi = 0;
		qaws_scalar best_d2 = QAWS_LITERAL(1e30);

		for (step = 0; step <= n_steps; step++)
		{
			qaws_surface_eval_result sr;
			qaws_scalar du2, dv2, d2;

			memset(&sr, 0, sizeof(sr));
			qaws_surface_evaluate(surface, qaws_clamp(u, QAWS_ZERO, QAWS_ONE),
				qaws_clamp(v, QAWS_ZERO, QAWS_ONE),
				QAWS_SURFACE_EVAL_POSITION, &sr);
			path_coords[actual_count * 3 + 0] = sr.position.x;
			path_coords[actual_count * 3 + 1] = sr.position.y;
			path_coords[actual_count * 3 + 2] = sr.position.z;

			/* Track closest approach in UV */
			du2 = u - end_u; dv2 = v - end_v;
			d2 = du2 * du2 + dv2 * dv2;
			if (d2 < best_d2) { best_d2 = d2; best_pi = actual_count; }

			actual_count++;
			if (step == n_steps) break;

			geodesic_rk4_step(surface, &u, &v, &udot, &vdot, dt);

			if (u < -QAWS_LITERAL(0.05) || u > QAWS_LITERAL(1.05) ||
				v < -QAWS_LITERAL(0.05) || v > QAWS_LITERAL(1.05))
			{
				u = qaws_clamp(u, QAWS_ZERO, QAWS_ONE);
				v = qaws_clamp(v, QAWS_ZERO, QAWS_ONE);
				{
					qaws_surface_eval_result sr2;
					memset(&sr2, 0, sizeof(sr2));
					qaws_surface_evaluate(surface, u, v,
						QAWS_SURFACE_EVAL_POSITION, &sr2);
					path_coords[actual_count * 3 + 0] = sr2.position.x;
					path_coords[actual_count * 3 + 1] = sr2.position.y;
					path_coords[actual_count * 3 + 2] = sr2.position.z;
					du2 = u - end_u; dv2 = v - end_v;
					d2 = du2 * du2 + dv2 * dv2;
					if (d2 < best_d2) { best_d2 = d2; best_pi = actual_count; }
					actual_count++;
				}
				break;
			}
		}

		/* Trim at closest approach to target */
		if (best_pi + 1 < actual_count)
			actual_count = best_pi + 1;

		/* Append exact target position */
		{
			qaws_surface_eval_result end_sr;
			memset(&end_sr, 0, sizeof(end_sr));
			qaws_surface_evaluate(surface, end_u, end_v,
				QAWS_SURFACE_EVAL_POSITION, &end_sr);
			path_coords[actual_count * 3 + 0] = end_sr.position.x;
			path_coords[actual_count * 3 + 1] = end_sr.position.y;
			path_coords[actual_count * 3 + 2] = end_sr.position.z;
			actual_count++;
		}

		/* Create clamped B-spline (exact endpoints) */
		status = create_path_bspline(
			path_coords, actual_count, 150, out_curve);
	}

	free(path_coords);
	return status;
}

/* ------------------------------------------------------------------ */
/*  Geodesic: tessellate surface and build mesh helper                 */
/* ------------------------------------------------------------------ */

static qaws_status geodesic_tessellate_and_build(
	qaws_surface const* surface,
	qaws_tessellation_desc const* tessellation_desc,
	qaws_tessellation_vertex** out_verts,
	unsigned int** out_indices,
	unsigned int* out_vert_count,
	unsigned int* out_idx_count,
	qaws_scalar** out_positions,
	qaws_scalar** out_uvs,
	qaws_internal_mesh* out_mesh)
{
	/* Use a regular grid tessellation for PDE quality (uniform triangles,
	   no degenerate elements from adaptive fan triangulation).
	   Grid resolution derived from tessellation_desc or defaults. */
	unsigned int grid_n;
	unsigned int vert_count, tri_count, idx_count;
	qaws_tessellation_vertex* verts = NULL;
	unsigned int* indices = NULL;
	qaws_scalar* positions = NULL;
	qaws_scalar* uvs = NULL;
	qaws_status status;
	unsigned int i, j;
	qaws_range u_range, v_range;
	qaws_scalar u_min, u_max, v_min, v_max;

	u_range = qaws_surface_get_u_range(surface);
	v_range = qaws_surface_get_v_range(surface);
	u_min = u_range.min_value; u_max = u_range.max_value;
	v_min = v_range.min_value; v_max = v_range.max_value;

	/* Determine grid resolution */
	if (tessellation_desc && tessellation_desc->max_depth > 0)
		grid_n = 1u << tessellation_desc->max_depth; /* 2^depth */
	else
		grid_n = 32; /* default 32x32 grid */

	/* Cap to avoid excessive memory */
	if (grid_n > 256) grid_n = 256;

	vert_count = (grid_n + 1) * (grid_n + 1);
	tri_count = grid_n * grid_n * 2;
	idx_count = tri_count * 3;

	verts = (qaws_tessellation_vertex*)malloc(vert_count * sizeof(qaws_tessellation_vertex));
	indices = (unsigned int*)malloc(idx_count * sizeof(unsigned int));
	positions = (qaws_scalar*)malloc(vert_count * 3 * sizeof(qaws_scalar));
	uvs = (qaws_scalar*)malloc(vert_count * 2 * sizeof(qaws_scalar));

	if (!verts || !indices || !positions || !uvs)
	{
		free(verts); free(indices); free(positions); free(uvs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Generate grid vertices */
	for (i = 0; i <= grid_n; i++)
	{
		qaws_scalar v = v_min + (v_max - v_min) * (qaws_scalar)i / (qaws_scalar)grid_n;
		for (j = 0; j <= grid_n; j++)
		{
			qaws_scalar u = u_min + (u_max - u_min) * (qaws_scalar)j / (qaws_scalar)grid_n;
			qaws_surface_eval_result sr;
			unsigned int idx = i * (grid_n + 1) + j;

			memset(&sr, 0, sizeof(sr));
			status = qaws_surface_evaluate(surface, u, v,
				QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &sr);
			if (status != QAWS_STATUS_OK)
			{
				free(verts); free(indices); free(positions); free(uvs);
				return status;
			}

			verts[idx].position = sr.position;
			verts[idx].normal = sr.normal;
			verts[idx].u = u;
			verts[idx].v = v;

			positions[idx * 3 + 0] = sr.position.x;
			positions[idx * 3 + 1] = sr.position.y;
			positions[idx * 3 + 2] = sr.position.z;
			uvs[idx * 2 + 0] = u;
			uvs[idx * 2 + 1] = v;
		}
	}

	/* Generate triangle indices (2 triangles per grid cell) */
	{
		unsigned int ii = 0;
		for (i = 0; i < grid_n; i++)
		{
			for (j = 0; j < grid_n; j++)
			{
				unsigned int v00 = i * (grid_n + 1) + j;
				unsigned int v10 = v00 + 1;
				unsigned int v01 = v00 + (grid_n + 1);
				unsigned int v11 = v01 + 1;

				/* Lower-left triangle */
				indices[ii++] = v00;
				indices[ii++] = v10;
				indices[ii++] = v01;

				/* Upper-right triangle */
				indices[ii++] = v10;
				indices[ii++] = v11;
				indices[ii++] = v01;
			}
		}
	}

	/* Build edge-adjacency mesh */
	status = qaws_internal_mesh_build(out_mesh, positions, uvs,
		vert_count, indices, idx_count);

	if (status != QAWS_STATUS_OK)
	{
		free(verts); free(indices); free(positions); free(uvs);
		return status;
	}

	*out_verts = verts;
	*out_indices = indices;
	*out_vert_count = vert_count;
	*out_idx_count = idx_count;
	*out_positions = positions;
	*out_uvs = uvs;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Geodesic via edge flipping (Sharp & Crane SGP 2020)                */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_compute_geodesic_flip(
	qaws_surface const* surface,
	qaws_scalar start_u, qaws_scalar start_v,
	qaws_scalar end_u, qaws_scalar end_v,
	qaws_tessellation_desc const* tessellation_desc,
	unsigned int max_iterations,
	qaws_curve** out_curve)
{
	qaws_tessellation_vertex* verts = NULL;
	unsigned int* indices = NULL;
	unsigned int vert_count = 0, idx_count = 0;
	qaws_scalar* positions = NULL;
	qaws_scalar* uvs = NULL;
	qaws_internal_mesh mesh;
	unsigned int src_v, dst_v;
	unsigned int* path = NULL;
	unsigned int path_count = 0;
	unsigned int max_iter;
	qaws_status status;

	if (!surface || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;
	max_iter = max_iterations > 0 ? max_iterations : 1000;

	memset(&mesh, 0, sizeof(mesh));

	status = geodesic_tessellate_and_build(surface, tessellation_desc,
		&verts, &indices, &vert_count, &idx_count,
		&positions, &uvs, &mesh);

	if (status != QAWS_STATUS_OK)
		return status;

	/* Find closest vertices to start/end UV */
	src_v = qaws_internal_mesh_closest_vertex_uv(&mesh, start_u, start_v);
	dst_v = qaws_internal_mesh_closest_vertex_uv(&mesh, end_u, end_v);

	if (src_v == dst_v)
	{
		qaws_internal_mesh_destroy(&mesh);
		free(verts); free(indices); free(positions); free(uvs);
		return QAWS_STATUS_DEGENERATE_CURVE;
	}

	/* Dijkstra shortest path */
	path = (unsigned int*)malloc(vert_count * sizeof(unsigned int));
	if (!path)
	{
		qaws_internal_mesh_destroy(&mesh);
		free(verts); free(indices); free(positions); free(uvs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	status = qaws_internal_mesh_dijkstra(&mesh, src_v, dst_v, path, &path_count);
	if (status != QAWS_STATUS_OK || path_count < 2)
	{
		free(path);
		qaws_internal_mesh_destroy(&mesh);
		free(verts); free(indices); free(positions); free(uvs);
		return status != QAWS_STATUS_OK ? status : QAWS_STATUS_NUMERICAL_FAILURE;
	}

	/* Iterative shortening via fan unfolding (simplified Sharp & Crane).
	   For each interior vertex, unfold the triangle fan between the
	   incoming and outgoing Dijkstra edges into 2D, then find where the
	   straight line from predecessor to successor crosses edges.
	   Multiple passes until converged. */
	{
		/* Build a refined point list from the Dijkstra vertex path.
		   Each point is a 3D position on the mesh.  Interior points
		   are iteratively shortened by fan unfolding. */
		unsigned int out_count;
		unsigned int pass;
		qaws_scalar* pts = NULL;   /* 3 * out_count */
		qaws_scalar* uv_pts = NULL; /* 2 * out_count */

		out_count = path_count;
		pts = (qaws_scalar*)malloc(out_count * 3 * sizeof(qaws_scalar));
		uv_pts = (qaws_scalar*)malloc(out_count * 2 * sizeof(qaws_scalar));
		if (!pts || !uv_pts)
		{
			free(pts); free(uv_pts);
			free(path);
			qaws_internal_mesh_destroy(&mesh);
			free(verts); free(indices); free(positions); free(uvs);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		{
			unsigned int pi;
			for (pi = 0; pi < out_count; pi++)
			{
				unsigned int vi = path[pi];
				pts[pi * 3 + 0] = positions[vi * 3 + 0];
				pts[pi * 3 + 1] = positions[vi * 3 + 1];
				pts[pi * 3 + 2] = positions[vi * 3 + 2];
				uv_pts[pi * 2 + 0] = uvs[vi * 2 + 0];
				uv_pts[pi * 2 + 1] = uvs[vi * 2 + 1];
			}
		}

		/* For each interior vertex, try to shortcut via fan unfolding.
		   Walk the triangle fan around path[i] from edge(path[i-1],path[i])
		   to edge(path[i],path[i+1]), unfolding into 2D as we go.
		   Then draw the straight line from p_{i-1} to p_{i+1} in the
		   unfolded picture and use the midpoint on that line as the
		   new position for path[i] (evaluated on the surface). */
		for (pass = 0; pass < max_iter; pass++)
		{
			int changed = 0;
			unsigned int pi;
			for (pi = 1; pi + 1 < out_count; pi++)
			{
				/* Find the edge from path[pi-1] to path[pi] */
				unsigned int v_prev = path[pi - 1];
				unsigned int v_cur  = path[pi];
				unsigned int v_next = path[pi + 1];
				unsigned int e_in = UINT_MAX, e_out = UINT_MAX;
				unsigned int ei;
				unsigned int tri_in = UINT_MAX;

				/* Find incoming/outgoing edges */
				for (ei = mesh.vert_edge_offset[v_cur];
					ei < mesh.vert_edge_offset[v_cur + 1]; ei++)
				{
					unsigned int eidx = mesh.vert_edges[ei];
					unsigned int ov = (mesh.edges[eidx].v[0] == v_cur) ?
						mesh.edges[eidx].v[1] : mesh.edges[eidx].v[0];
					if (ov == v_prev) e_in = eidx;
					if (ov == v_next) e_out = eidx;
				}
				if (e_in == UINT_MAX || e_out == UINT_MAX) continue;
				if (e_in == e_out) continue; /* consecutive edges are the same */

				/* Find the triangle on one side of e_in that contains v_cur.
				   We'll walk CW around v_cur toward e_out. */
				tri_in = mesh.edges[e_in].tri[0];
				if (tri_in == UINT_MAX) tri_in = mesh.edges[e_in].tri[1];

				/* Unfold triangle fan from e_in toward e_out.
				   Layout: v_cur at origin.
				   v_prev at (l_in, 0) where l_in = edge length of e_in. */
				{
					qaws_scalar l_in, prev_2d[2];
					qaws_scalar dx, dy, dz;
					unsigned int cur_edge = e_in;
					unsigned int cur_fan_tri = tri_in;
					unsigned int fan_steps = 0;
					unsigned int max_fan = 50;
					/* Accumulate angle around v_cur */
					qaws_scalar accum_angle = QAWS_ZERO;
					qaws_scalar next_2d[2] = {QAWS_ZERO, QAWS_ZERO};
					int found_out = 0;

					dx = positions[v_prev * 3] - positions[v_cur * 3];
					dy = positions[v_prev * 3 + 1] - positions[v_cur * 3 + 1];
					dz = positions[v_prev * 3 + 2] - positions[v_cur * 3 + 2];
					l_in = QAWS_SQRT(dx * dx + dy * dy + dz * dz);
					prev_2d[0] = l_in;
					prev_2d[1] = QAWS_ZERO;

					/* Walk the fan around v_cur.  In each triangle, find the
					   edge that is NOT cur_edge and NOT opposite v_cur, i.e.
					   the "next" edge in the fan. */
					while (fan_steps < max_fan)
					{
						int cur_local = -1; /* local index of v_cur */
						int prev_local = -1; /* local index of v on cur_edge */
						int third_local;
						unsigned int v_third;
						qaws_scalar l_cur_edge, l_third;
						qaws_scalar tri_angle, cos_angle;
						unsigned int next_edge;
						int mi;

						for (mi = 0; mi < 3; mi++)
						{
							if (mesh.tris[cur_fan_tri].v[mi] == v_cur) cur_local = mi;
						}
						if (cur_local < 0) break;

						/* The edge cur_edge connects v_cur and another vertex.
						   Find which local index that other vertex has. */
						{
							unsigned int ov = (mesh.edges[cur_edge].v[0] == v_cur) ?
								mesh.edges[cur_edge].v[1] : mesh.edges[cur_edge].v[0];
							for (mi = 0; mi < 3; mi++)
								if (mesh.tris[cur_fan_tri].v[mi] == ov) prev_local = mi;
						}
						if (prev_local < 0) break;

						/* The third vertex */
						third_local = 3 - cur_local - prev_local;
						v_third = mesh.tris[cur_fan_tri].v[third_local];

						/* Edge lengths from v_cur */
						dx = positions[mesh.tris[cur_fan_tri].v[prev_local] * 3] - positions[v_cur * 3];
						dy = positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 1] - positions[v_cur * 3 + 1];
						dz = positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 2] - positions[v_cur * 3 + 2];
						l_cur_edge = QAWS_SQRT(dx * dx + dy * dy + dz * dz);

						dx = positions[v_third * 3] - positions[v_cur * 3];
						dy = positions[v_third * 3 + 1] - positions[v_cur * 3 + 1];
						dz = positions[v_third * 3 + 2] - positions[v_cur * 3 + 2];
						l_third = QAWS_SQRT(dx * dx + dy * dy + dz * dz);

						/* Angle at v_cur between cur_edge and edge to v_third */
						{
							qaws_scalar l_opp; /* edge opposite v_cur */
							dx = positions[v_third * 3] - positions[mesh.tris[cur_fan_tri].v[prev_local] * 3];
							dy = positions[v_third * 3 + 1] - positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 1];
							dz = positions[v_third * 3 + 2] - positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 2];
							l_opp = QAWS_SQRT(dx * dx + dy * dy + dz * dz);

							if (l_cur_edge < QAWS_LITERAL(1e-15) || l_third < QAWS_LITERAL(1e-15))
								break;
							cos_angle = (l_cur_edge * l_cur_edge + l_third * l_third
								- l_opp * l_opp) / (QAWS_LITERAL(2.0) * l_cur_edge * l_third);
							if (cos_angle > QAWS_ONE) cos_angle = QAWS_ONE;
							if (cos_angle < -QAWS_ONE) cos_angle = -QAWS_ONE;
							tri_angle = QAWS_ATAN2(
								QAWS_SQRT(QAWS_ONE - cos_angle * cos_angle),
								cos_angle);
						}

						/* v_third in 2D: rotated by accum_angle + tri_angle from +x */
						{
							qaws_scalar total_a = accum_angle + tri_angle;
							next_2d[0] = l_third * QAWS_COS(total_a);
							next_2d[1] = l_third * QAWS_SIN(total_a);
						}

						/* Check if e_out connects v_cur and v_third */
						next_edge = mesh.tris[cur_fan_tri].edge[prev_local];
						/* edge[prev_local] is the edge opposite prev_local
						   = edge connecting v_cur and v_third */

						if (v_third == v_next)
						{
							found_out = 1;
							break;
						}

						/* Advance to next triangle in fan */
						accum_angle += tri_angle;
						cur_edge = next_edge;
						{
							unsigned int side = (mesh.edges[next_edge].tri[0] == cur_fan_tri) ? 1u : 0u;
							cur_fan_tri = mesh.edges[next_edge].tri[side];
							if (cur_fan_tri == UINT_MAX) break;
						}
						fan_steps++;
					}

					if (!found_out)
					{
						/* Try the other direction: walk the fan from e_in
						   in the opposite direction around v_cur */
						cur_edge = e_in;
						cur_fan_tri = mesh.edges[e_in].tri[1];
						if (cur_fan_tri == UINT_MAX || cur_fan_tri == tri_in)
							continue;
						accum_angle = QAWS_ZERO;
						found_out = 0;
						fan_steps = 0;

						while (fan_steps < max_fan)
						{
							int cur_local = -1;
							int prev_local = -1;
							int third_local;
							unsigned int v_third;
							qaws_scalar l_cur_edge, l_third;
							qaws_scalar tri_angle, cos_angle;
							unsigned int next_edge;
							int mi;

							for (mi = 0; mi < 3; mi++)
								if (mesh.tris[cur_fan_tri].v[mi] == v_cur) cur_local = mi;
							if (cur_local < 0) break;

							{
								unsigned int ov = (mesh.edges[cur_edge].v[0] == v_cur) ?
									mesh.edges[cur_edge].v[1] : mesh.edges[cur_edge].v[0];
								for (mi = 0; mi < 3; mi++)
									if (mesh.tris[cur_fan_tri].v[mi] == ov) prev_local = mi;
							}
							if (prev_local < 0) break;

							third_local = 3 - cur_local - prev_local;
							v_third = mesh.tris[cur_fan_tri].v[third_local];

							dx = positions[mesh.tris[cur_fan_tri].v[prev_local] * 3] - positions[v_cur * 3];
							dy = positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 1] - positions[v_cur * 3 + 1];
							dz = positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 2] - positions[v_cur * 3 + 2];
							l_cur_edge = QAWS_SQRT(dx * dx + dy * dy + dz * dz);

							dx = positions[v_third * 3] - positions[v_cur * 3];
							dy = positions[v_third * 3 + 1] - positions[v_cur * 3 + 1];
							dz = positions[v_third * 3 + 2] - positions[v_cur * 3 + 2];
							l_third = QAWS_SQRT(dx * dx + dy * dy + dz * dz);

							{
								qaws_scalar l_opp;
								dx = positions[v_third * 3] - positions[mesh.tris[cur_fan_tri].v[prev_local] * 3];
								dy = positions[v_third * 3 + 1] - positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 1];
								dz = positions[v_third * 3 + 2] - positions[mesh.tris[cur_fan_tri].v[prev_local] * 3 + 2];
								l_opp = QAWS_SQRT(dx * dx + dy * dy + dz * dz);
								if (l_cur_edge < QAWS_LITERAL(1e-15) || l_third < QAWS_LITERAL(1e-15))
									break;
								cos_angle = (l_cur_edge * l_cur_edge + l_third * l_third
									- l_opp * l_opp) / (QAWS_LITERAL(2.0) * l_cur_edge * l_third);
								if (cos_angle > QAWS_ONE) cos_angle = QAWS_ONE;
								if (cos_angle < -QAWS_ONE) cos_angle = -QAWS_ONE;
								tri_angle = QAWS_ATAN2(
								QAWS_SQRT(QAWS_ONE - cos_angle * cos_angle),
								cos_angle);
							}

							accum_angle -= tri_angle;
							next_2d[0] = l_third * QAWS_COS(accum_angle);
							next_2d[1] = l_third * QAWS_SIN(accum_angle);

							next_edge = mesh.tris[cur_fan_tri].edge[prev_local];

							if (v_third == v_next)
							{
								found_out = 1;
								break;
							}

							cur_edge = next_edge;
							{
								unsigned int side = (mesh.edges[next_edge].tri[0] == cur_fan_tri) ? 1u : 0u;
								cur_fan_tri = mesh.edges[next_edge].tri[side];
								if (cur_fan_tri == UINT_MAX) break;
							}
							fan_steps++;
						}
					}

					if (found_out)
					{
						/* In the unfolded picture:
						   v_cur at (0,0), v_prev at prev_2d, v_next at next_2d.
						   The geodesic shortcut is the midpoint of the straight
						   line from prev_2d to next_2d.  Map it back to UV by
						   interpolating the UV of the edge it lies on.
						   For simplicity, use the midpoint's distance ratio
						   to update the vertex position on the surface. */
						qaws_scalar mid_2d[2], t_along;
						qaws_scalar new_u, new_v;
						qaws_surface_eval_result sr;

						mid_2d[0] = (prev_2d[0] + next_2d[0]) * QAWS_LITERAL(0.5);
						mid_2d[1] = (prev_2d[1] + next_2d[1]) * QAWS_LITERAL(0.5);

						/* Express midpoint as fraction along v_cur's position:
						   project onto line from v_cur to midpoint.
						   UV = lerp between prev and next based on
						   distance ratio from unfolded picture. */
						{
							qaws_scalar dp = QAWS_SQRT(prev_2d[0] * prev_2d[0] + prev_2d[1] * prev_2d[1]);
							qaws_scalar dn = QAWS_SQRT(next_2d[0] * next_2d[0] + next_2d[1] * next_2d[1]);
							if (dp + dn > QAWS_LITERAL(1e-15))
								t_along = dp / (dp + dn);
							else
								t_along = QAWS_LITERAL(0.5);
						}
						new_u = (QAWS_ONE - t_along) * uv_pts[(pi - 1) * 2]
							+ t_along * uv_pts[(pi + 1) * 2];
						new_v = (QAWS_ONE - t_along) * uv_pts[(pi - 1) * 2 + 1]
							+ t_along * uv_pts[(pi + 1) * 2 + 1];

						memset(&sr, 0, sizeof(sr));
						qaws_surface_evaluate(surface, new_u, new_v,
							QAWS_SURFACE_EVAL_POSITION, &sr);

						/* Check if new position is closer to the line */
						{
							qaws_scalar old_len, new_len;
							dx = pts[pi * 3] - pts[(pi - 1) * 3];
							dy = pts[pi * 3 + 1] - pts[(pi - 1) * 3 + 1];
							dz = pts[pi * 3 + 2] - pts[(pi - 1) * 3 + 2];
							old_len = QAWS_SQRT(dx * dx + dy * dy + dz * dz);
							dx = pts[(pi + 1) * 3] - pts[pi * 3];
							dy = pts[(pi + 1) * 3 + 1] - pts[pi * 3 + 1];
							dz = pts[(pi + 1) * 3 + 2] - pts[pi * 3 + 2];
							old_len += QAWS_SQRT(dx * dx + dy * dy + dz * dz);

							dx = (qaws_scalar)sr.position.x - pts[(pi - 1) * 3];
							dy = (qaws_scalar)sr.position.y - pts[(pi - 1) * 3 + 1];
							dz = (qaws_scalar)sr.position.z - pts[(pi - 1) * 3 + 2];
							new_len = QAWS_SQRT(dx * dx + dy * dy + dz * dz);
							dx = pts[(pi + 1) * 3] - (qaws_scalar)sr.position.x;
							dy = pts[(pi + 1) * 3 + 1] - (qaws_scalar)sr.position.y;
							dz = pts[(pi + 1) * 3 + 2] - (qaws_scalar)sr.position.z;
							new_len += QAWS_SQRT(dx * dx + dy * dy + dz * dz);

							if (new_len < old_len - QAWS_LITERAL(1e-10))
							{
								pts[pi * 3 + 0] = sr.position.x;
								pts[pi * 3 + 1] = sr.position.y;
								pts[pi * 3 + 2] = sr.position.z;
								uv_pts[pi * 2 + 0] = new_u;
								uv_pts[pi * 2 + 1] = new_v;
								changed = 1;
							}
						}
					}
				}
			}
			if (!changed) break;
		}

		/* Set exact start/end positions */
		{
			qaws_surface_eval_result sr;
			memset(&sr, 0, sizeof(sr));
			qaws_surface_evaluate(surface, start_u, start_v,
				QAWS_SURFACE_EVAL_POSITION, &sr);
			pts[0] = sr.position.x;
			pts[1] = sr.position.y;
			pts[2] = sr.position.z;
			memset(&sr, 0, sizeof(sr));
			qaws_surface_evaluate(surface, end_u, end_v,
				QAWS_SURFACE_EVAL_POSITION, &sr);
			pts[(out_count - 1) * 3 + 0] = sr.position.x;
			pts[(out_count - 1) * 3 + 1] = sr.position.y;
			pts[(out_count - 1) * 3 + 2] = sr.position.z;
		}

		status = create_path_bspline(pts, out_count, 150, out_curve);

		free(pts);
		free(uv_pts);
	}

	free(path);
	qaws_internal_mesh_destroy(&mesh);
	free(verts);
	free(indices);
	free(positions);
	free(uvs);

	return status;
}

/* ------------------------------------------------------------------ */
/*  Heat method internal solver                                        */
/* ------------------------------------------------------------------ */

/* Solve heat method on a mesh to compute geodesic distances.
   Returns allocated distance array via out_distances (caller must free). */
static qaws_status heat_solve_distances(
	qaws_internal_mesh const* mesh,
	qaws_scalar const* positions,
	unsigned int vert_count,
	unsigned int src_v,
	qaws_scalar** out_distances)
{
	qaws_internal_csr L, A_minus_tL, L_reg;
	qaws_scalar* mass_diag = NULL;
	qaws_scalar* heat = NULL;
	qaws_scalar* rhs = NULL;
	qaws_scalar* div_X = NULL;
	qaws_scalar mean_edge_len, t_heat;
	unsigned int vi, ti;
	qaws_status status;

	*out_distances = NULL;
	memset(&L, 0, sizeof(L));
	memset(&A_minus_tL, 0, sizeof(A_minus_tL));
	memset(&L_reg, 0, sizeof(L_reg));

	/* Build Laplacian and mass matrix */
	mass_diag = (qaws_scalar*)malloc(vert_count * sizeof(qaws_scalar));
	if (!mass_diag) return QAWS_STATUS_ALLOCATION_FAILURE;

	status = qaws_internal_build_laplacian(mesh, &L, mass_diag);
	if (status != QAWS_STATUS_OK) { free(mass_diag); return status; }

	/* Compute mean edge length */
	{
		qaws_scalar sum = QAWS_ZERO;
		unsigned int ei;
		for (ei = 0; ei < mesh->edge_count; ei++)
			sum += qaws_internal_mesh_edge_length(mesh,
				mesh->edges[ei].v[0], mesh->edges[ei].v[1]);
		mean_edge_len = sum / (qaws_scalar)mesh->edge_count;
	}
	/* t = h^2 is the standard prescription (Crane et al. 2017) for small meshes.
	   For large meshes with many edges between source and farthest vertex,
	   the heat kernel decays as exp(-D^2/(4t)). With t=h^2, this is
	   exp(-(D/h)^2/4) which vanishes when D >> h.
	   Fix: ensure t >= D^2 / (4 * 30) so heat stays above ~exp(-30). */
	t_heat = mean_edge_len * mean_edge_len;
	{
		qaws_scalar max_r_sq = QAWS_ZERO;
		unsigned int vi2;
		for (vi2 = 0; vi2 < vert_count; vi2++)
		{
			qaws_scalar dx = positions[vi2 * 3]     - positions[src_v * 3];
			qaws_scalar dy = positions[vi2 * 3 + 1] - positions[src_v * 3 + 1];
			qaws_scalar dz = positions[vi2 * 3 + 2] - positions[src_v * 3 + 2];
			qaws_scalar r_sq = dx * dx + dy * dy + dz * dz;
			if (r_sq > max_r_sq) max_r_sq = r_sq;
		}
		{
			qaws_scalar t_min = max_r_sq / QAWS_LITERAL(120.0);
			if (t_min > t_heat) t_heat = t_min;
		}
	}

	/* Step I: Build (A - t*L) and solve (A - t*L) * heat = delta_source */
	{
		status = qaws_internal_csr_alloc(&A_minus_tL, vert_count, vert_count, L.nnz);
		if (status != QAWS_STATUS_OK) goto cleanup;

		memcpy(A_minus_tL.row_ptr, L.row_ptr, (vert_count + 1) * sizeof(unsigned int));
		memcpy(A_minus_tL.col_idx, L.col_idx, L.nnz * sizeof(unsigned int));

		for (vi = 0; vi < L.nnz; vi++)
			A_minus_tL.values[vi] = -t_heat * L.values[vi];

		for (vi = 0; vi < vert_count; vi++)
			A_minus_tL.values[A_minus_tL.row_ptr[vi]] += mass_diag[vi];
	}

	heat = (qaws_scalar*)calloc(vert_count, sizeof(qaws_scalar));
	rhs = (qaws_scalar*)calloc(vert_count, sizeof(qaws_scalar));
	if (!heat || !rhs) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	rhs[src_v] = QAWS_ONE;

	status = qaws_internal_cg_solve(&A_minus_tL, rhs, heat, 0, QAWS_ZERO);
	if (status != QAWS_STATUS_OK) goto cleanup;

	/* Step II: Per-triangle gradient X = -grad(u) / |grad(u)| */
	/* Step III: Integrated divergence div(X) at each vertex */
	div_X = (qaws_scalar*)calloc(vert_count, sizeof(qaws_scalar));
	if (!div_X) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	for (ti = 0; ti < mesh->tri_count; ti++)
	{
		unsigned int va = mesh->tris[ti].v[0];
		unsigned int vb = mesh->tris[ti].v[1];
		unsigned int vc = mesh->tris[ti].v[2];
		qaws_scalar const* pa = positions + va * 3;
		qaws_scalar const* pb = positions + vb * 3;
		qaws_scalar const* pc = positions + vc * 3;

		qaws_scalar e0x = pb[0] - pa[0], e0y = pb[1] - pa[1], e0z = pb[2] - pa[2];
		qaws_scalar e1x = pc[0] - pa[0], e1y = pc[1] - pa[1], e1z = pc[2] - pa[2];

		qaws_scalar nx = e0y * e1z - e0z * e1y;
		qaws_scalar ny = e0z * e1x - e0x * e1z;
		qaws_scalar nz = e0x * e1y - e0y * e1x;
		qaws_scalar area2 = QAWS_SQRT(nx * nx + ny * ny + nz * nz);

		if (area2 < QAWS_LITERAL(1e-15)) continue;

		{
			qaws_scalar ea_x = pc[0] - pb[0], ea_y = pc[1] - pb[1], ea_z = pc[2] - pb[2];
			qaws_scalar eb_x = pa[0] - pc[0], eb_y = pa[1] - pc[1], eb_z = pa[2] - pc[2];
			qaws_scalar ec_x = pb[0] - pa[0], ec_y = pb[1] - pa[1], ec_z = pb[2] - pa[2];

			qaws_scalar nxa_x = ny * ea_z - nz * ea_y;
			qaws_scalar nxa_y = nz * ea_x - nx * ea_z;
			qaws_scalar nxa_z = nx * ea_y - ny * ea_x;

			qaws_scalar nxb_x = ny * eb_z - nz * eb_y;
			qaws_scalar nxb_y = nz * eb_x - nx * eb_z;
			qaws_scalar nxb_z = nx * eb_y - ny * eb_x;

			qaws_scalar nxc_x = ny * ec_z - nz * ec_y;
			qaws_scalar nxc_y = nz * ec_x - nx * ec_z;
			qaws_scalar nxc_z = nx * ec_y - ny * ec_x;

			qaws_scalar inv_area2 = QAWS_ONE / (area2 * area2);

			qaws_scalar gx = (heat[va] * nxa_x + heat[vb] * nxb_x + heat[vc] * nxc_x) * inv_area2;
			qaws_scalar gy = (heat[va] * nxa_y + heat[vb] * nxb_y + heat[vc] * nxc_y) * inv_area2;
			qaws_scalar gz = (heat[va] * nxa_z + heat[vb] * nxb_z + heat[vc] * nxc_z) * inv_area2;

			qaws_scalar glen = QAWS_SQRT(gx * gx + gy * gy + gz * gz);
			qaws_scalar Xx, Xy, Xz;

			if (glen < QAWS_LITERAL(1e-15))
			{
				Xx = QAWS_ZERO; Xy = QAWS_ZERO; Xz = QAWS_ZERO;
			}
			else
			{
				Xx = -gx / glen;
				Xy = -gy / glen;
				Xz = -gz / glen;
			}

			{
				qaws_scalar dot_a = e0x * e1x + e0y * e1y + e0z * e1z;
				qaws_scalar cross_a_x = e0y * e1z - e0z * e1y;
				qaws_scalar cross_a_y = e0z * e1x - e0x * e1z;
				qaws_scalar cross_a_z = e0x * e1y - e0y * e1x;
				qaws_scalar cross_a_mag = QAWS_SQRT(cross_a_x * cross_a_x + cross_a_y * cross_a_y + cross_a_z * cross_a_z);
				qaws_scalar cot_a = (cross_a_mag > QAWS_LITERAL(1e-15)) ? dot_a / cross_a_mag : QAWS_ZERO;

				qaws_scalar f0x = pa[0] - pb[0], f0y = pa[1] - pb[1], f0z = pa[2] - pb[2];
				qaws_scalar f1x = pc[0] - pb[0], f1y = pc[1] - pb[1], f1z = pc[2] - pb[2];
				qaws_scalar dot_b = f0x * f1x + f0y * f1y + f0z * f1z;
				qaws_scalar cross_b_x = f0y * f1z - f0z * f1y;
				qaws_scalar cross_b_y = f0z * f1x - f0x * f1z;
				qaws_scalar cross_b_z = f0x * f1y - f0y * f1x;
				qaws_scalar cross_b_mag = QAWS_SQRT(cross_b_x * cross_b_x + cross_b_y * cross_b_y + cross_b_z * cross_b_z);
				qaws_scalar cot_b = (cross_b_mag > QAWS_LITERAL(1e-15)) ? dot_b / cross_b_mag : QAWS_ZERO;

				qaws_scalar g0x = pa[0] - pc[0], g0y = pa[1] - pc[1], g0z = pa[2] - pc[2];
				qaws_scalar g1x = pb[0] - pc[0], g1y = pb[1] - pc[1], g1z = pb[2] - pc[2];
				qaws_scalar dot_c = g0x * g1x + g0y * g1y + g0z * g1z;
				qaws_scalar cross_c_x = g0y * g1z - g0z * g1y;
				qaws_scalar cross_c_y = g0z * g1x - g0x * g1z;
				qaws_scalar cross_c_z = g0x * g1y - g0y * g1x;
				qaws_scalar cross_c_mag = QAWS_SQRT(cross_c_x * cross_c_x + cross_c_y * cross_c_y + cross_c_z * cross_c_z);
				qaws_scalar cot_c = (cross_c_mag > QAWS_LITERAL(1e-15)) ? dot_c / cross_c_mag : QAWS_ZERO;

				if (cot_a < QAWS_LITERAL(-1e6)) cot_a = QAWS_LITERAL(-1e6);
				if (cot_a > QAWS_LITERAL(1e6)) cot_a = QAWS_LITERAL(1e6);
				if (cot_b < QAWS_LITERAL(-1e6)) cot_b = QAWS_LITERAL(-1e6);
				if (cot_b > QAWS_LITERAL(1e6)) cot_b = QAWS_LITERAL(1e6);
				if (cot_c < QAWS_LITERAL(-1e6)) cot_c = QAWS_LITERAL(-1e6);
				if (cot_c > QAWS_LITERAL(1e6)) cot_c = QAWS_LITERAL(1e6);

				div_X[va] += QAWS_LITERAL(0.5) * (
					cot_c * (ec_x * Xx + ec_y * Xy + ec_z * Xz) +
					cot_b * (e1x * Xx + e1y * Xy + e1z * Xz));

				div_X[vb] += QAWS_LITERAL(0.5) * (
					cot_a * (ea_x * Xx + ea_y * Xy + ea_z * Xz) +
					cot_c * (-ec_x * Xx + -ec_y * Xy + -ec_z * Xz));

				div_X[vc] += QAWS_LITERAL(0.5) * (
					cot_b * (-e1x * Xx + -e1y * Xy + -e1z * Xz) +
					cot_a * (-ea_x * Xx + -ea_y * Xy + -ea_z * Xz));
			}
		}
	}

	/* Step IV: Solve Poisson equation for distance.
	   Our L has negative diagonal: L[i,i] = -Σw, L[i,j] = w.
	   The standard cotan Laplacian L_c (positive semidefinite) = -L.
	   Crane's paper: L_c * φ = -div(X).
	   So: (-L) * φ = -div(X), or (-L + εI) * φ = -div(X). */
	{
		qaws_scalar eps = QAWS_LITERAL(1e-6);
		qaws_scalar min_phi;

		status = qaws_internal_csr_alloc(&L_reg, vert_count, vert_count, L.nnz);
		if (status != QAWS_STATUS_OK) goto cleanup;

		memcpy(L_reg.row_ptr, L.row_ptr, (vert_count + 1) * sizeof(unsigned int));
		memcpy(L_reg.col_idx, L.col_idx, L.nnz * sizeof(unsigned int));

		/* -L is positive semidefinite */
		for (vi = 0; vi < L.nnz; vi++)
			L_reg.values[vi] = -L.values[vi];

		for (vi = 0; vi < vert_count; vi++)
			L_reg.values[L_reg.row_ptr[vi]] += eps;

		for (vi = 0; vi < vert_count; vi++)
			div_X[vi] = -div_X[vi];

		memset(heat, 0, vert_count * sizeof(qaws_scalar));
		status = qaws_internal_cg_solve(&L_reg, div_X, heat, 0, QAWS_ZERO);
		if (status != QAWS_STATUS_OK) goto cleanup;

		min_phi = heat[0];
		for (vi = 1; vi < vert_count; vi++)
		{
			if (heat[vi] < min_phi)
				min_phi = heat[vi];
		}
		for (vi = 0; vi < vert_count; vi++)
			heat[vi] -= min_phi;
	}

	/* Transfer ownership of distance array */
	*out_distances = heat;
	heat = NULL;

cleanup:
	qaws_internal_csr_free(&L_reg);
	qaws_internal_csr_free(&A_minus_tL);
	qaws_internal_csr_free(&L);
	free(div_X);
	free(heat);
	free(rhs);
	free(mass_diag);
	return status;
}

/* ------------------------------------------------------------------ */
/*  Geodesic distance via heat method (Crane et al. 2017)              */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_compute_geodesic_distance(
	qaws_surface const* surface,
	qaws_scalar source_u, qaws_scalar source_v,
	qaws_tessellation_desc const* tessellation_desc,
	qaws_scalar* out_u,
	qaws_scalar* out_v,
	qaws_scalar* out_distances,
	unsigned int capacity,
	unsigned int* out_count)
{
	qaws_tessellation_vertex* verts = NULL;
	unsigned int* indices = NULL;
	unsigned int vert_count = 0, idx_count = 0;
	qaws_scalar* positions = NULL;
	qaws_scalar* uvs_buf = NULL;
	qaws_internal_mesh mesh;
	unsigned int src_v, vi;
	qaws_scalar* distances = NULL;
	qaws_status status;

	if (!surface || !out_distances || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_count = 0;
	memset(&mesh, 0, sizeof(mesh));

	status = geodesic_tessellate_and_build(surface, tessellation_desc,
		&verts, &indices, &vert_count, &idx_count,
		&positions, &uvs_buf, &mesh);

	if (status != QAWS_STATUS_OK)
		return status;

	if (capacity < vert_count)
	{
		*out_count = vert_count;
		qaws_internal_mesh_destroy(&mesh);
		free(verts); free(indices); free(positions); free(uvs_buf);
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	src_v = qaws_internal_mesh_closest_vertex_uv(&mesh, source_u, source_v);

	status = heat_solve_distances(&mesh, positions, vert_count, src_v, &distances);
	if (status != QAWS_STATUS_OK) goto cleanup;

	for (vi = 0; vi < vert_count; vi++)
	{
		out_distances[vi] = distances[vi];
		if (out_u) out_u[vi] = uvs_buf[vi * 2 + 0];
		if (out_v) out_v[vi] = uvs_buf[vi * 2 + 1];
	}
	*out_count = vert_count;

cleanup:
	free(distances);
	qaws_internal_mesh_destroy(&mesh);
	free(verts);
	free(indices);
	free(positions);
	free(uvs_buf);

	return status;
}

/* ------------------------------------------------------------------ */
/*  Geodesic path via heat method + gradient descent                    */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_compute_geodesic_heat(
	qaws_surface const* surface,
	qaws_scalar start_u, qaws_scalar start_v,
	qaws_scalar end_u, qaws_scalar end_v,
	qaws_tessellation_desc const* tessellation_desc,
	qaws_curve** out_curve)
{
	qaws_tessellation_vertex* verts = NULL;
	unsigned int* indices = NULL;
	unsigned int vert_count = 0, idx_count = 0;
	qaws_scalar* positions = NULL;
	qaws_scalar* uvs_buf = NULL;
	qaws_internal_mesh mesh;
	unsigned int src_vi;
	qaws_scalar* distances = NULL;
	qaws_status status;

	if (!surface || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;
	memset(&mesh, 0, sizeof(mesh));

	status = geodesic_tessellate_and_build(surface, tessellation_desc,
		&verts, &indices, &vert_count, &idx_count,
		&positions, &uvs_buf, &mesh);
	if (status != QAWS_STATUS_OK) return status;

	src_vi = qaws_internal_mesh_closest_vertex_uv(&mesh, start_u, start_v);

	status = heat_solve_distances(&mesh, positions, vert_count, src_vi, &distances);
	if (status != QAWS_STATUS_OK) goto cleanup;

	/* Trace path from end toward source by continuous gradient descent
	   through triangle interiors (Crane et al. 2017).
	   In each triangle, lay out vertices in 2D from 3D edge lengths,
	   compute gradient of the piecewise-linear distance field, walk
	   along -gradient until hitting an edge, cross to the adjacent
	   triangle.  To prevent oscillation, never exit through the edge
	   we entered from; instead walk to the vertex with minimum
	   distance. */
	{
		unsigned int max_trace = 20000;
		qaws_scalar* trace_coords = NULL;
		qaws_scalar* path_coords = NULL;
		unsigned int trace_count = 0;
		unsigned int cur_tri;
		qaws_scalar bary[3];
		unsigned int entry_edge; /* edge we entered from (UINT_MAX = none) */
		unsigned int trace_iter;

		trace_coords = (qaws_scalar*)malloc(max_trace * 3 * sizeof(qaws_scalar));
		if (!trace_coords) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

		/* Find triangle containing (end_u, end_v) via UV barycentric test */
		cur_tri = UINT_MAX;
		{
			unsigned int ti;
			for (ti = 0; ti < mesh.tri_count; ti++)
			{
				unsigned int a = mesh.tris[ti].v[0], b = mesh.tris[ti].v[1], c = mesh.tris[ti].v[2];
				qaws_scalar u0 = uvs_buf[a*2], v0 = uvs_buf[a*2+1];
				qaws_scalar u1 = uvs_buf[b*2], v1 = uvs_buf[b*2+1];
				qaws_scalar u2 = uvs_buf[c*2], v2 = uvs_buf[c*2+1];
				qaws_scalar det = (u1-u0)*(v2-v0) - (u2-u0)*(v1-v0);
				qaws_scalar b1, b2, b0;
				if (QAWS_FABS(det) < QAWS_LITERAL(1e-15)) continue;
				b1 = ((end_u-u0)*(v2-v0) - (u2-u0)*(end_v-v0)) / det;
				b2 = ((u1-u0)*(end_v-v0) - (end_u-u0)*(v1-v0)) / det;
				b0 = QAWS_ONE - b1 - b2;
				if (b0 >= QAWS_LITERAL(-1e-4) && b1 >= QAWS_LITERAL(-1e-4) && b2 >= QAWS_LITERAL(-1e-4))
				{
					qaws_scalar s;
					cur_tri = ti;
					if (b0 < QAWS_ZERO) b0 = QAWS_ZERO;
					if (b1 < QAWS_ZERO) b1 = QAWS_ZERO;
					if (b2 < QAWS_ZERO) b2 = QAWS_ZERO;
					s = b0 + b1 + b2;
					bary[0] = b0/s; bary[1] = b1/s; bary[2] = b2/s;
					break;
				}
			}
		}
		/* Fallback: closest vertex */
		if (cur_tri == UINT_MAX)
		{
			unsigned int cv = qaws_internal_mesh_closest_vertex_uv(&mesh, end_u, end_v);
			if (cv < vert_count && mesh.vert_edge_offset[cv] < mesh.vert_edge_offset[cv+1])
			{
				unsigned int ei = mesh.vert_edges[mesh.vert_edge_offset[cv]];
				cur_tri = mesh.edges[ei].tri[0];
				if (cur_tri == UINT_MAX) cur_tri = mesh.edges[ei].tri[1];
			}
			if (cur_tri != UINT_MAX)
			{
				int mi;
				bary[0] = bary[1] = bary[2] = QAWS_ZERO;
				for (mi = 0; mi < 3; mi++)
					if (mesh.tris[cur_tri].v[mi] == cv) bary[mi] = QAWS_ONE;
			}
		}
		if (cur_tri == UINT_MAX)
		{
			free(trace_coords);
			status = QAWS_STATUS_NUMERICAL_FAILURE;
			goto cleanup;
		}
		entry_edge = UINT_MAX;

		/* Store initial position */
		{
			unsigned int a = mesh.tris[cur_tri].v[0], b = mesh.tris[cur_tri].v[1], c = mesh.tris[cur_tri].v[2];
			trace_coords[0] = bary[0]*positions[a*3] + bary[1]*positions[b*3] + bary[2]*positions[c*3];
			trace_coords[1] = bary[0]*positions[a*3+1] + bary[1]*positions[b*3+1] + bary[2]*positions[c*3+1];
			trace_coords[2] = bary[0]*positions[a*3+2] + bary[1]*positions[b*3+2] + bary[2]*positions[c*3+2];
			trace_count = 1;
		}

		for (trace_iter = 0; trace_iter < max_trace - 2; trace_iter++)
		{
			unsigned int tv[3];
			qaws_scalar d[3], q[3][2];
			qaws_scalar l01, l02, l12, cos_a, sin_a;
			qaws_scalar e0x, e0y, e1x, e1y, e2x, e2y;
			qaws_scalar area2, gx, gy, gmag;
			qaws_scalar db[3];
			qaws_scalar t_exit;
			int exit_vi, jj, kk;
			unsigned int edge_idx, adj_tri, adj_side, opp_v;
			qaws_scalar t_cross;

			tv[0] = mesh.tris[cur_tri].v[0];
			tv[1] = mesh.tris[cur_tri].v[1];
			tv[2] = mesh.tris[cur_tri].v[2];
			d[0] = distances[tv[0]]; d[1] = distances[tv[1]]; d[2] = distances[tv[2]];

			/* Stop if source vertex is in this triangle */
			if (tv[0] == src_vi || tv[1] == src_vi || tv[2] == src_vi)
			{
				trace_coords[trace_count*3+0] = positions[src_vi*3];
				trace_coords[trace_count*3+1] = positions[src_vi*3+1];
				trace_coords[trace_count*3+2] = positions[src_vi*3+2];
				trace_count++;
				break;
			}

			/* 2D layout from 3D edge lengths */
			{
				qaws_scalar dx, dy, dz;
				dx = positions[tv[1]*3]-positions[tv[0]*3];
				dy = positions[tv[1]*3+1]-positions[tv[0]*3+1];
				dz = positions[tv[1]*3+2]-positions[tv[0]*3+2];
				l01 = QAWS_SQRT(dx*dx+dy*dy+dz*dz);
				dx = positions[tv[2]*3]-positions[tv[0]*3];
				dy = positions[tv[2]*3+1]-positions[tv[0]*3+1];
				dz = positions[tv[2]*3+2]-positions[tv[0]*3+2];
				l02 = QAWS_SQRT(dx*dx+dy*dy+dz*dz);
				dx = positions[tv[2]*3]-positions[tv[1]*3];
				dy = positions[tv[2]*3+1]-positions[tv[1]*3+1];
				dz = positions[tv[2]*3+2]-positions[tv[1]*3+2];
				l12 = QAWS_SQRT(dx*dx+dy*dy+dz*dz);
			}
			if (l01 < QAWS_LITERAL(1e-15) || l02 < QAWS_LITERAL(1e-15)) break;
			q[0][0] = QAWS_ZERO; q[0][1] = QAWS_ZERO;
			q[1][0] = l01; q[1][1] = QAWS_ZERO;
			cos_a = (l01*l01 + l02*l02 - l12*l12) / (QAWS_LITERAL(2.0)*l01*l02);
			if (cos_a > QAWS_ONE) cos_a = QAWS_ONE;
			if (cos_a < -QAWS_ONE) cos_a = -QAWS_ONE;
			sin_a = QAWS_SQRT(QAWS_ONE - cos_a*cos_a);
			q[2][0] = l02*cos_a; q[2][1] = l02*sin_a;

			/* Gradient of distance in 2D */
			e0x = q[2][0]-q[1][0]; e0y = q[2][1]-q[1][1]; /* edge opp v0 */
			e1x = q[0][0]-q[2][0]; e1y = q[0][1]-q[2][1]; /* edge opp v1 */
			e2x = q[1][0]-q[0][0]; e2y = q[1][1]-q[0][1]; /* edge opp v2 */
			area2 = e2x*(q[2][1]-q[0][1]) - (q[2][0]-q[0][0])*e2y;
			if (QAWS_FABS(area2) < QAWS_LITERAL(1e-20)) break;

			gx = (d[0]*(-e0y) + d[1]*(-e1y) + d[2]*(-e2y)) / area2;
			gy = (d[0]*e0x + d[1]*e1x + d[2]*e2x) / area2;
			gmag = QAWS_SQRT(gx*gx + gy*gy);
			if (gmag < QAWS_LITERAL(1e-15)) break;

			/* Barycentric displacement for direction -gradient */
			{
				qaws_scalar ax = q[0][0]-q[2][0], ay = q[0][1]-q[2][1];
				qaws_scalar bx = q[1][0]-q[2][0], by = q[1][1]-q[2][1];
				qaws_scalar det = ax*by - ay*bx;
				if (QAWS_FABS(det) < QAWS_LITERAL(1e-20)) break;
				db[0] = (-gx*by - (-gy)*bx) / det;
				db[1] = (ax*(-gy) - ay*(-gx)) / det;
				db[2] = -db[0] - db[1];
			}

			/* Find exit edge (smallest t>0 where bary[i]+t*db[i]=0).
			   Skip the entry edge to prevent oscillation. */
			t_exit = QAWS_LITERAL(1e30);
			exit_vi = -1;
			{
				int ki;
				for (ki = 0; ki < 3; ki++)
				{
					if (db[ki] < QAWS_LITERAL(-1e-15))
					{
						qaws_scalar t = -bary[ki] / db[ki];
						if (t > QAWS_LITERAL(1e-12) && t < t_exit)
						{
							/* Check if this edge is the entry edge */
							if (entry_edge != UINT_MAX &&
								mesh.tris[cur_tri].edge[ki] == entry_edge)
								continue;
							t_exit = t;
							exit_vi = ki;
						}
					}
				}
			}

			/* If no forward exit, fall back: walk to vertex with min distance,
			   pick one of its adjacent triangles, and continue. */
			if (exit_vi < 0)
			{
				int best_li = 0;
				unsigned int best_v;
				if (d[1] < d[best_li]) best_li = 1;
				if (d[2] < d[best_li]) best_li = 2;
				best_v = tv[best_li];

				/* Add vertex position to trace */
				trace_coords[trace_count*3+0] = positions[best_v*3];
				trace_coords[trace_count*3+1] = positions[best_v*3+1];
				trace_coords[trace_count*3+2] = positions[best_v*3+2];
				trace_count++;
				if (best_v == src_vi) break;

				/* Pick adjacent triangle: the one with the steepest gradient
				   away from this vertex */
				{
					unsigned int ei_s = mesh.vert_edge_offset[best_v];
					unsigned int ei_e = mesh.vert_edge_offset[best_v + 1];
					unsigned int ei2, best_tri = UINT_MAX;
					qaws_scalar best_d_drop = QAWS_ZERO;
					for (ei2 = ei_s; ei2 < ei_e; ei2++)
					{
						unsigned int eidx = mesh.vert_edges[ei2];
						unsigned int s;
						for (s = 0; s < 2; s++)
						{
							unsigned int ti = mesh.edges[eidx].tri[s];
							int li;
							qaws_scalar d_min;
							if (ti == UINT_MAX || ti == cur_tri) continue;
							d_min = distances[best_v];
							for (li = 0; li < 3; li++)
								if (distances[mesh.tris[ti].v[li]] < d_min)
									d_min = distances[mesh.tris[ti].v[li]];
							if (distances[best_v] - d_min > best_d_drop)
							{
								best_d_drop = distances[best_v] - d_min;
								best_tri = ti;
							}
						}
					}
					if (best_tri == UINT_MAX) break;
					cur_tri = best_tri;
				}
				/* Set bary to be at best_v in the new triangle */
				{
					int mi;
					bary[0] = bary[1] = bary[2] = QAWS_ZERO;
					for (mi = 0; mi < 3; mi++)
						if (mesh.tris[cur_tri].v[mi] == best_v)
						{ bary[mi] = QAWS_ONE; break; }
				}
				entry_edge = UINT_MAX;
				continue;
			}

			/* Compute crossing point on exit edge */
			jj = (exit_vi + 1) % 3;
			kk = (exit_vi + 2) % 3;
			{
				qaws_scalar b_j = bary[jj] + t_exit*db[jj];
				qaws_scalar b_k = bary[kk] + t_exit*db[kk];
				qaws_scalar s;
				if (b_j < QAWS_ZERO) b_j = QAWS_ZERO;
				if (b_k < QAWS_ZERO) b_k = QAWS_ZERO;
				s = b_j + b_k;
				if (s < QAWS_LITERAL(1e-15)) break;
				t_cross = b_k / s;
			}

			/* 3D position at crossing */
			{
				unsigned int vj = tv[jj], vk = tv[kk];
				trace_coords[trace_count*3+0] = (QAWS_ONE-t_cross)*positions[vj*3]   + t_cross*positions[vk*3];
				trace_coords[trace_count*3+1] = (QAWS_ONE-t_cross)*positions[vj*3+1] + t_cross*positions[vk*3+1];
				trace_coords[trace_count*3+2] = (QAWS_ONE-t_cross)*positions[vj*3+2] + t_cross*positions[vk*3+2];
				trace_count++;
			}

			/* Cross to adjacent triangle */
			edge_idx = mesh.tris[cur_tri].edge[exit_vi];
			adj_side = (mesh.edges[edge_idx].tri[0] == cur_tri) ? 1u : 0u;
			adj_tri = mesh.edges[edge_idx].tri[adj_side];
			if (adj_tri == UINT_MAX) break; /* boundary */

			/* Barycentric coords in adjacent triangle */
			opp_v = mesh.edges[edge_idx].opposite[adj_side];
			{
				int m, opp_l = -1, a_l, b_l;
				qaws_scalar eps = QAWS_LITERAL(1e-6);
				for (m = 0; m < 3; m++)
					if (mesh.tris[adj_tri].v[m] == opp_v) { opp_l = m; break; }
				if (opp_l < 0) break;
				a_l = (opp_l + 1) % 3;
				b_l = (opp_l + 2) % 3;
				bary[opp_l] = eps;
				if (mesh.tris[adj_tri].v[a_l] == tv[jj])
				{ bary[a_l] = (QAWS_ONE-t_cross)*(QAWS_ONE-eps); bary[b_l] = t_cross*(QAWS_ONE-eps); }
				else
				{ bary[a_l] = t_cross*(QAWS_ONE-eps); bary[b_l] = (QAWS_ONE-t_cross)*(QAWS_ONE-eps); }
			}
			entry_edge = edge_idx;
			cur_tri = adj_tri;
		}

		if (trace_count < 2)
		{
			free(trace_coords);
			status = QAWS_STATUS_NUMERICAL_FAILURE;
			goto cleanup;
		}

		/* Reverse trace (end→start → start→end), fix exact endpoints */
		path_coords = (qaws_scalar*)malloc(trace_count * 3 * sizeof(qaws_scalar));
		if (!path_coords) { free(trace_coords); status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }
		{
			unsigned int pi;
			for (pi = 0; pi < trace_count; pi++)
			{
				path_coords[pi*3+0] = trace_coords[(trace_count-1-pi)*3+0];
				path_coords[pi*3+1] = trace_coords[(trace_count-1-pi)*3+1];
				path_coords[pi*3+2] = trace_coords[(trace_count-1-pi)*3+2];
			}
		}
		{
			qaws_surface_eval_result sr;
			memset(&sr, 0, sizeof(sr));
			qaws_surface_evaluate(surface, start_u, start_v, QAWS_SURFACE_EVAL_POSITION, &sr);
			path_coords[0] = sr.position.x; path_coords[1] = sr.position.y; path_coords[2] = sr.position.z;
		}
		{
			qaws_surface_eval_result sr;
			memset(&sr, 0, sizeof(sr));
			qaws_surface_evaluate(surface, end_u, end_v, QAWS_SURFACE_EVAL_POSITION, &sr);
			path_coords[(trace_count-1)*3+0] = sr.position.x;
			path_coords[(trace_count-1)*3+1] = sr.position.y;
			path_coords[(trace_count-1)*3+2] = sr.position.z;
		}

		status = create_path_bspline(path_coords, trace_count, 50, out_curve);
		free(trace_coords);
		free(path_coords);
	}

cleanup:
	free(distances);
	qaws_internal_mesh_destroy(&mesh);
	free(verts);
	free(indices);
	free(positions);
	free(uvs_buf);

	return status;
}
