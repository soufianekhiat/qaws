#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_PREPARE_H_F32
#define QAWS_PREPARE_H_F32
#define QAWS_PREPARE_H_BODY
#endif
#else
#ifndef QAWS_PREPARE_H_F64
#define QAWS_PREPARE_H_F64
#define QAWS_PREPARE_H_BODY
#endif
#endif
#ifdef QAWS_PREPARE_H_BODY
#undef QAWS_PREPARE_H_BODY
#define QAWS_PREPARE_H

#include "qaws_platform.h"

#if QAWS_BACKEND == QAWS_BACKEND_C

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Prepare functions — C only.
 *
 * Each takes raw geometry input and outputs pre-computed flat arrays
 * suitable for upload to any backend. The caller owns the output memory.
 */

/* ===================================================================
 * Hermite: points + tangents → a,b,c,d coefficients per span
 *
 * out_a..out_d: caller-allocated arrays of [span_count] vec2/vec3.
 * span_count = point_count - 1.
 * =================================================================== */

qaws_status qaws_hermite_prepare_2d(
	qaws_vec2* out_a,
	qaws_vec2* out_b,
	qaws_vec2* out_c,
	qaws_vec2* out_d,
	const qaws_vec2* points,
	const qaws_vec2* tangents,
	unsigned int point_count);

qaws_status qaws_hermite_prepare_3d(
	qaws_vec3* out_a,
	qaws_vec3* out_b,
	qaws_vec3* out_c,
	qaws_vec3* out_d,
	const qaws_vec3* points,
	const qaws_vec3* tangents,
	unsigned int point_count);

/* ===================================================================
 * Catmull-Rom: interpolation points → a,b,c,d + span_count
 * =================================================================== */

qaws_status qaws_catmull_rom_prepare_2d(
	qaws_vec2* out_a,
	qaws_vec2* out_b,
	qaws_vec2* out_c,
	qaws_vec2* out_d,
	unsigned int* out_span_count,
	const qaws_vec2* points,
	unsigned int point_count,
	qaws_parameterization param,
	int closed);

qaws_status qaws_catmull_rom_prepare_3d(
	qaws_vec3* out_a,
	qaws_vec3* out_b,
	qaws_vec3* out_c,
	qaws_vec3* out_d,
	unsigned int* out_span_count,
	const qaws_vec3* points,
	unsigned int point_count,
	qaws_parameterization param,
	int closed);

/* ===================================================================
 * Trajectory: positions + times → a,b,c,d + span_count
 * =================================================================== */

qaws_status qaws_trajectory_prepare_2d(
	qaws_vec2* out_a,
	qaws_vec2* out_b,
	qaws_vec2* out_c,
	qaws_vec2* out_d,
	unsigned int* out_span_count,
	const qaws_vec2* positions,
	const qaws_scalar* times,
	unsigned int key_count,
	unsigned int degree,
	int closed);

qaws_status qaws_trajectory_prepare_3d(
	qaws_vec3* out_a,
	qaws_vec3* out_b,
	qaws_vec3* out_c,
	qaws_vec3* out_d,
	unsigned int* out_span_count,
	const qaws_vec3* positions,
	const qaws_scalar* times,
	unsigned int key_count,
	unsigned int degree,
	int closed);

/* ===================================================================
 * Bezier: CPs → derivative CPs (D1, D2, D3)
 * =================================================================== */

qaws_status qaws_bezier_prepare_deriv_2d(
	qaws_vec2* out_d1_cp,
	qaws_vec2* out_d2_cp,
	qaws_vec2* out_d3_cp,
	const qaws_vec2* control_points,
	unsigned int degree);

qaws_status qaws_bezier_prepare_deriv_3d(
	qaws_vec3* out_d1_cp,
	qaws_vec3* out_d2_cp,
	qaws_vec3* out_d3_cp,
	const qaws_vec3* control_points,
	unsigned int degree);

/* ===================================================================
 * Rational Bezier: CPs + weights → homogeneous weighted CPs
 *
 * out_weighted_cp: (degree+1) * 3 scalars for 2D (wx, wy, w)
 *                  (degree+1) * 4 scalars for 3D (wx, wy, wz, w)
 * =================================================================== */

qaws_status qaws_rational_bezier_prepare_2d(
	qaws_scalar* out_weighted_cp,
	const qaws_vec2* control_points,
	const qaws_scalar* weights,
	unsigned int degree);

qaws_status qaws_rational_bezier_prepare_3d(
	qaws_scalar* out_weighted_cp,
	const qaws_vec3* control_points,
	const qaws_scalar* weights,
	unsigned int degree);

/* ===================================================================
 * Arc: center, radii, angles → axis vectors + theta range
 * =================================================================== */

qaws_status qaws_arc_prepare_2d(
	qaws_vec2* out_axis_u,
	qaws_vec2* out_axis_v,
	qaws_scalar* out_theta_start,
	qaws_scalar* out_theta_range,
	qaws_vec2 center,
	qaws_scalar radius_x,
	qaws_scalar radius_y,
	qaws_scalar start_angle,
	qaws_scalar end_angle);

qaws_status qaws_arc_prepare_3d(
	qaws_vec3* out_axis_u,
	qaws_vec3* out_axis_v,
	qaws_scalar* out_theta_start,
	qaws_scalar* out_theta_range,
	qaws_vec3 center,
	qaws_vec3 normal,
	qaws_scalar radius_x,
	qaws_scalar radius_y,
	qaws_scalar start_angle,
	qaws_scalar end_angle);

#endif /* QAWS_BACKEND == QAWS_BACKEND_C */

#endif /* QAWS_PREPARE_H */
