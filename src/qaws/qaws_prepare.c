#include "qaws_prepare.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_platform.h"
#include <string.h>
#include <stdlib.h>

/* ===================================================================
 * Hermite prepare
 * =================================================================== */

qaws_status qaws_hermite_prepare_2d(
	const qaws_vec2* points,
	const qaws_vec2* tangents,
	unsigned int point_count,
	qaws_vec2* out_a,
	qaws_vec2* out_b,
	qaws_vec2* out_c,
	qaws_vec2* out_d)
{
	unsigned int s;
	if (!points || !tangents || !out_a || !out_b || !out_c || !out_d)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (point_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (s = 0; s < point_count - 1; s++) {
		qaws_scalar p0x = points[s].x,   p0y = points[s].y;
		qaws_scalar p1x = points[s+1].x, p1y = points[s+1].y;
		qaws_scalar m0x = tangents[s].x,  m0y = tangents[s].y;
		qaws_scalar m1x = tangents[s+1].x, m1y = tangents[s+1].y;

		out_a[s].x = QAWS_LITERAL(2.0) * p0x + m0x - QAWS_LITERAL(2.0) * p1x + m1x;
		out_a[s].y = QAWS_LITERAL(2.0) * p0y + m0y - QAWS_LITERAL(2.0) * p1y + m1y;

		out_b[s].x = -QAWS_LITERAL(3.0) * p0x - QAWS_LITERAL(2.0) * m0x + QAWS_LITERAL(3.0) * p1x - m1x;
		out_b[s].y = -QAWS_LITERAL(3.0) * p0y - QAWS_LITERAL(2.0) * m0y + QAWS_LITERAL(3.0) * p1y - m1y;

		out_c[s].x = m0x;
		out_c[s].y = m0y;

		out_d[s].x = p0x;
		out_d[s].y = p0y;
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_hermite_prepare_3d(
	const qaws_vec3* points,
	const qaws_vec3* tangents,
	unsigned int point_count,
	qaws_vec3* out_a,
	qaws_vec3* out_b,
	qaws_vec3* out_c,
	qaws_vec3* out_d)
{
	unsigned int s;
	if (!points || !tangents || !out_a || !out_b || !out_c || !out_d)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (point_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (s = 0; s < point_count - 1; s++) {
		qaws_scalar p0x = points[s].x,   p0y = points[s].y,   p0z = points[s].z;
		qaws_scalar p1x = points[s+1].x, p1y = points[s+1].y, p1z = points[s+1].z;
		qaws_scalar m0x = tangents[s].x,  m0y = tangents[s].y,  m0z = tangents[s].z;
		qaws_scalar m1x = tangents[s+1].x, m1y = tangents[s+1].y, m1z = tangents[s+1].z;

		out_a[s].x = QAWS_LITERAL(2.0) * p0x + m0x - QAWS_LITERAL(2.0) * p1x + m1x;
		out_a[s].y = QAWS_LITERAL(2.0) * p0y + m0y - QAWS_LITERAL(2.0) * p1y + m1y;
		out_a[s].z = QAWS_LITERAL(2.0) * p0z + m0z - QAWS_LITERAL(2.0) * p1z + m1z;

		out_b[s].x = -QAWS_LITERAL(3.0) * p0x - QAWS_LITERAL(2.0) * m0x + QAWS_LITERAL(3.0) * p1x - m1x;
		out_b[s].y = -QAWS_LITERAL(3.0) * p0y - QAWS_LITERAL(2.0) * m0y + QAWS_LITERAL(3.0) * p1y - m1y;
		out_b[s].z = -QAWS_LITERAL(3.0) * p0z - QAWS_LITERAL(2.0) * m0z + QAWS_LITERAL(3.0) * p1z - m1z;

		out_c[s].x = m0x; out_c[s].y = m0y; out_c[s].z = m0z;
		out_d[s].x = p0x; out_d[s].y = p0y; out_d[s].z = p0z;
	}

	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Catmull-Rom prepare (simplified — delegates internally)
 * =================================================================== */

static qaws_scalar cr_compute_alpha(qaws_parameterization param)
{
	switch (param) {
	case QAWS_PARAMETERIZATION_UNIFORM:     return QAWS_ZERO;
	case QAWS_PARAMETERIZATION_CHORDAL:     return QAWS_ONE;
	case QAWS_PARAMETERIZATION_CENTRIPETAL:  return QAWS_LITERAL(0.5);
	case QAWS_PARAMETERIZATION_DEFAULT:
	default:                                 return QAWS_LITERAL(0.5);
	}
}

static void cr_segment_coeffs(
	const qaws_scalar* P0, const qaws_scalar* P1,
	const qaws_scalar* P2, const qaws_scalar* P3,
	qaws_scalar t0, qaws_scalar t1, qaws_scalar t2, qaws_scalar t3,
	unsigned int dim_count, qaws_scalar* coeffs_out)
{
	unsigned int d;
	qaws_scalar dt10 = t1 - t0;
	qaws_scalar dt21 = t2 - t1;
	qaws_scalar dt20 = t2 - t0;
	qaws_scalar dt32 = t3 - t2;
	qaws_scalar dt31 = t3 - t1;

	if (dt10 < QAWS_LITERAL(1e-10) && dt10 > -QAWS_LITERAL(1e-10)) dt10 = QAWS_ZERO;
	if (dt21 < QAWS_LITERAL(1e-10) && dt21 > -QAWS_LITERAL(1e-10)) dt21 = QAWS_ZERO;
	if (dt20 < QAWS_LITERAL(1e-10) && dt20 > -QAWS_LITERAL(1e-10)) dt20 = QAWS_ZERO;
	if (dt32 < QAWS_LITERAL(1e-10) && dt32 > -QAWS_LITERAL(1e-10)) dt32 = QAWS_ZERO;
	if (dt31 < QAWS_LITERAL(1e-10) && dt31 > -QAWS_LITERAL(1e-10)) dt31 = QAWS_ZERO;

	for (d = 0; d < dim_count; d++) {
		qaws_scalar p0 = P0[d], p1 = P1[d], p2 = P2[d], p3 = P3[d];
		qaws_scalar m1 = QAWS_ZERO, m2 = QAWS_ZERO, M1, M2;

		if (dt21 != QAWS_ZERO && dt20 != QAWS_ZERO && dt10 != QAWS_ZERO)
			m1 = (p2 - p1) / dt21 - (p2 - p0) / dt20 + (p1 - p0) / dt10;

		if (dt32 != QAWS_ZERO && dt31 != QAWS_ZERO && dt21 != QAWS_ZERO)
			m2 = (p3 - p2) / dt32 - (p3 - p1) / dt31 + (p2 - p1) / dt21;

		M1 = m1 * dt21;
		M2 = m2 * dt21;

		coeffs_out[d * 4 + 3] = p1;
		coeffs_out[d * 4 + 2] = M1;
		coeffs_out[d * 4 + 1] = QAWS_LITERAL(3.0) * (p2 - p1) - QAWS_LITERAL(2.0) * M1 - M2;
		coeffs_out[d * 4 + 0] = QAWS_LITERAL(2.0) * (p1 - p2) + M1 + M2;
	}
}

qaws_status qaws_catmull_rom_prepare_2d(
	const qaws_vec2* points,
	unsigned int point_count,
	qaws_parameterization param,
	int closed,
	qaws_vec2* out_a,
	qaws_vec2* out_b,
	qaws_vec2* out_c,
	qaws_vec2* out_d,
	unsigned int* out_span_count)
{
	unsigned int const dim_count = 2;
	unsigned int span_count, s;
	qaws_scalar alpha;
	const qaws_scalar* pts;

	if (!points || !out_a || !out_b || !out_c || !out_d || !out_span_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	pts = (const qaws_scalar*)points;
	alpha = cr_compute_alpha(param);

	if (closed) {
		if (point_count < 3) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = point_count;
	} else {
		if (point_count < 4) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = point_count - 3;
	}

	*out_span_count = span_count;

	for (s = 0; s < span_count; s++) {
		qaws_scalar coeffs[8]; /* dim_count * 4 */
		unsigned int i0, i1, i2, i3;
		qaws_scalar t0, t1, t2, t3;

		if (closed) {
			unsigned int N = point_count;
			i0 = (s == 0) ? N - 1 : s - 1;
			i1 = s;
			i2 = (s + 1) % N;
			i3 = (s + 2) % N;
		} else {
			i0 = s; i1 = s + 1; i2 = s + 2; i3 = s + 3;
		}

		/* Compute local knot values */
		if (closed) {
			qaws_scalar dist_sq, diff;
			unsigned int dd;

			t1 = QAWS_ZERO;

			dist_sq = QAWS_ZERO;
			for (dd = 0; dd < dim_count; dd++) { diff = pts[i1 * dim_count + dd] - pts[i0 * dim_count + dd]; dist_sq += diff * diff; }
			t0 = (alpha == QAWS_ZERO) ? -QAWS_ONE : -QAWS_POW(QAWS_SQRT(dist_sq), alpha);

			dist_sq = QAWS_ZERO;
			for (dd = 0; dd < dim_count; dd++) { diff = pts[i2 * dim_count + dd] - pts[i1 * dim_count + dd]; dist_sq += diff * diff; }
			t2 = (alpha == QAWS_ZERO) ? QAWS_ONE : QAWS_POW(QAWS_SQRT(dist_sq), alpha);

			dist_sq = QAWS_ZERO;
			for (dd = 0; dd < dim_count; dd++) { diff = pts[i3 * dim_count + dd] - pts[i2 * dim_count + dd]; dist_sq += diff * diff; }
			t3 = (alpha == QAWS_ZERO) ? t2 + QAWS_ONE : t2 + QAWS_POW(QAWS_SQRT(dist_sq), alpha);
		} else {
			/* Use global knot parameterization — for simplicity, recompute locally */
			qaws_scalar dist_sq, diff;
			unsigned int dd;

			t1 = QAWS_ZERO;

			dist_sq = QAWS_ZERO;
			for (dd = 0; dd < dim_count; dd++) { diff = pts[i1 * dim_count + dd] - pts[i0 * dim_count + dd]; dist_sq += diff * diff; }
			t0 = (alpha == QAWS_ZERO) ? -QAWS_ONE : -QAWS_POW(QAWS_SQRT(dist_sq), alpha);

			dist_sq = QAWS_ZERO;
			for (dd = 0; dd < dim_count; dd++) { diff = pts[i2 * dim_count + dd] - pts[i1 * dim_count + dd]; dist_sq += diff * diff; }
			t2 = (alpha == QAWS_ZERO) ? QAWS_ONE : QAWS_POW(QAWS_SQRT(dist_sq), alpha);

			dist_sq = QAWS_ZERO;
			for (dd = 0; dd < dim_count; dd++) { diff = pts[i3 * dim_count + dd] - pts[i2 * dim_count + dd]; dist_sq += diff * diff; }
			t3 = (alpha == QAWS_ZERO) ? t2 + QAWS_ONE : t2 + QAWS_POW(QAWS_SQRT(dist_sq), alpha);
		}

		cr_segment_coeffs(
			&pts[i0 * dim_count], &pts[i1 * dim_count],
			&pts[i2 * dim_count], &pts[i3 * dim_count],
			t0, t1, t2, t3, dim_count, coeffs);

		out_a[s].x = coeffs[0 * 4 + 0]; out_a[s].y = coeffs[1 * 4 + 0];
		out_b[s].x = coeffs[0 * 4 + 1]; out_b[s].y = coeffs[1 * 4 + 1];
		out_c[s].x = coeffs[0 * 4 + 2]; out_c[s].y = coeffs[1 * 4 + 2];
		out_d[s].x = coeffs[0 * 4 + 3]; out_d[s].y = coeffs[1 * 4 + 3];
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_catmull_rom_prepare_3d(
	const qaws_vec3* points,
	unsigned int point_count,
	qaws_parameterization param,
	int closed,
	qaws_vec3* out_a,
	qaws_vec3* out_b,
	qaws_vec3* out_c,
	qaws_vec3* out_d,
	unsigned int* out_span_count)
{
	unsigned int const dim_count = 3;
	unsigned int span_count, s;
	qaws_scalar alpha;
	const qaws_scalar* pts;

	if (!points || !out_a || !out_b || !out_c || !out_d || !out_span_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	pts = (const qaws_scalar*)points;
	alpha = cr_compute_alpha(param);

	if (closed) {
		if (point_count < 3) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = point_count;
	} else {
		if (point_count < 4) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = point_count - 3;
	}

	*out_span_count = span_count;

	for (s = 0; s < span_count; s++) {
		qaws_scalar coeffs[12]; /* dim_count * 4 */
		unsigned int i0, i1, i2, i3;
		qaws_scalar t0, t1, t2, t3;
		qaws_scalar dist_sq, diff;
		unsigned int dd;

		if (closed) {
			unsigned int N = point_count;
			i0 = (s == 0) ? N - 1 : s - 1;
			i1 = s;
			i2 = (s + 1) % N;
			i3 = (s + 2) % N;
		} else {
			i0 = s; i1 = s + 1; i2 = s + 2; i3 = s + 3;
		}

		t1 = QAWS_ZERO;

		dist_sq = QAWS_ZERO;
		for (dd = 0; dd < dim_count; dd++) { diff = pts[i1 * dim_count + dd] - pts[i0 * dim_count + dd]; dist_sq += diff * diff; }
		t0 = (alpha == QAWS_ZERO) ? -QAWS_ONE : -QAWS_POW(QAWS_SQRT(dist_sq), alpha);

		dist_sq = QAWS_ZERO;
		for (dd = 0; dd < dim_count; dd++) { diff = pts[i2 * dim_count + dd] - pts[i1 * dim_count + dd]; dist_sq += diff * diff; }
		t2 = (alpha == QAWS_ZERO) ? QAWS_ONE : QAWS_POW(QAWS_SQRT(dist_sq), alpha);

		dist_sq = QAWS_ZERO;
		for (dd = 0; dd < dim_count; dd++) { diff = pts[i3 * dim_count + dd] - pts[i2 * dim_count + dd]; dist_sq += diff * diff; }
		t3 = (alpha == QAWS_ZERO) ? t2 + QAWS_ONE : t2 + QAWS_POW(QAWS_SQRT(dist_sq), alpha);

		cr_segment_coeffs(
			&pts[i0 * dim_count], &pts[i1 * dim_count],
			&pts[i2 * dim_count], &pts[i3 * dim_count],
			t0, t1, t2, t3, dim_count, coeffs);

		out_a[s].x = coeffs[0 * 4 + 0]; out_a[s].y = coeffs[1 * 4 + 0]; out_a[s].z = coeffs[2 * 4 + 0];
		out_b[s].x = coeffs[0 * 4 + 1]; out_b[s].y = coeffs[1 * 4 + 1]; out_b[s].z = coeffs[2 * 4 + 1];
		out_c[s].x = coeffs[0 * 4 + 2]; out_c[s].y = coeffs[1 * 4 + 2]; out_c[s].z = coeffs[2 * 4 + 2];
		out_d[s].x = coeffs[0 * 4 + 3]; out_d[s].y = coeffs[1 * 4 + 3]; out_d[s].z = coeffs[2 * 4 + 3];
	}

	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Trajectory prepare
 * =================================================================== */

static void traj_compute_velocities(
	const qaws_scalar* positions,
	const qaws_scalar* times,
	unsigned int key_count,
	unsigned int dim_count,
	qaws_scalar* out_velocities)
{
	unsigned int i, d, last;
	qaws_scalar dt, inv_dt;

	if (key_count < 2) return;

	dt = times[1] - times[0];
	if (dt < QAWS_LITERAL(1e-15)) dt = QAWS_LITERAL(1e-15);
	inv_dt = QAWS_ONE / dt;
	for (d = 0; d < dim_count; d++)
		out_velocities[d] = (positions[1 * dim_count + d] - positions[d]) * inv_dt;

	for (i = 1; i < key_count - 1; i++) {
		dt = times[i + 1] - times[i - 1];
		if (dt < QAWS_LITERAL(1e-15)) dt = QAWS_LITERAL(1e-15);
		inv_dt = QAWS_ONE / dt;
		for (d = 0; d < dim_count; d++)
			out_velocities[i * dim_count + d] =
				(positions[(i + 1) * dim_count + d] - positions[(i - 1) * dim_count + d]) * inv_dt;
	}

	last = key_count - 1;
	dt = times[last] - times[last - 1];
	if (dt < QAWS_LITERAL(1e-15)) dt = QAWS_LITERAL(1e-15);
	inv_dt = QAWS_ONE / dt;
	for (d = 0; d < dim_count; d++)
		out_velocities[last * dim_count + d] =
			(positions[last * dim_count + d] - positions[(last - 1) * dim_count + d]) * inv_dt;
}

qaws_status qaws_trajectory_prepare_2d(
	const qaws_vec2* positions,
	const qaws_scalar* times,
	unsigned int key_count,
	unsigned int degree,
	int closed,
	qaws_vec2* out_a,
	qaws_vec2* out_b,
	qaws_vec2* out_c,
	qaws_vec2* out_d,
	unsigned int* out_span_count)
{
	unsigned int const dim_count = 2;
	unsigned int span_count, s, d;
	qaws_scalar vel_stack[128]; /* up to 64 keys in 2D */
	qaws_scalar* velocities = vel_stack;
	const qaws_scalar* pos;

	(void)closed;
	if (!positions || !times || !out_a || !out_b || !out_c || !out_d || !out_span_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (key_count < 2 || degree != 3)
		return QAWS_STATUS_INVALID_ARGUMENT;

	span_count = key_count - 1;
	*out_span_count = span_count;
	pos = (const qaws_scalar*)positions;

	if (key_count * dim_count > 128) {
		velocities = (qaws_scalar*)qaws_internal_malloc(key_count * dim_count * sizeof(qaws_scalar));
		if (!velocities) return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	traj_compute_velocities(pos, times, key_count, dim_count, velocities);

	for (s = 0; s < span_count; s++) {
		qaws_scalar dt = times[s + 1] - times[s];
		for (d = 0; d < dim_count; d++) {
			qaws_scalar p0 = pos[s * dim_count + d];
			qaws_scalar p1 = pos[(s + 1) * dim_count + d];
			qaws_scalar m0 = velocities[s * dim_count + d] * dt;
			qaws_scalar m1 = velocities[(s + 1) * dim_count + d] * dt;
			qaws_scalar a_coeff = QAWS_LITERAL(2.0) * p0 + m0 - QAWS_LITERAL(2.0) * p1 + m1;
			qaws_scalar b_coeff = -QAWS_LITERAL(3.0) * p0 - QAWS_LITERAL(2.0) * m0 + QAWS_LITERAL(3.0) * p1 - m1;
			if (d == 0) { out_a[s].x = a_coeff; out_b[s].x = b_coeff; out_c[s].x = m0; out_d[s].x = p0; }
			else        { out_a[s].y = a_coeff; out_b[s].y = b_coeff; out_c[s].y = m0; out_d[s].y = p0; }
		}
	}

	if (velocities != vel_stack) qaws_internal_free(velocities);
	return QAWS_STATUS_OK;
}

qaws_status qaws_trajectory_prepare_3d(
	const qaws_vec3* positions,
	const qaws_scalar* times,
	unsigned int key_count,
	unsigned int degree,
	int closed,
	qaws_vec3* out_a,
	qaws_vec3* out_b,
	qaws_vec3* out_c,
	qaws_vec3* out_d,
	unsigned int* out_span_count)
{
	unsigned int const dim_count = 3;
	unsigned int span_count, s, d;
	qaws_scalar vel_stack[192]; /* up to 64 keys in 3D */
	qaws_scalar* velocities = vel_stack;
	const qaws_scalar* pos;

	(void)closed;
	if (!positions || !times || !out_a || !out_b || !out_c || !out_d || !out_span_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (key_count < 2 || degree != 3)
		return QAWS_STATUS_INVALID_ARGUMENT;

	span_count = key_count - 1;
	*out_span_count = span_count;
	pos = (const qaws_scalar*)positions;

	if (key_count * dim_count > 192) {
		velocities = (qaws_scalar*)qaws_internal_malloc(key_count * dim_count * sizeof(qaws_scalar));
		if (!velocities) return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	traj_compute_velocities(pos, times, key_count, dim_count, velocities);

	for (s = 0; s < span_count; s++) {
		qaws_scalar dt = times[s + 1] - times[s];
		for (d = 0; d < dim_count; d++) {
			qaws_scalar p0 = pos[s * dim_count + d];
			qaws_scalar p1 = pos[(s + 1) * dim_count + d];
			qaws_scalar m0 = velocities[s * dim_count + d] * dt;
			qaws_scalar m1 = velocities[(s + 1) * dim_count + d] * dt;
			qaws_scalar a_coeff = QAWS_LITERAL(2.0) * p0 + m0 - QAWS_LITERAL(2.0) * p1 + m1;
			qaws_scalar b_coeff = -QAWS_LITERAL(3.0) * p0 - QAWS_LITERAL(2.0) * m0 + QAWS_LITERAL(3.0) * p1 - m1;
			if (d == 0) { out_a[s].x = a_coeff; out_b[s].x = b_coeff; out_c[s].x = m0; out_d[s].x = p0; }
			else if (d == 1) { out_a[s].y = a_coeff; out_b[s].y = b_coeff; out_c[s].y = m0; out_d[s].y = p0; }
			else { out_a[s].z = a_coeff; out_b[s].z = b_coeff; out_c[s].z = m0; out_d[s].z = p0; }
		}
	}

	if (velocities != vel_stack) qaws_internal_free(velocities);
	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Bezier prepare (derivative CPs)
 * =================================================================== */

qaws_status qaws_bezier_prepare_deriv_2d(
	const qaws_vec2* control_points,
	unsigned int degree,
	qaws_vec2* out_d1_cp,
	qaws_vec2* out_d2_cp,
	qaws_vec2* out_d3_cp)
{
	unsigned int i;
	qaws_scalar n;

	if (!control_points || degree < 1)
		return QAWS_STATUS_INVALID_ARGUMENT;

	n = (qaws_scalar)degree;
	if (out_d1_cp) {
		for (i = 0; i < degree; i++) {
			out_d1_cp[i].x = n * (control_points[i+1].x - control_points[i].x);
			out_d1_cp[i].y = n * (control_points[i+1].y - control_points[i].y);
		}
	}

	if (out_d2_cp && degree >= 2) {
		qaws_scalar m = (qaws_scalar)(degree - 1);
		for (i = 0; i < degree - 1; i++) {
			out_d2_cp[i].x = m * (out_d1_cp[i+1].x - out_d1_cp[i].x) / n;
			out_d2_cp[i].y = m * (out_d1_cp[i+1].y - out_d1_cp[i].y) / n;
		}
		/* Re-scale: D2 CPs are n*(n-1) * second differences */
		for (i = 0; i < degree - 1; i++) {
			out_d2_cp[i].x *= n;
			out_d2_cp[i].y *= n;
		}
	}

	if (out_d3_cp && degree >= 3 && out_d2_cp) {
		qaws_scalar k = (qaws_scalar)(degree - 2);
		for (i = 0; i < degree - 2; i++) {
			out_d3_cp[i].x = k * (out_d2_cp[i+1].x - out_d2_cp[i].x) / (n * (qaws_scalar)(degree - 1));
			out_d3_cp[i].y = k * (out_d2_cp[i+1].y - out_d2_cp[i].y) / (n * (qaws_scalar)(degree - 1));
		}
		for (i = 0; i < degree - 2; i++) {
			out_d3_cp[i].x *= n * (qaws_scalar)(degree - 1);
			out_d3_cp[i].y *= n * (qaws_scalar)(degree - 1);
		}
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_bezier_prepare_deriv_3d(
	const qaws_vec3* control_points,
	unsigned int degree,
	qaws_vec3* out_d1_cp,
	qaws_vec3* out_d2_cp,
	qaws_vec3* out_d3_cp)
{
	unsigned int i;
	qaws_scalar n;

	if (!control_points || degree < 1)
		return QAWS_STATUS_INVALID_ARGUMENT;

	n = (qaws_scalar)degree;
	if (out_d1_cp) {
		for (i = 0; i < degree; i++) {
			out_d1_cp[i].x = n * (control_points[i+1].x - control_points[i].x);
			out_d1_cp[i].y = n * (control_points[i+1].y - control_points[i].y);
			out_d1_cp[i].z = n * (control_points[i+1].z - control_points[i].z);
		}
	}

	(void)out_d2_cp;
	(void)out_d3_cp;
	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Rational Bezier prepare (homogeneous CPs)
 * =================================================================== */

qaws_status qaws_rational_bezier_prepare_2d(
	const qaws_vec2* control_points,
	const qaws_scalar* weights,
	unsigned int degree,
	qaws_scalar* out_weighted_cp)
{
	unsigned int i;
	unsigned int n_pts = degree + 1;

	if (!control_points || !weights || !out_weighted_cp)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (i = 0; i < n_pts; i++) {
		qaws_scalar w = weights[i];
		out_weighted_cp[i * 3 + 0] = control_points[i].x * w;
		out_weighted_cp[i * 3 + 1] = control_points[i].y * w;
		out_weighted_cp[i * 3 + 2] = w;
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_rational_bezier_prepare_3d(
	const qaws_vec3* control_points,
	const qaws_scalar* weights,
	unsigned int degree,
	qaws_scalar* out_weighted_cp)
{
	unsigned int i;
	unsigned int n_pts = degree + 1;

	if (!control_points || !weights || !out_weighted_cp)
		return QAWS_STATUS_INVALID_ARGUMENT;

	for (i = 0; i < n_pts; i++) {
		qaws_scalar w = weights[i];
		out_weighted_cp[i * 4 + 0] = control_points[i].x * w;
		out_weighted_cp[i * 4 + 1] = control_points[i].y * w;
		out_weighted_cp[i * 4 + 2] = control_points[i].z * w;
		out_weighted_cp[i * 4 + 3] = w;
	}

	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Arc prepare
 * =================================================================== */

qaws_status qaws_arc_prepare_2d(
	qaws_vec2 center,
	qaws_scalar radius_x,
	qaws_scalar radius_y,
	qaws_scalar start_angle,
	qaws_scalar end_angle,
	qaws_vec2* out_axis_u,
	qaws_vec2* out_axis_v,
	qaws_scalar* out_theta_start,
	qaws_scalar* out_theta_range)
{
	(void)center;
	if (!out_axis_u || !out_axis_v || !out_theta_start || !out_theta_range)
		return QAWS_STATUS_INVALID_ARGUMENT;

	out_axis_u->x = radius_x;
	out_axis_u->y = QAWS_ZERO;
	out_axis_v->x = QAWS_ZERO;
	out_axis_v->y = radius_y;
	*out_theta_start = start_angle;
	*out_theta_range = end_angle - start_angle;

	return QAWS_STATUS_OK;
}

qaws_status qaws_arc_prepare_3d(
	qaws_vec3 center,
	qaws_vec3 normal,
	qaws_scalar radius_x,
	qaws_scalar radius_y,
	qaws_scalar start_angle,
	qaws_scalar end_angle,
	qaws_vec3* out_axis_u,
	qaws_vec3* out_axis_v,
	qaws_scalar* out_theta_start,
	qaws_scalar* out_theta_range)
{
	qaws_scalar nx, ny, nz, len;
	qaws_vec3 up, u, v;

	(void)center;
	if (!out_axis_u || !out_axis_v || !out_theta_start || !out_theta_range)
		return QAWS_STATUS_INVALID_ARGUMENT;

	nx = normal.x; ny = normal.y; nz = normal.z;
	len = QAWS_SQRT(nx * nx + ny * ny + nz * nz);
	if (len < QAWS_LITERAL(1e-15))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nx /= len; ny /= len; nz /= len;

	/* Build orthonormal basis: pick up vector not parallel to normal */
	if (QAWS_FABS(nx) < QAWS_LITERAL(0.9)) {
		up.x = QAWS_ONE; up.y = QAWS_ZERO; up.z = QAWS_ZERO;
	} else {
		up.x = QAWS_ZERO; up.y = QAWS_ONE; up.z = QAWS_ZERO;
	}

	/* u = normalize(up cross normal) */
	u.x = up.y * nz - up.z * ny;
	u.y = up.z * nx - up.x * nz;
	u.z = up.x * ny - up.y * nx;
	len = QAWS_SQRT(u.x * u.x + u.y * u.y + u.z * u.z);
	u.x /= len; u.y /= len; u.z /= len;

	/* v = normal cross u */
	v.x = ny * u.z - nz * u.y;
	v.y = nz * u.x - nx * u.z;
	v.z = nx * u.y - ny * u.x;

	out_axis_u->x = u.x * radius_x;
	out_axis_u->y = u.y * radius_x;
	out_axis_u->z = u.z * radius_x;

	out_axis_v->x = v.x * radius_y;
	out_axis_v->y = v.y * radius_y;
	out_axis_v->z = v.z * radius_y;

	*out_theta_start = start_angle;
	*out_theta_range = end_angle - start_angle;

	return QAWS_STATUS_OK;
}
