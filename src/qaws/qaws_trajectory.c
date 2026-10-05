#include "qaws_trajectory.h"
#include "qaws_curve.h"
#include "qaws_prepare.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include "core/qaws_cubic_poly_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* ---------------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------------- */

static void trajectory_compute_velocities(
	qaws_scalar const *positions,
	qaws_scalar const *times,
	unsigned int key_count,
	unsigned int dim_count,
	qaws_scalar *out_velocities)
{
	unsigned int i, d;
	qaws_scalar dt, inv_dt;
	unsigned int last;

	if (key_count < 2)
		return;

	/* Endpoint 0: forward difference */
	dt = times[1] - times[0];
	if (dt < QAWS_LITERAL(1e-15)) dt = QAWS_LITERAL(1e-15);
	inv_dt = QAWS_ONE / dt;
	for (d = 0; d < dim_count; ++d)
	{
		out_velocities[0 * dim_count + d] =
			(positions[1 * dim_count + d] - positions[0 * dim_count + d]) * inv_dt;
	}

	/* Interior points: Catmull-Rom finite differences */
	for (i = 1; i < key_count - 1; ++i)
	{
		dt = times[i + 1] - times[i - 1];
		if (dt < QAWS_LITERAL(1e-15)) dt = QAWS_LITERAL(1e-15);
		inv_dt = QAWS_ONE / dt;
		for (d = 0; d < dim_count; ++d)
		{
			out_velocities[i * dim_count + d] =
				(positions[(i + 1) * dim_count + d] - positions[(i - 1) * dim_count + d]) * inv_dt;
		}
	}

	/* Endpoint N-1: backward difference */
	last = key_count - 1;
	dt = times[last] - times[last - 1];
	if (dt < QAWS_LITERAL(1e-15)) dt = QAWS_LITERAL(1e-15);
	inv_dt = QAWS_ONE / dt;
	for (d = 0; d < dim_count; ++d)
	{
		out_velocities[last * dim_count + d] =
			(positions[last * dim_count + d] - positions[(last - 1) * dim_count + d]) * inv_dt;
	}
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_2d
 * ------------------------------------------------------------------------- */

static qaws_status trajectory_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	unsigned int const dim_count = 2;
	unsigned int s = span_index;
	qaws_scalar dt;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar inv_dt, inv_dt2, inv_dt3;
	qaws_vec2 va, vb, vc, vd;
	int core_flags;
	qaws_eval_2d core;

	memset(out_result, 0, sizeof(*out_result));

	dt = impl->key_times[s + 1] - impl->key_times[s];
	if (dt < QAWS_LITERAL(1e-10))
		return QAWS_STATUS_DEGENERATE_CURVE;

	cx = &impl->span_coeffs[s * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[s * dim_count * 4 + 1 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];

	va.x = ax; va.y = ay;
	vb.x = bx; vb.y = by;
	vc.x = cx_coeff; vc.y = cy_coeff;
	vd.x = dx_val; vd.y = dy_val;

	/* Request d1 from core whenever any derivative is needed;
	   request d2 from core whenever D2 or D3 is needed. */
	core_flags = 0;
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		core_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		core_flags |= QAWS_EVAL_FLAG_D2;

	core = qaws_cubic_eval_2d(va, vb, vc, vd, local_t, core_flags);

	/* Position: P(u) — always computed */
	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position = core.position;

	/* D1: P'(u) / dt */
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) {
		inv_dt = QAWS_ONE / dt;
		out_result->d1.x = core.d1.x * inv_dt;
		out_result->d1.y = core.d1.y * inv_dt;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2: P''(u) / dt^2 */
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) {
		inv_dt2 = QAWS_ONE / (dt * dt);
		out_result->d2.x = core.d2.x * inv_dt2;
		out_result->d2.y = core.d2.y * inv_dt2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3: P'''(u) / dt^3 */
	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		inv_dt3 = QAWS_ONE / (dt * dt * dt);
		out_result->d3.x = (QAWS_LITERAL(6.0) * ax) * inv_dt3;
		out_result->d3.y = (QAWS_LITERAL(6.0) * ay) * inv_dt3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_3d
 * ------------------------------------------------------------------------- */

static qaws_status trajectory_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	unsigned int const dim_count = 3;
	unsigned int s = span_index;
	qaws_scalar dt;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar const *cz;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar az, bz, cz_coeff, dz_val;
	qaws_scalar inv_dt, inv_dt2, inv_dt3;
	qaws_vec3 va, vb, vc, vd;
	int core_flags;
	qaws_eval_3d core;

	memset(out_result, 0, sizeof(*out_result));

	dt = impl->key_times[s + 1] - impl->key_times[s];
	if (dt < QAWS_LITERAL(1e-10))
		return QAWS_STATUS_DEGENERATE_CURVE;

	cx = &impl->span_coeffs[s * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[s * dim_count * 4 + 1 * 4];
	cz = &impl->span_coeffs[s * dim_count * 4 + 2 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];
	az = cz[0]; bz = cz[1]; cz_coeff = cz[2]; dz_val = cz[3];

	va.x = ax; va.y = ay; va.z = az;
	vb.x = bx; vb.y = by; vb.z = bz;
	vc.x = cx_coeff; vc.y = cy_coeff; vc.z = cz_coeff;
	vd.x = dx_val; vd.y = dy_val; vd.z = dz_val;

	/* Request d1 from core whenever any derivative is needed;
	   request d2 from core whenever D2 or D3 is needed. */
	core_flags = 0;
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		core_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		core_flags |= QAWS_EVAL_FLAG_D2;

	core = qaws_cubic_eval_3d(va, vb, vc, vd, local_t, core_flags);

	/* Position: P(u) — always computed */
	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position = core.position;

	/* D1: P'(u) / dt */
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) {
		inv_dt = QAWS_ONE / dt;
		out_result->d1.x = core.d1.x * inv_dt;
		out_result->d1.y = core.d1.y * inv_dt;
		out_result->d1.z = core.d1.z * inv_dt;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2: P''(u) / dt^2 */
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) {
		inv_dt2 = QAWS_ONE / (dt * dt);
		out_result->d2.x = core.d2.x * inv_dt2;
		out_result->d2.y = core.d2.y * inv_dt2;
		out_result->d2.z = core.d2.z * inv_dt2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3: P'''(u) / dt^3 */
	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		inv_dt3 = QAWS_ONE / (dt * dt * dt);
		out_result->d3.x = (QAWS_LITERAL(6.0) * ax) * inv_dt3;
		out_result->d3.y = (QAWS_LITERAL(6.0) * ay) * inv_dt3;
		out_result->d3.z = (QAWS_LITERAL(6.0) * az) * inv_dt3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: destroy
 * ------------------------------------------------------------------------- */

static void trajectory_destroy_impl(void *impl, qaws_allocator const *allocator)
{
	qaws_trajectory_impl *ti = (qaws_trajectory_impl *)impl;
	if (ti)
	{
		qaws_internal_dealloc(allocator, ti->span_coeffs);
		qaws_internal_dealloc(allocator, ti->key_positions);
		qaws_internal_dealloc(allocator, ti->key_times);
		qaws_internal_dealloc(allocator, ti->key_velocities);
		qaws_internal_dealloc(allocator, ti->key_accelerations);
		qaws_internal_dealloc(allocator, ti);
	}
}

/* ---------------------------------------------------------------------------
 * Vtable: property queries
 * ------------------------------------------------------------------------- */

static int trajectory_is_closed(qaws_curve const *curve)
{
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	return impl->closed;
}

static int trajectory_is_periodic(qaws_curve const *curve)
{
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	return impl->closed;
}

static int trajectory_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity trajectory_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C1;
}

/* ---------------------------------------------------------------------------
 * Vtable definition
 * ------------------------------------------------------------------------- */

static qaws_curve_vtable const trajectory_vtable = {
	trajectory_eval_span_2d,
	trajectory_eval_span_3d,
	trajectory_destroy_impl,
	trajectory_is_closed,
	trajectory_is_periodic,
	trajectory_is_rational,
	trajectory_get_continuity,
	NULL /* diff */
};

/* ---------------------------------------------------------------------------
 * Creation
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_create_trajectory_ex(
	qaws_trajectory_desc const *desc,
	qaws_allocator const *allocator,
	qaws_curve **out_curve)
{
	unsigned int dim_count;
	unsigned int key_count;
	unsigned int span_count;
	unsigned int i;
	int user_velocities;
	qaws_range parameter_range;
	qaws_curve *curve;
	qaws_trajectory_impl *impl;
	size_t pos_size;
	size_t time_size;
	size_t vel_size;
	size_t acc_size;

	if (!desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	/* Validate dimension */
	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;

	dim_count = (unsigned int)desc->dimension;

	/* Only degree 3 (cubic Hermite) is supported */
	if (desc->degree != 3)
		return QAWS_STATUS_INVALID_DEGREE;

	key_count = desc->key_count;
	if (key_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (!desc->key_positions)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (!desc->key_times)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (desc->key_time_count != key_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate key_times are strictly increasing */
	for (i = 1; i < key_count; ++i)
	{
		if (desc->key_times[i] <= desc->key_times[i - 1])
			return QAWS_STATUS_INVALID_PARAMETER_RANGE;
	}

	span_count = key_count - 1;

	/* Parameter range */
	parameter_range.min_value = desc->key_times[0];
	parameter_range.max_value = desc->key_times[key_count - 1];

	/* Allocate curve */
	curve = qaws_internal_curve_alloc_ex(
		QAWS_CURVE_KIND_TRAJECTORY,
		desc->dimension,
		desc->degree,
		span_count,
		parameter_range,
		&trajectory_vtable,
		allocator);

	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Set span boundaries from key_times */
	for (i = 0; i < key_count; ++i)
		curve->span_boundaries[i] = desc->key_times[i];

	/* Allocate impl */
	impl = (qaws_trajectory_impl *)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_trajectory_impl));
	if (!impl)
	{
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(impl, 0, sizeof(qaws_trajectory_impl));

	impl->key_count = key_count;
	impl->closed = desc->closed;

	/* Copy key_positions */
	pos_size = (size_t)key_count * (size_t)dim_count * sizeof(qaws_scalar);
	impl->key_positions = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)pos_size);
	if (!impl->key_positions)
	{
		qaws_internal_dealloc(allocator, impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->key_positions, desc->key_positions, pos_size);

	/* Copy key_times */
	time_size = (size_t)key_count * sizeof(qaws_scalar);
	impl->key_times = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)time_size);
	if (!impl->key_times)
	{
		qaws_internal_dealloc(allocator, impl->key_positions);
		qaws_internal_dealloc(allocator, impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->key_times, desc->key_times, time_size);

	/* Handle velocities */
	vel_size = (size_t)key_count * (size_t)dim_count * sizeof(qaws_scalar);
	user_velocities = (desc->key_velocities && desc->key_velocity_count == key_count);

	if (user_velocities)
	{
		/* Copy user-provided velocities */
		impl->key_velocities = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)vel_size);
		if (!impl->key_velocities)
		{
			qaws_internal_dealloc(allocator, impl->key_times);
			qaws_internal_dealloc(allocator, impl->key_positions);
			qaws_internal_dealloc(allocator, impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		memcpy(impl->key_velocities, desc->key_velocities, vel_size);
		impl->key_velocity_count = key_count;
	}
	else
	{
		/* Compute velocities using finite differences */
		impl->key_velocities = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)vel_size);
		if (!impl->key_velocities)
		{
			qaws_internal_dealloc(allocator, impl->key_times);
			qaws_internal_dealloc(allocator, impl->key_positions);
			qaws_internal_dealloc(allocator, impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		trajectory_compute_velocities(
			impl->key_positions, impl->key_times,
			key_count, dim_count, impl->key_velocities);
		impl->key_velocity_count = key_count;
	}

	/* Copy key_accelerations if provided */
	impl->key_accelerations = NULL;
	impl->key_acceleration_count = 0;
	if (desc->key_accelerations && desc->key_acceleration_count > 0)
	{
		acc_size = (size_t)desc->key_acceleration_count * (size_t)dim_count * sizeof(qaws_scalar);
		impl->key_accelerations = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)acc_size);
		if (!impl->key_accelerations)
		{
			qaws_internal_dealloc(allocator, impl->key_velocities);
			qaws_internal_dealloc(allocator, impl->key_times);
			qaws_internal_dealloc(allocator, impl->key_positions);
			qaws_internal_dealloc(allocator, impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		memcpy(impl->key_accelerations, desc->key_accelerations, acc_size);
		impl->key_acceleration_count = desc->key_acceleration_count;
	}

	/* Precompute per-span polynomial coefficients */
	{
		unsigned int s;
		impl->span_coeffs = (qaws_scalar *)qaws_internal_alloc(allocator,
			(unsigned long)(sizeof(qaws_scalar) * (size_t)span_count * (size_t)dim_count * 4));
		if (!impl->span_coeffs) {
			qaws_internal_dealloc(allocator, impl->key_accelerations);
			qaws_internal_dealloc(allocator, impl->key_velocities);
			qaws_internal_dealloc(allocator, impl->key_times);
			qaws_internal_dealloc(allocator, impl->key_positions);
			qaws_internal_dealloc(allocator, impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		if (user_velocities) {
			/* User-provided velocities: compute coefficients inline */
			unsigned int d;
			for (s = 0; s < span_count; s++) {
				qaws_scalar dt = impl->key_times[s + 1] - impl->key_times[s];
				qaws_scalar const *p0 = impl->key_positions + s * dim_count;
				qaws_scalar const *p1 = impl->key_positions + (s + 1) * dim_count;
				qaws_scalar const *v0 = impl->key_velocities + s * dim_count;
				qaws_scalar const *v1 = impl->key_velocities + (s + 1) * dim_count;
				for (d = 0; d < dim_count; d++) {
					qaws_scalar m0 = v0[d] * dt;
					qaws_scalar m1 = v1[d] * dt;
					qaws_scalar a_coeff = QAWS_LITERAL(2.0) * p0[d] + m0 - QAWS_LITERAL(2.0) * p1[d] + m1;
					qaws_scalar b_coeff = -QAWS_LITERAL(3.0) * p0[d] - QAWS_LITERAL(2.0) * m0 + QAWS_LITERAL(3.0) * p1[d] - m1;
					qaws_scalar c_coeff = m0;
					qaws_scalar d_coeff = p0[d];
					impl->span_coeffs[(s * dim_count + d) * 4 + 0] = a_coeff;
					impl->span_coeffs[(s * dim_count + d) * 4 + 1] = b_coeff;
					impl->span_coeffs[(s * dim_count + d) * 4 + 2] = c_coeff;
					impl->span_coeffs[(s * dim_count + d) * 4 + 3] = d_coeff;
				}
			}
		} else if (dim_count == 2) {
			qaws_vec2 *tmp_a, *tmp_b, *tmp_c, *tmp_d;
			unsigned int prep_span_count = 0;
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec2);
			tmp_a = (qaws_vec2 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->span_coeffs);
				qaws_internal_dealloc(allocator, impl->key_accelerations);
				qaws_internal_dealloc(allocator, impl->key_velocities);
				qaws_internal_dealloc(allocator, impl->key_times);
				qaws_internal_dealloc(allocator, impl->key_positions);
				qaws_internal_dealloc(allocator, impl);
				qaws_curve_destroy(curve);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_trajectory_prepare_2d(
				(const qaws_vec2 *)impl->key_positions,
				impl->key_times, key_count, desc->degree,
				impl->closed,
				tmp_a, tmp_b, tmp_c, tmp_d, &prep_span_count);

			for (s = 0; s < span_count; s++) {
				impl->span_coeffs[s * 8 + 0] = tmp_a[s].x;
				impl->span_coeffs[s * 8 + 1] = tmp_b[s].x;
				impl->span_coeffs[s * 8 + 2] = tmp_c[s].x;
				impl->span_coeffs[s * 8 + 3] = tmp_d[s].x;
				impl->span_coeffs[s * 8 + 4] = tmp_a[s].y;
				impl->span_coeffs[s * 8 + 5] = tmp_b[s].y;
				impl->span_coeffs[s * 8 + 6] = tmp_c[s].y;
				impl->span_coeffs[s * 8 + 7] = tmp_d[s].y;
			}
			qaws_internal_dealloc(allocator, tmp_a);
		} else {
			qaws_vec3 *tmp_a, *tmp_b, *tmp_c, *tmp_d;
			unsigned int prep_span_count = 0;
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec3);
			tmp_a = (qaws_vec3 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->span_coeffs);
				qaws_internal_dealloc(allocator, impl->key_accelerations);
				qaws_internal_dealloc(allocator, impl->key_velocities);
				qaws_internal_dealloc(allocator, impl->key_times);
				qaws_internal_dealloc(allocator, impl->key_positions);
				qaws_internal_dealloc(allocator, impl);
				qaws_curve_destroy(curve);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_trajectory_prepare_3d(
				(const qaws_vec3 *)impl->key_positions,
				impl->key_times, key_count, desc->degree,
				impl->closed,
				tmp_a, tmp_b, tmp_c, tmp_d, &prep_span_count);

			for (s = 0; s < span_count; s++) {
				impl->span_coeffs[s * 12 + 0]  = tmp_a[s].x;
				impl->span_coeffs[s * 12 + 1]  = tmp_b[s].x;
				impl->span_coeffs[s * 12 + 2]  = tmp_c[s].x;
				impl->span_coeffs[s * 12 + 3]  = tmp_d[s].x;
				impl->span_coeffs[s * 12 + 4]  = tmp_a[s].y;
				impl->span_coeffs[s * 12 + 5]  = tmp_b[s].y;
				impl->span_coeffs[s * 12 + 6]  = tmp_c[s].y;
				impl->span_coeffs[s * 12 + 7]  = tmp_d[s].y;
				impl->span_coeffs[s * 12 + 8]  = tmp_a[s].z;
				impl->span_coeffs[s * 12 + 9]  = tmp_b[s].z;
				impl->span_coeffs[s * 12 + 10] = tmp_c[s].z;
				impl->span_coeffs[s * 12 + 11] = tmp_d[s].z;
			}
			qaws_internal_dealloc(allocator, tmp_a);
		}
	}

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_trajectory(
	qaws_trajectory_desc const *desc,
	qaws_curve **out_curve)
{
	return qaws_curve_create_trajectory_ex(desc, NULL, out_curve);
}
