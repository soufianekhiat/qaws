#include "qaws_catmull_rom.h"
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

/* -------------------------------------------------------------------------- */
/*  Helpers                                                                   */
/* -------------------------------------------------------------------------- */

static qaws_scalar compute_alpha(qaws_parameterization param)
{
	switch (param) {
	case QAWS_PARAMETERIZATION_UNIFORM:
		return QAWS_ZERO;
	case QAWS_PARAMETERIZATION_CHORDAL:
		return QAWS_ONE;
	case QAWS_PARAMETERIZATION_CENTRIPETAL:
		return QAWS_LITERAL(0.5);
	case QAWS_PARAMETERIZATION_DEFAULT:
	default:
		return QAWS_LITERAL(0.5); /* default is centripetal */
	}
}

static void compute_knot_params(
	qaws_scalar const *control_points,
	unsigned int control_point_count,
	unsigned int dim_count,
	qaws_scalar alpha,
	qaws_scalar *knot_params)
{
	unsigned int i, d;

	knot_params[0] = QAWS_ZERO;

	for (i = 1; i < control_point_count; i++) {
		qaws_scalar dist_sq = QAWS_ZERO;
		for (d = 0; d < dim_count; d++) {
			qaws_scalar diff = control_points[i * dim_count + d]
			                 - control_points[(i - 1) * dim_count + d];
			dist_sq += diff * diff;
		}

		if (alpha == QAWS_ZERO) {
			/* Uniform: each interval = 1 */
			knot_params[i] = knot_params[i - 1] + QAWS_ONE;
		} else {
			qaws_scalar dist = QAWS_SQRT(dist_sq);
			qaws_scalar interval = QAWS_POW(dist, alpha);
			knot_params[i] = knot_params[i - 1] + interval;
		}
	}
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (2D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status catmull_rom_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int const dim_count = 2;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_vec2 va, vb, vc, vd;
	qaws_eval_2d core;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];

	va.x = ax; va.y = ay;
	vb.x = bx; vb.y = by;
	vc.x = cx_coeff; vc.y = cy_coeff;
	vd.x = dx_val; vd.y = dy_val;

	core = qaws_cubic_eval_2d(va, vb, vc, vd, local_t, eval_flags);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position = core.position;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1 = core.d1;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2 = core.d2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = QAWS_LITERAL(6.0) * ax;
		out_result->d3.y = QAWS_LITERAL(6.0) * ay;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (3D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status catmull_rom_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int const dim_count = 3;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar const *cz;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar az, bz, cz_coeff, dz_val;
	qaws_vec3 va, vb, vc, vd;
	qaws_eval_3d core;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];
	cz = &impl->segment_coeffs[span_index * dim_count * 4 + 2 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];
	az = cz[0]; bz = cz[1]; cz_coeff = cz[2]; dz_val = cz[3];

	va.x = ax; va.y = ay; va.z = az;
	vb.x = bx; vb.y = by; vb.z = bz;
	vc.x = cx_coeff; vc.y = cy_coeff; vc.z = cz_coeff;
	vd.x = dx_val; vd.y = dy_val; vd.z = dz_val;

	core = qaws_cubic_eval_3d(va, vb, vc, vd, local_t, eval_flags);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position = core.position;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1 = core.d1;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2 = core.d2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = QAWS_LITERAL(6.0) * ax;
		out_result->d3.y = QAWS_LITERAL(6.0) * ay;
		out_result->d3.z = QAWS_LITERAL(6.0) * az;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: lifecycle and query functions                                     */
/* -------------------------------------------------------------------------- */

static void catmull_rom_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	qaws_catmull_rom_impl *cr = (qaws_catmull_rom_impl *)impl;
	if (cr) {
		qaws_internal_dealloc(allocator, cr->control_points);
		qaws_internal_dealloc(allocator, cr->knot_params);
		qaws_internal_dealloc(allocator, cr->segment_coeffs);
		qaws_internal_dealloc(allocator, cr);
	}
}

static int catmull_rom_is_closed(qaws_curve const *curve)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	return impl->closed;
}

static int catmull_rom_is_periodic(qaws_curve const *curve)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	return impl->closed;
}

static int catmull_rom_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity catmull_rom_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C1;
}

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                    */
/* -------------------------------------------------------------------------- */

static qaws_curve_vtable const catmull_rom_vtable = {
	catmull_rom_eval_span_2d,
	catmull_rom_eval_span_3d,
	catmull_rom_destroy_impl,
	catmull_rom_is_closed,
	catmull_rom_is_periodic,
	catmull_rom_is_rational,
	catmull_rom_get_continuity,
	NULL /* diff */
};

/* -------------------------------------------------------------------------- */
/*  Creation                                                                  */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_catmull_rom_ex(
	qaws_catmull_rom_desc const *desc,
	qaws_allocator const *allocator,
	qaws_curve **out_curve)
{
	qaws_catmull_rom_impl *impl = NULL;
	qaws_curve *curve = NULL;
	unsigned int dim_count;
	unsigned int span_count;
	unsigned int N;
	unsigned int i, s;
	qaws_scalar alpha;
	qaws_range range;

	/* --- Validation ---------------------------------------------------- */

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	if (desc->dimension != QAWS_DIMENSION_2D &&
	    desc->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (!desc->control_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	N = desc->control_point_count;

	if (desc->closed) {
		if (N < 3)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = N;
	} else {
		if (N < 4)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = N - 3;
	}

	/* --- Allocate implementation --------------------------------------- */

	impl = (qaws_catmull_rom_impl *)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_catmull_rom_impl));
	if (!impl)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(impl, 0, sizeof(qaws_catmull_rom_impl));

	impl->control_point_count = N;
	impl->parameterization = desc->parameterization;
	impl->closed = desc->closed ? 1 : 0;

	/* Copy control points */
	impl->control_points = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)(sizeof(qaws_scalar) * (size_t)(N * dim_count)));
	if (!impl->control_points) {
		qaws_internal_dealloc(allocator, impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points,
	       sizeof(qaws_scalar) * (size_t)(N * dim_count));

	/* Compute knot parameters */
	impl->knot_params = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)(sizeof(qaws_scalar) * (size_t)N));
	if (!impl->knot_params) {
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	alpha = compute_alpha(impl->parameterization);
	compute_knot_params(impl->control_points, N, dim_count,
	                    alpha, impl->knot_params);

	/* Compute segment coefficients via prepare functions */
	{
		size_t coeffs_size = (size_t)(span_count * dim_count * 4) * sizeof(qaws_scalar);
		impl->segment_coeffs = (qaws_scalar *)qaws_internal_alloc(allocator,
			(unsigned long)coeffs_size);
		if (!impl->segment_coeffs) {
			qaws_internal_dealloc(allocator, impl->knot_params);
			qaws_internal_dealloc(allocator, impl->control_points);
			qaws_internal_dealloc(allocator, impl);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		if (dim_count == 2) {
			qaws_vec2 *tmp_a, *tmp_b, *tmp_c, *tmp_d;
			unsigned int prep_span_count = 0;
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec2);
			tmp_a = (qaws_vec2 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->segment_coeffs);
				qaws_internal_dealloc(allocator, impl->knot_params);
				qaws_internal_dealloc(allocator, impl->control_points);
				qaws_internal_dealloc(allocator, impl);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_catmull_rom_prepare_2d(
				(const qaws_vec2 *)impl->control_points, N,
				impl->parameterization, impl->closed,
				tmp_a, tmp_b, tmp_c, tmp_d, &prep_span_count);

			for (s = 0; s < span_count; s++) {
				impl->segment_coeffs[s * 8 + 0] = tmp_a[s].x;
				impl->segment_coeffs[s * 8 + 1] = tmp_b[s].x;
				impl->segment_coeffs[s * 8 + 2] = tmp_c[s].x;
				impl->segment_coeffs[s * 8 + 3] = tmp_d[s].x;
				impl->segment_coeffs[s * 8 + 4] = tmp_a[s].y;
				impl->segment_coeffs[s * 8 + 5] = tmp_b[s].y;
				impl->segment_coeffs[s * 8 + 6] = tmp_c[s].y;
				impl->segment_coeffs[s * 8 + 7] = tmp_d[s].y;
			}
			qaws_internal_dealloc(allocator, tmp_a);
		} else {
			qaws_vec3 *tmp_a, *tmp_b, *tmp_c, *tmp_d;
			unsigned int prep_span_count = 0;
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec3);
			tmp_a = (qaws_vec3 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->segment_coeffs);
				qaws_internal_dealloc(allocator, impl->knot_params);
				qaws_internal_dealloc(allocator, impl->control_points);
				qaws_internal_dealloc(allocator, impl);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_catmull_rom_prepare_3d(
				(const qaws_vec3 *)impl->control_points, N,
				impl->parameterization, impl->closed,
				tmp_a, tmp_b, tmp_c, tmp_d, &prep_span_count);

			for (s = 0; s < span_count; s++) {
				impl->segment_coeffs[s * 12 + 0]  = tmp_a[s].x;
				impl->segment_coeffs[s * 12 + 1]  = tmp_b[s].x;
				impl->segment_coeffs[s * 12 + 2]  = tmp_c[s].x;
				impl->segment_coeffs[s * 12 + 3]  = tmp_d[s].x;
				impl->segment_coeffs[s * 12 + 4]  = tmp_a[s].y;
				impl->segment_coeffs[s * 12 + 5]  = tmp_b[s].y;
				impl->segment_coeffs[s * 12 + 6]  = tmp_c[s].y;
				impl->segment_coeffs[s * 12 + 7]  = tmp_d[s].y;
				impl->segment_coeffs[s * 12 + 8]  = tmp_a[s].z;
				impl->segment_coeffs[s * 12 + 9]  = tmp_b[s].z;
				impl->segment_coeffs[s * 12 + 10] = tmp_c[s].z;
				impl->segment_coeffs[s * 12 + 11] = tmp_d[s].z;
			}
			qaws_internal_dealloc(allocator, tmp_a);
		}
	}

	/* --- Allocate curve ------------------------------------------------ */

	range.min_value = QAWS_ZERO;
	range.max_value = (qaws_scalar)span_count;

	curve = qaws_internal_curve_alloc_ex(
		QAWS_CURVE_KIND_CATMULL_ROM,
		desc->dimension,
		3,            /* degree */
		span_count,
		range,
		&catmull_rom_vtable,
		allocator);
	if (!curve) {
		qaws_internal_dealloc(allocator, impl->segment_coeffs);
		qaws_internal_dealloc(allocator, impl->knot_params);
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Set span boundaries: uniform [0..span_count] */
	for (i = 0; i <= span_count; i++) {
		curve->span_boundaries[i] = (qaws_scalar)i;
	}

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_catmull_rom(
	qaws_catmull_rom_desc const *desc,
	qaws_curve **out_curve)
{
	return qaws_curve_create_catmull_rom_ex(desc, NULL, out_curve);
}
