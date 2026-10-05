#include "qaws_hermite.h"
#include "qaws_curve.h"
#include "qaws_prepare.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include "core/qaws_cubic_poly_core.h"
#include <stdlib.h>
#include <string.h>
#include "qaws_platform.h"

/* -------------------------------------------------------------------------- */
/*  Vtable functions                                                          */
/* -------------------------------------------------------------------------- */

static qaws_status hermite_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
	unsigned int const dim_count = 2;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_vec2 va, vb, vc, vd;
	qaws_eval_2d core;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->span_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[span_index * dim_count * 4 + 1 * 4];

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

static qaws_status hermite_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
	unsigned int const dim_count = 3;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar const *cz;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar az, bz, cz_coeff, dz_val;
	qaws_vec3 va, vb, vc, vd;
	qaws_eval_3d core;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->span_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[span_index * dim_count * 4 + 1 * 4];
	cz = &impl->span_coeffs[span_index * dim_count * 4 + 2 * 4];

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

static void hermite_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	qaws_hermite_impl *hi = (qaws_hermite_impl *)impl;
	if (hi) {
		qaws_internal_dealloc(allocator, hi->span_coeffs);
		qaws_internal_dealloc(allocator, hi->points);
		qaws_internal_dealloc(allocator, hi->tangents);
		qaws_internal_dealloc(allocator, hi);
	}
}

static int hermite_is_closed(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int hermite_is_periodic(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int hermite_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity hermite_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C1;
}

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                    */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*  Differential rules: linear in points and derivatives                      */
/* -------------------------------------------------------------------------- */

#define HERMITE_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT | QAWS_CAP_LINEAR)

static unsigned int hermite_describe_fields(qaws_curve const *curve, qaws_field_desc *out, unsigned int capacity)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_POINTS, (qaws_value_type)curve->dimension,
			impl->point_count, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, HERMITE_DIFF_CAPS);
	if (capacity >= 2)
		out[1] = qaws_internal_field_desc(QAWS_FIELD_DERIVATIVES, (qaws_value_type)curve->dimension,
			impl->point_count, QAWS_DOMAIN_DIRECTION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, HERMITE_DIFF_CAPS);
	return 2;
}

static qaws_status hermite_primal_field(qaws_curve const *curve, qaws_diff_field field,
	qaws_scalar const **out_data, unsigned int *out_count, unsigned int *out_components)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
	if (field == QAWS_FIELD_POINTS)
		*out_data = impl->points;
	else if (field == QAWS_FIELD_DERIVATIVES)
		*out_data = impl->tangents;
	else
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = impl->point_count;
	*out_components = (unsigned int)curve->dimension;
	return QAWS_STATUS_OK;
}

/* Cubic Hermite basis h00, h01 (points) and h10, h11 (derivatives) on a
   unit span: global and local parameters coincide. */
static qaws_status hermite_linear_support(qaws_curve const *curve, unsigned int span_index,
	qaws_scalar local_t, unsigned int order, qaws_local_support *out)
{
	qaws_scalar t = local_t, t2 = t * t, t3 = t2 * t;
	unsigned int r;
	(void)curve;
	if (order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	out->kind = QAWS_SUPPORT_LOCAL;
	out->range_count = 2;
	out->ranges[0].field = QAWS_FIELD_POINTS;
	out->ranges[0].first = span_index;
	out->ranges[0].count = 2;
	out->ranges[1].field = QAWS_FIELD_DERIVATIVES;
	out->ranges[1].first = span_index;
	out->ranges[1].count = 2;
	out->weight_field = QAWS_FIELD_NONE;
	out->has_weights = 1;
	out->order = order;

	for (r = 0; r <= order; r++) {
		qaws_scalar h00, h10, h01, h11;
		switch (r) {
		case 0:
			h00 = 2 * t3 - 3 * t2 + 1; h10 = t3 - 2 * t2 + t;
			h01 = -2 * t3 + 3 * t2;    h11 = t3 - t2;
			break;
		case 1:
			h00 = 6 * t2 - 6 * t;      h10 = 3 * t2 - 4 * t + 1;
			h01 = -6 * t2 + 6 * t;     h11 = 3 * t2 - 2 * t;
			break;
		case 2:
			h00 = 12 * t - 6;          h10 = 6 * t - 4;
			h01 = -12 * t + 6;         h11 = 6 * t - 2;
			break;
		case 3:
			h00 = 12; h10 = 6; h01 = -12; h11 = 6;
			break;
		default:
			h00 = h10 = h01 = h11 = 0;
			break;
		}
		out->weights[0][r][0] = h00;
		out->weights[0][r][1] = h01;
		out->weights[1][r][0] = h10;
		out->weights[1][r][1] = h11;
	}
	return QAWS_STATUS_OK;
}

static qaws_status hermite_rebuild(qaws_curve const *curve, qaws_diff_views const *values, qaws_curve **out_curve)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
	qaws_hermite_desc d;
	qaws_scalar const *pts = NULL, *ders = NULL;
	qaws_scalar *owned_pts = NULL, *owned_ders = NULL;
	unsigned int dim = (unsigned int)curve->dimension;
	qaws_status st = qaws_internal_field_override(values, QAWS_FIELD_POINTS, impl->points, impl->point_count, dim, &pts, &owned_pts);
	if (st != QAWS_STATUS_OK)
		return st;
	st = qaws_internal_field_override(values, QAWS_FIELD_DERIVATIVES, impl->tangents, impl->point_count, dim, &ders, &owned_ders);
	if (st == QAWS_STATUS_OK) {
		d.dimension = curve->dimension;
		d.degree = 3;
		d.points = pts;
		d.derivatives = ders;
		d.point_count = impl->point_count;
		d.derivative_count = impl->point_count;
		st = qaws_curve_create_hermite(&d, out_curve);
	}
	qaws_internal_dealloc(NULL, owned_pts);
	qaws_internal_dealloc(NULL, owned_ders);
	return st;
}

static qaws_curve_diff_vtable const hermite_diff_vtable = {
	HERMITE_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	hermite_describe_fields,
	hermite_primal_field,
	hermite_linear_support,
	NULL,
	NULL,
	NULL,
	hermite_rebuild
};

static qaws_curve_vtable const hermite_vtable = {
	hermite_eval_span_2d,
	hermite_eval_span_3d,
	hermite_destroy_impl,
	hermite_is_closed,
	hermite_is_periodic,
	hermite_is_rational,
	hermite_get_continuity,
	&hermite_diff_vtable
};

/* -------------------------------------------------------------------------- */
/*  Creation function                                                         */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_hermite_ex(
	qaws_hermite_desc const *desc,
	qaws_allocator const *allocator,
	qaws_curve **out_curve)
{
	qaws_status status;
	qaws_curve *curve;
	qaws_hermite_impl *impl;
	unsigned int dim_count;
	unsigned int span_count;
	size_t pts_size;
	qaws_range range;
	unsigned int i;

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	status = qaws_internal_validate_dimension(desc->dimension);
	if (status != QAWS_STATUS_OK)
		return status;

	if (desc->degree != 3)
		return QAWS_STATUS_INVALID_DEGREE;

	if (desc->point_count < 2)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (desc->derivative_count != desc->point_count)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (!desc->points || !desc->derivatives)
		return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	span_count = desc->point_count - 1;

	range.min_value = QAWS_ZERO;
	range.max_value = (qaws_scalar)span_count;

	curve = qaws_internal_curve_alloc_ex(
		QAWS_CURVE_KIND_HERMITE,
		desc->dimension,
		desc->degree,
		span_count,
		range,
		&hermite_vtable,
		allocator);
	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i <= span_count; i++) {
		curve->span_boundaries[i] = (qaws_scalar)i;
	}

	impl = (qaws_hermite_impl *)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_hermite_impl));
	if (!impl) {
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	pts_size = (size_t)desc->point_count * (size_t)dim_count * sizeof(qaws_scalar);

	impl->points = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)pts_size);
	if (!impl->points) {
		qaws_internal_dealloc(allocator, impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->points, desc->points, pts_size);

	impl->tangents = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)pts_size);
	if (!impl->tangents) {
		qaws_internal_dealloc(allocator, impl->points);
		qaws_internal_dealloc(allocator, impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->tangents, desc->derivatives, pts_size);

	impl->point_count = desc->point_count;

	/* Precompute per-span polynomial coefficients via prepare functions */
	{
		unsigned int s;
		impl->span_coeffs = (qaws_scalar *)qaws_internal_alloc(allocator,
			(unsigned long)(sizeof(qaws_scalar) * (size_t)span_count * (size_t)dim_count * 4));
		if (!impl->span_coeffs) {
			qaws_internal_dealloc(allocator, impl->tangents);
			qaws_internal_dealloc(allocator, impl->points);
			qaws_internal_dealloc(allocator, impl);
			qaws_internal_curve_free(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		if (dim_count == 2) {
			qaws_vec2 *tmp_a, *tmp_b, *tmp_c, *tmp_d;
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec2);
			tmp_a = (qaws_vec2 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->span_coeffs);
				qaws_internal_dealloc(allocator, impl->tangents);
				qaws_internal_dealloc(allocator, impl->points);
				qaws_internal_dealloc(allocator, impl);
				qaws_internal_curve_free(curve);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_hermite_prepare_2d(
				(const qaws_vec2 *)impl->points,
				(const qaws_vec2 *)impl->tangents,
				desc->point_count,
				tmp_a, tmp_b, tmp_c, tmp_d);

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
			size_t vec_size = (size_t)span_count * sizeof(qaws_vec3);
			tmp_a = (qaws_vec3 *)qaws_internal_alloc(allocator, (unsigned long)(vec_size * 4));
			if (!tmp_a) {
				qaws_internal_dealloc(allocator, impl->span_coeffs);
				qaws_internal_dealloc(allocator, impl->tangents);
				qaws_internal_dealloc(allocator, impl->points);
				qaws_internal_dealloc(allocator, impl);
				qaws_internal_curve_free(curve);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			tmp_b = tmp_a + span_count;
			tmp_c = tmp_b + span_count;
			tmp_d = tmp_c + span_count;

			qaws_hermite_prepare_3d(
				(const qaws_vec3 *)impl->points,
				(const qaws_vec3 *)impl->tangents,
				desc->point_count,
				tmp_a, tmp_b, tmp_c, tmp_d);

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

qaws_status qaws_curve_create_hermite(
	qaws_hermite_desc const *desc,
	qaws_curve **out_curve)
{
	return qaws_curve_create_hermite_ex(desc, NULL, out_curve);
}
