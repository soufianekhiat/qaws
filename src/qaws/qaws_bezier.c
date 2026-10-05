#include "qaws_bezier.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_decasteljau_core.h"
#include "qaws_platform.h"
#include <stdlib.h>
#include <string.h>

static qaws_status bezier_eval_span_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* cp = impl->control_points;
	qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 2];
	qaws_eval_2d core;
	unsigned int i, n;
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	n = (degree + 1) * 2;
	for (i = 0; i < n; i++)
		local_cp[i] = cp[i];

	core = qaws_decasteljau_eval_2d(local_cp, (int)degree, local_t, (int)eval_flags);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = core.position.x;
		out_result->position.y = core.position.y;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && degree >= 1) {
		out_result->d1.x = core.d1.x;
		out_result->d1.y = core.d1.y;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && degree >= 2) {
		out_result->d2.x = core.d2.x;
		out_result->d2.y = core.d2.y;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		qaws_scalar d1cp[QAWS_CORE_MAX_POINTS * 2], d2cp[QAWS_CORE_MAX_POINTS * 2], d3cp[QAWS_CORE_MAX_POINTS * 2];
		qaws_vec2 d3val;
		qaws_bezier_deriv_points_2d(local_cp, (int)degree, d1cp);
		qaws_bezier_deriv_points_2d(d1cp, (int)degree - 1, d2cp);
		qaws_bezier_deriv_points_2d(d2cp, (int)degree - 2, d3cp);
		d3val = qaws_decasteljau_2d(d3cp, (int)degree - 3, local_t);
		out_result->d3.x = d3val.x; out_result->d3.y = d3val.y;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}
	return QAWS_STATUS_OK;
}

static qaws_status bezier_eval_span_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* cp = impl->control_points;
	qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3];
	qaws_eval_3d core;
	unsigned int i, n;
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	n = (degree + 1) * 3;
	for (i = 0; i < n; i++)
		local_cp[i] = cp[i];

	core = qaws_decasteljau_eval_3d(local_cp, (int)degree, local_t, (int)eval_flags);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = core.position.x;
		out_result->position.y = core.position.y;
		out_result->position.z = core.position.z;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && degree >= 1) {
		out_result->d1.x = core.d1.x;
		out_result->d1.y = core.d1.y;
		out_result->d1.z = core.d1.z;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && degree >= 2) {
		out_result->d2.x = core.d2.x;
		out_result->d2.y = core.d2.y;
		out_result->d2.z = core.d2.z;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		qaws_scalar d1cp[QAWS_CORE_MAX_POINTS * 3], d2cp[QAWS_CORE_MAX_POINTS * 3], d3cp[QAWS_CORE_MAX_POINTS * 3];
		qaws_vec3 d3val;
		qaws_bezier_deriv_points_3d(local_cp, (int)degree, d1cp);
		qaws_bezier_deriv_points_3d(d1cp, (int)degree - 1, d2cp);
		qaws_bezier_deriv_points_3d(d2cp, (int)degree - 2, d3cp);
		d3val = qaws_decasteljau_3d(d3cp, (int)degree - 3, local_t);
		out_result->d3.x = d3val.x; out_result->d3.y = d3val.y; out_result->d3.z = d3val.z;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}
	return QAWS_STATUS_OK;
}

static void bezier_destroy_impl(void* impl, qaws_allocator const* allocator)
{
	qaws_bezier_impl* bi = (qaws_bezier_impl*)impl;
	if (bi) { qaws_internal_dealloc(allocator, bi->control_points); qaws_internal_dealloc(allocator, bi); }
}

static int bezier_is_closed(qaws_curve const* c)   { (void)c; return 0; }
static int bezier_is_periodic(qaws_curve const* c)  { (void)c; return 0; }
static int bezier_is_rational(qaws_curve const* c)  { (void)c; return 0; }
static qaws_continuity bezier_get_continuity(qaws_curve const* c) { (void)c; return QAWS_CONTINUITY_C3; }

/* -------------------------------------------------------------------------- */
/*  Differential rules: linear in the control points                          */
/* -------------------------------------------------------------------------- */

#define BEZIER_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT | QAWS_CAP_LINEAR)

static unsigned int bezier_describe_fields(qaws_curve const* curve, qaws_field_desc* out, unsigned int capacity)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CONTROL_POINTS, (qaws_value_type)curve->dimension,
			impl->control_point_count, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, BEZIER_DIFF_CAPS);
	return 1;
}

static qaws_status bezier_primal_field(qaws_curve const* curve, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	if (field != QAWS_FIELD_CONTROL_POINTS)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_data = impl->control_points;
	*out_count = impl->control_point_count;
	*out_components = (unsigned int)curve->dimension;
	return QAWS_STATUS_OK;
}

static qaws_status bezier_linear_support(qaws_curve const* curve, unsigned int span_index,
	qaws_scalar local_t, unsigned int order, qaws_local_support* out)
{
	qaws_scalar b[(QAWS_DIFF_MAX_ORDER + 1) * QAWS_DIFF_MAX_SUPPORT];
	unsigned int n = curve->degree, r, j;
	(void)span_index;
	if (n + 1 > QAWS_DIFF_MAX_SUPPORT || order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	qaws_internal_bernstein_derivs(n, local_t, order, b);

	out->kind = QAWS_SUPPORT_LOCAL;
	out->range_count = 1;
	out->ranges[0].field = QAWS_FIELD_CONTROL_POINTS;
	out->ranges[0].first = 0;
	out->ranges[0].count = n + 1;
	out->weight_field = QAWS_FIELD_NONE;
	out->has_weights = 1;
	out->order = order;
	for (r = 0; r <= order; r++)
		for (j = 0; j <= n; j++)
			out->weights[0][r][j] = b[r * (n + 1) + j];
	return QAWS_STATUS_OK;
}

static qaws_curve_diff_vtable const bezier_diff_vtable = {
	BEZIER_DIFF_CAPS,
	QAWS_DIFF_SMOOTH,
	bezier_describe_fields,
	bezier_primal_field,
	bezier_linear_support
};

static qaws_curve_vtable const bezier_vtable = {
	bezier_eval_span_2d, bezier_eval_span_3d, bezier_destroy_impl,
	bezier_is_closed, bezier_is_periodic, bezier_is_rational, bezier_get_continuity,
	&bezier_diff_vtable
};

qaws_status qaws_curve_create_bezier_ex(
	qaws_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve)
{
	unsigned int dim_count;
	size_t cp_size;
	qaws_range range;
	qaws_curve* curve;
	qaws_bezier_impl* impl;

	if (!desc || !out_curve) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (desc->degree < 1) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->control_point_count != desc->degree + 1) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points) return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (unsigned int)desc->dimension;
	range.min_value = 0; range.max_value = 1;

	curve = qaws_internal_curve_alloc_ex(QAWS_CURVE_KIND_BEZIER, desc->dimension, desc->degree, 1, range, &bezier_vtable, allocator);
	if (!curve) return QAWS_STATUS_ALLOCATION_FAILURE;
	curve->span_boundaries[0] = 0; curve->span_boundaries[1] = 1;

	impl = (qaws_bezier_impl*)qaws_internal_alloc(allocator, sizeof(qaws_bezier_impl));
	if (!impl) { qaws_curve_destroy(curve); return QAWS_STATUS_ALLOCATION_FAILURE; }

	cp_size = (size_t)dim_count * desc->control_point_count * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)cp_size);
	if (!impl->control_points) { qaws_internal_dealloc(allocator, impl); qaws_curve_destroy(curve); return QAWS_STATUS_ALLOCATION_FAILURE; }
	memcpy(impl->control_points, desc->control_points, cp_size);
	impl->control_point_count = desc->control_point_count;

	curve->impl = impl;
	*out_curve = curve;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_bezier(
	qaws_bezier_desc const* desc, qaws_curve** out_curve)
{
	return qaws_curve_create_bezier_ex(desc, NULL, out_curve);
}
