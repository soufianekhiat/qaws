#include "qaws_nurbs.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include "internal/qaws_internal_diff.h"
#include <stdlib.h>
#include <string.h>
#include "qaws_platform.h"
#include "core/qaws_nurbs_eval_core.h"

/* ---------------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------------- */

static int nurbs_compute_span_boundaries(
	qaws_scalar const *knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int control_point_count,
	qaws_allocator const *allocator,
	qaws_scalar **out_boundaries,
	unsigned int *out_span_count)
{
	unsigned int max_spans = knot_count;
	qaws_scalar *tmp = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)((size_t)max_spans * sizeof(qaws_scalar)));
	unsigned int count = 0;
	qaws_scalar prev;
	unsigned int i;

	(void)knot_count;

	if (!tmp)
		return -1;

	prev = knots[degree];
	tmp[count++] = prev;

	for (i = degree + 1; i <= control_point_count; ++i)
	{
		if (knots[i] > prev)
		{
			tmp[count++] = knots[i];
			prev = knots[i];
		}
	}

	*out_boundaries = tmp;
	*out_span_count = count - 1;
	return 0;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_2d
 * ------------------------------------------------------------------------- */

static qaws_status nurbs_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	unsigned int const dim_count = 2;
	qaws_nurbs_impl *impl = (qaws_nurbs_impl *)curve->impl;
	unsigned int degree = curve->degree;
	unsigned int order = degree + 1;
	qaws_scalar span_start;
	qaws_scalar span_end;
	qaws_scalar t;
	unsigned int knot_span;
	unsigned int max_deriv;
	qaws_scalar ders[4 * 64]; /* (max_deriv+1) * order, generous upper bound */
	qaws_scalar W_deriv[4] = {0, 0, 0, 0};
	qaws_scalar A_deriv[4][2] = {{0}};
	unsigned int k, j, idx, d;
	qaws_scalar Nkj, w, Nw;
	qaws_scalar inv_w0;
	qaws_scalar C[2];
	qaws_scalar C1[2] = {0, 0};
	qaws_scalar C2[2] = {0, 0};

	memset(out_result, 0, sizeof(*out_result));

	/* Map local_t [0,1] to global parameter t */
	span_start = curve->span_boundaries[span_index];
	span_end = curve->span_boundaries[span_index + 1];
	t = span_start + local_t * (span_end - span_start);

	/* Find knot span */
	knot_span = qaws_internal_find_knot_span(
		impl->knots, impl->knot_count, degree, impl->control_point_count, t);

	/* Determine max derivative order needed */
	max_deriv = 0;
	if (eval_flags & QAWS_EVAL_FLAG_D1) max_deriv = 1;
	if (eval_flags & QAWS_EVAL_FLAG_D2) max_deriv = 2;
	if (eval_flags & QAWS_EVAL_FLAG_D3) max_deriv = 3;

	/* Compute basis function derivatives */
	qaws_internal_bspline_basis_derivs(
		impl->knots, impl->knot_count, degree, knot_span, t, max_deriv, ders);

	/* Compute A^(k) and W^(k) for k = 0..max_deriv */
	for (k = 0; k <= max_deriv; ++k)
	{
		for (j = 0; j < order; ++j)
		{
			idx = knot_span - degree + j;
			Nkj = ders[k * order + j];
			w = impl->weights[idx];
			Nw = Nkj * w;
			W_deriv[k] += Nw;
			for (d = 0; d < dim_count; ++d)
				A_deriv[k][d] += Nw * impl->control_points[idx * dim_count + d];
		}
	}

	/* Check for degenerate weight */
	if (W_deriv[0] < QAWS_LITERAL(1e-15) && W_deriv[0] > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;

	inv_w0 = QAWS_ONE / W_deriv[0];

	/* Position: C = A[0] / W[0] */
	for (d = 0; d < dim_count; ++d)
		C[d] = A_deriv[0][d] * inv_w0;

	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position.x = C[0];
	out_result->position.y = C[1];

	/* First derivative: C' = (A' - C * W') / W[0] */
	if (max_deriv >= 1)
	{
		for (d = 0; d < dim_count; ++d)
			C1[d] = (A_deriv[1][d] - C[d] * W_deriv[1]) * inv_w0;
		out_result->d1.x = C1[0];
		out_result->d1.y = C1[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* Second derivative: C'' = (A'' - 2*C'*W' - C*W'') / W[0] */
	if (max_deriv >= 2)
	{
		for (d = 0; d < dim_count; ++d)
			C2[d] = (A_deriv[2][d] - QAWS_LITERAL(2.0) * C1[d] * W_deriv[1] - C[d] * W_deriv[2]) * inv_w0;
		out_result->d2.x = C2[0];
		out_result->d2.y = C2[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* Third derivative: C''' = (A''' - 3*C''*W' - 3*C'*W'' - C*W''') / W[0] */
	if (max_deriv >= 3)
	{
		out_result->d3.x = (A_deriv[3][0]
			- QAWS_LITERAL(3.0) * C2[0] * W_deriv[1]
			- QAWS_LITERAL(3.0) * C1[0] * W_deriv[2]
			- C[0] * W_deriv[3]) * inv_w0;
		out_result->d3.y = (A_deriv[3][1]
			- QAWS_LITERAL(3.0) * C2[1] * W_deriv[1]
			- QAWS_LITERAL(3.0) * C1[1] * W_deriv[2]
			- C[1] * W_deriv[3]) * inv_w0;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_3d
 * ------------------------------------------------------------------------- */

static qaws_status nurbs_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	unsigned int const dim_count = 3;
	qaws_nurbs_impl *impl = (qaws_nurbs_impl *)curve->impl;
	unsigned int degree = curve->degree;
	unsigned int order = degree + 1;
	qaws_scalar span_start;
	qaws_scalar span_end;
	qaws_scalar t;
	unsigned int knot_span;
	unsigned int max_deriv;
	qaws_scalar ders[4 * 64];
	qaws_scalar W_deriv[4] = {0, 0, 0, 0};
	qaws_scalar A_deriv[4][3] = {{0}};
	unsigned int k, j, idx, d;
	qaws_scalar Nkj, w, Nw;
	qaws_scalar inv_w0;
	qaws_scalar C[3];
	qaws_scalar C1[3] = {0, 0, 0};
	qaws_scalar C2[3] = {0, 0, 0};

	memset(out_result, 0, sizeof(*out_result));

	/* Map local_t [0,1] to global parameter t */
	span_start = curve->span_boundaries[span_index];
	span_end = curve->span_boundaries[span_index + 1];
	t = span_start + local_t * (span_end - span_start);

	/* Find knot span */
	knot_span = qaws_internal_find_knot_span(
		impl->knots, impl->knot_count, degree, impl->control_point_count, t);

	/* Determine max derivative order needed */
	max_deriv = 0;
	if (eval_flags & QAWS_EVAL_FLAG_D1) max_deriv = 1;
	if (eval_flags & QAWS_EVAL_FLAG_D2) max_deriv = 2;
	if (eval_flags & QAWS_EVAL_FLAG_D3) max_deriv = 3;

	/* Compute basis function derivatives */
	qaws_internal_bspline_basis_derivs(
		impl->knots, impl->knot_count, degree, knot_span, t, max_deriv, ders);

	/* Compute A^(k) and W^(k) for k = 0..max_deriv */
	for (k = 0; k <= max_deriv; ++k)
	{
		for (j = 0; j < order; ++j)
		{
			idx = knot_span - degree + j;
			Nkj = ders[k * order + j];
			w = impl->weights[idx];
			Nw = Nkj * w;
			W_deriv[k] += Nw;
			for (d = 0; d < dim_count; ++d)
				A_deriv[k][d] += Nw * impl->control_points[idx * dim_count + d];
		}
	}

	/* Check for degenerate weight */
	if (W_deriv[0] < QAWS_LITERAL(1e-15) && W_deriv[0] > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;

	inv_w0 = QAWS_ONE / W_deriv[0];

	/* Position: C = A[0] / W[0] */
	for (d = 0; d < dim_count; ++d)
		C[d] = A_deriv[0][d] * inv_w0;

	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position.x = C[0];
	out_result->position.y = C[1];
	out_result->position.z = C[2];

	/* First derivative: C' = (A' - C * W') / W[0] */
	if (max_deriv >= 1)
	{
		for (d = 0; d < dim_count; ++d)
			C1[d] = (A_deriv[1][d] - C[d] * W_deriv[1]) * inv_w0;
		out_result->d1.x = C1[0];
		out_result->d1.y = C1[1];
		out_result->d1.z = C1[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* Second derivative: C'' = (A'' - 2*C'*W' - C*W'') / W[0] */
	if (max_deriv >= 2)
	{
		for (d = 0; d < dim_count; ++d)
			C2[d] = (A_deriv[2][d] - QAWS_LITERAL(2.0) * C1[d] * W_deriv[1] - C[d] * W_deriv[2]) * inv_w0;
		out_result->d2.x = C2[0];
		out_result->d2.y = C2[1];
		out_result->d2.z = C2[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* Third derivative: C''' = (A''' - 3*C''*W' - 3*C'*W'' - C*W''') / W[0] */
	if (max_deriv >= 3)
	{
		out_result->d3.x = (A_deriv[3][0]
			- QAWS_LITERAL(3.0) * C2[0] * W_deriv[1]
			- QAWS_LITERAL(3.0) * C1[0] * W_deriv[2]
			- C[0] * W_deriv[3]) * inv_w0;
		out_result->d3.y = (A_deriv[3][1]
			- QAWS_LITERAL(3.0) * C2[1] * W_deriv[1]
			- QAWS_LITERAL(3.0) * C1[1] * W_deriv[2]
			- C[1] * W_deriv[3]) * inv_w0;
		out_result->d3.z = (A_deriv[3][2]
			- QAWS_LITERAL(3.0) * C2[2] * W_deriv[1]
			- QAWS_LITERAL(3.0) * C1[2] * W_deriv[2]
			- C[2] * W_deriv[3]) * inv_w0;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: destroy
 * ------------------------------------------------------------------------- */

static void nurbs_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	qaws_nurbs_impl *ni = (qaws_nurbs_impl *)impl;
	if (ni)
	{
		qaws_internal_dealloc(allocator, ni->control_points);
		qaws_internal_dealloc(allocator, ni->knots);
		qaws_internal_dealloc(allocator, ni->weights);
		qaws_internal_dealloc(allocator, ni);
	}
}

/* ---------------------------------------------------------------------------
 * Vtable: property queries
 * ------------------------------------------------------------------------- */

static int nurbs_is_closed(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int nurbs_is_periodic(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int nurbs_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 1;
}

static qaws_continuity nurbs_get_continuity(qaws_curve const *curve)
{
	switch (curve->degree)
	{
	case 1:  return QAWS_CONTINUITY_C0;
	case 2:  return QAWS_CONTINUITY_C1;
	case 3:  return QAWS_CONTINUITY_C2;
	default: return QAWS_CONTINUITY_C2;
	}
}

/* ---------------------------------------------------------------------------
 * Vtable definition
 * ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
 * Differential rules: rational in control points and weights
 * ------------------------------------------------------------------------- */

#define NURBS_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT)

static unsigned int nurbs_describe_fields(qaws_curve const *curve, qaws_field_desc *out, unsigned int capacity)
{
	qaws_nurbs_impl const *impl = (qaws_nurbs_impl const *)curve->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CONTROL_POINTS, (qaws_value_type)curve->dimension,
			impl->control_point_count, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, NURBS_DIFF_CAPS);
	if (capacity >= 2)
		out[1] = qaws_internal_field_desc(QAWS_FIELD_WEIGHTS, QAWS_VALUE_SCALAR,
			impl->control_point_count, QAWS_DOMAIN_WEIGHT, QAWS_CONSTRAINT_POSITIVE,
			QAWS_DIFF_SMOOTH, NURBS_DIFF_CAPS);
	if (capacity >= 3)
		out[2] = qaws_internal_field_desc(QAWS_FIELD_KNOTS, QAWS_VALUE_SCALAR,
			impl->knot_count, QAWS_DOMAIN_PARAMETRIC, QAWS_CONSTRAINT_MONOTONIC,
			QAWS_DIFF_UNSUPPORTED, 0u);
	return 3;
}

static qaws_status nurbs_primal_field(qaws_curve const *curve, qaws_diff_field field,
	qaws_scalar const **out_data, unsigned int *out_count, unsigned int *out_components)
{
	qaws_nurbs_impl const *impl = (qaws_nurbs_impl const *)curve->impl;
	switch (field)
	{
	case QAWS_FIELD_CONTROL_POINTS:
		*out_data = impl->control_points;
		*out_count = impl->control_point_count;
		*out_components = (unsigned int)curve->dimension;
		return QAWS_STATUS_OK;
	case QAWS_FIELD_WEIGHTS:
		*out_data = impl->weights;
		*out_count = impl->control_point_count;
		*out_components = 1;
		return QAWS_STATUS_OK;
	case QAWS_FIELD_KNOTS:
		*out_data = impl->knots;
		*out_count = impl->knot_count;
		*out_components = 1;
		return QAWS_STATUS_OK;
	default:
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
}

/* Homogeneous B-spline basis: the engine applies the quotient rule. */
static qaws_status nurbs_linear_support(qaws_curve const *curve, unsigned int span_index,
	qaws_scalar local_t, unsigned int order, qaws_local_support *out)
{
	qaws_nurbs_impl const *impl = (qaws_nurbs_impl const *)curve->impl;
	qaws_scalar ders[(QAWS_DIFF_MAX_ORDER + 1) * QAWS_DIFF_MAX_SUPPORT];
	unsigned int degree = curve->degree;
	unsigned int knot_span, r, j;
	qaws_scalar t;

	if (degree + 1 > QAWS_DIFF_MAX_SUPPORT || order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	/* Same parameter mapping and span search as nurbs_eval_span_*. */
	t = curve->span_boundaries[span_index]
	  + local_t * (curve->span_boundaries[span_index + 1]
	             - curve->span_boundaries[span_index]);
	knot_span = qaws_internal_find_knot_span(
		impl->knots, impl->knot_count, degree, impl->control_point_count, t);
	qaws_internal_bspline_basis_derivs_any(
		impl->knots, impl->knot_count, degree, knot_span, t, order, ders);

	out->kind = QAWS_SUPPORT_LOCAL;
	out->range_count = 1;
	out->ranges[0].field = QAWS_FIELD_CONTROL_POINTS;
	out->ranges[0].first = knot_span - degree;
	out->ranges[0].count = degree + 1;
	out->weight_field = QAWS_FIELD_WEIGHTS;
	out->has_weights = 1;
	out->order = order;
	for (r = 0; r <= order; r++)
		for (j = 0; j <= degree; j++)
			out->weights[0][r][j] = ders[r * (degree + 1) + j];
	return QAWS_STATUS_OK;
}

static qaws_status nurbs_rebuild(qaws_curve const *curve, qaws_diff_views const *values, qaws_curve **out_curve)
{
	qaws_nurbs_impl const *impl = (qaws_nurbs_impl const *)curve->impl;
	qaws_nurbs_desc d;
	qaws_scalar const *cps = NULL, *w = NULL, *knots = NULL;
	qaws_scalar *owned_cps = NULL, *owned_w = NULL, *owned_knots = NULL;
	qaws_status st = qaws_internal_field_override(values, QAWS_FIELD_CONTROL_POINTS, impl->control_points,
		impl->control_point_count, (unsigned int)curve->dimension, &cps, &owned_cps);
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_field_override(values, QAWS_FIELD_WEIGHTS, impl->weights, impl->control_point_count, 1, &w, &owned_w);
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_field_override(values, QAWS_FIELD_KNOTS, impl->knots, impl->knot_count, 1, &knots, &owned_knots);
	if (st == QAWS_STATUS_OK) {
		memset(&d, 0, sizeof(d));
		d.dimension = curve->dimension;
		d.degree = curve->degree;
		d.control_points = cps;
		d.control_point_count = impl->control_point_count;
		d.knots = knots;
		d.knot_count = impl->knot_count;
		d.weights = w;
		d.weight_count = impl->control_point_count;
		st = qaws_curve_create_nurbs(&d, out_curve);
	}
	qaws_internal_dealloc(NULL, owned_cps);
	qaws_internal_dealloc(NULL, owned_w);
	qaws_internal_dealloc(NULL, owned_knots);
	return st;
}

static qaws_curve_diff_vtable const nurbs_diff_vtable = {
	NURBS_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	nurbs_describe_fields,
	nurbs_primal_field,
	nurbs_linear_support,
	NULL,
	NULL,
	NULL,
	nurbs_rebuild
};

static qaws_curve_vtable const nurbs_vtable = {
	nurbs_eval_span_2d,
	nurbs_eval_span_3d,
	nurbs_destroy_impl,
	nurbs_is_closed,
	nurbs_is_periodic,
	nurbs_is_rational,
	nurbs_get_continuity,
	&nurbs_diff_vtable
};

/* ---------------------------------------------------------------------------
 * Creation
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_create_nurbs_ex(
	qaws_nurbs_desc const *desc,
	qaws_allocator const *allocator,
	qaws_curve **out_curve)
{
	unsigned int dim_count;
	unsigned int degree;
	unsigned int control_point_count;
	unsigned int knot_count;
	qaws_status status;
	qaws_scalar *span_boundaries = NULL;
	unsigned int span_count = 0;
	qaws_range parameter_range;
	qaws_curve *curve;
	qaws_nurbs_impl *impl;
	size_t cp_size;
	size_t knot_size;
	size_t weight_size;

	if (!desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	/* Validate dimension */
	status = qaws_internal_validate_dimension(desc->dimension);
	if (status != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;

	dim_count = (unsigned int)desc->dimension;
	degree = desc->degree;

	if (degree < 1)
		return QAWS_STATUS_INVALID_DEGREE;

	control_point_count = desc->control_point_count;
	if (control_point_count < degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (!desc->control_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (!desc->knots)
		return QAWS_STATUS_INVALID_ARGUMENT;

	knot_count = desc->knot_count;

	/* Validate knot vector */
	status = qaws_internal_validate_knot_vector(
		desc->knots, knot_count, degree, control_point_count);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Validate weights */
	if (!desc->weights)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (desc->weight_count != control_point_count)
		return QAWS_STATUS_INVALID_WEIGHT_COUNT;

	status = qaws_internal_validate_weights(desc->weights, desc->weight_count);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Compute span boundaries */
	if (nurbs_compute_span_boundaries(desc->knots, knot_count, degree,
		control_point_count, allocator, &span_boundaries, &span_count) != 0)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	if (span_count < 1)
	{
		qaws_internal_dealloc(allocator, span_boundaries);
		return QAWS_STATUS_INVALID_KNOT_VECTOR;
	}

	/* Parameter range */
	parameter_range.min_value = desc->knots[degree];
	parameter_range.max_value = desc->knots[control_point_count];

	/* Allocate curve */
	curve = qaws_internal_curve_alloc_ex(
		QAWS_CURVE_KIND_NURBS,
		desc->dimension,
		degree,
		span_count,
		parameter_range,
		&nurbs_vtable,
		allocator);

	if (!curve)
	{
		qaws_internal_dealloc(allocator, span_boundaries);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy span boundaries into the curve */
	memcpy(curve->span_boundaries, span_boundaries, (size_t)(span_count + 1) * sizeof(qaws_scalar));
	qaws_internal_dealloc(allocator, span_boundaries);

	/* Allocate impl */
	impl = (qaws_nurbs_impl *)qaws_internal_alloc(allocator,
		(unsigned long)sizeof(qaws_nurbs_impl));
	if (!impl)
	{
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(impl, 0, sizeof(qaws_nurbs_impl));

	/* Copy control points */
	cp_size = (size_t)control_point_count * (size_t)dim_count * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)cp_size);
	if (!impl->control_points)
	{
		qaws_internal_dealloc(allocator, impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points, cp_size);
	impl->control_point_count = control_point_count;

	/* Copy knots */
	knot_size = (size_t)knot_count * sizeof(qaws_scalar);
	impl->knots = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)knot_size);
	if (!impl->knots)
	{
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->knots, desc->knots, knot_size);
	impl->knot_count = knot_count;

	/* Copy weights */
	weight_size = (size_t)control_point_count * sizeof(qaws_scalar);
	impl->weights = (qaws_scalar *)qaws_internal_alloc(allocator,
		(unsigned long)weight_size);
	if (!impl->weights)
	{
		qaws_internal_dealloc(allocator, impl->knots);
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->weights, desc->weights, weight_size);
	impl->weight_count = control_point_count;

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_nurbs(
	qaws_nurbs_desc const *desc,
	qaws_curve **out_curve)
{
	return qaws_curve_create_nurbs_ex(desc, NULL, out_curve);
}
