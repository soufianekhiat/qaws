#include "qaws_polynomial.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_validation.h"
#include "internal/qaws_internal_diff.h"
#include <stdlib.h>
#include <string.h>
#include "qaws_platform.h"
#include "core/qaws_horner_core.h"

typedef struct qaws_polynomial_impl
{
	qaws_scalar* coefficients;         /* (degree+1) * dim_count scalars */
	unsigned int coefficient_count;
} qaws_polynomial_impl;

/* -------------------------------------------------------------------------- */
/*  Horner evaluation helpers                                                  */
/* -------------------------------------------------------------------------- */

/* Evaluate polynomial position using Horner's method for one component.
 * coeffs: pointer to the first scalar of the desired component (stride = dim_count)
 * degree: polynomial degree
 * dim_count: number of scalars per coefficient (2 or 3)
 * t: parameter value */
static qaws_scalar poly_horner(
	qaws_scalar const* coeffs, unsigned int degree,
	unsigned int dim_count, qaws_scalar t)
{
	int i;
	qaws_scalar result = coeffs[degree * dim_count];
	for (i = (int)degree - 1; i >= 0; i--)
		result = result * t + coeffs[(unsigned int)i * dim_count];
	return result;
}

/* Evaluate first derivative using Horner's method for one component.
 * P'(t) = c1 + 2*c2*t + 3*c3*t^2 + ... = sum_{i=1}^{degree} i * c_i * t^{i-1}
 * Horner form: c1 + t*(2*c2 + t*(3*c3 + ...)) */
static qaws_scalar poly_horner_d1(
	qaws_scalar const* coeffs, unsigned int degree,
	unsigned int dim_count, qaws_scalar t)
{
	int i;
	qaws_scalar result;
	if (degree < 1) return QAWS_ZERO;
	result = (qaws_scalar)degree * coeffs[degree * dim_count];
	for (i = (int)degree - 1; i >= 1; i--)
		result = result * t + (qaws_scalar)i * coeffs[(unsigned int)i * dim_count];
	return result;
}

/* Evaluate second derivative using Horner's method for one component.
 * P''(t) = 2*c2 + 6*c3*t + 12*c4*t^2 + ... = sum_{i=2}^{degree} i*(i-1) * c_i * t^{i-2}
 * Horner form: 2*c2 + t*(6*c3 + t*(12*c4 + ...)) */
static qaws_scalar poly_horner_d2(
	qaws_scalar const* coeffs, unsigned int degree,
	unsigned int dim_count, qaws_scalar t)
{
	int i;
	qaws_scalar result;
	if (degree < 2) return QAWS_ZERO;
	result = (qaws_scalar)(degree * (degree - 1)) * coeffs[degree * dim_count];
	for (i = (int)degree - 1; i >= 2; i--)
		result = result * t + (qaws_scalar)(i * (i - 1)) * coeffs[(unsigned int)i * dim_count];
	return result;
}

/* Evaluate third derivative using Horner's method for one component.
 * P'''(t) = 6*c3 + 24*c4*t + ... = sum_{i=3}^{degree} i*(i-1)*(i-2) * c_i * t^{i-3}
 * Horner form: 6*c3 + t*(24*c4 + ...) */
static qaws_scalar poly_horner_d3(
	qaws_scalar const* coeffs, unsigned int degree,
	unsigned int dim_count, qaws_scalar t)
{
	int i;
	qaws_scalar result;
	if (degree < 3) return QAWS_ZERO;
	result = (qaws_scalar)(degree * (degree - 1) * (degree - 2)) * coeffs[degree * dim_count];
	for (i = (int)degree - 1; i >= 3; i--)
		result = result * t + (qaws_scalar)(i * (i - 1) * (i - 2)) * coeffs[(unsigned int)i * dim_count];
	return result;
}

/* -------------------------------------------------------------------------- */
/*  Vtable functions                                                           */
/* -------------------------------------------------------------------------- */

static qaws_status polynomial_eval_span_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_polynomial_impl const* impl = (qaws_polynomial_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* coeffs = impl->coefficients;
	unsigned int const dim_count = 2;
	qaws_scalar t_min = curve->parameter_range.min_value;
	qaws_scalar t_max = curve->parameter_range.max_value;
	qaws_scalar t = t_min + local_t * (t_max - t_min);
	qaws_eval_2d core;
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	/* Core handles position, D1, D2 */
	core = qaws_horner_eval_2d((qaws_scalar *)coeffs, (int)degree, t, (int)eval_flags);

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

	/* D3 is not supported by the core; compute manually */
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		out_result->d3.x = poly_horner_d3(coeffs + 0, degree, dim_count, t);
		out_result->d3.y = poly_horner_d3(coeffs + 1, degree, dim_count, t);
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static qaws_status polynomial_eval_span_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_polynomial_impl const* impl = (qaws_polynomial_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* coeffs = impl->coefficients;
	unsigned int const dim_count = 3;
	qaws_scalar t_min = curve->parameter_range.min_value;
	qaws_scalar t_max = curve->parameter_range.max_value;
	qaws_scalar t = t_min + local_t * (t_max - t_min);
	qaws_eval_3d core;
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	/* Core handles position, D1, D2 */
	core = qaws_horner_eval_3d((qaws_scalar *)coeffs, (int)degree, t, (int)eval_flags);

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

	/* D3 is not supported by the core; compute manually */
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		out_result->d3.x = poly_horner_d3(coeffs + 0, degree, dim_count, t);
		out_result->d3.y = poly_horner_d3(coeffs + 1, degree, dim_count, t);
		out_result->d3.z = poly_horner_d3(coeffs + 2, degree, dim_count, t);
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static void polynomial_destroy_impl(void* impl, qaws_allocator const* allocator)
{
	qaws_polynomial_impl* pi = (qaws_polynomial_impl*)impl;
	if (pi) { qaws_internal_dealloc(allocator, pi->coefficients); qaws_internal_dealloc(allocator, pi); }
}

static int polynomial_is_closed(qaws_curve const* c)    { (void)c; return 0; }
static int polynomial_is_periodic(qaws_curve const* c)   { (void)c; return 0; }
static int polynomial_is_rational(qaws_curve const* c)   { (void)c; return 0; }
static qaws_continuity polynomial_get_continuity(qaws_curve const* c) { (void)c; return QAWS_CONTINUITY_C3; }

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                     */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*  Differential rules: linear in the coefficients                            */
/* -------------------------------------------------------------------------- */

#define POLYNOMIAL_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT | QAWS_CAP_LINEAR)

static unsigned int polynomial_describe_fields(qaws_curve const* curve, qaws_field_desc* out, unsigned int capacity)
{
	qaws_polynomial_impl const* impl = (qaws_polynomial_impl const*)curve->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_COEFFICIENTS, (qaws_value_type)curve->dimension,
			impl->coefficient_count, QAWS_DOMAIN_GENERIC, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, POLYNOMIAL_DIFF_CAPS);
	return 1;
}

static qaws_status polynomial_primal_field(qaws_curve const* curve, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_polynomial_impl const* impl = (qaws_polynomial_impl const*)curve->impl;
	if (field != QAWS_FIELD_COEFFICIENTS)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_data = impl->coefficients;
	*out_count = impl->coefficient_count;
	*out_components = (unsigned int)curve->dimension;
	return QAWS_STATUS_OK;
}

/* d^r/dt^r t^i = i!/(i-r)! t^(i-r) at the global parameter. */
static qaws_status polynomial_linear_support(qaws_curve const* curve, unsigned int span_index,
	qaws_scalar local_t, unsigned int order, qaws_local_support* out)
{
	qaws_polynomial_impl const* impl = (qaws_polynomial_impl const*)curve->impl;
	qaws_scalar t_min = curve->parameter_range.min_value;
	qaws_scalar t_max = curve->parameter_range.max_value;
	qaws_scalar t = t_min + local_t * (t_max - t_min);
	unsigned int n = impl->coefficient_count, r, i, q;
	(void)span_index;
	if (n > QAWS_DIFF_MAX_SUPPORT || order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	out->kind = QAWS_SUPPORT_GLOBAL;
	out->range_count = 1;
	out->ranges[0].field = QAWS_FIELD_COEFFICIENTS;
	out->ranges[0].first = 0;
	out->ranges[0].count = n;
	out->weight_field = QAWS_FIELD_NONE;
	out->has_weights = 1;
	out->order = order;
	for (r = 0; r <= order; r++) {
		for (i = 0; i < n; i++) {
			qaws_scalar w = QAWS_ZERO;
			if (i >= r) {
				w = QAWS_ONE;
				for (q = 0; q < r; q++)
					w *= (qaws_scalar)(i - q);
				for (q = 0; q < i - r; q++)
					w *= t;
			}
			out->weights[0][r][i] = w;
		}
	}
	return QAWS_STATUS_OK;
}

static qaws_curve_diff_vtable const polynomial_diff_vtable = {
	POLYNOMIAL_DIFF_CAPS,
	QAWS_DIFF_SMOOTH,
	polynomial_describe_fields,
	polynomial_primal_field,
	polynomial_linear_support
};

static qaws_curve_vtable const polynomial_vtable = {
	polynomial_eval_span_2d,
	polynomial_eval_span_3d,
	polynomial_destroy_impl,
	polynomial_is_closed,
	polynomial_is_periodic,
	polynomial_is_rational,
	polynomial_get_continuity,
	&polynomial_diff_vtable
};

/* -------------------------------------------------------------------------- */
/*  Creation function                                                          */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_polynomial(
	qaws_polynomial_desc const* desc,
	qaws_curve** out_curve)
{
	unsigned int dim_count;
	size_t coeff_size;
	qaws_range range;
	qaws_curve* curve;
	qaws_polynomial_impl* impl;

	if (!desc || !out_curve) return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;

	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (desc->degree < 1) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->coefficient_count != desc->degree + 1) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->coefficients) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->t_max <= desc->t_min) return QAWS_STATUS_INVALID_PARAMETER_RANGE;

	dim_count = (unsigned int)desc->dimension;
	range.min_value = desc->t_min;
	range.max_value = desc->t_max;

	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_POLYNOMIAL, desc->dimension, desc->degree,
		1, range, &polynomial_vtable);
	if (!curve) return QAWS_STATUS_ALLOCATION_FAILURE;
	curve->span_boundaries[0] = desc->t_min;
	curve->span_boundaries[1] = desc->t_max;

	impl = (qaws_polynomial_impl*)malloc(sizeof(qaws_polynomial_impl));
	if (!impl) { qaws_internal_curve_free(curve); return QAWS_STATUS_ALLOCATION_FAILURE; }

	coeff_size = (size_t)dim_count * desc->coefficient_count * sizeof(qaws_scalar);
	impl->coefficients = (qaws_scalar*)malloc(coeff_size);
	if (!impl->coefficients) { free(impl); qaws_internal_curve_free(curve); return QAWS_STATUS_ALLOCATION_FAILURE; }
	memcpy(impl->coefficients, desc->coefficients, coeff_size);
	impl->coefficient_count = desc->coefficient_count;

	curve->impl = impl;
	*out_curve = curve;
	return QAWS_STATUS_OK;
}
