#include "qaws_clothoid.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Impl struct
 * ------------------------------------------------------------------------- */

typedef struct qaws_clothoid_impl
{
	qaws_scalar origin_x, origin_y;
	qaws_scalar start_angle;
	qaws_scalar kappa_0;             /* start curvature */
	qaws_scalar kappa_1;             /* end curvature */
	qaws_scalar length;              /* total arc length */
} qaws_clothoid_impl;

/* ---------------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------------- */

/* Compute theta(s) = start_angle + kappa_0 * s + (kappa_1 - kappa_0) * s^2 / (2 * L) */
static qaws_scalar clothoid_theta(qaws_clothoid_impl const *impl, qaws_scalar s)
{
	qaws_scalar kd = (impl->kappa_1 - impl->kappa_0) / impl->length;
	return impl->start_angle + impl->kappa_0 * s + kd * s * s * (qaws_scalar)0.5;
}

/* Compute position at arc-length s using Simpson's rule integration */
static void clothoid_position(qaws_clothoid_impl const *impl, qaws_scalar s,
	qaws_scalar *out_x, qaws_scalar *out_y)
{
	unsigned int n = 64; /* Simpson panels (must be even) */
	qaws_scalar h = s / (qaws_scalar)n;
	qaws_scalar sum_x = 0, sum_y = 0;
	qaws_scalar kd = (impl->kappa_1 - impl->kappa_0) / impl->length;
	unsigned int i;

	for (i = 0; i <= n; i++)
	{
		qaws_scalar u = (qaws_scalar)i * h;
		qaws_scalar theta = impl->start_angle + impl->kappa_0 * u + kd * u * u * (qaws_scalar)0.5;
		qaws_scalar w = (i == 0 || i == n) ? (qaws_scalar)1 : (i % 2 == 1) ? (qaws_scalar)4 : (qaws_scalar)2;
		sum_x += w * (qaws_scalar)cos((double)theta);
		sum_y += w * (qaws_scalar)sin((double)theta);
	}
	*out_x = impl->origin_x + sum_x * h / (qaws_scalar)3;
	*out_y = impl->origin_y + sum_y * h / (qaws_scalar)3;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_2d
 * ------------------------------------------------------------------------- */

static qaws_status clothoid_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_clothoid_impl *impl = (qaws_clothoid_impl *)curve->impl;
	qaws_scalar L = impl->length;
	qaws_scalar s = local_t * L; /* local_t is [0,1]; map to arc-length */
	qaws_scalar theta;
	qaws_scalar cos_theta, sin_theta;
	qaws_scalar kappa, kappa_prime;

	(void)span_index;

	memset(out_result, 0, sizeof(*out_result));

	/* Position: numerical integration */
	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		if (s <= (qaws_scalar)0) {
			out_result->position.x = impl->origin_x;
			out_result->position.y = impl->origin_y;
		} else {
			clothoid_position(impl, s, &out_result->position.x, &out_result->position.y);
		}
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* Precompute theta and trig values for derivatives */
	theta = clothoid_theta(impl, s);
	cos_theta = (qaws_scalar)cos((double)theta);
	sin_theta = (qaws_scalar)sin((double)theta);

	/* Derivatives with respect to the curve parameter, the arc length s
	   (local_t only locates the point): D1 = (cos, sin) */
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		out_result->d1.x = cos_theta;
		out_result->d1.y = sin_theta;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2 = kappa (-sin, cos) */
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		kappa = impl->kappa_0 + (impl->kappa_1 - impl->kappa_0) * s / L;
		out_result->d2.x = kappa * (-sin_theta);
		out_result->d2.y = kappa * cos_theta;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3 = kappa' (-sin, cos) + kappa^2 (-cos, -sin) */
	if (eval_flags & QAWS_EVAL_FLAG_D3)
	{
		kappa = impl->kappa_0 + (impl->kappa_1 - impl->kappa_0) * s / L;
		kappa_prime = (impl->kappa_1 - impl->kappa_0) / L;
		out_result->d3.x = kappa_prime * (-sin_theta) + kappa * kappa * (-cos_theta);
		out_result->d3.y = kappa_prime * cos_theta + kappa * kappa * (-sin_theta);
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_3d
 * ------------------------------------------------------------------------- */

static qaws_status clothoid_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	(void)curve;
	(void)span_index;
	(void)local_t;
	(void)eval_flags;
	(void)out_result;
	return QAWS_STATUS_INVALID_DIMENSION;
}

/* ---------------------------------------------------------------------------
 * Vtable: destroy
 * ------------------------------------------------------------------------- */

static void clothoid_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	if (impl)
		qaws_internal_dealloc(allocator, impl);
}

/* ---------------------------------------------------------------------------
 * Vtable: property queries
 * ------------------------------------------------------------------------- */

static int clothoid_is_closed(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int clothoid_is_periodic(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int clothoid_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity clothoid_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C3;
}

/* ---------------------------------------------------------------------------
 * Vtable definition
 * ------------------------------------------------------------------------- */

static qaws_curve_vtable const clothoid_vtable = {
	clothoid_eval_span_2d,
	clothoid_eval_span_3d,
	clothoid_destroy_impl,
	clothoid_is_closed,
	clothoid_is_periodic,
	clothoid_is_rational,
	clothoid_get_continuity,
	NULL /* diff */
};

/* ---------------------------------------------------------------------------
 * Creation
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_create_clothoid(
	qaws_clothoid_desc const *desc,
	qaws_curve **out_curve)
{
	qaws_range parameter_range;
	qaws_curve *curve;
	qaws_clothoid_impl *impl;

	if (!desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	/* Length must be positive */
	if (desc->length <= (qaws_scalar)0)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Parameter range: arc-length [0, L] */
	parameter_range.min_value = (qaws_scalar)0;
	parameter_range.max_value = desc->length;

	/* Allocate curve: single span, degree 3 convention, 2D only */
	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_CLOTHOID,
		QAWS_DIMENSION_2D,
		3,
		1,
		parameter_range,
		&clothoid_vtable);

	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Set span boundaries: [0, L] */
	curve->span_boundaries[0] = (qaws_scalar)0;
	curve->span_boundaries[1] = desc->length;

	/* Allocate impl */
	impl = (qaws_clothoid_impl *)calloc(1, sizeof(qaws_clothoid_impl));
	if (!impl)
	{
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->origin_x = desc->origin_x;
	impl->origin_y = desc->origin_y;
	impl->start_angle = desc->start_angle;
	impl->kappa_0 = desc->start_curvature;
	impl->kappa_1 = desc->end_curvature;
	impl->length = desc->length;

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}
