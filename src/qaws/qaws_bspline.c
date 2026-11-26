#include "qaws_bspline.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include "internal/qaws_internal_diff.h"
#include <stdlib.h>
#include <string.h>
#include "qaws_platform.h"
#include "qaws_core_types.h"
#include "core/qaws_bspline_eval_core.h"
#include "internal/qaws_internal_span.h"
#include "simd/qaws_simd.h"

/* -------------------------------------------------------------------------- */
/*  Helpers                                                                   */
/* -------------------------------------------------------------------------- */

static qaws_status generate_uniform_knots(
	unsigned int degree,
	unsigned int control_point_count,
	qaws_allocator const *allocator,
	qaws_scalar **out_knots,
	unsigned int *out_knot_count)
{
	unsigned int knot_count = control_point_count + degree + 1;
	unsigned int n_internal = knot_count - 2 * (degree + 1);
	qaws_scalar *knots;
	unsigned int i;
	qaws_scalar max_val = (qaws_scalar)(control_point_count - degree);

	knots = (qaws_scalar *)qaws_internal_alloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)knot_count),
		allocator);
	if (!knots)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* First degree+1 knots = 0 (clamped) */
	for (i = 0; i <= degree; i++)
		knots[i] = QAWS_ZERO;

	/* Internal knots: uniformly spaced */
	for (i = 0; i < n_internal; i++)
		knots[degree + 1 + i] = (qaws_scalar)(i + 1);

	/* Last degree+1 knots = max_val (clamped) */
	for (i = 0; i <= degree; i++)
		knots[knot_count - 1 - i] = max_val;

	*out_knots = knots;
	*out_knot_count = knot_count;
	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (2D)                                                    */
/* -------------------------------------------------------------------------- */

/* the sums of a span: the basis functions and derivatives ders_buf ((max_deriv_order + 1) rows of degree + 1, es
   scalars apart) times the control points of knot span knot_span */
static qaws_status bspline_sum_2d(
	qaws_bspline_impl const *impl,
	unsigned int degree,
	unsigned int knot_span,
	qaws_scalar const *ders_buf,
	unsigned int es,
	unsigned int max_deriv_order,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	unsigned int const dim_count = 2;
	unsigned int const stride = degree + 1;
	unsigned int j;

	/* Evaluate position */
	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(0 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
		}
		out_result->position.x = vx;
		out_result->position.y = vy;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* Evaluate first derivative */
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && max_deriv_order >= 1) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(1 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
		}
		out_result->d1.x = vx;
		out_result->d1.y = vy;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* Evaluate second derivative */
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && max_deriv_order >= 2) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(2 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
		}
		out_result->d2.x = vx;
		out_result->d2.y = vy;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* Evaluate third derivative */
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && max_deriv_order >= 3) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(3 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
		}
		out_result->d3.x = vx;
		out_result->d3.y = vy;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static qaws_status bspline_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_bspline_impl const *impl =
		(qaws_bspline_impl const *)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar t;
	unsigned int knot_span;
	unsigned int max_deriv_order = 0;
	qaws_scalar ders_buf[4 * (QAWS_INTERNAL_BASIS_MAX_DEGREE + 1)]; /* (max_deriv+1) * (degree+1) */

	memset(out_result, 0, sizeof(*out_result));

	/* Map local_t [0,1] to global parameter */
	t = curve->span_boundaries[span_index]
	  + local_t * (curve->span_boundaries[span_index + 1]
	             - curve->span_boundaries[span_index]);

	/* Find knot span */
	knot_span = qaws_internal_find_knot_span(
		impl->knots, impl->knot_count,
		degree, impl->control_point_count, t);

	/* Determine max derivative order needed */
	if (eval_flags & QAWS_EVAL_FLAG_D3)
		max_deriv_order = 3;
	else if (eval_flags & QAWS_EVAL_FLAG_D2)
		max_deriv_order = 2;
	else if (eval_flags & QAWS_EVAL_FLAG_D1)
		max_deriv_order = 1;

	/* Limit to degree */
	if (max_deriv_order > degree)
		max_deriv_order = degree;

	/* Compute basis function derivatives */
	qaws_internal_bspline_basis_derivs(
		ders_buf, impl->knots,
		impl->knot_count, degree, knot_span,
		t, max_deriv_order);

	return bspline_sum_2d(impl, degree, knot_span, ders_buf, 1, max_deriv_order, eval_flags, out_result);
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (3D)                                                    */
/* -------------------------------------------------------------------------- */

/* the sums of a span: the basis functions and derivatives ders_buf ((max_deriv_order + 1) rows of degree + 1, es
   scalars apart) times the control points of knot span knot_span */
static qaws_status bspline_sum_3d(
	qaws_bspline_impl const *impl,
	unsigned int degree,
	unsigned int knot_span,
	qaws_scalar const *ders_buf,
	unsigned int es,
	unsigned int max_deriv_order,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	unsigned int const dim_count = 3;
	unsigned int const stride = degree + 1;
	unsigned int j;

	/* Evaluate position */
	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		qaws_scalar vz = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(0 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
			vz += basis * impl->control_points[cp_idx + 2];
		}
		out_result->position.x = vx;
		out_result->position.y = vy;
		out_result->position.z = vz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* Evaluate first derivative */
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && max_deriv_order >= 1) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		qaws_scalar vz = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(1 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
			vz += basis * impl->control_points[cp_idx + 2];
		}
		out_result->d1.x = vx;
		out_result->d1.y = vy;
		out_result->d1.z = vz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* Evaluate second derivative */
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && max_deriv_order >= 2) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		qaws_scalar vz = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(2 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
			vz += basis * impl->control_points[cp_idx + 2];
		}
		out_result->d2.x = vx;
		out_result->d2.y = vy;
		out_result->d2.z = vz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* Evaluate third derivative */
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && max_deriv_order >= 3) {
		qaws_scalar vx = QAWS_ZERO;
		qaws_scalar vy = QAWS_ZERO;
		qaws_scalar vz = QAWS_ZERO;
		for (j = 0; j <= degree; j++) {
			unsigned int cp_idx = (knot_span - degree + j) * dim_count;
			qaws_scalar basis = ders_buf[(3 * stride + j) * es];
			vx += basis * impl->control_points[cp_idx + 0];
			vy += basis * impl->control_points[cp_idx + 1];
			vz += basis * impl->control_points[cp_idx + 2];
		}
		out_result->d3.x = vx;
		out_result->d3.y = vy;
		out_result->d3.z = vz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static qaws_status bspline_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_bspline_impl const *impl =
		(qaws_bspline_impl const *)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar t;
	unsigned int knot_span;
	unsigned int max_deriv_order = 0;
	qaws_scalar ders_buf[4 * (QAWS_INTERNAL_BASIS_MAX_DEGREE + 1)]; /* (max_deriv+1) * (degree+1) */

	memset(out_result, 0, sizeof(*out_result));

	/* Map local_t [0,1] to global parameter */
	t = curve->span_boundaries[span_index]
	  + local_t * (curve->span_boundaries[span_index + 1]
	             - curve->span_boundaries[span_index]);

	/* Find knot span */
	knot_span = qaws_internal_find_knot_span(
		impl->knots, impl->knot_count,
		degree, impl->control_point_count, t);

	/* Determine max derivative order needed */
	if (eval_flags & QAWS_EVAL_FLAG_D3)
		max_deriv_order = 3;
	else if (eval_flags & QAWS_EVAL_FLAG_D2)
		max_deriv_order = 2;
	else if (eval_flags & QAWS_EVAL_FLAG_D1)
		max_deriv_order = 1;

	if (max_deriv_order > degree)
		max_deriv_order = degree;

	/* Compute basis function derivatives */
	qaws_internal_bspline_basis_derivs(
		ders_buf, impl->knots,
		impl->knot_count, degree, knot_span,
		t, max_deriv_order);

	return bspline_sum_3d(impl, degree, knot_span, ders_buf, 1, max_deriv_order, eval_flags, out_result);
}

/* -------------------------------------------------------------------------- */
/*  Batches: the basis functions a SIMD register's worth of parameters at once */
/* -------------------------------------------------------------------------- */

qaws_status qaws_internal_bspline_eval_batch(
	void *out_results,
	qaws_curve const *curve,
	qaws_scalar const *parameters,
	unsigned int count,
	unsigned int eval_flags)
{
	qaws_bspline_impl const *impl = (qaws_bspline_impl const *)curve->impl;
	unsigned int const degree = curve->degree, lanes = qaws_internal_bspline_basis_lanes(degree);
	unsigned int max_deriv_order = 0, i0, l, span_index = 0, ks = 0;
	int have = 0;
	QAWS_SIMD_ALIGN qaws_scalar ders[4 * 6 * QAWS_INTERNAL_BASIS_LANES_MAX];
	QAWS_SIMD_ALIGN qaws_scalar t[QAWS_INTERNAL_BASIS_LANES_MAX];
	unsigned int knot_span[QAWS_INTERNAL_BASIS_LANES_MAX];
	if (!lanes || curve->kind != QAWS_CURVE_KIND_BSPLINE)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (eval_flags & QAWS_EVAL_FLAG_D3)
		max_deriv_order = 3;
	else if (eval_flags & QAWS_EVAL_FLAG_D2)
		max_deriv_order = 2;
	else if (eval_flags & QAWS_EVAL_FLAG_D1)
		max_deriv_order = 1;
	if (max_deriv_order > degree)
		max_deriv_order = degree;
	for (i0 = 0; i0 < count; i0 += lanes)
	{
		unsigned int const n = count - i0 < lanes ? count - i0 : lanes;
		/* each parameter's span, global parameter and knot span as the one-point path finds them (the last repeated in
		   the lanes past the end) */
		for (l = 0; l < lanes; l++)
		{
			qaws_scalar const x = parameters[i0 + (l < n ? l : n - 1)];
			qaws_scalar const* bd = curve->span_boundaries;
			qaws_scalar local_t;
			/* a parameter strictly inside the last span (past the first boundary, short of the last) is in that span for
			   the search too, its local parameter computed as the search computes it */
			if (have && x > bd[0] && x < bd[curve->span_count] && bd[span_index] <= x && x < bd[span_index + 1])
			{
				qaws_scalar const span_len = bd[span_index + 1] - bd[span_index];
				local_t = span_len < (qaws_scalar)1e-10 ? (qaws_scalar)0 : (x - bd[span_index]) / span_len;
				if (local_t < (qaws_scalar)0) local_t = (qaws_scalar)0;
				else if (local_t > (qaws_scalar)1) local_t = (qaws_scalar)1;
			}
			else
				span_index = qaws_internal_find_span(&local_t, curve, x);
			t[l] = bd[span_index]
			  + local_t * (bd[span_index + 1]
			             - bd[span_index]);
			/* likewise the knot span, strictly inside the domain's ends */
			if (!(have && t[l] > impl->knots[degree] && t[l] < impl->knots[impl->control_point_count]
				&& impl->knots[ks] <= t[l] && t[l] < impl->knots[ks + 1]))
				ks = qaws_internal_find_knot_span(
					impl->knots, impl->knot_count,
					degree, impl->control_point_count, t[l]);
			knot_span[l] = ks;
			have = 1;
		}
		qaws_internal_bspline_basis_derivs_lanes(ders, impl->knots, knot_span, t, degree, max_deriv_order);
		for (l = 0; l < n; l++)
		{
			if (curve->dimension == QAWS_DIMENSION_2D)
			{
				qaws_eval_result_2d *r = (qaws_eval_result_2d *)out_results + i0 + l;
				memset(r, 0, sizeof(*r));
				bspline_sum_2d(impl, degree, knot_span[l], ders + l, QAWS_INTERNAL_BASIS_LANES_MAX, max_deriv_order, eval_flags, r);
			}
			else
			{
				qaws_eval_result_3d *r = (qaws_eval_result_3d *)out_results + i0 + l;
				memset(r, 0, sizeof(*r));
				bspline_sum_3d(impl, degree, knot_span[l], ders + l, QAWS_INTERNAL_BASIS_LANES_MAX, max_deriv_order, eval_flags, r);
			}
		}
	}
	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: lifecycle and query functions                                     */
/* -------------------------------------------------------------------------- */

static void bspline_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	qaws_bspline_impl *bi = (qaws_bspline_impl *)impl;
	if (bi) {
		qaws_internal_dealloc(bi->control_points, allocator);
		qaws_internal_dealloc(bi->knots, allocator);
		qaws_internal_dealloc(bi, allocator);
	}
}

static int bspline_is_closed(qaws_curve const *curve)
{
	return qaws_internal_curve_ends_meet(curve);
}

static int bspline_is_periodic(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int bspline_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity bspline_get_continuity(qaws_curve const *curve)
{
	unsigned int degree = curve->degree;
	if (degree >= 4)
		return QAWS_CONTINUITY_C3;
	if (degree == 3)
		return QAWS_CONTINUITY_C2;
	if (degree == 2)
		return QAWS_CONTINUITY_C1;
	/* degree 1 or 0 */
	return QAWS_CONTINUITY_C0;
}

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                    */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*  Differential rules: linear in the control points                          */
/* -------------------------------------------------------------------------- */

#define BSPLINE_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT | QAWS_CAP_LINEAR)

static unsigned int bspline_describe_fields(qaws_curve const *curve, qaws_field_desc *out, unsigned int capacity)
{
	qaws_bspline_impl const *impl = (qaws_bspline_impl const *)curve->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CONTROL_POINTS, (qaws_value_type)curve->dimension,
			impl->control_point_count, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, BSPLINE_DIFF_CAPS);
	if (capacity >= 2)
		out[1] = qaws_internal_field_desc(QAWS_FIELD_KNOTS, QAWS_VALUE_SCALAR,
			impl->knot_count, QAWS_DOMAIN_PARAMETRIC, QAWS_CONSTRAINT_MONOTONIC,
			QAWS_DIFF_PIECEWISE_SMOOTH, QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2);
	return 2;
}

static qaws_status bspline_primal_field(qaws_curve const *curve, qaws_diff_field field,
	qaws_scalar const **out_data, unsigned int *out_count, unsigned int *out_components)
{
	qaws_bspline_impl const *impl = (qaws_bspline_impl const *)curve->impl;
	if (field == QAWS_FIELD_CONTROL_POINTS) {
		*out_data = impl->control_points;
		*out_count = impl->control_point_count;
		*out_components = (unsigned int)curve->dimension;
		return QAWS_STATUS_OK;
	}
	if (field == QAWS_FIELD_KNOTS) {
		*out_data = impl->knots;
		*out_count = impl->knot_count;
		*out_components = 1;
		return QAWS_STATUS_OK;
	}
	return QAWS_STATUS_INVALID_ARGUMENT;
}

static qaws_status bspline_linear_support(qaws_curve const *curve, unsigned int span_index,
	qaws_scalar local_t, unsigned int order, qaws_local_support *out)
{
	qaws_bspline_impl const *impl = (qaws_bspline_impl const *)curve->impl;
	qaws_scalar ders[(QAWS_DIFF_MAX_ORDER + 1) * QAWS_DIFF_MAX_SUPPORT];
	unsigned int degree = curve->degree;
	unsigned int knot_span, r, j;
	qaws_scalar t;

	if (degree + 1 > QAWS_DIFF_MAX_SUPPORT || order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	/* Same parameter mapping and span search as bspline_eval_span_*. */
	t = curve->span_boundaries[span_index]
	  + local_t * (curve->span_boundaries[span_index + 1]
	             - curve->span_boundaries[span_index]);
	knot_span = qaws_internal_find_knot_span(
		impl->knots, impl->knot_count,
		degree, impl->control_point_count, t);

	qaws_internal_bspline_basis_derivs_any(
		ders, impl->knots, impl->knot_count, degree, knot_span, t, order);

	out->kind = QAWS_SUPPORT_LOCAL;
	out->range_count = 1;
	out->ranges[0].field = QAWS_FIELD_CONTROL_POINTS;
	out->ranges[0].first = knot_span - degree;
	out->ranges[0].count = degree + 1;
	out->weight_field = QAWS_FIELD_NONE;
	out->has_weights = 1;
	out->order = order;
	for (r = 0; r <= order; r++)
		for (j = 0; j <= degree; j++)
			out->weights[0][r][j] = ders[r * (degree + 1) + j];
	return QAWS_STATUS_OK;
}

static qaws_status bspline_rebuild(qaws_curve const *curve, qaws_diff_views const *values, qaws_curve **out_curve)
{
	qaws_bspline_impl const *impl = (qaws_bspline_impl const *)curve->impl;
	qaws_bspline_desc d;
	qaws_scalar const *cps = NULL, *knots = NULL;
	qaws_scalar *owned_cps = NULL, *owned_knots = NULL;
	qaws_status st = qaws_internal_field_override(&cps, &owned_cps, values,
		impl->control_points, impl->control_point_count, QAWS_FIELD_CONTROL_POINTS, (unsigned int)curve->dimension);
	if (st != QAWS_STATUS_OK)
		return st;
	st = qaws_internal_field_override(&knots, &owned_knots, values, impl->knots, impl->knot_count, QAWS_FIELD_KNOTS, 1);
	if (st == QAWS_STATUS_OK) {
		memset(&d, 0, sizeof(d));
		d.dimension = curve->dimension;
		d.degree = curve->degree;
		d.control_points = cps;
		d.control_point_count = impl->control_point_count;
		d.knots = knots;
		d.knot_count = impl->knot_count;
		st = qaws_curve_create_bspline(out_curve, &d);
	}
	qaws_internal_dealloc(owned_cps, NULL);
	qaws_internal_dealloc(owned_knots, NULL);
	return st;
}

static qaws_curve_diff_vtable const bspline_diff_vtable = {
	BSPLINE_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	bspline_describe_fields,
	bspline_primal_field,
	bspline_linear_support,
	NULL,
	NULL,
	NULL,
	bspline_rebuild
};

static qaws_curve_vtable const bspline_vtable = {
	bspline_eval_span_2d,
	bspline_eval_span_3d,
	bspline_destroy_impl,
	bspline_is_closed,
	bspline_is_periodic,
	bspline_is_rational,
	bspline_get_continuity,
	&bspline_diff_vtable
};

/* -------------------------------------------------------------------------- */
/*  Creation                                                                  */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_bspline_ex(
	qaws_curve **out_curve,
	qaws_bspline_desc const *desc,
	qaws_allocator const *allocator)
{
	qaws_bspline_impl *impl = NULL;
	qaws_curve *curve = NULL;
	unsigned int dim_count;
	unsigned int degree;
	unsigned int N;
	unsigned int knot_count;
	qaws_scalar *knots = NULL;
	int owns_knots = 0;
	unsigned int span_count;
	unsigned int i;
	qaws_status status;
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

	degree = desc->degree;
	if (degree < 1 || degree > QAWS_INTERNAL_BASIS_MAX_DEGREE)
		return QAWS_STATUS_INVALID_DEGREE;

	N = desc->control_point_count;
	if (N < degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	if (N > 0xFFFFFFFFu / 4)   /* N * (dim_count + 1) homogeneous scalars, N + knots, count in unsigned int */
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	/* --- Build or validate knot vector --------------------------------- */

	if (desc->is_uniform && !desc->knots) {
		/* Generate clamped uniform knot vector */
		status = generate_uniform_knots(degree, N, allocator, &knots, &knot_count);
		if (status != QAWS_STATUS_OK)
			return status;
		owns_knots = 1;
	} else if (desc->knots) {
		knot_count = desc->knot_count;

		/* Validate the knot vector */
		status = qaws_internal_validate_knot_vector(
			desc->knots, knot_count, degree, N);
		if (status != QAWS_STATUS_OK)
			return status;

		/* Copy the knot vector */
		knots = (qaws_scalar *)qaws_internal_alloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)knot_count),
			allocator);
		if (!knots)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		memcpy(knots, desc->knots,
		       sizeof(qaws_scalar) * (size_t)knot_count);
		owns_knots = 1;
	} else {
		/* No knots and not uniform: generate uniform by default */
		status = generate_uniform_knots(degree, N, allocator, &knots, &knot_count);
		if (status != QAWS_STATUS_OK)
			return status;
		owns_knots = 1;
	}

	/* --- Compute span boundaries --------------------------------------- */
	/*
	 * Find distinct knot intervals in [knots[degree], knots[N]].
	 * Each non-zero-length interval is one span.
	 */
	{
		qaws_scalar param_start = knots[degree];
		unsigned int max_spans = N - degree; /* upper bound on span count */
		qaws_scalar *temp_bounds;
		unsigned int count;

		temp_bounds = (qaws_scalar *)qaws_internal_alloc(qaws_internal_mul2(sizeof(qaws_scalar), ((size_t)(max_spans + 1))),
			allocator);
		if (!temp_bounds) {
			if (owns_knots) qaws_internal_dealloc(knots, allocator);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		count = 0;
		temp_bounds[0] = param_start;

		for (i = degree + 1; i <= N; i++) {
			if (knots[i] > temp_bounds[count]) {
				count++;
				temp_bounds[count] = knots[i];
			}
		}

		span_count = count;
		if (span_count < 1) {
			qaws_internal_dealloc(temp_bounds, allocator);
			if (owns_knots) qaws_internal_dealloc(knots, allocator);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}

		/* --- Allocate curve -------------------------------------------- */

		range.min_value = knots[degree];
		range.max_value = knots[N];

		curve = qaws_internal_curve_alloc_ex(
			QAWS_CURVE_KIND_BSPLINE,
			desc->dimension,
			degree,
			span_count,
			range,
			&bspline_vtable,
			allocator);
		if (!curve) {
			qaws_internal_dealloc(temp_bounds, allocator);
			if (owns_knots) qaws_internal_dealloc(knots, allocator);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		/* Copy span boundaries into the curve */
		memcpy(curve->span_boundaries, temp_bounds,
		       sizeof(qaws_scalar) * (size_t)(span_count + 1));
		qaws_internal_dealloc(temp_bounds, allocator);
	}

	/* --- Allocate and populate implementation -------------------------- */

	impl = (qaws_bspline_impl *)qaws_internal_alloc((unsigned long)sizeof(qaws_bspline_impl),
		allocator);
	if (!impl) {
		qaws_internal_curve_free(curve);
		if (owns_knots) qaws_internal_dealloc(knots, allocator);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(impl, 0, sizeof(qaws_bspline_impl));

	impl->control_point_count = N;
	impl->knots = knots;
	impl->knot_count = knot_count;
	owns_knots = 0; /* impl now owns the knots */

	/* Copy control points */
	impl->control_points = (qaws_scalar *)qaws_internal_alloc(qaws_internal_mul3(sizeof(qaws_scalar), N, dim_count),
		allocator);
	if (!impl->control_points) {
		qaws_internal_dealloc(impl->knots, allocator);
		qaws_internal_dealloc(impl, allocator);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points,
	       sizeof(qaws_scalar) * (size_t)(N * dim_count));

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_bspline(
	qaws_curve **out_curve,
	qaws_bspline_desc const *desc)
{
	return qaws_curve_create_bspline_ex(out_curve, desc, NULL);
}
