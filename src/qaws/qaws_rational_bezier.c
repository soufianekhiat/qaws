#include "qaws_rational_bezier.h"
#include "qaws_curve.h"
#include "qaws_prepare.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include "internal/qaws_internal_diff.h"
#include <stdlib.h>
#include <string.h>
#include "qaws_platform.h"
#include "core/qaws_rational_bezier_core.h"

/* ---------------------------------------------------------------------------
 * Helpers: evaluate homogeneous Bezier and its derivatives using De Casteljau
 *
 * The weighted control points are stored in homogeneous form:
 *   (w*x, w*y, [w*z], w)  -- (dim_count+1) components per point.
 *
 * We use qaws_internal_decasteljau and qaws_internal_bezier_derivative_points
 * which work on flat arrays with arbitrary dimension count, so we pass
 * hdim = dim_count + 1 to operate on the homogeneous coordinates.
 *
 * After evaluating the homogeneous curve and its derivatives, we apply the
 * quotient rule to recover the rational point and its derivatives:
 *   C(t) = P(t) / w(t)
 *   C'(t) = (P'(t) - C(t)*w'(t)) / w(t)
 *   C''(t) = (P''(t) - 2*C'(t)*w'(t) - C(t)*w''(t)) / w(t)
 *   C'''(t) = (P'''(t) - 3*C''(t)*w'(t) - 3*C'(t)*w''(t) - C(t)*w'''(t)) / w(t)
 * where P(t) is the spatial part of the homogeneous evaluation (the first
 * dim_count components) and w(t) is the last component.
 * ------------------------------------------------------------------------- */

/* Stack buffer size for De Casteljau scratch space.
   Must accommodate (degree+1) * (dim_count+1) scalars for the largest
   derivative level we compute. 256 scalars handles degree <= ~63 in 3D. */
#define RBEZ_STACK_BUF 256

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_2d
 * ------------------------------------------------------------------------- */

static qaws_status rbez_eval_span_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	unsigned int const dim_count = 2;
	unsigned int const hdim = dim_count + 1; /* 3: wx, wy, w */
	qaws_rational_bezier_impl const* impl =
		(qaws_rational_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* wp = impl->weighted_points;

	/* Homogeneous evaluation buffers */
	qaws_scalar hbuf[3]; /* hdim components */
	qaws_scalar d1p[RBEZ_STACK_BUF], d2p[RBEZ_STACK_BUF], d3p[RBEZ_STACK_BUF];
	qaws_scalar hd1[3], hd2[3], hd3[3];

	/* Rational results */
	qaws_scalar w0, inv_w0;
	qaws_scalar C[2], C1[2], C2[2];

	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	/* --- Position: evaluate homogeneous curve --- */
	qaws_internal_decasteljau(wp, degree, hdim, local_t, hbuf);

	w0 = hbuf[dim_count]; /* last component is weight */
	if (w0 < QAWS_LITERAL(1e-15) && w0 > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;
	inv_w0 = QAWS_ONE / w0;

	C[0] = hbuf[0] * inv_w0;
	C[1] = hbuf[1] * inv_w0;

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = C[0];
		out_result->position.y = C[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* --- First derivative --- */
	if ((eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		&& degree >= 1)
	{
		qaws_internal_bezier_derivative_points(wp, degree, hdim, d1p);
		qaws_internal_decasteljau(d1p, degree - 1, hdim, local_t, hd1);

		/* Quotient rule: C' = (P' - C * w') / w */
		C1[0] = (hd1[0] - C[0] * hd1[dim_count]) * inv_w0;
		C1[1] = (hd1[1] - C[1] * hd1[dim_count]) * inv_w0;

		if (eval_flags & QAWS_EVAL_FLAG_D1) {
			out_result->d1.x = C1[0];
			out_result->d1.y = C1[1];
			out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
		}
	}

	/* --- Second derivative --- */
	if ((eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		&& degree >= 2)
	{
		qaws_internal_bezier_derivative_points(d1p, degree - 1, hdim, d2p);
		qaws_internal_decasteljau(d2p, degree - 2, hdim, local_t, hd2);

		/* C'' = (P'' - 2*C'*w' - C*w'') / w */
		C2[0] = (hd2[0] - QAWS_LITERAL(2.0) * C1[0] * hd1[dim_count]
			- C[0] * hd2[dim_count]) * inv_w0;
		C2[1] = (hd2[1] - QAWS_LITERAL(2.0) * C1[1] * hd1[dim_count]
			- C[1] * hd2[dim_count]) * inv_w0;

		if (eval_flags & QAWS_EVAL_FLAG_D2) {
			out_result->d2.x = C2[0];
			out_result->d2.y = C2[1];
			out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
		}
	}

	/* --- Third derivative --- */
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3)
	{
		qaws_internal_bezier_derivative_points(d2p, degree - 2, hdim, d3p);
		qaws_internal_decasteljau(d3p, degree - 3, hdim, local_t, hd3);

		/* C''' = (P''' - 3*C''*w' - 3*C'*w'' - C*w''') / w */
		out_result->d3.x = (hd3[0]
			- QAWS_LITERAL(3.0) * C2[0] * hd1[dim_count]
			- QAWS_LITERAL(3.0) * C1[0] * hd2[dim_count]
			- C[0] * hd3[dim_count]) * inv_w0;
		out_result->d3.y = (hd3[1]
			- QAWS_LITERAL(3.0) * C2[1] * hd1[dim_count]
			- QAWS_LITERAL(3.0) * C1[1] * hd2[dim_count]
			- C[1] * hd3[dim_count]) * inv_w0;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_3d
 * ------------------------------------------------------------------------- */

static qaws_status rbez_eval_span_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	unsigned int const dim_count = 3;
	unsigned int const hdim = dim_count + 1; /* 4: wx, wy, wz, w */
	qaws_rational_bezier_impl const* impl =
		(qaws_rational_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* wp = impl->weighted_points;

	/* Homogeneous evaluation buffers */
	qaws_scalar hbuf[4]; /* hdim components */
	qaws_scalar d1p[RBEZ_STACK_BUF], d2p[RBEZ_STACK_BUF], d3p[RBEZ_STACK_BUF];
	qaws_scalar hd1[4], hd2[4], hd3[4];

	/* Rational results */
	qaws_scalar w0, inv_w0;
	qaws_scalar C[3], C1[3], C2[3];

	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	/* --- Position: evaluate homogeneous curve --- */
	qaws_internal_decasteljau(wp, degree, hdim, local_t, hbuf);

	w0 = hbuf[dim_count]; /* last component is weight */
	if (w0 < QAWS_LITERAL(1e-15) && w0 > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;
	inv_w0 = QAWS_ONE / w0;

	C[0] = hbuf[0] * inv_w0;
	C[1] = hbuf[1] * inv_w0;
	C[2] = hbuf[2] * inv_w0;

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = C[0];
		out_result->position.y = C[1];
		out_result->position.z = C[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* --- First derivative --- */
	if ((eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		&& degree >= 1)
	{
		qaws_internal_bezier_derivative_points(wp, degree, hdim, d1p);
		qaws_internal_decasteljau(d1p, degree - 1, hdim, local_t, hd1);

		/* Quotient rule: C' = (P' - C * w') / w */
		C1[0] = (hd1[0] - C[0] * hd1[dim_count]) * inv_w0;
		C1[1] = (hd1[1] - C[1] * hd1[dim_count]) * inv_w0;
		C1[2] = (hd1[2] - C[2] * hd1[dim_count]) * inv_w0;

		if (eval_flags & QAWS_EVAL_FLAG_D1) {
			out_result->d1.x = C1[0];
			out_result->d1.y = C1[1];
			out_result->d1.z = C1[2];
			out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
		}
	}

	/* --- Second derivative --- */
	if ((eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		&& degree >= 2)
	{
		qaws_internal_bezier_derivative_points(d1p, degree - 1, hdim, d2p);
		qaws_internal_decasteljau(d2p, degree - 2, hdim, local_t, hd2);

		/* C'' = (P'' - 2*C'*w' - C*w'') / w */
		C2[0] = (hd2[0] - QAWS_LITERAL(2.0) * C1[0] * hd1[dim_count]
			- C[0] * hd2[dim_count]) * inv_w0;
		C2[1] = (hd2[1] - QAWS_LITERAL(2.0) * C1[1] * hd1[dim_count]
			- C[1] * hd2[dim_count]) * inv_w0;
		C2[2] = (hd2[2] - QAWS_LITERAL(2.0) * C1[2] * hd1[dim_count]
			- C[2] * hd2[dim_count]) * inv_w0;

		if (eval_flags & QAWS_EVAL_FLAG_D2) {
			out_result->d2.x = C2[0];
			out_result->d2.y = C2[1];
			out_result->d2.z = C2[2];
			out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
		}
	}

	/* --- Third derivative --- */
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3)
	{
		qaws_internal_bezier_derivative_points(d2p, degree - 2, hdim, d3p);
		qaws_internal_decasteljau(d3p, degree - 3, hdim, local_t, hd3);

		/* C''' = (P''' - 3*C''*w' - 3*C'*w'' - C*w''') / w */
		out_result->d3.x = (hd3[0]
			- QAWS_LITERAL(3.0) * C2[0] * hd1[dim_count]
			- QAWS_LITERAL(3.0) * C1[0] * hd2[dim_count]
			- C[0] * hd3[dim_count]) * inv_w0;
		out_result->d3.y = (hd3[1]
			- QAWS_LITERAL(3.0) * C2[1] * hd1[dim_count]
			- QAWS_LITERAL(3.0) * C1[1] * hd2[dim_count]
			- C[1] * hd3[dim_count]) * inv_w0;
		out_result->d3.z = (hd3[2]
			- QAWS_LITERAL(3.0) * C2[2] * hd1[dim_count]
			- QAWS_LITERAL(3.0) * C1[2] * hd2[dim_count]
			- C[2] * hd3[dim_count]) * inv_w0;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: destroy
 * ------------------------------------------------------------------------- */

static void rbez_destroy_impl(void* impl, qaws_allocator const* allocator)
{
	qaws_rational_bezier_impl* ri = (qaws_rational_bezier_impl*)impl;
	if (ri) {
		qaws_internal_dealloc(allocator, ri->control_points);
		qaws_internal_dealloc(allocator, ri->weights);
		qaws_internal_dealloc(allocator, ri->weighted_points);
		qaws_internal_dealloc(allocator, ri);
	}
}

/* ---------------------------------------------------------------------------
 * Vtable: property queries
 * ------------------------------------------------------------------------- */

static int rbez_is_closed(qaws_curve const* c)   { (void)c; return 0; }
static int rbez_is_periodic(qaws_curve const* c)  { (void)c; return 0; }
static int rbez_is_rational(qaws_curve const* c)  { (void)c; return 1; }
static qaws_continuity rbez_get_continuity(qaws_curve const* c) { (void)c; return QAWS_CONTINUITY_C3; }

/* ---------------------------------------------------------------------------
 * Vtable definition
 * ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*  Differential rules: rational in control points and weights               */
/* -------------------------------------------------------------------------- */

#define RBEZ_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2 | \
	QAWS_CAP_LOCAL_SUPPORT)

static unsigned int rbez_describe_fields(qaws_curve const* curve, qaws_field_desc* out, unsigned int capacity)
{
	qaws_rational_bezier_impl const* impl = (qaws_rational_bezier_impl const*)curve->impl;
	if (capacity >= 1)
		out[0] = qaws_internal_field_desc(QAWS_FIELD_CONTROL_POINTS, (qaws_value_type)curve->dimension,
			impl->control_point_count, QAWS_DOMAIN_POSITION, QAWS_CONSTRAINT_NONE,
			QAWS_DIFF_SMOOTH, RBEZ_DIFF_CAPS);
	if (capacity >= 2)
		out[1] = qaws_internal_field_desc(QAWS_FIELD_WEIGHTS, QAWS_VALUE_SCALAR,
			impl->control_point_count, QAWS_DOMAIN_WEIGHT, QAWS_CONSTRAINT_POSITIVE,
			QAWS_DIFF_SMOOTH, RBEZ_DIFF_CAPS);
	return 2;
}

static qaws_status rbez_primal_field(qaws_curve const* curve, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	qaws_rational_bezier_impl const* impl = (qaws_rational_bezier_impl const*)curve->impl;
	if (field == QAWS_FIELD_CONTROL_POINTS) {
		*out_data = impl->control_points;
		*out_components = (unsigned int)curve->dimension;
	} else if (field == QAWS_FIELD_WEIGHTS) {
		*out_data = impl->weights;
		*out_components = 1;
	} else {
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	*out_count = impl->control_point_count;
	return QAWS_STATUS_OK;
}

/* Homogeneous Bernstein basis: the engine applies the quotient rule. */
static qaws_status rbez_linear_support(qaws_curve const* curve, unsigned int span_index,
	qaws_scalar local_t, unsigned int order, qaws_local_support* out)
{
	qaws_scalar b[(QAWS_DIFF_MAX_ORDER + 1) * QAWS_DIFF_MAX_SUPPORT];
	unsigned int n = curve->degree, r, j;
	(void)span_index;
	if (n + 1 > QAWS_DIFF_MAX_SUPPORT || order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	qaws_internal_bernstein_derivs(n, local_t, order, b);

	out->kind = QAWS_SUPPORT_GLOBAL;
	out->range_count = 1;
	out->ranges[0].field = QAWS_FIELD_CONTROL_POINTS;
	out->ranges[0].first = 0;
	out->ranges[0].count = n + 1;
	out->weight_field = QAWS_FIELD_WEIGHTS;
	out->has_weights = 1;
	out->order = order;
	for (r = 0; r <= order; r++)
		for (j = 0; j <= n; j++)
			out->weights[0][r][j] = b[r * (n + 1) + j];
	return QAWS_STATUS_OK;
}

static qaws_curve_diff_vtable const rbez_diff_vtable = {
	RBEZ_DIFF_CAPS,
	QAWS_DIFF_SMOOTH,
	rbez_describe_fields,
	rbez_primal_field,
	rbez_linear_support,
	NULL,
	NULL,
	NULL
};

static qaws_curve_vtable const rbez_vtable = {
	rbez_eval_span_2d,
	rbez_eval_span_3d,
	rbez_destroy_impl,
	rbez_is_closed,
	rbez_is_periodic,
	rbez_is_rational,
	rbez_get_continuity,
	&rbez_diff_vtable
};

/* ---------------------------------------------------------------------------
 * Creation
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_create_rational_bezier(
	qaws_rational_bezier_desc const* desc,
	qaws_curve** out_curve)
{
	unsigned int dim_count;
	unsigned int n_pts;
	unsigned int hdim;
	size_t cp_size, w_size, wp_size;
	qaws_range range;
	qaws_curve* curve;
	qaws_rational_bezier_impl* impl;
	qaws_status status;

	if (!desc || !out_curve) return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;

	/* Validate dimension */
	status = qaws_internal_validate_dimension(desc->dimension);
	if (status != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;

	dim_count = (unsigned int)desc->dimension;

	if (desc->degree < 1) return QAWS_STATUS_INVALID_DEGREE;

	n_pts = desc->degree + 1;
	if (desc->control_point_count != n_pts)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate weights */
	if (!desc->weights)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->weight_count != n_pts)
		return QAWS_STATUS_INVALID_WEIGHT_COUNT;

	status = qaws_internal_validate_weights(desc->weights, desc->weight_count);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Allocate curve */
	range.min_value = 0; range.max_value = 1;
	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_RATIONAL_BEZIER,
		desc->dimension, desc->degree, 1, range, &rbez_vtable);
	if (!curve) return QAWS_STATUS_ALLOCATION_FAILURE;
	curve->span_boundaries[0] = 0;
	curve->span_boundaries[1] = 1;

	/* Allocate impl */
	impl = (qaws_rational_bezier_impl*)calloc(1, sizeof(qaws_rational_bezier_impl));
	if (!impl) {
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy control points */
	cp_size = (size_t)dim_count * n_pts * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar*)malloc(cp_size);
	if (!impl->control_points) {
		free(impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points, cp_size);
	impl->control_point_count = n_pts;

	/* Copy weights */
	w_size = (size_t)n_pts * sizeof(qaws_scalar);
	impl->weights = (qaws_scalar*)malloc(w_size);
	if (!impl->weights) {
		free(impl->control_points);
		free(impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->weights, desc->weights, w_size);

	/* Pre-compute weighted (homogeneous) control points:
	   for each point i, store (w_i * x, w_i * y, [w_i * z], w_i) */
	hdim = dim_count + 1;
	wp_size = (size_t)hdim * n_pts * sizeof(qaws_scalar);
	impl->weighted_points = (qaws_scalar*)malloc(wp_size);
	if (!impl->weighted_points) {
		free(impl->weights);
		free(impl->control_points);
		free(impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	if (dim_count == 2) {
		qaws_rational_bezier_prepare_2d(
			(const qaws_vec2 *)desc->control_points,
			desc->weights, desc->degree,
			impl->weighted_points);
	} else {
		qaws_rational_bezier_prepare_3d(
			(const qaws_vec3 *)desc->control_points,
			desc->weights, desc->degree,
			impl->weighted_points);
	}

	curve->impl = impl;
	*out_curve = curve;
	return QAWS_STATUS_OK;
}
