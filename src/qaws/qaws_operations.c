#include "qaws_operations.h"
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_trajectory.h"
#include "qaws_yuksel.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_curve.h"
#include "qaws_composite.h"
#include "qaws_arc.h"
#include "qaws_surface.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_fit.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_span.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* ========================================================================== */
/*  Forward declarations of per-family split/join helpers                      */
/* ========================================================================== */

static qaws_status split_bezier(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right);

static qaws_status split_hermite(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right);

static qaws_status split_bspline(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right);

static qaws_status split_nurbs(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right);

static qaws_status split_catmull_rom(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right);

static qaws_status split_trajectory(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right);

static qaws_status join_hermite(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined);

static qaws_status join_bspline(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined);

static qaws_status join_nurbs(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined);

static qaws_status join_trajectory(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined);

static qaws_status join_catmull_rom(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined);

/* Split through qaws_curve_extract: two exact pieces, possibly of a simpler kind. */
static qaws_status split_by_extract(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right)
{
	qaws_status status = qaws_curve_extract(curve,
		curve->parameter_range.min_value, t, out_left);
	if (status != QAWS_STATUS_OK)
		return status;
	status = qaws_curve_extract(curve, t,
		curve->parameter_range.max_value, out_right);
	if (status != QAWS_STATUS_OK) {
		qaws_curve_destroy(*out_left);
		*out_left = NULL;
	}
	return status;
}

/* ========================================================================== */
/*  Public API: qaws_curve_split                                              */
/* ========================================================================== */

qaws_status qaws_curve_split(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_curve **out_left,
	qaws_curve **out_right)
{
	if (!curve || !out_left || !out_right)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_left = NULL;
	*out_right = NULL;

	/* Validate parameter is within range */
	if (parameter < curve->parameter_range.min_value ||
	    parameter > curve->parameter_range.max_value)
		return QAWS_STATUS_OUT_OF_RANGE;

	/* If parameter is at the boundary, splitting is degenerate */
	if (parameter <= curve->parameter_range.min_value ||
	    parameter >= curve->parameter_range.max_value)
		return QAWS_STATUS_INVALID_ARGUMENT;

	switch (curve->kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
		return split_bezier(curve, parameter, out_left, out_right);

	case QAWS_CURVE_KIND_HERMITE:
		/* a Hermite span cannot be cut inside and stay a uniform Hermite */
		if (parameter != (qaws_scalar)floor((double)parameter))
			return split_by_extract(curve, parameter, out_left, out_right);
		return split_hermite(curve, parameter, out_left, out_right);

	case QAWS_CURVE_KIND_BSPLINE:
		return split_bspline(curve, parameter, out_left, out_right);

	case QAWS_CURVE_KIND_NURBS:
		return split_nurbs(curve, parameter, out_left, out_right);

	case QAWS_CURVE_KIND_CATMULL_ROM:
	{
		qaws_status status = split_catmull_rom(curve, parameter, out_left, out_right);
		if (status == QAWS_STATUS_UNSUPPORTED_OPERATION)
			return split_by_extract(curve, parameter, out_left, out_right);
		return status;
	}

	case QAWS_CURVE_KIND_TRAJECTORY:
		return split_trajectory(curve, parameter, out_left, out_right);

	case QAWS_CURVE_KIND_YUKSEL:
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	case QAWS_CURVE_KIND_COMPOSITE:
	case QAWS_CURVE_KIND_ARC:
	case QAWS_CURVE_KIND_POLYNOMIAL:
	case QAWS_CURVE_KIND_CLOTHOID:
	case QAWS_CURVE_KIND_SUBDIVISION:
	case QAWS_CURVE_KIND_REPARAMETERIZED:
		return split_by_extract(curve, parameter, out_left, out_right);

	default:
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
}

/* ========================================================================== */
/*  Public API: qaws_curve_join                                               */
/* ========================================================================== */

qaws_status qaws_curve_join(
	qaws_curve const *curve_a,
	qaws_curve const *curve_b,
	qaws_curve **out_joined)
{
	if (!curve_a || !curve_b || !out_joined)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_joined = NULL;

	/* Both curves must have the same kind and dimension */
	if (curve_a->kind != curve_b->kind)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve_a->dimension != curve_b->dimension)
		return QAWS_STATUS_INVALID_ARGUMENT;

	switch (curve_a->kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	case QAWS_CURVE_KIND_HERMITE:
		return join_hermite(curve_a, curve_b, out_joined);

	case QAWS_CURVE_KIND_BSPLINE:
		return join_bspline(curve_a, curve_b, out_joined);

	case QAWS_CURVE_KIND_NURBS:
		return join_nurbs(curve_a, curve_b, out_joined);

	case QAWS_CURVE_KIND_CATMULL_ROM:
		return join_catmull_rom(curve_a, curve_b, out_joined);

	case QAWS_CURVE_KIND_TRAJECTORY:
		return join_trajectory(curve_a, curve_b, out_joined);

	case QAWS_CURVE_KIND_YUKSEL:
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	case QAWS_CURVE_KIND_COMPOSITE:
	case QAWS_CURVE_KIND_ARC:
	case QAWS_CURVE_KIND_POLYNOMIAL:
	case QAWS_CURVE_KIND_CLOTHOID:
	case QAWS_CURVE_KIND_SUBDIVISION:
	case QAWS_CURVE_KIND_REPARAMETERIZED:
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	default:
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
}

/* ========================================================================== */
/*  Bezier split -- De Casteljau                                              */
/* ========================================================================== */

static qaws_status split_bezier(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right)
{
	qaws_bezier_impl const *impl = (qaws_bezier_impl const *)curve->impl;
	unsigned int degree = curve->degree;
	unsigned int order = degree + 1;
	unsigned int dim_count = (unsigned int)curve->dimension;
	qaws_scalar *tri = NULL;
	qaws_scalar *left_cp = NULL;
	qaws_scalar *right_cp = NULL;
	qaws_scalar one_minus_t;
	unsigned int r, i, d;
	size_t tri_size;
	size_t cp_size;
	qaws_bezier_desc left_desc;
	qaws_bezier_desc right_desc;
	qaws_status status;

	/*
	 * Build the De Casteljau triangle as a flat 2D array:
	 *   tri[r * order * dim_count + i * dim_count + d]
	 * where r is the level (0..degree), i is the index within that level.
	 */
	tri_size = (size_t)order * (size_t)order * (size_t)dim_count * sizeof(qaws_scalar);
	tri = (qaws_scalar *)malloc(tri_size);
	if (!tri)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Level 0: copy the original control points */
	memcpy(tri, impl->control_points,
	       (size_t)order * (size_t)dim_count * sizeof(qaws_scalar));

	/* Build levels 1..degree */
	one_minus_t = (qaws_scalar)1.0 - t;
	for (r = 1; r <= degree; r++)
	{
		unsigned int prev_row = (r - 1) * order * dim_count;
		unsigned int curr_row = r * order * dim_count;
		for (i = 0; i <= degree - r; i++)
		{
			for (d = 0; d < dim_count; d++)
			{
				tri[curr_row + i * dim_count + d] =
					one_minus_t * tri[prev_row + i * dim_count + d] +
					t * tri[prev_row + (i + 1) * dim_count + d];
			}
		}
	}

	/* Extract left control points: tri[j][0] for j = 0..degree */
	cp_size = (size_t)order * (size_t)dim_count * sizeof(qaws_scalar);
	left_cp = (qaws_scalar *)malloc(cp_size);
	if (!left_cp)
	{
		free(tri);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (r = 0; r <= degree; r++)
	{
		for (d = 0; d < dim_count; d++)
			left_cp[r * dim_count + d] = tri[r * order * dim_count + d];
	}

	/* Extract right control points: tri[degree-j][j] for j = 0..degree */
	right_cp = (qaws_scalar *)malloc(cp_size);
	if (!right_cp)
	{
		free(left_cp);
		free(tri);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (r = 0; r <= degree; r++)
	{
		unsigned int level = degree - r;
		unsigned int idx = r;
		for (d = 0; d < dim_count; d++)
			right_cp[r * dim_count + d] = tri[level * order * dim_count + idx * dim_count + d];
	}

	free(tri);

	/* Create the left sub-curve */
	memset(&left_desc, 0, sizeof(left_desc));
	left_desc.dimension = curve->dimension;
	left_desc.degree = degree;
	left_desc.control_points = left_cp;
	left_desc.control_point_count = order;

	status = qaws_curve_create_bezier(&left_desc, out_left);
	if (status != QAWS_STATUS_OK)
	{
		free(left_cp);
		free(right_cp);
		return status;
	}

	/* Create the right sub-curve */
	memset(&right_desc, 0, sizeof(right_desc));
	right_desc.dimension = curve->dimension;
	right_desc.degree = degree;
	right_desc.control_points = right_cp;
	right_desc.control_point_count = order;

	status = qaws_curve_create_bezier(&right_desc, out_right);
	if (status != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(*out_left);
		*out_left = NULL;
		free(left_cp);
		free(right_cp);
		return status;
	}

	free(left_cp);
	free(right_cp);
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Hermite split                                                             */
/* ========================================================================== */

static qaws_status split_hermite(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
	unsigned int dim_count = (unsigned int)curve->dimension;
	unsigned int total_points = impl->point_count;
	unsigned int span_index;
	qaws_scalar local_t;
	qaws_scalar split_pos[3];
	qaws_scalar split_tan[3];
	unsigned int left_count;
	unsigned int right_count;
	qaws_scalar *left_pts = NULL;
	qaws_scalar *left_tan = NULL;
	qaws_scalar *right_pts = NULL;
	qaws_scalar *right_tan = NULL;
	size_t left_size;
	size_t right_size;
	qaws_hermite_desc left_desc;
	qaws_hermite_desc right_desc;
	qaws_status status;
	unsigned int i, d;

	/* Determine which span the parameter falls into */
	span_index = (unsigned int)floor((double)t);
	if (span_index >= total_points - 1)
		span_index = total_points - 2;

	local_t = t - (qaws_scalar)span_index;

	/* Check if splitting at an exact boundary */
	if (local_t < (qaws_scalar)1e-10)
	{
		/* Splitting at an integer boundary -- partition the arrays */
		if (span_index < 1)
			return QAWS_STATUS_INVALID_ARGUMENT;

		left_count = span_index + 1;
		right_count = total_points - span_index;
	}
	else
	{
		/* Splitting within a span -- evaluate position and D1 at the point */
		left_count = span_index + 2;
		right_count = total_points - span_index - 1 + 1;
	}

	/* Evaluate position and tangent at the split parameter */
	memset(split_pos, 0, sizeof(split_pos));
	memset(split_tan, 0, sizeof(split_tan));

	if (dim_count == 2)
	{
		qaws_eval_result_2d result;
		status = qaws_curve_evaluate_2d(curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);
		if (status != QAWS_STATUS_OK)
			return status;
		split_pos[0] = result.position.x;
		split_pos[1] = result.position.y;
		split_tan[0] = result.d1.x;
		split_tan[1] = result.d1.y;
	}
	else
	{
		qaws_eval_result_3d result;
		status = qaws_curve_evaluate_3d(curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);
		if (status != QAWS_STATUS_OK)
			return status;
		split_pos[0] = result.position.x;
		split_pos[1] = result.position.y;
		split_pos[2] = result.position.z;
		split_tan[0] = result.d1.x;
		split_tan[1] = result.d1.y;
		split_tan[2] = result.d1.z;
	}

	/* Allocate left curve data */
	left_size = (size_t)left_count * (size_t)dim_count * sizeof(qaws_scalar);
	left_pts = (qaws_scalar *)malloc(left_size);
	left_tan = (qaws_scalar *)malloc(left_size);
	if (!left_pts || !left_tan)
	{
		free(left_pts);
		free(left_tan);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Allocate right curve data */
	right_size = (size_t)right_count * (size_t)dim_count * sizeof(qaws_scalar);
	right_pts = (qaws_scalar *)malloc(right_size);
	right_tan = (qaws_scalar *)malloc(right_size);
	if (!right_pts || !right_tan)
	{
		free(left_pts);
		free(left_tan);
		free(right_pts);
		free(right_tan);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	if (local_t < (qaws_scalar)1e-10)
	{
		/* Boundary split: partition at span_index */
		memcpy(left_pts, impl->points,
		       (size_t)left_count * (size_t)dim_count * sizeof(qaws_scalar));
		memcpy(left_tan, impl->tangents,
		       (size_t)left_count * (size_t)dim_count * sizeof(qaws_scalar));

		memcpy(right_pts, &impl->points[span_index * dim_count],
		       (size_t)right_count * (size_t)dim_count * sizeof(qaws_scalar));
		memcpy(right_tan, &impl->tangents[span_index * dim_count],
		       (size_t)right_count * (size_t)dim_count * sizeof(qaws_scalar));
	}
	else
	{
		/* Within-span split */
		/* Left: points[0..span_index] + split_point */
		for (i = 0; i <= span_index; i++)
		{
			for (d = 0; d < dim_count; d++)
			{
				left_pts[i * dim_count + d] = impl->points[i * dim_count + d];
				left_tan[i * dim_count + d] = impl->tangents[i * dim_count + d];
			}
		}
		for (d = 0; d < dim_count; d++)
		{
			left_pts[(span_index + 1) * dim_count + d] = split_pos[d];
			left_tan[(span_index + 1) * dim_count + d] = split_tan[d];
		}

		/* Right: split_point + points[span_index+1..total-1] */
		for (d = 0; d < dim_count; d++)
		{
			right_pts[d] = split_pos[d];
			right_tan[d] = split_tan[d];
		}
		for (i = span_index + 1; i < total_points; i++)
		{
			unsigned int ri = i - span_index;
			for (d = 0; d < dim_count; d++)
			{
				right_pts[ri * dim_count + d] = impl->points[i * dim_count + d];
				right_tan[ri * dim_count + d] = impl->tangents[i * dim_count + d];
			}
		}
	}

	/* Create left Hermite curve */
	memset(&left_desc, 0, sizeof(left_desc));
	left_desc.dimension = curve->dimension;
	left_desc.degree = 3;
	left_desc.points = left_pts;
	left_desc.derivatives = left_tan;
	left_desc.point_count = left_count;
	left_desc.derivative_count = left_count;

	status = qaws_curve_create_hermite(&left_desc, out_left);
	if (status != QAWS_STATUS_OK)
	{
		free(left_pts);
		free(left_tan);
		free(right_pts);
		free(right_tan);
		return status;
	}

	/* Create right Hermite curve */
	memset(&right_desc, 0, sizeof(right_desc));
	right_desc.dimension = curve->dimension;
	right_desc.degree = 3;
	right_desc.points = right_pts;
	right_desc.derivatives = right_tan;
	right_desc.point_count = right_count;
	right_desc.derivative_count = right_count;

	status = qaws_curve_create_hermite(&right_desc, out_right);
	if (status != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(*out_left);
		*out_left = NULL;
		free(left_pts);
		free(left_tan);
		free(right_pts);
		free(right_tan);
		return status;
	}

	free(left_pts);
	free(left_tan);
	free(right_pts);
	free(right_tan);
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  B-Spline split -- Boehm knot insertion                                    */
/* ========================================================================== */

/*
 * Perform a single Boehm knot insertion of value t into the B-Spline.
 * Produces new_knots (knot_count+1) and new_cps (cp_count+1).
 * Caller must free the returned arrays.
 */
static qaws_status boehm_insert_knot(
	qaws_scalar const *knots, unsigned int knot_count,
	qaws_scalar const *cps, unsigned int cp_count,
	unsigned int degree, unsigned int dim_count,
	qaws_scalar t,
	qaws_scalar **out_new_knots, unsigned int *out_new_knot_count,
	qaws_scalar **out_new_cps, unsigned int *out_new_cp_count)
{
	unsigned int s;
	unsigned int new_knot_count;
	unsigned int new_cp_count;
	qaws_scalar *new_knots = NULL;
	qaws_scalar *new_cps = NULL;
	unsigned int i, d;
	size_t knot_size;
	size_t cp_size;

	/* Find knot span s: knots[s] <= t < knots[s+1] */
	s = qaws_internal_find_knot_span(knots, knot_count, degree, cp_count, t);

	new_knot_count = knot_count + 1;
	new_cp_count = cp_count + 1;

	knot_size = (size_t)new_knot_count * sizeof(qaws_scalar);
	new_knots = (qaws_scalar *)malloc(knot_size);
	if (!new_knots)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	cp_size = (size_t)new_cp_count * (size_t)dim_count * sizeof(qaws_scalar);
	new_cps = (qaws_scalar *)malloc(cp_size);
	if (!new_cps)
	{
		free(new_knots);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Insert knot into knot vector */
	for (i = 0; i <= s; i++)
		new_knots[i] = knots[i];
	new_knots[s + 1] = t;
	for (i = s + 1; i < knot_count; i++)
		new_knots[i + 1] = knots[i];

	/*
	 * New control points:
	 *   Q[i] = P[i]                                  for i <= s - degree
	 *   Q[i] = alpha_i * P[i] + (1-alpha_i) * P[i-1] for s-degree+1 <= i <= s
	 *   Q[i] = P[i-1]                                for i >= s+1
	 *
	 * where alpha_i = (t - knots[i]) / (knots[i+degree] - knots[i])
	 */
	for (i = 0; i < new_cp_count; i++)
	{
		if (i <= s - degree)
		{
			for (d = 0; d < dim_count; d++)
				new_cps[i * dim_count + d] = cps[i * dim_count + d];
		}
		else if (i >= s + 1)
		{
			for (d = 0; d < dim_count; d++)
				new_cps[i * dim_count + d] = cps[(i - 1) * dim_count + d];
		}
		else
		{
			/* s - degree + 1 <= i <= s */
			qaws_scalar denom = knots[i + degree] - knots[i];
			qaws_scalar alpha;

			if (denom < (qaws_scalar)1e-15 && denom > (qaws_scalar)-1e-15)
				alpha = (qaws_scalar)0.0;
			else
				alpha = (t - knots[i]) / denom;

			for (d = 0; d < dim_count; d++)
			{
				new_cps[i * dim_count + d] =
					alpha * cps[i * dim_count + d] +
					((qaws_scalar)1.0 - alpha) * cps[(i - 1) * dim_count + d];
			}
		}
	}

	*out_new_knots = new_knots;
	*out_new_knot_count = new_knot_count;
	*out_new_cps = new_cps;
	*out_new_cp_count = new_cp_count;
	return QAWS_STATUS_OK;
}

/*
 * Count the multiplicity of knot value t in the knot vector.
 */
static unsigned int count_knot_multiplicity(
	qaws_scalar const *knots, unsigned int knot_count, qaws_scalar t)
{
	unsigned int count = 0;
	unsigned int i;
	qaws_scalar tol = (qaws_scalar)1e-10;

	for (i = 0; i < knot_count; i++)
	{
		qaws_scalar diff = knots[i] - t;
		if (diff < (qaws_scalar)0.0)
			diff = -diff;
		if (diff <= tol)
			count++;
	}
	return count;
}

/*
 * Find the index in the knot vector where knot value t appears with
 * the required multiplicity (the last consecutive occurrence).
 * Returns the index of the first occurrence of t.
 */
static unsigned int find_knot_split_index(
	qaws_scalar const *knots, unsigned int knot_count, qaws_scalar t)
{
	unsigned int i;
	qaws_scalar tol = (qaws_scalar)1e-10;

	for (i = 0; i < knot_count; i++)
	{
		qaws_scalar diff = knots[i] - t;
		if (diff < (qaws_scalar)0.0)
			diff = -diff;
		if (diff <= tol)
			return i;
	}
	return 0;
}

static qaws_status split_bspline(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right)
{
	qaws_bspline_impl const *impl = (qaws_bspline_impl const *)curve->impl;
	unsigned int degree = curve->degree;
	unsigned int dim_count = (unsigned int)curve->dimension;
	unsigned int target_mult = degree + 1;
	unsigned int current_mult;
	unsigned int insertions_needed;
	qaws_scalar *work_knots = NULL;
	qaws_scalar *work_cps = NULL;
	unsigned int work_knot_count;
	unsigned int work_cp_count;
	unsigned int split_knot_idx;
	unsigned int left_cp_count;
	unsigned int right_cp_count;
	unsigned int left_knot_count;
	unsigned int right_knot_count;
	qaws_bspline_desc left_desc;
	qaws_bspline_desc right_desc;
	qaws_scalar *left_knots = NULL;
	qaws_scalar *left_cps = NULL;
	qaws_scalar *right_knots = NULL;
	qaws_scalar *right_cps = NULL;
	qaws_status status;
	unsigned int ins;
	unsigned int i, d;

	/* Start with a copy of the current knots and control points */
	work_knot_count = impl->knot_count;
	work_cp_count = impl->control_point_count;

	work_knots = (qaws_scalar *)malloc((size_t)work_knot_count * sizeof(qaws_scalar));
	if (!work_knots)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memcpy(work_knots, impl->knots, (size_t)work_knot_count * sizeof(qaws_scalar));

	work_cps = (qaws_scalar *)malloc((size_t)work_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	if (!work_cps)
	{
		free(work_knots);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(work_cps, impl->control_points,
	       (size_t)work_cp_count * (size_t)dim_count * sizeof(qaws_scalar));

	/* Insert knot t until multiplicity reaches degree+1 */
	current_mult = count_knot_multiplicity(work_knots, work_knot_count, t);
	insertions_needed = (target_mult > current_mult) ? (target_mult - current_mult) : 0;

	for (ins = 0; ins < insertions_needed; ins++)
	{
		qaws_scalar *new_knots = NULL;
		qaws_scalar *new_cps = NULL;
		unsigned int new_knot_count;
		unsigned int new_cp_count;

		status = boehm_insert_knot(
			work_knots, work_knot_count,
			work_cps, work_cp_count,
			degree, dim_count, t,
			&new_knots, &new_knot_count,
			&new_cps, &new_cp_count);

		if (status != QAWS_STATUS_OK)
		{
			free(work_knots);
			free(work_cps);
			return status;
		}

		free(work_knots);
		free(work_cps);
		work_knots = new_knots;
		work_cps = new_cps;
		work_knot_count = new_knot_count;
		work_cp_count = new_cp_count;
	}

	/*
	 * After insertion, t has multiplicity degree+1 in the knot vector.
	 * Find where the split occurs.
	 */
	split_knot_idx = find_knot_split_index(work_knots, work_knot_count, t);

	/*
	 * After knot insertion to multiplicity degree+1, t appears at indices
	 * [split_knot_idx .. split_knot_idx + degree].
	 *
	 * Left curve:
	 *   Knots: work_knots[0 .. split_knot_idx + degree]
	 *   CPs:   work_cps[0 .. split_knot_idx - 1]
	 *
	 * Right curve:
	 *   Knots: work_knots[split_knot_idx .. work_knot_count - 1]
	 *   CPs:   work_cps[split_knot_idx .. work_cp_count - 1]
	 */
	left_cp_count = split_knot_idx;
	left_knot_count = split_knot_idx + degree + 1;

	right_cp_count = work_cp_count - split_knot_idx;
	right_knot_count = work_knot_count - split_knot_idx;

	/* Validate counts */
	if (left_cp_count < degree + 1 || right_cp_count < degree + 1)
	{
		free(work_knots);
		free(work_cps);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}

	/* Allocate left arrays */
	left_knots = (qaws_scalar *)malloc((size_t)left_knot_count * sizeof(qaws_scalar));
	left_cps = (qaws_scalar *)malloc((size_t)left_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	if (!left_knots || !left_cps)
	{
		free(left_knots);
		free(left_cps);
		free(work_knots);
		free(work_cps);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Allocate right arrays */
	right_knots = (qaws_scalar *)malloc((size_t)right_knot_count * sizeof(qaws_scalar));
	right_cps = (qaws_scalar *)malloc((size_t)right_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	if (!right_knots || !right_cps)
	{
		free(left_knots);
		free(left_cps);
		free(right_knots);
		free(right_cps);
		free(work_knots);
		free(work_cps);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy left control points */
	for (i = 0; i < left_cp_count; i++)
	{
		for (d = 0; d < dim_count; d++)
			left_cps[i * dim_count + d] = work_cps[i * dim_count + d];
	}

	/* Copy left knots */
	for (i = 0; i < left_knot_count; i++)
		left_knots[i] = work_knots[i];

	/* Copy right control points */
	for (i = 0; i < right_cp_count; i++)
	{
		for (d = 0; d < dim_count; d++)
			right_cps[i * dim_count + d] = work_cps[(split_knot_idx + i) * dim_count + d];
	}

	/* Copy right knots */
	for (i = 0; i < right_knot_count; i++)
		right_knots[i] = work_knots[split_knot_idx + i];

	free(work_knots);
	free(work_cps);

	/* Create left B-Spline */
	memset(&left_desc, 0, sizeof(left_desc));
	left_desc.dimension = curve->dimension;
	left_desc.degree = degree;
	left_desc.control_points = left_cps;
	left_desc.control_point_count = left_cp_count;
	left_desc.knots = left_knots;
	left_desc.knot_count = left_knot_count;
	left_desc.is_uniform = 0;
	left_desc.is_closed = 0;

	status = qaws_curve_create_bspline(&left_desc, out_left);
	if (status != QAWS_STATUS_OK)
	{
		free(left_knots);
		free(left_cps);
		free(right_knots);
		free(right_cps);
		return status;
	}

	/* Create right B-Spline */
	memset(&right_desc, 0, sizeof(right_desc));
	right_desc.dimension = curve->dimension;
	right_desc.degree = degree;
	right_desc.control_points = right_cps;
	right_desc.control_point_count = right_cp_count;
	right_desc.knots = right_knots;
	right_desc.knot_count = right_knot_count;
	right_desc.is_uniform = 0;
	right_desc.is_closed = 0;

	status = qaws_curve_create_bspline(&right_desc, out_right);
	if (status != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(*out_left);
		*out_left = NULL;
		free(left_knots);
		free(left_cps);
		free(right_knots);
		free(right_cps);
		return status;
	}

	free(left_knots);
	free(left_cps);
	free(right_knots);
	free(right_cps);
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  NURBS split -- same as B-Spline but also split weights                    */
/* ========================================================================== */

/*
 * Perform a single Boehm knot insertion for NURBS (B-Spline CPs + weights).
 * Weights are treated as an extra 1D "control point" channel.
 */
static qaws_status boehm_insert_knot_nurbs(
	qaws_scalar const *knots, unsigned int knot_count,
	qaws_scalar const *cps, qaws_scalar const *weights, unsigned int cp_count,
	unsigned int degree, unsigned int dim_count,
	qaws_scalar t,
	qaws_scalar **out_new_knots, unsigned int *out_new_knot_count,
	qaws_scalar **out_new_cps, qaws_scalar **out_new_weights,
	unsigned int *out_new_cp_count)
{
	unsigned int s;
	unsigned int new_knot_count;
	unsigned int new_cp_count;
	qaws_scalar *new_knots = NULL;
	qaws_scalar *new_cps = NULL;
	qaws_scalar *new_weights = NULL;
	unsigned int i, d;

	s = qaws_internal_find_knot_span(knots, knot_count, degree, cp_count, t);

	new_knot_count = knot_count + 1;
	new_cp_count = cp_count + 1;

	new_knots = (qaws_scalar *)malloc((size_t)new_knot_count * sizeof(qaws_scalar));
	new_cps = (qaws_scalar *)malloc((size_t)new_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	new_weights = (qaws_scalar *)malloc((size_t)new_cp_count * sizeof(qaws_scalar));
	if (!new_knots || !new_cps || !new_weights)
	{
		free(new_knots);
		free(new_cps);
		free(new_weights);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Insert knot into knot vector */
	for (i = 0; i <= s; i++)
		new_knots[i] = knots[i];
	new_knots[s + 1] = t;
	for (i = s + 1; i < knot_count; i++)
		new_knots[i + 1] = knots[i];

	/* Compute new control points and weights */
	for (i = 0; i < new_cp_count; i++)
	{
		if (i <= s - degree)
		{
			for (d = 0; d < dim_count; d++)
				new_cps[i * dim_count + d] = cps[i * dim_count + d];
			new_weights[i] = weights[i];
		}
		else if (i >= s + 1)
		{
			for (d = 0; d < dim_count; d++)
				new_cps[i * dim_count + d] = cps[(i - 1) * dim_count + d];
			new_weights[i] = weights[i - 1];
		}
		else
		{
			qaws_scalar denom = knots[i + degree] - knots[i];
			qaws_scalar alpha;

			if (denom < (qaws_scalar)1e-15 && denom > (qaws_scalar)-1e-15)
				alpha = (qaws_scalar)0.0;
			else
				alpha = (t - knots[i]) / denom;

			/* Blend in homogeneous space (w P, w), then project back:
			   blending Euclidean points would change a rational curve. */
			{
				qaws_scalar wa = alpha * weights[i];
				qaws_scalar wb = ((qaws_scalar)1.0 - alpha) * weights[i - 1];
				new_weights[i] = wa + wb;
				for (d = 0; d < dim_count; d++)
				{
					new_cps[i * dim_count + d] =
						(wa * cps[i * dim_count + d] +
						 wb * cps[(i - 1) * dim_count + d]) / new_weights[i];
				}
			}
		}
	}

	*out_new_knots = new_knots;
	*out_new_knot_count = new_knot_count;
	*out_new_cps = new_cps;
	*out_new_weights = new_weights;
	*out_new_cp_count = new_cp_count;
	return QAWS_STATUS_OK;
}

static qaws_status split_nurbs(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right)
{
	qaws_nurbs_impl const *impl = (qaws_nurbs_impl const *)curve->impl;
	unsigned int degree = curve->degree;
	unsigned int dim_count = (unsigned int)curve->dimension;
	unsigned int target_mult = degree + 1;
	unsigned int current_mult;
	unsigned int insertions_needed;
	qaws_scalar *work_knots = NULL;
	qaws_scalar *work_cps = NULL;
	qaws_scalar *work_weights = NULL;
	unsigned int work_knot_count;
	unsigned int work_cp_count;
	unsigned int split_knot_idx;
	unsigned int left_cp_count;
	unsigned int right_cp_count;
	unsigned int left_knot_count;
	unsigned int right_knot_count;
	qaws_nurbs_desc left_desc;
	qaws_nurbs_desc right_desc;
	qaws_scalar *left_knots = NULL;
	qaws_scalar *left_cps = NULL;
	qaws_scalar *left_weights = NULL;
	qaws_scalar *right_knots = NULL;
	qaws_scalar *right_cps = NULL;
	qaws_scalar *right_weights = NULL;
	qaws_status status;
	unsigned int ins;
	unsigned int i, d;

	/* Copy working data */
	work_knot_count = impl->knot_count;
	work_cp_count = impl->control_point_count;

	work_knots = (qaws_scalar *)malloc((size_t)work_knot_count * sizeof(qaws_scalar));
	work_cps = (qaws_scalar *)malloc((size_t)work_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	work_weights = (qaws_scalar *)malloc((size_t)work_cp_count * sizeof(qaws_scalar));
	if (!work_knots || !work_cps || !work_weights)
	{
		free(work_knots);
		free(work_cps);
		free(work_weights);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(work_knots, impl->knots, (size_t)work_knot_count * sizeof(qaws_scalar));
	memcpy(work_cps, impl->control_points,
	       (size_t)work_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	memcpy(work_weights, impl->weights, (size_t)work_cp_count * sizeof(qaws_scalar));

	/* Insert knot t until multiplicity reaches degree+1 */
	current_mult = count_knot_multiplicity(work_knots, work_knot_count, t);
	insertions_needed = (target_mult > current_mult) ? (target_mult - current_mult) : 0;

	for (ins = 0; ins < insertions_needed; ins++)
	{
		qaws_scalar *new_knots = NULL;
		qaws_scalar *new_cps = NULL;
		qaws_scalar *new_weights = NULL;
		unsigned int new_knot_count;
		unsigned int new_cp_count;

		status = boehm_insert_knot_nurbs(
			work_knots, work_knot_count,
			work_cps, work_weights, work_cp_count,
			degree, dim_count, t,
			&new_knots, &new_knot_count,
			&new_cps, &new_weights, &new_cp_count);

		if (status != QAWS_STATUS_OK)
		{
			free(work_knots);
			free(work_cps);
			free(work_weights);
			return status;
		}

		free(work_knots);
		free(work_cps);
		free(work_weights);
		work_knots = new_knots;
		work_cps = new_cps;
		work_weights = new_weights;
		work_knot_count = new_knot_count;
		work_cp_count = new_cp_count;
	}

	/* Find the split point in the knot vector */
	split_knot_idx = find_knot_split_index(work_knots, work_knot_count, t);

	left_cp_count = split_knot_idx;
	left_knot_count = split_knot_idx + degree + 1;
	right_cp_count = work_cp_count - split_knot_idx;
	right_knot_count = work_knot_count - split_knot_idx;

	if (left_cp_count < degree + 1 || right_cp_count < degree + 1)
	{
		free(work_knots);
		free(work_cps);
		free(work_weights);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}

	/* Allocate output arrays */
	left_knots = (qaws_scalar *)malloc((size_t)left_knot_count * sizeof(qaws_scalar));
	left_cps = (qaws_scalar *)malloc((size_t)left_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	left_weights = (qaws_scalar *)malloc((size_t)left_cp_count * sizeof(qaws_scalar));
	right_knots = (qaws_scalar *)malloc((size_t)right_knot_count * sizeof(qaws_scalar));
	right_cps = (qaws_scalar *)malloc((size_t)right_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	right_weights = (qaws_scalar *)malloc((size_t)right_cp_count * sizeof(qaws_scalar));

	if (!left_knots || !left_cps || !left_weights ||
	    !right_knots || !right_cps || !right_weights)
	{
		free(left_knots);
		free(left_cps);
		free(left_weights);
		free(right_knots);
		free(right_cps);
		free(right_weights);
		free(work_knots);
		free(work_cps);
		free(work_weights);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Populate left arrays */
	for (i = 0; i < left_cp_count; i++)
	{
		for (d = 0; d < dim_count; d++)
			left_cps[i * dim_count + d] = work_cps[i * dim_count + d];
		left_weights[i] = work_weights[i];
	}
	for (i = 0; i < left_knot_count; i++)
		left_knots[i] = work_knots[i];

	/* Populate right arrays */
	for (i = 0; i < right_cp_count; i++)
	{
		for (d = 0; d < dim_count; d++)
			right_cps[i * dim_count + d] = work_cps[(split_knot_idx + i) * dim_count + d];
		right_weights[i] = work_weights[split_knot_idx + i];
	}
	for (i = 0; i < right_knot_count; i++)
		right_knots[i] = work_knots[split_knot_idx + i];

	free(work_knots);
	free(work_cps);
	free(work_weights);

	/* Create left NURBS */
	memset(&left_desc, 0, sizeof(left_desc));
	left_desc.dimension = curve->dimension;
	left_desc.degree = degree;
	left_desc.control_points = left_cps;
	left_desc.control_point_count = left_cp_count;
	left_desc.knots = left_knots;
	left_desc.knot_count = left_knot_count;
	left_desc.weights = left_weights;
	left_desc.weight_count = left_cp_count;
	left_desc.is_closed = 0;

	status = qaws_curve_create_nurbs(&left_desc, out_left);
	if (status != QAWS_STATUS_OK)
	{
		free(left_knots);
		free(left_cps);
		free(left_weights);
		free(right_knots);
		free(right_cps);
		free(right_weights);
		return status;
	}

	/* Create right NURBS */
	memset(&right_desc, 0, sizeof(right_desc));
	right_desc.dimension = curve->dimension;
	right_desc.degree = degree;
	right_desc.control_points = right_cps;
	right_desc.control_point_count = right_cp_count;
	right_desc.knots = right_knots;
	right_desc.knot_count = right_knot_count;
	right_desc.weights = right_weights;
	right_desc.weight_count = right_cp_count;
	right_desc.is_closed = 0;

	status = qaws_curve_create_nurbs(&right_desc, out_right);
	if (status != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(*out_left);
		*out_left = NULL;
		free(left_knots);
		free(left_cps);
		free(left_weights);
		free(right_knots);
		free(right_cps);
		free(right_weights);
		return status;
	}

	free(left_knots);
	free(left_cps);
	free(left_weights);
	free(right_knots);
	free(right_cps);
	free(right_weights);
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Catmull-Rom split -- boundary only, mid-span unsupported                  */
/* ========================================================================== */

static qaws_status split_catmull_rom(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int dim_count = (unsigned int)curve->dimension;
	unsigned int total_cp = impl->control_point_count;
	unsigned int span_count = curve->span_count;
	unsigned int split_span;
	qaws_scalar local_t;
	unsigned int left_cp_count;
	unsigned int right_cp_count;
	qaws_scalar *left_cp = NULL;
	qaws_scalar *right_cp = NULL;
	qaws_catmull_rom_desc left_desc;
	qaws_catmull_rom_desc right_desc;
	qaws_status status;

	/*
	 * Catmull-Rom spans are indexed from 0 to span_count-1.
	 * Each span uses 4 control points: [span-1, span, span+1, span+2]
	 * relative to the control point array (offset by 1 for the ghost point).
	 *
	 * For a clean split at a span boundary, partition the control points.
	 * For a mid-span split, we would need to convert to Bezier first,
	 * so we return UNSUPPORTED_OPERATION for that case.
	 */

	/* Find the span and local parameter */
	split_span = 0;
	local_t = (qaws_scalar)0.0;
	{
		unsigned int i;
		for (i = 0; i < span_count; i++)
		{
			if (t >= curve->span_boundaries[i] && t <= curve->span_boundaries[i + 1])
			{
				split_span = i;
				if (curve->span_boundaries[i + 1] - curve->span_boundaries[i] > (qaws_scalar)1e-15)
				{
					local_t = (t - curve->span_boundaries[i]) /
					          (curve->span_boundaries[i + 1] - curve->span_boundaries[i]);
				}
				break;
			}
		}
	}

	/* Only support splitting at span boundaries */
	if (local_t > (qaws_scalar)1e-10 && local_t < (qaws_scalar)(1.0 - 1e-10))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	/* Determine the boundary control point index */
	{
		unsigned int boundary_cp_idx;

		if (local_t < (qaws_scalar)1e-10)
			boundary_cp_idx = split_span + 1; /* CP index at span start (offset by 1 for ghost) */
		else
			boundary_cp_idx = split_span + 2; /* CP index at span end */

		/*
		 * Left curve gets control points [0 .. boundary_cp_idx + 1]
		 * (needs one extra ghost point on the right side)
		 * Right curve gets control points [boundary_cp_idx - 1 .. total_cp - 1]
		 * (needs one extra ghost point on the left side)
		 */
		left_cp_count = boundary_cp_idx + 2;
		if (left_cp_count > total_cp)
			left_cp_count = total_cp;

		if (boundary_cp_idx < 1)
			return QAWS_STATUS_INVALID_ARGUMENT;

		right_cp_count = total_cp - boundary_cp_idx + 1;
		if (right_cp_count > total_cp)
			right_cp_count = total_cp;

		/* Need at least 4 control points for a valid Catmull-Rom curve */
		if (left_cp_count < 4 || right_cp_count < 4)
			return QAWS_STATUS_INVALID_ARGUMENT;

		left_cp = (qaws_scalar *)malloc(
			(size_t)left_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
		right_cp = (qaws_scalar *)malloc(
			(size_t)right_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
		if (!left_cp || !right_cp)
		{
			free(left_cp);
			free(right_cp);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		memcpy(left_cp, impl->control_points,
		       (size_t)left_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
		memcpy(right_cp, &impl->control_points[(boundary_cp_idx - 1) * dim_count],
		       (size_t)right_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	}

	/* Create left Catmull-Rom */
	memset(&left_desc, 0, sizeof(left_desc));
	left_desc.dimension = curve->dimension;
	left_desc.control_points = left_cp;
	left_desc.control_point_count = left_cp_count;
	left_desc.parameterization = impl->parameterization;
	left_desc.closed = 0;

	status = qaws_curve_create_catmull_rom(&left_desc, out_left);
	if (status != QAWS_STATUS_OK)
	{
		free(left_cp);
		free(right_cp);
		return status;
	}

	/* Create right Catmull-Rom */
	memset(&right_desc, 0, sizeof(right_desc));
	right_desc.dimension = curve->dimension;
	right_desc.control_points = right_cp;
	right_desc.control_point_count = right_cp_count;
	right_desc.parameterization = impl->parameterization;
	right_desc.closed = 0;

	status = qaws_curve_create_catmull_rom(&right_desc, out_right);
	if (status != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(*out_left);
		*out_left = NULL;
		free(left_cp);
		free(right_cp);
		return status;
	}

	free(left_cp);
	free(right_cp);
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Trajectory split                                                          */
/* ========================================================================== */

static qaws_status split_trajectory(
	qaws_curve const *curve, qaws_scalar t,
	qaws_curve **out_left, qaws_curve **out_right)
{
	qaws_trajectory_impl const *impl =
		(qaws_trajectory_impl const *)curve->impl;
	unsigned int dim_count = (unsigned int)curve->dimension;
	unsigned int key_count = impl->key_count;
	unsigned int span_index;
	qaws_scalar split_pos[3];
	qaws_scalar split_vel[3];
	unsigned int left_count;
	unsigned int right_count;
	qaws_scalar *left_pos = NULL;
	qaws_scalar *left_times = NULL;
	qaws_scalar *left_vel = NULL;
	qaws_scalar *right_pos = NULL;
	qaws_scalar *right_times = NULL;
	qaws_scalar *right_vel = NULL;
	qaws_trajectory_desc left_desc;
	qaws_trajectory_desc right_desc;
	qaws_status status;
	unsigned int i, d;
	int at_boundary;

	/* Find the span */
	span_index = 0;
	for (i = 0; i < key_count - 1; i++)
	{
		if (t >= impl->key_times[i] && t <= impl->key_times[i + 1])
		{
			span_index = i;
			break;
		}
	}

	/* Determine if splitting at a key boundary */
	at_boundary = 0;
	{
		qaws_scalar diff = t - impl->key_times[span_index];
		if (diff < (qaws_scalar)0.0)
			diff = -diff;
		if (diff < (qaws_scalar)1e-10)
			at_boundary = 1;
	}

	if (at_boundary)
	{
		if (span_index < 1)
			return QAWS_STATUS_INVALID_ARGUMENT;

		left_count = span_index + 1;
		right_count = key_count - span_index;
	}
	else
	{
		left_count = span_index + 2;
		right_count = key_count - span_index;
	}

	/* Evaluate position and velocity at the split time */
	memset(split_pos, 0, sizeof(split_pos));
	memset(split_vel, 0, sizeof(split_vel));

	if (dim_count == 2)
	{
		qaws_eval_result_2d result;
		status = qaws_curve_evaluate_2d(curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);
		if (status != QAWS_STATUS_OK)
			return status;
		split_pos[0] = result.position.x;
		split_pos[1] = result.position.y;
		split_vel[0] = result.d1.x;
		split_vel[1] = result.d1.y;
	}
	else
	{
		qaws_eval_result_3d result;
		status = qaws_curve_evaluate_3d(curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);
		if (status != QAWS_STATUS_OK)
			return status;
		split_pos[0] = result.position.x;
		split_pos[1] = result.position.y;
		split_pos[2] = result.position.z;
		split_vel[0] = result.d1.x;
		split_vel[1] = result.d1.y;
		split_vel[2] = result.d1.z;
	}

	/* Allocate left arrays */
	left_pos = (qaws_scalar *)malloc((size_t)left_count * (size_t)dim_count * sizeof(qaws_scalar));
	left_times = (qaws_scalar *)malloc((size_t)left_count * sizeof(qaws_scalar));
	left_vel = (qaws_scalar *)malloc((size_t)left_count * (size_t)dim_count * sizeof(qaws_scalar));
	if (!left_pos || !left_times || !left_vel)
	{
		free(left_pos);
		free(left_times);
		free(left_vel);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Allocate right arrays */
	right_pos = (qaws_scalar *)malloc((size_t)right_count * (size_t)dim_count * sizeof(qaws_scalar));
	right_times = (qaws_scalar *)malloc((size_t)right_count * sizeof(qaws_scalar));
	right_vel = (qaws_scalar *)malloc((size_t)right_count * (size_t)dim_count * sizeof(qaws_scalar));
	if (!right_pos || !right_times || !right_vel)
	{
		free(left_pos);
		free(left_times);
		free(left_vel);
		free(right_pos);
		free(right_times);
		free(right_vel);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	if (at_boundary)
	{
		/* Boundary split: partition at span_index */
		for (i = 0; i < left_count; i++)
		{
			left_times[i] = impl->key_times[i];
			for (d = 0; d < dim_count; d++)
			{
				left_pos[i * dim_count + d] = impl->key_positions[i * dim_count + d];
				left_vel[i * dim_count + d] = impl->key_velocities[i * dim_count + d];
			}
		}

		for (i = 0; i < right_count; i++)
		{
			right_times[i] = impl->key_times[span_index + i];
			for (d = 0; d < dim_count; d++)
			{
				right_pos[i * dim_count + d] = impl->key_positions[(span_index + i) * dim_count + d];
				right_vel[i * dim_count + d] = impl->key_velocities[(span_index + i) * dim_count + d];
			}
		}
	}
	else
	{
		/* Within-span split */
		/* Left: keys[0..span_index] + split point */
		for (i = 0; i <= span_index; i++)
		{
			left_times[i] = impl->key_times[i];
			for (d = 0; d < dim_count; d++)
			{
				left_pos[i * dim_count + d] = impl->key_positions[i * dim_count + d];
				left_vel[i * dim_count + d] = impl->key_velocities[i * dim_count + d];
			}
		}
		left_times[span_index + 1] = t;
		for (d = 0; d < dim_count; d++)
		{
			left_pos[(span_index + 1) * dim_count + d] = split_pos[d];
			left_vel[(span_index + 1) * dim_count + d] = split_vel[d];
		}

		/* Right: split point + keys[span_index+1..n-1] */
		right_times[0] = t;
		for (d = 0; d < dim_count; d++)
		{
			right_pos[d] = split_pos[d];
			right_vel[d] = split_vel[d];
		}
		for (i = span_index + 1; i < key_count; i++)
		{
			unsigned int ri = i - span_index;
			right_times[ri] = impl->key_times[i];
			for (d = 0; d < dim_count; d++)
			{
				right_pos[ri * dim_count + d] = impl->key_positions[i * dim_count + d];
				right_vel[ri * dim_count + d] = impl->key_velocities[i * dim_count + d];
			}
		}
	}

	/* Create left trajectory */
	memset(&left_desc, 0, sizeof(left_desc));
	left_desc.dimension = curve->dimension;
	left_desc.degree = 3;
	left_desc.key_positions = left_pos;
	left_desc.key_count = left_count;
	left_desc.key_times = left_times;
	left_desc.key_time_count = left_count;
	left_desc.key_velocities = left_vel;
	left_desc.key_velocity_count = left_count;
	left_desc.key_accelerations = NULL;
	left_desc.key_acceleration_count = 0;
	left_desc.closed = 0;

	status = qaws_curve_create_trajectory(&left_desc, out_left);
	if (status != QAWS_STATUS_OK)
	{
		free(left_pos);
		free(left_times);
		free(left_vel);
		free(right_pos);
		free(right_times);
		free(right_vel);
		return status;
	}

	/* Create right trajectory */
	memset(&right_desc, 0, sizeof(right_desc));
	right_desc.dimension = curve->dimension;
	right_desc.degree = 3;
	right_desc.key_positions = right_pos;
	right_desc.key_count = right_count;
	right_desc.key_times = right_times;
	right_desc.key_time_count = right_count;
	right_desc.key_velocities = right_vel;
	right_desc.key_velocity_count = right_count;
	right_desc.key_accelerations = NULL;
	right_desc.key_acceleration_count = 0;
	right_desc.closed = 0;

	status = qaws_curve_create_trajectory(&right_desc, out_right);
	if (status != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(*out_left);
		*out_left = NULL;
		free(left_pos);
		free(left_times);
		free(left_vel);
		free(right_pos);
		free(right_times);
		free(right_vel);
		return status;
	}

	free(left_pos);
	free(left_times);
	free(left_vel);
	free(right_pos);
	free(right_times);
	free(right_vel);
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Hermite join                                                              */
/* ========================================================================== */

static qaws_status join_hermite(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined)
{
	qaws_hermite_impl const *impl_a = (qaws_hermite_impl const *)a->impl;
	qaws_hermite_impl const *impl_b = (qaws_hermite_impl const *)b->impl;
	unsigned int dim_count = (unsigned int)a->dimension;
	unsigned int count_a = impl_a->point_count;
	unsigned int count_b = impl_b->point_count;
	unsigned int joined_count;
	qaws_scalar *joined_pts = NULL;
	qaws_scalar *joined_tan = NULL;
	size_t joined_size;
	qaws_hermite_desc desc;
	qaws_status status;
	unsigned int i, d;

	/*
	 * Join by merging the shared endpoint:
	 * Use A's last point as the shared point.
	 * joined_count = count_a + count_b - 1
	 */
	joined_count = count_a + count_b - 1;

	joined_size = (size_t)joined_count * (size_t)dim_count * sizeof(qaws_scalar);
	joined_pts = (qaws_scalar *)malloc(joined_size);
	joined_tan = (qaws_scalar *)malloc(joined_size);
	if (!joined_pts || !joined_tan)
	{
		free(joined_pts);
		free(joined_tan);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy A's points and tangents */
	memcpy(joined_pts, impl_a->points,
	       (size_t)count_a * (size_t)dim_count * sizeof(qaws_scalar));
	memcpy(joined_tan, impl_a->tangents,
	       (size_t)count_a * (size_t)dim_count * sizeof(qaws_scalar));

	/* Copy B's points and tangents (skip B's first point -- it's the shared one) */
	for (i = 1; i < count_b; i++)
	{
		unsigned int dst = count_a + i - 1;
		for (d = 0; d < dim_count; d++)
		{
			joined_pts[dst * dim_count + d] = impl_b->points[i * dim_count + d];
			joined_tan[dst * dim_count + d] = impl_b->tangents[i * dim_count + d];
		}
	}

	/* Create the joined Hermite curve */
	memset(&desc, 0, sizeof(desc));
	desc.dimension = a->dimension;
	desc.degree = 3;
	desc.points = joined_pts;
	desc.derivatives = joined_tan;
	desc.point_count = joined_count;
	desc.derivative_count = joined_count;

	status = qaws_curve_create_hermite(&desc, out_joined);

	free(joined_pts);
	free(joined_tan);
	return status;
}

/* ========================================================================== */
/*  B-Spline join                                                             */
/* ========================================================================== */

static qaws_status join_bspline(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined)
{
	qaws_bspline_impl const *impl_a = (qaws_bspline_impl const *)a->impl;
	qaws_bspline_impl const *impl_b = (qaws_bspline_impl const *)b->impl;
	unsigned int degree = a->degree;
	unsigned int dim_count = (unsigned int)a->dimension;
	qaws_scalar junction;
	qaws_scalar b_shift;
	unsigned int joined_cp_count;
	unsigned int joined_knot_count;
	qaws_scalar *joined_cps = NULL;
	qaws_scalar *joined_knots = NULL;
	unsigned int skip_b_knots;
	qaws_bspline_desc desc;
	qaws_status status;
	unsigned int i, d;

	/* Both curves must have the same degree */
	if (b->degree != degree)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/*
	 * Join strategy:
	 *   junction = A's max parameter
	 *   Shift B's knots so B's min parameter aligns to junction
	 *   Concatenate knots (skip B's first degree+1 clamped knots)
	 *   Concatenate CPs (skip B's first CP which is the shared point)
	 */
	junction = a->parameter_range.max_value;
	b_shift = junction - b->parameter_range.min_value;

	/* CPs: A's CPs + B's CPs minus one shared point */
	joined_cp_count = impl_a->control_point_count + impl_b->control_point_count - 1;

	/*
	 * Knots: A's knots (minus last) + B's knots shifted (skip first degree+1).
	 * Dropping one knot from A's clamped end reduces the junction multiplicity
	 * from degree+1 to degree, giving C0 continuity at the shared point.
	 * Result: knot_count = cp_count + degree + 1.
	 */
	skip_b_knots = degree + 1;
	if (skip_b_knots > impl_b->knot_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	joined_knot_count = (impl_a->knot_count - 1) + (impl_b->knot_count - skip_b_knots);

	joined_cps = (qaws_scalar *)malloc(
		(size_t)joined_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	joined_knots = (qaws_scalar *)malloc(
		(size_t)joined_knot_count * sizeof(qaws_scalar));
	if (!joined_cps || !joined_knots)
	{
		free(joined_cps);
		free(joined_knots);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy A's control points */
	memcpy(joined_cps, impl_a->control_points,
	       (size_t)impl_a->control_point_count * (size_t)dim_count * sizeof(qaws_scalar));

	/* Copy B's control points (skip first, which is the shared point) */
	for (i = 1; i < impl_b->control_point_count; i++)
	{
		unsigned int dst = impl_a->control_point_count + i - 1;
		for (d = 0; d < dim_count; d++)
			joined_cps[dst * dim_count + d] = impl_b->control_points[i * dim_count + d];
	}

	/* Copy A's knots (all except the last) */
	memcpy(joined_knots, impl_a->knots,
	       (size_t)(impl_a->knot_count - 1) * sizeof(qaws_scalar));

	/* Copy B's knots (shifted), skipping first degree+1 */
	for (i = skip_b_knots; i < impl_b->knot_count; i++)
	{
		unsigned int dst = (impl_a->knot_count - 1) + i - skip_b_knots;
		joined_knots[dst] = impl_b->knots[i] + b_shift;
	}

	/* Create the joined B-Spline */
	memset(&desc, 0, sizeof(desc));
	desc.dimension = a->dimension;
	desc.degree = degree;
	desc.control_points = joined_cps;
	desc.control_point_count = joined_cp_count;
	desc.knots = joined_knots;
	desc.knot_count = joined_knot_count;
	desc.is_uniform = 0;
	desc.is_closed = 0;

	status = qaws_curve_create_bspline(&desc, out_joined);

	free(joined_cps);
	free(joined_knots);
	return status;
}

/* ========================================================================== */
/*  NURBS join                                                                */
/* ========================================================================== */

static qaws_status join_nurbs(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined)
{
	qaws_nurbs_impl const *impl_a = (qaws_nurbs_impl const *)a->impl;
	qaws_nurbs_impl const *impl_b = (qaws_nurbs_impl const *)b->impl;
	unsigned int degree = a->degree;
	unsigned int dim_count = (unsigned int)a->dimension;
	qaws_scalar junction;
	qaws_scalar b_shift;
	unsigned int joined_cp_count;
	unsigned int joined_knot_count;
	unsigned int skip_b_knots;
	qaws_scalar *joined_cps = NULL;
	qaws_scalar *joined_knots = NULL;
	qaws_scalar *joined_weights = NULL;
	qaws_nurbs_desc desc;
	qaws_status status;
	unsigned int i, d;

	if (b->degree != degree)
		return QAWS_STATUS_INVALID_ARGUMENT;

	junction = a->parameter_range.max_value;
	b_shift = junction - b->parameter_range.min_value;

	joined_cp_count = impl_a->control_point_count + impl_b->control_point_count - 1;

	skip_b_knots = degree + 1;
	if (skip_b_knots > impl_b->knot_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	joined_knot_count = (impl_a->knot_count - 1) + (impl_b->knot_count - skip_b_knots);

	joined_cps = (qaws_scalar *)malloc(
		(size_t)joined_cp_count * (size_t)dim_count * sizeof(qaws_scalar));
	joined_knots = (qaws_scalar *)malloc(
		(size_t)joined_knot_count * sizeof(qaws_scalar));
	joined_weights = (qaws_scalar *)malloc(
		(size_t)joined_cp_count * sizeof(qaws_scalar));
	if (!joined_cps || !joined_knots || !joined_weights)
	{
		free(joined_cps);
		free(joined_knots);
		free(joined_weights);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy A's control points and weights */
	memcpy(joined_cps, impl_a->control_points,
	       (size_t)impl_a->control_point_count * (size_t)dim_count * sizeof(qaws_scalar));
	memcpy(joined_weights, impl_a->weights,
	       (size_t)impl_a->control_point_count * sizeof(qaws_scalar));

	/* Copy B's control points and weights (skip shared first point) */
	for (i = 1; i < impl_b->control_point_count; i++)
	{
		unsigned int dst = impl_a->control_point_count + i - 1;
		for (d = 0; d < dim_count; d++)
			joined_cps[dst * dim_count + d] = impl_b->control_points[i * dim_count + d];
		joined_weights[dst] = impl_b->weights[i];
	}

	/* Copy A's knots (all except the last) */
	memcpy(joined_knots, impl_a->knots,
	       (size_t)(impl_a->knot_count - 1) * sizeof(qaws_scalar));

	/* Copy B's shifted knots (skip first degree+1) */
	for (i = skip_b_knots; i < impl_b->knot_count; i++)
	{
		unsigned int dst = (impl_a->knot_count - 1) + i - skip_b_knots;
		joined_knots[dst] = impl_b->knots[i] + b_shift;
	}

	/* Create the joined NURBS curve */
	memset(&desc, 0, sizeof(desc));
	desc.dimension = a->dimension;
	desc.degree = degree;
	desc.control_points = joined_cps;
	desc.control_point_count = joined_cp_count;
	desc.knots = joined_knots;
	desc.knot_count = joined_knot_count;
	desc.weights = joined_weights;
	desc.weight_count = joined_cp_count;
	desc.is_closed = 0;

	status = qaws_curve_create_nurbs(&desc, out_joined);

	free(joined_cps);
	free(joined_knots);
	free(joined_weights);
	return status;
}

/* ========================================================================== */
/*  Catmull-Rom join                                                          */
/* ========================================================================== */

static qaws_status join_catmull_rom(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined)
{
	qaws_catmull_rom_impl const *impl_a =
		(qaws_catmull_rom_impl const *)a->impl;
	qaws_catmull_rom_impl const *impl_b =
		(qaws_catmull_rom_impl const *)b->impl;
	unsigned int dim_count = (unsigned int)a->dimension;
	unsigned int count_a = impl_a->control_point_count;
	unsigned int count_b = impl_b->control_point_count;
	unsigned int joined_count;
	qaws_scalar *joined_cp = NULL;
	qaws_catmull_rom_desc desc;
	qaws_status status;
	unsigned int i, d;

	/*
	 * Concatenate control points, merging the shared boundary point.
	 * A's last point should coincide with B's first point.
	 * We use A's last point as the shared point.
	 */
	joined_count = count_a + count_b - 1;

	if (joined_count < 4)
		return QAWS_STATUS_INVALID_ARGUMENT;

	joined_cp = (qaws_scalar *)malloc(
		(size_t)joined_count * (size_t)dim_count * sizeof(qaws_scalar));
	if (!joined_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Copy A's control points */
	memcpy(joined_cp, impl_a->control_points,
	       (size_t)count_a * (size_t)dim_count * sizeof(qaws_scalar));

	/* Copy B's control points (skip the first shared point) */
	for (i = 1; i < count_b; i++)
	{
		unsigned int dst = count_a + i - 1;
		for (d = 0; d < dim_count; d++)
			joined_cp[dst * dim_count + d] = impl_b->control_points[i * dim_count + d];
	}

	memset(&desc, 0, sizeof(desc));
	desc.dimension = a->dimension;
	desc.control_points = joined_cp;
	desc.control_point_count = joined_count;
	desc.parameterization = impl_a->parameterization;
	desc.closed = 0;

	status = qaws_curve_create_catmull_rom(&desc, out_joined);

	free(joined_cp);
	return status;
}

/* ========================================================================== */
/*  Trajectory join                                                           */
/* ========================================================================== */

static qaws_status join_trajectory(
	qaws_curve const *a, qaws_curve const *b,
	qaws_curve **out_joined)
{
	qaws_trajectory_impl const *impl_a =
		(qaws_trajectory_impl const *)a->impl;
	qaws_trajectory_impl const *impl_b =
		(qaws_trajectory_impl const *)b->impl;
	unsigned int dim_count = (unsigned int)a->dimension;
	unsigned int count_a = impl_a->key_count;
	unsigned int count_b = impl_b->key_count;
	unsigned int joined_count;
	qaws_scalar time_shift;
	qaws_scalar *joined_pos = NULL;
	qaws_scalar *joined_times = NULL;
	qaws_scalar *joined_vel = NULL;
	qaws_trajectory_desc desc;
	qaws_status status;
	unsigned int i, d;

	/*
	 * Join by merging the shared key point.
	 * B's times are shifted so that B's start time equals A's end time.
	 */
	joined_count = count_a + count_b - 1;

	time_shift = impl_a->key_times[count_a - 1] - impl_b->key_times[0];

	joined_pos = (qaws_scalar *)malloc(
		(size_t)joined_count * (size_t)dim_count * sizeof(qaws_scalar));
	joined_times = (qaws_scalar *)malloc(
		(size_t)joined_count * sizeof(qaws_scalar));
	joined_vel = (qaws_scalar *)malloc(
		(size_t)joined_count * (size_t)dim_count * sizeof(qaws_scalar));
	if (!joined_pos || !joined_times || !joined_vel)
	{
		free(joined_pos);
		free(joined_times);
		free(joined_vel);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Copy A's data */
	memcpy(joined_pos, impl_a->key_positions,
	       (size_t)count_a * (size_t)dim_count * sizeof(qaws_scalar));
	memcpy(joined_times, impl_a->key_times,
	       (size_t)count_a * sizeof(qaws_scalar));
	memcpy(joined_vel, impl_a->key_velocities,
	       (size_t)count_a * (size_t)dim_count * sizeof(qaws_scalar));

	/* Copy B's data (skip the first shared key, shift times) */
	for (i = 1; i < count_b; i++)
	{
		unsigned int dst = count_a + i - 1;
		joined_times[dst] = impl_b->key_times[i] + time_shift;
		for (d = 0; d < dim_count; d++)
		{
			joined_pos[dst * dim_count + d] = impl_b->key_positions[i * dim_count + d];
			joined_vel[dst * dim_count + d] = impl_b->key_velocities[i * dim_count + d];
		}
	}

	/* Create the joined trajectory */
	memset(&desc, 0, sizeof(desc));
	desc.dimension = a->dimension;
	desc.degree = 3;
	desc.key_positions = joined_pos;
	desc.key_count = joined_count;
	desc.key_times = joined_times;
	desc.key_time_count = joined_count;
	desc.key_velocities = joined_vel;
	desc.key_velocity_count = joined_count;
	desc.key_accelerations = NULL;
	desc.key_acceleration_count = 0;
	desc.closed = 0;

	status = qaws_curve_create_trajectory(&desc, out_joined);

	free(joined_pos);
	free(joined_times);
	free(joined_vel);
	return status;
}

/* ========================================================================== */
/*  Curve offset (parallel curve)                                              */
/* ========================================================================== */

/*
 * Algorithm — "sweep" mental model:
 *
 * The offset of curve C(t) at signed distance d is O(t) = C(t) + d*N(t).
 * The forward factor f(t) = 1 - kappa(t)*d tells us when the offset moves
 * forward (f>0) or backward (f<0, cusp region).
 *
 * At cusps the raw offset polyline forms a backward loop that self-intersects.
 * We detect these self-intersections and trim the loop. The intersection point
 * lies on both forward branches so it IS at distance |d| from the curve.
 *
 * Optional trim step: for each output point, verify its distance to the
 * nearest curve point equals |d| (within tolerance). Points too close to
 * the curve are inside the swept region and get removed.
 */

#define OFFSET_SAMPLE_COUNT   512
#define OFFSET_MIN_PIECE_PTS  4
#define OFFSET_OUTPUT_PTS     48
#define OFFSET_TRIM_TOLERANCE ((qaws_scalar)0.10) /* 10% relative */

/* Segment-segment intersection for polyline self-crossing detection. */
static int offset_seg_isect(
	qaws_vec2 a0, qaws_vec2 a1, qaws_vec2 b0, qaws_vec2 b1,
	qaws_scalar *out_t, qaws_scalar *out_u)
{
	qaws_scalar dax = a1.x - a0.x, day = a1.y - a0.y;
	qaws_scalar dbx = b1.x - b0.x, dby = b1.y - b0.y;
	qaws_scalar denom = dax * dby - day * dbx;
	qaws_scalar abx, aby, t_val, u_val;

	if (denom > -(qaws_scalar)1e-12 && denom < (qaws_scalar)1e-12)
		return 0;

	abx = b0.x - a0.x;
	aby = b0.y - a0.y;
	t_val = (abx * dby - aby * dbx) / denom;
	u_val = (abx * day - aby * dax) / denom;

	if (t_val < (qaws_scalar)0 || t_val > (qaws_scalar)1) return 0;
	if (u_val < (qaws_scalar)0 || u_val > (qaws_scalar)1) return 0;

	*out_t = t_val;
	*out_u = u_val;
	return 1;
}

/* Truncate an offset curve to a normalised sub-range [frac0,frac1]
 * (0 = start, 1 = end) of the curve's actual parameter span.
 * Evaluates the old curve at n_trunc points within the mapped range,
 * builds a new open Catmull-Rom, and replaces out_curves[idx].       */
static void offset_trunc_curve(
	qaws_curve **out_curves, unsigned int idx,
	qaws_scalar frac0, qaws_scalar frac1)
{
	unsigned int n_trunc = 48, kk;
	qaws_scalar *flat;
	qaws_eval_result_2d er;
	qaws_catmull_rom_desc desc;
	qaws_curve *nc = NULL;
	qaws_status st;
	qaws_range rng = qaws_curve_get_parameter_range(out_curves[idx]);
	qaws_scalar rmin = rng.min_value, rlen = rng.max_value - rng.min_value;
	qaws_scalar t0 = rmin + frac0 * rlen;
	qaws_scalar t1 = rmin + frac1 * rlen;

	flat = (qaws_scalar *)malloc((n_trunc + 2) * 2 * sizeof(qaws_scalar));
	if (!flat) return;

	for (kk = 0; kk < n_trunc; ++kk)
	{
		qaws_scalar tt = t0 + (t1 - t0) * (qaws_scalar)kk / (qaws_scalar)(n_trunc - 1);
		qaws_curve_evaluate_2d(out_curves[idx], tt, QAWS_EVAL_FLAG_POSITION, &er);
		flat[(kk + 1) * 2]     = er.position.x;
		flat[(kk + 1) * 2 + 1] = er.position.y;
	}
	/* Phantom endpoint duplication */
	flat[0] = flat[2];  flat[1] = flat[3];
	flat[(n_trunc + 1) * 2]     = flat[n_trunc * 2];
	flat[(n_trunc + 1) * 2 + 1] = flat[n_trunc * 2 + 1];

	memset(&desc, 0, sizeof(desc));
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = flat;
	desc.control_point_count = n_trunc + 2;
	desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	desc.closed = 0;

	st = qaws_curve_create_catmull_rom(&desc, &nc);
	free(flat);

	if (st == QAWS_STATUS_OK && nc)
	{
		qaws_curve_destroy(out_curves[idx]);
		out_curves[idx] = nc;
	}
}

qaws_status qaws_curve_offset_2d(
	qaws_curve const *curve,
	qaws_scalar distance,
	int trim,
	qaws_curve **out_curves,
	unsigned int curve_capacity,
	unsigned int *out_count)
{
	unsigned int i, j;
	qaws_scalar range_min, range_max, t;
	qaws_eval_result_2d result;
	qaws_status status;

	unsigned int n_pts = OFFSET_SAMPLE_COUNT;
	qaws_vec2 *offset_pts = NULL;
	qaws_scalar *fwd = NULL;          /* forward factor per sample */

	/* Self-intersection bookkeeping */
	typedef struct { unsigned int a; qaws_scalar fa;
	                 unsigned int b; qaws_scalar fb; } xing_t;
	xing_t *xings = NULL;
	unsigned int n_xings = 0, xing_cap = 64;
	int *removed = NULL;

	unsigned int curve_count = 0;

	/* ---- Validate ---- */
	if (!curve || !out_curves || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	*out_count = 0;
	if (distance == (qaws_scalar)0)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* ---- Step 1: Sample offset + forward factor ---- */
	offset_pts = (qaws_vec2 *)malloc(n_pts * sizeof(qaws_vec2));
	fwd        = (qaws_scalar *)malloc(n_pts * sizeof(qaws_scalar));
	if (!offset_pts || !fwd) { free(offset_pts); free(fwd); return QAWS_STATUS_ALLOCATION_FAILURE; }

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	for (i = 0; i < n_pts; ++i)
	{
		qaws_scalar spd, cross, kappa, nx, ny;
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(n_pts - 1);

		status = qaws_curve_evaluate_2d(curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
			&result);
		if (status != QAWS_STATUS_OK) { free(offset_pts); free(fwd); return status; }

		spd = (qaws_scalar)sqrt((double)(
			result.d1.x * result.d1.x + result.d1.y * result.d1.y));

		if (spd < (qaws_scalar)1e-12)
		{
			fwd[i] = (qaws_scalar)1;
			if (i > 0) offset_pts[i] = offset_pts[i - 1];
			else { offset_pts[i] = result.position; offset_pts[i].x += distance; }
			continue;
		}

		cross = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
		kappa = cross / (spd * spd * spd);
		fwd[i] = (qaws_scalar)1 - kappa * distance;

		nx = -result.d1.y / spd;
		ny =  result.d1.x / spd;
		offset_pts[i].x = result.position.x + distance * nx;
		offset_pts[i].y = result.position.y + distance * ny;
	}

	/* ---- Step 2: Find self-intersections in offset polyline ---- */
	xings = (xing_t *)malloc(xing_cap * sizeof(xing_t));
	removed = (int *)calloc(n_pts, sizeof(int));
	if (!xings || !removed)
	{
		free(offset_pts); free(fwd); free(xings); free(removed);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	{
		/* min_gap: large enough to skip trivial adjacent noise,
		   small enough to catch cusp loops. */
		unsigned int min_gap = 8;
		for (i = 0; i + 1 < n_pts; ++i)
		{
			for (j = i + min_gap; j + 1 < n_pts; ++j)
			{
				qaws_scalar ft, fu;
				if (offset_seg_isect(offset_pts[i], offset_pts[i + 1],
					offset_pts[j], offset_pts[j + 1], &ft, &fu))
				{
					/* Only consider crossings that enclose backward samples
					   (cusp loops). Check if any sample in [i+1..j] has f<0. */
					int has_backward = 0;
					unsigned int k;
					for (k = i + 1; k <= j; ++k)
						if (fwd[k] < (qaws_scalar)0) { has_backward = 1; break; }
					if (!has_backward) continue;

					if (n_xings >= xing_cap)
					{
						xing_t *tmp;
						xing_cap *= 2;
						tmp = (xing_t *)realloc(xings, xing_cap * sizeof(xing_t));
						if (!tmp) break;
						xings = tmp;
					}
					xings[n_xings].a = i;
					xings[n_xings].fa = ft;
					xings[n_xings].b = j;
					xings[n_xings].fb = fu;
					n_xings++;
				}
			}
		}
	}

	/* ---- Step 3: Trim cusp loops ---- */
	for (i = 0; i < n_xings; ++i)
	{
		unsigned int a = xings[i].a, b = xings[i].b;
		qaws_scalar fa = xings[i].fa;
		qaws_vec2 cp;
		cp.x = offset_pts[a].x + fa * (offset_pts[a + 1].x - offset_pts[a].x);
		cp.y = offset_pts[a].y + fa * (offset_pts[a + 1].y - offset_pts[a].y);

		/* Mark loop interior for removal */
		for (j = a + 1; j <= b; ++j)
			removed[j] = 1;

		/* Bridge: place crossing point just after the trim */
		if (a + 1 < n_pts) { offset_pts[a + 1] = cp; removed[a + 1] = 0; }
	}

	free(xings);
	free(fwd);

	/* ---- Step 4: Distance-based trim (trim >= 1) ---- */
	if (trim >= 1)
	{
		qaws_scalar abs_d = distance < 0 ? -distance : distance;
		qaws_scalar lo = abs_d * ((qaws_scalar)1 - OFFSET_TRIM_TOLERANCE);

		for (i = 0; i < n_pts; ++i)
		{
			qaws_scalar closest_t, dx, dy, dist;
			if (removed[i]) continue;
			if (i == 0 || i == n_pts - 1) continue;

			status = qaws_curve_find_closest_parameter_2d(curve, offset_pts[i], &closest_t);
			if (status != QAWS_STATUS_OK) continue;
			qaws_curve_evaluate_2d(curve, closest_t, QAWS_EVAL_FLAG_POSITION, &result);
			dx = offset_pts[i].x - result.position.x;
			dy = offset_pts[i].y - result.position.y;
			dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));
			if (dist < lo)
				removed[i] = 1;
		}
	}

	/* ---- Step 5: Collect surviving points, split at gaps, fit curves ---- */
	/* For trim==2, each run is checked for self-intersection.
	 * If found, the run is split into two closed loops at the
	 * crossing point.  This is a true post-process of trim==1
	 * so trim==2 never produces more geometry than trim==1. */
	{
		unsigned int *surv_idx;
		unsigned int surv_len = 0;

		surv_idx = (unsigned int *)malloc(n_pts * sizeof(unsigned int));
		if (!surv_idx) { free(offset_pts); free(removed); return QAWS_STATUS_ALLOCATION_FAILURE; }

		for (i = 0; i < n_pts; ++i)
			if (!removed[i])
				surv_idx[surv_len++] = i;

		/* Identify run boundaries.  is_real_seg[k]=1 if surv_idx[k] and
		 * surv_idx[k+1] are within the same contiguous run (gap <= 3). */
		{
			int *is_real_seg = (int *)calloc(surv_len, sizeof(int));

			if (is_real_seg && surv_len > 1)
			{
				for (i = 0; i + 1 < surv_len; ++i)
					is_real_seg[i] = (surv_idx[i + 1] - surv_idx[i] <= 3) ? 1 : 0;
			}

			/* --- trim==2: remove intra-run self-intersection loops --- */
			if (trim == 2 && is_real_seg && surv_len >= 16)
			{
				int *surv_rm = (int *)calloc(surv_len, sizeof(int));
				if (surv_rm)
				{
					unsigned int run_s = 0;
					while (run_s < surv_len)
					{
						unsigned int run_e = run_s + 1, rl;
						while (run_e < surv_len && is_real_seg[run_e - 1])
							run_e++;
						rl = run_e - run_s;
						if (rl >= 16)
						{
							unsigned int ci = 0, cj = 0;
							qaws_scalar c_ft = 0;
							int found = 0;
							for (ci = run_s; ci + 1 < run_e && !found; ++ci)
								for (cj = ci + 8; cj + 1 < run_e; ++cj)
								{
									qaws_scalar ft, fu;
									if (offset_seg_isect(
										offset_pts[surv_idx[ci]], offset_pts[surv_idx[ci + 1]],
										offset_pts[surv_idx[cj]], offset_pts[surv_idx[cj + 1]],
										&ft, &fu))
									{ c_ft = ft; found = 1; break; }
								}
							if (found)
							{
								qaws_vec2 xp;
								unsigned int loop_n = cj - ci - 1;
								unsigned int tail1_n = ci - run_s;
								unsigned int tail2_n = run_e - 1 - cj;
								unsigned int m;
								xp.x = offset_pts[surv_idx[ci]].x + c_ft * (offset_pts[surv_idx[ci + 1]].x - offset_pts[surv_idx[ci]].x);
								xp.y = offset_pts[surv_idx[ci]].y + c_ft * (offset_pts[surv_idx[ci + 1]].y - offset_pts[surv_idx[ci]].y);
								if (loop_n >= tail1_n && loop_n >= tail2_n)
								{
									/* Loop is largest — bridge at crossing, keep ci..cj+1 */
									offset_pts[surv_idx[ci]] = xp;
									if (cj + 1 < run_e) offset_pts[surv_idx[cj + 1]] = xp;
									for (m = run_s; m < ci; ++m) surv_rm[m] = 1;
									for (m = cj + 2; m < run_e; ++m) surv_rm[m] = 1;
								}
								else if (tail1_n >= loop_n && tail1_n >= tail2_n)
								{
									/* First tail largest — keep run_s..ci+1, bridge ci+1=xp */
									if (ci + 1 < run_e) offset_pts[surv_idx[ci + 1]] = xp;
									for (m = ci + 2; m < run_e; ++m) surv_rm[m] = 1;
								}
								else
								{
									/* Second tail largest — keep cj..run_e-1, bridge cj=xp */
									offset_pts[surv_idx[cj]] = xp;
									for (m = run_s; m < cj; ++m) surv_rm[m] = 1;
								}
							}
						}
						run_s = run_e;
					}
					/* Compact surv_idx, removing marked entries */
					{
						unsigned int new_len = 0;
						for (i = 0; i < surv_len; ++i)
							if (!surv_rm[i])
								surv_idx[new_len++] = surv_idx[i];
						surv_len = new_len;
					}
					/* Rebuild is_real_seg for modified polyline */
					free(is_real_seg);
					is_real_seg = (int *)calloc(surv_len, sizeof(int));
					if (is_real_seg && surv_len > 1)
						for (i = 0; i + 1 < surv_len; ++i)
							is_real_seg[i] = (surv_idx[i + 1] - surv_idx[i] <= 3) ? 1 : 0;
					free(surv_rm);
				}
			}

			/* --- trim==2: inter-run crossing cleanup --- */
			if (trim == 2 && is_real_seg && surv_len >= 8)
			{
				int *surv_rm2 = (int *)calloc(surv_len, sizeof(int));
				if (surv_rm2)
				{
					unsigned int run_bd[128]; /* [start,end) pairs, max 64 runs */
					unsigned int n_runs = 0, ra, rb;

					/* Build run list */
					{
						unsigned int rs2 = 0;
						while (rs2 < surv_len && n_runs < 64)
						{
							unsigned int re2 = rs2 + 1;
							while (re2 < surv_len && is_real_seg[re2 - 1])
								re2++;
							if (re2 - rs2 >= OFFSET_MIN_PIECE_PTS)
							{
								run_bd[n_runs * 2] = rs2;
								run_bd[n_runs * 2 + 1] = re2;
								n_runs++;
							}
							rs2 = re2;
						}
					}

					/* Check each pair of runs for crossings */
					for (ra = 0; ra + 1 < n_runs; ++ra)
					{
						for (rb = ra + 1; rb < n_runs; ++rb)
						{
							int xfound = 0;
							unsigned int xai = 0, xbj = 0, ai2, bj2, m;
							qaws_scalar xft = 0;

							for (ai2 = run_bd[ra * 2]; ai2 + 1 < run_bd[ra * 2 + 1] && !xfound; ++ai2)
								for (bj2 = run_bd[rb * 2]; bj2 + 1 < run_bd[rb * 2 + 1]; ++bj2)
								{
									qaws_scalar ft, fu;
									if (offset_seg_isect(
										offset_pts[surv_idx[ai2]], offset_pts[surv_idx[ai2 + 1]],
										offset_pts[surv_idx[bj2]], offset_pts[surv_idx[bj2 + 1]],
										&ft, &fu))
									{
										xft = ft; xai = ai2; xbj = bj2;
										xfound = 1; break;
									}
								}

							if (xfound)
							{
								qaws_vec2 xp;
								unsigned int a_bef = xai - run_bd[ra * 2];
								unsigned int a_aft = run_bd[ra * 2 + 1] - 1 - xai;
								unsigned int b_bef = xbj - run_bd[rb * 2];
								unsigned int b_aft = run_bd[rb * 2 + 1] - 1 - xbj;

								xp.x = offset_pts[surv_idx[xai]].x + xft * (offset_pts[surv_idx[xai + 1]].x - offset_pts[surv_idx[xai]].x);
								xp.y = offset_pts[surv_idx[xai]].y + xft * (offset_pts[surv_idx[xai + 1]].y - offset_pts[surv_idx[xai]].y);

								/* Keep the combination that preserves the most
								 * geometry.  Only two combos avoid re-crossing:
								 *   A_before + B_after   or   A_after + B_before */
								if (a_bef + b_aft >= a_aft + b_bef)
								{
									/* Keep A before crossing, B after crossing */
									offset_pts[surv_idx[xai + 1]] = xp;
									for (m = xai + 2; m < run_bd[ra * 2 + 1]; ++m)
										surv_rm2[m] = 1;
									offset_pts[surv_idx[xbj]] = xp;
									for (m = run_bd[rb * 2]; m < xbj; ++m)
										surv_rm2[m] = 1;
								}
								else
								{
									/* Keep A after crossing, B before crossing */
									for (m = run_bd[ra * 2]; m < xai; ++m)
										surv_rm2[m] = 1;
									offset_pts[surv_idx[xai]] = xp;
									offset_pts[surv_idx[xbj + 1]] = xp;
									for (m = xbj + 2; m < run_bd[rb * 2 + 1]; ++m)
										surv_rm2[m] = 1;
								}
							}
						}
					}

					/* Compact */
					{
						unsigned int new_len = 0;
						for (i = 0; i < surv_len; ++i)
							if (!surv_rm2[i])
								surv_idx[new_len++] = surv_idx[i];
						surv_len = new_len;
					}

					free(surv_rm2);
				}
			}

			free(is_real_seg);

			/* --- Normal run-fitting --- */
			{
				unsigned int run_start = 0;

				while (run_start < surv_len)
				{
					unsigned int run_end = run_start + 1;
					qaws_scalar *flat_pts;
					unsigned int run_len, n_out, step_val;
					qaws_catmull_rom_desc cr_desc;
					qaws_curve *out_crv = NULL;

					while (run_end < surv_len &&
					       surv_idx[run_end] - surv_idx[run_end - 1] <= 3)
						run_end++;

					run_len = run_end - run_start;
					if (run_len < OFFSET_MIN_PIECE_PTS)
					{
						run_start = run_end;
						continue;
					}

					n_out = run_len < OFFSET_OUTPUT_PTS ? run_len : OFFSET_OUTPUT_PTS;
					if (n_out < 2) n_out = 2;

					flat_pts = (qaws_scalar *)malloc((n_out + 2) * 2 * sizeof(qaws_scalar));
					if (!flat_pts) { free(surv_idx); free(offset_pts); free(removed);
						return QAWS_STATUS_ALLOCATION_FAILURE; }

					flat_pts[0] = offset_pts[surv_idx[run_start]].x;
					flat_pts[1] = offset_pts[surv_idx[run_start]].y;

					step_val = run_len - 1;
					for (i = 0; i < n_out; ++i)
					{
						unsigned int src = (unsigned int)(
							(qaws_scalar)i * (qaws_scalar)step_val / (qaws_scalar)(n_out - 1)
							+ (qaws_scalar)0.5);
						if (src >= run_len) src = run_len - 1;
						flat_pts[(i + 1) * 2 + 0] = offset_pts[surv_idx[run_start + src]].x;
						flat_pts[(i + 1) * 2 + 1] = offset_pts[surv_idx[run_start + src]].y;
					}
					flat_pts[2] = offset_pts[surv_idx[run_start]].x;
					flat_pts[3] = offset_pts[surv_idx[run_start]].y;
					flat_pts[n_out * 2 + 0] = offset_pts[surv_idx[run_end - 1]].x;
					flat_pts[n_out * 2 + 1] = offset_pts[surv_idx[run_end - 1]].y;

					flat_pts[(n_out + 1) * 2 + 0] = offset_pts[surv_idx[run_end - 1]].x;
					flat_pts[(n_out + 1) * 2 + 1] = offset_pts[surv_idx[run_end - 1]].y;

					memset(&cr_desc, 0, sizeof(cr_desc));
					cr_desc.dimension = QAWS_DIMENSION_2D;
					cr_desc.control_points = flat_pts;
					cr_desc.control_point_count = n_out + 2;
					cr_desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
					cr_desc.closed = 0;

					status = qaws_curve_create_catmull_rom(&cr_desc, &out_crv);
					free(flat_pts);

					if (status == QAWS_STATUS_OK && out_crv)
					{
						if (curve_count < curve_capacity)
							out_curves[curve_count] = out_crv;
						else
							qaws_curve_destroy(out_crv);
						curve_count++;
					}

					run_start = run_end;
				}
			}
		}

		free(surv_idx);
		free(offset_pts);
		free(removed);
	}

	/* ---- trim==2 post-process: fix inter-curve crossings in fitted output ---- */
	if (trim == 2 && curve_count > 1)
	{
		unsigned int ca, cb;
		unsigned int eval_n = 128;
		int any_fix = 1;
		int iter_limit = 10;

		while (any_fix && iter_limit-- > 0)
		{
			any_fix = 0;
			for (ca = 0; ca + 1 < curve_count && ca < curve_capacity && !any_fix; ++ca)
			{
				for (cb = ca + 1; cb < curve_count && cb < curve_capacity; ++cb)
				{
					qaws_vec2 *pts_a, *pts_b;
					int xfound = 0;
					unsigned int xai = 0, xbj = 0;
					qaws_scalar xft = 0, xfu = 0;

					pts_a = (qaws_vec2 *)malloc(eval_n * sizeof(qaws_vec2));
					pts_b = (qaws_vec2 *)malloc(eval_n * sizeof(qaws_vec2));
					if (!pts_a || !pts_b) { free(pts_a); free(pts_b); continue; }

					{
						qaws_range ra2 = qaws_curve_get_parameter_range(out_curves[ca]);
						qaws_range rb2 = qaws_curve_get_parameter_range(out_curves[cb]);
						for (i = 0; i < eval_n; ++i)
						{
							qaws_scalar fa = (qaws_scalar)i / (qaws_scalar)(eval_n - 1);
							qaws_scalar fb = fa;
							qaws_curve_evaluate_2d(out_curves[ca],
								ra2.min_value + fa * (ra2.max_value - ra2.min_value),
								QAWS_EVAL_FLAG_POSITION, &result);
							pts_a[i] = result.position;
							qaws_curve_evaluate_2d(out_curves[cb],
								rb2.min_value + fb * (rb2.max_value - rb2.min_value),
								QAWS_EVAL_FLAG_POSITION, &result);
							pts_b[i] = result.position;
						}
					}

					for (i = 0; i + 1 < eval_n && !xfound; ++i)
						for (j = 0; j + 1 < eval_n; ++j)
						{
							qaws_scalar ft, fu;
							qaws_scalar ta_cand, tb_cand;
							if (offset_seg_isect(pts_a[i], pts_a[i + 1],
								pts_b[j], pts_b[j + 1], &ft, &fu))
							{
								/* Compute normalised parameters */
								ta_cand = ((qaws_scalar)i + ft) / (qaws_scalar)(eval_n - 1);
								tb_cand = ((qaws_scalar)j + fu) / (qaws_scalar)(eval_n - 1);
								/* Skip shared-endpoint crossings at curve extremes */
								if (ta_cand < (qaws_scalar)0.02 && tb_cand > (qaws_scalar)0.98)
									continue;
								if (ta_cand > (qaws_scalar)0.98 && tb_cand < (qaws_scalar)0.02)
									continue;
								xft = ft; xfu = fu;
								xai = i; xbj = j;
								xfound = 1; break;
							}
						}

					if (xfound)
					{
						unsigned int a_bef = xai, a_aft = eval_n - 2 - xai;
						unsigned int b_bef = xbj, b_aft = eval_n - 2 - xbj;
						qaws_scalar ta = ((qaws_scalar)xai + xft) / (qaws_scalar)(eval_n - 1);
						qaws_scalar tb = ((qaws_scalar)xbj + xfu) / (qaws_scalar)(eval_n - 1);

						if (a_bef + b_aft >= a_aft + b_bef)
						{
							offset_trunc_curve(out_curves, ca, (qaws_scalar)0, ta);
							offset_trunc_curve(out_curves, cb, tb, (qaws_scalar)1);
						}
						else
						{
							offset_trunc_curve(out_curves, ca, ta, (qaws_scalar)1);
							offset_trunc_curve(out_curves, cb, (qaws_scalar)0, tb);
						}
						any_fix = 1;
					}

					free(pts_a);
					free(pts_b);
					if (any_fix) break;
				}
			}
		}
	}

	*out_count = curve_count < curve_capacity ? curve_count : curve_capacity;
	return curve_count > 0 ? QAWS_STATUS_OK : QAWS_STATUS_INTERNAL_ERROR;
}

/* ========================================================================== */
/*  Arc-length reparameterization                                              */
/* ========================================================================== */

typedef struct qaws_arc_reparam_impl
{
	qaws_curve const* source;
	qaws_scalar* table_params;
	qaws_scalar* table_distances;
	unsigned int table_size;
	qaws_scalar total_arc_length;
} qaws_arc_reparam_impl;

/* ------------------------------------------------------------------ */
/*  Evaluation helpers                                                  */
/* ------------------------------------------------------------------ */

static qaws_status reparam_eval_span_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_arc_reparam_impl const* impl =
		(qaws_arc_reparam_impl const*)curve->impl;
	qaws_scalar s, src_t;
	qaws_eval_result_2d src;
	qaws_scalar speed, speed2, speed4;
	qaws_scalar L;
	unsigned int src_flags;
	qaws_status status;

	(void)span_index;

	L = impl->total_arc_length;
	s = local_t * L;
	src_t = qaws_internal_distance_to_parameter(
		impl->table_params, impl->table_distances,
		impl->table_size, s);

	/* We always need at least D1 from the source to compute speed for
	   the chain rule.  If D2 is requested we also need source D2. */
	src_flags = eval_flags | QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		src_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		src_flags |= QAWS_EVAL_FLAG_D2;
	if (eval_flags & QAWS_EVAL_FLAG_D3)
		src_flags |= QAWS_EVAL_FLAG_D3;

	status = qaws_curve_evaluate_2d(impl->source, src_t, src_flags, &src);
	if (status != QAWS_STATUS_OK)
		return status;

	memset(out_result, 0, sizeof(*out_result));
	out_result->valid_flags = 0;

	/* Position: pass through */
	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position = src.position;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* Compute speed = |D1| from source */
	speed = (qaws_scalar)sqrt((double)(
		src.d1.x * src.d1.x + src.d1.y * src.d1.y));
	if (speed < (qaws_scalar)1e-15)
		speed = (qaws_scalar)1e-15;

	/* D1: dC/d(local_t) = D1_src * L / speed */
	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		qaws_scalar scale = L / speed;
		out_result->d1.x = src.d1.x * scale;
		out_result->d1.y = src.d1.y * scale;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2: d2C/d(local_t)^2 = d2C/ds^2 * L^2
	   where d2C/ds^2 = (D2 * speed^2 - D1 * dot(D1,D2)) / speed^4 */
	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		qaws_scalar d1_dot_d2 = src.d1.x * src.d2.x + src.d1.y * src.d2.y;
		speed2 = speed * speed;
		speed4 = speed2 * speed2;
		qaws_scalar L2 = L * L;
		out_result->d2.x = (src.d2.x * speed2 - src.d1.x * d1_dot_d2) / speed4 * L2;
		out_result->d2.y = (src.d2.y * speed2 - src.d1.y * d1_dot_d2) / speed4 * L2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3: pass through with scaling d3C/ds^3 * L^3
	   Full formula is complex; use a simplified chain rule approximation.
	   For many practical uses (rendering, sampling), D3 is rarely needed.
	   We compute it exactly when source D3 is available:
	   d3C/ds^3 = (1/speed^3) * [D3/speed^3 - 3*(D1.D2)*D2/speed^5
	              - (|D2|^2 + D1.D3)*D1/speed^5 + 5*(D1.D2)^2*D1/speed^7]
	   Then multiply by L^3. */
	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		qaws_scalar d1_dot_d2 = src.d1.x * src.d2.x + src.d1.y * src.d2.y;
		qaws_scalar d1_dot_d3 = src.d1.x * src.d3.x + src.d1.y * src.d3.y;
		qaws_scalar d2_sq = src.d2.x * src.d2.x + src.d2.y * src.d2.y;
		speed2 = speed * speed;
		qaws_scalar speed3 = speed2 * speed;
		speed4 = speed2 * speed2;
		qaws_scalar speed5 = speed4 * speed;
		qaws_scalar speed7 = speed5 * speed2;
		qaws_scalar L3 = L * L * L;

		/* d3C/ds^3 components */
		qaws_scalar c3x = src.d3.x / speed3
			- (qaws_scalar)3 * d1_dot_d2 * src.d2.x / speed5
			- (d2_sq + d1_dot_d3) * src.d1.x / speed5
			+ (qaws_scalar)5 * d1_dot_d2 * d1_dot_d2 * src.d1.x / speed7;
		qaws_scalar c3y = src.d3.y / speed3
			- (qaws_scalar)3 * d1_dot_d2 * src.d2.y / speed5
			- (d2_sq + d1_dot_d3) * src.d1.y / speed5
			+ (qaws_scalar)5 * d1_dot_d2 * d1_dot_d2 * src.d1.y / speed7;

		out_result->d3.x = c3x * L3;
		out_result->d3.y = c3y * L3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static qaws_status reparam_eval_span_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_arc_reparam_impl const* impl =
		(qaws_arc_reparam_impl const*)curve->impl;
	qaws_scalar s, src_t;
	qaws_eval_result_3d src;
	qaws_scalar speed, speed2, speed4;
	qaws_scalar L;
	unsigned int src_flags;
	qaws_status status;

	(void)span_index;

	L = impl->total_arc_length;
	s = local_t * L;
	src_t = qaws_internal_distance_to_parameter(
		impl->table_params, impl->table_distances,
		impl->table_size, s);

	src_flags = eval_flags | QAWS_EVAL_FLAG_POSITION;
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		src_flags |= QAWS_EVAL_FLAG_D1;
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
		src_flags |= QAWS_EVAL_FLAG_D2;
	if (eval_flags & QAWS_EVAL_FLAG_D3)
		src_flags |= QAWS_EVAL_FLAG_D3;

	status = qaws_curve_evaluate_3d(impl->source, src_t, src_flags, &src);
	if (status != QAWS_STATUS_OK)
		return status;

	memset(out_result, 0, sizeof(*out_result));
	out_result->valid_flags = 0;

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position = src.position;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	speed = (qaws_scalar)sqrt((double)(
		src.d1.x * src.d1.x + src.d1.y * src.d1.y + src.d1.z * src.d1.z));
	if (speed < (qaws_scalar)1e-15)
		speed = (qaws_scalar)1e-15;

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		qaws_scalar scale = L / speed;
		out_result->d1.x = src.d1.x * scale;
		out_result->d1.y = src.d1.y * scale;
		out_result->d1.z = src.d1.z * scale;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		qaws_scalar d1_dot_d2 = src.d1.x * src.d2.x + src.d1.y * src.d2.y
			+ src.d1.z * src.d2.z;
		speed2 = speed * speed;
		speed4 = speed2 * speed2;
		qaws_scalar L2 = L * L;
		out_result->d2.x = (src.d2.x * speed2 - src.d1.x * d1_dot_d2) / speed4 * L2;
		out_result->d2.y = (src.d2.y * speed2 - src.d1.y * d1_dot_d2) / speed4 * L2;
		out_result->d2.z = (src.d2.z * speed2 - src.d1.z * d1_dot_d2) / speed4 * L2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		qaws_scalar d1_dot_d2 = src.d1.x * src.d2.x + src.d1.y * src.d2.y
			+ src.d1.z * src.d2.z;
		qaws_scalar d1_dot_d3 = src.d1.x * src.d3.x + src.d1.y * src.d3.y
			+ src.d1.z * src.d3.z;
		qaws_scalar d2_sq = src.d2.x * src.d2.x + src.d2.y * src.d2.y
			+ src.d2.z * src.d2.z;
		speed2 = speed * speed;
		qaws_scalar speed3 = speed2 * speed;
		speed4 = speed2 * speed2;
		qaws_scalar speed5 = speed4 * speed;
		qaws_scalar speed7 = speed5 * speed2;
		qaws_scalar L3 = L * L * L;

		qaws_scalar c3x = src.d3.x / speed3
			- (qaws_scalar)3 * d1_dot_d2 * src.d2.x / speed5
			- (d2_sq + d1_dot_d3) * src.d1.x / speed5
			+ (qaws_scalar)5 * d1_dot_d2 * d1_dot_d2 * src.d1.x / speed7;
		qaws_scalar c3y = src.d3.y / speed3
			- (qaws_scalar)3 * d1_dot_d2 * src.d2.y / speed5
			- (d2_sq + d1_dot_d3) * src.d1.y / speed5
			+ (qaws_scalar)5 * d1_dot_d2 * d1_dot_d2 * src.d1.y / speed7;
		qaws_scalar c3z = src.d3.z / speed3
			- (qaws_scalar)3 * d1_dot_d2 * src.d2.z / speed5
			- (d2_sq + d1_dot_d3) * src.d1.z / speed5
			+ (qaws_scalar)5 * d1_dot_d2 * d1_dot_d2 * src.d1.z / speed7;

		out_result->d3.x = c3x * L3;
		out_result->d3.y = c3y * L3;
		out_result->d3.z = c3z * L3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Destroy                                                            */
/* ------------------------------------------------------------------ */

static void reparam_destroy_impl(void* impl_ptr, qaws_allocator const* allocator)
{
	qaws_arc_reparam_impl* impl = (qaws_arc_reparam_impl*)impl_ptr;
	if (!impl) return;
	qaws_internal_dealloc(allocator, impl->table_params);
	qaws_internal_dealloc(allocator, impl->table_distances);
	qaws_internal_dealloc(allocator, impl);
}

/* ------------------------------------------------------------------ */
/*  Property queries                                                   */
/* ------------------------------------------------------------------ */

static int reparam_is_closed(qaws_curve const* curve)
{
	qaws_arc_reparam_impl const* impl =
		(qaws_arc_reparam_impl const*)curve->impl;
	if (impl->source->vtable && impl->source->vtable->is_closed)
		return impl->source->vtable->is_closed(impl->source);
	return 0;
}

static int reparam_is_periodic(qaws_curve const* curve)
{
	qaws_arc_reparam_impl const* impl =
		(qaws_arc_reparam_impl const*)curve->impl;
	if (impl->source->vtable && impl->source->vtable->is_periodic)
		return impl->source->vtable->is_periodic(impl->source);
	return 0;
}

static int reparam_is_rational(qaws_curve const* curve)
{
	qaws_arc_reparam_impl const* impl =
		(qaws_arc_reparam_impl const*)curve->impl;
	if (impl->source->vtable && impl->source->vtable->is_rational)
		return impl->source->vtable->is_rational(impl->source);
	return 0;
}

static qaws_continuity reparam_get_continuity(qaws_curve const* curve)
{
	qaws_arc_reparam_impl const* impl =
		(qaws_arc_reparam_impl const*)curve->impl;
	if (impl->source->vtable && impl->source->vtable->get_continuity)
		return impl->source->vtable->get_continuity(impl->source);
	return QAWS_CONTINUITY_NONE;
}

/* ------------------------------------------------------------------ */
/*  Vtable                                                             */
/* ------------------------------------------------------------------ */

static qaws_curve_vtable const reparam_vtable = {
	reparam_eval_span_2d,
	reparam_eval_span_3d,
	reparam_destroy_impl,
	reparam_is_closed,
	reparam_is_periodic,
	reparam_is_rational,
	reparam_get_continuity,
	NULL /* diff */
};

/* ------------------------------------------------------------------ */
/*  Public API                                                         */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_reparameterize_arc_length(
	qaws_curve const* curve,
	unsigned int table_resolution,
	qaws_curve** out_curve)
{
	qaws_arc_reparam_impl* impl = NULL;
	qaws_curve* result = NULL;
	qaws_range range;
	qaws_status status;
	unsigned int tbl_size;

	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	tbl_size = table_resolution > 0 ? table_resolution : 256;

	/* Allocate impl */
	impl = (qaws_arc_reparam_impl*)malloc(sizeof(qaws_arc_reparam_impl));
	if (!impl)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	memset(impl, 0, sizeof(*impl));
	impl->source = curve;
	impl->table_size = tbl_size;

	impl->table_params = (qaws_scalar*)malloc(
		sizeof(qaws_scalar) * (size_t)tbl_size);
	if (!impl->table_params) {
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->table_distances = (qaws_scalar*)malloc(
		sizeof(qaws_scalar) * (size_t)tbl_size);
	if (!impl->table_distances) {
		free(impl->table_params);
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Build arc-length table from the source curve */
	status = qaws_internal_build_arc_length_table(
		curve, tbl_size, impl->table_params, impl->table_distances);
	if (status != QAWS_STATUS_OK) {
		free(impl->table_distances);
		free(impl->table_params);
		free(impl);
		return status;
	}

	impl->total_arc_length = impl->table_distances[tbl_size - 1];

	/* Degenerate curve check */
	if (impl->total_arc_length <= (qaws_scalar)0) {
		free(impl->table_distances);
		free(impl->table_params);
		free(impl);
		return QAWS_STATUS_DEGENERATE_CURVE;
	}

	/* Allocate the wrapper curve: single span [0, total_arc_length] */
	range.min_value = (qaws_scalar)0;
	range.max_value = impl->total_arc_length;

	result = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_REPARAMETERIZED,
		curve->dimension,
		curve->degree,
		1,    /* span_count */
		range,
		&reparam_vtable);
	if (!result) {
		free(impl->table_distances);
		free(impl->table_params);
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	result->span_boundaries[0] = (qaws_scalar)0;
	result->span_boundaries[1] = impl->total_arc_length;
	result->impl = impl;

	*out_curve = result;
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Curve length matching (#23)                                               */
/* ========================================================================== */

/* A length-match wrapper reuses the reparameterization machinery:
   it builds arc-length tables for both curves, then maps from
   curve_a's arc-length domain to curve_b's parameter domain. */

typedef struct qaws_length_match_impl
{
	qaws_curve const* source_b;
	qaws_scalar* table_params_a;
	qaws_scalar* table_distances_a;
	unsigned int table_size_a;
	qaws_scalar total_arc_length_a;
	qaws_scalar* table_params_b;
	qaws_scalar* table_distances_b;
	unsigned int table_size_b;
	qaws_scalar total_arc_length_b;
} qaws_length_match_impl;

static qaws_status lm_eval_span_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_length_match_impl const* impl =
		(qaws_length_match_impl const*)curve->impl;
	qaws_scalar s, frac, dist_b, src_t;

	(void)span_index;

	/* s is in [0, total_arc_length_a] */
	s = curve->parameter_range.min_value +
		local_t * (curve->parameter_range.max_value - curve->parameter_range.min_value);

	/* Map s to fraction of A's length, then to same fraction of B's length */
	frac = (impl->total_arc_length_a > QAWS_ZERO) ?
		(s / impl->total_arc_length_a) : QAWS_ZERO;
	dist_b = frac * impl->total_arc_length_b;

	/* Map distance in B to source parameter */
	src_t = qaws_internal_distance_to_parameter(
		impl->table_params_b, impl->table_distances_b,
		impl->table_size_b, dist_b);

	return qaws_curve_evaluate_2d(impl->source_b, src_t, eval_flags, out_result);
}

static qaws_status lm_eval_span_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_length_match_impl const* impl =
		(qaws_length_match_impl const*)curve->impl;
	qaws_scalar s, frac, dist_b, src_t;

	(void)span_index;

	s = curve->parameter_range.min_value +
		local_t * (curve->parameter_range.max_value - curve->parameter_range.min_value);

	frac = (impl->total_arc_length_a > QAWS_ZERO) ?
		(s / impl->total_arc_length_a) : QAWS_ZERO;
	dist_b = frac * impl->total_arc_length_b;

	src_t = qaws_internal_distance_to_parameter(
		impl->table_params_b, impl->table_distances_b,
		impl->table_size_b, dist_b);

	return qaws_curve_evaluate_3d(impl->source_b, src_t, eval_flags, out_result);
}

static void lm_destroy_impl(void* impl_ptr, qaws_allocator const* allocator)
{
	qaws_length_match_impl* impl = (qaws_length_match_impl*)impl_ptr;
	if (!impl) return;
	qaws_internal_dealloc(allocator, impl->table_params_a);
	qaws_internal_dealloc(allocator, impl->table_distances_a);
	qaws_internal_dealloc(allocator, impl->table_params_b);
	qaws_internal_dealloc(allocator, impl->table_distances_b);
	qaws_internal_dealloc(allocator, impl);
}

static int lm_is_closed(qaws_curve const* curve)
{
	qaws_length_match_impl const* impl =
		(qaws_length_match_impl const*)curve->impl;
	if (impl->source_b->vtable && impl->source_b->vtable->is_closed)
		return impl->source_b->vtable->is_closed(impl->source_b);
	return 0;
}

static int lm_not(qaws_curve const* curve) { (void)curve; return 0; }
static qaws_continuity lm_cont(qaws_curve const* curve) { (void)curve; return QAWS_CONTINUITY_C2; }

static qaws_curve_vtable const lm_vtable = {
	lm_eval_span_2d,
	lm_eval_span_3d,
	lm_destroy_impl,
	lm_is_closed,
	lm_not, /* is_periodic */
	lm_not, /* is_rational */
	lm_cont,
	NULL /* diff */
};

int qaws_internal_reparam_source(
	qaws_curve const* curve,
	qaws_scalar t,
	qaws_curve const** out_source,
	qaws_scalar* out_source_t)
{
	qaws_scalar local_t;
	if (curve->vtable == &reparam_vtable)
	{
		qaws_arc_reparam_impl const* impl = (qaws_arc_reparam_impl const*)curve->impl;
		qaws_internal_find_span(curve, t, &local_t);
		*out_source = impl->source;
		*out_source_t = qaws_internal_distance_to_parameter(
			impl->table_params, impl->table_distances, impl->table_size,
			local_t * impl->total_arc_length);
		return 1;
	}
	if (curve->vtable == &lm_vtable)
	{
		qaws_length_match_impl const* impl = (qaws_length_match_impl const*)curve->impl;
		qaws_scalar s, frac;
		qaws_internal_find_span(curve, t, &local_t);
		s = curve->parameter_range.min_value +
			local_t * (curve->parameter_range.max_value - curve->parameter_range.min_value);
		frac = (impl->total_arc_length_a > QAWS_ZERO) ? (s / impl->total_arc_length_a) : QAWS_ZERO;
		*out_source = impl->source_b;
		*out_source_t = qaws_internal_distance_to_parameter(
			impl->table_params_b, impl->table_distances_b, impl->table_size_b,
			frac * impl->total_arc_length_b);
		return 1;
	}
	return 0;
}

qaws_status qaws_curve_match_arc_length(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	unsigned int table_resolution,
	qaws_curve** out_curve)
{
	qaws_length_match_impl* impl = NULL;
	qaws_curve* result = NULL;
	qaws_range range;
	qaws_status status;
	unsigned int tbl_size;

	if (!curve_a || !curve_b || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;
	tbl_size = table_resolution > 0 ? table_resolution : 256;

	impl = (qaws_length_match_impl*)malloc(sizeof(qaws_length_match_impl));
	if (!impl) return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(impl, 0, sizeof(*impl));

	impl->source_b = curve_b;
	impl->table_size_a = tbl_size;
	impl->table_size_b = tbl_size;

	/* Allocate tables */
	impl->table_params_a = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)tbl_size);
	impl->table_distances_a = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)tbl_size);
	impl->table_params_b = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)tbl_size);
	impl->table_distances_b = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)tbl_size);

	if (!impl->table_params_a || !impl->table_distances_a ||
		!impl->table_params_b || !impl->table_distances_b)
	{
		free(impl->table_params_a); free(impl->table_distances_a);
		free(impl->table_params_b); free(impl->table_distances_b);
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Build arc-length tables */
	status = qaws_internal_build_arc_length_table(
		curve_a, tbl_size, impl->table_params_a, impl->table_distances_a);
	if (status != QAWS_STATUS_OK) { lm_destroy_impl(impl, NULL); return status; }

	status = qaws_internal_build_arc_length_table(
		curve_b, tbl_size, impl->table_params_b, impl->table_distances_b);
	if (status != QAWS_STATUS_OK) { lm_destroy_impl(impl, NULL); return status; }

	impl->total_arc_length_a = impl->table_distances_a[tbl_size - 1];
	impl->total_arc_length_b = impl->table_distances_b[tbl_size - 1];

	if (impl->total_arc_length_a <= QAWS_ZERO || impl->total_arc_length_b <= QAWS_ZERO)
	{
		lm_destroy_impl(impl, NULL);
		return QAWS_STATUS_DEGENERATE_CURVE;
	}

	/* Create wrapper curve with parameter domain [0, arc_length_a] */
	range.min_value = QAWS_ZERO;
	range.max_value = impl->total_arc_length_a;

	result = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_REPARAMETERIZED,
		curve_b->dimension,
		curve_b->degree,
		1, range, &lm_vtable);
	if (!result) { lm_destroy_impl(impl, NULL); return QAWS_STATUS_ALLOCATION_FAILURE; }

	result->span_boundaries[0] = QAWS_ZERO;
	result->span_boundaries[1] = impl->total_arc_length_a;
	result->impl = impl;

	*out_curve = result;
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  3D curve offsetting (#21)                                                 */
/* ========================================================================== */

qaws_status qaws_curve_offset_3d(
	qaws_curve const* curve,
	qaws_scalar distance,
	int direction_mode,
	qaws_vec3 const* direction,
	qaws_surface const* surface,
	unsigned int sample_count,
	qaws_curve** out_curve)
{
	unsigned int n_samples;
	qaws_scalar* params = NULL;
	qaws_scalar* coords = NULL;
	qaws_range range;
	unsigned int i;
	unsigned int dim = 3;
	unsigned int n_cp;
	qaws_status status;

	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(curve) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (direction_mode == 0 && !direction)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (direction_mode == 2 && !surface)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;
	n_samples = sample_count > 0 ? sample_count : 256;
	range = qaws_curve_get_parameter_range(curve);

	params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samples);
	coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samples * dim);
	if (!params || !coords) { free(params); free(coords); return QAWS_STATUS_ALLOCATION_FAILURE; }

	for (i = 0; i < n_samples; i++)
	{
		qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
			(qaws_scalar)i / (qaws_scalar)(n_samples - 1);
		qaws_eval_result_3d er;
		qaws_vec3 offset_dir;
		qaws_scalar dir_len;

		memset(&er, 0, sizeof(er));
		status = qaws_curve_evaluate_3d(curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &er);
		if (status != QAWS_STATUS_OK) { free(params); free(coords); return status; }

		params[i] = t;

		if (direction_mode == 0)
		{
			/* Constant direction */
			offset_dir = *direction;
		}
		else if (direction_mode == 1)
		{
			/* Frenet normal */
			qaws_vec3 T_vec, N_vec, B_vec;
			status = qaws_curve_compute_frenet_frame_3d(curve, t, &T_vec, &N_vec, &B_vec);
			if (status != QAWS_STATUS_OK)
			{
				/* Fallback: use a perpendicular to the tangent */
				qaws_scalar tx = er.d1.x, ty = er.d1.y, tz = er.d1.z;
				if (QAWS_FABS(tx) < QAWS_FABS(tz))
				{ offset_dir.x = QAWS_ZERO; offset_dir.y = -tz; offset_dir.z = ty; }
				else
				{ offset_dir.x = -ty; offset_dir.y = tx; offset_dir.z = QAWS_ZERO; }
			}
			else
			{
				offset_dir = N_vec;
			}
		}
		else /* direction_mode == 2 */
		{
			/* Surface normal at closest point */
			qaws_scalar su, sv;
			qaws_vec3 closest;
			qaws_surface_eval_result sr;

			status = qaws_surface_find_closest_point(surface, er.position, &su, &sv, &closest);
			if (status != QAWS_STATUS_OK) { free(params); free(coords); return status; }

			memset(&sr, 0, sizeof(sr));
			qaws_surface_evaluate(surface, su, sv, QAWS_SURFACE_EVAL_NORMAL, &sr);
			offset_dir = sr.normal;
		}

		/* Normalize direction */
		dir_len = QAWS_SQRT(offset_dir.x * offset_dir.x +
			offset_dir.y * offset_dir.y + offset_dir.z * offset_dir.z);
		if (dir_len > QAWS_LITERAL(1e-12))
		{
			offset_dir.x /= dir_len;
			offset_dir.y /= dir_len;
			offset_dir.z /= dir_len;
		}

		coords[i * dim + 0] = er.position.x + distance * offset_dir.x;
		coords[i * dim + 1] = er.position.y + distance * offset_dir.y;
		coords[i * dim + 2] = er.position.z + distance * offset_dir.z;
	}

	/* Fit a B-spline through the offset samples */
	n_cp = n_samples / 4;
	if (n_cp < 8) n_cp = 8;
	if (n_cp > n_samples) n_cp = n_samples;

	status = qaws_internal_fit_bspline(
		QAWS_DIMENSION_3D, 3,
		params, coords,
		n_samples, n_cp,
		out_curve);

	free(params);
	free(coords);
	return status;
}

/* ========================================================================== */
/*  Approximate curve merging (#22)                                           */
/* ========================================================================== */

qaws_status qaws_curve_merge_chain(
	qaws_curve const* const* curves,
	unsigned int curve_count,
	unsigned int target_degree,
	qaws_scalar tolerance,
	qaws_curve** out_curves,
	unsigned int curve_capacity,
	unsigned int* out_count)
{
	unsigned int degree;
	qaws_dimension dim;
	unsigned int total_samples = 0;
	unsigned int samples_per_curve;
	qaws_scalar* params = NULL;
	qaws_scalar* coords = NULL;
	unsigned int dim_count;
	unsigned int ci, si, idx;
	qaws_scalar param_offset;
	unsigned int n_cp;
	qaws_status status;
	qaws_curve* merged = NULL;
	qaws_scalar max_err;

	if (!curves || curve_count == 0 || !out_curves || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve_capacity == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_count = 0;
	degree = target_degree > 0 ? target_degree : 3;
	dim = qaws_curve_get_dimension(curves[0]);
	dim_count = (unsigned int)dim;

	/* Sample all curves densely */
	samples_per_curve = 64;
	total_samples = samples_per_curve * curve_count;

	params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)total_samples);
	coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)total_samples * dim_count);
	if (!params || !coords) { free(params); free(coords); return QAWS_STATUS_ALLOCATION_FAILURE; }

	idx = 0;
	param_offset = QAWS_ZERO;
	for (ci = 0; ci < curve_count; ci++)
	{
		qaws_range range = qaws_curve_get_parameter_range(curves[ci]);
		qaws_scalar range_len = range.max_value - range.min_value;

		for (si = 0; si < samples_per_curve; si++)
		{
			qaws_scalar frac = (qaws_scalar)si / (qaws_scalar)(samples_per_curve - 1);
			qaws_scalar t = range.min_value + frac * range_len;

			params[idx] = param_offset + frac * range_len;

			if (dim == QAWS_DIMENSION_3D)
			{
				qaws_eval_result_3d er;
				memset(&er, 0, sizeof(er));
				qaws_curve_evaluate_3d(curves[ci], t, QAWS_EVAL_FLAG_POSITION, &er);
				coords[idx * dim_count + 0] = er.position.x;
				coords[idx * dim_count + 1] = er.position.y;
				coords[idx * dim_count + 2] = er.position.z;
			}
			else
			{
				qaws_eval_result_2d er;
				memset(&er, 0, sizeof(er));
				qaws_curve_evaluate_2d(curves[ci], t, QAWS_EVAL_FLAG_POSITION, &er);
				coords[idx * dim_count + 0] = er.position.x;
				coords[idx * dim_count + 1] = er.position.y;
			}
			idx++;
		}
		param_offset += range_len;
	}

	/* Fit a single B-spline through all samples */
	n_cp = total_samples / 4;
	if (n_cp < degree + 1) n_cp = degree + 1;
	if (n_cp > total_samples) n_cp = total_samples;

	status = qaws_internal_fit_bspline(dim, degree, params, coords,
		total_samples, n_cp, &merged);

	if (status != QAWS_STATUS_OK) { free(params); free(coords); return status; }

	/* Check maximum deviation */
	max_err = QAWS_ZERO;
	for (idx = 0; idx < total_samples; idx++)
	{
		qaws_scalar err;
		if (dim == QAWS_DIMENSION_3D)
		{
			qaws_eval_result_3d er;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(merged, params[idx], QAWS_EVAL_FLAG_POSITION, &er);
			err = QAWS_SQRT(
				(er.position.x - coords[idx * dim_count + 0]) * (er.position.x - coords[idx * dim_count + 0]) +
				(er.position.y - coords[idx * dim_count + 1]) * (er.position.y - coords[idx * dim_count + 1]) +
				(er.position.z - coords[idx * dim_count + 2]) * (er.position.z - coords[idx * dim_count + 2]));
		}
		else
		{
			qaws_eval_result_2d er;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(merged, params[idx], QAWS_EVAL_FLAG_POSITION, &er);
			err = QAWS_SQRT(
				(er.position.x - coords[idx * dim_count + 0]) * (er.position.x - coords[idx * dim_count + 0]) +
				(er.position.y - coords[idx * dim_count + 1]) * (er.position.y - coords[idx * dim_count + 1]));
		}
		if (err > max_err) max_err = err;
	}

	/* If tolerance exceeded and we have room, try with more control points */
	if (tolerance > QAWS_ZERO && max_err > tolerance && n_cp < total_samples)
	{
		qaws_curve_destroy(merged);
		merged = NULL;
		n_cp = total_samples / 2;
		if (n_cp < degree + 1) n_cp = degree + 1;
		status = qaws_internal_fit_bspline(dim, degree, params, coords,
			total_samples, n_cp, &merged);
		if (status != QAWS_STATUS_OK)
		{
			free(params); free(coords);
			return status;
		}
	}

	free(params);
	free(coords);

	out_curves[0] = merged;
	*out_count = 1;
	return QAWS_STATUS_OK;
}

/* ========================================================================== */
/*  Fillet and chamfer on 2D curves (#25)                                     */
/* ========================================================================== */

/* Helper: compute the angle between two 2D vectors */
static qaws_scalar vec2_angle(qaws_vec2 a, qaws_vec2 b)
{
	qaws_scalar dot = a.x * b.x + a.y * b.y;
	qaws_scalar cross = a.x * b.y - a.y * b.x;
	return QAWS_ATAN2(cross, dot);
}

qaws_status qaws_curve_fillet_2d(
	qaws_curve const* curve,
	qaws_scalar radius,
	qaws_curve** out_curve)
{
	qaws_curve_kind kind;
	unsigned int seg_count, total_segs;
	unsigned int seg_i, ji;
	qaws_curve** new_segments = NULL;
	unsigned int new_count = 0;
	qaws_composite_desc cdesc;
	qaws_status status = QAWS_STATUS_OK;

	/* Per-junction precomputed data */
	qaws_scalar* j_setback = NULL;
	qaws_scalar* j_angle = NULL;
	qaws_vec2* j_corner = NULL;
	qaws_vec2* j_tan_out = NULL;
	qaws_vec2* j_tan_in = NULL;
	qaws_scalar* j_d_out = NULL;
	qaws_scalar* j_d_in = NULL;
	int* j_active = NULL;
	unsigned int n_junctions;

	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (radius <= QAWS_ZERO)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(curve) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_curve = NULL;
	kind = qaws_curve_get_kind(curve);

	/* For non-composite curves, just clone */
	if (kind != QAWS_CURVE_KIND_COMPOSITE)
	{
		qaws_range range = qaws_curve_get_parameter_range(curve);
		unsigned int n_samp = 64;
		qaws_scalar* params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp);
		qaws_scalar* coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp * 2);
		unsigned int i;

		if (!params || !coords) { free(params); free(coords); return QAWS_STATUS_ALLOCATION_FAILURE; }

		for (i = 0; i < n_samp; i++)
		{
			qaws_eval_result_2d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)(n_samp - 1);
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &er);
			params[i] = t;
			coords[i * 2] = er.position.x;
			coords[i * 2 + 1] = er.position.y;
		}
		status = qaws_internal_fit_bspline(QAWS_DIMENSION_2D, 3, params, coords, n_samp, 16, out_curve);
		free(params); free(coords);
		return status;
	}

	/* Composite curve: two-pass approach */
	seg_count = qaws_curve_get_span_count(curve);
	n_junctions = (seg_count > 1) ? seg_count - 1 : 0;

	/* ---- Pass 1: precompute junction geometry ---- */
	if (n_junctions > 0)
	{
		j_setback = (qaws_scalar*)calloc((size_t)n_junctions, sizeof(qaws_scalar));
		j_angle   = (qaws_scalar*)calloc((size_t)n_junctions, sizeof(qaws_scalar));
		j_corner  = (qaws_vec2*)calloc((size_t)n_junctions, sizeof(qaws_vec2));
		j_tan_out = (qaws_vec2*)calloc((size_t)n_junctions, sizeof(qaws_vec2));
		j_tan_in  = (qaws_vec2*)calloc((size_t)n_junctions, sizeof(qaws_vec2));
		j_d_out   = (qaws_scalar*)calloc((size_t)n_junctions, sizeof(qaws_scalar));
		j_d_in    = (qaws_scalar*)calloc((size_t)n_junctions, sizeof(qaws_scalar));
		j_active  = (int*)calloc((size_t)n_junctions, sizeof(int));
		if (!j_setback || !j_angle || !j_corner || !j_tan_out || !j_tan_in
			|| !j_d_out || !j_d_in || !j_active)
		{
			free(j_setback); free(j_angle); free(j_corner);
			free(j_tan_out); free(j_tan_in); free(j_d_out); free(j_d_in); free(j_active);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		for (ji = 0; ji < n_junctions; ji++)
		{
			qaws_eval_result_2d er_end, er_start, er_corner;
			qaws_scalar t_junc = (qaws_scalar)(ji + 1);
			qaws_scalar half_angle;

			memset(&er_end, 0, sizeof(er_end));
			memset(&er_start, 0, sizeof(er_start));
			memset(&er_corner, 0, sizeof(er_corner));
			qaws_curve_evaluate_2d(curve, t_junc - QAWS_LITERAL(1e-6), QAWS_EVAL_FLAG_D1, &er_end);
			qaws_curve_evaluate_2d(curve, t_junc + QAWS_LITERAL(1e-6), QAWS_EVAL_FLAG_D1, &er_start);
			qaws_curve_evaluate_2d(curve, t_junc, QAWS_EVAL_FLAG_POSITION, &er_corner);

			j_tan_out[ji] = er_end.d1;
			j_tan_in[ji] = er_start.d1;
			j_corner[ji] = er_corner.position;
			j_angle[ji] = vec2_angle(er_end.d1, er_start.d1);

			j_d_out[ji] = QAWS_SQRT(j_tan_out[ji].x * j_tan_out[ji].x + j_tan_out[ji].y * j_tan_out[ji].y);
			j_d_in[ji] = QAWS_SQRT(j_tan_in[ji].x * j_tan_in[ji].x + j_tan_in[ji].y * j_tan_in[ji].y);

			if (QAWS_FABS(j_angle[ji]) > QAWS_LITERAL(0.087)
				&& j_d_out[ji] > QAWS_LITERAL(1e-12)
				&& j_d_in[ji] > QAWS_LITERAL(1e-12))
			{
				half_angle = QAWS_FABS(j_angle[ji]) * QAWS_LITERAL(0.5);
				if (QAWS_FABS(QAWS_SIN(half_angle)) > QAWS_LITERAL(1e-10))
				{
					j_setback[ji] = radius / QAWS_FABS(QAWS_SIN(half_angle))
						* QAWS_FABS(QAWS_COS(half_angle));
					j_active[ji] = 1;
				}
			}
		}
	}

	/* ---- Pass 2: build trimmed segments and insert arcs ---- */
	total_segs = seg_count * 3;
	new_segments = (qaws_curve**)calloc((size_t)total_segs, sizeof(qaws_curve*));
	if (!new_segments) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup_fillet; }

	for (seg_i = 0; seg_i < seg_count; seg_i++)
	{
		qaws_scalar t_lo = (qaws_scalar)seg_i;
		qaws_scalar t_hi = (qaws_scalar)(seg_i + 1);
		qaws_scalar seg_span = t_hi - t_lo;
		unsigned int n_samp = 32;
		unsigned int si;
		qaws_scalar* samp_params;
		qaws_scalar* samp_coords;

		/* Estimate segment arc length from endpoints */
		{
			qaws_eval_result_2d er_lo, er_hi;
			qaws_scalar seg_len, trim_lo, trim_hi;

			memset(&er_lo, 0, sizeof(er_lo));
			memset(&er_hi, 0, sizeof(er_hi));
			qaws_curve_evaluate_2d(curve, t_lo, QAWS_EVAL_FLAG_POSITION, &er_lo);
			qaws_curve_evaluate_2d(curve, t_hi, QAWS_EVAL_FLAG_POSITION, &er_hi);
			seg_len = QAWS_SQRT(
				(er_hi.position.x - er_lo.position.x) * (er_hi.position.x - er_lo.position.x) +
				(er_hi.position.y - er_lo.position.y) * (er_hi.position.y - er_lo.position.y));
			if (seg_len < QAWS_LITERAL(1e-12)) seg_len = QAWS_LITERAL(1e-12);

			/* Trim from start (previous junction) */
			trim_lo = QAWS_ZERO;
			if (seg_i > 0 && j_active[seg_i - 1])
				trim_lo = j_setback[seg_i - 1] / seg_len * seg_span;

			/* Trim from end (next junction) */
			trim_hi = QAWS_ZERO;
			if (seg_i + 1 < seg_count && j_active[seg_i])
				trim_hi = j_setback[seg_i] / seg_len * seg_span;

			/* Clamp trims so they don't overlap */
			if (trim_lo + trim_hi > seg_span * QAWS_LITERAL(0.9))
			{
				qaws_scalar scale = seg_span * QAWS_LITERAL(0.9) / (trim_lo + trim_hi);
				trim_lo *= scale;
				trim_hi *= scale;
			}

			t_lo += trim_lo;
			t_hi -= trim_hi;
		}

		/* Sample trimmed segment as B-spline copy */
		samp_params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp);
		samp_coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp * 2);
		if (!samp_params || !samp_coords)
		{
			free(samp_params); free(samp_coords);
			goto cleanup_fillet;
		}

		for (si = 0; si < n_samp; si++)
		{
			qaws_eval_result_2d er;
			qaws_scalar t = t_lo + (t_hi - t_lo) * (qaws_scalar)si / (qaws_scalar)(n_samp - 1);
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &er);
			samp_params[si] = (qaws_scalar)si / (qaws_scalar)(n_samp - 1);
			samp_coords[si * 2] = er.position.x;
			samp_coords[si * 2 + 1] = er.position.y;
		}

		{
			qaws_curve* seg_copy = NULL;
			status = qaws_internal_fit_bspline(QAWS_DIMENSION_2D, 3, samp_params, samp_coords, n_samp, 10, &seg_copy);
			free(samp_params); free(samp_coords);
			if (status != QAWS_STATUS_OK) goto cleanup_fillet;
			if (new_count >= total_segs) { qaws_curve_destroy(seg_copy); goto cleanup_fillet; }
			new_segments[new_count++] = seg_copy;
		}

		/* Insert fillet arc at this segment's end junction */
		if (seg_i < n_junctions && j_active[seg_i])
		{
			qaws_vec2 n_out, n_in, fillet_center;
			qaws_arc_segment arc_seg;
			qaws_arc_desc arc_desc;
			qaws_curve* arc_curve = NULL;
			qaws_scalar sign, half_angle;

			half_angle = QAWS_FABS(j_angle[seg_i]) * QAWS_LITERAL(0.5);
			n_out.x = -j_tan_out[seg_i].y / j_d_out[seg_i];
			n_out.y =  j_tan_out[seg_i].x / j_d_out[seg_i];
			n_in.x = -j_tan_in[seg_i].y / j_d_in[seg_i];
			n_in.y =  j_tan_in[seg_i].x / j_d_in[seg_i];

			sign = (j_angle[seg_i] > QAWS_ZERO) ? QAWS_ONE : -QAWS_ONE;
			{
				qaws_scalar bx = n_out.x + n_in.x;
				qaws_scalar by = n_out.y + n_in.y;
				qaws_scalar blen = QAWS_SQRT(bx * bx + by * by);
				qaws_scalar offset;
				if (blen < QAWS_LITERAL(1e-12)) blen = QAWS_LITERAL(1e-12);
				offset = radius / QAWS_FABS(QAWS_SIN(half_angle));
				fillet_center.x = j_corner[seg_i].x + sign * offset * bx / blen;
				fillet_center.y = j_corner[seg_i].y + sign * offset * by / blen;
			}

			memset(&arc_seg, 0, sizeof(arc_seg));
			arc_seg.center[0] = fillet_center.x;
			arc_seg.center[1] = fillet_center.y;
			arc_seg.radius = radius;

			{
				qaws_vec2 p_start, p_end;
				p_start.x = j_corner[seg_i].x - j_setback[seg_i] * j_tan_out[seg_i].x / j_d_out[seg_i];
				p_start.y = j_corner[seg_i].y - j_setback[seg_i] * j_tan_out[seg_i].y / j_d_out[seg_i];
				p_end.x = j_corner[seg_i].x + j_setback[seg_i] * j_tan_in[seg_i].x / j_d_in[seg_i];
				p_end.y = j_corner[seg_i].y + j_setback[seg_i] * j_tan_in[seg_i].y / j_d_in[seg_i];

				arc_seg.angle_start = QAWS_ATAN2(
					p_start.y - fillet_center.y, p_start.x - fillet_center.x);
				arc_seg.angle_end = QAWS_ATAN2(
					p_end.y - fillet_center.y, p_end.x - fillet_center.x);
			}

			memset(&arc_desc, 0, sizeof(arc_desc));
			arc_desc.dimension = QAWS_DIMENSION_2D;
			arc_desc.segments = &arc_seg;
			arc_desc.segment_count = 1;

			status = qaws_curve_create_arc(&arc_desc, &arc_curve);
			if (status == QAWS_STATUS_OK && arc_curve && new_count < total_segs)
				new_segments[new_count++] = arc_curve;
		}
	}

	if (new_count == 0)
	{
		free(new_segments);
		status = QAWS_STATUS_DEGENERATE_CURVE;
		goto cleanup_fillet;
	}

	/* Create composite from all segments */
	memset(&cdesc, 0, sizeof(cdesc));
	cdesc.dimension = QAWS_DIMENSION_2D;
	cdesc.segments = new_segments;
	cdesc.segment_count = new_count;
	status = qaws_curve_create_composite(&cdesc, out_curve);

	if (status != QAWS_STATUS_OK)
		goto cleanup_fillet;

	free(new_segments);
	free(j_setback); free(j_angle); free(j_corner);
	free(j_tan_out); free(j_tan_in); free(j_d_out); free(j_d_in); free(j_active);
	return QAWS_STATUS_OK;

cleanup_fillet:
	{
		unsigned int fi;
		for (fi = 0; fi < new_count; fi++)
			qaws_curve_destroy(new_segments[fi]);
		free(new_segments);
	}
	free(j_setback); free(j_angle); free(j_corner);
	free(j_tan_out); free(j_tan_in); free(j_d_out); free(j_d_in); free(j_active);
	return status;
}

qaws_status qaws_curve_chamfer_2d(
	qaws_curve const* curve,
	qaws_scalar distance,
	qaws_curve** out_curve)
{
	qaws_curve_kind kind;
	unsigned int seg_count, total_segs;
	unsigned int seg_i, ji;
	qaws_curve** new_segments = NULL;
	unsigned int new_count = 0;
	qaws_composite_desc cdesc;
	qaws_status status = QAWS_STATUS_OK;

	/* Per-junction precomputed data */
	qaws_scalar* j_angle = NULL;
	qaws_vec2* j_corner = NULL;
	qaws_vec2* j_tan_out = NULL;
	qaws_vec2* j_tan_in = NULL;
	qaws_scalar* j_d_out = NULL;
	qaws_scalar* j_d_in = NULL;
	int* j_active = NULL;
	unsigned int n_junctions;

	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (distance <= QAWS_ZERO)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(curve) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_curve = NULL;
	kind = qaws_curve_get_kind(curve);

	if (kind != QAWS_CURVE_KIND_COMPOSITE)
	{
		/* Non-composite: just copy */
		qaws_range range = qaws_curve_get_parameter_range(curve);
		unsigned int n_samp = 64;
		qaws_scalar* params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp);
		qaws_scalar* coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp * 2);
		unsigned int i;

		if (!params || !coords) { free(params); free(coords); return QAWS_STATUS_ALLOCATION_FAILURE; }

		for (i = 0; i < n_samp; i++)
		{
			qaws_eval_result_2d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)(n_samp - 1);
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &er);
			params[i] = t;
			coords[i * 2] = er.position.x;
			coords[i * 2 + 1] = er.position.y;
		}
		status = qaws_internal_fit_bspline(QAWS_DIMENSION_2D, 3, params, coords, n_samp, 16, out_curve);
		free(params); free(coords);
		return status;
	}

	/* Composite curve: two-pass approach */
	seg_count = qaws_curve_get_span_count(curve);
	n_junctions = (seg_count > 1) ? seg_count - 1 : 0;

	/* ---- Pass 1: precompute junction geometry ---- */
	if (n_junctions > 0)
	{
		j_angle   = (qaws_scalar*)calloc((size_t)n_junctions, sizeof(qaws_scalar));
		j_corner  = (qaws_vec2*)calloc((size_t)n_junctions, sizeof(qaws_vec2));
		j_tan_out = (qaws_vec2*)calloc((size_t)n_junctions, sizeof(qaws_vec2));
		j_tan_in  = (qaws_vec2*)calloc((size_t)n_junctions, sizeof(qaws_vec2));
		j_d_out   = (qaws_scalar*)calloc((size_t)n_junctions, sizeof(qaws_scalar));
		j_d_in    = (qaws_scalar*)calloc((size_t)n_junctions, sizeof(qaws_scalar));
		j_active  = (int*)calloc((size_t)n_junctions, sizeof(int));
		if (!j_angle || !j_corner || !j_tan_out || !j_tan_in
			|| !j_d_out || !j_d_in || !j_active)
		{
			free(j_angle); free(j_corner);
			free(j_tan_out); free(j_tan_in); free(j_d_out); free(j_d_in); free(j_active);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		for (ji = 0; ji < n_junctions; ji++)
		{
			qaws_eval_result_2d er_end, er_start, er_corner;
			qaws_scalar t_junc = (qaws_scalar)(ji + 1);

			memset(&er_end, 0, sizeof(er_end));
			memset(&er_start, 0, sizeof(er_start));
			memset(&er_corner, 0, sizeof(er_corner));
			qaws_curve_evaluate_2d(curve, t_junc - QAWS_LITERAL(1e-6), QAWS_EVAL_FLAG_D1, &er_end);
			qaws_curve_evaluate_2d(curve, t_junc + QAWS_LITERAL(1e-6), QAWS_EVAL_FLAG_D1, &er_start);
			qaws_curve_evaluate_2d(curve, t_junc, QAWS_EVAL_FLAG_POSITION, &er_corner);

			j_tan_out[ji] = er_end.d1;
			j_tan_in[ji] = er_start.d1;
			j_corner[ji] = er_corner.position;
			j_angle[ji] = vec2_angle(er_end.d1, er_start.d1);

			j_d_out[ji] = QAWS_SQRT(j_tan_out[ji].x * j_tan_out[ji].x + j_tan_out[ji].y * j_tan_out[ji].y);
			j_d_in[ji] = QAWS_SQRT(j_tan_in[ji].x * j_tan_in[ji].x + j_tan_in[ji].y * j_tan_in[ji].y);

			if (QAWS_FABS(j_angle[ji]) > QAWS_LITERAL(0.087)
				&& j_d_out[ji] > QAWS_LITERAL(1e-12)
				&& j_d_in[ji] > QAWS_LITERAL(1e-12))
			{
				j_active[ji] = 1;
			}
		}
	}

	/* ---- Pass 2: build trimmed segments and insert chamfer lines ---- */
	total_segs = seg_count * 2;
	new_segments = (qaws_curve**)calloc((size_t)total_segs, sizeof(qaws_curve*));
	if (!new_segments) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup_chamfer; }

	for (seg_i = 0; seg_i < seg_count; seg_i++)
	{
		qaws_scalar t_lo = (qaws_scalar)seg_i;
		qaws_scalar t_hi = (qaws_scalar)(seg_i + 1);
		qaws_scalar seg_span = t_hi - t_lo;
		unsigned int n_samp = 32;
		unsigned int si;
		qaws_scalar* samp_params;
		qaws_scalar* samp_coords;

		/* Trim segment parameter range based on chamfer distance */
		{
			qaws_eval_result_2d er_lo, er_hi;
			qaws_scalar seg_len, trim_lo, trim_hi;

			memset(&er_lo, 0, sizeof(er_lo));
			memset(&er_hi, 0, sizeof(er_hi));
			qaws_curve_evaluate_2d(curve, t_lo, QAWS_EVAL_FLAG_POSITION, &er_lo);
			qaws_curve_evaluate_2d(curve, t_hi, QAWS_EVAL_FLAG_POSITION, &er_hi);
			seg_len = QAWS_SQRT(
				(er_hi.position.x - er_lo.position.x) * (er_hi.position.x - er_lo.position.x) +
				(er_hi.position.y - er_lo.position.y) * (er_hi.position.y - er_lo.position.y));
			if (seg_len < QAWS_LITERAL(1e-12)) seg_len = QAWS_LITERAL(1e-12);

			trim_lo = QAWS_ZERO;
			if (seg_i > 0 && j_active[seg_i - 1])
				trim_lo = distance / seg_len * seg_span;

			trim_hi = QAWS_ZERO;
			if (seg_i + 1 < seg_count && j_active[seg_i])
				trim_hi = distance / seg_len * seg_span;

			if (trim_lo + trim_hi > seg_span * QAWS_LITERAL(0.9))
			{
				qaws_scalar scale = seg_span * QAWS_LITERAL(0.9) / (trim_lo + trim_hi);
				trim_lo *= scale;
				trim_hi *= scale;
			}

			t_lo += trim_lo;
			t_hi -= trim_hi;
		}

		/* Copy trimmed segment as B-spline */
		samp_params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp);
		samp_coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samp * 2);
		if (!samp_params || !samp_coords)
		{
			free(samp_params); free(samp_coords);
			goto cleanup_chamfer;
		}

		for (si = 0; si < n_samp; si++)
		{
			qaws_eval_result_2d er;
			qaws_scalar t = t_lo + (t_hi - t_lo) * (qaws_scalar)si / (qaws_scalar)(n_samp - 1);
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &er);
			samp_params[si] = (qaws_scalar)si / (qaws_scalar)(n_samp - 1);
			samp_coords[si * 2] = er.position.x;
			samp_coords[si * 2 + 1] = er.position.y;
		}

		{
			qaws_curve* seg_copy = NULL;
			status = qaws_internal_fit_bspline(QAWS_DIMENSION_2D, 3, samp_params, samp_coords, n_samp, 10, &seg_copy);
			free(samp_params); free(samp_coords);
			if (status != QAWS_STATUS_OK) goto cleanup_chamfer;
			if (new_count >= total_segs) { qaws_curve_destroy(seg_copy); goto cleanup_chamfer; }
			new_segments[new_count++] = seg_copy;
		}

		/* Insert chamfer line at junction */
		if (seg_i < n_junctions && j_active[seg_i])
		{
			qaws_vec2 p0, p1;
			qaws_scalar pts[4];
			qaws_bezier_desc bdesc;
			qaws_curve* line = NULL;

			p0.x = j_corner[seg_i].x - distance * j_tan_out[seg_i].x / j_d_out[seg_i];
			p0.y = j_corner[seg_i].y - distance * j_tan_out[seg_i].y / j_d_out[seg_i];
			p1.x = j_corner[seg_i].x + distance * j_tan_in[seg_i].x / j_d_in[seg_i];
			p1.y = j_corner[seg_i].y + distance * j_tan_in[seg_i].y / j_d_in[seg_i];

			pts[0] = p0.x; pts[1] = p0.y;
			pts[2] = p1.x; pts[3] = p1.y;

			memset(&bdesc, 0, sizeof(bdesc));
			bdesc.dimension = QAWS_DIMENSION_2D;
			bdesc.degree = 1;
			bdesc.control_points = pts;
			bdesc.control_point_count = 2;

			status = qaws_curve_create_bezier(&bdesc, &line);
			if (status == QAWS_STATUS_OK && line && new_count < total_segs)
				new_segments[new_count++] = line;
		}
	}

	if (new_count == 0)
	{
		free(new_segments);
		status = QAWS_STATUS_DEGENERATE_CURVE;
		goto cleanup_chamfer;
	}

	memset(&cdesc, 0, sizeof(cdesc));
	cdesc.dimension = QAWS_DIMENSION_2D;
	cdesc.segments = new_segments;
	cdesc.segment_count = new_count;
	status = qaws_curve_create_composite(&cdesc, out_curve);

	if (status != QAWS_STATUS_OK)
		goto cleanup_chamfer;

	free(new_segments);
	free(j_angle); free(j_corner);
	free(j_tan_out); free(j_tan_in); free(j_d_out); free(j_d_in); free(j_active);
	return QAWS_STATUS_OK;

cleanup_chamfer:
	{
		unsigned int fi;
		for (fi = 0; fi < new_count; fi++)
			qaws_curve_destroy(new_segments[fi]);
		free(new_segments);
	}
	free(j_angle); free(j_corner);
	free(j_tan_out); free(j_tan_in); free(j_d_out); free(j_d_in); free(j_active);
	return status;
}
