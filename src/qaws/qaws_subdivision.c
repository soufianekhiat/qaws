#include "qaws_subdivision.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* -------------------------------------------------------------------------- */
/*  Impl struct                                                               */
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
/*  Subdivision algorithms                                                    */
/* -------------------------------------------------------------------------- */

/*
 * Chaikin subdivision (one step).
 * Open curve: keeps endpoints, produces 2*(n-1)+1 = 2*n-1 points.
 * Closed curve: produces 2*n points.
 */
static qaws_scalar* chaikin_step(
	qaws_scalar const* pts,
	unsigned int count,
	unsigned int dim_count,
	int closed,
	unsigned int* out_count)
{
	qaws_scalar* result;
	unsigned int i, d;

	if (closed) {
		unsigned int new_count = 2 * count;
		result = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(new_count * dim_count));
		if (!result) return NULL;

		for (i = 0; i < count; i++) {
			unsigned int next = (i + 1) % count;
			for (d = 0; d < dim_count; d++) {
				result[(2 * i) * dim_count + d] =
					(qaws_scalar)0.75 * pts[i * dim_count + d] +
					(qaws_scalar)0.25 * pts[next * dim_count + d];
				result[(2 * i + 1) * dim_count + d] =
					(qaws_scalar)0.25 * pts[i * dim_count + d] +
					(qaws_scalar)0.75 * pts[next * dim_count + d];
			}
		}

		*out_count = new_count;
	} else {
		unsigned int new_count = 2 * (count - 1) + 1;
		result = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(new_count * dim_count));
		if (!result) return NULL;

		/* Keep first endpoint */
		for (d = 0; d < dim_count; d++)
			result[0 * dim_count + d] = pts[0 * dim_count + d];

		/* Interior cutting points */
		for (i = 0; i < count - 1; i++) {
			for (d = 0; d < dim_count; d++) {
				result[(2 * i + 1) * dim_count + d] =
					(qaws_scalar)0.75 * pts[i * dim_count + d] +
					(qaws_scalar)0.25 * pts[(i + 1) * dim_count + d];
				result[(2 * i + 2) * dim_count + d] =
					(qaws_scalar)0.25 * pts[i * dim_count + d] +
					(qaws_scalar)0.75 * pts[(i + 1) * dim_count + d];
			}
		}

		/* Keep last endpoint */
		for (d = 0; d < dim_count; d++)
			result[(new_count - 1) * dim_count + d] =
				pts[(count - 1) * dim_count + d];

		*out_count = new_count;
	}

	return result;
}

/*
 * Lane-Riesenfeld subdivision (one step, order k).
 *
 * Step 1: Point insertion (midpoint doubling).
 * Step 2: k rounds of averaging.
 *
 * For degree 3 (cubic B-spline limit): k = 2.
 * For degree 4 (quartic B-spline limit): k = 3.
 */
static qaws_scalar* lane_riesenfeld_step(
	qaws_scalar const* pts,
	unsigned int count,
	unsigned int dim_count,
	int closed,
	unsigned int averaging_rounds,
	unsigned int* out_count)
{
	qaws_scalar* doubled;
	unsigned int doubled_count;
	unsigned int r, i, d;

	/* --- Step 1: Point insertion (midpoint doubling) ------------------- */

	if (closed) {
		doubled_count = 2 * count;
		doubled = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(doubled_count * dim_count));
		if (!doubled) return NULL;

		for (i = 0; i < count; i++) {
			unsigned int next = (i + 1) % count;
			for (d = 0; d < dim_count; d++) {
				/* Keep original point */
				doubled[(2 * i) * dim_count + d] = pts[i * dim_count + d];
				/* Insert midpoint */
				doubled[(2 * i + 1) * dim_count + d] =
					(qaws_scalar)0.5 * (pts[i * dim_count + d] +
					                    pts[next * dim_count + d]);
			}
		}
	} else {
		doubled_count = 2 * count - 1;
		doubled = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(doubled_count * dim_count));
		if (!doubled) return NULL;

		for (i = 0; i < count - 1; i++) {
			for (d = 0; d < dim_count; d++) {
				doubled[(2 * i) * dim_count + d] = pts[i * dim_count + d];
				doubled[(2 * i + 1) * dim_count + d] =
					(qaws_scalar)0.5 * (pts[i * dim_count + d] +
					                    pts[(i + 1) * dim_count + d]);
			}
		}
		/* Last original point */
		for (d = 0; d < dim_count; d++)
			doubled[(doubled_count - 1) * dim_count + d] =
				pts[(count - 1) * dim_count + d];
	}

	/* --- Step 2: Averaging rounds -------------------------------------- */

	for (r = 0; r < averaging_rounds; r++) {
		if (closed) {
			/* Average all points with their next neighbor (wrapping) */
			qaws_scalar* temp = (qaws_scalar*)malloc(
				sizeof(qaws_scalar) * (size_t)(doubled_count * dim_count));
			if (!temp) { free(doubled); return NULL; }

			for (i = 0; i < doubled_count; i++) {
				unsigned int next = (i + 1) % doubled_count;
				for (d = 0; d < dim_count; d++) {
					temp[i * dim_count + d] =
						(qaws_scalar)0.5 * (doubled[i * dim_count + d] +
						                    doubled[next * dim_count + d]);
				}
			}

			free(doubled);
			doubled = temp;
		} else {
			/* Average interior points; keep endpoints fixed */
			qaws_scalar* temp = (qaws_scalar*)malloc(
				sizeof(qaws_scalar) * (size_t)(doubled_count * dim_count));
			if (!temp) { free(doubled); return NULL; }

			/* Keep first endpoint */
			for (d = 0; d < dim_count; d++)
				temp[0 * dim_count + d] = doubled[0 * dim_count + d];

			/* Average interior */
			for (i = 1; i < doubled_count - 1; i++) {
				for (d = 0; d < dim_count; d++) {
					temp[i * dim_count + d] =
						(qaws_scalar)0.5 * (doubled[i * dim_count + d] +
						                    doubled[(i + 1) * dim_count + d]);
				}
			}

			/* Keep last endpoint */
			for (d = 0; d < dim_count; d++)
				temp[(doubled_count - 1) * dim_count + d] =
					doubled[(doubled_count - 1) * dim_count + d];

			free(doubled);
			doubled = temp;
		}
	}

	*out_count = doubled_count;
	return doubled;
}

/*
 * Apply a subdivision scheme for the given number of levels.
 * Returns the refined point array (caller frees) and refined count.
 */
static qaws_scalar* apply_subdivision(
	qaws_scalar const* input_pts,
	unsigned int input_count,
	unsigned int dim_count,
	int closed,
	qaws_subdivision_scheme scheme,
	unsigned int levels,
	unsigned int* out_count)
{
	qaws_scalar* current;
	unsigned int current_count;
	unsigned int lev;

	/* Copy input so we can iterate */
	current_count = input_count;
	current = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(current_count * dim_count));
	if (!current) return NULL;
	memcpy(current, input_pts, sizeof(qaws_scalar) * (size_t)(current_count * dim_count));

	for (lev = 0; lev < levels; lev++) {
		qaws_scalar* next = NULL;
		unsigned int next_count = 0;

		switch (scheme) {
		case QAWS_SUBDIVISION_CHAIKIN:
			next = chaikin_step(current, current_count, dim_count,
			                    closed, &next_count);
			break;
		case QAWS_SUBDIVISION_LANE_RIESENFELD_3:
			next = lane_riesenfeld_step(current, current_count, dim_count,
			                            closed, 2, &next_count);
			break;
		case QAWS_SUBDIVISION_LANE_RIESENFELD_4:
			next = lane_riesenfeld_step(current, current_count, dim_count,
			                            closed, 3, &next_count);
			break;
		}

		free(current);
		if (!next) return NULL;

		current = next;
		current_count = next_count;
	}

	*out_count = current_count;
	return current;
}

/* -------------------------------------------------------------------------- */
/*  Uniform Catmull-Rom coefficient computation                               */
/* -------------------------------------------------------------------------- */

/*
 * Compute uniform Catmull-Rom cubic coefficients for a span between P1 and P2,
 * using the four-point stencil (P0, P1, P2, P3).
 *
 * P(t) = a*t^3 + b*t^2 + c*t + d,  t in [0,1]
 *
 *   a = -0.5*P0 + 1.5*P1 - 1.5*P2 + 0.5*P3
 *   b =      P0 - 2.5*P1 + 2.0*P2 - 0.5*P3
 *   c = -0.5*P0          + 0.5*P2
 *   d =               P1
 */
static void compute_uniform_cr_coeffs(
	qaws_scalar const* P0,
	qaws_scalar const* P1,
	qaws_scalar const* P2,
	qaws_scalar const* P3,
	unsigned int dim_count,
	qaws_scalar* coeffs_out)  /* dim_count * 4 scalars */
{
	unsigned int d;

	for (d = 0; d < dim_count; d++) {
		qaws_scalar p0 = P0[d];
		qaws_scalar p1 = P1[d];
		qaws_scalar p2 = P2[d];
		qaws_scalar p3 = P3[d];

		coeffs_out[d * 4 + 0] = (qaws_scalar)-0.5 * p0 + (qaws_scalar)1.5 * p1
		                       - (qaws_scalar)1.5 * p2 + (qaws_scalar)0.5 * p3;
		coeffs_out[d * 4 + 1] = p0 - (qaws_scalar)2.5 * p1
		                       + (qaws_scalar)2.0 * p2 - (qaws_scalar)0.5 * p3;
		coeffs_out[d * 4 + 2] = (qaws_scalar)-0.5 * p0 + (qaws_scalar)0.5 * p2;
		coeffs_out[d * 4 + 3] = p1;
	}
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (2D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status subdivision_eval_span_2d(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d* out_result)
{
	qaws_subdivision_impl const* impl =
		(qaws_subdivision_impl const*)curve->impl;
	unsigned int const dim_count = 2;
	qaws_scalar t = local_t;
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;
	qaws_scalar const* cx;
	qaws_scalar const* cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = ax * t3 + bx * t2 + cx_coeff * t + dx_val;
		out_result->position.y = ay * t3 + by * t2 + cy_coeff * t + dy_val;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = (qaws_scalar)3.0 * ax * t2 + (qaws_scalar)2.0 * bx * t + cx_coeff;
		out_result->d1.y = (qaws_scalar)3.0 * ay * t2 + (qaws_scalar)2.0 * by * t + cy_coeff;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = (qaws_scalar)6.0 * ax * t + (qaws_scalar)2.0 * bx;
		out_result->d2.y = (qaws_scalar)6.0 * ay * t + (qaws_scalar)2.0 * by;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = (qaws_scalar)6.0 * ax;
		out_result->d3.y = (qaws_scalar)6.0 * ay;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (3D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status subdivision_eval_span_3d(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d* out_result)
{
	qaws_subdivision_impl const* impl =
		(qaws_subdivision_impl const*)curve->impl;
	unsigned int const dim_count = 3;
	qaws_scalar t = local_t;
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;
	qaws_scalar const* cx;
	qaws_scalar const* cy;
	qaws_scalar const* cz;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar az, bz, cz_coeff, dz_val;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];
	cz = &impl->segment_coeffs[span_index * dim_count * 4 + 2 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];
	az = cz[0]; bz = cz[1]; cz_coeff = cz[2]; dz_val = cz[3];

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = ax * t3 + bx * t2 + cx_coeff * t + dx_val;
		out_result->position.y = ay * t3 + by * t2 + cy_coeff * t + dy_val;
		out_result->position.z = az * t3 + bz * t2 + cz_coeff * t + dz_val;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = (qaws_scalar)3.0 * ax * t2 + (qaws_scalar)2.0 * bx * t + cx_coeff;
		out_result->d1.y = (qaws_scalar)3.0 * ay * t2 + (qaws_scalar)2.0 * by * t + cy_coeff;
		out_result->d1.z = (qaws_scalar)3.0 * az * t2 + (qaws_scalar)2.0 * bz * t + cz_coeff;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = (qaws_scalar)6.0 * ax * t + (qaws_scalar)2.0 * bx;
		out_result->d2.y = (qaws_scalar)6.0 * ay * t + (qaws_scalar)2.0 * by;
		out_result->d2.z = (qaws_scalar)6.0 * az * t + (qaws_scalar)2.0 * bz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = (qaws_scalar)6.0 * ax;
		out_result->d3.y = (qaws_scalar)6.0 * ay;
		out_result->d3.z = (qaws_scalar)6.0 * az;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: lifecycle and query functions                                     */
/* -------------------------------------------------------------------------- */

static void subdivision_destroy_impl(void* impl, qaws_allocator const* allocator)
{
	qaws_subdivision_impl* si = (qaws_subdivision_impl*)impl;
	if (si) {
		qaws_internal_dealloc(allocator, si->refined_points);
		qaws_internal_dealloc(allocator, si->segment_coeffs);
		qaws_internal_dealloc(allocator, si);
	}
}

static int subdivision_is_closed(qaws_curve const* curve)
{
	qaws_subdivision_impl const* impl =
		(qaws_subdivision_impl const*)curve->impl;
	return impl->closed;
}

static int subdivision_is_periodic(qaws_curve const* curve)
{
	qaws_subdivision_impl const* impl =
		(qaws_subdivision_impl const*)curve->impl;
	return impl->closed;
}

static int subdivision_is_rational(qaws_curve const* curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity subdivision_get_continuity(qaws_curve const* curve)
{
	qaws_subdivision_impl const* impl =
		(qaws_subdivision_impl const*)curve->impl;

	switch (impl->scheme) {
	case QAWS_SUBDIVISION_CHAIKIN:
		return QAWS_CONTINUITY_C1;
	case QAWS_SUBDIVISION_LANE_RIESENFELD_3:
		return QAWS_CONTINUITY_C2;
	case QAWS_SUBDIVISION_LANE_RIESENFELD_4:
		return QAWS_CONTINUITY_C3;
	default:
		return QAWS_CONTINUITY_C1;
	}
}

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                    */
/* -------------------------------------------------------------------------- */

static qaws_curve_vtable const subdivision_vtable = {
	subdivision_eval_span_2d,
	subdivision_eval_span_3d,
	subdivision_destroy_impl,
	subdivision_is_closed,
	subdivision_is_periodic,
	subdivision_is_rational,
	subdivision_get_continuity,
	NULL /* diff */
};

/* -------------------------------------------------------------------------- */
/*  Helper: get scheme degree                                                 */
/* -------------------------------------------------------------------------- */

static unsigned int scheme_degree(qaws_subdivision_scheme scheme)
{
	switch (scheme) {
	case QAWS_SUBDIVISION_CHAIKIN:            return 2;
	case QAWS_SUBDIVISION_LANE_RIESENFELD_3: return 3;
	case QAWS_SUBDIVISION_LANE_RIESENFELD_4: return 4;
	default:                                  return 2;
	}
}

/* -------------------------------------------------------------------------- */
/*  Creation                                                                  */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_subdivision(
	qaws_subdivision_desc const* desc,
	qaws_curve** out_curve)
{
	qaws_subdivision_impl* impl = NULL;
	qaws_curve* curve = NULL;
	unsigned int dim_count;
	unsigned int span_count;
	unsigned int refined_count;
	unsigned int levels;
	unsigned int i, s;
	qaws_scalar* refined;
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

	if (desc->control_point_count < 3)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (desc->scheme != QAWS_SUBDIVISION_CHAIKIN &&
	    desc->scheme != QAWS_SUBDIVISION_LANE_RIESENFELD_3 &&
	    desc->scheme != QAWS_SUBDIVISION_LANE_RIESENFELD_4)
		return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	levels = desc->refinement_levels;
	if (levels == 0)
		levels = 6;  /* default */

	/* --- Apply subdivision --------------------------------------------- */

	refined = apply_subdivision(
		(qaws_scalar const*)desc->control_points,
		desc->control_point_count,
		dim_count,
		desc->closed ? 1 : 0,
		desc->scheme,
		levels,
		&refined_count);

	if (!refined)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Need at least 4 points for Catmull-Rom stencil (open) or 3 (closed) */
	if ((!desc->closed && refined_count < 4) ||
	    (desc->closed && refined_count < 3)) {
		free(refined);
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	}

	/* --- Compute span count -------------------------------------------- */

	if (desc->closed) {
		span_count = refined_count;
	} else {
		span_count = refined_count - 1;
	}

	/* --- Allocate implementation --------------------------------------- */

	impl = (qaws_subdivision_impl*)calloc(1, sizeof(qaws_subdivision_impl));
	if (!impl) {
		free(refined);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->refined_points = refined;
	impl->refined_count = refined_count;
	impl->closed = desc->closed ? 1 : 0;
	impl->scheme = desc->scheme;

	/* --- Compute Catmull-Rom coefficients for each span ---------------- */

	impl->segment_coeffs = (qaws_scalar*)calloc(
		(size_t)(span_count * dim_count * 4), sizeof(qaws_scalar));
	if (!impl->segment_coeffs) {
		free(impl->refined_points);
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	if (desc->closed) {
		/*
		 * Closed curve: span s goes from refined_points[s] to
		 * refined_points[(s+1)%N]. Stencil uses wrapping indices.
		 */
		unsigned int N = refined_count;
		for (s = 0; s < span_count; s++) {
			unsigned int i0 = (s + N - 1) % N;
			unsigned int i1 = s;
			unsigned int i2 = (s + 1) % N;
			unsigned int i3 = (s + 2) % N;

			compute_uniform_cr_coeffs(
				&refined[i0 * dim_count],
				&refined[i1 * dim_count],
				&refined[i2 * dim_count],
				&refined[i3 * dim_count],
				dim_count,
				&impl->segment_coeffs[s * dim_count * 4]);
		}
	} else {
		/*
		 * Open curve: span s goes from refined_points[s] to
		 * refined_points[s+1]. For the first and last spans we need
		 * phantom points. We reflect through the endpoint:
		 *   phantom_start = 2*P[0] - P[1]
		 *   phantom_end   = 2*P[N-1] - P[N-2]
		 */
		unsigned int N = refined_count;
		qaws_scalar* phantom_start = NULL;
		qaws_scalar* phantom_end = NULL;
		unsigned int d;

		phantom_start = (qaws_scalar*)malloc(sizeof(qaws_scalar) * dim_count);
		phantom_end = (qaws_scalar*)malloc(sizeof(qaws_scalar) * dim_count);
		if (!phantom_start || !phantom_end) {
			free(phantom_start);
			free(phantom_end);
			free(impl->segment_coeffs);
			free(impl->refined_points);
			free(impl);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		for (d = 0; d < dim_count; d++) {
			phantom_start[d] = (qaws_scalar)2.0 * refined[0 * dim_count + d]
			                 - refined[1 * dim_count + d];
			phantom_end[d] = (qaws_scalar)2.0 * refined[(N - 1) * dim_count + d]
			               - refined[(N - 2) * dim_count + d];
		}

		for (s = 0; s < span_count; s++) {
			qaws_scalar const* P0;
			qaws_scalar const* P1;
			qaws_scalar const* P2;
			qaws_scalar const* P3;

			/* P1 = refined[s], P2 = refined[s+1] */
			P1 = &refined[s * dim_count];
			P2 = &refined[(s + 1) * dim_count];

			/* P0: predecessor of P1 */
			if (s == 0)
				P0 = phantom_start;
			else
				P0 = &refined[(s - 1) * dim_count];

			/* P3: successor of P2 */
			if (s + 2 >= N)
				P3 = phantom_end;
			else
				P3 = &refined[(s + 2) * dim_count];

			compute_uniform_cr_coeffs(
				P0, P1, P2, P3,
				dim_count,
				&impl->segment_coeffs[s * dim_count * 4]);
		}

		free(phantom_start);
		free(phantom_end);
	}

	/* --- Allocate curve ------------------------------------------------ */

	range.min_value = (qaws_scalar)0.0;
	range.max_value = (qaws_scalar)span_count;

	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_SUBDIVISION,
		desc->dimension,
		scheme_degree(desc->scheme),
		span_count,
		range,
		&subdivision_vtable);
	if (!curve) {
		free(impl->segment_coeffs);
		free(impl->refined_points);
		free(impl);
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
