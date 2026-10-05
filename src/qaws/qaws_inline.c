#include "qaws_inline.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include <string.h>
#include <math.h>

/*
 * Inline curves store everything in a flat byte buffer (qaws_curve_inline.storage).
 *
 * Layout:
 *   [0]                       : qaws_curve struct
 *   [sizeof(qaws_curve)]      : span_bounds[2] (qaws_scalar)
 *   [after span_bounds]       : control point / coefficient data
 *   [after CP data]           : impl struct (qaws_bezier_impl or inline_poly_impl)
 *
 * The init functions placement-initialize everything. No heap allocation.
 */

/* Helper: get aligned offset */
static size_t align_up(size_t offset, size_t alignment)
{
	return (offset + alignment - 1) & ~(alignment - 1);
}

qaws_curve* qaws_curve_inline_get_curve(qaws_curve_inline* ic)
{
	if (!ic) return NULL;
	return (qaws_curve*)(void*)ic->storage;
}

/* ------------------------------------------------------------------ */
/*  Inline Bezier                                                      */
/* ------------------------------------------------------------------ */

static qaws_status inline_bezier_eval_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* cp = impl->control_points;
	qaws_scalar buf[2], d1p[128], d2p[128], d3p[128];
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_internal_decasteljau(cp, degree, 2, local_t, buf);
		out_result->position.x = buf[0]; out_result->position.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && degree >= 1) {
		qaws_internal_bezier_derivative_points(cp, degree, 2, d1p);
		qaws_internal_decasteljau(d1p, degree - 1, 2, local_t, buf);
		out_result->d1.x = buf[0]; out_result->d1.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && degree >= 2) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 2, d1p);
		qaws_internal_bezier_derivative_points(d1p, degree - 1, 2, d2p);
		qaws_internal_decasteljau(d2p, degree - 2, 2, local_t, buf);
		out_result->d2.x = buf[0]; out_result->d2.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 2, d1p);
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D2))
			qaws_internal_bezier_derivative_points(d1p, degree - 1, 2, d2p);
		qaws_internal_bezier_derivative_points(d2p, degree - 2, 2, d3p);
		qaws_internal_decasteljau(d3p, degree - 3, 2, local_t, buf);
		out_result->d3.x = buf[0]; out_result->d3.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}
	return QAWS_STATUS_OK;
}

static qaws_status inline_bezier_eval_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* cp = impl->control_points;
	qaws_scalar buf[3], d1p[192], d2p[192], d3p[192];
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_internal_decasteljau(cp, degree, 3, local_t, buf);
		out_result->position.x = buf[0]; out_result->position.y = buf[1]; out_result->position.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && degree >= 1) {
		qaws_internal_bezier_derivative_points(cp, degree, 3, d1p);
		qaws_internal_decasteljau(d1p, degree - 1, 3, local_t, buf);
		out_result->d1.x = buf[0]; out_result->d1.y = buf[1]; out_result->d1.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && degree >= 2) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 3, d1p);
		qaws_internal_bezier_derivative_points(d1p, degree - 1, 3, d2p);
		qaws_internal_decasteljau(d2p, degree - 2, 3, local_t, buf);
		out_result->d2.x = buf[0]; out_result->d2.y = buf[1]; out_result->d2.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 3, d1p);
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D2))
			qaws_internal_bezier_derivative_points(d1p, degree - 1, 3, d2p);
		qaws_internal_bezier_derivative_points(d2p, degree - 2, 3, d3p);
		qaws_internal_decasteljau(d3p, degree - 3, 3, local_t, buf);
		out_result->d3.x = buf[0]; out_result->d3.y = buf[1]; out_result->d3.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}
	return QAWS_STATUS_OK;
}

/* No-op destroy: inline curves have no heap allocation. */
static void inline_noop_destroy(void* impl, qaws_allocator const* allocator)
{
	(void)impl;
	(void)allocator;
}

static int inline_not_closed(qaws_curve const* c)    { (void)c; return 0; }
static int inline_not_periodic(qaws_curve const* c)   { (void)c; return 0; }
static int inline_not_rational(qaws_curve const* c)   { (void)c; return 0; }
static qaws_continuity inline_c3(qaws_curve const* c) { (void)c; return QAWS_CONTINUITY_C3; }

static qaws_curve_vtable const inline_bezier_vtable = {
	inline_bezier_eval_2d, inline_bezier_eval_3d, inline_noop_destroy,
	inline_not_closed, inline_not_periodic, inline_not_rational, inline_c3,
	NULL /* diff */
};

qaws_status qaws_curve_init_bezier_inline(
	qaws_bezier_desc const* desc,
	qaws_curve_inline* out)
{
	unsigned int dim_count, n_scalars;
	qaws_curve* curve;
	qaws_scalar* span_bounds;
	qaws_scalar* cp_data;
	qaws_bezier_impl* impl;
	size_t offset;

	if (!desc || !out) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (desc->degree < 1 || desc->degree > 7) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->control_point_count != desc->degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points) return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (unsigned int)desc->dimension;
	n_scalars = dim_count * desc->control_point_count;

	/* Compute layout offsets and check they fit */
	offset = align_up(sizeof(qaws_curve), sizeof(qaws_scalar));
	/* span_bounds: 2 scalars */
	{
		size_t sb_off = offset;
		offset += 2 * sizeof(qaws_scalar);
		/* CP data */
		offset = align_up(offset, sizeof(qaws_scalar));
		{
			size_t cp_off = offset;
			offset += (size_t)n_scalars * sizeof(qaws_scalar);
			/* impl struct */
			offset = align_up(offset, sizeof(void*));
			{
				size_t impl_off = offset;
				offset += sizeof(qaws_bezier_impl);

				if (offset > QAWS_INLINE_STORAGE_BYTES)
					return QAWS_STATUS_BUFFER_TOO_SMALL;

				/* Zero out and set up pointers */
				memset(out->storage, 0, QAWS_INLINE_STORAGE_BYTES);

				curve = (qaws_curve*)(void*)out->storage;
				span_bounds = (qaws_scalar*)(void*)&out->storage[sb_off];
				cp_data = (qaws_scalar*)(void*)&out->storage[cp_off];
				impl = (qaws_bezier_impl*)(void*)&out->storage[impl_off];

				/* Copy control points */
				memcpy(cp_data, desc->control_points,
					(size_t)n_scalars * sizeof(qaws_scalar));

				/* Set up impl */
				impl->control_points = cp_data;
				impl->control_point_count = desc->control_point_count;

				/* Set up span boundaries */
				span_bounds[0] = 0;
				span_bounds[1] = 1;

				/* Initialize curve struct */
				curve->kind = QAWS_CURVE_KIND_BEZIER;
				curve->dimension = desc->dimension;
				curve->degree = desc->degree;
				curve->span_count = 1;
				curve->parameter_range.min_value = 0;
				curve->parameter_range.max_value = 1;
				curve->span_boundaries = span_bounds;
				curve->vtable = &inline_bezier_vtable;
				curve->impl = impl;
				curve->allocator = NULL;
			}
		}
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Inline Polynomial                                                  */
/* ------------------------------------------------------------------ */

typedef struct inline_poly_impl
{
	qaws_scalar const* coefficients;
	unsigned int coefficient_count;
	unsigned int dim_count;
} inline_poly_impl;

static qaws_status inline_poly_eval_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	inline_poly_impl const* pi = (inline_poly_impl const*)curve->impl;
	unsigned int deg = curve->degree;
	qaws_scalar const* c = pi->coefficients;
	qaws_scalar t_min = curve->parameter_range.min_value;
	qaws_scalar t_max = curve->parameter_range.max_value;
	qaws_scalar t = t_min + local_t * (t_max - t_min);
	unsigned int i;
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_scalar px = c[deg * 2 + 0], py = c[deg * 2 + 1];
		for (i = deg; i > 0; i--) {
			px = px * t + c[(i - 1) * 2 + 0];
			py = py * t + c[(i - 1) * 2 + 1];
		}
		out_result->position.x = px; out_result->position.y = py;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && deg >= 1) {
		qaws_scalar dx = (qaws_scalar)deg * c[deg * 2 + 0];
		qaws_scalar dy = (qaws_scalar)deg * c[deg * 2 + 1];
		for (i = deg; i > 1; i--) {
			dx = dx * t + (qaws_scalar)(i - 1) * c[(i - 1) * 2 + 0];
			dy = dy * t + (qaws_scalar)(i - 1) * c[(i - 1) * 2 + 1];
		}
		out_result->d1.x = dx; out_result->d1.y = dy;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && deg >= 2) {
		qaws_scalar dx = (qaws_scalar)(deg * (deg - 1)) * c[deg * 2 + 0];
		qaws_scalar dy = (qaws_scalar)(deg * (deg - 1)) * c[deg * 2 + 1];
		for (i = deg; i > 2; i--) {
			dx = dx * t + (qaws_scalar)((i - 1) * (i - 2)) * c[(i - 1) * 2 + 0];
			dy = dy * t + (qaws_scalar)((i - 1) * (i - 2)) * c[(i - 1) * 2 + 1];
		}
		out_result->d2.x = dx; out_result->d2.y = dy;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	return QAWS_STATUS_OK;
}

static qaws_status inline_poly_eval_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	inline_poly_impl const* pi = (inline_poly_impl const*)curve->impl;
	unsigned int deg = curve->degree;
	qaws_scalar const* c = pi->coefficients;
	qaws_scalar t_min = curve->parameter_range.min_value;
	qaws_scalar t_max = curve->parameter_range.max_value;
	qaws_scalar t = t_min + local_t * (t_max - t_min);
	unsigned int i;
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_scalar px = c[deg * 3 + 0], py = c[deg * 3 + 1], pz = c[deg * 3 + 2];
		for (i = deg; i > 0; i--) {
			px = px * t + c[(i - 1) * 3 + 0];
			py = py * t + c[(i - 1) * 3 + 1];
			pz = pz * t + c[(i - 1) * 3 + 2];
		}
		out_result->position.x = px; out_result->position.y = py; out_result->position.z = pz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && deg >= 1) {
		qaws_scalar dx = (qaws_scalar)deg * c[deg * 3 + 0];
		qaws_scalar dy = (qaws_scalar)deg * c[deg * 3 + 1];
		qaws_scalar dz = (qaws_scalar)deg * c[deg * 3 + 2];
		for (i = deg; i > 1; i--) {
			dx = dx * t + (qaws_scalar)(i - 1) * c[(i - 1) * 3 + 0];
			dy = dy * t + (qaws_scalar)(i - 1) * c[(i - 1) * 3 + 1];
			dz = dz * t + (qaws_scalar)(i - 1) * c[(i - 1) * 3 + 2];
		}
		out_result->d1.x = dx; out_result->d1.y = dy; out_result->d1.z = dz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && deg >= 2) {
		qaws_scalar dx = (qaws_scalar)(deg * (deg - 1)) * c[deg * 3 + 0];
		qaws_scalar dy = (qaws_scalar)(deg * (deg - 1)) * c[deg * 3 + 1];
		qaws_scalar dz = (qaws_scalar)(deg * (deg - 1)) * c[deg * 3 + 2];
		for (i = deg; i > 2; i--) {
			dx = dx * t + (qaws_scalar)((i - 1) * (i - 2)) * c[(i - 1) * 3 + 0];
			dy = dy * t + (qaws_scalar)((i - 1) * (i - 2)) * c[(i - 1) * 3 + 1];
			dz = dz * t + (qaws_scalar)((i - 1) * (i - 2)) * c[(i - 1) * 3 + 2];
		}
		out_result->d2.x = dx; out_result->d2.y = dy; out_result->d2.z = dz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	return QAWS_STATUS_OK;
}

static qaws_curve_vtable const inline_poly_vtable = {
	inline_poly_eval_2d, inline_poly_eval_3d, inline_noop_destroy,
	inline_not_closed, inline_not_periodic, inline_not_rational, inline_c3,
	NULL /* diff */
};

qaws_status qaws_curve_init_polynomial_inline(
	qaws_polynomial_desc const* desc,
	qaws_curve_inline* out)
{
	unsigned int dim_count, n_scalars;
	qaws_curve* curve;
	qaws_scalar* span_bounds;
	qaws_scalar* coeff_data;
	inline_poly_impl* impl;
	size_t offset;
	qaws_scalar t_min, t_max;

	if (!desc || !out) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (desc->degree < 1 || desc->degree > 7) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->coefficient_count != desc->degree + 1)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->coefficients) return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (unsigned int)desc->dimension;
	n_scalars = dim_count * desc->coefficient_count;

	t_min = desc->t_min;
	t_max = desc->t_max;
	if (t_min == 0 && t_max == 0) { t_min = 0; t_max = 1; }

	/* Compute layout */
	offset = align_up(sizeof(qaws_curve), sizeof(qaws_scalar));
	{
		size_t sb_off = offset;
		offset += 2 * sizeof(qaws_scalar);
		offset = align_up(offset, sizeof(qaws_scalar));
		{
			size_t coeff_off = offset;
			offset += (size_t)n_scalars * sizeof(qaws_scalar);
			offset = align_up(offset, sizeof(void*));
			{
				size_t impl_off = offset;
				offset += sizeof(inline_poly_impl);

				if (offset > QAWS_INLINE_STORAGE_BYTES)
					return QAWS_STATUS_BUFFER_TOO_SMALL;

				memset(out->storage, 0, QAWS_INLINE_STORAGE_BYTES);

				curve = (qaws_curve*)(void*)out->storage;
				span_bounds = (qaws_scalar*)(void*)&out->storage[sb_off];
				coeff_data = (qaws_scalar*)(void*)&out->storage[coeff_off];
				impl = (inline_poly_impl*)(void*)&out->storage[impl_off];

				memcpy(coeff_data, desc->coefficients,
					(size_t)n_scalars * sizeof(qaws_scalar));

				impl->coefficients = coeff_data;
				impl->coefficient_count = desc->coefficient_count;
				impl->dim_count = dim_count;

				span_bounds[0] = t_min;
				span_bounds[1] = t_max;

				curve->kind = QAWS_CURVE_KIND_POLYNOMIAL;
				curve->dimension = desc->dimension;
				curve->degree = desc->degree;
				curve->span_count = 1;
				curve->parameter_range.min_value = t_min;
				curve->parameter_range.max_value = t_max;
				curve->span_boundaries = span_bounds;
				curve->vtable = &inline_poly_vtable;
				curve->impl = impl;
				curve->allocator = NULL;
			}
		}
	}

	return QAWS_STATUS_OK;
}
