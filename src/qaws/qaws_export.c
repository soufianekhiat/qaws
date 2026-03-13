#include "qaws_export.h"
#include "qaws_bezier.h"
#include "qaws_bspline.h"
#include "qaws_curve.h"
#include "qaws_convert.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_rational_bezier.h"
#include "qaws_sampling.h"
#include "qaws_surface.h"
#include "qaws_surface_types.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>

/* -------------------------------------------------------------------------- */
/*  Helpers                                                                    */
/* -------------------------------------------------------------------------- */

static qaws_scalar scalar_sqrt(qaws_scalar x)
{
#if QAWS_SCALAR_IS_FLOAT
	return sqrtf(x);
#else
	return sqrt(x);
#endif
}

static qaws_scalar scalar_floor(qaws_scalar x)
{
#if QAWS_SCALAR_IS_FLOAT
	return floorf(x);
#else
	return floor(x);
#endif
}

/* Solve A x = b via Gaussian elimination with partial pivoting.
   A is size x size (row-major, modified in place).
   b is size x dim_count (column-major per dimension: b[d * size + i]).
   Solution is written back into b.
   Returns QAWS_STATUS_OK on success, QAWS_STATUS_NUMERICAL_FAILURE if singular. */
static qaws_status solve_linear_system(
	qaws_scalar* A,
	qaws_scalar* b,
	unsigned int size,
	unsigned int dim_count)
{
	unsigned int i, j, k, d;

	for (i = 0; i < size; i++) {
		/* Partial pivoting: find row with largest absolute value in column i */
		unsigned int max_row = i;
		qaws_scalar max_val = A[i * size + i];
		qaws_scalar tmp;

		if (max_val < (qaws_scalar)0.0)
			max_val = -max_val;

		for (k = i + 1; k < size; k++) {
			qaws_scalar v = A[k * size + i];
			if (v < (qaws_scalar)0.0)
				v = -v;
			if (v > max_val) {
				max_val = v;
				max_row = k;
			}
		}

		if (max_val < (qaws_scalar)1e-12)
			return QAWS_STATUS_NUMERICAL_FAILURE;

		/* Swap rows i and max_row in A and b */
		if (max_row != i) {
			for (j = 0; j < size; j++) {
				tmp = A[i * size + j];
				A[i * size + j] = A[max_row * size + j];
				A[max_row * size + j] = tmp;
			}
			for (d = 0; d < dim_count; d++) {
				tmp = b[d * size + i];
				b[d * size + i] = b[d * size + max_row];
				b[d * size + max_row] = tmp;
			}
		}

		/* Eliminate below */
		for (k = i + 1; k < size; k++) {
			qaws_scalar factor = A[k * size + i] / A[i * size + i];
			A[k * size + i] = (qaws_scalar)0.0;
			for (j = i + 1; j < size; j++)
				A[k * size + j] -= factor * A[i * size + j];
			for (d = 0; d < dim_count; d++)
				b[d * size + k] -= factor * b[d * size + i];
		}
	}

	/* Back substitution */
	for (d = 0; d < dim_count; d++) {
		for (i = size; i-- > 0; ) {
			qaws_scalar sum = b[d * size + i];
			for (j = i + 1; j < size; j++)
				sum -= A[i * size + j] * b[d * size + j];
			b[d * size + i] = sum / A[i * size + i];
		}
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  B-spline least-squares fitting                                             */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_fit_bspline(
	qaws_bspline_fit_desc const* desc,
	qaws_curve** out_curve)
{
	unsigned int m;            /* data_point_count */
	unsigned int n;            /* control_point_count */
	unsigned int p;            /* degree */
	unsigned int dim_count;
	unsigned int knot_count;
	unsigned int sys_size;     /* n - 2: number of unknowns */
	qaws_scalar const* data;

	qaws_scalar* params  = NULL;  /* m parameter values */
	qaws_scalar* knots   = NULL;  /* n + p + 1 knots */
	qaws_scalar* NtN     = NULL;  /* sys_size x sys_size normal matrix */
	qaws_scalar* rhs     = NULL;  /* dim_count x sys_size right-hand side */
	qaws_scalar* basis   = NULL;  /* p + 1 basis function values */
	qaws_scalar* ctrl    = NULL;  /* n * dim_count control point scalars */
	qaws_status status = QAWS_STATUS_OK;

	unsigned int i, j, k, d;

	/* --- Validation ---------------------------------------------------- */

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	if (desc->dimension != QAWS_DIMENSION_2D &&
	    desc->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (!desc->data_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	m = desc->data_point_count;
	n = desc->control_point_count;
	p = desc->degree;
	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	data = (qaws_scalar const*)desc->data_points;

	if (p < 1)
		return QAWS_STATUS_INVALID_DEGREE;

	if (m < p + 1)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (n <= p)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (n > m)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	sys_size = n - 2;

	/* --- Step 1: Parameter assignment ---------------------------------- */

	params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)m);
	if (!params) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	if (desc->parameters) {
		memcpy(params, desc->parameters, sizeof(qaws_scalar) * (size_t)m);
	} else {
		/* Chord-length parameterization */
		qaws_scalar total_length = (qaws_scalar)0.0;
		params[0] = (qaws_scalar)0.0;
		for (i = 1; i < m; i++) {
			qaws_scalar dist_sq = (qaws_scalar)0.0;
			for (d = 0; d < dim_count; d++) {
				qaws_scalar diff = data[i * dim_count + d]
				                 - data[(i - 1) * dim_count + d];
				dist_sq += diff * diff;
			}
			total_length += scalar_sqrt(dist_sq);
			params[i] = total_length;
		}
		if (total_length > (qaws_scalar)0.0) {
			for (i = 1; i < m; i++)
				params[i] /= total_length;
		}
		params[m - 1] = (qaws_scalar)1.0; /* ensure exact endpoint */
	}

	/* --- Step 2: Knot placement (averaging) ---------------------------- */

	knot_count = n + p + 1;
	knots = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)knot_count);
	if (!knots) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	/* First p+1 knots = 0 */
	for (i = 0; i <= p; i++)
		knots[i] = (qaws_scalar)0.0;

	/* Last p+1 knots = 1 */
	for (i = 0; i <= p; i++)
		knots[knot_count - 1 - i] = (qaws_scalar)1.0;

	/* Internal knots via averaging */
	{
		qaws_scalar d_float = (qaws_scalar)(m - 1) / (qaws_scalar)(n - p);
		for (j = 1; j <= n - p - 1; j++) {
			qaws_scalar jd = (qaws_scalar)j * d_float;
			unsigned int ii = (unsigned int)scalar_floor(jd);
			qaws_scalar alpha = jd - (qaws_scalar)ii;
			/* Piegl & Tiller: knot[p+j] = (1-alpha)*t[i-1] + alpha*t[i] */
			knots[p + j] = ((qaws_scalar)1.0 - alpha) * params[ii - 1]
			             + alpha * params[ii];
		}
	}

	/* --- Step 3: Build normal equations N^T N Q = N^T R ---------------- */

	basis = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(p + 1));
	if (!basis) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	if (sys_size > 0) {
		NtN = (qaws_scalar*)calloc(
			(size_t)sys_size * (size_t)sys_size, sizeof(qaws_scalar));
		if (!NtN) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

		rhs = (qaws_scalar*)calloc(
			(size_t)dim_count * (size_t)sys_size, sizeof(qaws_scalar));
		if (!rhs) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

		/* For each interior data point k = 1..m-2, compute basis and
		   accumulate into the normal equations. */
		for (k = 1; k < m - 1; k++) {
			unsigned int span = qaws_internal_find_knot_span(
				knots, knot_count, p, n, params[k]);

			qaws_internal_bspline_basis(
				knots, knot_count, p, span, params[k], basis);

			/* Compute R_k = D_k - N_{0,p}(t_k)*D_0 - N_{n-1,p}(t_k)*D_{m-1}
			   for each dimension. */
			for (d = 0; d < dim_count; d++) {
				qaws_scalar r_kd = data[k * dim_count + d];

				/* Subtract endpoint contributions */
				for (j = 0; j <= p; j++) {
					int cp_idx = (int)span - (int)p + (int)j;
					if (cp_idx == 0) {
						r_kd -= basis[j] * data[0 * dim_count + d];
					} else if (cp_idx == (int)(n - 1)) {
						r_kd -= basis[j] * data[(m - 1) * dim_count + d];
					}
				}

				/* Accumulate into rhs: for each interior CP index i in [1..n-2] */
				for (j = 0; j <= p; j++) {
					int cp_idx = (int)span - (int)p + (int)j;
					if (cp_idx >= 1 && cp_idx <= (int)(n - 2)) {
						unsigned int ii = (unsigned int)(cp_idx - 1);
						rhs[d * sys_size + ii] += basis[j] * r_kd;
					}
				}
			}

			/* Accumulate into NtN: for each pair of interior CP indices */
			for (i = 0; i <= p; i++) {
				int ci = (int)span - (int)p + (int)i;
				if (ci < 1 || ci > (int)(n - 2))
					continue;
				for (j = 0; j <= p; j++) {
					int cj = (int)span - (int)p + (int)j;
					if (cj < 1 || cj > (int)(n - 2))
						continue;
					NtN[(unsigned int)(ci - 1) * sys_size + (unsigned int)(cj - 1)]
						+= basis[i] * basis[j];
				}
			}
		}

		/* --- Step 4: Solve ------------------------------------------------- */

		status = solve_linear_system(NtN, rhs, sys_size, dim_count);
		if (status != QAWS_STATUS_OK)
			goto cleanup;
	}

	/* --- Step 5: Assemble control points ------------------------------- */

	ctrl = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(n * dim_count));
	if (!ctrl) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	/* Q_0 = D_0 (first data point) */
	for (d = 0; d < dim_count; d++)
		ctrl[0 * dim_count + d] = data[0 * dim_count + d];

	/* Q_{n-1} = D_{m-1} (last data point) */
	for (d = 0; d < dim_count; d++)
		ctrl[(n - 1) * dim_count + d] = data[(m - 1) * dim_count + d];

	/* Interior control points from solved system */
	for (i = 0; i < sys_size; i++) {
		for (d = 0; d < dim_count; d++)
			ctrl[(i + 1) * dim_count + d] = rhs[d * sys_size + i];
	}

	/* --- Create B-spline curve ----------------------------------------- */
	{
		qaws_bspline_desc bspline_desc;
		memset(&bspline_desc, 0, sizeof(bspline_desc));
		bspline_desc.dimension = desc->dimension;
		bspline_desc.degree = p;
		bspline_desc.control_points = ctrl;
		bspline_desc.control_point_count = n;
		bspline_desc.knots = knots;
		bspline_desc.knot_count = knot_count;
		bspline_desc.is_uniform = 0;
		bspline_desc.is_closed = 0;

		status = qaws_curve_create_bspline(&bspline_desc, out_curve);
	}

cleanup:
	free(params);
	free(knots);
	free(NtN);
	free(rhs);
	free(basis);
	free(ctrl);
	return status;
}

/* ========================================================================== */
/*  SVG path export                                                            */
/* ========================================================================== */

/* Safe snprintf into a bounded buffer, advancing the write position. */
typedef struct svg_buf
{
	char* data;
	unsigned int capacity;
	unsigned int pos;
} svg_buf;

static void svgbuf_append(svg_buf* b, char const* fmt, ...)
{
	va_list ap;
	int n;
	unsigned int avail;
	if (b->pos + 1 >= b->capacity) return; /* need at least 1 byte + null */
	avail = b->capacity - b->pos;
	va_start(ap, fmt);
	n = vsnprintf(b->data + b->pos, (size_t)avail, fmt, ap);
	va_end(ap);
	if (n < 0)
	{
		/* MSVC: truncation returns -1. Mark buffer as full. */
		b->pos = b->capacity;
	}
	else if ((unsigned int)n >= avail)
	{
		/* C99: returns would-be length. Buffer was truncated. */
		b->pos = b->capacity;
	}
	else
	{
		b->pos += (unsigned int)n;
	}
}

/* Write M (move-to) command */
static void svg_emit_move(svg_buf* b, qaws_scalar x, qaws_scalar y)
{
	svgbuf_append(b, "M%.6g %.6g", (double)x, (double)y);
}

/* Write C (cubic bezier) command */
static void svg_emit_cubic(svg_buf* b,
	qaws_scalar x1, qaws_scalar y1,
	qaws_scalar x2, qaws_scalar y2,
	qaws_scalar x3, qaws_scalar y3)
{
	svgbuf_append(b, " C%.6g %.6g %.6g %.6g %.6g %.6g",
		(double)x1, (double)y1, (double)x2, (double)y2, (double)x3, (double)y3);
}

/* Write Q (quadratic bezier) command */
static void svg_emit_quadratic(svg_buf* b,
	qaws_scalar x1, qaws_scalar y1,
	qaws_scalar x2, qaws_scalar y2)
{
	svgbuf_append(b, " Q%.6g %.6g %.6g %.6g",
		(double)x1, (double)y1, (double)x2, (double)y2);
}

/* Write L (line-to) command */
static void svg_emit_line(svg_buf* b, qaws_scalar x, qaws_scalar y)
{
	svgbuf_append(b, " L%.6g %.6g", (double)x, (double)y);
}

/* ---- Export strategies per family ---- */

/* Bezier: direct mapping. Degree 1=L, 2=Q, 3=C, higher=approx */
static qaws_status svg_export_bezier(qaws_curve const* curve, svg_buf* b)
{
	unsigned int deg = qaws_curve_get_degree(curve);
	qaws_vec2 cp[64];
	unsigned int n = 0;
	qaws_status s;

	s = qaws_bezier_get_control_points(curve, cp, 64, &n);
	if (s != QAWS_STATUS_OK) return s;

	svg_emit_move(b, cp[0].x, cp[0].y);

	if (deg == 1)
	{
		svg_emit_line(b, cp[1].x, cp[1].y);
	}
	else if (deg == 2)
	{
		svg_emit_quadratic(b, cp[1].x, cp[1].y, cp[2].x, cp[2].y);
	}
	else if (deg == 3)
	{
		svg_emit_cubic(b, cp[1].x, cp[1].y, cp[2].x, cp[2].y, cp[3].x, cp[3].y);
	}
	else
	{
		/* Higher degree: sample and fit cubic segments */
		unsigned int npts = 4 * deg;
		unsigned int i;
		qaws_eval_result_2d r;
		qaws_vec2 prev;
		prev = cp[0];
		for (i = 1; i <= npts; i++)
		{
			qaws_scalar t = (qaws_scalar)i / (qaws_scalar)npts;
			qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
			svg_emit_line(b, r.position.x, r.position.y);
		}
	}
	return QAWS_STATUS_OK;
}

/* Hermite: convert each span to cubic Bezier */
static qaws_status svg_export_hermite(qaws_curve const* curve, svg_buf* b)
{
	unsigned int span_count = qaws_curve_get_span_count(curve);
	unsigned int span;
	int first = 1;

	for (span = 0; span < span_count; span++)
	{
		qaws_curve* bez = NULL;
		qaws_vec2 cp[4];
		unsigned int n = 0;
		qaws_status s;

		s = qaws_curve_convert_hermite_to_bezier(curve, span, &bez);
		if (s != QAWS_STATUS_OK) return s;

		qaws_bezier_get_control_points(bez, cp, 4, &n);

		if (first)
		{
			svg_emit_move(b, cp[0].x, cp[0].y);
			first = 0;
		}
		svg_emit_cubic(b, cp[1].x, cp[1].y, cp[2].x, cp[2].y, cp[3].x, cp[3].y);
		qaws_curve_destroy(bez);
	}
	return QAWS_STATUS_OK;
}

/* Catmull-Rom: convert each span to cubic Bezier */
static qaws_status svg_export_catmull_rom(qaws_curve const* curve, svg_buf* b)
{
	unsigned int span_count = qaws_curve_get_span_count(curve);
	unsigned int span;
	int first = 1;

	for (span = 0; span < span_count; span++)
	{
		qaws_curve* bez = NULL;
		qaws_vec2 cp[4];
		unsigned int n = 0;
		qaws_status s;

		s = qaws_curve_convert_catmull_rom_to_bezier(curve, span, &bez);
		if (s != QAWS_STATUS_OK) return s;

		qaws_bezier_get_control_points(bez, cp, 4, &n);

		if (first)
		{
			svg_emit_move(b, cp[0].x, cp[0].y);
			first = 0;
		}
		svg_emit_cubic(b, cp[1].x, cp[1].y, cp[2].x, cp[2].y, cp[3].x, cp[3].y);
		qaws_curve_destroy(bez);
	}
	return QAWS_STATUS_OK;
}

/* Rational Bezier: quadratic = exact Q command, otherwise approximate */
static qaws_status svg_export_rational_bezier(qaws_curve const* curve, svg_buf* b,
	unsigned int sample_count)
{
	unsigned int deg = qaws_curve_get_degree(curve);

	if (deg == 2)
	{
		/* Rational quadratic maps to SVG Q only if w1=1 (non-rational case).
		   For general rational, we must approximate. */
	}

	/* Approximate via sampling */
	{
		qaws_range range = qaws_curve_get_parameter_range(curve);
		qaws_scalar len = range.max_value - range.min_value;
		unsigned int i;
		qaws_eval_result_2d r;

		qaws_curve_evaluate_2d(curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &r);
		svg_emit_move(b, r.position.x, r.position.y);

		for (i = 1; i <= sample_count; i++)
		{
			qaws_scalar t = range.min_value + len * (qaws_scalar)i / (qaws_scalar)sample_count;
			qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
			svg_emit_line(b, r.position.x, r.position.y);
		}
	}
	return QAWS_STATUS_OK;
}

/*
 * Generic fallback: sample curve, then fit piecewise cubic Bezier segments.
 * Uses cubic Hermite interpolation: at each segment boundary, evaluate
 * position and tangent, then derive Bezier CPs from the Hermite form.
 */
static qaws_status svg_export_sampled(qaws_curve const* curve, svg_buf* b,
	unsigned int sample_count)
{
	qaws_range range = qaws_curve_get_parameter_range(curve);
	qaws_scalar len = range.max_value - range.min_value;
	unsigned int seg_count;
	unsigned int i;
	qaws_eval_result_2d prev, curr;
	qaws_status s;

	if (sample_count < 4) sample_count = 4;
	/* We'll create segments: every 3 sample intervals = 1 cubic segment.
	   Actually, the simplest high-quality approach: evaluate position+tangent
	   at evenly spaced points, then each pair of consecutive points defines
	   a cubic Hermite segment that we convert to Bezier. */
	seg_count = sample_count;

	s = qaws_curve_evaluate_2d(curve, range.min_value,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &prev);
	if (s != QAWS_STATUS_OK) return s;

	svg_emit_move(b, prev.position.x, prev.position.y);

	for (i = 1; i <= seg_count; i++)
	{
		qaws_scalar t = range.min_value + len * (qaws_scalar)i / (qaws_scalar)seg_count;
		qaws_scalar dt = len / (qaws_scalar)seg_count;
		qaws_vec2 b1, b2;

		s = qaws_curve_evaluate_2d(curve, t,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &curr);
		if (s != QAWS_STATUS_OK) return s;

		/* Hermite-to-Bezier: B1 = P0 + M0/(3), B2 = P1 - M1/(3)
		   where M0 = tangent * dt (scale to segment interval) */
		b1.x = prev.position.x + prev.d1.x * dt / (qaws_scalar)3;
		b1.y = prev.position.y + prev.d1.y * dt / (qaws_scalar)3;
		b2.x = curr.position.x - curr.d1.x * dt / (qaws_scalar)3;
		b2.y = curr.position.y - curr.d1.y * dt / (qaws_scalar)3;

		svg_emit_cubic(b, b1.x, b1.y, b2.x, b2.y,
			curr.position.x, curr.position.y);

		prev = curr;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_export_svg_path(
	qaws_curve const* curve,
	char* out_path_data,
	unsigned int capacity,
	unsigned int* out_length,
	unsigned int sample_count)
{
	svg_buf b;
	qaws_curve_kind kind;
	qaws_status s;

	if (!curve || !out_path_data || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(curve) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (capacity == 0)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	if (sample_count == 0)
		sample_count = 64;

	b.data = out_path_data;
	b.capacity = capacity;
	b.pos = 0;
	b.data[0] = '\0';

	kind = qaws_curve_get_kind(curve);

	switch (kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
		s = svg_export_bezier(curve, &b);
		break;
	case QAWS_CURVE_KIND_HERMITE:
		s = svg_export_hermite(curve, &b);
		break;
	case QAWS_CURVE_KIND_CATMULL_ROM:
		s = svg_export_catmull_rom(curve, &b);
		break;
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
		s = svg_export_rational_bezier(curve, &b, sample_count);
		break;
	default:
		/* All other families: sample + Hermite-to-Bezier fit */
		s = svg_export_sampled(curve, &b, sample_count);
		break;
	}

	/* Null-terminate */
	if (b.pos < b.capacity)
		b.data[b.pos] = '\0';
	else
		b.data[b.capacity - 1] = '\0';

	*out_length = b.pos;
	return s;
}

/* -------------------------------------------------------------------------- */
/*  Polyline export                                                           */
/* -------------------------------------------------------------------------- */

qaws_status qaws_polyline_export_2d(
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_polyline_sampling sampling,
	qaws_vec2* out_points,
	unsigned int* out_count)
{
	if (!curve || !out_points || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(curve) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (sampling == QAWS_POLYLINE_CURVATURE)
	{
		return qaws_curve_sample_curvature_weighted_2d(
			curve, sample_count, out_points, sample_count, out_count);
	}
	else
	{
		qaws_range range;
		qaws_scalar t, dt;
		unsigned int i;

		range = qaws_curve_get_parameter_range(curve);
		dt = (range.max_value - range.min_value) / (qaws_scalar)(sample_count - 1);

		for (i = 0; i < sample_count; i++)
		{
			qaws_eval_result_2d r;
			qaws_status s;

			t = (i == sample_count - 1) ? range.max_value : range.min_value + dt * (qaws_scalar)i;
			s = qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
			if (s != QAWS_STATUS_OK)
			{
				*out_count = i;
				return s;
			}
			out_points[i] = r.position;
		}

		*out_count = sample_count;
		return QAWS_STATUS_OK;
	}
}

qaws_status qaws_polyline_export_3d(
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_polyline_sampling sampling,
	qaws_vec3* out_points,
	unsigned int* out_count)
{
	if (!curve || !out_points || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (sample_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(curve) != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (sampling == QAWS_POLYLINE_CURVATURE)
	{
		return qaws_curve_sample_curvature_weighted_3d(
			curve, sample_count, out_points, sample_count, out_count);
	}
	else
	{
		qaws_range range;
		qaws_scalar t, dt;
		unsigned int i;

		range = qaws_curve_get_parameter_range(curve);
		dt = (range.max_value - range.min_value) / (qaws_scalar)(sample_count - 1);

		for (i = 0; i < sample_count; i++)
		{
			qaws_eval_result_3d r;
			qaws_status s;

			t = (i == sample_count - 1) ? range.max_value : range.min_value + dt * (qaws_scalar)i;
			s = qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
			if (s != QAWS_STATUS_OK)
			{
				*out_count = i;
				return s;
			}
			out_points[i] = r.position;
		}

		*out_count = sample_count;
		return QAWS_STATUS_OK;
	}
}

/* ========================================================================== */
/*  OBJ surface export                                                         */
/* ========================================================================== */

/* Safe snprintf into a bounded buffer, advancing the write position. */
typedef struct obj_buf
{
	char* data;
	unsigned int capacity;
	unsigned int pos;
} obj_buf;

static void objbuf_append(obj_buf* b, char const* fmt, ...)
{
	va_list ap;
	int n;
	unsigned int avail;
	if (b->pos + 1 >= b->capacity) return; /* need at least 1 byte + null */
	avail = b->capacity - b->pos;
	va_start(ap, fmt);
	n = vsnprintf(b->data + b->pos, (size_t)avail, fmt, ap);
	va_end(ap);
	if (n < 0)
	{
		/* MSVC: truncation returns -1. Mark buffer as full. */
		b->pos = b->capacity;
	}
	else if ((unsigned int)n >= avail)
	{
		/* C99: returns would-be length. Buffer was truncated. */
		b->pos = b->capacity;
	}
	else
	{
		b->pos += (unsigned int)n;
	}
}

qaws_status qaws_surface_export_obj(
	qaws_surface const* surface,
	unsigned int u_samples,
	unsigned int v_samples,
	int include_normals,
	char* out_obj_data,
	unsigned int capacity,
	unsigned int* out_length)
{
	obj_buf b;
	qaws_range ur, vr;
	qaws_scalar u_min, u_len, v_min, v_len;
	unsigned int ui, vi;
	unsigned int eval_flags;
	int flip = 0;

	/* --- Validation ---------------------------------------------------- */

	if (!surface || !out_obj_data || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (u_samples < 2 || v_samples < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (capacity == 0)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	b.data = out_obj_data;
	b.capacity = capacity;
	b.pos = 0;
	b.data[0] = '\0';

	ur = qaws_surface_get_u_range(surface);
	vr = qaws_surface_get_v_range(surface);
	u_min = ur.min_value;
	u_len = ur.max_value - ur.min_value;
	v_min = vr.min_value;
	v_len = vr.max_value - vr.min_value;

	eval_flags = QAWS_SURFACE_EVAL_POSITION;
	if (include_normals)
		eval_flags |= QAWS_SURFACE_EVAL_NORMAL;

	/* --- Auto-detect normal orientation -------------------------------- */
	if (include_normals)
	{
		qaws_surface_eval_result rc, r00, r10, r01, r11, rv0, rv1;
		qaws_vec3 centroid;
		qaws_scalar dx, dy, dz, dot_val, gap, diag;

		memset(&rc, 0, sizeof(rc));
		memset(&r00, 0, sizeof(r00));
		memset(&r10, 0, sizeof(r10));
		memset(&r01, 0, sizeof(r01));
		memset(&r11, 0, sizeof(r11));
		memset(&rv0, 0, sizeof(rv0));
		memset(&rv1, 0, sizeof(rv1));

		qaws_surface_evaluate(surface,
			u_min + u_len * (qaws_scalar)0.5,
			v_min + v_len * (qaws_scalar)0.5,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &rc);
		qaws_surface_evaluate(surface, u_min, v_min,
			QAWS_SURFACE_EVAL_POSITION, &r00);
		qaws_surface_evaluate(surface, u_min + u_len, v_min,
			QAWS_SURFACE_EVAL_POSITION, &r10);
		qaws_surface_evaluate(surface, u_min, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &r01);
		qaws_surface_evaluate(surface, u_min + u_len, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &r11);

		/* V-closure: pipe/swept surfaces wrap in v */
		qaws_surface_evaluate(surface,
			u_min + u_len * (qaws_scalar)0.5, v_min,
			QAWS_SURFACE_EVAL_POSITION, &rv0);
		qaws_surface_evaluate(surface,
			u_min + u_len * (qaws_scalar)0.5, v_min + v_len,
			QAWS_SURFACE_EVAL_POSITION, &rv1);

		dx = rv1.position.x - rv0.position.x;
		dy = rv1.position.y - rv0.position.y;
		dz = rv1.position.z - rv0.position.z;
		gap = dx * dx + dy * dy + dz * dz;

		dx = r11.position.x - r00.position.x;
		dy = r11.position.y - r00.position.y;
		dz = r11.position.z - r00.position.z;
		diag = dx * dx + dy * dy + dz * dz;

		if (gap < diag * (qaws_scalar)0.0001 || gap < (qaws_scalar)1e-10)
		{
			flip = 1;
		}

		/* Centroid heuristic for open surfaces */
		if (!flip)
		{
			centroid.x = (r00.position.x + r10.position.x + r01.position.x + r11.position.x) * (qaws_scalar)0.25;
			centroid.y = (r00.position.y + r10.position.y + r01.position.y + r11.position.y) * (qaws_scalar)0.25;
			centroid.z = (r00.position.z + r10.position.z + r01.position.z + r11.position.z) * (qaws_scalar)0.25;

			dx = centroid.x - rc.position.x;
			dy = centroid.y - rc.position.y;
			dz = centroid.z - rc.position.z;

			dot_val = rc.normal.x * dx + rc.normal.y * dy + rc.normal.z * dz;
			if (dot_val > (qaws_scalar)0.01)
				flip = 1;
		}
	}

	/* --- Emit vertices ------------------------------------------------- */

	for (ui = 0; ui < u_samples; ui++)
	{
		qaws_scalar u = u_min + u_len * (qaws_scalar)ui / (qaws_scalar)(u_samples - 1);
		for (vi = 0; vi < v_samples; vi++)
		{
			qaws_scalar v = v_min + v_len * (qaws_scalar)vi / (qaws_scalar)(v_samples - 1);
			qaws_surface_eval_result r;
			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surface, u, v, eval_flags, &r);
			objbuf_append(&b, "v %.6g %.6g %.6g\n",
				(double)r.position.x, (double)r.position.y, (double)r.position.z);
		}
	}

	/* --- Emit normals -------------------------------------------------- */

	if (include_normals)
	{
		qaws_scalar sign = flip ? (qaws_scalar)-1 : (qaws_scalar)1;
		for (ui = 0; ui < u_samples; ui++)
		{
			qaws_scalar u = u_min + u_len * (qaws_scalar)ui / (qaws_scalar)(u_samples - 1);
			for (vi = 0; vi < v_samples; vi++)
			{
				qaws_scalar v = v_min + v_len * (qaws_scalar)vi / (qaws_scalar)(v_samples - 1);
				qaws_surface_eval_result r;
				memset(&r, 0, sizeof(r));
				qaws_surface_evaluate(surface, u, v,
					QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
				objbuf_append(&b, "vn %.6g %.6g %.6g\n",
					(double)(r.normal.x * sign),
					(double)(r.normal.y * sign),
					(double)(r.normal.z * sign));
			}
		}
	}

	/* --- Emit triangle faces ------------------------------------------- */

	for (ui = 0; ui < u_samples - 1; ui++)
	{
		for (vi = 0; vi < v_samples - 1; vi++)
		{
			/* OBJ is 1-indexed */
			unsigned int v00 = 1 + ui * v_samples + vi;
			unsigned int v01 = v00 + 1;
			unsigned int v10 = v00 + v_samples;
			unsigned int v11 = v10 + 1;

			if (include_normals)
			{
				if (flip)
				{
					objbuf_append(&b, "f %u//%u %u//%u %u//%u\n",
						v00, v00, v01, v01, v10, v10);
					objbuf_append(&b, "f %u//%u %u//%u %u//%u\n",
						v01, v01, v11, v11, v10, v10);
				}
				else
				{
					objbuf_append(&b, "f %u//%u %u//%u %u//%u\n",
						v00, v00, v10, v10, v01, v01);
					objbuf_append(&b, "f %u//%u %u//%u %u//%u\n",
						v10, v10, v11, v11, v01, v01);
				}
			}
			else
			{
				if (flip)
				{
					objbuf_append(&b, "f %u %u %u\n", v00, v01, v10);
					objbuf_append(&b, "f %u %u %u\n", v01, v11, v10);
				}
				else
				{
					objbuf_append(&b, "f %u %u %u\n", v00, v10, v01);
					objbuf_append(&b, "f %u %u %u\n", v10, v11, v01);
				}
			}
		}
	}

	/* Null-terminate */
	if (b.pos < b.capacity)
		b.data[b.pos] = '\0';
	else
		b.data[b.capacity - 1] = '\0';

	*out_length = b.pos;
	return (b.pos >= b.capacity) ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}
