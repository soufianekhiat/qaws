#ifndef QAWS_INTERNAL_KINDS_H
#define QAWS_INTERNAL_KINDS_H

/* Impl structs of the curve kinds whose public headers carry their own types. */

#include "qaws_internal_types.h"
#include "../qaws_arc.h"
#include "../qaws_subdivision.h"

typedef struct qaws_arc_impl
{
	qaws_arc_segment* segments;
	unsigned int segment_count;
} qaws_arc_impl;

typedef struct qaws_polynomial_impl
{
	qaws_scalar* coefficients;         /* (degree+1) * dim_count scalars */
	unsigned int coefficient_count;
} qaws_polynomial_impl;

typedef struct qaws_clothoid_impl
{
	qaws_scalar origin_x, origin_y;
	qaws_scalar start_angle;
	qaws_scalar kappa_0;             /* start curvature */
	qaws_scalar kappa_1;             /* end curvature */
	qaws_scalar length;              /* total arc length */
	qaws_scalar origin[2];           /* origin, contiguous (CENTER field) */
	qaws_scalar rate;                /* (kappa_1 - kappa_0) / length (CURVATURE_RATE field) */
} qaws_clothoid_impl;

typedef struct qaws_subdivision_impl
{
	qaws_scalar* refined_points;      /* refined_count * dim_count */
	unsigned int refined_count;
	int closed;
	qaws_scalar* segment_coeffs;      /* uniform Catmull-Rom cubic coefficients */
	qaws_subdivision_scheme scheme;
} qaws_subdivision_impl;

/* Arc-length wrappers (QAWS_CURVE_KIND_REPARAMETERIZED): the source curve and
   its parameter at wrapper parameter t. Returns 0 when curve is not one. */
int qaws_internal_reparam_source(
	qaws_curve const** out_source,
	qaws_scalar* out_source_t,
	qaws_curve const* curve,
	qaws_scalar t);

/* The curve on [t0, t1] as a cubic B-spline (or Bezier) fitted span by span
   to 1e-9 of its extent (1e-5 in float builds); exact for cubic spans. */
qaws_status qaws_internal_curve_fit_cubic(qaws_curve** out_curve, qaws_curve const* curve, qaws_scalar t0, qaws_scalar t1);

/* 1 when the curve's two ends coincide, relative to the size of the curve. */
int qaws_internal_curve_ends_meet(qaws_curve const* curve);

/* A borrowed curve (QAWS_CURVE_KIND_VIEW) copied into an owned B-spline or NURBS curve, with its allocator. */
qaws_status qaws_internal_curve_view_copy(qaws_curve const* curve, qaws_curve** out_curve);

/* Span evaluators of owned families on per-span data, shared with the borrowed views (qaws_view.c), which rebuild that data from
   host memory at each evaluation: one arc segment; a trajectory span's dim * 4 coefficients and its duration; a subdivision span's
   dim * 4 uniform Catmull-Rom coefficients (from its four stencil points); the Yuksel blend of the two sub-curves of a span (first,
   last: an open curve's end spans, one sub-curve unblended) and the sub-curve at point i of a mode. */
qaws_status qaws_internal_arc_eval_2d(qaws_arc_segment const* seg, qaws_scalar local_t, unsigned int eval_flags, qaws_eval_result_2d* out_result);
qaws_status qaws_internal_arc_eval_3d(qaws_arc_segment const* seg, qaws_scalar local_t, unsigned int eval_flags, qaws_eval_result_3d* out_result);
int qaws_internal_arc_closed(qaws_dimension dimension, qaws_arc_segment const* first, qaws_arc_segment const* last);
qaws_status qaws_internal_trajectory_eval_2d(qaws_scalar const* span_coeffs, qaws_scalar dt, qaws_scalar local_t, unsigned int eval_flags,
	qaws_eval_result_2d* out_result);
qaws_status qaws_internal_trajectory_eval_3d(qaws_scalar const* span_coeffs, qaws_scalar dt, qaws_scalar local_t, unsigned int eval_flags,
	qaws_eval_result_3d* out_result);
void qaws_internal_subdivision_cr_coeffs(qaws_scalar const* P0, qaws_scalar const* P1, qaws_scalar const* P2, qaws_scalar const* P3,
	unsigned int dim_count, qaws_scalar* coeffs_out);
qaws_status qaws_internal_subdivision_eval_2d(qaws_scalar const* coeffs, qaws_scalar local_t, unsigned int eval_flags, qaws_eval_result_2d* out_result);
qaws_status qaws_internal_subdivision_eval_3d(qaws_scalar const* coeffs, qaws_scalar local_t, unsigned int eval_flags, qaws_eval_result_3d* out_result);
void qaws_internal_yuksel_subcurve(qaws_scalar const* pts, unsigned int n, unsigned int i, int closed, unsigned int dim, int mode,
	qaws_yuksel_subcurve* sc);
void qaws_internal_yuksel_blend(qaws_yuksel_subcurve const* sc1, qaws_yuksel_subcurve const* sc2, int is_first_span, int is_last_span,
	qaws_scalar t, unsigned int dim, unsigned int eval_flags, qaws_scalar* out_pos, qaws_scalar* out_d1, qaws_scalar* out_d2,
	qaws_scalar* out_d3);

#endif /* QAWS_INTERNAL_KINDS_H */
