#include "qaws_composite.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Impl struct                                                        */
/* ------------------------------------------------------------------ */

typedef struct qaws_composite_impl
{
	qaws_curve** segments;           /* owned segment curves */
	unsigned int segment_count;
} qaws_composite_impl;

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

/* Map a local_t in [0,1] (relative to a composite span) into the
   sub-curve's actual global parameter, then evaluate through the
   public API.  Also returns the derivative scale factor (b-a). */
static void composite_map_parameter(
	qaws_curve const* segment,
	qaws_scalar local_t,
	qaws_scalar* out_segment_param,
	qaws_scalar* out_deriv_scale)
{
	qaws_scalar a = segment->parameter_range.min_value;
	qaws_scalar b = segment->parameter_range.max_value;
	*out_segment_param = a + local_t * (b - a);
	*out_deriv_scale = b - a;
}

/* ------------------------------------------------------------------ */
/*  Evaluation                                                         */
/* ------------------------------------------------------------------ */

static qaws_status composite_eval_span_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	qaws_curve const* seg;
	qaws_scalar seg_param, scale;
	qaws_status status;

	if (span_index >= impl->segment_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	seg = impl->segments[span_index];
	composite_map_parameter(seg, local_t, &seg_param, &scale);

	status = qaws_curve_evaluate_2d(seg, seg_param, eval_flags, out_result);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Chain rule: d/dt_composite = d/dt_segment * dt_segment/dt_composite.
	   Each composite span has width 1 and the segment range has width
	   (b-a), so dt_segment/dt_local = scale = (b-a). */
	if (out_result->valid_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x *= scale;
		out_result->d1.y *= scale;
	}
	if (out_result->valid_flags & QAWS_EVAL_FLAG_D2) {
		qaws_scalar scale2 = scale * scale;
		out_result->d2.x *= scale2;
		out_result->d2.y *= scale2;
	}
	if (out_result->valid_flags & QAWS_EVAL_FLAG_D3) {
		qaws_scalar scale3 = scale * scale * scale;
		out_result->d3.x *= scale3;
		out_result->d3.y *= scale3;
	}

	return QAWS_STATUS_OK;
}

static qaws_status composite_eval_span_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	qaws_curve const* seg;
	qaws_scalar seg_param, scale;
	qaws_status status;

	if (span_index >= impl->segment_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	seg = impl->segments[span_index];
	composite_map_parameter(seg, local_t, &seg_param, &scale);

	status = qaws_curve_evaluate_3d(seg, seg_param, eval_flags, out_result);
	if (status != QAWS_STATUS_OK)
		return status;

	if (out_result->valid_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x *= scale;
		out_result->d1.y *= scale;
		out_result->d1.z *= scale;
	}
	if (out_result->valid_flags & QAWS_EVAL_FLAG_D2) {
		qaws_scalar scale2 = scale * scale;
		out_result->d2.x *= scale2;
		out_result->d2.y *= scale2;
		out_result->d2.z *= scale2;
	}
	if (out_result->valid_flags & QAWS_EVAL_FLAG_D3) {
		qaws_scalar scale3 = scale * scale * scale;
		out_result->d3.x *= scale3;
		out_result->d3.y *= scale3;
		out_result->d3.z *= scale3;
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Destroy                                                            */
/* ------------------------------------------------------------------ */

static void composite_destroy_impl(void* impl_ptr, qaws_allocator const* allocator)
{
	qaws_composite_impl* impl = (qaws_composite_impl*)impl_ptr;
	unsigned int i;
	if (!impl) return;

	if (impl->segments) {
		for (i = 0; i < impl->segment_count; i++)
			qaws_curve_destroy(impl->segments[i]);
		qaws_internal_dealloc(allocator, impl->segments);
	}
	qaws_internal_dealloc(allocator, impl);
}

/* ------------------------------------------------------------------ */
/*  Property queries                                                   */
/* ------------------------------------------------------------------ */

static int composite_is_closed(qaws_curve const* curve)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	qaws_curve const* first;
	qaws_curve const* last;
	qaws_eval_result_2d r0, r1;
	qaws_eval_result_3d r0_3d, r1_3d;
	qaws_scalar dist_sq;

	if (impl->segment_count == 0) return 0;

	first = impl->segments[0];
	last  = impl->segments[impl->segment_count - 1];

	if (curve->dimension == QAWS_DIMENSION_2D) {
		if (qaws_curve_evaluate_2d(first, first->parameter_range.min_value,
				QAWS_EVAL_FLAG_POSITION, &r0) != QAWS_STATUS_OK)
			return 0;
		if (qaws_curve_evaluate_2d(last, last->parameter_range.max_value,
				QAWS_EVAL_FLAG_POSITION, &r1) != QAWS_STATUS_OK)
			return 0;
		dist_sq = (r1.position.x - r0.position.x) * (r1.position.x - r0.position.x)
		        + (r1.position.y - r0.position.y) * (r1.position.y - r0.position.y);
	} else {
		if (qaws_curve_evaluate_3d(first, first->parameter_range.min_value,
				QAWS_EVAL_FLAG_POSITION, &r0_3d) != QAWS_STATUS_OK)
			return 0;
		if (qaws_curve_evaluate_3d(last, last->parameter_range.max_value,
				QAWS_EVAL_FLAG_POSITION, &r1_3d) != QAWS_STATUS_OK)
			return 0;
		dist_sq = (r1_3d.position.x - r0_3d.position.x) * (r1_3d.position.x - r0_3d.position.x)
		        + (r1_3d.position.y - r0_3d.position.y) * (r1_3d.position.y - r0_3d.position.y)
		        + (r1_3d.position.z - r0_3d.position.z) * (r1_3d.position.z - r0_3d.position.z);
	}

	return dist_sq < (qaws_scalar)1e-10 ? 1 : 0;
}

static int composite_is_periodic(qaws_curve const* curve)
{
	(void)curve;
	return 0;
}

static int composite_is_rational(qaws_curve const* curve)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	unsigned int i;
	for (i = 0; i < impl->segment_count; i++) {
		if (impl->segments[i]->vtable && impl->segments[i]->vtable->is_rational
			&& impl->segments[i]->vtable->is_rational(impl->segments[i]))
			return 1;
	}
	return 0;
}

static qaws_continuity composite_get_continuity(qaws_curve const* curve)
{
	(void)curve;
	/* Conservative: we can only guarantee C0 at segment junctions without
	   analysing derivative continuity between adjacent segments. */
	return QAWS_CONTINUITY_C0;
}

/* ------------------------------------------------------------------ */
/*  Vtable                                                             */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/*  Differential rules                                                 */
/*                                                                     */
/*  Span i is segment i (child i) on a unit-width span: global         */
/*  derivative k equals the segment derivative times (b - a)^k.        */
/* ------------------------------------------------------------------ */

#define COMPOSITE_DIFF_CAPS (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)

static unsigned int composite_describe_fields(qaws_curve const* curve, qaws_field_desc* out, unsigned int capacity)
{
	(void)curve; (void)out; (void)capacity;
	return 0;
}

static qaws_status composite_primal_field(qaws_curve const* curve, qaws_diff_field field,
	qaws_scalar const** out_data, unsigned int* out_count, unsigned int* out_components)
{
	(void)curve; (void)field; (void)out_data; (void)out_count; (void)out_components;
	return QAWS_STATUS_INVALID_ARGUMENT;
}

static unsigned int composite_children(qaws_curve const* curve, qaws_diff_child* out, unsigned int capacity)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	unsigned int i;
	for (i = 0; i < impl->segment_count && i < capacity; i++)
	{
		out[i].curve = impl->segments[i];
		out[i].surface = NULL;
	}
	return impl->segment_count;
}

static void composite_scale_jet(qaws_curve_jet_3d* j, qaws_scalar scale)
{
	qaws_scalar f = scale;
	unsigned int k;
	for (k = 1; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		j->d[k].x *= f;
		j->d[k].y *= f;
		j->d[k].z *= f;
		f *= scale;
	}
}

static qaws_status composite_tangent_span(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	unsigned int span_index, qaws_scalar local_t, qaws_scalar t_dot,
	unsigned int channels, qaws_diff_views const* views,
	qaws_curve_jet_3d* primal, qaws_curve_jet_3d* tangent, qaws_curve_jet_3d* tangent2)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	qaws_curve const* seg;
	qaws_scalar seg_t, scale;
	qaws_status st;

	if (span_index >= impl->segment_count)
		return QAWS_STATUS_OUT_OF_RANGE;
	seg = impl->segments[span_index];
	composite_map_parameter(seg, local_t, &seg_t, &scale);
	/* Each composite span has unit width, so dt_segment = scale * dt. */
	st = qaws_internal_curve_tangent_any(ctx, seg, seg_t, scale * t_dot, channels,
		qaws_internal_child_views(views, span_index), primal, tangent, tangent2);
	if (st != QAWS_STATUS_OK)
		return st;
	composite_scale_jet(primal, scale);
	composite_scale_jet(tangent, scale);
	if (tangent2)
		composite_scale_jet(tangent2, scale);
	return QAWS_STATUS_OK;
}

static qaws_status composite_adjoint_span(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	unsigned int span_index, qaws_scalar local_t, unsigned int channels,
	qaws_curve_jet_3d const* jet_adjoint, qaws_diff_views* views, qaws_scalar* t_adjoint)
{
	qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
	qaws_curve const* seg;
	qaws_scalar seg_t, scale, seg_bar = 0;
	qaws_curve_jet_3d j = *jet_adjoint;
	qaws_status st;

	if (span_index >= impl->segment_count)
		return QAWS_STATUS_OUT_OF_RANGE;
	seg = impl->segments[span_index];
	composite_map_parameter(seg, local_t, &seg_t, &scale);
	composite_scale_jet(&j, scale);
	st = qaws_internal_curve_adjoint_any(ctx, seg, seg_t, channels, &j,
		(qaws_diff_views*)qaws_internal_child_views(views, span_index), t_adjoint ? &seg_bar : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (t_adjoint)
		*t_adjoint += scale * seg_bar;
	return QAWS_STATUS_OK;
}

static qaws_curve_diff_vtable const composite_diff_vtable = {
	COMPOSITE_DIFF_CAPS,
	QAWS_DIFF_PIECEWISE_SMOOTH,
	composite_describe_fields,
	composite_primal_field,
	NULL,
	composite_tangent_span,
	composite_adjoint_span,
	composite_children
};

static qaws_curve_vtable const composite_vtable = {
	composite_eval_span_2d,
	composite_eval_span_3d,
	composite_destroy_impl,
	composite_is_closed,
	composite_is_periodic,
	composite_is_rational,
	composite_get_continuity,
	&composite_diff_vtable
};

/* ------------------------------------------------------------------ */
/*  Create                                                             */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_create_composite(
	qaws_composite_desc const* desc,
	qaws_curve** out_curve)
{
	unsigned int i;
	unsigned int max_degree;
	qaws_range range;
	qaws_curve* curve;
	qaws_composite_impl* impl;

	if (!desc || !out_curve) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->segments)    return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->segment_count < 1) return QAWS_STATUS_INVALID_ARGUMENT;

	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;

	/* Verify all segments share the same dimension */
	for (i = 0; i < desc->segment_count; i++) {
		if (!desc->segments[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (desc->segments[i]->dimension != desc->dimension)
			return QAWS_STATUS_INVALID_DIMENSION;
	}

	/* Compute the maximum degree across all segments */
	max_degree = 0;
	for (i = 0; i < desc->segment_count; i++) {
		if (desc->segments[i]->degree > max_degree)
			max_degree = desc->segments[i]->degree;
	}

	/* Parameter range: [0, segment_count] */
	range.min_value = (qaws_scalar)0;
	range.max_value = (qaws_scalar)desc->segment_count;

	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_COMPOSITE, desc->dimension,
		max_degree, desc->segment_count, range, &composite_vtable);
	if (!curve) return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Span boundaries: segment i occupies [i, i+1] */
	for (i = 0; i <= desc->segment_count; i++)
		curve->span_boundaries[i] = (qaws_scalar)i;

	/* Allocate impl */
	impl = (qaws_composite_impl*)malloc(sizeof(qaws_composite_impl));
	if (!impl) {
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->segment_count = desc->segment_count;
	impl->segments = (qaws_curve**)malloc(
		sizeof(qaws_curve*) * (size_t)desc->segment_count);
	if (!impl->segments) {
		free(impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Take ownership: copy the pointer array */
	memcpy(impl->segments, desc->segments,
		sizeof(qaws_curve*) * (size_t)desc->segment_count);

	curve->impl = impl;
	*out_curve = curve;
	return QAWS_STATUS_OK;
}
