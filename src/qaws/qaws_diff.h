#ifndef QAWS_DIFF_H
#define QAWS_DIFF_H

#include "qaws_diff_types.h"
#include "qaws_surface_types.h"

/*
 * Differentiable evaluation of curves and surfaces.
 *
 * Vocabulary:
 *   eval_tangent  : primal jet + tangent jet for a parameter/coordinate tangent
 *   eval_tangent2 : additionally the second directional derivative
 *   eval_adjoint  : accumulate parameter/coordinate adjoints (+=) from jet adjoints
 *
 * Batch entry points are the primary API (many samples, one shared object).
 * Single-sample functions are thin wrappers.
 *
 * ctx may be NULL (first order, scatter accumulation, malloc scratch,
 * no report). All functions are thread-safe on immutable objects; adjoint
 * functions write into caller storage and need external synchronization
 * when several threads share one adjoint view.
 */

/* ===================================================================
 * Context, report, views
 * =================================================================== */

void qaws_diff_context_init(qaws_diff_context* ctx);
void qaws_diff_report_reset(qaws_diff_report* report);

/* View over caller storage: element i, component c at data[i * components + c]. */
qaws_field_view qaws_field_view_make(
	qaws_diff_field field,
	qaws_scalar* data,
	unsigned int count,
	unsigned int components);

/* Returns the view of a field, or NULL when the field is absent (inactive). */
qaws_field_view* qaws_diff_views_find(
	qaws_diff_views const* views,
	qaws_diff_field field);

/* Zero every active element of every view (children included). */
void qaws_diff_views_clear(qaws_diff_views* views);

/* ===================================================================
 * Parameter keys
 * =================================================================== */

qaws_param_key qaws_param_key_make(
	qaws_diff_field field,
	unsigned int element,
	unsigned int component,
	unsigned int components);

/* Prefix the key with a child index (used by derived objects). */
qaws_status qaws_param_key_prepend_child(
	qaws_param_key* key,
	unsigned int child_index);

int qaws_param_key_compare(
	qaws_param_key const* a,
	qaws_param_key const* b);

/* "child[0]/control_points/4/z". Returns the length that the full string
   needs (excluding the terminator); writes at most capacity bytes. */
unsigned int qaws_param_key_to_string(
	qaws_param_key const* key,
	char* buffer,
	unsigned int capacity);

qaws_status qaws_param_key_parse(
	char const* text,
	qaws_param_key* out_key);

char const* qaws_diff_field_name(qaws_diff_field field);
qaws_diff_field qaws_diff_field_from_name(char const* name);

/* ===================================================================
 * Curve schema
 * =================================================================== */

unsigned int qaws_curve_get_diff_capabilities(qaws_curve const* curve);
qaws_diff_class qaws_curve_get_diff_class(qaws_curve const* curve);
qaws_coordinate_kind qaws_curve_get_coordinate_kind(qaws_curve const* curve);

qaws_status qaws_curve_describe_fields(
	qaws_curve const* curve,
	qaws_field_desc* out_fields,
	unsigned int capacity,
	unsigned int* out_count);

/* Copy the primal values of one field (count * components scalars). */
qaws_status qaws_curve_read_field(
	qaws_curve const* curve,
	qaws_diff_field field,
	qaws_scalar* out_values,
	unsigned int capacity,
	unsigned int* out_scalar_count);

/* Parameters influencing the evaluation at t. For QAWS_CAP_LINEAR curves
   the weights are exact basis derivatives for orders 0..order. */
qaws_status qaws_curve_local_support(
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int order,
	qaws_local_support* out_support);

/*
 * Support index for gather-style adjoint accumulation (CPU or GPU):
 * for every element e of `field`, the samples influenced by it are
 * out_samples[out_offsets[e] .. out_offsets[e+1]).
 * out_offsets needs element_count + 1 entries. When capacity is too small
 * the function returns QAWS_STATUS_BUFFER_TOO_SMALL and *out_entry_count
 * holds the required size.
 */
qaws_status qaws_curve_build_support_index(
	qaws_curve const* curve,
	qaws_scalar const* t,
	unsigned int count,
	qaws_diff_field field,
	unsigned int* out_offsets,
	unsigned int offset_capacity,
	unsigned int* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_entry_count);

/* ===================================================================
 * Curve evaluation: batch
 *
 * t_tangent may be NULL (coordinate held fixed); param_tangent may be NULL
 * (parameters held fixed). out_primal may be NULL.
 * =================================================================== */

qaws_status qaws_curve_eval_batch_tangent_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal,
	qaws_curve_jet_2d* out_tangent);

qaws_status qaws_curve_eval_batch_tangent_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal,
	qaws_curve_jet_3d* out_tangent);

qaws_status qaws_curve_eval_batch_tangent2_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal,
	qaws_curve_jet_2d* out_tangent,
	qaws_curve_jet_2d* out_tangent2);

qaws_status qaws_curve_eval_batch_tangent2_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal,
	qaws_curve_jet_3d* out_tangent,
	qaws_curve_jet_3d* out_tangent2);

/* Accumulates (+=) into param_adjoint and t_adjoint (either may be NULL).
   Only the channels in `channels` of each out_adjoint jet are read. */
qaws_status qaws_curve_eval_batch_adjoint_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	unsigned int count,
	unsigned int channels,
	qaws_curve_jet_2d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);

qaws_status qaws_curve_eval_batch_adjoint_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	unsigned int count,
	unsigned int channels,
	qaws_curve_jet_3d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);

/* ===================================================================
 * Curve evaluation: single sample
 * =================================================================== */

qaws_status qaws_curve_eval_tangent_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	qaws_scalar t_tangent,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal,
	qaws_curve_jet_2d* out_tangent);

qaws_status qaws_curve_eval_tangent_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	qaws_scalar t_tangent,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal,
	qaws_curve_jet_3d* out_tangent);

qaws_status qaws_curve_eval_adjoint_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int channels,
	qaws_curve_jet_2d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);

qaws_status qaws_curve_eval_adjoint_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int channels,
	qaws_curve_jet_3d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);

/* ===================================================================
 * Surface schema
 * =================================================================== */

unsigned int qaws_surface_get_diff_capabilities(qaws_surface const* surface);
qaws_diff_class qaws_surface_get_diff_class(qaws_surface const* surface);

qaws_status qaws_surface_describe_fields(
	qaws_surface const* surface,
	qaws_field_desc* out_fields,
	unsigned int capacity,
	unsigned int* out_count);

qaws_status qaws_surface_read_field(
	qaws_surface const* surface,
	qaws_diff_field field,
	qaws_scalar* out_values,
	unsigned int capacity,
	unsigned int* out_scalar_count);

/* Objects the surface rules chain into; child i matches views->children[i]. */
qaws_status qaws_surface_diff_children(
	qaws_surface const* surface,
	qaws_diff_child* out_children,
	unsigned int capacity,
	unsigned int* out_count);

qaws_status qaws_surface_local_support(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int order,
	qaws_surface_support* out_support);

qaws_status qaws_surface_build_support_index(
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	unsigned int count,
	qaws_diff_field field,
	unsigned int* out_offsets,
	unsigned int offset_capacity,
	unsigned int* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_entry_count);

/* ===================================================================
 * Surface jets and evaluation
 *
 * Channels are qaws_surface_jet_channel bits (all partials up to total
 * order 3). (u, v) are clamped to the surface domain like
 * qaws_surface_evaluate.
 * =================================================================== */

/* Primal spatial jet (analytic, no finite differences). */
qaws_status qaws_surface_eval_jet(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int channels,
	qaws_surface_jet* out_jet);

qaws_status qaws_surface_eval_batch_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	qaws_scalar const* u_tangent,
	qaws_scalar const* v_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal,
	qaws_surface_jet* out_tangent);

qaws_status qaws_surface_eval_batch_tangent2(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	qaws_scalar const* u_tangent,
	qaws_scalar const* v_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal,
	qaws_surface_jet* out_tangent,
	qaws_surface_jet* out_tangent2);

qaws_status qaws_surface_eval_batch_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	unsigned int count,
	unsigned int channels,
	qaws_surface_jet const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* u_adjoint,
	qaws_scalar* v_adjoint);

qaws_status qaws_surface_eval_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_scalar u_tangent,
	qaws_scalar v_tangent,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal,
	qaws_surface_jet* out_tangent);

qaws_status qaws_surface_eval_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int channels,
	qaws_surface_jet const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* u_adjoint,
	qaws_scalar* v_adjoint);

#endif /* QAWS_DIFF_H */
