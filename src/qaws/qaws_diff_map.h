#ifndef QAWS_DIFF_MAP_H
#define QAWS_DIFF_MAP_H

#include "qaws_diff_types.h"

/*
 * Differential maps of geometry-building operations.
 *
 * An operation that builds new geometry (split, join, conversion, degree
 * change, fitting, ...) returns, next to its output objects, a map from
 * the parameters of its inputs to the parameters of its outputs:
 *
 *   qaws_diff_map_tangent : out = J in        (outputs are overwritten)
 *   qaws_diff_map_adjoint : in += J^T out
 *
 * Inputs and outputs are numbered objects, each addressed through its own
 * qaws_diff_views (NULL entries are inactive). Operations document their
 * numbering, e.g. split: inputs { curve, parameter }, outputs { left, right }.
 *
 * Maps are either materialized sparse Jacobians (exact, evaluated at the
 * operation's inputs) or operators that apply J without forming it.
 */

typedef struct qaws_diff_map qaws_diff_map;

typedef enum qaws_diff_map_kind
{
	QAWS_DIFF_MAP_LINEAR_SPARSE = 0,
	QAWS_DIFF_MAP_OPERATOR
} qaws_diff_map_kind;

/* Component value meaning "the same weight on every component". */
#define QAWS_DIFF_MAP_ALL_COMPONENTS 0xFFu

typedef struct qaws_diff_map_entry
{
	unsigned int out_object, out_element;
	unsigned int in_object, in_element;
	unsigned char out_field, in_field;         /* qaws_diff_field */
	unsigned char out_component, in_component; /* or QAWS_DIFF_MAP_ALL_COMPONENTS */
	qaws_scalar weight;
} qaws_diff_map_entry;

qaws_diff_map_kind qaws_diff_map_get_kind(qaws_diff_map const* map);
unsigned int qaws_diff_map_input_count(qaws_diff_map const* map);
unsigned int qaws_diff_map_output_count(qaws_diff_map const* map);

qaws_status qaws_diff_map_tangent(
	qaws_diff_map const* map,
	qaws_diff_context const* ctx,
	qaws_diff_views const* const* in_tangents,
	unsigned int in_count,
	qaws_diff_views* const* out_tangents,
	unsigned int out_count);

qaws_status qaws_diff_map_adjoint(
	qaws_diff_map const* map,
	qaws_diff_context const* ctx,
	qaws_diff_views const* const* out_adjoints,
	unsigned int out_count,
	qaws_diff_views* const* in_adjoints,
	unsigned int in_count);

/* Sparse entries of a LINEAR_SPARSE map (BUFFER_TOO_SMALL reports the size). */
qaws_status qaws_diff_map_get_entries(
	qaws_diff_map const* map,
	qaws_diff_map_entry* out_entries,
	unsigned int capacity,
	unsigned int* out_count);

void qaws_diff_map_destroy(qaws_diff_map* map);

/* ===================================================================
 * Operations returning maps
 *
 * Polynomial families are exact linear maps of their fields; NURBS inputs
 * are linear in homogeneous coordinates and their maps include the
 * projection. The split parameter is differentiable for Bezier curves;
 * for other families a split parameter tangent is refused.
 * =================================================================== */

/* inputs: 0 = curve, 1 = split parameter; outputs: 0 = left, 1 = right */
qaws_status qaws_curve_split_diff(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_curve** out_left,
	qaws_curve** out_right,
	qaws_diff_map** out_map);

/* inputs: 0 = curve_a, 1 = curve_b; outputs: 0 = joined */
qaws_status qaws_curve_join_diff(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_curve** out_joined,
	qaws_diff_map** out_map);

/* inputs: 0 = curve; outputs: 0 = converted curve */
qaws_status qaws_curve_convert_hermite_to_bezier_diff(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_curve** out_bezier,
	qaws_diff_map** out_map);

qaws_status qaws_curve_convert_bezier_to_bspline_diff(
	qaws_curve const* curve,
	qaws_curve** out_bspline,
	qaws_diff_map** out_map);

qaws_status qaws_curve_convert_bspline_to_nurbs_diff(
	qaws_curve const* curve,
	qaws_curve** out_nurbs,
	qaws_diff_map** out_map);

qaws_status qaws_curve_elevate_degree_diff(
	qaws_curve const* curve,
	qaws_curve** out_elevated,
	qaws_diff_map** out_map);

qaws_status qaws_curve_reduce_degree_diff(
	qaws_curve const* curve,
	qaws_curve** out_reduced,
	qaws_diff_map** out_map);

#endif /* QAWS_DIFF_MAP_H */
