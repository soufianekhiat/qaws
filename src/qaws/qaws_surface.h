#ifndef QAWS_SURFACE_H
#define QAWS_SURFACE_H

#include "qaws_surface_types.h"

/* Generic surface API. All functions are thread-safe on immutable surfaces. */

qaws_surface_kind qaws_surface_get_kind(qaws_surface const* surface);
unsigned int qaws_surface_get_u_degree(qaws_surface const* surface);
unsigned int qaws_surface_get_v_degree(qaws_surface const* surface);
qaws_range qaws_surface_get_u_range(qaws_surface const* surface);
qaws_range qaws_surface_get_v_range(qaws_surface const* surface);
int qaws_surface_is_rational(qaws_surface const* surface);

qaws_status qaws_surface_evaluate(
	qaws_surface_eval_result* out_result,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags);

qaws_status qaws_surface_get_normal(
	qaws_vec3* out_normal,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v);

/* Sample a grid of positions on the surface.
   out_positions must hold u_count * v_count elements (row-major: v varies fastest). */
qaws_status qaws_surface_sample(
	qaws_vec3* out_positions,
	unsigned int* out_position_count,
	qaws_surface const* surface,
	unsigned int u_count,
	unsigned int v_count,
	unsigned int position_capacity);

void qaws_surface_destroy(qaws_surface* surface);

#endif /* QAWS_SURFACE_H */
