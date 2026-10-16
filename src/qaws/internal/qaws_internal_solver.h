#ifndef QAWS_INTERNAL_SOLVER_H
#define QAWS_INTERNAL_SOLVER_H

#include "../qaws_types.h"

/* Newton-Raphson iteration for closest-point projection.
   Given a curve evaluation callback, find the parameter t that minimizes
   distance from the curve to a target point.
   Returns QAWS_STATUS_OK on convergence. */
typedef void (*qaws_eval_pos_d1_fn)(
	void const* ctx,
	qaws_scalar t,
	qaws_scalar* out_pos,
	qaws_scalar* out_d1,
	unsigned int dim_count);

qaws_scalar qaws_internal_closest_point_newton(
	qaws_scalar const* target,
	unsigned int dim_count,
	qaws_eval_pos_d1_fn eval_fn,
	qaws_scalar t_min,
	qaws_scalar t_max,
	qaws_scalar t_initial,
	unsigned int max_iterations,
	void const* ctx);

#endif /* QAWS_INTERNAL_SOLVER_H */
