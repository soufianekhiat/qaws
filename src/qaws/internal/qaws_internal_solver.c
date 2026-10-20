#include "qaws_internal_solver.h"
#include <math.h>

qaws_scalar qaws_internal_closest_point_newton(
	qaws_scalar const* target,
	unsigned int dim_count,
	qaws_eval_pos_d1_fn eval_fn,
	qaws_scalar t_min,
	qaws_scalar t_max,
	qaws_scalar t_initial,
	unsigned int max_iterations,
	void const* ctx)
{
	qaws_scalar t = t_initial;
	qaws_scalar tolerance = (qaws_scalar)1e-7;

	qaws_scalar pos[3];
	qaws_scalar d1[3];

	for (unsigned int iter = 0; iter < max_iterations; iter++) {
		eval_fn(ctx, t, pos, d1, dim_count);

		qaws_scalar dot_diff_d1 = (qaws_scalar)0;
		qaws_scalar dot_d1_d1 = (qaws_scalar)0;

		for (unsigned int d = 0; d < dim_count; d++) {
			qaws_scalar diff = pos[d] - target[d];
			dot_diff_d1 += diff * d1[d];
			dot_d1_d1 += d1[d] * d1[d];
		}

		if (dot_d1_d1 < (qaws_scalar)1e-14)
			break;

		qaws_scalar dt = dot_diff_d1 / dot_d1_d1;
		t = t - dt;

		if (t < t_min) t = t_min;
		if (t > t_max) t = t_max;

		if ((qaws_scalar)fabs((double)dt) < tolerance)
			break;
	}

	return t;
}
