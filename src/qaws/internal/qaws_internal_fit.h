#ifndef QAWS_INTERNAL_FIT_H
#define QAWS_INTERNAL_FIT_H

#include "qaws_internal_types.h"

/*
 * Least-squares B-spline fitting from sample points.
 *
 * Given N sample points in 2D or 3D, fit a B-spline of the specified degree
 * with a given number of control points.  The knot vector is computed from
 * chord-length parameterization.
 *
 * The caller provides sample_params[sample_count] (parameter values for each
 * sample, in ascending order) and sample_coords[sample_count * dim_count]
 * (interleaved coordinate values).
 *
 * The output is a newly created qaws_curve (B-spline).
 */

qaws_status qaws_internal_fit_bspline(
	qaws_curve** out_curve,
	qaws_scalar const* sample_params,
	qaws_scalar const* sample_coords,
	qaws_dimension dimension,
	unsigned int sample_count,
	unsigned int control_point_count,
	unsigned int degree);

/*
 * Fit a B-spline through sample points with specified parameter range.
 * Same as above but the knot vector endpoints are set to [t_min, t_max].
 */
qaws_status qaws_internal_fit_bspline_range(
	qaws_curve** out_curve,
	qaws_scalar const* sample_params,
	qaws_scalar const* sample_coords,
	qaws_dimension dimension,
	unsigned int sample_count,
	unsigned int control_point_count,
	unsigned int degree,
	qaws_scalar t_min,
	qaws_scalar t_max);

#endif /* QAWS_INTERNAL_FIT_H */
