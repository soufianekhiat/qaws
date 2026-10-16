#ifndef QAWS_INTERNAL_ARC_LENGTH_H
#define QAWS_INTERNAL_ARC_LENGTH_H

#include "qaws_internal_types.h"

/* Compute arc length of a curve segment between two global parameter values
   using Gaussian quadrature. */
qaws_status qaws_internal_arc_length(
	qaws_scalar* out_length,
	qaws_curve const* curve,
	qaws_scalar t_min,
	qaws_scalar t_max);

/* Build a cumulative arc-length table for a curve.
   table_params and table_distances must each hold table_size entries.
   table_params[i] are uniformly spaced parameters from min to max.
   table_distances[i] are cumulative arc lengths. */
qaws_status qaws_internal_build_arc_length_table(
	qaws_scalar* table_params,
	qaws_scalar* table_distances,
	qaws_curve const* curve,
	unsigned int table_size);

/* Map a distance value to a parameter using an arc-length table.
   Uses binary search + linear interpolation. */
qaws_scalar qaws_internal_distance_to_parameter(
	qaws_scalar const* table_params,
	qaws_scalar const* table_distances,
	unsigned int table_size,
	qaws_scalar distance);

/* Map a parameter to a cumulative distance using an arc-length table. */
qaws_scalar qaws_internal_parameter_to_distance(
	qaws_scalar const* table_params,
	qaws_scalar const* table_distances,
	unsigned int table_size,
	qaws_scalar parameter);

#endif /* QAWS_INTERNAL_ARC_LENGTH_H */
