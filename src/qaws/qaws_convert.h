#ifndef QAWS_CONVERT_H
#define QAWS_CONVERT_H

#include "qaws_types.h"
#include "qaws_status.h"

/* Thread-safe: operates on immutable input curves. */

/* Family conversions */
qaws_status qaws_curve_convert_hermite_to_bezier(
	qaws_curve** out_bezier,
	qaws_curve const* curve,
	unsigned int span_index);

qaws_status qaws_curve_convert_catmull_rom_to_bezier(
	qaws_curve** out_bezier,
	qaws_curve const* curve,
	unsigned int span_index);

qaws_status qaws_curve_convert_bezier_to_bspline(
	qaws_curve** out_bspline,
	qaws_curve const* curve);

qaws_status qaws_curve_convert_bspline_to_nurbs(
	qaws_curve** out_nurbs,
	qaws_curve const* curve);

/* Degree elevation */
qaws_status qaws_curve_elevate_degree(
	qaws_curve** out_elevated,
	qaws_curve const* curve);

/* Degree reduction (Bezier only, n -> n-1) */
qaws_status qaws_curve_reduce_degree(
	qaws_curve** out_reduced,
	qaws_curve const* curve);

#endif /* QAWS_CONVERT_H */
