#ifndef QAWS_IMPORT_H
#define QAWS_IMPORT_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Polyline import.
 *
 * Create a curve from a raw point array (polyline).
 *
 * qaws_polyline_import_catmull_rom:
 *   Creates a Catmull-Rom spline that interpolates through the given points.
 *   points:           flat array of qaws_vec2 or qaws_vec3.
 *   point_count:      number of points (>= 2).
 *   parameterization: knot spacing (0 = centripetal default).
 *   closed:           non-zero for a closed loop.
 *
 * qaws_polyline_import_trajectory:
 *   Creates a cubic trajectory that interpolates through the given points
 *   with uniform timing.
 *   points:           flat array of qaws_vec2 or qaws_vec3.
 *   point_count:      number of points (>= 2).
 *   closed:           non-zero for a closed loop.
 */
qaws_status qaws_polyline_import_catmull_rom(
	qaws_curve** out_curve,
	void const* points,
	qaws_dimension dimension,
	unsigned int point_count,
	qaws_parameterization parameterization,
	int closed);

qaws_status qaws_polyline_import_trajectory(
	qaws_curve** out_curve,
	void const* points,
	qaws_dimension dimension,
	unsigned int point_count,
	int closed);

#endif /* QAWS_IMPORT_H */
