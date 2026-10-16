#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_IMPORT_H_F32
#define QAWS_IMPORT_H_F32
#define QAWS_IMPORT_H_BODY
#endif
#else
#ifndef QAWS_IMPORT_H_F64
#define QAWS_IMPORT_H_F64
#define QAWS_IMPORT_H_BODY
#endif
#endif
#ifdef QAWS_IMPORT_H_BODY
#undef QAWS_IMPORT_H_BODY
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
	qaws_dimension dimension,
	void const* points,
	unsigned int point_count,
	qaws_parameterization parameterization,
	int closed,
	qaws_curve** out_curve);

qaws_status qaws_polyline_import_trajectory(
	qaws_dimension dimension,
	void const* points,
	unsigned int point_count,
	int closed,
	qaws_curve** out_curve);

#endif /* QAWS_IMPORT_H */
