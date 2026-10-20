#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_EXPORT_H_F32
#define QAWS_EXPORT_H_F32
#define QAWS_EXPORT_H_BODY
#endif
#else
#ifndef QAWS_EXPORT_H_F64
#define QAWS_EXPORT_H_F64
#define QAWS_EXPORT_H_BODY
#endif
#endif
#ifdef QAWS_EXPORT_H_BODY
#undef QAWS_EXPORT_H_BODY
#define QAWS_EXPORT_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"

typedef struct qaws_bspline_fit_desc
{
	qaws_dimension dimension;
	void const* data_points;           /* m * dim_count scalars (flat array of qaws_vec2 or qaws_vec3) */
	unsigned int data_point_count;     /* m >= degree + 1 */
	unsigned int degree;               /* desired B-spline degree (typically 3) */
	unsigned int control_point_count;  /* n: desired number of CPs (degree < n <= m) */
	qaws_scalar const* parameters;     /* optional: m parameter values in [0,1], NULL = from the points as below */
	qaws_parameterization parameterization;  /* without parameters: DEFAULT and CHORDAL chord lengths, UNIFORM equal steps,
	                                            CENTRIPETAL square roots of the chord lengths */
} qaws_bspline_fit_desc;

qaws_status qaws_curve_fit_bspline(
	qaws_curve** out_curve,
	qaws_bspline_fit_desc const* desc);

/*
 * SVG path export.
 *
 * Convert a 2D curve to an SVG path data string (the "d" attribute).
 *
 * For Bezier: maps directly to M/C/Q/L commands (exact).
 * For Hermite/Catmull-Rom: converts each span to cubic Bezier (exact).
 * For all other families: samples the curve and fits piecewise cubic
 * Bezier segments (approximate, controlled by sample_count).
 *
 * out_path_data:   caller-allocated buffer for the SVG path string.
 * capacity:        buffer size in bytes (including null terminator).
 * out_length:      actual string length written (excluding null terminator).
 *                  If out_length >= capacity, the string was truncated.
 * sample_count:    number of sample points for non-Bezier families (0 = default 64).
 *
 * The curve must be 2D.
 */
qaws_status qaws_curve_export_svg_path(
	char* out_path_data,
	unsigned int* out_length,
	qaws_curve const* curve,
	unsigned int capacity,
	unsigned int sample_count);

/*
 * Polyline export sampling mode.
 *
 * QAWS_POLYLINE_UNIFORM:    uniform parameter spacing (default).
 * QAWS_POLYLINE_CURVATURE:  curvature-weighted — more points in high-curvature
 *                           regions, fewer in straight segments.
 */
#ifndef QAWS_ENUM_QAWS_POLYLINE_SAMPLING
#define QAWS_ENUM_QAWS_POLYLINE_SAMPLING
typedef enum qaws_polyline_sampling
{
	QAWS_POLYLINE_UNIFORM   = 0,
	QAWS_POLYLINE_CURVATURE = 1
} qaws_polyline_sampling;
#endif

/*
 * Polyline export.
 *
 * Sample any curve and output positions as a flat point array (polyline).
 *
 * curve:         any curve (2D or 3D, any family).
 * sample_count:  number of sample points (>= 2).
 * sampling:      QAWS_POLYLINE_UNIFORM or QAWS_POLYLINE_CURVATURE.
 * out_points:    caller-allocated buffer for qaws_vec2 (2D) or qaws_vec3 (3D).
 *                Must hold at least sample_count elements.
 * out_count:     receives the number of points actually written.
 *                For uniform sampling this equals sample_count.
 *                For curvature sampling this may be <= sample_count.
 */
qaws_status qaws_polyline_export_2d(
	qaws_vec2* out_points,
	unsigned int* out_count,
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_polyline_sampling sampling);

qaws_status qaws_polyline_export_3d(
	qaws_vec3* out_points,
	unsigned int* out_count,
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_polyline_sampling sampling);

/*
 * OBJ surface export.
 *
 * Export a surface as a Wavefront OBJ triangle mesh string.
 *
 * The surface is sampled on a u_samples x v_samples grid and triangulated.
 * Vertex normals are included when include_normals is non-zero.
 *
 * surface:         any surface.
 * u_samples:       grid density in u direction (>= 2).
 * v_samples:       grid density in v direction (>= 2).
 * include_normals: non-zero to emit "vn" lines and face normal indices.
 * out_obj_data:    caller-allocated buffer for the OBJ string.
 * capacity:        buffer size in bytes (including null terminator).
 * out_length:      actual string length written (excluding null).
 */
qaws_status qaws_surface_export_obj(
	char* out_obj_data,
	unsigned int* out_length,
	qaws_surface const* surface,
	unsigned int capacity,
	unsigned int u_samples,
	unsigned int v_samples,
	int include_normals);

#endif /* QAWS_EXPORT_H */
