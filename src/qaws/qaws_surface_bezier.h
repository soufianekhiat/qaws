#ifndef QAWS_SURFACE_BEZIER_H
#define QAWS_SURFACE_BEZIER_H

#include "qaws_surface_types.h"

typedef struct qaws_surface_bezier_desc
{
	unsigned int u_degree;
	unsigned int v_degree;
	/* Control point grid: (u_degree+1) x (v_degree+1) qaws_vec3 values.
	   Row-major: row i contains points P[i][0..v_degree], where i is the u index. */
	qaws_vec3 const* control_points;
	unsigned int u_point_count;   /* must be u_degree + 1 */
	unsigned int v_point_count;   /* must be v_degree + 1 */
} qaws_surface_bezier_desc;

qaws_status qaws_surface_create_bezier(
	qaws_surface** out_surface,
	qaws_surface_bezier_desc const* desc);

#endif /* QAWS_SURFACE_BEZIER_H */
