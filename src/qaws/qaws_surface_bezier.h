#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_BEZIER_H_F32
#define QAWS_SURFACE_BEZIER_H_F32
#define QAWS_SURFACE_BEZIER_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_BEZIER_H_F64
#define QAWS_SURFACE_BEZIER_H_F64
#define QAWS_SURFACE_BEZIER_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_BEZIER_H_BODY
#undef QAWS_SURFACE_BEZIER_H_BODY
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
	qaws_surface_bezier_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_BEZIER_H */
