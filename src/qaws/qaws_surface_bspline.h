#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_BSPLINE_H_F32
#define QAWS_SURFACE_BSPLINE_H_F32
#define QAWS_SURFACE_BSPLINE_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_BSPLINE_H_F64
#define QAWS_SURFACE_BSPLINE_H_F64
#define QAWS_SURFACE_BSPLINE_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_BSPLINE_H_BODY
#undef QAWS_SURFACE_BSPLINE_H_BODY
#define QAWS_SURFACE_BSPLINE_H

#include "qaws_surface_types.h"

typedef struct qaws_surface_bspline_desc
{
	unsigned int u_degree;
	unsigned int v_degree;
	/* Control point grid: u_point_count x v_point_count qaws_vec3 values.
	   Row-major: row i contains points P[i][0..v_point_count-1]. */
	qaws_vec3 const* control_points;
	unsigned int u_point_count;
	unsigned int v_point_count;
	qaws_scalar const* u_knots;
	unsigned int u_knot_count;   /* 0 for auto uniform clamped */
	qaws_scalar const* v_knots;
	unsigned int v_knot_count;   /* 0 for auto uniform clamped */
} qaws_surface_bspline_desc;

qaws_status qaws_surface_create_bspline(
	qaws_surface** out_surface,
	qaws_surface_bspline_desc const* desc);

#endif /* QAWS_SURFACE_BSPLINE_H */
