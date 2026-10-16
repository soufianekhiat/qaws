#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_NURBS_H_F32
#define QAWS_SURFACE_NURBS_H_F32
#define QAWS_SURFACE_NURBS_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_NURBS_H_F64
#define QAWS_SURFACE_NURBS_H_F64
#define QAWS_SURFACE_NURBS_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_NURBS_H_BODY
#undef QAWS_SURFACE_NURBS_H_BODY
#define QAWS_SURFACE_NURBS_H

#include "qaws_surface_types.h"

typedef struct qaws_surface_nurbs_desc
{
	unsigned int u_degree;
	unsigned int v_degree;
	/* Control point grid: u_point_count x v_point_count qaws_vec3 values. */
	qaws_vec3 const* control_points;
	unsigned int u_point_count;
	unsigned int v_point_count;
	/* Weight grid: u_point_count x v_point_count scalars, same layout as CPs. */
	qaws_scalar const* weights;
	qaws_scalar const* u_knots;
	unsigned int u_knot_count;   /* 0 for auto uniform clamped */
	qaws_scalar const* v_knots;
	unsigned int v_knot_count;   /* 0 for auto uniform clamped */
} qaws_surface_nurbs_desc;

qaws_status qaws_surface_create_nurbs(
	qaws_surface_nurbs_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_NURBS_H */
