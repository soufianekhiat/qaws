#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_NURBS_H_F32
#define QAWS_NURBS_H_F32
#define QAWS_NURBS_H_BODY
#endif
#else
#ifndef QAWS_NURBS_H_F64
#define QAWS_NURBS_H_F64
#define QAWS_NURBS_H_BODY
#endif
#endif
#ifdef QAWS_NURBS_H_BODY
#undef QAWS_NURBS_H_BODY
#define QAWS_NURBS_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_nurbs_desc
{
	qaws_dimension dimension;
	unsigned int degree;
	void const* control_points;
	unsigned int control_point_count;
	qaws_scalar const* knots;
	unsigned int knot_count;
	qaws_scalar const* weights;
	unsigned int weight_count;
	int is_closed;
} qaws_nurbs_desc;

qaws_status qaws_curve_create_nurbs(
	qaws_nurbs_desc const* desc,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_nurbs_ex(
	qaws_nurbs_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

#endif /* QAWS_NURBS_H */
