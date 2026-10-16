#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_BEZIER_H_F32
#define QAWS_BEZIER_H_F32
#define QAWS_BEZIER_H_BODY
#endif
#else
#ifndef QAWS_BEZIER_H_F64
#define QAWS_BEZIER_H_F64
#define QAWS_BEZIER_H_BODY
#endif
#endif
#ifdef QAWS_BEZIER_H_BODY
#undef QAWS_BEZIER_H_BODY
#define QAWS_BEZIER_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_bezier_desc
{
	qaws_dimension dimension;
	unsigned int degree;
	void const* control_points;
	unsigned int control_point_count;
} qaws_bezier_desc;

qaws_status qaws_curve_create_bezier(
	qaws_bezier_desc const* desc,
	qaws_curve** out_curve);

#endif /* QAWS_BEZIER_H */
