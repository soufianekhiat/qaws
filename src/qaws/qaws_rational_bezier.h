#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_RATIONAL_BEZIER_H_F32
#define QAWS_RATIONAL_BEZIER_H_F32
#define QAWS_RATIONAL_BEZIER_H_BODY
#endif
#else
#ifndef QAWS_RATIONAL_BEZIER_H_F64
#define QAWS_RATIONAL_BEZIER_H_F64
#define QAWS_RATIONAL_BEZIER_H_BODY
#endif
#endif
#ifdef QAWS_RATIONAL_BEZIER_H_BODY
#undef QAWS_RATIONAL_BEZIER_H_BODY
#define QAWS_RATIONAL_BEZIER_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_rational_bezier_desc
{
	qaws_dimension dimension;
	unsigned int degree;
	void const* control_points;      /* dim_count * (degree+1) scalars */
	unsigned int control_point_count; /* must be degree+1 */
	qaws_scalar const* weights;      /* degree+1 positive weights */
	unsigned int weight_count;        /* must be degree+1 */
} qaws_rational_bezier_desc;

qaws_status qaws_curve_create_rational_bezier(
	qaws_rational_bezier_desc const* desc,
	qaws_curve** out_curve);

#endif /* QAWS_RATIONAL_BEZIER_H */
