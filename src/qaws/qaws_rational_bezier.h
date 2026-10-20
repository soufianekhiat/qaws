#ifndef QAWS_RATIONAL_BEZIER_H
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
	qaws_curve** out_curve,
	qaws_rational_bezier_desc const* desc);

#endif /* QAWS_RATIONAL_BEZIER_H */
