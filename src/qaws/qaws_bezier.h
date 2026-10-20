#ifndef QAWS_BEZIER_H
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
	qaws_curve** out_curve,
	qaws_bezier_desc const* desc);

#endif /* QAWS_BEZIER_H */
