#ifndef QAWS_NURBS_H
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
	qaws_curve** out_curve,
	qaws_nurbs_desc const* desc);

qaws_status qaws_curve_create_nurbs_ex(
	qaws_curve** out_curve,
	qaws_nurbs_desc const* desc,
	qaws_allocator const* allocator);

#endif /* QAWS_NURBS_H */
