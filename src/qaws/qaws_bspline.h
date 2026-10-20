#ifndef QAWS_BSPLINE_H
#define QAWS_BSPLINE_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_bspline_desc
{
	qaws_dimension dimension;
	unsigned int degree;
	void const* control_points;
	unsigned int control_point_count;
	qaws_scalar const* knots;
	unsigned int knot_count;
	int is_uniform;
	int is_closed;
} qaws_bspline_desc;

qaws_status qaws_curve_create_bspline(
	qaws_curve** out_curve,
	qaws_bspline_desc const* desc);

qaws_status qaws_curve_create_bspline_ex(
	qaws_curve** out_curve,
	qaws_bspline_desc const* desc,
	qaws_allocator const* allocator);

#endif /* QAWS_BSPLINE_H */
