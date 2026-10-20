#ifndef QAWS_HERMITE_H
#define QAWS_HERMITE_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_hermite_desc
{
	qaws_dimension dimension;
	unsigned int degree;
	void const* points;
	void const* derivatives;
	unsigned int point_count;
	unsigned int derivative_count;
} qaws_hermite_desc;

qaws_status qaws_curve_create_hermite(
	qaws_curve** out_curve,
	qaws_hermite_desc const* desc);

#endif /* QAWS_HERMITE_H */
