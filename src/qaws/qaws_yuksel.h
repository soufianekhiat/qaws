#ifndef QAWS_YUKSEL_H
#define QAWS_YUKSEL_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef enum qaws_yuksel_mode
{
	QAWS_YUKSEL_MODE_BEZIER = 0,
	QAWS_YUKSEL_MODE_CIRCULAR,
	QAWS_YUKSEL_MODE_ELLIPTICAL,
	QAWS_YUKSEL_MODE_HYBRID
} qaws_yuksel_mode;

typedef struct qaws_yuksel_desc
{
	qaws_dimension dimension;
	void const* control_points;
	unsigned int control_point_count;
	qaws_yuksel_mode mode;
	int closed;
} qaws_yuksel_desc;

qaws_status qaws_curve_create_yuksel(
	qaws_curve** out_curve,
	qaws_yuksel_desc const* desc);

#endif /* QAWS_YUKSEL_H */
