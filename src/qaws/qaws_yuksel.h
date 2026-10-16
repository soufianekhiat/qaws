#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_YUKSEL_H_F32
#define QAWS_YUKSEL_H_F32
#define QAWS_YUKSEL_H_BODY
#endif
#else
#ifndef QAWS_YUKSEL_H_F64
#define QAWS_YUKSEL_H_F64
#define QAWS_YUKSEL_H_BODY
#endif
#endif
#ifdef QAWS_YUKSEL_H_BODY
#undef QAWS_YUKSEL_H_BODY
#define QAWS_YUKSEL_H

#include "qaws_types.h"
#include "qaws_status.h"

#ifndef QAWS_ENUM_QAWS_YUKSEL_MODE
#define QAWS_ENUM_QAWS_YUKSEL_MODE
typedef enum qaws_yuksel_mode
{
	QAWS_YUKSEL_MODE_BEZIER = 0,
	QAWS_YUKSEL_MODE_CIRCULAR,
	QAWS_YUKSEL_MODE_ELLIPTICAL,
	QAWS_YUKSEL_MODE_HYBRID
} qaws_yuksel_mode;
#endif

typedef struct qaws_yuksel_desc
{
	qaws_dimension dimension;
	void const* control_points;
	unsigned int control_point_count;
	qaws_yuksel_mode mode;
	int closed;
} qaws_yuksel_desc;

qaws_status qaws_curve_create_yuksel(
	qaws_yuksel_desc const* desc,
	qaws_curve** out_curve);

#endif /* QAWS_YUKSEL_H */
