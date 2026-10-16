#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_CATMULL_ROM_H_F32
#define QAWS_CATMULL_ROM_H_F32
#define QAWS_CATMULL_ROM_H_BODY
#endif
#else
#ifndef QAWS_CATMULL_ROM_H_F64
#define QAWS_CATMULL_ROM_H_F64
#define QAWS_CATMULL_ROM_H_BODY
#endif
#endif
#ifdef QAWS_CATMULL_ROM_H_BODY
#undef QAWS_CATMULL_ROM_H_BODY
#define QAWS_CATMULL_ROM_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_catmull_rom_desc
{
	qaws_dimension dimension;
	void const* control_points;
	unsigned int control_point_count;
	qaws_parameterization parameterization;
	int closed;
} qaws_catmull_rom_desc;

qaws_status qaws_curve_create_catmull_rom(
	qaws_curve** out_curve,
	qaws_catmull_rom_desc const* desc);

#endif /* QAWS_CATMULL_ROM_H */
