#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_HERMITE_H_F32
#define QAWS_HERMITE_H_F32
#define QAWS_HERMITE_H_BODY
#endif
#else
#ifndef QAWS_HERMITE_H_F64
#define QAWS_HERMITE_H_F64
#define QAWS_HERMITE_H_BODY
#endif
#endif
#ifdef QAWS_HERMITE_H_BODY
#undef QAWS_HERMITE_H_BODY
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
	qaws_hermite_desc const* desc,
	qaws_curve** out_curve);

#endif /* QAWS_HERMITE_H */
