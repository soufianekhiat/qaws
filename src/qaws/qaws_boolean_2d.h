#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_BOOLEAN_2D_H_F32
#define QAWS_BOOLEAN_2D_H_F32
#define QAWS_BOOLEAN_2D_H_BODY
#endif
#else
#ifndef QAWS_BOOLEAN_2D_H_F64
#define QAWS_BOOLEAN_2D_H_F64
#define QAWS_BOOLEAN_2D_H_BODY
#endif
#endif
#ifdef QAWS_BOOLEAN_2D_H_BODY
#undef QAWS_BOOLEAN_2D_H_BODY
#define QAWS_BOOLEAN_2D_H

#include "qaws_types.h"
#include "qaws_status.h"

#ifndef QAWS_ENUM_QAWS_BOOLEAN_OP
#define QAWS_ENUM_QAWS_BOOLEAN_OP
typedef enum qaws_boolean_op {
	QAWS_BOOLEAN_UNION = 0,
	QAWS_BOOLEAN_INTERSECTION,
	QAWS_BOOLEAN_DIFFERENCE
} qaws_boolean_op;
#endif

#endif /* QAWS_BOOLEAN_2D_H */
