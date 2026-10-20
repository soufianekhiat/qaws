#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SUBDIVISION_H_F32
#define QAWS_SUBDIVISION_H_F32
#define QAWS_SUBDIVISION_H_BODY
#endif
#else
#ifndef QAWS_SUBDIVISION_H_F64
#define QAWS_SUBDIVISION_H_F64
#define QAWS_SUBDIVISION_H_BODY
#endif
#endif
#ifdef QAWS_SUBDIVISION_H_BODY
#undef QAWS_SUBDIVISION_H_BODY
#define QAWS_SUBDIVISION_H

#include "qaws_types.h"
#include "qaws_status.h"

#ifndef QAWS_ENUM_QAWS_SUBDIVISION_SCHEME
#define QAWS_ENUM_QAWS_SUBDIVISION_SCHEME
typedef enum qaws_subdivision_scheme
{
	QAWS_SUBDIVISION_CHAIKIN = 0,        /* Quadratic B-spline limit (C1) */
	QAWS_SUBDIVISION_LANE_RIESENFELD_3,  /* Cubic B-spline limit (C2) */
	QAWS_SUBDIVISION_LANE_RIESENFELD_4   /* Quartic B-spline limit (C3) */
} qaws_subdivision_scheme;
#endif

typedef struct qaws_subdivision_desc
{
	qaws_dimension dimension;
	qaws_subdivision_scheme scheme;
	void const* control_points;          /* initial control polygon */
	unsigned int control_point_count;    /* >= 3 */
	int closed;                          /* 1 for closed polygon */
	unsigned int refinement_levels;      /* number of subdivision iterations (default 6) */
} qaws_subdivision_desc;

qaws_status qaws_curve_create_subdivision(
	qaws_curve** out_curve,
	qaws_subdivision_desc const* desc);

#endif /* QAWS_SUBDIVISION_H */
