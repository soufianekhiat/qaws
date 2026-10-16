#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_ARC_H_F32
#define QAWS_ARC_H_F32
#define QAWS_ARC_H_BODY
#endif
#else
#ifndef QAWS_ARC_H_F64
#define QAWS_ARC_H_F64
#define QAWS_ARC_H_BODY
#endif
#endif
#ifdef QAWS_ARC_H_BODY
#undef QAWS_ARC_H_BODY
#define QAWS_ARC_H

#include "qaws_types.h"
#include "qaws_status.h"

/* A single arc segment defined by center, radius, start angle, end angle.
 * For 2D: angles are in radians, CCW from +x axis.
 * For 3D: arcs lie in a plane defined by axis_u and axis_v vectors.
 * Radius must be > 0 (degenerate zero-radius segments are not supported). */
typedef struct qaws_arc_segment
{
	qaws_scalar center[3];
	qaws_scalar radius;
	qaws_scalar angle_start;    /* radians */
	qaws_scalar angle_end;      /* radians */
	/* 3D only: plane basis vectors (ignored for 2D) */
	qaws_scalar axis_u[3];     /* unit vector: cos direction */
	qaws_scalar axis_v[3];     /* unit vector: sin direction */
} qaws_arc_segment;

typedef struct qaws_arc_desc
{
	qaws_dimension dimension;
	qaws_arc_segment const* segments;
	unsigned int segment_count;    /* >= 1 */
} qaws_arc_desc;

qaws_status qaws_curve_create_arc(
	qaws_curve** out_curve,
	qaws_arc_desc const* desc);

#endif /* QAWS_ARC_H */
