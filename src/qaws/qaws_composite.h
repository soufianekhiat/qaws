#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_COMPOSITE_H_F32
#define QAWS_COMPOSITE_H_F32
#define QAWS_COMPOSITE_H_BODY
#endif
#else
#ifndef QAWS_COMPOSITE_H_F64
#define QAWS_COMPOSITE_H_F64
#define QAWS_COMPOSITE_H_BODY
#endif
#endif
#ifdef QAWS_COMPOSITE_H_BODY
#undef QAWS_COMPOSITE_H_BODY
#define QAWS_COMPOSITE_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Composite / multi-segment curve.
 *
 * Chains multiple heterogeneous qaws_curve objects end-to-end into a
 * single qaws_curve.  Each segment occupies one unit in the composite
 * parameter space: segment i maps to [i, i+1].  The total parameter
 * range is [0, segment_count].
 *
 * Ownership: the create function takes ownership of every curve pointer
 * in the segments array.  The caller must NOT destroy them afterwards;
 * they will be destroyed when the composite curve is destroyed.
 */

typedef struct qaws_composite_desc
{
	qaws_dimension dimension;
	qaws_curve** segments;           /* array of segment curves; ownership transfers to composite */
	unsigned int segment_count;      /* >= 1 */
} qaws_composite_desc;

qaws_status qaws_curve_create_composite(
	qaws_curve** out_curve,
	qaws_composite_desc const* desc);

#endif /* QAWS_COMPOSITE_H */
