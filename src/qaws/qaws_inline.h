#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_INLINE_H_F32
#define QAWS_INLINE_H_F32
#define QAWS_INLINE_H_BODY
#endif
#else
#ifndef QAWS_INLINE_H_F64
#define QAWS_INLINE_H_F64
#define QAWS_INLINE_H_BODY
#endif
#endif
#ifdef QAWS_INLINE_H_BODY
#undef QAWS_INLINE_H_BODY
#define QAWS_INLINE_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_bezier.h"
#include "qaws_polynomial.h"

/*
 * Stack-allocated (inline) curve initialization.
 *
 * These functions initialize a qaws_curve_inline struct without any
 * heap allocation. The resulting curve is usable with all evaluation
 * and inspection functions via &inline_curve.curve.
 *
 * Constraints:
 *   - Single span only (degree+1 control points for Bezier)
 *   - Maximum degree 7 (8 CPs x 3 dimensions = 24 scalars for CPs)
 *   - 2D or 3D
 *   - Bezier and Polynomial families only
 *
 * The qaws_curve_inline struct is self-contained. It can be copied
 * with memcpy or assigned. Do NOT call qaws_curve_destroy on it.
 *
 * Example:
 *   qaws_curve_inline ic;
 *   qaws_vec3 pts[] = {{0,0,0}, {1,2,3}, {4,5,6}, {7,8,9}};
 *   qaws_bezier_desc desc = { QAWS_DIMENSION_3D, 3, pts, 4 };
 *   qaws_curve_init_bezier_inline(&ic, &desc);
 *   // Use &ic.curve with qaws_curve_evaluate_3d, etc.
 */

qaws_status qaws_curve_init_bezier_inline(
	qaws_curve_inline* out,
	qaws_bezier_desc const* desc);

qaws_status qaws_curve_init_polynomial_inline(
	qaws_curve_inline* out,
	qaws_polynomial_desc const* desc);

#endif /* QAWS_INLINE_H */
