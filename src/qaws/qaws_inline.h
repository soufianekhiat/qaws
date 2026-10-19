#ifndef QAWS_INLINE_H
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
 * and inspection functions via qaws_curve_inline_get_curve(&inline_curve).
 *
 * Constraints:
 *   - Single span only (degree+1 control points for Bezier)
 *   - Maximum degree 7 (8 CPs x 3 dimensions = 24 scalars for CPs)
 *   - 2D or 3D
 *   - Bezier and Polynomial families only
 *
 * The qaws_curve_inline is self-contained: it may be copied with memcpy or
 * assigned, and the copy used through qaws_curve_inline_get_curve(&copy),
 * which re-points it at its own storage. Do NOT call qaws_curve_destroy on it.
 *
 * Example:
 *   qaws_curve_inline ic;
 *   qaws_vec3 pts[] = {{0,0,0}, {1,2,3}, {4,5,6}, {7,8,9}};
 *   qaws_bezier_desc desc = { QAWS_DIMENSION_3D, 3, pts, 4 };
 *   qaws_curve_init_bezier_inline(&desc, &ic);
 *   // qaws_curve_evaluate_3d(qaws_curve_inline_get_curve(&ic), ...)
 */

qaws_status qaws_curve_init_bezier_inline(
	qaws_bezier_desc const* desc,
	qaws_curve_inline* out);

qaws_status qaws_curve_init_polynomial_inline(
	qaws_polynomial_desc const* desc,
	qaws_curve_inline* out);

#endif /* QAWS_INLINE_H */
