#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_CLOTHOID_H_F32
#define QAWS_CLOTHOID_H_F32
#define QAWS_CLOTHOID_H_BODY
#endif
#else
#ifndef QAWS_CLOTHOID_H_F64
#define QAWS_CLOTHOID_H_F64
#define QAWS_CLOTHOID_H_BODY
#endif
#endif
#ifdef QAWS_CLOTHOID_H_BODY
#undef QAWS_CLOTHOID_H_BODY
#define QAWS_CLOTHOID_H

#include "qaws_types.h"
#include "qaws_status.h"

/* Clothoid (Euler spiral / Cornu spiral) defined by:
 *   - origin: starting point
 *   - start_angle: initial heading (radians)
 *   - start_curvature: curvature at s=0 (kappa_0)
 *   - end_curvature: curvature at s=L (kappa_1)
 *   - length: total arc length L
 *
 * The curvature varies linearly with arc length:
 *   kappa(s) = kappa_0 + (kappa_1 - kappa_0) * s / L
 *
 * Position is computed by integrating:
 *   x(s) = origin.x + integral_0^s cos(theta(u)) du
 *   y(s) = origin.y + integral_0^s sin(theta(u)) du
 * where theta(s) = start_angle + kappa_0*s + (kappa_1-kappa_0)*s^2/(2*L)
 *
 * 2D only (3D clothoids are not standard). */
typedef struct qaws_clothoid_desc
{
	qaws_scalar origin_x;
	qaws_scalar origin_y;
	qaws_scalar start_angle;        /* radians */
	qaws_scalar start_curvature;    /* kappa_0 */
	qaws_scalar end_curvature;      /* kappa_1 */
	qaws_scalar length;             /* arc length L > 0 */
} qaws_clothoid_desc;

qaws_status qaws_curve_create_clothoid(
	qaws_curve** out_curve,
	qaws_clothoid_desc const* desc);

#endif /* QAWS_CLOTHOID_H */
