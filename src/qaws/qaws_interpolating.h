#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_INTERPOLATING_H_F32
#define QAWS_INTERPOLATING_H_F32
#define QAWS_INTERPOLATING_H_BODY
#endif
#else
#ifndef QAWS_INTERPOLATING_H_F64
#define QAWS_INTERPOLATING_H_F64
#define QAWS_INTERPOLATING_H_BODY
#endif
#endif
#ifdef QAWS_INTERPOLATING_H_BODY
#undef QAWS_INTERPOLATING_H_BODY
#define QAWS_INTERPOLATING_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Smooth 2D splines solved from points alone: the tangents (and, for the
 * curvature-matching methods, the shapes) come from a global solve over all
 * the points, so moving one point moves the whole curve a little.
 *
 * The result is a composite of Bezier segments, its parameter domain
 * [0, segment_count]. Consecutive duplicate points are dropped first.
 *
 *   CUBIC_C2  Cubic segments, equal parameter intervals, C2 at every point;
 *             an open curve has zero second derivative at its ends. Point i
 *             is at parameter i.
 *   HOBBY     Cubic segments with Hobby's tangent angles and velocities: G1
 *             at every point, close to circular arcs where the points allow
 *             it. Point i is at parameter i.
 *   KAPPA     Quadratic segments, one per interior point, G2 where they join
 *             (G1 where the curve inflects at a join: there only the size of
 *             the curvature matches), each with its curvature maximum at its
 *             point. Segments join
 *             between points, so a point is inside a segment, not at an
 *             integer parameter; an open curve needs at least 3 points.
 *   CUBIC_G2  Cubic segments with the tangent directions and lengths solved
 *             so the curvature matches at every point: G2. Point i is at
 *             parameter i.
 *   S9        Segments whose curvature is a function of arc length, set by
 *             their end angles, the angles solved so the curvature matches at
 *             every point: rounder than a single cubic allows, with low
 *             curvature variation. The result is those segments as cubics,
 *             so G2 holds to the accuracy of that approximation, and a point
 *             is at a join between cubics, several per point: not at an
 *             integer parameter matching its index. A solve uses about
 *             145 KB of stack.
 *
 * An open curve of two distinct points is the straight segment for every
 * method. A closed curve needs at least 3 distinct points.
 *
 * Returns INVALID_ARGUMENT for a NULL or non-finite input or an unknown
 * method, DEGENERATE_CURVE when too few distinct points remain, and
 * NUMERICAL_FAILURE when the solve gives no finite result.
 */

#ifndef QAWS_ENUM_QAWS_INTERPOLATING_METHOD
#define QAWS_ENUM_QAWS_INTERPOLATING_METHOD
typedef enum qaws_interpolating_method
{
	QAWS_INTERPOLATING_CUBIC_C2 = 0,
	QAWS_INTERPOLATING_HOBBY,
	QAWS_INTERPOLATING_KAPPA,
	QAWS_INTERPOLATING_CUBIC_G2,
	QAWS_INTERPOLATING_S9
} qaws_interpolating_method;
#endif

typedef struct qaws_interpolating_spline_desc
{
	qaws_vec2 const* points;
	unsigned int point_count;
	qaws_interpolating_method method;
	int closed;                        /* non-zero: the last point joins the first */
} qaws_interpolating_spline_desc;

qaws_status qaws_curve_create_interpolating_spline(
	qaws_curve** out_curve,
	qaws_interpolating_spline_desc const* desc);

#endif /* QAWS_INTERPOLATING_H */
