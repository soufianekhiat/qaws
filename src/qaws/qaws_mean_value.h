#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MEAN_VALUE_H_F32
#define QAWS_MEAN_VALUE_H_F32
#define QAWS_MEAN_VALUE_H_BODY
#endif
#else
#ifndef QAWS_MEAN_VALUE_H_F64
#define QAWS_MEAN_VALUE_H_F64
#define QAWS_MEAN_VALUE_H_BODY
#endif
#endif
#ifdef QAWS_MEAN_VALUE_H_BODY
#undef QAWS_MEAN_VALUE_H_BODY
#define QAWS_MEAN_VALUE_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_batch_executor.h"

/*
 * Mean value coordinates over regions bounded by curves of any kind
 * (concave, with holes): generalized barycentric coordinates whose cage is
 * the boundary itself, its control points the handles.
 *
 * The transfinite mean value interpolant (Ju, Schaefer and Warren 2005;
 * Dyken and Floater 2009) of a function f given on the boundary is
 *
 *     g(x) = integral f(c) k ds / integral k ds,
 *     k = (c - x) . n / |c - x|^3,
 *
 * over every boundary loop, n the outward normal. It is the average of f
 * over the directions round x, each boundary point seen weighted by
 * 1 / distance; along a ray crossing the boundary several times the
 * crossings alternate in sign, so concave regions and holes need nothing
 * more. g interpolates f on the boundary and reproduces linear functions.
 *
 * Each boundary curve is a sum of its control points times its basis
 * functions (Bernstein, B-spline, rational for weighted kinds), so the
 * handles' coordinates are the integrals of their basis functions:
 *
 *     lambda_i(x) = integral B_i k ds / integral k ds,
 *
 * summing to 1, with sum lambda_i P_i = x. On a polygon they are Floater's
 * mean value coordinates. Moving the handles moves the boundary curves, and
 * x -> sum lambda_i(x) P_i' carries the inside along: a cage deformation
 * whose cage is the region's own curves.
 */

typedef struct qaws_mean_value qaws_mean_value;

typedef struct qaws_mean_value_desc
{
	qaws_path_2d const* paths;      /* closed loops, simple and not crossing one another; nested
	                                   loops alternate between outside and hole, whatever their
	                                   orientation */
	unsigned int path_count;
	qaws_scalar tolerance;          /* relative accuracy of the integrals; 0 = 1e-10 (1e-5 in float) */
	qaws_batch_executor const* executor;   /* optional: the many-point calls in parallel */
} qaws_mean_value_desc;

qaws_status qaws_mean_value_create(qaws_mean_value_desc const* desc, qaws_mean_value** out_coords);
void qaws_mean_value_destroy(qaws_mean_value* coords);

/* The handles: the boundary curves' control points, in loop and curve
   order, a point shared by consecutive curves (or a loop's two ends) taken
   once. Lines and polylines have their vertices; Beziers, B-splines and
   NURBS their control points; arcs the control points of rational
   quadratic pieces of at most 90 degrees; other kinds those of a cubic
   B-spline fitted to them (exact for cubic Hermite and Catmull-Rom). */
unsigned int qaws_mean_value_get_handle_count(qaws_mean_value const* coords);
qaws_status qaws_mean_value_get_handle(qaws_mean_value const* coords, unsigned int index, qaws_vec2* out_point);

/* The coordinates of a point: handle_count weights. Inside the region they
   are smooth; on the boundary they are its basis functions there. Outside,
   where the denominator can vanish, NUMERICAL_FAILURE when it does.
   Thread-safe. */
qaws_status qaws_mean_value_coordinates(qaws_mean_value const* coords, qaws_vec2 point, qaws_scalar* out_weights);

/* The cage deformation: the points mapped by sum lambda_i(x) P_i', with
   handle_count moved handles (the handles themselves give the points
   back). */
qaws_status qaws_mean_value_deform(qaws_mean_value const* coords, qaws_vec2 const* moved_handles,
	qaws_vec2 const* points, unsigned int count, qaws_vec2* out_points);

/* The transfinite interpolant of a function given on the boundary: dim
   values f(path, curve, parameter, point) per boundary point, interpolated
   at point. f is called from the calling thread only. */
typedef void (*qaws_mean_value_boundary_fn)(void* user, unsigned int path, unsigned int curve,
	qaws_scalar parameter, qaws_vec2 point, qaws_scalar* out_values);
qaws_status qaws_mean_value_interpolate(qaws_mean_value const* coords, qaws_vec2 point,
	qaws_mean_value_boundary_fn fn, void* user, unsigned int dim, qaws_scalar* out_values);

#endif /* QAWS_MEAN_VALUE_H */
