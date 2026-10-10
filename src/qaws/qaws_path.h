#ifndef QAWS_PATH_H
#define QAWS_PATH_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * 2D paths and regions of curves.
 *
 * A path is a list of 2D curves joined end to end: each curve ends where the
 * next one starts, and a closed path's last curve ends where its first one
 * starts. A polygon is a path of one degree-1 B-spline. A region is a list of
 * closed paths read with a fill rule, as in Clipper2: a point is inside when
 * the winding number of all the paths around it passes the rule.
 *
 * Thread-safe on immutable curves.
 */

typedef struct qaws_path_2d
{
	qaws_curve const* const* curves;
	unsigned int curve_count;
	int closed;
} qaws_path_2d;

typedef enum qaws_fill_rule
{
	QAWS_FILL_EVEN_ODD = 0,     /* odd winding numbers are inside */
	QAWS_FILL_NON_ZERO,         /* non-zero winding numbers are inside */
	QAWS_FILL_POSITIVE,         /* positive winding numbers are inside */
	QAWS_FILL_NEGATIVE          /* negative winding numbers are inside */
} qaws_fill_rule;

typedef enum qaws_point_location
{
	QAWS_POINT_OUTSIDE = 0,
	QAWS_POINT_INSIDE,
	QAWS_POINT_ON               /* within the tolerance of a boundary */
} qaws_point_location;

/* 1 when winding number w is inside under the rule. */
int qaws_fill_rule_inside(qaws_fill_rule rule, int winding);

/* ------------------------------------------------------------------------ */
/*  Construction                                                            */
/* ------------------------------------------------------------------------ */

/* A polyline as a degree-1 B-spline through the points (2 scalars each),
   with parameter i at point i. closed = 1 adds the edge back to the first
   point (do not repeat it). */
qaws_status qaws_curve_create_polyline_2d(
	qaws_scalar const* points,
	unsigned int point_count,
	int closed,
	qaws_curve** out_curve);

/* An exact ellipse as a closed rational quadratic NURBS (four quarter arcs),
   counter-clockwise from center + rx (cos rot, sin rot). rx == ry gives a
   circle. */
qaws_status qaws_curve_create_ellipse_2d(
	qaws_scalar center_x, qaws_scalar center_y,
	qaws_scalar radius_x, qaws_scalar radius_y,
	qaws_scalar rotation,
	qaws_curve** out_curve);

/* The curve under x' = m[0] x + m[1] y + m[2], y' = m[3] x + m[4] y + m[5].
   Bezier, rational Bezier, B-spline, NURBS and polynomial curves keep their
   kind; arcs and clothoids keep theirs under a similarity (rotation, uniform
   scale, reflection, translation). Under other maps an arc becomes an exact
   NURBS and a clothoid a cubic B-spline fitted to 1e-9 of its extent. Other
   kinds become their exact qaws_curve_extract form first. */
qaws_status qaws_curve_transform_2d(
	qaws_curve const* curve,
	qaws_scalar const m[6],
	qaws_curve** out_curve);

/* ------------------------------------------------------------------------ */
/*  Measures                                                                */
/* ------------------------------------------------------------------------ */

/* Signed area enclosed by a closed path, positive counter-clockwise
   (x right, y up): the integral of x dy around it. Gauss-Legendre per span,
   exact for polynomial spans, adaptive for the others. */
qaws_status qaws_path_compute_area_2d(qaws_path_2d const* path, qaws_scalar* out_area);

/* Sum of the signed areas of the paths of a region. */
qaws_status qaws_region_compute_area_2d(qaws_path_2d const* paths, unsigned int path_count, qaws_scalar* out_area);

/* 1 when the closed path's signed area is >= 0 (Clipper2's IsPositive). */
int qaws_path_is_positive_2d(qaws_path_2d const* path);

/* Arc length of a path. */
qaws_status qaws_path_compute_length_2d(qaws_path_2d const* path, qaws_scalar* out_length);

/* Tight bounds of a path: the end points and every point where a
   coordinate's derivative vanishes. */
qaws_status qaws_path_compute_bounds_2d(qaws_path_2d const* path, qaws_vec2* out_min, qaws_vec2* out_max);

/* Winding number of a closed path around a point; QAWS_STATUS_INVALID_ARGUMENT
   when the point lies on the path (within tolerance). */
qaws_status qaws_path_compute_winding_2d(qaws_path_2d const* path, qaws_vec2 point, int* out_winding);

/* Where a point lies in a region: ON within `tolerance` of a boundary
   (0 = 1e-10 of the region's extent, 1e-5 in float builds), else INSIDE or
   OUTSIDE by the fill rule on the total winding number. */
qaws_status qaws_region_locate_point_2d(
	qaws_path_2d const* paths,
	unsigned int path_count,
	qaws_fill_rule fill_rule,
	qaws_vec2 point,
	qaws_scalar tolerance,
	qaws_point_location* out_location);

/* ------------------------------------------------------------------------ */
/*  Polylines (Clipper2's path utilities on point arrays)                   */
/* ------------------------------------------------------------------------ */

/* Each takes `count` points (2 scalars each) and writes at most `count`
   points to `out` (which may equal `points`), with the kept count in
   *out_count. */

/* Consecutive repeated points removed; closed also drops a last point equal
   to the first. */
qaws_status qaws_polyline_strip_duplicates_2d(qaws_scalar const* points, unsigned int count, int closed,
	qaws_scalar* out, unsigned int* out_count);

/* Points on the line through their neighbours removed, and spikes (a point
   whose neighbours lie on the same side along the line) too, by an exact
   orientation test. Fewer than 3 points left (2 for open) gives 0. */
qaws_status qaws_polyline_trim_collinear_2d(qaws_scalar const* points, unsigned int count, int closed,
	qaws_scalar* out, unsigned int* out_count);

/* Clipper2's SimplifyPath: repeatedly removes the point nearest the line
   through its neighbours while that distance is at most epsilon. */
qaws_status qaws_polyline_simplify_2d(qaws_scalar const* points, unsigned int count, int closed, qaws_scalar epsilon,
	qaws_scalar* out, unsigned int* out_count);

/* Ramer-Douglas-Peucker: keeps the ends, recursively keeps the point
   farthest from the chord while it is farther than epsilon. */
qaws_status qaws_polyline_rdp_2d(qaws_scalar const* points, unsigned int count, qaws_scalar epsilon,
	qaws_scalar* out, unsigned int* out_count);

#endif /* QAWS_PATH_H */
