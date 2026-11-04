#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_INSPECT_H_F32
#define QAWS_INSPECT_H_F32
#define QAWS_INSPECT_H_BODY
#endif
#else
#ifndef QAWS_INSPECT_H_F64
#define QAWS_INSPECT_H_F64
#define QAWS_INSPECT_H_BODY
#endif
#endif
#ifdef QAWS_INSPECT_H_BODY
#undef QAWS_INSPECT_H_BODY
#define QAWS_INSPECT_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"

/* All inspection functions are thread-safe on immutable curves. */

/* Generic inspection */
qaws_curve_kind qaws_curve_get_kind(qaws_curve const* curve);
qaws_dimension  qaws_curve_get_dimension(qaws_curve const* curve);
unsigned int    qaws_curve_get_degree(qaws_curve const* curve);
unsigned int    qaws_curve_get_span_count(qaws_curve const* curve);
qaws_range      qaws_curve_get_parameter_range(qaws_curve const* curve);
int             qaws_curve_is_closed(qaws_curve const* curve);
int             qaws_curve_is_periodic(qaws_curve const* curve);
int             qaws_curve_is_rational(qaws_curve const* curve);
qaws_continuity qaws_curve_get_continuity(qaws_curve const* curve);

/* Analysis helpers */
qaws_status qaws_curve_compute_arc_length(
	qaws_scalar* out_length,
	qaws_curve const* curve,
	qaws_scalar parameter_min,
	qaws_scalar parameter_max);
/* The parameter at arc length `length` after parameter_min (safeguarded Newton on the measure above); OUT_OF_RANGE for
   a negative length or one past the end of the range. */
qaws_status qaws_curve_compute_parameter_at_arc_length(
	qaws_scalar* out_parameter,
	qaws_curve const* curve,
	qaws_scalar parameter_min,
	qaws_scalar length);

qaws_status qaws_curve_compute_bounds_2d(
	qaws_vec2* out_min,
	qaws_vec2* out_max,
	qaws_curve const* curve);

qaws_status qaws_curve_compute_bounds_3d(
	qaws_vec3* out_min,
	qaws_vec3* out_max,
	qaws_curve const* curve);

qaws_status qaws_curve_find_closest_parameter_2d(
	qaws_scalar* out_parameter,
	qaws_curve const* curve,
	qaws_vec2 point);

qaws_status qaws_curve_find_closest_parameter_3d(
	qaws_scalar* out_parameter,
	qaws_curve const* curve,
	qaws_vec3 point);

/* The closest point with a certificate, for curves with Bezier spans (Bezier, rational Bezier, B-spline, NURBS and borrowed
   splines; the two functions above use it too): branch and bound over the spans' Bezier forms with convex-hull and fat-line
   bounds. The global minimum distance lies in [*out_distance - *out_gap, *out_distance]; the gap is within 1e-9 of the
   curve's extent (QAWS_EPSILON in float) unless the work bound (4096 halvings) is reached, as for a point at a circle's
   center. QAWS_STATUS_UNSUPPORTED_OPERATION for the other families and for non-positive weights. */
qaws_status qaws_curve_find_closest_point_certified_2d(
	qaws_scalar* out_parameter,
	qaws_scalar* out_distance,
	qaws_scalar* out_gap,
	qaws_curve const* curve,
	qaws_vec2 point);

qaws_status qaws_curve_find_closest_point_certified_3d(
	qaws_scalar* out_parameter,
	qaws_scalar* out_distance,
	qaws_scalar* out_gap,
	qaws_curve const* curve,
	qaws_vec3 point);

qaws_status qaws_curve_get_span_continuity(
	qaws_continuity* out_continuity,
	qaws_curve const* curve,
	unsigned int boundary_index);

/* Derived geometric helpers */
qaws_status qaws_curve_compute_tangent_2d(
	qaws_vec2* out_tangent,
	qaws_curve const* curve,
	qaws_scalar parameter);

qaws_status qaws_curve_compute_tangent_3d(
	qaws_vec3* out_tangent,
	qaws_curve const* curve,
	qaws_scalar parameter);

qaws_status qaws_curve_compute_curvature_2d(
	qaws_scalar* out_curvature,
	qaws_curve const* curve,
	qaws_scalar parameter);

qaws_status qaws_curve_compute_curvature_3d(
	qaws_scalar* out_curvature,
	qaws_curve const* curve,
	qaws_scalar parameter);

qaws_status qaws_curve_compute_torsion_3d(
	qaws_scalar* out_torsion,
	qaws_curve const* curve,
	qaws_scalar parameter);

qaws_status qaws_curve_compute_speed(
	qaws_scalar* out_speed,
	qaws_curve const* curve,
	qaws_scalar parameter);

/* Frenet frame */
qaws_status qaws_curve_compute_normal_2d(
	qaws_vec2* out_normal,
	qaws_curve const* curve,
	qaws_scalar parameter);

qaws_status qaws_curve_compute_frenet_frame_3d(
	qaws_vec3* out_tangent,
	qaws_vec3* out_normal,
	qaws_vec3* out_binormal,
	qaws_curve const* curve,
	qaws_scalar parameter);

/* Rotation-minimizing (parallel-transport, Bishop) frames at count non-decreasing parameters: out_frames[9 i ..] = the unit tangent,
   the normal and the binormal (tangent x normal), row i for parameters[i]. Defined where the Frenet frame is not (straight pieces,
   inflections) and free of its twist: each normal follows from the previous by double reflection (core/qaws_frame_core.h; Wang et
   al. 2008), fourth-order accurate in the parameter step (test 147: halving the step divides the error by 16). The first normal is
   initial_normal made perpendicular to the tangent (NULL, or parallel to it: a perpendicular chosen from the tangent alone). A
   closed curve's frame need not close (its holonomy). QAWS_STATUS_DEGENERATE_CURVE where the speed vanishes. No allocation. */
qaws_status qaws_curve_compute_rotation_minimizing_frames_3d(
	qaws_scalar* out_frames,
	qaws_curve const* curve,
	qaws_scalar const* parameters,
	qaws_vec3 const* initial_normal,
	unsigned int count);

/* Inflection point detection */
qaws_status qaws_curve_find_inflection_points(
	qaws_scalar* out_parameters,
	unsigned int* out_count,
	qaws_curve const* curve,
	unsigned int parameter_capacity);

/* Extrema detection */
qaws_status qaws_curve_find_extrema(
	qaws_scalar* out_parameters,
	unsigned int* out_count,
	qaws_curve const* curve,
	unsigned int parameter_capacity,
	unsigned int axis);

/* Curvature comb data */
typedef struct qaws_curvature_sample_2d {
	qaws_vec2 position;
	qaws_scalar curvature;
	qaws_vec2 normal;
} qaws_curvature_sample_2d;

typedef struct qaws_curvature_sample_3d {
	qaws_vec3 position;
	qaws_scalar curvature;
	qaws_vec3 normal;
} qaws_curvature_sample_3d;

qaws_status qaws_curve_compute_curvature_comb_2d(
	qaws_curvature_sample_2d* out_samples,
	qaws_curve const* curve,
	unsigned int sample_count,
	unsigned int sample_capacity);

qaws_status qaws_curve_compute_curvature_comb_3d(
	qaws_curvature_sample_3d* out_samples,
	qaws_curve const* curve,
	unsigned int sample_count,
	unsigned int sample_capacity);

/* Winding number */
qaws_status qaws_curve_compute_winding_number_2d(
	int* out_winding_number,
	qaws_curve const* curve,
	qaws_vec2 point);

/* Intersection results */
typedef struct qaws_intersection_2d {
	qaws_scalar parameter_a;
	qaws_scalar parameter_b;
	qaws_vec2 position;
} qaws_intersection_2d;

typedef struct qaws_intersection_3d {
	qaws_scalar parameter_a;
	qaws_scalar parameter_b;
	qaws_vec3 position;
} qaws_intersection_3d;

/* Family-specific inspection */
qaws_status qaws_bezier_get_control_points(
	void* out_control_points,
	unsigned int* out_point_count,
	qaws_curve const* curve,
	unsigned int point_capacity);

/* Control points of a Bezier, rational Bezier, B-spline or NURBS curve
   (Cartesian, dimension scalars each; weights through qaws_nurbs_get_weights
   or the rational Bezier desc). Other kinds: UNSUPPORTED_OPERATION. */
qaws_status qaws_curve_get_control_points(
	void* out_control_points,
	unsigned int* out_point_count,
	qaws_curve const* curve,
	unsigned int point_capacity);

qaws_status qaws_bspline_get_knots(
	qaws_scalar* out_knots,
	unsigned int* out_knot_count,
	qaws_curve const* curve,
	unsigned int knot_capacity);

qaws_status qaws_nurbs_get_weights(
	qaws_scalar* out_weights,
	unsigned int* out_weight_count,
	qaws_curve const* curve,
	unsigned int weight_capacity);

/* Surface inspection */

qaws_status qaws_surface_compute_bounds(
	qaws_vec3* out_min,
	qaws_vec3* out_max,
	qaws_surface const* surface);

qaws_status qaws_surface_compute_area(
	qaws_scalar* out_area,
	qaws_surface const* surface);

qaws_status qaws_surface_compute_gaussian_curvature(
	qaws_scalar* out_curvature,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v);

qaws_status qaws_surface_compute_mean_curvature(
	qaws_scalar* out_curvature,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v);

typedef struct qaws_surface_curvature_result {
	qaws_scalar gaussian;
	qaws_scalar mean;
	qaws_scalar kappa1;
	qaws_scalar kappa2;
} qaws_surface_curvature_result;

qaws_status qaws_surface_compute_principal_curvatures(
	qaws_surface_curvature_result* out_result,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v);

/* Hausdorff distance */

qaws_status qaws_curve_compute_hausdorff_distance_2d(
	qaws_scalar* out_distance,
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	unsigned int sample_count);

qaws_status qaws_curve_compute_hausdorff_distance_3d(
	qaws_scalar* out_distance,
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	unsigned int sample_count);

/* Closest point on surface */

/* Find the (u,v) parameters of the point on a surface closest to a given 3D point.
   An 8x8 grid and Newton on the distance function gradient, then certified as
   below where the surface has Bezier patches. */
qaws_status qaws_surface_find_closest_point(
	qaws_scalar* out_u,
	qaws_scalar* out_v,
	qaws_vec3* out_closest_point,
	qaws_surface const* surface,
	qaws_vec3 point);

/* The closest point with a certificate, for Bezier, B-spline, NURBS and borrowed surfaces of degree at most 15: branch and
   bound over the Bezier patches with bounds by their nets' boxes, axis-aligned and in each patch's own frame. The global
   minimum distance lies in [*out_distance - *out_gap, *out_distance]; the gap is within 1e-9 of the extent (QAWS_EPSILON in
   float) unless the work bound (8192 splits) is reached. QAWS_STATUS_UNSUPPORTED_OPERATION for other surfaces. */
qaws_status qaws_surface_find_closest_point_certified(
	qaws_scalar* out_u,
	qaws_scalar* out_v,
	qaws_vec3* out_closest_point,
	qaws_scalar* out_distance,
	qaws_scalar* out_gap,
	qaws_surface const* surface,
	qaws_vec3 point);

/* Guaranteed ranges over a parameter box: for every (u, v) in [u0, u1] x [v0, v1], S(u, v) lies in
   [out_lo[0..2], out_hi[0..2]]; with order 1, S_u and S_v in out_lo / out_hi [3..5] and [6..8]; with order 2, S_uu, S_uv and
   S_vv in [9..11], [12..14], [15..17]. From the Bezier nets of the box - each knot span's B-spline window (its derivatives'
   windows for the derivatives) blossomed over the part of the box in the span - every operation rounded outward:
   Bernstein's hull bounds the homogeneous parts, the quotient rule in interval arithmetic the rational ones (S_u = (X_u -
   S W_u) / W ...). Wider than the true ranges, never narrower; shrinking with the box. Bezier, B-spline, NURBS (positive
   weights) and borrowed spline surfaces of degree at most 15; QAWS_STATUS_UNSUPPORTED_OPERATION for other surfaces,
   QAWS_STATUS_NUMERICAL_FAILURE where a weight's range reaches 0, QAWS_STATUS_INVALID_ARGUMENT for a box empty, not finite or
   past the domain, an order past 2. */
qaws_status qaws_surface_enclose(
	double* out_lo,
	double* out_hi,
	qaws_surface const* surface,
	double u0,
	double u1,
	double v0,
	double v1,
	unsigned int order);

/* Guaranteed ranges over a parameter interval: for every t in [t0, t1], C(t) lies in [out_lo[0 .. dim - 1], out_hi[...]];
   with order 1, C' in [dim .. 2 dim - 1]; with order 2, C'' in [2 dim .. 3 dim - 1]. As qaws_surface_enclose: each span's
   window - a Hermite or uniform Catmull-Rom span its Bezier points, in intervals - blossomed over the interval's part in
   the span, every operation rounded outward, the quotient rule for rational curves. Bezier, rational Bezier, B-spline,
   NURBS (positive weights), Hermite, uniform Catmull-Rom and borrowed spline curves of degree at most 15;
   QAWS_STATUS_UNSUPPORTED_OPERATION for other curves, QAWS_STATUS_NUMERICAL_FAILURE where a weight's range reaches 0,
   QAWS_STATUS_INVALID_ARGUMENT for an interval empty, not finite or past the range, an order past 2. */
qaws_status qaws_curve_enclose(
	double* out_lo,
	double* out_hi,
	qaws_curve const* curve,
	double t0,
	double t1,
	unsigned int order);

/* Curve-plane intersection */

/* Plane defined by a point and a normal vector. */
typedef struct qaws_plane {
	qaws_vec3 point;
	qaws_vec3 normal;
} qaws_plane;

/* Surface-curve intersection */

/* Find points where a 3D curve pierces a surface.
   Returns (u,v,t) parameter triples for each intersection.
   Uses Newton iteration on the system S(u,v) - C(t) = 0. */
typedef struct qaws_surface_curve_intersection {
	qaws_scalar u;
	qaws_scalar v;
	qaws_scalar t;
	qaws_vec3 position;
} qaws_surface_curve_intersection;

/* Adaptive tessellation */

/* Adaptive tessellation result: vertex + index buffers for triangle mesh.
   Caller allocates buffers. out_vertex_count and out_index_count receive
   the actual counts written.

   max_depth: maximum subdivision depth (0 = default 5).
   curvature_threshold: subdivide quads where curvature exceeds this (0 = default 0.1).
   max_edge_length: subdivide edges longer than this in 3D space (0 = disabled). */

typedef struct qaws_tessellation_desc {
	unsigned int max_depth;
	qaws_scalar curvature_threshold;
	qaws_scalar max_edge_length;
} qaws_tessellation_desc;

typedef struct qaws_tessellation_vertex {
	qaws_vec3 position;
	qaws_vec3 normal;
	qaws_scalar u;
	qaws_scalar v;
} qaws_tessellation_vertex;

qaws_status qaws_surface_tessellate(
	qaws_tessellation_vertex* out_vertices,
	unsigned int* out_vertex_count,
	unsigned int* out_indices,
	unsigned int* out_index_count,
	qaws_surface const* surface,
	unsigned int vertex_capacity,
	unsigned int index_capacity,
	qaws_tessellation_desc const* desc);

/* Adaptive tessellation by tolerances into caller buffers, without cracks and without allocation: a quadtree over (u, v)
   whose quads split until the surface at their 3 x 3 test points is within the chord tolerance of the bilinear patch of
   their corners and the normal turns by at most the angle tolerance across them. Deterministic: call with NULL buffers
   for the counts, then with buffers of those sizes (BUFFER_TOO_SMALL with the counts otherwise). Triangles are
   counter-clockwise in (u, v). Vertices are written per quad: a point shared by quads repeats, with the same bits. */
typedef struct qaws_surface_tessellation_desc
{
	qaws_scalar chord_tolerance;        /* distance; 0: unused */
	qaws_scalar angle_tolerance;        /* radians; 0: unused (one of the two is needed) */
	unsigned int min_depth;             /* quarterings everywhere first; 0: 2 */
	unsigned int max_depth;             /* quarterings at most; 0: 10, at most 20 */
} qaws_surface_tessellation_desc;

qaws_status qaws_surface_tessellate_adaptive(
	qaws_tessellation_vertex* out_vertices,
	unsigned int* out_vertex_count,
	unsigned int* out_indices,
	unsigned int* out_index_count,
	qaws_surface const* surface,
	unsigned int vertex_capacity,
	unsigned int index_capacity,
	qaws_surface_tessellation_desc const* desc);

/* Curve projection onto surface */

/* Project a 3D curve onto a surface, producing a 2D curve in (u,v) parameter space.
   Each point C(t) is mapped to the closest (u,v) on the surface.
   sample_count controls the projection resolution (0 = default 128).
   The output is a B-spline curve in 2D with x=u, y=v. */
qaws_status qaws_surface_project_curve(
	qaws_curve** out_uv_curve,
	qaws_surface const* surface,
	qaws_curve const* curve,
	unsigned int sample_count);

#endif /* QAWS_INSPECT_H */
