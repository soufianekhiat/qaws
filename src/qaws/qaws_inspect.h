#ifndef QAWS_INSPECT_H
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
	qaws_curve const* curve,
	qaws_scalar parameter_min,
	qaws_scalar parameter_max,
	qaws_scalar* out_length);

qaws_status qaws_curve_compute_bounds_2d(
	qaws_curve const* curve,
	qaws_vec2* out_min,
	qaws_vec2* out_max);

qaws_status qaws_curve_compute_bounds_3d(
	qaws_curve const* curve,
	qaws_vec3* out_min,
	qaws_vec3* out_max);

qaws_status qaws_curve_find_closest_parameter_2d(
	qaws_curve const* curve,
	qaws_vec2 point,
	qaws_scalar* out_parameter);

qaws_status qaws_curve_find_closest_parameter_3d(
	qaws_curve const* curve,
	qaws_vec3 point,
	qaws_scalar* out_parameter);

qaws_status qaws_curve_get_span_continuity(
	qaws_curve const* curve,
	unsigned int boundary_index,
	qaws_continuity* out_continuity);

/* Derived geometric helpers */
qaws_status qaws_curve_compute_tangent_2d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec2* out_tangent);

qaws_status qaws_curve_compute_tangent_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec3* out_tangent);

qaws_status qaws_curve_compute_curvature_2d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_curvature);

qaws_status qaws_curve_compute_curvature_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_curvature);

qaws_status qaws_curve_compute_torsion_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_torsion);

qaws_status qaws_curve_compute_speed(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_speed);

/* Frenet frame */
qaws_status qaws_curve_compute_normal_2d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec2* out_normal);

qaws_status qaws_curve_compute_frenet_frame_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec3* out_tangent,
	qaws_vec3* out_normal,
	qaws_vec3* out_binormal);

/* Inflection point detection */
qaws_status qaws_curve_find_inflection_points(
	qaws_curve const* curve,
	qaws_scalar* out_parameters,
	unsigned int parameter_capacity,
	unsigned int* out_count);

/* Extrema detection */
qaws_status qaws_curve_find_extrema(
	qaws_curve const* curve,
	unsigned int axis,
	qaws_scalar* out_parameters,
	unsigned int parameter_capacity,
	unsigned int* out_count);

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
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_curvature_sample_2d* out_samples,
	unsigned int sample_capacity);

qaws_status qaws_curve_compute_curvature_comb_3d(
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_curvature_sample_3d* out_samples,
	unsigned int sample_capacity);

/* Winding number */
qaws_status qaws_curve_compute_winding_number_2d(
	qaws_curve const* curve,
	qaws_vec2 point,
	int* out_winding_number);

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

/* Self-intersection detection */
qaws_status qaws_curve_find_self_intersections_2d(
	qaws_curve const* curve,
	qaws_intersection_2d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

qaws_status qaws_curve_find_self_intersections_3d(
	qaws_curve const* curve,
	qaws_intersection_3d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

/* Curve-curve intersection */
qaws_status qaws_curve_find_intersections_2d(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_intersection_2d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

qaws_status qaws_curve_find_intersections_3d(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_intersection_3d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

/* Family-specific inspection */
qaws_status qaws_bezier_get_control_points(
	qaws_curve const* curve,
	void* out_control_points,
	unsigned int point_capacity,
	unsigned int* out_point_count);

qaws_status qaws_bspline_get_knots(
	qaws_curve const* curve,
	qaws_scalar* out_knots,
	unsigned int knot_capacity,
	unsigned int* out_knot_count);

qaws_status qaws_nurbs_get_weights(
	qaws_curve const* curve,
	qaws_scalar* out_weights,
	unsigned int weight_capacity,
	unsigned int* out_weight_count);

/* Surface inspection */

qaws_status qaws_surface_compute_bounds(
	qaws_surface const* surface,
	qaws_vec3* out_min,
	qaws_vec3* out_max);

qaws_status qaws_surface_compute_area(
	qaws_surface const* surface,
	qaws_scalar* out_area);

qaws_status qaws_surface_compute_gaussian_curvature(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_scalar* out_curvature);

qaws_status qaws_surface_compute_mean_curvature(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_scalar* out_curvature);

typedef struct qaws_surface_curvature_result {
	qaws_scalar gaussian;
	qaws_scalar mean;
	qaws_scalar kappa1;
	qaws_scalar kappa2;
} qaws_surface_curvature_result;

qaws_status qaws_surface_compute_principal_curvatures(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_surface_curvature_result* out_result);

/* Hausdorff distance */

qaws_status qaws_curve_compute_hausdorff_distance_2d(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	unsigned int sample_count,
	qaws_scalar* out_distance);

qaws_status qaws_curve_compute_hausdorff_distance_3d(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	unsigned int sample_count,
	qaws_scalar* out_distance);

/* Closest point on surface */

/* Find the (u,v) parameters of the point on a surface closest to a given 3D point.
   Uses Newton iteration on the distance function gradient. */
qaws_status qaws_surface_find_closest_point(
	qaws_surface const* surface,
	qaws_vec3 point,
	qaws_scalar* out_u,
	qaws_scalar* out_v,
	qaws_vec3* out_closest_point);

/* Curve-plane intersection */

/* Plane defined by a point and a normal vector. */
typedef struct qaws_plane {
	qaws_vec3 point;
	qaws_vec3 normal;
} qaws_plane;

/* Find parameter values where a 3D curve crosses a plane.
   Uses bisection root-finding on the signed distance function dot(C(t) - plane.point, plane.normal). */
qaws_status qaws_curve_find_plane_intersections(
	qaws_curve const* curve,
	qaws_plane const* plane,
	qaws_scalar* out_parameters,
	qaws_vec3* out_positions,
	unsigned int capacity,
	unsigned int* out_count);

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

qaws_status qaws_surface_find_curve_intersections(
	qaws_surface const* surface,
	qaws_curve const* curve,
	qaws_surface_curve_intersection* out_intersections,
	unsigned int capacity,
	unsigned int* out_count);

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
	qaws_surface const* surface,
	qaws_tessellation_desc const* desc,
	qaws_tessellation_vertex* out_vertices,
	unsigned int vertex_capacity,
	unsigned int* out_vertex_count,
	unsigned int* out_indices,
	unsigned int index_capacity,
	unsigned int* out_index_count);

/* Curve projection onto surface */

/* Project a 3D curve onto a surface, producing a 2D curve in (u,v) parameter space.
   Each point C(t) is mapped to the closest (u,v) on the surface.
   sample_count controls the projection resolution (0 = default 128).
   The output is a B-spline curve in 2D with x=u, y=v. */
qaws_status qaws_surface_project_curve(
	qaws_surface const* surface,
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_curve** out_uv_curve);

/* Geodesic curve on surface */

/* Compute a geodesic curve (shortest path) between two points on a surface.
   Uses a shooting method with RK4 integration of the geodesic ODE.
   max_iterations controls the shooting refinement (0 = default 20).
   step_count controls the integration resolution (0 = default 200).
   The output is a 3D B-spline curve lying on the surface. */
qaws_status qaws_surface_compute_geodesic(
	qaws_surface const* surface,
	qaws_scalar start_u,
	qaws_scalar start_v,
	qaws_scalar end_u,
	qaws_scalar end_v,
	unsigned int max_iterations,
	unsigned int step_count,
	qaws_curve** out_curve);

#endif /* QAWS_INSPECT_H */
