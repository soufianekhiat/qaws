#ifndef QAWS_INTERNAL_TYPES_H
#define QAWS_INTERNAL_TYPES_H

#include "../qaws_types.h"
#include "../qaws_status.h"

typedef struct qaws_curve_vtable qaws_curve_vtable;
typedef struct qaws_curve_diff_vtable qaws_curve_diff_vtable;

struct qaws_curve
{
	qaws_curve_kind kind;
	qaws_dimension dimension;
	unsigned int degree;
	unsigned int span_count;
	qaws_range parameter_range;
	qaws_scalar* span_boundaries;
	qaws_curve_vtable const* vtable;
	void* impl;
	qaws_allocator const* allocator; /* NULL = use malloc/free */
};

struct qaws_curve_vtable
{
	qaws_status (*eval_span_2d)(
		qaws_curve const* curve,
		unsigned int span_index,
		qaws_scalar local_t,
		unsigned int eval_flags,
		qaws_eval_result_2d* out_result);

	qaws_status (*eval_span_3d)(
		qaws_curve const* curve,
		unsigned int span_index,
		qaws_scalar local_t,
		unsigned int eval_flags,
		qaws_eval_result_3d* out_result);

	void (*destroy_impl)(void* impl, qaws_allocator const* allocator);

	int (*is_closed)(qaws_curve const* curve);
	int (*is_periodic)(qaws_curve const* curve);
	int (*is_rational)(qaws_curve const* curve);
	qaws_continuity (*get_continuity)(qaws_curve const* curve);

	/* Differential rules, NULL when the family is not differentiable yet. */
	qaws_curve_diff_vtable const* diff;
};

struct qaws_traversal
{
	qaws_curve const* curve;
	qaws_traversal_desc desc;
	qaws_scalar total_arc_length;
	unsigned int table_size;
	qaws_scalar* table_params;
	qaws_scalar* table_distances;
	qaws_scalar current_distance;

	/* Multi-curve chain (NULL for single-curve) */
	qaws_curve const** chain_curves;
	unsigned int chain_count;
	qaws_scalar* chain_param_offsets; /* chain_count + 1: cumulative param offsets */

	/* Custom speed function time-to-distance table */
	qaws_scalar* custom_time_table;    /* time values */
	qaws_scalar* custom_dist_table;    /* corresponding distance values */
	unsigned int custom_table_size;    /* 0 if not used */
};

/* Family-specific impl structs */

typedef struct qaws_bezier_impl
{
	qaws_scalar* control_points;
	unsigned int control_point_count;
} qaws_bezier_impl;

typedef struct qaws_hermite_impl
{
	qaws_scalar* points;
	qaws_scalar* tangents;
	unsigned int point_count;
	qaws_scalar* span_coeffs; /* dim_count * 4 scalars per span: a,b,c,d */
} qaws_hermite_impl;

typedef struct qaws_catmull_rom_impl
{
	qaws_scalar* control_points;
	unsigned int control_point_count;
	qaws_parameterization parameterization;
	int closed;
	qaws_scalar* knot_params;
	qaws_scalar* segment_coeffs;
} qaws_catmull_rom_impl;

typedef struct qaws_bspline_impl
{
	qaws_scalar* control_points;
	unsigned int control_point_count;
	qaws_scalar* knots;
	unsigned int knot_count;
} qaws_bspline_impl;

typedef struct qaws_nurbs_impl
{
	qaws_scalar* control_points;
	unsigned int control_point_count;
	qaws_scalar* knots;
	unsigned int knot_count;
	qaws_scalar* weights;
	unsigned int weight_count;
} qaws_nurbs_impl;

typedef struct qaws_trajectory_impl
{
	qaws_scalar* key_positions;
	unsigned int key_count;
	qaws_scalar* key_times;
	qaws_scalar* key_velocities;
	unsigned int key_velocity_count;
	qaws_scalar* key_accelerations;
	unsigned int key_acceleration_count;
	int closed;
	qaws_scalar* span_coeffs; /* dim_count * 4 scalars per span: a,b,c,d (in normalized [0,1] local_t) */
} qaws_trajectory_impl;

typedef struct qaws_yuksel_subcurve
{
	/* Quadratic Bezier sub-curve: p0, p1, p2, and the optimal t parameter */
	qaws_scalar p0[3];
	qaws_scalar p1[3];
	qaws_scalar p2[3];
	qaws_scalar t1;
	/* Circular/elliptical: center, radius, angle limits */
	qaws_scalar center[3];
	qaws_scalar radius;
	qaws_scalar angle_start;
	qaws_scalar angle_end;
	qaws_scalar axis_u[3];
	qaws_scalar axis_v[3];
	/* For elliptical: second radius */
	qaws_scalar radius_b;
	int use_arc;
} qaws_yuksel_subcurve;

typedef struct qaws_yuksel_impl
{
	qaws_scalar* control_points;
	unsigned int control_point_count;
	int closed;
	int mode;
	qaws_yuksel_subcurve* subcurves;
} qaws_yuksel_impl;

typedef struct qaws_rational_bezier_impl
{
	qaws_scalar* control_points;
	unsigned int control_point_count;
	qaws_scalar* weights;
	qaws_scalar* weighted_points; /* homogeneous: (w*x, w*y, [w*z], w) per point */
} qaws_rational_bezier_impl;

#endif /* QAWS_INTERNAL_TYPES_H */
