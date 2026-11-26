#ifndef QAWS_INTERNAL_SURFACE_H
#define QAWS_INTERNAL_SURFACE_H

#include "../qaws_platform.h"
#include "../qaws_surface_types.h"

typedef struct qaws_surface_vtable qaws_surface_vtable;
typedef struct qaws_surface_diff_vtable qaws_surface_diff_vtable;

struct qaws_surface
{
	qaws_surface_kind kind;
	unsigned int u_degree;
	unsigned int v_degree;
	qaws_range u_range;
	qaws_range v_range;
	qaws_surface_vtable const* vtable;
	void* impl;
	qaws_allocator const* allocator; /* NULL = use malloc/free */
};

struct qaws_surface_vtable
{
	qaws_status (*eval_3d)(
		qaws_surface const* surface,
		qaws_scalar u,
		qaws_scalar v,
		unsigned int eval_flags,
		qaws_surface_eval_result* out_result);

	void (*destroy_impl)(void* impl, qaws_allocator const* allocator);

	int (*is_rational)(qaws_surface const* surface);

	/* Differential rules, NULL when the family is not differentiable yet. */
	qaws_surface_diff_vtable const* diff;
};

/* Bezier surface impl */
typedef struct qaws_surface_bezier_impl
{
	qaws_scalar* control_points;  /* 3 scalars per point, row-major */
	unsigned int u_count;
	unsigned int v_count;
} qaws_surface_bezier_impl;

/* B-spline surface impl */
typedef struct qaws_surface_bspline_impl
{
	qaws_scalar* control_points;
	unsigned int u_count;
	unsigned int v_count;
	qaws_scalar* u_knots;
	unsigned int u_knot_count;
	qaws_scalar* v_knots;
	unsigned int v_knot_count;
} qaws_surface_bspline_impl;

/* NURBS surface impl */
typedef struct qaws_surface_nurbs_impl
{
	qaws_scalar* control_points;
	unsigned int u_count;
	unsigned int v_count;
	qaws_scalar* weights;
	qaws_scalar* u_knots;
	unsigned int u_knot_count;
	qaws_scalar* v_knots;
	unsigned int v_knot_count;
} qaws_surface_nurbs_impl;

/* Surface of revolution impl: S(u, v) = origin + r(t) (cos(theta) X + sin(theta) Y) + h(t) Z, theta = u angle,
   t = prof_range.min + v (prof_range.max - prof_range.min), (r, h) the 2D profile */
typedef struct qaws_surface_revolution_impl
{
	qaws_curve const* profile;
	qaws_range prof_range;
	qaws_vec3 origin;
	qaws_vec3 X_axis;   /* perpendicular to rotation axis */
	qaws_vec3 Y_axis;   /* perpendicular to rotation axis and X_axis */
	qaws_vec3 Z_axis;   /* normalized rotation axis direction */
	qaws_scalar angle;
} qaws_surface_revolution_impl;

qaws_surface* qaws_internal_surface_alloc(
	qaws_surface_kind kind,
	unsigned int u_degree,
	unsigned int v_degree,
	qaws_range u_range,
	qaws_range v_range,
	qaws_surface_vtable const* vtable);

qaws_surface* qaws_internal_surface_alloc_ex(
	qaws_surface_kind kind,
	unsigned int u_degree,
	unsigned int v_degree,
	qaws_range u_range,
	qaws_range v_range,
	qaws_surface_vtable const* vtable,
	qaws_allocator const* allocator);

void qaws_internal_surface_free(qaws_surface* surface);

/* Generate uniform clamped knot vector into out_knots.
   knot_count = num_cp + degree + 1. Returns knot_count. */
unsigned int qaws_internal_surface_uniform_knots(
	qaws_scalar* out_knots,
	unsigned int capacity,
	unsigned int degree,
	unsigned int num_cp);

/* Normalize the cross product du x dv into out.
   Falls back to +Z when the cross product is degenerate. */
QAWS_INLINE void qaws_internal_surface_normal(qaws_vec3* out, qaws_vec3 du, qaws_vec3 dv)
{
	qaws_scalar nx = du.y * dv.z - du.z * dv.y;
	qaws_scalar ny = du.z * dv.x - du.x * dv.z;
	qaws_scalar nz = du.x * dv.y - du.y * dv.x;
	qaws_scalar len = QAWS_SQRT(nx * nx + ny * ny + nz * nz);
	if (len > QAWS_LITERAL(1e-12))
	{
		out->x = nx / len; out->y = ny / len; out->z = nz / len;
	}
	else
	{
		out->x = 0; out->y = 0; out->z = 1;
	}
}

/* Rows of qaws_internal_catmull_rom_weights: derivative orders 0..5. */
#define QAWS_INTERNAL_CR_ROWS 6

/* The Catmull-Rom blend as explicit weights: blend(t) = sum_k out_w[r][k] *
   pts[out_index[k]] for the r-th derivative with respect to t (indices may
   repeat at the ends). Matches qaws_internal_surface_catmull_rom_blend. */
unsigned int qaws_internal_catmull_rom_weights(
	unsigned int* out_index,
	qaws_scalar const* params,
	unsigned int n_pts,
	qaws_scalar t,
	qaws_scalar (*out_w)[4]);

/* Non-uniform Catmull-Rom blend through n_pts points with the given
   parameter values. out_deriv may be NULL. */
void qaws_internal_surface_catmull_rom_blend(
	qaws_vec3* out_pos,
	qaws_vec3* out_deriv,
	qaws_vec3 const* pts,
	qaws_scalar const* params,
	unsigned int n_pts,
	qaws_scalar t);

/* A trimmed surface's base (NULL when the surface is not trimmed) and its trim loops: loop index's curves (2D, in the
   base's (u, v)), their count and whether it is an outer loop (qaws_shatter.c) */
qaws_surface const* qaws_internal_trim_base(qaws_surface const* surface);
unsigned int qaws_internal_trim_loop_count(qaws_surface const* surface);
struct qaws_curve;
struct qaws_curve const* const* qaws_internal_trim_loop(unsigned int* out_count, int* out_outer, qaws_surface const* surface, unsigned int index);

/* qaws_surface_evaluate_batch on B-spline and NURBS surfaces of degrees 1 to 5: the u and v basis functions of a SIMD
   register's worth of points at once (qaws_internal_bspline_basis_derivs_lanes), each result as qaws_surface_evaluate
   makes it, bit for bit (the first failing status); QAWS_STATUS_UNSUPPORTED_OPERATION for the other surfaces and degrees
   and without SIMD (the caller loops). */
qaws_status qaws_internal_bspline_surface_eval_batch(
	qaws_surface_eval_result* out_results,
	qaws_surface const* surface,
	qaws_scalar const* u_params,
	qaws_scalar const* v_params,
	unsigned int count,
	unsigned int eval_flags);
qaws_status qaws_internal_nurbs_surface_eval_batch(
	qaws_surface_eval_result* out_results,
	qaws_surface const* surface,
	qaws_scalar const* u_params,
	qaws_scalar const* v_params,
	unsigned int count,
	unsigned int eval_flags);

#endif /* QAWS_INTERNAL_SURFACE_H */
