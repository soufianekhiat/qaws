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
	unsigned int degree,
	unsigned int num_cp,
	qaws_scalar* out_knots,
	unsigned int capacity);

/* Normalize the cross product du x dv into out.
   Falls back to +Z when the cross product is degenerate. */
QAWS_INLINE void qaws_internal_surface_normal(qaws_vec3 du, qaws_vec3 dv, qaws_vec3* out)
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
	qaws_scalar const* params,
	unsigned int n_pts,
	qaws_scalar t,
	unsigned int* out_index,
	qaws_scalar (*out_w)[4]);

/* Non-uniform Catmull-Rom blend through n_pts points with the given
   parameter values. out_deriv may be NULL. */
void qaws_internal_surface_catmull_rom_blend(
	qaws_vec3 const* pts,
	qaws_scalar const* params,
	unsigned int n_pts,
	qaws_scalar t,
	qaws_vec3* out_pos,
	qaws_vec3* out_deriv);

#endif /* QAWS_INTERNAL_SURFACE_H */
