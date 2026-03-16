#ifndef QAWS_SURFACE_TSPLINE_H
#define QAWS_SURFACE_TSPLINE_H

#include "qaws_surface_types.h"

/* T-spline surface: locally refinable surface using a T-mesh.

   Unlike tensor-product B-splines where all control points share global
   knot vectors, T-splines assign local knot vectors to each control point.
   This allows local refinement: adding detail in one area without inserting
   entire rows/columns of control points.

   Each control point has:
   - A 3D position
   - An optional weight (1.0 = non-rational)
   - A local u-knot vector of length (2*u_degree+2)
   - A local v-knot vector of length (2*v_degree+2)

   The basis function for control point i is the product of two
   univariate B-spline basis functions defined by its local knot vectors:
   B_i(u,v) = N[u_knots_i](u) * N[v_knots_i](v)

   S(u,v) = sum_i w_i * P_i * B_i(u,v) / sum_i w_i * B_i(u,v)

   For non-rational T-splines (all weights = 1), this simplifies to:
   S(u,v) = sum_i P_i * B_i(u,v) / sum_i B_i(u,v)

   The local knot vectors are inferred from the T-mesh connectivity,
   but in this API they are provided explicitly by the user. */

typedef struct qaws_tspline_control_point
{
	qaws_vec3 position;
	qaws_scalar weight;         /* 1.0 for non-rational */
	qaws_scalar const* u_knots; /* local u-knot vector, length = 2*u_degree+2 */
	qaws_scalar const* v_knots; /* local v-knot vector, length = 2*v_degree+2 */
} qaws_tspline_control_point;

typedef struct qaws_surface_tspline_desc
{
	qaws_tspline_control_point const* control_points;
	unsigned int control_point_count;
	unsigned int u_degree;      /* typically 3 (cubic) */
	unsigned int v_degree;      /* typically 3 (cubic) */
	qaws_range u_range;         /* parameter domain in u */
	qaws_range v_range;         /* parameter domain in v */
} qaws_surface_tspline_desc;

qaws_status qaws_surface_create_tspline(
	qaws_surface_tspline_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_TSPLINE_H */
