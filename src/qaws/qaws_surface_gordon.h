#ifndef QAWS_SURFACE_GORDON_H
#define QAWS_SURFACE_GORDON_H

#include "qaws_surface_types.h"
#include "qaws_curve.h"

/* Gordon surface: interpolates a network of u-curves and v-curves.
   Generalization of Coons patch to arbitrary grid density.

   S(u,v) = L_u(u,v) + L_v(u,v) - T(u,v)

   L_u: loft through u-curves in v direction (each u-curve evaluated at u, blended at v)
   L_v: loft through v-curves in u direction (each v-curve evaluated at v, blended at u)
   T:   tensor product through intersection points Q[i][j] = u_curve[i](v_param[j])

   Compatibility: u_curve[i] evaluated at v_param[j] must equal v_curve[j] evaluated at u_param[i]. */

typedef struct qaws_surface_gordon_desc
{
	qaws_curve const* const* u_curves;  /* M u-isoparameter curves (3D, borrowed) */
	unsigned int u_curve_count;          /* M >= 2 */
	qaws_scalar const* v_params;         /* M parameter values where u-curves sit in v */
	qaws_curve const* const* v_curves;  /* N v-isoparameter curves (3D, borrowed) */
	unsigned int v_curve_count;          /* N >= 2 */
	qaws_scalar const* u_params;         /* N parameter values where v-curves sit in u */
} qaws_surface_gordon_desc;

qaws_status qaws_surface_create_gordon(
	qaws_surface_gordon_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_GORDON_H */
