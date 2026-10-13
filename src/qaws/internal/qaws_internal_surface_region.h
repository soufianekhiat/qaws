#ifndef QAWS_INTERNAL_SURFACE_REGION_H
#define QAWS_INTERNAL_SURFACE_REGION_H

#include "../qaws_types.h"
#include "../qaws_status.h"
#include "../qaws_surface_types.h"

/* A surface's parameter domain, its closed directions and the tolerance in
   (u, v), shared by the surface curve and region code. */
typedef struct qaws_internal_sg_ctx
{
	qaws_surface const* s;
	double u0, u1, v0, v1;
	double pu, pv;              /* periods of closed directions (0: open) */
	double diag, tol;           /* tolerance in (u, v) */
} qaws_internal_sg_ctx;

qaws_status qaws_internal_sg_init(qaws_internal_sg_ctx* c, qaws_surface const* s, double tol);

/* g = {E, F, G}; gam (may be NULL) the Christoffel symbols; (u, v) wrapped
   in closed directions */
qaws_status qaws_internal_sg_metric(qaws_internal_sg_ctx const* c, double u, double v, double* g, double* gam);

/* the surface at (u, v), wrapped in closed directions */
qaws_status qaws_internal_sg_eval(qaws_internal_sg_ctx const* c, double u, double v, unsigned int flags,
	qaws_surface_eval_result* r);

/* E, F, G without the singular check; sqrt(E G - F^2), zero at a pole */
qaws_status qaws_internal_sg_first_form(qaws_internal_sg_ctx const* c, double u, double v, double* g);
qaws_status qaws_internal_sg_area_element(qaws_internal_sg_ctx const* c, double u, double v, double* out);

/* y = (u, v, u', v'): the base curve's point at t and its unit conormal
   (side > 0 left of travel); y advanced by len along its geodesic */
qaws_status qaws_internal_sg_conormal(qaws_internal_sg_ctx const* c, qaws_curve const* base, double t, double side, double* y);
qaws_status qaws_internal_sg_integrate(qaws_internal_sg_ctx const* c, double* y, double len);

/* the geodesic offset of the (u, v) curve base by dist (> 0 left), with the
   conormal geodesics from its ends to the offset's (exactly its ends) */
qaws_status qaws_internal_sg_offset(qaws_internal_sg_ctx const* c, qaws_curve const* base, double dist,
	qaws_curve** out_offset, qaws_curve** out_ray0, qaws_curve** out_ray1);

#endif /* QAWS_INTERNAL_SURFACE_REGION_H */
