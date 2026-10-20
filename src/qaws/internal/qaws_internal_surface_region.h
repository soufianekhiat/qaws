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

qaws_status qaws_internal_sg_init(qaws_surface const* s, double tol, qaws_internal_sg_ctx* c);

/* g = {E, F, G}; gam (may be NULL) the Christoffel symbols; (u, v) wrapped
   in closed directions */
qaws_status qaws_internal_sg_metric(double* g, double* gam, double u, double v, qaws_internal_sg_ctx const* c);

/* the surface at (u, v), wrapped in closed directions */
qaws_status qaws_internal_sg_eval(qaws_surface_eval_result* r, double u, double v, unsigned int flags,
	qaws_internal_sg_ctx const* c);

/* E, F, G without the singular check; sqrt(E G - F^2), zero at a pole */
qaws_status qaws_internal_sg_first_form(double* g, double u, double v, qaws_internal_sg_ctx const* c);
qaws_status qaws_internal_sg_area_element(double* out, double u, double v, qaws_internal_sg_ctx const* c);

/* y = (u, v, u', v'): the base curve's point at t and its unit conormal
   (side > 0 left of travel); y advanced by len along its geodesic */
qaws_status qaws_internal_sg_conormal(double* y, qaws_curve const* base, double t, double side, qaws_internal_sg_ctx const* c);
qaws_status qaws_internal_sg_integrate(double* y, double len, qaws_internal_sg_ctx const* c);

/* the geodesic offset of the (u, v) curve base by dist (> 0 left), with the
   conormal geodesics from its ends to the offset's (exactly its ends) */
/* the boundary foot of x: (s, t) with the conormal geodesic from base(s) on
   the side (> 0 left) reaching x after t (x's period copy nearest it),
   Newton from the given (s, t), s kept in [lo, hi]; 1 when it lands */
int qaws_internal_sg_foot(double* s, double* t, qaws_curve const* base, double const* x, double side,
	double lo, double hi, qaws_internal_sg_ctx const* c);

/* the geodesic from a reaching x (its nearest period copy): angle th (as
   sg_dir) and length L by shooting from the given ones; 1 when it lands */
int qaws_internal_sg_reach(double* th, double* L, double const* a, double const* x, qaws_internal_sg_ctx const* c);

/* the angle of a (u, v) direction (from S_u, about S_u x S_v); the geodesic
   leaving a at angle th run L: its end and velocity (vel may be NULL) */
qaws_status qaws_internal_sg_angle(double* out, double const* d, double u, double v, qaws_internal_sg_ctx const* c);
qaws_status qaws_internal_sg_shoot(double* end, double* vel, double const* a, double th, double L, qaws_internal_sg_ctx const* c);

qaws_status qaws_internal_sg_offset(qaws_curve** out_offset, qaws_curve** out_ray0, qaws_curve** out_ray1,
	qaws_curve const* base, double dist, qaws_internal_sg_ctx const* c);

#endif /* QAWS_INTERNAL_SURFACE_REGION_H */
