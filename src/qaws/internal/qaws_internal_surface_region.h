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

#endif /* QAWS_INTERNAL_SURFACE_REGION_H */
