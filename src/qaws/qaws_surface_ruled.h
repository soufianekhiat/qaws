#ifndef QAWS_SURFACE_RULED_H
#define QAWS_SURFACE_RULED_H

#include "qaws_surface_types.h"

/* Ruled surface: linear interpolation between two 3D boundary curves.
   S(u,v) = (1-v) * curve_a(u) + v * curve_b(u)

   u parameter is mapped to each curve's own domain [0,1].
   v parameter blends linearly from curve_a (v=0) to curve_b (v=1).

   Both curves must be 3D. The surface borrows (does not own) the curves.
   The caller must keep both curves alive for the lifetime of the surface. */

typedef struct qaws_surface_ruled_desc
{
	qaws_curve const* curve_a;    /* first boundary curve, 3D (borrowed) */
	qaws_curve const* curve_b;    /* second boundary curve, 3D (borrowed) */
} qaws_surface_ruled_desc;

qaws_status qaws_surface_create_ruled(
	qaws_surface** out_surface,
	qaws_surface_ruled_desc const* desc);

#endif /* QAWS_SURFACE_RULED_H */
