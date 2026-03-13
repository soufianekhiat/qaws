#ifndef QAWS_SURFACE_COONS_H
#define QAWS_SURFACE_COONS_H

#include "qaws_surface_types.h"

/* Coons patch: bilinearly blended surface from four boundary curves.
   S(u,v) = Lc(u,v) + Ld(u,v) - B(u,v)
   where Lc is the ruled surface interpolating c0(u) and c1(u),
   Ld is the ruled surface interpolating d0(v) and d1(v),
   and B is the bilinear interpolation of the four corner points.

   Boundary layout:
     c0(u) : bottom edge (v=0), from u=0 to u=1
     c1(u) : top edge (v=1), from u=0 to u=1
     d0(v) : left edge (u=0), from v=0 to v=1
     d1(v) : right edge (u=1), from v=0 to v=1

   Corner consistency: c0(0)=d0(0), c0(1)=d1(0), c1(0)=d0(1), c1(1)=d1(1).
   All curves must be 3D. The surface borrows the curves. */

typedef struct qaws_surface_coons_desc
{
	qaws_curve const* c0;    /* bottom edge (v=0), 3D (borrowed) */
	qaws_curve const* c1;    /* top edge (v=1), 3D (borrowed) */
	qaws_curve const* d0;    /* left edge (u=0), 3D (borrowed) */
	qaws_curve const* d1;    /* right edge (u=1), 3D (borrowed) */
} qaws_surface_coons_desc;

qaws_status qaws_surface_create_coons(
	qaws_surface_coons_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_COONS_H */
