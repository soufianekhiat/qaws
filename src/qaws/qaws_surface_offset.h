#ifndef QAWS_SURFACE_OFFSET_H
#define QAWS_SURFACE_OFFSET_H

#include "qaws_surface_types.h"

/* Offset surface: displaces a base surface by a constant distance along its normal.
   S_off(u,v) = S_base(u,v) + distance * N_base(u,v)

   Positive distance offsets along the surface normal direction.
   The base surface is borrowed (not owned). */

typedef struct qaws_surface_offset_desc
{
	qaws_surface const* base;   /* base surface (borrowed) */
	qaws_scalar distance;       /* offset distance (positive = along normal) */
} qaws_surface_offset_desc;

qaws_status qaws_surface_create_offset(
	qaws_surface_offset_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_OFFSET_H */
