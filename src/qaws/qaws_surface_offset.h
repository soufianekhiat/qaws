#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_OFFSET_H_F32
#define QAWS_SURFACE_OFFSET_H_F32
#define QAWS_SURFACE_OFFSET_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_OFFSET_H_F64
#define QAWS_SURFACE_OFFSET_H_F64
#define QAWS_SURFACE_OFFSET_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_OFFSET_H_BODY
#undef QAWS_SURFACE_OFFSET_H_BODY
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
	qaws_surface** out_surface,
	qaws_surface_offset_desc const* desc);

#endif /* QAWS_SURFACE_OFFSET_H */
