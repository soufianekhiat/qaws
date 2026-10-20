#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_TRIM_H_F32
#define QAWS_SURFACE_TRIM_H_F32
#define QAWS_SURFACE_TRIM_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_TRIM_H_F64
#define QAWS_SURFACE_TRIM_H_F64
#define QAWS_SURFACE_TRIM_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_TRIM_H_BODY
#undef QAWS_SURFACE_TRIM_H_BODY
#define QAWS_SURFACE_TRIM_H

#include "qaws_surface_types.h"
#include "qaws_curve.h"

/* Trimmed surface: restricts a base surface to a region defined by trim loops
   in (u,v) parameter space.

   A trim loop is a closed sequence of 2D curves in parameter space.
   The outer boundary defines the kept region (points inside are kept).
   Inner loops define holes (points inside are removed).

   If no outer loop is provided, the full [0,1]x[0,1] domain is the outer boundary.
   The base surface is borrowed (not owned). Trim curves are borrowed. */

typedef struct qaws_trim_loop
{
	qaws_curve const* const* curves;  /* array of 2D curves forming closed loop (borrowed) */
	unsigned int curve_count;
	int is_outer;                      /* 1 = outer boundary, 0 = hole */
} qaws_trim_loop;

typedef struct qaws_surface_trim_desc
{
	qaws_surface const* base;          /* base surface (borrowed) */
	qaws_trim_loop const* loops;       /* array of trim loops */
	unsigned int loop_count;
} qaws_surface_trim_desc;

qaws_status qaws_surface_create_trimmed(
	qaws_surface** out_surface,
	qaws_surface_trim_desc const* desc);

/* Query whether a (u,v) point is inside the trimmed region.
   Returns QAWS_STATUS_OK and sets *out_inside to 1 (inside) or 0 (outside). */
qaws_status qaws_surface_trim_contains(
	int* out_inside,
	qaws_surface const* trimmed_surface, qaws_scalar u,
	qaws_scalar v);

#endif /* QAWS_SURFACE_TRIM_H */
