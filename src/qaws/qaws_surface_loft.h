#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_LOFT_H_F32
#define QAWS_SURFACE_LOFT_H_F32
#define QAWS_SURFACE_LOFT_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_LOFT_H_F64
#define QAWS_SURFACE_LOFT_H_F64
#define QAWS_SURFACE_LOFT_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_LOFT_H_BODY
#undef QAWS_SURFACE_LOFT_H_BODY
#define QAWS_SURFACE_LOFT_H

#include "qaws_surface_types.h"
#include "qaws_curve.h"

/* Lofted (skinning) surface: interpolates a sequence of 3D cross-section curves.
   u parameterizes along each cross-section [0,1].
   v parameterizes between cross-sections [0,1].

   At v = v_parameters[i], the surface passes through section_curves[i](u).
   Between sections, cubic B-spline blending interpolates the cross-section positions. */

typedef struct qaws_surface_loft_desc
{
	qaws_curve const* const* sections;  /* array of N 3D cross-section curves (borrowed) */
	unsigned int section_count;          /* N >= 2 */
	qaws_scalar const* v_parameters;    /* N parameter values in [0,1] for each section, NULL = uniform spacing */
} qaws_surface_loft_desc;

qaws_status qaws_surface_create_loft(
	qaws_surface** out_surface,
	qaws_surface_loft_desc const* desc);

#endif /* QAWS_SURFACE_LOFT_H */
