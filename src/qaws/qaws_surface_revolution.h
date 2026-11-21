#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_REVOLUTION_H_F32
#define QAWS_SURFACE_REVOLUTION_H_F32
#define QAWS_SURFACE_REVOLUTION_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_REVOLUTION_H_F64
#define QAWS_SURFACE_REVOLUTION_H_F64
#define QAWS_SURFACE_REVOLUTION_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_REVOLUTION_H_BODY
#undef QAWS_SURFACE_REVOLUTION_H_BODY
#define QAWS_SURFACE_REVOLUTION_H

#include "qaws_surface_types.h"

/* Surface of revolution: a 2D profile curve (r, h) rotated about an axis.
   The profile's x is the distance from the axis, its y the height along it.
   u in [0, 1] is the turn (theta = u * angle); v in [0, 1] runs over the
   profile's parameter range.

   S(u,v) = axis_origin + r(v)*cos(theta)*X + r(v)*sin(theta)*Y + h(v)*Z

   Z is axis_direction normalized (INVALID_ARGUMENT when it is zero); X is
   the coordinate axis least aligned with Z made perpendicular to it; Y is
   Z x X.
   angle: total rotation in radians (0: 2*PI, a full revolution).

   The profile must be 2D (INVALID_DIMENSION otherwise). The surface borrows
   the profile curve. */

typedef struct qaws_surface_revolution_desc
{
	qaws_curve const* profile;        /* 2D profile curve (borrowed) */
	qaws_vec3 axis_origin;            /* point on rotation axis */
	qaws_vec3 axis_direction;         /* rotation axis direction (will be normalized) */
	qaws_scalar angle;                /* rotation angle in radians (0 = 2*PI) */
} qaws_surface_revolution_desc;

qaws_status qaws_surface_create_revolution(
	qaws_surface** out_surface,
	qaws_surface_revolution_desc const* desc);

#endif /* QAWS_SURFACE_REVOLUTION_H */
