#ifndef QAWS_SURFACE_SWEPT_H
#define QAWS_SURFACE_SWEPT_H

#include "qaws_surface_types.h"

/* Swept surface: cross-section profile swept along a 3D path curve.
   S(u,v) = path(u) + profile_x(v) * N(u) + profile_y(v) * B(u)
   where N, B are the Frenet normal and binormal of the path.

   u parameter follows the path curve domain.
   v parameter follows the profile curve domain.

   The path must be a 3D curve. The profile must be a 2D curve whose
   x,y coordinates map to the normal and binormal directions of the frame.

   The surface borrows (does not own) the path and profile curves.
   The caller must keep both curves alive for the lifetime of the surface. */

typedef struct qaws_surface_swept_desc
{
	qaws_curve const* path;       /* 3D path curve (borrowed) */
	qaws_curve const* profile;    /* 2D cross-section curve (borrowed) */
	qaws_scalar scale;            /* uniform scale of cross-section (0 = default 1.0) */
} qaws_surface_swept_desc;

qaws_status qaws_surface_create_swept(
	qaws_surface** out_surface,
	qaws_surface_swept_desc const* desc);

#endif /* QAWS_SURFACE_SWEPT_H */
