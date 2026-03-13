#ifndef QAWS_SURFACE_EXTRUSION_H
#define QAWS_SURFACE_EXTRUSION_H

#include "qaws_surface_types.h"

/* Extrusion surface: a profile curve translated along a direction vector.
   S(u,v) = profile(u) + v * direction

   u parameter follows the profile curve.
   v parameter controls distance along the extrusion direction (0 to length).

   The profile can be 2D (lifted to z=0 plane) or 3D.
   The surface borrows the profile curve. */

typedef struct qaws_surface_extrusion_desc
{
	qaws_curve const* profile;   /* profile curve, 2D or 3D (borrowed) */
	qaws_vec3 direction;         /* extrusion direction vector */
	qaws_scalar length;          /* extrusion length (0 = use |direction| as unit) */
} qaws_surface_extrusion_desc;

qaws_status qaws_surface_create_extrusion(
	qaws_surface_extrusion_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_EXTRUSION_H */
