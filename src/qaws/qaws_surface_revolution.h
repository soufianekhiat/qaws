#ifndef QAWS_SURFACE_REVOLUTION_H
#define QAWS_SURFACE_REVOLUTION_H

#include "qaws_surface_types.h"

/* Surface of revolution: a 2D profile curve rotated around an axis.
   The profile curve lies in the XZ plane (x = radius, z = height).
   u parameter is the rotation angle [0, angle_range] mapped to [0,1].
   v parameter follows the profile curve domain.

   S(u,v) = axis_origin + r(v)*cos(theta(u))*X + r(v)*sin(theta(u))*Y + h(v)*Z
   where r(v) = profile_x(v) and h(v) = profile_y(v), theta = u * angle.

   axis_origin: point on the rotation axis.
   axis_direction: unit direction of the rotation axis (default Y-up: (0,1,0)).
   angle: total rotation in radians (default 2*PI for full revolution).

   The profile must be 2D. The surface borrows the profile curve. */

typedef struct qaws_surface_revolution_desc
{
	qaws_curve const* profile;        /* 2D profile curve (borrowed) */
	qaws_vec3 axis_origin;            /* point on rotation axis */
	qaws_vec3 axis_direction;         /* rotation axis direction (will be normalized) */
	qaws_scalar angle;                /* rotation angle in radians (0 = 2*PI) */
} qaws_surface_revolution_desc;

qaws_status qaws_surface_create_revolution(
	qaws_surface_revolution_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_REVOLUTION_H */
