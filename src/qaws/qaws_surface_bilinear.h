#ifndef QAWS_SURFACE_BILINEAR_H
#define QAWS_SURFACE_BILINEAR_H

#include "qaws_surface_types.h"

/* Bilinear patch: degree (1,1) surface defined by four corner points.
   S(u,v) = (1-u)(1-v)P00 + u(1-v)P10 + (1-u)vP01 + uvP11

   Corner layout:
     P01 --- P11
      |       |
     P00 --- P10
   u increases left to right, v increases bottom to top. */

typedef struct qaws_surface_bilinear_desc
{
	qaws_vec3 p00;   /* corner at (u=0, v=0) */
	qaws_vec3 p10;   /* corner at (u=1, v=0) */
	qaws_vec3 p01;   /* corner at (u=0, v=1) */
	qaws_vec3 p11;   /* corner at (u=1, v=1) */
} qaws_surface_bilinear_desc;

qaws_status qaws_surface_create_bilinear(
	qaws_surface_bilinear_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_BILINEAR_H */
