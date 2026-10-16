#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_BIQUADRATIC_H_F32
#define QAWS_SURFACE_BIQUADRATIC_H_F32
#define QAWS_SURFACE_BIQUADRATIC_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_BIQUADRATIC_H_F64
#define QAWS_SURFACE_BIQUADRATIC_H_F64
#define QAWS_SURFACE_BIQUADRATIC_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_BIQUADRATIC_H_BODY
#undef QAWS_SURFACE_BIQUADRATIC_H_BODY
#define QAWS_SURFACE_BIQUADRATIC_H

#include "qaws_surface_types.h"

/* Biquadratic patch: degree (2,2) surface with 3x3 control points.
   Evaluated using Bernstein basis B_i^2(t) = {(1-t)^2, 2t(1-t), t^2}.

   Control point layout (row-major, 3 rows of 3):
     cp[0..2] = v=0 row (bottom), cp[3..5] = v=1/2 row, cp[6..8] = v=1 row
   Within each row, u increases: cp[0]=u=0, cp[1]=u=1/2, cp[2]=u=1 */

typedef struct qaws_surface_biquadratic_desc
{
	qaws_vec3 control_points[9];  /* 3x3 grid, row-major */
} qaws_surface_biquadratic_desc;

qaws_status qaws_surface_create_biquadratic(
	qaws_surface** out_surface,
	qaws_surface_biquadratic_desc const* desc);

#endif /* QAWS_SURFACE_BIQUADRATIC_H */
