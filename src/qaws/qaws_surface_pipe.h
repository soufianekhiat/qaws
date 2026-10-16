#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_PIPE_H_F32
#define QAWS_SURFACE_PIPE_H_F32
#define QAWS_SURFACE_PIPE_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_PIPE_H_F64
#define QAWS_SURFACE_PIPE_H_F64
#define QAWS_SURFACE_PIPE_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_PIPE_H_BODY
#undef QAWS_SURFACE_PIPE_H_BODY
#define QAWS_SURFACE_PIPE_H

#include "qaws_surface_types.h"
#include "qaws_curve.h"

/* Pipe surface: sweeps a circular or elliptical cross-section along a 3D path curve.
   S(u,v) = path(u) + radius_x * cos(2*pi*v) * N(u) + radius_y * sin(2*pi*v) * B(u)

   u parameterizes along the path [0,1], v parameterizes around the cross-section [0,1].
   N(u), B(u) are from a parallel transport (Bishop) frame along the path curve.
   For circular pipe, set radius_x = radius_y. */

typedef struct qaws_surface_pipe_desc
{
	qaws_curve const* path;        /* 3D path curve (borrowed) */
	qaws_scalar radius_x;          /* cross-section radius in normal direction */
	qaws_scalar radius_y;          /* cross-section radius in binormal direction (0 = same as radius_x) */
} qaws_surface_pipe_desc;

qaws_status qaws_surface_create_pipe(
	qaws_surface_pipe_desc const* desc,
	qaws_surface** out_surface);

#endif /* QAWS_SURFACE_PIPE_H */
