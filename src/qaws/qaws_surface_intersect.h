#ifndef QAWS_SURFACE_INTERSECT_H
#define QAWS_SURFACE_INTERSECT_H

#include "qaws_surface_types.h"

/* Surface-surface intersection: finds curves where two surfaces meet.
   Uses a marching method with Newton refinement.

   Results are returned as arrays of intersection points.
   Each point contains parameters on both surfaces and the 3D position. */

typedef struct qaws_ssi_point
{
	qaws_scalar u1, v1;   /* parameters on surface A */
	qaws_scalar u2, v2;   /* parameters on surface B */
	qaws_vec3 position;    /* 3D intersection point */
} qaws_ssi_point;

typedef struct qaws_ssi_curve
{
	qaws_ssi_point* points;
	unsigned int point_count;
} qaws_ssi_curve;

typedef struct qaws_ssi_desc
{
	qaws_surface const* surface_a;
	qaws_surface const* surface_b;
	qaws_scalar tolerance;          /* convergence tolerance (0 = default 1e-6) */
	unsigned int grid_samples;      /* seed grid density per axis (0 = default 20) */
	unsigned int max_march_steps;   /* max steps per curve (0 = default 500) */
	qaws_scalar step_size;          /* march step size (0 = auto) */
} qaws_ssi_desc;

/* Find intersection curves between two surfaces.
   out_curves: caller-provided array of qaws_ssi_curve (capacity = curve_capacity).
   point_buffer: caller-provided buffer for all intersection points.
   Each curve in out_curves will point into point_buffer.
   Returns the number of curves found in *out_curve_count. */
qaws_status qaws_surface_intersect(
	qaws_ssi_desc const* desc,
	qaws_ssi_curve* out_curves,
	unsigned int curve_capacity,
	unsigned int* out_curve_count,
	qaws_ssi_point* point_buffer,
	unsigned int point_capacity);

#endif /* QAWS_SURFACE_INTERSECT_H */
