#ifndef QAWS_SURFACE_FILLET_H
#define QAWS_SURFACE_FILLET_H

#include "qaws_surface_types.h"

/* Surface fillet: generates a smooth blend surface connecting two adjacent surfaces.

   The fillet is computed using the rolling-ball method:
   1. Find the intersection curve between the two surfaces.
   2. At each point on the intersection curve, compute surface normals.
   3. The fillet ball center lies along the normal bisector at distance
      radius / cos(half_angle) from the intersection point.
   4. For each spine point, construct the fillet cross-section as a
      circular arc from the contact point on surface_a to the contact
      point on surface_b.

   The result is a smooth G1 blend surface connecting the two original
   surfaces.

   Both surfaces are borrowed (not owned). */

typedef struct qaws_surface_fillet_desc
{
	qaws_surface const* surface_a;  /* first surface (borrowed) */
	qaws_surface const* surface_b;  /* second surface (borrowed) */
	qaws_scalar radius;             /* fillet radius */
	unsigned int sample_count;      /* spine sampling resolution (0 = default 64) */
	unsigned int arc_segments;      /* circular arc segments (0 = default 8) */
} qaws_surface_fillet_desc;

/* Create a fillet surface between two adjacent surfaces.
   The fillet is a swept circular arc along the intersection spine.
   Returns a new surface representing the fillet blend. */
qaws_status qaws_surface_create_fillet(
	qaws_surface_fillet_desc const* desc,
	qaws_surface** out_surface);

/* Convenience: compute fillet and also return the trim curves on each surface.
   out_trim_a: 2D curve in surface_a parameter space where fillet meets surface_a.
   out_trim_b: 2D curve in surface_b parameter space where fillet meets surface_b.
   These curves can be used to trim the original surfaces. */
qaws_status qaws_surface_create_fillet_with_trims(
	qaws_surface_fillet_desc const* desc,
	qaws_surface** out_fillet,
	qaws_curve** out_trim_a,
	qaws_curve** out_trim_b);

#endif /* QAWS_SURFACE_FILLET_H */
