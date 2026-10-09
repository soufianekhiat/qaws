#ifndef QAWS_INTERNAL_FLATTEN_H
#define QAWS_INTERNAL_FLATTEN_H

#include "../qaws_types.h"
#include "../qaws_status.h"
#include "../qaws_surface_types.h"

/*
 * Flattening for the batched operations.
 *
 * Curves: span by span into chord segments; a piece is split while its
 * midpoint or a quarter point lies farther than the flatness bound from
 * its chord. Surfaces: an adaptive (u, v) quadtree; a patch is split while
 * its centre or an edge midpoint lies farther than the bound from the
 * bilinear interpolation of its corners. Every piece carries a box
 * inflated by twice the measured deviation (r).
 */

typedef struct qaws_flat_seg
{
	qaws_scalar t0, t1;
	qaws_scalar p0[3], p1[3];
	qaws_scalar lo[3], hi[3];   /* inflated box */
	qaws_scalar r;              /* inflation */
	unsigned int owner;         /* curve index */
	unsigned int index;         /* position along its curve */
} qaws_flat_seg;

typedef struct qaws_flat_patch
{
	qaws_scalar u0, u1, v0, v1;
	qaws_scalar p[4][3];        /* corners (u0,v0) (u1,v0) (u0,v1) (u1,v1) */
	qaws_scalar lo[3], hi[3];
	qaws_scalar r;
	unsigned int owner;         /* surface index */
} qaws_flat_patch;

/* position (and d1 when d != NULL) of a 2D or 3D curve as 3 scalars (z = 0 in 2D) */
qaws_status qaws_internal_curve_point(qaws_curve const* c, unsigned int dim, qaws_scalar t, unsigned int flags, qaws_scalar* p, qaws_scalar* d);

/* appends the segments of curve `owner` to *segs (grown with realloc) */
qaws_status qaws_internal_flatten_curve(qaws_curve const* c, unsigned int dim, qaws_scalar flatness, unsigned int owner,
	qaws_flat_seg** segs, unsigned int* count, unsigned int* capacity);

/* appends the patches of surface `owner` to *patches (grown with realloc) */
qaws_status qaws_internal_flatten_surface(qaws_surface const* s, qaws_scalar flatness, unsigned int owner,
	qaws_flat_patch** patches, unsigned int* count, unsigned int* capacity);

/* largest axis extent of a coarse sampling of curves and surfaces (1 when empty) */
qaws_scalar qaws_internal_flatten_extent(qaws_curve const* const* curves, unsigned int curve_count, unsigned int dim,
	qaws_surface const* const* surfaces, unsigned int surface_count);

#endif /* QAWS_INTERNAL_FLATTEN_H */
