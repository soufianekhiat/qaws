#ifndef QAWS_INTERNAL_BATCH_H
#define QAWS_INTERNAL_BATCH_H

#include "../qaws_curve_batch.h"
#include "../qaws_surface_batch.h"
#include "qaws_internal_flatten.h"

/* Prepared sets, shared by the curve and the surface batches. */

struct qaws_curve_set
{
	qaws_curve_batch_desc desc;     /* curves and families owned by the set */
	unsigned int dim;
	qaws_flat_seg* segs;
	unsigned int nseg;
	unsigned int* seg_count;
	qaws_scalar flat, ext, pos_tol;
};

struct qaws_surface_set
{
	qaws_surface_batch_desc desc;   /* surfaces and families owned by the set */
	qaws_flat_patch* patches;
	unsigned int npatch;
	qaws_scalar flat, ext;
};

#endif /* QAWS_INTERNAL_BATCH_H */
