#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_DECIMATE_H_F32
#define QAWS_MESH_DECIMATE_H_F32
#define QAWS_MESH_DECIMATE_H_BODY
#endif
#else
#ifndef QAWS_MESH_DECIMATE_H_F64
#define QAWS_MESH_DECIMATE_H_F64
#define QAWS_MESH_DECIMATE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_DECIMATE_H_BODY
#undef QAWS_MESH_DECIMATE_H_BODY
#define QAWS_MESH_DECIMATE_H

#include "qaws_mesh_measure.h"

/*
 * Decimation by quadric error metrics (Garland, Heckbert 1997), in place:
 * edges collapsed cheapest first, each merged vertex placed where the sum
 * of its two ends' quadrics is least (the planes of the faces they had,
 * weighted by area). The error is read as a distance: the root of the
 * area-weighted mean squared distance to those planes.
 *
 * What stays: feature corners, and the mesh's feature lines and border,
 * held near their place by planes standing on their edges (the border's
 * vertices held still with QAWS_DECIMATE_LOCK_BORDER); collapses that would
 * turn a face over or move a feature vertex off its line are refused.
 * Attributes follow the collapses by their rules (qaws_mesh.h).
 *
 * It stops at the target face count, before a collapse whose error passes
 * max_error, or when no collapse is left. With QAWS_DECIMATE_MEASURE the
 * report holds the result's certified deviation from the input
 * (qaws_mesh_3d_measure_deviation against the input as a LINEAR
 * reference).
 */

#define QAWS_DECIMATE_LOCK_BORDER 1u        /* the border's vertices are neither removed nor moved */
#define QAWS_DECIMATE_MEASURE 2u            /* the report's deviation measured */

typedef struct qaws_decimate_desc
{
	unsigned int target_faces;          /* stop at this many faces; 0: no count */
	qaws_scalar max_error;              /* stop before an error past this distance; 0: no bound */
	unsigned int flags;                 /* QAWS_DECIMATE_* */
} qaws_decimate_desc;

typedef struct qaws_decimate_report
{
	unsigned int faces_before, faces_after;
	unsigned int collapses;
	qaws_scalar max_error;              /* the largest error of a collapse made */
	qaws_mesh_deviation deviation;      /* with QAWS_DECIMATE_MEASURE */
} qaws_decimate_report;

/* report may be NULL; desc NULL: half the faces */
qaws_status qaws_mesh_3d_decimate(qaws_decimate_report* report, qaws_mesh_3d* mesh, qaws_decimate_desc const* desc);

#endif /* QAWS_MESH_DECIMATE_H */
