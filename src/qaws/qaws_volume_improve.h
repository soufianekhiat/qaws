#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_VOLUME_IMPROVE_H_F32
#define QAWS_VOLUME_IMPROVE_H_F32
#define QAWS_VOLUME_IMPROVE_H_BODY
#endif
#else
#ifndef QAWS_VOLUME_IMPROVE_H_F64
#define QAWS_VOLUME_IMPROVE_H_F64
#define QAWS_VOLUME_IMPROVE_H_BODY
#endif
#endif
#ifdef QAWS_VOLUME_IMPROVE_H_BODY
#undef QAWS_VOLUME_IMPROVE_H_BODY
#define QAWS_VOLUME_IMPROVE_H

#include "qaws_volume_mesh.h"

/*
 * Improvement of tetrahedral meshes: slivers and other poor tetrahedra
 * removed by smoothing and flips (Klingner and Shewchuk 2007, Aggressive
 * tetrahedral mesh improvement), the border kept exactly.
 *
 * A tetrahedron's quality is the smallest sine of its six dihedral angles
 * (negative when it is inverted): small for both too sharp and too flat
 * angles, slivers included. Each pass, on the tetrahedra below the target:
 *   edge removal  an inside edge's ring of n tetrahedra (3 to 10) replaced
 *                 by the 2n - 4 of the ring polygon's best triangulation
 *                 (Shewchuk 2002's dynamic program: the triangulation whose
 *                 worst tetrahedron is best; n = 3 is the 3-2 flip);
 *   face removal  a facet between two tetrahedra replaced by the edge
 *                 between their apexes (the 2-3 flip);
 *   smoothing     an inside vertex moved up the gradient of the worst
 *                 quality round it (the steepest ascent of the smallest of
 *                 its active qualities, the active gradients' convex hull's
 *                 nearest point to zero, a line search); each move kept
 *                 only when it raises the worst.
 * An operation is kept only when it raises the worst quality of the cells it
 * replaces or moves. Passes stop when one changes nothing.
 *
 * The border (facets with no tetrahedron across) and its vertices stay, as
 * do the other cells (hexahedra, prisms, pyramids) and their vertices: the
 * result's border is the input's, triangle for triangle, and its volume the
 * same. The result: the other cells first, in their order, then the
 * tetrahedra; the vertices unchanged in number and order.
 */

typedef struct qaws_volume_improve_desc
{
	unsigned int iterations;            /* passes; 0: 10 */
	qaws_scalar target_dihedral;        /* degrees: tetrahedra whose dihedral angles all lie within [t, 180 - t] are left alone; 0: 40 */
	int no_smoothing;                   /* 1: flips only */
	int no_flips;                       /* 1: smoothing only */
} qaws_volume_improve_desc;

typedef struct qaws_volume_improve_report
{
	unsigned int iterations;            /* passes made */
	unsigned int edge_removals, face_removals, moves;
	qaws_scalar min_quality_before, min_quality_after;      /* the smallest sine of a dihedral angle */
	qaws_volume_quality quality;        /* the result's */
} qaws_volume_improve_report;

/* desc and report may be NULL. INVALID_ARGUMENT for a mesh without tetrahedra. */
qaws_status qaws_volume_mesh_improve(qaws_volume_mesh** out_mesh, qaws_volume_improve_report* report, qaws_volume_mesh const* mesh,
	qaws_volume_improve_desc const* desc);

#endif /* QAWS_VOLUME_IMPROVE_H */
