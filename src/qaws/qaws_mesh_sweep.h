#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_SWEEP_H_F32
#define QAWS_MESH_SWEEP_H_F32
#define QAWS_MESH_SWEEP_H_BODY
#endif
#else
#ifndef QAWS_MESH_SWEEP_H_F64
#define QAWS_MESH_SWEEP_H_F64
#define QAWS_MESH_SWEEP_H_BODY
#endif
#endif
#ifdef QAWS_MESH_SWEEP_H_BODY
#undef QAWS_MESH_SWEEP_H_BODY
#define QAWS_MESH_SWEEP_H

#include "qaws_mesh.h"
#include "qaws_mesh_boolean.h"
#include "qaws_motion.h"

/*
 * Swept solids: the region a closed mesh covers while it moves (a
 * qaws_motion over t in [0, 1]).
 *
 *   - Time is cut into slabs, each halved until every vertex's trajectory
 *     lies within tolerance of its chord over it (checked at the slab's
 *     middle and quarters; a straight motion is one slab).
 *   - Over a slab the vertices move along their chords. A point the solid
 *     covers during it and not at its start was crossed by a face moving
 *     toward it: a front face, its normal at the start against its
 *     vertices' mean displacement. The slab's sweep is the solid at its start
 *     with its front faces extruded: those faces at the start (turned) and
 *     at the end, and walls along their region's border (each border edge's
 *     quad between its two positions, split along a diagonal).
 *   - The swept solid is the solid at 0 united with every slab's extrusion:
 *     the self-union of them all (qaws_mesh_boolean), its winding-number
 *     rule resolving their overlaps and an extrusion's folds.
 * A translation is exact: a cube moved along a segment is its Minkowski sum
 * with it. Curved motions are within the tolerance (the chords), rotating
 * walls within their twist (a bilinear quad split in two).
 *
 * A vertex motion must have the mesh's vertex count; the mesh's positions
 * are then ignored (the motion's at each time are used).
 */

typedef struct qaws_sweep_desc
{
	qaws_scalar tolerance;              /* the chords' largest deviation; 0: 1e-3 of the mesh's box diagonal */
	unsigned int max_slabs;             /* 0: 256 */
	qaws_batch_executor const* executor;   /* the Boolean's pair tests and samples through it; NULL: this thread */
} qaws_sweep_desc;

typedef struct qaws_sweep_report
{
	unsigned int slabs;
	qaws_scalar deviation;              /* the largest chord deviation measured (at the slabs' middles and quarters) */
	unsigned int front_faces;           /* summed over the slabs */
	qaws_boolean_report boolean;
	int closed;
} qaws_sweep_report;

/* desc and report may be NULL. INVALID_ARGUMENT for a vertex motion of another vertex count. An empty result: OK and NULL. */
qaws_status qaws_mesh_3d_sweep(qaws_mesh_3d** out_mesh, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion,
	qaws_sweep_desc const* desc);

#endif /* QAWS_MESH_SWEEP_H */
