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
 *     quad between its two positions, split along the diagonal that folds it
 *     outward: two slabs' turning walls overlap, never leave a sliver).
 *   - The swept solid is the solid at 0 united with every slab's extrusion:
 *     the self-union of them all (qaws_mesh_boolean), its winding-number
 *     rule resolving their overlaps and an extrusion's folds.
 * A translation is exact: a cube moved along a segment is its Minkowski sum
 * with it. Curved motions are within the tolerance (the chords), rotating
 * walls within their twist (a bilinear quad split in two).
 *
 * A vertex motion must have the mesh's vertex count; the mesh's positions
 * are then ignored (the motion's at each time are used).
 *
 * The mesh's attributes are carried as qaws_mesh_3d_boolean carries them
 * (vertex ones come out per corner): each layer repeats the mesh's values, a
 * wall along its edge the front face's corners' there; directions (NORMAL,
 * TANGENT, VECTOR, and NORMALIZED of 3 components) turn with a rigid
 * motion's rotation, a QUATERNION is composed with it, and a NORMAL face
 * attribute on a wall is the wall's own normal. Features: the mesh's copies,
 * the paths its feature corners trace and the cuts, kept where the result
 * turns by more than half the mesh's softest feature angle (the mesh's own
 * at the first and last poses wherever it turns); an edge sliding along
 * itself (an envelope) traces no feature.
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

/*
 * Continuous Booleans: a tool swept by a sequence of motions (a tool path,
 * each motion one move) against a stock: out = sweep op stock, op one of
 *   DIFFERENCE    carving: the stock less every point the tool covers
 *                 (CNC simulation; the report's removed volume);
 *   UNION         sculpting: the stroke added;
 *   INTERSECTION  contact: the part of the stock the tool passes through.
 * Each move's sweep (qaws_mesh_3d_sweep), then all of them as one operand
 * (overlapping closed parts: their winding numbers add, a union) against the
 * stock by qaws_mesh_boolean. The volumes are the meshes' signed volumes.
 * The tool's and the stock's attributes are carried, matched by name
 * (qaws_mesh_3d_boolean).
 */
typedef struct qaws_sweep_boolean_report
{
	unsigned int moves;
	unsigned int slabs;                  /* over all the moves */
	qaws_scalar stock_volume;
	qaws_scalar result_volume;
	qaws_scalar changed_volume;          /* DIFFERENCE: removed (stock - result); UNION: added (result - stock); INTERSECTION: the result */
	qaws_boolean_report boolean;
	int closed;
} qaws_sweep_boolean_report;

/* op: QAWS_BOOLEAN_DIFFERENCE, _UNION or _INTERSECTION (sweep op stock, the stock taken from for DIFFERENCE); report and desc
   may be NULL. INVALID_ARGUMENT for no move, a NULL move, an unknown op. An empty result: OK and NULL. */
qaws_status qaws_mesh_sweep_boolean(qaws_mesh_3d** out_mesh, qaws_sweep_boolean_report* report, qaws_mesh_3d const* tool,
	qaws_motion const* const* moves, qaws_mesh_3d const* stock, unsigned int move_count, qaws_boolean_op op, qaws_sweep_desc const* desc);

#endif /* QAWS_MESH_SWEEP_H */
