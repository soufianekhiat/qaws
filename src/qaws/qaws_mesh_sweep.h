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
 *     covers during it and not at its start was crossed by a face moving out
 *     toward it. For each face with a corner moving out of the solid by more
 *     than the tolerance (its normal at the start against its displacement)
 *     the slab adds the convex hull of the face's two positions: it holds
 *     every position the face passes through (a turning face moving out at
 *     some corners and in at others included), and never turns inside out
 *     as an extruded face's twisted walls would.
 *   - The swept solid is the solid at 0 united with every slab's hulls: the
 *     self-union of them all (qaws_mesh_boolean). Caps that coincide turned
 *     the other way (a hull's end and the next one's start) cancel and are
 *     left out, as are any two triangles of the soup on the same vertices
 *     turned opposite ways.
 * A translation is exact: a cube moved along a segment is its Minkowski sum
 * with it. Curved motions are within the tolerance (the chords); where an
 * edge slides along itself while turning, the hull fills the notch of its
 * quad, so the result may pass the true sweep by about the slabs' turn times
 * the edge's length.
 *
 * Spins: a rigid motion turning a whole turn or more about one axis (its
 * angular velocity keeping its direction and sense) passes over the same
 * places again and again: a grinding wheel, a drill, an end mill. When it
 * moves at most the tolerance along any turn, or the body is itself symmetric
 * about the axis within the tolerance (every vertex, edge middle and face
 * centre on the border of its profile: any drift), Each turn reaches every angle, so its sweep is
 * the body's solid of revolution about the axis (through the body's origin
 * at 0) swept along the drift (the translation, as above). The solid of
 * revolution (the profile's points on a straight line between their
 * neighbours left out) comes from the body's profile in the half plane of the
 * axis:
 * sliced across the axis between its vertex heights (a slice's topology
 * holds between them), each slice's pieces give the distance range they
 * cover (from 0 when the axis runs through them, out to their farthest
 * point), joined into the profile (the union of those quads, exact integer
 * clipping), then revolved with chords within the tolerance. Attributes
 * come from the body's point the turn first brings to each place (the
 * least angle behind it on its circle: the passes coincide, the first
 * stands; a flat face across the circle its own point), directions turned
 * by that angle; a NORMAL straying past 30 degrees from its face (a corner or
 * edge of the body tracing the surface has no normal of it) the face's own;
 * edge values do not survive the revolution (zeros). The report's turns says
 * a spin was found. An asymmetric body spinning while it drifts past the
 * tolerance a turn is swept by slabs as any motion (its hulls overlap turn
 * after turn: slow), or with revolve_spins revolved anyway: the result then
 * passes the true sweep by up to the report's spin_drift (a CNC simulator's
 * spindle tool).
 *
 * A vertex motion must have the mesh's vertex count; the mesh's positions
 * are then ignored (the motion's at each time are used).
 *
 * The mesh should be closed: an open one (a box without its top) is swept as
 * its faces are, its volume right, the result open where it is (repair it
 * first, qaws_mesh_repair.h). A capped max_slabs leaves chords past the
 * tolerance: the report's deviation says by how much.
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
	int simplify;                       /* nonzero: the vertices in flat regions and on straight creases the overlapping slabs
	                                       leave removed (each collapse within a thousandth of the tolerance), the flat regions
	                                       then turned toward Delaunay: far fewer faces, the shape, features and attributes kept */
	qaws_batch_executor const* executor;   /* the Boolean's pair tests and samples through it; NULL: this thread */
	int revolve_spins;                  /* nonzero: a spin drifting past the tolerance a turn is revolved even when the body is not
	                                       symmetric about its axis (a one-flute cutter: the result may pass the true sweep by up
	                                       to the report's spin_drift); 0: such a spin is swept by slabs (slow) */
} qaws_sweep_desc;

typedef struct qaws_sweep_report
{
	unsigned int slabs;
	qaws_scalar deviation;              /* the largest chord deviation measured (at the slabs' middles and quarters) */
	unsigned int front_faces;           /* summed over the slabs */
	qaws_boolean_report boolean;
	int closed;
	qaws_scalar turns;                  /* a spin's turns, swept as its solid of revolution along its drift; 0 otherwise */
	qaws_scalar spin_drift;             /* a spin's largest drift along one turn */
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

/*
 * Continuous collision detection: the first time a moving mesh meets an
 * obstacle, standing still or moving too (its own motion), or none. Time is
 * cut into the sweep's slabs (both bodies' chords within the tolerance);
 * slab by slab, in time order, each triangle's swept prism (its two
 * positions, its sides split in two, its three tetrahedra) is tested against
 * the obstacle's prisms its box meets (a BVH over the obstacle's prisms,
 * refit per step; exact triangle predicates), a corner of either inside the
 * other's prism counting too: a thin wall is not passed through at any
 * speed. The first slab that touches is halved, its touching halves kept,
 * until a step moves no vertex of either more than the tolerance: the
 * contact lies between time and time_upper (no contact before time, up to
 * the tolerance: two moving prisms may meet with the triangles a step apart
 * in time).
 * Straight moves (both bodies' rotations fixed, their translations linear:
 * every vertex at q + T0 + t (T1 - T0)) are solved exactly: a vertex
 * reaching a face's plane and two edges coplanar are linear in t, so among
 * the pairs touching about the bracket each such event's time is a
 * rational, its point checked inside the face or both edges, all on the
 * exact kernel's integers; the earliest is the contact (exact set, time and
 * time_upper the double nearest it). Pairs only within the tolerance meet
 * no event, and the search goes on past them.
 * Meshes meeting at 0 (triangles crossing, or one inside the other by the
 * winding number) give time 0. The overlap's volume over the whole motion
 * is qaws_mesh_sweep_boolean's INTERSECTION.
 */
typedef struct qaws_contact_report
{
	int contact;                        /* the meshes meet during the motion */
	qaws_scalar time, time_upper;       /* the first contact between them (1 and 1 with none) */
	int exact;                          /* the time exact (straight moves, or 0): the double nearest the rational */
	unsigned int face;                  /* the moving mesh's triangle and the obstacle's that meet first (QAWS_MESH_NONE when
	                                       one mesh is inside the other at 0, or with none) */
	unsigned int obstacle_face;
	qaws_vec3 point;                    /* near the contact: the obstacle triangle's centre (or a vertex inside the other) */
	unsigned int slabs;
	unsigned int tests;                 /* triangle pairs tested */
} qaws_contact_report;

/* obstacle_motion NULL: the obstacle stands still. desc: the tolerance and max_slabs (as the sweep's); report may be NULL.
   INVALID_ARGUMENT for a vertex motion of another vertex count, a mesh without faces */
qaws_status qaws_mesh_3d_first_contact(qaws_contact_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_mesh_3d const* obstacle,
	qaws_motion const* obstacle_motion, qaws_sweep_desc const* desc);

#endif /* QAWS_MESH_SWEEP_H */
