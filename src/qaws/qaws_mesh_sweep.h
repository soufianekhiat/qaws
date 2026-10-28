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
 * the edge's length (desc.strict bounds it by the tolerance, on more slabs).
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
 * an edge of the revolution whose two ends come from one edge of the body
 * takes its values, as does a ring a body vertex or edge traces from the
 * body's edge there lying along it (within 30 degrees), other edges zeros.
 * The report's turns says
 * a spin was found. An asymmetric body spinning while it drifts past the
 * tolerance a turn is swept implicitly (qaws_mesh_3d_sweep_implicit, the
 * report's implicit set: within a cell of desc.resolution's grid, its
 * attributes carried; by slabs its hulls would pile up turn after turn), or
 * with revolve_spins revolved anyway: the result then
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
 * wall along its edge the front face's corners' there (a wall lying in the
 * plane of the face across that edge, a face sliding along itself: that
 * face's); directions (NORMAL,
 * TANGENT, VECTOR, and NORMALIZED of 3 components) turn with a rigid
 * motion's rotation, a QUATERNION is composed with it, and a NORMAL face
 * attribute on a wall is the wall's own normal. Features: the mesh's copies,
 * the paths the vertices of its feature edges trace (a box's corners, a
 * rim's envelope where it runs along the motion), its feature edges' chords
 * across two poses and the cuts, kept where the result turns by more than
 * half the mesh's softest feature angle (the mesh's own at the first and
 * last poses wherever it turns).
 */

/* Attributes a sweep can add (qaws_sweep_desc.generate, ORed), per corner of the result, from the time its place was swept
   (the first touch: the implicit sweep's least time, the revolution's first pass) and the body's point that swept it:
     TIME         "time" (LINEAR): that time in [0, 1];
     INDEX        "index" (NEAREST): the whole turns the body had made by then (a spin's pass);
     ANGLE        "angle" (LINEAR): the angle it had turned by then (radians, summed over turns);
     ORIENTATION  "orientation" (QUATERNION, w x y z): its rotation then;
     POSITION     "pose_position" (LINEAR, 3): its translation then;
     VELOCITY     "velocity" (VECTOR, 3): the point's velocity then;
     DISTANCE     "distance" (LINEAR): the length its origin had travelled by then;
     SOURCE       "source_face" (NEAREST): the body's face of that point, and "source_uv" (LINEAR, 2): the point's weights on
                  the face's triangle's second and third corners (any of the body's values can be looked up from them);
     PASSES       "passes" (NEAREST): how many times the body passed there (see qaws_sweep_combine; by slabs 1);
     DWELL        "dwell" (LINEAR): the time it stayed within reach over them (the implicit sweep: within 1.5 cells; else 0).
   A vertex motion has no turn, pose or origin: INDEX, ANGLE and DISTANCE 0, ORIENTATION the identity, POSITION 0 (its
   VELOCITY the point's). A generated name replaces an attribute of that name on the body. */
#ifndef QAWS_ENUM_QAWS_SWEEP_GENERATE
#define QAWS_ENUM_QAWS_SWEEP_GENERATE
typedef enum qaws_sweep_generate
{
	QAWS_SWEEP_GENERATE_TIME = 1,
	QAWS_SWEEP_GENERATE_INDEX = 2,
	QAWS_SWEEP_GENERATE_ANGLE = 4,
	QAWS_SWEEP_GENERATE_ORIENTATION = 8,
	QAWS_SWEEP_GENERATE_POSITION = 16,
	QAWS_SWEEP_GENERATE_VELOCITY = 32,
	QAWS_SWEEP_GENERATE_DISTANCE = 64,
	QAWS_SWEEP_GENERATE_SOURCE = 128,
	QAWS_SWEEP_GENERATE_PASSES = 256,
	QAWS_SWEEP_GENERATE_DWELL = 512
} qaws_sweep_generate;
#endif

/* How the values of the passes over one place combine (qaws_sweep_desc.combines, per attribute by name: the body's, a
   track's or a generated one). A place is passed again when the body comes back to it: a spin each turn, the implicit
   sweep each run of time it comes within 1.5 cells; by slabs only the first pass is known. FIRST (the default) the first
   touch's; LAST the last pass's; MINIMUM, MAXIMUM and SUM per component; AVERAGE by the attribute's rule. */
#ifndef QAWS_ENUM_QAWS_SWEEP_COMBINE
#define QAWS_ENUM_QAWS_SWEEP_COMBINE
typedef enum qaws_sweep_combine
{
	QAWS_SWEEP_COMBINE_FIRST = 0,
	QAWS_SWEEP_COMBINE_LAST,
	QAWS_SWEEP_COMBINE_MINIMUM,
	QAWS_SWEEP_COMBINE_MAXIMUM,
	QAWS_SWEEP_COMBINE_AVERAGE,
	QAWS_SWEEP_COMBINE_SUM
} qaws_sweep_combine;
#endif

typedef struct qaws_attribute_combine
{
	char const* name;
	qaws_sweep_combine combine;
} qaws_attribute_combine;

/* An attribute of the body animated over the motion: keys at increasing times in [0, 1], each holding every value of the
   attribute (its element count times its components, as qaws_mesh_3d_get_attribute gives them). The sweep reads it at the
   time each place was swept (as the generated attributes), between the two keys about it blended by the attribute's rule
   (colours in their space, quaternions on the short arc, NEAREST the nearer key); before the first key the first, after
   the last the last. Borrowed for the call. */
typedef struct qaws_attribute_track
{
	char const* name;                   /* an attribute of the body (vertex, corner or face) */
	qaws_scalar const* times;           /* key_count, increasing */
	qaws_scalar const* values;          /* key_count times the attribute's values */
	unsigned int key_count;
} qaws_attribute_track;

typedef struct qaws_sweep_desc
{
	qaws_scalar tolerance;              /* the chords' largest deviation; 0: 1e-3 of the mesh's box diagonal */
	unsigned int max_slabs;             /* 0: 256 */
	int simplify;                       /* nonzero: the vertices in flat regions and on straight creases the overlapping slabs
	                                       leave removed (each collapse within a thousandth of the tolerance), the flat regions
	                                       then turned toward Delaunay: far fewer faces, the shape, features and attributes kept */
	qaws_batch_executor const* executor;   /* the Boolean's pair tests and samples, contact detection's shares of triangles and
	                                          the implicit sweep's queries through it; NULL: this thread */
	int revolve_spins;                  /* nonzero: a spin drifting past the tolerance a turn is revolved even when the body is not
	                                       symmetric about its axis (a one-flute cutter: the result may pass the true sweep by up
	                                       to the report's spin_drift); 0: such a spin is swept by slabs (slow) */
	unsigned int resolution;            /* qaws_mesh_3d_sweep_implicit: cells along the result's box's longest side; 0: 64 */
	unsigned int generate;              /* qaws_sweep_generate flags ORed: the attributes to add */
	qaws_attribute_track const* tracks; /* attributes animated over the motion (NULL: none) */
	unsigned int track_count;
	qaws_attribute_combine const* combines;   /* per attribute how its passes combine (NULL: the first touch for all) */
	unsigned int combine_count;
	int strict;                         /* nonzero: slabs also halved until no edge twists by more than 4 tolerances over one (an
	                                       edge sliding along itself while turning: its hull's notch, a quarter of the twist,
	                                       within the tolerance), so the result passes the true sweep by about the tolerance
	                                       at most; many more slabs where long edges turn (max_slabs still caps them) */
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
	int implicit;                       /* swept on a grid (qaws_mesh_3d_sweep_implicit): an asymmetric spin drifting past the
	                                       tolerance a turn without revolve_spins */
} qaws_sweep_report;

/* desc and report may be NULL. INVALID_ARGUMENT for a vertex motion of another vertex count. An empty result: OK and NULL. */
qaws_status qaws_mesh_3d_sweep(qaws_mesh_3d** out_mesh, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion,
	qaws_sweep_desc const* desc);

/*
 * The implicit route: the swept solid as the level set phi = 0 of
 * phi(x) = min over t of the moving shape's signed distance at x (the
 * BVH's, negative inside by the winding number), contoured on a grid by dual
 * contouring as qaws_mesh_3d_offset's CONTOUR (crossings by regula falsi,
 * their normals the closest point's direction, or on the surface phi's
 * differences: sharp edges and corners kept where the cells resolve them). phi is sampled in time with no vertex moving
 * more than two cells between samples (a vertex motion: a quarter of one), each
 * point's samples pruned by their bounding spheres; a rigid motion's best
 * time is refined by golden sections between the samples beside it (no
 * staircase where a face slides along the sweep). Smooth, never inside out,
 * any motion (spins, vertex motions, self-overlapping paths) at the cost of
 * the grid: detail below the cell is lost, edges between cells are bevelled.
 * For sculpting and as a cross-check of qaws_mesh_3d_sweep. Attributes:
 * each result vertex's source is the shape's closest point at its best time
 * (a face touching it over a span of time: the first touch, as the spins'
 * revolution; a rigid motion's bisected to it); a
 * corner takes it on that triangle or, where it lies on an edge or corner of
 * the shape (a wall an edge sweeps, a seam), on the triangle there whose
 * normal at that time is nearest the corner's own triangle's: values sharp
 * across seams, continuous where they are. Rules as qaws_mesh_3d_sweep
 * (directions turned with a rigid motion, vertex values per corner, a
 * NORMAL straying past 30 degrees its triangle's own); edge values zeros,
 * features none. desc: resolution, max_slabs (the most time
 * samples; 0: 4096) and executor are read. The report: slabs the time
 * samples, deviation half the most a vertex moves between two (the sampling's
 * bound before refinement), closed.
 */
qaws_status qaws_mesh_3d_sweep_implicit(qaws_mesh_3d** out_mesh, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion,
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
 * no event, and the search goes on past them. A vertex motion moving each
 * vertex on its chord at a constant rate is straight too: the events are
 * roots of cubics in t on exact integers, isolated by exact halving in the
 * Bernstein basis to 2^-62 and checked there (time within 2^-62, rounded).
 * Otherwise a vertex motion's bracket is conservative: a triangle whose
 * corners move apart sweeps a curved surface its prisms hold, so they may
 * touch before it does (none comes before time).
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

/* A mesh against itself (a vertex motion: a cloth folding, a strip curling, parts of one mesh meeting): the same search
   over the pairs of its triangles sharing no vertex (vertices at one position at 0 count as one), each pair once;
   report.face and obstacle_face are both its triangles. Triangles crossing at 0 give time 0; a rigid motion keeps the
   mesh apart from itself (no contact past 0). Neighbours (a shared vertex) are not tested: a fold within one vertex's
   ring is missed. INVALID_ARGUMENT as qaws_mesh_3d_first_contact */
qaws_status qaws_mesh_3d_first_self_contact(qaws_contact_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_sweep_desc const* desc);

#endif /* QAWS_MESH_SWEEP_H */
