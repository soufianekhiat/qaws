#ifndef QAWS_CLIP_H
#define QAWS_CLIP_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_batch_executor.h"

/*
 * Boolean operations on regions of curves (Clipper2's Clipper64 / ClipperD,
 * on curves).
 *
 * Subjects and clips are regions: closed paths read with the fill rule. Open
 * subjects are clipped by the closed regions. The result is made of new
 * curves that lie exactly on the inputs (qaws_curve_extract pieces), in closed
 * paths with their nesting (Clipper2's PolyTree) and open paths.
 *
 * How it works: every input curve goes through one batched intersection
 * (crossings, touches and shared stretches); hits within tolerance are one
 * vertex; curves are cut at every vertex into edges, and edges that run
 * along each other are merged into one edge that carries the winding steps
 * of all of them. The edges form a planar graph; its faces get their
 * subject and clip winding numbers by walking across edges from the
 * unbounded face, and a face is inside the result by the clip type on the
 * fill rule of both:
 *
 *   intersection  in subject and in clip       union  in subject or in clip
 *   difference    in subject and not in clip   xor    in exactly one
 *
 * The boundary between inside and outside faces, walked with the inside on
 * the left, gives the result: outer paths counter-clockwise, holes
 * clockwise (reversed with QAWS_CLIP_REVERSE_SOLUTION). Two parts touching
 * at a point come out as two paths. Nesting comes from the faces, not from
 * point-in-polygon tests.
 *
 * Thread-safe on immutable curves.
 */

typedef enum qaws_clip_type
{
	QAWS_CLIP_NONE = 0,
	QAWS_CLIP_INTERSECTION,
	QAWS_CLIP_UNION,
	QAWS_CLIP_DIFFERENCE,
	QAWS_CLIP_XOR
} qaws_clip_type;

#define QAWS_CLIP_PRESERVE_COLLINEAR 1u  /* keep vertices between collinear line pieces */
#define QAWS_CLIP_REVERSE_SOLUTION   2u  /* outer paths clockwise, holes counter-clockwise */

#define QAWS_CLIP_NONE_INDEX 0xFFFFFFFFu

/* Where an output vertex comes from: the input curve it lies on, and the
   second input curve at a crossing (operand 0 subject, 1 clip, 2 open
   subject; path and curve index into that operand's array). */
typedef struct qaws_clip_vertex
{
	qaws_vec2 position;
	unsigned int operand_a, path_a, curve_a;
	qaws_scalar parameter_a;
	unsigned int operand_b, path_b, curve_b;   /* QAWS_CLIP_NONE_INDEX when on one curve only */
	qaws_scalar parameter_b;
	qaws_scalar z;                            /* from z_fn, else 0 */
} qaws_clip_vertex;

/* Called for every output vertex where two input curves meet (Clipper2's Z
   callback); returns the vertex's z. */
typedef qaws_scalar (*qaws_clip_z_fn)(void* user, qaws_clip_vertex const* vertex);

typedef struct qaws_clip_desc
{
	qaws_path_2d const* subjects;
	unsigned int subject_count;
	qaws_path_2d const* open_subjects;
	unsigned int open_subject_count;
	qaws_path_2d const* clips;
	unsigned int clip_count;
	qaws_clip_type clip_type;
	qaws_fill_rule fill_rule;
	unsigned int flags;                     /* QAWS_CLIP_* */
	qaws_scalar tolerance;                  /* points this close are one vertex; 0 = 6.4e-9 of the extent (1e-4 in float) */
	qaws_batch_executor const* executor;    /* optional, for the batched intersection */
	qaws_clip_z_fn z_fn;                    /* optional */
	void* z_user;
} qaws_clip_desc;

typedef struct qaws_clip_result qaws_clip_result;

qaws_status qaws_clip_execute(qaws_clip_desc const* desc, qaws_clip_result** out_result);
void qaws_clip_result_destroy(qaws_clip_result* result);

/* Closed paths. The views stay valid until the result is destroyed. */
unsigned int qaws_clip_result_get_path_count(qaws_clip_result const* result);
qaws_status qaws_clip_result_get_path(qaws_clip_result const* result, unsigned int index, qaws_path_2d* out_path);
/* Nesting (Clipper2's PolyTree): the enclosing path, QAWS_CLIP_NONE_INDEX at
   the top; holes are the paths at odd depth. */
unsigned int qaws_clip_result_get_parent(qaws_clip_result const* result, unsigned int index);
int qaws_clip_result_is_hole(qaws_clip_result const* result, unsigned int index);
unsigned int qaws_clip_result_get_depth(qaws_clip_result const* result, unsigned int index);

/* Open paths (pieces of the open subjects). */
unsigned int qaws_clip_result_get_open_path_count(qaws_clip_result const* result);
qaws_status qaws_clip_result_get_open_path(qaws_clip_result const* result, unsigned int index, qaws_path_2d* out_path);

/* The vertex at the start of each curve of a path (closed: open = 0). */
qaws_status qaws_clip_result_get_vertices(qaws_clip_result const* result, int open, unsigned int index,
	qaws_clip_vertex const** out_vertices, unsigned int* out_count);

/* One-call helpers (Clipper2's BooleanOp / Union / Intersect / Difference / Xor). */
qaws_status qaws_clip_boolean(qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_path_2d const* subjects, unsigned int subject_count,
	qaws_path_2d const* clips, unsigned int clip_count,
	qaws_clip_result** out_result);

#endif /* QAWS_CLIP_H */
