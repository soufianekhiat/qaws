#ifndef QAWS_OFFSET_H
#define QAWS_OFFSET_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_clip.h"

/*
 * Offsetting paths of curves (Clipper2's ClipperOffset / InflatePaths, on
 * curves).
 *
 * Each path is offset piece by piece: a line by translation, an arc by a
 * concentric arc (both exact), any other curve by G1 cubics fitted to the
 * true offset C(t) + d N(t) within the tolerance. Where the offset turns
 * back (1 + d k < 0, k the curvature) its cusps are found exactly and the
 * backward run is replaced by the evolute (the centres of curvature), after
 * "Fast GPU stroke expansion" (Levien and Uguray, HPG 2024); a variable
 * delta (delta_fn) fits cubics through offset points instead.
 * Corners get the join type on their outer side; on the inner side the two
 * offsets are linked through the corner (as Clipper2 does) and the loops
 * that makes are removed by the final union. Open paths get the end type.
 * Every raw offset loop then goes through one union with the positive fill
 * rule (qaws_clip), which removes the self-overlaps, cusps and inner loops.
 *
 * Positive deltas grow regions: the offset is to the right of travel for a
 * group whose outermost path (the one with the lowest point) runs
 * counter-clockwise, to the left otherwise. Open paths are offset by |delta|.
 *
 * Thread-safe on immutable curves.
 */

typedef enum qaws_join_type
{
	QAWS_JOIN_SQUARE = 0,   /* corners cut square at exactly |delta| from the corner */
	QAWS_JOIN_BEVEL,        /* corners cut by the line between the two offsets */
	QAWS_JOIN_ROUND,        /* a circular arc of radius |delta| round the corner */
	QAWS_JOIN_MITER         /* the offsets extended to meet, squared beyond the miter limit */
} qaws_join_type;

typedef enum qaws_end_type
{
	QAWS_END_POLYGON = 0,   /* closed path, offset on one side */
	QAWS_END_JOINED,        /* open path offset on both sides, ends joined with the join type */
	QAWS_END_BUTT,          /* open path, ends cut square at the end points */
	QAWS_END_SQUARE,        /* open path, ends extended by |delta| and cut square */
	QAWS_END_ROUND          /* open path, half circles at the ends */
} qaws_end_type;

typedef struct qaws_offset_group
{
	qaws_path_2d const* paths;
	unsigned int path_count;
	qaws_join_type join_type;
	qaws_end_type end_type;
} qaws_offset_group;

/* Variable offset (Clipper2's delta callback): the delta at a point of path
   `path` of group `group`, with the unit normal of the offset side. */
typedef qaws_scalar (*qaws_offset_delta_fn)(void* user, unsigned int group, unsigned int path,
	qaws_vec2 point, qaws_vec2 normal);

typedef struct qaws_offset_desc
{
	qaws_offset_group const* groups;
	unsigned int group_count;
	qaws_scalar delta;
	qaws_scalar miter_limit;          /* in multiples of |delta|; 0 = 2 (Clipper2's default) */
	qaws_scalar tolerance;            /* curve offsets' fit; 0 = 1e-6 of max(|delta|, extent) */
	unsigned int flags;               /* QAWS_CLIP_PRESERVE_COLLINEAR, QAWS_CLIP_REVERSE_SOLUTION */
	qaws_offset_delta_fn delta_fn;    /* optional: replaces delta */
	void* delta_user;
	qaws_batch_executor const* executor;
} qaws_offset_desc;

/* The offset region as a qaws_clip result (paths, nesting). */
qaws_status qaws_offset_execute(qaws_offset_desc const* desc, qaws_clip_result** out_result);

/* One group, one delta (Clipper2's InflatePaths). */
qaws_status qaws_offset_paths(qaws_path_2d const* paths, unsigned int path_count, qaws_scalar delta,
	qaws_join_type join_type, qaws_end_type end_type, qaws_clip_result** out_result);

#endif /* QAWS_OFFSET_H */
