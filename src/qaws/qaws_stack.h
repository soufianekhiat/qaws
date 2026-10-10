#ifndef QAWS_STACK_H
#define QAWS_STACK_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_clip.h"
#include "qaws_batch_executor.h"

/*
 * 2.5D shapes as stacks of parallel contours.
 *
 * A stack is a list of levels at increasing heights z (any spacing), each a
 * region (closed paths read with the fill rule). Between two levels the
 * shape is interpolated; outside the first and last level it is empty. A
 * cone is a circle at z = 0 and a point at the apex; a sphere, 64 circles; a
 * height map, its contour regions { h >= z }.
 *
 * Interpolation between levels k and k + 1 (s = (z - z_k) / (z_k+1 - z_k)):
 *
 *   MATCH     contours paired one to one between the levels and blended:
 *             C(t) = (1 - s) C_k(t) + s C_k+1(t). Contours of the same
 *             structure (kind, degree, knots) blend control points and
 *             weights, exactly: two concentric NURBS circles give the circle
 *             of the blended radius. A point contour (one degenerate curve)
 *             takes the other contour's structure, so a circle and a point
 *             give a cone exactly. Others are resampled by arc length (start
 *             and direction aligned) into closed Catmull-Rom curves.
 *   DISTANCE  the section is where (1 - s) d_k + s d_k+1 <= 0, d the signed
 *             distance to each level's region (negative inside): handles
 *             contours that split, merge or appear between levels. Traced
 *             by marching squares on a grid and refined onto the zero set;
 *             the result is polylines.
 *   AUTO      MATCH when both levels have the same number of contours, each
 *             one curve, else DISTANCE.
 *
 * LINEAR interpolation blends the two levels around z; CUBIC blends four
 * with non-uniform Catmull-Rom weights (tangents from the neighbouring
 * levels), so a sphere given by circles follows its curved profile.
 *
 * Thread-safe on immutable curves.
 */

typedef enum qaws_stack_interp
{
	QAWS_STACK_LINEAR = 0,
	QAWS_STACK_CUBIC
} qaws_stack_interp;

typedef enum qaws_stack_blend
{
	QAWS_STACK_BLEND_AUTO = 0,
	QAWS_STACK_BLEND_MATCH,
	QAWS_STACK_BLEND_DISTANCE
} qaws_stack_blend;

typedef struct qaws_stack_level
{
	qaws_scalar z;
	qaws_path_2d const* paths;
	unsigned int path_count;
} qaws_stack_level;

typedef struct qaws_stack_2d
{
	qaws_stack_level const* levels;     /* increasing z */
	unsigned int level_count;
	qaws_fill_rule fill_rule;
	qaws_stack_interp interp;
	qaws_stack_blend blend;
	unsigned int grid;                  /* DISTANCE: marching squares cells across; 0 = 192 */
} qaws_stack_2d;

/* The region of a stack at height z (nothing outside its levels). */
qaws_status qaws_stack_section_2d(qaws_stack_2d const* stack, qaws_scalar z, qaws_clip_result** out_result);

/* A Boolean of two stacks: one level at every height either stack has a
   level at (sorted, repeats once), each the clip type of the two sections
   there, solved level by level (in parallel through the executor). */
typedef struct qaws_stack_result qaws_stack_result;

qaws_status qaws_stack_boolean_2d(qaws_clip_type clip_type, qaws_stack_2d const* a, qaws_stack_2d const* b,
	qaws_batch_executor const* executor, qaws_stack_result** out_result);
void qaws_stack_result_destroy(qaws_stack_result* result);
unsigned int qaws_stack_result_get_level_count(qaws_stack_result const* result);
qaws_scalar qaws_stack_result_get_z(qaws_stack_result const* result, unsigned int level);
qaws_clip_result const* qaws_stack_result_get_level(qaws_stack_result const* result, unsigned int level);

#endif /* QAWS_STACK_H */
