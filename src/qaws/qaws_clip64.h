#ifndef QAWS_CLIP64_H
#define QAWS_CLIP64_H

#include <stdint.h>
#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_clip.h"

/*
 * Exact Boolean operations on int64 polygons (Clipper2's Clipper64).
 *
 * Coordinates are integers with |x|, |y| <= QAWS_CLIP64_MAX_COORD (2^61 - 1,
 * Clipper2's range). Every predicate is exact: crossings are rational points
 * (the numerator and denominator of a segment / segment crossing fit in 256
 * and 128 bits), vertices are equal only when they are the same point, and
 * the order of edges round a vertex, the orientation of faces and the
 * containment of parts come from exact integer signs. Collinear overlaps,
 * touching vertices and vertices on edges are resolved exactly; nothing
 * depends on a tolerance.
 *
 * The result's topology is exact. Its vertices are given both rounded to
 * the nearest integer (Clipper2's output, cleaned of the repeated and
 * collinear points rounding makes) and as the double nearest each exact
 * rational point.
 *
 * Same clip types, fill rules, flags (QAWS_CLIP_PRESERVE_COLLINEAR,
 * QAWS_CLIP_REVERSE_SOLUTION), nesting and open subjects as qaws_clip.
 *
 * Thread-safe.
 */

#define QAWS_CLIP64_MAX_COORD ((int64_t)0x1FFFFFFFFFFFFFFF)

typedef struct qaws_path64
{
	int64_t const* points;      /* x, y pairs */
	unsigned int point_count;
} qaws_path64;

typedef struct qaws_clip64_desc
{
	qaws_path64 const* subjects;
	unsigned int subject_count;
	qaws_path64 const* open_subjects;
	unsigned int open_subject_count;
	qaws_path64 const* clips;
	unsigned int clip_count;
	qaws_clip_type clip_type;
	qaws_fill_rule fill_rule;
	unsigned int flags;         /* QAWS_CLIP_* */
} qaws_clip64_desc;

typedef struct qaws_clip64_result qaws_clip64_result;

qaws_status qaws_clip64_execute(qaws_clip64_desc const* desc, qaws_clip64_result** out_result);
void qaws_clip64_result_destroy(qaws_clip64_result* result);

/* Closed paths, rounded to integers (x, y pairs). */
unsigned int qaws_clip64_result_get_path_count(qaws_clip64_result const* result);
qaws_status qaws_clip64_result_get_path(qaws_clip64_result const* result, unsigned int index,
	int64_t const** out_points, unsigned int* out_count);
/* The same path's exact vertices before rounding (as the nearest doubles):
   every point where the exact boundary turns, so it may have more vertices
   than the rounded path. Paths that round to nothing (no area left) are
   dropped from both. */
qaws_status qaws_clip64_result_get_path_exact(qaws_clip64_result const* result, unsigned int index,
	double const** out_points, unsigned int* out_count);
unsigned int qaws_clip64_result_get_parent(qaws_clip64_result const* result, unsigned int index);
int qaws_clip64_result_is_hole(qaws_clip64_result const* result, unsigned int index);
unsigned int qaws_clip64_result_get_depth(qaws_clip64_result const* result, unsigned int index);

/* Open paths, rounded. */
unsigned int qaws_clip64_result_get_open_path_count(qaws_clip64_result const* result);
qaws_status qaws_clip64_result_get_open_path(qaws_clip64_result const* result, unsigned int index,
	int64_t const** out_points, unsigned int* out_count);

/* Twice the signed area of an int64 path, exactly as a double rounded once
   (Clipper2's Area times 2). */
double qaws_path64_area2(int64_t const* points, unsigned int count);

#endif /* QAWS_CLIP64_H */
