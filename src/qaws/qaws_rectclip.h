#ifndef QAWS_RECTCLIP_H
#define QAWS_RECTCLIP_H

#include <stdint.h>
#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_clip.h"
#include "qaws_clip64.h"

/*
 * Rectangle clipping and Minkowski sums (Clipper2's RectClip,
 * RectClipLines, MinkowskiSum, MinkowskiDiff).
 *
 * Rectangles: rect = { left, bottom, right, top }. Closed paths are clipped
 * to the rectangle each on its own (as Clipper2's RectClip: a path's
 * self-overlaps are not merged with other paths', and the even-odd rule is
 * not applied across paths); open paths become the pieces inside.
 *
 * Minkowski: the pattern (a closed polygon) swept along the path: for each
 * edge of the path and each edge of the pattern, the quad between the two
 * translated copies, and the union of all quads (non-zero). Sum translates
 * by +pattern, difference by -pattern. A closed path's sweep is the band
 * round its boundary (Clipper2's convention).
 *
 * Thread-safe.
 */

/* Each closed path clipped to the rectangle; results appended in path
   order into one result. */
qaws_status qaws_rect_clip_2d(qaws_scalar const rect[4], qaws_path_2d const* paths, unsigned int path_count,
	qaws_clip_result** out_result);

/* Open paths: the pieces inside the rectangle (open paths of the result). */
qaws_status qaws_rect_clip_lines_2d(qaws_scalar const rect[4], qaws_path_2d const* paths, unsigned int path_count,
	qaws_clip_result** out_result);

qaws_status qaws_rect_clip64(int64_t const rect[4], qaws_path64 const* paths, unsigned int path_count,
	qaws_clip64_result** out_result);
qaws_status qaws_rect_clip_lines64(int64_t const rect[4], qaws_path64 const* paths, unsigned int path_count,
	qaws_clip64_result** out_result);

/* Minkowski sum / difference of a pattern polygon (pattern_count points) and
   a polyline path (path_count points, closed or open). */
qaws_status qaws_minkowski_sum_2d(qaws_scalar const* pattern, unsigned int pattern_count,
	qaws_scalar const* path, unsigned int path_count, int closed, qaws_clip_result** out_result);
qaws_status qaws_minkowski_diff_2d(qaws_scalar const* pattern, unsigned int pattern_count,
	qaws_scalar const* path, unsigned int path_count, int closed, qaws_clip_result** out_result);
qaws_status qaws_minkowski_sum64(int64_t const* pattern, unsigned int pattern_count,
	int64_t const* path, unsigned int path_count, int closed, qaws_clip64_result** out_result);
qaws_status qaws_minkowski_diff64(int64_t const* pattern, unsigned int pattern_count,
	int64_t const* path, unsigned int path_count, int closed, qaws_clip64_result** out_result);

#endif /* QAWS_RECTCLIP_H */
