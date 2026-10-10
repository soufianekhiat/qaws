# API Reference: qaws_rectclip.h

Rectangle clipping and Minkowski sums. These are Clipper2's `RectClip`,
`RectClipLines`, `MinkowskiSum` and `MinkowskiDiff`, on both Boolean
engines.

```c
/* rect = { left, bottom, right, top } */
qaws_status qaws_rect_clip_2d(qaws_scalar const rect[4], qaws_path_2d const* paths, unsigned int n, qaws_clip_result** out);
qaws_status qaws_rect_clip_lines_2d(qaws_scalar const rect[4], qaws_path_2d const* paths, unsigned int n, qaws_clip_result** out);
qaws_status qaws_rect_clip64(int64_t const rect[4], qaws_path64 const* paths, unsigned int n, qaws_clip64_result** out);
qaws_status qaws_rect_clip_lines64(int64_t const rect[4], qaws_path64 const* paths, unsigned int n, qaws_clip64_result** out);

qaws_status qaws_minkowski_sum_2d(pattern, pattern_count, path, path_count, closed, qaws_clip_result** out);
qaws_status qaws_minkowski_diff_2d(...);
qaws_status qaws_minkowski_sum64(int64_t const* pattern, ..., qaws_clip64_result** out);
qaws_status qaws_minkowski_diff64(...);
```

## Rectangle clipping

- **Closed paths** are clipped to the rectangle one by one, as Clipper2's
  `RectClip` does. Two overlapping paths stay two, and the results come out
  in path order in one result.
- **Open paths** are cut to the pieces inside the rectangle, which become
  the open paths of the result.
- **Paths of any curve kind** can be clipped in the float version. A
  circle clipped at a corner gives exactly a quarter of its arc, for
  example.
- **A path that crosses itself** has its overlap resolved by the non-zero
  rule. Clipper2's `RectClip` leaves self-intersections in place.

## Minkowski

The pattern is a closed polygon, swept along a polyline path. For each
edge of the path and each edge of the pattern, the quad between the two
translated copies of the pattern edge is made counter-clockwise. The result
is the union of all the quads under the non-zero rule:

- the sum translates by +pattern, the difference by −pattern;
- a closed path's sweep is the band around its boundary, as in Clipper2.

The int64 versions build the quads exactly, so sums of coordinates must
stay within `QAWS_CLIP64_MAX_COORD`.

## Numbers

Test 88 checks:

- **A square, a square and a circle clipped to [0, 10] × [0, 5]:** three
  paths with area 25 + 18 + π.
- **A zig-zag clipped by the rectangle's lines:** one piece inside.
- **int64:** area 43, exactly.
- **A square swept along a square path:** the band 12² − 8² = 80, in
  float and exactly in int64.
- **A square swept along a segment:** 12 × 2.
- **A triangle's sum and difference along a segment:** equal areas.
