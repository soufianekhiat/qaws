# API Reference: qaws_offset.h

Offsetting paths of curves. This is Clipper2's `ClipperOffset` and
`InflatePaths`, working on curves as well as polygons.

```c
typedef enum qaws_join_type { QAWS_JOIN_SQUARE, QAWS_JOIN_BEVEL, QAWS_JOIN_ROUND, QAWS_JOIN_MITER } qaws_join_type;
typedef enum qaws_end_type  { QAWS_END_POLYGON, QAWS_END_JOINED, QAWS_END_BUTT, QAWS_END_SQUARE, QAWS_END_ROUND } qaws_end_type;

typedef struct qaws_offset_group {
	qaws_path_2d const* paths; unsigned int path_count;
	qaws_join_type join_type;  qaws_end_type end_type;
} qaws_offset_group;

typedef struct qaws_offset_desc {
	qaws_offset_group const* groups; unsigned int group_count;
	qaws_scalar delta;
	qaws_scalar miter_limit;          /* x |delta|; 0 = 2 */
	qaws_scalar tolerance;            /* curve offsets' fit; 0 = 1e-6 of max(|delta|, extent) */
	unsigned int flags;               /* QAWS_CLIP_PRESERVE_COLLINEAR, QAWS_CLIP_REVERSE_SOLUTION */
	qaws_offset_delta_fn delta_fn;    /* variable offset: delta at a point and its normal */
	void* delta_user;
	qaws_batch_executor const* executor;
} qaws_offset_desc;

qaws_status qaws_offset_execute(qaws_offset_desc const* desc, qaws_clip_result** out_result);
qaws_status qaws_offset_paths(paths, path_count, delta, join_type, end_type, out_result);   /* InflatePaths */
```

The result is a `qaws_clip_result`, with the same paths, nesting and
vertices as the Boolean operations (see [clip.md](clip.md)).

## How it works

1. **Pieces.** Each path is cut into pieces: polylines by segment,
   composites by segment, arcs by arc segment.
   - A line is translated by delta along its normal.
   - An arc becomes the concentric arc of radius r ± delta.
   - Any other curve becomes G1 cubics fitted to the true offset
     C(t) + delta N(t). See [Curves](#curves) below. Only the direction of
     C' is used, so every kind of parameterization works.
   - Under a delta callback, the offset varies along the curve. Cubics are
     then put through four offset points, and each is split until three
     probes lie within the tolerance.
2. **Corners.**
   - On the outer side of a corner, the join is placed between the two
     offsets:

     | Join | Shape |
     |---|---|
     | round | an exact arc of radius \|delta\| |
     | miter | the two offsets extended to meet; squared beyond the miter limit (Clipper2's rule `1 + cos a > 2 / ml²`) |
     | square | cut perpendicular to the bisector at exactly \|delta\| from the corner |
     | bevel | a straight line between the two offsets |

   - On the inner side, the two offsets are linked through the corner, as
     Clipper2 does. The loops this makes are removed later.
3. **Open paths.** These are offset by |delta| and travel the path forward
   and back. Butt, square and round ends become the caps at the two
   turnarounds. A joined end treats the two turnarounds as corners with
   the join type.
4. **Clean-up.** Every raw loop goes through one union with the positive
   fill rule (`qaws_clip`). That removes self-overlaps, the cusps of
   curve offsets, and anything that shrank past nothing.

**Orientation.** Positive deltas grow regions. A closed group whose
outermost path, the one with the lowest point, runs clockwise is offset
the other way. Its loops are turned round for the union, and the result
keeps the input's orientation, as Clipper2 does.

## Curves

This follows "Fast GPU stroke expansion" (Levien and Uguray, HPG 2024),
which Dear Widgets ports from Rust. The paper outputs line segments for a
GPU rasterizer; qaws outputs curves.

1. **Euler cuts.** Each smooth span is cut, in travel order, wherever it
   stops looking like an Euler spiral (clothoid). The test is:
   - the paper's error estimate on the end points and derivatives;
   - plus the gap at the middle to the Hermite cubic of those ends, for
     curve kinds that are not cubics;
   - and a check that the tangent at the middle does not turn back from
     either end.

   Between two cuts, the curvature k is then close to linear.
2. **Cusps of the offset.** The offset O = C + d n has
   O' = |C'| (1 + d k) T. It runs forward where 1 + d k > 0 and backward
   where 1 + d k < 0, with a cusp at each root. With k close to linear
   between cuts, each root lies between two cuts of opposite sign. It is
   then polished on the true curve by Illinois steps. Test 90 checks the
   cusps against Mathematica to 1e-10.
3. **Forward runs** become cubics with their ends and end tangents held,
   so neighbouring cubics meet G1. The two arm lengths come from the
   closer of two starting guesses:
   - least squares on samples;
   - the Hermite arms from the true speed, which vanishes at a cusp, so
     the cubic's end is a cusp too.

   Gauss-Newton on the samples' distances then refines the arms. A cubic
   is split in two while it misses the offset by more than the tolerance.
   A NURBS circle grown by 1 takes 16 cubics, the minimum for that
   tolerance. The old fit took 64.
4. **Backward runs** become the evolute E = C − n / k, the centres of
   curvature. The offset meets the evolute at its cusps. The chain runs:
   - along O up to the cusp,
   - along E to the next cusp,
   - then along O again.

   This is the smooth form of Clipper's link through an inner corner, and
   in the limit of a sharp corner it is exactly that link. For one side,
   the backward run alone already winds with the sign of delta, so the
   positive union would be right without the evolute. The evolute is used
   because it is smooth where the backward run ends in two cusps: a
   parabola offset toward its focus at 1e-9 takes 30 cubics, against 54
   with the backward run. E has its own cusps where k peaks, and its
   cubics are cut there too.
5. **Cusps of the curve** (C' = 0, where the tangent turns back between
   two of the finest cuts) split the run. There the offset jumps from
   C + d n to C − d n, and the true offset is the half circle between
   them, which a round join supplies. Collinear control points that double
   back, the paper's `do_linear` case, are handled the same way.
6. **Output.** Straight pieces become lines and runs of curved pieces
   become cubic B-splines. Coincident straight edges then merge exactly in
   the union.

Figure 6 of `qaws_better2d_showcase` shows the raw chains, the cusps and
the evolutes. Figure 7 shows the hard cases from the paper:
- a cusp, at two widths;
- a straight line doubling back;
- a loop;
- a near-cusp of radius 0.01;
- a tight parabola;
- a tiny S under a stroke three times its size;
- a wide inflection.

Test 90 checks each of these against a distance oracle on a grid: a point
lies in the stroke if and only if its distance to the curve is below
delta.

## Numbers

Test 87 checks these closed forms:

- **A 10 × 10 square grown by 1:**
  - miter: 144
  - square: 132 + 8√2
  - bevel: 142
  - round: 140 + π
- **The same square shrunk:** by 1 it is 64; by 6 nothing is left.
- **A square with a hole, grown by 1:** 140. The outer edge grows and the
  hole shrinks.
- **A segment of length 10, delta 1:**
  - butt: 20
  - square: 24
  - round: 20 + π
  - joined with round joins: 20 + π
- **Circles:**
  - An arc-kind circle of radius 5 grown by 1 gives 36π exactly, from
    concentric arcs.
  - A NURBS circle grown by 1 gives 36π, and shrunk by 2 gives 9π, both
    within the fit tolerance.
- **An open Bézier with round ends:** 2dL + πd².
- **A thin triangle with miter joins:** the apex is cut within 2 deltas.
- **A constant delta callback** equals the plain offset.

Test 90 checks:

- the offset's cusps against Mathematica (a parabola and both sides of a
  cubic S), each one a joint of the chain;
- the chain lies on the true offset or the evolute;
- the chain is G1 except at the cusps;
- stars grown and shrunk past their radii of curvature, and the hard cases
  of figure 7, against the distance oracle, in f64 and f32.

`qaws_better2d_showcase` (figure 4) runs 12 offset series in 0.04 s.

## int64 paths (Clipper2's ClipperOffset)

```c
typedef struct qaws_offset64_group { qaws_path64 const* paths; unsigned int path_count; qaws_join_type join_type; qaws_end_type end_type; } qaws_offset64_group;
typedef struct qaws_offset64_desc {
	qaws_offset64_group const* groups; unsigned int group_count;
	double delta;
	double miter_limit;       /* x |delta|; 0 = 2 */
	double arc_tolerance;     /* 0 = |delta| / 500 */
	unsigned int flags;       /* QAWS_CLIP_PRESERVE_COLLINEAR, QAWS_CLIP_REVERSE_SOLUTION */
} qaws_offset64_desc;
qaws_status qaws_offset64_execute(desc, qaws_clip64_result** out_result);
qaws_status qaws_offset64_paths(paths, path_count, delta, join_type, end_type, out_result);   /* InflatePaths */
```

This is Clipper2's `ClipperOffset` algorithm, ported to C (Angus Johnson,
Boost licence). The result matches Clipper2's own:
- every point is rounded to the integer grid as Clipper2 rounds it;
- round joins and ends are arcs of straight steps within the arc
  tolerance;
- at a concave corner the path runs through the vertex itself;
- miter, square and bevel joins follow Clipper2's rules, including its
  case for nearly straight joins.

The clean-up is one exact union (`qaws_clip64`): positive, or negative for
a reversed group. The result keeps the input's orientation unless the
reverse flag is set.

Test 93 runs Clipper2's TestOffsets.cpp and TestOffsetOrientation.cpp in C,
with Clipper2's own expectations:
- the Offsets.txt corpus: exactly one outer path, on int64 paths and on
  curves (`qaws_offset`);
- the regression issues #405, #424, #448, #456, #482, #593, #715, #724,
  #733, #873 and #965: vertex counts, distances within the arc tolerance,
  orientation, over-shrunk paths vanishing, and an area of exactly 122.

All 15 checks pass.
