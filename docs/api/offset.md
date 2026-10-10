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
   - Any other curve becomes a cubic B-spline through the true offset
     C(t) + delta N(t), in pieces. Each piece is the cubic through four
     offset points, split until three probes lie within the tolerance.
     Only the direction of C' is used, so every kind of parameterization
     works.
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

`qaws_better2d_showcase` (figure 4) runs 12 offset series in 0.04 s.
