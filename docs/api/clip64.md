# API Reference: qaws_clip64.h

Exact Boolean operations on int64 polygons. This is Clipper2's `Clipper64`,
with exact topology.

```c
#define QAWS_CLIP64_MAX_COORD ((int64_t)0x1FFFFFFFFFFFFFFF)   /* 2^61 - 1, Clipper2's range */

typedef struct qaws_path64 { int64_t const* points; unsigned int point_count; } qaws_path64;

typedef struct qaws_clip64_desc {
	qaws_path64 const* subjects;      unsigned int subject_count;
	qaws_path64 const* open_subjects; unsigned int open_subject_count;
	qaws_path64 const* clips;         unsigned int clip_count;
	qaws_clip_type clip_type;         /* as qaws_clip */
	qaws_fill_rule fill_rule;
	unsigned int flags;               /* QAWS_CLIP_PRESERVE_COLLINEAR, QAWS_CLIP_REVERSE_SOLUTION */
} qaws_clip64_desc;

qaws_status qaws_clip64_execute(qaws_clip64_desc const* desc, qaws_clip64_result** out_result);
```

Closed paths are closed implicitly; don't repeat the first point.
Coordinates outside the range return `QAWS_STATUS_OUT_OF_RANGE`.

## Exactness

Clipper2's orientation tests are exact, but it rounds every new vertex to
the integer grid and then repairs the damage with tolerances. Here nothing
is rounded before the topology is fixed:

- **Crossings** are rational points `(a·den + r·tn) / den`. Here `den` and
  `tn` are the exact cross products of the segments, below 2^126. The
  numerators stay below 2^190.
- **Collinear segments** cut each other at the end points that lie inside
  the other one, and straight edges between the same two vertices are one
  edge. A shared stretch is therefore one edge, with the summed winding
  steps of every segment along it.
- **Vertices** are equal only when they are the same rational point.
  Points are sorted by exact comparison, with a double prefilter.
- **The edges around a vertex** are ordered by half-plane, then by the
  exact sign of the cross product of the integer segment directions.
- **A cycle's orientation** is its turn at its lowest vertex.
- **Containment of separate parts** is an exact crossing count.
- **Arithmetic** is on 384-bit integers (`internal/qaws_internal_wide.h`).
  Cross-product signs, rounding and the segment-pair filter decide in
  double when the error bound allows and fall back to the exact value
  otherwise.

Faces, windings, clip types, fill rules, loops cut at touch points, nesting
and open subjects all work as in [clip.md](clip.md).

## Result

```c
qaws_clip64_result_get_path(r, i, &points, &count);         /* rounded to the nearest integer */
qaws_clip64_result_get_path_exact(r, i, &xy, &count);       /* exact vertices as the nearest doubles */
qaws_clip64_result_get_parent / is_hole / get_depth / get_open_path ...
double qaws_path64_area2(int64_t const* points, unsigned int count);   /* exact 2 x area, rounded once */
```

Output vertices are rounded to the nearest integer. Clipper2 truncates
instead. Rounding is the only inexact step, and it can leave repeated
points and new collinear corners in a rounded path:

- Repeated points are dropped.
- Collinear corners are dropped unless `QAWS_CLIP_PRESERVE_COLLINEAR` is
  set. Spikes are always dropped.
- A path that rounds to no area is dropped.

The exact vertices, before any of this, are available alongside.

## Parity and numbers

Test 85 runs Clipper2's `Polygons.txt`, 195 records, through this engine:

- **Areas** match Clipper2 to within its own tolerances, or within the
  perimeter of the result. That bound covers the different rounding.
- **Path counts** lie between ours with touching pieces joined and the
  exact engine's count. Pieces that touch at a point and slivers under one
  unit are where Clipper2's rounding and ours part ways.
- **An oracle** samples a grid of points per record and checks the result
  against the inputs' fill rule and clip type. It finds 0 wrong points out
  of 110,016.
- **Record 62** is the one exception. Clipper2 stores 2 paths there, while
  the exact union is one polygon.

Test 86 checks:

- every clip type;
- shared edges, corner touches and identical inputs;
- rectangles at 2^60 overlapping on half their edges;
- Clipper2's #831 bow tie at 0x4000000000000;
- fill rules, nesting and open subjects.

Timings come from `qaws_clipper2_compare` (`bench/clipper2_compare.cpp`,
built with `-DQAWS_BUILD_CLIPPER2_BENCH=ON -DQAWS_CLIPPER2_DIR=...`). They
are the best of 3 runs in a Release build on one thread:

| Case | Clipper2 | qaws_clip64 | qaws_clip (polylines) |
|---|---|---|---|
| 2 random polygons of 100 vertices, intersection | 0.0007 s | 0.012 s | 0.014 s |
| 2 random polygons of 400 vertices | 0.015 s | 0.35 s | 0.34 s |
| union of 2000 squares | 0.005 s | 0.046 s | 0.068 s |
| union of 8000 squares | 0.021 s | 0.34 s | 0.32 s |
| Clipper2's 195 test records | 0.026 s | 0.33 s | 0.54 s |

Both of our engines are 9 to 23 times slower than Clipper2. They give the
same areas, and their topology is exact rather than repaired.

The exact engine builds every crossing's rational point and every cut's
exact parameter as soon as it finds them. Most of the time goes there,
and into sorting 384-bit records. Computing the exact values only when a
double comparison cannot decide is the next step.
