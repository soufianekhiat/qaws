# API Reference: qaws_curve_batch.h

Intersections of N curves at once.

---

## Why a batch

Intersecting N curves pair by pair re-samples both curves for every pair and
costs about N² / 2 pair tests. The batch shares the work:

1. **Flatten once.** Every curve, knot span by knot span, becomes chord
   segments. A piece is split while its midpoint or a quarter point lies
   farther than the flatness bound from its chord. Each segment's box is its
   chord box inflated by twice the measured deviation.
2. **One grid.** A uniform grid over every segment of every curve, sized to
   about one segment per cell and filled by counting sort.
3. **Candidates.** Inside each cell, segment pairs from curves allowed to meet
   whose boxes overlap. Each pair is handled once, in the cell holding the
   low corner of the two boxes' overlap.
4. **Refine.** Chords closer than their two inflations seed Newton on the true
   curves (Gauss-Newton in 3D), polished to rounding level.
5. **Merge.** Coincident hits of one curve pair are kept once. Output is sorted
   by `(curve_a, curve_b, parameter_a)`.

The cost grows with the total segment count plus the number of nearby segment
pairs, not with the number of curve pairs.

## Families

`families` gives each curve an id, and curves with the same id are never
intersected with each other. For a heightfield, give the contour lines one
family and the gradient lines another:

- contours of different levels never cross;
- gradient lines meet only at critical points, where they are tangent.

Only the contour × gradient pairs are then searched.

## Types

```c
#define QAWS_CURVE_BATCH_SELF 1u   /* also intersect each curve with itself */

typedef struct qaws_curve_batch_desc {
	qaws_curve const* const* curves;
	unsigned int curve_count;
	unsigned int const* families;   /* optional, one id per curve */
	unsigned int flags;             /* QAWS_CURVE_BATCH_* */
	qaws_scalar flatness;           /* chord deviation bound; 0 = 2^-10 of the scene extent */
} qaws_curve_batch_desc;

typedef struct qaws_curve_batch_hit_2d {   /* _3d: qaws_vec3 position */
	unsigned int curve_a, curve_b;          /* curve_a < curve_b, or equal for a self-intersection */
	qaws_scalar parameter_a, parameter_b;   /* parameter_a < parameter_b when curve_a == curve_b */
	qaws_vec2 position;
} qaws_curve_batch_hit_2d;

typedef struct qaws_curve_batch_stats {
	unsigned int segment_count, cell_count, candidate_count, newton_count, hit_count;
} qaws_curve_batch_stats;
```

## Functions

```c
qaws_status qaws_curve_batch_find_intersections_2d(
	qaws_curve_batch_desc const* desc,
	qaws_curve_batch_hit_2d* out_hits, unsigned int hit_capacity,
	unsigned int* out_count, qaws_curve_batch_stats* out_stats);

qaws_status qaws_curve_batch_find_intersections_3d(
	qaws_curve_batch_desc const* desc,
	qaws_curve_batch_hit_3d* out_hits, unsigned int hit_capacity,
	unsigned int* out_count, qaws_curve_batch_stats* out_stats);
```

The functions write `min(count, hit_capacity)` hits. `*out_count` is the total
found, so compare it with the capacity and call again with a larger buffer.
`out_stats` may be NULL.

Every curve must have the function's dimension, or the call returns
`QAWS_STATUS_INVALID_DIMENSION`. The functions are thread-safe on immutable
curves.

Tangential contacts are not reported, because Newton's Jacobian is singular
there. The pairwise `qaws_curve_find_intersections_*` calls have the same
limit.

## Numbers

These are f64 timings from test 79 and `qaws_batch_showcase` (Release build).

| Scene | Batch | All pairs | Speedup |
|---|---|---|---|
| 16 × 16 conic contours × polynomial gradient lines | 0.002 s, 512 hits | 0.025 s, 512 hits | ×12 |
| 128 × 128 of the same | 0.09 s, 32768 hits | 1.3 s, 27498 hits (misses) | ×14 |
| Gaussian heightfield, 72 B-spline curves | 0.010 s, 375 hits | 3.2 s, 375 hits | ×317 |
| Gaussian heightfield, 143 B-spline curves | 0.032 s, 1446 hits | 15.6 s, 1446 hits | ×489 |

## Certified version

`qaws_exact_curve_batch_hits` (see [exact.md](exact.md)) does the same
for exact curves. It uses the same grid over sound span boxes, then
certifies every candidate span pair. On the 72-curve heightfield it takes
0.23 s, where pairwise exact calls take 2.1 s.
