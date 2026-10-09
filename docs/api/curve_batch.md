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

## Prepared sets

```c
qaws_status qaws_curve_set_create(qaws_curve_batch_desc const* desc, qaws_curve_set** out_set);
void qaws_curve_set_destroy(qaws_curve_set* set);
unsigned int qaws_curve_set_get_segment_count(qaws_curve_set const* set);

qaws_status qaws_curve_set_find_intersections_2d(qaws_curve_set const* set, qaws_curve_set const* other,
	qaws_curve_batch_hit_2d* out_hits, unsigned int hit_capacity, unsigned int* out_count, qaws_curve_batch_stats* out_stats);
/* and _3d */
```

A set flattens its curves once. A query re-runs only the grid and the
narrow phase:

- `other == NULL` intersects the set with itself, using its families and
  flags, exactly as the one-shot call.
- Otherwise every curve of `set` is intersected with every curve of
  `other`. `curve_a` indexes `set` and `curve_b` indexes `other`.

The set keeps pointers to its curves, so they must outlive it unchanged.

Example: 64 fixed contours queried against 6 frames of 64 gradient lines
give the same 49152 hits as the one-shot calls, in 0.085 s instead of
0.097 s at equal flatness. Newton on the crossings dominates there. The
saving grows with the size and span count of the fixed set.

## Level crossings

```c
typedef qaws_scalar (*qaws_scalar_field_fn)(void* user, qaws_scalar const* point, qaws_scalar* gradient);

typedef struct qaws_level_crossing_desc {
	qaws_curve const* const* curves; unsigned int curve_count;   /* all 2D or all 3D */
	qaws_scalar_field_fn field; void* user;
	qaws_scalar const* levels; unsigned int level_count;         /* strictly ascending */
	qaws_scalar flatness;
} qaws_level_crossing_desc;

typedef struct qaws_level_crossing {
	unsigned int curve, level;
	qaws_scalar parameter;
	qaws_vec3 position;
} qaws_level_crossing;            /* sorted by (curve, parameter) */

qaws_status qaws_curve_batch_find_level_crossings(qaws_level_crossing_desc const* desc,
	qaws_level_crossing* out_crossings, unsigned int capacity, unsigned int* out_count);
```

This finds where curves cross the level sets h = L of a scalar field, for
every level at once, without building the level curves. For a heightfield,
the curves are the gradient lines and the levels are the contour heights,
so each crossing of a gradient line with a contour is found directly.

1. Each curve is flattened once.
2. The field values at a segment's ends bracket the levels it crosses (a
   binary search in the sorted levels). A segment is split while h is not
   close to linear along it, so a level crossed twice inside it is not
   missed. A value exactly on a level belongs to the piece where it is the
   lower end, so each crossing is counted once.
3. Each bracket is solved for h(C(t)) = L. With the gradient this uses
   Newton (gradient · C'); without it, Illinois secant. Both are kept inside
   the bracket.

The field may leave `gradient` untouched when it has none. A curve that
touches a level without crossing it is not reported.

Example: 128 parabola gradient lines against 128 levels of x² + 2y² give
the 32768 crossings of the batch with contour curves, in 0.010 s instead of
0.074 s. On the Gaussian heightfield of `qaws_batch_showcase` they are
3 to 5 times faster than the batch against the traced contour curves.

## Closest points

```c
typedef struct qaws_closest_desc {
	qaws_curve const* const* curves; unsigned int curve_count;   /* all 2D or all 3D */
	qaws_scalar const* points; unsigned int point_count;         /* 2 or 3 scalars per point */
	qaws_scalar max_distance;                                    /* 0 = no limit */
	qaws_scalar flatness;
} qaws_closest_desc;

typedef struct qaws_closest_point {
	unsigned int curve;               /* QAWS_CURVE_BATCH_NONE: none within max_distance */
	qaws_scalar parameter, distance;
	qaws_vec3 position;
} qaws_closest_point;

qaws_status qaws_curve_batch_find_closest(qaws_closest_desc const* desc, qaws_closest_point* out_points,
	qaws_curve_batch_stats* out_stats);
qaws_status qaws_curve_set_find_closest(qaws_curve_set const* set, qaws_scalar const* points, unsigned int point_count,
	qaws_scalar max_distance, qaws_closest_point* out_points, qaws_curve_batch_stats* out_stats);
```

This returns, for every query point, the nearest point over all curves.

1. One grid over the flattened segments of every curve is searched ring by
   ring outward from the point.
2. A segment is kept while its chord distance minus its inflation can beat
   the best chord distance plus inflation seen so far. The search stops
   when the next ring is farther than that bound.
3. The surviving segments are refined by Newton on (C(t) − p) · C'(t) = 0,
   clamped to the curve's domain. The nearest true point wins, and the
   lower curve index wins a tie.

With `max_distance`, points farther than it from every curve get
`QAWS_CURVE_BATCH_NONE`.

Example: 2000 points against 40 concentric NURBS circles take 0.004 s,
against about 0.31 s for `qaws_curve_find_closest_parameter_2d` on every
curve. The distances match the closed form to about 1e-15 in 2D and 3D.

## Parallel execution

```c
typedef void (*qaws_batch_task_fn)(void* ctx, unsigned int begin, unsigned int end);
typedef struct qaws_batch_executor {
	void (*parallel_for)(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx);
	void* user;
} qaws_batch_executor;
```

qaws has no threads of its own. Every batch desc ends with an optional
`qaws_batch_executor const* executor`. Prepared sets copy it from the desc
they are created with.

- The batch splits its independent work into at most 256 chunks it owns.
  The work is query points, rays, or segments for level crossings.
- The executor's `parallel_for` must call `task(ctx, begin, end)` over
  ranges covering `[0, count)` exactly once, on any threads in any order,
  and return when all are done.
- Each chunk has its own scratch memory and statistics, so the results are
  identical whatever the executor does.
- Without an executor, the work runs on the calling thread.

It currently covers closest points on curves and surfaces, ray casting
and level crossings.

Test 81 runs 20000 surface closest points, 20000 rays, 20000 curve closest
points and level crossings. A scrambled single-threaded executor and a
4-thread executor (Win32 threads, or pthreads) both give the serial
results bit for bit, the threads 3.3 times faster.
