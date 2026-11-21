# API Reference: qaws_curve_batch.h

Queries over N curves at once: closest points for many query points, on
one-shot descriptors or on prepared sets that keep the flattening. The
header also declares the types of batched curve intersections and level
crossings, shared by code built on qaws.

Precision-bearing names exist as `name_f32` and `name_f64`; this page uses
the unsuffixed names (see [precision.md](precision.md)). Thread-safe on
immutable curves. Parallel execution: [batch_executor.md](batch_executor.md).

---

## Flattening

Every curve is flattened once into chord segments whose boxes enclose the
curve, inflated by the measured chord deviation. One uniform grid over all
segments of all curves gives the candidates, and only those are refined by
Newton on the true curves. The cost grows with the total segment count plus
the number of nearby candidates, not with the number of curve pairs.

`flatness` is the chord deviation bound; 0 means 2^-10 of the scene extent.

## Batch descriptor and set

```c
#define QAWS_CURVE_BATCH_SELF 1u   /* also intersect each curve with itself */

typedef struct qaws_curve_batch_desc {
    qaws_curve const* const* curves;
    unsigned int curve_count;
    unsigned int const* families;          /* optional, one id per curve */
    unsigned int flags;                    /* QAWS_CURVE_BATCH_* */
    qaws_scalar flatness;                  /* 0 = 2^-10 of the scene extent */
    qaws_batch_executor const* executor;   /* optional */
    qaws_progress const* progress;         /* optional: cancel and progress reports */
} qaws_curve_batch_desc;
```

`families` and `flags` concern intersections: two curves with the same
family id are never intersected with each other (for a heightfield, the
contour lines in one family and the gradient lines in another), and
`QAWS_CURVE_BATCH_SELF` also intersects each curve with itself.

### qaws_curve_set_create

```c
typedef struct qaws_curve_set qaws_curve_set;

qaws_status qaws_curve_set_create(qaws_curve_set** out_set, qaws_curve_batch_desc const* desc);
```

Flattens the curves once (the costly part: curve evaluations) for curves
queried again and again. The curves must be all 2D or all 3D
(`QAWS_STATUS_INVALID_DIMENSION` otherwise); a NULL desc, curve array entry
or `out_set` gives `QAWS_STATUS_INVALID_ARGUMENT`. The set copies the desc,
its curve pointer array and its families, and keeps the executor. It keeps
pointers to the curves: they must outlive the set and stay unchanged.

### qaws_curve_set_destroy

```c
void qaws_curve_set_destroy(qaws_curve_set* set);
```

### qaws_curve_set_get_segment_count

```c
unsigned int qaws_curve_set_get_segment_count(qaws_curve_set const* set);
```

The chord segments over all curves of the set.

## Statistics

```c
typedef struct qaws_curve_batch_stats {
    unsigned int segment_count;     /* chord segments over all curves */
    unsigned int cell_count;        /* grid cells */
    unsigned int candidate_count;   /* candidates whose boxes overlap */
    unsigned int newton_count;      /* candidates close enough to refine */
    unsigned int hit_count;         /* distinct results */
} qaws_curve_batch_stats;
```

For closest points, `candidate_count` counts the segments refined and
`hit_count` the points that found a curve.

## Closest points

```c
#define QAWS_CURVE_BATCH_NONE 0xFFFFFFFFu

typedef struct qaws_closest_desc {
    qaws_curve const* const* curves;       /* all 2D or all 3D */
    unsigned int curve_count;
    qaws_scalar const* points;             /* point_count points of 2 or 3 scalars */
    unsigned int point_count;
    qaws_scalar max_distance;              /* 0 = no limit */
    qaws_scalar flatness;                  /* 0 = 2^-10 of the extent */
    qaws_batch_executor const* executor;   /* optional */
    qaws_progress const* progress;         /* optional: cancel and progress reports */
} qaws_closest_desc;

typedef struct qaws_closest_point {
    unsigned int curve;                    /* QAWS_CURVE_BATCH_NONE: no curve within max_distance */
    qaws_scalar parameter;
    qaws_scalar distance;
    qaws_vec3 position;                    /* z = 0 for 2D */
} qaws_closest_point;
```

### qaws_curve_batch_find_closest

```c
qaws_status qaws_curve_batch_find_closest(qaws_closest_point* out_points, qaws_curve_batch_stats* out_stats,
    qaws_closest_desc const* desc);
```

For every query point, the nearest point over all curves. `out_points`
receives `point_count` results; `out_stats` may be NULL.

1. One grid over the flattened segments of every curve is searched ring by
   ring outward from the point.
2. A segment is kept while its chord distance minus its inflation can beat
   the best chord distance plus inflation seen.
3. The survivors are refined by Newton on (C(t) - p) . C'(t) = 0, clamped
   to the curve's domain. The nearest true point wins; the lower curve
   index wins a tie.

With `max_distance`, points farther than it from every curve get
`curve = QAWS_CURVE_BATCH_NONE`. Returns `QAWS_STATUS_INVALID_ARGUMENT` for
a NULL desc or missing arrays, `QAWS_STATUS_INVALID_DIMENSION` for mixed
dimensions.

### qaws_curve_set_find_closest

```c
qaws_status qaws_curve_set_find_closest(qaws_closest_point* out_points, qaws_curve_batch_stats* out_stats,
    qaws_curve_set const* set, qaws_scalar const* points, unsigned int point_count, qaws_scalar max_distance);
```

The same on a prepared set, with the set's flatness and executor.

```c
qaws_status qaws_curve_set_find_closest_ex(qaws_closest_point* out_points, qaws_curve_batch_stats* out_stats,
    qaws_curve_set const* set, qaws_scalar const* points, unsigned int point_count, qaws_scalar max_distance,
    qaws_progress const* progress);
```

The same with a progress ([progress.md](progress.md); may be NULL): stage `query`. A set does not keep the progress of
its desc; `qaws_curve_batch_find_closest` reports `flatten` then `query`, `qaws_curve_set_create` `flatten`. Cancelled:
`QAWS_STATUS_CANCELLED`, no set made, the closest points unspecified.

## Intersection and level-crossing types

These types describe batched curve intersections and level crossings. No
function of this header computes them.

### Hits

```c
#define QAWS_CURVE_HIT_CROSSING 0u
#define QAWS_CURVE_HIT_TOUCH    1u
#define QAWS_CURVE_HIT_OVERLAP  2u

typedef struct qaws_curve_batch_hit_2d {   /* qaws_curve_batch_hit_3d: qaws_vec3 position */
    unsigned int curve_a;
    unsigned int curve_b;
    qaws_scalar parameter_a;
    qaws_scalar parameter_b;
    qaws_vec2 position;
    unsigned int kind;                     /* QAWS_CURVE_HIT_* */
    qaws_scalar parameter_a_end;           /* OVERLAP: other end of the shared stretch; else parameter_a */
    qaws_scalar parameter_b_end;           /* OVERLAP: other end on curve_b; else parameter_b */
} qaws_curve_batch_hit_2d;
```

- **CROSSING** passes from one side to the other (transversal, or tangent
  with an odd contact such as at an inflection).
- **TOUCH** meets without crossing (an even tangency).
- **OVERLAP** is a stretch the two curves share: `[parameter_a,
  parameter_a_end]` on `curve_a` (increasing) and `[parameter_b,
  parameter_b_end]` on `curve_b` (decreasing when they run opposite ways),
  from `position` to the end point.
- 3D hits are CROSSING or OVERLAP.
- Ordering: `curve_a < curve_b`, or `curve_a == curve_b` with
  `parameter_a < parameter_b` for a self-intersection; sorted by
  `(curve_a, curve_b, parameter_a)`.

### Level crossings

Where N curves cross the level sets h = L of a scalar field, for every
level at once. Each segment's field values at its ends bracket the levels
it crosses (a segment is split where h is not close to linear along it);
each bracket is solved for h(C(t)) = L by Newton with bisection safeguards
(secant without a gradient). A curve touching a level without crossing it
is not reported.

```c
typedef qaws_scalar (*qaws_scalar_field_fn)(
    void* user,
    qaws_scalar const* point,       /* x, y (, z) */
    qaws_scalar* gradient);         /* NULL, or receives dh/dx, dh/dy (, dh/dz) */

typedef struct qaws_level_crossing_desc {
    qaws_curve const* const* curves;       /* all 2D or all 3D */
    unsigned int curve_count;
    qaws_scalar_field_fn field;
    void* user;
    qaws_scalar const* levels;             /* strictly ascending */
    unsigned int level_count;
    qaws_scalar flatness;                  /* 0 = 2^-10 of the extent */
    qaws_batch_executor const* executor;   /* optional */
    qaws_progress const* progress;         /* optional: cancel and progress reports */
} qaws_level_crossing_desc;

typedef struct qaws_level_crossing {       /* sorted by (curve, parameter) */
    unsigned int curve;
    unsigned int level;
    qaws_scalar parameter;
    qaws_vec3 position;                    /* z = 0 for 2D */
} qaws_level_crossing;
```
