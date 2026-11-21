# API Reference: qaws_surface_batch.h

Closest points over M surfaces at once, prepared surface sets, and the
descriptor and result types of surface batches.

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)). All functions are thread-safe on
immutable surfaces. The optional `executor` fields take a
`qaws_batch_executor` ([batch_executor.md](batch_executor.md)) that runs the
work in parallel; NULL runs it on the calling thread.
`QAWS_CURVE_BATCH_NONE` comes from qaws_curve_batch.h
([curve_batch.md](curve_batch.md)).

---

## Patches and grid

- Every surface is flattened into (u, v) patches, each with a box inflated
  by its measured deviation from the surface.
- `flatness` is the deviation bound; 0 means 2^-9 of the scene extent.
- One uniform grid over all patches answers the queries.

## Prepared surface sets

```c
typedef struct qaws_surface_batch_desc {
	qaws_surface const* const* surfaces;
	unsigned int surface_count;
	unsigned int const* families;          /* optional, one id per surface */
	qaws_scalar flatness;                  /* deviation bound; 0 = 2^-9 of the scene extent */
	qaws_batch_executor const* executor;   /* optional */
} qaws_surface_batch_desc;

typedef struct qaws_surface_set qaws_surface_set;   /* opaque */

qaws_status  qaws_surface_set_create(qaws_surface_set** out_set, qaws_surface_batch_desc const* desc);
void         qaws_surface_set_destroy(qaws_surface_set* set);
unsigned int qaws_surface_set_get_patch_count(qaws_surface_set const* set);
```

A set does the flattening once, for surfaces queried again and again (a
terrain against changing points). It copies the `surfaces` and `families`
arrays but keeps pointers to the surfaces themselves: they must outlive the
set and stay unchanged. Its `executor` is used by every query on the set.
`families` groups surfaces for pair queries over a set against itself
(surfaces of one family are not paired); against another set, `surface_a`
indexes the first set and `surface_b` the other.

- `qaws_surface_set_create`: `QAWS_STATUS_INVALID_ARGUMENT` on a NULL
  `desc`/`out_set` or a NULL surface, `QAWS_STATUS_ALLOCATION_FAILURE` when
  memory runs out. `*out_set` is NULL on failure.
- `qaws_surface_set_destroy`: NULL is accepted.
- `qaws_surface_set_get_patch_count`: the number of patches; 0 for NULL.

## Closest points on surfaces

```c
typedef struct qaws_surface_batch_closest_desc {
	qaws_surface const* const* surfaces;
	unsigned int surface_count;
	qaws_scalar const* points;             /* point_count points of 3 scalars */
	unsigned int point_count;
	qaws_scalar max_distance;              /* 0 = no limit */
	qaws_scalar flatness;                  /* 0 = 2^-9 of the extent */
	qaws_batch_executor const* executor;   /* optional */
} qaws_surface_batch_closest_desc;

typedef struct qaws_surface_batch_closest {
	unsigned int surface;     /* QAWS_CURVE_BATCH_NONE: none within max_distance */
	qaws_scalar u, v;
	qaws_scalar distance;
	qaws_vec3 position;
} qaws_surface_batch_closest;

typedef struct qaws_surface_batch_stats {
	unsigned int segment_count;     /* curve chord segments */
	unsigned int patch_count;       /* surface patches */
	unsigned int cell_count;        /* grid cells */
	unsigned int candidate_count;   /* piece pairs whose boxes overlap */
	unsigned int newton_count;      /* candidates refined */
	unsigned int hit_count;
} qaws_surface_batch_stats;

qaws_status qaws_surface_batch_find_closest(
	qaws_surface_batch_closest* out_points,
	qaws_surface_batch_stats* out_stats,
	qaws_surface_batch_closest_desc const* desc);

qaws_status qaws_surface_set_find_closest(
	qaws_surface_batch_closest* out_points,
	qaws_surface_batch_stats* out_stats,
	qaws_surface_set const* set,
	qaws_scalar const* points,
	unsigned int point_count,
	qaws_scalar max_distance);
```

For every query point, the nearest point over all surfaces.
`out_points` receives `point_count` results. `out_stats` may be NULL; here
`candidate_count` counts the patches refined and `hit_count` the points that
found a surface.

1. The grid is searched ring by ring outward from the point.
2. A patch is refined while its box can beat the nearest true surface point
   seen so far (patch corners, then refined points).
3. Refinement is Newton on (S - p) . S_u = (S - p) . S_v = 0, clamped to the
   domain, seeded at the nearest point of the patch's two triangles.
4. The nearest wins; on a tie, the lower surface index.

Local minima inside one patch that the seed does not reach are not
reported; the flatness bound keeps patches close to planar. For a single
surface, see `qaws_surface_find_closest_point` in [inspect.md](inspect.md).

`QAWS_STATUS_INVALID_ARGUMENT`: a NULL `desc` or `set`, NULL `surfaces`
with a nonzero count, a NULL surface, or NULL `points`/`out_points` with a
nonzero `point_count`.

## Other batch types

These descriptor and result types are declared in this header for batched
curve-surface intersections, surface-surface intersection curves and ray
casts. No function in qaws_surface_batch.h takes them.

```c
typedef struct qaws_curve_surface_batch_desc {
	qaws_curve const* const* curves;       /* 3D */
	unsigned int curve_count;
	qaws_surface const* const* surfaces;
	unsigned int surface_count;
	qaws_scalar flatness;                  /* 0 = 2^-9 of the scene extent */
	qaws_batch_executor const* executor;   /* optional */
} qaws_curve_surface_batch_desc;

/* sorted by (curve, surface, t) */
typedef struct qaws_curve_surface_batch_hit {
	unsigned int curve;
	unsigned int surface;
	qaws_scalar t, u, v;
	qaws_vec3 position;
} qaws_curve_surface_batch_hit;

/* points[first .. first + count - 1] of a point buffer; closed when the last
   point joins the first; sorted by (surface_a, surface_b) */
typedef struct qaws_surface_batch_curve {
	unsigned int surface_a;   /* surface_a < surface_b */
	unsigned int surface_b;
	unsigned int first, count;
	int closed;
} qaws_surface_batch_curve;

/* rays o + t d, 0 <= t <= max_t */
typedef struct qaws_surface_ray_desc {
	qaws_surface const* const* surfaces;
	unsigned int surface_count;
	qaws_scalar const* origins;            /* 3 scalars per ray */
	qaws_scalar const* directions;         /* 3 scalars per ray; t in units of each */
	unsigned int ray_count;
	qaws_scalar max_t;                     /* 0 = no limit */
	qaws_scalar flatness;                  /* 0 = 2^-9 of the extent */
	qaws_batch_executor const* executor;   /* optional */
} qaws_surface_ray_desc;

typedef struct qaws_surface_ray_hit {
	unsigned int surface;     /* QAWS_CURVE_BATCH_NONE: no hit */
	qaws_scalar t, u, v;
	qaws_vec3 position;
} qaws_surface_ray_hit;
```

The header describes the intended methods: curve-surface pairs seeded by a
chord segment against a patch's two triangles and refined by Newton on
S(u, v) = C(t); surface-surface segments moved onto both surfaces by a
minimum-norm Newton on S_a(u1, v1) = S_b(u2, v2) and chained into
polylines; rays walking the grid cell by cell (Amanatides-Woo) and refined
by Newton on S(u, v) = o + t d. Tangential contacts and grazes are not
reported.
