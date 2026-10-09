# API Reference: qaws_surface_batch.h

Batched intersections with surfaces. These calls share the flattening and
the grid of [curve_batch.md](curve_batch.md).

---

## Pieces and grid

- **Curves** become chord segments, split span by span until the midpoint
  and both quarter points lie within the flatness bound of the chord.
- **Surfaces** become an adaptive (u, v) quadtree of patches. A patch is
  split until its centre and edge midpoints lie within the bound of the
  bilinear interpolation of its corners. Every patch is split at least
  twice per direction.
- Every piece carries a box inflated by twice its measured deviation.
- One uniform grid over all pieces visits each overlapping pair once.

## Curves x surfaces

```c
typedef struct qaws_curve_surface_batch_desc {
	qaws_curve const* const* curves;   unsigned int curve_count;   /* 3D */
	qaws_surface const* const* surfaces; unsigned int surface_count;
	qaws_scalar flatness;              /* 0 = 2^-9 of the scene extent */
} qaws_curve_surface_batch_desc;

typedef struct qaws_curve_surface_batch_hit {
	unsigned int curve, surface;
	qaws_scalar t, u, v;
	qaws_vec3 position;
} qaws_curve_surface_batch_hit;      /* sorted by (curve, surface, t) */

qaws_status qaws_curve_surface_batch_find_intersections(
	qaws_curve_surface_batch_desc const* desc,
	qaws_curve_surface_batch_hit* out_hits, unsigned int hit_capacity,
	unsigned int* out_count, qaws_surface_batch_stats* out_stats);
```

Only segment / patch pairs are kept from the grid. Each pair is processed
in three steps:

1. The segment is intersected with the patch's two triangles. The
   barycentric and segment slack is scaled by both inflations.
2. The result seeds Newton on S(u, v) = C(t), which is then polished to
   rounding level.
3. Coincident hits of one curve / surface pair are kept once.

## Surfaces x surfaces

```c
typedef struct qaws_surface_batch_desc {
	qaws_surface const* const* surfaces; unsigned int surface_count;
	unsigned int const* families;      /* optional: same id, never intersected */
	qaws_scalar flatness;
} qaws_surface_batch_desc;

typedef struct qaws_surface_batch_curve {
	unsigned int surface_a, surface_b; /* a < b */
	unsigned int first, count;         /* into the qaws_ssi_point buffer */
	int closed;
} qaws_surface_batch_curve;

qaws_status qaws_surface_batch_find_intersections(
	qaws_surface_batch_desc const* desc,
	qaws_surface_batch_curve* out_curves, unsigned int curve_capacity, unsigned int* out_curve_count,
	qaws_ssi_point* out_points, unsigned int point_capacity, unsigned int* out_point_count,
	qaws_surface_batch_stats* out_stats);
```

Patch pairs of different surfaces (and different families, when given)
are processed in four steps:

1. Their triangles are intersected. A vertex lying exactly on the other
   triangle's plane counts as a crossing.
2. Each small segment's ends are moved onto both true surfaces by a
   minimum-norm Newton on S_a(u1, v1) = S_b(u2, v2):
   dx = Jᵀ (J Jᵀ)⁻¹ F.
3. End points that coincide are joined. Each dangling end is then bridged
   to its nearest dangling partner within four flatness bounds, one to one,
   so runs of short segments never collapse.
4. Ends still dangling away from the surface boundaries mark real gaps,
   where the triangles missed at a shallow crossing angle. From each one,
   the call marches along the curve tangent n_a × n_b in steps of two
   flatness bounds. Each step is corrected onto both surfaces, until
   another dangling end is in reach.
5. The segments are walked into open or closed polylines.

A curve or point that does not fit in the buffers is not written, and the
counts are the totals.

Tangential contacts and branches closer than the flatness bound are not
resolved.

## Certified versions

[exact.md](exact.md) documents `qaws_exact_curve_surface_batch_hits`. It
works on exact curves and surfaces: it puts sound control boxes of all
spans and patches in one grid, then certifies every candidate span / patch
pair as `qaws_exact_curve_surface_hits` does. That call now skips disjoint
span / patch boxes too.

Every desc here, float and exact, takes an optional `executor` that runs the
grid cells (or the query chunks) on the caller's threads; see
[curve_batch.md](curve_batch.md#parallel-execution).

## Numbers

These are f64 timings from test 81 and `qaws_batch_showcase` (Release build).

| Scene | Batch | All pairs |
|---|---|---|
| 300 lines × 3 paraboloids | 383 hits (closed form), 0.007 s | 385 (2 duplicates), 0.44 s |
| 3 paraboloids × 6 planes | 19 curves, every circle closed / every arc, 0.009 s | marching: 6 curves, 0.18 s |
| B-spline terrain × 16 planes (contours) | 31 curves (the certified count), 0.014 s | marching: 5 curves, 0.17 s |
| 150 arcs × terrain | 129 impacts, 0.003 s | 129, 0.073 s |

`qaws_exact_surface_batch_hits` gives certified intersection curves of N
exact surfaces. It grids sound patch boxes, then runs each surface pair's
certified solve on its overlapping patch pairs only. On the B-spline
terrain against 16 planes it certifies 31 contour branches in 0.5 s. The
float batch gives the same branch count on all 16 levels.

## Prepared sets

```c
qaws_status qaws_surface_set_create(qaws_surface_batch_desc const* desc, qaws_surface_set** out_set);
void qaws_surface_set_destroy(qaws_surface_set* set);
unsigned int qaws_surface_set_get_patch_count(qaws_surface_set const* set);

qaws_status qaws_surface_set_find_intersections(qaws_surface_set const* set, qaws_surface_set const* other,
	qaws_surface_batch_curve* out_curves, unsigned int curve_capacity, unsigned int* out_curve_count,
	qaws_ssi_point* out_points, unsigned int point_capacity, unsigned int* out_point_count, qaws_surface_batch_stats* out_stats);

qaws_status qaws_curve_set_find_surface_intersections(qaws_curve_set const* curves, qaws_surface_set const* surfaces,
	qaws_curve_surface_batch_hit* out_hits, unsigned int hit_capacity, unsigned int* out_count, qaws_surface_batch_stats* out_stats);
```

A surface set flattens its surfaces once, for example a terrain queried
against changing planes, curves or other surfaces. A query re-runs only
the grid, the narrow phase and the chaining:

- `other == NULL` intersects the set with itself, using its families.
- Otherwise every surface of `set` is intersected with every surface of
  `other`. `surface_a` indexes `set` and `surface_b` indexes `other`.
- A prepared 3D curve set (`qaws_curve_set`, see
  [curve_batch.md](curve_batch.md)) can be queried against a prepared
  surface set.

Sets keep pointers to their surfaces, so the surfaces must outlive them
unchanged. At equal flatness, the results equal the one-shot calls curve
for curve and point for point (test 81).

## Closest points on surfaces

```c
typedef struct qaws_surface_batch_closest_desc {
	qaws_surface const* const* surfaces; unsigned int surface_count;
	qaws_scalar const* points; unsigned int point_count;   /* 3 scalars per point */
	qaws_scalar max_distance;                              /* 0 = no limit */
	qaws_scalar flatness;
} qaws_surface_batch_closest_desc;

typedef struct qaws_surface_batch_closest {
	unsigned int surface;             /* QAWS_CURVE_BATCH_NONE: none within max_distance */
	qaws_scalar u, v, distance;
	qaws_vec3 position;
} qaws_surface_batch_closest;

qaws_status qaws_surface_batch_find_closest(qaws_surface_batch_closest_desc const* desc,
	qaws_surface_batch_closest* out_points, qaws_surface_batch_stats* out_stats);
qaws_status qaws_surface_set_find_closest(qaws_surface_set const* set, qaws_scalar const* points, unsigned int point_count,
	qaws_scalar max_distance, qaws_surface_batch_closest* out_points, qaws_surface_batch_stats* out_stats);
```

This returns, for every query point, the nearest point over all surfaces.

1. One grid over the flattened patches is searched ring by ring outward
   from the point. Patch corners lie on the surface and bound the search.
2. Patches are refined nearest box first, while a box can beat the best
   point. Each is seeded at the nearest point of its two triangles.
3. Newton runs on (S − p) · S_u = (S − p) · S_v = 0 with the full
   Hessian. Gauss-Newton is used where the Hessian is not positive.
4. A parameter on its bound with the descent pointing out of the domain
   stays there while the other takes the 1D step along that edge; a
   parameter the 2D step would take out stops on its bound the same way.

Test 81 runs 400 points against three paraboloids and two planes. No
answer is farther than a 160 × 160 sampling of the surfaces, while the
pairwise `qaws_surface_find_closest_point` lands farther on 32 of them.

## Ray casting

```c
typedef struct qaws_surface_ray_desc {
	qaws_surface const* const* surfaces; unsigned int surface_count;
	qaws_scalar const* origins;          /* 3 scalars per ray */
	qaws_scalar const* directions;       /* 3 scalars per ray; t in units of each */
	unsigned int ray_count;
	qaws_scalar max_t;                   /* 0 = no limit */
	qaws_scalar flatness;
} qaws_surface_ray_desc;

typedef struct qaws_surface_ray_hit {
	unsigned int surface;                /* QAWS_CURVE_BATCH_NONE: no hit */
	qaws_scalar t, u, v;
	qaws_vec3 position;
} qaws_surface_ray_hit;

qaws_status qaws_surface_batch_raycast(qaws_surface_ray_desc const* desc, qaws_surface_ray_hit* out_hits,
	qaws_surface_batch_stats* out_stats);
qaws_status qaws_surface_set_raycast(qaws_surface_set const* set, qaws_scalar const* origins, qaws_scalar const* directions,
	unsigned int ray_count, qaws_scalar max_t, qaws_surface_ray_hit* out_hits, qaws_surface_batch_stats* out_stats);
```

This returns the first hit of every ray o + t d, with 0 ≤ t ≤ max_t.

1. Each ray walks the grid over the flattened patches cell by cell
   (Amanatides–Woo).
2. A patch met is tested against its two triangles, with slack from its
   inflation, and the result seeds Newton on S(u, v) = o + t d.
3. The walk stops once the best hit lies before the current cell's exit,
   because any patch hit earlier along the ray overlaps a cell already
   visited.

Tangential grazes are not reported.

Test 81 casts 1000 rays against three paraboloids and a plane. Every
first hit matches the closed form, with t error under 4e-13, in 0.007 s.
The pairwise curve / surface call on segment curves would take about
1.8 s, and it misses one ray in 100.
