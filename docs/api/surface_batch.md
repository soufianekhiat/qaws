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
4. The segments are walked into open or closed polylines.

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

## Numbers

These are f64 timings from test 81 and `qaws_batch_showcase` (Release build).

| Scene | Batch | All pairs |
|---|---|---|
| 300 lines × 3 paraboloids | 383 hits (closed form), 0.007 s | 385 (2 duplicates), 0.44 s |
| 3 paraboloids × 6 planes | 19 curves, every circle closed / every arc, 0.009 s | marching: 6 curves, 0.18 s |
| B-spline terrain × 16 planes (contours) | 32 curves, 0.013 s | marching: 5 curves, 0.16 s |
| 150 arcs × terrain | 129 impacts, 0.003 s | 129, 0.073 s |
