# API Reference: qaws_surface_batch.h

Closest points over N surfaces at once, and the types of surface batches.
The hit and ray types (`qaws_curve_surface_batch_*`, `qaws_surface_batch_curve`,
`qaws_surface_ray_*`) belong to batched intersections and ray casts, which
are not part of qaws.

---

## Pieces and grid

- **Surfaces** become an adaptive (u, v) quadtree of patches. A patch is
  split until its centre and edge midpoints lie within the bound of the
  bilinear interpolation of its corners. Every patch is split at least
  twice per direction.
- Every piece carries a box inflated by twice its measured deviation.
- One uniform grid over all pieces is searched ring by ring from each query.

## Prepared sets

```c
qaws_status qaws_surface_set_create(qaws_surface_set** out_set, qaws_surface_batch_desc const* desc);
void qaws_surface_set_destroy(qaws_surface_set* set);
unsigned int qaws_surface_set_get_patch_count(qaws_surface_set const* set);

typedef struct qaws_surface_batch_desc {
	qaws_surface const* const* surfaces; unsigned int surface_count;
	unsigned int const* families;      /* used by batched intersections */
	qaws_scalar flatness;
} qaws_surface_batch_desc;
```

A surface set flattens its surfaces once, and every query on it reuses the
patches; `qaws_surface_set_find_closest` below is one. Sets keep pointers
to their surfaces, so the surfaces must outlive them unchanged.

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

qaws_status qaws_surface_batch_find_closest(qaws_surface_batch_closest* out_points,
	qaws_surface_batch_stats* out_stats, qaws_surface_batch_closest_desc const* desc);
qaws_status qaws_surface_set_find_closest(qaws_surface_batch_closest* out_points, qaws_surface_batch_stats* out_stats, qaws_surface_set const* set,
	qaws_scalar const* points, unsigned int point_count, qaws_scalar max_distance);
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
answer is farther than a 160 × 160 sampling of the surfaces. The pairwise
`qaws_surface_find_closest_point` (certified, test 149) is farther on none;
its 8 × 8 grid and Newton alone were farther on 32.
