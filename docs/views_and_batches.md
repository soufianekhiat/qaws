# Views, Prepared Evaluation and Batches

This page covers three ways to evaluate geometry that differ from the usual owned `qaws_curve` and `qaws_surface`:

- **Views** read geometry from the host's own arrays without copying them.
- **Prepared objects** turn a view, or raw input, into data that is fast to evaluate.
- **Batches** answer one query over many curves or surfaces at once.

API pages: [view](api/view.md), [prepare](api/prepare.md), [curve_batch](api/curve_batch.md),
[surface_batch](api/surface_batch.md), [batch_executor](api/batch_executor.md).

> Names are given without their precision suffix. Each one exists as `name_f32` and `name_f64`; see
> [precision.md](precision.md).

## Borrowed views

A `qaws_array_view` describes a strided array that the host owns. Its elements are f32 or f64, with 1 to 4
components. The array is either:

- **packed** (`data` set): element i, component c is at `data + i * stride + c * size`. This covers interleaved
  records of any stride.
- **planar** (`data` NULL): component c of element i is at `planes[c] + i * stride`.

A stride of 0 means tightly packed. Values are read byte by byte, so no alignment is required. A view's element type
is independent of the library's `qaws_scalar`: each value is converted as it is read, and no converted copy is made.

A `qaws_curve_view` describes a curve through array views. It covers splines (B-spline, NURBS and Bezier: points,
optional weights, optional knots), Hermite, Catmull-Rom, polynomial, arc, trajectory, Yuksel and subdivision curves.
Each family evaluates as the owned family it names. A `qaws_surface_view` describes a tensor-product B-spline or NURBS
grid in row-major order.

`qaws_array_view_validate`, `qaws_curve_view_validate` and `qaws_surface_view_validate` check counts, layouts, knot
order and weight signs. Every function that takes a view runs these checks. `qaws_array_view_check_extent` also
confirms that every byte a view reads lies inside a given host buffer.

### What you can make from a view

| Object | Owns | Use |
|---|---|---|
| `qaws_curve_prepared` / `qaws_surface_prepared` | derived data only: knots in the computing precision, the span table | fastest evaluation, one parameter or a range of them |
| borrowed curve or surface (`qaws_curve_create_borrowed`, `qaws_surface_create_borrowed`) | a prepared object inside | a normal `qaws_curve` / `qaws_surface` (kind `VIEW`) for sampling, arc length, closest points, frames, surface jets to third order |
| owned copy (`qaws_curve_create_from_view`, `qaws_surface_create_from_view`) | everything | a normal owned B-spline, NURBS or family curve; the view is read during the call only |

A borrowed curve works with every operation that only evaluates. Extract, split and reverse return owned curves,
which they make through an owned copy. Operations that read a family's own arrays refuse a borrowed curve, because the
host holds those arrays. These include control-point and knot inspection, transforms, and derivatives with respect to
control points.

### Lifetime contract

- qaws never writes to, frees or keeps the host memory a view points to, except as stated below.
- A function that takes a view reads it during that call only.
- A prepared object or a borrowed curve or surface reads the view's arrays for its whole life. The arrays must stay
  valid (not freed, not moved) until the object is destroyed, and must not change while an evaluation runs.
- Between evaluations, the host may change anything the object does not depend on.
  `qaws_curve_prepared_dependencies` returns `QAWS_DEPENDS_*` flags: `POINTS`, `WEIGHTS`, `KNOTS` and `LAYOUT`
  (pointers, strides, counts, types, degree). An object never depends on the values of the points or weights, so a
  solver or an animation may rewrite them in place without preparing again. Any change to something the object depends
  on requires destroying it and preparing again.
- The prepared object's own memory comes from the given allocator (or the default chain). Alternatively,
  `qaws_curve_prepare_into` places it in caller memory of `qaws_curve_prepared_size` bytes, aligned to 16, with no
  allocation at all. Destroying such an object does nothing.
- Evaluation never allocates.

### Range evaluation

`qaws_curve_prepared_evaluate_range` evaluates parameters `begin` to `end - 1` of a parameter view into the output
views (`position`, `d1`, `d2`, `d3`, each f32 or f64, each optional). `qaws_surface_prepared_evaluate_range` does the
same for surfaces, with outputs `position`, `du`, `dv`, `normal`, `duu`, `duv` and `dvv`. Disjoint ranges write
disjoint elements, so the host's workers can split one call's range between them. On failure the call returns the
first failing status, and the outputs are written up to that element.

## Backend preparation (qaws_prepare.h)

`qaws_prepare.h` is a separate, C-only set of functions. It turns raw geometry into the flat arrays that the shader and
array backends evaluate (see [backends.md](backends.md)):

- Hermite, Catmull-Rom and trajectory curves become per-span `a, b, c, d` coefficients.
- A Bezier curve gets the control points of its derivatives.
- A rational Bezier curve gets homogeneous weighted control points.
- An arc gets its axis vectors and angle range.

The caller allocates every output array. These functions keep no state.

## Batches

A batch answers one query over N curves or M surfaces at once. Each curve is flattened into chord segments, and each
surface into an adaptive quadtree of (u, v) patches. Every piece's box is inflated by its measured deviation, and one
uniform grid over all pieces finds the candidates. Newton's method then refines each candidate on the true geometry.
The cost therefore grows with the number of pieces and their nearby pairs, not with the number of curve or surface
pairs.

- **Closest points.** `qaws_curve_batch_find_closest` and `qaws_surface_batch_find_closest` return, for every query
  point, the nearest point over all inputs. On a tie, the lower index wins. With `max_distance`, a point that has no
  input within that distance gets `QAWS_CURVE_BATCH_NONE`.
- **Flatness.** This is the bound on the deviation of a piece from the true geometry. A value of 0 means 2^-10 of the
  scene extent for curves and 2^-9 for surfaces. A smaller bound gives more pieces and fewer missed local minima. Within
  one surface patch, a local minimum that the seed does not reach is not reported.
- **Prepared sets.** `qaws_curve_set` and `qaws_surface_set` flatten their inputs once, which is the costly part. Every
  `*_set_find_closest` call then reuses that work. A set keeps pointers to its curves or surfaces: they must outlive it
  and stay unchanged.
- **Shared types.** The intersection, level-crossing, ray and surface-curve hit types in these headers
  (`qaws_curve_batch_hit_*`, `qaws_level_crossing*`, `qaws_curve_surface_batch_*`, `qaws_surface_batch_curve`,
  `qaws_surface_ray_*`, `families`, `QAWS_CURVE_BATCH_SELF`) describe batched intersections and ray casts. qaws does
  not provide those operations. The types are declared here so that libraries built on qaws can share them.

### Parallel execution

Every batch desc has an optional `qaws_batch_executor`, and prepared sets copy it from their desc. The batch splits its
independent work (the query points) into at most 256 chunks, and the executor runs those chunks on the host's threads.
Each chunk has its own scratch memory and statistics, and these are merged in chunk order, so the results do not
depend on the executor. See [threading.md](threading.md).
