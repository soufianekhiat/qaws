# API Reference: qaws_operations.h

Curve operations that build new curves from existing ones: split, extract, join, offsets (2D trimmed, 2D cubic,
3D), arc-length reparameterization and length matching, chain merging, fillets and chamfers. All are thread-safe:
they only read the input curves.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Ownership: every `qaws_curve**` output is a new curve owned by the caller (destroy it with `qaws_curve_destroy`,
see [curve.md](curve.md)). Outputs are allocated with the source curve's allocator and set to NULL on failure.
Wrappers (`reparameterize_arc_length`, `match_arc_length`) keep a non-owning reference to their source, which
must outlive them.

Related: [inspect.md](inspect.md), [composite.md](composite.md), concepts in [../operations.md](../operations.md).

---

## Split and extract

```c
qaws_status qaws_curve_split(
    qaws_curve** out_left,
    qaws_curve** out_right,
    qaws_curve const* curve,
    qaws_scalar parameter);

qaws_status qaws_curve_extract(
    qaws_curve** out_curve,
    qaws_curve const* curve,
    qaws_scalar t0,
    qaws_scalar t1);
```

`qaws_curve_split` cuts a curve at `parameter` into the piece before (`out_left`) and after (`out_right`).

- `QAWS_STATUS_OUT_OF_RANGE` outside the domain; `QAWS_STATUS_INVALID_ARGUMENT` at either end of the domain (an
  empty piece) or for a null pointer.
- Every kind can be split. Bezier, B-spline, NURBS and trajectories split into two curves of their own kind, as do
  Hermite and Catmull-Rom curves at a span boundary. Everywhere else (including a cut inside a Hermite span) the
  halves come from `qaws_curve_extract`.

`qaws_curve_extract` returns the piece between `t0` and `t1` as a new curve that lies exactly on the source; with
`t0 > t1` the piece is reversed. Parameters within rounding of the domain are clamped onto it.

| Source kind | Piece |
|---|---|
| Bezier, rational Bezier | same kind (de Casteljau, homogeneous for rational) |
| B-spline, NURBS | same kind by knot insertion; the source parameters are kept |
| Polynomial | same kind |
| Arc | arc with new start and end angles |
| Clothoid | clothoid starting at C(t0) with the heading and curvatures there |
| Hermite, Catmull-Rom, trajectory, subdivision | a cubic Bezier (one span) or a cubic B-spline whose knots are the source parameters |
| Composite | the piece of one segment, or a composite of pieces |
| Reparameterized | the piece of its source |
| Yuksel | no polynomial form: cubic B-spline fitted to 1e-9 of the piece's extent (1e-5 in float builds) |

---

## Join

```c
qaws_status qaws_curve_join(
    qaws_curve** out_joined,
    qaws_curve const* curve_a,
    qaws_curve const* curve_b);
```

Joins two curves end to end into one curve of the same kind.

- Both curves must have the same kind and dimension (`QAWS_STATUS_INVALID_ARGUMENT` otherwise).
- Supported kinds: Hermite, B-spline, NURBS (same degree; B's parameters are shifted to start where A's end, and
  B's first control point is taken as the shared point, so the join is C0), Catmull-Rom and trajectory.
- `QAWS_STATUS_UNSUPPORTED_OPERATION` for every other kind (Bezier, rational Bezier, composite, arc, polynomial,
  clothoid, Yuksel, subdivision, wrappers): build a composite instead ([composite.md](composite.md)).

---

## 2D offset with trimming

```c
qaws_status qaws_curve_offset_2d(
    qaws_curve** out_curves,
    unsigned int* out_count,
    qaws_curve const* curve,
    unsigned int curve_capacity,
    qaws_scalar distance,
    int trim);
```

The offset (parallel) curve at signed `distance`, positive to the left of the travel direction. 2D only.

- The offset is sampled at 512 points. Where the curvature radius is below `abs(distance)` the raw offset reverses;
  the backward loops are trimmed at their self-intersections.
- `trim` levels:

| trim | Cleanup |
|---|---|
| 0 | raw offset, cusp loop removal only |
| 1 | + distance trim: removes points closer than `abs(distance)` to the curve |
| 2 | + self-intersection cleanup: each remaining loop becomes a separate closed curve (closed at the crossing), tails are discarded, and crossings between output curves are trimmed |

- Each output piece is a centripetal Catmull-Rom curve through the surviving samples, so the result approximates
  the offset. At most `curve_capacity` curves are written into `out_curves`; extra pieces are destroyed and
  `*out_count` is the number written.
- Returns `QAWS_STATUS_INVALID_ARGUMENT` (null pointer, `distance == 0`), `QAWS_STATUS_INVALID_DIMENSION`,
  `QAWS_STATUS_ALLOCATION_FAILURE`, or `QAWS_STATUS_INTERNAL_ERROR` when no piece survives.

---

## 2D offset as cubic Beziers

```c
qaws_status qaws_curve_offset_cubic_2d(
    qaws_curve** out_curve,
    qaws_curve const* curve,
    qaws_scalar distance,
    qaws_scalar tolerance);
```

The parallel curve at signed `distance` (positive to the left) as a composite of cubic Bezier segments, within
`tolerance` (> 0, in the curve's units) of the true offset. 2D only, any family.

- Each span is first made cubics (a cubic polynomial span is one exactly). Each cubic is split at its inflections
  and curvature maxima, and each piece is tried as a single offset cubic, then as circular arcs, then halved.
- A cusp becomes an arc about it, a corner between spans a line.
- Where the curvature radius is below `abs(distance)` the result keeps its loops: use `qaws_curve_offset_2d` to
  have them trimmed.
- `QAWS_STATUS_DEGENERATE_CURVE` when the curve has no extent.

---

## 3D offset along a direction field

```c
qaws_status qaws_curve_offset_3d(
    qaws_curve** out_curve,
    qaws_curve const* curve,
    qaws_surface const* surface,
    qaws_vec3 const* direction,
    unsigned int sample_count,
    qaws_scalar distance,
    int direction_mode);
```

Moves `sample_count` samples (0: 256) of a 3D curve by `distance` along a unit direction and fits a cubic B-spline
through them (`sample_count / 4` control points, at least 8). The result approximates the offset.

| direction_mode | Direction | Requires |
|---|---|---|
| 0 | the constant `*direction` (normalized) | `direction` non-NULL |
| 1 | the curve's Frenet normal | - |
| 2 | the surface normal at the closest surface point | `surface` non-NULL |

Returns `QAWS_STATUS_INVALID_ARGUMENT` (null pointer or a missing `direction` / `surface`),
`QAWS_STATUS_INVALID_DIMENSION` (2D curve), `QAWS_STATUS_ALLOCATION_FAILURE`, or an evaluation or fitting status.

---

## Arc-length reparameterization

```c
qaws_status qaws_curve_reparameterize_arc_length(
    qaws_curve** out_curve,
    qaws_curve const* curve,
    unsigned int table_resolution);

qaws_status qaws_curve_match_arc_length(
    qaws_curve** out_curve,
    qaws_curve const* curve_a,
    qaws_curve const* curve_b,
    unsigned int table_resolution);
```

- `reparameterize_arc_length`: a wrapper curve (kind `QAWS_CURVE_KIND_REPARAMETERIZED`, one span) with domain
  `[0, total arc length]` whose parameter is the arc length along the source. Derivatives are returned with respect
  to that arc length.
- `match_arc_length`: a wrapper over `curve_b` with domain `[0, arc_length(curve_a)]`; parameter `s` evaluates
  `curve_b` at the same fraction of its length as `s` is of `curve_a`'s, so both curves are traversed over the same
  total length.
- `table_resolution`: arc-length lookup table size (0: 256).
- The wrapper references its source without owning it: `curve` / `curve_b` must outlive it. `curve_a` is only
  measured.
- `QAWS_STATUS_DEGENERATE_CURVE` for a curve of zero length; `QAWS_STATUS_INVALID_ARGUMENT`,
  `QAWS_STATUS_ALLOCATION_FAILURE`.

---

## Merge a chain

```c
qaws_status qaws_curve_merge_chain(
    qaws_curve** out_curves,
    unsigned int* out_count,
    qaws_curve const* const* curves,
    unsigned int curve_count,
    unsigned int curve_capacity,
    unsigned int target_degree,
    qaws_scalar tolerance);
```

Replaces a chain of segments (which should share endpoints, C0) by a single approximating B-spline of degree
`target_degree` (0: 3), fitted through 64 samples per input curve; the parameter runs over the chain, each curve
contributing its domain length.

- `tolerance`: when positive and the first fit (`samples / 4` control points) deviates more than it at the samples,
  one refit with `samples / 2` control points is made. The tolerance is not guaranteed after that refit.
- Always writes one curve: `out_curves[0]`, `*out_count = 1`.
- `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer, `curve_count == 0` or `curve_capacity == 0`.

---

## Fillet and chamfer

```c
qaws_status qaws_curve_fillet_2d(
    qaws_curve** out_curve,
    qaws_curve const* curve,
    qaws_scalar radius);

qaws_status qaws_curve_chamfer_2d(
    qaws_curve** out_curve,
    qaws_curve const* curve,
    qaws_scalar distance);
```

Round (`fillet`: a circular arc of `radius`) or cut (`chamfer`: a straight line set back `distance` from the
corner) the sharp corners of a 2D composite curve; the result is a new composite.

- Corners are the junctions between composite segments where the tangent turns by more than about 5 degrees
  (0.087 radians). The segments next to a corner are trimmed back and refitted as cubic B-splines, and the arc or
  line is inserted between them.
- A curve that is not a composite has no corners: the result is a cubic B-spline fitted to 64 samples of it (an
  approximation, not a copy).
- `radius` / `distance` must be positive (`QAWS_STATUS_INVALID_ARGUMENT`); `QAWS_STATUS_INVALID_DIMENSION` for a 3D
  curve.
