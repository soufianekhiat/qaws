# API Reference: qaws_operations.h

Curve operations (split, join, offset, reparameterization).

Thread-safe: operates on immutable input curves.

---

## qaws_curve_split

```c
qaws_status qaws_curve_split(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_curve** out_left,
    qaws_curve** out_right);
```

Splits a curve at the given parameter into two sub-curves.

**Parameters:**
- `curve` -- The curve to split.
- `parameter` -- Parameter value at which to split. Must be within the curve's domain.
- `out_left` -- Receives the left sub-curve (from domain start to `parameter`).
- `out_right` -- Receives the right sub-curve (from `parameter` to domain end).

**Returns:** `QAWS_STATUS_OK` on success.

Every kind can be split. Bezier, B-spline, NURBS and cubic trajectories
split into two curves of their own kind, as does a Hermite or Catmull-Rom
curve at a span boundary. Everywhere else the halves come from
`qaws_curve_extract` below, so a cut inside a Hermite span gives two exact
cubic B-splines.

---

## qaws_curve_extract

```c
qaws_status qaws_curve_extract(
    qaws_curve const* curve,
    qaws_scalar t0,
    qaws_scalar t1,
    qaws_curve** out_curve);
```

Returns the piece of `curve` between `t0` and `t1` as a new curve that lies
exactly on the source. With `t0 > t1` the piece is reversed. Parameters
within rounding of the domain are clamped onto it.

| Source kind | Piece |
|---|---|
| Bezier, rational Bezier | same kind, de Casteljau (homogeneous for rational), domain [0, 1] |
| B-spline, NURBS | same kind, knot insertion; the source parameters are kept |
| Polynomial | same kind, domain narrowed to [t0, t1] |
| Arc | arc with new start and end angles |
| Clothoid | clothoid starting at C(t0) with the heading and curvatures there |
| Hermite, Catmull-Rom (any parameterization), trajectory, subdivision | each cubic span rebuilt as a Bezier span: one span gives a cubic Bezier, more give a cubic B-spline whose knots are the source parameters |
| Composite | the piece of one segment, or a composite of pieces |
| Reparameterized | the piece of its source |
| Yuksel | no polynomial form: cubic spans fitted to 1e-9 of the piece's extent (1e-5 in float builds) |

Test 82 extracts 8 parameter pairs from a curve of each kind. Every piece
starts at C(t0), ends at C(t1), and stays within 1e-9 of the extent from the
source, both ways (f64; 1e-5 in float builds).

---

## qaws_curve_join

```c
qaws_status qaws_curve_join(
    qaws_curve const* curve_a,
    qaws_curve const* curve_b,
    qaws_curve** out_joined);
```

Joins two curves end-to-end into a single curve.

**Parameters:**
- `curve_a` -- First curve.
- `curve_b` -- Second curve. Must have the same dimension as `curve_a`.
- `out_joined` -- Receives the joined curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_offset_2d

```c
qaws_status qaws_curve_offset_2d(
    qaws_curve const* curve,
    qaws_scalar distance,
    int trim,
    qaws_curve** out_curves,
    unsigned int curve_capacity,
    unsigned int* out_count);
```

Computes the offset (parallel) curve at a signed distance. 2D only.

Positive distance offsets to the left of the travel direction. At cusps where the curvature radius is less than the absolute distance, backward loops are trimmed at their self-intersection points to produce a clean boundary.

**Parameters:**
- `curve` -- The 2D source curve.
- `distance` -- Signed offset distance. Positive = left, negative = right.
- `trim` -- Trim level controlling cleanup:

| Trim | Description |
|---|---|
| 0 | Raw offset -- cusp loop removal only. |
| 1 | + Distance-based trim: removes points closer than `|distance|` to the source curve. |
| 2 | + Self-intersection cleanup: detects remaining self-intersections, extracts each loop as a separate closed curve, and discards tails. Produces clean, non-self-intersecting closed loops. |

- `out_curves` -- Caller-allocated array to receive output curve pointers.
- `curve_capacity` -- Maximum number of curves that `out_curves` can hold.
- `out_count` -- Receives the number of curves actually written.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_reparameterize_arc_length

```c
qaws_status qaws_curve_reparameterize_arc_length(
    qaws_curve const* curve,
    unsigned int table_resolution,
    qaws_curve** out_curve);
```

Creates an arc-length reparameterized wrapper curve.

The resulting curve has parameter domain [0, total_arc_length] and maps uniformly to arc length along the source curve.

**Parameters:**
- `curve` -- The source curve. Must outlive the returned wrapper (non-owning reference).
- `table_resolution` -- Size of the arc-length lookup table. Pass 0 for the default (256).
- `out_curve` -- Receives the reparameterized wrapper curve.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The source curve must remain valid for the lifetime of the returned wrapper curve.
- The wrapper does not take ownership of the source curve.
