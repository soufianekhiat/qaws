# Exact / certified geometry (`qaws_exact.h`)

C-only, built with `QAWS_ENABLE_EXACT` (default ON), independent of
`QAWS_SCALAR_IS_FLOAT` and of the GPU backends. Plan and status:
`docs/exact_kernel_plan.md`.

## Contract

- **Predicates on doubles** are certified: a finite double is an exact dyadic
  rational, so `qaws_exact_orient2d`, `qaws_exact_orient3d` and
  `qaws_exact_compare_ratio` return the true sign of the polynomial in the
  given doubles. An f64 evaluation with Shewchuk's error bound decides when
  it can; otherwise the inputs are decoded exactly and evaluated in
  multi-limb integers. The f64 filter never decides a zero.
- **Exact curves** are exact relative to their lattice inputs: coordinates
  are rounded once onto `x = i 2^space_exp2`, weights scaled by a common
  power of two, knots and parameters onto `t = T 2^-shift`. After that
  boundary nothing rounds. Every value returned as a double is the
  correctly rounded value of an exact rational.
- **Overflow is never silent**: past the 2048-bit integer budget an
  operation returns `QAWS_STATUS_EXACT_RANGE_EXCEEDED`.
- **Undecidable is reported, not guessed**: `QAWS_STATUS_CERTIFICATION_FAILED`
  (for example a winding number queried exactly on the curve).

## Predicates

```c
qaws_status qaws_exact_orient2d(qaws_exact_sign* out_sign, qaws_exact_path* out_path, double const a[2], double const b[2], double const c[2]);
qaws_status qaws_exact_orient3d(qaws_exact_sign* out_sign, qaws_exact_path* out_path, double const a[3], double const b[3], double const c[3], double const d[3]);
qaws_status qaws_exact_compare_ratio(qaws_exact_sign* out_sign, qaws_exact_path* out_path, double a, double b, double c, double d);
```

`out_path` (may be NULL) tells whether the f64 filter or the exact path
decided. Inputs spanning more than about 900 binary orders of magnitude
return `QAWS_STATUS_EXACT_RANGE_EXCEEDED`.

## Exact curves

```c
void qaws_exact_desc_default(qaws_exact_desc* desc);   /* 2^-20 lattice, 26 coordinate bits, 24 parameter and weight bits */
qaws_status qaws_exact_curve_prepare(qaws_exact_curve** out_curve, qaws_exact_report* out_report, qaws_curve const* curve, qaws_exact_desc const* desc);
void qaws_exact_curve_destroy(qaws_exact_curve* curve);
qaws_status qaws_exact_curve_evaluate(double* out, qaws_exact_report* out_report, qaws_exact_curve const* curve, double t, unsigned int order);
unsigned int qaws_exact_curve_span_count(qaws_exact_curve const* curve);
qaws_status qaws_exact_curve_span_bezier(unsigned int* out_degree, double* out_t0, double* out_t1, double* out_points, double* out_weights, qaws_exact_curve const* curve, unsigned int s);
```

A prepared curve holds one integer homogeneous Bezier per span. Bezier and
rational Bezier curves are one span. B-spline and NURBS spans come from
exact blossoms (de Boor triangles with the span ends as arguments), reduced
by gcd and cleared to one common denominator. Evaluation is division-free
De Casteljau on the k-th integer derivative polygon. The exact quotient rule
gives `C^(k)` as one rational per component. Span selection: `a <= T < b`,
the last span at the domain end.

Cubic Hermite and uniform Catmull-Rom spans become integer cubic Beziers
with one common factor (3 and 6). Polynomial coefficients are decoded
exactly as dyadics, re-expanded around the quantized domain start and
converted to Bernstein form over one integer denominator; only the domain
ends are quantized. Chordal and centripetal Catmull-Rom need square roots
and return `QAWS_STATUS_EXACT_UNSUPPORTED`.

The report gives the quality (`QAWS_NUMERIC_EXACT_RATIONAL`), whether
inputs were quantized, the widest stored integer and the quantization
errors.

Integer budget at full precision (24-bit knots and weights, 26-bit
coordinates), measured by test 67:

| degree | span storage (bits) | highest derivative within 2048 bits |
|---|---|---|
| 1..5 | 49..341 | C''' |
| 6..7 | 427..521 | C'' |
| 8..9 | 618..711 | C' |
| 10..16 | 801..1343 | position |

Coarser lattices raise every row.

## Certified winding numbers

```c
qaws_status qaws_exact_winding_2d(int* out_winding, qaws_exact_curve const* const* pieces, double const p[2], unsigned int count);
```

The ray crossings `x > p_x` are decided on the exact integer Bernstein
coefficients of `Y - p_y W` and `X - p_x W`, with exact halving
subdivision. The loop must close exactly and share one exact space.

## Exact surfaces

```c
qaws_status qaws_exact_surface_prepare(qaws_exact_surface** out_surface, qaws_exact_report* out_report, qaws_surface const* surface, qaws_exact_desc const* desc);
void qaws_exact_surface_destroy(qaws_exact_surface* surface);
qaws_status qaws_exact_surface_evaluate(double* out, double* out_normal, qaws_exact_surface const* surface, double u, double v, unsigned int order);
void qaws_exact_surface_patch_count(unsigned int* out_u_count, unsigned int* out_v_count, qaws_exact_surface const* surface);
qaws_status qaws_exact_surface_patch_bezier(unsigned int* out_p, unsigned int* out_q, double out_rect[4], double* out_points, double* out_weights, qaws_exact_surface const* surface, unsigned int iu, unsigned int iv);
```

Bezier, B-spline and NURBS surfaces become one integer homogeneous Bezier
patch per non-empty knot rectangle: exact blossoms along u on each control
column, then along v on each row, one lcm for the net. `order` 0, 1, 2
returns S; then Su, Sv; then Suu, Suv, Svv. They come from integer
derivative nets (rows, then the column) and the bivariate quotient rule,
for example `Suv = (Xuv W^2 - Wuv X W - Wu Nv - Wv Nu) / W^3` with
`Nu = Xu W - X Wu`. `out_normal` is `Su x Sv` with exact components, so
its sign along any lattice direction is certified. The unit normal needs a
square root and is not exact.

## Certified closest points

The exact batch descs (`qaws_exact_batch_desc`, `qaws_exact_ssi_batch_desc`)
give N exact curves or surfaces in one exact space, with an optional
`executor`.

### Over N curves

```c
qaws_status qaws_exact_curve_batch_closest(qaws_exact_closest_point* out_points, qaws_exact_batch_stats* out_stats, double const* points,
	unsigned int point_count, qaws_exact_batch_desc const* desc);
```

For every query point (exact doubles, world units), the call returns the
nearest point over the exact curves, with its curve, a parameter enclosure,
a distance enclosure and a `certified` flag.

- The query point becomes P / Q on the lattice, with Q a power of two.
- On a span C = X / W, the squared distance is stationary at the roots of
  the integer polynomial g(s) = Σ (Q X_c − P_c W)(X_c' W − X_c W'), of
  degree 3n − 2.
- Its roots in (0, 1), isolated exactly, and the span ends are the
  candidates. A knot counts once, and the start of a closed curve once.
- Each candidate's squared distance E / V, with E = Σ (Q X_c − P_c W)²
  and V = (Q W)², is enclosed by the ratios of their Bernstein coefficients
  restricted together to its interval.
- Spans come from one grid over their sound control boxes, nearest box
  first, while a box can still win.
- The winner is certified when its enclosure lies below every other
  candidate's. Overlapping root intervals are refined first.

A span whose roots do not isolate (a double root, or g ≡ 0) enters as one
candidate covering the whole span. An exact tie, such as a point midway
between two lines, is reported with `certified = 0`.

300 dyadic points against dyadic parabolas and segments are all certified,
with distance enclosures at most 1.1e-16 wide, in 0.04 s. Every answer of
the float batch lies inside its enclosure.

### Over M surfaces

```c
qaws_status qaws_exact_surface_batch_closest(qaws_exact_surface_closest_point* out_points, qaws_exact_batch_stats* out_stats, double const* points,
	unsigned int point_count, qaws_exact_ssi_batch_desc const* desc);
```

For every query point, the call returns the nearest point over the exact
surfaces (up to bicubic patches): its surface, (u, v) enclosures, a
distance enclosure and a `certified` flag. On a patch S = X / W the
squared distance to p = P / Q is stationary where G_u = G_v = 0, with
G_u = Σ (Q X_c − P_c W)(X_c,u W − X_c W_u) and G_v likewise.

- **Localising the interior.** Branch and bound on E / V over quarter
  boxes keeps only the boxes that can beat the best distance. The corners
  of each box are exact surface points, so their distances tighten that
  bound as the search goes.
- **Solving.** On each surviving box, G_u and G_v are restricted exactly
  (integer de Casteljau to any dyadic interval). The system with
  3w − 1 = 0 goes to the certified 3D solver, refined only to 2⁻²² in u
  and v.
- **Edges and corners.** These are exact curves, solved by the curve
  method. Each edge and corner counts once over the surface.
- **Certification.** Every candidate gets a distance enclosure, and the
  nearest is certified when its enclosure lies below every other
  candidate's. A pruned box is proven farther than an actual surface
  point, so the nearest point is always among the candidates.

120 dyadic points against three paraboloids and two planes.
all are certified, every answer of the float batch lies inside its
enclosure, and the run takes 1.0 s. A point midway between two planes is
reported as a tie.

### Parallel execution

The exact batch descs (`qaws_exact_batch_desc`, `qaws_exact_ssi_batch_desc`)
take the optional `executor` of the float batches
([curve_batch.md](curve_batch.md#parallel-execution)). The closest-point
batches run query chunks of 4. The kernel has no mutable globals, and the
results are bit-identical whatever the executor does.

## Capability matrix

| quantity | status |
|---|---|
| orient2d / orient3d / ratio compare on doubles | certified |
| Bezier, rational Bezier: position, C', C'', C''' | exact rational |
| B-spline, NURBS: position, C', C'', C''' | exact rational (budget above) |
| span selection | exact (integer knots) |
| winding number of a closed loop | certified |
| Hermite (cubic), uniform Catmull-Rom (open, closed) | exact rational (lattice points and tangents) |
| polynomial | exact rational (dyadic coefficients taken as given, no quantization) |
| composites | exact rational (segment spans moved to [i, i + 1]; non-dyadic bounds reported as a parameter quantization) |
| chordal / centripetal Catmull-Rom | exact relative to the runtime's frozen preparation (QAWS_EXACT_FLAG_PREP_QUANTIZED) |
| Yuksel curves | not rational (trigonometric blend of sub-curves) |
| arcs, clothoids (sin/cos, Fresnel) | not rational |
| unit normals, curvature, arc length (sqrt) | not rational |
| Bezier, B-spline, NURBS surfaces: S, Su, Sv, Suu, Suv, Svv, Su x Sv | exact rational |
| unit surface normal, curvatures (sqrt) | not rational |
| other surface families (sweeps, lofts, offsets, ...) | planned or not rational |
| closest points over N curves | certified: roots of the distance derivative isolated exactly, distance enclosures, unique nearest proven or a tie reported |
| closest points over M surfaces | certified: interior stationary points by the 3D solver on boxes localised by branch and bound, edges as exact curves, unique nearest proven or a tie reported |

Exact intersections, ray casts and 2D Booleans are not part of qaws; they
are built on these exact curves and surfaces.

Tests: 63 (integers), 64 (predicates), 65 (Bezier), 66 (winding), 67
(B-spline / NURBS), 68 (Hermite, Catmull-Rom, polynomial), 69 (surfaces),
76 (composites, frozen Catmull-Rom), all against Mathematica exact
references. Figures: `examples/exact_showcase.c`.
