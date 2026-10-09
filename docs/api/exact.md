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
qaws_status qaws_exact_orient2d(double const a[2], double const b[2], double const c[2], qaws_exact_sign* out_sign, qaws_exact_path* out_path);
qaws_status qaws_exact_orient3d(double const a[3], double const b[3], double const c[3], double const d[3], qaws_exact_sign* out_sign, qaws_exact_path* out_path);
qaws_status qaws_exact_compare_ratio(double a, double b, double c, double d, qaws_exact_sign* out_sign, qaws_exact_path* out_path);
```

`out_path` (may be NULL) tells whether the f64 filter or the exact path
decided. Inputs spanning more than about 900 binary orders of magnitude
return `QAWS_STATUS_EXACT_RANGE_EXCEEDED`.

## Exact curves

```c
void qaws_exact_desc_default(qaws_exact_desc* desc);   /* 2^-20 lattice, 26 coordinate bits, 24 parameter and weight bits */
qaws_status qaws_exact_curve_prepare(qaws_exact_desc const* desc, qaws_curve const* curve, qaws_exact_curve** out_curve, qaws_exact_report* out_report);
void qaws_exact_curve_destroy(qaws_exact_curve* curve);
qaws_status qaws_exact_curve_evaluate(qaws_exact_curve const* curve, double t, unsigned int order, double* out, qaws_exact_report* out_report);
unsigned int qaws_exact_curve_span_count(qaws_exact_curve const* curve);
qaws_status qaws_exact_curve_span_bezier(qaws_exact_curve const* curve, unsigned int s, unsigned int* out_degree, double* out_t0, double* out_t1, double* out_points, double* out_weights);
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
qaws_status qaws_exact_winding_2d(qaws_exact_curve const* const* pieces, unsigned int count, double const p[2], int* out_winding);
```

The ray crossings `x > p_x` are decided on the exact integer Bernstein
coefficients of `Y - p_y W` and `X - p_x W`, with exact halving
subdivision. The loop must close exactly and share one exact space.

## Exact surfaces

```c
qaws_status qaws_exact_surface_prepare(qaws_exact_desc const* desc, qaws_surface const* surface, qaws_exact_surface** out_surface, qaws_exact_report* out_report);
void qaws_exact_surface_destroy(qaws_exact_surface* surface);
qaws_status qaws_exact_surface_evaluate(qaws_exact_surface const* surface, double u, double v, unsigned int order, double* out, double* out_normal);
void qaws_exact_surface_patch_count(qaws_exact_surface const* surface, unsigned int* out_u_count, unsigned int* out_v_count);
qaws_status qaws_exact_surface_patch_bezier(qaws_exact_surface const* surface, unsigned int iu, unsigned int iv, unsigned int* out_p, unsigned int* out_q, double out_rect[4], double* out_points, double* out_weights);
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

## Certified intersections

```c
qaws_status qaws_exact_curve_line_hits(qaws_exact_curve const* curve, double const p0[2], double const p1[2], double width, qaws_exact_hit* out_hits, unsigned int capacity, unsigned int* out_count);
qaws_status qaws_exact_curve_plane_hits(qaws_exact_curve const* curve, double const point[3], double const normal[3], double width, qaws_exact_hit* out_hits, unsigned int capacity, unsigned int* out_count);
```

On each span, `W N . (C - P)` is an integer Bernstein polynomial (the line
or plane is taken exactly from its doubles). Its roots are isolated by
subdivision at midpoints, counting sign variations of the coefficients:
zero variations, no root; one, exactly one simple root. Hits are sorted
and come in three kinds:

- `QAWS_EXACT_HIT_POINT`: an exact intersection parameter (a dyadic root,
  found when a subdivision or refinement point lands on it);
- `QAWS_EXACT_HIT_CROSSING`: `(t_lo, t_hi)` holds exactly one transversal
  intersection, refined by sign bisection down to `width` (0: about 2^-60
  of a span), with the bounds rounded outward;
- `QAWS_EXACT_HIT_OVERLAP`: the curve lies in the line or plane on
  `[t_lo, t_hi]`, merged across spans.

No intersection lies outside the hits, and a hit on a knot is reported
once. A tangency at a non-dyadic parameter is a multiple root that signs
cannot tell from two close crossings, so it returns
`QAWS_STATUS_CERTIFICATION_FAILED` instead of a guess. Two crossings
2^-40 apart are still separated.

```c
qaws_status qaws_exact_curve_curve_hits(qaws_exact_curve const* a, qaws_exact_curve const* b, qaws_exact_pair* out_pairs, unsigned int capacity, unsigned int* out_count);
```

Two 2D curves are compared span against span:

1. The lower-degree span I (degree n from 2 to 6) is implicitized
   exactly. With `Bez(P, Q)` the Bezout matrix of two of its homogeneous
   polynomials,
   `M(x, y, w) = x Bez(Y, W) + y Bez(W, X) + w Bez(X, Y)` is singular
   exactly on its curve.
2. The other span S is substituted: `g(r) = det M(S(r))` is an integer
   polynomial of degree `n m <= 64`, and its roots are isolated as above.
3. At a root, the kernel of the symmetric `M` is `(1, s, s^2, ...)`, and
   I's own parameter is `s = C01 / C00` (cofactors, polynomials in r). A
   root is kept only when `0 <= s <= 1` is proven by the signs of `C00`,
   `C01` and `C00 - C01` on its interval, refining until they settle.
   The bound on s is the range of the coefficient ratios. When s lies
   exactly on the span end at an irrational root (a gradient line through
   a knot of a contour), `C01` (or `C00 - C01`) never settles: then
   `gcd(R, C01)` (or `gcd(R, C00 - C01)`) decides it: one sign variation
   of the gcd on the root's interval means the root is the gcd's, and s is
   exactly 0 (or 1); none means it is not; otherwise the interval is refined.
4. Segments against segments are solved directly. A degenerate span (a
   line traced by a quadratic, for example) is retried with the roles
   swapped.

`QAWS_STATUS_CERTIFICATION_FAILED` is returned for:
- tangencies;
- a common component (overlapping curves);
- crossings through a singular point.

All curves share one exact space: the coordinate lattice decides what
"touching" means, and a 2^-k gap below the lattice step is quantized
away.

In 3D, the system is solved in a coordinate plane where the projected
curves do not share a component (curves in a vertical plane use another
one). Each projected root must also satisfy the remaining coordinate's
equation, `H(s(r), r) = 0`, which is an exact zero test at an algebraic
number:
- the gcd of R and H is computed by Brown's modular algorithm (gcds
  modulo 31-bit primes, Chinese remaindering, then exact trial division
  of both polynomials as the proof);
- a constant gcd proves that no projected crossing is a 3D intersection;
- otherwise the gcd's sign variations on the isolating interval decide.

```c
qaws_status qaws_exact_curve_self_hits(qaws_exact_curve const* curve, qaws_exact_pair* out_pairs, unsigned int capacity, unsigned int* out_count);
```

Self-intersections are pairs a < b with C(a) = C(b), in 2D or 3D:

- **Inside a span:** the divided differences
  `(X(s) W(t) - X(t) W(s)) / (s - t)` are eliminated with a Bezout matrix
  whose entries are polynomials in t. The constraint s > t is the sign of
  `C01 - t C00`. Spans proven monotone along a coordinate (or the sum or
  difference of two) are skipped.
- **Across spans:** the span pairs are compared as above.
- **Not reported:** shared knots and the closing point of a closed curve.
- **Spans on one algebraic curve** (the arcs of a NURBS conic) give a
  zero resultant. They are separated by inverting their endpoints
  exactly, including the conic parameter's point at infinity. Arcs that
  touch at their ends are fine; an overlap is refused.
- **Refused:** an irrational cusp, and a span folding back on itself.

## Certified batched curve intersections

```c
qaws_status qaws_exact_curve_batch_hits(qaws_exact_batch_desc const* desc, qaws_exact_batch_hit* out_hits, unsigned int capacity,
	unsigned int* out_count, qaws_exact_batch_stats* out_stats);
```

N exact curves at once (all 2D or all 3D, one exact space):

1. Every span gets a box from its control points. When all its weights are
   positive the span lies in the hull, and the correctly rounded ratios
   widened by two ulps make the box sound. A span with a non-positive weight
   gets an unbounded box.
2. One uniform grid over all span boxes (the broad phase shared with the
   float batch) gives the span pairs of different curves whose boxes
   overlap.
3. Only those pairs are certified, exactly as `qaws_exact_curve_curve_hits`
   does. That call now skips span pairs with disjoint boxes too.

`families` skips the pairs of curves with the same id, and
`QAWS_EXACT_BATCH_SELF` adds each curve's self-intersections. A curve pair
that fails certification contributes no hits and is counted in
`stats.uncertified_count`. Every other pair is still reported, and the call
returns `QAWS_STATUS_CERTIFICATION_FAILED`.

On a closed curve, a hit at its start is reported once, by its last span.

Example: on the 72-curve Gaussian heightfield of `qaws_batch_showcase`
(8267 spans), all 375 contour / gradient crossings are certified in 0.23 s.
The pairwise exact call over every contour / gradient pair takes 2.1 s, and
the float batch takes 0.014 s.

## Certified line / surface intersections (exact ray casting)

```c
qaws_status qaws_exact_surface_line_hits(qaws_exact_surface const* surface, double const p0[3], double const p1[3], qaws_exact_surface_hit* out_hits, unsigned int capacity, unsigned int* out_count);
```

The line `x = p0 + t (p1 - p0)` is computed per patch:

1. The line is the meet of two planes, with normals `d x e_k` and
   `d x (d x e_k)`, exact in dyadics.
2. On a patch, each plane gives an integer polynomial in (u, v).
3. One parameter is eliminated with a Bezout matrix whose entries are
   polynomials in the other: the lower degree first, the other when the
   resultant vanishes. A degree-1 elimination is solved directly.
4. The roots are isolated, and the eliminated parameter is recovered from
   the cofactors and proven in [0, 1].
5. The line parameter t is enclosed over the exact sub-patch: the patch is
   restricted to the (u, v) box by integer De Casteljau at rational cut
   points, and t is the range of `d . (X - p0 W) / (W |d|^2)` over its
   control net (the weights are positive).

Hits are sorted by t, and a hit on a patch edge is reported once. A
tangent ray or a line lying on the surface returns
`QAWS_STATUS_CERTIFICATION_FAILED`.

## Certified curve / surface intersections

```c
qaws_status qaws_exact_curve_surface_hits(qaws_exact_curve const* curve, qaws_exact_surface const* surface, qaws_exact_curve_surface_hit* out_hits, unsigned int capacity, unsigned int* out_count);
```

Resultant elimination would be far past the integer budget here (degree
around 650 for a cubic curve against a bicubic patch), so this uses
certified subdivision. Per curve span and surface patch, the three
equations `X_c(u, v) W_C(t) - C_c(t) W(u, v) = 0` have exact integer
coefficients in the product Bernstein basis over the (u, v, t) box. The
boxes are subdivided exactly (integer De Casteljau, gcd-normalized), and
each box is tested:

- **Exclusion:** some equation, or some preconditioned combination
  `G_i = sum_c Y_ic F_c`, has one strict sign on the box, so it holds no
  root.
- **Exactly one root:** G passes the Poincare-Miranda sign test on the
  box faces (existence), and its interval Jacobian, bounded by Bernstein
  derivative coefficients, is strictly diagonally dominant (G is
  injective, so the root is unique). Y is an approximate inverse Jacobian
  rounded to integers; any Y keeps the test exact.

Certified boxes are shrunk by bisection to about 2^-44. Each step keeps
the half that is proven to hold the root, either because the other half
is excluded or because Miranda holds. When the root sits on a cut, the
middle half [1/4, 3/4] is tried, and a stuck direction waits until the
others have shrunk. A root on a shared face is reported once.

`QAWS_STATUS_CERTIFICATION_FAILED` is returned when the box budget runs
out: a tangency, the curve lying on the surface, or crossings closer than
the subdivision can separate.

### Curves x surfaces at once

```c
qaws_status qaws_exact_curve_surface_batch_hits(qaws_exact_surface_batch_desc const* desc, qaws_exact_curve_surface_batch_hit* out_hits,
	unsigned int capacity, unsigned int* out_count, qaws_exact_batch_stats* out_stats);
```

The call takes N 3D exact curves and M exact surfaces in one exact space.
Every span and patch gets a sound box from its control points: the
correctly rounded ratios widened by two ulps, or an unbounded box when a
weight is not positive. One grid gives the span / patch pairs whose boxes
overlap, and each is certified as by `qaws_exact_curve_surface_hits`.

A root on a span or patch edge is reported once. A curve / surface pair
that fails certification is counted in `stats.uncertified_count`, and the
call returns `QAWS_STATUS_CERTIFICATION_FAILED` after reporting every
other pair.

## Certified surface / surface intersection curves

```c
qaws_status qaws_exact_surface_surface_hits(qaws_exact_surface const* a, qaws_exact_surface const* b, unsigned int min_depth, qaws_exact_ssi_point* out_points, unsigned int point_capacity, unsigned int* out_point_count, qaws_exact_ssi_branch* out_branches, unsigned int branch_capacity, unsigned int* out_branch_count);
```

Per patch pair, the three equations `X^a_c W^b - X^b_c W^a = 0` in
(u1, v1, u2, v2) have exact integer Bernstein coefficients. The 4D box is
subdivided exactly. The cuts fall at 7/16 of a box rather than its
middle, so symmetric and dyadic positions of a curve do not fall on box
edges.

- **Excluded boxes:** one equation keeps a strict sign.
- **Regular boxes:** a 3 x 3 minor of the interval Jacobian has an exact
  interval determinant that excludes 0. By the implicit function theorem,
  every solution arc in the box is then a monotone graph over the
  remaining variable, so it has two distinct ends on the box boundary and
  cannot close inside.
- **Face points:** the arcs' ends are the solutions of three equations in
  three unknowns on each face, certified by the same solver as curve /
  surface intersections.
- **Topology:** a regular box with no face point holds no arc. With one,
  the curve only touches the box. With two, it holds exactly one arc
  joining them. Anything else is subdivided again.

Branches chain these points across boxes and patch pairs. Points on a
closed surface's seam (a NURBS conic tube, for example) are identified at
both parameter ends. Between consecutive points of a branch lies exactly
one smooth intersection arc, and every intersection curve is covered.
`min_depth` subdivides regular boxes further to give more points.
Tangential contact, overlapping surfaces and singular intersection points
return `QAWS_STATUS_CERTIFICATION_FAILED`.

While building this, two soundness details of the 3-equation solver were
fixed:
- **Invertible preconditioner:** the integer preconditioner Y is now scaled
  row by row and must have a non-zero exact determinant before its
  Miranda test counts. A singular Y could certify a box without a root.
- **Roots on a cut:** a root sitting exactly on a cut is certified in the
  middle [1/4, 3/4]^3 of the box, and that region is then marked as
  covered.

## Certified 2D Booleans

```c
qaws_status qaws_exact_boolean_2d(qaws_exact_curve const* a, qaws_exact_curve const* b, unsigned int operation, qaws_exact_piece* out_pieces, unsigned int piece_capacity, unsigned int* out_piece_count, qaws_exact_loop* out_loops, unsigned int loop_capacity, unsigned int* out_loop_count);
```

Regions are bounded by closed, simple 2D exact curves; the operation is
`QAWS_BOOLEAN_UNION`, `_INTERSECTION` or `_DIFFERENCE` (a minus b). The
result is made of pieces of the inputs, never refitted. Each piece is
`(region, start, end, reversed)`, with its start and end as parameter
enclosures.

1. Both boundaries are checked to be simple (certified self-intersections).
2. They are cut at their certified crossings.
3. Each piece is classified by the certified winding number of the other
   boundary around an exact rational point inside it (integer De
   Casteljau at a rational parameter).
4. Pieces are kept by the operation and linked into loops. Two kept
   pieces meet at each crossing; at a touching point the same curve goes
   on.

Tangencies and overlaps return `QAWS_STATUS_CERTIFICATION_FAILED`.

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
| curve / line (2D), curve / plane (3D) intersections | certified: exact points, one-root intervals, overlaps; non-dyadic tangencies refused |
| curve / curve intersections (2D, 3D) | certified: implicitization of the lower-degree span, root isolation, exact inversion (span ends at irrational roots by a gcd test); 3D: exact zero test of the third coordinate |
| N curves at once (2D, 3D) | certified: one grid over sound span boxes, only overlapping span pairs certified, families, uncertified pairs counted |
| N curves x M surfaces at once | certified: one grid over sound span and patch boxes, only overlapping span / patch pairs certified |
| self-intersections (2D, 3D) | certified: divided differences inside spans, span pairs, knots and closing points excluded, spans on one conic told apart |
| line / surface intersections (ray casting) | certified: u, v and t enclosed, per patch Bezout elimination |
| 2D Boolean operations (union, intersection, difference) | certified: pieces of the source curves, cut at certified crossings, classified by exact winding |
| curve / surface intersections | certified: exact Bernstein subdivision, Miranda existence, diagonal-dominance uniqueness |
| surface / surface intersections | certified topology: branches of certified points, one smooth arc between consecutive points (seams identified) |

Tests: 63 (integers), 64 (predicates), 65 (Bezier), 66 (winding), 67
(B-spline / NURBS), 68 (Hermite, Catmull-Rom, polynomial), 69 (surfaces), 70 (line / plane hits), 71 (curve / curve hits), 72 (3D curve / curve), 73 (self-intersections), 74 (line / surface), 75 (Booleans), 76 (composites, frozen Catmull-Rom), 77 (curve / surface), 78 (surface / surface), 80 (batched curves, against the float batch), all against
Mathematica exact references. Figures:
`examples/exact_showcase.c`.
