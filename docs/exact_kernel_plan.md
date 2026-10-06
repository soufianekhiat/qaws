# Exact rational geometry kernel: development plan

Working plan for branch `exact`, derived from
`qaws_exact_rational_kernel_design.md` (the specification) with the review
changes agreed on 2026-10-06:

1. B-spline/NURBS extraction width bound made concrete before coding it, and
   the supported degree stated from it (not "0..16" for extracted splines).
2. Predicates also run on raw IEEE doubles (exact dyadic decode, adaptive
   fallback), next to the lattice path: callers get robust signs without a
   quantization policy.
3. A thin vertical slice first (integers, exact Bezier, filtered predicates,
   one real consumer), extraction and surfaces after.
4. Test numbers follow the current tree (63 and up); span selection follows
   the current runtime, including the repeated-end-knot rule of `9045cdb`.

Rules that hold throughout (from the specification):

- `qaws_scalar`, `QAWS_BACKEND` and every GPU core header stay untouched by
  the exact path; the exact code is C-only, under `src/qaws/exact/`.
- One rounding boundary: input quantization (and the final conversion back
  to f32/f64 when asked). Never multiply-and-round inside an exact kernel.
- Overflow is impossible by proof, or reported (`QAWS_STATUS_EXACT_RANGE_EXCEEDED`);
  never wrapped, never silently approximated.
- Epsilons may trigger the exact fallback; they never decide a sign.
- Ground truth: Mathematica (`wolframscript`) exact rational vectors under
  `tests/reference/`, plus in-process identities for fuzzing.

## Slice A: arithmetic foundation

Files: `src/qaws/exact/qaws_exact_int.h/.c`, `src/qaws/qaws_exact.h`.

- `qaws_exact_int`: sign-magnitude, 32-bit limbs (portable C11 products in
  `uint64_t`, no compiler intrinsics), fixed capacity `QAWS_EXACT_MAX_BITS`
  (2048) with a used-limb count so small values stay cheap.
- Operations: set from int64 / from double (exact when integral) / from
  decimal text, add, sub, mul (checked), mul by small, exact division by a
  small divisor, compare, sign, shifts, bit length, to double (round to
  nearest, documented error), to decimal text.
- Every operation returns a status; capacity overflow is
  `QAWS_STATUS_EXACT_RANGE_EXCEEDED`.
- Status codes appended to `qaws_status`: `EXACT_UNSUPPORTED`,
  `EXACT_RANGE_EXCEEDED`, `EXACT_INCOMPATIBLE_SPACE`, `CERTIFICATION_FAILED`.
- Dyadic decode: a finite double is `m * 2^e` exactly (53-bit `m`).

Tests (63_exact_int): Mathematica vectors for add/sub/mul/compare/shift on
values up to ~1900 bits, carries across every limb, signs, zero; randomized
identities `(a+b)-b = a`, `(a*b) / b`, distributivity; overflow reported at
the capacity edge.

## Slice B: exact predicates on doubles

Files: `src/qaws/exact/qaws_exact_predicates.c`.

- `qaws_exact_orient2d(a, b, c)`, `qaws_exact_compare_rational`,
  `qaws_exact_point_equal`: f64 evaluation with a proven error bound first
  (Shewchuk-style static/dynamic filter); when the bound does not prove the
  sign, exact evaluation on the dyadic-decoded inputs (all coordinates
  scaled to the smallest exponent; bounded by the width budget, otherwise
  `RANGE_EXCEEDED`).
- The filter never returns zero; zero is always decided exactly.
- Statistics hook (filter hits / fallbacks) for the tests.

Tests (64_exact_predicates): Mathematica signs on adversarial near-collinear
sets (points along a line perturbed by single ulps); millions of random and
near-degenerate cases checked against the exact path (no wrong sign); the
naive f64 determinant shown to fail on the same set.

Showcase: the classic "orientation of a point near a line" grid (naive f64
sign map with its wrong-sign speckle next to the exact map).

## Slice C: exact Bezier and rational Bezier evaluation

Files: `src/qaws/exact/qaws_exact_bezier.c`, `qaws_exact_quantize.c`,
`qaws_exact_curve.c`.

- Exact space descriptor: power-of-two lattice step, coordinate bits
  (default 26), rounding (nearest-even default); parameter lattice per curve
  (default 24 bits); weights normalized and quantized (24 bits).
- `qaws_exact_curve_prepare` for Bezier and rational Bezier (2D/3D): integer
  homogeneous control points, report with quantization errors.
- Division-free scaled De Casteljau at `s = a/b`; endpoints returned exactly.
- Derivatives D1..D3 from integer derivative polygons, projected with the
  quotient rule, scaled to the public parameter.
- `qaws_exact_curve_evaluate_2d/3d` (round once to `qaws_scalar`) and an
  exact text accessor (numerator / denominator per component) for tests.
- Width planner: `H + d (T + 1)` for evaluation, `H + k (1 + log2 d)` for
  derivative polygons, asserted in debug builds.

Tests (65_exact_bezier): Mathematica exact rationals for positions and D1..D3
of Bezier and rational Bezier, degrees 1..16, 2D and 3D, at lattice
parameters including endpoints; scale invariance; agreement with the f64
evaluator to rounding.

Showcase: exact vs f32 vs f64 evaluation error along a high-degree rational
curve (log scale), and points the f32 evaluator places on the wrong side of
a line that the exact kernel classifies correctly.

## Slice D: a real consumer

Wire the exact predicates into one existing decision that today uses a
tolerance (candidates: the 2D Boolean crossing classification, span
selection at knot values, curve/line side tests in inspection). Report the
numeric quality next to the existing result; the sampled Boolean output
stays labelled approximate.

## Slice E: Hermite, polynomial, uniform Catmull-Rom, composites

Exact conversion to Bernstein form (Hermite: common factor 3; polynomial:
monomial to Bernstein with one common denominator; uniform Catmull-Rom:
rational cubic coefficients). Chordal/centripetal: prepared-quantized mode
with `QAWS_EXACT_FLAG_PREP_QUANTIZED`.

## Slice F: B-spline and NURBS extraction

- Width bound first: local extraction by knot insertion divides by knot
  differences at every step; with a per-span common denominator the growth
  is up to about `p^2 T` bits (degree p, knot bits T). Derive the exact
  bound, choose LCM-based clearing, and state the supported degree at the
  default lattices (expected: p <= 7 within 2048 bits).
- Local extraction per non-empty span, homogeneous for NURBS, integer knot
  spans, span selection identical to the runtime.

Tests: Mathematica exact NURBS values on random and repeated knot vectors;
range-exceeded status beyond the bound.

## Slice G: tensor-product surfaces

Bezier / B-spline / NURBS patches, exact partials and unnormalized normal;
unit normal reported approximate.

## Slice H: certified intersections and topology

Bernstein subdivision, interval root isolation, exact predicates in the
Boolean pipeline keeping source span pieces instead of refitting.

## Bookkeeping

- Build option `QAWS_ENABLE_EXACT` (default ON), independent of
  `QAWS_SCALAR_IS_FLOAT`; f32 and f64 builds run every exact test.
- Docs: `docs/api/exact.md` with the capability matrix (what is and is not
  exact), kept in sync slice by slice.
- Each slice ends with: tests green in f32 and f64, a showcase figure, a
  commit.

## Status (2026-10-06)

| slice | state | commits |
|---|---|---|
| A integers | done: add/sub/mul/shift/div/gcd, nearest double, 2048 bits | cda00e7, 95bd5c9 |
| B predicates on doubles | done: orient2d, orient3d, ratio compare | 5c3402e |
| C exact Bezier | done: rational Bezier, degrees 1..16, C..C''' | a78d309 (fix 9d451c8) |
| D consumer | done: certified winding numbers | c3b8199 |
| F B-spline / NURBS | done: exact blossom extraction; budget measured | 95bd5c9, c3226b1 |
| E other families | done: Hermite, uniform Catmull-Rom, polynomial (composites later) | |
| G surfaces | planned | |
| H intersections | planned | |

Extraction width turned out linear in the degree (about 85 bits per
degree at full precision, gcd/lcm reduction), not p^2 T: every degree up to
16 extracts; the derivative order is what the 2048-bit budget limits (see
docs/api/exact.md).
