# Exact kernel

The exact kernel answers geometric questions whose answer must never be
wrong: the sign of an orientation, the winding number of a loop around a
point, the nearest point on a set of curves or surfaces. It sits next to the
float API, prepares its own objects from float curves and surfaces, and
either returns a proven answer or says it cannot. API reference:
[api/exact.md](api/exact.md).

It is C only, runs on the CPU, and is compiled when `QAWS_ENABLE_EXACT` is
on (default ON).

## What it models

### Doubles are exact rationals

A finite double is the dyadic rational `m 2^e`. The sign of a polynomial in
double inputs is therefore a well defined exact quantity, and the predicates
(`qaws_exact_orient2d`, `qaws_exact_orient3d`, `qaws_exact_compare_ratio`)
return that sign for the inputs exactly as given: nothing is rounded on the
way in.

### Curves and surfaces on lattices

Curves and surfaces go through one rounding boundary, the preparation
(`qaws_exact_curve_prepare`, `qaws_exact_surface_prepare`), driven by a
`qaws_exact_desc`:

- coordinates onto the lattice `x = i 2^space_exp2` (round to nearest, ties
  to even), `|i| < 2^coord_bits`;
- weights scaled by a common power of two to `weight_bits` (a common factor
  cancels in a rational curve, so integer weights stay exact);
- knots and parameters onto `t = T 2^-shift`, `shift` derived from the
  domain and `param_bits`.

The result is integer homogeneous Bezier data: one integer `(X, W)` control
polygon per curve span, one integer control net per surface patch. B-spline
and NURBS spans come from exact blossoms cleared to one common denominator;
Hermite and uniform Catmull-Rom spans are integer cubic Beziers with a common
factor; polynomial coefficients are taken exactly as dyadics. After that
boundary nothing rounds. The report (`qaws_exact_report`) states what the
boundary did: `QAWS_EXACT_FLAG_INPUT_QUANTIZED` and the largest position,
weight and parameter roundings.

So an exact curve is exact relative to its lattice inputs, not to the float
object it came from. The difference is bounded and reported. Chordal and
centripetal Catmull-Rom need square roots for their knot spacing; the
runtime's own float preparation is frozen and taken exactly, flagged
`QAWS_EXACT_FLAG_PREP_QUANTIZED`.

The rational domain covers Bezier, rational Bezier, B-spline, NURBS, cubic
Hermite, Catmull-Rom, polynomial and composite curves (2D and 3D), and
Bezier, B-spline and NURBS surfaces. Families whose geometry is not
rational (arcs and clothoids through sin/cos and Fresnel integrals, Yuksel
curves, sweeps, offsets, ...) return `QAWS_STATUS_EXACT_UNSUPPORTED`.
Quantities needing a square root (unit normals, curvature, arc length) are
not offered; the surface normal is the exact unnormalized `Su x Sv`.

### Exact spaces

`space_exp2` names the exact space. Operations combining several objects
(a winding loop, a closest-point batch) require one space and return
`QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE` otherwise: integers on different
lattices are not comparable without a rescale the kernel does not make
silently.

## Arithmetic and budgets

Values are fixed-capacity signed integers: sign and magnitude in 32-bit
limbs, at most 2048 bits. Limb products fit in `uint64_t`, so the
arithmetic is portable C11 with no compiler intrinsics, and a value lives on
the stack with no allocation per operation. Every operation that can grow a
value is checked; a result that does not fit returns
`QAWS_STATUS_EXACT_RANGE_EXCEEDED`. Nothing wraps, and overflow is never
silent.

The budget limits what can be computed exactly:

- predicates: inputs spanning more than about 900 binary orders of
  magnitude;
- curves: degree up to 16, and derivatives whose exact quotient-rule
  numerators fit. High degrees reach the budget for `C'''` or `C''` before
  they reach it for the position; coarser lattices (fewer `coord_bits`,
  `param_bits`, `weight_bits`) leave more room;
- surfaces: degree up to 16, up to 64 control points per direction.

Past the budget the call returns `QAWS_STATUS_EXACT_RANGE_EXCEEDED`, never
an approximate value.

## How results are computed

- **Predicates** evaluate in f64 first with Shewchuk's proven error bound.
  When the bound proves the sign, that path decides
  (`QAWS_EXACT_PATH_FILTER`). Otherwise the inputs are decoded into
  integers scaled to the smallest exponent among them (a positive factor
  that keeps the sign) and evaluated exactly (`QAWS_EXACT_PATH_EXACT`). The
  filter never decides a zero. Points sharing a coordinate give a zero
  column and are decided at once.
- **Evaluation** is division-free De Casteljau on multi-limb integers;
  derivatives come from integer derivative polygons and nets with the exact
  (bivariate) quotient rule. Each component is one rational, converted to
  double once, correctly rounded.
- **Winding numbers** count the signed crossings of the ray `x > p_x` on the
  exact integer Bernstein coefficients of `Y - p_y W` and `X - p_x W`,
  subdividing at 1/2 exactly until each piece is decided.
- **Closest points** isolate the roots of the integer polynomial that makes
  the squared distance stationary (on curves), or solve the two-equation
  system with a certified 3D solver (on surfaces), and enclose every
  candidate's distance. A broad phase (a grid over sound control-point
  boxes) limits the work to spans or patches that can still win.

## Guarantees

- A sign returned by a predicate is the true sign of the exact polynomial in
  the given doubles.
- A value returned by an evaluation is the correctly rounded value of the
  exact rational of the prepared object at the lattice parameter.
- An enclosure (`*_lo`, `*_hi`) contains the true value and is rounded
  outward; `lo == hi` means the value is exact.
- A winding number returned with `QAWS_STATUS_OK` is the exact winding
  number of the prepared loop.
- A closest point with `certified = 1` is proven to be the unique nearest
  candidate; with `certified = 0` it is still the nearest candidate found,
  but a tie or a too-close competitor could not be ruled out.

The float API gives no such guarantees: it rounds at every step. Every
answer of the float closest-point batches is expected to fall inside the
exact enclosure; the exact kernel is the reference to check float results
against, and the place to decide what float arithmetic cannot.

## Statuses

| status | meaning |
|---|---|
| `QAWS_STATUS_EXACT_UNSUPPORTED` | the operation or family is outside the rational domain (non-rational geometry, a solver past its degree limit) |
| `QAWS_STATUS_EXACT_RANGE_EXCEEDED` | an integer would outgrow the 2048-bit budget, an input does not fit the lattice, or the inputs span too many exponents |
| `QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE` | objects combined from different exact spaces (`space_exp2`) |
| `QAWS_STATUS_CERTIFICATION_FAILED` | the answer is undefined or cannot be proven within the subdivision budget; for a winding number, the point lies on the curve (or within 200 halvings of it) |

`QAWS_STATUS_INVALID_ARGUMENT` covers non-finite inputs, NULL required
pointers, out-of-range descriptor fields and loops that do not close
exactly. Undecidable answers are reported, never guessed. See
[error_handling.md](error_handling.md) and [api/status.md](api/status.md).

## Threading and memory

The kernel has no mutable globals. Prepared curves and surfaces are
read-only after preparation and may be shared between threads. The
closest-point batches take the optional `qaws_batch_executor` of the float
batches and run query points in chunks of 4; the results are bit-identical
whatever the executor does (see [threading.md](threading.md)). Prepared
objects are allocated by the library and freed with their `*_destroy`
function (see [memory.md](memory.md)).

## Precision

The exact kernel works on `double` in both precision builds. Its names
still carry the precision suffix (`_f32`, `_f64`) because preparation reads
the float objects of that precision (see [precision.md](precision.md)).

## Design reasons

- **Dyadic inputs, integer storage.** Doubles and lattice values are
  dyadic, so homogeneous integer coordinates represent them without loss
  and a common power of two is free to remove. Rationals with arbitrary
  denominators would need gcd work at every step.
- **One rounding boundary.** Quantizing once at preparation and reporting
  the error keeps the exact part closed: every later operation is exact
  relative to the prepared object, and the user knows how far that object is
  from the float one.
- **Division-free evaluation.** Homogeneous De Casteljau and the quotient
  rule keep numerators and one denominator, so a value is rounded exactly
  once, at the end.
- **Fixed capacity, checked growth.** A hard 2048-bit bound keeps memory
  predictable and allocation-free in inner loops; a checked overflow turns
  "too big" into a status rather than a wrong answer.
- **Filter first.** Most predicate calls are far from degenerate; the f64
  filter with a proven bound decides them at float cost, and only the
  ambiguous ones pay for exact integers.
- **Enclosures instead of values.** Intersection and stationary-point
  parameters are algebraic, not rational. The kernel returns certified
  intervals around them rather than pretending to an exact double.
- **Separate from the float API.** The exact objects are prepared from the
  float ones and never change them: the float API keeps its speed and GPU
  portability, and the exact kernel is used where proof matters.
