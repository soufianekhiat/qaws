# QAWS Exact Rational Geometry Kernel
## Integer-Lattice Inputs, Fixed-Width Arithmetic, Certified Predicates, and Exact Rational Span Evaluation

**Status:** Implementation design / developer handoff  
**Repository:** `https://github.com/soufianekhiat/qaws`  
**Validated baseline:** branch `differentiable`, commit `923415d81dfae064ca64742ab8d4bf2a676eee33` (2026-10-05)  
**Scope:** QAWS only. This document is self-contained and does not depend on any external design discussion.

---

## 1. Executive decision

QAWS should gain an optional CPU exact-geometry subsystem, but **it must not be implemented by replacing `qaws_scalar` with a fixed-point type**.

The current QAWS architecture uses `qaws_scalar` as a field-like numerical type across C, HLSL, GLSL, Halide, SIMD, runtime geometry, inspection, and differentiation. A fixed-point `qaws_scalar` would still round after multiplication/division, would not make B-splines or NURBS exact, would break transcendental operations, and would couple exact geometry to shader backends that cannot efficiently support the required wide integers.

Instead, implement exact geometry as an independent C-only arithmetic layer with this model:

1. Input positions, parameters/knots, and weights are quantized once onto explicit integer lattices.
2. After that quantization boundary, exact-capable algorithms perform **no rounding**.
3. Intermediate values widen into fixed-capacity multi-limb signed integers.
4. Exact-capable curve and surface spans are canonicalized into **integer homogeneous Bezier spans**.
5. Runtime exact evaluation uses division-free, wide-integer De Casteljau.
6. Cartesian values are interpreted as exact rational values only at the boundary (`X/W`, `Y/W`, `Z/W`).
7. Predicates run through f64 certified filters first and fall back to exact integer evaluation when the sign is uncertain.
8. Operations involving `sqrt`, trigonometric functions, Fresnel integrals, generic iterative roots, normalization, or fitting remain approximate unless a separate certification algorithm is explicitly implemented.
9. Existing f32/f64, HLSL, GLSL, Halide, SIMD, and public differentiability paths remain intact.

The resulting architecture is not a second copy of QAWS. It is an exact/certified side path used for robust evaluation, predicates, and eventually topology-sensitive operations.

The most important invariant is:

> **Fixed lattice at input, widened integer/rational arithmetic after input. Never fixed-point multiply-and-round inside an exact kernel.**

---

## 2. Why this design fits the current QAWS repository

This proposal was checked against the current `differentiable` branch rather than `main`. `main` is the older v1.0.0 snapshot; the active branch contains the newer multi-backend core and differentiation architecture.

The relevant current code properties are:

### 2.1 `qaws_scalar` is currently an execution/backend scalar

`src/qaws/qaws_platform.h` selects:

- C + f32/f64 through `QAWS_SCALAR_IS_FLOAT`
- HLSL + f32
- GLSL + f32
- Halide + f32/f64 symbolic `Expr`

It also maps `sqrt`, `sin`, `cos`, `atan2`, `pow`, `fmod`, comparisons, literals, and select operations to the active backend.

Therefore exact arithmetic is **not** another value for `QAWS_BACKEND`, and it is **not** another value of `QAWS_SCALAR_IS_FLOAT`.

Execution backend and exact arithmetic are different axes:

```text
Execution target                    Robust/exact numeric path
----------------                    -------------------------
C f32              +                ordinary floating
C f64              +                ordinary floating
C f32/f64          +                exact/certified sidecar
HLSL f32           +                ordinary floating
GLSL f32           +                ordinary floating
Halide f32/f64     +                ordinary floating
```

The exact sidecar is initially C/CPU only.

### 2.2 QAWS already has a useful `prepare + core` split

The branch contains backend-neutral core evaluators under `src/qaws/core/`:

- `qaws_decasteljau_core.h`
- `qaws_rational_bezier_core.h`
- `qaws_bspline_basis_core.h`
- `qaws_bspline_eval_core.h`
- `qaws_nurbs_eval_core.h`
- `qaws_cubic_poly_core.h`
- `qaws_horner_core.h`
- `qaws_arc_core.h`
- `qaws_surface_core.h`
- `qaws_dual_core.h`

and C-only preparation in `qaws_prepare.h/.c`.

The exact subsystem should follow the same architectural idea but use a different canonical prepared representation: exact homogeneous Bezier spans.

### 2.3 Current core arithmetic assumes ordinary scalar division

For example, De Casteljau directly evaluates expressions of the form:

```c
u * p0 + t * p1
```

and the B-spline basis code directly evaluates:

```c
(t - knot_i) / (knot_j - knot_i)
```

Replacing `qaws_scalar` with a conventional Q-format fixed point would round those operations and therefore would not be exact.

### 2.4 Rational Bezier is already conceptually homogeneous

`qaws_rational_bezier_core.h` evaluates weighted homogeneous control points and projects using `1 / w` at the end. The exact path should preserve the homogeneous value instead of projecting early.

This is the key bridge to a unified exact representation.

### 2.5 QAWS runtime objects currently store floating fields

`src/qaws/internal/qaws_internal_types.h` stores control points, knots, weights, coefficients, span boundaries, etc. as `qaws_scalar` arrays. `qaws_curve` is opaque to public callers and dispatches through a family vtable.

Do not replace this storage. Create a separate immutable `qaws_exact_curve` / `qaws_exact_surface` prepared object.

### 2.6 The current differentiation API also uses `qaws_scalar`

`qaws_diff_types.h` represents primal fields, tangents, adjoints, local support weights, reports, condition numbers, residuals, etc. with `qaws_scalar`.

The exact subsystem must therefore be orthogonal to the current differentiation ABI. Exact internal derivatives may be added without changing the existing public JVP/VJP API.

### 2.7 Several current high-level algorithms are explicitly approximate

Examples in the present branch include:

- `qaws_curve_compute_bounds_*`: sampled bounds for the families without Bezier spans.
- `qaws_boolean_2d.c`: sampled arcs and B-spline refitting.
- generic intersection/inspection paths: tolerance-based candidate generation and Newton-style refinement.
- normal computation: cross product followed by `sqrt` normalization.
- arcs: `sin` / `cos`.

Adding an exact scalar type underneath these functions would **not** automatically make them exact. Exactness must be attached to algorithms whose mathematical operations are closed under the exact representation.

---

## 3. Exactness contract

The exact subsystem needs precise terminology. Avoid marketing a result as simply "exact" without saying relative to what.

### 3.1 Exact after quantization

The primary production contract is:

> QAWS exact results are mathematically exact relative to the integer-lattice inputs produced by the selected quantization descriptor.

If an input position `x` is mapped to lattice coordinate `i`, then the exact kernel works with `i`, not with the pre-quantization floating value.

No later arithmetic may silently requantize or round.

### 3.2 Four numeric qualities

Use a small common vocabulary in reports:

```c
typedef enum qaws_numeric_quality
{
    QAWS_NUMERIC_APPROXIMATE = 0,
    QAWS_NUMERIC_CERTIFIED,
    QAWS_NUMERIC_EXACT_RATIONAL,
    QAWS_NUMERIC_UNSUPPORTED
} qaws_numeric_quality;
```

Meaning:

- `APPROXIMATE`: ordinary floating computation; no proof bound.
- `CERTIFIED`: returned classification/sign/bound is mathematically guaranteed although an approximate filter may have been used internally.
- `EXACT_RATIONAL`: all reported exact quantities belong to the rational field of the quantized inputs.
- `UNSUPPORTED`: the requested exact operation leaves the supported domain or exceeds the configured arithmetic range.

### 3.3 Exactness is independent of differentiability validity

Do not overload `qaws_diff_validity`.

These are independent questions:

- Is the numerical value exact/certified/approximate?
- Is the derivative smooth, one-sided, locally valid under a frozen span/active set, ill-conditioned, etc.?

For example, a B-spline point can be `EXACT_RATIONAL` while its derivative at a repeated knot is `QAWS_DIFF_AT_BOUNDARY`.

---

## 4. Goals and non-goals

### 4.1 Goals

The first production implementation should provide:

1. Deterministic integer quantization of QAWS geometry.
2. Portable fixed-capacity multi-limb signed arithmetic in C11.
3. Exact rational evaluation for polynomial/rational spline families.
4. Exact first/second/third spatial derivatives where the family is rational.
5. Certified comparisons and geometric predicates with f64 filters and exact fallback.
6. Exact span selection for quantized parameter values.
7. Exact-capable curve preparation without modifying current f32/f64 storage.
8. Stable behavior across f32 and f64 builds once geometry has entered the exact lattice.
9. Explicit reporting when an operation leaves the rational domain.
10. A multiprecision reference path used for verification, not production hot loops.

### 4.2 Non-goals for V1

Do not attempt these in the first implementation:

- Replace the ordinary QAWS evaluator.
- Run wide exact arithmetic in HLSL/GLSL/Halide.
- Make trigonometric arcs rational.
- Make clothoids/Fresnel integrals rational.
- Make normalized vectors rational in general.
- Build a general algebraic-number library.
- Claim exact generic curve/curve or surface/surface roots merely because the input curves are rational.
- Claim the current sampled/refitted 2D Boolean implementation is exact.
- Use arbitrary precision in production hot paths.
- Flatten all QAWS values into one giant generic rational type.

---

## 5. High-level architecture

Add a C-only exact subsystem beside the existing runtime:

```text
                         Existing QAWS
                  +------------------------+
                  | qaws_curve/qaws_surface|
                  | f32/f64 runtime        |
                  +-----------+------------+
                              |
                              | prepare/quantize once
                              v
                  +------------------------+
                  | qaws_exact_curve       |
                  | qaws_exact_surface     |
                  | immutable prepared data|
                  +-----------+------------+
                              |
             +----------------+----------------+
             |                                 |
             v                                 v
 +-------------------------+       +-------------------------+
 | exact span evaluation   |       | certified predicates    |
 | integer homogeneous     |       | f64 filter -> exact     |
 | Bezier                  |       | fallback                |
 +-------------------------+       +-------------------------+
```

The existing `qaws_curve` remains the authoring/runtime object. The exact object is a prepared representation for robust CPU work.

### 5.1 Why a separate object instead of an internal pointer in `qaws_curve`

Prefer a separate object initially:

- Exact preparation has a quantization policy; one floating curve may be prepared in several exact spaces.
- It avoids increasing every ordinary curve allocation.
- It preserves `qaws_curve_inline`'s fixed 512-byte storage assumptions.
- It keeps GPU/backend code isolated.
- It allows exact objects to be cached by callers according to their own lifetime policy.
- It avoids hidden lazy allocation inside otherwise immutable/thread-safe curve objects.

Both objects should remain immutable and safe for concurrent read-only evaluation.

---

## 6. Input lattice model

### 6.1 Spatial lattice

All exact curves participating in the same cross-curve geometric operation should share one spatial exact context.

Represent a spatial coordinate as:

```text
x_exact = i * 2^space_exp2
```

where `i` is a signed integer lattice coordinate.

The physical power-of-two scale is outside the core integer predicates. Translation and scaling back to world units happen only at API boundaries.

Recommended production storage:

```c
typedef int32_t qaws_exact_coord;
```

with a configurable usable bit count. A reasonable initial default is 26 signed magnitude bits within `int32_t`; this leaves headroom for coordinate differences while still providing about 67 million lattice positions across the signed domain. The correctness of the exact kernel must **not** depend on the default being 26.

Proposed descriptor:

```c
typedef struct qaws_exact_space_desc
{
    int space_exp2;              /* lattice step = 2^space_exp2 world units */
    unsigned int coord_bits;     /* usable signed coordinate bits, default 26 */
    int64_t origin_lattice[3];   /* optional local frame origin in lattice units */
} qaws_exact_space_desc;
```

Internally, geometry may be stored relative to `origin_lattice`. All curves used together must have compatible exact spaces or be explicitly converted to a common one.

### 6.2 Deterministic position quantization

Quantization must specify rounding exactly. Use round-to-nearest, ties-to-even unless the caller requests another documented policy.

```c
typedef enum qaws_exact_rounding
{
    QAWS_EXACT_ROUND_NEAREST_EVEN = 0,
    QAWS_EXACT_ROUND_FLOOR,
    QAWS_EXACT_ROUND_CEIL,
    QAWS_EXACT_ROUND_TOWARD_ZERO
} qaws_exact_rounding;
```

Every exact prepared object stores the quantization descriptor required to reconstruct its physical coordinates.

### 6.3 Parameter and knot lattice

For B-splines/NURBS, `t` and knots should share an integer parameter lattice per curve:

```text
t_exact = T * 2^param_exp2
U_i     = K_i * 2^param_exp2
```

Then basis ratios become integer ratios:

```text
(t - U_i) / (U_j - U_i)
    ==
(T - K_i) / (K_j - K_i)
```

The common power-of-two factor cancels exactly.

Recommended initial storage:

```c
typedef int32_t qaws_exact_param;
```

with `param_bits` defaulting to 24 and an automatically selected `param_exp2` based on the curve's parameter range. The exact width planner, not the nominal storage type, determines whether a prepared kernel is safe.

### 6.4 Weights

NURBS weights should also be quantized once to positive integers with a shared scale:

```text
w_i_exact = W_i * 2^weight_exp2
```

A common scale applied to all weights cancels from the NURBS quotient, so the exact evaluator can operate directly on integer `W_i`.

Recommended initial `weight_bits`: 24. Normalize the floating weights before quantization when safe, for example relative to the maximum weight, because common scaling does not change the rational curve.

Reject non-positive weights where the existing QAWS family requires positivity.

### 6.5 Exact interpretation of existing f32/f64 values

A finite IEEE-754 f32/f64 value is itself an exact dyadic rational. For debugging/reference conversion, QAWS should provide a lossless decoder from `qaws_scalar` bits to `(integer significand, power-of-two exponent)`.

Production exact preparation should normally use the explicit lattice policy above because it gives bounded ranges and predictable widths.

### 6.6 Direct integer input

Preparing from an existing `qaws_curve` is sufficient for V1 functionality, but a direct integer-input API should be added before exact QAWS is advertised as preserving externally supplied integer geometry.

Reason: in an f32 build, integers above `2^24` cannot be represented exactly if they first pass through `qaws_scalar`.

A future low-level API may accept integer control points/knots/weights directly and build the same exact prepared representation without any floating round trip.

---

## 7. Arithmetic model

### 7.1 Do not implement ordinary Q-format arithmetic

This is incorrect for the exact path:

```text
Q26 * Q26 -> round back to Q26
```

because every multiplication destroys exactness.

Instead:

```text
int32 lattice input
      |
      v
wide integer products/sums
      |
      +--> exact rational interpretation when division is mathematically required
```

### 7.2 Fixed-capacity multi-limb signed integers

Add an internal two's-complement or sign+magnitude integer implementation built from 64-bit limbs.

Suggested maximum production capacity:

```text
128, 256, 512, 1024, and optionally 2048 bits
```

The implementation may use one fixed maximum-capacity structure with an active-width field, or generated width-specific structs/functions. Correctness requirements are the same:

- no heap allocation in primitive arithmetic;
- deterministic two's-complement/sign behavior;
- checked addition/subtraction/multiplication;
- exact comparison;
- exact sign;
- left/right power-of-two shifts;
- exact small-integer multiply/divide where divisibility is proven;
- conversion to f64 with a documented error bound;
- debug assertions that discarded upper limbs are only sign extension when a proven narrower result is used.

Example internal shape:

```c
#define QAWS_EXACT_MAX_LIMBS 32 /* 2048 bits */

typedef struct qaws_exact_int
{
    uint64_t limb[QAWS_EXACT_MAX_LIMBS];
} qaws_exact_int;
```

A more specialized generated implementation is acceptable if benchmarks justify it.

### 7.3 Rational values

Use rational values only where they are truly needed:

```c
typedef struct qaws_exact_rat
{
    qaws_exact_int num;
    qaws_exact_int den; /* always positive, non-zero */
} qaws_exact_rat;
```

Rules:

- canonicalize the sign into `num`;
- denominator is always positive;
- reduce by GCD only when useful for controlling growth, not automatically after every operation;
- cross-cancel before multiplication/division where possible;
- compare `a/b` and `c/d` by the exact sign of `a*d - c*b`;
- use width proofs before every cross product.

The exact evaluation hot path should avoid general rational objects whenever possible. That is why homogeneous Bezier spans are the canonical representation.

### 7.4 Homogeneous exact points

Use integer homogeneous coordinates:

```text
2D: (X, Y, W)       -> (X/W, Y/W)
3D: (X, Y, Z, W)    -> (X/W, Y/W, Z/W)
```

A common non-zero scale applied to all components is irrelevant.

This allows polynomial and rational curves to use the same exact evaluation kernel.

---

## 8. Canonical prepared representation: integer homogeneous Bezier spans

This is the central design decision.

Every curve family that can be expressed rationally on a span should be prepared into one or more Bezier spans whose homogeneous control points are integers.

Proposed internal representation:

```c
typedef struct qaws_exact_bezier_span
{
    unsigned int degree;
    qaws_dimension dimension;

    /* Global parameter interval in this curve's integer parameter lattice. */
    qaws_exact_param parameter_begin;
    qaws_exact_param parameter_end;

    /* (degree + 1) * (dimension + 1) signed wide integers. */
    qaws_exact_int* homogeneous_cp;

    /* Width selected/proven for this span. */
    unsigned int storage_bits;
} qaws_exact_bezier_span;

struct qaws_exact_curve
{
    qaws_curve_kind source_kind;
    qaws_dimension dimension;
    unsigned int span_count;
    qaws_exact_bezier_span* spans;
    qaws_exact_space_desc space;
    int param_exp2;
    qaws_numeric_quality quality;
    unsigned int flags;
};
```

`homogeneous_cp` may be stored inline/contiguously in one allocation for the full exact curve.

### 8.1 Why one representation is better than exact versions of every evaluator

With canonical exact Bezier spans:

- Bezier is direct.
- Rational Bezier is direct homogeneous data.
- Hermite becomes cubic Bezier exactly.
- uniform Catmull-Rom becomes cubic Bezier exactly.
- polynomial segments can be converted to Bernstein form exactly.
- B-splines become Bezier spans by exact local extraction.
- NURBS become rational Bezier spans by exact homogeneous local extraction.
- tensor-product spline surfaces use the same idea in two dimensions.

This avoids implementing an exact Cox-de Boor runtime, exact NURBS quotient recurrence, exact Horner path, etc. The current floating paths remain optimized for their respective backends; the exact path has one canonical arithmetic kernel.

---

## 9. Exact curve preparation by family

### 9.1 Bezier

**Quality:** `EXACT_RATIONAL` after position quantization.

For an ordinary polynomial Bezier with integer lattice control points:

```text
H_i = (x_i, y_i, 1)
```

or

```text
H_i = (x_i, y_i, z_i, 1)
```

If the source control points are rational after some conversion, clear one common denominator across the entire span and give every control point the same homogeneous `W`.

### 9.2 Rational Bezier

**Quality:** `EXACT_RATIONAL` after position/weight quantization.

Construct integer homogeneous points from integer position and integer weight:

```text
2D: H_i = (W_i * X_i, W_i * Y_i, W_i)
3D: H_i = (W_i * X_i, W_i * Y_i, W_i * Z_i, W_i)
```

If additional rational coefficients appear during preparation, clear one common denominator across the span. Common scaling of all homogeneous control points does not change the curve.

### 9.3 Hermite

**Quality:** `EXACT_RATIONAL` if points and tangent vectors are quantized rationally.

For one unit span:

```text
B0 = P0
B1 = P0 + M0/3
B2 = P1 - M1/3
B3 = P1
```

Multiply all homogeneous control points by a common factor of 3 to avoid division.

The current QAWS Hermite implementation uses unit-length public spans (`[i, i+1]`), so this conversion matches the present parameterization directly.

### 9.4 Polynomial curve

**Quality:** `EXACT_RATIONAL` when coefficients are quantized rationally.

Convert the monomial basis to Bernstein basis with exact integer/rational combinatorial coefficients. Clear one common denominator for the resulting span.

Do not evaluate exact polynomial curves through the floating Horner core.

### 9.5 Uniform Catmull-Rom

**Quality:** `EXACT_RATIONAL` after point quantization.

Uniform Catmull-Rom has rational cubic coefficients and can be converted directly to homogeneous cubic Bezier spans.

### 9.6 Chordal and centripetal Catmull-Rom

The parameterization depends on Euclidean distances and therefore introduces `sqrt` in general.

Two modes are possible:

- `STRICT_RATIONAL`: unsupported.
- `EXACT_AFTER_PREP_QUANTIZATION`: run the existing floating parameterization preparation, quantize the resulting span parameters/control coefficients once, then treat the prepared curve exactly thereafter.

The report must distinguish this from a curve derived exactly from integer inputs.

### 9.7 B-spline

**Quality:** `EXACT_RATIONAL` after control-point/knot quantization.

Do not implement exact runtime Cox-de Boor as the primary path.

Instead, for each non-empty knot span:

1. Read the local `degree + 1` control points.
2. Compute the local Bezier extraction matrix exactly from integer knot differences.
3. Apply the rational extraction matrix to the local homogeneous control points.
4. Clear one common denominator across all extracted homogeneous control points/components in that span.
5. Store the resulting integer homogeneous Bezier span.

Prefer local Bezier extraction over repeatedly inserting knots into a global copy of the entire curve. Local extraction gives tighter, span-local width bounds and avoids unnecessary growth across unrelated spans.

Span selection uses the quantized integer knot vector, not floating comparisons.

### 9.8 NURBS

**Quality:** `EXACT_RATIONAL` after position/knot/weight quantization.

1. Convert each input point to integer homogeneous form `(w*x, w*y, [w*z], w)`.
2. Apply the same exact local Bezier extraction used by B-splines, but to homogeneous points.
3. Clear a common denominator across the complete span.
4. Store integer homogeneous Bezier control points.

At runtime NURBS then uses exactly the same evaluator as Bezier.

### 9.9 Trajectory

Trajectory exactness depends on how coefficients are produced.

- Explicit positions/times/velocities/accelerations quantized to rational lattices can be rationally prepared.
- Any automatically generated quantity that uses only `+ - * /` remains exact.
- If preparation adds a non-rational normalization or iterative fit, the curve becomes `EXACT_AFTER_PREP_QUANTIZATION` or unsupported in strict mode.

Audit `qaws_trajectory.c`/`qaws_prepare.c` during implementation and classify each preparation branch explicitly.

### 9.10 Composite

A composite is exact if:

- every referenced child span has an exact representation;
- the composite parameter mapping is affine/rational;
- no resampling/refitting occurs.

The exact prepared object may simply concatenate child exact spans with exact global parameter intervals.

### 9.11 Arc

**Strict rational exactness:** unsupported in general.

Current arc evaluation calls `sin` and `cos`. Generic values are transcendental and cannot be represented by the rational kernel.

Do not special-case a few rational/quadrantal angles in V1; that creates surprising family semantics.

### 9.12 Clothoid

**Strict rational exactness:** unsupported.

Fresnel/trigonometric evaluation leaves the rational field.

### 9.13 Yuksel curves

The current differentiable branch includes root solving and trigonometric blending in the Yuksel evaluator/preparation path. Treat strict exact rational preparation as unsupported in V1.

A prepared-exact mode may quantize an already generated Bezier approximation/subcurve, but its report must say that exactness begins only after that preparation boundary.

### 9.14 Reparameterized/derived curves

Classify by the parameter map:

- affine/rational map: exact-capable;
- map using arc length, normalization, roots, trig, etc.: not strict rational.

Never infer exactness solely from the child curve type.

---

## 10. Division-free exact De Casteljau

Let the local Bezier parameter be the exact rational:

```text
s = a / b
```

with `0 <= a <= b` and `b > 0`.

Ordinary De Casteljau is:

```text
Q_i^r = (1-s) Q_i^(r-1) + s Q_(i+1)^(r-1)
```

Multiply the whole level by the common denominator `b`:

```text
Q_i^r_scaled = (b-a) Q_i^(r-1)_scaled + a Q_(i+1)^(r-1)_scaled
```

At level `r`, every control point has the same omitted scale `b^r`. Homogeneous projection is invariant to a common scale, so the evaluator never has to divide by `b^r`.

This gives an exact evaluator containing only:

- wide integer subtraction (`b-a`);
- wide integer multiply-by-parameter-integer;
- wide integer addition;
- exact sign/zero tests.

No general rational object is required in the main evaluation loop.

### 10.1 Local parameter from global parameter

For a prepared span with integer parameter bounds `u0` and `u1`, and integer global parameter `t`:

```text
a = t - u0
b = u1 - u0
s = a / b
```

This is exact. Clamp behavior and boundary span selection must match the current public QAWS semantics.

### 10.2 Endpoint behavior

At `a = 0`, return the first homogeneous control point exactly.

At `a = b`, return the last homogeneous control point exactly.

Do not pass endpoints through a generic slow path.

---

## 11. Exact spatial derivatives

### 11.1 Polynomial homogeneous derivative control polygons

For degree `n` homogeneous Bezier control points `H_i`, the first derivative control polygon is:

```text
D_i = n * (H_(i+1) - H_i)
```

Second and third derivatives are formed recursively in the same way.

All operations are integer and exact.

Evaluate each derivative polygon with the same scaled De Casteljau kernel.

### 11.2 Cartesian derivative of a rational curve

Let the evaluated homogeneous curve be:

```text
H = (X, W)
```

where `X` represents one Cartesian numerator component. Let `X'`, `W'`, `X''`, `W''` be exact homogeneous derivatives.

Then:

```text
C  = X / W
C' = (X' W - X W') / W^2
```

and:

```text
C'' = [X'' W^2 - X W'' W - 2 X' W W' + 2 X (W')^2] / W^3
```

Use equivalent algebra chosen to minimize intermediate width, but keep the derivation documented and tested symbolically.

Third derivative may be added similarly if required by current `QAWS_EVAL_FLAG_D3` support. Derive it once, test it against the multiprecision reference, and do not hand-expand multiple inconsistent versions.

### 11.3 Derivative parameter convention

The exact evaluator must preserve QAWS's existing public parameter semantics.

If a canonical Bezier span uses normalized `s in [0,1]` but the public curve parameter is `t in [u0,u1]`, then:

```text
ds/dt = 1 / (u1-u0)
```

so:

```text
dC/dt   = dC/ds   / (u1-u0)
d2C/dt2 = d2C/ds2 / (u1-u0)^2
d3C/dt3 = d3C/ds3 / (u1-u0)^3
```

Apply these scale factors exactly in the rational result.

This is especially important for B-spline/NURBS spans because the current runtime maps normalized span-local parameters back to the global knot parameter before basis evaluation.

---

## 12. Bit-width planning and overflow policy

Exact arithmetic is only exact if overflow is impossible.

### 12.1 Widths are part of the algorithm contract

Every exact kernel must have a conservative mathematical bit bound.

Do not rely on "probably enough" widths.

Preparation computes the required width before allocating/storing a span. Runtime exact evaluation then dispatches to a width class that is already proven sufficient.

### 12.2 Bezier evaluation bound

Assume:

- initial homogeneous component magnitude needs `H` bits;
- local parameter numerator/denominator need at most `T` bits;
- degree is `d`.

One scaled De Casteljau level computes a sum of two `T`-bit by `B`-bit products. A conservative recurrence is:

```text
B_0 = H
B_r <= B_(r-1) + T + 1
```

Therefore:

```text
B_eval <= H + d * (T + 1)
```

This is a useful first-order planner bound.

Example only, not a hard-coded invariant:

```text
H = 50 bits
T = 24 bits
d = 16

B_eval <= 50 + 16*25 = 450 bits
```

A 512-bit evaluation kernel is therefore plausible for a degree-16 rational span under those input bounds. Preparation/extraction and later predicates may require more.

### 12.3 Derivative polygon bound

Each derivative stage adds:

- one bit for subtraction;
- approximately `ceil(log2(degree))` bits for multiplication by degree.

A conservative bound for derivative order `k` is:

```text
H_k <= H + k * (1 + ceil(log2(d)))
```

followed by evaluation of degree `d-k`:

```text
B_eval_k <= H_k + (d-k) * (T+1)
```

Use a tighter generated bound in production if worthwhile.

### 12.4 Rational comparison doubles width

Comparing:

```text
a/b < c/d
```

requires the sign of:

```text
a*d - c*b
```

so a predicate may require roughly twice the bit width of either rational operand.

This means a curve evaluation that fits in 512 bits may need a 1024-bit comparison kernel.

The width planner must distinguish:

- storage width;
- evaluation width;
- predicate/cross-product width;
- derivative width.

### 12.5 Bezier extraction bound

Exact B-spline/NURBS extraction coefficients are products/ratios of integer knot differences. The bound depends on:

- degree;
- knot integer bit width;
- local knot multiplicities;
- number of extraction steps;
- control-point/weight widths.

Derive and document the local extraction bound independently of the runtime De Casteljau bound.

Prefer local span extraction precisely because it limits this growth.

### 12.6 Surfaces

For tensor-product Bezier evaluation, a first conservative bound is additive in the two parameter directions:

```text
B_surface <= H + du*(Tu+1) + dv*(Tv+1)
```

B-spline/NURBS surface extraction adds its own preparation growth before that evaluation bound.

High degree in both axes can exceed 1024 bits. The implementation must be allowed to select a 2048-bit kernel or return `QAWS_STATUS_EXACT_RANGE_EXCEEDED` rather than silently round.

### 12.7 Production failure policy

Add status codes by appending to `qaws_status`:

```c
QAWS_STATUS_EXACT_UNSUPPORTED,
QAWS_STATUS_EXACT_RANGE_EXCEEDED,
QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE,
QAWS_STATUS_CERTIFICATION_FAILED
```

Never map an overflow/range failure to `QAWS_STATUS_NUMERICAL_FAILURE`; callers need to know that the exact algorithm was valid but the configured integer budget was insufficient.

Never silently fall back from exact to approximate unless the caller explicitly selected such a policy.

---

## 13. Exact/certified predicates

Exact predicates are often more valuable than exact display coordinates because they stabilize topological decisions.

### 13.1 Initial predicate set

Implement exact/certified forms of:

- wide integer sign/zero;
- rational compare;
- homogeneous coordinate compare (`X1/W1` vs `X2/W2`);
- exact point equality in the same exact space;
- 2D orientation for rational points;
- 3D orientation/volume determinant for rational points where needed;
- interval containment in parameter space;
- exact knot/span boundary comparison;
- side-of-line / side-of-plane when line/plane coefficients are exact rational/integer;
- dot/cross sign predicates used by inspection/topology code.

### 13.2 f64 certified filters

Do not run multi-limb exact arithmetic for every easy case.

Each hot predicate should have this shape:

```text
exact integer/rational inputs
          |
          v
f64 approximation + conservative error bound
          |
     +----+----+
     |         |
 sign proven?  uncertain / near zero
     |         |
     v         v
 return     exact wide-int fallback
```

Important rules:

- A filter may return positive or negative only when the error bound proves the sign.
- A filter must never declare exact zero.
- Degenerate/near-degenerate cases always fall back.
- Filter constants/error derivations belong in source comments/tests, not in undocumented magic epsilons.

An interval-style error bound is preferable to one universal epsilon.

### 13.3 Certified bounds

The current sampled bounds APIs are approximate. Do not silently change their semantics.

Add new certified-bound APIs later, based on the convex hull property of Bezier/homogeneous spans and recursive subdivision where rational projection complicates bounds.

---

## 14. Intersections and roots: exact input does not imply rational roots

This section is critical.

Two rational curves can intersect at parameter values that are algebraic and not rational. A rational surface intersection is even more general.

Therefore a fixed-rational backend alone cannot promise an exact Cartesian intersection point for all QAWS intersections.

### 14.1 V1 policy

For intersection/topology code, distinguish:

1. exact/certified predicates on rational control geometry;
2. certified parameter intervals containing a root;
3. approximate root coordinates within a certified interval;
4. exact rational roots when the root happens to be rational and is detected as such.

### 14.2 Recommended future rational-Bezier intersection strategy

A robust path can use:

- exact/certified homogeneous Bezier bounds;
- recursive subdivision;
- Bernstein sign variation / interval exclusion;
- exact orientation predicates;
- interval root isolation;
- Newton only as an accelerator inside a certified bracket, never as the proof of existence/uniqueness.

Return a certification report describing whether the root is isolated and whether the returned coordinate is exact or merely enclosed.

### 14.3 No tolerance-based topology decisions in certified mode

A predicate like:

```text
abs(f) < 1e-5
```

must not decide equality/intersection topology in certified mode.

It may decide whether to invoke the exact fallback.

---

## 15. 2D Boolean implications

The current `qaws_boolean_2d.c` is not suitable as an exact output algorithm because it samples kept arcs and fits a new B-spline.

Even if intersection classification becomes exact, sampling and fitting change the curve.

Therefore:

- exact predicates may be integrated into the existing Boolean as a robustness improvement;
- the output must still be reported as approximate while the refit remains;
- a future exact/certified Boolean must retain trimmed source curve spans (or exact split spans) instead of resampling/refitting them.

A future exact Boolean pipeline should be conceptually:

```text
exact prepared spans
    -> certified intersections/root intervals
    -> exact/certified span split parameters
    -> winding/topology classification with exact predicates
    -> retain/reverse source span pieces
    -> composite result without geometric refit
```

Do not label the existing Boolean result `EXACT_RATIONAL` simply because exact predicates were introduced.

---

## 16. Surface support

### 16.1 Tensor-product Bezier surfaces

**Strict rational:** yes.

Quantized control points become homogeneous integer control points. Evaluate by nested scaled De Casteljau in `u` and `v`.

Exact partial derivatives `Su`, `Sv`, `Suu`, `Suv`, `Svv`, and third-order homogeneous jets are rational and can be derived from derivative control nets.

### 16.2 B-spline/NURBS surfaces

**Strict rational:** yes, subject to width limits.

Prepare each non-empty `(u,v)` patch into an exact tensor-product homogeneous Bezier patch using exact local extraction in both axes.

Because width growth can be substantial for high degree in both axes, surfaces should use the same width planner and may require 1024/2048-bit kernels.

### 16.3 Bilinear/biquadratic patches

**Strict rational:** yes after position quantization.

They are polynomial and map naturally to the same exact patch representation.

### 16.4 Normal vectors

The unnormalized normal:

```text
N = Su x Sv
```

is rational if `Su` and `Sv` are rational.

The unit normal:

```text
n = N / sqrt(N.N)
```

is generally not rational.

Therefore expose/report these separately:

```text
unnormalized normal: EXACT_RATIONAL
unit normal:         APPROXIMATE (or future algebraic/certified)
```

The current `qaws_surface_core.h` explicitly normalizes using `sqrt`; the exact subsystem must not hide that domain transition.

### 16.5 Derived surfaces

Classify by construction, not by name.

Likely strict-rational candidates when all inputs are exact and only rational blending is used:

- ruled surface;
- simple extrusion by a rational vector;
- some Coons constructions;
- polynomial/rational lofts with exact coefficients.

Not strict-rational in general:

- revolution using trigonometric parameterization;
- pipe frames involving normalization;
- offset surfaces using unit normals;
- many fillet constructions;
- geodesic/closest-point constructions;
- algorithms containing iterative fitting or normalization.

Audit each implementation before enabling an exact capability bit.

---

## 17. Interaction with QAWS differentiation

### 17.1 Keep the existing public differentiation ABI

Do not replace `qaws_scalar` in:

- `qaws_field_view`;
- tangent/adjoint buffers;
- `qaws_curve_jet_*`;
- `qaws_surface_jet`;
- reports and conditioning metrics.

This would break the current C/GPU architecture and is unnecessary.

### 17.2 Exact spatial jets are still useful

For rational exact spans, the following are exact rational quantities at rational evaluation parameters:

- position;
- D1/D2/D3;
- surface partial derivatives before normalization;
- homogeneous basis/support coefficients;
- rational projection derivatives.

The exact subsystem should compute these internally and only round when returning through the existing floating API.

### 17.3 Exact JVP/VJP is mathematically possible

The current derivative rules for `+`, `-`, `*`, `/`, dot, cross, rational projection, etc. remain in the rational field if primal values and tangent/cotangent seeds are rational.

Finite f32/f64 tangent values can be decoded exactly as dyadic rationals.

Therefore a later phase may implement exact-internal:

- JVP;
- VJP;
- second directional derivatives;
- mixed second-order contractions.

The public output may still be f32/f64 unless a separate exact-value API is introduced.

### 17.4 Operations containing `sqrt` remain outside the rational exact domain

Current `qaws_dual_core.h` contains exact analytic derivative formulas for `sqrt` and normalization, but their numerical values are not generally rational.

Analytic differentiability and rational exactness are different properties.

Do not mark `qaws_dual1_sqrt`, normalized normals, curvature magnitudes, etc. as exact rational merely because their derivative formulas are analytic.

### 17.5 Exactness report + diff report

When exact/certified APIs and differentiation meet, return both concepts:

```text
numeric quality: EXACT_RATIONAL / CERTIFIED / APPROXIMATE

differentiability:
    class
    validity
    frozen state
    condition number
    residual
    branch gap
```

This is especially important at knots, active-set solves, and topology-dependent operations.

---

## 18. Proposed public API

Names below are proposed. Semantics are more important than exact spelling.

### 18.1 Header

Add:

```text
src/qaws/qaws_exact.h
```

and include it from `qaws.h` when `QAWS_ENABLE_EXACT` is enabled, or always expose declarations and compile stubs according to existing QAWS style.

### 18.2 Context and policy

```c
typedef struct qaws_exact_context qaws_exact_context;
typedef struct qaws_exact_curve qaws_exact_curve;
typedef struct qaws_exact_surface qaws_exact_surface;

typedef enum qaws_exact_fallback
{
    QAWS_EXACT_FALLBACK_FAIL = 0,
    QAWS_EXACT_FALLBACK_APPROXIMATE
} qaws_exact_fallback;

typedef struct qaws_exact_desc
{
    qaws_exact_space_desc space;

    unsigned int param_bits;       /* suggested default 24 */
    unsigned int weight_bits;      /* suggested default 24 */
    unsigned int max_integer_bits; /* e.g. 1024 or 2048 */

    qaws_exact_rounding rounding;
    qaws_exact_fallback fallback;
} qaws_exact_desc;
```

Creation:

```c
qaws_status qaws_exact_context_create(
    qaws_exact_desc const* desc,
    qaws_allocator const* allocator,
    qaws_exact_context** out_context);

void qaws_exact_context_destroy(qaws_exact_context* context);
```

### 18.3 Prepare an existing curve

```c
qaws_status qaws_exact_curve_prepare(
    qaws_exact_context const* context,
    qaws_curve const* curve,
    qaws_exact_curve** out_curve);

void qaws_exact_curve_destroy(qaws_exact_curve* curve);
```

Likewise for surfaces:

```c
qaws_status qaws_exact_surface_prepare(
    qaws_exact_context const* context,
    qaws_surface const* surface,
    qaws_exact_surface** out_surface);
```

Preparation either succeeds with a declared numeric quality or returns a specific exact status.

### 18.4 Report

```c
typedef enum qaws_exact_flag
{
    QAWS_EXACT_FLAG_NONE                    = 0,
    QAWS_EXACT_FLAG_INPUT_QUANTIZED         = 1 << 0,
    QAWS_EXACT_FLAG_PREP_QUANTIZED          = 1 << 1,
    QAWS_EXACT_FLAG_TRANSCENDENTAL_BOUNDARY = 1 << 2,
    QAWS_EXACT_FLAG_REFERENCE_FALLBACK      = 1 << 3
} qaws_exact_flag;

typedef struct qaws_exact_report
{
    qaws_numeric_quality quality;
    unsigned int flags;

    unsigned int storage_bits;
    unsigned int evaluation_bits;
    unsigned int predicate_bits;

    qaws_scalar max_position_quantization_error;
    qaws_scalar max_parameter_quantization_error;
    qaws_scalar max_weight_quantization_error;
} qaws_exact_report;
```

The report is descriptive. Correctness must not depend on callers inspecting it after the fact.

### 18.5 Floating result from exact internal evaluation

Provide an easy bridge first:

```c
qaws_status qaws_exact_curve_evaluate_2d(
    qaws_exact_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags,
    qaws_eval_result_2d* out_result,
    qaws_exact_report* out_report);

qaws_status qaws_exact_curve_evaluate_3d(
    qaws_exact_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags,
    qaws_eval_result_3d* out_result,
    qaws_exact_report* out_report);
```

Semantics:

- quantize/decode the input parameter according to the exact object's policy;
- choose the span with exact integer comparisons;
- evaluate exactly internally;
- round once when writing `qaws_scalar` output.

The returned floating coordinate is not itself an exact representation, but it is the correctly rounded/defined conversion of an exact internal rational result.

### 18.6 Low-level exact value access

Do not force a large fixed public `qaws_exact_int` ABI into V1 unless required by a caller.

For testing/debugging, provide formatting/query helpers such as:

```c
qaws_status qaws_exact_curve_evaluate_fraction_text(...);
```

or an opaque exact-value handle API later.

The main production value of the subsystem is robust geometry and certified decisions, not exposing 1024-bit fractions to every caller.

### 18.7 Certified predicates API

Examples:

```c
typedef enum qaws_exact_sign
{
    QAWS_EXACT_NEGATIVE = -1,
    QAWS_EXACT_ZERO = 0,
    QAWS_EXACT_POSITIVE = 1
} qaws_exact_sign;

qaws_status qaws_exact_compare_curve_coordinates(...);
qaws_status qaws_exact_orient2d(...);
qaws_status qaws_exact_point_equal(...);
```

High-level inspect/operation APIs should consume these internally rather than exposing every primitive publicly unless useful to users.

---

## 19. Proposed internal file layout

```text
src/qaws/
    qaws_exact.h
    qaws_exact.c

    exact/
        qaws_exact_int.h
        qaws_exact_int.c
        qaws_exact_convert.h
        qaws_exact_convert.c
        qaws_exact_filter.h
        qaws_exact_filter.c
        qaws_exact_quantize.h
        qaws_exact_quantize.c
        qaws_exact_types.h

        qaws_exact_bezier.h
        qaws_exact_bezier.c
        qaws_exact_extract.h
        qaws_exact_extract.c
        qaws_exact_curve_prepare.c
        qaws_exact_surface_prepare.c
        qaws_exact_predicates.c
        qaws_exact_derivatives.c

        qaws_exact_bounds_generated.h

tests/
    reference/
        qaws_bigint_ref.c
        qaws_bigint_ref.h
        qaws_rational_ref.c
        qaws_rational_ref.h

    58_exact_int.c
    59_exact_quantize.c
    60_exact_bezier.c
    61_exact_rational_bezier.c
    62_exact_bspline.c
    63_exact_nurbs.c
    64_exact_surfaces.c
    65_exact_predicates.c
    66_exact_diff.c
```

Numbers are illustrative; follow the repository's current test numbering when implementation begins.

Do not place the exact integer implementation inside `qaws_platform.h`. It is a different subsystem from C/HLSL/GLSL/Halide scalar abstraction.

---

## 20. Build-system integration

Add an independent option:

```cmake
option(QAWS_ENABLE_EXACT "Enable exact/certified CPU geometry" ON)
```

`QAWS_ENABLE_EXACT` must be independent of:

```text
QAWS_SCALAR_IS_FLOAT
QAWS_BACKEND
SIMD selection
```

Recommended CI additions:

- Linux f32 + exact
- Linux f64 + exact
- Windows f32 + exact
- Windows f64 + exact
- macOS at least one precision + exact
- existing shader/Halide compile tests unchanged
- a permutation that builds with `QAWS_ENABLE_EXACT=OFF`

Exact integer unit tests must produce identical integer results across compilers/platforms.

---

## 21. Reference arithmetic and verification

Production fixed-width arithmetic needs a mathematically trustworthy oracle.

Implement a test-only dynamically sized signed integer/rational reference under `tests/reference/`, or generate oracle vectors offline. Do not make a heavyweight external multiprecision dependency part of the public QAWS runtime.

Every fixed-width kernel should be tested against the reference for:

- random values within proven bounds;
- maximum positive/negative values;
- carry/borrow across every limb;
- signed multiplication;
- exact zero;
- equality after different common homogeneous scales;
- rational comparison;
- quantization edge cases;
- degenerate knot spans;
- repeated knots;
- degree 0 through `QAWS_CORE_MAX_DEGREE` (currently 16);
- derivative orders 0 through 3;
- 2D and 3D.

Debug builds should assert all theoretical width assumptions at runtime.

---

## 22. Tests required before enabling each family

A family must not claim exact support until all of these pass.

### 22.1 Exact value identity

For rational test data, compare exact numerator/denominator or homogeneous equality against the reference backend, not an epsilon.

### 22.2 Floating cross-check

Convert the exact result to f64 and compare against the existing f64 QAWS evaluator. This catches preparation/parameter-semantic mistakes but is not the exactness proof.

### 22.3 Endpoint and knot tests

Test:

- exact endpoints;
- every interior knot;
- just-left / just-right lattice neighbors;
- repeated knots;
- zero-length/invalid spans rejected consistently;
- span selection identical to documented QAWS conventions.

### 22.4 Scale invariance

Multiply every homogeneous control point in a span by a non-zero integer and prove all Cartesian results/predicates remain identical.

### 22.5 Width boundary tests

Construct geometry that reaches each proven bound and verify:

- the selected width succeeds;
- one configuration beyond the supported bound returns `QAWS_STATUS_EXACT_RANGE_EXCEEDED`;
- no wraparound is possible in release builds.

### 22.6 Filter correctness

For every certified floating filter:

- compare millions of random cases with exact sign;
- heavily sample near-zero cases;
- ensure every false-confidence case is impossible;
- measure fallback rate separately from correctness.

### 22.7 Differential identities

For exact-capable derivatives:

- derivative of endpoint formulas;
- exact quotient-rule identities;
- JVP/VJP dot-product identity when exact differential support is added;
- second-order contraction identities;
- finite differences may be used as a smoke test only, never as the oracle.

---

## 23. Performance model

The exact subsystem should be designed around the assumption that exact arithmetic is **not** the default dense sampling path.

Use ordinary f32/f64/GPU evaluation for:

- drawing/tessellation;
- dense preview;
- large sample batches;
- shader execution.

Use exact/certified arithmetic for:

- topological decisions;
- degeneracy handling;
- equality/order predicates;
- span/knot boundary decisions;
- robust trimming/Boolean/intersection certification;
- reference validation;
- explicit exact evaluation requested by a caller.

### 23.1 Prepare once

B-spline/NURBS Bezier extraction is a preparation cost. Cache prepared exact spans in the immutable exact object.

### 23.2 Filter first

Certified f64 filters should resolve the overwhelming majority of non-degenerate predicates without multi-limb fallback.

### 23.3 Width specialize only if benchmarks justify it

Begin with clear, correct fixed-capacity arithmetic. Add hand-specialized 128/256/512-bit paths only after profiling shows a meaningful bottleneck.

---

## 24. Implementation phases

### Phase 0 - Arithmetic foundation

Deliver:

- exact integer types and checked operations;
- f64 conversion with error bounds;
- rational comparison/cross-cancellation;
- quantization descriptors;
- test-only multiprecision reference;
- width-planning utilities;
- new exact status codes;
- exact build option.

Exit gate: fixed-width arithmetic is bit-exact against the reference under randomized and boundary testing.

### Phase 1 - Exact Bezier kernel

Deliver:

- integer homogeneous Bezier span representation;
- exact position evaluation 2D/3D;
- exact D1/D2/D3;
- ordinary Bezier preparation;
- rational Bezier preparation;
- Hermite exact preparation;
- polynomial exact preparation;
- exact span/parameter selection;
- exact-to-f32/f64 conversion.

Exit gate: exact Bezier/rational Bezier/Hermite tests pass degree-by-degree and derivative-by-derivative.

### Phase 2 - B-spline and NURBS

Deliver:

- integer knot quantization;
- local exact Bezier extraction matrix;
- exact B-spline span preparation;
- exact homogeneous NURBS span preparation;
- repeated-knot/boundary handling;
- extraction width proofs.

Exit gate: exact results match multiprecision reference for randomized knot vectors within the supported contract.

### Phase 3 - Certified predicates

Deliver:

- exact rational comparisons;
- orient2d and core geometric signs;
- f64 filter infrastructure;
- exact fallback;
- filter proof comments/tests;
- integration into selected inspection predicates.

Exit gate: no wrong sign in randomized/adversarial test corpus; performance fallback rate measured.

### Phase 4 - Tensor-product surfaces

Deliver:

- exact Bezier patch evaluation;
- exact partial derivatives;
- exact B-spline/NURBS patch extraction;
- exact unnormalized normal;
- explicit approximate boundary for unit normal;
- surface width planning up to configured max width.

### Phase 5 - Robust high-level operations

Incrementally replace tolerance-only topology decisions with certified predicates/root intervals.

Do **not** begin by claiming the current 2D Boolean is exact. First remove/refactor its sampled-refit output path where exact output is desired.

### Phase 6 - Exact-internal differentiation

Optional after the geometry kernel is stable:

- exact local support weights;
- exact internal JVP/VJP for rational rules;
- exact second-order directional rules;
- conversion to current `qaws_scalar` public derivative outputs;
- combined exactness + differentiability reporting.

---

## 25. Capability matrix for the first implementation

Legend:

- `YES`: strict rational exact after input lattice quantization.
- `PREP`: exact only after an approximate preparation result is explicitly quantized/frozen.
- `NO`: leaves rational domain in general.
- `LATER`: mathematically possible but not V1 priority.

| QAWS family / quantity | Strict rational | Notes |
|---|---:|---|
| Bezier curve | YES | Canonical direct path |
| Rational Bezier curve | YES | Integer homogeneous |
| Hermite curve | YES | Exact cubic conversion |
| Polynomial curve | YES | Exact Bernstein conversion |
| Uniform Catmull-Rom | YES | Rational cubic coefficients |
| Chordal Catmull-Rom | PREP | Distance uses sqrt |
| Centripetal Catmull-Rom | PREP | Parameterization leaves rational field |
| B-spline curve | YES | Exact local Bezier extraction |
| NURBS curve | YES | Homogeneous extraction |
| Composite | YES/LATER | If child maps are rational |
| Trajectory | YES/PREP | Depends on preparation branch |
| Arc | NO | sin/cos |
| Clothoid | NO | Fresnel/trig |
| Yuksel | PREP/NO | Current path includes root/trig operations |
| Bezier surface | YES | Tensor-product exact |
| B-spline surface | YES | Width may be large |
| NURBS surface | YES | Width may be large |
| Bilinear / biquadratic surface | YES | Polynomial |
| Unnormalized surface normal | YES | Cross product of rational derivatives |
| Unit normal | NO | sqrt normalization |
| Arc length | NO in general | Integral includes sqrt |
| Curvature scalar | NO in general | Norm/sqrt denominators |
| Span selection | YES | Integer parameter comparisons |
| Rational coordinate comparison | YES | Cross multiply |
| 2D orientation | YES | Exact determinant |
| Generic curve root/intersection coordinate | NO in general | Root may be algebraic |
| Certified root interval | LATER | Interval/Bernstein isolation |
| Current sampled 2D Boolean output | NO | Samples + refits |

---

## 26. Changes that must NOT be made

These are explicit anti-patterns for the implementation.

### 26.1 Do not add `QAWS_BACKEND_FIXED`

`QAWS_BACKEND` selects C/HLSL/GLSL/Halide execution semantics. Exact integer geometry is not another shader/backend target.

### 26.2 Do not redefine `qaws_scalar`

Do not add:

```c
typedef q26 qaws_scalar;
```

or a C++ rational class behind `qaws_scalar`.

It would break current scalar assumptions and would not make transcendental/high-level algorithms exact.

### 26.3 Do not round after every arithmetic operation

The only intended rounding boundary is input quantization (and final conversion back to f32/f64 when requested).

### 26.4 Do not use epsilon as equality in certified mode

Epsilon can trigger fallback, not decide exact equality.

### 26.5 Do not silently widen beyond the proven kernel contract

If the planner says a value needs more than the configured width, fail or use an explicitly selected reference fallback. Never truncate.

### 26.6 Do not label prepared approximations as exact source geometry

If a chordal/centripetal/trigonometric/fitted preparation is quantized, later evaluation can be exact relative to that prepared representation, but the report must contain `QAWS_EXACT_FLAG_PREP_QUANTIZED`.

### 26.7 Do not make existing GPU core headers depend on multi-limb exact types

Keep exact arithmetic C-only and isolated.

---

## 27. Repository-specific implementation touchpoints

The following files were inspected on branch `differentiable` and should guide implementation.

### `src/qaws/qaws_platform.h`

Current C/HLSL/GLSL/Halide scalar abstraction. Leave the backend model intact. Add no wide integer arithmetic here.

### `src/qaws/qaws_core_types.h`

Currently defines `QAWS_CORE_MAX_DEGREE = 16` and core vector/eval types. The exact width planner/tests should cover degrees through this maximum, although a requested exact configuration may legitimately exceed the configured max integer width and return a range status.

### `src/qaws/core/qaws_decasteljau_core.h`

Reference for current Bezier semantics and derivative-control-point construction. Implement a separate integer homogeneous De Casteljau kernel rather than templating this header into a giant multi-limb abstraction.

### `src/qaws/core/qaws_rational_bezier_core.h`

Confirms current rational Bezier already uses homogeneous weighted points and quotient-rule derivatives. Exact preparation should mirror this representation but keep homogeneous coordinates integer and delay projection.

### `src/qaws/core/qaws_bspline_basis_core.h`

Confirms runtime basis evaluation divides by knot differences. The exact path should canonicalize spans instead of running this scalar division recurrence with fixed-point rounding.

### `src/qaws/core/qaws_nurbs_eval_core.h`

Confirms current NURBS implementation evaluates B-spline basis derivatives, accumulates weighted numerator/denominator derivatives, and projects with quotient rules. Exact homogeneous Bezier span preparation replaces this entire runtime path for exact evaluation.

### `src/qaws/core/qaws_surface_core.h`

Confirms Bezier/B-spline surface evaluation is already tensor-product based and that normal output explicitly uses `sqrt`. Preserve rational partials separately from approximate normalization.

### `src/qaws/core/qaws_dual_core.h`

Confirms second-order analytic rules for multiplication/division/dot/cross/sqrt/normalization. Rational arithmetic can support the first group exactly; sqrt/normalize mark a rational-domain boundary.

### `src/qaws/internal/qaws_internal_types.h`

Current `qaws_curve` and family implementations store floating `qaws_scalar` arrays. Do not replace them. Exact data belongs in a separate prepared object.

### `src/qaws/qaws_eval.c` and `src/qaws/internal/qaws_internal_span.c`

Current public evaluation locates a span and computes normalized local parameters. The exact evaluator must reproduce the same clamping/boundary semantics with integer comparisons and exact local rational parameters.

### `src/qaws/qaws_bezier.c`, `qaws_hermite.c`, `qaws_bspline.c`

Useful family-specific reference for current parameter ranges, span conventions, derivative flags, and differential capabilities.

### `src/qaws/qaws_diff_types.h`

Keep its structured field/tangent/adjoint API. Add exactness metadata separately rather than replacing `qaws_scalar` differential views.

### `src/qaws/qaws_prepare.h/.c`

Existing precedent for C-only preparation feeding portable core evaluation. Exact preparation should follow the same philosophy with a different prepared representation.

### `src/qaws/qaws_boolean_2d.c`

Current Boolean samples kept arcs and B-spline-fits the result. This is a hard boundary: exact predicates alone cannot make this output exact.

---

## 28. Definition of done for the first mergeable milestone

The first exact milestone is complete when all of the following are true:

1. `QAWS_ENABLE_EXACT` builds cleanly in f32 and f64 C configurations.
2. Existing QAWS tests remain unchanged/green.
3. Exact arithmetic has no production dependency outside QAWS/C runtime.
4. Position/parameter/weight quantization is deterministic and documented.
5. Bezier and rational Bezier exact preparation/evaluation support degree 0..16 subject only to a documented max-width configuration.
6. D1/D2/D3 exact rational derivatives pass reference tests.
7. Hermite exact conversion passes identity tests against its source definition.
8. B-spline local Bezier extraction is exact for repeated and non-uniform knots.
9. NURBS homogeneous extraction/evaluation matches the reference rational implementation exactly.
10. Exact span boundary behavior matches QAWS's public parameter conventions.
11. At least rational compare, point equality, and orient2d have certified f64 filters with exact fallback.
12. Every exact kernel has a written bit-bound derivation and debug bound assertions.
13. Overflow/range exhaustion returns a specific status; no silent wrap or approximate fallback occurs in strict mode.
14. The exact report clearly distinguishes input quantization from approximate preparation quantization.
15. No HLSL/GLSL/Halide core API is changed to accommodate the exact subsystem.
16. Documentation contains a feature matrix stating what is and is not exact.

---

## 29. Recommended implementation order inside the codebase

For the developer, the shortest low-risk path is:

```text
1. qaws_exact_int + reference tests
2. quantization + exact space
3. integer homogeneous Bezier representation
4. scaled exact De Casteljau
5. exact rational projection + D1/D2/D3
6. qaws_exact_curve_prepare for Bezier/rational Bezier/Hermite
7. width planner + debug proofs
8. local exact B-spline Bezier extraction
9. NURBS homogeneous extraction
10. certified rational comparisons/orient2d filters
11. tensor-product surfaces
12. inspection/topology integration
13. optional exact-internal differentiation
```

Do not begin with Boolean/intersection rewrites. First establish a small exact kernel whose mathematical contract is easy to prove and exhaustively test.

---

## 30. Final architectural summary

The intended QAWS numeric architecture after this work is:

```text
                           QAWS geometry API
                                 |
               +-----------------+------------------+
               |                                    |
               v                                    v
      ordinary runtime                        exact preparation
   qaws_scalar f32/f64                   quantize once to lattices
               |                                    |
       +-------+-------+                    canonicalize spans
       |               |                             |
       v               v                             v
      CPU       HLSL/GLSL/Halide        integer homogeneous Bezier
       |               |                             |
       +-------+-------+                 +-----------+-----------+
               |                         |                       |
               v                         v                       v
       fast approximate eval        exact rational eval     certified predicates
                                                      f64 filter -> exact fallback
```

The key separation is deliberate:

- Floating QAWS remains fast, portable, SIMD/GPU-friendly, and differentiable.
- Exact QAWS is CPU-side, fixed-width, deterministic, and used where numerical ambiguity can change geometry/topology.
- Rational families meet in one canonical exact representation instead of each getting a bespoke fixed-point evaluator.
- Non-rational operations are reported explicitly rather than being given a false "exact" label.

This provides QAWS with a robust exact foundation without compromising the architecture already present in the current repository.
