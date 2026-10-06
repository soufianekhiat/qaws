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

## Capability matrix

| quantity | status |
|---|---|
| orient2d / orient3d / ratio compare on doubles | certified |
| Bezier, rational Bezier: position, C', C'', C''' | exact rational |
| B-spline, NURBS: position, C', C'', C''' | exact rational (budget above) |
| span selection | exact (integer knots) |
| winding number of a closed loop | certified |
| Hermite, polynomial, uniform Catmull-Rom, composites | planned |
| chordal / centripetal Catmull-Rom, Yuksel | exact only after a frozen preparation (planned) |
| arcs, clothoids (sin/cos, Fresnel) | not rational |
| unit normals, curvature, arc length (sqrt) | not rational |
| tensor-product surfaces | planned |
| curve intersections | certified root intervals planned (roots are algebraic in general) |

Tests: 63 (integers), 64 (predicates), 65 (Bezier), 66 (winding), 67
(B-spline / NURBS), all against Mathematica exact references. Figures:
`examples/exact_showcase.c`.
