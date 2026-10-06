# Differentiation

Qaws computes exact derivatives of curve and surface evaluation, of differential geometry, of geometric solves (closest points, intersections, extrema), of geometry-building operations (split, join, conversions, fitting) and of integral functionals (length, bending, area, thin plate, Willmore). Derivatives are formed analytically from the same basis functions the evaluator uses; the library never differentiates by finite steps.

All declarations live in `qaws_diff_types.h`, `qaws_diff.h`, `qaws_diff_geometry.h`, `qaws_diff_ops.h`, `qaws_diff_map.h` and `qaws_diff_functionals.h`, all included by `qaws.h`. The per-function reference is in [api/diff.md](api/diff.md).

## Vocabulary

| Term | Meaning |
|---|---|
| **tangent** | forward derivative along one direction: `out = J * in` |
| **adjoint** | reverse derivative: `in += J^T * out`, accumulated into caller buffers |
| **tangent2** | second derivative along the same direction: `d^2 out / de^2` |
| **HVP** | Hessian-vector product of a scalar functional, provided directly where the output is linear in the fields, otherwise composed by the caller |

Tangents overwrite their outputs. Adjoints always accumulate (`+=`), so several contributions can be pulled back into the same buffers.

## Parameters are fields

A curve or surface exposes its parameters as named **fields**. `qaws_curve_describe_fields` / `qaws_surface_describe_fields` list them with a `qaws_field_desc`:

- `field`: `QAWS_FIELD_CONTROL_POINTS`, `QAWS_FIELD_WEIGHTS`, `QAWS_FIELD_KNOTS`, `QAWS_FIELD_U_KNOTS`, `QAWS_FIELD_V_KNOTS`, `QAWS_FIELD_POINTS`, `QAWS_FIELD_DERIVATIVES`, `QAWS_FIELD_COEFFICIENTS`, `QAWS_FIELD_CENTER`, `QAWS_FIELD_OFFSET_DISTANCE`, `QAWS_FIELD_DIRECTION`, ...
- `value_type` and `count`: scalar, vec2 or vec3, and the number of elements.
- `domain` and `constraint`: what the values mean (position, weight, parametric) and their admissible set (knots are `QAWS_CONSTRAINT_MONOTONIC`, weights are positive).
- `diff_class`: `QAWS_DIFF_SMOOTH`, `QAWS_DIFF_PIECEWISE_SMOOTH` (smooth inside knot spans), `QAWS_DIFF_ACTIVE_SET`, ... or `QAWS_DIFF_UNSUPPORTED`.
- `capabilities`: `QAWS_CAP_TANGENT`, `QAWS_CAP_ADJOINT`, `QAWS_CAP_TANGENT2`, `QAWS_CAP_LOCAL_SUPPORT`, `QAWS_CAP_LINEAR`.

`qaws_curve_read_field` copies the primal values; `qaws_curve_clone_with_fields` builds a new object of the same family with some fields replaced.

### Views

Derivatives are exchanged through **views** over caller storage. A `qaws_field_view` points at one field's tangent or adjoint values; a `qaws_diff_views` groups the views of one object. A field without a view is inactive (its tangent is zero and its adjoint is not written).

```c
qaws_scalar cp_dot[7 * 3] = { 0 };
qaws_field_view fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, cp_dot, 7, 3);
qaws_diff_views views = { &fv, 1, NULL, 0 };
```

Two masks select parameters without copying:

- `fv.active`: one byte per element, 0 = inactive. Used for "boundary fixed, interior free".
- `fv.component_mask`: bit `c` set = component `c` active. Used for "only heights move".

Inactive entries read as zero in tangents and receive nothing in adjoints. `stride` supports interleaved storage.

### Parameter keys

`qaws_param_key` names a single scalar parameter through derived objects: `child[0]/control_points/3/y`. `qaws_param_key_to_string` and `qaws_param_key_parse` convert keys to and from that text, so optimizers and editors can store stable parameter paths.

## Jets

Evaluations are returned as **jets**.

- Curves: `qaws_curve_jet_2d` / `qaws_curve_jet_3d` hold `d[0..3]` = position and the first three derivatives. `channels` uses the `QAWS_EVAL_FLAG_*` bits.
- Surfaces: `qaws_surface_jet` holds the ten partials up to total order three, at index `t(t+1)/2 + b` for `d^(a+b) S / du^a dv^b` with `t = a + b`: `S, Su, Sv, Suu, Suv, Svv, Suuu, Suuv, Suvv, Svvv`. `channels` uses `QAWS_SJET_*` bits (`QAWS_SJET_ORDER1/2/3` for full orders).

A tangent jet is the derivative of every requested partial along one direction; an adjoint jet holds one adjoint per partial.

## Tangents

```c
qaws_curve_jet_3d p, dp;
qaws_scalar cp_dot[7 * 3] = { 0 };
qaws_field_view fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, cp_dot, 7, 3);
qaws_diff_views views = { &fv, 1, NULL, 0 };

cp_dot[2 * 3 + 1] = 1;   /* move control point 2 along y */
qaws_curve_eval_tangent_3d(NULL, curve, t, 0, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1,
    &views, &p, &dp);
/* dp.d[0]: how C(t) moves, dp.d[1]: how C'(t) moves */
```

`t_tangent` adds motion of the evaluation coordinate itself; its contribution is the next spatial derivative. Batch variants (`qaws_curve_eval_batch_tangent_*`, `qaws_surface_eval_batch_tangent`) take arrays of coordinates and coordinate tangents.

## Adjoints

```c
qaws_curve_jet_3d bar[N];        /* d loss / d jet for every sample */
qaws_scalar cp_bar[7 * 3] = { 0 };
qaws_scalar t_bar[N] = { 0 };
qaws_field_view fv = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, cp_bar, 7, 3);
qaws_diff_views grad = { &fv, 1, NULL, 0 };

qaws_curve_eval_batch_adjoint_3d(&ctx, curve, ts, N, QAWS_EVAL_FLAG_POSITION, bar, &grad, t_bar);
```

Adjoints are batch first. `qaws_diff_context.accumulation` selects how sample contributions reach the parameter buffers; all three give the same result:

| Strategy | Behaviour |
|---|---|
| `QAWS_ACCUMULATE_SCATTER` | every sample adds directly into the parameter adjoints |
| `QAWS_ACCUMULATE_TILED` | tiles of `tile_size` samples accumulate into a dense local buffer, flushed once per tile (the CPU analogue of workgroup shared memory) |
| `QAWS_ACCUMULATE_GATHER` | contributions are bucketed per parameter (counting sort) and each parameter sums its own bucket: deterministic, conflict free |

`qaws_curve_local_support` / `qaws_surface_local_support` return the exact basis weights of one evaluation, and `qaws_*_build_support_index` builds the sample-to-parameter index for custom kernels.

## Second order

`qaws_curve_eval_batch_tangent2_*` and `qaws_surface_eval_batch_tangent2` return the second derivative along the combined direction (parameters, coordinates). They are available where `QAWS_CAP_TANGENT2` is set, which includes rational families and knots. Integral functionals provide a direct HVP (`qaws_*_functional_hvp`) for families that are linear in their fields (`QAWS_CAP_LINEAR`); for the others, compose HVPs from tangents and adjoints. Second order is selective by design: only primitives where it is cheap and exact provide it.

## Families

| Family | Fields with derivatives |
|---|---|
| Bezier, Hermite, polynomial | control points, points and derivatives, coefficients (linear) |
| B-spline curve | control points (linear), knots |
| NURBS curve, rational Bezier | control points, weights (quotient rule), knots (NURBS) |
| Catmull-Rom (uniform, chordal, centripetal), Yuksel C2 (Bezier mode) | interpolated points (non-linear: parameterization and sub-curve parameter) |
| Clothoid | origin (`CENTER`), start angle, start curvature and curvature rate; the position integral is differentiated exactly in dual numbers |
| Composite | its segments as children |
| Bezier, B-spline, bilinear, biquadratic surfaces | control points (linear); B-spline also U and V knots |
| NURBS surface | control points, weights, U and V knots |
| Extrusion, ruled, revolution, Coons, loft, Gordon, offset surfaces | their own fields (direction, center, angle, distance, ...) plus their input curves and surfaces as children |

Rational families use homogeneous jets and the quotient recurrence `C^(k) = (A^(k) - sum binom(k,i) W^(i) C^(k-i)) / W`, forward and reverse. Derived objects chain into their inputs: the tangent views of input `i` live in `views->children[i]`, and parameter keys carry the matching `child[i]` prefix (`qaws_curve_diff_children` / `qaws_surface_diff_children` list them).

### Knots

Knots are differentiable fields of B-spline and NURBS curves and surfaces. A knot direction re-evaluates the basis rows in dual numbers (first and second order), so tangents, second tangents and adjoints include knots. The knot span containing each sample is frozen (`QAWS_FREEZE_SPAN`): moving a knot past a sample is a discrete event. Curve functionals differentiate their quadrature exactly as a function of the knots, including the motion of span boundaries; surface functionals integrate on a cell grid over the domain, so interior knots change the integrand at fixed nodes while the end knots of the domain move every node and scale the weights (value, tangent, second tangent and gradient; no direct HVP for knots). Operations returning differential maps hold the knots of their inputs fixed.

## Geometry

`qaws_diff_geometry.h` turns jets into differential geometry with derivatives: speed, unit tangent, normal (and binormal), curvature and torsion for curves; normal, first and second fundamental forms, Gaussian, mean and principal curvatures for surfaces. Each has an eval (value, tangent and tangent2 from primal, tangent and tangent2 jets) and an adjoint (pullback onto the jet), so it composes with any evaluation. Degenerate configurations (zero speed, umbilics for principal curvatures) are reported through the validity.

## Implicit operations

Operations defined by a solve are differentiated with the implicit function theorem, never by differentiating solver iterations: for `g(x, theta) = 0`, `dx/dtheta = -(dg/dx)^-1 dg/dtheta`.

| Operation | Condition |
|---|---|
| closest point on a curve / surface | `C' . (C - q) = 0`, `Su . (S - q) = Sv . (S - q) = 0` |
| curve pair (closest approach, intersection) | `A' . (A - B) = 0`, `B' . (A - B) = 0` |
| curve crossing a plane | `n . (C - p) = 0` |
| curve piercing a surface | `S(u, v) - C(t) = 0` |
| extremum along a direction | `e . C'(t) = 0` |
| inflection (xy) | `x'y'' - y'x'' = 0` |

Closest points search globally (sampled profile, refined by Newton). The other operations differentiate one solution: pass a seed from the discrete finders (`qaws_curve_find_intersections_*`, `qaws_curve_find_plane_intersections`, `qaws_surface_find_curve_intersections`, `qaws_curve_find_extrema`, `qaws_curve_find_inflection_points`); the solution is refined by Newton and differentiated. Which solutions exist is a discrete choice (`QAWS_FREEZE_ACTIVE_SET`). Inputs that are not objects (query point, plane, direction) have their own tangent and adjoint arguments.

Each solve reports:

- `condition_number`: conditioning of `dg/dx`;
- `residual`: relative residual at the returned solution;
- `branch_gap`: distance gap to the best competing solution (closest points);
- `validity`: `QAWS_DIFF_ILL_CONDITIONED` near singular configurations (tangential intersections, flat extrema), `QAWS_DIFF_AMBIGUOUS` when a competing solution is within 0.5% of the distance, `QAWS_DIFF_VALID_LOCALLY` when a coordinate sits on its domain boundary and is held fixed.

## Differential maps

Operations that build new geometry return, next to their outputs, a `qaws_diff_map` from the parameters of their inputs to the parameters of their outputs:

```c
qaws_curve *left, *right;
qaws_diff_map* map;
qaws_curve_split_diff(curve, 0.4f, &left, &right, &map);
/* inputs: 0 = curve, 1 = split parameter; outputs: 0 = left, 1 = right */
qaws_diff_map_tangent(map, NULL, in_tangents, 2, out_tangents, 2);  /* out = J in */
qaws_diff_map_adjoint(map, NULL, out_adjoints, 2, in_adjoints, 2);  /* in += J^T out */
qaws_diff_map_destroy(map);
```

Available: split, join, Hermite to Bezier, Bezier to B-spline, B-spline to NURBS, degree elevation and reduction, the least-squares B-spline fit and the 3D offset (constant direction or Frenet normal; the sample set is frozen). Polynomial constructions are exact linear maps; NURBS maps include the homogeneous projection. `qaws_curve_fit_bspline_diff` differentiates the whole fit: chord-length parameters, averaged knots, basis functions and normal equations, with respect to the data points and (when given) the sample parameters. `qaws_diff_map_get_entries` exposes the sparse Jacobian.

## Integral functionals

`qaws_diff_functionals.h` integrates with Gauss-Legendre quadrature and differentiates the quadrature exactly:

- curves: `QAWS_FUNCTIONAL_LENGTH`, `QAWS_FUNCTIONAL_BENDING` (`int |C''|^2 dt`), `QAWS_FUNCTIONAL_CURVATURE_SQUARED` (`int kappa^2 ds`);
- surfaces: `QAWS_FUNCTIONAL_AREA`, `QAWS_FUNCTIONAL_THIN_PLATE`, `QAWS_FUNCTIONAL_WILLMORE`.

Each has eval (value, tangent, tangent2), gradient (accumulated into views) and HVP.

## Inverse-CDF sampling

`qaws_diff_sampling.h` samples a curve by inverting a measure `M(t) = int rho |C'|`: arc length (constant speed), curvature-weighted (`rho = sqrt(floor^2 + kappa^2)`) or a user density field in space. The sample at `sigma = distance + fraction * M_total` has a parameter `t` defined by `M(t) = sigma`, differentiated implicitly to first and second order (`tangent`, `tangent2`), pulled back as an adjoint (one quadrature pass for a whole batch, plus distance adjoints) and as a Hessian-vector product for linear families.

Time traversal maps onto the same machinery: `qaws_traversal_cdf_targets` turns times (easing, motion profile, wrap mode) into arc-length targets with their speed and acceleration, so samples of a traversal have exact first and second derivatives in time and in the curve parameters.

Surfaces warp points of the unit square by the marginal and conditional inverse CDFs of `rho(S) |S_u x S_v|` (area or a density field): stratified or blue-noise points stay stratified on the surface. The forward pass solves the discrete equations in dual numbers (first and second order), the adjoint pulls back through the conditional then the marginal equation, and the HVP polarizes the second order forward pass.

## Context and reports

`qaws_diff_context` carries the accumulation strategy, tile size, scratch allocator and an optional `qaws_diff_report`. Initialize with `qaws_diff_context_init` and reset the report with `qaws_diff_report_reset`. The report collects, over all samples of a call, the worst `diff_class` and `validity`, the frozen discrete states the result relies on (`frozen_used`), the largest condition number and residual, the smallest branch gap and the index of the worst sample. A result that relies on a frozen state is exact for that state; the report tells the caller when the state matters.

## Verification

The test suites `49_diff_model` to `59_diff_surface_sampling` check every rule in single and double precision:

- tangents against central differences of rebuilt objects (finite differences appear only in tests);
- adjoints through the identity `<ybar, J xdot> = <J^T ybar, xdot>`, for all three accumulation strategies;
- second tangents against differences of tangents;
- implicit solves against re-solved perturbed problems.
- analytic ground truth from Mathematica where a closed derivation exists (`tests/reference/*.wls`, 30 digits): the inverse-CDF samples (arc length, curvature, density), their first and second tangents, gradient and HVP match to 1e-11.

## Examples

`examples/diff_showcase.c` (CMake option `QAWS_BUILD_EXAMPLES`, target `qaws_diff_showcase`) writes SVG figures to `showcase/`:

1. B-spline fit to data with adjoints
2. Sensitivity of a curve to each control point
3. Fairing with bending energy
4. NURBS weight optimization
5. Surface fit and Gaussian curvature sensitivity
6. Vase from a revolved profile
7. Coons patch fairing through its boundary curves
8. Gordon and loft surfaces through curve networks
9. Projection fit with closest points and a validity map
10. Soap film: area minimization with Newton steps on HVPs
11. Parameter correction through the differentiable fit, with knot derivatives
12. Arch design from intersection and extremum adjoints
13. Knot placement: fitting a crease by moving surface knots
14. Constant-speed sampling: first and second order sample tangents, and a fit of the samples by Newton-CG with exact Hessian-vector products
15. Inverse-CDF sampling under arc length, curvature and a density field, with the sample tangents
16. Surfaces: stratified points warped onto a patch by the area and density inverse CDFs
17. Blue noise on a patch: Gaussian repulsion of the warped samples descended through the xi adjoints of the warp

## GPU kernels

`core/qaws_bspline_diff_core.h` provides the B-spline tangent and adjoint kernels for every backend (C, HLSL, GLSL, Halide): `qaws_bspline_tangent_3d`, `qaws_bspline_adjoint_cp_3d` (the per-element term of a gather pass) and `qaws_bspline_adjoint_t_3d`. `examples/diff_bspline_adjoint.hlsl` runs them as compute shaders: one thread per sample for tangents and parameter adjoints, and one thread per control point for the control point adjoints. The gather pass loops over the samples listed by `qaws_curve_build_support_index`, so it needs no atomics and its result is deterministic. Test `57_diff_core` checks the kernels against the C runtime.

`core/qaws_bspline_sampling_core.h` does the same for inverse-CDF sampling and integral functionals of B-spline curves. The integrands are the arc-length and curvature measures and the length, bending and curvature-squared functionals, all on second order dual numbers along the control point tangents. The kernels are:

- `qaws_bspline_integrate`: integral over one span, with its first and second rates (the default composite rule).
- `qaws_bspline_cdf_solve`: safeguarded Newton with a fixed iteration count, branch free.
- `qaws_bspline_cdf_tangent`: first and second sample tangents by the implicit function theorem.
- `qaws_bspline_cdf_lambda`: the adjoint multiplier.
- `qaws_bspline_node_adjoint_cp`: the per-node term of a gather over control points.

`examples/diff_cdf_sampling.hlsl` chains them as compute passes:

1. One thread per span computes the measures.
2. A prefix sum combines them.
3. One thread per sample solves it and computes its first and second tangents.
4. One thread per sample computes its multiplier.
5. One thread per control point gathers the adjoint over the nodes and samples of the spans it supports.

Test `60_diff_sampling_core` drives the kernels the same way on the C backend. It reproduces the runtime samples, tangents, adjoints, functional values and gradients to round-off.

`core/qaws_bspline_surface_functional_core.h` covers the area, thin-plate and Willmore functionals of B-spline surfaces. The runtime integrates on a uniform grid of cells with 4 x 4 Gauss points each, and the cells need not align with the knots, so every node is evaluated on the knot window of its own spans:

- `qaws_surface_cell_node` gives a node's parameters and weight.
- `qaws_bspline_surface_node` gives the weighted integrand with its first and second rates along the control point tangents. One thread per node, then a sum.
- `qaws_bspline_surface_node_adjoint_cp` is the per-node term of a gather over control points.

`examples/diff_surface_functionals.hlsl` runs them as compute passes. Test `61_diff_surface_functional_core` reproduces the runtime values, tangents and gradients to round-off.

`core/qaws_bspline_surface_sampling_core.h` warps points of the unit square onto a B-spline surface by its area measure. It is the GPU form of `qaws_surface_cdf_sample_*` with `QAWS_MEASURE_AREA`. The curvature measure needs fourth derivatives, beyond the kernels' `QAWS_CORE_MAX_DERIV`; the density measure needs a callback. The whole control grid is passed, since a v-line crosses several knot spans. The kernels are:

- `qaws_patch_cell_mass`: the u cells, followed by a prefix sum.
- `qaws_patch_cdf_solve`: the samples, with a fixed iteration count.
- `qaws_patch_cdf_tangent`: first and second tangents along the control point and xi tangents, by dual Newton steps on the discrete equations.
- `qaws_patch_cdf_multipliers`: lambda, mu and the xi adjoints.
- `qaws_patch_line_adjoint_cp`: the per-line term of the control point gather. The lines are the sample's conditional lines, the lines of its partial u cell and those of the full cells, with the area gradient `dw/dS_u = S_v x n`, `dw/dS_v = n x S_u`.

`examples/diff_surface_sampling.hlsl` runs them as compute passes. Test `62_diff_surface_sampling_core` reproduces the runtime samples, tangents and adjoints to round-off.

`examples/diff_applications.c` (target `qaws_diff_applications`) applies the API to real data:

1. Hair strands from a photo: structure-tensor orientation field, evenly spaced streamlines, B-spline strands aligned through unit-tangent adjoints.
2. Photo vectorization: isocontours become centripetal Catmull-Rom splines whose interpolation points are optimized onto the contours.
3. Triangulated mesh to six bicubic patches with ADMM: exact per-patch least squares at foot points (local support weights, thin-plate HVPs), seam control points in consensus.
4. Non-rigid registration: CMA-ES basin hopping on the pose with gradient refinements inside the fitness.
5. Hair grooming: 3D strands with fixed roots under length and bending functionals, gravity, head collision and unit-tangent alignment to a combing field.
6. Hair from a photo (single-view hair modeling): the projected unit tangents of the 3D strands follow the photo orientation field and stay on its hair mask.
7. Rational patches: NURBS weights freed after the ADMM fit, seams kept closed.
8. Real haircuts: vector hair strands from six portrait photos (seeded hair mask, orientation field, B-spline strands).
9. 3D hair from frontal portraits: head placed from the face, visible strands follow the photo, hidden strands follow the 3D priors.
10. Single-view hair modeling on plain-background portraits (long, curly, ponytail, profile bob), after Chai et al. 2012/2013 and Hu et al. 2015: background flooded from the border; Gabor orientation (32 angles, refinement pass); image strands lifted onto a hair volume (head plus inflated silhouette) and fitted as 3D B-splines through the batch adjoint; their tangents constrain a voxel orientation-tensor field diffused through the volume; 14k strands grown from the scalp with RK2 (gathered toward a tie for the ponytail), fill strands in empty cells, per-photo curl as helical offsets; a small rasterizer renders them with Kajiya-Kay shading and density self-shadowing.
11. Font outline fitting: glyphs rendered from an installed font (`examples/glyph_to_ppm.ps1`) are traced at the half-ink isoline; each contour becomes a periodic cubic B-spline with triple knots at detected corners and arc-length knots elsewhere, fitted through the batch adjoint with parameter correction.
12. Shape optimization, 2D vesicles: a closed cubic B-spline minimizes its bending energy at fixed length and reduced area (augmented Lagrangian) by Newton-CG on exact Hessian-vector products of the curvature-squared, bending and length functionals and of the enclosed area.
13. CAD reverse engineering with G1 seams: the six-patch ADMM fit of the blob scan is refined as one model, with shared seam control points (C0 exact) and `1 - (n_a . n_b)^2` at seam samples pulled back through the surface adjoint on the `S_u`, `S_v` channels. The largest crease drops from 20 degrees to 0.24 (cube corners, where three patches meet, excluded). Zebra stripes show the result.

Photos are read as PPM; `examples/photo_to_ppm.ps1` converts any image (`-Crop "x,y,w,h"` in fractions).
