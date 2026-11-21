# Differentiation

Qaws computes exact derivatives of curve and surface evaluation, of local differential geometry, of geometric solves (closest points, intersection points, extrema, inflections), of geometry-building operations (split, join, conversions, degree change, fitting, offset) and of inverse-CDF sampling. Derivatives are formed analytically from the same basis functions the evaluator uses; no numerical differentiation is ever returned. Finite differences appear only in the tests that check the rules.

The declarations live in `qaws_diff_types.h`, `qaws_diff.h`, `qaws_diff_geometry.h`, `qaws_diff_ops.h`, `qaws_diff_map.h` and `qaws_diff_sampling.h`, all included by `qaws.h`. The per-function reference is [api/diff.md](api/diff.md). Like the rest of the library, every function exists in single and double precision ([precision](precision.md)).

## What is differentiated

Two independent axes:

- the **evaluation coordinate** (`t` for curves, `u`, `v` for surfaces): handled by spatial jets, which already carry the derivatives in the coordinate;
- the **parameters** of the object (control points, weights, knots, radius, ...): handled by tangents and adjoints.

A tangent jet is therefore the derivative, along a parameter (and optionally coordinate) direction, of every spatial derivative in the jet.

## Vocabulary

| Term | Meaning |
|---|---|
| **tangent** | forward derivative along one direction: `y' = J x'` |
| **adjoint** | reverse derivative: `x_bar += J^T y_bar`, accumulated into caller buffers |
| **tangent2** | second directional derivative along the same direction |
| **HVP** | Hessian-vector product of a scalar objective, provided for inverse-CDF sampling; elsewhere composed by the caller from tangents and adjoints |

Tangents overwrite their outputs. Adjoints always accumulate (`+=`), so several contributions can be pulled back into the same buffers; `qaws_diff_views_clear` zeroes them first when needed.

Forward mode suits few directions and many outputs (sensitivity of a whole curve to one control point); reverse mode suits one scalar loss and many parameters (fitting). Both are exposed for every rule, so the caller chooses per problem.

## Parameters are fields

An object exposes its parameters as typed **fields** rather than one flat vector. `qaws_curve_describe_fields` / `qaws_surface_describe_fields` list them as `qaws_field_desc`:

- `field`: `QAWS_FIELD_CONTROL_POINTS`, `QAWS_FIELD_WEIGHTS`, `QAWS_FIELD_KNOTS`, `QAWS_FIELD_U_KNOTS`, `QAWS_FIELD_V_KNOTS`, `QAWS_FIELD_POINTS`, `QAWS_FIELD_DERIVATIVES`, `QAWS_FIELD_COEFFICIENTS`, `QAWS_FIELD_CENTER`, `QAWS_FIELD_ANGLE_START`, `QAWS_FIELD_ANGLE_END`, `QAWS_FIELD_OFFSET_DISTANCE`, `QAWS_FIELD_CURVATURE`, `QAWS_FIELD_CURVATURE_RATE`, `QAWS_FIELD_DIRECTION`, `QAWS_FIELD_PARAMETER`, ... (numeric values are stable and append only).
- `value_type` and `count`: scalar, vec2 or vec3, and the number of elements.
- `domain` and `constraint`: what the values mean (position, weight, parametric, ...) and their admissible set (knots are `QAWS_CONSTRAINT_MONOTONIC`, weights `QAWS_CONSTRAINT_POSITIVE`), with suggested bounds `default_lower` / `default_upper`.
- `diff_class`: `QAWS_DIFF_SMOOTH`, `QAWS_DIFF_PIECEWISE_SMOOTH` (smooth inside knot spans or segments), `QAWS_DIFF_ACTIVE_SET`, `QAWS_DIFF_SUBGRADIENT`, `QAWS_DIFF_DISCRETE` or `QAWS_DIFF_UNSUPPORTED`.
- `capabilities`: `QAWS_CAP_TANGENT`, `QAWS_CAP_ADJOINT`, `QAWS_CAP_TANGENT2`, `QAWS_CAP_LOCAL_SUPPORT`, `QAWS_CAP_LINEAR` (output linear in all fields).

`qaws_curve_read_field` / `qaws_surface_read_field` copy the primal values; `qaws_curve_clone_with_fields` builds a new curve of the same family with some fields replaced (an optimizer step).

Why fields: the caller never packs a parameter vector and the library never needs to know the optimizer's layout. Flattening, bounds and optimizer variables belong to the caller, who can map fields to any solver.

### Views

Derivatives are exchanged through **views** over caller storage. A `qaws_field_view` addresses one field's tangent or adjoint values (element `i`, component `c` at `data[i * stride + c]`); a `qaws_diff_views` groups the views of one object and, for derived objects, the views of its children. A field without a view is inactive: its tangent is zero and its adjoint is not written.

```c
qaws_scalar cp_dot[7 * 3] = { 0 };
qaws_field_view fv = qaws_field_view_make(cp_dot, 7, QAWS_FIELD_CONTROL_POINTS, 3);
qaws_diff_views views = { &fv, 1, NULL, 0 };
```

Two masks select parameters without copying:

- `fv.active`: one byte per element, 0 = inactive ("boundary fixed, interior free").
- `fv.component_mask`: bit `c` set = component `c` active ("only heights move").

Inactive entries read as zero in tangents and receive nothing in adjoints. `stride` supports interleaved storage.

### Parameter keys

`qaws_param_key` names one scalar parameter through derived objects: `child[0]/control_points/3/y`. `qaws_param_key_to_string` and `qaws_param_key_parse` convert keys to and from that text, so optimizers and editors can store stable parameter paths.

### Derived objects

A derived object (composite curve, extrusion, ruled, revolution, Coons, loft, Gordon, offset surface, ...) chains its rules into the objects it references. `qaws_curve_diff_children` / `qaws_surface_diff_children` list them; child `i` is `views->children[i]` in tangent and adjoint storage and `child[i]` in parameter keys. Derivatives with respect to a profile curve of a revolved surface are thus reached through the surface's own calls.

## Jets

Evaluations are returned as **jets**:

- Curves: `qaws_curve_jet_2d` / `qaws_curve_jet_3d` hold `d[0..3]`, the position and the first three derivatives. `channels` uses the `QAWS_EVAL_FLAG_POSITION/D1/D2/D3` bits.
- Surfaces: `qaws_surface_jet` holds the ten partials up to total order three, `d[0..9] = S, Su, Sv, Suu, Suv, Svv, Suuu, Suuv, Suvv, Svvv`. `channels` uses the `QAWS_SJET_*` bits (`QAWS_SJET_ORDER1/2/3` for full orders). `qaws_surface_eval_jet` evaluates a primal jet analytically.

A tangent jet holds the derivative of every requested partial along one direction; an adjoint jet holds one adjoint per partial. Jets are the common currency: geometry rules read jets, and sampling and implicit operations build on them.

## Tangents

```c
qaws_curve_jet_3d p, dp;
qaws_scalar cp_dot[7 * 3] = { 0 };
qaws_field_view fv = qaws_field_view_make(cp_dot, 7, QAWS_FIELD_CONTROL_POINTS, 3);
qaws_diff_views views = { &fv, 1, NULL, 0 };

cp_dot[2 * 3 + 1] = 1;   /* move control point 2 along y */
qaws_curve_eval_tangent_3d(&p, &dp, curve, &views, t,
    0, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, NULL);
/* dp.d[0]: how C(t) moves, dp.d[1]: how C'(t) moves */
```

`t_tangent` adds motion of the evaluation coordinate itself; its contribution is the next spatial derivative. The batch functions (`qaws_curve_eval_batch_tangent_2d/3d`, `qaws_surface_eval_batch_tangent`) take arrays of coordinates and coordinate tangents and are the primary API; the single-sample functions wrap them.

## Adjoints

```c
qaws_curve_jet_3d bar[N];        /* d loss / d jet for every sample */
qaws_scalar cp_bar[7 * 3] = { 0 };
qaws_scalar t_bar[N] = { 0 };
qaws_field_view fv = qaws_field_view_make(cp_bar, 7, QAWS_FIELD_CONTROL_POINTS, 3);
qaws_diff_views grad = { &fv, 1, NULL, 0 };

qaws_curve_eval_batch_adjoint_3d(&grad, t_bar, curve, ts, bar, N, QAWS_EVAL_FLAG_POSITION, &ctx);
```

Adjoints are batch first: one shared object, many samples, one pass. `qaws_diff_context.accumulation` selects how sample contributions reach the parameter buffers; all three give the same sum up to floating-point ordering:

| Strategy | Behaviour |
|---|---|
| `QAWS_ACCUMULATE_SCATTER` | every sample adds directly into the parameter adjoints |
| `QAWS_ACCUMULATE_TILED` | tiles of `tile_size` samples accumulate into a dense local buffer, flushed once per tile (the CPU analogue of workgroup shared memory) |
| `QAWS_ACCUMULATE_GATHER` | contributions are bucketed per parameter and each parameter sums its own bucket in sample order: deterministic, free of write conflicts |

`qaws_curve_local_support` / `qaws_surface_local_support` return the basis weights of one evaluation (exact for `QAWS_CAP_LINEAR` families; the homogeneous basis for rational ones), and `qaws_curve_build_support_index` / `qaws_surface_build_support_index` build the parameter-to-samples index used by gather kernels.

## Second order

`qaws_curve_eval_batch_tangent2_2d/3d` and `qaws_surface_eval_batch_tangent2` return the second derivative along the combined direction (parameters and coordinates). They are available where `QAWS_CAP_TANGENT2` is set, which includes rational families and knots. Hessian-vector products are provided directly by the inverse-CDF sampling functions (`qaws_curve_cdf_sample_hvp`, `qaws_surface_cdf_sample_hvp`); for other objectives the caller composes them from tangents, second tangents and adjoints. Second order is selective by design: only primitives where it is cheap and exact provide it. `QAWS_CAP_DIRECT_HVP` exists in the capability enum but no object family reports it.

## Families

| Family | Fields with derivatives |
|---|---|
| Bezier, Hermite, polynomial | control points, points and derivatives, coefficients (linear) |
| B-spline curve | control points (linear), knots |
| NURBS curve, rational Bezier | control points, weights (quotient rule), knots (NURBS) |
| Catmull-Rom (uniform, chordal, centripetal), Yuksel C2 | interpolated points (non-linear through the parameterization) |
| Clothoid | origin (`CENTER`), start angle, start curvature and curvature rate; the position integral is differentiated exactly in dual numbers |
| Composite | its segments as children |
| Bezier, B-spline, bilinear, biquadratic surfaces | control points (linear); B-spline also U and V knots |
| NURBS surface | control points, weights, U and V knots |
| Extrusion, ruled, revolution, Coons, loft, Gordon, offset surfaces | their own fields (direction, angle, offset distance, ...) plus their input curves and surfaces as children |

Rational families use homogeneous jets and the quotient recurrence `C^(k) = (A^(k) - sum binom(k,i) W^(i) C^(k-i)) / W`, forward and reverse. `qaws_curve_get_diff_capabilities`, `qaws_curve_get_diff_class` and the surface equivalents report what a given object supports, so callers query instead of hard-coding the table.

### Knots

Knots are differentiable fields of B-spline and NURBS curves and surfaces. A knot direction re-evaluates the basis rows in dual numbers (first and second order), so tangents, second tangents and adjoints include knots. The knot span containing each sample is frozen (`QAWS_FREEZE_SPAN`): moving a knot past a sample is a discrete event. In inverse-CDF sampling knots also move the quadrature spans and the start of the domain, and both are differentiated. Split, join, conversion and degree maps hold the knots of their inputs fixed; the B-spline fit returns its knots as functions of the data.

## Geometry

`qaws_diff_geometry.h` turns jets into differential geometry with derivatives: speed, unit tangent, normal (and binormal), curvature and torsion for curves; unit normal, first and second fundamental forms, Gaussian, mean and principal curvatures for surfaces. Every quantity is a closed-form function of a jet, so each has an eval (value, tangent and tangent2 from primal, tangent and tangent2 jets) and an adjoint (pullback onto the jet), and composes with any evaluation:

```
parameters --eval_tangent--> jet tangent --geometry_eval--> quantity tangent
quantity adjoint --geometry_adjoint--> jet adjoint --eval_adjoint--> parameters
```

Degenerate configurations are reported through the validity: `QAWS_DIFF_ILL_CONDITIONED` when a quantity is undefined (zero curvature leaves the normal, binormal and torsion undefined; umbilics for principal curvatures), `QAWS_DIFF_INVALID` for a singular jet (zero speed, degenerate normal). Undefined quantities and their derivatives are returned as zero.

## Implicit operations

Operations defined by a solve are differentiated with the implicit function theorem, never by differentiating solver iterations: for `g(x, theta) = 0`, `dx/dtheta = -(dg/dx)^-1 dg/dtheta`. This makes the derivative independent of the solver's iteration count and starting point.

| Operation | Condition |
|---|---|
| closest point on a curve / surface | `C' . (C - q) = 0`; `Su . (S - q) = Sv . (S - q) = 0` |
| curve pair (closest approach, intersection) | `A' . (A - B) = 0`, `B' . (A - B) = 0` |
| curve crossing a plane | `n . (C - p) = 0` |
| curve piercing a surface | `S(u, v) - C(t) = 0` |
| extremum along a direction | `e . C'(t) = 0` |
| inflection (xy) | `x'y'' - y'x'' = 0` |

Closest points solve the primal problem themselves. The other operations differentiate one solution: the caller passes a seed from a discrete finder (an intersection search, `qaws_curve_find_extrema`, `qaws_curve_find_inflection_points`); the solution is refined by Newton and differentiated. Which solutions exist is a discrete choice (`QAWS_FREEZE_ACTIVE_SET`). Inputs that are not objects (query point, plane, direction) have their own tangent and adjoint arguments.

Each solve reports through `ctx->report`:

- `condition_number`: conditioning of `dg/dx`;
- `residual`: relative residual at the returned solution;
- `branch_gap`: distance gap to the best competing solution;
- `validity`: `QAWS_DIFF_ILL_CONDITIONED` near singular configurations (tangential intersections, flat extrema), `QAWS_DIFF_AMBIGUOUS` when a competing solution is within 0.5% of the distance, `QAWS_DIFF_VALID_LOCALLY` when a coordinate sits on its domain boundary and is held fixed.

## Differential maps

Operations that build new geometry return, next to their outputs, a `qaws_diff_map` from the parameters of their inputs to the parameters of their outputs. Inputs and outputs are numbered objects, each with its own `qaws_diff_views`:

```c
qaws_curve *left, *right;
qaws_diff_map* map;
qaws_curve_split_diff(&left, &right, &map, curve, 0.4f);
/* inputs: 0 = curve, 1 = split parameter; outputs: 0 = left, 1 = right */
qaws_diff_map_tangent(out_tangents, map, in_tangents, 2, 2, NULL);  /* out = J in */
qaws_diff_map_adjoint(in_adjoints, map, out_adjoints, 2, 2, NULL);  /* in += J^T out */
qaws_diff_map_destroy(map);
```

Available: split, join, Hermite to Bezier, Bezier to B-spline, B-spline to NURBS, degree elevation and reduction, the least-squares B-spline fit and the 3D offset (constant direction, Frenet normal or the normal of a surface at the closest point; the sample set is frozen). Polynomial constructions are exact linear maps; NURBS maps include the homogeneous projection. `qaws_curve_fit_bspline_diff` differentiates the whole fit: chord-length parameters, averaged knots, basis functions and normal equations, with respect to the data points and (when given) the sample parameters. The maps are materialized sparse Jacobians; `qaws_diff_map_get_entries` exposes them.

Why maps: a chain of modelling operations stays differentiable without the caller re-deriving each one. Pulling an adjoint on the output back through the map gives the adjoint on the inputs, which then flows into the input objects' own rules.

## Inverse-CDF sampling

`qaws_diff_sampling.h` samples a curve by inverting a measure `M(t) = int rho |C'|`: arc length (constant speed), curvature-weighted (`rho = sqrt(floor^2 + kappa^2)`) or a user density field in space. The sample at `sigma = distance + fraction * M_total` has a parameter `t` defined by `M(t) = sigma`, differentiated implicitly to first and second order (`qaws_curve_cdf_sample_tangent`), pulled back as an adjoint (one quadrature pass for a whole batch, plus distance adjoints, `qaws_curve_cdf_sample_adjoint`) and as a Hessian-vector product (`qaws_curve_cdf_sample_hvp`). The forward value and every derivative use the same composite Gauss-Legendre rule, so the derivatives are exact for the computed samples.

Time traversal maps onto the same machinery: `qaws_traversal_cdf_targets` turns traversal inputs (easing, motion profile, wrap mode) into arc-length targets with their first and second rates, so traversal samples have exact first and second derivatives in time and in the curve parameters.

Surfaces warp points of the unit square by the marginal and conditional inverse CDFs of `rho(S) |S_u x S_v|` (area, density or curvature measure): stratified or blue-noise points stay stratified on the surface. The forward pass solves the discrete equations in dual numbers (first and second order), the adjoint pulls back through the conditional then the marginal equation, and the HVP differentiates the backward pass along the direction or polarizes the second order forward pass.

## Context and reports

`qaws_diff_context` carries the derivative order, the frozen states the caller accepts, the accumulation strategy, the tile size, a scratch allocator and an optional `qaws_diff_report`. Initialize it with `qaws_diff_context_init`; `ctx = NULL` means first order, scatter accumulation, `malloc` scratch and no report. Reset the report with `qaws_diff_report_reset`.

The report collects, over all samples of a call, the worst `diff_class` and `validity`, the frozen discrete states the result relies on (`frozen_used`), the largest condition number and residual, the smallest branch gap and the index of the worst sample. A result that relies on a frozen state is exact for that state; the report tells the caller when the state matters.

## Threading

All functions are thread-safe on immutable objects. Adjoint functions write into caller storage, so threads sharing one adjoint view (or one report) need external synchronization; giving each thread its own views and summing them is the usual pattern. See [threading](threading.md).

## GPU kernels

The installed portable headers in `src/qaws/core/` compile for C, HLSL, GLSL and Halide (see [backends](backends.md)):

- `qaws_bspline_diff_core.h`: B-spline tangent and adjoint kernels `qaws_bspline_tangent_3d`, `qaws_bspline_adjoint_cp_3d` (the per-element term of a gather pass) and `qaws_bspline_adjoint_t_3d`. `examples/diff_bspline_adjoint.hlsl` runs them as compute shaders; the gather pass loops over the samples listed by `qaws_curve_build_support_index`, so it needs no atomics and its result is deterministic.
- `qaws_bspline_sampling_core.h`: inverse-CDF sampling of B-spline curves on second order dual numbers: `qaws_bspline_integrate` (one span with its first and second rates), `qaws_bspline_cdf_solve` (safeguarded Newton with a fixed iteration count, branch free), `qaws_bspline_cdf_tangent` (first and second sample tangents), `qaws_bspline_cdf_lambda` (the adjoint multiplier) and `qaws_bspline_node_adjoint_cp` (the per-node term of a control point gather). `examples/diff_cdf_sampling.hlsl` chains them as compute passes: spans, prefix sum, samples, multipliers, gather.
- `qaws_bspline_surface_sampling_core.h`: area-measure warping of the unit square onto a B-spline surface, the GPU form of `qaws_surface_cdf_sample_*` with `QAWS_MEASURE_AREA`: `qaws_patch_cell_mass`, `qaws_patch_cdf_solve`, `qaws_patch_cdf_tangent`, `qaws_patch_cdf_multipliers` and `qaws_patch_line_adjoint_cp`. `examples/diff_surface_sampling.hlsl` runs them as compute passes.

Tests `57_diff_core`, `60_diff_sampling_core` and `62_diff_surface_sampling_core` drive the kernels on the C backend and compare them with the runtime.

## Verification

The diff test suites (`49_diff_model` to `62_diff_surface_sampling_core`) run in single and double precision and check:

- tangents against central differences of rebuilt objects;
- adjoints through the identity `<ybar, J xdot> = <J^T ybar, xdot>`, for all three accumulation strategies;
- second tangents against differences of tangents;
- analytic ground truth from Mathematica where a closed derivation exists (`tests/reference/*.wls`, 30 digits): curve and surface jets, differential geometry, differential maps and surface sampling.

## Examples

`examples/diff_showcase.c` (CMake option `QAWS_BUILD_EXAMPLES`, target `qaws_diff_showcase`) writes SVG figures to `shots/` in the working directory:

- 1: B-spline fit to data by adjoint gradient descent
- 2: Sensitivity maps: what a point is affected by (adjoint) and what a control point affects (tangent)
- 3: Curvature fairing through the geometry rules
- 4: NURBS weight optimization: a parabola becomes a circle
- 5: Height-field surface fit (z only, by component mask) and Gaussian curvature sensitivity
- 6: Surface of revolution fitted to scan points through its profile curve (child views)
- 7: Coons patch faired by editing its boundary curves
- 8: Curve networks (Gordon and loft surfaces)
- 9: Projection fit with closest points
- 13: Knot placement: fitting a crease by moving surface knots
- 14: Constant-speed sampling: first and second order sample tangents, and a fit by Newton-CG on exact HVPs
- 15: Inverse-CDF samples under arc length, curvature and a density field, with their tangents
- 16: Stratified points warped onto a patch by the area, density and curvature inverse CDFs
- 17: Blue noise on a patch: repulsion of the warped samples descended through the xi adjoints
- 18: Knot tangents of arc-length samples, and interior knots moved by the knot adjoint

`examples/diff_applications.c` (target `qaws_diff_applications`) applies the API to real data: non-rigid registration of a template curve to a scan, CMA-ES on the pose with gradient refinements inside the fitness. Photos are read as binary PPM; `examples/photo_to_ppm.ps1` converts any image.
