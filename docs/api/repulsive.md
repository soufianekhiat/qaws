# API Reference: qaws_repulsive.h

Tangent-point (repulsive) energies of polygonal curves and triangle
meshes, with exact first and second derivatives. After Yu, Schumacher and
Crane, *Repulsive Curves* (2021), and Yu, Brakensiek, Schumacher and Crane,
*Repulsive Surfaces* (2021). The discrete models, the reconciliation of the
papers with their reference code, and the contracts are in
[`docs/repulsive_plan.md`](../repulsive_plan.md).

Milestones R1 to R6: the energies and their derivatives, evaluated
directly over all pairs or by Barnes-Hut; the constrained fractional
Sobolev solvers for curves (with adapters to qaws curves) and for triangle
meshes, each with a collision-checked mode, dense or hierarchical (no
|V|^2 storage); the surfaces' Poisson-based preconditioner; curves on
parametric surfaces, native control-point optimization, the derivatives
of converged solutions; and GPU kernels for the energy, its gradient and
the metric product.

```c
qaws_repulsive_curves_desc d = { 3, vertex_count, edges, edge_count, 0, 0 };   /* (alpha, beta) = (3, 6) */
qaws_repulsive_curves* c = NULL;
qaws_repulsive_curves_create(&c, &d);
qaws_repulsive_curves_eval(&e, &de, &d2e, c, positions, direction);   /* E, dE[v], d2E[v, v] */
qaws_repulsive_curves_gradient(gradient, &e, c, positions, 1);        /* gradient += seed dE/dx */
qaws_repulsive_curves_hvp(hv, c, positions, direction);               /* hv += H v */
qaws_repulsive_curves_destroy(c);

qaws_repulsive_mesh_desc m = { vertex_count, triangles, triangle_count, 0 };   /* p = 6 */
/* qaws_repulsive_mesh_create / _eval / _gradient / _hvp / _destroy alike */
```

## The energies

- **Curves** (the paper's Eq. 17): for every ordered pair of edges I, J with
  no common vertex, l_I l_J / 4 times the sum over their endpoints of
  |T_I x (x_i - x_j)|^alpha / |x_i - x_j|^beta, T_I the edge's unit
  tangent. The input is a graph (vertices and edges), so it can be one
  curve, several (they repel each other), open or closed, or a network.
  2D or 3D. Edges meeting at a vertex never interact with each other.
- **Meshes** (the paper's Eq. 10): for every ordered pair of distinct faces,
  a_S a_T |N_S . (X_S - X_T)|^p / |X_S - X_T|^(2p), with X the barycenter,
  a the area and N the unit normal. The input is a positions/index view: 3
  coordinates per vertex and 3 indices per face. Connectivity is fixed for
  the object's lifetime.

| | Default | Allowed | Scaling of E under x -> c x |
|---|---|---|---|
| curves | (alpha, beta) = (3, 6) | alpha >= 2, beta in [alpha + 2, 2 alpha + 1] | c^(2 + alpha - beta) |
| meshes | p = 6 | p > 4 | c^(4 - p) |

The curve discretization is the paper's, not the reference code's: TPESC
works on vertices with vertex tangents and dual lengths, which is a
different discrete energy. See the plan, section 2.

## Derivatives

Every derivative is exact for the discrete sum, including how the
tangents, normals, lengths and areas depend on the positions. Each pair
term is evaluated once in a second-order forward jet. No finite
differences are used.

- `_eval`: the value, and with a direction (else NULL) the first and second
  derivatives along it: a JVP and the second directional derivative.
- `_gradient`: `gradient += seed * dE/dx`, a VJP of the scalar energy. This
  is the ordinary coordinate gradient. It is not a descent direction: the
  Sobolev-preconditioned direction is a solver's business (R2).
- `_hvp`: `hv += H v`.

Where a point lies on another edge's line (|T x d| = 0, as on any straight
curve), the energy is twice differentiable for alpha >= 2 and the
derivatives are their limits. For alpha > 2 those are zero; for alpha = 2
they are the square's.

## Contracts

- Positions arrive with every call, so one object serves many
  configurations. Objects are immutable and may be used from several
  threads at once. There is no global state.
- Arithmetic is in double in both builds; f32 builds round the outputs.
- On failure the outputs are untouched: nothing is partially accumulated.
- `INVALID_ARGUMENT`: bad indices, or exponents outside the ranges above.
- `DEGENERATE_CURVE`: a zero-length edge or a zero-area face, which has no
  tangent or normal.
- `NUMERICAL_FAILURE`: two interacting points coincide. That is a
  self-intersection with infinite energy. It is never softened with an
  epsilon into a finite value.
- Cost: O(n^2) pairs, direct. "Direct" means the discrete sum without
  approximation, not the continuum integral. Hierarchical acceleration is
  R4.
- The energy's weight is not a clearance: it is not a strand radius or a
  shell thickness.

## Checks (test 112)

| Check | f64 | f32 |
|---|---|---|
| Against Mathematica (the same discrete sums differentiated symbolically, 30 digits; `tests/reference/112_repulsive.wls`): energy, dE[v], d2E[v, v], gradient and Hv for a closed 3D curve, an open 2D curve (alpha = 2), a network (a junction, a triangle, an open polyline), an octahedron and three patches (p = 5.5) | 1.1e-14 (test bound 1e-12) | 1.4e-6 (bound 1e-5) |
| Translation and rotation invariance; the scaling degree for E and the gradient | 2.3e-15 | 1.6e-6 |
| Adjoint identity, Hessian symmetry | 2e-14 | 5.3e-7 |
| Hv against central differences of gradients (best of h = 1e-3 to 1e-7) | 1.8e-10 | |
| A straight curve: E, gradient and Hv exactly 0 (alpha = 3); Hv against differences (alpha = 2) | | |
| Refusals, with outputs untouched | | |

## Curve optimization (R2)

```c
qaws_repulsive_curves_solve_desc sd = { 0 };
sd.hold_edge_lengths = 1;          /* or hold_total_length; pins, held edge directions */
sd.collision_checked = 1;          /* validate every accepted motion */
qaws_repulsive_curves_solve(&report, positions, c, &sd);              /* in place */

qaws_repulsive_curves_solver* s = NULL;                               /* or step by step */
qaws_repulsive_curves_solver_create(&s, c, positions, &sd);
qaws_repulsive_curves_solver_step(&report, positions, s);             /* report.reason: RUNNING or why it stopped */
qaws_repulsive_curves_solver_direction(g, s, positions);              /* the descent direction, not the gradient */
qaws_repulsive_curves_solver_destroy(s);
```

The method is the paper's (Sec. 4.2 and 5.2 to 5.4), direct and dense.

- **Metric.** A = B + B0 acts on vertex functions:
  - The high-order term acts on edge differences D_I u = (u_i2 - u_i1) T_I / l_I, weighted
    w_IJ = l_I l_J / 4 sum 1 / |x_i - x_j|^(2 sigma + 1), with
    sigma = (beta - 1) / alpha - 1.
  - The low-order term acts on edge averages, weighted by
    k_4^2 / |x_i - x_j|^(2 sigma + 1).
  - Both are assembled over the same pairs as the energy: edges with no
    common vertex.
- **Direction.** The solver solves [A~ C^T; C 0] [g; l] = [dE; 0], where A~
  is the metric on each coordinate and C is the constraints' Jacobian. It
  is solved by blocks (not CG), which A~ = A on each coordinate allows:
  pinned vertices are eliminated, A's free block is inverted (|V|^3), then
  Y = A~^-1 C^T and the Schur complement (C Y) l = C A~^-1 dE are solved
  by LU (k^3). That is |V|^3 + k^3 rather than (dim |V| + k)^3, which
  matters when every edge's length is held (k ~ |V|): the paper's
  600-vertex scenes take about 1 s a step, not 10. With nothing pinned, A
  (whose kernel is the constants) is lifted by tau 1 1^T / |V| first; the
  barycenter rows keep g's mean 0, so the solution is unchanged. The
  result g is the descent direction in the metric, and
  `_direction` returns it. It is not the coordinate gradient
  (`qaws_repulsive_curves_gradient`); neither stands in for the other.
- **Step.** Backtracking on the largest vertex move h:
  - h starts at twice the last accepted move (a quarter of the mean edge
    length at first);
  - the trial point is y = x - (h / |g|_inf) g;
  - y is projected back onto the constraints in the same metric, along
    Y: dx = -Y m. First chord steps with the direction's factors (a solve
    each), then Newton steps (C(y) Y) m = phi(y) if those stop
    contracting. It aims at 1e-3 of the tolerance: a residual near the
    tolerance moves the energy more than the last steps lower it;
  - the step is accepted when Armijo holds (sufficient decrease 0.1 of the
    linear prediction: with less, a step nearly twice the best one passes,
    and the iterates zigzag across a narrow valley without converging) on the projected y and, in
    collision-checked mode, the motion x -> y passes the check.
- **Constraints.** Each is held at its value when the solver is created:
  - pinned vertices (an obstacle is a pinned component);
  - held edge directions: the products of the edge with the direction's
    normals stay 0 (one row in 2D, two in 3D);
  - the total length, or every edge's length;
  - the vertex barycenter, held by default when nothing is pinned, since
    the metric does not see translations.

  Rows that would repeat pins are dropped: the length or direction of an
  edge whose ends are both pinned. Total and per-edge lengths together are
  refused. Without barycenter or pins the system is singular, and the step
  reports NUMERICAL_FAILURE rather than a silent result.
- **Stopping.** The stationarity is <dE, g> l / (E |g|_inf): the energy's
  relative first-order fall for a move of one mean edge length l along g.
  The solve is CONVERGED when it drops under the tolerance (1e-4).

| Termination | Meaning |
|---|---|
| CONVERGED | stationarity under the tolerance |
| ITERATION_LIMIT | the budget spent (solve) |
| INFEASIBLE | no trial step could be projected back onto the constraints |
| COLLISION_LIMITED | every trial step collided, or the accepted ones vanish (under 1e-8 l) with collisions rejected |
| STALLED | no trial step lowered the energy, or the accepted ones vanish (under 1e-12 l) |
| NUMERICAL_FAILURE | a singular system or an unusable energy |
| CANCELLED | the host's cancel callback stopped it (status QAWS_STATUS_CANCELLED) |

The report also counts the trial steps rejected by the energy, by the
collision check and by the projection. It gives the constraint residual and
the size of the last step.

### Host control (test 145)

```c
sd.progress.cancel = my_cancel;      /* int (*)(void* user): non-zero stops */
sd.progress.on_step = my_on_step;    /* void (*)(report, positions, user), after each accepted step of _solve */
sd.progress.user_data = &my_state;
qaws_repulsive_curves_solver_objective(s, positions, &E, gradient);   /* what the solver lowers, and its gradient */
```

The same `progress` member is in the mesh desc. It applies to the curves
on a surface and the control-point solvers through their curves desc.

- **Cancellation.** The callback is called on the caller's thread, never
  from the executor's. It is polled at the start of a step, between the
  step's stages, and before each line-search trial's projection and
  energy. The stages are the energy and gradient, the metric's assembly,
  and its factorization (or the hierarchical direction's solves). A step
  therefore stops within one stage; a dense factorization is the longest.
  On cancellation the step returns QAWS_STATUS_CANCELLED with the reason
  CANCELLED. The positions are untouched: a step writes them only when it
  accepts. `_solve` stops with the positions of the last accepted step. A
  time budget is a cancel callback that reads a clock. The solver keeps
  stepping when called again.
- **Progress.** `_solve` calls `on_step` after each accepted step. It
  passes the report (energy, stationarity, residual, step size) and the
  caller's positions array, read-only during the call.
- **Objective.** `qaws_repulsive_*_solver_objective` evaluates the
  function the solver lowers: the energy plus every term the desc adds
  (length, room; area, volume, obstacles for meshes), as `report.energy`
  shows it. It also returns its gradient with respect to the solver's
  variables (positions, uv, control points).

### Line search and GMRES settings (test 145)

```c
sd.armijo = 0.2;                     /* sufficient decrease, in (0, 0.5); 0: 0.1 */
sd.backtrack = 0.25;                 /* step factor per refused trial, in (0, 1); 0: 0.5 */
sd.max_trials = 10;                  /* trials per step; 0: 40 */
sd.gmres_max_iterations = 500;       /* 0: 1000 (hierarchical), 2400 (derivatives) */
sd.gmres_restart = 30;               /* 0: 40 (hierarchical), 80 (derivatives) */
sd.gmres_tolerance = 1e-6;           /* relative residual; 0: 1e-8 (hierarchical), 1e-11 (derivatives) */
```

The same members are in the mesh desc, with the same defaults; the curves
on a surface and the control-point solvers take them from their curves
desc. Zero keeps the previous behaviour bit for bit. A trial is accepted
when E(y) <= E - armijo tau <dE, g>; each refused trial shrinks the step by
`backtrack`; after `max_trials` refused trials the step stops (STALLED,
COLLISION_LIMITED or INFEASIBLE, by the last refusal). The GMRES limits
apply to the hierarchical mode's scalar solves and to the solution
derivatives' KKT solve, each with its own defaults; a value set applies to
both. A solve whose restarts stop gaining is still accepted under 1e-6
(hierarchical) or 1e-7 (derivatives). Out-of-range values (armijo outside
[0, 0.5), backtrack outside [0, 1), a negative or NaN tolerance) make the
solver's creation return INVALID_ARGUMENT.

| Check (test 145) | Result |
|---|---|
| The defaults given explicitly against zeros: 8 steps, dense and hierarchical | positions and reports bit for bit |
| Hierarchical, `gmres_max_iterations` 1 | the step's reason NUMERICAL_FAILURE, nothing moved |
| Hierarchical, `gmres_tolerance` 1e-3 | the energy still falls (2.96796 against 2.96800 by default after 8 steps), other positions |
| A mesh step from a first trial 40 mean edge lengths long | backtrack 0.5: 7 refused trials, then accepted; 0.1: 2; `max_trials` 2: STALLED after 2, positions untouched |
| armijo 0.5, -0.1, 0.7; backtrack 1, -0.5; `gmres_tolerance` -1 | INVALID_ARGUMENT from the curve and the mesh solver |

### Allocators (test 145)

```c
qaws_repulsive_curves_create_ex(&cd, &allocator, &curves);
qaws_repulsive_curves_solver_create_ex(curves, &sd, positions, &allocator, &solver);
```

Every creator has an `_ex` variant (curves, meshes, control points, and the
four solvers). The object keeps a copy of the allocator. Its blocks come
from it, and so do the blocks of every call on it: eval, gradient, hvp, bh,
the curves-on-surface evaluations, a solver's steps, direction, objective
and derivatives, and `*_solve`. Each call installs the allocator as the
thread allocator on entry and restores the previous one on every return.
Per-call temporaries still come from the thread's scratch allocator when
one is installed (`qaws_exchange_thread_scratch_allocator`). A NULL
allocator, and the non-`_ex` creators, keep the usual order. A solver made
without `_ex` follows the usual order whatever its curves' or mesh's
allocator; `*_solve` uses the allocator of the object it is given. With an
executor, the temporaries of jobs on other threads follow those threads'
order.

| Check (test 145) | Result |
|---|---|
| Curves and a mesh made with `_ex` and a counting allocator; the curves' eval, gradient, hvp, bh; a dense collision-checked and a hierarchical solver (`_ex`) stepping, their direction and objective; `_solve`; a mesh solver converged, its JVP (checked) and VJP; all destroyed. A second counting allocator installed as the thread allocator | 1756 blocks from the objects' allocator, none left after destroy, none from the thread allocator |
| The same with a scratch arena installed | no call leaves a block in the arena (`live` 0 after each), peak 23 KB |

### The two modes

- **Energy-based.** Every accepted step lowers the energy and holds the
  constraints, and nothing more. The discrete energy samples only the
  vertices, so on a coarse polygon a strand can pass through another
  between vertices. A 48-vertex trefoil with its total length held ends as
  a circle (figure 32). At 96 vertices the same trefoil stays knotted, but
  that is the energy's doing, not a guarantee.
- **Collision-checked.** Every accepted motion is validated: the straight
  path from the old positions to the projected new ones, so the
  projection's own movement is included.
  - **Check:** conservative advancement on every pair of edges sharing no
    vertex. A pair's distance is Lipschitz in time, with the sum of its
    largest endpoint moves as the constant. Each pair must stay at least
    min(clearance, half its distance at the step's start) apart; the
    clearance is 1e-3 l by default.
  - **Guarantee:** no two such edges ever touch, so the topology is kept: a
    knot keeps its type and a link its linking number.
  - **Input:** if edges already touch, it is refused with INVALID_ARGUMENT.
    There is no repair and no promise of untangling.
  - **Exemption:** edges that share a vertex are not checked against each
    other, since they meet there by construction. At a junction that means
    every pair of its edges, and only those.

### Adapters

`qaws_repulsive_polyline_from_curves` turns qaws curves (2D or 3D) into the
working polyline.
- Each curve's parameter range is split until every edge is at most
  edge_length long and its middle lies within tolerance of its chord.
- Each vertex keeps its source curve and parameter.
- A curve whose ends meet becomes a loop.

`qaws_repulsive_polyline_to_curve` refits one curve's vertices as a
degree-1 B-spline with those parameters as knots. The result is exactly the
working polyline, parametrized as the source, and nothing is smoothed, so a
collision-checked result stays valid for the refitted curve. Smooth refits,
and optimizing a curve's own control points, are R5.

### Checks (test 113)

| Check | f64 | f32 |
|---|---|---|
| The metric against Mathematica (`tests/reference/113_repulsive_metric.wls`: from the paper's bilinear form, not its increment rules), a closed 3D heptagon and an open 2D curve | 2.2e-16 | 8e-8 |
| The metric symmetric, constants in its kernel, positive semidefinite | exact, 3.4e-16 | |
| The direction tangent to the held edge lengths, descending | 1.3e-16 | 5.8e-9 |
| A trefoil (48 vertices), edge lengths held: converged, the energy falling at every step | 17 steps, lengths within 9.6e-9 | within 1e-5 |
| An open 2D curve, ends pinned, end tangents and total length held, collision-checked | pins 1e-21, tangents 2.8e-15, length 7.8e-15 | 4.5e-6 |
| A Hopf link, edge lengths held, collision-checked: linking number (an independent Gauss integral) | -1.0000 before and after, 16 steps | same |
| The path check: a sweep through another segment with both ends clear is caught; sliding past passes; closing under the clearance is caught | | |
| A circle and a cubic -> working polyline (provenance exact, the circle a loop, edges at most 0.15) -> collision-checked solve -> refit through every vertex at its parameter | 0 / 1e-14 | 1.7e-7 |
| An obstacle loop pinned with its edge lengths held: solves, the loop unmoved | | |
| Refusals: touching input (collision-checked), total and edge lengths together; free translations reported as NUMERICAL_FAILURE; the iteration budget | | |

Cost: each step assembles the dense metric (O(E^2)) and factors a dense
(dim V + k)^2 system, so the solver is meant for hundreds of vertices; a
48-vertex trefoil takes about 25 ms a step. Hierarchical evaluation and
iterative solves are R4.

Figure 32 (`examples/better2d_showcase/32_repulsive_solve.c`) shows the two
modes on a trefoil, pinned strands round a fixed loop, and a Hopf link.

## Mesh optimization (R3)

```c
qaws_repulsive_mesh_solve_desc sd = { 0 };
sd.hold_area = 1; sd.hold_volume = 1;      /* volume: closed meshes only */
sd.collision_checked = 1;                  /* vertex-face and edge-edge checks on every accepted motion */
qaws_repulsive_mesh_solve(&report, positions, mesh, &sd);
/* or qaws_repulsive_mesh_solver_create / _step / _direction / _destroy, as for curves */
```

The surfaces' fractional Sobolev descent follows Repulsive Surfaces,
Sec. 3.2 and 3.3. It runs on the same core as the curves: the constrained
direction, the backtracking step, the Newton projection, the two modes and
the report.

- **Metric.** A = B + B0 acts on vertex functions, with s = 2 - 2/p and
  sigma = s - 1. It is summed over every ordered pair of distinct faces
  with weight a_S a_T / |X_S - X_T|^(2 sigma + 2):
  - the high-order term (Eq. 12) acts on the per-face gradients of the hat
    functions, D u(S) = N_S x (u_i e_jk + u_j e_ki + u_k e_ij) / (2 a_S);
  - the low-order term (Eq. 13) acts on face averages, weighted by the
    p = 2 kernel (N_S . d)^2 / |d|^4.
- **Constraints.** Each is held at its value when the solver is created:
  - pinned vertices (an obstacle is a pinned component);
  - the total area;
  - the enclosed volume (1/6) sum x_i . (x_j x x_k), for closed meshes only
    (every edge on exactly two faces, oppositely oriented), refused
    otherwise;
  - the barycenter, by default when nothing is pinned.
- **Area weight.** `area_weight` w > 0 adds w times the area to the
  energy: the paper's "minimize area" (Fig. 9, centre). The descent then
  runs without the energy's Hessian-vector product, and the solution's
  derivatives are refused (unsupported). A negative weight is refused.
- **Initial step.** `initial_step` sets the first trial's largest vertex
  move, over the mean edge length (0: 0.25). A solver made again on a
  remeshed surface resumes with twice its predecessor's `last_step`, where
  that solver would have tried next. A negative one is refused.
- **Translations.** The fractional metric barely resists a component
  moving as a whole, so its solve can let one rigid motion dominate the
  direction. With `translations` 2, each connected component without a pin
  is solved with its barycenter held, then moved by its dual-area mean of
  the differential (Repulsive Surfaces, Fig. 8; the paper's published
  program declares this constraint but never uses it). The overall
  barycenter is still held when `hold_barycenter` asks. With 1 the whole
  mesh does so and its barycenter drifts (the program's
  `allow_barycenter_shift`; refused with pins). Each group takes 3 rows,
  so the hierarchical mode (32 rows) takes about 10 components. No
  solution derivatives with either.
- **Moving targets.** `qaws_repulsive_mesh_solver_set_targets` changes the
  held area and volume between steps, as the paper's program's growth
  schedules do (`constrain area 8 200`: the target raised a little each
  step until 8 times the input's after 200). The next projections head
  there, so the target should move by little at a time.
- **Volume deviation.** `volume_weight` w > 0 adds w ((V - Vt) / Vt)^2 to
  the energy, Vt = `volume_target` (0: the input's volume), the program's
  `minimize volume_deviation`. A closed mesh; no solution derivatives.
- **Obstacles.** `obstacles` (`qaws_repulsive_obstacle`: positions,
  triangles or none, a weight) are fixed shapes the surface is kept from,
  the program's `obstacle` and `point_cloud_obstacle`. The weight times
  their tangent-point energy with the surface is added: with triangles,
  both orders of each pair of faces; a point cloud, each point weighing 1
  against the surface's faces and normals. Each obstacle's tree is built
  once; a surface face visits it as in the Barnes-Hut energy, a far node
  by its forward and reverse terms (its averaged projector), a near leaf
  item by item (then exactly). They add nothing to the metric and have no
  solution derivatives. `collision_checked` checks them as it checks the
  mesh against itself: a mesh vertex against an obstacle's faces, a mesh
  face against its vertices or points, edges against edges, each motion
  keeping at least min(clearance, half the distance at its start). Their
  energy alone repels, but a step lands where it lands: unchecked, a shrink
  wrap passed between a sparse cloud's points, and an envelope through the
  level inside it.
- **Step cap.** `max_step` caps a step's largest vertex move, in mean edge
  lengths (0: none; the trial size otherwise doubles after each accepted
  step).

  The area and volume rows are divided by the mean edge length l and l^2,
  so that one tolerance (in lengths) serves every row.
- **Collision-checked.** Conservative advancement checks every vertex
  against every face it is not a corner of, and every two edges that share
  no vertex. The rule is the curves': each pair stays at least
  min(clearance, half its start distance) apart. Contact between triangle
  meshes always begins as such a pair, so no two triangles touch beyond the
  vertices and edges they share. The input is also checked for an edge
  passing through a face it shares no vertex with, since motion from a
  clear state cannot reach that without a vertex-face or edge-edge contact
  first. Touching input is refused.

### The preconditioner (Sec. 5.2)

The paper solves the metric iteratively, by GMRES, preconditioned by
(-L)^+ L^(2 - s) (-L)^+:
- L is the cotan Laplacian, so (-L)^+ is its pseudo-inverse on mean-zero
  functions (a bordered dense LU here);
- L^(2 - s) is the fractional operator of order 2 (2 - s), so the
  composition has order -2s, A's inverse order.

This is built and validated here
(`qaws_internal_repulsive_mesh_gmres`, restarted GMRES with right
preconditioning). The dense mode solves directly (A inverted), which is
exact at these sizes; the hierarchical mode (R4, below) uses the preconditioner,
with the Laplacian factored sparsely.

One deliberate change from the paper concerns L^(2 - s), which is
discretized by vertex quadrature: the sum over pairs of distinct vertices
of (u_i - u_j)(v_i - v_j) m_i m_j / r_ij^(2(2 - s) + 2), with m the lumped
areas. The paper's form (Eq. 11) works on face averages. On a
3-colourable mesh (every subdivided octahedron or icosahedron), a vertex
function that differs between two colour classes averages to zero on
every face. That form then has a kernel beyond the constants, and the
preconditioned system is singular there. On a lumpy 258-vertex sphere it
had two extra zero eigenvalues, and GMRES stalled at 3000 iterations. The
vertex form's kernel is the constants only.

### Checks (test 114)

| Check | f64 | f32 |
|---|---|---|
| The metric against Mathematica (`tests/reference/114_repulsive_mesh_metric.wls`: from the paper's bilinear forms, Eqs. 12 and 13), an open patch (p 6) and an octahedron (p 5.5) | 4.8e-16 | 1e-8 |
| The metric symmetric, constants in its kernel, positive semidefinite | exact, 3.4e-16 | |
| GMRES to 1e-10 on a lumpy sphere (258 vertices): unpreconditioned against (-L)^+ L^(2-s) (-L)^+ | 510 against 13 iterations, solutions 9.8e-10 apart | same |
| The preconditioned spectrum (Mathematica, on the same sphere): condition of P A, excluding the constants | 2.43 (A alone 2750, (-L)^+ A 36) | |
| A bumpy sphere (128 faces), area and volume held | E 16.43 -> 5.41 in 25 steps, area 8.8e-14, volume 2.1e-13 | 1.6e-9, 4.2e-9 |
| A sheet under a pinned one, area held, collision-checked: the pinned sheet unmoved | 0 | 0 |
| The path check: a vertex through a face and an edge through an edge, both ends clear, caught | | |
| Refusals: volume on an open mesh; an edge through a face in the collision-checked mode; free translations reported | | |
| The area weight (2): the energy at the input against the unweighted energy + 2 A; 20 steps, the area 12.59 -> 4.16 (unweighted it inflates to 861); solution derivatives and a negative weight refused | within 1e-10 relative | within 1e-4 |
| The initial step: the first move within a cap of 0.01 mean edges, longer with a cap of 1; a negative one refused | | |
| Translations (2), two spheres: each one's mean direction against its dual-area mean of the differential, less their mean; 10 steps part them with the overall barycenter kept; a pinned sphere stays while the other moves; the drift (1) refused with pins | within 1e-10 relative; the barycenter 1e-12 | 1e-4, 1e-5 |
| Moving targets: the held area raised 5 % a step for 10 steps reaches 1.5 times its input; a volume not held and a negative area refused. The volume's deviation (3, toward twice the volume): the energy at the input against the unweighted + 3/4, the volume growing; a negative weight refused | 1e-8; 1e-10 relative | 1e-5; 1e-4 |
| Obstacles: a sphere's mesh obstacle (weight 2, theta 0) against twice the two meshes' cross terms (the union's energy less each one's), and its gradient; hierarchical (0.25) within 5 %; a point cloud (weight 3) against the brute-force sum; drifting free (translations 1), the sphere moves away from it; a weight of 0 and an index out of range refused | 1e-10 relative, the gradient 1e-9 | same (in double) |
| An obstacle collision-checked: a sphere's volume raised 1.6x toward a sheet just above it, weighted 1e-6 (too light to hold it by its energy): grown past 1.4x, its top stays under the sheet | | |
| The step cap (0.005): no move over 5 steps beyond it, whatever the trial size grew to | | |

Cost: the metric and the systems are dense, O(F^2) to assemble and
(3V + k)^3 to factor. A 128-face sphere takes well under a second per
solve.

Figure 33 (`examples/better2d_showcase/33_repulsive_mesh.c`) shows a bumpy
sphere before and after, with its area and volume held, and a membrane
with its rim pinned, pushed away from a fixed sheet in the
collision-checked mode.

## Hierarchical evaluation (R4): Barnes-Hut

```c
qaws_repulsive_bh_stats stats;
qaws_repulsive_curves_bh(gradient, &energy, &stats, curves, positions, theta, 1);   /* theta 0: the direct sum */
qaws_repulsive_mesh_bh(gradient, &energy, &stats, mesh, positions, theta, 1);
```

Every ordered pair's term uses its first element's own tangent or normal.
For a target element I, a whole tree node of far elements is therefore
taken at once, through its total length or area W and its weighted
centroid X.

- **Tree.** A bounding volume hierarchy over the elements' boxes, split
  at the median centre of the widest axis (quickselect, no globals) down
  to 8 per leaf.
- **Far.** A node is far for I when
  max(r_I, r_node) <= theta dist(box_I, box_node), with r the boxes' half
  diagonals and dist the gap between the boxes. A node holding an element
  that shares a vertex with I is at distance 0, so it is never far: the
  adjacency rules hold under approximation.
- **Energy.** Near pairs are the direct evaluator's exact terms, in closed
  form (value and first derivatives; no jets on this path). A far
  node contributes w_I W k(x_I, X, T_I). The I side keeps the direct
  term's quadrature (the trapezoid over an edge's two endpoints, a face's
  barycentre), so the error comes from the node's size alone.
- **Gradient.** This follows the paper's scheme (Surfaces, Sec. 4.2):
  - each element's vertices get, from its own traversal, the derivative of
    both orders of each pair, e(I, J) + e(J, I), taken in I's vertices;
  - this derivative is exact for the near pairs;
  - for a far node, the reverse term replaces each member's tangent or
    normal by the node's averaged projector (sum w (I - T T^T), or
    sum w N N^T, over W).

  It approximates the direct gradient and converges to it as theta falls;
  it is not the derivative of the approximated energy. (The exact
  derivative of the approximation was tried first. It turns each far pair
  into a monopole force on the whole node, and on a 300-edge trefoil at
  theta 0.25 it came out 14% off the direct gradient. This scheme gives
  6%.)
- **Not provided:** hierarchical HVPs.

| Check (test 115) | f64 | f32 |
|---|---|---|
| theta 0 against the direct evaluator: trefoil (300 edges), lumpy sphere (512 faces) | energy 3.5e-16 and 1.6e-14, gradient 2.7e-13 and 3.2e-15 | 1e-6 |
| Trefoil, theta 0.5 / 0.25 / 0.1: energy off | 3.1e-3 / 5.2e-4 / 3.8e-5 | |
| Trefoil, theta 0.5 / 0.25 / 0.1: gradient off (max norm) | 0.16 / 0.061 / 0.0016 | |
| Lumpy sphere, theta 0.5 / 0.25 / 0.1: energy off | 0.037 / 0.0015 / 1.7e-14 | |
| Lumpy sphere, theta 0.5 / 0.25 / 0.1: gradient off | 0.0026 / 1.3e-4 / 2.2e-15 | |
| A sphere of 2048 faces, theta 0.25, energy and gradient | 0.81 s against 5.5 s direct, energy within 1.6e-2 (times reported, not asserted) | |

The curves' gradient converges more slowly than their energy. The
kernel's derivative falls like 1/d^7, so a monopole far field at theta
0.25 is coarse there; theta 0.1 is within 0.2%. Figure 34 plots the errors
against theta and the time against size: a 2000-edge trefoil takes 0.70 s at
theta 0.25 against 12 to 21 s direct (the times vary with the machine's load).

## Hierarchical solve (R4)

```c
qaws_repulsive_curves_solve_desc sd = { 0 };
sd.hold_total_length = 1;
sd.theta = 0.25;                           /* 0: the dense mode above */
qaws_repulsive_curves_solve(&report, positions, curves, &sd);
/* meshes alike: qaws_repulsive_mesh_solve_desc.theta */
```

With `theta > 0` the solvers store nothing of size |V|^2:

- **Energy and gradient.** These are the Barnes-Hut ones above.
- **The metric.** A u is applied from Barnes-Hut interaction lists, built
  once per step. A near pair carries its exact weights, the same ones as
  the dense assembly. A far node carries one weight per unit length or
  area, together with its sums of w g and w u. The low-order term's
  reverse kernel uses the node's averaged projector, as the gradient does.
  At theta 0 the product is the dense metric's.
- **The direction.** [A~ C^T; C 0] [g; l] = [dE; 0] is solved with:
  - pinned vertices eliminated: their rows and columns are dropped, and g
    is 0 there;
  - the remaining rows (barycenter, lengths, held directions, area,
    volume; at most 32) solved by their Schur complement;
  - z = A~^-1 dE and Y = A~^-1 C^T by flexible GMRES (restart 40,
    relative tolerance 1e-8), one scalar solve per coordinate that a
    right-hand side has;
  - (C Y) l = C z, then g = z - Y l.

  With nothing pinned, A's kernel (the constants) is lifted by
  A + tau 1 1^T / |V|. The barycenter rows keep g's mean at 0, so the
  solution is the same.
- **The preconditioner** is the paper's (-L)^-1 L^(2 - s) (-L)^-1. The
  Laplacian uses 1/l weights for curves and cotan weights for meshes. It
  is shifted by a small multiple of the lumped mass, so that every
  component is invertible, and solved by a profile Cholesky in a reverse
  Cuthill-McKee order, factored once per step. L^(2 - s) is the
  fractional operator by vertex quadrature, applied hierarchically too.
  The constants on the free vertices are passed through by 1/tau.
- **The projection** after a step is Newton's in the same metric. It
  reuses Y: dx = -Y m with (C(y) Y) m = Phi(y).
- **The collision-checked mode** is the dense one's: the same pairs, the
  same rule, through a broad phase (below).
- **Refused here:** per-edge lengths (one row per edge is too many for the
  Schur complement).

| Check (test 116) | f64 | f32 |
|---|---|---|
| The metric from the lists at theta 0 against the dense metric: two circles in 2D, a trefoil (200 edges each), a lumpy sphere (512 faces) | 3.9e-15, 4.3e-15, 2.1e-15 | same |
| The same at theta 0.5 / 0.25 / 0.1 | curves 0.0048 to 0.0069 / 4.4e-4 to 1.2e-3 / 2.4e-5; the sphere 1.1e-4 / 6.9e-6 / 2e-15 | |
| The iterative direction with every pair near, against the dense direction. Curves: barycenter; two pins and the length; 2D with a held direction. Meshes: area and volume; pins and area | 5.6e-10 to 4.4e-8 | 9e-9 to 4.6e-8 |
| The same at theta 0.25: 1 - cos against the dense direction | 2.3e-6 at most | same |
| A 600-edge trefoil, length held, theta 0.25 | E 40.9 -> 34.4, converged in 12 steps, length within 3.2e-13 | 9.5e-13 |
| A 2048-face sphere, area and volume held, theta 0.25 | E 3.70 -> 3.17 in 8 steps, held within 3.9e-12 | 2.9e-6 |

Time per step (f64, one core, the profile above; the times vary with the
machine's load):

| Size | Dense | Hierarchical (theta 0.25) |
|---|---|---|
| Curves, 250 / 500 / 1000 edges | 0.37 / 2.25 / 15.8 s | 0.14 / 0.47 / 0.79 s |
| Curves, 2000 / 4000 edges | | 1.96 / 5.4 s |
| Meshes, 512 / 2048 faces | 0.91 / 32.7 s | 0.56 / 4.3 s |
| Meshes, 8192 faces | | 24.7 s |

On meshes the step is dominated by the metric products. At 2048 faces,
theta 0.25 keeps most pairs near. On curves, GMRES takes about 45
iterations per solve; on meshes, about 11.

Figure 35 shows a 600-edge (2, 5) torus knot relaxed in this mode, and the
time per step against size.

## Curves on a surface (R5)

```c
qaws_repulsive_surface_desc sd = { surface, 1, 1 };     /* periodic in u and v: declared */
qaws_repulsive_curves_solve_desc d = { 0 };
d.hold_edge_lengths = 1;
qaws_repulsive_surface_curves_solve(&report, uv, curves, &sd, &d);   /* uv: two per vertex */
qaws_repulsive_surface_positions(positions, uv, vertex_count, &sd);
/* and _eval / _gradient / _hvp in uv */
```

The unknowns are the vertices' surface parameters, with x_i = S(u_i, v_i).
The energy is the curves' own, measured in the embedding: no UV or
geodesic distance enters the published energy.

- **Derivatives in uv.** The gradient is J^T g. The HVP is
  J^T H J v + sum g . S''[v, .], and the second directional derivative
  d2E[Jv, Jv] + g . S''[v, v]. S's first and second derivatives come from
  qaws_surface_evaluate.
- **The descent.** It runs on the shared core with the metric pulled back,
  the blocks A_ij J_i^T J_j (an N x N metric, N = 2|V|). A 1e-10 multiple
  of the diagonal blocks lifts what A cannot see: a curve sliding along a
  cylinder's axis is a symmetry of the energy. Steps move the parameters,
  so the curves stay on S exactly and no projection onto the surface is
  needed.
- **Constraints.**
  - Pins hold a vertex's parameters.
  - Lengths and held directions are the curves' own, in 3D: the rows of an
    internal x-space curve solver times J.
  - No barycenter is held.
- **The domain** is explicit:
  - a declared periodic parameter wraps, so seams are crossed freely, and
    uv stays unwrapped (continuous);
  - a parameter that is not periodic stays in its range;
  - a trimmed surface's vertices stay in its region;
  - a vertex at a singular point (S_u x S_v about 0: a pole) is refused.
  - A step that would leave the domain is rejected like one that raises
    the energy. Input outside the domain is refused (OUT_OF_RANGE).
- **Not covered.** Only vertices are kept in the domain: an edge's chord
  is straight in 3D. One surface is one chart: crossing between the
  patches of a multi-patch model is not supported. The dense mode only.
- **Collision-checked.** The straight motion of each vertex from x to x',
  a chord of its path on S, is checked with the curves' rule.

| Check (test 117) | f64 | f32 |
|---|---|---|
| Against Mathematica (`tests/reference/117_repulsive_surface.wls`): E(S(uv)) for a closed hexagon, and an open polyline with a segment, on a bicubic Bezier patch. E, dE and d2E along a direction are composed in one variable; the gradient and Hv come from the symbolic x-space derivatives and S's exact ones | 1.3e-14 | 3.2e-6 |
| Torus, two loops across the u seam: E with uv shifted by a period; the uv gradient against central differences | 0; 1.2e-9 | 0; 1.1e-3 |
| Torus, three loops, length held: E 534 -> 40.6 in 30 steps, length 2.2e-14, off the torus | 1.1e-16 | 7.5e-8 |
| Torus, collision-checked, two pins: E falls, the pins' parameters held | 2.6e-24 | 3e-21 |
| Patch (not periodic), an open polyline growing: the energy falls, every parameter stays in [0, 1] | yes | yes |
| Refusals: outside the range, in a trimmed hole, at a sphere's pole | OUT_OF_RANGE | OUT_OF_RANGE |

## Native control points (R5)

```c
qaws_repulsive_controls_desc cd = { 0 };
cd.polyline.curves = curves; cd.polyline.curve_count = n; cd.polyline.dimension = 3; cd.polyline.edge_length = h;
qaws_repulsive_controls_create(&ctl, &cd);
qaws_repulsive_controls_read(P, ctl);                      /* every curve's control points */
qaws_repulsive_controls_solve(&report, P, ctl, &solve_desc);
qaws_repulsive_controls_to_curve(&curve_i, ctl, P, i);    /* the same family, new control points */
```

The working polyline of the R2 adapter is a linear function of the curves'
control points: x_s = sum_k B_sk P_k, the same weights for every
coordinate. Its samples (curve, parameter) are chosen at creation and then
held, so the derivatives are those of that fixed sampling.

- **B** comes from the curves' own tangent evaluation (qaws_diff.h), one
  control point at a time, and is checked twice:
  - a control's tangent in one coordinate moves that coordinate only;
  - B P reproduces the samples.

  A curve that fails is refused (UNSUPPORTED_OPERATION). Bezier, B-spline
  and NURBS (weights held) qualify.
- **Derivatives in P** are B^T those in x, exactly, since the map is
  linear; the HVP is B^T H B v.
- **The solver** descends in P with the metric B^T A B per coordinate.
  With theta > 0, the energy is Barnes-Hut's and B^T A B is assembled
  column by column from the hierarchical metric (|P|^2 stored, never
  |V|^2).
  - Pins name control points.
  - The barycenter is the samples'.
  - Lengths and held directions are the polyline's.
- **Collision-checked.** The samples move straight when the controls do,
  and that motion is checked. The curve between the samples is not.

## Derivatives of a converged solution (R5)

```c
qaws_repulsive_targets d = { 0 };
d.pinned = dpins;                     /* the pinned positions' changes */
d.edge_lengths = dl;                  /* curves; or d.length; meshes: d.area, d.volume */
qaws_repulsive_curves_solver_solution_jvp(dx, &gauge, solver, &d);
qaws_repulsive_curves_solver_solution_vjp(&targets_bar, &gauge, solver, xbar);
```

At a converged x*, dE = C^T lambda and every row holds its target. Moving
the targets moves the solution by the KKT system of the Lagrangian:
[K C^T; C 0] [dx; mu] = [0; dtargets], with K = Hess E - sum lambda_r
Hess c_r. The descent's preconditioner is not differentiated.

- **The terms.** lambda is C's least-squares multiplier of dE. K is applied
  exactly: the energy's HVP and the constraints' Hessians (edge lengths;
  area and volume per face).
- **The solve.** Flexible GMRES, preconditioned by the descent's dense
  Sobolev system (the same principal part). The VJP is the same solve with
  [xbar; 0], the system being symmetric.
- **Symmetries.** The rigid motions every constraint row allows are found
  as the null space of (C M)^T (C M), over the translations and rotations
  M, and gauged out: dx is orthogonal to them, and *out_gauge counts them.
  A loop held by its barycenter keeps three rotations; a strand pinned at
  both ends keeps the rotation about the line through its pins.
- **Status.**
  - NOT_CONVERGED unless the solver's last reason is CONVERGED.
  - NUMERICAL_FAILURE for dependent rows, or a system GMRES cannot reduce
    to 1e-7 (not a strict minimum on the constraint set), or not within
    the desc's `gmres_*` limits.
  - UNSUPPORTED_OPERATION in the hierarchical mode, and for the surface
    and control-point solvers.
- **Validity report.** The `_checked` variants of the four functions take
  a `qaws_repulsive_derivative_report` (the unchecked ones pass NULL). It
  is filled when the KKT solve ran (OK, or NUMERICAL_FAILURE because GMRES
  did not converge), zeroed when the call was refused before it:
  - `residual`: |b - K u| / |b| of the returned solution, recomputed with
    the exact operator;
  - `iterations`: GMRES iterations over all restarts;
  - `converged`: 1 when the solve was accepted (under `gmres_tolerance`,
    or under 1e-7 once restarts stopped gaining), else 0 and the call
    returns NUMERICAL_FAILURE with the outputs untouched;
  - `condition_estimate`: the ratio of the extreme singular values of
    GMRES's Hessenberg matrices (the largest over the restart cycles).
    These lie within the singular values of the preconditioned operator
    K P^-1, so it is a lower bound of that operator's condition number,
    not of K's. It is >= 1, HUGE_VAL for a singular cycle, 1 when no
    iteration ran;
  - `gauge`: as `*out_gauge`.
- **Not provided.** The derivative of a fixed sequence of steps, and
  derivatives in the exponents.
- **A practical note.** The descent converges well when the problem is
  well posed: every edge's length held for curves, no vertex free to slide
  within a surface for meshes. Holding only the total length leaves a
  curve's vertices free to bunch, and the descent crawls. The paper uses
  per-edge lengths, and remeshing for surfaces, for this reason.

| Check (test 118) | f64 | f32 |
|---|---|---|
| Control points (a B-spline, a Bezier, a NURBS): the gradient against central differences; the HVP against differences of the gradient; t and t2 against both | 3e-10; 3.7e-9; 2.9e-15 | 4.6e-5; 2.7e-4; 1.3e-7 |
| Control-point solves, length held, dense and Barnes-Hut: E 26 -> 1.5e-4 and 0.17 in 15 steps; the curves rebuilt reproduce the solved samples | 1.7e-13 | 7.6e-5 |
| A strand pinned at both ends, edge lengths held, a pin moved: the JVP against re-solved differences (less the free rotation about the pins' line, the one gauge) | 1.9e-5 (limited by the re-solves' convergence) | 1.2e-3 |
| The adjoint identity xbar . JVP(d) = VJP(xbar) . d | to 12 digits | 1.4e-7 |
| A loop, barycenter and edge lengths held, scaled: three rotations gauged, the JVP orthogonal to them, against the re-solved difference less its rotation | 1.7e-18; 2.3e-8 | 1.6e-10; 2.4e-5 |
| A lumpy octahedron, area and volume changed: against the re-solved difference less its rotation (the differences' floor: 5e-4) | 2.8e-4 | 6.9e-3 |
| Refusals: before convergence; in the hierarchical mode | NOT_CONVERGED; UNSUPPORTED | same |
| The report on the converged strand (default limits): residual, iterations, converged, condition estimate, gauge; the unchecked JVP | 4.9e-12, 16, 1, 191, 1; bit for bit | not measured |
| `gmres_tolerance` 1e-4; `gmres_max_iterations` 2 | 1.1e-5 in 8 iterations, converged; NUMERICAL_FAILURE, converged 0, residual 0.995 | not measured |
| A mesh's VJP, a pin, barycenter, area and volume held (one gauge: the rotation about the pin-barycenter line): the adjoint identity over three random pairs over every target | 6.4e-13 relative | not measured |
| The same mesh, the area target alone moved by +-4e-3 and +-2e-3 from the solution (re-solved to convergence, Richardson extrapolated): the JVP against it, less the free rotation | 6.1e-4 (bound 2e-3: the re-solves stop at stationarity 1e-6) | not measured |

Figure 36 (`examples/better2d_showcase/36_repulsive_surface.c`) shows five
crowded loops on a torus spreading over it with every edge's length held,
and two B-splines, their ends pinned, relaxed through their control
points.

## On the GPU (R6)

`core/qaws_repulsive_core.h` holds the pair terms as backend-neutral
kernels (C, HLSL, GLSL):

- the curve and mesh pair energies, with their gradients in either
  element's vertices, in closed form;
- the metric's pair weights (kB, kB0), and the per-element derivative
  coefficients.

The loops over pairs belong to the caller. `examples/repulsive.hlsl`
drives them in two passes:

- **CSElement**, one thread per element I, looping over every element J:
  - I's energy share;
  - its gradient block sum_J d/dI [e(I, J) + e(J, I)];
  - its row of the metric product.
- **CSGather**, one thread per vertex: the gradient and A u from the
  elements holding the vertex, through a CSR adjacency.

No atomics, so results do not depend on scheduling. Storage is
O(elements), never N^2. The energy is the shares' sum, taken on the host
in double.

Two numerical choices matter in single precision:

- **The cross product itself.** |T x d|^2 is computed from it, not as
  |T|^2 |d|^2 - (T . d)^2, which loses every digit for neighbours nearly
  in line.
- **Powers by products.** Powers whose doubled exponent is a whole number
  (alpha / 2, beta / 2, p) are computed by products and a square root, not
  by the GPU's pow (exp2 of log2, about 5e-6 off).

**What runs on the GPU:** the direct energy, its gradient, and the
metric's product.

**What does not:** HVPs, the Barnes-Hut traversal, the constraint rows,
and the solves themselves. The host orchestrates (line search, GMRES,
projection) with these products. This is not a full GPU solver.

`examples/repulsive_gpu_parity.c` (Direct3D 11, Windows only, outside the
library) runs the kernels on the GPU and compares them with the CPU
library in double, given the same (float-rounded) positions:

| Case (RTX 3060 Laptop) | Energy | Gradient | Metric product | GPU, one pass (all pairs, with the metric product) | CPU direct | CPU Barnes-Hut (theta 0.25) |
|---|---|---|---|---|---|---|
| Trefoil, 300 edges | 1.7e-9 | 1.4e-5 | 7.7e-7 | 1 to 4 ms | 0.24 s | 0.03 s |
| Trefoil, 2000 edges | 3.1e-8 | 3.8e-4 | 1.8e-6 | 6 to 25 ms | 11.9 s | 0.38 s |
| Trefoil, 8000 edges | 2.6e-9 | 1.1e-4 | | 0.10 s | 187 s | 1.7 s |
| Lumpy sphere, 512 faces | 7.6e-8 | 2.3e-6 | 9.3e-7 | 3 to 13 ms | 0.34 s | 0.07 s |
| Lumpy sphere, 2048 faces | 3.4e-8 | 9.4e-6 | 1.1e-6 | 10 to 50 ms | 5.7 s | 0.83 s |
| Lumpy sphere, 8192 faces | 2.5e-7 | 3e-5 | | 40 to 70 ms | 87 s | 5.0 s |

GPU times vary from run to run with the machine's load; ranges are shown.

**Single-precision positions.** For fine curves, positions in single
precision limit the gradient more than the arithmetic does. Rounding a
trefoil's coordinates (about 2 in size) to float moves its exact gradient
by 2e-3 at 300 edges, 0.65 at 2000, and 40 at 8000. The gradient of a
finely sampled curve is a small difference of large neighbour terms,
sensitive to each edge's position. The same holds for the f32 library's
inputs. The kernels match the double CPU to 4e-4 at worst when given the
same positions. Meshes are far less sensitive (1.2e-4 at 8192 faces).

| Check (test 119: the core on the C backend, driven as the kernels drive it) | f64 | f32 |
|---|---|---|
| Against the library: energy, gradient and metric product for a trefoil with a ring (3D), an open wave (2D, z = 0), and a lumpy sphere | 2.7e-15 at most | 1.2e-6 at most |

Figure 37 (`examples/better2d_showcase/37_repulsive_packing.c`) shows
curves packed into a region. A container is a pinned loop; collision
checking keeps the curves inside. One folded strand fills a disc, and
three strands fill a square.

### The showcase scenes

| Scene (from the specification) | Figure |
|---|---|
| Pinned strands around obstacles | 32 |
| Curves packed into a region | 37 |
| Curves distributed over a surface | 36 |
| A knot relaxed without crossings | 32 (a collision-checked trefoil), 35 (a 600-edge torus knot) |
| A constrained surface optimized near another surface | 33 |

## The paper's figures

`examples/repulsive_paper.c` (target `qaws_repulsive_paper`) reproduces
figures of Repulsive Curves from the paper's own scenes, copied under
their MIT license into `examples/data/repulsive_curves`. Each figure uses
the scene's curve, its (alpha, beta), and its constraints: every edge's
length and the barycenter. The dense mode runs until it converges
(stationarity 1e-4), on a 20-thread executor:

| Figure | Scene | Vertices | Steps | Energy | Time |
|---|---|---|---|---|---|
| 2 | Freedman's unknot to the circle (alpha 2, beta 4.5) | 592 | 450 | 2281 -> 6.08 | 247 s |
| 27 | K5 embedded in space (alpha 3, beta 6) | 635 | 235 | 1332 -> 1326 | 123 s |
| 27 | K3,3 embedded in space | 573 | 394 | 420 -> 359 | 136 s |
| 33 | an earbud cable untangled | 529 | 585 | -> 104.7 | 159 s |
| 26 | a tree drawn in the plane (2D) | 793 | 464 | 1.4e5 -> 1.98e4 | 337 s |
| 31 | four agents' paths planned in space-time (below) | 772 | 400 (not converged) | 1433 -> 593 | 531 s |

Figure 31 stops at 400 steps. By then the paths are set, and the agents'
closest approach in a time slice has grown from 0.08 to 0.84 in a
crossbar 2 wide. The tail to stationarity 1e-4 is slow: 3e-4 at 750
steps, the energy still falling by 0.5 every 50.

`qaws_repulsive_paper data-directory [figure [most steps]]` runs one
figure (1 to 10), or caps every solve for a preview.

Not reproduced: the interpolation scenes (TangentInterpolation and
CurveInterpolation). They need objectives on the lengths' differences.


### Repulsive Surfaces

The same program reproduces figures of Repulsive Surfaces. The Hopf link,
the Borromean rings and the genus-2 scenes are the paper's meshes, copied
under their MIT license into `examples/data/repulsive_surfaces`. The cylinders are
generated, in the paper's proportions. Each runs the hierarchical mode
(theta 0.25, p 6) until it stalls at the hierarchical stationarity floor
(about 6e-3). The surfaces are rendered by a small z-buffered rasterizer
(smooth shading from the vertex normals, supersampled), one PNG per panel
beside the figure's SVG:

| Figure | Scene | Faces | Steps | Energy |
|---|---|---|---|---|
| 9 (left) | a cylinder twice as tall as wide, both rims pinned, ballooning to a near-sphere | 3360 | 82 | 2.244 -> 0.397 |
| 9 (centre) | two pinned rings 3 apart (radius 1), the area weighted 1, 10 and 100: a neck kept in each | 2592 | 30, 29, 38 | 18.6, 134.9, 1077 |
| (a scene of the paper's code) | a Hopf link, the area held | 9216 | 40 | 842.5 -> 303.3 (256 s) |
| 1 | linked handcuffs (genus 2) to the unlinked double torus, the area held, remeshed after every step | 10506 | 700 | 376 -> 7.551 (75 min) |

The links hold the area and the barycenter, as the paper's own scene
files do (with limits of 200 and 350 iterations): unconstrained, they
inflate (the energy falls toward 0 as the surfaces grow).

The paper's program remeshes the surface after every step by default (a
setting of its interface, not of its scene files). Without it, its
Fig. 1 never unlinks: after the scene's 700 iterations it is still
interlocked (E 76.1 exact; this solver 77.1), and 300 more steps leave
it so (71.2). The
figures of closed surfaces therefore remesh as it does, without its edge
splits and collapses (the reference run without them unlinks the same
way): five rounds of smoothing (each vertex half way to the area-weighted
circumcentres of its faces, in its tangent plane) and Delaunay flips,
the area scaled back, and the solver made again, resuming its step. A
face the smoothing would turn over keeps its vertices where they were:
unguarded, one round on a mesh that had not been remeshed folded 81 face
pairs through each other. Remeshed, Fig. 1 unlinks in about 300 steps;
after 700 its energy is 7.551, the paper's program 7.548 (exact).

Against the reference implementation (the paper's code, scored by its
exact all-pairs energy): on the Fig. 9 (left) cylinder it reaches 0.3823
and this solver 0.3852; on the Hopf link 303.25 against 303.69. On the
Borromean rings it descends much faster (4686 against 7033 after 350
steps), but only by passing the rings through each other (24 to 28
intersecting face pairs between rings over its first 9 steps; the linking
numbers come back 0); this solver's step, held by its largest vertex
move, keeps them apart.

Not reproduced: implicit obstacles, growth schedules, remeshing of surfaces
with a boundary, or p <= 4 (the solver needs p > 4).

## Space-time paths and rooms (Sec. 8.5)

The paper plans the paths of agents moving in a 2D room as curves in
(x, time, z). If the paths are kept apart in space-time, they are kept
apart in every time slice. The environment is the room extruded along
time. Three settings of `qaws_repulsive_curves_solve_desc` provide this:

- **`hold_coordinate`** (1 to the dimension) holds that coordinate of
  every vertex, as one constraint row each. Time then stays exact: each
  path remains a graph over time, and no agent moves back in time. The
  paper does not enforce this; its length term keeps the braid stretched
  in practice. The barycenter's row along it is dropped, since the held
  rows already fix it. These rows are too many for the hierarchical mode,
  so only the dense mode accepts it.
- **`length_weight`** adds the weight times the total length to the
  energy. This is the paper's `optimize_length`.
- **`room`** is a simple polygon that the vertices stay inside. In 2D it
  lies in the plane. In 3D it lies in the plane of the other two
  coordinates, extruded along the held one (or z). Each vertex is
  repelled by room_weight times the sum over walls of the integral of
  rho^-2 along the wall, where rho is the distance in the plane. Per
  wall, in closed form:
  - the integral is atan2(|h| (b - a), h^2 + a b) / |h|, where h is the
    vertex's height over the wall's line and [a, b] is the wall along it;
  - near the line and beside the wall, a series takes over;
  - this is half the integral of r^-3 over the wall extruded without end.
  The paper's obstacle is a mesh's potential, p = beta - alpha = 3 in its
  scene. On the extruded room, that potential is twice this term, so the
  paper's weight 0.05 is room_weight 0.1 here. A vertex must start
  strictly inside: creation refuses one that does not. A trial step that
  leaves the room has no energy, so it is refused.

These terms have no Hessian products, so a solver using them refuses
solution derivatives (UNSUPPORTED_OPERATION).

Test 120 checks the room term against Mathematica integrating the
definition (`tests/reference/120_repulsive_room.wls`). On the paper's H,
at five points, the energy agrees to 7.6e-16 and the gradient to 1.6e-15
relative. One of the points lies 1e-6 off a wall's line, beside the wall,
so the series is checked too. The 2D integral also matches half the 3D
integral of r^-3. Finally, four agents swap corners through a square room
(65 vertices each). The solve converges in 32 steps with times held to
5e-15 and every vertex inside after every step. The closest agents in a
time slice go from 0.12 to 2.07 apart: they share a roundabout. With 25
vertices each, the time steps were too coarse and the solve settled in a
pinwheel instead. There, the four agents trade places within a single
edge, 0.59 apart.

## CPU performance

The ideas come from the acceleration structures of TinyBVH, nanoflann and
Farq (Farq's `farq_bvh.hpp`: median and Morton builds, 8-wide nodes with
outward-rounded float bounds, AVX box tests, a parallel build). They are
applied where the profile pointed: the per-interaction cost first, then
threads, then the collision checks' broad phase.

- **Closed forms.** The Barnes-Hut far terms are closed forms, value and
  gradient, like the near pairs. No jets remain on the hierarchical path.
  Their accuracy is unchanged to the digit (test 115).
- **The dense mode** takes its energy from Barnes-Hut with theta 0, which
  opens every node, so it computes the direct terms in closed form. The
  jets' evaluator, kept for the Hessian products, is 6 times slower. The
  dense mode's constrained system is solved by blocks (see Direction
  above).
- **Powers.** Every pair term takes powers. When twice the exponent is a
  whole number (alpha / 2, beta / 2 and p usually are), the power is
  computed by products and one square root.
- **An executor** (`qaws_repulsive_curves_desc.executor`, and the mesh's)
  runs the work in parallel, following Qaws's convention: Qaws has no
  threads of its own, and results do not depend on the executor. Test 115
  runs the chunks backwards and checks bit-identical results.
  - **Barnes-Hut:** each element writes only its own energy share and
    gradient block; a gather in element order follows.
  - **The interaction lists:** built in chunks of elements and joined in
    element order.
  - **The hierarchical direction:** its scalar GMRES solves run at once,
    each with its own workspace (the products' scratch is per solve).
    Identical right-hand sides, such as the barycenter rows, are solved
    once.
- **The collision checks' broad phase** is a self traversal of the same
  median-split tree over each element's swept box, grown by half the
  clearance. Pairs of boxes farther apart than the clearance cannot come
  closer than it, since a straight motion stays in the box of its two
  positions. Qaws's uniform grid sizes its cells for boxes spread through a
  volume; along a curve its cells hold hundreds of edges each.
- **Fixes found on the way:**
  - the hierarchical GMRES now accepts the rounding floor (1e-6) once its
    restarts stop gaining; at 8000 edges the floor sat just above 1e-8 and
    every step failed;
  - a solver's creation now takes its initial energy in its own mode; the
    direct one cost minutes at 8000 edges.

| Lumpy trefoil (f64, 20 threads) | Before | After, serial | After, threads |
|---|---|---|---|
| Barnes-Hut energy and gradient, 8000 edges | 1.63 s | 0.24 s | 0.035 s |
| Barnes-Hut energy and gradient, 32000 edges | 10.9 s | 1.06 s | 0.15 s |
| Hierarchical step, 2000 edges | 1.7 s | 1.5 s | 0.36 s |
| Hierarchical step, 8000 edges | failed (GMRES) | 10.6 s | 1.9 s |
| Collision check of one step, 8000 / 32000 edges | 2.6 s / 61 s (the grid; all pairs before it) | 4.4 ms / 18 ms | |
| Dense step, Freedman unknot (592 vertices, every edge's length held: k = 595) | 10 to 15 s | about 1 s | |

Not done here:

- **Wide nodes with SIMD box tests** (Farq's 8-wide AVX traversal). After
  the closed forms, the traversal is no longer what dominates.
- **A Morton build.** The tree's build is a small share.
- **Parallelism inside one metric product.** With fewer solves than
  threads, the extra threads stay idle; at 8000 edges a step's 7 solves use
  7 of 20.
