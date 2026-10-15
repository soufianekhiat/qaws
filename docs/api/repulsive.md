# API Reference: qaws_repulsive.h

Tangent-point (repulsive) energies of polygonal curves and triangle
meshes, with exact first and second derivatives. After Yu, Schumacher and
Crane, *Repulsive Curves* (2021), and Yu, Brakensiek, Schumacher and Crane,
*Repulsive Surfaces* (2021). The discrete models, the reconciliation of the
papers with their reference code, and the contracts are in
[`docs/repulsive_plan.md`](../repulsive_plan.md).

Milestones R1 and R2: the energies and their derivatives, evaluated
directly over all pairs, and the curves' constrained fractional Sobolev
solver with a collision-checked mode and adapters to qaws curves. Meshes
have no solver yet (R3).

```c
qaws_repulsive_curves_desc d = { 3, vertex_count, edges, edge_count, 0, 0 };   /* (alpha, beta) = (3, 6) */
qaws_repulsive_curves* c = NULL;
qaws_repulsive_curves_create(&d, &c);
qaws_repulsive_curves_eval(c, positions, direction, &e, &de, &d2e);   /* E, dE[v], d2E[v, v] */
qaws_repulsive_curves_gradient(c, positions, 1, gradient, &e);        /* gradient += seed dE/dx */
qaws_repulsive_curves_hvp(c, positions, direction, hv);               /* hv += H v */
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
qaws_repulsive_curves_solve(c, &sd, positions, &report);              /* in place */

qaws_repulsive_curves_solver* s = NULL;                               /* or step by step */
qaws_repulsive_curves_solver_create(c, &sd, positions, &s);
qaws_repulsive_curves_solver_step(s, positions, &report);             /* report.reason: RUNNING or why it stopped */
qaws_repulsive_curves_solver_direction(s, positions, g);              /* the descent direction, not the gradient */
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
  is the metric on each coordinate and C is the constraints' Jacobian. The
  saddle-point system is solved as such, by dense LU with partial pivoting
  (not CG). The result g is the descent direction in the metric, and
  `_direction` returns it. It is not the coordinate gradient
  (`qaws_repulsive_curves_gradient`); neither stands in for the other.
- **Step.** Backtracking on the largest vertex move h:
  - h starts at twice the last accepted move (a quarter of the mean edge
    length at first);
  - the trial point is y = x - (h / |g|_inf) g;
  - y is projected back onto the constraints by Newton steps in the same
    metric, [A~ C(y)^T; C(y) 0] [dx; m] = [0; -phi(y)];
  - the step is accepted when Armijo holds on the projected y and, in
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

The report also counts the trial steps rejected by the energy, by the
collision check and by the projection. It gives the constraint residual and
the size of the last step.

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
