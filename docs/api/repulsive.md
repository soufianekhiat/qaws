# API Reference: qaws_repulsive.h

Tangent-point (repulsive) energies of polygonal curves and triangle
meshes, with exact first and second derivatives. After Yu, Schumacher and
Crane, *Repulsive Curves* (2021), and Yu, Brakensiek, Schumacher and Crane,
*Repulsive Surfaces* (2021). The discrete models, the reconciliation of the
papers with their reference code, and the contracts are in
[`docs/repulsive_plan.md`](../repulsive_plan.md).

This is milestone R1: the energies and their derivatives, evaluated
directly over all pairs. There are no solvers yet (R2 for curves, R3 for
meshes), so nothing here moves geometry or promises collision-free motion.

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
