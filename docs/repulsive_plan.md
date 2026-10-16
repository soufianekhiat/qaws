# Repulsive curves and surfaces: specification

Tangent-point energies for curves and surfaces, their derivatives, and
later their constrained solvers, as differentiable geometry operations of
qaws. This document fixes the discrete models, the conventions, the
contracts and the order of work. It is the reference the code and the
tests are checked against.

## 1. References (pinned)

| Work | Source | Revision used |
|---|---|---|
| Repulsive Curves (Yu, Schumacher, Crane, TOG 2021) | paper, supplemental: cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves | the published PDF |
| Repulsive Curves, reference implementation | github.com/ythea/repulsive-curves (MIT, (c) 2019 Christopher Yu) | `78a0c15da6edfc19a3f6607aebc59d769b1e59da` |
| Repulsive Surfaces (Yu, Brakensiek, Schumacher, Crane, TOG 2021) | paper: cs.cmu.edu/~kmcrane/Projects/RepulsiveSurfaces | the published PDF |
| Repulsive Surfaces, reference implementation | github.com/ythea/repulsive-surfaces (MIT, (c) 2019 Nick Sharp) | `c51b4643ba3886e711daa136d32935eba9603bac` |

No code is copied from the reference implementations: qaws's evaluator is
written from the papers' formulas, in C, and checked against independent
Mathematica fixtures. Should any code be ported later (the hierarchical
trees, the remeshing), its file keeps the MIT notice of its origin. The
reference applications and their dependencies (polyscope, Eigen,
geometry-central) never enter the qaws runtime.

## 2. Discrete curve energy: the paper's, not TPESC's

The smooth energy (paper, Eq. 3) is

    E(gamma) = integral over M x M of k(gamma(x), gamma(y), T(x)) dx dy,
    k(p, q, T) = |T x (p - q)|^alpha / |p - q|^beta,

with |T x (p - q)| = |(I - T T^T)(p - q)|, both measures arc length.

The paper and the reference code discretize it differently:

| | Paper, Sec. 5.1 (Eqs. 16 to 18) | Reference `TPESC` (`src/tpe_energy_sc.cpp`) |
|---|---|---|
| Points | edge endpoints, trapezoidal rule on each edge pair | vertices |
| Tangent | the edge's, T_I = (gamma_i2 - gamma_i1) / l_I | the vertex's: the normalized sum of its edges' unit tangents |
| Weights | l_I l_J / 4 per endpoint pair | the vertex's mean incident edge length (an end vertex gets its whole edge, an inner one the mean of two) |
| Excluded | edge pairs sharing a vertex | only i = j; a vertex of more than two edges contributes nothing as the tangent point |
| Pairs | ordered (I, J) and (J, I) | ordered (i, j) and (j, i) |

qaws implements the paper's discretization:

    E = sum over ordered pairs (I, J) of edges with no common vertex of
        l_I l_J (1/4) sum over i in I, j in J of k(gamma_i, gamma_j, T_I).

Reasons: it is the published model (Fig. 8, Eq. 17), it has a tangent
defined everywhere (an edge always has one; a vertex of a network has
none), and the exclusion is by topology, the same rule for one curve,
several curves and networks. Comparisons are made against this exact
formula only: qaws's C against Mathematica's evaluation of Eq. 17 on the
same polygons. The TPESC vertex model is not compared against, since it is
a different discrete energy (it converges to the same smooth energy under
refinement). A vertex variant may be added later for parity with the
reference code; it would be a separate, named model.

Defaults (alpha, beta) = (3, 6), the paper's choice outside its Figs. 14 to
19 (which use (2, 4.5)). The energy is defined for alpha > 1 and
beta in [alpha + 2, 2 alpha + 1] (paper, Sec. 3.2); qaws requires
alpha >= 2 as well, so that the energy is twice differentiable where a
point lies on another edge's line (|T x d| = 0; see 5.3), and refuses other
exponents with INVALID_ARGUMENT.

Scaling: every term is homogeneous of degree 2 + alpha - beta in the
coordinates, so E(c X) = c^(2 + alpha - beta) E(X) exactly (c^-1 for the
defaults). It is invariant under rigid motions and under reversing an
edge's orientation (T enters only through |T x d|).

The input is a graph: vertices and oriented edges. One curve, several
curves (mutual repulsion), open and closed, in 2D (|T x d| as the 2D
cross product) or 3D: all the same formula. A vertex of three or more
edges is a junction: the edges meeting there never interact with each
other (they share it), and edges that do not touch it interact as usual.
That is the only junction rule; explicit junction semantics for networks
(R2+) will refine it, never widen it to a neighbourhood.

## 3. Discrete surface energy: the paper's, which the code shares

The smooth energy (paper, Eqs. 2 and 3, the constant 2^p omitted) is

    E^p(f) = integral over M x M of |P(x)(f(x) - f(y))|^p / |f(x) - f(y)|^(2p) dx dy,

P(x) = N N^T the projector on the normal. The paper's discretization
(Eqs. 9 and 10) uses one-point (barycenter) quadrature on each face:

    E = sum over ordered pairs (S, T) of faces, S != T, of
        a(S) a(T) |N_S . (X_S - X_T)|^p / |X_S - X_T|^(2p),

X the barycenter, a the area, N the unit normal. The reference code
(`src/energy/tpe_all_pairs.cpp`) computes the same sum, as unordered pairs
i < j with the kernel symmetrized (|N_i . v|^alpha + |N_j . v|^alpha) /
|v|^beta, alpha = p, beta = 2p; no pair is excluded (faces sharing a
vertex or an edge have distinct barycenters). qaws implements this. Default
p = 6 (the paper's throughout); p > 4 is required (paper, Sec. 2.1), and
qaws takes p as a real number > 4. Scaling: degree 4 - p (c^-2 for p = 6).

The input is a triangle mesh as a positions/index view: positions
(3 per vertex) and triangles (3 indices each). Connectivity is fixed for
the lifetime of a prepared object. Several components repel mutually.

## 4. Conventions

- Ordered pairs: both energies sum every ordered pair once (equivalently
  each unordered pair with the symmetrized kernel). No factor 1/2.
- The tangents, normals, lengths and areas are functions of the vertex
  positions, and every derivative includes them.
- Arithmetic: double internally in both builds; inputs and outputs are
  qaws_scalar. f32 builds round the outputs.
- No softening. A pair whose points coincide (a self-intersection) has
  infinite energy: the call fails with NUMERICAL_FAILURE, never returns a
  finite number. A zero-length edge or a zero-area face has no tangent or
  normal: DEGENERATE_CURVE. A regularized model, if ever useful, will be a
  separate, named energy.
- Repulsion is not clearance: the energy's weight is not a strand radius or
  a shell thickness. Clearance belongs to the collision-checked mode (R2).

## 5. Derivatives (R1)

Exact derivatives of the discrete energy above, on the same branch (the
same pairs, exclusions and quadrature):

- value, tangent (first directional derivative, a JVP) and tangent2 (the
  second directional derivative along the same direction);
- gradient, accumulated (+=) and scaled by a seed (a VJP of the scalar);
- Hessian-vector product, accumulated (+=).

5.1 Method. Each pair term is a function of 4 vertices (curves: the two
edges' endpoints) or 6 (surfaces: the two triangles' corners). It is
evaluated once in a second-order forward jet: value, the gradient in all
its variables, the derivative along the direction, its second derivative
along it, and the derivative of the gradient along it (the HVP). The jet's
rules are those of hyper-dual numbers with a vector first slot. The pair's
results are scattered to its vertices. No finite differences.

5.2 Validation: Mathematica (symbolic derivatives of the same discrete
sums, evaluated at 30 digits) for energy, gradient, tangent, tangent2 and
HVP; adjoint identity <grad, v> = tangent; Hessian symmetry
<u, H v> = <v, H u>; finite-difference sweeps of the gradient;
translation and rotation invariance; the scaling degree.

5.3 Where the energy is not smooth. |T x d|^alpha with T x d = 0 (a point
on another edge's line, e.g. any straight open curve) is C^2 for
alpha >= 2: the jet's power rule takes value, gradient and Hessian 0 there
(alpha > 2) or the square's (alpha = 2), the limits of the smooth formulas.
(N . d)^p likewise for p > 4.

5.4 What R1 does not differentiate: the solver (R2+). Later contracts keep
three things apart: the energy's derivatives (here), derivatives of a
fixed sequence of solver steps, and derivatives of a converged constrained
solution (through the KKT system, with its rank and nullspaces checked).
The Sobolev-preconditioned descent direction is not the gradient and is
never returned as one. An unavailable derivative is an error status, never
a zero.

## 6. API (R1)

    qaws_repulsive_curves_create(&curves, desc);        topology and exponents, validated
    qaws_repulsive_curves_eval(&e, &t, &t2, curves, positions, direction);
    qaws_repulsive_curves_gradient(gradient, &e, curves, positions, seed);
    qaws_repulsive_curves_hvp(hv, curves, positions, direction);
    qaws_repulsive_curves_destroy(curves);

    qaws_repulsive_mesh_create(&mesh, desc);            positions/index view
    qaws_repulsive_mesh_eval / _gradient / _hvp / _destroy

A prepared object holds the topology and the exponents only. Positions
come with each call, so an optimizer evaluates many configurations against
one object. Objects are immutable after creation; calls may run
concurrently. No global state.

## 7. Guarantees and non-guarantees

| Operation | Kind | Contract |
|---|---|---|
| energy, tangent, tangent2, gradient, HVP (R1) | exact for the discrete model | the discrete sums of sections 2 and 3, in double, direct all-pairs (O(n^2)): no hierarchical approximation |
| energy-based optimization (R2/R3) | numerical | reduces the energy; does not promise collision-free motion |
| collision-checked optimization (R2/R3) | checked | accepted updates (after constraint projection) validated by continuous collision tests: segment-segment for curves, vertex-face and edge-edge for meshes; invalid starts refused or sent to a separate repair; no claim of global untangling |
| topology preservation | checked, when requested | a crossing in an accepted step is a rejection, so a knot cannot pass through itself numerically |

"Direct" means the discrete sum without approximation, not the continuum
integral: the discrete energy differs from the smooth one by the
quadrature error.

## 8. Milestones

| Milestone | Deliverable | State |
|---|---|---|
| R0 | this specification: models, reconciliation, pinned references, contracts | done |
| R1 | direct curve and surface energies with first- and second-order derivatives; Mathematica fixtures; invariance, scaling, adjoint, symmetry and refusal tests in f64 and f32 | done: `qaws_repulsive.h`, test 112 (Mathematica within 1.1e-14 in f64, 1.4e-6 in f32), figure 31 |
| R2 | constrained curve solver (fractional Sobolev metric with its low-order term, constraints, translation nullspace), collision-checked stepping, polyline and qaws curve adapters | done: `qaws_repulsive_solve.c`, `qaws_repulsive_adapt.c`, test 113 (metric against Mathematica within 2.2e-16), figure 32 |
| R3 | surface solver (fixed connectivity, pins, global area/volume, the Poisson-based preconditioner), collision checks | done: `qaws_repulsive_mesh_solve.c` on the shared core, test 114 (metric against Mathematica within 4.8e-16; preconditioned GMRES 13 iterations against 510), figure 33. L^(2-s) by vertex quadrature: the face-average form has a spurious kernel on 3-colourable meshes |
| R4 | hierarchical energy/gradient and matrix-free operators, convergence to the direct evaluator, benchmarks | done: Barnes-Hut energy and gradient (test 115, figure 34); the metric from Barnes-Hut lists, the iterative constrained solve (pins eliminated, Schur complement, GMRES with the Laplacian sandwich, profile Cholesky), dense against hierarchical benchmarks (test 116, figure 35); collision checks through a tree broad phase |
| R5 | surface-constrained curves (repulsion in the embedded 3D geometry, seams and charts handled), native control-point optimization, differentiated solves | done: curves on parametric surfaces in their parameters (seams declared, ranges, trims and poles refused; one chart), control-point optimization through the curves' own tangents, derivatives of converged solutions by the KKT system with symmetries gauged (tests 117, 118, figure 36); not provided: derivatives of step sequences; surface curves and controls dense only |
| R6 | HLSL kernels, backend parity, docs, showcases | done: backend-neutral pair kernels (core/qaws_repulsive_core.h), HLSL energy, gradient and metric product (examples/repulsive.hlsl), run on a GPU against the CPU (examples/repulsive_gpu_parity.c), test 119, figure 37; on the GPU: energy, gradient and metric products only, the solve orchestrated by the host |
