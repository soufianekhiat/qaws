#ifndef QAWS_REPULSIVE_H
#define QAWS_REPULSIVE_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"
#include "qaws_batch_executor.h"

/*
 * Tangent-point (repulsive) energies of polygonal curves and triangle
 * meshes, after Yu, Schumacher and Crane, "Repulsive Curves" (2021), and
 * Yu, Brakensiek, Schumacher and Crane, "Repulsive Surfaces" (2021). The
 * discrete models, their conventions and the contracts are specified in
 * docs/repulsive_plan.md; in short:
 *
 *   curves:   E = sum over ordered pairs (I, J) of edges with no common
 *             vertex of l_I l_J (1/4) sum_{i in I, j in J}
 *             |T_I x (x_i - x_j)|^alpha / |x_i - x_j|^beta
 *             (the paper's Eq. 17: trapezoidal rule on each edge pair,
 *             the edge's tangent T_I)
 *   surfaces: E = sum over ordered pairs (S, T), S != T, of faces of
 *             a_S a_T |N_S . (X_S - X_T)|^p / |X_S - X_T|^(2p)
 *             (the paper's Eq. 10: barycenters X, areas a, unit normals N)
 *
 * Direct all-pairs evaluation, O(n^2), in double in both builds: the
 * discrete sums themselves, no hierarchical approximation. Derivatives are
 * exact for these sums (a second-order forward jet per pair), including
 * the tangents, normals, lengths and areas' dependence on the positions.
 *
 * Positions come with each call (dimension per vertex, vertex major); a
 * prepared object holds only the topology and the exponents, is immutable
 * and may be used from several threads at once. On any failure the
 * outputs are left untouched.
 *
 * Statuses: INVALID_ARGUMENT for bad input or exponents outside the
 * energy's range; DEGENERATE_CURVE for a zero-length edge or a zero-area
 * face (no tangent or normal); NUMERICAL_FAILURE when two interacting
 * points coincide (a self-intersection: the energy is infinite, never
 * softened into a finite value).
 */

/* ---- curves ------------------------------------------------------------ */

typedef struct qaws_repulsive_curves_desc
{
	unsigned int dimension;            /* 2 or 3 */
	unsigned int vertex_count;
	unsigned int const* edges;         /* 2 vertex indices per edge; one curve, several, open or closed, or a network */
	unsigned int edge_count;
	qaws_scalar alpha;                 /* 0 and 0: (3, 6). alpha >= 2, beta in [alpha + 2, 2 alpha + 1] */
	qaws_scalar beta;
	qaws_batch_executor const* executor;   /* optional: runs the evaluations in parallel (the results do not depend on it) */
} qaws_repulsive_curves_desc;

typedef struct qaws_repulsive_curves qaws_repulsive_curves;

qaws_status qaws_repulsive_curves_create(qaws_repulsive_curves_desc const* desc, qaws_repulsive_curves** out_curves);
void qaws_repulsive_curves_destroy(qaws_repulsive_curves* curves);

/* The energy at positions; with a direction (else NULL, and the derivative
   outputs 0), its first and second derivatives along it (a JVP and the
   second directional derivative). Any output may be NULL. */
qaws_status qaws_repulsive_curves_eval(qaws_repulsive_curves const* curves, qaws_scalar const* positions,
	qaws_scalar const* direction, qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2);

/* gradient += seed * dE/dpositions (a VJP of the scalar energy); the
   energy in out_value (may be NULL). */
qaws_status qaws_repulsive_curves_gradient(qaws_repulsive_curves const* curves, qaws_scalar const* positions, qaws_scalar seed,
	qaws_scalar* gradient, qaws_scalar* out_value);

/* out_hv += H direction, H the Hessian of the energy at positions. */
qaws_status qaws_repulsive_curves_hvp(qaws_repulsive_curves const* curves, qaws_scalar const* positions, qaws_scalar const* direction,
	qaws_scalar* out_hv);

/* ---- curve optimization (R2) ------------------------------------------- */

/*
 * Fractional Sobolev descent (Repulsive Curves, Sec. 4.2 and 5.2 to 5.4):
 * the metric A = B + B0 (the high-order term on edge differences and the
 * low-order term, exponent 2 sigma + 1, sigma = (beta - 1) / alpha - 1),
 * assembled densely; the descent direction g from the constrained system
 * [A C^T; C 0] [g; l] = [dE; 0] (dense LU with partial pivoting, the
 * saddle-point system solved as such); a step along -g, its size
 * backtracked (Armijo); the result projected back onto the constraints by
 * Newton steps in the same metric, [A C^T; C 0] [x; m] = [0; -Phi].
 *
 * Constraints are held at their values at the solver's creation: pinned
 * vertices (an obstacle is a pinned component), edges whose direction is
 * held, the total length or every edge's length, and the vertex
 * barycenter (by default when nothing is pinned: the metric does not see
 * translations). Redundant sets (total and every edge length) are refused.
 *
 * Two modes. Energy-based: each accepted step lowers the energy and holds
 * the constraints; nothing promises the curve does not pass through itself
 * between or within steps. Collision-checked: every accepted motion (the
 * straight path from the old positions to the projected new ones, so the
 * projection's movement included) is validated by conservative advancement
 * on every pair of edges sharing no vertex: along it each pair stays at
 * least min(clearance, half its distance at the step's start) apart, so no
 * two edges ever touch and the topology (a knot's type, a link's linking)
 * is kept. An input whose edges already touch is refused. There is no
 * repair and no promise of untangling.
 *
 * Hierarchical mode (theta > 0, R4; also for meshes): nothing of size
 * |V|^2 is stored. The energy and gradient are the Barnes-Hut ones below;
 * the metric is applied from Barnes-Hut interaction lists, and the
 * direction is solved by flexible GMRES. The preconditioner is the
 * Laplacian sandwich (-L)^-1 L^(2 - s) (-L)^-1 (Repulsive Surfaces, Sec.
 * 5.2), with the Laplacian factored by a profile Cholesky. Pinned vertices
 * are eliminated, and the other rows are solved by their Schur complement:
 * at most 32 rows, so per-edge lengths are refused here. The Newton
 * projection reuses those solves. The step is an approximation of the
 * dense mode's, converging to it as theta falls. The collision check is
 * the dense mode's, through the same tree broad phase.
 */

typedef enum qaws_repulsive_termination
{
	QAWS_REPULSIVE_RUNNING = 0,            /* a step taken, not stopped */
	QAWS_REPULSIVE_CONVERGED,              /* stationarity under the tolerance */
	QAWS_REPULSIVE_ITERATION_LIMIT,        /* the iteration budget spent */
	QAWS_REPULSIVE_INFEASIBLE,             /* the constraints could not be met again */
	QAWS_REPULSIVE_COLLISION_LIMITED,      /* every shorter step collided (collision-checked) */
	QAWS_REPULSIVE_STALLED,                /* no step lowered the energy */
	QAWS_REPULSIVE_NUMERICAL_FAILURE       /* a singular system or an unusable energy */
} qaws_repulsive_termination;

typedef struct qaws_repulsive_curves_solve_desc
{
	unsigned int const* pinned;            /* vertices held where they are */
	unsigned int pinned_count;
	unsigned int const* tangent_edges;     /* edges whose direction is held */
	unsigned int tangent_edge_count;
	int hold_total_length;
	int hold_edge_lengths;
	int hold_barycenter;                   /* 0: when nothing is pinned; 1: always; -1: never */
	int collision_checked;
	qaws_scalar clearance;                 /* collision-checked; 0: 1e-3 of the mean edge length */
	unsigned int max_iterations;           /* solve: 0 = 100 */
	qaws_scalar tolerance;                 /* stationarity (see the report); 0 = 1e-4 */
	qaws_scalar constraint_tolerance;      /* |Phi| over the mean edge length; 0 = 1e-8 (1e-5 in f32) */
	qaws_scalar theta;                     /* 0: dense (exact metric, LU); > 0: hierarchical (Barnes-Hut energy, gradient and metric at this theta, GMRES; see below) */
} qaws_repulsive_curves_solve_desc;

typedef struct qaws_repulsive_report
{
	qaws_repulsive_termination reason;
	unsigned int iterations;               /* accepted steps */
	qaws_scalar initial_energy, energy;
	qaws_scalar stationarity;              /* <dE, g> l / (E |g|_inf): the energy's relative first-order fall
	                                          for a move of one mean edge length l along the direction */
	qaws_scalar constraint_residual;       /* |Phi|_inf / l */
	qaws_scalar last_step;                 /* the last accepted step's largest vertex move / l */
	unsigned int rejected_energy;          /* trial steps rejected: no Armijo decrease or no energy */
	unsigned int rejected_collision;       /* trial steps rejected by the collision check */
	unsigned int rejected_projection;      /* trial steps whose projection failed */
} qaws_repulsive_report;

typedef struct qaws_repulsive_curves_solver qaws_repulsive_curves_solver;

/* A solver for one configuration: the constraints' targets taken from
   positions; checks the input (collision-checked: no touching edges). */
qaws_status qaws_repulsive_curves_solver_create(qaws_repulsive_curves const* curves, qaws_repulsive_curves_solve_desc const* desc,
	qaws_scalar const* positions, qaws_repulsive_curves_solver** out_solver);
void qaws_repulsive_curves_solver_destroy(qaws_repulsive_curves_solver* solver);

/* One step: positions updated in place when a step is accepted; the report
   (may be NULL) says what happened (RUNNING, or why it stopped). */
qaws_status qaws_repulsive_curves_solver_step(qaws_repulsive_curves_solver* solver, qaws_scalar* positions, qaws_repulsive_report* out_report);

/* Steps until a termination reason (create, step, destroy). */
qaws_status qaws_repulsive_curves_solve(qaws_repulsive_curves const* curves, qaws_repulsive_curves_solve_desc const* desc,
	qaws_scalar* positions, qaws_repulsive_report* out_report);

/* The descent direction g at positions (the projected fractional Sobolev
   gradient: A g + C^T l = dE, C g = 0), not the coordinate gradient. */
qaws_status qaws_repulsive_curves_solver_direction(qaws_repulsive_curves_solver const* solver, qaws_scalar const* positions,
	qaws_scalar* out_direction);

/* ---- adapters: qaws curves <-> working polylines ------------------------ */

/*
 * A working polyline for qaws curves: each curve's parameter range split
 * until every edge is at most edge_length long and its middle within
 * tolerance of its chord. Every vertex keeps where it came from: its curve
 * and parameter. A curve whose ends coincide (within 1e-9 of its extent;
 * 1e-5 in f32) closes into a loop. All curves share the dimension given (2 or 3).
 */
typedef struct qaws_repulsive_polyline_desc
{
	qaws_curve const* const* curves;
	unsigned int curve_count;
	unsigned int dimension;
	qaws_scalar edge_length;           /* the longest edge; must be > 0 */
	qaws_scalar tolerance;             /* chord deviation; 0: edge_length / 100 */
} qaws_repulsive_polyline_desc;

/* Positions (dimension per vertex), edges (2 per edge), and per vertex its
   source curve and parameter (either may be NULL). BUFFER_TOO_SMALL with
   the counts set when a capacity is short (call with 0 to size). */
qaws_status qaws_repulsive_polyline_from_curves(qaws_repulsive_polyline_desc const* desc,
	qaws_scalar* out_positions, unsigned int vertex_capacity, unsigned int* out_edges, unsigned int edge_capacity,
	unsigned int* out_source_curve, qaws_scalar* out_source_t, unsigned int* out_vertex_count, unsigned int* out_edge_count);

/* Curve `curve`'s vertices (as produced above, with their source parameters)
   at new positions, as a degree-1 B-spline whose knots are those
   parameters: exactly the working polyline, parametrized as the source
   (a closed curve's loop closed at its end parameter). No smoothing: what
   was checked on the polyline holds for this curve. */
qaws_status qaws_repulsive_polyline_to_curve(unsigned int dimension, qaws_scalar const* positions, unsigned int vertex_count,
	unsigned int const* source_curve, qaws_scalar const* source_t, unsigned int curve, qaws_curve const* source,
	qaws_curve** out_curve);

/* ---- triangle meshes --------------------------------------------------- */

typedef struct qaws_repulsive_mesh_desc
{
	unsigned int vertex_count;         /* positions are 3 per vertex */
	unsigned int const* triangles;     /* 3 vertex indices per face; connectivity fixed */
	unsigned int triangle_count;
	qaws_scalar p;                     /* 0: 6. p > 4 */
	qaws_batch_executor const* executor;   /* optional, as for curves */
} qaws_repulsive_mesh_desc;

typedef struct qaws_repulsive_mesh qaws_repulsive_mesh;

qaws_status qaws_repulsive_mesh_create(qaws_repulsive_mesh_desc const* desc, qaws_repulsive_mesh** out_mesh);
void qaws_repulsive_mesh_destroy(qaws_repulsive_mesh* mesh);

qaws_status qaws_repulsive_mesh_eval(qaws_repulsive_mesh const* mesh, qaws_scalar const* positions,
	qaws_scalar const* direction, qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2);
qaws_status qaws_repulsive_mesh_gradient(qaws_repulsive_mesh const* mesh, qaws_scalar const* positions, qaws_scalar seed,
	qaws_scalar* gradient, qaws_scalar* out_value);
qaws_status qaws_repulsive_mesh_hvp(qaws_repulsive_mesh const* mesh, qaws_scalar const* positions, qaws_scalar const* direction,
	qaws_scalar* out_hv);

/* ---- hierarchical evaluation (R4) --------------------------------------- */

/*
 * Barnes-Hut energy and gradient (Repulsive Curves, Sec. 6.1; Repulsive
 * Surfaces, Sec. 4). Each ordered pair's term uses its first element's own
 * tangent or normal, so for an element I the elements of a tree node far
 * from it are taken at once, by their total length (area) and weighted
 * centroid with one-point quadrature: w_I W k(x_I, X, T_I). A node is far
 * when max(r_I, r_node) <= theta dist(box_I, box_node), r the boxes' half
 * diagonals; a node holding an element sharing a vertex with I is never
 * far. The rest are the direct evaluator's exact terms. theta = 0 is the
 * direct sum; larger is faster and coarser (0.25 to 0.5 typical).
 *
 * The gradient (gradient += seed dE/dx, may be NULL) follows the paper's
 * scheme: each element's vertices get, from its own traversal, the
 * derivative of both orders of each pair in those vertices, exactly for the
 * pairs near it and, for a far node, with the node's averaged projector
 * standing for its members' tangents or normals in the reverse term. It
 * approximates the direct gradient (converging to it as theta falls); it
 * is not the derivative of the approximated energy. No HVP here:
 * hierarchical second derivatives are not provided.
 */
typedef struct qaws_repulsive_bh_stats
{
	unsigned long long near_pairs;         /* exact element pairs */
	unsigned long long far_pairs;          /* element-node pairs taken at once */
	unsigned int nodes;
} qaws_repulsive_bh_stats;

qaws_status qaws_repulsive_curves_bh(qaws_repulsive_curves const* curves, qaws_scalar const* positions, qaws_scalar theta, qaws_scalar seed,
	qaws_scalar* gradient, qaws_scalar* out_value, qaws_repulsive_bh_stats* out_stats);
qaws_status qaws_repulsive_mesh_bh(qaws_repulsive_mesh const* mesh, qaws_scalar const* positions, qaws_scalar theta, qaws_scalar seed,
	qaws_scalar* gradient, qaws_scalar* out_value, qaws_repulsive_bh_stats* out_stats);

/* ---- mesh optimization (R3) --------------------------------------------- */

/*
 * The surfaces' fractional Sobolev descent (Repulsive Surfaces, Sec. 3.2
 * and 3.3): the metric A = B + B0 on vertex functions (the high-order term
 * on the per-face gradients of the hat functions, the low-order term on
 * face averages weighted by the p = 2 kernel, exponent 2 sigma + 2,
 * sigma = s - 1, s = 2 - 2 / p), over every ordered pair of distinct
 * faces; then the same constrained direction, step, projection, modes and
 * report as the curves. Connectivity is fixed.
 *
 * Constraints, held at their values at the solver's creation: pinned
 * vertices (an obstacle is a pinned component), the total area, the
 * enclosed volume (a closed mesh only: every edge on exactly two faces),
 * and the barycenter (by default when nothing is pinned).
 *
 * Collision-checked: every accepted motion keeps each vertex and face it
 * does not belong to, and each two edges sharing no vertex, at least
 * min(clearance, half their distance at the step's start) apart. Contact
 * between triangle meshes always begins at such a pair, so no two
 * triangles ever touch beyond their shared vertices and edges.
 */
typedef struct qaws_repulsive_mesh_solve_desc
{
	unsigned int const* pinned;
	unsigned int pinned_count;
	int hold_area;
	int hold_volume;
	int hold_barycenter;                   /* 0: when nothing is pinned; 1: always; -1: never */
	int collision_checked;
	qaws_scalar clearance;                 /* 0: 1e-3 of the mean edge length */
	unsigned int max_iterations;           /* solve: 0 = 100 */
	qaws_scalar tolerance;                 /* stationarity; 0 = 1e-4 */
	qaws_scalar constraint_tolerance;      /* |Phi| over the mean edge length l (area rows over l, volume rows over l^2); 0 = 1e-8 (1e-5 in f32) */
	qaws_scalar theta;                     /* 0: dense; > 0: hierarchical, as for curves */
} qaws_repulsive_mesh_solve_desc;

typedef struct qaws_repulsive_mesh_solver qaws_repulsive_mesh_solver;

qaws_status qaws_repulsive_mesh_solver_create(qaws_repulsive_mesh const* mesh, qaws_repulsive_mesh_solve_desc const* desc,
	qaws_scalar const* positions, qaws_repulsive_mesh_solver** out_solver);
void qaws_repulsive_mesh_solver_destroy(qaws_repulsive_mesh_solver* solver);
qaws_status qaws_repulsive_mesh_solver_step(qaws_repulsive_mesh_solver* solver, qaws_scalar* positions, qaws_repulsive_report* out_report);
qaws_status qaws_repulsive_mesh_solver_direction(qaws_repulsive_mesh_solver const* solver, qaws_scalar const* positions,
	qaws_scalar* out_direction);
qaws_status qaws_repulsive_mesh_solve(qaws_repulsive_mesh const* mesh, qaws_repulsive_mesh_solve_desc const* desc,
	qaws_scalar* positions, qaws_repulsive_report* out_report);

/* ---- derivatives of a converged solution (R5) ---------------------------- */

/*
 * At a solver whose last step reported CONVERGED (dense mode), the
 * derivative of the solution x* with respect to the constraints' targets.
 * It comes from the KKT system of the Lagrangian E - lambda . c, not from
 * the descent's preconditioner:
 *   [K C^T; C 0] [dx; mu] = [0; dtargets],  K = Hess E - sum lambda_r Hess c_r,
 * K applied exactly (the energy's HVP and the constraints' Hessians),
 * solved by GMRES preconditioned by the descent's Sobolev system. JVP:
 * targets' changes to dx. VJP: xbar to the targets' adjoints (the same
 * symmetric system).
 *
 * Contracts:
 *   - NOT_CONVERGED unless the solver's last reason is CONVERGED: the
 *     identity holds at a stationary point only (as close as the
 *     tolerance; tighten it for accurate derivatives);
 *   - rotations about the barycenter that every constraint allows (nothing
 *     pinned, no held direction) leave the solution non-isolated. They are
 *     gauged out (dx orthogonal to them), and *out_gauge says how many (0
 *     to 3);
 *   - NUMERICAL_FAILURE for dependent constraint rows, or a KKT system
 *     GMRES cannot solve (not a strict minimum on the constraint set);
 *   - UNSUPPORTED_OPERATION in the hierarchical mode;
 *   - not provided: the derivative of a fixed sequence of steps, and
 *     derivatives in the energy's exponents.
 * Held directions have no target here (their rows' changes are 0).
 */
typedef struct qaws_repulsive_targets
{
	qaws_scalar* pinned;                   /* the pinned positions (dimension each, the desc's order); NULL: none */
	qaws_scalar barycenter[3];             /* the held barycenter */
	qaws_scalar length;                    /* curves: the total length */
	qaws_scalar* edge_lengths;             /* curves: each edge's (hold_edge_lengths), NULL: none */
	qaws_scalar area, volume;              /* meshes */
} qaws_repulsive_targets;

qaws_status qaws_repulsive_curves_solver_solution_jvp(qaws_repulsive_curves_solver const* solver, qaws_repulsive_targets const* dtargets,
	qaws_scalar* out_dx, unsigned int* out_gauge);
qaws_status qaws_repulsive_curves_solver_solution_vjp(qaws_repulsive_curves_solver const* solver, qaws_scalar const* xbar,
	qaws_repulsive_targets* out_targets_bar, unsigned int* out_gauge);
qaws_status qaws_repulsive_mesh_solver_solution_jvp(qaws_repulsive_mesh_solver const* solver, qaws_repulsive_targets const* dtargets,
	qaws_scalar* out_dx, unsigned int* out_gauge);
qaws_status qaws_repulsive_mesh_solver_solution_vjp(qaws_repulsive_mesh_solver const* solver, qaws_scalar const* xbar,
	qaws_repulsive_targets* out_targets_bar, unsigned int* out_gauge);

/* ---- curves on a surface (R5) ------------------------------------------- */

/*
 * A 3D curve set whose vertices live on a parametric surface S: vertex i at
 * x_i = S(u_i, v_i), its parameters the unknowns (uv: two per vertex). The
 * energy is the curves' own, measured in the embedding (no UV or geodesic
 * distance enters it), and its derivatives in uv carry S's first and
 * second derivatives (the chain rule, exact).
 *
 * The descent (dense mode only) pulls the metric back through S, the
 * blocks A_ij J_i^T J_j, J_i = [S_u S_v] at vertex i; a 1e-10 multiple of
 * the diagonal blocks lifts what A cannot see (a curve sliding along a
 * cylinder's axis). Steps move the parameters, so the curve stays on S
 * exactly; no projection onto the surface is needed. Pins hold a vertex's
 * parameters; lengths and held directions are the curves' (measured in 3D);
 * no barycenter is held.
 *
 * The domain is explicit:
 *   - a parameter declared periodic wraps (a seam is crossed freely; the
 *     uv returned are unwrapped, continuous across it); it must be periodic
 *     in S, which is not checked;
 *   - a parameter not declared periodic must stay in its range;
 *   - a trimmed surface's vertices must stay in the trimmed region;
 *   - a vertex at a singular point (S_u x S_v about 0, a pole) is refused.
 * A step that would leave the domain is rejected like one that raises the
 * energy (the line search backtracks); an input outside it is refused.
 * Only the vertices are kept in the domain: an edge's chord is straight in
 * 3D and may cross a trim's hole. One surface, one chart: curves crossing
 * between patches of a multi-patch model are not supported.
 *
 * Collision-checked: every accepted motion is checked as the straight
 * motion of each vertex in 3D from x to x' (a chord of its path on S), with
 * the curves' rule.
 */
typedef struct qaws_repulsive_surface_desc
{
	qaws_surface const* surface;           /* borrowed */
	int periodic_u, periodic_v;            /* declared seams */
} qaws_repulsive_surface_desc;

/* the positions S(uv) of count vertices (3 each) */
qaws_status qaws_repulsive_surface_positions(qaws_repulsive_surface_desc const* surface, qaws_scalar const* uv, unsigned int count,
	qaws_scalar* out_positions);

/* the energy E(S(uv)) and, along duv, its first and second directional
   derivatives; gradient += seed dE/duv; out_hv += H duv (H the Hessian in
   uv, S's second derivatives included) */
qaws_status qaws_repulsive_curves_on_surface_eval(qaws_repulsive_curves const* curves, qaws_repulsive_surface_desc const* surface,
	qaws_scalar const* uv, qaws_scalar const* duv, qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2);
qaws_status qaws_repulsive_curves_on_surface_gradient(qaws_repulsive_curves const* curves, qaws_repulsive_surface_desc const* surface,
	qaws_scalar const* uv, qaws_scalar seed, qaws_scalar* gradient, qaws_scalar* out_value);
qaws_status qaws_repulsive_curves_on_surface_hvp(qaws_repulsive_curves const* curves, qaws_repulsive_surface_desc const* surface,
	qaws_scalar const* uv, qaws_scalar const* duv, qaws_scalar* out_hv);

typedef struct qaws_repulsive_surface_curves_solver qaws_repulsive_surface_curves_solver;

/* desc as for curves; hold_barycenter and theta must be 0 */
qaws_status qaws_repulsive_surface_curves_solver_create(qaws_repulsive_curves const* curves, qaws_repulsive_surface_desc const* surface,
	qaws_repulsive_curves_solve_desc const* desc, qaws_scalar const* uv, qaws_repulsive_surface_curves_solver** out_solver);
void qaws_repulsive_surface_curves_solver_destroy(qaws_repulsive_surface_curves_solver* solver);
qaws_status qaws_repulsive_surface_curves_solver_step(qaws_repulsive_surface_curves_solver* solver, qaws_scalar* uv,
	qaws_repulsive_report* out_report);
qaws_status qaws_repulsive_surface_curves_solver_direction(qaws_repulsive_surface_curves_solver const* solver, qaws_scalar const* uv,
	qaws_scalar* out_direction);
qaws_status qaws_repulsive_surface_curves_solve(qaws_repulsive_curves const* curves, qaws_repulsive_surface_desc const* surface,
	qaws_repulsive_curves_solve_desc const* desc, qaws_scalar* uv, qaws_repulsive_report* out_report);

/* ---- native control points (R5) ---------------------------------------- */

/*
 * The working polyline of qaws_repulsive_polyline_from_curves as a linear
 * function of the curves' control points: x_s = sum_k B_sk P_k, the same
 * weights for every coordinate. The samples (curve, parameter) are chosen
 * once, at creation, and held: the derivatives are those of that fixed
 * sampling. B comes from the curves' own tangent evaluation (qaws_diff.h)
 * and is checked: a control's tangent in one coordinate moves that
 * coordinate only, and B P reproduces the samples. A curve that fails the
 * checks (one not linear in its control points with its other fields held)
 * is refused with UNSUPPORTED_OPERATION. Bezier, B-spline and NURBS curves
 * (weights held) qualify.
 *
 * Energy derivatives in P are B^T those in x: exact, the map being linear
 * (the HVP is B^T H B v). The solver descends in P with the metric B^T A B
 * per coordinate. With theta > 0 the energy is the Barnes-Hut one and
 * B^T A B is assembled column by column from the hierarchical metric, so
 * nothing of size |V|^2 is stored, only |P|^2. Pins name control points.
 * The barycenter is the samples'. Lengths and held directions are the
 * polyline's. Collision-checked: the samples move straight when the
 * controls do, and that motion is checked. The curve between the samples
 * is not checked.
 */
typedef struct qaws_repulsive_controls_desc
{
	qaws_repulsive_polyline_desc polyline; /* the curves and their sampling */
	qaws_scalar alpha, beta;               /* 0: (3, 6) */
} qaws_repulsive_controls_desc;

typedef struct qaws_repulsive_controls qaws_repulsive_controls;

qaws_status qaws_repulsive_controls_create(qaws_repulsive_controls_desc const* desc, qaws_repulsive_controls** out_controls);
void qaws_repulsive_controls_destroy(qaws_repulsive_controls* controls);
/* the control points (all curves', in order, dimension each), the samples, the edges */
void qaws_repulsive_controls_get_counts(qaws_repulsive_controls const* controls, unsigned int* out_control_count, unsigned int* out_vertex_count,
	unsigned int* out_edge_count);
/* the energy over the polyline (for evaluating positions directly) */
qaws_repulsive_curves const* qaws_repulsive_controls_get_curves(qaws_repulsive_controls const* controls);
/* the control points as read from the curves at creation */
qaws_status qaws_repulsive_controls_read(qaws_repulsive_controls const* controls, qaws_scalar* out_control_points);
/* the samples B P */
qaws_status qaws_repulsive_controls_positions(qaws_repulsive_controls const* controls, qaws_scalar const* control_points, qaws_scalar* out_positions);
/* as the curves' eval / gradient / hvp, in the control points */
qaws_status qaws_repulsive_controls_eval(qaws_repulsive_controls const* controls, qaws_scalar const* control_points, qaws_scalar const* direction,
	qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2);
qaws_status qaws_repulsive_controls_gradient(qaws_repulsive_controls const* controls, qaws_scalar const* control_points, qaws_scalar seed,
	qaws_scalar* gradient, qaws_scalar* out_value);
qaws_status qaws_repulsive_controls_hvp(qaws_repulsive_controls const* controls, qaws_scalar const* control_points, qaws_scalar const* direction,
	qaws_scalar* out_hv);
/* curve `curve` with these control points: the same family, its other
   fields unchanged (qaws_curve_clone_with_fields) */
qaws_status qaws_repulsive_controls_to_curve(qaws_repulsive_controls const* controls, qaws_scalar const* control_points, unsigned int curve,
	qaws_curve** out_curve);

typedef struct qaws_repulsive_controls_solver qaws_repulsive_controls_solver;

qaws_status qaws_repulsive_controls_solver_create(qaws_repulsive_controls const* controls, qaws_repulsive_curves_solve_desc const* desc,
	qaws_scalar const* control_points, qaws_repulsive_controls_solver** out_solver);
void qaws_repulsive_controls_solver_destroy(qaws_repulsive_controls_solver* solver);
qaws_status qaws_repulsive_controls_solver_step(qaws_repulsive_controls_solver* solver, qaws_scalar* control_points, qaws_repulsive_report* out_report);
qaws_status qaws_repulsive_controls_solver_direction(qaws_repulsive_controls_solver const* solver, qaws_scalar const* control_points,
	qaws_scalar* out_direction);
qaws_status qaws_repulsive_controls_solve(qaws_repulsive_controls const* controls, qaws_repulsive_curves_solve_desc const* desc,
	qaws_scalar* control_points, qaws_repulsive_report* out_report);

#endif /* QAWS_REPULSIVE_H */
