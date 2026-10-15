#ifndef QAWS_REPULSIVE_H
#define QAWS_REPULSIVE_H

#include "qaws_types.h"
#include "qaws_status.h"

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

#endif /* QAWS_REPULSIVE_H */
