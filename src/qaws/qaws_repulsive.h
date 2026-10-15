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
