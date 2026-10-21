#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_REPULSIVE_MESH_CURVES_H_F32
#define QAWS_REPULSIVE_MESH_CURVES_H_F32
#define QAWS_REPULSIVE_MESH_CURVES_H_BODY
#endif
#else
#ifndef QAWS_REPULSIVE_MESH_CURVES_H_F64
#define QAWS_REPULSIVE_MESH_CURVES_H_F64
#define QAWS_REPULSIVE_MESH_CURVES_H_BODY
#endif
#endif
#ifdef QAWS_REPULSIVE_MESH_CURVES_H_BODY
#undef QAWS_REPULSIVE_MESH_CURVES_H_BODY
#define QAWS_REPULSIVE_MESH_CURVES_H

#include "qaws_repulsive.h"
#include "qaws_mesh.h"

/*
 * Repulsive curves on a triangle mesh, their distances measured along it.
 *
 * A 3D curve set (qaws_repulsive_curves, dimension 3) whose vertices are
 * locations on a mesh's triangles (polygons fanned from their first vertex,
 * as qaws_mesh_3d_get_triangle). The energy is the curves' (Eq. 17 of
 * qaws_repulsive.h) with each pair's distance r taken along the surface:
 *
 *   |T_I x (x_i - x_j)|^alpha / |x_i - x_j|^beta  becomes
 *   (r^2 - (T_I . (x_i - x_j))^2)^(alpha / 2) / r^beta,
 *
 * r the larger of the chord |x_i - x_j| and the geodesic distance, by the
 * heat method (qaws_mesh_geodesic.h) from each vertex's location, the two
 * ways averaged. The geodesic is never shorter than the chord, so where the
 * heat method falls short (near its source) the chord is kept; where the
 * surface is flat the energy is the curves' own. Curves on two sheets close
 * in space but far apart along the surface (the two sides of a thin shell,
 * the folds of a sheet) hardly repel each other. With euclidean set r is
 * the chord: the curves' energy, the vertices kept on the mesh.
 *
 * The gradient: exact for the chord's terms and the edges' tangents and
 * lengths; a geodesic distance's gradient at x_i is the gradient of the
 * distance field from x_j in x_i's triangle (the field from x_i moving
 * with x_i taken by symmetry: the same), so it carries the heat method's
 * errors.
 *
 * Moves: each step is taken in each vertex's tangent plane (its triangle's)
 * and the vertex walked there along the triangles, straight across each
 * edge (qaws_mesh_reference_walk): curves cross faces and stay on the mesh
 * without projecting. A walk stops at a border. The descent is the curves'
 * dense one (the fractional Sobolev metric of the chords, pulled back to
 * the tangent planes), with their pins, lengths and held directions
 * (measured in 3D); hold_barycenter and theta must be 0. Collision-checked:
 * each vertex's motion checked as the chord from x to x'.
 *
 * Cost: per energy, one geodesic distance per vertex (two sparse solves on
 * the mesh each, prefactored) and the all-pairs sum; the n x n distances
 * are stored.
 *
 * Patches and shells reach this through a triangulation of them.
 */

typedef struct qaws_repulsive_mesh_curves_desc
{
	qaws_mesh_3d const* mesh;           /* read at creation (copied) */
	int euclidean;                      /* the chord's length only (no geodesic distances) */
	qaws_scalar geodesic_time;          /* the heat method's time factor (qaws_mesh_geodesic_desc.time); 0: 1 */
} qaws_repulsive_mesh_curves_desc;

typedef struct qaws_repulsive_mesh_curves_solver qaws_repulsive_mesh_curves_solver;

/* locations: one per curve vertex. INVALID_ARGUMENT for a location off the mesh's triangles, a theta or barycenter set */
qaws_status qaws_repulsive_mesh_curves_solver_create(qaws_repulsive_mesh_curves_solver** out_solver, qaws_repulsive_curves const* curves,
	qaws_mesh_location const* locations, qaws_repulsive_mesh_curves_desc const* mesh, qaws_repulsive_curves_solve_desc const* desc);
qaws_status qaws_repulsive_mesh_curves_solver_create_ex(qaws_repulsive_mesh_curves_solver** out_solver, qaws_repulsive_curves const* curves,
	qaws_mesh_location const* locations, qaws_repulsive_mesh_curves_desc const* mesh, qaws_repulsive_curves_solve_desc const* desc,
	qaws_allocator const* allocator);
void qaws_repulsive_mesh_curves_solver_destroy(qaws_repulsive_mesh_curves_solver* solver);

/* one step from locations (the solver's own, or the caller's when they differ), written back */
qaws_status qaws_repulsive_mesh_curves_solver_step(qaws_repulsive_report* out_report, qaws_mesh_location* locations,
	qaws_repulsive_mesh_curves_solver* solver);

/* the energy at the solver's locations and (may be NULL) its gradient in 3D (3 per vertex, in the vertices' triangles' planes) */
qaws_status qaws_repulsive_mesh_curves_solver_energy(qaws_scalar* out_energy, qaws_scalar* out_gradient,
	qaws_repulsive_mesh_curves_solver const* solver);

/* the positions (3 per vertex) of the solver's locations */
qaws_status qaws_repulsive_mesh_curves_solver_get_positions(qaws_scalar* out_positions, qaws_repulsive_mesh_curves_solver const* solver);

/* steps until a termination reason or max_iterations; on_step is given the positions (3 per vertex) */
qaws_status qaws_repulsive_mesh_curves_solve(qaws_repulsive_report* out_report, qaws_mesh_location* locations, qaws_repulsive_curves const* curves,
	qaws_repulsive_mesh_curves_desc const* mesh, qaws_repulsive_curves_solve_desc const* desc);

#endif /* QAWS_REPULSIVE_MESH_CURVES_H */
