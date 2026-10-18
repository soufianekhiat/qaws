#ifndef QAWS_QUAD_GLOBAL_H
#define QAWS_QUAD_GLOBAL_H

/*
 * GLOBAL quads: an integer-grid map (Bommes, Zimmer, Kobbelt 2009; Bommes
 * et al. 2013) of INSTANT's dense mesh, read off as QEx does (Ebke,
 * Bommes, Chen, Kobbelt 2013). Per vertex two unknowns, its place in its
 * lattice cell in its own frame; every edge's transition a quarter turn and
 * a whole shift, so no cut graph is needed. The integers: INSTANT's steps
 * made consistent by the min-cost flow (qaws_quad_flow.h) first; when they
 * give no map, the map's own: the cones (singular faces, corners all the
 * way round) moved by whole steps along the paths of a spanning tree
 * between them and to the border, the map's energy over those steps a
 * small quadratic (one sparse solve per step on one factor), rounded
 * greedily as MIQ does and refined by whole steps; each cone then onto its
 * nearest free vertex. The places minimize sum A |grad u - X / rho|^2 +
 * |grad v - Y / rho|^2 over the faces (X, Y: the field in the face) in one
 * sparse LDL^T solve, vertices on feature lines and the border held on
 * their whole isoline, corners and cones on a lattice point (a corner's
 * place carried round its triangles from the fixed point of the turn round
 * it); folds left are undone vertex by vertex (kernel centres), a vertex
 * wound twice round is no map. A singular face whose cone sits on a
 * vertex's lattice point makes that vertex the cone; else it is split
 * round its cone (three triangles). Extraction: the lattice points in each
 * triangle are the vertices; from each, its four (at a cone three or five)
 * unit isolines traced across the triangles through the transitions to the
 * next point, each nudged a little to its left so it never runs along an
 * edge or through a vertex; faces follow the isolines, turning left by a
 * quarter at each point.
 */

#include "qaws_mesh_internal.h"
#include "../qaws_quadrangulate.h"

/* A quad mesh read off: vertices x (3 each) with footpoints, kind (2: pinned, 1: on a line), faces. */
typedef struct qd_out
{
	double* x;
	qaws_mesh_location* foot;
	unsigned char* kind;
	unsigned int nv;
	unsigned int *fv, *fs, nf, nc;
} qd_out;

void qd_out_free(qd_out* o);

/* the map as a triangle soup: per triangle its corners' positions (9) and uv (6, in lattice steps) */
typedef struct qg_map
{
	unsigned int nt;
	double *x, *uv;
} qg_map;

void qg_map_free(qg_map* m);

/* The map of dm (level 0 of the fields, L, at lattice spacing rho) from its consistent steps T, K and the quads read off
   it into o: positions on dm's faces (not yet on the reference), as footpoint the nearest dense vertex's (a hint; exact
   at a corner), kinds. map (may be NULL): the map kept. INVALID_TOPOLOGY, o empty, when neither those steps nor the
   quantized ones give a map (a face whose steps do not close, a cone off the lattice, a fold left, a vertex wound twice
   round, isolines that do not close into quads). */
qaws_status qg_quads(qaws_hm const* dm, mf_level const* L, double rho, int const* T, unsigned char const* K, qd_out* o, qg_map* map);

/* qaws_mesh_quadrangulate's GLOBAL map of mesh (desc->method GLOBAL), for figures and tests; INVALID_TOPOLOGY when it
   fell back on INSTANT (no map) */
qaws_status qg_parametrization(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qg_map* out);

#endif /* QAWS_QUAD_GLOBAL_H */
