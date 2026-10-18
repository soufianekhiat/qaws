#ifndef QAWS_MESH_INTERNAL_H
#define QAWS_MESH_INTERNAL_H

/*
 * The meshing module's internals shared between its files: the
 * interpolants' patches and the mesh reference's layout.
 */

#include "qaws_mesh_core.h"
#include "../qaws_mesh_reference.h"
#include "../internal/qaws_internal_broadphase.h"

/* ---- interpolants ------------------------------------------------------
 *
 * One triangle's patch as a cubic Bezier triangle (control points indexed
 * by the powers of the barycentrics b0, b1, b2: 300 030 003 210 120 021 012
 * 102 201 111) with a quadratic normal field (200 020 002 110 011 101).
 * PN builds the curved net; LINEAR is the flat triangle degree-elevated, so
 * one evaluation serves both. A location (a, b) is b1 = a, b2 = b.
 */

typedef struct mi_tri
{
	double p[10][3];
	double n[6][3];
} mi_tri;

/* P: the corners; N: their unit normals; crease (per edge i, from corner i to
   i + 1): NULL or the edge's unit tangents at its two ends, 6 doubles (the
   curve leaves each end along them) */
void mi_build(mi_tri* t, qaws_interpolant kind, double const* const* P, double const* const* N, double const* const* crease);

/* position and derivatives in (a, b), the unit normal and its derivatives;
   any output may be NULL */
void mi_eval(mi_tri const* t, double a, double b, double* x, double* xa, double* xb, double* xaa, double* xab, double* xbb,
	double* n, double* na, double* nb);

/* principal curvatures k1 >= k2 and unit directions from the tangents and
   the normal field's derivatives (the shape operator) */
void mi_curvature(double const* xa, double const* xb, double const* n, double const* na, double const* nb, double* k1, double* k2,
	double* d1, double* d2);

/* ---- the mesh reference ------------------------------------------------ */

struct qaws_mesh_reference
{
	qaws_interpolant kind;
	qaws_mesh_3d* mesh;             /* the input, copied */
	unsigned int nt;
	unsigned int* tv;               /* per triangle: 3 vertices */
	unsigned int* tc;               /* per triangle: its 3 corner slots in the copy's store */
	unsigned int* tf;               /* per triangle: its face */
	unsigned int* ta;               /* per triangle edge (3): the neighbour's triangle * 3 + edge, or HM_NONE */
	unsigned char* tcrease;         /* per triangle edge: on a feature */
	mi_tri* patch;
	qaws_bp_box* box;               /* per triangle: its control net's box */
	qaws_bp_grid* grid;
	unsigned int* vloc;             /* per vertex: triangle * 3 + corner */
	double extent;                  /* the mesh's size */
};

/* the closest point of patch t to p: its location (a, b) and squared distance */
double mr_closest_on_patch(qaws_mesh_reference const* r, unsigned int t, double const* p, double* out_a, double* out_b);

#endif /* QAWS_MESH_INTERNAL_H */
