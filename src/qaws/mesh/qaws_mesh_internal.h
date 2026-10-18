#ifndef QAWS_MESH_INTERNAL_H
#define QAWS_MESH_INTERNAL_H

/*
 * The meshing module's internals shared between its files: the
 * interpolants' patches and the mesh reference's layout.
 */

#include "qaws_mesh_core.h"
#include "../qaws_mesh_reference.h"
#include "../qaws_mesh_field.h"
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
	unsigned int* treg;             /* per triangle: its region (connected across edges that are not creases) */
	unsigned int nreg;
	mi_tri* patch;
	qaws_bp_box* box;               /* per triangle: its control net's box */
	qaws_bp_grid* grid;
	unsigned int* vloc;             /* per vertex: triangle * 3 + corner */
	double extent;                  /* the mesh's size */
};

/* the closest point of patch t to p: its location (a, b) and squared distance */
double mr_closest_on_patch(qaws_mesh_reference const* r, unsigned int t, double const* p, double* out_a, double* out_b);

/* scratch for repeated projections (stamps per triangle, no clearing) */
typedef struct mr_scratch
{
	unsigned int* tried;
	unsigned int* queued;
	unsigned int* queue;
	unsigned int gen;
} mr_scratch;

int mr_scratch_init(mr_scratch* s, unsigned int nt);
void mr_scratch_free(mr_scratch* s);
/* the closest location to p: reachable from triangle `hint` (HM_NONE: anywhere); its squared distance */
double mr_project(qaws_mesh_reference const* r, double const* p, unsigned int hint, mr_scratch* s, qaws_mesh_location* out);
/* the closest location to p on the creases and the border reachable from `hint`; HUGE_VAL when none */
double mr_project_crease(qaws_mesh_reference const* r, double const* p, unsigned int hint, mr_scratch* s, qaws_mesh_location* out);
/* as mr_project, within a region (the triangles between creases): a hint outside it is ignored */
double mr_project_region(qaws_mesh_reference const* r, double const* p, unsigned int hint, unsigned int region, mr_scratch* s,
	qaws_mesh_location* out);
/* the straightest walk by a tangent vector; with stop_at_creases it ends just short of a crease it meets, and returns 1 */
int mr_walk(qaws_mesh_reference const* r, qaws_mesh_location from, double const* tangent, int stop_at_creases, qaws_mesh_location* out_to);
void mr_eval(qaws_mesh_reference const* r, qaws_mesh_location const* l, double* x, double* n);

/* ---- size fields -------------------------------------------------------- */

double sf_eval(qaws_size_field const* f, qaws_mesh_location const* l, double const* x);

/* the reference mesh's attributes onto m, at its vertices' footpoints and its faces' centre locations fl (remeshing, quads) */
qaws_status mx_transfer(qaws_hm* m, qaws_mesh_reference const* r, qaws_mesh_location const* fl);

/* ---- direction and position fields (Instant Meshes) ---------------------
 *
 * A hierarchy of a mesh's vertices, pairs merged level by level; per vertex
 * a 4-symmetric orientation q (a unit tangent: q and n x q up to sign) and a
 * position o (a lattice point near the vertex: the lattice o + scale (i q +
 * j n x q)). Solved coarse to fine by Gauss-Seidel, constraints from feature
 * edges and the border (the orientation along them, the position on them;
 * corners pin the lattice).
 */

typedef struct mf_level
{
	unsigned int n;
	double *x, *nrm, *area;
	unsigned int *adj0, *adj;       /* neighbours (CSR) */
	unsigned int* parent;           /* at the next level */
	unsigned int *child0, *child;   /* at the previous level (CSR) */
	double *q, *o;
	double* cq;                     /* constraint direction (zero: none) */
	unsigned char* con;             /* MF_ALONG: q along cq; MF_ON_LINE: o on the line x + t cq; MF_PINNED: o at x */
} mf_level;

#define MF_ALONG 1u
#define MF_ON_LINE 2u
#define MF_PINNED 4u

typedef struct mf_hier
{
	mf_level* lv;
	unsigned int nl;
} mf_hier;

/* level 0: the vertices of a triangle store, its edges, normals (NULL: angle-weighted), constraints from its features */
qaws_status mf_build(mf_hier* h, qaws_hm const* m, double const* normals);
void mf_free(mf_hier* h);
/* init (3 per level-0 vertex, zero: none) may be NULL: a preferred start, carried up the levels */
void mf_orientations(mf_hier* h, unsigned int sweeps, double const* init);
void mf_positions(mf_hier* h, double scale, unsigned int sweeps);
/* the 4-symmetric copies of q0 and q1 closest to each other */
void mf_compat_orientation(double const* q0, double const* n0, double const* q1, double const* n1, double* out0, double* out1);
/* the lattice points of (q0, o0) and (q1, o1) closest to each other near x0 and x1 */
void mf_compat_position(double const* x0, double const* n0, double const* q0, double const* o0, double const* x1, double const* n1,
	double const* q1, double const* o1, double scale, double* out0, double* out1);
/* the quarter turns taking q1 onto q0 (0 to 3) */
int mf_rotation(double const* q0, double const* n0, double const* q1, double const* n1);
/* p, a point of x1's tangent plane, turned into x0's about the point near both (the rotation taking n1 onto n0): unfolded across a crease */
void mf_unfold(double const* x0, double const* n0, double const* x1, double const* n1, double const* p, double* out);

#endif /* QAWS_MESH_INTERNAL_H */
