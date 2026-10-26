#ifndef QAWS_MESH_INTERNAL_H
#define QAWS_MESH_INTERNAL_H

/*
 * The meshing module's internals shared between its files: the
 * interpolants' patches and the mesh reference's layout.
 */

#include "qaws_mesh_core.h"
#include "../qaws_mesh_reference.h"
#include "../qaws_mesh_field.h"
#include "../qaws_volume_mesh.h"
#include "../internal/qaws_internal_broadphase.h"
#include "../internal/qaws_internal_ldlt.h"
#include "../qaws_mesh_boolean.h"

/* qaws_mesh_boolean; one: edges inside flat regions turned only between pieces of one source face; out_turned (may be
   NULL, freed by the caller): per result face, 1 where it was turned over (qaws_internal_free; NULL on failure) */
qaws_status sx_boolean(qaws_mesh_3d** out, qaws_boolean_report* report, qaws_mesh_3d_desc const* a, qaws_mesh_3d_desc const* b,
	qaws_boolean_op op, qaws_intersect_desc const* desc, int one, unsigned char** out_turned);
/* edges (2 vertices each) keyed and sorted with their indices (2 per edge; qaws_internal_free), and edge (a, b)'s index
   among them (QAWS_MESH_NONE: none) */
unsigned long long* cy_keys(unsigned int const* e, unsigned int ne);
unsigned int cy_find(unsigned long long const* k, unsigned int ne, unsigned int a, unsigned int b);
/* qaws_mesh_3d_boolean of soups of meshes: a's (na of them) as one operand, b's (nb; 0: the self-union) as the other */
qaws_status cy_boolean(qaws_mesh_3d** out, qaws_boolean_report* report, qaws_mesh_3d const* const* a, unsigned int na,
	qaws_mesh_3d const* const* b, unsigned int nb, qaws_boolean_op op, qaws_intersect_desc const* desc);

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
/* a times factor (a must outlive it) */
qaws_status sf_create_scaled(qaws_size_field const* a, double factor, qaws_size_field** out);

/* ---- metric fields ------------------------------------------------------ */

struct qaws_metric_field
{
	qaws_mesh_reference const* ref;
	double* m;                      /* six per vertex of the reference's mesh: xx, xy, xz, yy, yz, zz */
};

/* the metric at l (six values) */
void mf_eval(qaws_metric_field const* f, qaws_mesh_location const* l, double* M);
/* a symmetric 2x2 g (g00, g01, g11) in the tangent basis (t1, t2) of normal n, its eigenvalues made positive and clamped to
   [lo, hi], as a 3x3 metric M (six values) that measures the normal by the smaller one; returns their product */
double mf_lift(double* M, double const* g, double const* t1, double const* t2, double const* n, double lo, double hi);
/* a tangent basis of unit normal n */
void mf_basis(double const* n, double* t1, double* t2);

/* the reference mesh's attributes onto m, at its vertices' footpoints and its faces' centre locations fl (remeshing, quads) */
qaws_status mx_transfer(qaws_hm* m, qaws_mesh_reference const* r, qaws_mesh_location const* fl);

/* qaws_mesh_repair.c: triangles P and Q (vertex ids ip, iq) intersect beyond their shared corners and edge, by exact
   predicates */
int rp_tri_tri_at(double const (*P)[3], double const (*Q)[3], unsigned int const* ip, unsigned int const* iq);

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
	double* sc;                     /* the lattice's local spacing over the scale (NULL: 1 everywhere; mf_scales) */
	double* wt;                     /* per vertex its weight in the orientations' smoothing, an edge the mean of its ends' (NULL: 1; mf_weights) */
} mf_level;

#define MF_ALONG 1u
#define MF_ON_LINE 2u
#define MF_PINNED 4u
#define MF_PREFER 8u                    /* q drawn toward cq, as much as |cq| times its neighbours together */
#define MF_GUIDE 16u                    /* on a user's guide line (with MF_ALONG and MF_ON_LINE), not a crease */
#define MF_CONE 32u                     /* a placed cone's vertex: the steps round its star kept by the flow (qaws_quad_flow.c) */

typedef struct mf_hier
{
	mf_level* lv;
	unsigned int nl;
} mf_hier;

/* level 0: the vertices of a triangle store, its edges, normals (NULL: angle-weighted), constraints from its features */
qaws_status mf_build(mf_hier* h, qaws_hm const* m, double const* normals);
void mf_free(mf_hier* h);
/* init (3 per level-0 vertex, zero: none) may be NULL: a preferred direction carried up the levels, the start; with pull also
   a pull on every level as strong as its length (1: as much as the neighbours together) */
void mf_orientations(mf_hier* h, unsigned int sweeps, double const* init, int pull);
void mf_positions(mf_hier* h, double scale, unsigned int sweeps);
/* user guidance on level 0, carried up the levels: per vertex a direction (3) and a weight w (> 0: MF_PREFER, cq the
   direction made tangent times w; < 0 or infinite: MF_ALONG; 0: none); vertices already MF_ALONG keep theirs */
void mf_guide(mf_hier* h, double const* dir, double const* w);
/* per level-0 vertex the lattice's spacing over the scale (all positive), carried up area-weighted; 0 on allocation failure */
int mf_scales(mf_hier* h, double const* s0);
/* the spacing between vertices i and j of L at the scale: scale times their mean sc */
double mf_pair_scale(mf_level const* L, unsigned int i, unsigned int j, double scale);
/* per level-0 vertex a weight for the orientations' smoothing (all positive), carried up area-weighted; 0 on allocation failure */
int mf_weights(mf_hier* h, double const* w0);
/* singularity attractors, points p_a of strength s_a (> 0 draws singularities, < 0 pushes them off) reaching r, at the n
   points x of areas area: out[0, n), LOCAL's smoothing weight exp(-sum s_a exp(-|x - p_a|^2 / r^2)) (singularities cost
   less where it is low); out[n, 2 n), SMOOTHEST's share of cones (s_a spread by area times that Gaussian, cf_steer) */
void mf_attract(double const* x, double const* area, unsigned int n, qaws_scalar const* points, qaws_scalar const* strengths, unsigned int count,
	double r, double* out);
/* the 4-symmetric copies of q0 and q1 closest to each other */
void mf_compat_orientation(double const* q0, double const* n0, double const* q1, double const* n1, double* out0, double* out1);
/* the lattice points of (q0, o0) and (q1, o1) closest to each other near x0 and x1 */
void mf_compat_position(double const* x0, double const* n0, double const* q0, double const* o0, double const* x1, double const* n1,
	double const* q1, double const* o1, double scale, double* out0, double* out1);
/* the quarter turns taking q1 onto q0 (0 to 3) */
int mf_rotation(double const* q0, double const* n0, double const* q1, double const* n1);
/* p, a point of x1's tangent plane, turned into x0's about the point near both (the rotation taking n1 onto n0): unfolded across a crease */
void mf_unfold(double const* x0, double const* n0, double const* x1, double const* n1, double const* p, double* out);

/* ---- volume meshes ------------------------------------------------------ */

struct qaws_volume_mesh
{
	double* x;                      /* 3 per vertex */
	unsigned int nv;
	unsigned char* type;            /* per cell, a qaws_cell_type */
	unsigned int *c0, *cv;          /* per cell its first corner (nc + 1 of them); per corner its vertex */
	unsigned int *f0, *fo;          /* per cell its first facet (nc + 1); per facet the facet across, HM_NONE on the border */
	unsigned int nc;
};

/* corners and facets per cell type; facet f's corners, outward (*n of them) */
extern unsigned char const vm_corners[4], vm_facets[4];
unsigned char const* vm_facet(unsigned int type, unsigned int f, unsigned int* n);
/* a mesh from its arrays (the facets matched here); INVALID_ARGUMENT for a bad index or type, INVALID_TOPOLOGY for unmatched facets */
qaws_status vm_create(double const* x, unsigned int nv, unsigned int const* cv, unsigned char const* type, unsigned int nc, qaws_volume_mesh** out);
/* the border points the tetrahedra's recovery added (steiner, per vertex) moved inside where they can be, their border
   triangles merged back: out_moved of them, out_left not (qaws_volume_improve.c) */
qaws_status vim_unsplit_border(qaws_volume_mesh** out_mesh, unsigned int* out_moved, unsigned int* out_left, qaws_volume_mesh const* mesh,
	unsigned char const* steiner);
/* a hexahedron's smallest scaled Jacobian over its corners (1 for a cube, negative when a corner turns over) */
double vm_hex_quality(double const* x, unsigned int const* v);

/* frame fields on tetrahedra (qaws_frame_field.c): per cell its axes (9 doubles, R allocated by the caller), the
   singular edges' masks (may be NULL); mbo: 0, or MBO rounds of that many diffusion steps after the sweeps (the singular
   curves parted to about sqrt(mbo) cells) */
qaws_status ff_frames(qaws_volume_mesh const* m, unsigned int sweeps, unsigned int mbo, double* R, unsigned char* singular);
/* the frames combed into one gauge (each cell's axes permuted onto its neighbour's, breadth first); 1 when every facet
   then agrees without a permutation (no singular edge, no turn round a hole), 0 else or on allocation failure */
int ff_comb(qaws_volume_mesh const* m, double* R);
/* hexahedra of tetrahedra t's frames R (combed in place) with singular edges: an integer-grid map with transitions across a
   cut (CubeCover), its lattice read off (as HexEx) (qaws_hex_seamless.c); UNSUPPORTED_OPERATION for a field no hexahedra
   follow (an edge turning by other than a quarter turn, a cell turned over in the map) */
qaws_status hs_hexes(qaws_volume_mesh const* t, double* R, double h, qaws_volume_mesh** out_mesh, qaws_mesh_report* rep);

/* ---- sparse Laplacians (qaws_mesh_smooth.c), shared by smoothing, fairing, deformation and parameterization ---- */

/* a sparse matrix by rows, columns sorted in each */
typedef struct sm_csr
{
	unsigned int n;
	unsigned int *rp, *ci;
	double* v;
} sm_csr;

/* a symmetric system restricted to the rows that move (the others held: their columns go to the right-hand side), factored */
typedef struct sm_system
{
	sm_csr const* A;
	unsigned int* row;          /* per row of A its row here, 0xFFFFFFFF: held */
	unsigned int m;
	qaws_internal_ldlt L;
} sm_system;

void sm_csr_free(sm_csr* A);
/* K (K_ii the sum of the weights, K_ij minus the weight: half the cotangents facing the edge, clamped at 0 with clamp; 1 with
   uniform) over triangles t (3 each) of positions x (dim 3), and the lumped masses (1 with uniform); 0 on allocation failure */
int sm_laplacian(sm_csr* K, double* mass, double const* x, unsigned int nv, unsigned int const* t, unsigned int nt, int uniform, int clamp);
/* C = A diag(d) B (d NULL: the identity) */
int sm_mul(sm_csr* C, sm_csr const* A, double const* d, sm_csr const* B);
/* the system over the rows where move is nonzero; xyz: 3 per row of A for the ordering (NULL: natural) */
qaws_status sm_system_factor(sm_system* S, sm_csr const* A, double const* xyz, unsigned char const* move);
/* x (dim per row of A) solved on the moving rows from rhs, the held rows' x read; its own workspace */
qaws_status sm_system_solve(sm_system const* S, double* x, double const* rhs, unsigned int dim);
void sm_system_free(sm_system* S);
/* labels by alpha-expansion graph cuts (Boykov, Veksler, Zabih 2001): n sites, L labels (at most 255); data: L costs per site;
   nbr: 3 neighbours per site (0xFFFFFFFF: none; symmetric), wgt: those pairs' weights; a pair's cost its weight times V[a L + b]
   (a metric: V[a L + a] = 0, the triangle inequality). lab: the start in, the labels out; cycles over the labels until none
   lowers the energy (rounds at most). 0 on allocation failure */
int sm_label_cut(unsigned char* lab, double const* data, unsigned int n, unsigned int L, unsigned int const* nbr, double const* wgt,
	double const* V, unsigned int rounds);

/* the symmetric 3 x 3 A (row major) diagonalized in place by Jacobi rotations, its eigenvectors the columns of V (qaws_mesh_watertight.c) */
void wt_sym_eigen3(double* A, double* V);

/* a grid edge the contoured boundary crosses (qaws_mesh_watertight.c, qaws_mesh_offset.c) */
typedef struct wt_cross
{
	size_t p;                  /* the edge's lower grid point */
	unsigned int axis;
	double x[3], nrm[3];       /* the crossing and its normal, from inside to outside */
} wt_cross;

/* a field to contour (wt_contour_field, qaws_mesh_offset.c): phi at n points (3 each), negative inside, and when grad is
   not NULL its gradient (3 each) */
typedef qaws_status (*wt_field_fn)(void* ctx, double const* x, size_t n, double* phi, double* grad);

typedef struct wt_contour_info
{
	unsigned int grid[3];      /* grid points per axis */
	double cell;
	unsigned int crossings, sharp;
	int closed;
	size_t evaluated;          /* grid points the field was evaluated at */
} wt_contour_info;

/* phi = 0 contoured on a grid over lo..hi grown by margin, res cells along the longest side (*out NULL when empty). lip > 0:
   phi changes by at most lip per unit length (a distance: 1), so it is evaluated on a coarse lattice first and exactly only
   at the points of blocks the level may cross, the others given their block's sign: the same contour, fewer evaluations */
qaws_status wt_contour_field(qaws_mesh_3d** out, wt_contour_info* info, wt_field_fn f, void* ctx, double const* lo, double const* hi, double margin,
	unsigned int res, double lip);

/* dual contouring of the crossings X on a grid (lo, cell h, n points per axis; in per point: 1 inside) into a mesh through
   qaws_mesh_repair (non-manifold places split); *out_sharp += its vertices on edges or corners */
qaws_status wt_dual_contour(qaws_mesh_3d** out, unsigned int* out_sharp, int* out_closed, double const* lo, double h, unsigned int const* n,
	unsigned char const* in, wt_cross const* X, size_t nc);

#endif /* QAWS_MESH_INTERNAL_H */
