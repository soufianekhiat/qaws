#ifndef QAWS_MESH_CORE_H
#define QAWS_MESH_CORE_H

/*
 * The corner store of the meshing module, under the 2D and 3D mesh types
 * (geogram's facet corners, the corner table of Rossignac et al.). A face is
 * a run of consecutive corners; each corner holds its vertex and its
 * opposite: the corner of the neighbouring face at the far end of the
 * corner's edge (the edge from its vertex to the next corner's). Faces'
 * adjacency, the dual graph, is stored with the faces.
 *
 * A half-edge is a corner: next and prev within its face's run, the twin is
 * the opposite corner (HM_NONE on the border: border edges have no stored
 * half-edge). A vertex keeps one of its corners: on the border, the one
 * whose own edge is on the border, so turning counter-clockwise from it
 * (the opposite of the previous corner) meets every face round the vertex.
 *
 * Local operations replace a patch of faces by new ones and re-stitch them
 * to the faces round the patch. Removed vertices and faces are holes until
 * hm_compact moves the last ones into them (indices dense again); freed
 * corner slots are garbage, the pool rebuilt when half of it is.
 */

#include "../qaws_mesh.h"

#define HM_NONE 0xFFFFFFFFu

#define HM_DELETED 1u
#define HM_FEATURE 2u               /* corners: the corner's edge is a feature (on both its corners) */
#define HM_CORNER 4u                /* vertices: a feature corner */
#define HM_MARK 8u                  /* faces: in the patch being replaced */

typedef struct hm_attr
{
	char name[32];
	unsigned int element;           /* qaws_mesh_element */
	unsigned int rule;              /* qaws_attribute_rule */
	unsigned int n;                 /* components, at most 16 */
	double* v;                      /* n per slot of its element (corner slots for corners) */
} hm_attr;

/* the patch an operation builds: old faces, new faces of corners with
   attribute sources (positions in the old faces' corners, weighted) */
typedef struct hm_patch
{
	unsigned int *old, *off;        /* old faces; their first corner position */
	unsigned int nold, cold, total;
	unsigned int *fs, *fa;          /* new faces' sizes; the old face (index in old) they take attributes from */
	unsigned int nnew, cnew;
	unsigned int *fv, *s0, *sn;     /* new corners' vertices; their sources */
	unsigned char* hint;            /* flags or'ed into new corners */
	unsigned int ncv, ccv;
	unsigned int* sp;               /* sources: old corner positions */
	double* sw;                     /* and weights */
	unsigned int nsrc, csrc;
	int failed;
} hm_patch;

typedef struct qaws_hm
{
	unsigned int dim;
	double* x;                      /* dim per vertex */
	unsigned int* vc;               /* a corner of the vertex; HM_NONE: isolated */
	unsigned int* vsrc;             /* the given vertex a copy comes from */
	unsigned char* vflag;
	unsigned int nv, capv;
	unsigned int* cv;               /* corner: vertex */
	unsigned int* co;               /* opposite corner; HM_NONE: border */
	unsigned int* cf;               /* face; HM_NONE: free slot */
	unsigned char* cflag;
	unsigned int nc, capc, cfree;   /* slots used, capacity, free slots below nc */
	unsigned int* fc;               /* face: first corner */
	unsigned int* fn;               /* corners */
	unsigned char* fflag;
	unsigned int nf, capf;
	hm_attr* attr;
	unsigned int na;
	qaws_mesh_location* foot;       /* per vertex; NULL: none */
	unsigned int* hole[2];          /* removed vertices, faces */
	unsigned int nhole[2];
	hm_patch pb;                    /* scratch of the operations */
} qaws_hm;

struct qaws_mesh_3d
{
	qaws_hm hm;
};

/* ---- building and checking ------------------------------------------- */

qaws_status hm_build(qaws_hm* m, unsigned int dim, double const* x, unsigned int nv, unsigned int const* fv,
	unsigned int const* fs, unsigned int nf);
void hm_free(qaws_hm* m);
qaws_status hm_copy(qaws_hm* dst, qaws_hm const* src);
/* room for dv more vertices, dc corners, df faces; 0 on allocation failure (nothing changed) */
int hm_reserve(qaws_hm* m, unsigned int dv, unsigned int dc, unsigned int df);
qaws_status hm_check(qaws_hm const* m);

/* ---- navigation ------------------------------------------------------- */

unsigned int hm_next(qaws_hm const* m, unsigned int c);
unsigned int hm_prev(qaws_hm const* m, unsigned int c);
unsigned int hm_to(qaws_hm const* m, unsigned int c);      /* the vertex at the far end of c's edge */
unsigned int hm_ccw(qaws_hm const* m, unsigned int c);     /* the next corner of c's vertex counter-clockwise; HM_NONE at the border */
unsigned int hm_cw(qaws_hm const* m, unsigned int c);
/* the corner after c round vertex v counter-clockwise, HM_NONE past the last (from vc[v]) */
unsigned int hm_round(qaws_hm const* m, unsigned int v, unsigned int c);
unsigned int hm_find(qaws_hm const* m, unsigned int a, unsigned int b);        /* the corner of half-edge a -> b, or HM_NONE */
unsigned int hm_find_edge(qaws_hm const* m, unsigned int a, unsigned int b);   /* a corner of edge {a, b}, either way */
int hm_is_border_vertex(qaws_hm const* m, unsigned int v);
unsigned int hm_valence(qaws_hm const* m, unsigned int v);                     /* edges */
/* v's neighbours counter-clockwise; returns their count (out written up to cap) */
unsigned int hm_ring(qaws_hm const* m, unsigned int v, unsigned int* out, unsigned int cap);
void hm_fix_vertex(qaws_hm* m, unsigned int v, unsigned int c);                /* vc[v] from c, one of its corners */
double hm_length(qaws_hm const* m, unsigned int c);
void hm_face_normal(qaws_hm const* m, unsigned int f, double* out);           /* area-weighted (Newell), 3D */
unsigned int hm_edge_count(qaws_hm const* m);

/* ---- operations: 0 or HM_NONE on allocation failure, nothing changed -- */

/* a vertex at x on c's edge, attributes at t from cv[c] to hm_to(c); with
   triangulate the triangles either side are split in two, other faces take
   the vertex. Returns it. */
unsigned int hm_split_edge(qaws_hm* m, unsigned int c, double t, double const* x, int triangulate);
/* a vertex at x in face f (NULL: its centroid, weights exact) joined to its corners; the first triangle keeps f */
unsigned int hm_split_face(qaws_hm* m, unsigned int f, double const* x);
/* face cf[c0] cut by a new edge from cv[c0] to cv[c1] (corners of one face, not neighbours) */
int hm_cut_face(qaws_hm* m, unsigned int c0, unsigned int c1);
/* a, one end of c's edge, merged into the other end b, which moves to x,
   its attributes at t from a's (t = 1: b's) */
int hm_can_collapse(qaws_hm const* m, unsigned int c, unsigned int a);
int hm_collapse(qaws_hm* m, unsigned int c, unsigned int a, double t, double const* x);
int hm_can_flip(qaws_hm const* m, unsigned int c);
int hm_flip(qaws_hm* m, unsigned int c);
void hm_compact(qaws_hm* m);

/* ---- attributes and footpoints --------------------------------------- */

/* index, or HM_NONE on allocation failure; values zeroed */
unsigned int hm_attr_add(qaws_hm* m, char const* name, unsigned int element, unsigned int rule, unsigned int n);
void hm_attr_remove(qaws_hm* m, unsigned int a);
double* hm_attr_at(qaws_hm const* m, unsigned int a, unsigned int i);
void hm_lerp(double* dst, double const* a, double const* b, double t, unsigned int n, unsigned int rule);
int hm_foot_enable(qaws_hm* m);
unsigned int hm_new_vertex(qaws_hm* m, double const* x);   /* capacity reserved by the caller */
void hm_kill_vertex(qaws_hm* m, unsigned int v);
unsigned int hm_feature_degree(qaws_hm const* m, unsigned int v);   /* feature edges at v */

/* a 3D mesh from double arrays (the operators' outputs) */
qaws_status hm_mesh_3d_create(double const* x, unsigned int nv, unsigned int const* fv, unsigned int const* fs, unsigned int nf,
	qaws_mesh_3d** out_mesh);

#endif /* QAWS_MESH_CORE_H */
